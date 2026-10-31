#!/usr/bin/env python3
"""windowcalc 依赖自举：让 agent 替用户把环境装好，用户无感。

设计目标
--------
任何 agent（dsh / opencode / claude-code / codex …）拿到这个技能目录后，
**只需要调用 `python3 run.py <子命令>`**，不必让用户先执行 install.sh、
也不必预先 `pip install numpy pyyaml`。第一次调用时自动补齐依赖。

三级策略（按顺序，命中即停）
----------------------------
1. **当前解释器已能 import numpy + yaml + openpyxl** → 直接用，零副作用。
   这是最常见的情况（科学计算机器/conda 环境多半已有 numpy），
   不建 venv、不写任何文件、耗时 < 50 ms。
2. **技能目录内 `.venv/`** → 建并安装（自包含，随技能目录走，不污染 home）。
3. **`~/.local/share/fluctuation-analysis/venv`** → 技能目录只读时回退。
   这两个位置也正是 install.sh 使用的，所以不会出现两份 venv。

逃生舱
------
`WINDOWCALC_PYTHON=/path/to/python`  指定一个已装好依赖的解释器，跳过自举。
`WINDOWCALC_VENV=/path/to/venv`      指定 venv 位置。
`WINDOWCALC_PIP_INDEX=https://…`     pip 源（内网/镜像）。

命令行
------
    python3 bootstrap.py            # 确保依赖，成功打印解释器路径
    python3 bootstrap.py --doctor   # 只诊断，不安装（agent 排查用）
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

__all__ = ["ensure", "doctor", "deps_ok", "missing_deps", "BootstrapError", "skill_dir"]

DEPS = ("numpy", "yaml", "openpyxl")       # import 名（openpyxl 让用户能直接给 .xlsx）
PKGS = {"numpy": "numpy", "yaml": "pyyaml", "openpyxl": "openpyxl"}  # import 名 → pip 名
MARKER = ".windowcalc-deps.json"             # 记录装过什么，便于 doctor
PIP_TIMEOUT = 900


class BootstrapError(RuntimeError):
    """自举失败。message 面向用户，hint 给下一步该怎么做。"""

    def __init__(self, message: str, hint: str = ""):
        super().__init__(message)
        self.hint = hint


# --------------------------------------------------------------------------
# 基础探测
# --------------------------------------------------------------------------
def skill_dir() -> Path:
    return Path(__file__).resolve().parent


def missing_deps() -> list[str]:
    """当前解释器缺哪些依赖（import 名）。"""
    missing = []
    for mod in DEPS:
        try:
            __import__(mod)
        except Exception:
            missing.append(mod)
    return missing


def deps_ok() -> bool:
    return not missing_deps()


def _venv_python(venv: Path) -> Path:
    if os.name == "nt":
        return venv / "Scripts" / "python.exe"
    return venv / "bin" / "python"


def _python_has_deps(py: str | Path) -> bool:
    """用一个子进程验证某个解释器能否 import 全部依赖（不抛异常，只返回真假）。"""
    try:
        out = subprocess.run(
            [str(py), "-c", "import " + ", ".join(DEPS)],
            capture_output=True, timeout=120,
        )
        return out.returncode == 0
    except Exception:
        return False


def _writable(path: Path) -> bool:
    """试探目录是否可写（os.access 在 ACL/只读挂载下不可靠）。"""
    try:
        path.mkdir(parents=True, exist_ok=True)
        probe = path / f".w-{os.getpid()}.tmp"
        probe.write_text("x", encoding="utf-8")
        probe.unlink()
        return True
    except Exception:
        return False


def candidate_venvs() -> list[Path]:
    env = os.environ.get("WINDOWCALC_VENV")
    if env:
        return [Path(env).expanduser()]
    local = skill_dir() / ".venv"
    home_venv = Path.home() / ".local" / "share" / "fluctuation-analysis" / "venv"
    if os.name == "nt" and len(str(local)) > 100:
        # Windows 有 260 字符路径上限，numpy 展开后很深；技能目录太深时改用短路径
        short = Path(os.environ.get("LOCALAPPDATA") or Path.home()) / "fluctuation-analysis" / "venv"
        return [short, local, home_venv]
    return [local, home_venv]


# --------------------------------------------------------------------------
# 并发保护：agent 可能同时并发调几条命令，两个 pip 写同一个 venv 会装坏
# --------------------------------------------------------------------------
class _Lock:
    def __init__(self, venv: Path):
        self.path = venv.with_suffix(".lock")

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        for _ in range(120):                      # 最多等 60 s
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, str(os.getpid()).encode())
                os.close(fd)
                return self
            except FileExistsError:
                # 陈旧锁（进程已死）直接清掉
                try:
                    age = time.time() - self.path.stat().st_mtime
                    if age > 300:
                        self.path.unlink(missing_ok=True)
                        continue
                except Exception:
                    pass
                time.sleep(0.5)
        return self                                # 拿不到锁也继续，不阻塞使用

    def __exit__(self, *exc):
        try:
            self.path.unlink(missing_ok=True)
        except Exception:
            pass


# --------------------------------------------------------------------------
# 建 venv + 装依赖
# --------------------------------------------------------------------------
def _build_venv(venv: Path, log) -> None:
    py = _venv_python(venv)
    if not py.exists():
        proc = subprocess.run(
            [sys.executable, "-m", "venv", str(venv)],
            capture_output=True, text=True, encoding="utf-8", errors="replace",
        )
        if proc.returncode != 0:
            raise _VenvBuildFailed((proc.stderr or proc.stdout or "").strip())

    pkgs = [PKGS[m] for m in DEPS if not _python_has_deps_single(py, m)]
    if not pkgs:                                   # 已有依赖
        return

    log(f"正在安装依赖（{', '.join(pkgs)}）…")
    last_err = ""
    for i, src in enumerate(_pip_sources()):
        cmd = [str(py), "-m", "pip", "install", "--quiet",
               "--disable-pip-version-check", "--no-input"]
        if src:
            cmd += ["--index-url", src]
        label = src or "默认源"
        if i:                                      # 第一个源失败时才提示换源
            log(f"默认源不可用，改用镜像：{label}")
        try:
            proc = subprocess.run(cmd + pkgs, capture_output=True, text=True,
                                  encoding="utf-8", errors="replace", timeout=PIP_TIMEOUT)
        except subprocess.TimeoutExpired:
            last_err = f"{label} 超时（{PIP_TIMEOUT}s）"
            continue
        if proc.returncode == 0:
            break
        last_err = (proc.stderr or proc.stdout or "").strip()
    else:
        raise _PipFailed(last_err)

    if not _python_has_deps(py):
        raise _PipFailed("pip 报告成功，但依赖仍无法 import")


# 国内可达性更好的镜像：默认源拉不动时自动兜底，让用户什么都不用配
MIRRORS = (
    "https://pypi.tuna.tsinghua.edu.cn/simple",
    "https://mirrors.aliyun.com/pypi/simple/",
)


def _pip_sources() -> list[str | None]:
    """按顺序尝试的 pip 源：显式指定 → 默认 → 内置镜像。"""
    out: list[str | None] = []
    env = os.environ.get("WINDOWCALC_PIP_INDEX")
    if env:
        out.append(env)
    out.append(None)
    out.extend(m for m in MIRRORS if m not in out)
    return out


def _python_has_deps_single(py: Path, mod: str) -> bool:
    try:
        return subprocess.run([str(py), "-c", f"import {mod}"],
                              capture_output=True, timeout=60).returncode == 0
    except Exception:
        return False


class _VenvBuildFailed(Exception):
    pass


class _PipFailed(Exception):
    pass


def _write_marker(venv: Path) -> None:
    try:
        (venv / MARKER).write_text(json.dumps({
            "installed": [PKGS[m] for m in DEPS],
            "python": sys.version.split()[0],
            "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }, ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception:
        pass


# --------------------------------------------------------------------------
# 主入口
# --------------------------------------------------------------------------
def ensure(force: bool = False, log=None) -> str | None:
    """确保依赖可用。

    返回 None  = 当前解释器就够用，不要换。
    返回 "路径" = 请用这个解释器重新执行（调用方负责 exec）。
    """
    log = log or _default_log

    if not force and deps_ok():
        return None

    # 已经自举过一轮还缺依赖 → 不要无限重入
    if os.environ.get("WINDOWCALC_BOOTSTRAP_DONE"):
        raise BootstrapError(
            "依赖自举后仍无法 import：" + ", ".join(missing_deps()),
            hint="可能 venv 被破坏。删掉下面两个目录后重试：\n"
                 f"  {skill_dir() / '.venv'}\n"
                 f"  {Path.home() / '.local/share/fluctuation-analysis/venv'}",
        )

    # 1) 逃生舱：指定解释器
    alt = os.environ.get("WINDOWCALC_PYTHON")
    if alt:
        if _python_has_deps(alt):
            return alt
        raise BootstrapError(
            f"WINDOWCALC_PYTHON={alt} 不可用，或它没有 numpy/pyyaml",
            hint="换一个已装好依赖的解释器，或取消该环境变量让程序自动安装。",
        )

    # 2) 已有可用 venv
    for venv in candidate_venvs():
        py = _venv_python(venv)
        if py.exists() and _python_has_deps(py):
            return str(py)

    # 3) 建 venv（逐候选；目录不可写或建不出来就换下一个）
    errors: list[str] = []
    for venv in candidate_venvs():
        if not _writable(venv.parent):
            errors.append(f"{venv}：目录不可写")
            continue
        try:
            with _Lock(venv):
                _build_venv(venv, log)
            _write_marker(venv)
            return str(_venv_python(venv))
        except _VenvBuildFailed as exc:
            errors.append(f"{venv}：建 venv 失败（{_summarize(str(exc))}）")
            continue
        except _PipFailed as exc:
            # pip 失败通常与目录无关（网络/权限/磁盘），不再换目录重试
            raise BootstrapError(
                "无法自动安装 numpy / pyyaml：" + _summarize(str(exc)),
                hint=_pip_hint(),
            )
        except subprocess.TimeoutExpired:
            raise BootstrapError(
                f"pip install 超时（{PIP_TIMEOUT}s）",
                hint="网络太慢或被墙。可设置镜像后重试：\n"
                     "  export WINDOWCALC_PIP_INDEX=https://pypi.tuna.tsinghua.edu.cn/simple\n"
                     "或指定已有解释器：\n"
                     "  export WINDOWCALC_PYTHON=$(which python3)",
            )

    raise BootstrapError(
        "无法自动准备 Python 依赖：" + "；".join(errors),
        hint=_pip_hint(),
    )


def _summarize(text: str, limit: int = 400) -> str:
    """从 pip 输出里提炼一句有用的话（首行常是没营养的 "ERROR: Exception:"）。"""
    lines = [l.strip() for l in (text or "").splitlines() if l.strip()]
    if not lines:
        return "（无输出）"
    keep = [l for l in lines if not l.startswith(("WARNING:", "warning:"))] or lines
    tail = keep[-1]
    head = next((l for l in keep if l.startswith("ERROR:") and "Exception:" not in l), "")
    msg = f"{head}｜{tail}" if head and head != tail else tail
    return msg[:limit]


def _pip_hint() -> str:
    return (
        "下一步（任选其一）：\n"
        "  1. 有网：确认能访问 PyPI，必要时设镜像\n"
        "       export WINDOWCALC_PIP_INDEX=https://pypi.tuna.tsinghua.edu.cn/simple\n"
        "     然后重跑同一条命令，会自动安装。\n"
        "  2. 已装好 numpy/pyyaml 的 Python（如 conda 环境）：\n"
        "       export WINDOWCALC_PYTHON=/path/to/python\n"
        "  3. 手动装：pip install numpy pyyaml\n"
        "  （本程序不会修改系统 Python，也不会写用户 site-packages）"
    )


def _default_log(msg: str) -> None:
    print(f"[windowcalc] {msg}", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# 诊断
# --------------------------------------------------------------------------
def doctor() -> int:
    """只诊断不安装。0=可直接用；2=能自动解决；1=需要人工介入。"""
    sd = skill_dir()
    print(f"技能目录 : {sd}")
    try:
        ver = (sd / "windowcalc" / "__init__.py").read_text(encoding="utf-8")
        for line in ver.splitlines():
            if "__version__" in line:
                print("技能版本 : " + line.split("=")[1].strip().strip('"'))
                break
    except Exception:
        pass
    print(f"解释器   : {sys.executable} ({sys.version.split()[0]})")

    missing = missing_deps()
    if not missing:
        import numpy, yaml, openpyxl
        print(f"依赖     : numpy {numpy.__version__} / pyyaml {yaml.__version__} / openpyxl {openpyxl.__version__}  ✅")
    else:
        print(f"依赖     : 缺 {', '.join(missing)}")

    print("venv 候选:")
    for venv in candidate_venvs():
        py = _venv_python(venv)
        if py.exists():
            ok = _python_has_deps(py)
            print(f"  {venv}  {'✅ 可用' if ok else '⚠️ 存在但缺依赖'}")
        else:
            print(f"  {venv}  （不存在）")

    alt = os.environ.get("WINDOWCALC_PYTHON")
    if alt:
        print(f"WINDOWCALC_PYTHON : {alt} -> {'✅' if _python_has_deps(alt) else '❌'}")
    if os.environ.get("WINDOWCALC_PIP_INDEX"):
        print(f"pip 源    : {os.environ['WINDOWCALC_PIP_INDEX']}")

    if not missing:
        print("\n结论: 可直接运行，无需安装。")
        return 0
    for venv in candidate_venvs():
        if _venv_python(venv).exists() and _python_has_deps(_venv_python(venv)):
            print(f"\n结论: 首次调用会自动切换到已有 venv：{_venv_python(venv)}")
            return 2
    print("\n结论: 缺依赖，首次调用会自动创建 venv 并安装（需能访问 PyPI）。")
    return 2


def main(argv: list[str]) -> int:
    if "--doctor" in argv:
        return doctor()
    try:
        target = ensure(force="--force" in argv,
                        log=lambda m: print(f"[windowcalc] {m}", file=sys.stderr, flush=True))
    except BootstrapError as exc:
        print(f"[windowcalc] 安装失败：{exc}", file=sys.stderr)
        if exc.hint:
            print(exc.hint, file=sys.stderr)
        return 1
    print(target or sys.executable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
