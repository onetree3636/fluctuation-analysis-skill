#!/usr/bin/env python3
"""windowcalc 跨平台入口（**免安装**：依赖缺了会自己装）。

    python3 run.py rate --input 计数.csv --nts 3e8     # Linux / macOS
    py -3 run.py rate --input 计数.csv --nts 3e8       # Windows

这是**给 agent 用的默认入口**。第一次调用时若缺 numpy / pyyaml，
会自动在技能目录内建 `.venv` 并装好（见 bootstrap.py），用户不需要做任何事。
已装好依赖的解释器则零副作用直接运行。

诊断：`python3 run.py --doctor`
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

# Windows 重定向输出时默认走本地代码页（GBK），"⚠" 之类会 UnicodeEncodeError
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8")
    except Exception:
        pass


def _relaunch(python: str) -> int:
    """用另一个解释器重跑本脚本（bootstrap 刚建好的 venv）。"""
    argv = [python, str(Path(__file__).resolve()), *sys.argv[1:]]
    env = dict(os.environ, WINDOWCALC_BOOTSTRAP_DONE="1", PYTHONUTF8="1")
    if os.name == "nt":
        import subprocess
        return subprocess.run(argv, env=env).returncode
    os.execve(python, argv, env)          # 不返回
    return 1


def main() -> int:
    import bootstrap

    if "--doctor" in sys.argv[1:]:
        return bootstrap.doctor()

    try:
        target = bootstrap.ensure()
    except bootstrap.BootstrapError as exc:
        print(f"[windowcalc] 无法准备运行环境：{exc}", file=sys.stderr)
        if exc.hint:
            print(exc.hint, file=sys.stderr)
        return 1

    if target:
        return _relaunch(target)

    from windowcalc.cli import main as cli_main
    return cli_main(sys.argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
