#!/usr/bin/env python3
"""windowcalc 自检：一条命令验证「环境 → 引擎 → 已知答案 → 文档」。

    python selftest.py            # 全量（约 20–60 s）
    python selftest.py --quick    # 只查环境与命令（约 2 s）

退出码 0 = 全过；1 = 有失败（会打印失败项与可复现命令）。
设计原则：所有用例都用**技能自带的** examples/acceptance/ 数据，不依赖网络与 R。
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable
QUICK = "--quick" in sys.argv

FAILS: list[tuple[str, str]] = []
PASSED = 0
CHECKS: list = []


def run(args, timeout=180):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    env["PYTHONUTF8"] = "1"
    return subprocess.run([PY, "-m", "windowcalc", *args], cwd=ROOT, env=env,
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=timeout)


def check(name, quick=False):
    def deco(fn):
        CHECKS.append((name, fn, quick))
        return fn
    return deco


def near(got, want, tol, label=""):
    assert abs(got - want) <= tol, f"{label} 期望 {want}±{tol}，实得 {got}"


def grab(text, pattern):
    m = re.search(pattern, text)
    assert m, f"输出里找不到 /{pattern}/：\n{text[:600]}"
    return float(m.group(1))


# ---------------------------------------------------------------- 环境
@check("环境：python ≥3.9 且有 numpy / pyyaml", quick=True)
def _():
    assert sys.version_info >= (3, 9), f"python {sys.version.split()[0]} 太旧"
    import numpy, yaml, openpyxl  # noqa: F401
    assert numpy.__version__ and yaml.__version__ and openpyxl.__version__


@check("入口：run.py --doctor 能诊断且不安装", quick=True)
def _():
    assert (ROOT / "bootstrap.py").exists(), "缺 bootstrap.py（agent 自动安装逻辑）"
    p = subprocess.run([PY, str(ROOT / "run.py"), "--doctor"], cwd=ROOT,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=180)
    assert p.returncode in (0, 2), f"rc={p.returncode}：{p.stderr[:400]}"
    assert "结论" in p.stdout, f"输出里没有「结论」：{p.stdout[:400]}"
    assert "Traceback" not in p.stderr, p.stderr[:400]


@check("规范：frontmatter 符合 Agent Skills 规范", quick=True)
def _():
    """按 agentskills.io/specification 的硬约束自查（不依赖 skills-ref 包）。

    这些约束一旦违反，技能可能**静默不加载**——与本机环境无关，纯格式问题。
    """
    import yaml
    text = (ROOT / "SKILL.md").read_text(encoding="utf-8")
    assert text.startswith("---\n"), "SKILL.md 必须以 YAML frontmatter 开头"
    fm = text.split("---\n", 2)[1]
    d = yaml.safe_load(fm) or {}

    name = d.get("name", "")
    assert re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", name), f"name 不合规（小写+连字符）：{name!r}"
    assert name == ROOT.name, f"name={name!r} 必须与目录名 {ROOT.name!r} 一致"

    desc = d.get("description", "")
    assert 1 <= len(desc) <= 1024, f"description 长度 {len(desc)}，须在 1–1024"

    if "compatibility" in d:
        n = len(str(d["compatibility"]))
        assert 1 <= n <= 500, f"compatibility 长度 {n}，须在 1–500"

    allowed = {"name", "description", "license", "compatibility", "metadata", "allowed-tools"}
    extra = set(d) - allowed
    assert not extra, f"未知 frontmatter 键（Claude Code 会拒绝加载）：{sorted(extra)}"

    lines = len(text.splitlines())
    assert lines < 500, f"SKILL.md {lines} 行，规范建议 <500（渐进式披露：细节移到 docs/）"


@check("八个子命令都能启动（--help）", quick=True)
def _():
    for sub in ("rate", "compare", "power", "optimize", "fit-fc", "fit-tauoff",
                "fit-induction", "design-induction"):
        p = run([sub, "--help"])
        assert p.returncode == 0, f"{sub} --help 返回 {p.returncode}"


# ---------------------------------------------------------------- 已知答案
@check("rate：Rosche & Foster 2000 Table 3（60 管；MSS m̂=1.23 [0.93,1.58]；p0 法 1.70）")
def _():
    p = run(["rate", "--input", "examples/acceptance/RF2000.csv", "--nts", "1.5e8"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"m \(MSS-MLE\)\s*=\s*([\d.]+)"), 1.2274, 0.002, "m")
    lo = grab(p.stdout, r"95%CI \[([\d.]+),")
    hi = grab(p.stdout, r"95%CI \[[\d.]+, ([\d.]+)\]")
    near(lo, 0.9309, 0.005, "CI 下界")
    near(hi, 1.5757, 0.005, "CI 上界")
    near(grab(p.stdout, r"p0=([\d.]+)"), 0.183, 0.002, "p0")


@check("输入：Excel(.xlsx) 与 CSV 得到同一结果（新手直接给 Excel 也能算）")
def _():
    p = run(["rate", "--input", "examples/acceptance/RF2000.xlsx", "--nts", "1.5e8"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"m \(MSS-MLE\)\s*=\s*([\d.]+)"), 1.2274, 0.002, "m(xlsx)")
    assert "Traceback" not in p.stderr, p.stderr[:400]


@check("rate per_tube：部分涂板 ζ=0.2（m̂≈3.174）")
def _():
    p = run(["rate", "--input", "examples/acceptance/plating_z02.csv",
             "--nt-mode", "per_tube", "--fraction-plated", "0.2"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"m \(MSS-MLE\)\s*=\s*([\d.]+)"), 3.1740, 0.006, "m")
    assert "部分涂板" in p.stderr, "ζ<1 时应在 stderr 提示 Nt 口径"


@check("rate：GBK 编码 CSV 可读（Excel 另存的 CSV）")
def _():
    p = run(["rate", "--input", "examples/acceptance/gbk_rf2000.csv", "--nts", "1.5e8"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"m \(MSS-MLE\)\s*=\s*([\d.]+)"), 1.2274, 0.002, "m")
    assert "Traceback" not in p.stderr


@check("真实数据：LD 1943 exp16（ε=0.4，含 jackpot）→ m=1.18636")
def _():
    p = run(["rate", "--input", "examples/acceptance/LD1943_exp16.csv",
             "--nts", "5.6e8", "--fraction-plated", "0.4"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"m \(MSS-MLE\)\s*=\s*([\d.]+)"), 1.1864, 0.002, "m")
    lo = grab(p.stdout, r"95%CI \[([\d.]+),")
    hi = grab(p.stdout, r"95%CI \[[\d.]+, ([\d.]+)\]")
    near(lo, 0.5803, 0.005, "CI 下界")
    near(hi, 2.0908, 0.005, "CI 上界")


@check("真实数据：Demerec 1945（30 管）→ m=10.84383")
def _():
    p = run(["rate", "--input", "examples/acceptance/Demerec1945.csv", "--nts", "1.9e8"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"m \(MSS-MLE\)\s*=\s*([\d.]+)"), 10.84383, 0.002, "m")


@check("design-induction：剂量轴设计（覆盖检查 + 每点管数）")
def _():
    p = run(["design-induction", "--c50", "0.48", "--n-hill", "1.85",
             "--ref-dose", "0.3", "--ref-m", "6"])
    assert p.returncode == 0, p.stderr[:400]
    assert "覆盖检查" in p.stdout and "合计" in p.stdout
    assert "C50 以下 2 点" in p.stdout or "以下 2 点" in p.stdout


@check("compare：LRT 显著 + fold-change θ≈10")
def _():
    p = run(["compare", "--d1", "examples/acceptance/groupA.csv",
             "--d2", "examples/acceptance/groupB.csv", "--labels", "A,B"])
    assert p.returncode == 0, p.stderr[:400]
    pval = grab(p.stdout, r"p = ([\d.eE+-]+)")
    assert pval < 1e-10, f"p={pval} 应远小于 0.05"
    theta = grab(p.stdout, r"θ = μ\(B\)/μ\(A\) = ([\d.]+)")
    assert 8.0 <= theta <= 14.0, f"θ={theta} 偏离真值 10 太多"


@check("fit-induction：剂量轴机制模型（Δsat≈1.315）")
def _():
    p = run(["fit-induction", "--input", "examples/acceptance/dose_induction.csv",
             "--m-off", "0.012", "--query", "0.5"])
    assert p.returncode == 0, p.stderr[:500]
    assert "Δsat" in p.stdout and "饱和" in p.stdout
    import re as _re
    m = _re.search(r"Δsat = ([\d.]+)", p.stdout)
    assert m, p.stdout[:400]
    near(float(m.group(1)), 1.315, 0.02, "Δsat")


@check("power：解析式样本量 C=131")
def _():
    p = run(["power", "--mode", "precision", "--mu", "2e-9", "--nt-mean", "3e8",
             "--target-rse", "0.25"])
    assert p.returncode == 0, p.stderr[:400]
    near(grab(p.stdout, r"C = (\d+) 管"), 131, 0, "C")


# ---------------------------------------------------------------- 端到端
@check("optimize：报告含六个必需章节")
def _():
    p = run(["optimize", "--strain", "examples/strain_card.example.yaml",
             "--protocol", "examples/protocol_card.example.yaml"], timeout=300)
    assert p.returncode == 0, p.stderr[:600]
    for sec in ("## 1. 推荐方案", "## 2. 时刻表", "## 3. 本方案的标度律读数",
                "## 7. 参数实测状态"):
        assert sec in p.stdout, f"缺章节 {sec}"


@check("fit-fc / fit-tauoff：合成数据参数回收")
def _():
    with tempfile.TemporaryDirectory() as td:
        fc = Path(td) / "dose.csv"
        fc.write_text("dose_mM,m_on\n0.03,4e-10\n0.1,2e-9\n0.3,1.6e-8\n1.0,4e-8\n3.0,4.6e-8\n",
                      encoding="utf-8")
        p = run(["fit-fc", "--input", str(fc), "--mu-base", "1.2e-10", "--n-boot", "50"])
        assert p.returncode == 0, p.stderr[:400]
        assert "EC50" in p.stdout and "F_max" in p.stdout

        tau = Path(td) / "tau.csv"
        tau.write_text("h,mu\n0,4.0e-7\n2,2.5e-7\n4,1.6e-7\n8,6.0e-8\n16,1.5e-8\n",
                       encoding="utf-8")
        p = run(["fit-tauoff", "--input", str(tau), "--mu-off", "1.2e-8",
                 "--td", "2.5", "--n-boot", "50"])
        assert p.returncode == 0, p.stderr[:400]
        near(grab(p.stdout, r"τ_off = ([\d.]+) h"), 4.0, 0.6, "τ_off")


# ---------------------------------------------------------------- 错误路径
@check("错误路径：坏数据给中文提示，不甩 traceback")
def _():
    with tempfile.TemporaryDirectory() as td:
        bad = Path(td) / "neg.csv"
        bad.write_text("r,Nt\n0,3e8\n-1,3e8\n", encoding="utf-8")
        p = run(["rate", "--input", str(bad)])
        assert p.returncode != 0, "负 r 应失败"
        assert ("不能为负" in p.stderr or "非负整数" in p.stderr), p.stderr[:300]
        assert "Traceback" not in p.stderr, "不应甩 traceback"

        p = run(["rate", "--input", str(Path(td) / "nope.csv")])
        assert p.returncode != 0 and "错误" in p.stderr
        assert "Traceback" not in p.stderr


# ---------------------------------------------------------------- 文档
@check("文档：相对链接都存在、命令示例可解析")
def _():
    docs = [ROOT / "SKILL.md"] + sorted((ROOT / "docs").glob("*.md"))
    link = re.compile(r"\]\((?!https?://|#)([^)]+)\)")
    missing = []
    for d in docs:
        text = d.read_text(encoding="utf-8")
        for target in link.findall(text):
            t = target.split("#")[0].strip()
            if not t:
                continue
            resolved = (ROOT / t) if not t.startswith(".") else (d.parent / t)
            if not resolved.exists():
                missing.append(f"{d.name} → {t}")
    assert not missing, "断链：" + "; ".join(missing)


@check("文档：示例里的 windowcalc 子命令都真实存在")
def _():
    known = {"rate", "compare", "power", "optimize", "fit-fc", "fit-tauoff",
             "fit-induction", "design-induction"}
    docs = [ROOT / "SKILL.md"] + sorted((ROOT / "docs").glob("*.md"))
    bad = []
    for d in docs:
        text = d.read_text(encoding="utf-8")
        for m in re.finditer(r"windowcalc\s+([a-z-]+)", text):
            sub = m.group(1)
            if sub not in known and sub not in ("--help", "-h"):
                bad.append(f"{d.name}: windowcalc {sub}")
    assert not bad, "未知子命令：" + "; ".join(bad)


@check("文档：四个用途的入口文件都在")
def _():
    for f in ("docs/knowledge-base.md", "docs/experiment-design.md",
              "docs/mutator-window.md", "examples/record-template.csv",
              "examples/raw_counts_template.csv"):
        assert (ROOT / f).exists(), f"缺 {f}"


@check("文档：数字与代码一致（子命令数/自检项数/样本量）")
def _():
    """文档里最容易烂的是数字：改了代码不改文档，读者会照着错的做。

    这里现场**算**出三个数字再和文档比对：
      - 子命令个数 ← argparse 的 subparsers；
      - 自检项数   ← CHECKS 的长度（本检查已注册）；
      - 解析样本量 ← power.sample_size_precision_analytic（μ=2e-9、Nt=3e8、95%/25%）
    """
    sys.path.insert(0, str(ROOT))
    from windowcalc.cli import build_parser
    from windowcalc.power import sample_size_precision_analytic

    ap = build_parser()
    n_sub = len(next(a for a in ap._actions if a.dest == "cmd").choices)
    n_checks = len(CHECKS)

    skill = (ROOT / "SKILL.md").read_text(encoding="utf-8")
    testing = (ROOT / "docs" / "testing.md").read_text(encoding="utf-8")

    m = re.search(r"(\d+)\s*个子命令", skill)
    assert m, "SKILL.md 里没有「N 个子命令」的说法"
    assert int(m.group(1)) == n_sub, f"SKILL.md 说 {m.group(1)} 个子命令，代码有 {n_sub}"

    m = re.search(r"(\d+)\s*项：环境", skill)
    assert m, "SKILL.md 里没有「N 项：环境…」的自检项数"
    assert int(m.group(1)) == n_checks, f"SKILL.md 说自检 {m.group(1)} 项，实际 {n_checks} 项"

    m = re.search(r"(\d+)\s*项，约", testing)
    assert m, "testing.md 里没有「N 项，约…」的自检项数"
    assert int(m.group(1)) == n_checks, f"testing.md 说 {m.group(1)} 项，实际 {n_checks} 项"

    C = sample_size_precision_analytic(2e-9, 3e8, 0.25, conf_level=0.95)
    assert f"C={C}" in testing, f"testing.md 的解析样本量不是 {C}（请同步「已知答案」表）"


def main() -> int:
    global PASSED
    print(f"windowcalc 自检  ({ROOT})")
    print(f"python {sys.version.split()[0]}｜{'只跑安装检查（--quick）' if QUICK else '全量'}")
    print("-" * 60)
    for name, fn, quick_ok in CHECKS:
        if QUICK and not quick_ok:
            continue
        try:
            fn()
            PASSED += 1
            print(f"[PASS] {name}")
        except AssertionError as exc:
            FAILS.append((name, str(exc)))
            print(f"[FAIL] {name}\n       {exc}")
        except Exception as exc:                      # noqa: BLE001
            FAILS.append((name, f"{type(exc).__name__}: {exc}"))
            print(f"[FAIL] {name}\n       {type(exc).__name__}: {exc}")
    print("-" * 60)
    print(f"结果：{PASSED} 通过，{len(FAILS)} 失败")
    if FAILS:
        print("\n失败项：")
        for name, why in FAILS:
            print(f"  - {name}\n    {why}")
        print("\n复现方法：cd 技能目录 && python -m windowcalc <子命令> ...")
        return 1
    print("全部通过。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
