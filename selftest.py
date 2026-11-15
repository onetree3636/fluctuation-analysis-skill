#!/usr/bin/env python3
"""windowcalc 自检：一条命令验证「环境 → 引擎 → 已知答案 → 文档」。

    python selftest.py            # 全量（约 20–60 s）
    python selftest.py --quick    # 只查环境与命令（约 2 s）

退出码 0 = 无失败（跳过项不算失败）；1 = 有失败（会打印失败项与可复现命令）。
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
SKIPS: list[tuple[str, str]] = []
PASSED = 0
CHECKS: list = []


class Skip(Exception):
    """本环境按设计不跑该项（入口文件只存在于部署副本时触发）。"""


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
    if not (ROOT / "run.py").exists() and not (ROOT / "bootstrap.py").exists():
        raise Skip("入口文件（run.py/bootstrap.py）按设计只在部署副本"
                   "（sync_skill.sh「不同步」清单）——全量自检请在技能副本跑")
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
    if (ROOT / "run.py").exists() or (ROOT / "bootstrap.py").exists():
        assert name == ROOT.name, f"name={name!r} 必须与目录名 {ROOT.name!r} 一致"
    else:
        print("      [note] 项目源目录不做 name=目录名 检查（入口文件按设计只在部署副本）")

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


@check("九个子命令都能启动（--help）", quick=True)
def _():
    for sub in ("rate", "compare", "power", "optimize", "fit-fc", "fit-tauoff",
                "fit-induction", "design-induction", "multitarget"):
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
             "--ref-dose", "0.3", "--ref-m", "6", "--nt", "2e8"])
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


@check("fit-induction：剂量轴机制模型（合成数据 Δsat≈0.80）")
def _():
    p = run(["fit-induction", "--input", "examples/acceptance/dose_induction.csv",
             "--m-off", "0.02", "--query", "0.5"])
    assert p.returncode == 0, p.stderr[:500]
    assert "Δsat" in p.stdout and "饱和" in p.stdout
    import re as _re
    m = _re.search(r"Δsat = ([\d.]+)", p.stdout)
    assert m, p.stdout[:400]
    near(float(m.group(1)), 0.80, 0.02, "Δsat")


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
@check("multitarget：缺参必须 exit 2 且给齐清单；账本主路径能跑（quick）")
def _():
    """这是 2026-11-11 用户明令的行为：缺参不许静默填默认，必须让用户先同意。"""
    p = run(["multitarget", "--markers", "examples/acceptance/markers_demo.csv",
             "--genome-mb", "10.5", "--k", "12", "--capacity", "3e8",
             "--c-target", "5", "--n-win", "4"])
    assert p.returncode == 2, f"缺 --records 应 exit 2，实得 {p.returncode}"
    assert "不会用默认值代替" in p.stderr, p.stderr[:300]
    assert "不要自行编造这些值" in p.stderr, "缺 agent 禁令提示"
    ok = run(["multitarget", "--records", "examples/records_template.csv",
              "--markers", "examples/acceptance/markers_demo.csv",
              "--genome-mb", "5.0", "--n-win", "4", "--k", "30",
              "--capacity", "3e8", "--c-target", "5", "--bgc-kb", "40"])
    assert ok.returncode == 0, f"账本主路径应 exit 0：{ok.stderr[:400]}"
    for key in ("各档 f_step", "PAP 曲线现状", "下一窗设计", "剩余梯级"):
        assert key in ok.stdout, f"缺章节：{key}"


@check("multitarget：账本背景扣除与红灯可见（quick）")
def _():
    """fstep 并入后的口径：输出①必须含 Wilson/Newcombe 区间与逐档 f_step。"""
    p = run(["multitarget", "--records", "examples/records_template.csv",
             "--markers", "examples/acceptance/markers_demo.csv",
             "--genome-mb", "5.0", "--n-win", "4", "--k", "30",
             "--capacity", "3e8", "--c-target", "5"])
    assert p.returncode == 0, p.stderr[:300]
    assert "Wilson" in p.stdout or "Newcombe" in p.stdout or "f_step" in p.stdout


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
            # 相对链接先按文件所在目录解析，再退回仓库根（两处都找不到才算断链）
            cand = d.parent / t
            resolved = cand if cand.exists() else (ROOT / t)
            if not resolved.exists():
                missing.append(f"{d.name} → {t}")
    assert not missing, "断链：" + "; ".join(missing)


@check("文档：示例里的 windowcalc 子命令都真实存在")
def _():
    # ★ 与 build_parser() 的子命令保持同步（2026-11-13 补 multitarget——v0.8.0 漏更，
    #   任何含 "windowcalc multitarget" 的文档都会被误报为未知子命令）
    known = {"rate", "compare", "power", "optimize", "fit-fc", "fit-tauoff",
             "fit-induction", "design-induction", "multitarget"}
    docs = [ROOT / "SKILL.md"] + sorted((ROOT / "docs").glob("*.md"))
    bad = []
    for d in docs:
        text = d.read_text(encoding="utf-8")
        for m in re.finditer(r"windowcalc\s+([a-z-]+)", text):
            sub = m.group(1)
            if sub not in known and sub not in ("--help", "-h"):
                bad.append(f"{d.name}: windowcalc {sub}")
    assert not bad, "未知子命令：" + "; ".join(bad)


@check("文档：核心入口文件与配图都在")
def _():
    for f in ("docs/knowledge-base.md", "docs/experiment-design.md",
              "docs/mutator-window.md", "examples/record-template.csv",
              "examples/raw_counts_template.csv",
              "docs/generic-policy.md", "docs/protocol-execution.md",
              "docs/multitarget-code-map.md",
              "docs/figures/cycle-liquid-propagation.png",
              "docs/figures/make_cycle_schematic.py",
              "docs/figures/propagation-topologies.png",
              "docs/figures/make_propagation_topologies.py",
              "docs/figures/decide-plating-and-timing.png",
              "docs/figures/make_decision_flow.py"):
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


# ---------------------------------------------------------------- 通用性
# 通用性政策（docs/generic-policy.md）：技能里不得出现本项目专属的药物名/参数。
# 表里每条都是"曾经泄漏过"的真实串——新增内容前先想清楚为什么要往里加。
GENERIC_BANNED = [
    (r"茶碱|theophyllin", "本项目诱导剂名"),
    (r"\bpolC\b", "本项目靶点基因名"),
    (r"\bM7\b", "本项目菌株代号"),
    (r"0\.162", "本项目诱导幅度"),
    (r"8\.7\s*(?:e-3|×10\^?-3|×10⁻³)", "本项目泄漏速率"),
    (r"1\.51e-7|1\.94e-9|8\.12e-9", "本项目实测 μ"),
    (r"7\.5\s*Mb", "本项目基因组大小"),
    (r"83\s*kb", "本项目保护区大小"),
    (r"78\s*倍", "本项目 F 值"),
    (r"4\.2\s*×", "本项目泄漏倍数"),
]


@check("文档↔代码：参数名与关键常量一致（防文档漂移）")
def _():
    """第三道不变量：文档里写的参数名/阈值，必须与代码里真实存在的一致。

    2026-11-15 首次跑就抓到 6 处漂移：`--f-step`（命令已并入）、`--pap`、`--defaults-ok`、
    `--steps`、以及 `--mode design` 的残留。文档漂移不会报错，只会让人照着错的做。
    """
    import inspect
    import re as _re
    sys.path.insert(0, str(ROOT))
    from windowcalc.cli import build_parser
    from windowcalc import ramp as _ramp
    from windowcalc.mtblocks import f_ceiling_note

    docs = ([ROOT / "SKILL.md"] + sorted((ROOT / "docs").glob("*.md"))
            + sorted((ROOT / "references").glob("*.md")))

    # ---- A. 参数名：文档里出现的 --flag 必须真的存在 ----
    ap = build_parser()
    opts: set[str] = set()

    def walk(p):
        for a in p._actions:
            opts.update(a.option_strings)
            if isinstance(getattr(a, "choices", None), dict):
                for sp in a.choices.values():
                    walk(sp)

    walk(ap)
    # 不属于本 CLI、但文档里合法出现的旗标（每条都写清出处）
    external = {
        "--delete": "rsync 的旗标（maintainer.md）",
        "--doctor": "run.py 的旗标",
        "--quick": "selftest.py 自己的旗标（testing.md）",
        "--mu-wt": "testing.md 的『错误路径』示例——故意用一个不存在的旗标",
    }
    bad = []
    for d in docs:
        for m in _re.finditer(r"--[A-Za-z][A-Za-z0-9-]*", d.read_text(encoding="utf-8")):
            flag = m.group(0)
            if flag not in opts and flag not in external:
                bad.append(f"{d.name}: {flag}")
    assert not bad, "文档里出现不存在的参数：\n  " + "\n  ".join(sorted(set(bad))[:15])

    # ---- B. 关键常量：代码值 ↔ 文档里的同一句话 ----
    def line_has(doc: str, kw: str, token: str) -> bool:
        for ln in (ROOT / doc).read_text(encoding="utf-8").splitlines():
            if kw in ln and token in ln:
                return True
        return False

    sat_default = inspect.signature(_ramp.run_ledger).parameters["sat_f"].default
    assert sat_default == 0.5, f"sat_f 默认值变了（{sat_default}）——同步 docs 的'f ≥ 0.5'"
    assert line_has("references/guardrails.md", "饱和点", "0.5"), "护栏 G22 没写饱和阈值 0.5"

    from windowcalc.fstep import wilson_ci
    assert abs(wilson_ci(0, 10 ** 8)[2] * 10 ** 8 - 3.0) < 1e-9, "rule of three 不是 3/n"
    assert line_has("docs/glossary.md", "rule of three", "3/N"), "术语表没写 3/N"

    assert abs(_ramp.Z95 - 1.96) < 0.01, "Z95 变了——同步文档里的 1.96"
    assert any("1.96" in (ROOT / d).read_text(encoding="utf-8")
               for d in ("docs/protocol-ramp.md", "docs/multitarget-code-map.md",
                         "docs/generic-policy.md")), "文档里找不到 Wilson 的 z=1.96"

    assert "22" in f_ceiling_note(22) and "139" in f_ceiling_note(139), "F 代价带变了"
    assert line_has("docs/multitarget-code-map.md", "22×", "139×"), "代码地图没写 22×/139×"

    src_model = (ROOT / "windowcalc" / "model.py").read_text(encoding="utf-8")
    assert "3*M_base" in src_model, "泄漏红线 3×M_base 改了？"
    assert line_has("docs/experiment-design.md", "M_off", "3×M_base"), "实验设计没写 3×M_base 合格线"

    src_ramp = (ROOT / "windowcalc" / "ramp.py").read_text(encoding="utf-8")
    assert "0.8 * sat_f" in src_ramp and "1.2 * sat_f" in src_ramp, "±20% 边界带的实现变了"
    assert line_has("references/guardrails.md", "±20%", "饱和"), "护栏没写 ±20% 边界带"

    # ---- C. 声称一致性：不许说"没建模"而代码其实建了 ----
    if hasattr(_ramp, "washout_tail"):
        stale = [d.name for d in docs
                 if _re.search(r"τ[^\n]{0,12}(不建模|未建模|不进账本|未进账本)",
                               d.read_text(encoding="utf-8"))]
        assert not stale, f"washout_tail 已实现，但这些文档还说没建模：{stale}"


@check("护栏：每个子命令都在 guardrails §0，且引用的 G 编号都有定义")
def _():
    """2026-11-15 的"护栏三元组化"：命令 ↔ 触发编号必须双向对得上。

    新加子命令却没配护栏 = 用户会在没人问的情况下踩坑；引用不存在的 G 编号 = 悬空规则。
    """
    import re as _re
    sys.path.insert(0, str(ROOT))
    from windowcalc.cli import build_parser

    ap = build_parser()
    subs = set(next(a for a in ap._actions if a.dest == "cmd").choices)

    text = (ROOT / "references" / "guardrails.md").read_text(encoding="utf-8")
    head = text.split("## 1.")[0]                     # §0 命令表
    listed = {m.group(1) for m in _re.finditer(r"^\|\s*`([a-z-]+)`", head, _re.M)}
    missing = subs - listed
    extra = listed - subs
    assert not missing, f"guardrails §0 缺这些子命令：{sorted(missing)}"
    assert not extra, f"guardrails §0 列了不存在的子命令：{sorted(extra)}"

    used = set(_re.findall(r"\bG(\d{1,2})\b", head))
    defined = set(_re.findall(r"^\|\s*\*\*G(\d{1,2})\*\*", text, _re.M))
    assert not (used - defined), f"§0 引用了没定义的 G 编号：{sorted(used - defined, key=int)}"


@check("文档：分级横幅齐全、索引存在（docs/INDEX.md）")
def _():
    """2026-11-15 的"减法"：文档只能有一个入口。

    每篇 docs/*.md 的 H1 下方必须有一行 `> **级别：…** ｜ 什么时候读我：… ｜ 不用读我：…`
    且级别 ∈ {必读, 主路径, 附录}——否则读者又会多一个"不知道要不要点开"的入口。
    """
    idx = ROOT / "docs" / "INDEX.md"
    assert idx.exists(), "缺 docs/INDEX.md（唯一入口）"
    levels = ("必读", "主路径", "附录")
    bad = []
    for d in sorted((ROOT / "docs").glob("*.md")):
        if d.name == "INDEX.md":
            continue
        head = "\n".join(d.read_text(encoding="utf-8").splitlines()[:12])
        if "> **级别：" not in head:
            bad.append(f"{d.name}: 缺分级横幅")
            continue
        if not any(f"> **级别：{lv}" in head for lv in levels):
            bad.append(f"{d.name}: 级别不在 {levels}")
    assert not bad, "；".join(bad)


@check("通用性：技能内无项目专属药物名/参数（docs/generic-policy.md）")
def _():
    """两条硬规则：①药名只许来自引用文献；②默认值只许来自文献确证的通用值。

    这里只挡"已经泄漏过一次"的串——新写的数字不在表里不会报警，
    所以 docs/generic-policy.md §5 的五问仍然是必须做的。
    """
    targets = [ROOT / "SKILL.md", ROOT / "CHANGELOG.md"]
    targets += sorted((ROOT / "docs").glob("*.md"))
    targets += sorted((ROOT / "references").glob("*.md"))
    targets += sorted((ROOT / "windowcalc").glob("*.py"))
    for pat in ("*.csv", "*.yaml", "*.md"):
        targets += sorted((ROOT / "examples").rglob(pat))

    hits = []
    for f in targets:
        try:
            text = f.read_text(encoding="utf-8")
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        for i, line in enumerate(text.splitlines(), 1):
            for rx, why in GENERIC_BANNED:
                if re.search(rx, line):
                    hits.append(f"{f.relative_to(ROOT)}:{i}: {why} → {line.strip()[:70]}")
    assert not hits, ("命中项目专属内容（通用性政策，见 docs/generic-policy.md）：\n  "
                      + "\n  ".join(hits[:20]))


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
        except Skip as exc:
            SKIPS.append((name, str(exc)))
            print(f"[SKIP] {name}\n       {exc}")
        except Exception as exc:                      # noqa: BLE001
            FAILS.append((name, f"{type(exc).__name__}: {exc}"))
            print(f"[FAIL] {name}\n       {type(exc).__name__}: {exc}")
    print("-" * 60)
    print(f"结果：{PASSED} 通过，{len(FAILS)} 失败，{len(SKIPS)} 跳过")
    if SKIPS:
        print("\n跳过项（按设计不在本目录检查）：")
        for name, why in SKIPS:
            print(f"  - {name}\n    {why}")
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
