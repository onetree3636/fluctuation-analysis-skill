"""对拍测试：**同一批公式的两套实现必须给出同一个数**。

背景（2026-11-15 代码审查的发现）：

  · `mtblocks.window_table_f_step()` / `minimal_n0()` 是文档里的"主路径积木"
    （单窗：给 f_step 与 Nt，解出最大 N₀ 与乘客）；
  · CLI 真正跑的是 `ramp.run_ledger()` → `plan_ramp()`，里面**重写了一遍**
    同样的账（N₀ = Nt − k/f、乘客 = 代数 × M、失守 = 1 − exp(−E[P]·L/G)）。

两套实现的公式现在一致，但**没有任何测试锁住**——改一处、忘另一处就会静默分叉
（历史上真发生过同类事故：BGC 分母写成 1e6、Wilson 参数顺序写反）。

本文件把三件事钉死：
  1. `minimal_n0` 的代数与账本里 N₀ 设计行（ramp.py `nt − k/f`）同源；
  2. `window_table_f_step` 的乘客 == 账本 `base_passengers`（同一 m_on、同一 n_win）；
  3. `window_table_f_step` 的失守概率 == `1 − exp(−E[P] × L/G)`，
     且梯级的累计乘客严格等于「基数 + i × 每档」，`n_post` 只加在最后一档。

改公式时：先跑本文件。它红 = 两处已经分叉。
"""

import math

import pytest

from windowcalc.mtblocks import Marker, minimal_n0, window_table_f_step
from windowcalc.ramp import run_ledger

# ---- 测试用的**通用**参数（不是任何项目实测值）-----------------------
MK = [Marker("A", 2.0e-8 * 7, 7)]      # μ_perbp = 2e-8
GENOME_MB = 10.0                        # 基因组大小（Mb）
F_STEP = 1.0e-6                         # 供给（克隆数/CFU）
K = 12.0                                # 目标克隆数
NT = 1.0e9                              # 饱和密度
BGC_KB = 100.0                          # 保护区大小


def _records_one_cycle() -> list[tuple[str, float, str, float, int]]:
    """一个周期的账本：1 个饱和点 + 2 个尾点（够 fit_pap 的 2 点下限）。"""
    return [
        ("1", 1.0, "post", 1.0e9, int(8.0e8)),   # f ≈ 0.8 → 饱和，用来定位 c_now
        ("1", 2.0, "post", 1.0e9, 1000),         # f = 1e-6
        ("1", 3.0, "post", 1.0e9, 100),          # f = 1e-7
    ]


# ---------------------------------------------------------------- 1. N₀ 同源

@pytest.mark.parametrize("f,k,nt", [
    (1e-6, 12.0, 1e9),
    (3e-7, 12.0, 2e9),
    (1e-5, 50.0, 5e8),
])
def test_minimal_n0_matches_ledger_n0_row(f, k, nt):
    """`minimal_n0` 必须等于解析式 `Nt − k/f`（积木版自身自洽）。"""
    r = minimal_n0(f, k, nt)
    assert r["feasible"]
    assert r["n0"] == pytest.approx(nt - k / f, rel=1e-12)
    # 代数口径自洽：ΔN = k/f 且 n_gen = log2(Nt/N0)
    assert r["n_gen"] == pytest.approx(math.log2(nt / r["n0"]), rel=1e-12)
    # 供给回头：rate × (Nt − N0) 必须刚好等于 k（取等号）
    assert r["events"] == pytest.approx(k, rel=1e-12)


@pytest.mark.parametrize("f,k,nt", [
    (1e-6, 12.0, 1e9),
    (1e-5, 50.0, 5e8),
])
def test_rendered_n0_design_line_is_correct(f, k, nt):
    """**实跑路径**的 N₀ 设计行（`render_ledger`）必须与解析式一致。

    2026-11-15 补：此前这一行只有 `minimal_n0` 一侧被断言，而真正打印给用户的
    `ramp.render_ledger` 那行**零测试覆盖**——改坏了 19 项 pin + 55 项其余测试全绿
    （见 docs/multitarget-code-map.md §6 的变异测试记录）。
    """
    from windowcalc.ramp import render_ledger

    plan = run_ledger(_records_one_cycle(), m_on=0.1, m_off=0.01,
                      genome_mb=GENOME_MB, n_win=4.0, n_plate=0.0, n_post=0.0,
                      k=k, capacity=1e9, c_target=6.0, bgc_kb=None, margin=3.0,
                      c_start=1.0, nt=nt)
    out = render_ledger(plan)
    # 该行用的是"尾部里最大的 f"（升档用），本账本里 = 1e-6
    r1_f = plan["ladder"][0]["f_step"] if plan["status"] == "可行" else f
    expect = nt - k / r1_f
    assert f"N₀ ≤ Nt_sat − k/f_step" in out, "报告少了 N₀ 设计行"
    assert f"= {expect:.3g}（" in out, (
        f"N₀ 设计行算错：应出现 '= {expect:.3g}（'；实际报告片段：\n"
        + "\n".join(ln for ln in out.splitlines() if "N₀" in ln))


def test_minimal_n0_generations_are_consistent_with_divisions():
    """n_gen 必须满足 ΔN = N0·(2^n − 1)——这是"代数 vs 总分裂数"的换算桥。"""
    r = minimal_n0(1e-6, 12.0, 1e9)
    n, n0 = r["n_gen"], r["n0"]
    assert n0 * (2.0 ** n - 1.0) == pytest.approx(r["delta_n"], rel=1e-9)


# ---------------------------------------------------------------- 2. 乘客同源

def test_window_table_passengers_equal_ledger_base():
    """mtblocks 的"单窗乘客" == 账本的"已完成档基数"（同一 m_on、同一 n_win）。"""
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    assert w["row"]["feasible"]

    plan = run_ledger(_records_one_cycle(),
                      m_on=w["m_genome"], m_off=None, genome_mb=GENOME_MB,
                      n_win=w["row"]["n_gen"],      # ← 用 mtblocks 解出的代数
                      n_plate=0.0, n_post=0.0,
                      k=K, capacity=1.0e9, c_target=10.0, bgc_kb=BGC_KB,
                      margin=3.0)
    assert plan["cycles"] == ["1"]
    assert plan["base_passengers"] == pytest.approx(w["row"]["passengers"], rel=1e-12)

    # 并且两边都等于 n_gen × (μ_perbp × G)——这是文档里的三段账第一项
    expect = w["row"]["n_gen"] * w["m_genome"]
    assert w["row"]["passengers"] == pytest.approx(expect, rel=1e-12)
    assert plan["base_passengers"] == pytest.approx(expect, rel=1e-12)


# ---------------------------------------------------------------- 3. 失守概率同源

def test_window_table_bgc_equals_poisson_single_shot():
    """失守 = 1 − exp(−E[P] × L/G)：单次代入，绝不做每档概率相加。"""
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT,
                            genome_mb=GENOME_MB, bgc_kb=BGC_KB)
    p = w["row"]["passengers"]
    frac = BGC_KB * 1000.0 / (GENOME_MB * 1e6)
    assert w["bgc_frac"] == pytest.approx(frac, rel=1e-12)
    assert w["row"]["p_bgc"] == pytest.approx(1.0 - math.exp(-p * frac), rel=1e-12)


def test_ladder_cumulative_matches_three_segment_ledger():
    """梯级累计 = 基数 + i ×（n_win·M_on + n_plate·M_off）；n_post 只加在最后一档。"""
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    n_win, n_plate, n_post = 4.0, 2.0, 3.0
    m_on, m_off = w["m_genome"], 0.01

    plan = run_ledger(_records_one_cycle(), m_on=m_on, m_off=m_off,
                      genome_mb=GENOME_MB, n_win=n_win, n_plate=n_plate,
                      n_post=n_post, k=K, capacity=1.0e9, c_target=10.0,
                      bgc_kb=None, margin=3.0)
    assert plan["status"] == "可行"
    base = plan["base_passengers"]
    assert base == pytest.approx(n_win * m_on + n_plate * m_off, rel=1e-12)

    per_rung = n_win * m_on + n_plate * m_off
    ladder = plan["ladder"]
    for i, r in enumerate(ladder, start=1):
        assert r["p_win"] == pytest.approx(n_win * m_on, rel=1e-12)
        assert r["p_plate"] == pytest.approx(n_plate * m_off, rel=1e-12)
        assert r["passengers"] == pytest.approx(per_rung, rel=1e-12)
        cum = base + i * per_rung
        if i == len(ladder):
            cum += n_post * m_off        # n_post 只加在**最后一档**的累计里
        assert r["cum_passengers"] == pytest.approx(cum, rel=1e-12)

    # 总乘客：基数 + 剩余梯级 + n_post（只一次）
    assert plan["p_total"] == pytest.approx(
        base + len(ladder) * per_rung + n_post * m_off, rel=1e-12)


# ---------------------------------------------------------------- 4. 洗脱尾（2026-11-15 新增）

def test_washout_tail_closed_form():
    """tail = (M_on − M_off)/(e^(td/τ) − 1)：按"每代 td 小时"的离散求和。"""
    from windowcalc.ramp import washout_tail
    m_on, m_off, tau, td = 0.2, 0.01, 5.0, 3.0
    assert washout_tail(m_on, m_off, tau, td) == pytest.approx(
        (m_on - m_off) / math.expm1(td / tau), rel=1e-12)
    # 都不给 → 0（报告里会注明"未计入"）
    assert washout_tail(m_on, m_off, None, None) == 0.0


def test_washout_tail_continuous_limit():
    """td ≪ τ 时退化为 (M_on−M_off)·τ/td，与 fit_tauoff 的尾部积分同口径。"""
    from windowcalc.ramp import washout_tail
    m_on, m_off, tau, td = 0.2, 0.01, 1.0e4, 1.0
    # 离散和与连续积分差一个 e^x/(e^x−1) 因子，x=td/τ≪1 时 ≈ 1 + x/2
    assert washout_tail(m_on, m_off, tau, td) == pytest.approx(
        (m_on - m_off) * tau / td, rel=1e-3)
    assert washout_tail(m_on, m_off, tau, td) < (m_on - m_off) * tau / td


def test_washout_tail_requires_both_parameters():
    from windowcalc.ramp import washout_tail
    with pytest.raises(ValueError, match="同时"):
        washout_tail(0.2, 0.01, 5.0, None)
    with pytest.raises(ValueError, match="同时"):
        washout_tail(0.2, 0.01, None, 3.0)


def test_ledger_counts_tail_once_per_washout():
    """每次洗脱计一次：已完成档基数 + 每剩余档 + 终点后。"""
    from windowcalc.ramp import washout_tail
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    m_on, m_off = w["m_genome"], 0.01
    n_win, n_plate, n_post = 4.0, 2.0, 3.0
    tail = washout_tail(m_on, m_off, 5.0, 3.0)

    plan = run_ledger(_records_one_cycle(), m_on=m_on, m_off=m_off,
                      genome_mb=GENOME_MB, n_win=n_win, n_plate=n_plate,
                      n_post=n_post, k=K, capacity=1.0e9, c_target=10.0,
                      bgc_kb=None, margin=3.0, tau_off_h=5.0, td_h=3.0)
    assert plan["p_tail"] == pytest.approx(tail, rel=1e-12)
    # 基数含尾
    assert plan["base_passengers"] == pytest.approx(
        n_win * m_on + n_plate * m_off + tail, rel=1e-12)
    # 每档含尾
    for r in plan["ladder"]:
        assert r["p_tail"] == pytest.approx(tail, rel=1e-12)
        assert r["passengers"] == pytest.approx(
            n_win * m_on + n_plate * m_off + tail, rel=1e-12)
    # 总账 = 基数 + 剩余档 + 终点后段（含尾）
    assert plan["p_total"] == pytest.approx(
        plan["base_passengers"]
        + len(plan["ladder"]) * (n_win * m_on + n_plate * m_off + tail)
        + n_post * m_off + tail, rel=1e-12)


# ---------------------------------------------------------------- 5. 饱和阈值

def _records_with_a_borderline_point():
    """一个点正好落在 0.6——阈值 0.5 时是饱和、0.7 时进尾部拟合。"""
    return [
        ("1", 1.0, "post", 1.0e9, int(6.0e8)),   # f = 0.6（边界点）
        ("1", 2.0, "post", 1.0e9, 1000),         # f = 1e-6
        ("1", 3.0, "post", 1.0e9, 100),          # f = 1e-7
    ]


def test_saturation_threshold_is_a_parameter_not_a_hardcode():
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    common = dict(m_on=w["m_genome"], m_off=None, genome_mb=GENOME_MB,
                  n_win=4.0, n_plate=0.0, n_post=0.0, k=K, capacity=1.0e9,
                  c_target=10.0, bgc_kb=None, margin=3.0)

    p5 = run_ledger(_records_with_a_borderline_point(), sat_f=0.5, **common)
    assert p5["sat_f"] == 0.5
    assert [c for c, _, _ in p5["saturated_excluded"]] == [1.0]   # 0.6 被当饱和
    assert p5["fit"]["n_pos"] == 2                                # 只剩 2 个尾点

    p7 = run_ledger(_records_with_a_borderline_point(), sat_f=0.7, **common)
    assert p7["sat_f"] == 0.7
    assert p7["saturated_excluded"] == []                          # 0.6 不再是饱和
    assert p7["fit"]["n_pos"] == 3                                 # 进拟合，且有 CI
    assert p7["fit"]["delta_ci"] is not None
    # 边界带标注：0.6 落在 0.7 的 ±20% 带内
    assert any(abs(c - 1.0) < 1e-12 for c, _ in p7["near_sat"])


def test_saturation_threshold_validation():
    with pytest.raises(ValueError, match="sat_f"):
        run_ledger(_records_one_cycle(), m_on=0.1, m_off=None, genome_mb=GENOME_MB,
                   n_win=4.0, n_plate=0.0, n_post=0.0, k=K, capacity=1.0e9,
                   c_target=10.0, bgc_kb=None, margin=3.0, sat_f=1.0)


def test_delta_ci_brackets_point_estimate():
    """δ 的 95%CI 必须满足 lo < δ < hi（2026-11-15 修过一处上下颠倒）。"""
    from windowcalc.ramp import fit_pap
    fit = fit_pap([(2, 3000, 3e8), (3, 300, 3e8), (4, 45, 3e8), (5, 3, 3e8)])
    lo, hi = fit["delta_ci"]
    assert lo < fit["delta"] < hi
    assert fit["delta"] == pytest.approx(-1.0 / fit["slope"], rel=1e-12)


# ---------------------------------------------------------------- 6. 传代拓扑（2026-11-15 新增）

def test_propagation_manual_uses_n_plate_as_is():
    from windowcalc.ramp import inter_rung_generations
    n, src = inter_rung_generations(mode="manual", n_plate=12.0)
    assert n == 12.0 and "manual" in src


def test_propagation_plate_is_log2_cells_per_colony():
    from windowcalc.ramp import inter_rung_generations
    n, src = inter_rung_generations(mode="plate", cells_per_colony=1024.0)
    assert n == pytest.approx(10.0, rel=1e-12)          # log2(1024)
    assert "plate" in src
    with pytest.raises(ValueError, match="cells-per-colony"):
        inter_rung_generations(mode="plate")
    with pytest.raises(ValueError):
        inter_rung_generations(mode="plate", cells_per_colony=1.0)


def test_propagation_liquid_is_log2_n0_over_survivors():
    from windowcalc.ramp import inter_rung_generations
    nt, k, f = 1.0e9, 12.0, 1.0e-6
    n, src = inter_rung_generations(mode="liquid", nt=nt, k=k, f=f)
    n0 = nt - k / f                        # 缺省：刚好够 k 的最大接种量
    assert n == pytest.approx(math.log2(n0 / (f * nt)), rel=1e-12)
    assert "liquid" in src
    # 显式给 N₀ 时按它算
    n2, _ = inter_rung_generations(mode="liquid", nt=nt, k=k, f=f, n0=5.0e7)
    assert n2 == pytest.approx(math.log2(5.0e7 / (f * nt)), rel=1e-12)
    # 缺 --nt 直接报错（不猜）
    with pytest.raises(ValueError, match="--nt"):
        inter_rung_generations(mode="liquid", k=k, f=f)
    # N₀ 给不出（f 太小 / k 太大）
    with pytest.raises(ValueError, match="N₀"):
        inter_rung_generations(mode="liquid", nt=1.0e6, k=12.0, f=1.0e-6)


def test_ledger_base_uses_each_cycles_own_f():
    """liquid 模式下，已完成档的档间代数按**该周期自己的 f** 算，不是一刀切。"""
    from windowcalc.ramp import inter_rung_generations
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    m_on, m_off, nt = w["m_genome"], 0.01, 1.0e9
    n_win = 4.0
    recs = [
        ("1", 1.0, "post", 1e9, int(8e8)),   # 饱和
        ("1", 2.0, "post", 1e9, 1000),       # f=1e-6
        ("1", 3.0, "post", 1e9, 100),        # f=1e-7
        ("2", 1.0, "post", 1e9, int(8e8)),   # 饱和
        ("2", 2.0, "post", 1e9, 10000),      # f=1e-5 ← 与 cycle 1 不同
        ("2", 3.0, "post", 1e9, 1000),       # f=1e-6
    ]
    plan = run_ledger(recs, m_on=m_on, m_off=m_off, genome_mb=GENOME_MB,
                      n_win=n_win, n_plate=0.0, n_post=0.0, k=K,
                      capacity=1.0e9, c_target=10.0, bgc_kb=None, margin=3.0,
                      nt=nt, propagation="liquid")
    n1, _ = inter_rung_generations(mode="liquid", nt=nt, k=K, f=1e-6)
    n2, _ = inter_rung_generations(mode="liquid", nt=nt, k=K, f=1e-5)
    assert plan["base_passengers"] == pytest.approx(
        2 * n_win * m_on + (n1 + n2) * m_off, rel=1e-12)
    # 两个周期的档间代数不同（证明不是一刀切）
    assert n1 != n2


def test_cli_liquid_mode_without_nt_refuses():
    from windowcalc.cli import cmd_multitarget
    import argparse
    args = argparse.Namespace(
        records=str(_TMP_RECORDS_PATH()), markers=str(_TMP_MARKERS_PATH()),
        m_on=None, genome_mb=GENOME_MB, m_off=0.01, n_win=4.0, n_plate=0.0, n_post=0.0,
        k=K, capacity=1.0e9, c_target=10.0, bgc_kb=None, margin=3.0,
        c_start=None, nt=None, sat_f=0.5, tau_off_h=None, td=None,
        propagation="liquid", cells_per_colony=None, n0=None,
        json=False, out=None)
    with pytest.raises(ValueError, match="--nt"):
        cmd_multitarget(args)


def _TMP_RECORDS_PATH():
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write("cycle,c,when,cfu,clones\n"
                "1,1,post,1e9,800000000\n1,2,post,1e9,1000\n1,3,post,1e9,100\n")
    return path


def _TMP_MARKERS_PATH():
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w", encoding="utf-8") as f:
        f.write("marker,mu,teff_bp\nmarkerA,1.0e-7,10\n")
    return path
