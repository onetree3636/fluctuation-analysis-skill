"""引擎（model）与优化器（optimize）测试。

回代锚点：gkab583 WGS——*S. lividans* 四个 mutator 株传一轮（20–27 代），
全强度 M7（≈3–4 突变/基因组/代）额外乘客 66/74/96 个（第 4 株丢突变只有 6 个，剂量内对照）。
静态公式 E[m]=n·M 在该锚点上应给出 ~72（24 代 × 3.0）。
"""

import math

import pytest

from windowcalc.model import (Phase, ProtocolResult, RateParam, StrainCard,
                              avg_activity, evaluate, monte_carlo, phase_generations)
from windowcalc.optimize import Decision, optimize, poisson_sf


def make_card(**over):
    d = dict(
        name="test",
        M_base=RateParam(0.004, (0.002, 0.008), "literature"),
        M_on_ref=RateParam(0.4, (0.2, 0.9), "assumed"),
        M_off=RateParam(0.012, (0.006, 0.03), "assumed"),
        mu_supply=RateParam(2.0e-8, (1.0e-8, 4.0e-8), "assumed"),
        background_false=RateParam(1.0e-9),
        tau_on_h=RateParam(2.0, (1.0, 4.0)),
        tau_off_h=RateParam(4.0, (2.0, 8.0)),
        F_doses=[{"dose_mM": 0.03, "F": 10}, {"dose_mM": 0.3, "F": 100},
                 {"dose_mM": 3.0, "F": 400}],
        ref_dose_mM=0.3,
        td_h_off=2.0, td_h_on=2.5, Nt_saturation=3.0e8,
    )
    d.update(over)
    return StrainCard(**d)


def make_phases():
    return [
        Phase("恢复预培养", "pre", duration_h=24, N0=1.0e5, Nt=1.0e8),
        Phase("诱变脉冲", "pulse", duration_h=18, dose_mM=0.3, N0=1.0e7),
        Phase("洗脱操作", "decay", duration_h=2, generations=0),
        Phase("菌落长出", "off", generations=15),
        Phase("治愈传代", "off", generations=3),
    ]


def test_avg_activity_shapes():
    p = Phase("pulse", "pulse", duration_h=10)
    assert avg_activity(p, tau_on=0.0, tau_off=4.0) == pytest.approx(1.0)
    assert avg_activity(p, tau_on=4.0, tau_off=4.0) == pytest.approx(0.8)
    d = Phase("decay", "decay", duration_h=4.0)
    assert avg_activity(d, 0, 4.0) == pytest.approx((1 - math.exp(-1.0)), rel=1e-9)
    assert avg_activity(Phase("pre", "pre"), 0, 4) == 0.0


def test_scaling_law_supply_vs_T_pulse():
    """供给 ∝ 总分裂：达到饱和后，脉冲时长不再增加 events；负荷仍随时长线性增加。"""
    card = make_card()
    # 14 h（td_on=2.5 → 5.6 代，1e7×2^5.6 > 3e8 饱和）与 18 h 均饱和
    short = evaluate(card, [Phase("脉冲", "pulse", duration_h=14, N0=1.0e7, dose_mM=0.3)])
    long_ = evaluate(card, [Phase("脉冲", "pulse", duration_h=18, N0=1.0e7, dose_mM=0.3)])
    assert long_.events == pytest.approx(short.events, rel=0.02)   # 都已饱和
    assert long_.E_passengers > short.E_passengers                 # 负荷仍上升


def test_scaling_law_load_vs_N0():
    """负荷 ∝ 每世系代数：高密度接种（N₀ 大）→ 更少乘客、相近供给。"""
    card = make_card()
    dense = evaluate(card, [Phase("脉冲", "pulse", duration_h=18, N0=3.0e7, dose_mM=0.3)])
    dilute = evaluate(card, [Phase("脉冲", "pulse", duration_h=18, N0=1.0e6, dose_mM=0.3)])
    assert dense.E_passengers < dilute.E_passengers
    assert dense.events > 0.5 * dilute.events  # 供给同量级（略低，饱和差）


def test_off_phases_leak_into_load():
    """off 相位泄漏计入负荷：同脉冲、不同 off 代数 → 负荷差 = Δn·M_off。"""
    card = make_card()
    base = [Phase("脉冲", "pulse", duration_h=6, N0=1.0e7, dose_mM=0.3)]
    r1 = evaluate(card, base)
    r2 = evaluate(card, base + [Phase("治愈", "off", generations=10)])
    assert r2.E_passengers - r1.E_passengers == pytest.approx(10 * card.M_off.value, rel=1e-9)


def test_gkab583_anchor():
    """回代 gkab583：24 代全强度（M_on=3.0）→ E[P]≈72，落在实测 66–96 邻域。"""
    card = make_card(M_on_ref=RateParam(3.0, (2.0, 4.0), "literature"),
                     F_doses=[{"dose_mM": 0.3, "F": 750.0}],   # 形状比锚点：3.0/0.004
                     ref_dose_mM=0.3)
    r = evaluate(card, [Phase("gkab583 一轮", "pulse", duration_h=0, generations=24,
                              a_bar_override=1.0, dose_mM=0.3)])
    assert r.E_passengers == pytest.approx(24 * 3.0, rel=0.01)  # 静态公式复现
    assert 60 <= r.E_passengers * 0.9 <= 100                    # 与 WGS 锚点相容


def test_mc_percentiles_ordered():
    card = make_card()
    mc = monte_carlo(card, make_phases(), dose_mM=3.0, T_pulse_h=18.0, N0=5.0e7,
                     n_samples=2000, seed=1)
    for key in ("E_passengers", "events", "p_success"):
        p5, p50, p95 = mc[key]
        assert p5 <= p50 <= p95


def test_poisson_sf():
    assert poisson_sf(0.0, 1) == 0.0
    assert poisson_sf(3.0, 1) == pytest.approx(1 - math.exp(-3.0))
    assert poisson_sf(19.6, 12) > 0.94   # k=12 @0.95 所需 λ
    assert poisson_sf(3.0, 12) < 1e-3


def test_optimizer_prefers_dense_N0_and_min_dose():
    """三旋钮：可行域内选最小剂量（降 F 即降负荷）与最大接种密度（降 n）。"""
    card = make_card()
    phases = make_phases()
    dec = Decision(k_target_clones=5, P_success_min=0.95,
                   T_pulse_range_h=(3, 24), N0_min=1.0e6, N0_max=5.0e7,
                   doses_available=[0.03, 0.3, 3.0], mc_samples=800, seed=0)
    cands, diag = optimize(card, phases, dec)
    assert diag == [] and cands
    best = cands[0]
    assert best.N0 == pytest.approx(5.0e7)          # 标度律：密度拉满
    # 供给约束下的最小可行剂量（0.03 mM F=10 供给不足）
    feas_doses = {c.dose_mM for c in cands}
    assert 0.03 not in feas_doses
    # 乘客数落在合理区间（≈n_pulse·M_on + off 项）
    assert 1.0 <= best.point.E_passengers <= 15.0


def test_optimizer_reports_infeasible():
    card = make_card(mu_supply=RateParam(1.0e-11, (5e-12, 2e-11), "assumed"))
    dec = Decision(k_target_clones=12, P_success_min=0.95,
                   T_pulse_range_h=(3, 24), N0_min=1.0e6, N0_max=5.0e7,
                   doses_available=[0.03, 0.3], mc_samples=200)
    cands, diag = optimize(card, make_phases(), dec)
    assert cands == [] and diag


def test_background_redlight():
    """背景假抗性淹没供给 → 红灯。"""
    card = make_card(background_false=RateParam(1.0e-3))
    r = evaluate(card, make_phases(), dose_mM=3.0)
    assert any("红灯" in f for f in r.flags)


def test_leak_redlight():
    card = make_card(M_off=RateParam(0.05))  # >3×M_base=0.012
    r = evaluate(card, make_phases(), dose_mM=0.3)
    assert any("红灯" in f for f in r.flags)


# ---------------------------------------------------------------- v0.4 诱导剂机制模型
def _card_with_inducer(**kw):
    from windowcalc.model import InducerSpec, M_at, RateParam, StrainCard
    card = StrainCard(
        name="t", M_base=RateParam(0.004), M_on_ref=RateParam(0.4),
        M_off=RateParam(0.012), mu_supply=RateParam(2e-8),
        tau_on_h=RateParam(2.0), tau_off_h=RateParam(4.0),
        inducer=InducerSpec(**kw))
    return card, M_at


def test_inducer_model_saturates_and_normalizes():
    """M(c) 饱和到 M_off+Δsat；activity 以参考剂量归一（ref 处 =1）。"""
    card, M_at = _card_with_inducer(amplitude=1.315, ec50_mM=0.48, n_hill=1.85)
    m_low, a_low, _ = M_at(card, 0.03)
    m_mid, a_mid, _ = M_at(card, 0.48)
    m_high, a_high, _ = M_at(card, 10.0)
    assert m_low < m_mid < m_high
    assert m_high - card.M_off.value > 0.98 * 1.315          # 10 mM 已接近饱和
    assert m_high - card.M_off.value <= 1.315 * (1 + 1e-12)  # 且永不越界
    _, a_ref, _ = M_at(card, card.ref_dose_mM)
    assert a_ref == pytest.approx(1.0, rel=1e-9)             # 参考剂量处 activity=1


def test_inducer_dose_response_is_bounded():
    """饱和意味着 M 有上限：c→大 时 M 不再增长（旧 F 表会无限外推）。"""
    card, M_at = _card_with_inducer(amplitude=1.315, ec50_mM=0.48, n_hill=1.85)
    m1, m10, m100 = M_at(card, 1.0)[0], M_at(card, 10.0)[0], M_at(card, 100.0)[0]
    assert m100 - card.M_off.value <= 1.315 * (1 + 1e-12)     # 有上限
    assert (m100 - m10) < 0.1 * (m10 - m1)                    # 增量随剂量迅速衰减


# ---------------------------------------------------------------- 无剂量轴（组成型）
def _constitutive_fixture():
    from windowcalc.model import Phase, RateParam, StrainCard
    from windowcalc.optimize import Decision
    card = StrainCard(
        name="const", M_base=RateParam(0.004), M_on_ref=RateParam(0.4),
        M_off=RateParam(0.4), mu_supply=RateParam(5e-8),
        tau_on_h=RateParam(2.0), tau_off_h=RateParam(4.0),
        td_h_on=2.5, Nt_saturation=3.0e8)
    phases = [Phase(name="生长+诱变（全程）", type="pulse", duration_h=24, N0=1.0e6)]
    dec = Decision(k_target_clones=4, P_success_min=0.95,
                   T_pulse_range_h=(3.0, 24.0), N0_min=1.0e6, N0_max=5.0e7,
                   doses_available=[], mc_samples=300)
    return card, phases, dec


def test_no_dose_axis_degrades_to_two_knobs():
    """无诱导剂：报告不显示剂量列、标注'剂量轴：无'，并给出比值固定的退化提示。"""
    from windowcalc.optimize import optimize
    from windowcalc.report import render_optimize_report
    card, phases, dec = _constitutive_fixture()
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, diag, dec)
    assert "剂量轴：无（恒定表达）" in md
    assert "剂量 (mM)" not in md
    assert "无开关窗口" in md
    assert "每事件乘客" in md and "不是** M/μ_supply" in md


def test_no_dose_axis_failure_advice_does_not_mention_inducer():
    """无剂量轴且无解时，不能建议'加诱导剂'。"""
    from windowcalc.optimize import optimize
    card, phases, dec = _constitutive_fixture()
    dec.k_target_clones = 200              # 造一个必然无解的场景
    cands, diag = optimize(card, phases, dec)
    assert not cands and diag
    text = " ".join(diag)
    assert "加诱导剂" not in text
    assert "无剂量轴" in text


def test_plasmid_loss_reduces_off_phase_load():
    """温敏质粒治愈：off 相位的泄漏随代数衰减，乘客负荷应低于不衰减。"""
    from windowcalc.model import (InducerSpec, Phase, RateParam, StrainCard,
                                  evaluate)
    card = StrainCard(
        name="ts", M_base=RateParam(0.004), M_on_ref=RateParam(0.4),
        M_off=RateParam(0.12), mu_supply=RateParam(2e-8),
        tau_on_h=RateParam(2.0), tau_off_h=RateParam(4.0),
        inducer=InducerSpec(amplitude=1.315, ec50_mM=0.48, n_hill=1.85))
    base = [Phase(name="pulse", type="pulse", duration_h=6, dose_mM=0.3, N0=1e7)]
    no_cure = base + [Phase(name="off", type="off", generations=4)]
    with_cure = base + [Phase(name="off", type="off", generations=4,
                              plasmid_loss_per_gen=0.5)]
    a = evaluate(card, no_cure)
    b = evaluate(card, with_cure)
    assert b.E_passengers < a.E_passengers          # 治愈降低泄漏负荷
    assert b.events == pytest.approx(a.events, rel=1e-12)   # 供给不受影响
    off_a = next(x for x in a.phases if x.name == "off").passengers
    off_b = next(x for x in b.phases if x.name == "off").passengers
    assert off_b == pytest.approx(off_a * 0.338, rel=0.02)   # 4 代每代减半 → 平均 33.8%


# ---------------------------------------------------------------- v0.5 保护区账本
def test_protected_hits_math():
    """f = 拷贝×长度/总 DNA；E = 乘客×f；P = 1−exp(−E)。"""
    import math as _m
    from windowcalc.model import RateParam, StrainCard, protected_hits
    card = StrainCard(
        name="t", M_base=RateParam(0.004), M_on_ref=RateParam(0.4),
        M_off=RateParam(0.012), mu_supply=RateParam(2e-8),
        tau_on_h=RateParam(2.0), tau_off_h=RateParam(4.0),
        genome_mb=9.0, plasmids=[{"name": "p", "kb": 12.7, "copies": 30}],
        protected=[{"name": "BGC", "kb": 100, "copies": 1}])
    total = 9.0 * 1000 + 12.7 * 30
    f = 100 / total
    prot, p_any = protected_hits(card, 6.8)
    assert prot[0]["f"] == pytest.approx(f, rel=1e-12)
    assert prot[0]["expected"] == pytest.approx(6.8 * f, rel=1e-12)
    assert prot[0]["prob"] == pytest.approx(1 - _m.exp(-6.8 * f), rel=1e-12)
    assert p_any == pytest.approx(1 - _m.exp(-6.8 * f), rel=1e-12)


def test_protected_absent_when_not_configured():
    from windowcalc.model import RateParam, StrainCard, protected_hits
    card = StrainCard(name="t", M_base=RateParam(0.004), M_on_ref=RateParam(0.4),
                      M_off=RateParam(0.012), mu_supply=RateParam(2e-8),
                      tau_on_h=RateParam(2.0), tau_off_h=RateParam(4.0))
    assert protected_hits(card, 6.8) == ([], 0.0)


def test_report_has_protected_sections():
    """配置保护区后，报告应出现 1b 账本与 1c 预算表。"""
    from windowcalc.cli import load_strain_card
    from windowcalc.optimize import optimize
    from windowcalc.report import render_optimize_report
    from windowcalc.cli import load_protocol
    card = load_strain_card("examples/strain_card.inducer.example.yaml")
    phases, dec = load_protocol("examples/protocol_card.dual_switch.example.yaml")
    dec.mc_samples = 200
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, diag, dec)
    assert "1b. 保护区账本" in md and "1c. 不同误伤预算" in md
    assert ("产物基因簇" in md) or ("BGC" in md)


def test_poisson_sf_semantics():
    """P0 回归：p_success 必须是 Poisson 上尾 P(≥k)，k=1 才退化为 1−exp(−λ)。"""
    from windowcalc.model import poisson_sf
    assert poisson_sf(2.0, 1) == pytest.approx(1 - math.exp(-2.0))
    assert poisson_sf(1.0, 2) == pytest.approx(1 - 2 * math.exp(-1.0))
    assert poisson_sf(0.0, 3) == 0.0                       # λ=0 时不可能有 ≥1 个事件
    assert poisson_sf(5.0, 3) > poisson_sf(5.0, 4) > poisson_sf(5.0, 5)
    assert poisson_sf(50.0, 5) > 0.9999                    # λ≫k → 几乎必然
    assert poisson_sf(0.5, 40) < 1e-10                     # k≫λ → 几乎不可能
    arr = poisson_sf(__import__("numpy").array([1.0, 2.0]), 2)
    assert arr.shape == (2,) and arr[1] > arr[0]           # ndarray 透传


def test_evaluate_p_success_uses_k():
    """P0 回归：evaluate 的 p_success 必须随 k 单调下降（此前恒为 P(≥1)≈100%）。"""
    card = make_card()
    phases = [Phase("脉冲", "pulse", duration_h=18, N0=1.0e7, dose_mM=0.3)]
    base = evaluate(card, phases)
    p1 = evaluate(card, phases, k=1).p_success
    p2 = evaluate(card, phases, k=2).p_success
    p5 = evaluate(card, phases, k=5).p_success
    assert p1 > p2 > p5
    assert p1 == pytest.approx(1 - math.exp(-base.events))


def test_monte_carlo_p_success_uses_k():
    """P0 回归：MC 的 p_success_p5 也要按 P(≥k) 算，否则稳健性门槛形同虚设。"""
    card = make_card()
    phases = make_phases()
    mc1 = monte_carlo(card, phases, 0.3, 18.0, 1.0e7, n_samples=200, k=1)
    mc5 = monte_carlo(card, phases, 0.3, 18.0, 1.0e7, n_samples=200, k=5)
    assert mc1["p_success_p5"] >= mc5["p_success_p5"]
