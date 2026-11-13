"""不变量与回归测试（2026-10-29 复核修复的守护线）。

这些测试断言的是**跨模块的一致性**，不是某个函数的输出值——T1/T2 就是当初
让 P0-1/P0-2 溜过去的缺口：`evaluate` 与 `monte_carlo` 各写了一遍负荷公式，
而旧测试只检查"分位数有序"。
"""

from __future__ import annotations

import math

import pytest

from windowcalc.model import (InducerSpec, Phase, RateParam, StrainCard,
                              card_problems, evaluate, monte_carlo,
                              phase_problems, phase_terms, tornado, val)
from windowcalc.optimize import Decision, optimize
from windowcalc.report import render_optimize_report


def make_card(**kw):
    base = dict(
        name="t",
        M_base=RateParam(0.004, (0.002, 0.008), "literature"),
        M_on_ref=RateParam(0.4, (0.2, 0.9), "assumed"),
        M_off=RateParam(0.012, (0.006, 0.03), "assumed"),
        mu_supply=RateParam(2e-8, (1e-8, 4e-8), "assumed"),
        tau_on_h=RateParam(2.0, (1.0, 4.0), "assumed"),
        tau_off_h=RateParam(5.0, (0.5, 5.0), "assumed"),
        ref_dose_mM=0.3, td_h_off=2.0, td_h_on=2.5, Nt_saturation=3e8,
        F_doses=[{"dose_mM": 0.3, "F": 100}],
    )
    base.update(kw)
    return StrainCard(**base)


# ---------------------------------------------------------------- T1：MC 与点估计同口径

def test_mc_median_matches_point_estimate_with_decay_phase():
    """P0-1：monte_carlo 曾把 decay 相位的活性丢掉，中位数比点估计低 40%。

    这里用"短脉冲 + 长衰变尾"（decay 负荷占比 >30%），要求两者相差 <15%。
    """
    card = make_card(tau_off_h=RateParam(8.0, (4.0, 12.0), "assumed"))
    phases = [Phase("脉冲", "pulse", duration_h=3.0, N0=1e7, dose_mM=0.3),
              Phase("衰变尾", "decay", duration_h=8.0, N0=1e7),
              Phase("选择板", "off", generations=15)]
    pt = evaluate(card, phases, dose_mM=0.3, T_pulse_h=3.0, N0=1e7)
    mc = monte_carlo(card, phases, 0.3, 3.0, 1e7, n_samples=2000, seed=0)
    p50 = mc["E_passengers"][1]
    assert pt.E_passengers > 0
    assert abs(p50 / pt.E_passengers - 1.0) < 0.15, (
        f"MC P50={p50:.3f} 与点估计 {pt.E_passengers:.3f} 差 {abs(p50/pt.E_passengers-1):.0%}"
        "——两者口径又不一致了")


def test_decay_phase_load_is_not_just_leak():
    """decay 相位的负荷必须含 (M_on−M_off)·ā，而不只是 M_off。"""
    card = make_card()
    only_leak = [Phase("衰变尾", "decay", duration_h=8.0, generations=3)]
    r = evaluate(card, [Phase("脉冲", "pulse", duration_h=1, N0=1e7)] + only_leak,
                 dose_mM=0.3, T_pulse_h=1.0, N0=1e7)
    decay_row = [p for p in r.phases if p.type == "decay"][0]
    assert decay_row.a_bar > 0.4
    assert decay_row.M_eff > 5 * val(card.M_off)


# ---------------------------------------------------------------- T2：供给按全相位累计

def test_supply_accrues_in_decay_phase():
    """P0-2：文档公式是 Σ_i，代码曾只算 pulse。给 decay 相位 N0 后必须贡献供给。"""
    card = make_card()
    base = [Phase("脉冲", "pulse", duration_h=6.0, N0=1e7, dose_mM=0.3)]
    r1 = evaluate(card, base, dose_mM=0.3, T_pulse_h=6.0, N0=1e7)
    r2 = evaluate(card, base + [Phase("衰变尾", "decay", duration_h=6.0, N0=1e7)],
                  dose_mM=0.3, T_pulse_h=6.0, N0=1e7)
    assert r2.events > r1.events * 1.1


def test_phase_terms_is_shared_scalar_and_vector():
    """同一份 phase_terms 吃标量与数组，且两者一致（单一事实源）。"""
    import numpy as np
    card = make_card()
    ph = Phase("脉冲", "pulse", duration_h=6.0, N0=1e7)
    s = phase_terms(ph, td_on=2.5, td_off=2.0, nt_sat=3e8,
                    tau_on=2.0, tau_off=5.0, m_on=0.4, m_off=0.012)
    v = phase_terms(ph, td_on=np.full(3, 2.5), td_off=np.full(3, 2.0), nt_sat=np.full(3, 3e8),
                    tau_on=np.full(3, 2.0), tau_off=np.full(3, 5.0),
                    m_on=np.full(3, 0.4), m_off=np.full(3, 0.012))
    assert np.allclose(v["passengers"], s["passengers"])
    assert np.allclose(v["supply"], s["supply"])


# ---------------------------------------------------------------- T3：报告恒等式

def test_report_no_switch_ratio_is_not_claimed_as_M_over_mu():
    """P0-3：报告曾把"每事件乘客"标成 = M/μ_supply（差 7 个数量级）。"""
    card = make_card()
    phases = [Phase("脉冲", "pulse", duration_h=6.0, N0=1e7, dose_mM=0.3)]
    dec = Decision(k_target_clones=2, P_success_min=0.3, T_pulse_range_h=(6.0, 6.0),
                   N0_min=1e7, N0_max=1e7, doses_available=[0.3], mc_samples=200)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, [], dec)
    assert "无开关窗口" in md
    assert "不是** M/μ_supply" in md or "**不是** M/μ_supply" in md
    # 打印的比值必须等于 E/events
    ratio = cands[0].point.E_passengers / cands[0].point.events
    assert f"{ratio:.3g}" in md


# ---------------------------------------------------------------- T4/T5：坏输入必须报错

def test_fraction_plated_above_one_is_rejected():
    from windowcalc.mss_mle import analyze
    with pytest.raises(ValueError, match="不在 \\(0,1\\]"):
        analyze([0, 0, 1, 2, 3], 3e8, fraction_plated=1.5)
    with pytest.raises(ValueError):
        analyze([0, 0, 1, 2, 3], 3e8, fraction_plated=0.0)


def test_negative_doubling_time_is_caught_before_computing():
    card = make_card(td_h_on=-2.5)
    problems = card_problems(card)
    assert any("td_h_on" in p for p in problems), problems


def test_card_problems_checks_order_and_physics():
    bad = make_card(M_off=RateParam(0.5), ref_dose_mM=0.0)
    problems = card_problems(bad)
    assert any("M_off" in p and "M_on_ref" in p for p in problems)
    assert any("ref_dose_mM" in p for p in problems)


def test_phase_problems_rejects_impossible_phases():
    phases = [Phase("p", "pulse", duration_h=-1.0, plasmid_loss_per_gen=1.0)]
    problems = phase_problems(phases)
    assert any("duration_h" in p for p in problems)
    assert any("plasmid_loss_per_gen" in p for p in problems)


def test_cli_refuses_invalid_card(tmp_path):
    """端到端：负 td 必须变成"错误：卡片不合法"，而不是 events=0 的假诊断。"""
    from windowcalc.cli import main
    strain = tmp_path / "s.yaml"
    strain.write_text(
        "strain: x\n"
        "M_base: {value: 0.004}\n"
        "M_on_ref: {value: 0.4}\n"
        "M_off: {value: 0.012}\n"
        "mu_supply: {value: 2.0e-8}\n"
        "tau_on_h: {value: 2.0}\n"
        "tau_off_h: {value: 5.0}\n"
        "ref_dose_mM: 0.3\n"
        "td_h_off: 2.0\n"
        "td_h_on: -2.5\n"
        "Nt_saturation: 3.0e8\n", encoding="utf-8")
    proto = tmp_path / "p.yaml"
    proto.write_text(
        "phases:\n"
        '  - {name: p, type: "pulse", duration_h: 6, N0: 1.0e7, dose_mM: 0.3}\n'
        "decision:\n"
        "  k_target_clones: 3\n  P_success_min: 0.5\n"
        "  T_pulse_range_h: [6, 6]\n  N0_min: 1.0e7\n  N0_max: 1.0e7\n",
        encoding="utf-8")
    rc = main(["optimize", "--strain", str(strain), "--protocol", str(proto)])
    assert rc == 2


# ---------------------------------------------------------------- T6：CLI 与库同一口径

def test_cli_precision_equals_library(capsys):
    """P1-6：CLI 曾把 1σ 管数先 ceil 再乘 z²（135），与库函数（131）不一致。"""
    from windowcalc.cli import main
    from windowcalc.power import sample_size_precision_analytic
    C = sample_size_precision_analytic(2e-9, 3e8, 0.25, conf_level=0.95)
    assert C == 131
    rc = main(["power", "--mode", "precision", "--mu", "2e-9",
               "--nt-mean", "3e8", "--target-rse", "0.25"])
    assert rc == 0
    out = capsys.readouterr().out
    assert f"C = {C} 管" in out, out
    assert "135" not in out.split("（")[0]


# ---------------------------------------------------------------- T7：结构参数 dict 形式

def test_structural_params_accept_ci_dict(tmp_path):
    """P2-3：文档说"每个参数可带 value/ci"，td_h_on 写成 dict 曾直接 TypeError。"""
    from windowcalc.cli import load_strain_card
    card_yaml = tmp_path / "s.yaml"
    card_yaml.write_text(
        "strain: x\n"
        "M_base: {value: 0.004, ci: [0.002, 0.008]}\n"
        "M_on_ref: {value: 0.4, ci: [0.2, 0.9]}\n"
        "M_off: {value: 0.012, ci: [0.006, 0.03]}\n"
        "mu_supply: {value: 2.0e-8, ci: [1.0e-8, 4.0e-8]}\n"
        "tau_on_h: {value: 2.0, ci: [1.0, 4.0]}\n"
        "tau_off_h: {value: 5.0, ci: [0.5, 5.0]}\n"
        "ref_dose_mM: 0.3\n"
        "td_h_off: 2.0\n"
        "td_h_on: {value: 2.5, ci: [2.0, 3.0], status: measured}\n"
        "Nt_saturation: {value: 3.0e8, ci: [2.0e8, 4.0e8]}\n", encoding="utf-8")
    card = load_strain_card(str(card_yaml))
    assert val(card.td_h_on) == pytest.approx(2.5)
    assert card.td_h_on.ci == (2.0, 3.0)
    assert not card_problems(card)


# ---------------------------------------------------------------- T8：示例卡都能跑

def test_all_example_cards_load_and_render():
    """示例卡是文档的门面：坏一张就等于文档给了不能用的模板。"""
    import pathlib
    from windowcalc.cli import load_protocol, load_strain_card
    from windowcalc.optimize import optimize
    ex = pathlib.Path(__file__).resolve().parent.parent / "examples"
    for sc in sorted(ex.glob("strain_card*.yaml")):
        card = load_strain_card(str(sc))
        assert not card_problems(card), (sc.name, card_problems(card))
        for pc in sorted(ex.glob("protocol_card*.yaml")):
            phases, dec = load_protocol(str(pc))
            assert not phase_problems(phases), (pc.name, phase_problems(phases))
            dec.mc_samples = 200
            cands, diag = optimize(card, phases, dec)
            assert cands, (sc.name, pc.name, diag)
            md = render_optimize_report(card, phases, cands, diag, dec)
            assert "附录 A" in md and md.count("#") > 5


# ---------------------------------------------------------------- T9：selection 必须显式

def test_protected_requires_explicit_selection(tmp_path):
    """P2-9：selection 缺省 silent 是静默默认——现在直接报错。"""
    from windowcalc.cli import load_strain_card
    bad = tmp_path / "s.yaml"
    bad.write_text(
        "strain: x\nM_base: {value: 0.004}\nM_on_ref: {value: 0.4}\n"
        "M_off: {value: 0.012}\nmu_supply: {value: 2e-8}\n"
        "tau_on_h: {value: 2}\ntau_off_h: {value: 5}\nref_dose_mM: 0.3\n"
        "td_h_off: 2\ntd_h_on: 2.5\nNt_saturation: 3e8\n"
        "protected:\n  - {name: BGC, kb: 22, copies: 1}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="selection"):
        load_strain_card(str(bad))


# ---------------------------------------------------------------- T10：候选池多样性

def test_candidate_pool_is_diverse():
    """P1-10：旧实现只取全局 E 最小的 20 个 → 全挤在同一剂量，MC 等于没做。"""
    card = make_card(mu_supply=RateParam(2e-6, (1e-6, 4e-6)),  # 供给充足 → 多个剂量都可行
                     F_doses=[{"dose_mM": 0.03, "F": 10}, {"dose_mM": 0.3, "F": 100},
                              {"dose_mM": 3.0, "F": 400}])
    phases = [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=15)]
    dec = Decision(k_target_clones=3, P_success_min=0.9, T_pulse_range_h=(3.0, 24.0),
                   N0_min=1e6, N0_max=1e8, doses_available=[0.03, 0.3, 3.0],
                   mc_samples=100, top_k_mc=12)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    assert cands.meta["n_doses_in_pool"] >= 2, cands.meta
    assert cands.meta["n_pool"] > len(cands), cands.meta


# ---------------------------------------------------------------- T11：边界解

def test_boundary_solution_is_flagged():
    card = make_card()
    phases = [Phase("脉冲", "pulse", duration_h=6.0, N0=1e7, dose_mM=0.3)]
    dec = Decision(k_target_clones=2, P_success_min=0.3, T_pulse_range_h=(6.0, 9.0),
                   N0_min=1e7, N0_max=1e7, doses_available=[0.3], mc_samples=100)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, [], dec)
    assert "边界解" in md


# ---------------------------------------------------------------- T12：YAML 错误要友好

def test_bad_yaml_gives_friendly_error(tmp_path, capsys):
    """P2-4：坏 YAML 曾直接甩 yaml.scanner.ScannerError 栈。"""
    from windowcalc.cli import main
    bad = tmp_path / "bad.yaml"
    bad.write_text("strain: x\n  M_base: 1\n", encoding="utf-8")
    proto = tmp_path / "p.yaml"
    proto.write_text("phases: []\n", encoding="utf-8")
    rc = main(["optimize", "--strain", str(bad), "--protocol", str(proto)])
    assert rc == 2
    err = capsys.readouterr().err
    assert "YAML 语法错误" in err and "Traceback" not in err


# ---------------------------------------------------------------- T13：剂量轴 CI 进 MC

def test_inducer_ci_widens_parameter_interval():
    """P1-3：fit-induction 的 CI 写进 inducer.ci 后，MC 的参数区间必须变宽。"""
    base = make_card(
        inducer=InducerSpec(name="c", amplitude=1.315, ec50_mM=0.48, n_hill=1.85),
        F_doses=[])
    wide = make_card(
        inducer=InducerSpec(name="c", amplitude=1.315, ec50_mM=0.48, n_hill=1.85,
                            ci={"amplitude": (1.0, 1.7), "c50": (0.4, 0.6),
                                "n": (1.5, 2.2)}),
        F_doses=[])
    phases = [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7, dose_mM=0.3)]
    a = monte_carlo(base, phases, 0.3, 12.0, 1e7, n_samples=1500, seed=0)["E_passengers"]
    b = monte_carlo(wide, phases, 0.3, 12.0, 1e7, n_samples=1500, seed=0)["E_passengers"]
    assert (b[2] / b[0]) > 1.05 * (a[2] / a[0])


# ---------------------------------------------------------------- T14：逐管 μ 搜索域

def test_per_tube_mu_domain_keeps_every_tube_in_range():
    """P1-11：下界应由最小 N 决定、上界由最大 N 决定（旧实现 max/min 对调）。

    判据：域端点处 μ·N_min == M_MIN、μ·N_max == M_MAX——即"每一管"都落在 MSS 适用域内。
    """
    import numpy as np
    from windowcalc.mss_mle import M_MAX, M_MIN, _mu_domain
    nt = np.array([1e8, 2e8, 3e8])
    lo, hi = _mu_domain(nt)
    assert math.isclose(math.exp(lo) * nt.min(), M_MIN, rel_tol=1e-12)
    assert math.isclose(math.exp(hi) * nt.max(), M_MAX, rel_tol=1e-12)
    # 域内任一点都必须让所有管满足 M_MIN ≤ m_i ≤ M_MAX
    for x in (lo, 0.5 * (lo + hi), hi):
        m = math.exp(x) * nt
        assert m.min() >= M_MIN * (1 - 1e-9) and m.max() <= M_MAX * (1 + 1e-9)


# ---------------------------------------------------------------- T15–T19：小白视角的护栏

def test_rate_snippet_warns_per_locus():
    """N0-1：μ 是单位点速率，片段必须标 per_locus 并写"不能直接粘"。"""
    from windowcalc.mss_mle import analyze
    res = analyze([0, 0, 1, 2, 3], 3e8)
    snip = res.as_yaml_snippet("M_base")
    assert "per_locus_per_passage" in snip
    assert "不能直接粘" in snip
    assert "全基因组" in snip


def test_per_locus_magnitude_triggers_warning_and_redlight():
    """N0-1/N0-3：M 填成 1e-9 量级时要报警，而不是给出"零乘客"的推荐。"""
    from windowcalc.model import card_warnings
    card = make_card(M_base=RateParam(2.1e-9), M_on_ref=RateParam(2.1e-7),
                     M_off=RateParam(3e-9))
    assert card_warnings(card), "量级可疑却没告警"
    r = evaluate(card, [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7)], 0.3, 12.0, 1e7)
    assert any("单位点" in f and "红灯" in f for f in r.flags), r.flags


def test_report_degenerate_load_sentence():
    """N0-3：两项都是 0 时不能输出"off 项不可忽略"。"""
    card = make_card(M_on_ref=RateParam(2.1e-7), M_off=RateParam(3e-9))
    phases = [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=15)]
    res = evaluate(card, phases, 0.3, 12.0, 1e7)
    assert res.E_passengers < 1e-3
    from windowcalc.optimize import Decision, optimize
    dec = Decision(k_target_clones=2, P_success_min=0.3, T_pulse_range_h=(12.0, 12.0),
                   N0_min=1e7, N0_max=1e7, doses_available=[0.3], mc_samples=100)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, [], dec)
    assert "负荷几乎为 0" in md
    assert "off 项不可忽略" not in md


def test_nts_without_unit_prints_warning(tmp_path, capsys):
    """N0-2：--nts 未给 --nt-unit 时必须提示单位（μ 差一个体积倍）。"""
    from windowcalc.cli import main
    csv = tmp_path / "d.csv"
    csv.write_text("r,Nt\n0,3e8\n1,3e8\n0,3e8\n2,3e8\n", encoding="utf-8")
    assert main(["rate", "--input", str(csv), "--nts", "3e8"]) == 0
    out = capsys.readouterr().out
    assert "CFU/管" in out and "per_mL" in out


def test_compare_accepts_constant_nt(tmp_path, capsys):
    """N1-2：只有均值 Nt 的用户也能做两批比较。"""
    from windowcalc.cli import main
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    a.write_text("r\n0\n1\n0\n2\n0\n", encoding="utf-8")
    b.write_text("r\n0\n0\n1\n0\n3\n", encoding="utf-8")
    rc = main(["compare", "--d1", str(a), "--d2", str(b),
               "--nt1", "3e8", "--nt2", "3e8", "--labels", "WT,M7"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "fold-change" in out and "常数" in out


# ---------------------------------------------------------------- T20–T26：小白视角（第二轮）

def test_snippet_key_can_be_a_card_field():
    """N1-3：--label WT 的键名不是卡片字段，直接粘会报"缺 M_base"。"""
    from windowcalc.mss_mle import analyze
    r = analyze([0, 0, 1, 2, 3], 3e8)
    assert "不是菌株卡字段" in r.as_yaml_snippet("WT")          # 提醒改名
    assert "\nM_base:" in r.as_yaml_snippet("WT", key="M_base")     # 显式给字段
    assert "不是菌株卡字段" not in r.as_yaml_snippet("WT", key="M_base")


def test_tube_index_column_gets_specific_hint(tmp_path):
    """N1-4：多一列管号时，提示要说"管号"，而不是"频率/列顺序反了"。"""
    from windowcalc.cli import main
    csv = tmp_path / "t.csv"
    csv.write_text("管号,r,Nt\n1,0,3.0e8\n2,1,3.0e8\n", encoding="utf-8")
    assert main(["rate", "--input", str(csv)]) == 2


def test_tube_index_hint_text(tmp_path, capsys):
    from windowcalc.cli import main
    csv = tmp_path / "t.csv"
    csv.write_text("管号,r,Nt\n1,0,3.0e8\n2,1,3.0e8\n", encoding="utf-8")
    main(["rate", "--input", str(csv)])
    err = capsys.readouterr().err
    assert "管号" in err and "删掉第一列" in err


def test_fit_csv_skips_non_numeric_header(tmp_path, capsys):
    """N1-5：中文表头不能被当成数据（旧实现抛 could not convert string to float）。"""
    from windowcalc.cli import main
    csv = tmp_path / "dose.csv"
    csv.write_text("浓度,m\n0.1,1.2\n0.3,4.5\n1.0,9.0\n", encoding="utf-8")
    rc = main(["fit-induction", "--input", str(csv), "--m-off", "0.1", "--n-boot", "20"])
    assert rc == 0
    assert "剂量响应" in capsys.readouterr().out


def test_bare_run_gives_chinese_menu(capsys):
    """N2-1：裸跑不应该是 argparse 的英文报错。"""
    from windowcalc.cli import main
    assert main([]) == 0
    out = capsys.readouterr().out
    assert "你想做什么" in out and "windowcalc rate" in out


def test_report_has_reading_guide_and_anchor():
    """N2-2/N2-3/N2-4：报告要有"怎么读"、乘客量级锚点、events 的白话解释。"""
    card = make_card()
    phases = [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=15)]
    dec = Decision(k_target_clones=3, P_success_min=0.5, T_pulse_range_h=(6.0, 12.0),
                   N0_min=1e7, N0_max=1e7, doses_available=[0.3], mc_samples=100)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, [], dec)
    assert "怎么读" in md
    assert "不开 mutator" in md and "66–96" in md
    assert "不同来源" in md


def test_rate_error_hint_is_context_aware(tmp_path, capsys):
    """CSV 报错不该提"卡片写法"（那是 optimize 的文档）。"""
    from windowcalc.cli import main
    bad = tmp_path / "bad.csv"
    bad.write_text("r,Nt\nx,3e8\n", encoding="utf-8")
    main(["rate", "--input", str(bad)])
    err = capsys.readouterr().err
    assert "卡片写法" not in err
    assert "tutorial.md" in err


# ---------------------------------------------------------------- T27–T31：passenger_scale

def _per_locus_card(**kw):
    base = dict(M_base=RateParam(2.1e-9, (5.4e-10, 5.3e-9)),
                M_on_ref=RateParam(2.1e-7, (1e-7, 4e-7)),
                M_off=RateParam(3.0e-9, (1e-9, 9e-9)))
    base.update(kw)
    return make_card(**base)


def test_passenger_scale_scales_load_only():
    """passenger_scale 只乘负荷：E 线性放大，events 与红灯比值不变。"""
    phases = [Phase("脉冲", "pulse", duration_h=24.0, N0=1e7, dose_mM=0.3)]
    a = evaluate(_per_locus_card(), phases, 0.3, 24.0, 1e7)
    b = evaluate(_per_locus_card(passenger_scale=RateParam(1e6, (3e5, 3e6))),
                 phases, 0.3, 24.0, 1e7)
    assert b.E_passengers == pytest.approx(a.E_passengers * 1e6, rel=1e-9)
    assert b.events == pytest.approx(a.events, rel=1e-12)
    # 红灯判据是 M_off/M_base 的比值 → 与因子无关
    assert any("泄漏" in f for f in a.flags) == any("泄漏" in f for f in b.flags)


def test_passenger_scale_does_not_change_the_recommendation():
    """★ 尺度不变性：推荐窗口与 passenger_scale 无关（所以它只是"报告旋钮"）。"""
    phases = [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=15)]
    dec = Decision(k_target_clones=3, P_success_min=0.9, T_pulse_range_h=(3.0, 24.0),
                   N0_min=1e6, N0_max=1e8, doses_available=[0.3, 3.0], mc_samples=100)
    a, _ = optimize(_per_locus_card(F_doses=[{"dose_mM": 0.3, "F": 100},
                                             {"dose_mM": 3.0, "F": 400}]), phases, dec)
    b, _ = optimize(_per_locus_card(F_doses=[{"dose_mM": 0.3, "F": 100},
                                             {"dose_mM": 3.0, "F": 400}],
                                    passenger_scale=RateParam(1e6, (3e5, 3e6))), phases, dec)
    assert [(c.dose_mM, round(c.T_pulse_h, 6), round(c.N0, 6)) for c in a] == \
           [(c.dose_mM, round(c.T_pulse_h, 6), round(c.N0, 6)) for c in b]


def test_passenger_scale_enters_tornado_and_mc():
    """因子带 CI 时：tornado 有它一行，MC 的参数区间随之变宽。"""
    card = _per_locus_card(passenger_scale=RateParam(1e6, (3e5, 3e6)))
    phases = [Phase("脉冲", "pulse", duration_h=24.0, N0=1e7, dose_mM=0.3)]
    rows = tornado(card, phases, 0.3, 24.0, 1e7)
    assert rows and rows[0]["param"] == "passenger_scale"
    mc = monte_carlo(card, phases, 0.3, 24.0, 1e7, n_samples=600, seed=0)
    p5, _, p95 = mc["E_passengers"]
    assert p95 / p5 > 3          # 因子 CI（10 倍）+ 其它参数 → 区间明显变宽


def test_passenger_scale_must_be_positive():
    assert any("passenger_scale" in p for p in card_problems(_per_locus_card(
        passenger_scale=RateParam(0.0))))


def test_report_shows_passenger_scale_and_its_neutrality():
    card = _per_locus_card(passenger_scale=RateParam(1e6, (3e5, 3e6)))
    phases = [Phase("脉冲", "pulse", duration_h=12.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=15)]
    dec = Decision(k_target_clones=2, P_success_min=0.3, T_pulse_range_h=(12.0, 12.0),
                   N0_min=1e7, N0_max=1e7, doses_available=[0.3], mc_samples=100)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    md = render_optimize_report(card, phases, cands, [], dec)
    assert "passenger_scale" in md
    assert "窗口推荐与该因子无关" in md
    assert "本节 M_eff 与乘客贡献已乘" in md


def test_protection_cap_binds_and_reports_shortfall():
    """保护优先（硬约束）：k 要求的窗口超过保护区上限时，按上限截断并报告缺口。"""
    card = make_card(
        genome_mb=8.4,
        protected=[{"name": "BGC", "kb": 200.0, "copies": 1, "selection": "silent"}],
    )
    phases = [Phase("脉冲", "pulse", duration_h=6.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=10)]
    dec = Decision(k_target_clones=3, P_success_min=0.9,
                   T_pulse_range_h=(3.0, 24.0), N0_min=1e6, N0_max=1e8,
                   doses_available=[0.3], P_protected_max=0.012, mc_samples=50)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    m = cands.meta
    assert m["protection_binding"] is True
    assert m["protection_impossible"] is False
    assert m["actual_events"] < 3
    assert cands[0].point.p_other <= 0.012 + 1e-12
    assert m["protected_binding"] is True
    md = render_optimize_report(card, phases, cands, diag, dec)
    assert "保护上限 binding" in md
    assert "不要为了 k 超窗" in md


def test_mutator_and_bgc_caps_are_independent():
    """两条上限独立：mutator 自失活先 binding 时，不应把 BGC 也算作 binding。"""
    card = make_card(
        genome_mb=8.4,
        protected=[
            {"name": "mutator", "kb": 50.0, "copies": 1, "selection": "silent", "role": "mutator"},
            {"name": "BGC", "kb": 20.0, "copies": 1, "selection": "silent", "role": "bgc"},
        ],
    )
    phases = [Phase("脉冲", "pulse", duration_h=6.0, N0=1e7, dose_mM=0.3),
              Phase("选择板", "off", generations=10)]
    dec = Decision(k_target_clones=3, P_success_min=0.9,
                   T_pulse_range_h=(3.0, 24.0), N0_min=1e6, N0_max=1e8,
                   doses_available=[0.3], P_protected_max=0.90, P_mutator_max=0.003,
                   mc_samples=50, T_grid=6, N0_grid=5)
    cands, diag = optimize(card, phases, dec)
    assert cands, diag
    m = cands.meta
    assert m["protection_binding"] is True
    assert m["mutator_binding"] is True
    assert m["protected_binding"] is False
    assert m["protection_impossible"] is False
    assert m["actual_events"] < 3
    assert cands[0].point.p_mutator <= 0.003 + 1e-12


def test_p0_crosscheck_uses_dilution_model():
    """ζ<1 时 p0 必须解 L-D 稀释模型 G(1−ζ)=p0，而不是直接 −ln(p0)。"""
    from windowcalc.mss_mle import analyze, analyze_per_tube, p0_with_plating
    r = [0, 1, 2, 0, 3]          # 2/5 零管
    res = analyze(r, nts=3e8, fraction_plated=0.5)
    assert res.m_p0_valid is True
    assert abs(p0_with_plating(res.m_p0, 0.5) - 0.4) < 1e-4
    res2 = analyze_per_tube(r, [3e8] * 5, fraction_plated=0.5)
    assert abs(res2.m_p0 - res.m_p0) < 1e-9


def test_low_m_warns():
    """m<0.3 应在 stderr 给出“不建议”提示，而不是静默输出。"""
    from windowcalc.mss_mle import analyze
    res = analyze([0, 0, 0, 1, 0], nts=3e8)
    assert any("<0.3" in n for n in res.notes)
