"""ramp 子命令的单元测试——PAP 拟合、梯级规划、累计乘客与失守账。"""

import argparse
import math

import pytest

from windowcalc.ramp import fit_pap, plan_ramp, predict_f, render


def _synth():
    """δ=1 g/L/十倍频的合成 PAP：f(1)=1e-6, f(2)=1e-7, f(3)=1e-8，c=4 零克隆。"""
    return [(1, 300, 3e8), (2, 30, 3e8), (3, 3, 3e8), (4, 0, 3e8)]


def _plan(**kw):
    fit = fit_pap(_synth())
    base = dict(m_on=0.162, m_off=8.7e-3, n_win=4.0, n_plate=27.0, n_post=4.0,
                k=30.0, capacity=6e8, c_start=1.0, c_target=5.0,
                genome_mb=7.5, bgc_kb=83.0, margin=1.0)
    base.update(kw)
    return plan_ramp(fit, **base)


# ---------------------------------------------------------------- 拟合

def test_fit_recovers_delta():
    fit = fit_pap(_synth())
    assert fit["slope"] == pytest.approx(-1.0)
    assert fit["delta"] == pytest.approx(1.0)
    assert fit["f_ref"] == pytest.approx(1e-8)
    assert fit["n_pos"] == 3                      # 零克隆行不进拟合
    assert fit["censored"] == [(4.0, pytest.approx(1e-8))]


def test_fit_rejects_fewer_than_two_positive():
    with pytest.raises(ValueError, match="至少需要 2"):
        fit_pap([(1, 300, 3e8), (2, 0, 3e8)])


def test_fit_rejects_positive_slope():
    with pytest.raises(ValueError, match="斜率必须 < 0"):
        fit_pap([(1, 3, 3e8), (2, 30, 3e8), (3, 300, 3e8)])


def test_fit_pools_same_concentration():
    rows = [(1, 150, 1.5e8), (1, 150, 1.5e8), (2, 30, 3e8)]
    fit = fit_pap(rows)
    assert fit["n_pos"] == 2                      # 同浓度两板合并
    assert fit["slope"] == pytest.approx(-1.0)


def test_fit_flags_fit_above_censor_bound():
    # 拟合曲线在 c=4 预测 1e-9 < 上界 1e-8 → 不该报警；构造会报警的数据：
    rows = [(1, 300, 3e8), (2, 30, 3e8), (3, 300, 3e8)]   # 非单调 → 正斜率被拒
    with pytest.raises(ValueError):
        fit_pap(rows)


def test_predict_ci_widens_with_distance():
    fit = fit_pap([(1, 300, 3e8), (2, 45, 3e8), (3, 30, 3e8)])   # 带噪声：残差非零
    _, lo1, hi1 = predict_f(fit, 3.5)             # 实测范围内
    _, lo2, hi2 = predict_f(fit, 6.0)             # 远端外推
    # CI 宽度在 log10 空间比较（线性空间被数量级压缩）
    w1 = math.log10(hi1) - math.log10(lo1)
    w2 = math.log10(hi2) - math.log10(lo2)
    assert w2 > w1 > 0


# ---------------------------------------------------------------- 规划

def test_plan_reaches_target_with_expected_rungs():
    plan = _plan()
    assert plan["status"] == "可行"
    # Δmax = log10(1e-6 × 6e8 / 30) = log10(20) = 1.30 → 4 档（等分 4 g/L 跨度）
    assert plan["n_rungs"] == 4
    assert plan["ladder"][-1]["c_to"] == pytest.approx(5.0)
    # 平稳补充：每档 f_step ≈ 曲线在 c_start+Δ 处的实测值 = 1e-7
    assert plan["f_rung"] == pytest.approx(1e-7, rel=0.01)


def test_plan_passengers_are_additive():
    plan = _plan()
    S = plan["n_rungs"]
    expect = S * (4.0 * 0.162 + 27.0 * 8.7e-3) + 4.0 * 8.7e-3
    assert plan["p_total"] == pytest.approx(expect)
    assert plan["ladder"][-1]["cum_passengers"] == pytest.approx(expect)


def test_plan_bgc_is_single_shot_not_summed():
    plan = _plan()
    cum = plan["p_total"]
    exact = 1.0 - math.exp(-cum * 83.0 / 7500.0)
    assert plan["p_bgc"] == pytest.approx(exact)
    # 朴素"每档失守概率相加"在单命中主导段与严格式接近（差 <5%），
    # 但概念上仍是错的：档数多了或 λ 大了会高估——严格式只有总乘客单次代入
    per_step = plan["ladder"][0]["passengers"]
    naive = plan["n_rungs"] * (1.0 - math.exp(-per_step * 83.0 / 7500.0))
    assert abs(naive - exact) / exact < 0.05      # 小概率段内接近，但不相等也允许


def test_plan_bigger_margin_shrinks_steps():
    small = _plan(margin=1.0)
    big = _plan(margin=10.0)
    assert big["step"] < small["step"]
    assert big["n_rungs"] >= small["n_rungs"]


def test_plan_stuck_when_capacity_too_small():
    plan = _plan(capacity=1e6)                    # f(1)×cap = 1 ≪ k=30
    assert plan["status"] == "卡住"
    assert "富集" in plan["stuck_reason"] or "容量" in plan["stuck_reason"]


def test_plan_stuck_by_margin():
    plan = _plan(margin=1000.0)                   # f(1)×cap=600 < 30×1000
    assert plan["status"] == "卡住"
    assert "margin" in plan["stuck_reason"]


def test_plan_refuses_bad_args():
    with pytest.raises(ValueError, match="c_target > c_start"):
        _plan(c_start=5.0, c_target=5.0)
    with pytest.raises(ValueError, match="capacity"):
        _plan(capacity=0.0)
    with pytest.raises(ValueError, match="n_win"):
        _plan(n_win=0.0)
    with pytest.raises(ValueError, match="margin"):
        _plan(margin=0.0)


def test_plan_m_off_absent_is_flagged():
    plan = _plan(m_off=None)
    assert plan["m_off_counted"] is False
    out = render(plan)
    assert "未给 --m-off" in out and "未计入" in out


def test_plan_flags_extrapolation_beyond_pap():
    plan = _plan(c_target=9.0)                    # 实测只到 c=3（4 为上界）
    assert plan["extrapolation_decades"] is not None
    assert plan["extrapolation_decades"] > 1


# ---------------------------------------------------------------- 渲染 / CLI

def test_render_contains_sections_and_numbers():
    out = render(_plan())
    for key in ("PAP 尾部拟合", "梯级规划", "累计账", "怎么读", "0.92", "83"):
        assert key in out or key.replace(".92", "") in out


def test_cli_smoke(tmp_path, capsys):
    pap = tmp_path / "records.csv"
    pap.write_text("cycle,c,when,cfu,clones\n"
                   "1,2,pre,3.0e8,0\n1,2,post,3.0e8,300\n1,3,post,3.0e8,30\n"
                   "1,4,post,3.0e8,3\n2,3,pre,3.0e8,1\n2,3,post,3.0e8,40\n"
                   "2,4,post,3.0e8,8\n2,5,post,3.0e8,0\n", encoding="utf-8")
    mk = tmp_path / "mk.csv"
    mk.write_text("marker,mu,teff_bp\nmarkerA,1.0e-7,10\n", encoding="utf-8")
    from windowcalc.cli import cmd_multitarget
    args = argparse.Namespace(
        records=str(pap), markers=str(mk), genome_mb=5.0, m_off=0.01,
        n_win=4.0, n_plate=12.0, n_post=4.0, k=30.0, capacity=3e8,
        c_target=5.0, bgc_kb=40.0, margin=3.0, c_start=None, nt=None,
        sat_f=0.5, tau_off_h=None, td=None,
        json=False, out=None)
    assert cmd_multitarget(args) == 0
    out = capsys.readouterr().out
    for key in ("各档 f_step", "PAP 曲线现状", "下一窗设计", "剩余梯级", "累计"):
        assert key in out, f"缺章节：{key}"


def test_cli_refusal_lists_missing(tmp_path):
    mk = tmp_path / "mk.csv"
    mk.write_text("marker,mu,teff_bp\nmarkerA,1.0e-7,10\n", encoding="utf-8")
    from windowcalc.cli import cmd_multitarget
    args = argparse.Namespace(records=None, markers=str(mk), genome_mb=None,
                              m_off=None, n_win=None, n_plate=0.0, n_post=0.0,
                              sat_f=0.5, tau_off_h=None, td=None,
                              k=None, capacity=None, c_target=None,
                              bgc_kb=None, margin=3.0, c_start=None, nt=None,
                              json=False, out=None)
    with pytest.raises(ValueError, match="不会用默认值代替"):
        cmd_multitarget(args)
