"""MSS-MLE 验收测试。

主验收数据：Rosche & Foster 2000, Methods 20:4–17, Table 3（Exp.2 全样本，52 管）：
    r:    0   1   2   3   4   5   6   7   8   9
    n:   11  17  12   3   4   1   1   2   0   1
该文 Table 4 给 MSS-MLE：m = 1.11，95% CL [0.77, 1.49]。
"""

import math

import numpy as np
import pytest

from windowcalc import mss_mle
from windowcalc.mss_mle import (analyze, analyze_per_tube, ld_pmf, log_likelihood,
                                mle_m, simulate)

# R&F 2000 Table 3
RF_EXP2 = {0: 11, 1: 17, 2: 12, 3: 3, 4: 4, 5: 1, 6: 1, 7: 2, 8: 0, 9: 1}
RF_EXP2_M, RF_EXP2_CI = 1.11, (0.77, 1.49)


def test_pmf_matches_compound_poisson_first_principles():
    """递推 p1、p2 与复合泊松一阶原理严格一致。"""
    m = 1.11
    p = ld_pmf(m, 2)
    assert p[0] == pytest.approx(np.exp(-m), rel=1e-12)
    # p1 = m·e^{-m}·g1, g1 = 1/(1·2) = 1/2
    assert p[1] == pytest.approx(m * np.exp(-m) / 2.0, rel=1e-12)
    # p2 = m·e^{-m}·g2 + (m²/2)e^{-m}·g1², g2 = 1/(2·3) = 1/6
    assert p[2] == pytest.approx(
        m * np.exp(-m) / 6.0 + (m**2 / 2.0) * np.exp(-m) / 4.0, rel=1e-12
    )


@pytest.mark.parametrize(
    "m,r_max",
    [(0.5, 300), (1.11, 500), (3.0, 2000), (10.0, 5000)],
)
def test_pmf_tail_mass_within_bound(m, r_max):
    """有限 r_max 的质量和 <1 是重尾模型性质（缺失 ≈ m/(r_max+1)），不得归一化。"""
    p = ld_pmf(m, r_max)
    missing = 1.0 - p.sum()
    # 渐近界 m/(r_max+1) 是一跳近似，有限 r_max 下留 10% + 绝对项余量
    assert 0.0 <= missing <= 1.1 * m / (r_max + 1) + 1e-4


def test_rf_exp2_m_and_ci():
    """主验收：复现 R&F Exp.2 的 m̂=1.11 与 CI。

    CI 容差 0.25：本实现用 profile likelihood，R&F/FALCOR 用 Stewart 1994 的
    ln(m) 正态近似——小 m 下后者更宽（下界 0.77 vs LR 约 0.81），差异是方法差
    不是错；最终验收以 FALCOR 网页对数为准（README §4.1）。
    """
    res = analyze(RF_EXP2, nts=3.0e8)  # Nt 只作标度，m 不受影响
    assert res.m == pytest.approx(RF_EXP2_M, abs=0.03)
    assert res.m_ci[0] == pytest.approx(RF_EXP2_CI[0], abs=0.25)
    assert res.m_ci[1] == pytest.approx(RF_EXP2_CI[1], abs=0.25)
    assert res.p0 == pytest.approx(11 / 52, abs=1e-9)
    assert res.n_cultures == 52
    # p₀ 交叉核验：该文给出 m=1.55，且落在有效区间
    assert res.m_p0_valid
    assert res.m_p0 == pytest.approx(1.55, abs=0.01)


def test_simulation_recovery():
    """模拟回收：已知 m 的复合泊松数据应被收回 ±15%（C=40）。"""
    for m_true, seed in [(0.8, 1), (2.0, 2), (6.0, 3)]:
        data = simulate(m_true, n_cultures=40, seed=seed)
        res = analyze(list(data), nts=1e8)
        assert res.m == pytest.approx(m_true, rel=0.15), f"m_true={m_true}"


def test_mean_nt_scaling():
    res1 = analyze(RF_EXP2, nts=1.0e8)
    res2 = analyze(RF_EXP2, nts=4.0e8)
    assert res2.mu == pytest.approx(res1.mu / 4.0)
    assert res1.m == res2.m


def test_partial_plating_thinning():
    """部分涂板用二项稀释似然（rSalvador 式）：应收回真实 m。

    旧版"先 r/z 再 ML"对有缺口的观测直方图有偏，已弃用（模块 docstring）。
    """
    data = simulate(2.0, n_cultures=60, seed=7)  # 整管"真"计数
    full = analyze(list(data), nts=1e8)
    assert full.m == pytest.approx(2.0, rel=0.15)
    rng = np.random.default_rng(99)
    thinned = [int(rng.binomial(int(r), 0.1)) for r in data]
    part = analyze(thinned, nts=1e8 * 0.1, fraction_plated=0.1)
    assert part.m == pytest.approx(2.0, rel=0.35)
    assert part.n_cultures == 60


def test_lik_unimodal_and_mle_at_peak():
    rs, ns = mss_mle._expand_counts(RF_EXP2)
    xs = np.linspace(np.log(0.05), np.log(20), 400)
    lls = [log_likelihood(np.exp(x), rs, ns) for x in xs]
    i = int(np.argmax(lls))
    assert 0 < i < len(lls) - 1  # 单峰（非边界）
    m_hat, ll_hat = mle_m(RF_EXP2)
    assert ll_hat >= lls[i] - 1e-9  # 精修不劣于粗网格


def test_all_zero_cultures():
    res = analyze({0: 30}, nts=1e8)
    assert res.m <= 0.01
    assert res.mu_ci[1] < 1e-6


# ---------------------------------------------------------------- 部分涂板（z<1）
# 回归背景：多项式系数法重构时把 p_obs(r) 的次数误当作 r（截断到 :r+1），
# 而二项混合后 p_obs(r) 是 j_max 次多项式——p_0 被算小三倍，μ̂ 偏 −61%。
# 以下三条锁死该缺陷。

def _direct_thinned_pmf(m, z, r_max, j_max=8000):
    """独立参考：显式二项混合 Σ_j p_j·C(j,r)z^r(1−z)^(j−r)。"""
    p = ld_pmf(m, j_max)
    lf = np.array([math.lgamma(n + 1) for n in range(j_max + 2)])
    out = np.zeros(r_max + 1)
    for r in range(r_max + 1):
        jj = np.arange(r, j_max + 1)
        logC = lf[jj] - lf[r] - lf[jj - r]
        out[r] = float(np.sum(np.exp(logC) * z ** r * (1 - z) ** (jj - r) * p[jj]))
    return out


@pytest.mark.parametrize("z,r_max", [(0.1, 20), (0.2, 60)])
def test_partial_plating_pmf_matches_direct_mixture(z, r_max):
    """两种模式的 pmf 都必须等于显式二项混合（相对误差 <1e-8）。"""
    from windowcalc.mss_mle import _make_pmf_fn, _pmf_from_coeffs, _poly_coeffs
    m = 1.5
    ref = _direct_thinned_pmf(m, z, r_max)
    mean_mode = _make_pmf_fn(z, r_max)(m)
    tube_mode = _pmf_from_coeffs(_poly_coeffs(r_max, z), m, r_max)
    assert np.max(np.abs(mean_mode - ref) / np.maximum(ref, 1e-300)) < 1e-8
    assert np.max(np.abs(tube_mode - ref) / np.maximum(ref, 1e-300)) < 1e-8
    assert mean_mode[0] > 0.5 * ref[0]        # p0 被截断时会掉到 1/3
    assert 0.95 < mean_mode.sum() <= 1.0


def test_partial_plating_large_rmax_is_finite():
    """r_max=200 的二项混合不得出现 inf/nan（旧实现此处溢出，LLR 变 nan）。"""
    from windowcalc.mss_mle import _make_pmf_fn
    pm = _make_pmf_fn(0.2, 200)(3.0)
    assert np.all(np.isfinite(pm))
    assert np.all(pm >= 0)
    assert pm.sum() > 0.95


def test_partial_plating_recovers_mu_both_modes():
    """部分涂板下 μ 回收（4 次重复均值 ±25%），且恒定 Nt 时两模式一致。"""
    from windowcalc.compare import simulate_assay
    mu_true, nt, z = 1e-8, 3e8, 0.2
    ratios, same = [], []
    for seed in range(4):
        cnt = simulate_assay(mu_true, [nt] * 40, fraction_plated=z, seed=seed)
        b = analyze_per_tube(cnt, [nt] * 40, fraction_plated=z)
        a = analyze(cnt, nt, fraction_plated=z)
        ratios.append(b.mu / mu_true)
        same.append(b.m == pytest.approx(a.m, rel=1e-6))
    assert 0.75 < float(np.mean(ratios)) < 1.30, ratios
    assert all(same)


def test_negative_or_fractional_counts_rejected():
    """负 r / 非整数 r 必须报错——否则 r=-1 会把系数阵索引成末行，静默给出错值。"""
    with pytest.raises(ValueError, match="不能为负"):
        analyze([0, 1, -1, 2], nts=1e8)
    with pytest.raises(ValueError, match="整数"):
        analyze([0, 1, 1.5, 2], nts=1e8)
    with pytest.raises(ValueError, match="不能为负"):
        analyze_per_tube([0, -1, 2], [1e8, 1e8, 1e8])


def test_partial_plating_matches_rsalvador_golden():
    """与 rSalvador 1.8 `prob.LD.plating` 的逐点对拍（2026-10-28 实测值）。

    这两组值同时与 flan `dflan(..., plateff=)` 一致（差 ≤1e-16）——说明三个实现
    用的是同一个二项稀释模型，不是"另一种微菌落设计"。
    """
    from windowcalc.mss_mle import _make_pmf_fn
    golden = {
        (0.1, 1.5): [0.68129206905796125, 0.17695742593254926,
                     0.060093685068583098, 0.026430383878710644],
        (0.8, 5.0): [0.011529215046068471, 0.026682698003124485,
                     0.039437977970030451],
    }
    for (z, m), vals in golden.items():
        p = _make_pmf_fn(z, 20)(m)
        assert p[:len(vals)] == pytest.approx(vals, rel=1e-12), (z, m)
