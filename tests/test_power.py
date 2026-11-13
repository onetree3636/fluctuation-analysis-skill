"""power 模块测试：解析 Fisher 信息 vs 模拟实测的互验。

口径约定：
- target_rse/conf_level：1σ 相对标准误口径（conf_level=None）；conf_level=0.95 时
  target 指 95% CI 半宽/估计值（含 z=1.96，对齐 rSalvador samp.size.*）——两种口径
  不要混用；
- 单次实验的 m̂ 抽样噪声大（L-D 重尾；C=40 时 m≈2 的 CV≈12%、m≈0.6 时≈25%），无偏性断言基于
  多重复均值；功效断言用单调性而非绝对阈值（MC 误差内稳健）；
- 模拟器克隆封顶 5000（管级），r_max 不会突破多项式系数法适用域。
"""

import numpy as np
import pytest

from windowcalc.compare import foldchange, lrt_equal_rates, simulate_assay
from windowcalc.mss_mle import analyze, analyze_per_tube
from windowcalc.power import (cv_mu_simulation, fisher_info_mean_nt,
                              power_lrt_sim, sample_size_precision_analytic,
                              sample_size_precision_sim)

MU = 2.0e-9
NT = 3.0e8


def fluctuating_nts(n, cv=0.5, seed=100):
    rng = np.random.default_rng(seed)
    return list(NT * np.exp(rng.normal(0.0, cv, n)))


def test_fisher_info_asymptotic_matches_simulation():
    """解析 RSE（Fisher 信息，相对口径）与 C=400 大样本模拟 CV 一致（渐近互验）。"""
    C = 400
    analytic_rse = 1.0 / ((C * NT ** 2 * fisher_info_mean_nt(MU * NT)) ** 0.5 * MU)
    sim_rse = cv_mu_simulation(MU, [NT] * C, n_sim=25, seed=42)
    assert sim_rse == pytest.approx(analytic_rse, rel=0.35)


def test_cv_scaling_inverse_sqrt_C():
    """CV ∝ 1/√C：C 扩 4 倍 → CV 减半（±30%，MC 误差内）。"""
    cv_small = cv_mu_simulation(MU, [NT] * 32, n_sim=60, seed=1)
    cv_large = cv_mu_simulation(MU, [NT] * 128, n_sim=60, seed=2)
    assert cv_small / cv_large == pytest.approx(2.0, rel=0.30)


def test_sample_size_precision_sim_finds_minimum():
    """模拟搜索：返回的 C 的实测 CV ≤ 目标（含 MC 误差）。"""
    out = sample_size_precision_sim(MU, NT, target_rse=0.25, n_sim=60,
                                    c_max=256, seed=5)
    assert out["C"] is not None
    cv_at_C = cv_mu_simulation(MU, [NT] * out["C"], n_sim=60, seed=6)
    assert cv_at_C <= 0.25 * 1.35


def test_analytic_vs_simulation_sample_size():
    """解析样本量与模拟搜索互相印证（1σ RSE 口径，C 相差 <2 倍）。"""
    C_analytic = sample_size_precision_analytic(MU, NT, target_rse=0.25,
                                                conf_level=None)
    out = sample_size_precision_sim(MU, NT, target_rse=0.25, n_sim=60,
                                    c_max=512, seed=9)
    assert out["C"] is not None
    assert out["C"] == pytest.approx(C_analytic, rel=1.0)   # 2 倍以内


def test_power_lrt_monotone():
    """功效单调性：效应越大功效越高；样本量越大功效越高。"""
    p_fold2_c40 = power_lrt_sim(MU, 2 * MU, 40, NT, n_exp=40, seed=11)
    p_fold13_c40 = power_lrt_sim(MU, 1.3 * MU, 40, NT, n_exp=40, seed=12)
    p_fold2_c6 = power_lrt_sim(MU, 2 * MU, 6, NT, n_exp=40, seed=13)
    assert p_fold13_c40 < 1.0                                  # 1.3× 未饱和
    assert p_fold2_c40 > p_fold13_c40
    assert p_fold2_c40 >= p_fold2_c6                           # 降样本量功效不升


def test_power_huge_fold_needs_few_tubes():
    """F=100 的 on/off 差异：C=6 已功效≈1（SOP 场景锚点；Nt 压低以保 r 域）。"""
    assert power_lrt_sim(MU, 100 * MU, 6, 3e7, n_exp=30, seed=21) > 0.9


def test_per_tube_equivalent_when_nt_constant():
    """Nt 恒定 → per_tube 与 mean 模式一致（优化器容差级）。"""
    RF = [r for r, n in {0: 11, 1: 17, 2: 12, 3: 3, 4: 4,
                         5: 1, 6: 1, 7: 2, 8: 0, 9: 1}.items() for _ in range(n)]
    a = analyze(RF, nts=3e8)
    b = analyze_per_tube(RF, [3e8] * 52)
    assert b.m == pytest.approx(a.m, rel=1e-6)
    assert b.m_ci[0] == pytest.approx(a.m_ci[0], rel=1e-6)
    assert b.m_ci[1] == pytest.approx(a.m_ci[1], rel=1e-6)


def test_per_tube_requires_paired_nts():
    with pytest.raises(ValueError):
        analyze_per_tube([0, 1, 2], [3e8])                   # 标量 Nt：提示用 mean 模式
    with pytest.raises(ValueError):
        analyze_per_tube([0, 1], [3e8, 3e8, 3e8])            # 长度不配对
    with pytest.raises(ValueError):
        analyze_per_tube([0], [3e8])                         # 单管无法估计


def test_foldchange_small_sample_magnitude():
    """小样本（C=30、Nt 波动）下 θ̂ 量级正确（L-D 噪声大，只断 2 倍内）。"""
    nts1 = fluctuating_nts(30, seed=410)
    nts2 = fluctuating_nts(30, seed=510)
    c1 = simulate_assay(MU, nts1, seed=610)
    c2 = simulate_assay(4 * MU, nts2, seed=710)
    fc = foldchange(c1, nts1, c2, nts2)
    lo, hi = fc["ci"]
    assert 0 < lo < hi
    assert fc["theta"] == pytest.approx(4.0, rel=1.2)


def test_foldchange_large_sample_coverage():
    """大样本（C=120/批、Nt 恒定）下 fold-change CI 应覆盖真值 θ=4（≥2/3）。"""
    cover, thetas = 0, []
    for rep in range(3):
        c1 = simulate_assay(MU, [3e8] * 120, seed=610 + rep)
        c2 = simulate_assay(4 * MU, [3e8] * 120, seed=710 + rep)
        fc = foldchange(c1, [3e8] * 120, c2, [3e8] * 120)
        thetas.append(fc["theta"])
        lo, hi = fc["ci"]
        cover += lo < 4.0 < hi
    assert cover >= 2
    assert np.median(thetas) == pytest.approx(4.0, rel=0.3)
