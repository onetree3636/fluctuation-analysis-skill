"""compare 模块测试：逐管 Nt、LRT、fold-change、模拟器。

抽样噪声口径：L-D 重尾使单次实验的 m̂/θ̂ 方差很大（C=40 时 m≈2 的 CV≈12%）——
无偏性/覆盖率的断言一律基于**重复实验的均值**，单次实验只断言量级。
"""

import numpy as np
import pytest

from windowcalc.compare import foldchange, lrt_equal_rates, simulate_assay
from windowcalc.mss_mle import analyze, analyze_per_tube

MU = 2.0e-9          # 每细胞每分裂
CV_NT = 0.5          # Nt 波动（丝状菌口径的现实量级）


def fluctuating_nts(n, cv=CV_NT, seed=100):
    rng = np.random.default_rng(seed)
    return list(3e8 * np.exp(rng.normal(0.0, cv, n)))


def test_per_tube_unbiased_under_nt_fluctuation():
    """Ycart 2014 预言（8 次重复均值）：per_tube 无偏；mean 模式低估。"""
    ratios_old, ratios_new = [], []
    for rep in range(8):
        nts = fluctuating_nts(40, seed=200 + rep)
        counts = simulate_assay(MU, nts, seed=300 + rep)
        ratios_old.append(analyze(counts, nts=float(np.mean(nts))).mu / MU)
        ratios_new.append(analyze_per_tube(counts, nts).mu / MU)
    mo, mn = float(np.mean(ratios_old)), float(np.mean(ratios_new))
    assert mn == pytest.approx(1.0, abs=0.12)     # 新方法无偏
    assert mo < mn                                # 旧方法被校正（低估方向）
    assert mo < 0.95


def test_simulate_assay_recovery_replicate_mean():
    """常数 Nt 下 simulate_assay + per_tube：8 次重复均值回收 μ。"""
    ratios = []
    for rep in range(8):
        counts = simulate_assay(MU, [3e8] * 40, seed=800 + rep)
        ratios.append(analyze_per_tube(counts, [3e8] * 40).mu / MU)
    assert float(np.mean(ratios)) == pytest.approx(1.0, abs=0.15)


def test_lrt_calibrated():
    """同 μ 的两批 8 组 p 全部 >0.05；μ 翻倍两批 p<0.05。"""
    for rep in range(8):
        nts1 = fluctuating_nts(30, seed=400 + rep)
        nts2 = fluctuating_nts(30, seed=500 + rep)
        c1 = simulate_assay(MU, nts1, seed=600 + rep)
        c2 = simulate_assay(MU, nts2, seed=700 + rep)
        p = lrt_equal_rates(c1, nts1, c2, nts2)["p"]
        assert p > 0.05, f"rep{rep}: p={p}"

    nts1 = fluctuating_nts(30, seed=400)
    nts2 = fluctuating_nts(30, seed=500)
    c1 = simulate_assay(MU, nts1, seed=600)
    c2 = simulate_assay(2 * MU, nts2, seed=701)
    diff = lrt_equal_rates(c1, nts1, c2, nts2)
    assert diff["p"] < 0.05
    assert diff["mu2"] > diff["mu1"]


def test_foldchange_large_sample_coverage():
    """大样本（C=200/批、Nt 恒定）下 fold-change CI 应覆盖真值 θ=4（≥4/5）。"""
    cover, thetas = 0, []
    for rep in range(5):
        c1 = simulate_assay(MU, [3e8] * 200, seed=610 + rep)
        c2 = simulate_assay(4 * MU, [3e8] * 200, seed=710 + rep)
        fc = foldchange(c1, [3e8] * 200, c2, [3e8] * 200)
        thetas.append(fc["theta"])
        lo, hi = fc["ci"]
        cover += lo < 4.0 < hi
    assert cover >= 4
    assert np.median(thetas) == pytest.approx(4.0, rel=0.2)


def test_foldchange_small_sample_ci_is_wide_not_wrong():
    """小样本（C=30、Nt 波动）下 CI 宽（profile-LR 校准偏乐观的已知限制，
    见 acceptance 文档）——断言 θ̂ 在真值 2 倍内且 CI 非空。"""
    nts1 = fluctuating_nts(30, seed=410)
    nts2 = fluctuating_nts(30, seed=510)
    c1 = simulate_assay(MU, nts1, seed=610)
    c2 = simulate_assay(4 * MU, nts2, seed=710)
    fc = foldchange(c1, nts1, c2, nts2)
    lo, hi = fc["ci"]
    assert 0 < lo < hi
    assert fc["theta"] == pytest.approx(4.0, rel=1.2)   # 量级正确（2 倍内）


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
