"""fitting 模块测试：合成数据参数回收 + 退化保护。"""

import math

import numpy as np
import pytest

from windowcalc.fitting import fit_hill_fc, fit_tau_off, hill_fc


def test_hill_synthetic_recovery():
    """已知 Hill 参数的合成剂量响应 → 三参数回收。"""
    fmax_t, ec50_t, n_t = 400.0, 0.5, 2.0
    doses = np.array([0.03, 0.1, 0.3, 1.0, 3.0])
    rng = np.random.default_rng(3)
    folds = hill_fc(doses, fmax_t, ec50_t, n_t) * (1 + rng.normal(0, 0.05, len(doses)))
    fit = fit_hill_fc(doses, folds * 0.004, m_base=0.004, n_boot=50, seed=1)
    assert fit.fmax == pytest.approx(fmax_t, rel=0.5)
    assert fit.ec50 == pytest.approx(ec50_t, rel=0.6)
    assert fit.n_exp == pytest.approx(n_t, rel=0.6)
    assert fit.r2 > 0.8
    # 内插：中间剂量处的预测应接近真曲线
    assert fit.predict(0.5) == pytest.approx(
        float(hill_fc(0.5, fmax_t, ec50_t, n_t)), rel=0.5)


def test_hill_degenerate_inputs():
    """单点 → DEGENERATE；非单调 → 警告。"""
    fit = fit_hill_fc([0.3], [0.4 * 0.004], m_base=0.004)
    assert any("DEGENERATE" in w for w in fit.warnings)
    fit2 = fit_hill_fc([0.03, 0.3, 3.0], [200 * 0.004, 300 * 0.004, 100 * 0.004],
                       m_base=0.004, n_boot=30, seed=2)
    assert any("非单调" in w for w in fit2.warnings)


def test_tau_off_synthetic_recovery():
    """已知 τ=4h 的衰变系列 → τ 回收、CI 非空、尾部负荷有值。"""
    tau_t, mu_off_t, mu_on_t = 4.0, 0.012, 0.4
    t = np.array([0.0, 2.0, 4.0, 8.0, 16.0])
    rng = np.random.default_rng(5)
    mus = mu_off_t + (mu_on_t - mu_off_t) * np.exp(-t / tau_t) \
        * (1 + rng.normal(0, 0.06, len(t)))
    fit = fit_tau_off(t, mus, mu_off=mu_off_t, td_h=2.5, m_on=mu_on_t,
                      n_boot=100, seed=6)
    assert fit.tau_off_h == pytest.approx(tau_t, rel=0.5)
    assert fit.ci is not None and fit.ci[0] < fit.ci[1]
    assert fit.tail_passengers == pytest.approx(
        (mu_on_t - mu_off_t) * tau_t / 2.5, rel=0.5)


def test_tau_off_too_few_points():
    fit = fit_tau_off([0.0, 4.0], [0.4, 0.2], mu_off=0.012)
    assert any("DEGENERATE" in w for w in fit.warnings)


def test_tau_off_long_tau_warning():
    """τ>12h → 长尾警示。"""
    t = np.array([0.0, 4.0, 8.0, 24.0, 48.0])
    mus = 0.012 + (0.4 - 0.012) * np.exp(-t / 20.0)
    fit = fit_tau_off(t, mus, mu_off=0.012, n_boot=50, seed=7)
    assert any("τ_off" in w and "长" in w for w in fit.warnings)


def test_hill_degenerate_inputs_do_not_crash():
    """三类退化数据必须返回 DEGENERATE，而不是 TypeError/ZeroDivisionError。"""
    # 全部 F<1（抑制型或 M_on 单位填错）
    f1 = fit_hill_fc([0.1, 0.3], [1e-9, 8e-10], m_base=1e-9)
    assert any("DEGENERATE" in w for w in f1.warnings)
    assert math.isnan(f1.ec50)
    # 全部剂量为 0
    f2 = fit_hill_fc([0.0, 0.0], [1e-9, 1e-9], m_base=1e-9)
    assert any("DEGENERATE" in w for w in f2.warnings)
    # 完全平坦（TSS=0）
    f3 = fit_hill_fc([0.03, 0.1, 0.3], [1e-7, 1e-7, 1e-7], m_base=1e-9, n_boot=20)
    assert math.isnan(f3.r2)
    assert any("平坦" in w for w in f3.warnings)


def test_tau_tail_uses_m_on():
    """尾部负荷的参考 μ_on：给了 m_on 就用它，否则用拟合值；1τ 版 = 63.2% 全量。"""
    t = np.array([0.0, 2.0, 4.0, 8.0, 16.0])
    mus = 0.012 + (0.4 - 0.012) * np.exp(-t / 4.0)
    a = fit_tau_off(t, mus, mu_off=0.012, td_h=2.5, n_boot=20, seed=1)
    b = fit_tau_off(t, mus, mu_off=0.012, td_h=2.5, m_on=0.8, n_boot=20, seed=1)
    assert a.tail_passengers == pytest.approx((a.mu_on_fit - 0.012) * a.tau_off_h / 2.5, rel=1e-9)
    assert b.tail_passengers == pytest.approx((0.8 - 0.012) * b.tau_off_h / 2.5, rel=1e-9)
    assert a.tail_passengers != pytest.approx(b.tail_passengers, rel=1e-6)
    assert a.tail_1tau == pytest.approx(a.tail_passengers * (1 - math.exp(-1)), rel=1e-12)
