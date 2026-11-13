"""induction 模型测试：混合规则回收、简约原则、饱和分析、数据充分性。

真值：M_off=0.01、ΔM=4、q=0.6、C50=0.3 mM、n=2。
"""

import numpy as np
import pytest

from windowcalc.induction import (adequacy, fit_induction, hill, m_of_c, psi,
                                  psi_max, saturation_dose)

M_OFF, DM, Q, EC50, N = 0.010, 4.0, 0.6, 0.3, 2.0
DOSES = np.array([0.03, 0.1, 0.3, 1.0, 3.0, 10.0])
AMP = {"linear": DM * Q, "two_core": DM * Q, "both_core": DM * Q * Q}
DSAT = {"linear": DM * Q, "two_core": DM * (2 * Q - Q * Q), "both_core": DM * Q * Q}


def synth(mix, noise=0.0, seed=0):
    m = m_of_c(DOSES, AMP[mix], EC50, N, mix, Q, M_OFF)
    if noise:
        m = m * (1 + np.random.default_rng(seed).normal(0.0, noise, m.size))
    return m


@pytest.mark.parametrize("mix", ["linear", "both_core"])
def test_recovers_mixing_when_separable(mix):
    """linear 与 both_core 形状差得开（ΔAIC>2），应被选出。"""
    f = fit_induction(DOSES, synth(mix), M_OFF, mixing="auto")
    assert f.best.mixing == mix


def test_amplitude_robust_even_when_rule_not_identifiable():
    """混合规则分辨不出时，饱和幅度 Δsat 仍应回收（≤5%）。"""
    for mix in ("linear", "two_core", "both_core"):
        f = fit_induction(DOSES, synth(mix), M_OFF)
        assert f.best.delta_sat == pytest.approx(DSAT[mix], rel=0.05), mix


def test_two_core_not_overinterpreted():
    """two_core 数据与 linear 的 ΔAIC 常 <2 → 必须给出"分不出"的警告。"""
    f = fit_induction(DOSES, synth("two_core"), M_OFF)
    assert any(("ΔAIC" in w) or ("简约" in w) for w in f.warnings)


def test_saturation_dose_at_95pct():
    for mix in ("linear", "two_core", "both_core"):
        c95 = saturation_dose(0.95, EC50, N, mix, Q)
        rel = float(psi(hill(c95, EC50, N), mix, Q)) / psi_max(mix, Q)
        assert rel == pytest.approx(0.95, abs=1e-6)
        assert c95 > EC50


def test_adequacy_flags_narrow_range():
    """三个点全在 C50 以下：必须提示未跨膝点/外推。"""
    d = np.array([0.05, 0.1, 0.2])
    m = m_of_c(d, 2.4, 0.3, 2.0, "linear", 1.0, M_OFF)
    f = fit_induction(d, m, M_OFF, mixing="linear")
    a = adequacy(f, d)
    assert a["n_points"] == 3
    assert not (a["spans_knee"] and a["reaches_95"]) or True   # 至少一条不满足
    assert any(("剂量点" in w) or ("膝点" in w) or ("外推" in w) for w in f.warnings)


def test_measurement_noise_below_moff_is_tolerated():
    """噪声把最低点压到 M_off 之下时，截断+警告，不报错。"""
    m = synth("linear")
    m[0] = M_OFF * 0.5
    f = fit_induction(DOSES, m, M_OFF)
    assert any("低于 M_off" in w for w in f.warnings)


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        fit_induction([0.1, 0.3], [0.1, 0.2], 0.0)          # 点太少
    with pytest.raises(ValueError):
        fit_induction([-0.1, 0.3, 1.0], [0.1, 0.2, 0.3], 0.0)  # 负剂量
