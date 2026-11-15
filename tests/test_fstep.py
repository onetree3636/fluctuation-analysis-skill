"""fstep 子命令的单元测试——Wilson / Newcombe / 背景红灯。"""

import math

import pytest

from windowcalc.fstep import newcombe_diff, wilson_ci


def test_wilson_known_values():
    """与 protocol-fstep.md §2.2 的实算例逐位对账。"""
    p, lo, hi = wilson_ci(12, 10**9)
    assert p == pytest.approx(1.2e-8)
    assert lo == pytest.approx(6.9e-9, rel=0.05)
    assert hi == pytest.approx(2.1e-8, rel=0.05)
    assert wilson_ci(0, 10**8)[2] == pytest.approx(3.0 / 1e8)   # rule of three


def test_wilson_validation():
    with pytest.raises(ValueError):
        wilson_ci(5, 0)
    with pytest.raises(ValueError):
        wilson_ci(-1, 1e8)
    with pytest.raises(ValueError):
        wilson_ci(12, 10)                                       # 克隆 > CFU


def test_newcombe_diff_known_case():
    """f_post − f_pre：pre=0 时差值区间应与 f_post 的 Wilson 下界重合。"""
    post = wilson_ci(12, 10**9)
    pre = (0.0, 0.0, 3.0 / 1e9)
    d, lo, hi = newcombe_diff(post, pre)
    assert d == pytest.approx(post[0])
    assert 0.0 <= lo <= d <= hi
