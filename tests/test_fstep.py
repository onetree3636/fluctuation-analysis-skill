"""fstep 子命令的单元测试——Wilson / Newcombe / 背景红灯。"""

import math

import pytest

from windowcalc.fstep import newcombe_diff, render, run_step, wilson_ci


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


def test_run_step_red_on_background():
    """f_pre ≥ 0.5×f_post → 背景淹没红灯（protocol-fstep.md §3）。"""
    r = run_step("1", [("pre", 1e9, 10), ("post", 1e9, 12)])
    assert r["red"]
    assert r["f_step"] == pytest.approx(max(0.0, r["f_post"][0] - r["f_pre"][0]))


def test_run_step_normal_and_paste_line():
    r = run_step("1", [("pre", 3e8, 0), ("post", 3e8, 12)])
    assert not r["red"]
    assert r["f_step"] == pytest.approx(4e-8, rel=1e-6)
    text = render([r])
    assert "--f-step 4e-08" in text
    assert "背景淹没" not in text


def test_run_step_requires_both_sides():
    with pytest.raises(ValueError):
        run_step("1", [("post", 3e8, 12)])
