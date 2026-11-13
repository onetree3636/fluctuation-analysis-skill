"""积木层（mtblocks）的单元测试——B 律、f_step 实测供给。

几何供给路线（window_table/target_set/first_step_events）已于 2026-11-13 随
design/fit/classify 移除；这里只测：过原点回归、解析解 N₀、f_step 主路径、提示带。
"""

import math

import pytest

from windowcalc.mtblocks import (Marker, f_ceiling_note, minimal_n0, mu_perbp,
                                 window_table_f_step)


# --------------------------------------------------------------- [B1] 过原点回归

def test_mu_perbp_recovers_slope():
    per = 2.0e-8
    mk = [Marker("A", per * 7, 7), Marker("B", per * 81, 81)]
    r = mu_perbp(mk)
    assert r["mu_perbp"] == pytest.approx(per, rel=1e-12)
    assert r["spread"] == pytest.approx(1.0, rel=1e-12)
    assert r["n"] == 2


def test_mu_perbp_flags_inconsistent_markers():
    mk = [Marker("A", 2.0e-8 * 7, 7), Marker("B", 1.0e-8 * 81, 81)]
    assert mu_perbp(mk)["spread"] == pytest.approx(2.0, rel=1e-9)


def test_single_marker_degenerates_to_division():
    r = mu_perbp([Marker("A", 1.4e-7, 7)])
    assert r["mu_perbp"] == pytest.approx(2.0e-8, rel=1e-9)
    assert r["spread"] == 1.0


def test_marker_input_validation():
    with pytest.raises(ValueError):
        Marker("A", 0.0, 7)
    with pytest.raises(ValueError):
        Marker("A", 1e-7, 0.0)
    with pytest.raises(ValueError):
        mu_perbp([])


# --------------------------------------------------------------- 解析解 N₀

def test_minimal_n0_is_closed_form():
    """解出来必须刚好够 k：供给(N0) == k。"""
    rate, k, nt = 3e-7, 12.0, 1e9
    r = minimal_n0(rate, k, nt)
    assert r["feasible"]
    assert r["events"] == pytest.approx(k, rel=1e-12)
    assert r["delta_n"] == pytest.approx(k / rate, rel=1e-12)
    assert r["n_gen"] == pytest.approx(math.log2(nt / r["n0"]))


def test_minimal_n0_reports_shortfall_when_infeasible():
    r = minimal_n0(1e-12, 12.0, 1e9)
    assert not r["feasible"]
    assert r["shortfall"] > 0
    assert r["n0"] is None


def test_minimal_n0_validation():
    with pytest.raises(ValueError):
        minimal_n0(3e-7, 12, 1e9, n0_min=2e9)
    with pytest.raises(ValueError):
        minimal_n0(3e-7, 0, 1e9)


def test_f_ceiling_note_bands():
    assert "最优档" in f_ceiling_note(10)
    assert "越过" in f_ceiling_note(78)
    assert "139" in f_ceiling_note(200)


# --------------------------------------------------------------- f_step 主路径

def test_f_step_path_has_no_target_geometry():
    """实测 f_step 路线不依赖任何靶点几何：结果里没有 target 键。"""
    mk = [Marker("rpoB", 1.5e-7, 7)]
    res = window_table_f_step(mk, f_step=3e-7, k=12, nt=1e9, genome_mb=10.5, bgc_kb=100)
    assert "target" not in res
    assert res["supply_source"] == "measured"
    assert res["row"]["feasible"]


def test_f_step_linear_in_supply():
    mk = [Marker("rpoB", 1.5e-7, 7)]
    a = window_table_f_step(mk, 1e-7, k=12, nt=2e9, genome_mb=10.5)["row"]
    b = window_table_f_step(mk, 2e-7, k=12, nt=2e9, genome_mb=10.5)["row"]
    assert b["delta_n"] == pytest.approx(a["delta_n"] / 2, rel=1e-12)


def test_f_step_infeasible_and_validation():
    mk = [Marker("rpoB", 1.5e-7, 7)]
    r = window_table_f_step(mk, 1e-15, k=12, nt=1e9, genome_mb=10.5)["row"]
    assert not r["feasible"] and r["shortfall"] > 0
    with pytest.raises(ValueError):
        window_table_f_step(mk, 0.0, k=12, nt=1e9, genome_mb=10.5)
