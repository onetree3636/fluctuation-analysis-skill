"""对拍测试：**同一批公式的两套实现必须给出同一个数**。

背景（2026-11-15 代码审查的发现）：

  · `mtblocks.window_table_f_step()` / `minimal_n0()` 是文档里的"主路径积木"
    （单窗：给 f_step 与 Nt，解出最大 N₀ 与乘客）；
  · CLI 真正跑的是 `ramp.run_ledger()` → `plan_ramp()`，里面**重写了一遍**
    同样的账（N₀ = Nt − k/f、乘客 = 代数 × M、失守 = 1 − exp(−E[P]·L/G)）。

两套实现的公式现在一致，但**没有任何测试锁住**——改一处、忘另一处就会静默分叉
（历史上真发生过同类事故：BGC 分母写成 1e6、Wilson 参数顺序写反）。

本文件把三件事钉死：
  1. `minimal_n0` 的代数与账本里 N₀ 设计行（ramp.py `nt − k/f`）同源；
  2. `window_table_f_step` 的乘客 == 账本 `base_passengers`（同一 m_on、同一 n_win）；
  3. `window_table_f_step` 的失守概率 == `1 − exp(−E[P] × L/G)`，
     且梯级的累计乘客严格等于「基数 + i × 每档」，`n_post` 只加在最后一档。

改公式时：先跑本文件。它红 = 两处已经分叉。
"""

import math

import pytest

from windowcalc.mtblocks import Marker, minimal_n0, window_table_f_step
from windowcalc.ramp import run_ledger

# ---- 测试用的**通用**参数（不是任何项目实测值）-----------------------
MK = [Marker("A", 2.0e-8 * 7, 7)]      # μ_perbp = 2e-8
GENOME_MB = 10.0                        # 基因组大小（Mb）
F_STEP = 1.0e-6                         # 供给（克隆数/CFU）
K = 12.0                                # 目标克隆数
NT = 1.0e9                              # 饱和密度
BGC_KB = 100.0                          # 保护区大小


def _records_one_cycle() -> list[tuple[str, float, str, float, int]]:
    """一个周期的账本：1 个饱和点 + 2 个尾点（够 fit_pap 的 2 点下限）。"""
    return [
        ("1", 1.0, "post", 1.0e9, int(8.0e8)),   # f ≈ 0.8 → 饱和，用来定位 c_now
        ("1", 2.0, "post", 1.0e9, 1000),         # f = 1e-6
        ("1", 3.0, "post", 1.0e9, 100),          # f = 1e-7
    ]


# ---------------------------------------------------------------- 1. N₀ 同源

@pytest.mark.parametrize("f,k,nt", [
    (1e-6, 12.0, 1e9),
    (3e-7, 12.0, 2e9),
    (1e-5, 50.0, 5e8),
])
def test_minimal_n0_matches_ledger_n0_row(f, k, nt):
    """账本第三节的 N₀ 设计行写的是 `Nt − k/f`——必须与 minimal_n0 完全一致。"""
    r = minimal_n0(f, k, nt)
    assert r["feasible"]
    assert r["n0"] == pytest.approx(nt - k / f, rel=1e-12)
    # 代数口径自洽：ΔN = k/f 且 n_gen = log2(Nt/N0)
    assert r["n_gen"] == pytest.approx(math.log2(nt / r["n0"]), rel=1e-12)
    # 供给回头：rate × (Nt − N0) 必须刚好等于 k（取等号）
    assert r["events"] == pytest.approx(k, rel=1e-12)


def test_minimal_n0_generations_are_consistent_with_divisions():
    """n_gen 必须满足 ΔN = N0·(2^n − 1)——这是"代数 vs 总分裂数"的换算桥。"""
    r = minimal_n0(1e-6, 12.0, 1e9)
    n, n0 = r["n_gen"], r["n0"]
    assert n0 * (2.0 ** n - 1.0) == pytest.approx(r["delta_n"], rel=1e-9)


# ---------------------------------------------------------------- 2. 乘客同源

def test_window_table_passengers_equal_ledger_base():
    """mtblocks 的"单窗乘客" == 账本的"已完成档基数"（同一 m_on、同一 n_win）。"""
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    assert w["row"]["feasible"]

    plan = run_ledger(_records_one_cycle(),
                      m_on=w["m_genome"], m_off=None, genome_mb=GENOME_MB,
                      n_win=w["row"]["n_gen"],      # ← 用 mtblocks 解出的代数
                      n_plate=0.0, n_post=0.0,
                      k=K, capacity=1.0e9, c_target=10.0, bgc_kb=BGC_KB,
                      margin=3.0)
    assert plan["cycles"] == ["1"]
    assert plan["base_passengers"] == pytest.approx(w["row"]["passengers"], rel=1e-12)

    # 并且两边都等于 n_gen × (μ_perbp × G)——这是文档里的三段账第一项
    expect = w["row"]["n_gen"] * w["m_genome"]
    assert w["row"]["passengers"] == pytest.approx(expect, rel=1e-12)
    assert plan["base_passengers"] == pytest.approx(expect, rel=1e-12)


# ---------------------------------------------------------------- 3. 失守概率同源

def test_window_table_bgc_equals_poisson_single_shot():
    """失守 = 1 − exp(−E[P] × L/G)：单次代入，绝不做每档概率相加。"""
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT,
                            genome_mb=GENOME_MB, bgc_kb=BGC_KB)
    p = w["row"]["passengers"]
    frac = BGC_KB * 1000.0 / (GENOME_MB * 1e6)
    assert w["bgc_frac"] == pytest.approx(frac, rel=1e-12)
    assert w["row"]["p_bgc"] == pytest.approx(1.0 - math.exp(-p * frac), rel=1e-12)


def test_ladder_cumulative_matches_three_segment_ledger():
    """梯级累计 = 基数 + i ×（n_win·M_on + n_plate·M_off）；n_post 只加在最后一档。"""
    w = window_table_f_step(MK, f_step=F_STEP, k=K, nt=NT, genome_mb=GENOME_MB)
    n_win, n_plate, n_post = 4.0, 2.0, 3.0
    m_on, m_off = w["m_genome"], 0.01

    plan = run_ledger(_records_one_cycle(), m_on=m_on, m_off=m_off,
                      genome_mb=GENOME_MB, n_win=n_win, n_plate=n_plate,
                      n_post=n_post, k=K, capacity=1.0e9, c_target=10.0,
                      bgc_kb=None, margin=3.0)
    assert plan["status"] == "可行"
    base = plan["base_passengers"]
    assert base == pytest.approx(n_win * m_on + n_plate * m_off, rel=1e-12)

    per_rung = n_win * m_on + n_plate * m_off
    ladder = plan["ladder"]
    for i, r in enumerate(ladder, start=1):
        assert r["p_win"] == pytest.approx(n_win * m_on, rel=1e-12)
        assert r["p_plate"] == pytest.approx(n_plate * m_off, rel=1e-12)
        assert r["passengers"] == pytest.approx(per_rung, rel=1e-12)
        cum = base + i * per_rung
        if i == len(ladder):
            cum += n_post * m_off        # n_post 只加在**最后一档**的累计里
        assert r["cum_passengers"] == pytest.approx(cum, rel=1e-12)

    # 总乘客：基数 + 剩余梯级 + n_post（只一次）
    assert plan["p_total"] == pytest.approx(
        base + len(ladder) * per_rung + n_post * m_off, rel=1e-12)
