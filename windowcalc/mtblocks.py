"""多靶点窗口的"积木层"——只借用成熟模型，不发明新理论。

**本模块只建 B 律（串联、中间态可被选择）下的窗口：供给用铺板实测的 f_step，
μ_perbp 只进乘客（M = μ_perbp × G）。**

2026-11-13 收缩：几何供给路线（`window_table` / `target_set` / `first_step_events`，
即 μ_perbp × T_set × entry_fraction）随 `--mode design/fit/classify` 一并移除——
多靶点/耐受药的 T_set 不可枚举，属假精度（CHANGELOG 2026-11-12）。C 律（谷型）
更早（2026-11-12）移除，公式与推导记录在 docs/multitarget-model.md §2。

积木与出处：

  [B1] 过原点回归  `μ_i = μ_perbp · T_i`            → `mu_perbp()`
       最小二乘；把 ≥2 个单靶点标定折成同一个"每碱基突变率"。
       文献原型：Lang & Murray 2008（逐标记 ÷T_eff）；均匀性边界：Jee 2016（位点间 ≥10×）。

  [B4] 相位积分——**不重写**，乘客仍用 `model` 的 `代数 n × M` 口径。

  [B5] 突变负荷上限（只作提示）                      → `f_ceiling_note()`
       Sprouffske 2018, PLoS Genet：22× 最优、139× 显著变差。

  [B6] 供给        `f_step`（铺板实测，克隆数/CFU）  → `window_table_f_step()`
       f_step 一次性吸收「靶点个数 × 长度 × 有益比例 × μ_perbp」；
       `minimal_n0()` 解析解 `N0 ≤ Nt − k/rate`。

明确不做（都不在本项目协议里出现）：
  - 几何供给路线（--targets / T_set / entry_fraction）——2026-11-13 随三个模式移除
  - C 律 / 谷型（2026-11-12 移除）；选择富集（发生在涂板/传代之后）；
    克隆干扰（`S = N·μ ≫ 1` 的判定在 `rate`/`power` 侧，属于校准实验设计）
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence


# --------------------------------------------------------------------------- 标定

@dataclass(frozen=True)
class Marker:
    """一个单靶点标定标记（一次波动实验的结果）。

    mu       : 每 locus 每代速率（`rate` 子命令的输出，per_locus_per_passage）
    teff_bp  : 该标记的**有效靶标**（能产生该表型的碱基当量，**不是基因长度**）
    f        : 可选，该标记的 F = μ_mut/μ_WT（用于跨标记一致性检查）
    """
    name: str
    mu: float
    teff_bp: float
    f: float | None = None

    def __post_init__(self):
        if not (self.mu > 0):
            raise ValueError(f"标记 {self.name!r} 的 mu={self.mu!r} 必须 > 0")
        if not (self.teff_bp > 0):
            raise ValueError(
                f"标记 {self.name!r} 的 teff_bp={self.teff_bp!r} 必须 > 0——"
                "T_eff 是\"能产生该表型的碱基当量\"，不是基因长度；"
                "填不出来的标记不能用来标定（见 references/markers.md §1）")
        if self.f is not None and self.f <= 0:
            raise ValueError(f"标记 {self.name!r} 的 f 必须 > 0")


def mu_perbp(markers: Sequence[Marker]) -> dict:
    """[B1] 过原点回归：μ_i = μ_perbp · T_i  ⇒  μ̂ = Σ(T_i μ_i) / Σ(T_i²)。

    为什么不是"对 μ 取平均"：T_eff 因标记而异（7 bp 与 81 bp 的 μ 不可平均）。
    产出的是**标记位点平均**的每碱基率——当基因组平均用的前提是各标记 λ̄≈1
    （未测；L&M 2008 的两个报告基因 λ=0.32/0.55）。误差棒见 docs/multitarget-model.md §8.4。

    返回：
        mu_perbp   斜率
        markers    逐标记的隐含 perbp 与相对斜率的比
        spread     max(比)/min(比)；>2 说明"每碱基率均匀"这条假设在你的菌里不成立
    """
    if not markers:
        raise ValueError("至少需要 1 个标定标记（0 个时不能外推，也不该用默认值顶替）")
    slope = (sum(m.teff_bp * m.mu for m in markers)
             / sum(m.teff_bp ** 2 for m in markers))
    rows = []
    for m in markers:
        implied = m.mu / m.teff_bp
        rows.append(dict(name=m.name, mu=m.mu, teff_bp=m.teff_bp,
                         perbp=implied, ratio=implied / slope))
    ratios = [r["ratio"] for r in rows]
    return dict(mu_perbp=float(slope), markers=rows,
                spread=float(max(ratios) / min(ratios)), n=len(rows),
                teff_total=float(sum(m.teff_bp for m in markers)))


# --------------------------------------------------------------------------- 窗口

def minimal_n0(rate: float, k: float, nt: float,
               n0_min: float = 1.0) -> dict:
    """给定每代供给速率（实测 f_step），解出**最大**的 N₀（= 最少代数 = 最少乘客）。

    一阶供给是解析的：`rate · (Nt − N0) ≥ k` ⇒ `N0 ≤ Nt − k/rate`。
    取等号即"刚好够 k"，所以不需要数值搜索。
    """
    if not (nt > n0_min > 0):
        raise ValueError(f"需要 Nt > n0_min > 0（当前 {n0_min!r}, {nt!r}）")
    if k <= 0:
        raise ValueError("k 必须 > 0")
    if rate <= 0:
        raise ValueError("供给速率必须 > 0")
    need = float(k) / float(rate)
    n0 = float(nt) - need
    if n0 < n0_min:
        ev_max = float(rate) * (float(nt) - float(n0_min))
        return dict(feasible=False, n0=None, events=ev_max, n_gen=None,
                    delta_n=None, shortfall=float(k - ev_max))
    return dict(feasible=True, n0=n0, events=float(k),
                n_gen=float(math.log2(nt / n0)), delta_n=need, shortfall=0.0)


def window_table_f_step(markers: Sequence[Marker], f_step: float,
                        k: float, nt: float, genome_mb: float,
                        bgc_kb: float | None = None, n0_min: float = 1.0) -> dict:
    """**主路径**：供给用实测 f_step，乘客用标定链的 M。

    `f_step` = 每细胞每代"变得能在下一档浓度生长"的概率（铺板测：克隆数/CFU）。
    它一次性吸收了 `T_set × 有益比例 × μ_perbp`，所以**不需要**知道靶点个数与长度。

    标定标记在这里只用于一条支路：`M = μ_perbp × G`（乘客预算）。
    """
    if not (f_step > 0):
        raise ValueError(f"f_step={f_step!r} 必须 > 0（克隆数/CFU，铺板实测）")
    if not (genome_mb > 0):
        raise ValueError("genome_mb 必须 > 0（用来把 μ_perbp 折成全基因组负荷）")
    cal = mu_perbp(markers)
    m_genome = cal["mu_perbp"] * genome_mb * 1e6
    bgc_frac = (bgc_kb * 1000.0) / (genome_mb * 1e6) if bgc_kb else None

    r = minimal_n0(f_step, k, nt, n0_min=n0_min)
    if r["feasible"]:
        r["passengers"] = r["n_gen"] * m_genome
        r["p_bgc"] = (None if bgc_frac is None
                      else 1.0 - math.exp(-r["passengers"] * bgc_frac))
    else:
        r["passengers"] = None
        r["p_bgc"] = None

    return dict(calibration=cal, f_step=float(f_step), m_genome=m_genome,
                k=float(k), nt=float(nt), genome_mb=float(genome_mb),
                bgc_kb=bgc_kb, bgc_frac=bgc_frac, row=r, supply_source="measured")


# --------------------------------------------------------------------------- 提示

def f_ceiling_note(fold: float) -> str:
    """[B5] 只作提示：Sprouffske 2018 的 22× 最优 / 139× 变差，不是硬阈值。"""
    if fold <= 0:
        raise ValueError("fold 必须 > 0")
    if fold <= 25:
        return f"F={fold:g}× 在 Sprouffske 2018 的\"最优档\"附近（22×），长期代价小。"
    if fold < 100:
        return (f"F={fold:g}× 已越过 22× 的收益峰值、未到 139× 的显著变差档——"
                "短期筛克隆有利，长期维持产量需监控。")
    return (f"F={fold:g}× ≥139×：Sprouffske 2018 里该档适应度显著最低（p=0.02）、"
            "90 种新环境全败；只适合极短脉冲。")
