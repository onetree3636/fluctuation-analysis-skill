"""样本量与功效分析——rSalvador `samp.size.*` 的 Python 对应实现 + 模拟实测版。

两种口径（都围绕 μ，每细胞每分裂突变概率）：

1. **解析版**（快，常数 Nt）：Fisher 信息 I(m) = Σ_r (∂p_r/∂m)²/p_r。
   利用复合泊松恒等式 ∂p_r/∂m = (p⊛g)_r − p_r（g 为骨架 1/(k(k+1))），
   与 rSalvador `fisher.LD.plating` 的 U = Σ (p⊛h)²/p 同构（h=[−1,g]）。
   样本量 n = z²/(ψ²·N²·I(μN))，ψ = 目标相对标准误 RSE（z 由置信水平换算）。

Nt 波动一律参数化为对数正态的 **σ_ln**（不是 CV；CV→σ 用 `sigma_ln_from_cv`）。

2. **模拟版**（一般，覆盖 Nt 波动/部分涂板/小 m）：用已验证的
   `simulate_assay` + warm-start MLE 蒙特卡洛实测 CV 与 LRT 功效。
   warm-start 只改搜索路径不改最优点（似然单峰已验证）——与真实分析的 MLE 同一估计量。
"""

from __future__ import annotations

import math
from statistics import NormalDist

import numpy as np

from .mss_mle import _ll_per_tube_factory, _maximize_1d, ld_pmf
from .compare import simulate_assay

Z_95 = 1.959963984540054


def sigma_ln_from_cv(cv: float) -> float:
    """把 Nt 的真 CV（std/mean）换算成对数正态的 σ_ln：σ = √ln(1+CV²)。

    内部一律用 σ_ln（等比波动）；CLI 的 --nt-cv 走这个换算，别把 CV 直接当 σ。
    """
    return math.sqrt(math.log1p(cv * cv)) if cv > 0 else 0.0


def fisher_info_mean_nt(m: float, n_max: int = 300) -> float:
    """MSS 模型下 m 的每管期望 Fisher 信息 I(m) = Σ_r (∂p_r/∂m)²/p_r。

    恒等式 ∂p_r/∂m = (p⊛g)_r − p_r（g_k = 1/(k(k+1))，g_0=0）——
    与 rSalvador `LD.score.info` 的 h-卷积同构（phi=1）。
    """
    p = ld_pmf(m, n_max)
    g = np.concatenate([[0.0], 1.0 / (np.arange(1, n_max + 1) * (np.arange(1, n_max + 1) + 1.0))])
    pg = np.convolve(p, g)[: n_max + 1]
    score = pg / p - 1.0
    return float(np.sum(p * score ** 2))


def sample_size_precision_analytic(mu: float, nt_mean: float, target_rse: float = 0.25,
                                   conf_level: float = 0.95) -> int:
    """常数 Nt 下，把 μ̂ 的置信半宽/估计值压到 target_rse 所需的平行培养物数。

    每管对 μ 的 Fisher 信息 = N²·I(μN)（链式法则）；C 管：总信息 = C·N²·I；
    相对 RSE(μ̂) = 1/(μ·√(C·N²·I)) ⇒ C = 1/(μ²·ψ²·N²·I)。
    conf_level=0.95 时 z=1.96（对齐 rSalvador samp.size.* 的口径；conf_level=None/0 按 1σ RSE）。

    ★ 2026-10-29：**只在这里取整一次**。旧实现先返回 1σ 管数的 ceil（34→35）再由 CLI 乘 z²，
    把 95% 口径的答案抬高了 3%（135 vs 正确的 131）——testing.md 曾把 135 冻成验收值。
    """
    z = (NormalDist().inv_cdf(0.5 + conf_level / 2.0) if conf_level else 1.0)
    info = fisher_info_mean_nt(mu * nt_mean)
    return int(math.ceil((z / (mu * target_rse)) ** 2 / (nt_mean ** 2 * info)))


def _mu_hat_warm(counts, nts, mu_guess, fraction_plated: float = 1.0) -> float:
    """warm-start MLE：在真值 ±ln100 的括号内黄金分割（似然单峰，已验证）。"""
    rs = np.asarray(counts, dtype=int)
    nt_arr = np.asarray(nts, dtype=float)
    ll = _ll_per_tube_factory(rs, nt_arr, fraction_plated)
    x0 = math.log(mu_guess)
    x_hat, _ = ll.warm(x0)
    return math.exp(x_hat)


def cv_mu_simulation(mu: float, nts, fraction_plated: float = 1.0,
                     n_sim: int = 200, seed: int = 0) -> float:
    """模拟实测 μ̂ 的相对标准误（CV）。逐管 Nt 模式。"""
    vals = []
    for i in range(n_sim):
        counts = simulate_assay(mu, nts, fraction_plated=fraction_plated, seed=seed + i)
        vals.append(_mu_hat_warm(counts, nts, mu, fraction_plated))
    vals = np.array(vals)
    return float(vals.std(ddof=1) / max(vals.mean(), 1e-300))


def sample_size_precision_sim(mu: float, nt_mean: float, target_rse: float = 0.25,
                              sigma_ln_nt: float = 0.0, fraction_plated: float = 1.0,
                              n_sim: int = 150, c_max: int = 256, seed: int = 0) -> dict:
    """搜索把模拟实测 CV(μ̂) 压到 target_rse 的最小 C（倍增 + 二分收缩）。

    sigma_ln_nt：Nt 对数正态的 **σ_ln**（不是 CV；CV→σ 用 `sigma_ln_from_cv`）。
    """
    def nts_at(C, rng):
        if sigma_ln_nt > 0:
            return list(nt_mean * np.exp(rng.normal(0.0, sigma_ln_nt, C)))
        return [nt_mean] * C

    def cv_at(C, seed):
        rng = np.random.default_rng(seed)
        return cv_mu_simulation(mu, nts_at(C, rng), fraction_plated, n_sim, seed)

    C, hist, rounds = 8, [], 0
    while C <= c_max:
        rounds += 1
        v = cv_at(C, seed + C)
        hist.append((C, v))
        if v <= target_rse:
            break
        C *= 2
    else:
        return dict(C=None, history=hist, note=f"C>{c_max} 仍未达标")
    lo, hi = C // 2, C
    for _ in range(3):
        if hi - lo <= 2:
            break
        mid = (lo + hi) // 2
        rounds += 1
        v = cv_at(mid, seed + mid)
        hist.append((mid, v))
        if v <= target_rse:
            hi = mid
        else:
            lo = mid
    return dict(C=hi, history=hist, rounds=rounds,
                note=f"{rounds} 轮 × {n_sim} 次模拟/轮")


def _p_lrt_warm(c1, n1, c2, n2, mu1_g, mu2_g, fraction_plated):
    """warm 版 LRT p 值：H0 下 μ2 = μ1（pooled 括号取两真值的几何均值邻域）。"""
    rs1, rs2 = np.asarray(c1, dtype=int), np.asarray(c2, dtype=int)
    n1, n2 = np.asarray(n1, float), np.asarray(n2, float)
    ll1 = _ll_per_tube_factory(rs1, n1, fraction_plated)
    ll2 = _ll_per_tube_factory(rs2, n2, fraction_plated)

    x1, _ = ll1.warm(math.log(mu1_g))
    x2, _ = ll2.warm(math.log(mu2_g))
    mu1_h, mu2_h = math.exp(x1), math.exp(x2)

    def pooled(x):
        return ll1(x) + ll2(x)

    x_geo = 0.5 * (x1 + x2)
    lo = min(x1, x2) - math.log(100.0)
    hi = max(x1, x2) + math.log(100.0)
    _, ll0 = _maximize_local(pooled, x_geo, lo, hi)
    _, ll1_sep = ll1.warm(x1)          # 注意解包顺序：(x, ll)
    _, ll2_sep = ll2.warm(x2)
    llr = max(2.0 * ((ll1_sep + ll2_sep) - ll0), 0.0)
    return math.erfc(math.sqrt(llr / 2.0))


def _maximize_local(ll_fn, x_guess, ln_lo, ln_hi):
    lo = max(x_guess - math.log(100.0), ln_lo)
    hi = min(x_guess + math.log(100.0), ln_hi)
    if hi <= lo:
        return x_guess, -math.inf
    return _maximize_1d(ll_fn, lo, hi)


def power_lrt_sim(mu1: float, mu2: float, C: int, nt_mean: float, sigma_ln_nt: float = 0.0,
                  fraction_plated: float = 1.0, alpha: float = 0.05,
                  n_exp: int = 100, seed: int = 0) -> float:
    """LRT（compare 的速率齐性检验）在给定样本量下的功效（模拟）。

    sigma_ln_nt：Nt 对数正态的 **σ_ln**（不是 CV）。
    """
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_exp):
        if sigma_ln_nt > 0:
            n1 = list(nt_mean * np.exp(rng.normal(0.0, sigma_ln_nt, C)))
            n2 = list(nt_mean * np.exp(rng.normal(0.0, sigma_ln_nt, C)))
        else:
            n1 = [nt_mean] * C
            n2 = [nt_mean] * C
        c1 = simulate_assay(mu1, n1, fraction_plated, seed=int(rng.integers(1 << 30)))
        c2 = simulate_assay(mu2, n2, fraction_plated, seed=int(rng.integers(1 << 30)))
        if _p_lrt_warm(c1, n1, c2, n2, mu1, mu2, fraction_plated) < alpha:
            hits += 1
    return hits / n_exp


def sample_size_power_sim(mu1: float, fold: float, target_power: float = 0.9,
                          nt_mean: float = 3e8, sigma_ln_nt: float = 0.0,
                          fraction_plated: float = 1.0,
                          alpha: float = 0.05, n_exp: int = 60,
                          c_max: int = 200, seed: int = 0) -> dict:
    """搜索检出 fold 倍差异达到 target_power 的最小 C（倍增 + 二分收缩）。

    sigma_ln_nt：Nt 对数正态的 **σ_ln**（不是 CV）；fraction_plated 与实测同口径。
    """
    mu2 = mu1 * fold
    hist, C, rounds = [], 8, 0
    while C <= c_max:
        rounds += 1
        p = power_lrt_sim(mu1, mu2, C, nt_mean, sigma_ln_nt, fraction_plated=fraction_plated,
                          alpha=alpha, n_exp=n_exp, seed=seed + C)
        hist.append((C, p))
        if p >= target_power:
            break
        C *= 2
    else:
        return dict(C=None, history=hist, note=f"C>{c_max} 仍未达标")
    lo, hi = C // 2, C
    for _ in range(3):
        if hi - lo <= 2:
            break
        mid = (lo + hi) // 2
        rounds += 1
        p = power_lrt_sim(mu1, mu2, mid, nt_mean, sigma_ln_nt, fraction_plated=fraction_plated,
                          alpha=alpha, n_exp=n_exp, seed=seed + mid)
        hist.append((mid, p))
        if p >= target_power:
            hi = mid
        else:
            lo = mid
    return dict(C=hi, history=hist, rounds=rounds,
                note=f"{rounds} 轮 × {n_exp} 实验/轮")
