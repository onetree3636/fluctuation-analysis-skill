"""两批波动实验的比较与模拟——rSalvador LRT.LD/newton.foldchange 与 flan 模拟器的 Python 对应实现。

- `simulate_assay`：已知 μ、逐管 Nt、涂板比例 z 的波动实验模拟（事件 Poisson(μ·N_i)
  → LD 骨架克隆 → 二项稀释）。
- `lrt_equal_rates`：两批数据速率齐性的似然比检验（H0: μ1=μ2，LLR ~ χ²₁）。
- `foldchange`：μ2/μ1 的 profile-likelihood 点估计与 95% CI（F = M_on/M_base 的直接需求）。
- fitness/death 的完整模型不在此实现（需生成函数数值反演，见 flan）——该场景请用
  zrway fluctuation 环境的 flan 复核。
"""

from __future__ import annotations

import math

import numpy as np

from .mss_mle import _ll_per_tube_factory, _maximize_1d, _profile_ci_1d

CHI2_1_95_HALF = 3.841458820694124 / 2.0


def _skeleton_clone_sum(k: int, rng, cap: int | None = None) -> int:
    """k 个突变事件的克隆大小合计（骨架 P(size)=1/(size(size+1))，逆变换采样）。

    cap：单个克隆的大小上限。理想 LD 骨架均值无穷（P(size≥j)≈1/j），偶发天文数字
    的克隆会拖垮下游多项式系数（O(r_max²)）——真实有限群体中克隆不可能超过
    终末细胞数，故默认 cap=N_i（生物学正确的截断）。
    """
    if k <= 0:
        return 0
    u = rng.uniform(1e-12, 1.0, size=k)
    sizes = np.ceil(1.0 / u - 1.0).astype(np.int64)
    if cap is not None:
        sizes = np.minimum(sizes, cap)
    return int(sizes.sum())


def simulate_assay(mu: float, nts, fraction_plated: float = 1.0, seed: int = 0) -> list:
    """模拟一份波动实验：逐管事件 Poisson(μ·N_i)，克隆按 LD 骨架生长（封顶 5000），
    涂板二项稀释。返回逐管 r（观测突变体数）。

    封顶 5000 的理由：理想骨架 P(size≥j)≈1/j 无上界，偶发大克隆会把分析端
    `_poly_coeffs` 顶到 O(r²) 爆内存；截断质量 ~K/5000 对 μ 估计的偏倚可忽略
    （真实实验中克隆本就被终末细胞数 N_i 自然截断）。
    """
    rng = np.random.default_rng(seed)
    out = []
    for N in nts:
        k = int(rng.poisson(mu * N))
        r = _skeleton_clone_sum(k, rng, cap=5000)
        r = min(r, 5000)          # 管级截断：多克隆合计也可能破 poly 上限
        if fraction_plated < 1.0:
            r = int(rng.binomial(r, fraction_plated))
        out.append(r)
    return out


def _ll_pair_factory(counts1, nts1, counts2, nts2, z):
    rs1 = np.asarray(counts1, dtype=int)
    rs2 = np.asarray(counts2, dtype=int)
    nt1 = np.asarray(nts1, dtype=float)
    nt2 = np.asarray(nts2, dtype=float)
    ll1 = _ll_per_tube_factory(rs1, nt1, z)
    ll2 = _ll_per_tube_factory(rs2, nt2, z)
    return ll1, ll2, nt1.min(), nt1.max(), nt2.min(), nt2.max()


def _mu_mle(ll_fn, nt_min, nt_max):
    ln_lo, ln_hi = math.log(1e-6 / nt_max), math.log(1e3 / max(nt_min, 1.0))
    x, ll = _maximize_1d(ll_fn, ln_lo, ln_hi)
    return math.exp(x), ll, ln_lo, ln_hi


def lrt_equal_rates(counts1, nts1, counts2, nts2, fraction_plated: float = 1.0) -> dict:
    """H0: μ1 = μ2 的似然比检验（df=1）。

    返回 mu1/mu2 的点估计、LLR、p 值（p = erfc(sqrt(LLR/2))，χ²₁ 生存函数）。
    """
    ll1, ll2, n1lo, n1hi, n2lo, n2hi = _ll_pair_factory(counts1, nts1, counts2, nts2, fraction_plated)
    mu1, ll1m, _, _ = _mu_mle(ll1, n1lo, n1hi)
    mu2, ll2m, _, _ = _mu_mle(ll2, n2lo, n2hi)

    def pooled(x):
        return ll1(x) + ll2(x)

    lo = math.log(1e-6 / max(n1hi, n2hi))
    hi = math.log(1e3 / max(min(n1lo, n2lo), 1.0))
    _, ll_pooled = _maximize_1d(pooled, lo, hi)

    llr = 2.0 * ((ll1m + ll2m) - ll_pooled)
    llr = max(llr, 0.0)
    p = math.erfc(math.sqrt(llr / 2.0))
    return dict(mu1=mu1, mu2=mu2, llr=llr, p=p)


def foldchange(counts1, nts1, counts2, nts2, fraction_plated: float = 1.0) -> dict:
    """θ = μ2/μ1 的 profile-likelihood 点估计与 95% CI。

    内层对 ln μ1（给定 θ）最大化 ll1(μ1)+ll2(θ·μ1)；外层对 ln θ 游走+二分
    （阈值 ΔLL = 1.921）。θ>1 即 μ2 相对更高。
    """
    ll1, ll2, n1lo, n1hi, n2lo, n2hi = _ll_pair_factory(counts1, nts1, counts2, nts2, fraction_plated)
    mu1, _, _, _ = _mu_mle(ll1, n1lo, n1hi)
    mu2, _, _, _ = _mu_mle(ll2, n2lo, n2hi)
    theta_hat = mu2 / mu1

    ln1_lo = math.log(1e-6 / n1hi)
    ln1_hi = math.log(1e3 / max(n1lo, 1.0))
    last_x1 = [math.log(mu1)]                   # warm 链：上一个 θ 的最优 μ1 是下一个的初值

    def ll_prof(ln_theta: float) -> float:
        # 内层 ln μ1 的可行域：μ1 在 D1 有效域内，且 θ·μ1 在 D2 有效域内
        lo = max(ln1_lo, math.log(1e-6 / n2hi) - ln_theta)
        hi = min(ln1_hi, math.log(1e3 / max(n2lo, 1.0)) - ln_theta)
        if hi <= lo:
            return -math.inf

        def inner(x):
            return ll1(x) + ll2(x + ln_theta)

        # warm 链：内层最优点随 θ 平滑移动，用上一次解 ±ln50 做括号，纯 golden（括号内单峰）
        xg = min(max(last_x1[0], lo + 1e-9), hi - 1e-9)
        x1_opt, v = _maximize_1d(inner, max(xg - math.log(50.0), lo),
                                 min(xg + math.log(50.0), hi), n_grid=0)
        last_x1[0] = x1_opt
        return v

    ll_max = ll_prof(math.log(theta_hat))
    lo_ln, hi_ln = _profile_ci_1d(ll_prof, math.log(theta_hat), ll_max,
                                  math.log(1e-8), math.log(1e8))
    return dict(theta=theta_hat,
                ci=(math.exp(lo_ln), math.exp(hi_ln)),
                mu1=mu1, mu2=mu2)
