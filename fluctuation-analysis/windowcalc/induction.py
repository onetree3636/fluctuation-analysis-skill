"""诱导剂剂量轴的机制模型：表达 → 竞争占位 → 突变率。

三个方程（README/技能 mutator-window 的 v0.4 模型）：

    E(c)  = Emax · cⁿ / (EC50ⁿ + cⁿ)          表达（Hill）
    p(c)  = E(c) / (E(c) + K) = q · h(c)      占位（Langmuir 竞争）
            h(c) = cⁿ/(cⁿ+C50ⁿ)，q = Emax/(Emax+K) ≤ 1
    M(c)  = M_off + ΔM · φ(p)                 速率

混合规则 φ（复制的两种 Pol III 核心怎么装）：
    linear     φ = p                 平均效应
    two_core   φ = 2p − p²           至少一个核心是 M7（ψ = 2h − q h²）
    both_core  φ = p²                两个核心都要 M7（ψ = h²）

可辨识性（重要，别把不可辨识的量当结论）：
    - 从 (c, M) 数据只能稳健定出 **A = ΔM·q**（有效占位幅度）、C50、n；
    - q 只在 two_core 下由曲率可分（因为 ψ=2h−q h² 里 2A 与 Aq 可分别定）；
    - both_core 下 q 与 A 退化（只能定 B = A·q）；
    - linear 下 q 完全不可辨识。
    要拆开 q 与 ΔM，必须有独立表达读出（报告基因/定量）或"纯 M7"锚点。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

MIXINGS = ("linear", "two_core", "both_core")
_MIX_LABEL = {"linear": "线性 φ=p", "two_core": "至少一个核心 φ=2p−p²",
              "both_core": "两个核心 φ=p²"}


def hill(c, ec50, n):
    """h(c) = cⁿ/(cⁿ + C50ⁿ)；c=0 → 0。

    ec50 / n 可以是数组（MC 逐样本采样剂量轴参数时），返回广播后的形状。
    """
    c = np.maximum(np.asarray(c, dtype=float), 0.0)
    e = np.maximum(np.asarray(ec50, dtype=float), 1e-12)
    nn = np.asarray(n, dtype=float)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        cn = np.power(c, nn)
        out = np.where(c > 0, cn / (np.power(e, nn) + cn), 0.0)
    return float(out) if np.ndim(out) == 0 else out


def psi(h, mixing: str, q: float = 1.0):
    """归一化混合基函数 ψ：M = M_off + A·ψ。"""
    h = np.asarray(h, dtype=float)
    if mixing == "linear":
        return h
    if mixing == "two_core":
        q = min(max(float(q), 1e-6), 1.0)
        return 2.0 * h - q * h * h
    if mixing == "both_core":
        return h * h
    raise ValueError(f"mixing={mixing!r} 不在 {MIXINGS}")


def psi_max(mixing: str, q: float = 1.0) -> float:
    """ψ 在 h=1 的值（= 饱和时的相对幅度）。"""
    if mixing == "linear":
        return 1.0
    if mixing == "two_core":
        q = min(max(float(q), 1e-6), 1.0)
        return 2.0 - q
    return 1.0


def m_of_c(c, amplitude: float, ec50: float, n: float, mixing: str,
           q: float = 1.0, m_off: float = 0.0):
    """M(c) = M_off + A·ψ(h(c))。"""
    return m_off + amplitude * psi(hill(c, ec50, n), mixing, q)


def saturation_dose(frac: float, ec50: float, n: float, mixing: str,
                    q: float = 1.0) -> float:
    """达到 ψ(1) 的 frac（如 0.95）所需的剂量 c。"""
    target = frac * psi_max(mixing, q)
    lo, hi = math.log(max(ec50, 1e-9)) - 8.0, math.log(max(ec50, 1e-9)) + 8.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if float(psi(hill(math.exp(mid), ec50, n), mixing, q)) >= target:
            hi = mid
        else:
            lo = mid
    return math.exp(0.5 * (lo + hi))


# ---------------------------------------------------------------- 拟合

@dataclass
class MixFit:
    mixing: str
    amplitude: float          # A（M_off 之上的系数）
    ec50_mM: float
    n_hill: float
    q: float = 1.0            # 仅 two_core 有意义
    ss: float = math.inf
    n: int = 0
    k: int = 3
    warnings: list = field(default_factory=list)

    @property
    def aic(self) -> float:
        if not np.isfinite(self.ss) or self.n <= 0:
            return math.inf
        return self.n * math.log(max(self.ss / self.n, 1e-300)) + 2 * self.k

    @property
    def delta_sat(self) -> float:
        """饱和时在 M_off 之上的可观测增幅（= A·ψ(1)）。"""
        return self.amplitude * psi_max(self.mixing, self.q)


@dataclass
class InductionFit:
    best: MixFit
    models: list                      # 所有参与比较的模型
    m_off: float
    ci: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def predict(self, c) -> float:
        return float(m_of_c(c, self.best.amplitude, self.best.ec50_mM,
                            self.best.n_hill, self.best.mixing, self.best.q,
                            self.m_off))


def _solve_amplitude(doses, y, ec50, n, mixing, q):
    """给定 (ec50, n, q)，A 的最小二乘解（模型对 A 线性）。"""
    h = hill(doses, ec50, n)
    w = psi(h, mixing, q)
    den = float(np.dot(w, w))
    if den <= 0:
        return None
    return float(np.dot(w, y) / den)


def _fit_one(doses, y, mixing):
    """网格 + 坐标下降。y = M − M_off（≥0 的部分）。"""
    d = np.asarray(doses, float)
    y = np.asarray(y, float)
    lo_d, hi_d = max(d[d > 0].min(), 1e-6), d.max()
    ec_grid = np.geomspace(lo_d / 10.0, hi_d * 10.0, 40)
    n_grid = np.geomspace(0.5, 6.0, 20)
    q_grid = np.linspace(0.2, 1.0, 9) if mixing == "two_core" else [1.0]

    best = None
    for q in q_grid:
        for ec in ec_grid:
            for nn in n_grid:
                a = _solve_amplitude(d, y, ec, nn, mixing, q)
                if a is None or a <= 0:
                    continue
                pred = a * psi(hill(d, ec, nn), mixing, q)
                ss = float(np.sum((y - pred) ** 2))
                if best is None or ss < best[0]:
                    best = (ss, a, ec, nn, q)
    if best is None:
        return None
    ss, a, ec, nn, q = best

    # 坐标下降精修
    for _ in range(3):
        for idx in range(3):
            if idx == 0:
                cand = np.geomspace(max(ec / 3, lo_d / 100), ec * 3, 25)
                for c_ec in cand:
                    aa = _solve_amplitude(d, y, c_ec, nn, mixing, q)
                    if aa and aa > 0:
                        s = float(np.sum((y - aa * psi(hill(d, c_ec, nn), mixing, q)) ** 2))
                        if s < ss:
                            ss, a, ec = s, aa, c_ec
            elif idx == 1:
                for c_n in np.geomspace(max(nn / 3, 0.2), nn * 3, 25):
                    aa = _solve_amplitude(d, y, ec, c_n, mixing, q)
                    if aa and aa > 0:
                        s = float(np.sum((y - aa * psi(hill(d, ec, c_n), mixing, q)) ** 2))
                        if s < ss:
                            ss, a, nn = s, aa, c_n
            else:
                if mixing == "two_core":
                    for c_q in np.linspace(max(0.05, q - 0.2), min(1.0, q + 0.2), 17):
                        aa = _solve_amplitude(d, y, ec, nn, mixing, c_q)
                        if aa and aa > 0:
                            s = float(np.sum((y - aa * psi(hill(d, ec, nn), mixing, c_q)) ** 2))
                            if s < ss:
                                ss, a, q = s, aa, c_q
    k = 4 if mixing == "two_core" else 3
    return MixFit(mixing=mixing, amplitude=a, ec50_mM=ec, n_hill=nn, q=q,
                  ss=ss, n=len(d), k=k)


def fit_induction(doses, m_values, m_off: float, mixing: str = "auto",
                  n_boot: int = 200, seed: int = 0) -> InductionFit:
    """拟合剂量-速率曲线；mixing='auto' 时用 AIC 在三种混合规则里选。"""
    d = np.asarray(doses, float)
    m = np.asarray(m_values, float)
    if d.size != m.size or d.size < 3:
        raise ValueError("至少 3 个剂量点，且 doses 与 m 等长")
    if np.any(d < 0):
        raise ValueError("剂量不能为负")
    below = m < m_off
    y = np.maximum(m - m_off, 0.0)
    if np.any(below):
        warnings_pre = [f"{int(below.sum())} 个点的实测 M 低于 M_off（测量噪声或 M_off 偏大）：已按 0 处理"]
    else:
        warnings_pre = []

    cand = MIXINGS if mixing == "auto" else (mixing,)
    fits = [f for f in (_fit_one(d, y, mx) for mx in cand) if f is not None]
    if not fits:
        raise ValueError("拟合失败：检查剂量是否全为 0 或数据是否单调")
    full = sorted(fits, key=lambda f: f.aic)
    _order = {"linear": 0, "both_core": 1, "two_core": 2}
    best = full[0]
    demoted = None
    for f in full[1:]:
        if f.aic - best.aic < 2.0 and _order[f.mixing] < _order[best.mixing]:
            demoted, best = best, f
    fits = [best] + [f for f in full if f is not best]

    warns = list(warnings_pre)
    if demoted is not None:
        warns.append(f"AIC 最优是 {demoted.mixing}（{demoted.aic:.1f}），但与 {best.mixing}"
                     f"（{best.aic:.1f}）差 <2：按简约原则取 {best.mixing}——数据分不出两者")
    if d.size < 4:
        warns.append("只有 3 个剂量点：n 与 q 不稳，建议 ≥5 点")
    c50 = best.ec50_mM
    if d.max() < c50:
        warns.append(f"最大剂量 {d.max():g} < C50 {c50:g}：数据未跨膝点，C50/饱和幅度不可信")
    else:
        c95 = saturation_dose(0.95, best.ec50_mM, best.n_hill, best.mixing, best.q)
        if d.max() < c95:
            warns.append(f"最大剂量 {d.max():g} < 95% 饱和剂量 {c95:g}：饱和幅度 A 仍属外推")
    second = next((f for f in full if f is not best), None)
    if second is not None and second.aic - best.aic < 2.0:
        warns.append(f"次佳模型 {second.mixing} 的 ΔAIC={second.aic - best.aic:.2f} <2："
                     f"数据分不出两种混合规则，别过度解读 q")

    # 残差 bootstrap CI（参数 + 饱和幅度）
    rng = np.random.default_rng(seed)
    boot = []
    fitted = best.amplitude * psi(hill(d, best.ec50_mM, best.n_hill), best.mixing, best.q)
    resid = y - fitted
    resid = resid - resid.mean()
    if n_boot and np.std(resid) > 1e-15 * max(abs(y).max(), 1.0):
        for _ in range(n_boot):
            yb = fitted + rng.choice(resid, size=len(y), replace=True)
            f = _fit_one(d, yb, best.mixing)
            if f is not None:
                boot.append((f.amplitude, f.ec50_mM, f.n_hill, f.q, f.delta_sat))
    ci = {}
    if len(boot) >= 20:
        arr = np.array(boot)
        for j, key in enumerate(("amplitude", "ec50_mM", "n_hill", "q", "delta_sat")):
            if best.mixing != "two_core" and key in ("q",):
                continue
            ci[key] = (float(np.percentile(arr[:, j], 2.5)),
                       float(np.percentile(arr[:, j], 97.5)))
    else:
        warns.append("残差过小或 bootstrap 样本不足，CI 不可用")

    return InductionFit(best=best, models=fits, m_off=m_off, ci=ci, warnings=warns)


def saturation_table(fit: InductionFit, doses) -> list:
    """每个实测剂量处：h、ψ/ψ(1)、M 预测、距饱和还差多少。"""
    b = fit.best
    rows = []
    for c in doses:
        h = float(hill(c, b.ec50_mM, b.n_hill))
        rel = float(psi(h, b.mixing, b.q)) / psi_max(b.mixing, b.q)
        m_pred = fit.predict(c)
        act = rel                      # 以最大活性归一（参考剂量处见 activity_ref）
        rows.append(dict(dose_mM=float(c), h=h, fraction_of_Amax=rel,
                         activity=act, supply_per_load=(act / m_pred if m_pred > 0 else 0.0),
                         M_pred=m_pred,
                         M_off=fit.m_off,
                         delta=m_pred - fit.m_off))
    return rows


def adequacy(fit: InductionFit, doses) -> dict:
    """数据充分性：是否跨膝点、是否到饱和、能否分辨混合规则。"""
    b = fit.best
    d = np.asarray(doses, float)
    c95 = saturation_dose(0.95, b.ec50_mM, b.n_hill, b.mixing, b.q)
    second = fit.models[1] if len(fit.models) > 1 else None
    return dict(
        spans_knee=bool(d.max() >= b.ec50_mM),
        reaches_95=bool(d.max() >= c95),
        c50_mM=b.ec50_mM,
        c95_mM=c95,
        max_dose=float(d.max()),
        n_points=int(d.size),
        model_gap_aic=(None if second is None else second.aic - b.aic),
    )
