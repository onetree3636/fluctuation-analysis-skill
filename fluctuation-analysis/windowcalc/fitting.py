"""v0.3 拟合模块——F(c) Hill 剂量响应 与 τ_off 洗脱衰变。

数据到位（波动分析规程梯队二）后，把测定结果按 CSV 喂进来即可：

- F(c)：各剂量的 M_on 实测 → Hill 拟合 → 任意剂量的 F 内插/CI（三旋钮优化的
  剂量轴从"离散点"升级为"连续曲线"）。
- τ_off：洗脱后 0/2/4/8/16 h 的系列速率 → 指数衰变拟合 → 有效窗口的尾部负荷。

两个拟合都用 bootstrap 给 CI（重采样数据点，非参数），并带退化保护
（数据点不足/非单调 → 明确警告，不硬拟合）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# Hill 拟合的搜索域（对数剂量轴）
_EC50_LO, _EC50_HI = 1e-4, 1e4     # mM
_N_LO, _N_HI = 0.3, 8.0
_FMAX_LO, _FMAX_HI = 1.01, 1e5


# ---------------------------------------------------------------- F(c) Hill

def hill_fc(c: float | np.ndarray, fmax: float, ec50: float, n_exp: float,
            fmin: float = 1.0) -> float | np.ndarray:
    """F(c) = fmin + (fmax−fmin)·cⁿ/(EC50ⁿ+cⁿ)。c=0 → fmin（泄漏已并入 fmin）。"""
    c = np.maximum(np.asarray(c, dtype=float), 0.0)
    ec50 = max(ec50, 1e-12)
    frac = np.where(c > 0, c ** n_exp / (ec50 ** n_exp + c ** n_exp), 0.0)
    return fmin + (fmax - fmin) * frac


def _hill_ss(doses, folds, fmax, ec50, n_exp, fmin):
    return float(np.sum((np.asarray(folds) - hill_fc(doses, fmax, ec50, n_exp, fmin)) ** 2))


def _coord_hill(doses, folds, fmin=1.0, fmax_init=None):
    """三参数 (fmax, ec50, n) 的网格+坐标下降（无 scipy 依赖）。

    返回 (fmax, ec50, n, ss)；数据无法支撑 Hill 形状时返回 None（由调用方判退化）。
    """
    doses = np.asarray(doses, float)
    folds = np.asarray(folds, float)
    fmax0 = fmax_init or max(float(folds.max()) / max(fmin, 1e-9) * 1.2, 2.0)
    fmax0 = float(min(max(fmax0, _FMAX_LO), min(_FMAX_HI, fmax0 * 50)))

    # 粗网格：(ec50, n) 全扫，fmax 每行解析最优（线性最小二乘）
    ec_grid = np.geomspace(_EC50_LO, _EC50_HI, 120)
    n_grid = np.geomspace(_N_LO, _N_HI, 40)
    best = (None, None, None, math.inf)
    for ec in ec_grid:
        for nn in n_grid:
            frac = np.where(doses > 0, (doses / ec) ** nn / (1.0 + (doses / ec) ** nn), 0.0)
            # fmax 线性解：min Σ (folds − fmin − fmax·frac)² → fmax* = Σ frac·(folds−fmin)/Σ frac²
            num = float(np.sum(frac * (folds - fmin)))
            den = float(np.sum(frac ** 2))
            if den <= 0:
                continue
            fmax_star = num / den
            if not (_FMAX_LO <= fmax_star <= _FMAX_HI):
                continue
            ss = _hill_ss(doses, folds, fmax_star, ec, nn, fmin)
            if ss < best[3]:
                best = (fmax_star, ec, nn, ss)
    if best[0] is None:
        return None
    fmax_b, ec_b, n_b, ss_b = best

    # 坐标下降精修（3 轮 × 每参数对数括号收缩）
    spans = [(fmax_b, fmax_b / 10, fmax_b * 10), (ec_b, ec_b / 20, ec_b * 20),
             (n_b, max(n_b / 5, _N_LO), min(n_b * 5, _N_HI))]
    cur = [fmax_b, ec_b, n_b]
    for _ in range(3):
        for pi, (val, lo, hi) in enumerate(spans):
            lo = max(lo, (_FMAX_LO, _EC50_LO, _N_LO)[pi])
            hi = min(hi, (_FMAX_HI, _EC50_HI, _N_HI)[pi])
            trial = np.geomspace(max(lo, 1e-9), hi, 25)
            sss = [math.inf] * len(trial)
            for t_i, t in enumerate(trial):
                c3 = cur.copy()
                c3[pi] = t
                sss[t_i] = _hill_ss(doses, folds, c3[0], c3[1], c3[2], fmin)
            j = int(np.argmin(sss))
            cur[pi] = float(trial[j])
            spans[pi] = (cur[pi], max(trial[max(j - 1, 0)], lo), min(trial[min(j + 1, len(trial) - 1)], hi))
    ss = _hill_ss(doses, folds, cur[0], cur[1], cur[2], fmin)
    return cur[0], cur[1], cur[2], ss


@dataclass
class HillFit:
    fmax: float
    ec50: float
    n_exp: float
    fmin: float
    r2: float
    ci: dict = field(default_factory=dict)     # fmax/ec50/n → (lo, hi) bootstrap
    warnings: list = field(default_factory=list)

    def predict(self, dose_mM: float) -> float:
        return float(hill_fc(dose_mM, self.fmax, self.ec50, self.n_exp, self.fmin))


def fit_hill_fc(doses, m_on, m_base: float, ci_pairs=None, fmin: float = 1.0,
                n_boot: int = 200, seed: int = 0) -> HillFit:
    """剂量响应拟合。

    doses : 各剂量（mM，须 ≥1 个 >0 点）
    m_on  : 该剂量的 M_on 实测（每基因组每代，与 m_base 同口径）
    m_base  : 基线速率（F = m_on/m_base）
    ci_pairs : 每点 (ci_lo, ci_hi)——给了就按它的尺度加噪做 bootstrap（可选）；
               不给则用残差重采样。bootstrap 只重采误差，**不重采剂量点**。
    """
    doses = np.asarray(doses, float)
    folds = np.asarray(m_on, float) / float(m_base)
    warnings = []
    if ci_pairs is not None and len(ci_pairs) != len(doses):
        raise ValueError(
            f"ci_pairs 有 {len(ci_pairs)} 个但剂量点有 {len(doses)} 个——须逐点配对")
    if np.any(folds < 0) or not np.all(np.isfinite(folds)):
        raise ValueError("m_on/m_base 必须是非负有限数（检查 m_base 是否填错/为 0）")
    if len(doses) < 2:
        warnings.append("剂量点 <2：无法拟合 Hill，返回常数 F=last/first")
        f = float(folds[-1]) if len(folds) else 1.0
        return HillFit(fmax=f, ec50=math.nan, n_exp=math.nan, fmin=fmin,
                       r2=math.nan, warnings=warnings + ["DEGENERATE"])
    if np.any(np.diff(folds) < 0) and len(doses) > 2:
        warnings.append("F(c) 非单调下降——检查剂量排序与测定噪声")

    if np.all(doses <= 0):
        warnings.append("所有剂量为 0：无法定 EC50（只能给常数 F）")
        return HillFit(fmax=float(folds.mean()), ec50=math.nan, n_exp=math.nan,
                       fmin=fmin, r2=math.nan, warnings=warnings + ["DEGENERATE"])

    fitted = _coord_hill(doses, folds, fmin)
    if fitted is None:
        warnings.append(
            "Hill 形状拟合失败（数据全在基线以下 / 单调方向与预期相反）："
            f"只返回常数 F={float(folds.max()):.3g}，不给 EC50——请检查单位与剂量排序")
        return HillFit(fmax=float(max(folds.max(), fmin)), ec50=math.nan, n_exp=math.nan,
                       fmin=fmin, r2=math.nan, warnings=warnings + ["DEGENERATE"])
    fmax_b, ec_b, n_b, _ = fitted
    ss0 = _hill_ss(doses, folds, fmax_b, ec_b, n_b, fmin)
    tss = float(np.sum((folds - folds.mean()) ** 2))
    r2 = 1.0 - ss0 / tss if tss > 0 else math.nan      # 完全平坦 → 无残差可解释
    if tss <= 0:
        warnings.append("F 完全平坦（TSS=0）：Hill 参数不可辨识，R² 记为 nan")
    fit0 = HillFit(fmax=fmax_b, ec50=ec_b, n_exp=n_b, fmin=fmin, r2=r2)

    # bootstrap CI：**残差重采样**（剂量设计固定）。
    # 不能重采剂量点：3–5 个点重采会造出重复剂量，等价于给某个点加权，
    # fmax/EC50 的 CI 会被撑到几个数量级（实测 fmax [357, 1.45e4]）。
    rng = np.random.default_rng(seed)
    fitted_vals = hill_fc(doses, fmax_b, ec_b, n_b, fmin)
    resid = folds - fitted_vals
    resid = resid - resid.mean()
    scale = float(np.max(np.abs(fitted_vals))) or 1.0
    boots = []
    if np.std(resid) <= 1e-12 * scale:
        warnings.append("残差≈0（数据由 Hill 曲线精确生成）：bootstrap CI 无意义，已跳过")
    else:
        if len(np.unique(doses)) == len(doses):
            warnings.append("每剂量只测 1 次：CI 来自残差重采样，包含模型失拟——只作量级参考；"
                            "要真正的测量不确定度需每剂量 ≥2 次重复")
        if ci_pairs is not None:
            sig = np.array([max((hi - lo) / (2 * 1.96), 1e-15) for lo, hi in ci_pairs], float)
        for _ in range(n_boot):
            if ci_pairs is not None:                     # 有逐点 CI → 按它的尺度加噪
                y_b = fitted_vals + rng.normal(0.0, 1.0, len(doses)) * sig
            else:
                y_b = fitted_vals + rng.choice(resid, size=len(resid), replace=True)
            r = _coord_hill(doses, y_b, fmin, fmax_init=fmax_b)
            if r:
                boots.append(r)
    if len(boots) >= 20:
        arr = np.array(boots)
        for j, key in enumerate(("fmax", "ec50", "n_exp")):
            fit0.ci[key] = (float(np.percentile(arr[:, j], 2.5)),
                            float(np.percentile(arr[:, j], 97.5)))
    else:
        warnings.append("bootstrap 有效样本不足，CI 不可用")
    fit0.warnings = warnings
    return fit0


# ---------------------------------------------------------------- τ_off 衰变

@dataclass
class TauOffFit:
    tau_off_h: float
    ci: tuple = None
    halflife_h: float = None
    mu_off: float = None
    mu_on_fit: float = None
    tail_passengers: float = None     # 洗脱后尾部全积分 ∫₀^∞ 的额外乘客（每基因组）
    tail_1tau: float = None           # 只再培养 1 个 τ 的额外乘客（=63.2% 全量）
    warnings: list = field(default_factory=list)


def fit_tau_off(hours, mus, mu_off: float, td_h: float = 2.0,
                m_on: float | None = None, n_boot: int = 200, seed: int = 0) -> TauOffFit:
    """洗脱后系列速率的指数衰变拟合：μ(t) = μ_off + (μ(0)−μ_off)·e^(−t/τ)。

    hours : 洗脱后取样时间（h，含 0）
    mus   : 各时间点的 μ̂ 实测（per cell per division；mean 模式或 per_tube 均可）
    mu_off: off 态速率（strain_card 的 M_off，视为已知）
    m_on  : 诱导态速率（尾部负荷的参考值）；缺省用拟合出的 μ_on。注意它只影响
            尾部负荷的折算，不参与 τ 的拟合（τ 由曲线的衰减形状决定）。
    """
    t = np.asarray(hours, float)
    y = np.asarray(mus, float)
    warnings = []
    if len(t) < 3:
        warnings.append("时间点 <3：τ 不可辨识，返回 NA")
        return TauOffFit(tau_off_h=math.nan, warnings=warnings + ["DEGENERATE"])
    if np.any(y <= 0):
        raise ValueError("μ̂ 必须为正（对数/速率口径）；请检查输入是否为负或 0")

    def model(tau, mu_on):
        return mu_off + (mu_on - mu_off) * np.exp(-t / max(tau, 1e-6))

    def ss(tau, mu_on):
        r = y - model(tau, mu_on)
        return float(np.sum(r ** 2))

    # τ 一维搜索（μ_on 内层线性解：给定 τ，μ_on* = Σ w·(y−μ_off)/Σ w²，w=e^(−t/τ)）
    taus = np.geomspace(0.05, 200, 400)
    best = (None, None, math.inf)
    for tau in taus:
        w = np.exp(-t / tau)
        num = float(np.sum(w * (y - mu_off)))
        den = float(np.sum(w ** 2))
        if den <= 0:
            continue
        mu_on = num / den
        s = ss(tau, mu_on)
        if s < best[2]:
            best = (tau, mu_on, s)
    tau_b, mu_on_b, ss_b = best

    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        noise = rng.normal(0.0, max(math.sqrt(ss_b / len(t)), 1e-12), len(t))
        r = y + noise
        bb = (None, None, math.inf)
        for tau in np.geomspace(0.05, 200, 200):
            w = np.exp(-t / tau)
            num = float(np.sum(w * (r - mu_off)))
            den = float(np.sum(w ** 2))
            mu_on_b2 = num / den if den > 0 else mu_on_b
            s2 = float(np.sum((r - (mu_off + (mu_on_b2 - mu_off) * np.exp(-t / tau))) ** 2))
            if s2 < bb[2]:
                bb = (tau, mu_on_b2, s2)
        if bb[0] and math.isfinite(bb[0]):
            boots.append(bb[0])
    ci = ((float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))
          if len(boots) >= 20 else None)
    if ci is None:
        warnings.append("bootstrap 有效样本不足，CI 不可用")

    # 尾部负荷：洗脱后残留活性带来的一辈子额外乘客 = ∫₀^∞(μ(t)−μ_off)dt / td_h
    # （单位：每基因组；m_on 给定时以它为准——t=0 点噪声大时用它替换拟合值）
    mu_on_ref = float(m_on) if m_on is not None else mu_on_b
    tail = (mu_on_ref - mu_off) * tau_b / td_h
    tail_1tau = tail * (1.0 - math.exp(-1.0))      # 只再培养 1 个 τ 的量（63.2%）

    out = TauOffFit(tau_off_h=tau_b, ci=ci,
                    halflife_h=tau_b * math.log(2.0),
                    mu_off=mu_off, mu_on_fit=mu_on_b, tail_passengers=tail,
                    tail_1tau=tail_1tau, warnings=warnings)
    if tau_b > 12:
        out.warnings.append(f"τ_off={tau_b:.1f} h 很长——有效窗口远大于名义脉冲，"
                            "考虑洗脱后低温/换培终止表达")
    return out
