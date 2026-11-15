"""ramp：多靶点渐升的梯级规划器——把"预测"锚在 PAP 尾部实测曲线上。

预测的三层（全部锚在可测/已测的量上，每一步被下一步的实测审计）：

  [R1] PAP 尾部拟合    f(c) = f_ref × 10^(−(c − c_ref)/δ)
       输入 (c, 克隆数, CFU)——**同一个窗口后群体**在多个浓度上的铺板计数
       （其中一块本来就是下一步的工作板，多涂几块即得）。
       δ = 十倍频浓度宽：频率每掉一个数量级，浓度能抬多少。
       零克隆行不进拟合，按 rule of three 作检测上界记录（拟合高于上界要报警）。
  [R2] 每档 f_step     平稳补充假设：等步长下，每档的升档频率 ≈ 首档实测值
       （每个诱导窗口的突变供给把尾部重新填满）。实证：PLoS Genet 2025
       （pgen.1011627）的梯度孔自适应爬升与其"每次突变抬升 log-MIC"的
       加性模型；临床异质性耐药的 PAP 尾部 log-线性也是同一观察。
  [R3] 梯级累计        乘客三段账按窗可加（无性谱系单调、驱动与乘客独立泊松流）；
       失守概率 = 1 − exp(−E[P_总] × L_prot/G)——**总乘客单次代入**，
       不把每档的失守概率相加（概率相加会重复计入双命中）。

规划逻辑：Δ 越大 → 档数越少 → 乘客越少，但 f_step 沿曲线下落越深、越不可行。
规划器解出**可行（含 margin）的最大步长**，给出最省乘客的梯级，
并附 2 个更保守的备选档距供权衡。

明确不预测：突变身份、MIC 轨迹的分子机制、CNV 何时发生（WGS 对账的事）。
明确不做：μ×T_set / 几何供给外推（2026-11-13 已移除的旧路线）——
本模块的预测全部锚在 PAP 实测曲线上，曲线之外的每一步都要重新实测。
"""

from __future__ import annotations

import math
from typing import Sequence

from .fstep import newcombe_diff, wilson_ci

Z95 = 1.959963985


# --------------------------------------------------------------------------- PAP 拟合

def fit_pap(points: Sequence[tuple[float, int, float]],
            c_ref: float | None = None) -> dict:
    """PAP 尾部拟合：log10 f = a + s·c（s < 0），δ = −1/s。

    points: (c, clones, cfu)——同一窗口后群体在各浓度的铺板计数；
            同浓度多板自动合并（克隆数、CFU 分别求和，等价 run_step 的 pool）。
    零克隆行：不进拟合，按 rule of three 记 3/CFU 上界（censored），
    拟合值高于上界的点会报警。
    返回含 slope/intercept/δ/CI/预测函数所需全部量的 dict。
    """
    if not points:
        raise ValueError("PAP 点集为空——至少需要 2 个有克隆的浓度点才能定 δ")
    pool: dict[float, list[tuple[int, float]]] = {}
    for c, k, n in points:
        c = float(c)
        if c < 0:
            raise ValueError(f"浓度 c={c!r} 不能为负")
        if n <= 0:
            raise ValueError(f"c={c:g} 的涂板 CFU 必须 > 0")
        k = int(k)
        if k < 0 or k > n:
            raise ValueError(f"c={c:g} 的克隆数 {k} 非法（须 0 ≤ k ≤ CFU={n:g}）")
        pool.setdefault(c, [0, 0.0])
        pool[c][0] += k
        pool[c][1] += n
    pos = sorted((c, k, n) for c, (k, n) in pool.items() if k > 0)
    censored = sorted((c, 3.0 / n) for c, (k, n) in pool.items() if k == 0)
    if len(pos) < 2:
        raise ValueError(
            "PAP 至少需要 2 个**有克隆**的浓度点才能定尾部斜率 δ"
            f"（当前 {len(pos)} 个；零克隆行只能当检测上界，不进拟合）")
    xs = [p[0] for p in pos]
    ys = [math.log10(p[1] / p[2]) for p in pos]
    n = len(pos)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx <= 0:
        raise ValueError("PAP 各点的浓度相同——至少要 2 个不同浓度才能定 δ")
    sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    slope = sxy / sxx
    if slope >= 0:
        raise ValueError(
            "PAP 尾部斜率必须 < 0（浓度越高频率越低）；拟合出正斜率说明"
            "数据有误或该表型不是耐受（检查克隆/CFU 是否填反）")
    intercept = ybar - slope * xbar
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    delta = -1.0 / slope
    sigma2 = None
    se_slope = None
    delta_ci = None
    ci_includes_zero = False
    pred_se = None
    if n >= 3:
        sigma2 = ss_res / (n - 2)
        se_slope = math.sqrt(sigma2 / sxx)
        # slope 95%CI → δ CI（s<0 时 1/s 单调减，端点交换）
        s_lo, s_hi = slope - Z95 * se_slope, slope + Z95 * se_slope
        # δ = −1/s：s<0 时 −1/s 随 s 增大而增大，所以**下端点由 s_lo、上端点由 s_hi**。
        # （2026-11-15 修：此前写成 (−1/s_hi, −1/s_lo)，区间上下颠倒，打印成 [大, 小]。）
        delta_ci = (-1.0 / s_lo if s_lo < 0 else float("inf"),
                    -1.0 / s_hi if s_hi < 0 else float("inf"))
        ci_includes_zero = (s_hi >= 0)
        pred_se = lambda c: math.sqrt(sigma2 * (1 + (c - xbar) ** 2 / sxx))  # noqa: E731
    c_ref = float(c_ref) if c_ref is not None else pos[-1][0]
    f_ref = 10.0 ** (intercept + slope * c_ref)
    if f_ref <= 0 or not (f_ref < float("inf")):
        raise ValueError("PAP 拟合出非正的 f_ref——检查输入数量级")
    flags = []
    if ci_includes_zero:
        flags.append("斜率的 95%CI 含 0（s≥0）：δ 上界不可定——该尾段几乎是一条水平线，"
                     "补更高浓度的点，或确认表型确实是耐受")
    for c, ub in censored:
        if (intercept + slope * c) > math.log10(ub):
            flags.append(f"c={c:g} 实测 0 克隆（上界 {ub:.1e}），但拟合值 "
                         f"{10 ** (intercept + slope * c):.1e} 高于上界——"
                         "该浓度附近曲线可能偏乐观")
    return dict(points=points, pooled=pos, censored=censored, n_pos=n,
                slope=slope, intercept=intercept, r2=r2,
                delta=delta, delta_ci=delta_ci, sigma2=sigma2,
                xbar=xbar, sxx=sxx,
                c_ref=c_ref, f_ref=f_ref, flags=flags)


def predict_f(fit: dict, c: float) -> tuple[float, float | None, float | None]:
    """拟合曲线在 c 处的频率点估计与 95%CI（n<3 时 CI 为 None 并须由调用方注明）。"""
    lg = fit["intercept"] + fit["slope"] * float(c)
    f = 10.0 ** lg
    lo = hi = None
    if fit.get("sigma2") is not None:
        se = math.sqrt(fit["sigma2"] * (1 + (c - fit["xbar"]) ** 2 / fit["sxx"]))
        lo = 10.0 ** (lg - Z95 * se)
        hi = 10.0 ** (lg + Z95 * se)
    return f, lo, hi


# --------------------------------------------------------------------------- 洗脱尾

PROPAGATION_MODES = ("manual", "plate", "liquid")


def inter_rung_generations(*, mode: str, n_plate: float = 0.0,
                           cells_per_colony: float | None = None,
                           nt: float | None = None, k: float | None = None,
                           n0: float | None = None, f: float | None = None,
                           where: str = "") -> tuple[float, str]:
    """档间扩增代数——按传代拓扑**自己算**，不再让用户手工折算。

    manual : 直接用 `n_plate`（旧行为，向后兼容）
    plate  : 板传代，1 个细胞 → 可见菌落需要 log2(cells_per_colony) 代
    liquid : 液传代，存活者 `f·Nt` 长到下一窗的 N₀ 需要 log2(N₀/(f·Nt)) 代，
             N₀ 缺省取"刚好够 k"的最大接种量 `Nt − k/f`（与账本第三节同一设计）

    返回 (代数, 来源说明)。参数不够时**直接报错**，不猜默认值。
    """
    if mode not in PROPAGATION_MODES:
        raise ValueError(f"propagation={mode!r} 不认识；可选 {PROPAGATION_MODES}")
    tag = f"（{where}）" if where else ""

    if mode == "manual":
        return float(n_plate), f"manual：n_plate={n_plate:g}"

    if mode == "plate":
        if cells_per_colony is None:
            raise ValueError(f"plate 传代需要 --cells-per-colony{tag}——"
                             "它是菌落大小的实测/假设值，本技能不内置")
        if not (cells_per_colony >= 2):
            raise ValueError(f"cells_per_colony={cells_per_colony!r} 必须 ≥ 2{tag}")
        return math.log2(cells_per_colony), \
            f"plate：log2({cells_per_colony:g}) = 每克隆细胞数"

    # liquid
    if nt is None:
        raise ValueError(f"liquid 传代需要 --nt（周期末饱和密度）{tag}——"
                         "没有它算不出存活者要长几代")
    if f is None or not (f > 0):
        return 0.0, "liquid：该档无实测 f（按 0 代计，须人工核对）"
    n0_use = float(n0) if n0 is not None else float(nt) - float(k) / f
    if not (n0_use > 0):
        raise ValueError(
            f"liquid 传代给不出 N₀：Nt − k/f = {n0_use:.3g} ≤ 0{tag}——"
            "该档 f 太小/ k 太大，先加大容量或降低 k")
    if not (f * nt > 0):
        raise ValueError(f"liquid 传代需要 f·Nt > 0（当前 f={f!r}, nt={nt!r}）{tag}")
    return math.log2(n0_use / (f * nt)), \
        f"liquid：log2({n0_use:.3g}/({f:.2e}×{nt:.3g}))"


def washout_tail(m_on: float, m_off: float | None,
                 tau_off_h: float | None, td_h: float | None) -> float:
    """洗脱后残留活性的额外乘客（每克隆、**每次洗脱**）。

    洗脱时 mutase/效应蛋白不是瞬间消失，而是按 τ_off 衰减；板/液扩的前几代
    仍带着 > M_off 的活性。按"每代 td_h 小时"离散求和：

        tail = (M_on − M_off) · Σ_{g≥1} exp(−g·td/τ)
             = (M_on − M_off) / (exp(td/τ) − 1)              [推导]

    连续极限（td ≪ τ）退化为 (M_on−M_off)·τ/td，与 `fit_tauoff` 的
    ∫₀^∞(μ−μ_off)dt / td 一致。

    两个参数**必须同时给或同时不给**：只给一个时无法折算，直接报错
    （本技能不内置任何默认值——见 docs/generic-policy.md）。
    """
    if tau_off_h is None and td_h is None:
        return 0.0
    if tau_off_h is None or td_h is None:
        raise ValueError(
            "洗脱尾需要 --tau-off-h 与 --td 同时给出（只给一个无法折算）："
            f"当前 tau_off_h={tau_off_h!r}、td_h={td_h!r}")
    if not (tau_off_h > 0):
        raise ValueError(f"tau_off_h={tau_off_h!r} 必须 > 0")
    if not (td_h > 0):
        raise ValueError(f"td_h={td_h!r} 必须 > 0")
    if m_off is None:
        return 0.0                     # 没有泄漏基线就无从定义"超出部分"
    excess = float(m_on) - float(m_off)
    if excess <= 0:
        return 0.0
    return excess / math.expm1(td_h / tau_off_h)


# --------------------------------------------------------------------------- 规划

def plan_ramp(fit: dict, *, m_on: float, k: float, capacity: float,
              c_start: float, c_target: float, genome_mb: float,
              m_off: float | None = None,
              n_win: float = 1.0, n_plate: float = 0.0, n_post: float = 0.0,
              bgc_kb: float | None = None, margin: float = 3.0,
              max_rungs: int = 200, cum_base: float = 0.0,
              tau_off_h: float | None = None, td_h: float | None = None,
              propagation: str = "manual",
              cells_per_colony: float | None = None,
              nt: float | None = None, n0: float | None = None) -> dict:
    """解可行（含 margin）的**最大步长** Δ，输出最省乘客的均匀梯级与累计账。

    平稳补充假设：等步长下每档的 f_step ≈ PAP 曲线在 c_start+Δ 处的实测值
    （窗口的突变供给把尾部重新填满）。该假设在机制天花板/谱系效应介入处失效，
    因此规划结果**每档都要被实测审计**（回填 --pap 重算）。

    m_off 为空 → 窗口外段（板内扩张/终点后）不计入乘客，报告中显式注明。
    tau_off_h/td_h 同时给出时，每档再计一次**洗脱尾**（见 `washout_tail`）；
    只给一个直接报错，都不给则不计并在报告中注明。
    """
    if not (m_on > 0):
        raise ValueError(f"m_on={m_on!r} 必须 > 0（μ_perbp × G）")
    if m_off is not None and m_off < 0:
        raise ValueError(f"m_off={m_off!r} 不能为负")
    if not (k > 0 and capacity > 0):
        raise ValueError("k 与 capacity 都必须 > 0")
    if not (c_target > c_start >= 0):
        raise ValueError(f"需要 c_target > c_start（当前 {c_start!r} → {c_target!r}）")
    if not (margin > 0):
        raise ValueError(f"margin={margin!r} 必须 > 0")
    if not (n_win > 0):
        raise ValueError(f"n_win={n_win!r} 必须 > 0（每档诱导窗口的分裂代数）")
    if not (genome_mb > 0):
        raise ValueError("genome_mb 必须 > 0（保护区占比 L_prot/G 的分母）")
    span = float(c_target) - float(c_start)
    bgc_frac = (bgc_kb * 1000.0) / (genome_mb * 1e6) if bgc_kb else None
    p_tail = washout_tail(m_on, m_off, tau_off_h, td_h)

    # 可行最大步长：f(c_start+Δ) = k·margin/capacity
    thr = (k * margin) / capacity
    lg_thr = math.log10(thr)
    # log10 f(c_start+Δ) = intercept + slope·(c_start+Δ) = lg_thr
    delta_step = (lg_thr - fit["intercept"]) / fit["slope"] - float(c_start)
    f_start, f_start_lo, f_start_hi = predict_f(fit, c_start)
    if f_start * capacity < k:
        return dict(status="卡住", ladder=[], span=span, delta=None, k=float(k),
                    capacity=float(capacity), margin=float(margin),
                    c_start=float(c_start), c_target=float(c_target),
                    m_on=m_on, m_off=m_off, n_win=n_win, n_plate=n_plate,
                    n_post=n_post, bgc_kb=bgc_kb, genome_mb=genome_mb,
                    bgc_frac=bgc_frac, f_start=f_start, stuck_c=c_start,
                    stuck_reason=(
                        f"起步即不可行：当前群体在 c_start 处的频率 {f_start:.2e} × "
                        f"容量 {capacity:.3g} = {f_start * capacity:.3g} < k={k:g}。"
                        "先富集（液体无药/低药预培养扩大耐受群）、加大容量或降低 k。"),
                    p_total=None, p_bgc=None, flags=list(fit.get("flags", [])),
                    p_tail=p_tail, tau_off_h=tau_off_h, td_h=td_h)
    if delta_step <= 0:
        # 起点 f 满足 k，但满足不了 margin——margin 下卡住
        return dict(status="卡住", ladder=[], span=span, delta=None, k=float(k),
                    capacity=float(capacity), margin=float(margin),
                    c_start=float(c_start), c_target=float(c_target),
                    m_on=m_on, m_off=m_off, n_win=n_win, n_plate=n_plate,
                    n_post=n_post, bgc_kb=bgc_kb, genome_mb=genome_mb,
                    bgc_frac=bgc_frac, f_start=f_start, stuck_c=c_start,
                    stuck_reason=(
                        f"margin={margin:g} 下起步即卡住：f(c_start)·capacity = "
                        f"{f_start * capacity:.3g} < k×margin = {k * margin:.3g}。"
                        "降低 margin（接受更少克隆）或加大容量。"),
                    p_total=None, p_bgc=None, flags=list(fit.get("flags", [])),
                    p_tail=p_tail, tau_off_h=tau_off_h, td_h=td_h)
    delta_step = min(delta_step, span)
    n_rungs = max(1, int(math.ceil(span / delta_step - 1e-12)))
    n_rungs = min(n_rungs, max_rungs)
    step = span / n_rungs

    m_off_eff = m_off if m_off is not None else 0.0
    ladder = []
    cum = float(cum_base)          # 已完成档的乘客基数（账本解释器传入）
    c_from = float(c_start)
    f_rung, f_rung_lo, f_rung_hi = predict_f(fit, c_start + step)
    for i in range(1, n_rungs + 1):
        c_to = c_start + i * step
        pw = n_win * m_on
        n_inter, n_inter_src = inter_rung_generations(
            mode=propagation, n_plate=n_plate, cells_per_colony=cells_per_colony,
            nt=nt, k=k, n0=n0, f=f_rung, where=f"第 {i} 档")
        pp = n_inter * m_off_eff
        p_step = pw + pp + p_tail
        cum += p_step
        p_bgc = None if bgc_frac is None else 1.0 - math.exp(-cum * bgc_frac)
        ladder.append(dict(
            i=i, c_from=c_from, c_to=(c_target if i == n_rungs else round(c_to, 12)),
            f_step=f_rung, f_lo=f_rung_lo, f_hi=f_rung_hi,
            need_cfu=k / f_rung if f_rung > 0 else float("inf"),
            surplus=f_rung * capacity / k,
            p_win=pw, p_plate=pp, p_tail=p_tail, passengers=p_step,
            n_inter=n_inter, n_inter_source=n_inter_src,
            cum_passengers=cum, p_bgc=p_bgc))
        c_from = c_to
    # 终点后
    if n_post > 0:
        cum += n_post * m_off_eff + p_tail        # 末次洗脱同样有尾
        if ladder:
            ladder[-1]["cum_passengers"] = cum
            ladder[-1]["p_bgc"] = (None if bgc_frac is None else
                                   1.0 - math.exp(-cum * bgc_frac))
    max_meas = max(p[0] for p in fit["pooled"])
    extrap = (max(c_start + step, c_target) - max_meas) / fit["delta"] \
        if fit["delta"] > 0 else None
    return dict(status="可行", ladder=ladder, span=span,
                step=step, n_rungs=n_rungs, k=float(k), capacity=float(capacity),
                cum_base=float(cum_base),
                margin=float(margin), c_start=float(c_start),
                c_target=float(c_target), m_on=m_on, m_off=m_off,
                genome_mb=genome_mb, bgc_frac=bgc_frac,
                n_win=n_win, n_plate=n_plate, n_post=n_post, bgc_kb=bgc_kb,
                p_tail=p_tail, tau_off_h=tau_off_h, td_h=td_h,
                propagation=propagation, cells_per_colony=cells_per_colony,
                f_start=f_start, f_rung=f_rung,
                p_total=cum, p_bgc=(None if bgc_frac is None else
                                    1.0 - math.exp(-cum * bgc_frac)),
                horizon_max_c=max_meas,
                extrapolation_decades=extrap,
                stuck_c=None, stuck_reason=None,
                flags=list(fit.get("flags", [])),
                m_off_counted=(m_off is not None))


# --------------------------------------------------------------------------- 渲染

def render(plan: dict) -> str:
    L: list[str] = []
    fit = plan.get("fit", {})
    L.append("一、PAP 尾部拟合（预测的锚点）")
    if fit:
        L.append(f"  f(c) = f_ref × 10^(−(c−c_ref)/δ)；正克隆点 {fit['n_pos']} 个，"
                 f"零克隆上界行 {len(fit['censored'])} 条")
        L.append(f"  f_ref(c={fit['c_ref']:g}) = {fit['f_ref']:.3g}"
                 + (f"   δ = {fit['delta']:.3g} /十倍频" if fit.get("delta") else ""))
        if fit.get("delta_ci"):
            dlo, dhi = fit["delta_ci"]
            L.append(f"  δ 95%CI [{dlo:.3g}, {dhi:.3g}]   R² = {fit['r2']:.3f}")
        elif fit.get("n_pos", 0) < 3:
            L.append("  ⚠ 正克隆点 <3：δ 的 CI 不可估——预测只有点值，"
                     "下档实测后务必回填")
        for fl in fit.get("flags", []):
            L.append(f"  ⚠ {fl}")
    else:
        L.append("  （规划未带拟合明细）")
    L.append("")
    L.append("二、梯级规划")
    if plan["status"] != "可行":
        L.append(f"  ⛔ {plan['stuck_reason']}")
        L.append("")
        return "\n".join(L)
    L.append(f"  可行最大步长 Δ = {plan['step']:.3g}"
             f"（margin {plan['margin']:g}）→ {plan['n_rungs']} 档："
             f"{plan['c_start']:g} → {plan['c_target']:g}")
    L.append(f"  每档 f_step ≈ {plan['f_rung']:.3g}（平稳补充假设；每档实测后回填）")
    L.append("")
    L.append("  档   c_from→c_to        预测f_step     需CFU    乘客   累计    "
             + ("BGC失守" if plan["bgc_kb"] else ""))
    for r in plan["ladder"]:
        bgc = f"{r['p_bgc']:.2%}" if r["p_bgc"] is not None else ""
        L.append(f"  {r['i']:>2}   {r['c_from']:g}→{r['c_to']:<g}   "
                 f"{r['f_step']:.2e}   {r['need_cfu']:8.3g}  {r['passengers']:.2f}  "
                 f"{r['cum_passengers']:5.2f}   {bgc}")
    L.append("")
    L.append("三、累计账")
    tot = plan["p_total"]
    L.append(f"  总乘客 E[P] = {tot:.2f}（窗口 Σ{plan['n_win'] * plan['m_on']:.2f}"
             f" + 板内 Σ{sum(r['p_plate'] for r in plan['ladder']):.2f}"
             + (f" + 终点后 {plan['n_post'] * plan['m_off']:.2f}"
                if plan["n_post"] > 0 and plan["m_off"] else "") + ")")
    if plan["p_bgc"] is not None:
        L.append(f"  保护区（{plan['bgc_kb']:g} kb）失守 = 1 − exp(−E[P]×L/G) "
                 f"= {plan['p_bgc']:.2%}（总乘客单次代入，不是每档相加）")
    if not plan["m_off_counted"]:
        L.append("  ⚠ 未给 --m-off：窗口外段（板内扩张/终点后）**未计入**——"
                 "本会话已证该项与窗口同量级，务必补测")
    L.append("")
    L.append("四、怎么读")
    L.append("  - 预测锚在 PAP 实测曲线上，平稳补充假设会随机制天花板/谱系效应失效——")
    L.append("    每档实测 f_step 后回填 --pap 重算；预测偏出 CI → 曲线重拟合")
    L.append("  - 富余（surplus）大的档可评估**关诱导**：省 (M_on−M_off)×n_win/档"
             "（先修泄漏，否则\"关\"也要付 M_off/M_base 的地板价）")
    if plan.get("extrapolation_decades") is not None and plan["extrapolation_decades"] > 1:
        L.append(f"  ⚠ 规划外推超出 PAP 实测范围约 {plan['extrapolation_decades']:.1f} 个 δ"
                 "——远端档的预测最脆弱，优先实测")
    return "\n".join(L)


# --------------------------------------------------------------------------- 账本解释器

def run_ledger(records: Sequence[tuple[str, float, str, float, int]], *,
               m_on: float, m_off: float | None, genome_mb: float,
               n_win: float, n_plate: float, n_post: float,
               k: float, capacity: float, c_target: float,
               bgc_kb: float | None, margin: float,
               c_start: float | None = None,
               nt: float | None = None,
               sat_f: float = 0.5,
               tau_off_h: float | None = None,
               td_h: float | None = None,
               propagation: str = "manual",
               cells_per_colony: float | None = None,
               n0: float | None = None) -> dict:
    """账本解释器：records（cycle, c, when, cfu, clones）→ 现状 + 下一窗设计 + 剩余梯级 + 累计账。

    语义（对应 protocol-fstep / protocol-ramp）：
      - 每个 (cycle, c) 的 pre/post 板配对做背景扣除（Newcombe 差值；
        红灯 = f_pre ≥ 0.5×f_post；无 pre 则 f_step = f_post 并注明“未扣背景”）；
      - PAP 曲线：全部 post 行（跨档合并、按 c 池化）拟合一条尾部——
        拟合残差本身就是平稳补充假设的审计（出 CI 的点报警）；
        频率 ≥ sat_f 的饱和行不进拟合（对尾部无信息）[约定]；
        sat_f 是显式参数（默认 0.5 = "长到与未处理相当"），写进报告；
        落在阈值 ±20% 相对带内的点单独标出（分类对阈值敏感）；
      - 当前位置 c_now = 有实测 f_post ≥ 0.5 的最高浓度（--c-start 可覆盖）；
      - 剩余梯级：从 c_now 到 c_target 按可行最大步长均匀分档；
        累计乘客 = 已完成档基数（每档 n_win×M_on + n_plate×M_off + 洗脱尾，[约定]按设计值）
                   + 剩余档同式 + n_post×M_off + 洗脱尾（只计一次）；
        洗脱尾见 `washout_tail`（需 --tau-off-h 与 --td 同时给）；
      - 失守概率 = 1 − exp(−E[P_总] × L_prot/G)（总乘客单次代入）。
    """
    if not records:
        raise ValueError("账本没有任何记录——先完成一个窗口并铺板（pre/post）")
    if not (0.0 < sat_f < 1.0):
        raise ValueError(f"sat_f={sat_f!r} 必须在 (0, 1) 之间——"
                         "它是\"频率多高算饱和\"的阈值")
    cycles = sorted({str(r[0]) for r in records},
                    key=lambda s: (0, int(s)) if s.isdigit() else (1, s))

    by_cc: dict[tuple[str, float], list[tuple[str, float, int]]] = {}
    for cyc, c_plate, w, cfu, clones in records:
        if w not in ("pre", "post"):
            raise ValueError(f"when 必须是 pre/post（当前 {w!r}，cycle {cyc}）")
        by_cc.setdefault((str(cyc), float(c_plate)),
                         []).append((w, float(cfu), int(clones)))

    # ---- ① 逐 (cycle, c) 的 f_step 表 ----
    def _pool(rows: Sequence[tuple[str, float, int]]) -> tuple[int, int]:
        return (sum(n for _, n, _ in rows), sum(k for _, _, k in rows))

    cycle_table = []
    f_by_cc: dict[tuple[str, float], dict] = {}
    for cyc in cycles:
        for (cc, c) in sorted(by_cc, key=lambda t: t[1]):
            if cc != cyc:
                continue
            rows = by_cc[(cc, c)]
            has_pre = any(r[0] == "pre" for r in rows)
            has_post = any(r[0] == "post" for r in rows)
            entry = dict(cycle=cyc, c=c,
                         cfu_pre=_pool([r for r in rows if r[0] == "pre"])[0] if has_pre else 0,
                         clones_pre=_pool([r for r in rows if r[0] == "pre"])[1] if has_pre else 0,
                         cfu_post=_pool([r for r in rows if r[0] == "post"])[0] if has_post else 0,
                         clones_post=_pool([r for r in rows if r[0] == "post"])[1] if has_post else 0)
            if not has_post:
                entry.update(f_step=None, red=False, note="只有 pre 板，无 post")
            else:
                pk, pn = _pool([r for r in rows if r[0] == "post"])
                f_post = wilson_ci(int(pn), pk)          # (克隆数, CFU)
                entry.update(f_post=f_post, saturated=(f_post[0] >= sat_f))
                if has_pre:
                    qk, qn = _pool([r for r in rows if r[0] == "pre"])
                    f_pre = wilson_ci(int(qn), qk)       # (克隆数, CFU)
                    d, dlo, dhi = newcombe_diff(f_post, f_pre)
                    red = (f_pre[0] >= 0.5 * f_post[0]) if f_post[0] > 0 else False
                    entry.update(f_pre=f_pre, f_step=max(0.0, d),
                                 ci=(dlo, dhi), red=red)
                else:
                    entry.update(f_pre=None, f_step=f_post[0],
                                 ci=(f_post[1], f_post[2]), red=False,
                                 note="未扣背景（无 pre 板）")
                f_by_cc[(cyc, c)] = entry
            cycle_table.append(entry)

    # ---- ② 当前位置 + PAP 拟合（平稳性审计）----
    post_pts = [(e["c"], e["clones_post"], e["cfu_post"])
                for e in f_by_cc.values()]
    saturated = [(c, k, n) for c, k, n in post_pts if k / n >= sat_f]
    tail = [(c, k, n) for c, k, n in post_pts if k / n < sat_f]
    # 阈值边界带：分类对 sat_f 敏感的点单独标出（±20% 相对）
    near_sat = [(c, k / n) for c, k, n in post_pts
                if 0.8 * sat_f <= k / n < 1.2 * sat_f]
    if c_start is not None:
        c_now, c_now_source = float(c_start), "--c-start 显式指定"
    else:
        cand = [c for c, k, n in post_pts if k / n >= sat_f]
        if cand:
            c_now, c_now_source = max(cand), f"账本中 f_post ≥ {sat_f:g} 的最高浓度"
        elif tail:
            c_now, c_now_source = min(c for c, _, _ in tail), "最低尾部分档（无饱和点，弱证据）"
        else:
            raise ValueError(
                "无法从账本判定当前位置（没有 f_post ≥ 0.5 的记录）——"
                "请用 --c-start 显式指定，或补一块生长饱和的板")
    if len(tail) < 2:
        raise ValueError(
            "PAP 尾部点不足（<2 个 f<0.5 的浓度点）——无法拟合 δ、无法规划剩余梯级。"
            "下一窗请多涂几块浓度板补点（每块都是下一步的工作板）。")
    fit = fit_pap(tail)

    # ---- ③ 剩余梯级（含已完成档基数）----
    done_cycles = sorted({e["cycle"] for e in f_by_cc.values()},
                         key=lambda s: (0, int(s)) if s.isdigit() else (1, s))
    # ⚠ 与 mtblocks.window_table_f_step 的乘客口径同源——改这里要同步改那里，
    #   并跑 tests/test_ledger_formula_pin.py（docs/multitarget-code-map.md §6）
    p_tail = washout_tail(m_on, m_off, tau_off_h, td_h)
    m_off_eff = m_off if m_off is not None else 0.0

    def _cycle_f(cyc: str) -> float | None:
        """该周期"升档"用的 f：尾部（f < sat_f）里**最大**的那个 f。

        对应"你实际升到的那一档浓度"——它是可数区里最低的浓度、也是 f 最大的尾点。
        不能用"最高浓度的 f"：那一点是纯检测点，f 更小，会把 N₀ 算成负数。
        """
        rows = [e for e in f_by_cc.values()
                if e["cycle"] == cyc and e.get("f_step") is not None
                and 0 < e["f_step"] < sat_f]
        if not rows:
            rows = [e for e in f_by_cc.values()
                    if e["cycle"] == cyc and e.get("f_step")]
        if not rows:
            return None
        return float(max(e["f_step"] for e in rows))

    base, base_detail = 0.0, []
    for cyc in done_cycles:
        n_inter, src = inter_rung_generations(
            mode=propagation, n_plate=n_plate, cells_per_colony=cells_per_colony,
            nt=nt, k=k, n0=n0, f=_cycle_f(cyc), where=f"cycle {cyc}")
        base += n_win * m_on + n_inter * m_off_eff + p_tail
        base_detail.append(dict(cycle=cyc, n_inter=n_inter, source=src))
    plan = plan_ramp(fit, m_on=m_on, k=k, capacity=capacity,
                     c_start=c_now, c_target=c_target, genome_mb=genome_mb,
                     m_off=m_off, n_win=n_win, n_plate=n_plate, n_post=n_post,
                     bgc_kb=bgc_kb, margin=margin, cum_base=base,
                     tau_off_h=tau_off_h, td_h=td_h,
                     propagation=propagation, cells_per_colony=cells_per_colony,
                     nt=nt, n0=n0)
    plan["fit"] = fit
    plan["c_now"] = c_now
    plan["c_now_source"] = c_now_source
    plan["cycles"] = done_cycles
    plan["base_passengers"] = base
    plan["cycle_table"] = cycle_table
    plan["saturated_excluded"] = saturated
    plan["sat_f"] = float(sat_f)
    plan["near_sat"] = near_sat
    plan["p_tail"] = p_tail
    plan["base_detail"] = base_detail
    plan["nt"] = nt

    # ---- ④ 下一窗设计（用最新实测 f_step 优先，曲线预测并列）----
    f_latest = None
    if cycle_table:
        done_entries = [e for e in cycle_table
                        if e["cycle"] == done_cycles[-1] and e.get("f_step") is not None]
        if done_entries:
            f_latest = max(done_entries, key=lambda e: e["c"])
    plan["f_latest"] = f_latest
    return plan


def render_ledger(plan: dict) -> str:
    L: list[str] = []
    cyc_tb = plan["cycle_table"]
    L.append("一、各档 f_step（账本 records）")
    cur = None
    for e in cyc_tb:
        if e["cycle"] != cur:
            cur = e["cycle"]
            L.append(f"  == cycle {cur} ==")
        tag = f"c={e['c']:g} g/L"
        if not e.get("f_step") and e.get("f_step") != 0.0:
            L.append(f"    {tag}: 只有 pre 板，无 post——无法计 f_step")
            continue
        fp = e.get("f_post")
        line = (f"    {tag}: f_post = {e['clones_post']}/{e['cfu_post']:.3g} = {fp[0]:.2e} "
                f"[{fp[1]:.1e}, {fp[2]:.1e}]")
        if e.get("f_pre") is not None:
            fq = e["f_pre"]
            line += (f"  | f_pre = {e['clones_pre']}/{e['cfu_pre']:.3g} = {fq[0]:.2e}  "
                     f"f_step = {e['f_step']:.2e} [{e['ci'][0]:.1e}, {e['ci'][1]:.1e}]")
            if e["red"]:
                line += "  ⛔背景淹没"
        else:
            line += "  | ⚠ 未扣背景（无 pre 板）"
        if e.get("saturated"):
            line += "  [饱和点，未进拟合]"
        L.append(line)
    L.append("")

    L.append("二、PAP 曲线现状（平稳补充假设的审计）")
    fit = plan["fit"]
    L.append(f"  正克隆点 {fit['n_pos']}（饱和剔除 {len(plan['saturated_excluded'])}，"
             f"阈值 f ≥ {plan.get('sat_f', 0.5):g}）"
             f" | f_ref(c={fit['c_ref']:g}) = {fit['f_ref']:.2e}"
             f" | δ = {fit['delta']:.3g} /十倍频")
    if fit.get("delta_ci"):
        dlo, dhi = fit["delta_ci"]
        L.append(f"  δ 95%CI [{dlo:.3g}, {dhi:.3g}]   R² = {fit['r2']:.3f}")
    else:
        L.append("  ⚠ 正克隆点 <3：δ 无 CI——下档多涂浓度板补点")
    for fl in fit.get("flags", []):
        L.append(f"  ⚠ {fl}")
    for c, fr in plan.get("near_sat", []):
        L.append(f"  ⚠ c={c:g} 的 f={fr:.2f} 落在饱和阈值 ±20% 带内——"
                 "换个阈值它就会进出拟合，本档分类不稳，建议补点")
    L.append(f"  当前位置 c_now = {plan['c_now']:g}（{plan['c_now_source']}）")
    L.append("")

    if plan["status"] != "可行":
        L.append("三、剩余梯级")
        L.append(f"  ⛔ {plan['stuck_reason']}")
        L.append("  （保守规划视野到此为止；每档实测回填账本重跑——实际富集后通常优于保守预测）")
        return "\n".join(L)
    L.append("三、下一窗设计（第 " + str(len(plan["cycles"]) + 1) + " 档）")
    r1 = plan["ladder"][0]
    L.append(f"  下一档浓度 c = {r1['c_to']:g}（剩余梯级首档）")
    if plan.get("propagation"):
        L.append(f"  档间扩增：**{plan['propagation']}** 模式 → "
                 f"{r1.get('n_inter', 0):.2f} 代/档（{r1.get('n_inter_source', '')}）")
    L.append(f"  预测 f_step（PAP 曲线）= {r1['f_step']:.3g}"
             + (f" [{r1['f_lo']:.1e}, {r1['f_hi']:.1e}]" if r1["f_lo"] is not None else "（点值，<3 正克隆点无 CI）"))
    if plan["f_latest"] is not None:
        e = plan["f_latest"]
        if e["f_step"] > 0:
            L.append(f"  最新实测（cycle {e['cycle']}, c={e['c']:g}）= {e['f_step']:.3g}——与预测对照，偏出 CI 即重拟合")
        else:
            L.append("  ⚠ 上一档实测 f_step = 0：本档无升档供给（背景淹没或未长）——先查原因再前进")
    if plan.get("nt") is not None:
        n0 = plan["nt"] - plan["k"] / r1["f_step"]
        if n0 > 0:
            L.append(f"  N₀ ≤ Nt_sat − k/f_step = {plan['nt']:.3g} − {plan['k']:g}/"
                     f"{r1['f_step']:.2e} = {n0:.3g}（高密度接种 = 最强降荷旋钮）")
        else:
            L.append(f"  ⛔ 下一窗给不出 k：Nt_sat − k/f = {n0:.3g} < 0"
                     "（f_step 太小 / Nt_sat 太小 / k 太大）")
    else:
        L.append("  （未给 --nt：跳过 N₀ 设计行）")
    L.append("")

    L.append("四、剩余梯级与累计账")
    L.append(f"  可行最大步长 Δ = {plan['step']:.3g}（margin {plan['margin']:g}）→ "
             f"剩余 {plan['n_rungs']} 档：{plan['c_start']:g} → {plan['c_target']:g}"
             f"（当前 c_now = {plan['c_now']:g}）")
    L.append("")
    L.append("  档   c_from→c_to        预测f_step     需CFU    乘客   累计    "
             + ("BGC失守" if plan["bgc_kb"] else ""))
    for r in plan["ladder"]:
        bgc = f"{r['p_bgc']:.2%}" if r["p_bgc"] is not None else ""
        L.append(f"  {r['i']:>2}   {r['c_from']:g}→{r['c_to']:<g}   "
                 f"{r['f_step']:.2e}   {r['need_cfu']:8.3g}  {r['passengers']:.2f}  "
                 f"{r['cum_passengers']:5.2f}   {bgc}")
    L.append("")
    tot = plan["p_total"]
    L.append(f"  总乘客 E[P] = {tot:.2f}（已完成档基数 {plan['base_passengers']:.2f}"
             f" + 剩余梯级 {tot - plan['base_passengers']:.2f}）")
    p_tail = plan.get("p_tail", 0.0)
    if p_tail:
        L.append(f"  其中洗脱尾 = {p_tail:.2f}/档×次（τ_off = {plan['tau_off_h']:g} h，"
                 f"td = {plan['td_h']:g} h；每次洗脱计一次：已完成档 + 每剩余档 + 终点后）")
    elif plan.get("tau_off_h") is None:
        L.append("  ⚠ 未给 --tau-off-h/--td：洗脱尾**未计入**——它与板内段同量级，"
                 "漏掉会低估每档乘客（见 docs/protocol-ramp.md §3.4）")
    if plan["p_bgc"] is not None:
        L.append(f"  保护区（{plan['bgc_kb']:g} kb）失守 = 1 − exp(−E[P]×L/G) "
                 f"= {plan['p_bgc']:.2%}（总乘客单次代入，不是每档相加）")
    if not plan["m_off_counted"]:
        L.append("  ⚠ 未给 --m-off：窗口外段（板内扩张/终点后）**未计入**——"
                 "板内扩张可与窗口同量级，务必补测")
    L.append("")
    L.append("五、怎么读")
    L.append("  - 预测锚在 PAP 实测曲线上；平稳补充假设的审计就是本报告第二节——")
    L.append("    每档实测回填账本重跑，拟合残差出 CI → 曲线重拟合/拆档")
    L.append("  - 富余大的档可评估**关诱导**（省 (M_on−M_off)×n_win/档）——先修泄漏，")
    L.append("    否则“关”也有泄漏地板价")
    L.append("  - 不预测突变身份/CNV；保护区完整性每档抽检（PCR/产物检定）")
    return "\n".join(L)
