"""三旋钮优化：脉冲时长 T_pulse × 诱导剂浓度 c × 接种密度 N₀。

两阶段：
  1) 点估计网格全扫 → 用成功率约束过滤 → 取"每个剂量每个 T 的 E[乘客] 最优"
     作为多样化候选池（**不是**只取全局最低 E 的 N 个——那样 20 个候选会挤在同一个角落，
     第二阶段等于没做，见 review_2026-10-29 §P1-10）；
  2) 对候选池做蒙特卡洛（参数 CI + 克隆间离散）→ 先按"供给稳健性"过滤
     （P(≥k 事件) 的 P5 仍达标），再按 E[乘客] 排序。

约束：P(独立事件数 ≥ k_target) ≥ P_min；T_pulse ≥ τ_on + 0.5 h；N₀ ≤ N₀_max。
无可行解时输出诊断（供给瓶颈在哪：F 不够 / 涂板量不够 / 背景淹没）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .model import ProtocolResult, StrainCard, evaluate, monte_carlo, replace, val


def poisson_sf(lam: float, k: int) -> float:
    """P(Poisson(λ) ≥ k)——至少 k 个独立事件的概率。"""
    if lam <= 0:
        return 1.0 if k <= 0 else 0.0
    if k <= 0:
        return 1.0
    term = math.exp(-lam)
    cdf = term
    for i in range(1, k):
        term *= lam / i
        cdf += term
    return max(1.0 - cdf, 0.0)


@dataclass
class Decision:
    k_target_clones: int = 12          # 想要的独立 hit 数（进初筛的多样性）
    P_success_min: float = 0.95
    T_pulse_range_h: tuple = (3.0, 24.0)
    N0_min: float = 1.0e6
    N0_max: float = 1.0e8
    doses_available: list = field(default_factory=list)   # 空则从模型/F 表生成
    P_protected_max: float = 0.10      # BGC/其他保护区的误伤硬上限（保护优先）
    P_mutator_max: float = 0.05        # mutator 自失活的硬上限
    N0_grid: int = 8
    T_grid: int = 10
    mc_samples: int = 3000
    seed: int = 0
    top_k_mc: int = 20
    top_k_report: int = 5


class CandidateList(list):
    """候选列表（附带"不同误伤预算下的最优解"表与本次运行的网格/设置，供报告用）。"""
    budget_table: list = []
    meta: dict = {}


@dataclass
class Candidate:
    dose_mM: float
    T_pulse_h: float
    N0: float
    point: ProtocolResult
    mc: dict | None = None

    @property
    def E_passengers(self) -> float:
        return self.point.E_passengers


def effective_doses(card: StrainCard, dec: Decision) -> tuple[list, str]:
    """本次要扫的剂量轴 → (doses, 来源说明)。

    优先级：协议卡 doses_available > 诱导剂模型自动布点 > F 表实测点 > 参考剂量单点。
    **配了 inducer 却没给 doses_available 时不再静默退化成单点**——那会让"三旋钮"
    悄悄变成"两旋钮"，而报告还写着"剂量轴：无"。
    """
    if dec.doses_available:
        ds = sorted({float(d) for d in dec.doses_available})
        src = "按协议卡 doses_available"
        if card.inducer is not None and len(ds) == 1:
            src += ("（⚠ 只有 1 个剂量点，虽然配了 inducer 模型——剂量旋钮实际没有被扫；"
                    "删掉 doses_available 可由模型自动布点）")
        return ds, src
    if card.inducer is not None:
        c50 = max(float(card.inducer.ec50_mM), 1e-9)
        ds = sorted({round(float(x), 6) for x in np.geomspace(c50 / 10.0, c50 * 10.0, 7)})
        return ds, (f"由诱导剂模型自动生成 7 点（C50/10–10×C50 = {ds[0]:g}–{ds[-1]:g} mM；"
                    "若你的可用剂量有限，请在协议卡里写 doses_available）")
    if card.F_doses:
        ds = sorted({float(d["dose_mM"]) for d in card.F_doses})
        return ds, "按 F_doses 的实测剂量点"
    return [val(card.ref_dose_mM)], "无（恒定表达）"


def _diversify(cands: list, k: int) -> tuple[list, int]:
    """按剂量分箱、箱内按 T 去重，各取 E 最小的若干个 → 候选池。

    返回 (picked, 池中出现的剂量数)。旧实现直接取全局 E 最小的 top_k，
    实测前 5 名全是同一个剂量、相邻 T，MC 阶段形同虚设。
    """
    by_dose = {}
    for c in cands:
        by_dose.setdefault(round(float(c.dose_mM), 9), []).append(c)
    per_dose = max(2, k // max(len(by_dose), 1)) if k else 2
    picked = []
    for group in by_dose.values():
        group.sort(key=lambda c: (c.point.E_passengers, -c.point.events))
        seen_T = set()
        cnt = 0
        for c in group:
            tkey = round(float(c.T_pulse_h), 6)
            if tkey in seen_T:          # 同一 T 的多个 N₀ 只留最好的一个
                continue
            seen_T.add(tkey)
            picked.append(c)
            cnt += 1
            if cnt >= per_dose:
                break
    picked.sort(key=lambda c: (c.point.E_passengers, -c.point.events))
    if k and len(picked) > k:
        picked = picked[:k]
    return picked, len(by_dose)


def optimize(card: StrainCard, phases: list, dec: Decision) -> tuple[list, list]:
    doses, dose_source = effective_doses(card, dec)
    tau_on = val(card.tau_on_h)
    t_lo = max(dec.T_pulse_range_h[0], tau_on + 0.5)
    t_hi = max(dec.T_pulse_range_h[1], t_lo * 1.5)
    T_grid = np.geomspace(t_lo, t_hi, dec.T_grid)
    N0_grid = np.geomspace(min(dec.N0_min, dec.N0_max), dec.N0_max, dec.N0_grid)

    all_c = []
    for dose in doses:
        for T in T_grid:
            for N0 in N0_grid:
                pt = evaluate(card, phases, dose_mM=dose, T_pulse_h=float(T), N0=float(N0))
                all_c.append(Candidate(float(dose), float(T), float(N0), pt))
    if not all_c:
        return [], _diagnose_failure(card, phases, dec, doses, T_grid, N0_grid, dose_source)

    def _p_succ(c: Candidate) -> float:
        return poisson_sf(c.point.events, dec.k_target_clones)

    cap_prot = getattr(dec, "P_protected_max", None)
    cap_mut = getattr(dec, "P_mutator_max", None)
    has_prot = any(c.point.protected for c in all_c)
    use_cap = has_prot and (cap_prot is not None or cap_mut is not None)

    def _under_cap(c: Candidate) -> bool:
        if not use_cap:
            return True
        if cap_prot is not None and c.point.p_other > cap_prot + 1e-15:
            return False
        if cap_mut is not None and c.point.p_mutator > cap_mut + 1e-15:
            return False
        return True

    event_ok = [c for c in all_c if _p_succ(c) >= dec.P_success_min]
    feasible = [c for c in event_ok if _under_cap(c)]

    # 连"忽略保护上限"都拿不到 k → 是供给问题，不是保护问题
    if not event_ok:
        return [], _diagnose_failure(card, phases, dec, doses, T_grid, N0_grid, dose_source)

    # 保护上限优先（硬约束）：k 要求的窗口不得越过 mutator/保护区允许的上限。
    # 若二者冲突，按上限截断窗口，并报告实际 events 与缺口。
    if not feasible and use_cap:
        under_cap = [c for c in all_c if _under_cap(c)]
        common = dict(
            doses=list(doses), dose_source=dose_source,
            T_range=(float(T_grid[0]), float(T_grid[-1]), int(dec.T_grid)),
            N0_range=(float(N0_grid[0]), float(N0_grid[-1]), int(dec.N0_grid)),
            n_feasible=0, n_pool=1, n_dropped_by_mc=0,
            mc_samples=int(dec.mc_samples), seed=int(dec.seed),
            k=int(dec.k_target_clones), P_min=float(dec.P_success_min),
            protection_binding=True,
            protection_cap=float(cap_prot) if cap_prot is not None else None,
            mutator_cap=float(cap_mut) if cap_mut is not None else None,
        )

        # 哪条上限在"挡路"：看全网格里 events 最大的方案会违反哪条
        _top = max(event_ok, key=lambda c: c.point.events)
        _flags = dict(
            protected_binding=bool(cap_prot is not None and _top.point.p_other > cap_prot + 1e-15),
            mutator_binding=bool(cap_mut is not None and _top.point.p_mutator > cap_mut + 1e-15),
        )

        if under_cap:
            best = max(under_cap, key=lambda c: (c.point.events, -c.point.p_any_protected))
            best.mc = monte_carlo(card, phases, best.dose_mM, best.T_pulse_h, best.N0,
                                  n_samples=dec.mc_samples, seed=dec.seed)
            out = CandidateList([best])
            out.budget_table = _budget_table(under_cap)
            out.meta = dict(common,
                            n_doses_in_pool=len({c.dose_mM for c in under_cap}),
                            min_p_any=min((c.point.p_any_protected for c in under_cap),
                                          default=float("nan")),
                            min_p_mutator=min((c.point.p_mutator for c in under_cap),
                                              default=float("nan")),
                            min_p_other=min((c.point.p_other for c in under_cap),
                                            default=float("nan")),
                            protection_impossible=False,
                            p_mutator=float(best.point.p_mutator),
                            p_other=float(best.point.p_other),
                            actual_events=float(best.point.events),
                            actual_p_success=float(_p_succ(best)),
                            shortfall=max(0.0, float(dec.k_target_clones - best.point.events)),
                            **_flags)
            return out, []
        # 连最短窗口都超上限：没有安全方案；取"最接近上限"的候选
        def _violation(c: Candidate) -> float:
            v = 0.0
            if cap_prot:
                v = max(v, c.point.p_other / cap_prot)
            if cap_mut:
                v = max(v, c.point.p_mutator / cap_mut)
            return v

        best = min(all_c, key=_violation)
        out = CandidateList([best])
        out.meta = dict(common,
                        n_doses_in_pool=len({c.dose_mM for c in all_c}),
                        min_p_any=float(best.point.p_any_protected),
                        min_p_mutator=float(best.point.p_mutator),
                        min_p_other=float(best.point.p_other),
                        protection_impossible=True,
                        p_mutator=float(best.point.p_mutator),
                        p_other=float(best.point.p_other),
                        actual_events=float(best.point.events),
                        actual_p_success=float(_p_succ(best)),
                        shortfall=max(0.0, float(dec.k_target_clones - best.point.events)),
                        **_flags)
        return out, []
    if not feasible:
        return [], _diagnose_failure(card, phases, dec, doses, T_grid, N0_grid, dose_source)

    cands = feasible

    def _primary(c: Candidate) -> float:
        """主目标：配置了保护区就按"保护区期望命中"排，否则按总乘客。"""
        if c.point.protected:
            return sum(r["expected"] for r in c.point.protected
                       if r.get("selection", "silent") == "silent")
        return c.point.E_passengers

    cands.sort(key=lambda c: (_primary(c), c.point.E_passengers, -c.point.events))
    diverse, n_dose = _diversify(cands, dec.top_k_mc)
    for c in diverse:
        c.mc = monte_carlo(card, phases, c.dose_mM, c.T_pulse_h, c.N0,
                           n_samples=dec.mc_samples, seed=dec.seed)
    feasible_mc = [c for c in diverse
                   if c.mc is None or c.mc["p_success_p5"] >= dec.P_success_min - 0.05]
    dropped_by_mc = len(diverse) - len(feasible_mc)
    pool = feasible_mc or diverse
    pool.sort(key=lambda c: (_primary(c),
                             c.mc["E_passengers"][2] if c.mc else c.point.E_passengers))
    out = CandidateList(pool[: dec.top_k_report])
    out.budget_table = _budget_table(pool)
    out.meta = dict(
        doses=list(doses), dose_source=dose_source,
        T_range=(float(T_grid[0]), float(T_grid[-1]), int(dec.T_grid)),
        N0_range=(float(N0_grid[0]), float(N0_grid[-1]), int(dec.N0_grid)),
        n_feasible=len(cands), n_pool=len(diverse), n_doses_in_pool=n_dose,
        n_dropped_by_mc=dropped_by_mc,
        mc_samples=int(dec.mc_samples), seed=int(dec.seed),
        k=int(dec.k_target_clones), P_min=float(dec.P_success_min),
        min_p_any=min((c.point.p_any_protected for c in pool), default=float("nan")),
        min_p_mutator=min((c.point.p_mutator for c in pool), default=float("nan")),
        min_p_other=min((c.point.p_other for c in pool), default=float("nan")),
        protection_binding=False,
        protection_cap=float(cap_prot) if cap_prot is not None else None,
        mutator_cap=float(cap_mut) if cap_mut is not None else None,
        protected_binding=False, mutator_binding=False,
        protection_impossible=False,
    )
    return out, []


def _budget_table(pool: list, budgets=(0.01, 0.02, 0.05, 0.10, 0.20)) -> list:
    """不同"任一保护区被命中概率"预算下，最优可行解是谁。"""
    rows = []
    for b in budgets:
        pick = next((c for c in pool
                     if c.point.protected and c.point.p_any_protected <= b), None)
        rows.append(dict(budget=b, cand=pick))
    return rows


def _diagnose_failure(card: StrainCard, phases: list, dec: Decision,
                      doses, T_grid, N0_grid, dose_source: str = "") -> list:
    """无可行解：找供给瓶颈。"""
    best = None
    for dose in doses:
        for T in T_grid:
            for N0 in N0_grid:
                pt = evaluate(card, phases, dose_mM=dose, T_pulse_h=float(T), N0=float(N0))
                if best is None or pt.events > best[0].events:
                    best = (pt, dose, float(T), float(N0))
    if best is None:
        return ["无可行解且无法评估任何组合：检查协议卡片（缺 pulse 相位或 N0 无效）"]
    pt, dose, T, N0 = best
    need = -math.log(1.0 - dec.P_success_min) if dec.k_target_clones == 1 else None
    has_dose_axis = len({float(d) for d in doses}) > 1
    where = (f"剂量 {dose} mM, T={T:.1f} h, N0={N0:.1e}" if has_dose_axis
             else f"T={T:.1f} h, N0={N0:.1e}")
    msg = [
        f"无可行解：全网格内最大供给 events={pt.events:.2f}（{where}），"
        f"仍达不到 k={dec.k_target_clones} 个独立 hit @ P≥{dec.P_success_min:.2f}。"
    ]
    if dose_source:
        msg.append(f"本次剂量轴来源：{dose_source}")
    if has_dose_axis:
        msg.append("瓶颈排查顺序：① μ_supply 太低 → 提高筛选浓度重测 / 换标记（rsmG 低抗档）；"
                   "② 剂量轴不够 → 测更高诱导剂剂量（参数 7）；"
                   "③ 涂板量 → N_plated 上限是否已用满；"
                   "④ 若背景红灯在前，先解决假抗性再谈供给。")
    else:
        msg.append("无剂量轴（恒定表达）——没有可调的浓度旋钮，排查顺序："
                   "① μ_supply 太低 → 提高筛选浓度重测 / 换标记；"
                   "② N₀/Nt 太小 → 增加总分裂数（更大的接种量或培养体积）；"
                   "③ 涂板量 → N_plated 上限是否已用满；"
                   "④ 若背景红灯在前，先解决假抗性再谈供给。")
    if need:
        msg.append(f"(k=1 时所需 events≈{need:.2f}）")
    return msg
