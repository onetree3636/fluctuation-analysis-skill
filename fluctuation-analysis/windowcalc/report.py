"""Markdown 报告生成——推荐时刻表 / 乘客瀑布 / 蒙特卡洛 / tornado / 红灯 / 待测清单。"""

from __future__ import annotations

import datetime as _dt

from . import __version__
from .model import (StrainCard, avg_activity, params_without_ci,
                    phase_generations, scale_of, val)


def _fmt(x: float, digits: int = 2) -> str:
    if x == 0:
        return "0"
    if 0.01 <= abs(x) < 1000:
        return f"{x:.{digits}f}"
    return f"{x:.2e}"


def _provenance_table(card: StrainCard) -> list[tuple[str, object]]:
    rows = [
        ("M_base（基线，每基因组每代）", card.M_base),
        ("M_on_ref（参考剂量诱导态）", card.M_on_ref),
        ("M_off（off 态泄漏）", card.M_off),
        ("μ_supply（供给端标记速率，选定浓度）", card.mu_supply),
        ("背景假抗性频率", card.background_false),
        ("τ_on（诱导上升沿，h）", card.tau_on_h),
        ("τ_off（洗脱衰变尾，h）★参数5", card.tau_off_h),
    ]
    return rows


def render_optimize_report(card: StrainCard, phases: list, candidates: list,
                           diagnostics: list, decision) -> str:
    lines = []
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    n_measured = 0
    n_total = 0
    prov_lines = []
    for label, p in _provenance_table(card):
        if p is None:
            continue
        n_total += 1
        mark = {"measured": "实测", "literature": "文献", "assumed": "★未测"}.get(p.status, p.status)
        if p.status == "measured":
            n_measured += 1
        ci = f", 95%CI {p.ci[0]:.2e}–{p.ci[1]:.2e}" if p.ci else ""
        prov_lines.append(f"| {label} | {_fmt(p.value)}{ci} | {mark} | {p.note} |")
    for d in card.F_doses:
        st = {"measured": "实测", "literature": "文献", "assumed": "★未测"}.get(d.get("status", "assumed"))
        prov_lines.append(
            f"| F(c) @ {d['dose_mM']} mM | {d['F']}× | {st} | 剂量响应（参数 7） |")

    meta = getattr(candidates, "meta", {}) or {}
    lines.append("# mutator 表达窗口优化报告")
    lines.append("")
    lines.append(f"> 菌株：{card.name}｜引擎 v{__version__}｜生成：{now}")
    doses_used = meta.get("doses") or list(
        getattr(decision, "doses_available", []) or [val(card.ref_dose_mM)])
    has_dose_axis = len({float(d) for d in doses_used}) > 1
    dose_src = meta.get("dose_source")
    if not dose_src:
        if getattr(card, "inducer", None):
            dose_src = "诱导剂模型"
        elif card.F_doses:
            dose_src = "F(c)"
        else:
            dose_src = "无剂量轴（恒定表达）"
    lines.append(f"> 参数实测率：{n_measured}/{n_total}｜剂量轴：{dose_src}")
    base_pulse = next((p for p in phases if p.type == "pulse"), None)
    if base_pulse is not None and candidates:
        b0 = candidates[0]
        d0 = (base_pulse.dose_mM if base_pulse.dose_mM is not None else val(card.ref_dose_mM))
        lines.append(f"> 输入协议：pulse {base_pulse.duration_h:g} h @ {d0:g} mM → "
                     f"推荐 {b0.T_pulse_h:g} h @ {b0.dose_mM:g} mM")
    ps = scale_of(card)
    if ps != 1.0:
        p_sc = card.passenger_scale
        ci_txt = (f"（CI {p_sc.ci[0]:g}–{p_sc.ci[1]:g}）"
                  if hasattr(p_sc, "ci") and p_sc.ci else "")
        lines.append(f"> 绝对负荷：M 按单位点速率填入，已乘 `passenger_scale` = **{ps:g}**{ci_txt}"
                     "——该因子是假设值（需 WGS 校准），绝对乘客数/误伤概率的可信度取决于它；"
                     "**窗口推荐与该因子无关**（尺度不变）。")
    if candidates and not diagnostics:
        lines.append("")
        lines.append(f"> **怎么读**：① 第 1 节选方案（E[乘客] 越小越好，P(≥k) 要不低于 "
                     f"{decision.P_success_min:.0%}）；② 第 6 节看有没有红灯，有就先解决；"
                     "③ 第 4/5/7 节是「还缺什么参数、先测哪个」；④ 附录 A 是复现信息。")
    lines.append("")
    if diagnostics:
        lines.append("## ⛔ 无可行解")
        lines.append("")
        for d in diagnostics:
            lines.append(f"- {d}")
        lines.append("")
    if not candidates:
        return "\n".join(lines) + "\n"

    best = candidates[0]
    lines.append("## 1. 推荐方案")
    lines.append("")
    has_prot = bool(candidates and getattr(candidates[0].point, "protected", None))
    tail = " P(误伤) |" if has_prot else ""
    tail_sep = "---|" if has_prot else ""
    if has_dose_axis:
        lines.append("| # | 剂量 (mM) | T_pulse (h) | N₀ (CFU) | events | P(≥k) | E[乘客/克隆] |" + tail + " MC P5–P95（仅参数） |")
        lines.append("|---|---|---|---|---|---|---|" + tail_sep + "---|")
    else:
        lines.append("| # | T_pulse (h) | N₀ (CFU) | events | P(≥k) | E[乘客/克隆] |" + tail + " MC P5–P95（仅参数） |")
        lines.append("|---|---|---|---|---|---|" + tail_sep + "---|")
    for i, c in enumerate(candidates, 1):
        mc = c.mc
        mc_s = f"{mc['E_passengers'][0]:.1f}–{mc['E_passengers'][2]:.1f}" if mc else "—"
        ps = "≈100%" if c.point.p_success >= 0.9995 else f"{c.point.p_success:.1%}"
        if mc:
            ps = f"{ps}（MC P5 {mc['p_success_p5']:.1%}）"
        prot_s = f" {c.point.p_any_protected:.1%} |" if has_prot else ""
        if has_dose_axis:
            lines.append(
                f"| {i} | {c.dose_mM} | {c.T_pulse_h:.1f} | {c.N0:.1e} | "
                f"{c.point.events:.1f} | {ps} | {c.point.E_passengers:.1f} |" + prot_s + f" {mc_s} |")
        else:
            lines.append(
                f"| {i} | {c.T_pulse_h:.1f} | {c.N0:.1e} | "
                f"{c.point.events:.1f} | {ps} | {c.point.E_passengers:.1f} |" + prot_s + f" {mc_s} |")
    lines.append("")

    # 边界解（P1-9）：两个旋钮都顶在网格端点时，必须说出来，否则读者以为这是内点最优
    bounds = []
    if has_dose_axis and meta.get("doses"):
        if best.dose_mM >= max(doses_used) - 1e-9:
            bounds.append(f"剂量顶到网格上限 {max(doses_used):g} mM")
        elif best.dose_mM <= min(doses_used) + 1e-9:
            bounds.append(f"剂量顶到网格下限 {min(doses_used):g} mM")
    if meta.get("T_range"):
        _, t_hi = meta["T_range"][0], meta["T_range"][1]
        if best.T_pulse_h >= t_hi - 1e-9:
            bounds.append(f"T_pulse 顶到上限 {t_hi:g} h")
        elif best.T_pulse_h <= meta["T_range"][0] + 1e-9:
            bounds.append(f"T_pulse 顶到下限 {meta['T_range'][0]:g} h（可能是 τ_on+0.5 的约束，而非你的设计意图）")
    if meta.get("N0_range"):
        n_lo, n_hi = meta["N0_range"][0], meta["N0_range"][1]
        if best.N0 >= n_hi * (1 - 1e-9):
            bounds.append(f"N₀ 顶到上限 {n_hi:.1e}")
        elif best.N0 <= n_lo * (1 + 1e-9):
            bounds.append(f"N₀ 顶到下限 {n_lo:.1e}")
    n_total_gen = sum(p.n for p in best.point.phases)
    baseline = n_total_gen * val(card.M_base)
    occ = best.point.E_passengers
    if occ < 0.05:
        lines.append(f"> 参考量级：按你卡的 M_base 算，「不开 mutator」的基线 ≈ {baseline:.2g} 个乘客/克隆；"
                     f"本方案 {occ:.2g} ≈ 基线——等于没开 mutator，先查 M_on 的量纲/单位。")
    else:
        lines.append(f"> 参考量级：按你卡的 M_base 算，「不开 mutator」的基线 ≈ {baseline:.2g} 个乘客/克隆；"
                     f"全强度 M7 一轮 20–27 代 ≈ **66–96**（gkab583 WGS）——本方案 {occ:.1f} ≈ 全强度的 {occ / 81:.0%}。")
    lines.append("")
    if meta.get("protection_binding"):
        cap_p = meta.get("protection_cap")
        cap_m = meta.get("mutator_cap")
        binds = []
        if meta.get("protected_binding") and cap_p is not None:
            binds.append(f"BGC/其他保护区上限 P≤{cap_p:.1%}")
        if meta.get("mutator_binding") and cap_m is not None:
            binds.append(f"mutator 自失活上限 P≤{cap_m:.1%}")
        if not binds:
            binds.append("保护上限")
        if meta.get("protection_impossible"):
            lines.append(f"> ⛔ **保护上限不可行**：即使最短窗口也超过 {'；'.join(binds)}"
                         f"（网格内最低：mutator {meta.get('min_p_mutator', float('nan')):.1%}、"
                         f"其他 {meta.get('min_p_other', float('nan')):.1%}）。"
                         "当前没有安全方案——先缩短必须分裂的阶段 / 降低 M_off / 换更安全的构型，"
                         "而不是加大窗口。")
        else:
            lines.append(f"> ⛔ **保护上限 binding（硬约束优先）**：目标 k={meta.get('k')} 需要的窗口会被 "
                         f"{'；'.join(binds)} 截断。已按上限截断窗口："
                         f"实际 E[events]={meta.get('actual_events', float('nan')):.1f}"
                         f"（缺口 {meta.get('shortfall', float('nan')):.1f}，"
                         f"P(≥k)={meta.get('actual_p_success', float('nan')):.1%}）。"
                         "**不要为了 k 超窗**；补缺口：降 k / 提高 N₀ / 换更大选择靶标 / 改用局部编辑。")
        lines.append("")
    if bounds:
        lines.append("> ⚠ **边界解**：" + "；".join(bounds)
                     + "。真实最优点可能在网格外——放宽 `T_pulse_range_h` / `doses_available` / "
                       "`N0_max` 再算一次；若放宽后仍顶在边界，说明是**约束**（k 与 P_min）决定了答案。")
        lines.append("")

    if has_prot:
        lines.append("## 1b. 保护区账本（最优解）")
        lines.append("")
        lines.append("> 假设突变在总 DNA 上均匀分布：`f = 拷贝数×长度/总 DNA`，"
                     "`P(≥1) = 1−exp(−乘客×f)`。热点/偏好不在模型内。")
        lines.append("")
        lines.append("| 保护区 | 角色 | 拷贝 | 占比 f | 期望命中 | P(≥1 拷贝) | P(全部拷贝) | 性质 |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for r in best.point.protected:
            sel = "静默（计入克隆质量）" if r.get("selection", "silent") == "silent" \
                else "选择淘汰（计入产量）"
            lines.append(f"| {r['name']} | {r.get('role', 'other')} | {r['copies']:g} | {r['f']:.2%} | "
                         f"{r['expected']:.3f} | {r['prob']:.1%} | {r.get('prob_all_copies', 0):.1e} | {sel} |")
        lam_mut = sum(r['expected'] for r in best.point.protected
                      if r.get('selection', 'silent') == 'silent' and r.get('role') == 'mutator')
        lam_other = sum(r['expected'] for r in best.point.protected
                        if r.get('selection', 'silent') == 'silent' and r.get('role') != 'mutator')
        lines.append(f"| **mutator 自失活（role=mutator 合计）** | mutator | | | {lam_mut:.3f} "
                     f"| **{best.point.p_mutator:.1%}** | | |")
        lines.append(f"| **BGC/其他误伤（role≠mutator 合计）** | bgc/other | | | {lam_other:.3f} "
                     f"| **{best.point.p_other:.1%}** | | |")
        lines.append(f"| **克隆质量风险（静默区合计）** | | | | {lam_mut + lam_other:.3f} "
                     f"| **{best.point.p_any_protected:.1%}** | | |")
        lines.append("")
        bt = getattr(candidates, "budget_table", None)
        if bt:
            lines.append("## 1c. 不同误伤预算下的可选方案")
            lines.append("")
            lines.append("| 误伤预算 P(任一保护区) | 可达方案（T_pulse / 剂量 / N₀） | events | E[乘客] |")
            lines.append("|---|---|---|---|")
            for row in bt:
                c = row["cand"]
                if c is None:
                    lines.append(
                        f"| ≤{row['budget']:.0%} | 候选池里没有这么低负荷的方案"
                        f"（池内最低 {meta.get('min_p_any', float('nan')):.1%}）"
                        "——不是「拿不到 k 个事件」，而是满足 k 的方案负荷都超预算 | — | — |")
                else:
                    lines.append(
                        f"| ≤{row['budget']:.0%} | T={c.T_pulse_h:.1f} h / {c.dose_mM:g} mM / "
                        f"N₀={c.N0:.1e} | {c.point.events:.1f} | {c.point.E_passengers:.1f} |")
            lines.append("")

    lines.append("## 2. 时刻表（最优解逐相位分解）")
    lines.append("")
    lines.append("| 相位 | 类型 | 时长 (h) | ā（平均活性） | n（每世系代数） | M_eff（/基因组/代） | 乘客贡献 |")
    lines.append("|---|---|---|---|---|---|---|")
    missing_dur = False
    for pr in best.point.phases:
        if pr.duration_h > 0:
            dur = f"{pr.duration_h:.1f}"
        else:
            dur, missing_dur = "—", True
        lines.append(
            f"| {pr.name} | {pr.type} | {dur} | {pr.a_bar:.2f} | {pr.n:.1f} | "
            f"{_fmt(pr.M_eff, 3)} | {pr.passengers:.2f} |")
    lines.append(f"| **合计** | | | | | | **{best.point.E_passengers:.1f}** |")
    lines.append("")
    if missing_dur:
        lines.append("> 「—」= 协议卡片只给了代数没给时长；换算式 h = n × td（td 见菌株卡的 td_h_on/off）。")
        lines.append("")
    if scale_of(card) != 1.0:
        lines.append(f"> 本节 M_eff 与乘客贡献已乘 `passenger_scale` = {scale_of(card):g}"
                     "（即「折算后」的绝对口径，不是卡片里填的原始数字）。")
        lines.append("")

    ind = getattr(card, "inducer", None)
    if ind is not None:
        from .induction import hill, psi, psi_max, saturation_dose
        sp = ind
        d0 = best.dose_mM
        h0 = float(hill(d0, sp.ec50_mM, sp.n_hill))
        rel0 = float(psi(h0, sp.mixing, sp.q)) / psi_max(sp.mixing, sp.q)
        c95 = saturation_dose(0.95, sp.ec50_mM, sp.n_hill, sp.mixing, sp.q)
        lines.append("## 2b. 剂量轴（诱导剂机制模型）")
        lines.append("")
        lines.append(f"- 模型：{sp.mixing}；C50 = {sp.ec50_mM:g} mM，n = {sp.n_hill:.2f}，"
                     f"饱和幅度 Δsat = {sp.amplitude * psi_max(sp.mixing, sp.q):.3g}"
                     "（有效占位幅度，不是占比 p）")
        lines.append(f"- 参考剂量 {val(card.ref_dose_mM):g} mM、选定剂量 {d0:g} mM；"
                     f"95% 饱和剂量 {c95:.3g} mM → 当前剂量处于 **{rel0:.0%}** 饱和"
                     + ("（已饱和：再加诱导剂只增负荷）" if rel0 >= 0.95 else ""))
        lines.append("")
        lines.append("| 剂量 (mM) | h | ψ/ψ(1) | M(c)（每基因组每代） |")
        lines.append("|---|---|---|---|")
        for c in sorted({val(card.ref_dose_mM), d0, c95}):
            h = float(hill(c, sp.ec50_mM, sp.n_hill))
            rel = float(psi(h, sp.mixing, sp.q)) / psi_max(sp.mixing, sp.q)
            m = card.M_off.value + sp.amplitude * float(psi(h, sp.mixing, sp.q))
            lines.append(f"| {c:g} | {h:.3f} | {rel:.2f} | {m:.4g} |")
        lines.append("")

    lines.append("## 3. 本方案的标度律读数")
    lines.append("")
    div_total = sum(p.divisions for p in best.point.phases)
    off_load = sum(p.passengers for p in best.point.phases if p.type != "pulse")
    pulse_load = sum(p.passengers for p in best.point.phases if p.type == "pulse")
    lines.append(f"- 供给 ∝ 窗口内总分裂：本方案 ΔN = {div_total:.2e}（涂板量即供给上限）；")
    if max(pulse_load, off_load) < 1e-3:
        # 两个数都≈0 时不能说"谁主导"——那多半是 M 的量纲/单位错了
        lines.append(f"- 负荷几乎为 0（脉冲项 {pulse_load:.2g}、off 项 {off_load:.2g}）："
                     "先查 M 的量级/单位（全基因组每代 vs 单位点 μ），此刻讨论「谁主导」没有意义。")
    else:
        lines.append(f"- 负荷 ∝ 每世系代数：脉冲项 {pulse_load:.1f} vs off 泄漏项 {off_load:.1f}——"
                     + ("off 项不可忽略，静态公式 E[m]=n·M 低估了负荷。" if off_load > 0.2 * pulse_load
                        else "脉冲项主导，静态公式近似可用。"))
    lines.append(f"- 成功率按独立事件计：events={best.point.events:.1f} ≈ 大约能筛到 "
                 f"**{best.point.events:.0f} 个不同来源**的克隆（不是菌落数：一个 jackpot 克隆"
                 f"贡献几百个菌落，也只算 1 个事件）→ P(≥{decision.k_target_clones} 独立 hit)="
                 f"{best.point.p_success:.1%}。")
    div_pulse = sum(p.divisions for p in best.point.phases if p.type == "pulse")
    div_other = div_total - div_pulse
    if div_other <= 1e-12 * max(div_total, 1.0) and best.point.events > 0:
        ratio = best.point.E_passengers / best.point.events
        lines.append(
            f"- ⚠️ **无开关窗口**：所有分裂都在诱变期内，每事件乘客 ≈ {ratio:.3g}"
            "（= 每克隆乘客 ÷ 独立事件数 = (M/μ_supply)·(n/ΔN) ≈ M/(μ_supply·Nt)；**不是** M/μ_supply）。"
            "此时「优化」只剩一条路：在满足 k 的前提下把总分裂数压到最低，无法靠调窗口改善比值。")
    lines.append("")

    if best.mc:
        lines.append("## 4. 蒙特卡洛（两套区间，别混着读）")
        lines.append("")
        e5, e50, e95 = best.mc["E_passengers"]
        c5, c50, c95 = best.mc.get("E_passengers_clone", best.mc["E_passengers"])
        v5, v50, v95 = best.mc["events"]
        lines.append(f"- **只含参数不确定性**（与第 1 节点估计同口径）："
                     f"每克隆乘客 P5 {e5:.1f}｜P50 {e50:.1f}｜P95 {e95:.1f}")
        sigma = best.mc.get("dispersion_sigma_ln", 0.0)
        if sigma > 0:
            lines.append(f"- **再叠加克隆间离散**（σ_ln={sigma:g}，生物变异——参数测得再准也不会收窄）："
                         f"P5 {c5:.1f}｜P95 {c95:.1f}")
        lines.append(f"- events：P5 {v5:.1f}｜P50 {v50:.1f}｜P95 {v95:.1f}")
        lines.append(f"- P(≥{decision.k_target_clones} hit) 的 P5：{best.mc['p_success_p5']:.1%}")
        miss = params_without_ci(card)
        if miss:
            lines.append("- ⚠ 以下参数**没有 CI**：RateParam 类按 ±2× 对数正态假设（不是测量值），"
                         "inducer.* 与 td/Nt_saturation 按常数处理："
                         + "、".join(miss))
        lines.append("")

    lines.append("## 5. 敏感性（tornado，CI 端到端）")
    lines.append("")
    from .model import tornado as _tornado
    rows = _tornado(card, phases, best.dose_mM, best.T_pulse_h, best.N0)
    if rows:
        lines.append("| 参数 | E[P] @CI下界 | E[P] @CI上界 | events 变化 |")
        lines.append("|---|---|---|---|")
        for r in rows[:6]:
            lines.append(
                f"| {r['param']} | {r['E_low']:.1f} | {r['E_high']:.1f} | "
                f"{r['ev_low']:.1f} ↔ {r['ev_high']:.1f} |")
        lines.append("")

    flags = best.point.flags
    if flags:
        lines.append("## 6. 红灯与警示")
        lines.append("")
        for f in flags:
            lines.append(f"- {f}")
        lines.append("")

    lines.append("## 7. 参数实测状态（★ = 软件正在用假设值，结果可信度受限）")
    lines.append("")
    lines.append("| 参数 | 值 | 状态 | 备注 |")
    lines.append("|---|---|---|---|")
    for pl in prov_lines:
        lines.append(pl)
    lines.append("")
    lines.append("**待测优先级（按第 5 节 tornado 的 CI 端到端变化排序，不是事先规定）**")
    lines.append("")
    by_e = [r for r in rows if r.get("range_E", 0.0) > 0]
    by_ev = sorted([r for r in rows if r.get("range_ev", 0.0) > 0],
                   key=lambda r: r["range_ev"], reverse=True)
    if by_e:
        lines.append("- 对**负荷**（乘客预算）影响最大：" + " > ".join(
            f"{r['param']}（{r['E_low']:.1f}→{r['E_high']:.1f}）" for r in by_e[:3]))
    if by_ev:
        lines.append("- 对**供给**（能否拿到 k 个独立克隆）影响最大：" + " > ".join(
            f"{r['param']}（{r['ev_low']:.1f}→{r['ev_high']:.1f}）" for r in by_ev[:3]))
    if not by_e and not by_ev:
        lines.append("- 没有任何参数带 CI，tornado 为空——先给参数补 `ci: [lo, hi]` 才能排序。")
    else:
        lines.append("- 判据：负荷决定乘客预算、供给决定克隆数，谁先卡住就先测谁（见第 1/3 节）。")
    lines.append("")

    # 附录 A：复现信息（P2-7）——没有这一段，报告换个 --mc-samples 就不可复原
    lines.append("## 附录 A. 复现信息")
    lines.append("")
    lines.append("| 项 | 值 |")
    lines.append("|---|---|")
    lines.append(f"| 引擎版本 | windowcalc v{__version__} |")
    lines.append(f"| 生成时间 | {now} |")
    lines.append(f"| 菌株名 | {card.name} |")
    if meta:
        tr = meta.get("T_range")
        nr = meta.get("N0_range")
        lines.append(f"| 剂量网格 | {len(meta.get('doses', []))} 点："
                     + "、".join(f"{d:g}" for d in meta.get("doses", [])) + " mM"
                     f"（{meta.get('dose_source', '')}） |")
        if tr:
            lines.append(f"| T_pulse 网格 | {tr[0]:.2f}–{tr[1]:.2f} h，{tr[2]} 点（几何布点） |")
        if nr:
            lines.append(f"| N₀ 网格 | {nr[0]:.2e}–{nr[1]:.2e} CFU，{nr[2]} 点 |")
        lines.append(f"| 可行候选 / 候选池 / MC 剔除 | {meta.get('n_feasible')} / "
                     f"{meta.get('n_pool')}（覆盖 {meta.get('n_doses_in_pool')} 个剂量） / "
                     f"{meta.get('n_dropped_by_mc')} |")
        lines.append(f"| MC 抽样 | {meta.get('mc_samples')} 次，seed={meta.get('seed')} |")
        lines.append(f"| 目标 | k={meta.get('k')} 个独立事件，P_min={meta.get('P_min'):.2f} |")
    lines.append(f"| 结构参数 | td_on={val(card.td_h_on):g} h，td_off={val(card.td_h_off):g} h，"
                 f"Nt_saturation={val(card.Nt_saturation):.2e} CFU，"
                 f"ref_dose={val(card.ref_dose_mM):g} mM |")
    if getattr(card, "inducer", None):
        sp = card.inducer
        lines.append(f"| 诱导剂模型 | {sp.name}：{sp.mixing}，A={sp.amplitude:.4g}，"
                     f"C50={sp.ec50_mM:.4g} mM，n={sp.n_hill:.3g}，q={sp.q:.3g}"
                     + (f"，CI={sp.ci}" if sp.ci else "，无 CI") + " |")
    if card.passenger_scale is not None:
        p_sc = card.passenger_scale
        lines.append(f"| 乘客折算因子 | passenger_scale={val(p_sc):g}"
                     + (f"，CI {p_sc.ci[0]:g}–{p_sc.ci[1]:g}" if hasattr(p_sc, "ci") and p_sc.ci else "")
                     + "（只影响绝对负荷，不影响推荐） |")
    lines.append(f"| 克隆间离散 | clone_dispersion_sigma_ln="
                 f"{card.clone_dispersion_sigma_ln:g} |")
    lines.append("")
    lines.append("> 复现：把本节的值填回 strain_card / protocol_card（或直接保留卡片文件），"
                 "用同一 `--mc-samples` 与 seed 重跑即可得到同一份报告。")
    lines.append("")
    return "\n".join(lines) + "\n"
