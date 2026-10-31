"""相位积分引擎——供给（独立事件数）/ 负荷（每克隆期望乘客）/ 成功率。

模型（README §2）：
    开关活性 a(t)：0→1 上升沿（τ_on，线性近似），1→e^(−t/τ_off) 衰减
      M(t) = M_off + (M_on(c) − M_off)·a(t)        每基因组每代乘客突变数
      M_on(c) = M_base × F(c)
    供给  events = Σ_相位 μ_supply(c)·ā_i·ΔN_i           ΔN = 相位内总分裂（≈ Nt − N₀）
    负荷  E[P]  = Σ_相位 n_i · M(ā_i)                n = 每世系代数
    成功  P = 1 − exp(−events)                       events 按独立事件计（jackpot 是假供给）

★ 单一事实源（2026-10-29）：**单个相位的 (ā, n, M_eff, 乘客, 供给) 只在 `phase_terms()` 里算一次**，
  `evaluate`（点估计）与 `monte_carlo`（逐样本）都调它。此前两处各写一遍，导致
  monte_carlo 漏掉 decay 相位的活性、供给只算 pulse 相位——见 review_2026-10-29_professional.md §1。

两条标度律（引擎必须原样呈现给用户）：
    供给 ∝ 窗口内总分裂次数；负荷 ∝ 每世系在窗口内经历的代数 n = log2(Nt/N₀)。

诚实简化（README §4.2）：菌丝体口径随波动分析（per CFU per passage）；表达负担经
td_on/td_off 之差进入；衰变尾窗口内默认不分裂（generations 可显式覆盖）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

PHASE_TYPES = ("pre", "pulse", "decay", "off")


# ----------------------------------------------------------------------------- 输入结构

@dataclass
class RateParam:
    """带 95% CI 与出处的标量参数；CI 用于蒙特卡洛与 tornado。"""
    value: float
    ci: tuple[float, float] = None
    status: str = "assumed"        # measured | assumed | literature
    note: str = ""

    @property
    def sigma_ln(self) -> float:
        """对数正态 σ：由非对称 CI 反推（取两侧均值），无 CI 时默认 2 倍。

        ⚠ 无 CI 时这个 2 倍是**假设**，不是测量——报告 §4/§7 会把它标出来
        （见 `params_without_ci()`），别把它当成已测得的不确定度。
        """
        if self.ci and self.value > 0 and self.ci[0] > 0 and self.ci[1] > self.ci[0]:
            up = math.log(self.ci[1] / self.value) / 1.96
            dn = math.log(self.value / self.ci[0]) / 1.96
            return max((up + dn) / 2.0, 1e-6)
        return math.log(2.0) / 1.96


@dataclass
class InducerSpec:
    """诱导剂剂量轴的机制模型参数（见 induction.py）。

    ec50_mM 字段存的是 **C50（速率半饱和剂量）**——不是启动子的 EC50；
    两者关系 C50 = EC50·(1−q)^(1/n)，q→1 时才近似相等。YAML 里推荐写 `c50_mM`（旧名也认）。

    只有 amplitude（= ΔM·q，有效占位幅度）是稳健可辨识的；q 仅在 two_core 下有弱辨识。

    ci：可选 {参数名: (lo, hi)}——`{"amplitude": (...), "c50": (...), "n": (...)}`。
    把 `fit-induction` 的 bootstrap 区间填进来，MC 才会传播剂量轴的不确定性。
    """
    name: str = "诱导剂"
    amplitude: float = 0.0
    ec50_mM: float = 0.3      # C50：速率半饱和剂量（mM）
    n_hill: float = 2.0
    mixing: str = "linear"          # linear | two_core | both_core
    q: float = 1.0
    ci: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.mixing not in ("linear", "two_core", "both_core"):
            raise ValueError(f"inducer.mixing={self.mixing!r} 不在 linear/two_core/both_core")
        if self.amplitude < 0:
            raise ValueError("inducer.amplitude 不能为负")
        if self.ec50_mM <= 0:
            raise ValueError(f"inducer.C50 必须为正（当前 {self.ec50_mM}）")
        if self.n_hill <= 0:
            raise ValueError(f"inducer.n_hill 必须为正（当前 {self.n_hill}）")
        if not (0 < self.q <= 1):
            raise ValueError(f"inducer.q 必须在 (0,1]（当前 {self.q}）")


@dataclass
class StrainCard:
    """菌株卡——全部速率/动力学参数（windowcalc rate 子命令的输出直接粘贴）。

    td_h_on / td_h_off / Nt_saturation / ref_dose_mM 既接受裸数值，也接受
    `{value, ci, status}`（与 M_* 同一套写法）——它们决定代数 n，是负荷的一阶驱动量，
    能带 CI 才谈得上"传播不确定性"。
    """
    name: str = "strain"
    M_base: RateParam = None            # 每基因组每代，基线
    M_on_ref: RateParam = None        # 参考剂量下的每基因组每代（=M_base×F_ref）
    M_off: RateParam = None           # off 态泄漏
    mu_supply: RateParam = None           # 供给端标记速率（每分裂，选定浓度实测）
    mu_supply_source: str | None = None   # μ_supply 是用哪个标记/药物测的（元数据，用于外推提示）
    selection_marker: str | None = None   # 窗口实际用哪个标记筛选（元数据）
    background_false: RateParam = None  # 孢子直接涂选择板的假抗性频率
    tau_on_h: RateParam = None
    tau_off_h: RateParam = None       # ★ 关键待测
    F_doses: list = field(default_factory=list)   # [{"dose_mM":0.3,"F":100,"status":"measured"}]
    ref_dose_mM: float = 0.3
    td_h_off: float = 2.0
    td_h_on: float = 2.5
    Nt_saturation: float = 3.0e8      # CFU（单管/单板当量）
    clone_dispersion_sigma_ln: float = 1.0  # gkab583：克隆间频率差 15×
    passenger_scale: RateParam | float | None = None  # μ(单位点)→M(全基因组) 的折算因子
    inducer: InducerSpec | None = None      # 有则用机制模型（表达×竞争占位）替代 F 表
    genome_mb: float | None = None          # 基因组大小（Mb）；不填则不做保护区分析
    plasmids: list = field(default_factory=list)   # [{name, kb, copies}]：用于总 DNA 折算
    protected: list = field(default_factory=list)  # [{name, kb, copies, note}]：不能误伤的区域


def val(x) -> float:
    """RateParam 或裸数值 → float（结构参数也允许携带 CI）。"""
    return x.value if isinstance(x, RateParam) else float(x)


def scale_of(card) -> float:
    """乘客折算因子 passenger_scale（缺省 1.0）。

    作用：把卡片里 M_base/M_on/M_off 的**输入口径**换算到"全基因组每代"的绝对负荷。
    只影响乘客数与误伤概率，**不影响**窗口推荐（目标函数与约束对该因子是尺度不变的）、
    也不影响供给 events（供给与 M 无关）与 M_off/M_base 的红灯比值。
    """
    x = getattr(card, "passenger_scale", None)
    return 1.0 if x is None else val(x)


def ci_of(x):
    """RateParam / 裸数值 → CI 或 None。"""
    return x.ci if isinstance(x, RateParam) else None


@dataclass
class Phase:
    name: str
    type: str = "off"                 # pre | pulse | decay | off
    duration_h: float = 0.0
    generations: float | None = None  # 显式代数（优先）
    N0: float | None = None
    Nt: float | None = None
    dose_mM: float | None = None      # pulse 相位的诱导剂浓度
    td_h: float | None = None         # 覆盖 strain card
    a_bar_override: float | None = None  # 显式平均活性（回代文献锚点用，如 gkab583）
    plasmid_loss_per_gen: float = 0.0    # 温敏质粒每代丢失/稀释比例（0.5=每代减半；0=不衰减）

    def __post_init__(self):
        if self.type not in PHASE_TYPES:
            raise ValueError(f"相位类型 {self.type!r} 不在 {PHASE_TYPES}")


# ----------------------------------------------------------------------------- 活性函数

def _avg_geometric_decay(rate: float, n):
    """(1−rate)^x 在 x∈[0,n] 上的平均（= 相位内平均剩余质粒比例）。

    n 可以是数组（td 逐样本采样时 n 也是数组）。
    """
    if rate is None or rate <= 0:
        return 1.0
    r = min(max(float(rate), 0.0), 0.999999)
    if r <= 1e-9:
        return 1.0
    nn = np.asarray(n, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(nn > 0, ((1.0 - r) ** nn - 1.0) / (nn * math.log(1.0 - r)), 1.0)
    return float(out) if np.ndim(out) == 0 else out


def _avg_activity_vec(phase: Phase, tau_on, tau_off):
    """`avg_activity` 的向量化内核（tau_on/tau_off 可为逐样本数组）。

    pulse：ā = max(T − τ_on/2, 0)/T
    decay：ā = τ_off·(1−e^(−T/τ_off))/T
    """
    if phase.a_bar_override is not None:
        return float(phase.a_bar_override)
    if phase.type in ("pre", "off"):
        return 0.0
    T = max(phase.duration_h, 0.0)
    if T == 0:
        return 0.0 if phase.type == "pulse" else 1.0
    if phase.type == "pulse":
        return np.maximum(T - np.asarray(tau_on, dtype=float) / 2.0, 0.0) / T
    tau = np.asarray(tau_off, dtype=float)
    with np.errstate(over="ignore"):
        return tau * (1.0 - np.exp(-T / tau)) / T


def avg_activity(phase: Phase, tau_on: float, tau_off: float) -> float:
    """相位内平均开关活性 ā ∈ [0,1]（标量版；数组版见 `phase_terms`）。"""
    return float(_avg_activity_vec(phase, tau_on, tau_off))


def _generations_vec(phase: Phase, td_off, td_on, nt_sat):
    """`phase_generations` 的向量化内核（td / nt_sat 可为数组）。"""
    if phase.generations is not None:
        return float(phase.generations)
    td = phase.td_h if phase.td_h else (td_on if phase.type == "pulse" else td_off)
    td = np.asarray(td, dtype=float)
    if phase.N0 is not None:
        if phase.Nt is not None:
            n_raw = math.log2(phase.Nt / phase.N0) if phase.Nt > phase.N0 else 0.0
        elif phase.duration_h and phase.duration_h > 0:
            n_raw = phase.duration_h / td
        else:
            n_raw = 0.0
        if nt_sat is not None and phase.N0 > 0:
            reached = np.minimum(phase.N0 * np.power(2.0, n_raw), np.asarray(nt_sat, float))
            return np.log2(reached / phase.N0)
        return n_raw
    if phase.duration_h and phase.duration_h > 0:
        return phase.duration_h / td
    return 0.0


def phase_generations(phase: Phase, td_off: float, td_on: float,
                      Nt_sat: float | None = None) -> float:
    """相位内每世系经历的代数 n。

    优先级：显式 generations > N0/Nt（含饱和封顶） > duration/td。
    只有 N0（无 Nt）时按 duration/td 增长并用 Nt_sat 封顶——与 divisions 同一口径。
    """
    return float(_generations_vec(phase, td_off, td_on, Nt_sat))


def _divisions_vec(phase: Phase, td_on, td_off, nt_sat):
    """`phase_divisions` 的向量化内核。"""
    if phase.N0 is None:
        return 0.0
    td = phase.td_h if phase.td_h else (td_on if phase.type == "pulse" else td_off)
    td = np.asarray(td, dtype=float)
    if phase.generations is not None:
        n_reached = float(phase.generations)
    elif phase.Nt is not None:
        n_reached = math.log2(phase.Nt / phase.N0) if phase.Nt > phase.N0 else 0.0
    elif phase.duration_h and phase.duration_h > 0:
        n_reached = phase.duration_h / td
    else:
        n_reached = 0.0
    reached = np.minimum(phase.N0 * np.power(2.0, n_reached), np.asarray(nt_sat, float))
    return np.maximum(reached - phase.N0, 0.0)


def phase_divisions(phase: Phase, td_on: float, Nt_sat: float,
                    td_off: float | None = None) -> float:
    """相位内总分裂次数 ΔN（只对给了 N0 的相位有意义；饱和封顶）。

    td_off 缺省时与旧行为一致（非 pulse 相位也用 td_on）；显式给出则按相位类型选。
    """
    return float(_divisions_vec(phase, td_on, td_on if td_off is None else td_off, Nt_sat))


def phase_terms(phase: Phase, *, td_on, td_off, nt_sat, tau_on, tau_off,
                m_on, m_off, activity=1.0) -> dict:
    """★ 单个相位的全部量——**负荷与供给的唯一实现**。

    `evaluate` 传标量、`monte_carlo` 传逐样本数组（numpy 广播），两者不可能再算错口径。

    返回 dict：
        a_bar       相位内平均活性 ā
        n           每世系代数
        leak        该相位的有效泄漏（含质粒稀释衰减）
        M_eff       有效速率 = leak + (M_on − leak)·ā
        passengers  n·M_eff（每克隆期望乘客）
        divisions   ΔN（相位内总分裂，只对给了 N0 的相位非零）
        supply      activity·ā·ΔN（**不含 μ_supply**，由调用方乘）
    """
    abar = _avg_activity_vec(phase, tau_on, tau_off)
    n = _generations_vec(phase, td_off, td_on, nt_sat)
    leak = m_off
    if phase.plasmid_loss_per_gen and phase.plasmid_loss_per_gen > 0:
        leak = m_off * _avg_geometric_decay(phase.plasmid_loss_per_gen, n)
    m_eff = leak + (m_on - leak) * abar
    div = _divisions_vec(phase, td_on, td_off, nt_sat)
    return dict(a_bar=abar, n=n, leak=leak, M_eff=m_eff,
                passengers=n * m_eff, divisions=div,
                supply=activity * abar * div)


def F_at(card: StrainCard, dose_mM: float) -> tuple[float, str]:
    """剂量 c → 倍数 F(c)。无剂量响应数据时任何剂量都返回参考 F（单旋钮退化）。"""
    m_on_ref, m_base = val(card.M_on_ref), val(card.M_base)
    if not card.F_doses:
        return (m_on_ref / m_base if m_base > 0 else 1.0), "参考剂量（无剂量响应数据）"
    doses = sorted(card.F_doses, key=lambda d: d["dose_mM"])
    for d in doses:
        if abs(d["dose_mM"] - dose_mM) < 1e-9:
            return float(d["F"]), d.get("status", "measured")
    # 线性内插（对数剂量轴）；超出范围按最近端点（不外推）
    xs = [d["dose_mM"] for d in doses]
    ys = [float(d["F"]) for d in doses]
    x = min(max(dose_mM, xs[0]), xs[-1])
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            t = (math.log(x / xs[i])) / (math.log(xs[i + 1] / xs[i])) if xs[i] > 0 else 0.0
            return ys[i] + t * (ys[i + 1] - ys[i]), "interpolated"
    return ys[-1], "edge"


def M_at(card: "StrainCard", dose_mM: float) -> tuple[float, float, str]:
    """剂量 → (M(c), activity, status)。

    - 有 `inducer`：M(c) = M_off + A·ψ(h(c))，activity = ψ(h(c))/ψ(h(ref))（供给按同一活性缩放）；
    - 无 `inducer`：沿用 F 表（M(c) = M_on_ref·F/F(ref)，activity 同 F/F(ref)），行为与旧版一致。
    """
    if card.inducer is not None:
        from .induction import hill, psi

        sp = card.inducer
        h = float(hill(dose_mM, sp.ec50_mM, sp.n_hill))
        h_ref = float(hill(val(card.ref_dose_mM), sp.ec50_mM, sp.n_hill))
        w = float(psi(h, sp.mixing, sp.q))
        w_ref = float(psi(h_ref, sp.mixing, sp.q))
        activity = (w / w_ref) if w_ref > 0 else 1.0
        return val(card.M_off) + sp.amplitude * w, activity, "inducer-model"

    F, status = F_at(card, dose_mM)
    F_ref = F_at(card, val(card.ref_dose_mM))[0] or 1.0
    shape = F / F_ref
    return val(card.M_on_ref) * shape, shape, status


# ----------------------------------------------------------------------------- 核心计算

@dataclass
class PhaseResult:
    name: str
    type: str
    a_bar: float
    n: float                      # 每世系代数
    M_eff: float                  # 该相位有效突变率（每基因组每代）
    passengers: float             # 该相位负荷贡献
    divisions: float              # ΔN
    events: float                 # 该相位供给贡献
    duration_h: float = 0.0       # 相位时长（报告时刻表用）
    td_used: float | None = None


@dataclass
class ProtocolResult:
    phases: list
    events: float
    E_passengers: float
    p_success: float
    flags: list = field(default_factory=list)
    protected: list = field(default_factory=list)   # [{name, f, expected, prob, role}]
    p_any_protected: float = 0.0                     # 任一保护区被打中的概率（全部 silent）
    p_mutator: float = 0.0                           # mutator 自失活概率（role=mutator 的 silent 区）
    p_other: float = 0.0                             # 非 mutator 保护区误伤概率（BGC/其他 silent 区）


def evaluate(card: StrainCard, phases: list, dose_mM: float | None = None,
             T_pulse_h: float | None = None, N0: float | None = None) -> ProtocolResult:
    """对一份（可被决策变量覆盖的）协议做点估计。

    覆盖规则：dose_mM/T_pulse_h/N0 只替换第一个 pulse 相位的对应字段。

    **decay 相位沿用触发它的脉冲的 (M_on, activity)**——去诱导那一刻的速率是脉冲末端的
    速率，不是参考剂量处的速率；否则脉冲跑在 3 mM 而衰变尾按参考剂量算，会把尾巴
    低估一个 activity 因子。
    """
    phases = [replace(p) for p in phases]
    pulse_seen = 0
    for p in phases:
        if p.type == "pulse":
            pulse_seen += 1
            if dose_mM is not None and pulse_seen == 1:
                p.dose_mM = dose_mM
            if T_pulse_h is not None and pulse_seen == 1:
                p.duration_h = T_pulse_h
            if N0 is not None and pulse_seen == 1:
                p.N0 = N0
    if pulse_seen == 0:
        raise ValueError("协议中没有 pulse 相位")

    M_base, M_off, mu_supply = val(card.M_base), val(card.M_off), val(card.mu_supply)
    ps = scale_of(card)                 # 乘客折算因子（只作用在负荷侧）
    tau_on, tau_off = val(card.tau_on_h), val(card.tau_off_h)
    td_on, td_off, nt_sat = val(card.td_h_on), val(card.td_h_off), val(card.Nt_saturation)
    flags = []

    results, events_total, load_total = [], 0.0, 0.0
    # 一致性检查：参考剂量处的模型/表值与直接测量的 M_on_ref 应相符
    if card.inducer is not None:
        m_ref_model = M_at(card, val(card.ref_dose_mM))[0]
        if val(card.M_on_ref) > 0 and abs(m_ref_model - val(card.M_on_ref)) > 0.2 * max(
                m_ref_model, val(card.M_on_ref), 1e-12):
            flags_check = (f"inducer 模型在参考剂量给 M={m_ref_model:.4g}，与 "
                           f"M_on_ref={val(card.M_on_ref):.4g} 差 >20%——两组数据应出自同批测定")
        else:
            flags_check = None
    elif card.F_doses and val(card.M_base) > 0:
        f_shape_ref = F_at(card, val(card.ref_dose_mM))[0] or 1.0
        f_rate_ref = val(card.M_on_ref) / val(card.M_base)
        if abs(f_shape_ref - f_rate_ref) > 0.2 * max(f_shape_ref, f_rate_ref, 1e-12):
            flags_check = (f"F_doses@{val(card.ref_dose_mM)} mM（{f_shape_ref:.0f}×）与 "
                           f"M_on_ref/M_base（{f_rate_ref:.0f}×）差 >20%——两组数据应出自同批测定")
        else:
            flags_check = None
    else:
        flags_check = None

    last_on = None          # 最近一个 pulse 的 (M_on, activity)
    for p in phases:
        if p.type == "decay" and last_on is not None:
            m_on, activity = last_on
            F_status = "decay（沿用脉冲末端活性）"
        else:
            m_on, activity, F_status = M_at(
                card, p.dose_mM if p.dose_mM is not None else val(card.ref_dose_mM))
            if p.type == "pulse":
                last_on = (m_on, activity)
        t = phase_terms(p, td_on=td_on, td_off=td_off, nt_sat=nt_sat,
                        tau_on=tau_on, tau_off=tau_off,
                        m_on=m_on * ps, m_off=M_off * ps, activity=activity)
        ev = mu_supply * float(t["supply"])
        events_total += ev
        load_total += float(t["passengers"])
        results.append(PhaseResult(p.name, p.type, float(t["a_bar"]), float(t["n"]),
                                   float(t["M_eff"]), float(t["passengers"]),
                                   float(t["divisions"]), ev, float(p.duration_h),
                                   p.td_h or (td_on if p.type == "pulse" else td_off)))
        if p.type == "pulse" and F_status == "参考剂量（无剂量响应数据）" and dose_mM is not None \
                and abs(dose_mM - val(card.ref_dose_mM)) > 1e-9:
            flags.append(f"剂量 {dose_mM} mM 无 F(c) 数据，按参考 F 计算（参数 7 未测）")

    # 红灯（README §4 / 规程 §6.4.9）
    if flags_check:
        flags.append("警示：" + flags_check)
    if M_off > 3.0 * M_base:
        flags.append(f"红灯：M_off={M_off:.4g} > 3×M_base={3*M_base:.4g}（泄漏超标，规程 §6.4.9 合格线）")
    if card.background_false and val(card.background_false) > 0:
        # 涂板量 = 脉冲结束时的细胞数（N0 + ΔN）
        N_plated = sum(p.N0 + phase_divisions(p, td_on, nt_sat, td_off)
                       for p in phases if p.type == "pulse" and p.N0 is not None)
        if N_plated > 0 and val(card.background_false) * N_plated > events_total:
            flags.append(
                f"红灯：背景假抗性（b·N_plated={val(card.background_false)*N_plated:.1f}）"
                f"将淹没真供给（events={events_total:.1f}）——该药/该浓度不可用")
    if tau_off > 6.0:
        flags.append(f"警示：τ_off={tau_off:.1f} h 偏长，有效窗口显著大于名义脉冲（参数 5 待实测收紧）")
    tiny = [n for n in ("M_base", "M_on_ref", "M_off")
            if getattr(card, n) is not None and 0 < val(getattr(card, n)) < M_LOCI_LIKE]
    if tiny:
        if ps > 1.0 + 1e-12:
            flags.append(f"警示：{'、'.join(tiny)} 是按单位点速率（μ）填的，已按 "
                         f"passenger_scale={ps:g} 折算成绝对负荷——该因子是假设值，"
                         "绝对乘客数/误伤概率的可信度取决于它（不是 WGS 校准值就当区间读）")
        else:
            flags.append(f"红灯：{'、'.join(tiny)} 的数值看起来是**单位点**速率（μ），不是全基因组每代——"
                         "乘客会被低估到近似 0。若确实是 μ，请加 `passenger_scale`（10⁵–10⁶，"
                         "用 WGS 校准）折算后再看绝对风险")

    prot, p_any = protected_hits(card, load_total)
    lam_mut = sum(r["expected"] for r in prot
                  if r["selection"] == "silent" and r.get("role") == "mutator")
    lam_other = sum(r["expected"] for r in prot
                    if r["selection"] == "silent" and r.get("role") != "mutator")
    return ProtocolResult(results, events_total, load_total,
                          1.0 - math.exp(-events_total), flags, prot, p_any,
                          p_mutator=1.0 - math.exp(-lam_mut),
                          p_other=1.0 - math.exp(-lam_other))


def protected_hits(card: "StrainCard", passengers: float) -> tuple[list, float]:
    """把"总乘客"折算成各保护区的期望命中数与命中概率。

    假设突变在总 DNA 上均匀分布：f = 拷贝数×长度 / 总 DNA；
    E[hits] = 乘客 × f；P(≥1) = 1 − exp(−E[hits])（泊松近似）。
    未配置 genome_mb 或 protected 时返回 ([], 0.0)。
    """
    if not card.genome_mb or not card.protected:
        return [], 0.0
    total_kb = float(card.genome_mb) * 1000.0
    for pl in card.plasmids:
        total_kb += float(pl.get("kb", 0.0)) * float(pl.get("copies", 1.0))
    if total_kb <= 0:
        return [], 0.0
    out, lam_silent = [], 0.0
    for r in card.protected:
        copies = max(float(r.get("copies", 1.0)), 1.0)
        f = float(r["kb"]) * copies / total_kb
        lam = max(passengers, 0.0) * f
        per_copy = lam / copies
        # 至少一个拷贝被打中 vs 全部拷贝被打中（后者才可能功能丧失）
        p_any = 1.0 - math.exp(-lam)
        p_all = (1.0 - math.exp(-per_copy)) ** copies
        sel = str(r.get("selection", "silent"))
        out.append(dict(name=str(r.get("name", "region")), f=f, copies=copies,
                        selection=sel, expected=lam, prob=p_any, prob_all_copies=p_all,
                        role=str(r.get("role", "other"))))
        if sel == "silent":
            lam_silent += lam
    return out, 1.0 - math.exp(-lam_silent)


# ----------------------------------------------------------------------------- 参数清单与校验

#: 结构参数（不是 RateParam，但同样决定结果）：字段名 → 说明
STRUCT_PARAMS = {
    "td_h_on": "诱导态倍增时间（h）",
    "td_h_off": "off 态倍增时间（h）",
    "Nt_saturation": "饱和密度（CFU/管）",
    "ref_dose_mM": "参考剂量（mM）",
}


#: M 的量级提示：低于这个数几乎只能是"单位点速率"，不可能是全基因组每代
M_LOCI_LIKE = 1e-5
_M_SCALE_HINT = ("模型的 M 是**全基因组每代**的乘客数（示例卡 M_base=0.004 是 Drake 量级）；"
                 "Rif/St 波动分析给的是**单位点**速率 μ（~1e-9）。两者差 10³–10⁶ 倍，"
                 "直接填会让乘客负荷被低估到近似为 0——只有比值（如 M_off/M_base）可以跨量纲比较。")


def card_warnings(card: "StrainCard") -> list:
    """不阻断计算、但会让结论失真的"量级/单位"警告（报告与 stderr 都会显示）。"""
    out = []
    tiny = []
    for name in ("M_base", "M_on_ref", "M_off"):
        x = getattr(card, name, None)
        if x is None:
            continue
        v = val(x)
        if 0 < v < M_LOCI_LIKE:
            tiny.append(f"{name}={v:.3g}")
    if tiny:
        ps = scale_of(card)
        if ps > 1.0 + 1e-12:
            out.append("、".join(tiny) + f" 是单位点速率（μ），已按 passenger_scale={ps:g} 折算；"
                       "该因子要拿自家 WGS 对子校准（同一株：测 Rif μ + WGS 数乘客，两者之商）")
        else:
            out.append("、".join(tiny) + " 看起来是**单位点**速率（μ）。" + _M_SCALE_HINT)
    return out


def params_without_ci(card: "StrainCard") -> list:
    """列出"没有 CI、MC 只能按 ±2 倍（或按常数）处理"的参数——报告必须标出来。"""
    out = []
    for name in ("M_base", "M_on_ref", "M_off", "mu_supply", "tau_on_h", "tau_off_h"):
        p = getattr(card, name, None)
        if p is None:
            out.append(name)
        elif isinstance(p, RateParam):
            if p.ci is None:
                out.append(name)
        else:
            out.append(name)
    if card.inducer is not None:
        have = card.inducer.ci or {}
        for name in ("amplitude", "c50", "n"):
            if name not in have:
                out.append(f"inducer.{name}")
    return out


def _check_range(problems, label, x, lo=None, hi=None):
    try:
        v = val(x)
    except (TypeError, ValueError):
        problems.append(f"{label} 不是数值（读到 {x!r}）")
        return
    if not math.isfinite(v):
        problems.append(f"{label} 不是有限数（{v}）")
        return
    if lo is not None and v <= lo:
        problems.append(f"{label} = {v:g}，必须 > {lo:g}")
    if hi is not None and v > hi:
        problems.append(f"{label} = {v:g}，必须 ≤ {hi:g}")


def card_problems(card: "StrainCard", decision=None) -> list:
    """菌株卡/协议卡的**合法性**检查（不只检查"字段在不在"）。

    返回人类可读的问题列表；`cli` 有则报错并一次列全，**不允许带着负的 td 一路算下去**
    （否则会以"μ_supply 太低"的假诊断出现在报告里）。
    """
    p = []
    for name, lo in (("M_base", 0.0), ("M_on_ref", 0.0), ("mu_supply", 0.0),
                     ("tau_on_h", 0.0), ("tau_off_h", 0.0)):
        x = getattr(card, name, None)
        if x is None:
            p.append(f"{name} 缺失")
        else:
            _check_range(p, name, x, lo=lo)
    x = getattr(card, "M_off", None)
    if x is None:
        p.append("M_off 缺失")
    else:
        _check_range(p, "M_off", x, lo=-1e-300)      # 允许 0（完全关不掉的泄漏才 >0）
    for name in ("td_h_on", "td_h_off", "Nt_saturation", "ref_dose_mM"):
        x = getattr(card, name, None)
        if x is None:
            p.append(f"{name} 缺失")
        else:
            _check_range(p, name, x, lo=0.0)
    if card.M_on_ref is not None and card.M_off is not None:
        if val(card.M_off) > val(card.M_on_ref) * 1.000001:
            p.append(f"M_off={val(card.M_off):g} > M_on_ref={val(card.M_on_ref):g}"
                     "——off 泄漏不可能高于诱导态，检查两列是否写反")
    if card.passenger_scale is not None:
        _check_range(p, "passenger_scale", card.passenger_scale, lo=0.0)
    if card.clone_dispersion_sigma_ln is not None and card.clone_dispersion_sigma_ln < 0:
        p.append(f"clone_dispersion_sigma_ln = {card.clone_dispersion_sigma_ln:g}，必须 ≥ 0")
    if card.genome_mb is not None:
        _check_range(p, "genome_mb", card.genome_mb, lo=0.0)
    for i, pl in enumerate(card.plasmids or []):
        _check_range(p, f"plasmids[{i}].kb", pl.get("kb"), lo=0.0)
        _check_range(p, f"plasmids[{i}].copies", pl.get("copies", 1.0), lo=0.0)
    for i, r in enumerate(card.protected or []):
        nm = r.get("name", f"protected[{i}]")
        if r.get("kb") is None:
            p.append(f"{nm} 缺 kb")
        else:
            _check_range(p, f"{nm}.kb", r.get("kb"), lo=0.0)
        _check_range(p, f"{nm}.copies", r.get("copies", 1.0), lo=0.0)
        if str(r.get("selection", "silent")) not in ("silent", "lethal"):
            p.append(f"{nm}.selection={r.get('selection')!r} 不在 silent/lethal")
        role = str(r.get("role", "other"))
        if role not in ("mutator", "bgc", "other"):
            p.append(f"{nm}.role={role!r} 不在 mutator/bgc/other")
    for i, fd in enumerate(card.F_doses or []):
        _check_range(p, f"F_doses[{i}].dose_mM", fd.get("dose_mM"), lo=-1e-300)
        _check_range(p, f"F_doses[{i}].F", fd.get("F"), lo=0.0)
    if card.inducer is not None:
        sp = card.inducer
        _check_range(p, "inducer.amplitude", sp.amplitude, lo=-1e-300)
        _check_range(p, "inducer.c50_mM", sp.ec50_mM, lo=0.0)
        _check_range(p, "inducer.n_hill", sp.n_hill, lo=0.0)
        if not (0 < sp.q <= 1):
            p.append(f"inducer.q = {sp.q:g}，必须在 (0,1]")
    if decision is not None:
        if decision.k_target_clones < 1:
            p.append(f"k_target_clones = {decision.k_target_clones}，必须 ≥ 1")
        if not (0 < decision.P_success_min < 1):
            p.append(f"P_success_min = {decision.P_success_min:g}，必须在 (0,1)")
        cap = getattr(decision, "P_protected_max", None)
        if cap is not None and not (0 < cap <= 1):
            p.append(f"P_protected_max = {cap:g}，必须在 (0,1]")
        cap_mut = getattr(decision, "P_mutator_max", None)
        if cap_mut is not None and not (0 < cap_mut <= 1):
            p.append(f"P_mutator_max = {cap_mut:g}，必须在 (0,1]")
        lo, hi = decision.T_pulse_range_h
        if not (0 < lo <= hi):
            p.append(f"T_pulse_range_h = [{lo:g}, {hi:g}]，必须满足 0 < 下限 ≤ 上限")
        if not (0 < decision.N0_min <= decision.N0_max):
            p.append(f"N0_min/N0_max = {decision.N0_min:g}/{decision.N0_max:g}，"
                     "必须满足 0 < 下限 ≤ 上限")
        for d in decision.doses_available or []:
            if d <= 0:
                p.append(f"doses_available 含 ≤0 的剂量（{d:g} mM）")
    return p


def phase_problems(phases: list) -> list:
    """协议相位的基本合法性（负时长、负代数、质粒丢失率越界等）。"""
    p = []
    if not any(x.type == "pulse" for x in phases):
        p.append("协议里没有 pulse 相位")
    for ph in phases:
        if ph.duration_h < 0:
            p.append(f"相位 {ph.name!r} 的 duration_h = {ph.duration_h:g}，不能为负")
        if ph.generations is not None and ph.generations < 0:
            p.append(f"相位 {ph.name!r} 的 generations = {ph.generations:g}，不能为负")
        for k in ("N0", "Nt"):
            v = getattr(ph, k)
            if v is not None and v < 0:
                p.append(f"相位 {ph.name!r} 的 {k} = {v:g}，不能为负")
        if ph.dose_mM is not None and ph.dose_mM < 0:
            p.append(f"相位 {ph.name!r} 的 dose_mM = {ph.dose_mM:g}，不能为负")
        if not (0.0 <= ph.plasmid_loss_per_gen < 1.0):
            p.append(f"相位 {ph.name!r} 的 plasmid_loss_per_gen = {ph.plasmid_loss_per_gen:g}，"
                     "必须在 [0,1)（1 表示一代全丢，用 0.999 近似）")
        if ph.td_h is not None and ph.td_h <= 0:
            p.append(f"相位 {ph.name!r} 的 td_h = {ph.td_h:g}，必须 > 0")
        if ph.a_bar_override is not None and not (0 <= ph.a_bar_override <= 1):
            p.append(f"相位 {ph.name!r} 的 a_bar_override = {ph.a_bar_override:g}，必须在 [0,1]")
    return p


# ----------------------------------------------------------------------------- 蒙特卡洛

def _sample_rate(rng: np.random.Generator, p: RateParam, n: int) -> np.ndarray:
    return p.value * np.exp(rng.normal(0.0, p.sigma_ln, n))


def _sample_ci(rng, value: float, ci, n: int) -> np.ndarray:
    """有 CI 按对数正态采样；无 CI 时返回常数（**不**制造 ±2 倍）。"""
    v = float(value)
    if ci is None or v <= 0 or ci[0] <= 0 or ci[1] <= ci[0]:
        return np.full(n, v)
    up = math.log(ci[1] / v) / 1.96
    dn = math.log(v / ci[0]) / 1.96
    return v * np.exp(rng.normal(0.0, max((up + dn) / 2.0, 1e-6), n))


def _sample_struct(rng, x, n: int) -> np.ndarray:
    """结构参数（td / Nt_saturation / ref_dose）：带 CI 才采样，否则广播常数。"""
    return _sample_ci(rng, val(x), ci_of(x), n)


def _dose_terms(card: "StrainCard", dose: float, n: int, rng,
                s_M_off: np.ndarray, s_M_on_ref: np.ndarray):
    """逐样本 (M_on, activity)。剂量轴参数带 CI 时会一起采样。"""
    if card.inducer is not None:
        from .induction import hill, psi

        sp = card.inducer
        A = _sample_ci(rng, sp.amplitude, (sp.ci or {}).get("amplitude"), n)
        c50 = _sample_ci(rng, sp.ec50_mM, (sp.ci or {}).get("c50"), n)
        nh = _sample_ci(rng, sp.n_hill, (sp.ci or {}).get("n"), n)
        h = hill(dose, c50, nh)
        h_ref = hill(val(card.ref_dose_mM), c50, nh)
        w = psi(h, sp.mixing, sp.q)
        w_ref = psi(h_ref, sp.mixing, sp.q)
        activity = np.where(w_ref > 0, w / np.maximum(w_ref, 1e-300), 1.0)
        return s_M_off + A * w, activity

    F, _ = F_at(card, dose)
    F_ref = F_at(card, val(card.ref_dose_mM))[0] or 1.0
    shape = F / F_ref
    return s_M_on_ref * shape, np.full(n, shape)


def monte_carlo(card: StrainCard, phases: list, dose_mM: float, T_pulse_h: float,
                N0: float, n_samples: int = 4000, seed: int = 0) -> dict:
    """参数不确定性 + 克隆间离散 → **两套分位数，分开给**。

    - `E_passengers`：只含参数不确定性（与点估计同口径，可直接和 `evaluate` 比）；
    - `E_passengers_clone`：再乘克隆间离散 exp(N(0, σ_ln))（生物变异，测量再准也不收窄）；
    - `events` / `p_success`：供给侧（不乘克隆离散——events 是培养物级计数）。
    """
    rng = np.random.default_rng(seed)
    n = n_samples
    s_M_off = _sample_rate(rng, card.M_off, n)
    s_M_on_ref = _sample_rate(rng, card.M_on_ref, n)
    s_mu_supply = _sample_rate(rng, card.mu_supply, n)
    s_tau_on = _sample_rate(rng, card.tau_on_h, n)
    s_tau_off = _sample_rate(rng, card.tau_off_h, n)
    s_scale = _sample_struct(rng, card.passenger_scale if card.passenger_scale is not None else 1.0, n)
    s_td_on = _sample_struct(rng, card.td_h_on, n)
    s_td_off = _sample_struct(rng, card.td_h_off, n)
    s_nsat = _sample_struct(rng, card.Nt_saturation, n)

    pulse_idx = next(i for i, p in enumerate(phases) if p.type == "pulse")
    pp = replace(phases[pulse_idx], dose_mM=dose_mM, duration_h=T_pulse_h, N0=N0)

    m_pulse, act_pulse = _dose_terms(card, dose_mM, n, rng, s_M_off, s_M_on_ref)
    m_ref, act_ref = _dose_terms(card, val(card.ref_dose_mM), n, rng, s_M_off, s_M_on_ref)

    load = np.zeros(n)
    supply = np.zeros(n)
    for i, p in enumerate(phases):
        q = pp if i == pulse_idx else p
        if p.type in ("pulse", "decay"):
            m_on, act = m_pulse, act_pulse
        else:
            m_on, act = m_ref, act_ref
        t = phase_terms(q, td_on=s_td_on, td_off=s_td_off, nt_sat=s_nsat,
                        tau_on=s_tau_on, tau_off=s_tau_off,
                        m_on=m_on * s_scale, m_off=s_M_off * s_scale, activity=act)
        load = load + t["passengers"]
        supply = supply + t["supply"]

    events = s_mu_supply * supply
    disp = np.exp(rng.normal(0.0, card.clone_dispersion_sigma_ln, n))
    p_succ = 1.0 - np.exp(-events)

    def _pct(a):
        return (float(np.percentile(a, 5)), float(np.percentile(a, 50)),
                float(np.percentile(a, 95)))

    return {
        "E_passengers": _pct(load),
        "E_passengers_clone": _pct(load * disp),
        "events": _pct(events),
        "p_success": _pct(p_succ),
        "p_success_p5": float(np.percentile(p_succ, 5)),
        "dispersion_sigma_ln": float(card.clone_dispersion_sigma_ln),
        "n_samples": int(n),
    }


def tornado(card: StrainCard, phases: list, dose_mM: float, T_pulse_h: float,
            N0: float) -> list:
    """单参数敏感性：各参数取 CI 下/上界时 E[P] 与 events 的变化（点估计引擎）。

    只列**这个模型真的在用**的参数——否则会给出"零敏感度"的假象：
    inducer 模型下 M_on_ref / F_doses 不参与计算，它们不会出现在表里。
    """
    rows = []

    def add(label, ci, apply_fn):
        row = {"param": label}
        for tag, v in (("low", ci[0]), ("high", ci[1])):
            restore = apply_fn(v)
            try:
                r = evaluate(card, phases, dose_mM, T_pulse_h, N0)
            finally:
                restore()
            row[f"E_{tag}"] = r.E_passengers
            row[f"ev_{tag}"] = r.events
        row["range_E"] = abs(row["E_high"] - row["E_low"])
        row["range_ev"] = abs(row["ev_high"] - row["ev_low"])
        rows.append(row)

    for label in ("M_off", "mu_supply", "tau_on_h", "tau_off_h", "M_on_ref"):
        p = getattr(card, label)
        if p is None or ci_of(p) is None:
            continue
        if label == "M_on_ref" and card.inducer is not None:
            continue                      # 该模型不用 M_on_ref
        add(label, p.ci, lambda v, p=p: _set_attr(p, "value", v))

    if card.inducer is not None:
        sp = card.inducer
        for key, attr, label in (("amplitude", "amplitude", "inducer.amplitude"),
                                 ("c50", "ec50_mM", "inducer.c50_mM"),
                                 ("n", "n_hill", "inducer.n_hill")):
            ci = (sp.ci or {}).get(key)
            if ci is None:
                continue
            add(label, ci, lambda v, sp=sp, attr=attr: _set_attr(sp, attr, v))

    if card.passenger_scale is not None and ci_of(card.passenger_scale) is not None:
        add("passenger_scale", ci_of(card.passenger_scale),
            lambda v: _set_card(card, "passenger_scale", v))

    for label in STRUCT_PARAMS:
        obj = getattr(card, label, None)
        if obj is None or ci_of(obj) is None:
            continue
        add(label, ci_of(obj), lambda v, label=label: _set_card(card, label, v))

    rows.sort(key=lambda r: r["range_E"], reverse=True)
    return rows


def _set_attr(obj, name, v):
    """临时设置对象属性 name=v，返回恢复函数（RateParam.value 与 InducerSpec.* 通用）。"""
    saved = getattr(obj, name)
    setattr(obj, name, v)
    return lambda: setattr(obj, name, saved)


def _set_card(card, label, v):
    """临时替换菌株卡上的结构参数（可能是 RateParam 也可能是裸数值），返回恢复函数。"""
    saved = getattr(card, label)
    setattr(card, label, v)
    return lambda: setattr(card, label, saved)
