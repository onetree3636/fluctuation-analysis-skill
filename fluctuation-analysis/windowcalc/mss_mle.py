"""MSS-MLE（Ma–Sandri–Sarkar 最大似然）突变率估计——内置实现。

出处：
- 递推式：Ma, Sandri & Sarkar, J Appl Prob 1992;29:255–267；
  以 Rosche & Foster, Methods 2000;20:4–17 的 Table 2 Eq.[16] 转写：
      p0 = exp(-m)
      p_r = (m/r) * sum_{i=0}^{r-1} p_i / (r - i + 1)
  递推给出每个 p_r 的**精确值**（复合泊松表示），不做尾部折叠、不做重归一——
  与 FALCOR 的实现口径一致（对观测 r 逐点取 p_r 进似然）。
  已用复合泊松一阶原理核对 p1、p2（严格相等），并以该文 Table 3 全样本
  （52 管，MSS-MLE m=1.11，95% CL 0.77–1.49）验收——见 tests/test_mss_mle.py。
- 速率定义：mu = m / mean(Nt)（FALCOR 原文口径，Hall 2009）。
- 置信区间：profile likelihood（Delta_LL = chi2_{1,0.95}/2 = 1.9207）。
  FALCOR 用 Stewart 1994 的 ln(m) 正态近似（经 Foster 2006）——小 m 下该近似
  比 LR 区间宽（R&F Exp.2：Stewart 0.77 vs 本实现 0.81），最终以 FALCOR 网页
  对数为准（浏览器可用，curl 403）。
- 部分涂板：二项稀释似然——每个突变细胞以概率 z=涂板比例被观测，
  P_obs(r) = Σ_{j≥r} p_j·C(j,r)·z^r·(1-z)^(j-r)。（R&F 的"先 r/z 再 ML"对
  有缺口的观测直方图有偏，弃用。）
  与 rSalvador `newton.LD.plating(e=z)`、flan `mutestim(plateff=z)` **是同一模型**：
  pmf 逐点相同（≤1e-16），同一数据 m̂ 相对差 2e-7（求解器容差）；见
  `source_code_audit.md` §4 与 tests/test_mss_mle.py。

与 FALCOR 的分工：本模块只覆盖 MSS-MLE + p0 交叉核验；LC 中位数法、频率法见
《波动分析_突变率测定操作规程.md》§5，本软件不重复实现。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# MSS-MLE 的 m 适用域（规程 §0/§5.2；p0=exp(-m) 在 m>745 下溢）
M_MIN, M_MAX = 1e-4, 300.0
#: 菌株卡里可以由 rate 输出直接填的字段名（片段键名必须是其中之一）
KNOWN_CARD_FIELDS = ("M_base", "M_on_ref", "M_off", "mu_supply", "background_false")
# profile-likelihood 阈值：2*(LLmax-LL) = chi2_{1,0.95}=3.8415
CHI2_1_95_HALF = 3.841458820694124 / 2.0


def ld_pmf(m: float, r_max: int) -> np.ndarray:
    """Luria–Delbrück 分布 p_r（r=0..r_max），MSS 递推，逐点精确。

    注意：理想 LD 分布重尾（P(clone≥k)≈1/k），有限 r_max 的质量和 <1 是模型性质
    （缺失质量 ≈ m/(r_max+1)），不是数值误差——**不要**按有限和归一化。
    """
    if m <= 0:
        p = np.zeros(r_max + 1)
        p[0] = 1.0
        return p
    if m > M_MAX:
        raise ValueError(
            f"m={m:.1f} 超出 MSS-MLE 适用域（>{M_MAX:.0f}）：改实验设计压 m（规程 §3.2）"
        )
    p = np.empty(r_max + 1, dtype=float)
    p[0] = math.exp(-m)
    for r in range(1, r_max + 1):
        # 分母 = r-i+1, i=0..r-1  →  r+1, r, ..., 2
        p[r] = (m / r) * np.dot(p[:r], 1.0 / np.arange(r + 1, 1, -1, dtype=float))
    return p


def _make_pmf_fn(fraction_plated: float, r_max_obs: int):
    """返回 pmf_fn(m) -> p_obs[0..r_max_obs]。

    整管（z=1）：恒等截断。部分涂板（z<1）：二项稀释——走多项式系数法
    （`_poly_coeffs`），与逐管模式共用同一份精确实现。
    """
    if r_max_obs > _POLY_RMAX_CAP:
        raise ValueError(
            f"max(r)={r_max_obs} 超出适用域（>{_POLY_RMAX_CAP}）：疑似 jackpot / 计数异常——"
            "先确认不是污染或板上饱和；确属真值时可 winsor（文献建议 r>150 并为一类）后重算")
    z = float(fraction_plated)
    if not (0.0 < z <= 1.0):
        # 1.2 之类的值以前会被当成"≥1 → 整管"静默处理，等于悄悄换成另一个模型
        raise ValueError(
            f"涂板比例 fraction_plated={z:g} 不在 (0,1]：必须是 (0,1] 的小数"
            "（填了 50 这类百分数请改成 0.5）")
    if z >= 1.0:
        return lambda m: ld_pmf(m, r_max_obs)
    B = _poly_coeffs(r_max_obs, z)

    def pmf_fn(m: float) -> np.ndarray:
        return _pmf_from_coeffs(B, m, r_max_obs)

    def pmf_rows(m: float, rows) -> np.ndarray:
        """只算指定 r 的概率（对数似然用，省掉无关行）。"""
        return _pmf_from_coeffs(B, m, r_max_obs, rows=rows)

    pmf_fn.rows = pmf_rows
    return pmf_fn


def _expand_counts(counts) -> tuple[np.ndarray, np.ndarray]:
    """接受 {r: n} 字典 / 逐管 r 一维数组，返回 (r_values, n_cultures)。"""
    if isinstance(counts, dict):
        items = sorted((int(k), int(v)) for k, v in counts.items())
        rs = np.array([k for k, _ in items], dtype=int)
        ns = np.array([v for _, v in items], dtype=float)
    else:
        raw = np.asarray(counts)
        if raw.dtype.kind == "f" and not np.all(raw == np.floor(raw)):
            raise ValueError("逐管计数 r 必须是整数（发现非整数值）——检查数据录入")
        arr = raw.astype(int)
        if arr.ndim != 1 or arr.size == 0:
            raise ValueError("counts 应为逐管 r 的一维数组或 {r: 管数} 字典")
        if np.any(arr < 0):
            raise ValueError(
                f"计数 r 不能为负（发现 r={int(arr.min())}）——负值会被直方图静默丢掉，"
                "C 与似然都不再是你以为的那份数据；请检查数据录入")
        rs = np.arange(arr.max() + 1)
        ns = np.array([(arr == r).sum() for r in rs], dtype=float)
    _check_rs_ns(rs, ns)
    if ns.sum() < 1:
        raise ValueError("没有培养物数据")
    return rs, ns


def _check_rs_ns(rs: np.ndarray, ns: np.ndarray) -> None:
    """校验计数与管数：r 必须是非负整数（负号会把系数阵索引成末行，静默出错）。"""
    if rs.size and int(rs.min()) < 0:
        raise ValueError(
            f"计数 r 不能为负（发现 r={int(rs.min())}）——负数会被当成索引从末尾取值，"
            "静默给出错误的似然；请检查数据录入")
    if np.any(~np.isfinite(ns)) or np.any(ns < 0):
        raise ValueError("每类管数 n_r 必须是非负有限数")


def log_likelihood(m: float, rs: np.ndarray, ns: np.ndarray, pmf_fn=None) -> float:
    """对数似然 LL(m) = Σ n_r · ln(p_obs_r)。"""
    if pmf_fn is None:
        pmf_fn = _make_pmf_fn(1.0, int(rs.max()))
    rows_fn = getattr(pmf_fn, "rows", None)
    p = rows_fn(m, rs) if rows_fn is not None else pmf_fn(m)[rs]
    if np.any(p <= 0) or not np.all(np.isfinite(p)):
        return -math.inf
    return float(np.dot(ns, np.log(p)))


def _mle_from_rs_ns(rs, ns, pmf_fn) -> tuple[float, float]:
    """ML 估计 m：粗网格 + 黄金分割精修（对 ln m 一维搜索）。"""
    xs = np.linspace(math.log(M_MIN), math.log(M_MAX), 200)
    lls = np.array([log_likelihood(math.exp(x), rs, ns, pmf_fn) for x in xs])
    i = int(np.nanargmax(lls))
    lo, hi = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - inv_phi * (b - a), a + inv_phi * (b - a)
    fc, fd = log_likelihood(math.exp(c), rs, ns, pmf_fn), log_likelihood(math.exp(d), rs, ns, pmf_fn)
    for _ in range(80):
        if b - a < 1e-6:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = log_likelihood(math.exp(c), rs, ns, pmf_fn)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = log_likelihood(math.exp(d), rs, ns, pmf_fn)
    x_hat = (a + b) / 2.0
    return math.exp(x_hat), log_likelihood(math.exp(x_hat), rs, ns, pmf_fn)


def mle_m(counts) -> tuple[float, float]:
    rs, ns = _expand_counts(counts)
    pmf_fn = _make_pmf_fn(1.0, int(rs.max()))
    return _mle_from_rs_ns(rs, ns, pmf_fn)


def profile_ci(rs, ns, m_hat: float, ll_hat: float, pmf_fn) -> tuple[float, float]:
    """profile-likelihood 95% CI：在两侧解 LL = LL_hat − 1.921。

    LL 对 ln m 近似单峰；从 m̂ 向两侧游走找跨阈点，再二分收紧。
    游走到边界仍不跨阈（如全零管数据）→ 返回边界值（CI 开口）。
    直接吃已展开的 (rs, ns)——勿传原始 counts（会被误读为逐管观测）。
    """
    target = ll_hat - CHI2_1_95_HALF
    log_min, log_max = math.log(M_MIN), math.log(M_MAX)

    def solve_side(direction: int) -> float:
        x, step = math.log(m_hat), 0.2
        x_last_ge = x  # 最后一个 LL ≥ target 的点
        x_cross = None
        for _ in range(60):
            xn = x + direction * step
            if xn <= log_min or xn >= log_max:
                return math.exp(log_min if direction < 0 else log_max)  # 边界
            if log_likelihood(math.exp(xn), rs, ns, pmf_fn) < target:
                x_cross = xn
                break
            x_last_ge, x, step = xn, xn, step * 1.6
        if x_cross is None:  # 保险：游走用尽仍未跨阈
            return math.exp(log_min if direction < 0 else log_max)
        a, b = x_last_ge, x_cross  # LL(a) ≥ target > LL(b)
        for _ in range(80):
            mid = 0.5 * (a + b)
            if log_likelihood(math.exp(mid), rs, ns, pmf_fn) >= target:
                a = mid
            else:
                b = mid
            if abs(b - a) < 1e-6:
                break
        return math.exp(a)

    return solve_side(-1), solve_side(+1)


def p0_with_plating(m: float, zeta: float, j_max: int | None = None) -> float:
    """部分涂板下的零管概率：P0 = G(1−ζ)，G 为 L-D 分布的 PGF。

    ζ=1 时退化为 e^(−m)。j_max 只需让 (1−ζ)^j 足够小：ζ≥0.25 取 400，否则 800。
    """
    if zeta >= 1.0:
        return math.exp(-m)
    if j_max is None:
        j_max = 400 if zeta >= 0.25 else 800
    p = ld_pmf(m, j_max)
    q = 1.0 - zeta
    return float(np.dot(p, q ** np.arange(j_max + 1)))


def p0_crosscheck(rs: np.ndarray, ns: np.ndarray,
                  fraction_plated: float = 1.0) -> tuple[float, bool]:
    """p₀ 法：整管 m=−ln(p₀)；部分涂板解 G(1−ζ)=p₀（L-D 稀释模型）。仅 m∈[0.3,2.3] 有效。"""
    p0 = float(ns[rs == 0].sum()) / float(ns.sum())
    if p0 <= 0:
        return math.inf, False
    if fraction_plated >= 1.0:
        m = -math.log(p0)
    else:
        lo, hi = 1e-4, 50.0
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            if p0_with_plating(mid, fraction_plated) > p0:
                lo = mid
            else:
                hi = mid
        m = math.sqrt(lo * hi)
    return m, 0.3 <= m <= 2.3


@dataclass
class RateResult:
    m: float
    m_ci: tuple[float, float]
    mu: float                      # per CFU per passage（= m / mean Nt）
    mu_ci: tuple[float, float]
    mean_nt: float
    n_cultures: int
    p0: float
    m_p0: float                    # p₀ 法交叉核验
    m_p0_valid: bool
    mean_r: float
    median_r: float
    notes: list = field(default_factory=list)

    def as_yaml_snippet(self, label: str = "group", suppress_mu: bool = False,
                        key: str = None) -> str:
        """粘进 strain_card 的 YAML 片段。

        key：写进 YAML 的字段名。缺省用 label——但 label 常是人起的名字（如 "WT"），
        不是卡片字段（M_base/M_on_ref/M_off/mu_supply），直接粘会报"缺 M_base"。
        所以只要 key 不是已知字段名，就在片段里显式提示改名。

        suppress_mu=True（未提供 Nt 时）：**不能**输出 μ——否则用户会把 m/1 当成速率
        粘进卡片（差 8–9 个数量级）。改成注释模板，提示补 Nt 后自己算。
        """
        key = key or label
        if key not in KNOWN_CARD_FIELDS:
            caution_key = (f"# ⚠ 键名 \"{key}\" 不是菌株卡字段——粘进去会报\"缺 M_base\"。"
                           f"请改成 {' / '.join(KNOWN_CARD_FIELDS[:4])} 之一"
                           "（或下次用 `--field M_base` 直接生成）\n")
        else:
            caution_key = ""
        lo, hi = self.mu_ci
        p0_note = "" if self.m_p0_valid else "（无效区间，仅参考）"
        head = (f"# {label}: C={self.n_cultures}, mean_r={self.mean_r:.2f}, "
                f"p0={self.p0:.3f}, m={self.m:.3f} [CI {self.m_ci[0]:.2f}–{self.m_ci[1]:.2f}]"
                f"｜p0 法 m={self.m_p0:.2f}{p0_note}")
        # ★ 2026-10-29：μ 是**单位点**（如 rpoB/Rif）的每细胞每次分裂速率；
        #   窗口模型里的 M_base/M_on/M_off 是**全基因组每代**的乘客数（示例卡 0.004 = Drake 量级）。
        #   两者差 10³–10⁶ 倍。不写清楚，新手会把 μ 直接粘进 M_base，得到"零乘客"的假结论。
        caution = ("# ⚠ 这是**单位点**速率（per locus，如 Rif/rpoB），不是模型的 M：\n"
                   "#   M_base / M_on / M_off 要的是**全基因组每代**的乘客数（示例卡 0.004 = Drake 量级），\n"
                   "#   两者差 10³–10⁶ 倍，**不能直接粘**；只有比值（如 M_off/M_base）可以跨量纲比较。\n")
        if suppress_mu:
            return (
                f"{head}\n"
                f"{caution_key}"
                f"{caution}"
                f"# ⚠ 未提供 Nt：μ 未计算（m 本身有效）。补上各管 Nt 均值后填：\n"
                f"#   value = {self.m:.6g} / mean(Nt)   ，ci 同理除以 mean(Nt)\n"
                f"# {key}:\n"
                f"#   value: <m/mean_Nt>\n"
                f"#   ci: [<lo/mean_Nt>, <hi/mean_Nt>]\n"
                f"#   unit: per_locus_per_passage\n"
            )
        return (
            f"{head}\n"
            f"{caution_key}"
            f"{caution}"
            f"{key}:\n"
            f"  value: {self.mu:.3e}\n"
            f"  ci: [{lo:.3e}, {hi:.3e}]\n"
            f"  unit: per_locus_per_passage\n"
        )


def analyze(counts, nts, fraction_plated: float = 1.0, label: str = "group") -> RateResult:
    """主入口。

    counts : {r: 管数} 或逐管 r 一维数组（部分涂板时给**观测值**，fraction_plated<1）
    nts    : 逐管 Nt（CFU/管）；或标量 = 各管均值
    fraction_plated : 每管涂板比例（整管=1.0），二项稀释模型
    """
    rs, ns = _expand_counts(counts)
    pmf_fn = _make_pmf_fn(fraction_plated, int(rs.max()))
    nts_arr = np.atleast_1d(np.asarray(nts, dtype=float))
    if nts_arr.size == 1:
        mean_nt = float(nts_arr[0])
    else:
        if nts_arr.size != int(ns.sum()):
            raise ValueError(
                f"nts 有 {nts_arr.size} 个但 counts 共 {int(ns.sum())} 管；"
                "逐管配对时 counts 请传逐管一维数组"
            )
        mean_nt = float(nts_arr.mean())

    m_hat, ll_hat = _mle_from_rs_ns(rs, ns, pmf_fn)
    lo, hi = profile_ci(rs, ns, m_hat, ll_hat, pmf_fn)
    m_p0, m_p0_valid = p0_crosscheck(rs, ns, fraction_plated)
    notes = []
    if fraction_plated < 1.0:
        notes.append(f"部分涂板 ζ={fraction_plated:g}：二项稀释似然；Nt 必须是整管总 CFU"
                     "（不是涂板那部分的），且不要手工把 r 除以 ζ")
    if ns.sum() < 20:
        notes.append(f"C={int(ns.sum())}<20：管数偏少——1σ=20% 精度所需 C 随 m 变"
                     "（m≈1 约 36 管、m≈2 约 23、m≈0.6 约 54）；用 power 子命令按你的 m 算")
    if m_hat > 15:
        notes.append(f"m={m_hat:.1f}>15：超出设计目标区间，按保守处理")
    if m_hat < 0.3:
        notes.append(f"m={m_hat:.2f}<0.3：R&F 不建议（除非大量培养物）；精度很差，"
                     "建议提高 Nt 或增加 C")
    if m_p0_valid and max(m_p0, m_hat) / min(m_p0, m_hat) > 2.0:
        notes.append("p₀ 法与 MSS-MLE 差 >2×：查数据（规程 §6）")
    if lo <= M_MIN * 1.001 or hi >= M_MAX * 0.999:
        notes.append("CI 触及适用域边界：m 过低（数据不足）或过高（设计需压 m）")
    mu = m_hat / mean_nt
    return RateResult(
        m=m_hat, m_ci=(lo, hi), mu=mu, mu_ci=(lo / mean_nt, hi / mean_nt),
        mean_nt=mean_nt, n_cultures=int(ns.sum()),
        p0=float(ns[rs == 0].sum()) / float(ns.sum()),
        m_p0=m_p0, m_p0_valid=m_p0_valid,
        mean_r=float((rs * ns).sum() / ns.sum()),
        median_r=float(rs[np.searchsorted(np.cumsum(ns), ns.sum() / 2.0)]),
        notes=notes,
    )


def simulate(m_expected: float, n_cultures: int, seed: int = 0) -> np.ndarray:
    """复合泊松模拟 LD 数据（**仅性质测试用**）：事件数~Poisson(m)，克隆大小 P(k)=1/(k(k+1))。

    不做克隆封顶——理想骨架均值无穷，偶尔会吐出一个 1e9 级的 jackpot，足以把下游
    系数法顶到适用域外。要跑估计器的模拟回收请用 `compare.simulate_assay`（有封顶）。
    """
    rng = np.random.default_rng(seed)
    out = np.zeros(n_cultures, dtype=int)
    for i in range(n_cultures):
        n_events = rng.poisson(m_expected)
        if n_events:
            u = rng.uniform(1e-12, 1.0, size=n_events)
            k = np.ceil(1.0 / u - 1.0).astype(int)  # 逆变换采样，F(k)=1-1/(k+1)
            out[i] = int(k.sum())
    return out


# ================= 逐管 Nt 模式（Ycart 2014 无偏估计） =================
#
# 经典 MSS-MLE 假设各管终末细胞数 Nt 相同；Nt 波动时速率被系统性低估（Ycart & Veziris
# 2014, PLoS ONE）。逐管模式把参数改为 μ（每细胞每分裂突变概率），每管事件数
# K_i ~ Poisson(μ·N_i)，观测 r_i ~ MSS-pmf(m_i=μ·N_i)（部分涂板再叠二项稀释）。
#
# 数值实现：MSS pmf 可写成 m 的多项式 p_r(m) = e^{−m}·Σ_{k≤r} a_{r,k}·m^k，
# 系数 a_{r,k} 与 m 无关（a_{0,0}=1；a_{r,k}=(1/r)Σ_{i<r} a_{i,k−1}/(r−i+1)）。
# 预计算三角阵后任意 (r, m) 用 Horner 求值 O(r)。
# ⚠️ 部分涂板（z<1）时 p_obs(r) 是 m 的 **j_max** 次多项式（二项混合把次数抬上去），
# 必须用整行系数求值；系数全为正、值 ≤ e^m，Horner 天然无溢出。

_POLY_CACHE = {}
_POLY_RMAX_CAP = 5000       # 观测 r 上限（系数阵 O(r²) 内存）
_BINOM_JMAX_CAP = 3000      # 二项混合的 j（未涂板真值）上限（j_max³/3 建表成本 ~2.5 s）
_POLY_CACHE_BUDGET = 32_000_000   # 缓存预算（元素数，≈256 MB）；超了先淘汰最大的那张
_POLY_CACHE_ELEMS = 0


def _partial_jmax(r_max_obs: int, z: float) -> int:
    """部分涂板二项混合需要的 j 上限。

    给定观测 r，C(j,r)·z^r·(1−z)^(j−r) 的质量峰在 j≈r/z，标准差
    √(r(1−z))/z——取 +10 sd +10 的余量。不要用 30·r_max 这类系数：z 小时
    矩阵会大到爆内存，而且先 exp(logC) 再乘会溢出成 inf/nan。
    """
    if z >= 1.0:
        return int(r_max_obs)
    j = r_max_obs / z + 10.0 * math.sqrt(max(r_max_obs, 1) * (1.0 - z)) / z + 10.0
    if j > _BINOM_JMAX_CAP:
        raise ValueError(
            f"部分涂板 + 大 jackpot 超出精确二项混合的适用域：max(r)={r_max_obs}, "
            f"fraction_plated={z:g} 需要 j≈{j:.0f}（> {_BINOM_JMAX_CAP}）。"
            "建议整管涂板，或先确认该管不是污染/计数错误；确属真 jackpot 时"
            "请在报告中注明该管并按 winsor 处理后重算（不要用截断后的混合硬算）")
    return int(math.ceil(j))


def _pmf_from_coeffs(B: np.ndarray, m: float, r_max_obs: int, rows=None) -> np.ndarray:
    """p_obs(r|m) = e^{−m}·Σ_k B[r,k]·m^k（Horner，对 k 向量化）。

    B 的每一行是**完整**的 j_max+1 个系数：部分涂板下 p_obs(r) 关于 m 的次数是
    j_max（不是 r）。截断到 :r+1 会把 p_0 算小三倍（回归测试见 tests）。
    """
    if rows is not None:
        rows = np.asarray(rows, dtype=int)
        out = np.zeros(rows.size)
        for k in range(B.shape[1] - 1, -1, -1):
            out *= m
            out += B[rows, k]
        return out * math.exp(-m)
    out = np.zeros(r_max_obs + 1)
    for k in range(B.shape[1] - 1, -1, -1):
        out *= m
        out += B[: r_max_obs + 1, k]
    return out * math.exp(-m)


def _poly_coeffs(r_max_obs: int, z: float = 1.0):
    """返回 B[r,k]：p_obs(r|m) = e^{−m}·Σ_k B[r,k]·m^k（z<1 时含二项稀释）。

    成本 O(j_max³/3)（j_max = r_max 或 r_max/z + 余量）——j_max 过大（数据含
    极端 jackpot）会又慢又占内存；观测 r 上限 5000（再大请检查数据/winsor，
    或用 mean 模式；power.py 的模拟器已对克隆封顶于 N_i，正常不会触发）。
    """
    if r_max_obs > _POLY_RMAX_CAP:
        raise ValueError(
            f"max(r)={r_max_obs} 超出多项式系数法适用域（>{_POLY_RMAX_CAP}）：数据含极端 jackpot——"
            "请检查数据或用 mean 模式/winsor 处理")
    z = float(z)
    if not (0.0 < z <= 1.0):
        raise ValueError(f"fraction_plated={z} 不在 (0,1]：涂板比例必须为正且 ≤1")
    key = (int(r_max_obs), round(z, 9))
    if key in _POLY_CACHE:
        return _POLY_CACHE[key]
    j_max = _partial_jmax(int(r_max_obs), z)
    A = np.zeros((j_max + 1, j_max + 1))
    A[0, 0] = 1.0
    for r in range(1, j_max + 1):
        w = 1.0 / np.arange(r + 1, 1, -1)              # 1/(r−i+1), i=0..r−1
        A[r, 1:r + 1] = (w @ A[0:r, 0:r]) / r
    if z >= 1.0:
        B = A
    else:
        # 二项转移阵在**对数域**构造：先 exp(logC) 会溢出（C(6000,200)≈1e380）
        lgam = np.vectorize(math.lgamma)
        log_fact = lgam(np.arange(1, j_max + 2))       # log_fact[n] = ln(n!)
        log_z, log_1z = math.log(z), math.log1p(-z)
        T = np.zeros((r_max_obs + 1, j_max + 1))
        for r in range(r_max_obs + 1):
            jj = np.arange(r, j_max + 1)
            logC = log_fact[jj] - log_fact[r] - log_fact[jj - r]
            logT = logC + r * log_z + (jj - r) * log_1z
            T[r, r:] = np.exp(logT)                    # logT ≤ 0 → 不会溢出
        B = T @ A
    _cache_put(key, B)
    return B


def _cache_put(key, B: np.ndarray) -> None:
    """带预算的缓存：大 jackpot 的 r_max 差异很大，无界缓存会吃光内存；
    但一律不缓存又会让同一个大 r_max 反复建表（O(r³)，秒级）。折中：超预算先淘汰最大的。"""
    global _POLY_CACHE_ELEMS
    if B.size > _POLY_CACHE_BUDGET:          # 单张就超预算：不缓存，也別把别的挤掉
        return
    while _POLY_CACHE and _POLY_CACHE_ELEMS + B.size > _POLY_CACHE_BUDGET:
        big = max(_POLY_CACHE, key=lambda k: _POLY_CACHE[k].size)
        _POLY_CACHE_ELEMS -= _POLY_CACHE.pop(big).size
    _POLY_CACHE[key] = B
    _POLY_CACHE_ELEMS += B.size


def _ll_per_tube_factory(rs: np.ndarray, nt_vec: np.ndarray, z: float = 1.0):
    """逐管模式对数似然工厂：ll(ln μ) = Σ_i log p_obs(r_i; μ·N_i)。

    性能：按唯一 r 值分组的向量化 Horner（np.polyval 批量），比逐管循环快 ~10×；
    附 warm 属性（模拟研究用：真值已知时在真值 ±ln100 内黄金分割，跳过全局网格）。
    部分涂板（z<1）用**整行**系数（次数 j_max，不是 r）。
    """
    B = _poly_coeffs(int(rs.max()), z)
    full_row = float(z) < 1.0
    groups = {}
    for r, N in zip(rs, nt_vec):
        groups.setdefault(int(r), []).append(float(N))
    groups = {r: np.array(Ns) for r, Ns in groups.items()}
    coeffs = {r: np.ascontiguousarray(B[r, ::-1] if full_row else B[r, :r + 1][::-1])
              for r in groups}

    def _total(mu: float) -> float:
        total = 0.0
        with np.errstate(over="ignore"):    # 极端 μ 下 Horner 溢出 → s=inf → −inf，行为正确
            for r, Ns in groups.items():
                ms = mu * Ns
                s = np.polyval(coeffs[r], ms)      # 向量化 Horner
                if np.any(s <= 0.0) or not np.all(np.isfinite(s)):
                    return -math.inf
                total += float(np.sum(-ms + np.log(s)))
        return total

    def ll(ln_mu: float) -> float:
        return _total(math.exp(ln_mu))

    ll.total = _total
    ll.warm = lambda ln_mu: _golden(ll, ln_mu - math.log(100.0),
                                    ln_mu + math.log(100.0))
    return ll


def _maximize_1d(ll_fn, ln_lo: float, ln_hi: float, n_grid: int = 200):
    """一维 ML：粗网格 + 黄金分割精修。n_grid=0 → 纯 golden（要求括号内单峰，
    且解不贴边——贴边自动退回全网格重搜）。"""
    if n_grid > 0:
        xs = np.linspace(ln_lo, ln_hi, n_grid)
        lls = np.array([ll_fn(x) for x in xs])
        i = int(np.nanargmax(lls))
        lo, hi = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
        x, ll = _golden(ll_fn, lo, hi)
        return x, ll
    x, ll = _golden(ll_fn, ln_lo, ln_hi)
    if ln_lo + 1e-6 < x < ln_hi - 1e-6:
        return x, ll
    xs = np.linspace(ln_lo, ln_hi, 200)
    lls = np.array([ll_fn(x) for x in xs])
    i = int(np.nanargmax(lls))
    lo, hi = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    return _golden(ll_fn, lo, hi)


def _golden(ll_fn, a: float, b: float, tol: float = 1e-8):
    """纯黄金分割（假设 [a,b] 内单峰）——warm 路径用，省去全局网格。"""
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - inv_phi * (b - a), a + inv_phi * (b - a)
    fc, fd = ll_fn(c), ll_fn(d)
    for _ in range(80):
        if b - a < tol:
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = ll_fn(c)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = ll_fn(d)
    x = (a + b) / 2.0
    return x, ll_fn(x)


def _profile_ci_1d(ll_fn, x_hat: float, ll_hat: float, ln_lo: float, ln_hi: float):
    """profile-likelihood CI（ln 参数）：两侧游走找跨阈再二分。阈值 chi2_{1,0.95}/2。"""
    target = ll_hat - CHI2_1_95_HALF

    def solve_side(direction: int) -> float:
        x, step = x_hat, 0.2
        x_last_ge = x
        x_cross = None
        for _ in range(60):
            xn = x + direction * step
            if xn <= ln_lo or xn >= ln_hi:
                return ln_lo if direction < 0 else ln_hi
            if ll_fn(xn) < target:
                x_cross = xn
                break
            x_last_ge, x, step = xn, xn, step * 1.6
        if x_cross is None:
            return ln_lo if direction < 0 else ln_hi
        a, b = x_last_ge, x_cross
        for _ in range(80):
            mid = 0.5 * (a + b)
            if ll_fn(mid) >= target:
                a = mid
            else:
                b = mid
            if abs(b - a) < 1e-7:
                break
        return a

    return solve_side(-1), solve_side(+1)


def _mu_domain(nt_arr: np.ndarray) -> tuple[float, float]:
    """ln μ 的搜索域：保证**每一管**的 m_i = μ·N_i 都落在 [M_MIN, M_MAX]。

    下界由**最小** N 决定（μ·N_min ≥ M_MIN），上界由**最大** N 决定（μ·N_max ≤ M_MAX）。
    2026-10-29 修：旧实现把 max/min 对调，搜索域被放宽 (N_max/N_min)² 倍，
    200 点粗网格被摊薄（最优解是内点，所以数值影响小，但属纯逻辑错误）。
    """
    return (math.log(M_MIN / nt_arr.min()), math.log(M_MAX / nt_arr.max()))


def analyze_per_tube(counts, nts, fraction_plated: float = 1.0,
                     label: str = "group") -> RateResult:
    """逐管 Nt 模式主入口（Ycart 2014 无偏估计）。

    counts : 逐管 r 的一维数组（须与 nts 逐管配对）
    nts    : 逐管终末细胞数 N_i（实测）
    fraction_plated : 涂板比例 z（二项稀释）
    返回 RateResult：mu = μ̂（每细胞每分裂），m 字段为 μ̂·mean(Nt) 的等效 m。
    """
    raw = np.asarray(counts)
    if raw.dtype.kind == "f" and not np.all(raw == np.floor(raw)):
        raise ValueError("逐管计数 r 必须是整数（发现非整数值）——检查数据录入")
    arr = raw.astype(int)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("per_tube 模式要求逐管 r 的一维数组")
    if np.any(arr < 0):
        raise ValueError(
            f"计数 r 不能为负（发现 r={int(arr.min())}）——负数会被当成系数阵的行索引，"
            "静默给出错误的似然；请检查数据录入")
    nt_arr = np.atleast_1d(np.asarray(nts, dtype=float))
    if nt_arr.size == 1 and arr.size > 1:
        raise ValueError(
            "per_tube 模式需要逐管 Nt（与 counts 一一配对）；常数 Nt 请用 mean 模式（analyze）")
    if nt_arr.size != arr.size:
        raise ValueError(f"nts 有 {nt_arr.size} 个但 counts 共 {arr.size} 管——逐管模式须一一配对")
    if arr.size < 2:
        raise ValueError("至少需要 2 管才能估计速率")
    rs = arr
    mean_nt = float(nt_arr.mean())

    ll = _ll_per_tube_factory(rs, nt_arr, fraction_plated)
    ln_lo, ln_hi = _mu_domain(nt_arr)
    x_hat, ll_hat = _maximize_1d(ll, ln_lo, ln_hi)
    lo_ln, hi_ln = _profile_ci_1d(ll, x_hat, ll_hat, ln_lo, ln_hi)
    mu_hat, mu_lo, mu_hi = math.exp(x_hat), math.exp(lo_ln), math.exp(hi_ln)

    m_p0, m_p0_valid = p0_crosscheck(rs, np.ones_like(rs, dtype=float), fraction_plated)
    notes = [f"nt_mode=per_tube（Ycart 2014 逐管校正，CV_Nt={nt_arr.std()/mean_nt:.2f}）"]
    if fraction_plated < 1.0:
        notes.append(f"部分涂板 ζ={fraction_plated:g}：二项稀释似然；Nt 必须是整管总 CFU"
                     "（不是涂板那部分的），且不要手工把 r 除以 ζ")
    if arr.size < 20:
        notes.append(f"C={arr.size}<20：MSS-MLE 精度不足（规程 §3.2）")
    if mu_hat * mean_nt > 15:
        notes.append("等效 m>15：超出设计目标区间")
    if mu_hat * mean_nt < 0.3:
        notes.append(f"等效 m={mu_hat * mean_nt:.2f}<0.3：R&F 不建议（除非大量培养物）；"
                     "精度很差，建议提高 Nt 或增加 C")

    return RateResult(
        m=mu_hat * mean_nt,
        m_ci=(mu_lo * mean_nt, mu_hi * mean_nt),
        mu=mu_hat, mu_ci=(mu_lo, mu_hi),
        mean_nt=mean_nt, n_cultures=int(arr.size),
        p0=float((rs == 0).sum()) / arr.size,
        m_p0=m_p0, m_p0_valid=m_p0_valid,
        mean_r=float(rs.mean()),
        median_r=float(np.median(rs)),
        notes=notes,
    )
