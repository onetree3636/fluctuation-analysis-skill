# 多靶点路径的代码地图——模型、算法、数据处理

> 给**维护者与审稿人**看：这条路径用了哪些模型、每个公式落在哪个函数里、
> 数据在进入模型前被怎么处理、以及**哪些东西明确不建模**。
> 标注体系：[文献] 有出处 / [推导] 算术或从文献推出 / [约定] 本项目口径。
> 通用性政策：本文档不含任何项目实测数值（见 `generic-policy.md`）。

---

## 0. 一句话

多靶点路径 = **一个账本解释器 + 六个成熟模型的组合**。核心分工是一句话：
**供给只认实测 `f_step`；μ 只进乘客（`M = μ_perbp × G`）。**

```
records.csv ─┐
markers.csv ─┼→ cmd_multitarget (cli.py:1157)
             │      ├─ mu_perbp()      ← 过原点回归
             │      └─ run_ledger()    ← ramp.py:312  账本解释器
             │             ├─ wilson_ci / newcombe_diff   （逐档 f_step）
             │             ├─ fit_pap()                   （尾部曲线）
             │             └─ plan_ramp() → render_ledger()
             └─ 输出五节：各档 f_step｜PAP 现状｜下一窗设计｜剩余梯级｜怎么读
```

---

## 1. 代码地图

| 文件 | 行数 | 角色 | 关键入口 |
|---|---|---|---|
| `windowcalc/ramp.py` | 600+ | **主引擎**：账本 → 现状 → 规划 | `run_ledger` · `fit_pap` · `plan_ramp` · `washout_tail` (139) · `render_ledger` |
| `windowcalc/fstep.py` | 101 | 比例估计的数学底座 | `wilson_ci` (23) · `newcombe_diff` (43) |
| `windowcalc/mtblocks.py` | 162 | 标定积木（B1 / B5 / B6） | `mu_perbp` (66) · `minimal_n0` (95) · `window_table_f_step` (118) |
| `windowcalc/induction.py` | 301 | 剂量轴机制模型 → `M_on(c)` | `fit_induction` (197) · `psi` (49) |
| `windowcalc/fitting.py` | 279 | 剂量-响应 Hill、洗脱衰变 | `fit_tau_off` (204) · `hill_fc` (28) |
| `windowcalc/model.py` | 879 | 窗口引擎（j=1 的 `optimize` 用） | `phase_terms` (261) · `monte_carlo` (754) · `protected_hits` (495) |
| `windowcalc/optimize.py` | 322 | 三旋钮网格 + MC 选优 | `optimize` (117) |
| `windowcalc/cli.py` | 1456 | 参数校验、CSV 读入、拼装 | `cmd_multitarget` (1157) |

---

## 2. 模型

### 2.1 比例估计（f_step）

```
Wilson 分数区间（fstep.py:23）
  den = 1 + z²/n
  center = (p + z²/2n)/den
  half   = z/den · √(p(1−p)/n + z²/4n²)        z = 1.959963985
  k = 0 特例 → 上界 3/n（rule of three）——"零克隆"是"≤ 3/n"，不是"没有"   [文献]
Newcombe 差值（fstep.py:43，1998 法 10）  [文献]
  d = p_post − p_pre；区间由两个 Wilson 区间混合构造
```

**口径** [约定]：`f_step = f_post − f_pre`；`f_pre ≥ 0.5 × f_post` → 背景淹没红灯，本档作废。
没铺 pre 板 → 不扣背景，且报告必须注明。

### 2.2 耐受限尾部（PAP 曲线）

```
log10 f(c) = a + s·c ,  s < 0            δ = −1/s（十倍频浓度宽）   [推导]
OLS：s = Sxy/Sxx；R²；σ² = SSres/(n−2)；se(s) = √(σ²/Sxx)
δ 的 CI：s 的 CI 经 −1/s 单调变换（端点交换）
预测：se(c) = √(σ²[1 + (c−x̄)²/Sxx])
```

三点数据处理规则（都在 `fit_pap` 里）：
- **同浓度多板池化**：克隆数与 CFU 分别求和（不是对比例取平均）；
- **零克隆行删失**：不进拟合，只作 `3/n` 上界；拟合值高于上界时报警；
- **饱和行剔除**：`f ≥ 0.5` 视为饱和（`run_ledger` 里判定），对尾部无信息。

### 2.3 阶梯规划

```
可行最大步长（ramp.py:163）  f(c_start + Δ) = k·margin / capacity
  ⇒ Δ = (log10(threshold) − intercept)/slope − c_start      ← 解析解，无搜索
档数 n = ⌈span/Δ⌉；步长 = span/n（均匀分档）
每档需 CFU = k / f；富余 = f·capacity/k
```

**平稳补充假设** [约定，`ramp.py:140`]：等步长下每档 `f_step` 都取曲线在
`c_start + Δ` 处的值（窗口把尾部重新填满）。该假设**每档都要被实测审计**：
回填账本重跑，拟合残差出 CI → 曲线重拟合或拆档。
`margin`（默认 3）是冗余系数：预测克隆数 ≥ `margin × k` 才允许开档。

### 2.4 乘客账（三段可加）

```
M_on = μ_perbp × G                       ← [B1] 过原点回归 μ̂ = Σ(T_i μ_i)/Σ(T_i²)   [文献]
每档 = n_win·M_on + n_plate·M_off         ← 窗口段 + 板内扩张段
总   = 已完成档基数 + 剩余梯级 Σ + n_post·M_off（只计一次）
```

`μ_perbp` 用**过原点回归**而不是平均：`T_eff` 因标记而异，μ 不可直接平均。
同时输出 `spread = max/min`，>2 说明"每碱基率均匀"这条假设在本菌里不成立。

**档间代数怎么来**（2026-11-15 接入 `--propagation`）：`inter_rung_generations()`
（ramp.py）按拓扑算——`manual` 用用户给的 `--n-plate`；`plate` 用
`log₂(cells_per_colony)`；`liquid` 用 `log₂(N₀/(f·Nt))`，N₀ 缺省取 `Nt − k/f`。
已完成档按**各自周期的 f** 算（不是一刀切），来源写进报告的"档间扩增"一行。

**第四段：洗脱尾**（2026-11-15 接入）——`washout_tail()`（ramp.py:139）：

```
tail = (M_on − M_off)/(e^(td/τ_off) − 1)     每次洗脱计一次（档 + 终点后）
```

由 `--tau-off-h` + `--td` 驱动；**只给一个直接报错**，都不给则报告明写"未计入"。
连续极限退化为 `(M_on−M_off)·τ/td`，与 `fit_tauoff` 的尾部积分同口径（有测试锁住）。

### 2.5 N₀ 的解析解（降荷的最强旋钮）

```
rate·(Nt − N0) ≥ k   ⇒   N0 ≤ Nt − k/rate ,   n_gen = log2(Nt/N0)   [推导]
```

一阶供给是线性的 ⇒ 取等号即"刚好够 k、代数最少、乘客最少"，**不需要数值搜索**。
账本第三节的 "N₀ ≤ Nt_sat − k/f" 是同一式子（`ramp.py:495`）。

### 2.6 保护区命中（泊松近似，含多拷贝）

```
f = kb·copies / (G + Σ质粒 kb·copies)
E = 乘客 × f
P(≥1) = 1 − exp(−E)                      ← 至少一个拷贝被打中
P(全部) = (1 − exp(−E/copies))^copies     ← 全部拷贝被打中才可能功能丧失
失守 = 1 − exp(−E[P]_总 × L_prot/G)        ← 总乘客单次代入，绝不按档相加
```

### 2.7 剂量轴机制模型（induction.py）

```
E(c) = Emax·cⁿ/(EC50ⁿ + cⁿ)                     Hill（表达量）
p(c) = E/(E+K) = q·h(c)                          Langmuir 竞争占位
M(c) = M_off + ΔM·φ(p)，  φ ∈ {p, 2p−p², p²}     三种复制叉装配模式
```

对 `A = ΔM·q` **线性** ⇒ 给定 (C50, n, q) 时 A 有闭式最小二乘解；
外层 40×20×9 网格 + 3 轮坐标下降；**AIC** 在三种混合规则里选（two_core 多一个参数）。

### 2.8 洗脱衰变（fitting.py）

```
μ(t) = μ_off + (μ_on − μ_off)·e^(−t/τ)
τ：一维网格（对数，400 点），μ_on 内层闭式解（profile 最小二乘）
CI：残差噪声的参数自助法（B 次），取 2.5/97.5 分位
尾部负荷 = (μ_on − μ_off)·τ/td     ← ∫₀^∞(μ−μ_off)dt / td
```

### 2.9 配套：窗口引擎（j=1 的 `optimize`，共享 M 口径）

```
每相位：ā（平均活性）；n（每世系代数）；M_eff = leak + (M_on − leak)·ā
乘客 = Σ n_i·M_eff,i；  供给 events = Σ μ_supply·ā_i·ΔN_i
成功 P = P(Poisson(events) ≥ k)      [文献]（jackpot 不算独立供体）
```

MC：参数 CI → 对数正态 σ（由 CI 反推；无 CI 时按 2× 假设并标出）；
**克隆间离散** `exp(N(0, σ_ln))` 只乘乘客，不乘 events。

---

## 3. 算法（数值方法）

| 方法 | 用在哪 | 规模 |
|---|---|---|
| 闭式解 | Wilson、Newcombe、μ_perbp、N₀、步长 Δ、给定 φ 的 A | O(1) |
| OLS + 标准误 | PAP 尾部 | O(n) |
| 网格搜索 | 剂量轴 (C50, n, q) | 40×20×9 |
| 坐标下降 | 剂量轴精修 | 3 轮 × 75 点 |
| 一维 profile | τ_off（τ 网格 × 内层闭式 μ_on） | 400 点 |
| 参数自助法 | τ_off 的 CI | B=200 |
| 蒙特卡洛 | 窗口引擎参数不确定性 | 4000 样本 |
| 泊松尾（递推） | P(≥k) | O(k) |
| 泊松近似 | 保护区命中 | O(1) |
| 闭式离散和 | 洗脱尾 `1/(e^(td/τ)−1)` | O(1) |

---

## 4. 数据处理（进入模型之前）

| 处理 | 位置 | 为什么 |
|---|---|---|
| 编码探测 + Excel 直读 | `cli.py:_read_text` | 中文 Windows 的 CSV 常是 GBK；不让用户"另存为" |
| 注释与表头跳过 | `cli.py:_read_xy_csv` / `_read_markers_csv` | 模板带 `#` 注释 |
| 多板池化 | `ramp.py:_pool` | 克隆数、CFU 分别求和 = 按 CFU 加权，不是平均比例 |
| 配对扣除 | `ramp.py:369` | 背景（假抗性 + 残留耐受）必须减掉 |
| 饱和切分 | `ramp.py:386` | 对 δ 无信息，还会压平曲线 |
| 删失 | `ramp.py` `fit_pap` 的 `censored` | 零克隆保留为"上界"，不丢信息 |
| 饱和阈值 | `run_ledger(sat_f=…)` | 显式参数（默认 0.5）；±20% 带内的点单独标出 |
| 周期排序 | `ramp.py:335` | 按数值排序，不是字典序（否则 "10" < "2"） |
| 缺参拒绝 | `cli.py:1164` | 不用默认值顶替（通用性政策） |

---

## 5. 明确不建模（边界，写在代码注释里）

1. **不用 μ×T_set 外推供给**（2026-11-13 删除整条几何路线）；
2. **没有 F^j**：j 步不是 F 的幂；
3. **无 C 律（谷型）**（2026-11-12 移除）；
4. **无谱系/上位效应**：每档独立、乘客可加，靠"平稳补充假设 + 残差审计"兜底；
5. **选择富集不建模**（它发生在涂板/传代之后）；
6. **板内突变混入克隆计数**：只在 `protocol-fstep.md` 里提示（双时间点读），模型里没有对应项；
7. **洗脱尾只在给全 τ_off+td 时才计**——不给就明确不计，不猜默认值。

---

## 6. 维护者须知：同一批公式有**两套**实现

| 公式 | 积木版（文档里的"主路径"） | 实跑版（CLI 真正调用） |
|---|---|---|
| N₀ = Nt − k/f | `mtblocks.minimal_n0` (95) | 账本第三节设计行 `ramp.py:495` |
| 乘客 = n·M_on | `mtblocks.window_table_f_step` (118) | `ramp.run_ledger` 的 `base`（409）+ `plan_ramp`（206） |
| 失守 = 1 − exp(−E·L/G) | 同上（`row["p_bgc"]`） | `plan_ramp:210`（按梯级累计） |

**现状** [2026-11-15 核查]：CLI 只用了 `mu_perbp`；`minimal_n0` 与
`window_table_f_step` 仅被测试调用——它们是**文档积木**，不是实跑路径。
两套公式现在数值一致，但**不是同一段代码**。

⇒ **改公式时两处都要改，并跑 `tests/test_ledger_formula_pin.py`。**
该测试用变异测试验证过：只改 mtblocks 或只改 ramp，都会立刻红
（mtblocks 侧 4 项失败、ramp 侧 2 项失败）。

```bash
python -m pytest tests/test_ledger_formula_pin.py -q     # 对拍
python -m pytest -q                                      # 全量
```

> 位置说明：这两个测试在**开发仓库**的 `tests/` 下；技能分发副本不含 `tests/`
> （见 `sync_skill.sh` 的"不同步"清单）。所以"改公式跑对拍"是**维护者在本仓库**的动作，
> 技能使用者不需要（也无法）跑它。

---

## 7. 数据流（一次完整调用）

```
rate  ──→ μ_i（逐标记，MSS-MLE + CI）
          └─→ markers.csv ──→ mu_perbp ──→ μ_perbp ──→ M_on = μ_perbp × G
records.csv（cycle, c, when, cfu, clones）
   └─→ run_ledger ──→ 逐档 f_step（Wilson/Newcombe）
                    ──→ 饱和点定 c_now
                    ──→ 尾部点 ──→ fit_pap ──→ δ, f_ref
                    ──→ plan_ramp ──→ 梯级、需 CFU、乘客、失守
                    ──→ render_ledger ──→ 五节报告
下一档实测 ──→ 回填 records.csv ──→ 重跑（唯一状态文件，逐行累积、永不改写）
```
