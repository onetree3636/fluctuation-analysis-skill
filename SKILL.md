---
name: fluctuation-analysis
description: 把波动实验的菌落计数变成带 95%CI 的突变率（MSS-MLE，支持逐管 Nt 与部分涂板），并做两批比较、样本量/功效、mutator 表达窗口、多靶点窗口（f_step 实测供给）、剂量响应与 τ_off 拟合。用户给出 r,Nt 计数、问“这两批有没有差 / 要铺多少管 / 窗口怎么定 / 多靶点药怎么算 / 剂量响应怎么拟合”，或要设计波动实验、问波动分析原理时使用。
license: MIT
compatibility: Python 3.9+；Linux/macOS/Windows；首次运行自动建 .venv 装 numpy/pyyaml/openpyxl，不需要 R。
allowed-tools: Bash Read Write Edit Glob Grep
---

# 波动分析与 mutator 窗口

入口只有一条：`python3 run.py <子命令>`（Windows：`run.cmd <子命令>`）。
第一次运行会自动准备环境；**不要让用户 `pip install`，不要运行 `install.sh`**。

## 设计目标（先读这个，再读铁律）

> **把最便宜的实验数据（菌落计数），变成最贵的工程决策（开 mutator 做耐受进化）的依据——
> 每个数字带误差棒、带来源、可被推翻；缺数据时拒绝给数。**

技能回答三类问题：

1. **测量**——"我的突变率是多少？可信吗？两组一样吗？要铺多少管？"（`rate` / `compare` / `power`）
2. **决策**——"mutator 窗口开多久/多强，能拿到 k 个耐受克隆，乘客压到最低？"（`optimize`；多靶点走 `multitarget --records`——账本解释器）
3. **机制**——"剂量怎么影响突变？去诱导后多久关掉？"（`fit-fc` / `fit-tauoff` / `design-induction` / `fit-induction`）

五条原则（下面的铁律是它们的可执行形式）：

1. **能实测的一律实测**；必须外推的，链条要短、要带误差棒（μ 只进乘客预算，供给一律实测）；
2. **每个数字带来源标签**（measured / literature / guess；报告 assumed 全标 ★）；
3. **缺数据时拒绝给数**（缺参五要素、U 类不给数、S≫1 不给绝对时间）；
4. **红灯必须转达**（泄漏/背景/边界解/spread>2——它们是结论的一部分，不是附录）；
5. **文档数字与代码同源**（selftest 现场算给你看）。

**通用性是硬约束（2026-11-15 用户明令，细则见 `docs/generic-policy.md`）**

- **药物名 / 诱导剂名只允许来自被引用的文献**；本项目用的药、诱导剂、靶点基因、菌株代号**一律不得出现**。
- **默认值只允许来自文献确证的通用值**；本项目实测或推导出来的数（速率、泄漏、诱导幅度、
  基因组大小、窗口代数、保护区、剂量、F）**不得作为默认值、示例值或"典型值"出现**。
- 落点：代码缺参就拒绝（不回落内置数）；示例只用文献数据集或明确标注的合成数据；文档算例用符号。
- `selftest.py` 有「通用性」检查会扫描全文并判失败——新增内容前先跑它。

**边界（本技能不答，被问到就明说）**："耐药会快几倍"（研究问题）、"临床算不算耐药"
（断点问题）、"乘客里哪个是有害的"（WGS 功能注释，交给基因组管线）。

## 铁律

- 在技能目录里执行 `run.py`；不要假设 `windowcalc` 在 PATH 上。
- 首次慢（20–60 s，装依赖）是正常的；只有打印“无法准备运行环境”才算失败。
- 失败时原样转达它给的三条救援（镜像 / `WINDOWCALC_PYTHON` / 手动安装），不要自创方案。
- 不在回答里汇报“我装了依赖 / 建了 venv”，除非失败需要用户决策。
- 排查：`python3 run.py --doctor`（0=可用，2=首次自动装，1=需人工）；沙箱不能执行 shell 时让用户跑一次并贴回输出。
- **用户没有菌株卡时，不要让用户写 YAML**：按 `references/intake.md` 用 6 个带默认值的问题问齐，agent 自己写卡；协议卡默认用内置模板（`optimize` 不带 `--protocol`），只在用户有明确要求时改 k/成功率/上限。
- **μ_supply 没测时不要只填默认**：先给测量方案（12 管试点 + MIC/背景，见 `docs/experiment-design.md`），或明确说“默认值下窗口只作方向参考”；不能把 assumed 当实测。
- **跑任何子命令前，先过 `references/guardrails.md` 的必问清单**：多板/整管涂板、Nt 单位、频率 vs 速率、质粒 mutator 角色、保护上限取值——用户不会主动说，但都会改变结论。
- **主动查外推，不等用户提**：用户只要说“我测过突变率/我有 μ”，就先问两句——**这个速率是用哪个标记/药物测的？窗口里打算用哪个标记筛？** 两者不同（哪怕用户没意识到）→ 按 `references/markers.md` 走 4 问检查清单；不得直接把 μ_X 填成 `mu_supply`。
- **外推必须给结论**：明确说“可以给先验 / 只能当量级参考 / 不可外推”三选一，并给出下一步（试点或实测）；未确认目标标记类型前，不得继续算窗口。
- **数值参数必填**：`power` 的 `--mu/--mu1/--nt-mean`、`fit-tauoff` 的 `--td`、`design-induction` 的 `--nt`
  都是必填——本技能不内置任何项目实测默认值（通用性政策），缺了就报缺，不许猜。
- **缺参不得自行填默认**：必须把“缺什么 / 默认值多少 / 用默认的后果 / 怎么拿到”原样呈现，**等用户明确同意**后才加 `--fill-defaults`；agent 不得替用户同意。
- **多靶点只有一条计算路径**：`multitarget --records …`（供给=铺板实测）+ 标定链只算乘客（M = μ_perbp × G）。**禁止用 μ×T_set 外推供给**——几何路线与 `--mode design/fit/classify` 已于 2026-11-13 删除；F 不进窗口计算，只做 QC 与 22×/139× 代价判断；谷型（C 律）已于 2026-11-12 移除。

## 先判断用户要哪一类

| 用户说 | 做 | 细节 |
|---|---|---|
| **「看不懂 / 太复杂 / 我是新手」** | **不要解释术语**；直接把 `docs/start-here.md` 的三个问题问出来 | 该文档 |
| 给了 r/Nt 计数表 | `rate` | `references/commands.md` |
| 两批速率有没有差 | `compare` | 同上 |
| 要铺多少管 / 功效 | `power`（先确认 1σ 还是 95% CI 口径） | 同上 |
| 测什么 / 实验怎么设计 | 按 `docs/experiment-design.md` 给完整设计 | 该文档 |
| 窗口 / 表达多久 / 乘客预算 | **先过 j 三问**（表后）→ j=1：`references/intake.md` 问齐 → `optimize`；j≥2：`references/multitarget-intake.md` 问齐 → `multitarget --records`（每档结束重跑重算 j） | `docs/mutator-window.md` / `docs/protocol-ramp.md` |
| 剂量响应 / τ_off | `fit-fc` / `fit-tauoff` | `references/commands.md` |
| 概念 / 原理 / 设计逻辑 / 选估计器 | 读 `docs/knowledge-base.md` 后回答（开篇有定义与“七个为什么”） | 该文档 |
| 用一种标记测速、想外推到另一种 | 先读 `references/markers.md`，走 4 问检查清单；**不得直接把 μ_X 填成 μ_supply** | 该文档 |
| **多靶点/膜活性药：渐升怎么规划、乘客多少** | **先过 j 三问**（表后）→ 按 `references/multitarget-intake.md` 问齐实测值（**PAP 板 / 标定标记 / k·capacity·Nt**）→ `multitarget --records …`（乘客速率：`--markers` 标记链 或 `--m-on` 直给）；先过 `references/guardrails.md` §4（G16–G20）。**没有账本时不要直接要 `--records`**——先按 `docs/protocol-execution.md` §1.5 的 S0/S1 清单（MIC → 基线梯度 → c₁ → 首窗 → 7 块板） | `docs/protocol-ramp.md` |
| **多靶点：试管与平板之间怎么倒（涂多少 / 剩的去哪 / 板还是液传代 / 能不能不等板）** | 按 `docs/protocol-execution.md` 的决策树给：涂板比例由克隆数反算、剩余五去处、拓扑二选一、并行四条护栏 | `docs/protocol-execution.md` |
| **PAP 是什么 / 怎么测 / 有文献吗 / 板没出结果怎么预测下一档 / 不涂板直接液接液行不行** | 按 `references/faq-pap.md` 五问逐条答：Q1 定义与 δ 含义 · Q2 三个测量时机与板量 · Q3 文献锚点（含"尾部 log-线性未钉 DOI"的坦白）· Q4 预测/转接/审计三时刻拆解 · Q5 液接液✓与零涂板的五项损失 + 三档配置 | `references/faq-pap.md` |


### j 三问——任何窗口设计前必过（决定实验形状）

> 对新手只说一句："**先别挑命令——先回答：你要跨几档？**" 这决定实验是
> 一次诱导脉冲（j=1）还是"每档 pre/post 配对板 + 逐档审计"的渐升（j≥2），
> 两者的人力、板量、可审计性完全不同。

| 问 | 不知道时 |
|---|---|
| ① 你想把**有效浓度**提高几倍？（目标浓度 ÷ 现在能压住的浓度） | 必答 |
| ② 这类药「**升一档**」一般能提高几倍？（查同源药文献；脂肽常见 20–30 倍） | 查不到 → **第一档当预实验**，测出来再算 |
| ③ 要跨几档：j = ⌈log(①)/log(②)⌉ | j=1 → 一个窗口搞定（`optimize`）；j≥2 → 逐档账本（`multitarget`）；**每档结束重算 j** |

②取文献值要打折：可及路径可能 <20%（ncomms13002 五突变仅 22 条可行）、
mutator 会改道（pbio.3001878）——单步效应按实测走，文献只定方向。

## 9 个子命令（速查）

| 子命令 | 用途 | 关键旗标 |
|---|---|---|
| rate | 计数 → μ | `--field`（M_base/M_on_ref/M_off/mu_supply）、`--nt-mode per_tube`、`--fraction-plated` |
| compare | 两批比较 | `--d1/--d2`、`--labels`、`--nt1/--nt2` |
| power | 样本量 / 功效 | `--mode`（precision/power）、`--conf-level 0` = 1σ |
| optimize | 窗口三旋钮 | `--strain`（必填）+ `--protocol`（可省）；`--fill-defaults` 只补 4 项约定量 |
| design-induction | 剂量轴布点 | `--strain` / `--c50` / `--doses` |
| fit-induction | 剂量机制 | `--m-off`、`--query` |
| fit-fc | 剂量响应 | `--mu-base` |
| fit-tauoff | 衰变 | `--mu-off`、`--td` |
| multitarget | 多靶点窗口：f_step 实测供给 + 乘客标定链（M = μ_perbp × G） | `--records --markers --k --nt --genome-mb [--bgc-kb]`；缺参**不跑**、不得编造 |

**rate 输入**：每行一管 `r[,Nt[,volume_mL]]`；可读 `.csv/.txt/.xlsx`（第一个工作表）；
表头可留；第 1 列 r、第 2 列 Nt；Nt 是 **CFU/管**（CFU/mL 用 `--nt-unit per_mL --volume-mL`）。
输出先回显输入，再给 m̂±95%CI、μ、p₀ 交叉核验与 YAML 片段。

⚠ 两条必须转达：

1. 用 `--nts` 未给 `--nt-unit` 时按 CFU/管 解释——CFU/mL 差一个体积倍，是 μ 最常见的错；
2. YAML 片段是**单位点**速率，不能直接粘进 `M_base`（模型 M 是全基因组每代，差 10³–10⁶ 倍）。

## 报告口径（转述结果时必须遵守）

- 报**速率**不报频率：r/Nt 受 L-D 重尾支配，不可复现。
- **标记无关**：任何“可数、低背景”的选择标记都行（Rif^R、Str^R…）；标记只影响生物学解释（靶标大小、T_eff、背景），不影响 m/μ 的算法。
- `rate` 的 μ 是**单位点**；`M_base/M_on/M_off` 是**全基因组每代**；只有比值可跨量纲。
- 样本量两套口径不可混：`--conf-level 0` 是 1σ 相对标准误；默认 0.95 是“95% CI 半宽”，同样写 20% 差 3.84 倍。
- 丝状菌（如放线菌）：μ 是 “per CFU per passage” 近似；跨文献只比倍数 F，不比绝对值。
- 单次实验噪声大（C=40、m≈2 时 CV≈14%）；判差异看 `compare` 的 LRT p 与 fold-change CI，不看点估计是否接近。
- 红灯（泄漏超标、背景淹没、τ_off 过长、边界解）必须转达，不得省略。
- **保护优先（硬约束）**：两条独立上限——`P_mutator_max`（默认 0.05，mutator 自失活）与 `P_protected_max`（默认 0.10，BGC/其他保护区）；k 若需要更长的窗口，**按上限截断**，报告实际 events、缺口与哪条 binding，**不得为了 k 超窗**。
- 用户问“要铺多少管”时给完整设计（C、N₀、板浓度、计数时机、对照、记录表）；**没有 m 先验就用默认 C=30**（m 1–3 够用），或先铺 12 管试点。

## 常见坑（Top 10）

1. `--nt-mode per_tube` 配标量 Nt 会报错（故意的）；CSV 每行与 Nt 必须一一配对。
2. 部分涂板只用 `--fraction-plated`（Nt 填整管总 CFU）；手工 r/ζ 会系统性偏低。
3. m<0.3 不可靠（R&F 建议 ≥0.3）；m 搜索域上限 300（触及给警示）；r>5000 直接报错。
4. jackpot 管（r 上千）先排饱和/污染；部分涂板叠加 jackpot 会报错。
5. `target_rse` 与 95% CI 半宽别混（差 z²）。
6. MC 两个区间别混：参数不确定性 vs 再叠加克隆离散（后者补测不会收窄）。
7. 无 CI 的参数按 ±2× 对数正态假设；要真区间就在卡片写 `ci`。
8. 边界解（剂量/T/N₀ 顶到网格端点）先放宽范围再算。
9. 复现看报告附录 A（版本/网格/seed）。
10. 对外模拟（`power --sim` / `compare`）克隆封顶 5000；内部性质测试的 `simulate` 不封顶（已注明仅测试用）。
11. **频率数据别喂 `rate`**：单次铺板的克隆频率（f_step 升档板）不是平行培养计数，MSS-MLE 不适用——把铺板记录写进账本（cycle,c,when,cfu,clones），`multitarget` 自动做背景扣除（含 Wilson/Newcombe 区间）——规程：`docs/protocol-fstep.md`。

## 细节在哪

**入口只有一个：`docs/INDEX.md`**——14 篇文档按"我现在要干什么"路由；
每篇开头都有一行「级别 ｜ 什么时候读我 ｜ 不用读我」，按那行判断要不要继续看。

- **必读 3 篇**：`docs/tutorial.md`（算突变率）· `docs/protocol-execution.md`（走一个浓度周期）·
  `docs/multitarget-code-map.md`（维护者：改公式/审代码 → 改一处要改两处）
- **主路径 4 篇**：`docs/start-here.md`（多靶点 5 分钟判断）· `docs/experiment-design.md`（设计波动实验）·
  `docs/protocol-fstep.md`（测 f_step）· `docs/protocol-ramp.md`（排阶梯与乘客账）
- **附录 7 篇**：`knowledge-base.md`（含 §7 μ→M 换算、§8 计数诊断）·
  `multitarget-model.md`（含 §9 分类轴、§10 证据分级）· `mutator-window.md` · `glossary.md` ·
  `multitarget-tutorial.md`（**命令已撤，只读概念**）· `generic-policy.md` · `testing.md`（含历史核对）
- `docs/figures/`：三张符号化配图（**不含任何数值**，生成脚本同目录、可复现）——
  周期结构 / 两种传代拓扑 / 两个决策。用户问"一个浓度周期长什么样"时，展示
  `figures/cycle-liquid-propagation.png` 并按 `protocol-execution.md` 解读。
- **`references/`（agent 专用）**：`commands.md`（9 个子命令逐条）· `guardrails.md`（必问清单）·
  `markers.md`（标记外推）· `multitarget-intake.md`（多靶点接入）· `faq-pap.md`（PAP 五问：定义/测量/文献/时序/不涂板的代价）·
  `guardrails.md` §0（命令→该查哪几条）·
  `install.md` / `maintainer.md`（部署与发布）
- 维护者自检 **28 项：环境**、入口自举、规范校验、**9 个子命令**、Excel、已知答案、端到端报告、文档一致性、错误路径、断链、多靶点缺参拦截、**通用性（无项目专属药物名/参数）**、**文档分级（横幅 + INDEX）**、**护栏覆盖（每个子命令都在 guardrails §0）**、**护栏反查（每个 G 编号都挂在 §0，无孤儿规则）**、**文档↔代码（参数名/常量/声称一致）**。
