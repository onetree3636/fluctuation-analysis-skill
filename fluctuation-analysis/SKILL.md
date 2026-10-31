---
name: fluctuation-analysis
description: 把波动实验的菌落计数变成带 95%CI 的突变率（MSS-MLE，支持逐管 Nt 与部分涂板），并做两批比较、样本量/功效、mutator 表达窗口、剂量响应与 τ_off 拟合。用户给出 r,Nt 计数、问“这两批有没有差 / 要铺多少管 / 窗口怎么定 / 剂量响应怎么拟合”，或要设计波动实验、问波动分析原理时使用。
license: MIT
compatibility: Python 3.9+；Linux/macOS/Windows；首次运行自动建 .venv 装 numpy/pyyaml/openpyxl，不需要 R。
allowed-tools: Bash Read Write Edit Glob Grep
---

# 波动分析与 mutator 窗口

入口只有一条：`python3 run.py <子命令>`（Windows：`run.cmd <子命令>`）。
第一次运行会自动准备环境；**不要让用户 `pip install`，不要运行 `install.sh`**。

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

## 先判断用户要哪一类

| 用户说 | 做 | 细节 |
|---|---|---|
| 给了 r/Nt 计数表 | `rate` | `references/commands.md` |
| 两批速率有没有差 | `compare` | 同上 |
| 要铺多少管 / 功效 | `power`（先确认 1σ 还是 95% CI 口径） | 同上 |
| 测什么 / 实验怎么设计 | 按 `docs/experiment-design.md` 给完整设计 | 该文档 |
| 窗口 / 表达多久 / 乘客预算 | 先按 `references/intake.md` 问齐菌株参数 → `design-induction` → `fit-induction` → `optimize` | `docs/mutator-window.md` |
| 剂量响应 / τ_off | `fit-fc` / `fit-tauoff` | `references/commands.md` |
| 概念 / 原理 / 设计逻辑 / 选估计器 | 读 `docs/knowledge-base.md` 后回答（开篇有定义与“七个为什么”） | 该文档 |
| 用一种标记测速、想外推到另一种 | 先读 `references/markers.md`，走 4 问检查清单；**不得直接把 μ_X 填成 μ_supply** | 该文档 |

## 8 个子命令（速查）

| 子命令 | 用途 | 关键旗标 |
|---|---|---|
| rate | 计数 → μ | `--field`（M_base/M_on_ref/M_off/mu_supply）、`--nt-mode per_tube`、`--fraction-plated` |
| compare | 两批比较 | `--d1/--d2`、`--labels`、`--nt1/--nt2` |
| power | 样本量 / 功效 | `--mode`（precision/power）、`--conf-level 0` = 1σ |
| optimize | 窗口三旋钮 | `--strain`（必填）+ `--protocol`（可省）；`--fill-defaults` 补缺项 |
| design-induction | 剂量轴布点 | `--strain` / `--c50` / `--doses` |
| fit-induction | 剂量机制 | `--m-off`、`--query` |
| fit-fc | 剂量响应 | `--mu-base` |
| fit-tauoff | 衰变 | `--mu-off`、`--td` |

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

## 细节在哪

- `references/commands.md`：子命令完整用法、输入格式、输出与边界、实测数字。
- `references/install.md`：首次运行、venv/镜像、Windows、导入其它 agent、故障救援。
- `references/maintainer.md`：selftest、sync_skill、与 rSalvador/flan 对拍、发布检查单。
- `references/guardrails.md`：新手盲区排查表（哪些话用户不会主动说、必须问）。
- `references/markers.md`：标记外推（哪种能按 T_eff 换算、哪种必须实测）。
- `docs/experiment-design.md`、`docs/mutator-window.md`、`docs/knowledge-base.md`、`docs/tutorial.md`、`docs/glossary.md`。
- `docs/viable-count-problems.md`：活菌计数不准（随机误差/低回收/团块/VBNC/超范围）的诊断与校正/替代计数。
- 维护者自检 **21 项：环境**、入口自举、规范校验、**8 个子命令**、Excel、已知答案、端到端报告、文档一致性、错误路径、断链。
