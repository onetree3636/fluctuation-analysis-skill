# 文档索引——先看这一页

> 19 篇文档不是 19 个入口。**入口只有这里。**
> 每篇文档开头都有一行「级别 / 什么时候读我 / 不用读我」，按那行判断要不要继续看。

---

## 一、按"我现在要干什么"选

| 我现在的处境 | 读这几篇（按顺序） |
|---|---|
| **手上有 r/Nt 计数表，要算突变率** | `tutorial.md` → 需要原理再看 `knowledge-base.md` |
| **要设计一次波动实验（多少管、几个剂量）** | `experiment-design.md` |
| **计数不对劲（回收低、成团、对不上）** | `viable-count-problems.md` |
| **要开 mutator 拿耐受克隆（单步）** | `start-here.md` → `mutator-window.md`（原理）→ 跑 `optimize` |
| **要多靶点渐升（j≥2）** | `start-here.md` → `protocol-ramp.md`（排阶梯）→ `protocol-execution.md`（怎么走一个周期） |
| **要测 f_step（配对铺板）** | `protocol-fstep.md` |
| **某个药能不能外推到另一个药** | `multitarget-extrapolation.md` → `resistance-evidence-tiers.md` |
| **单位点 μ 怎么换成全基因组 M** | `mu-genome-derivation.md` |
| **术语看不懂** | `glossary.md` |
| **我要改代码 / 审模型 / 加内容** | `multitarget-code-map.md` → `generic-policy.md` → `testing.md` |

---

## 二、级别一览（三档）

### 必读（3 篇）

| 文档 | 一句话 | 什么时候**不用**读 |
|---|---|---|
| `tutorial.md` | 从 40 管计数到突变率，10 分钟 | 你已经会跑 `rate` |
| `protocol-execution.md` | 一个浓度周期怎么走：涂多少、传什么、何时能往下走 | 你只做单步窗口、不做渐升 |
| `multitarget-code-map.md` | 模型/公式/算法落在哪个函数（维护者必读） | 你只是使用者 |

### 主路径（4 篇）

`start-here.md`（多靶点 5 分钟判断）· `experiment-design.md`（实验设计）·
`protocol-fstep.md`（测 f_step）· `protocol-ramp.md`（排阶梯与乘客账）

### 附录（12 篇）

原理：`knowledge-base.md` · `mutator-window.md` · `multitarget-model.md` · `mu-genome-derivation.md`
判断：`multitarget-extrapolation.md` · `resistance-evidence-tiers.md`
诊断：`viable-count-problems.md`
工具书：`glossary.md`
概念演示（命令已撤，只读概念）：`multitarget-tutorial.md`
维护者：`generic-policy.md` · `testing.md`
历史：`software_audit_2026-10-30.md`

---

## 三、两条纪律（写文档前先读）

1. **能塞进已有文档的一节，就不新开一篇。** 19 篇已是上限。
2. **新增/改动前先问三句**：它让哪一次实验决策更省或更准？没有它用户会犯什么具体错？
   能不能放进已有文档？——答不上就别加（详见 `generic-policy.md` §5）。

---

## 四、代码与图在哪

- 命令用法：`references/commands.md`（9 个子命令逐条）
- 领域护栏（agent 必问清单）：`references/guardrails.md`
- 多靶点接入：`references/multitarget-intake.md` · `references/multitarget-guardrails.md`
- 三张符号化配图（不含任何数值）：`figures/`——周期结构、两种传代拓扑、两个决策
