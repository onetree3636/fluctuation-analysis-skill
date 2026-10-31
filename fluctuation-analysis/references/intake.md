# 菌株卡/协议卡：怎么问用户（agent intake）

> 目标：用户**不用写 YAML、不用看协议卡**。agent 用几个带默认值的问题问齐，自己写卡、自己跑。

## 0. 三条原则

1. **不要让用户编辑 YAML**——问问题，agent 写卡。
2. **"不知道"是合法答案**——用默认值，报告会标 ★未测；不要因为缺参数就卡住。
3. **协议卡默认用内置模板**：`optimize --strain card.yaml`（不带 `--protocol`）直接跑；只问 4 个决策值。

## 1. 必问（6 个，带默认）

| # | 问用户 | 字段 | 默认（不知道时） |
|---|---|---|---|
| 1 | 基因组多大（Mb）？ | `genome_mb` | 必填（没有就不算保护区账本） |
| 2 | 要保护的产物基因簇多大、几拷贝？ | `protected[{role: bgc, kb, copies}]` | 可跳过 |
| 3 | mutator 在质粒上吗？质粒多大、拷贝数多少？ | `plasmids[{kb, copies}]` + `protected[{role: mutator}]` | 可跳过 |
| 4 | 这个速率是**用哪个标记/药物**测的？窗口里**打算用哪个标记**筛？ | `mu_supply` + `mu_supply_source` / `selection_marker` | 两个都写进卡片；**不同 → 按 `references/markers.md` 处理**；没测 → 先给测量方案（见 1.1）；只想要量级才用 2×10⁻⁸ |
| 5 | WT 的全基因组乘客速率有吗（M_base）？ | `M_base` | **0.004**（Drake 量级，assumed） |
| 6 | 诱导/关闭时间常数、倍增时间？ | `tau_on_h` / `tau_off_h` / `td_h_on` / `td_h_off` | **2 h / 4 h / 2.5 h / 2.0 h** |

### 1.1 μ_supply 没测怎么办（不要直接填默认）

μ_supply 是窗口的**供给端主参数**：`events ∝ μ_supply`，达到 k 所需的窗口长度 ≈ 1/μ_supply。猜错 10× → 窗口、事件数、保护上限命中都会偏。所以：

- **优先先测**：按 `docs/experiment-design.md` 做一轮波动分析（12 管试点 → 正式 C），选择板就是实际要用的标记与浓度；
- 测完 `rate --field mu_supply` 回填卡片，再跑 `optimize`；
- **暂时不测**：可以用默认 2×10⁻⁸ 做“量级演练”，但结论必须写“μ_supply 未测（assumed），窗口只作方向参考”；
- 不要把它当成和 `M_base` 一样的“整体 scale”——M 的整体缩放不改变窗口 argmin，μ_supply 会改变。
- 如果用户是用**另一种标记**测的（如四环素），先看 `references/markers.md`：目标标记是单拷贝蛋白才可按 T_eff 比例给先验；rRNA 多拷贝/酶灭活/外排必须重测。

## 2. 选问（用户主动提供或已实测才问）

| 问 | 字段 | 默认 |
|---|---|---|
| mutator 全开时的速率？ | `M_on_ref` | 100×M_base |
| off 态泄漏？ | `M_off` | 3×M_base |
| 饱和密度？ | `Nt_saturation` | 3×10⁸ CFU/管 |
| 参考诱导剂量？ | `ref_dose_mM` | 0.3 mM |
| 假抗性背景？ | `background_false` | 不填（不做红灯检查） |
| 有 CI 吗？ | 各字段 `ci: [lo, hi]` | 不填（MC 按 ±2× 对数正态） |

## 3. 协议卡只问 4 个

| 问 | 字段 | 默认 |
|---|---|---|
| 想拿几个独立事件？ | `k_target_clones` | 12 |
| 成功率下限？ | `P_success_min` | 0.95 |
| 两条保护上限？ | `P_protected_max` / `P_mutator_max` | 0.10 / 0.05 |
| 脉冲时长 / 接种量范围？ | `T_pulse_range_h` / `N0_min` / `N0_max` | 3–24 h / 1e6–1e8 |

其余（相位表、doses、MC 次数）用内置模板；用户不需要知道。

## 4. 写卡与运行

- 回答齐 → agent 写 `<strain>.yaml`；缺的必填项可以两个办法：
  - 写全（推荐，`status: assumed` 标清楚）；
  - 或直接 `optimize --strain card.yaml --fill-defaults`（缺项自动补默认，stderr 会列出补了哪些）。
- 协议卡：`optimize --strain card.yaml` 即可（不带 `--protocol`）。
- 跑完把 `保护上限 binding` / `红灯` / `★未测` 三件事如实转告。

## 5. 对话示例

> **用户**：帮我设计一个 mutator 窗口，拿 12 个克隆。
> **agent**：先确认 4 件事（不知道就说不知道）：① 基因组多大？② 产物基因簇多大？③ mutator 在质粒上吗、拷贝数多少？④ 筛选标记的速率测过吗？
> **用户**：基因组 8.4 Mb；簇 100 kb；质粒 16 kb、20 拷贝；速率没测。
> **agent**：好，μ_supply 先用 2×10⁻⁸（标 assumed），M_base 先用 0.004，其余用默认；协议卡用内置的（k=12、成功率 0.95、两条上限 10%/5%）。跑完告诉你实际事件数、保护上限有没有 binding、哪些是未测参数。
