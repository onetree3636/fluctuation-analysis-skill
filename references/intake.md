# 菌株卡/协议卡：怎么问用户（agent intake）

> 目标：用户**不用写 YAML、不用看协议卡**。agent 用几个带默认值的问题问齐，自己写卡、自己跑。

## 0. 三条原则

1. **不要让用户编辑 YAML**——问问题，agent 写卡。
2. **"不知道"是合法答案，但只有 4 项能落到默认值**（`M_base`/`M_on_ref`/`M_off`/`Nt_saturation`）。
   另外 6 项（`mu_supply`、τ_on/τ_off、td_on/td_off、`ref_dose`）是**你体系的实测物理量**——
   不知道就明说，让用户给个数或标 `status: assumed`，**不要替他编**。
3. **协议卡默认用内置模板**：`optimize --strain card.yaml`（不带 `--protocol`）直接跑；只问 4 个决策值。

## 1. 必问（6 个，带默认）

| # | 问用户 | 字段 | 默认（不知道时） |
|---|---|---|---|
| 1 | 你的菌基因组大约多少 Mb？（查已发表基因组即可） | `genome_mb` | 必填（没有就不算保护区账本） |
| 2 | 要保住的产物基因簇有多大（kb）？在染色体上几拷贝？ | `protected[{role: bgc, kb, copies}]` | 可跳过 |
| 3 | 你的 mutator 是**质粒上**还是染色体上？质粒的话多大、几拷贝？ | `plasmids[{kb, copies}]` + `protected[{role: mutator}]` | 可跳过 |
| 4 | 这个速率是用**哪个药**筛出来的？窗口里打算用**哪个药**筛？（两件事，也要问两遍） | `mu_supply` + `mu_supply_source` / `selection_marker` | 两个都写进卡片；**不同 → 按 `references/markers.md` 处理**；没测 → 先给测量方案（见 1.1）；只想要量级 → 引一条**文献**的同类体系数值，并在卡片里写出处（标 assumed） |
| 5 | WT 的全基因组乘客速率有吗（M_base）？ | `M_base` | **0.004**（Drake 量级，assumed） |
| 6 | 细菌在你手里**几小时分裂一次**？去掉诱导剂后，mutator 大约**多久关掉**？ | `tau_on_h` / `tau_off_h` / `td_h_on` / `td_h_off` | **必填**：没测过倍增时间 → 用你自己生长曲线的对数期；τ_off 没测 → 先按 1.44×倍增时间 起，标 assumed |

### 1.1 μ_supply 没测怎么办（不要直接填默认）

μ_supply 是窗口的**供给端主参数**：`events ∝ μ_supply`，达到 k 所需的窗口长度 ≈ 1/μ_supply。猜错 10× → 窗口、事件数、保护上限命中都会偏。所以：

- **优先先测**：按 `docs/experiment-design.md` 做一轮波动分析（12 管试点 → 正式 C），选择板就是实际要用的标记与浓度；
- 测完 `rate --field mu_supply` 回填卡片，再跑 `optimize`；
- **暂时不测**：可以引一条**文献**里同类体系的速率做“量级演练”（卡片里写清出处、标 assumed），但结论必须写“μ_supply 未测（assumed），窗口只作方向参考”。**不要用任何内置数字**——本技能不内置它；
- 不要把它当成和 `M_base` 一样的“整体 scale”——M 的整体缩放不改变窗口 argmin，μ_supply 会改变。
- 如果用户是用**另一种标记**测的（如四环素），先看 `references/markers.md`：目标标记是单拷贝蛋白才可按 T_eff 比例给先验；rRNA 多拷贝/酶灭活/外排必须重测。

## 2. 选问（用户主动提供或已实测才问）

| 问 | 字段 | 默认 |
|---|---|---|
| mutator 全开时的速率？ | `M_on_ref` | 100×M_base |
| off 态泄漏？ | `M_off` | 3×M_base |
| 饱和密度？ | `Nt_saturation` | 3×10⁸ CFU/管 |
| 参考诱导剂量？ | `ref_dose_mM` | **必须你给**（诱导剂不同差几百倍，不能内置） |
| 假抗性背景？ | `background_false` | 不填（不做红灯检查） |
| 有 CI 吗？ | 各字段 `ci: [lo, hi]` | 不填（MC 按 ±2× 对数正态） |

## 3. 协议卡只问 4 个

| 问 | 字段 | 默认 |
|---|---|---|
| 这一轮想要几个**独立的抗性克隆**（不是重复验证同一个）？ | `k_target_clones` | 12 |
| 有多大把握拿到这么多？（把握越大，需要的细胞越多） | `P_success_min` | 0.95 |
| 允许误伤多少？（产物基因簇 / mutator 自身） | `P_protected_max` / `P_mutator_max` | 0.10 / 0.05 |
| 一次诱导最多开多久 / 起始细胞数大致范围？ | `T_pulse_range_h` / `N0_min` / `N0_max` | 3–24 h / 1e6–1e8 |

其余（相位表、doses、MC 次数）用内置模板；用户不需要知道。

## 4. 写卡与运行

- 回答齐 → agent 写 `<strain>.yaml`；缺的必填项可以两个办法：
  - 写全（推荐，`status: assumed` 标清楚）；
  - 或 `optimize --strain card.yaml --fill-defaults`——**只补 4 项与体系无关的量**
    （M_base / M_on_ref / M_off / Nt_saturation）；`mu_supply`、τ_on/τ_off、td_on/td_off、
    `ref_dose_mM` 这 6 项是你体系的实测物理量，它**不会**替你编，缺就报错（政策规则二）。
- 协议卡：`optimize --strain card.yaml` 即可（不带 `--protocol`）。
- 跑完把 `保护上限 binding` / `红灯` / `★未测` 三件事如实转告。

## 5. 对话示例

> **用户**：帮我设计一个 mutator 窗口，拿 12 个克隆。
> **agent**：先确认 4 件事（不知道就说不知道）：① 基因组多大？② 产物基因簇多大？③ mutator 在质粒上吗、拷贝数多少？④ 筛选标记的速率测过吗？
> **用户**：基因组 8.4 Mb；簇 100 kb；质粒 16 kb、20 拷贝；速率没测。
> **agent**：好。μ_supply 你没有实测——要么先做一轮 12 管试点，要么给我一条**文献**里同类体系的速率（我标 assumed 并写出处）。M_base 先用 0.004（Drake 量级），其余 4 项用默认；协议卡用内置的（k=12、成功率 0.95、两条上限 10%/5%）。跑完告诉你实际拿到几个事件、保护上限有没有被顶住、哪些参数还是未测。
