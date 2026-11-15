# Changelog

## 未发布（工作区，2026-11-14 起）

### Added（2026-11-15，把 F 归位、把 M_on 的入口打开）

设计问题：多靶点路径里 F 到底要不要测？查代码发现 **F 一次都没进计算**，
而 `μ` 是结构上必填（只能走标记链）。两个改动让口径与代码一致：

- **`--m-on`：直接给全基因组乘客速率**（`/基因组/代`），跳过标记链。
  与 `--markers` **二选一**（同时给 / 都不给都直接报错，并说明两者区别）。
  给 WGS 突变累积测出 M_on 的实验室一条路——不必为了凑标记去补波动实验，
  也不会把 per-bp 的量误乘一遍 G。
- **F 归位到"代价"**：markers 若带 `f` 列，`multitarget` 报告末加一行
  "代价参考（F 只做 QC/表述，**不进供给**）：markerA F=30× → …"，
  调用此前**没有调用者**的 `f_ceiling_note()`（22× 最优 / 139× 变差）。
  报告开头新增"乘客速率 M_on = …（来源：标记链 μ_perbp×G / 直接给定）"。
- 测试 192 → 197：`--m-on` 直给、与 `--markers` 互斥、两者都缺、非正值拒绝、
  markers 带 f 时出代价提示。
- 文档：`protocol-ramp.md` 命令示例与参数、`multitarget-code-map.md` §2.4、
  SKILL.md 路由行同步。

### Changed（2026-11-15，文档合并：19 → 14 篇）

按主题合并 5 篇，**标题、公式、结论一条未删**（并入后逐篇核对：标题 100%、公式/数值行 100%）：

| 原文件 | 并入 | 位置 |
|---|---|---|
| `multitarget-extrapolation.md` | `multitarget-model.md` | §9 分类轴 |
| `resistance-evidence-tiers.md` | `multitarget-model.md` | §10 证据分级 |
| `mu-genome-derivation.md` | `knowledge-base.md` | §7 μ→M 换算 |
| `viable-count-problems.md` | `knowledge-base.md` | §8 计数诊断 |
| `software_audit_2026-10-30.md` | `testing.md` | §9 历史核对 |

- 合并方式：源文档**整篇**挂在目标的新 H2 下、标题降一级，正文逐字保留；每段开头注明
  "2026-11-15 由 `X.md` 整篇并入，内容未删"。
- 目标文档的标题与分级横幅同步（`multitarget-model.md` → "数学模型、分类轴与证据分级"；
  `knowledge-base.md` → "含 μ→M 换算与计数诊断"）。
- 全仓引用改写（README / SKILL / INDEX / 各文档 / references），指向合并后的新位置与小节号。
- `INDEX.md` 重写：路由表按"我要干什么"给到小节（如"计数不对劲 → knowledge-base §8"），
  分级改为必读 3 / 主路径 4 / 附录 7。

### Changed（2026-11-15，文档减法）

- **新增 `docs/INDEX.md`：文档的唯一入口**。19 篇文档按"我现在要干什么"路由
  （算突变率 / 设计实验 / 开 mutator / 多靶点渐升 / 测 f_step / 外推 / 改代码…），
  并给出三档分级与"什么时候**不用**读"。
- **每篇文档 H1 下方加分级横幅**：`> **级别：必读/主路径/附录** ｜ 什么时候读我 ｜ 不用读我`。
  必读 3 篇（tutorial / protocol-execution / code-map）、主路径 4 篇、附录 12 篇。
- **SKILL.md 的"细节在哪"从 20 条平铺改为三级结构**（必读 3 / 主路径 4 / 附录 12 + references 一组），
  并明确"入口只有一个：docs/INDEX.md"。
- **selftest 新增第 25 项**：每篇 docs/*.md 必须有横幅且级别合法、INDEX 必须存在
  ——防止以后又长出平行入口。SKILL.md / testing.md 的项数同步。
- 内容一篇未删：只做索引与分级。

### Added（2026-11-15，传代拓扑参数化）

- **`--propagation {manual,plate,liquid}`：档间代数由拓扑自己算**，不再让用户手工折算。
  - `manual`（默认，向后兼容）= 沿用 `--n-plate`；
  - `plate` = `log₂(--cells-per-colony)`（菌落大小是实测/假设值，技能不内置）；
  - `liquid` = `log₂(N₀/(f·Nt))`，需 `--nt`；N₀ 缺省取"刚好够 k"的 `Nt − k/f`，
    `--n0` 可覆盖（给得太大时该档拿不到 k）。
  - 缺参直接报错（如 liquid 无 `--nt`、plate 无 `--cells-per-colony`），不猜默认值。
- **已完成档按各自周期的 f 算档间代数**（不是一刀切）：`_cycle_f()` 取该周期尾部
  （`f < sat_f`）里**最大**的 f —— 对应"你实际升到的那一档浓度"；
  用"最高浓度的 f"会把 N₀ 算成负数（纯检测点的 f 更小）。
- 报告新增一行："档间扩增：liquid 模式 → N 代/档（来源公式）"。
- 测试 187 → 192：manual 直通、plate 的 log₂、liquid 的公式/`--n0`/缺 `--nt`/N₀ 不可行、
  **两周期用各自 f** 的 base 累加、CLI `liquid` 缺 `--nt` 拒绝。

### Fixed / Added（2026-11-15，收紧三项）

- **P0（数字会错）：`fit_pap` 的 δ 95%CI 上下颠倒**。`δ = −1/s` 在 s<0 时随 s 单调增，
  下端点应由 `s_lo` 算，此前写成 `(−1/s_hi, −1/s_lo)`，报告里会打印成 `δ 95%CI [大, 小]`
  （实测例 `[1.26, 1.06]`）。已改正并加断言 `lo < δ < hi`；斜率 CI 含 0 时另给一条 flag。
- **洗脱尾进账本**：新增 `ramp.washout_tail()` 与 `multitarget --tau-off-h/--td`。
  按 `(M_on−M_off)/(e^(td/τ)−1)` 计，**每次洗脱一次**（已完成档基数 + 每剩余档 + 终点后）；
  只给一个参数直接报错（不猜默认值），都不给则报告明写"未计入"。
  连续极限退化为 `(M_on−M_off)·τ/td`，与 `fit-tauoff` 的尾部积分同口径。
- **饱和阈值从硬编码改为参数**：`--sat-f`（默认 0.5 = "长到与未处理相当"），
  报告第二节显示所用阈值；落在阈值 **±20% 带**内的点单独标出（分类对阈值敏感）。
  此前 `f ≥ 0.5` 散落在 `run_ledger` 两处，边界附近的结论会静默跳变。
- **测试**：对拍文件从 7 项扩到 14 项——洗脱尾闭式解、连续极限、参数必须成对、
  每档计尾一次、饱和阈值可调（同一条记录在 0.5 / 0.7 下分别"饱和"/"进拟合"）、
  阈值校验、δ CI 包含点估计。
- **通用性**：清掉 `ramp.py` 报告串里残留的项目数值（把"泄漏倍数的具体数字"改成"M_off/M_base 的地板价"）。
- 文档同步：`protocol-ramp.md` §3.4 改四段账、命令示例补两个新参数；
  `protocol-execution.md` §7.2/§7.4 更新（尾项已建模）；`multitarget-code-map.md` 同步。

### Added（2026-11-15，代码审查产物）

- **`tests/test_ledger_formula_pin.py`：把"同一批公式的两套实现"钉在一起**。
  审查发现 `mtblocks.minimal_n0 / window_table_f_step`（文档积木）与
  `ramp.run_ledger / plan_ramp`（实跑）各写了一遍同样的账，公式一致但无测试锁住。
  7 项对拍：N₀ 同源、代数↔分裂数换算、乘客 = n·(μ_perbp×G)、失守 = 1−exp(−E·L/G)、
  梯级累计 = 基数 + i×每档、`n_post` 只在最后一档。
  **已用变异测试验证**：只改 mtblocks 或只改 ramp，测试立刻红（4 项 / 2 项失败）。
- **`docs/multitarget-code-map.md`：多靶点路径的代码地图**——九个模型的公式与落点、
  算法清单、八项数据处理、七条明确不建模的边界、维护者须知（改公式要改两处 + 对拍命令）、
  一次完整调用的数据流。
- `mtblocks.minimal_n0 / window_table_f_step` 与 `ramp.run_ledger` 的注释里加回指，
  写明"同一公式还有一份实跑实现"。

### Changed（2026-11-15）——通用性政策：技能不得携带项目专属药物名与参数

**用户定的两条硬规则**（细则 `docs/generic-policy.md`）：
① 药物名/诱导剂名只允许来自**被引用的文献**；② 默认值只允许来自**文献确证的通用值**，
本地项目推导的数不得作为默认值、示例值或"典型值"出现。

- **代码：项目推导的默认值改为必填**——`power --mu/--mu1/--nt-mean`、`fit-tauoff --td`、
  `design-induction --nt` 均为必填并给缺参提示；`--n-plate` 帮助文本去掉"1e8≈27 代"这类项目推导。
- **示例数据：换成文献数据集或明确标注的合成数据**——`dose_induction.csv` 改为合成的 Hill 曲线
  （文件头写生成参数）、`markers_demo.csv` 改为虚构演示值（只演示列格式）；示例菌株卡/记录模板去项目株代号。
- **文档：算例一律改为符号式**——`protocol-execution.md` 全篇改写为 `M_on/M_off/f/Nt/n_win/n_plate`
  的符号算例（并新增"把稀释倍数当档间代数"这一易错点的说明）；`protocol-ramp.md` §3.4/§3.6/§7、
  `protocol-fstep.md` §7、`experiment-design.md`、`multitarget-tutorial.md`、`start-here.md`、
  `glossary.md`、`references/intake.md`、`references/markers.md` 的项目数值与项目株代号一并移除/改占位符。
- **新增自动守卫**：`selftest.py` 加「通用性」检查（第 24 项），扫描
  `SKILL.md/CHANGELOG/docs/references/windowcalc/examples`，命中项目专属名字或数值指纹即失败并打印 `文件:行`。
- **SKILL.md**：把两条规则写进"设计目标"后的硬约束段，并在铁律里加"数值参数必填"。

### Added（2026-11-15）

- **`docs/protocol-execution.md`：浓度周期的执行层**——补上 `protocol-ramp`（阶梯怎么排）与
  `protocol-fstep`（f_step 怎么测）之间那段"试管与平板之间怎么倒"的空白，教学 + 选型两用：
  - **涂多少**：不是"全涂"，按 `N_plated = C_target/f_step` 反算；附 f_step → 涂板比例表
    （10⁻⁵ → 1%，10⁻⁶ → 10%）与"第一档三档夹逼"预实验法；放线菌 10⁷ CFU/板、团块打散、孢子口径三个上限。
  - **剩的去哪**：五个去处（传代主线 90–95% / 涂板测量 5–10% / 甘油 ×3 / 测序沉淀 / 备份活菌）+ 板量预算 6 块/档
    （利用"一块板两用"：c_{N+1} 板同时是本档 f_post 与下一档 f_pre）。
  - **传代拓扑二选一**：板传代（严格但 20–27 代、C 个谱系瓶颈、洗脱偏差）vs 液传代
    （无瓶颈、代数 `log2(N0/(f·Nt))`、但需防"选择不严"）；含"什么情况必须板传代"的四种情形。
  - **桥接验证**：`f_liq ≥ f_plate` 是严格不等式而非误差（板问"能否走完 27 代"、液问"能否启动生长"），
    给出 P/L/C 三臂设计与判读阈值（比值 ≥10 = 液体选择漏）。
  - **并行时序**：为什么不等板（4 条理由）、四条护栏（先冻菌 / 24–48 h 浊度闸 / 步长反推 / 板臂涂 t=0 同管）、
    并行的固有代价与三条折中；与 `protocol-fstep.md` §1.2 互补（该节给时间线，本文给护栏）。
  - **乘客账随拓扑变**：档间代数 `n_grow = log2(N₀/(f·Nt))`；给出板/液两套符号式，
    并**明确更正**本会话早期"液传代省 0.4"的错误估算（错因：按稀释倍数算代数，忽略了只有 f·Nt 个细胞在分裂）。
  - 一页操作卡（T0→T+8 时间线）、七条红线、待实测清单（τ_off / td / Nt / 每克隆细胞数 / f_liq-f_plate 比值 / 洗菌回收率）。
- **传代拓扑推荐"液传代"（板只做审计）**——理由：文献五篇渐升协议全是液体群体传代，本规程早期的
  "板上传代"是为绕开"异质性生长使 OD 不可靠"的局部偏离，代价是每档一次单谱系瓶颈 + 20–27 代结构下限。
  改法：判读不靠 OD，改靠"48 h CFU 增长倍数"。`docs/protocol-ramp.md` §0/§1.1 同步改写（消除与
  `protocol-fstep.md` §1.2 的相互矛盾），SKILL.md 路由表与 README 文档地图新增本文档入口。
- **洗菌与携带**：涂板前离心洗一次（去诱导剂），附残留浓度换算公式（C_ind × 涂板体积/板体积）；**洗菌回收率**作为与 f_step 同级的必录项（<50% 本档降级为方向性）。

### Added
- **新子命令 `ramp`：渐升梯级规划器——多靶点场景的"预测"落点**。11-13 路线收缩
  删掉的是"用不可测的量（μ×T_set）做预测"；本命令把预测重新锚在**可测的量**上：
  每档铺板时多涂几块浓度板，即得 PAP 尾部曲线（f(c) = f_ref × 10^(−(c−c_ref)/δ)，
  δ = 十倍频浓度宽），规划器据此输出——可行最大步长 Δ（含 margin）、推荐梯级表
  （每档预测 f_step ±CI、需铺 CFU、富余倍数）、**累计乘客三段账与失守概率**
  （总乘客单次代入，不做每档概率相加）。
  - 平稳补充假设（等步长下每档 f_step ≈ 首档实测）显式声明、每档被实测审计
    （回填 --pap 重算，预测偏出 CI 报警重拟合）；
  - 供给富余的档输出"可评估关诱导"提示（按需诱导门控——泄漏修好后才有意义）；
  - 外推超出 PAP 实测范围 >1 个 δ 时显式警告；
  - 零克隆行按 rule of three 作检测上界、不进拟合（拟合高于上界报警）；
  - 缺 --m-off 时窗口外段不计入并显式注明（本会话已证该项与窗口同量级）。
  18 项新测试（合成 δ 恢复、梯级规划、卡住/边界、渲染、CLI 缺参拒绝）。
  配套规程文档 `docs/protocol-ramp.md`（文献两大流派、预测边界精读、
  规划器数学推导、记录字段表——教学级详述）；SKILL.md 速查表同步。
  文献锚点：梯度孔自适应爬升（pgen.1011627）、PAP 尾部 log-线性（异质性耐药
  标准方法）、"每次突变抬升 log-MIC"的加性模型（同前）。

### Fixed
- **P1（数字会错）：`plan_ramp` 的 BGC 失守概率曾把保护区占比分母写成 1e6**，
  未用实际基因组大小（风险被按基因组大小放大）。已加 `genome_mb` 参数并
  按 bgc_kb/(genome_mb×1e6) 计算，测试锁定严格式。

## 0.9.0 — 2026-11-13

（定稿说明：本版合并原「未发布 2026-11-12」与「未发布 2026-11-13」两个工作区区块；
条目内保留原始日期。11-12 条目中因 11-13 删文件而失效的指针已就地加后注修正。）

### Removed
- **`multitarget --mode design/fit/classify` 整体移除（2026-11-13，用户拍板）**：
  三个模式服务的都是"用 μ 外推多靶点供给"的旧路线——供给已改为 `--f-step` 铺板实测，
  α 不在主路径任何环节被消费。连带删除：`windowcalc/multitarget.py`（421 行，
  含标度律/α 拟合/S-C-U 判定；C 律恢复公式已在 docs/multitarget-model.md §2 留档）、
  `--targets/--entry-fraction` 旗标、`mtblocks` 的 `window_table/target_set/first_step_events`
  （几何三件套）、`tests/test_multitarget.py`（21 项）与 test_mtblocks 的 8 项几何测试。
  `multitarget` 子命令保留，**等价于 window 主路径**
  （`--f-step --markers --k --nt --genome-mb [--bgc-kb]`；`--mode` 旗标删除）。
  S/C/U 概念保留在 `docs/multitarget-extrapolation.md`（规划步数预期用，不进窗口计算）。
  连带修复 selftest 文档检查的 `known` 子命令集合——v0.8.0 加第 9 个子命令时漏更
  （硬编码 8 个，任何含 "windowcalc multitarget" 的文档都会被误报未知子命令；
  现已补齐并加"与 build_parser 保持同步"注释）。
  路由层同步：SKILL.md（铁律/用户行/子命令表/description）、intake.md（路由表 + 必问改为
  f_step/标定/k·G·Nt 三样实测）、start-here.md（问题 2 重写）、guardrails.md 重写
  （f_step 三问 + M 误差棒与 indel 下限）、model/tutorial/extrapolation 顶部加路线收缩注。
- **C 律（谷型）整体移除（2026-11-12）**。理由：本项目全部是"逐步升浓度"的连续性状协议
  （如 1 g/L → 5 g/L），中间态在结构上必然可被选择，谷型不会发生；保留它只会让人
  误以为报告里的 `F^j` 是要认真读的输出。移除范围：
  `multitarget.t_valley`、`design` 的"C 谷型"行、`F^j` 最坏情况与相关字段、
  `mtblocks` 的多命中核（`events_multi_hit` / `poisson_sf`）、CLI 的 `--js`。
  如需恢复：`T_j = s^(j−2)/(2·N·μ^j)`（Iwasa/Michor/Nowak 2004；Armitage & Doll 1954）；
  公式与推导留档于 `docs/multitarget-model.md` §2（原在 `multitarget.py` 模块头，
  该文件已随上条移除）。

### Fixed
- **P0（数字会错）：`p_success` 从 P(≥1) 改为真正的 P(≥k)**。此前 `evaluate`/`monte_carlo`
  算的是 `1−exp(−events)`（至少 1 个事件的概率），报告与 MC 稳健性门槛却把它标成
  "P(≥k)"并拿去过滤——大 events 下恒≈100%，`optimize` 的 MC 门槛形同虚设。
  现在 `evaluate(..., k=)` / `monte_carlo(..., k=)` 接收 k（optimize 传 `k_target_clones`），
  `p_success` = Poisson 上尾 `poisson_sf(events, k)`——该实现从 optimize.py 提升为
  model.py 的共享函数（标量/ndarray 皆可；k=1 时逐位退化为旧值，初筛 `_p_succ` 语义不变）。
  meta 新增 `mc_gate_all_dropped`；报告 §4 在门槛拦下全部候选时给出 ⚠️ 提示。
  回归测试 ×3（上尾语义 / evaluate 随 k 单调 / MC 随 k 单调）；全量 **179 通过**。

### Added
- **新子命令 `fstep`**：f_step 铺板记录（step,when,cfu_plated,clones）→ 背景扣除后的
  可粘贴 `--f-step` 值。单比例 95%CI 用 Wilson（Brown/Cai/DasGupta 2001 推荐；
  k=0 用 rule of three 上界 3/n，Hanley & Lippman-Hand 1983）；两比例之差用
  Newcombe (1998) 法 10；f_pre ≥ 0.5×f_post 打背景淹没红灯（protocol-fstep.md §3）。
  配套：`examples/fstep_template.csv`、回归测试 ×6、selftest 第 23 项、SKILL.md 速查表。
- **`docs/protocol-fstep.md`**：f_step 铺板规程（配对设计、3/f 板量倒推、Wilson CI 数例、
  背景扣除、与 μ_supply 的口径关系、七个坑）——全部数字标 [文献]/[推导]/[约定] 三级。
  §1.1 新增**开跑前基线体检（微型 PAP）**：f_pre 沿浓度档 = 群体 MIC 分布上尾的实测，
  专抓异质性耐药/持续细胞（f_pre 不随浓度下降 → 停，先做单克隆复测）。
- **glossary §七**（f_step 时代 9 词条）+ **rate 口径提醒**（频率数据别喂 MSS-MLE，
  stderr 固定提示）+ SKILL.md 常见坑第 11 条。
- **文档路由收缩注**：model/tutorial/extrapolation 顶部注明 design/fit/classify 已删除、
  现行主路径为 `multitarget --f-step`（配套 Removed 条目）。
- **provenance 逐条标出处（2026-11-12）**：原 `multitarget.py` 模块头逐条标出处与状态；
  `docs/multitarget-model.md` 新增 §8（经典 Haldane/Iwasa/Armitage–Doll vs 作者综合
  vs 启发式阈值）+ §8.3 本地可核文献。〔11-13 后注：`multitarget.py` 已随 Removed 条目
  整体删除，出处记录现由 `windowcalc/mtblocks.py` 模块头与 docs 承接。〕
- **`windowcalc/mtblocks.py`（2026-11-12）**：多靶点窗口积木层（过原点回归、入口折扣、
  一阶供给）。发布时配套 `multitarget --mode window`；〔11-13 后注：`--mode` 已删除，
  现为 `multitarget --f-step` 主路径。〕

### Changed
- **文档出处补链（仅文档，无代码改动）**：`docs/multitarget-model.md` 新增 §8.4
  「位点→基因组换算的出处（target-size 法）」——Drake 1998（*Genetics* 148:1667）、
  Lang & Murray 2008（*Genetics* 178:67，多标记 ÷T_eff 原型）、Lee 2012
  （*PNAS* 109:E2774，WGS 真值）、Takemoto 2018（`10.1093/nar/gky481`，
  `T_eff(Rif)≈7 bp` 的同株配对实测出处）、Jee 2016 / Dillon 2016（均匀性 ≥10× 边界）；
  附外部文献对照（Tanneur 2025 的 MA-WGS 实测，与换算结果同量级）
  与适用范围提醒（M 只进乘客支路；rRNA 靶标不可外推）。
  项目内多靶点外推可行性论证文档的引用表补同行四条，
  §4 的 `T_eff(Rif)≈7 bp` 由"项目口径，需核对"升级为"gky481 同株内部对实测"。
- **§8.4 追加核对（读 PMC2206112 全文后）**：Lang & Murray 2008 的贡献是
  「T_eff 实测（207+227 个抗性克隆突变谱 + bootstrap CI）+ 每标记各自除」，它**没有**
  多标记联合回归——两 loci 每碱基率差 1.7×（URA3 3.80 vs CAN1 6.44 ×10⁻¹⁰/bp/代）
  且显著，恰是 `mu_perbp` 的 `spread` 检查要抓的对象。合并成单斜率是本工具在
  "均匀性一阶成立"下的扩展；`spread>2` 报警线与 L&M 自家数据（1.69）相容。
  已知差距两条：① 本工具 T_eff 靠文献/WGS 自校准，代码未实现"从突变谱测 T_eff"；
  ② 回归未用 μ_i 的 CI 加权、M 无 CI 传播；速率模块无涂板后生长两参数修正
  （L&M 实测 1× canavanine/α-factor 因此违反 L-D，靠 10× 浓度规避——本工具靠协议规避）。
- **L3 对拍环境常驻化（自建 R 参考环境）**：确认 `fluctuation` env（R 4.4.3 + rsalvador 1.9 + flan，
  由内部部署脚本部署）并实跑复验三组验收数据——Demerec1945
  **10.84383**（逐位一致）、LD1943_exp16 `newton.LD.plating(e=0.4)` **1.18636**（逐位一致）、
  plating_z02 **3.173983** vs 本文 3.174（六位一致）。`testing.md` §4/§8 已记录，
  含"部分涂板必须用 plating 估计器（整板估计器给 0.694）"的现场教训、
  悬空指针更正（数值记录实际在 `dev/source_code_audit.md` §4），
  以及 rsalvador 无多标记合并功能、其 foldchange/LRT 系列可供 `compare` 对拍的备注。
- **垃圾清理（用户确认后执行，2026-11-13）**：删除 `报告_示例.md`（2026-10-26 旧示例，
  通篇展示未修复的 P0 "P(≥k) 100%" 标签 bug——修复后可重新生成）、
  `SKILL.md.bak-20261030`、孤儿数据 `data/example_M7on_rif.csv`（全项目无引用）；
  清理 `__pycache__` 与 `.pytest_cache`（项目源 + 技能副本，运行时产物，自动再生）。
  `data/` 空目录保留（README 约定的用户 CSV 位置）。
  【同日已执行】`.bak-20261111/` 移至工作区 `archive/`（保留回滚能力，未删除）；
  13 个维护者文件移入根目录 `dev/`——**有意偏离"docs/dev"提议**：`docs/` 会被整目录同步进
  技能副本，`l4_answer_key.md`（L4 盲测答案）绝不能进技能副本，故统一放不同步的 `dev/`；
  残留引用（README 3 个链接、testing.md、maintainer.md、引擎 docstring）已全部改指 `dev/`，
  selftest 断链检查通过。
  【再更新，同日】`dev/` 已按用户指示**整体删除**；`l4_answer_key.md` 移至工作区
  `archive/l4_answer_key.md`（L4 盲测答案，保留）；README/testing.md 的引用已改为
  指向 CHANGELOG 与 archive。
  **另**：技能副本仍是 11-12 午间旧快照（缺 `mtblocks.py` 与 `--f-step`），
  需用户在可写 shell 跑 `bash sync_skill.sh && python selftest.py`。
  【2026-11-13 晚更新】同步已完成并验证（7 项逐字节一致；技能副本 selftest **22/22**）。
- **selftest 增加 [SKIP] 语义**：`run.py`/`bootstrap.py` 按设计只在部署副本
  （sync_skill.sh「不同步」清单），因此在**项目源**目录跑 selftest 时「入口」一项显示
  [SKIP]、规范项的 name=目录名 检查自动跳过并打印 note——退出码只由失败项决定。
  实测：项目源 21 通过 + 1 跳过；部署副本 22/22。此前在项目源跑出的 2 个 FAIL
  均为此设计差异，非回归。
- **`--mode window` 主路径改为 `--f-step`（实测，2026-11-12）**：`f_step`（克隆数/CFU）
  一次性吸收「靶点个数 × 靶点长度 × 有益比例」，因此**不需要**知道多靶点几何。
  `--targets` 降级为「仅已知靶点的蛋白药（Rif 式）」的可选路线，`--help` 已标注
  不适用于膜/耐受体系。理由：对"耐受"这种连续性状，靶点个数与长度不可枚举，
  `T_eff`/`T_set` 属于假精度。〔11-13 后注：`--targets`/`--mode` 均已整体删除，见 Removed。〕
- **`design` 输出删除"预期代次"表（2026-11-12）**，改为「绝对代次不预测」：`S≫1` 时
  稀有突变前提失效，该表给出的数字会被误读。`t_wt`/`t_mut` 仍保留在 JSON 里，
  只是不再渲染。〔11-13 后注：该输出随 `--mode design` 移除而不复存在。〕
- `mtblocks`（2026-11-12）：`accessible_fraction` → **`entry_fraction`**（明确为
  "靶点里能当第一步的比例"，膜靶药文献例（ncomms13002）= 2/5 ≈ 0.4）；窗口求解改成解析式
  `N0 ≤ Nt − k/(μ·activity)`，不再需要数值搜索。
- `SKILL.md` / `docs/multitarget-model.md` / `tutorial` / `glossary` / `extrapolation` /
  `references/multitarget-guardrails.md` 同步（2026-11-12）：**加速只报线性 ≈ F**；
  `F^j` 仅保留为"不要当预测"的警示。

## 0.8.0 — 2026-11-11

新增**多靶点外推**（第 9 个子命令 `multitarget`），并收紧"缺参填默认"的行为。

- **新子命令** `multitarget --mode design|fit|classify`：
  - `design`：算 `N·μ` 饱和 → 三种标度律（单步/串联可选/谷型）预期代次 → 可分辨的 α 区间 → 96 孔排布表。
    饱和时**不给 α 预测值**（模型前提失效），改为要求直接测 `T_L`。
  - `fit`：由三档 `(F, T)` 走 log-log 斜率，给 `α` 与 CI；2 档时明确标"CI 不可估"。
  - `classify`：按抗性证据定 S/C/U（可判 P=表型终止），**缺项降级为 U 并给下一步**，不猜。
- **默认结论改为线性**：`design` 的输出以「加速 ≈ F 倍」开头，`F^j` 降级为「最坏情况」，并标注"不是预测"。依据：独立多 locus 模型校验（`papers/多靶点抗性-理论/alpha-validation/`）实测 α 在 0.85–1.6，从未接近 j。
- **关键概念修正**：加速比不是 `F^j`，而是 `F^α`；`F^j` 只在"谷型（中间态不可选）"成立。
  串联但中间态可选时，等待时间是若干 `1/μ` 的**和**，加速仍是 `F`。
  并联靶点数 `k` 在 mutator/WT 的**比值**里约掉，只影响绝对时间。
- **行为收紧（用户明令）**：缺必填参数时**不得静默填默认**。`_require` 的报错升级为"决策级"——
  讲清默认值/后果/出路，并写明"agent 不得自行添加 `--fill-defaults`"。`multitarget` 用 `--defaults-ok` 显式同意。
- 新增文档：`docs/multitarget-model.md`（数学）、`docs/multitarget-extrapolation.md`（S/C/U 分类轴）、
  `docs/resistance-evidence-tiers.md`（证据四级链）、`references/multitarget-intake.md`（给 agent 的 L2 路由）。
- 自检 21 → **22 项**（新增"缺参 exit 2 + agent 禁令"检查）；pytest 135 → **156 项**。

## 0.7.5 — 2026-10-31

- 新增 `docs/viable-count-problems.md`：活菌计数不准的六种情况（随机误差、低回收、WT/突变体回收不同、团块/链、VBNC、超范围）及诊断与校正/替代计数（显微、超声曲线、spike-in、MPN、qPCR、流式）。
- 明确 ζ（涂板比例）与 e（回收率）的 Nt 口径区别，避免 `--fraction-plated` 重复校正。
- `SKILL.md` 与 `docs/knowledge-base.md` 加入口。

## 0.7.4 — 2026-10-30

- 实现 p0 部分涂板的精确稀释模型：`P0 = G(1−ζ)`（L-D PGF，bisection 解 m）；LD16（ζ=0.4）得 0.9787，与 rSalvador `p0.LD.plating=0.97868` 一致。
- 重跑并更新：ζ 手工 r/ζ 偏差（m=1/C=60/8 reps → 21%/39%/72%；winsor 13%/30%/68%；稀释似然 1.01–1.05）；`power --sim` 功效（C=40 68–72%、60 87%、75 93%、90 96.5%）。
- 更新 `docs/software_audit_2026-10-30.md` 与 `docs/knowledge-base.md`。

## 0.7.3 — 2026-10-30

文献核实与软件行为修正：

- **Axioms 2023 的 6%/98%**：原文已核（Devin & Couce 2023, Axioms 12:1100；引 Couce & Blázquez 2011），数字逐字一致。
- **FALCOR CI**：原文确认按 Stewart (1994) 的 ln(m) 正态近似（Foster 2006 近似）。
- **rSalvador CI**：原文确认所有 CI 均为 likelihood-ratio（profile）法，`confint.LD` 同法。
- **LD1943 / Demerec 验收数据**：本工具输出与 rSalvador 原文示例逐位一致（Demerec 10.8438 [8.6505,13.1948]；LD16 1.1864 [0.5803,2.0908]）。
- **p0 交叉核验**：部分涂板（ζ<1）下不再判为有效（p0 未按稀释模型校正），避免误导；以 MSS-MLE 为准。
- **m<0.3 提示**：新增 stderr 警告（R&F 不建议、精度差）——原来静默输出。
- 新增 `docs/software_audit_2026-10-30.md`（CLI 行为逐条核对）。

## 0.7.2 — 2026-10-30

- 去质粒化：删掉文档与示例里的具体质粒名；拷贝数 `N` 与限制温度一律由用户提供，治愈代数只保留公式 `ceil(log2 N)`。
- `T_eff≈7 bp` 标注为“文献口径，需自行核对”。

## 0.7.1 — 2026-10-30

文献核实后的更正（无 API 变更）：

- **四环素靶标**：由误写的 23S 更正为 **30S/16S**；抗性以efflux/核糖体保护为主，rRNA 点突变少见。
- **rrn 拷贝数**：由"放线菌常有 5–7"改为按属区分（链霉菌 5–7，S. coelicolor 6；分枝杆菌 1–2）。
- **m>15 估计器**：不再写"MSS 唯一可用"，改为"中位数/p0 不适用；用最大似然类（MSS/GF/NR-MLE）"。
- **"均值发散"**改为"重尾、样本均值被 jackpot 主导"（有限 Nt 下均值有限）。
- **m≥0.3** 明确为 R&F 的建议（除非大量培养物），不再写"任何估计器都不可靠"。
- **rSalvador 引用页码**：7:3661 → **7:3849–3856**。
- **RF2000 验收数据**：旧文件是 52 管截断版且把中位数法 1.11 当成 MSS；已替换为完整 60 管 Table 3，MSS=1.2274（原文 1.23）。
- 同步更新 `docs/testing.md`、`examples/README.md`。

## 0.7.0 — 2026-10-30

### Breaking
- 菌株卡字段中性化：`M_WT` → **`M_base`**，`mu_St` → **`mu_supply`**（旧名不再兼容，旧卡需手动改键）。
- `fit-fc` 的 `--mu-wt` → `--mu-base`（旧名删除）。

### Added
- 两条独立保护上限：`P_mutator_max`（默认 0.05，`role: mutator`）与 `P_protected_max`（默认 0.10，`role: bgc/other`）；冲突时按上限截断窗口并报告实际 events、缺口与 binding 来源。
- `protected[].role`：`mutator` / `bgc` / `other`。
- `optimize --protocol` 可省略（用内置默认协议卡）；`--fill-defaults` 用默认假设值补菌株卡缺项（报告标 ★未测）。
- 外推元数据：`mu_supply_source` / `selection_marker`；两者不同时自动提示按 `references/markers.md` 处理。
- 新手护栏：`references/intake.md`（提问清单）、`references/guardrails.md`（盲区排查表）、`references/markers.md`（标记外推）。
- `rate` 结束提示整管/多板；`compare` 结束提示 p 值解释。

### Changed
- 文档去项目化：示例菌株/协议改为通用占位；文献菌株只作"文献示例"标注。
- `SKILL.md` 精简为路由 + 铁律 + 速查（细节移入 `references/`）。
- 新增 `README.md`。

### Migration
旧卡把 `M_WT:` 改成 `M_base:`、`mu_St:` 改成 `mu_supply:`；命令 `--mu-wt` 改成 `--mu-base`。
