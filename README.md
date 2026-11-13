# window-calculator —— mutator 表达窗口计算器

> 建立：2026-10-26。定位：为 mutator 表达窗口的瞬变设计提供**相位分辨的突变预算优化**——
> 在"筛到足够多抗性克隆"与"每克隆乘客负荷 ≤ 预算"之间，算出最优的表达窗口
> （脉冲时长 × 诱导剂浓度 × 接种密度/涂板量）。

## 0. 设计决策记录（2026-10-26，与用户确认）

| # | 决策 | 选择 |
|---|---|---|
| 1 | 决策变量（旋钮） | **三旋钮**：脉冲时长 T_pulse × 诱导剂浓度 c × 接种密度 N₀/涂板量 |
| 2 | 突变率数据 | **内置 MSS-MLE**（吃波动分析原始计数，须与 FALCOR/rSalvador 对数验证） |
| 3 | 软件形态 | **CLI + Markdown 报告**，代码与文档同目录 |

## 1. 软件定位：不是计时器，是相位分辨的预算优化器

设计规程里的静态公式 `E[m] = n × M_on` 只算脉冲项。实际流程中乘客负荷**逐相位累积**：

| 相位 | 诱导状态 | 负荷贡献 | 静态公式 |
|---|---|---|---|
| 0 接合后恢复/预培养 | off（泄漏） | n₀ × M_off | ✗ |
| 1 脉冲（含上升沿 τ_on） | on | n_pulse × M_on | ✓（部分） |
| 2 洗脱衰变尾 | on→off（τ_off） | 蛋白不降解就持续突变 | ✗ |
| 3 选择板菌落长出 | off（泄漏） | ~10–20 代 × M_off | ✗ |
| 4 治愈/孢子传代 | off（泄漏） | 2–4 代 × M_off | ✗ |

全强度 M7（F≈1000×）时脉冲项绝对主导，静态公式够用；**弱 mutator（10–100×）＋多筛克隆**时
off 相位贡献与脉冲同量级——这正是原设计规程自己推荐（"反直觉建议"）的策略，静态公式在此系统性低估负荷。

## 2. 模型

### 2.1 三个方程

```
开关活性 a(t)：0→1 上升（τ_on），1→e^(−t/τ_off) 衰减（τ_off）
  M_on(c) = M_on_ref × [F(c)/F(ref_dose)]         锚定直接测量的 M_on_ref；剂量只给形状比
  M(t) = M_off + (M_on(c) − M_off)·a(t)           每基因组每代乘客突变数

供给（独立事件数）   events = Σ_i μ_supply · [F(c_i)/F(ref)] · ā_i · ΔN_i
负荷（每克隆期望乘客）E[P]  = Σ_i M(ā_i) · n_i        n_i = 该相位每世系代数
成功率              P = P(≥k)（Poisson 上尾）          （独立事件，非抗性菌落数——jackpot 是假供给；v0.9.0 起为真 P(≥k)，k=目标克隆数）
```

参数化要点（v0.1 教训，已进测试）：诱导态速率锚定 **M_on_ref**（M7-on 波动分析的直接
测量值），M_base 只用于泄漏红线（≤3×M_base）与 F(c) 缺省回退；μ_supply 与 M_on_ref 的不确定性
各自独立进入，不得通过 F_ref 交叉污染（tornado 会暴露这种耦合）。

### 2.2 两条标度律（软件最值钱的输出）

1. **供给 ∝ 窗口内总分裂次数**（≈ Nt − N₀，即最终涂板 CFU），与脉冲长短基本无关——
   长到饱和，事件数就定了；
2. **负荷 ∝ 每世系在窗口内经历的代数** n = log₂(Nt/N₀)——**高密度接种是最强降荷旋钮**。

推论：把 F 拉满不增加"有效"供给（超出筛选通量的事件数是浪费的负荷），只增加乘客。
优化器据此在 (T_pulse, c, N₀) 空间里找 min E[P] s.t. P(events≥k) ≥ 0.95。

### 2.3 不确定性传播

- 波动分析的 95% CI → 参数按对数正态采样（σ = ln(CI_high/点估计)/1.96）；
- 克隆间离散（gkab583：Rif 频率差 15×）→ 额外乘性 lognormal（σ_ln 默认 1.0，strain card 可调）；
- 输出 E[P] 与 events 的 P5/P50/P95，不给假点值。

## 3. 需要提前测定的参数

### A. 已有规程覆盖（波动分析规程 §2/§3.1 已排期）

| # | 参数 | 测法 | 软件用途 |
|---|---|---|---|
| 1 | M_base | 波动分析，Rif 板 | 负荷基线 |
| 2 | M_on @0.3 mM | 同上 | 脉冲项 |
| 3 | M_off（泄漏，合格 ≤2–3× WT） | 同上 | 相位 0/3/4 |
| 4 | td（±纤维二糖、±Apra） | 生长曲线 | 代数换算；表达负担 |

### B. 必须新增（现有文档缺口）

| # | 参数 | 测法 | 为什么关键 |
|---|---|---|---|
| 5 | **τ_off 洗脱衰变尾** | 加诱导 12 h → 洗脱 → 0/2/4/8/16 h 小规模波动分析，拟合 e^(−t/τ) | ★ 最可能推翻设计：关启动子≠关 mutator，DnaE1-M7 (~130 kDa) 的功能半衰期决定有效窗口。报告基因只能作辅证，**以 Rif 速率衰减为准** |
| 6 | **μ_supply @选定浓度** | 波动分析选择板换成 St（5–10× MIC），WT/M7-on/off 同批 | 供给端主参数。不可按每碱基折算（低浓度富集 rsmG 类、高浓度 K88 类；位点速率差 420×）。顺带测背景假抗性率 b |
| 7 | **F(c) 剂量响应** | ≥3 个纤维二糖浓度各一轮波动分析 → Hill 拟合 | 没有它只有 T_pulse 单旋钮；弱 mutator 策略下降 F 恰是主要优化空间 |

### C. 选配

- 首脉冲后 WGS 3–5 克隆直接数乘客（先导实验）→ 软件的**验收实验**：预测→对账→校准；
- 治愈动力学复核（30 °C 无抗 60%/两代是否在本菌成立）。

### D. 不测不管

- per-base 理论折算（T_eff≈7 bp 陷阱，规程已写明）——软件直接吃实测 μ_supply；
- 产量命中率（下游筛选漏斗，软件只回答"需要多少抗性克隆"）。

## 4. 软件架构

```
window-calculator/
├── windowcalc/
│   ├── mss_mle.py     # MSS-MLE：MSS 递推 + Brent 最大化 + profile-likelihood CI + p₀ 交叉核验
│   ├── model.py       # 相位积分引擎：events / E[P] / P_succ
│   ├── optimize.py    # 三旋钮网格+局部精修；蒙特卡洛（对数正态参数 + 克隆间离散）
│   ├── report.py      # Markdown 报告：时刻表 / 乘客瀑布 / tornado / 红灯
│   └── cli.py         # rate / optimize 两个子命令
├── examples/          # strain_card / protocol_card 样例
├── tests/             # R&F 2000 Exp.2 验证 + 性质测试
└── data/              # 波动分析原始计数 CSV（r,Nt 两列）
```

**用法**：

```bash
cd window-calculator
# ① 波动分析速率（内置 MSS-MLE；--nt-mode per_tube 为 Ycart 逐管校正）
python -m windowcalc rate --input data/批次.csv --label M7_on --nts 3e8
python -m windowcalc rate --input data/批次.csv --nt-mode per_tube        # CSV 需带逐管 Nt 列

# ② 两批速率比较（LRT + fold-change CI）
python -m windowcalc compare --d1 data/wt.csv --d2 data/m7on.csv --labels WT,M7-on

# ③ 样本量/功效分析（precision=估速精度；power=两批比较检出功效）
python -m windowcalc power --mode precision --mu 2e-9 --nt-mean 3e8 --nt-cv 0.5 --sim
python -m windowcalc power --mode power --mu1 2e-9 --fold 2 --target-power 0.9 --nt-cv 0.5

# ④ 表达窗口三旋钮优化（strain/protocol 卡片见 examples/）
python -m windowcalc optimize --strain examples/strain_card.example.yaml \
    --protocol examples/protocol_card.example.yaml --out 报告.md
```

运行环境：任何带 numpy + pyyaml 的 python3（无 R、无网络依赖）。
自建 R 参考环境的 rSalvador/flan 仅作独立复核，非运行依赖。

### 4.1 MSS-MLE 实现要点（内置，与 FALCOR 的关系）

- 分布递推（Ma–Sandri–Sarkar 1992；Rosche & Foster 2000 Eq.16）：
  `p₀ = e^(−m)`；`p_r = (m/r) · Σ_{i=0}^{r−1} p_i/(r−i+1)`。
  已用复合泊松一阶原理核对 p₁、p₂（严格相等），并对 Table 3 全样本复算（见 tests）。
- μ = m/Ñt（Ñt = 各管 Nt 算术平均，FALCOR 原文口径）。
- 95% CI 用 **profile likelihood**（ΔLL = χ²₁/2 = 1.921），比 FALCOR 的 Stewart/Foster 正态近似
  更自洽；验收时与 FALCOR 网页对数（浏览器可用，curl 403）。
- 部分涂板（ζ=涂板比例<1）：**二项稀释似然** P_obs(r)=Σ_j p_j·C(j,r)ζ^r(1−ζ)^(j−r)，
  与逐管模式共用同一套多项式系数（2026-10-27 复核修复：旧版把 p_obs 的次数误当 r 截断，
  μ̂ 偏 −61%；详见 CHANGELOG 对应条目）。
  ✅ 与参照实现**同一模型、数值等价**：rSalvador `newton.LD.plating(e=ζ)` 与
  flan `mutestim(plateff=ζ)` 给出同一 pmf（差 ≤1e-16）；同一份数据 m̂ 相对差 2e-7
  （求解器容差）。旧版注记"rSalvador 是 η 微菌落模型、设计不同"**已作废**（见
  CHANGELOG 对应条目（原审计文件 2026-11-13 已清理））。
  ⚠️ **不要手工用 r/ζ 校正**：稀释抬高 p₀，除以 ζ 补不回来——实测（真值 1.000、C=60、10 次重复）
  ζ=0.5 → 0.78、ζ=0.3 → 0.59、ζ=0.1 → 0.29；二项稀释似然 0.99/1.04/0.99。
- 交叉核验：p₀ 法 m = −ln(p₀)（仅 m≈0.3–2.3 有效），报告并列输出。
- m 适用域：**m ≥ 0.3**（R&F 2000：低于 0.3 任何估计器都不可靠）；1.5–15 只是**中位数法**的适用域，MSS-MLE 全域有效（Foster 2006）；实现支持到 ~300（p₀=e^(−m) 下溢前），>300 报错。

### 4.2 引擎的诚实简化（v0.1）

- 菌丝体"每代"按波动分析同口径（per CFU per passage），n 由实测生长曲线换算——丝状菌口径
  局限与规程 §8 相同，比较只认相对量级；
- 表达负担以 td_on 与 td_off 之差进入脉冲相位代数（负反馈）；
- 衰变尾窗口内默认按 duration_h/td_off 折算代数（洗后即涂板请显式写 generations: 0）；
- 谱偏（PHP 缺陷产 indel，Rif 板测不到）不在模型内——先导 WGS 对账时校准 M_on。

## 5. 版本规划

| v0.2 | **rSalvador/flan 功能整合**：逐管 Nt 无偏估计（Ycart 2014，多项式系数法）、LRT 速率齐性、fold-change CI、模拟器——`compare` 子命令 + `rate --nt-mode per_tube`；与自建 R 参考 flan 交叉对数通过（见 CHANGELOG 对应条目） | **已实现**（29 测试全绿） |
| v0.2.1 | **样本量/功效分析**（`power` 子命令，rSalvador `samp.size.*` 对应实现）：解析 Fisher 信息版（常数 Nt，快）＋模拟实测版（一般，含 Nt 波动/部分涂板）；Fisher 信息经恒等式 ∂p_r/∂m=(p⊛g)_r−p_r 用自家 pmf 推导，与模拟 CV 渐近互验 | **已实现**（39 测试全绿） |
| v0.4 | **诱导剂剂量轴的机制模型**：表达（Hill）× 竞争占位（Langmuir）→ 速率；三种混合规则（linear / two_core / both_core）AIC 比较；`fit-induction` 子命令 + `strain_card.inducer` 块；报告新增「剂量轴饱和分析」。可辨识性：只报 Δsat = A·ψ(1)（linear 下 = ΔM·q），不报占比 p（q 与 ΔM 不可分，需报告基因或纯 M7 锚点才能拆） | **已实现** |
| v0.3.1 | **全量代码复核**（2026-10-27）：修复部分涂板两模式（数学错误 + 溢出）、单列 CSV 崩溃、无 Nt 时 YAML 输出假 μ、负/非整数计数静默错、Hill 退化数据崩溃、τ_off 死参数等 9 项；新增 16 个回归测试 | **已实现**（60 测试全绿；见 CHANGELOG 对应条目） |
| v0.3 | **拟合模块**（`fitting.py` + `fit-fc`/`fit-tauoff` 子命令）：F(c) Hill 剂量响应三参数拟合（fmax/EC50/n，bootstrap CI，剂量查询直出 strain_card.F_doses）；τ_off 洗脱衰变单参数拟合（μ_off 固定，CI + 尾部负荷折算）。合成数据参数回收验证；退化保护（点数不足/非单调/长 τ 警示） | **已实现**（44 测试全绿）——数据一到即可用 |
| **v0.6.0** | **专业复核 5 批修复**（2026-10-29）：P0 级 3 项会算错数字的缺陷（MC 漏 decay 活性、供给只算 pulse、"乘客/事件比=M/μ_supply"恒等式错误）、11 项统计口径（MC 两套区间分离、CI 缺失告警、tornado 剔除未用参数、结构参数 CI、样本量双重取整）、14 项接口（traceback→中文错误、边界解标注、候选池多样性、复现附录 A、`--strain` 读卡）；`model.phase_terms()` 成为负荷/供给的单一事实源 | **已实现**（124 测试全绿；见 CHANGELOG 对应条目） |
| **v0.6.3** | **`passenger_scale`（μ 单位点 → M 全基因组 的折算因子）**：只乘在绝对负荷上，窗口推荐对其尺度不变（有测试守着）；带 CI 时自动进 MC 与 tornado 第一行——把「绝对风险取决于这个未测因子到什么程度」摆到明面上；配套修正文档里「按基因组/位点长度折算」的错误建议（长度比 ~2e3 vs 速率比 1e5–1e6，按长度算会低估 ~500 倍） | **已实现**（124 测试全绿） |
| **v0.6.4** | **agent 零配置自举 + Agent Skills 规范合规（分发层）**：`run.py` 首次运行即检查依赖 → 技能目录内建 `.venv` → 默认源不通时自动换镜像 → 装好重跑原命令；`--doctor` 只诊断不安装；`install.sh` 与 `run.py` 共用同一套 bootstrap；Windows 加 `run.cmd`（自动探测 PATH / py / conda 默认位置，跳过 WindowsApps stub）；**frontmatter 按 [agentskills.io 规范](https://agentskills.io/specification)修正**（`allowed-tools` 改空格分隔字符串、`compatibility` 改为环境声明 ≤500 字符、无未知键），`agentskills validate` 通过；自检新增规范校验项。实测：**无 numpy 的新机器、隔离 HOME、默认 PyPI 不可达 → 40 s 自动装好并算出 m=1.1067** | **已实现**（124 测试全绿；自检 18→20 项） |
| **v0.6.5** | **Windows 实机验收 + 可分发**：`run.cmd` 在 conda 环境跑通（`--doctor` + `rate`，结果与 Linux 一致）；补 `LICENSE`（MIT 全文）与 `requirements.txt`；分发打包 `zip` + `tar.gz`（42 文件，排除 `__pycache__`/`.venv`/`.pyc`），**两个包解压后独立跑自检均 20/20** | **已实现**（2026-10-29 用户实机实测） |
| **v0.6.6** | **Excel 直接输入**：`--input` 接受 `.xlsx` / `.xlsm`（读第一个工作表，表头可留，**不需要另存为 CSV**）；`.xls`（97-2003）报错并提示另存。实现放在 `_read_text` 层（Excel → CSV 文本），下游解析逻辑**零改动**——表头跳过、列检查、行号提示全部复用，含逗号/引号/换行的单元格由 `csv.writer` 正确转义。`openpyxl` 加入运行时依赖（bootstrap 自动装）。新增 4 个单测（xlsx≡csv、中文表头+空行、`.xls` 拒绝提示、随包验收文件一致性）+ 1 项自检 | **已实现**（128 测试全绿；自检 20→21 项） |
| v1.0 | 先导脉冲 WGS 对账闭环（预测→实测→贝叶斯更新 M_on） | 待先导实验 |

## 6. 验收标准

1. `tests/` 全绿：R&F 2000 Table 3 复算 m̂ = 1.11 ± 0.03、CI 与 0.77–1.49 相容；模拟回收 ±15%；
2. 回代 gkab583：输入 20–27 代、全强度 M7，输出乘客 66–96 区间；
3. 先导实验前出预测，WGS 3–5 克隆实测对账（乘客数落 P5–P95）。

## 技能形态（分发，2026-10-29 更新为 agent 自举）

本软件已打包为 agent 技能 **`fluctuation-analysis`**，可导入 dsh / opencode / Claude Code / codex 等：
实验室成员说"帮我算这批波动分析/这两批速率有没有差/要铺多少管"，agent 加载技能后**自己把环境装好**
（首次在技能目录建 `.venv`、装 numpy/pyyaml；默认源不通时自动换镜像）再计算并解读——
**用户不执行任何安装命令**。

- **agent 入口**：`python3 run.py <子命令>`。`run.py` 内含自举逻辑（`bootstrap.py`），
  已装好依赖的解释器则零副作用直接跑；`python3 run.py --doctor` 只诊断不安装。
- **人类想敲短命令**：`bash install.sh` → `~/.local/bin/windowcalc`（与 `run.py` 共用同一套 bootstrap）。
- **引擎同步**：本目录为开发源，改完跑 `bash sync_skill.sh`（单向覆盖 `windowcalc/` 与 `examples/`，
  并自动 `diff` 校验），然后在技能目录跑 `python3 selftest.py`（19 项）。
- SKILL.md 内含领域口径（速率 vs 频率、丝状菌声明、两套样本量口径、红灯转达义务）——
  agent 的通用知识里没有这些，技能的存在意义就在于此。

**技能会回答"测什么数据、怎么设计实验"吗？**——会。技能自带自包含参考
`docs/experiment-design.md`（参数测定三梯队、波动实验设计速查、记录表格式、典型问答），
SKILL.md 要求 agent 在用户提出新测定需求时**主动给出完整设计**（管数、接种量、板浓度、
计数时机、对照组、模板），而非只给一条命令。记录表模板：`examples/record-template.csv`。
