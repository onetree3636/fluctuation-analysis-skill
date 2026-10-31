# 子命令完整手册

> 由 `SKILL.md` 按需加载。示例统一写 `python3 run.py <子命令>`（Windows 用 `run.cmd`）；
> 装了用户级 `windowcalc` 时与之等价。**不要让用户 `pip install`，不要运行 `install.sh`。**

## 0. 通用输入

- 每行一管：`r[,Nt[,volume_mL]]`（r = 菌落数；Nt = **CFU/管**；可省略 Nt，但 μ 无法算）。
- `--input` 直接接受 Excel：`.xlsx` / `.xlsm`（读第一个工作表，不用另存为 CSV），
  或 `.csv` / `.txt`（UTF-8、UTF-8 BOM、GBK 都认）。`.xls`（97-2003）不支持，会提示另存为。
- 表头行可以留（中英文都会跳过）；列的位置必须对：第 1 列 r、第 2 列 Nt。
- 多一列“管号”、Nt 用了 CFU/mL 这两种最常见的错，报错/回显里都有针对性提示。
- 记录表模板：`examples/record-template.csv`；原始计数模板：`examples/raw_counts_template.csv`。

## 1. rate —— 原始计数 → 速率（最常用）

```bash
python3 run.py rate --input 计数.csv --label mutator_on --field M_on_ref --nts 3e8   # 片段键名=卡片字段
python3 run.py rate --input 计数.csv --label mutator_on --nts 3e8        # 经典口径
python3 run.py rate --input 计数.csv --nt-mode per_tube             # 逐管 Nt 校正（CSV 须带 Nt 列）
python3 run.py rate --input 计数.csv --fraction-plated 0.1          # 部分涂板（Nt 用整管总 CFU）
python3 run.py rate --input 计数.csv --nt-unit per_mL --volume-mL 2 # Nt 以 CFU/mL 给出
```

- `--label` 只是这组数据的名字；想直接得到能粘进菌株卡的片段，用
  `--field M_base|M_on_ref|M_off|mu_supply`（不给 `--field` 时片段键名 = label，粘进卡片会报“缺 M_base”，
  片段里会提醒改名）。
- 字段名是菌株卡的固定键：`mu_supply` 指“供给端/第二标记”的速率，**与具体药物无关**。
- 若实验室测的是 CFU/mL，加第三列 `volume_mL`（自动换算），或对 `--nts` 用
  `--nt-unit per_mL --volume-mL 2`。
- 输出先打印**输入回显**（C、r 合计、Nt 均值），再给 m̂±95%CI、μ、p₀ 交叉核验、YAML 片段。

⚠ **两条必须转达给用户的提示**：

1. 用 `--nts` 又没给 `--nt-unit` 时，输出会提示“按 CFU/管 解释”——CFU/mL 与 CFU/管 差一个体积倍，
   这是 μ 最常见的错（文档自己说占九成）；
2. YAML 片段里的 μ 是**单位点**速率，`unit: per_locus_per_passage`，**不能直接粘进 M_base**
   （模型的 M 是全基因组每代，差 10³–10⁶ 倍）；片段里已写警告，别省略。

## 2. compare —— 两批速率比较（LRT + fold-change）

```bash
python3 run.py compare --d1 WT.csv --d2 mutator_on.csv --labels WT,mutator-on
# 只有均值 Nt 时（不用逐管列）：
python3 run.py compare --d1 WT.csv --d2 mutator_on.csv --nt1 3e8 --nt2 3e8
```

有逐管 `r,Nt` 两列时用第一种。输出：两批 μ̂、LRT p 值（速率差异是否显著）、fold-change θ 与 95%CI。
`--nt1/--nt2` 是“常数 Nt”口径：方便，但各管 Nt 有波动时会系统性低估速率。

## 3. power —— 样本量/功效分析

```bash
python3 run.py power --mode precision --mu 2e-9 --nt-mean 3e8 --nt-cv 0.5 --sim  # 精度型：要多少管
python3 run.py power --mode power --mu1 2e-9 --fold 2 --target-power 0.9 --C 40  # 检出型：功效够不够
```

- 精度型两套口径：`--conf-level 0` = 1σ 相对标准误；默认 `0.95` = “95% CI 半宽/估计值”。
  同样写 20%，后者约为前者的 3.84 倍（m≈1 时 36 管 vs 138 管）。**不要混用**。
- 检出型是两批 LRT 比较的功效；C 不够时加大每批管数或提高 fold。

## 4. optimize —— mutator 表达窗口三旋钮优化

```bash
python3 run.py optimize --strain strain_card.yaml --out 报告.md            # 协议卡用内置默认
python3 run.py optimize --strain strain_card.yaml --protocol my.yaml --out 报告.md
python3 run.py optimize --strain minimal.yaml --fill-defaults --out 报告.md  # 缺项用默认假设值
```

样例卡片在技能 `examples/` 与项目目录。**这是专家功能**：strain_card 有 **10 项必填**
（M_base / M_on_ref / M_off / μ_supply / τ_on / τ_off / td_on / td_off / Nt_saturation / ref_dose）
与协议卡 5 项必填（k / P_min / T 区间 / N₀ 区间）；缺项默认**直接报错**（不会静默替代）。
只能先用假设值时显式写 `status: assumed`，报告会把它们标成“★未测”——如实转告用户。

**无感用法**：协议卡可省略（用内置默认）；菌株卡缺项可用 `--fill-defaults` 自动补默认假设值（stderr 会列出补了哪些，报告标 ★未测）。**不要让用户手写 YAML**——按 `references/intake.md` 用几个带默认值的问题问齐即可。

⚠ **μ_supply 是 events 端主参数**：没测就先做波动分析（12 管试点），不要只填默认——它和 `M_base` 的整体 scale 不同，会改变窗口 argmin；默认值下窗口只能当方向参考。
卡片可写 `mu_supply_source`（这个速率是用哪个标记测的）与 `selection_marker`（窗口用哪个标记筛）；两者不同时 `optimize` 会自动提示外推（见 `references/markers.md`）。

**保护优先**：协议卡有两条独立硬上限——`P_mutator_max`（默认 0.05，`role: mutator` 的保护区）与 `P_protected_max`（默认 0.10，`role: bgc/other`）。k 要求的窗口若超过任一条，optimize 按上限截断，报告 `保护上限 binding` + 实际 events + 缺口 + 哪条 binding；不要为了 k 超窗。连最短窗口都超限时报 `保护上限不可行`。

## 5. fit-fc —— 剂量响应

```bash
python3 run.py fit-fc --input 剂量.csv --mu-base 0.004 --query 0.5 12     # CSV: dose_mM,m_on
```

## 6. fit-tauoff —— 衰变

```bash
python3 run.py fit-tauoff --input 衰变.csv --mu-off 0.012 --td 2.5      # CSV: h,mu
```

## 7. fit-induction —— 诱导剂剂量轴的机制模型（v0.4）

把“诱导剂浓度 → 表达量 → 竞争占位 → 突变率”显式建模，并回答**饱和点在哪**：

```bash
python3 run.py fit-induction --input 剂量.csv --m-off 0.012 --query 0.05 0.5 5
# CSV: dose_mM,M（每个剂量一轮波动分析测得的 M，与 M_off 同单位）
```

- 三个混合规则（`linear` / `two_core` / `both_core`）用 AIC 比较，ΔAIC<2 时按简约取更简单的；
- 输出 **Δsat（有效占位幅度）、C50、n**，以及“95% 饱和剂量”和每个剂量的饱和百分比；
- ⚠️ 从突变率数据只能稳健定出 `Δsat = A·ψ(1)`（linear 下才等于 ΔM·q），**不能**据此宣称
  “突变型聚合酶占比 p”（q 与 ΔM 不可辨识；要拆开需报告基因或纯 mutator 锚点）；
- 拟合结果用 `strain_card` 的 `inducer:` 块接进 `optimize`：

```yaml
inducer:
  name: 纤维二糖      # 任意诱导剂
  mixing: linear      # linear 下 amplitude 才等于 Δsat；two_core/both_core 见 docs §2.1
  amplitude: 1.315    # = A（linear 时 = Δsat）
  c50_mM: 0.48        # 速率半饱和剂量（旧名 ec50_mM 也认）
  n_hill: 1.85
  ci:                 # 可选：把 fit-induction 的 bootstrap 区间粘进来，MC 才会传播剂量轴
    amplitude: [1.31, 1.32]
    c50: [0.45, 0.51]
    n: [1.60, 2.10]
```

## 8. design-induction —— 剂量轴方案（用户只说几个浓度）

用户说“帮我设计 N 个诱导剂浓度”时，**直接调用，不要让用户先给 C50**：

```bash
python3 run.py design-induction --strain 菌株卡.yaml --ref-m 6        # 卡里有 inducer.c50_mM 就直接读
python3 run.py design-induction --c50 0.48 --ref-m 6 --n-doses 5       # 只有预期 C50
python3 run.py design-induction --n-doses 4 --range 0.03 3             # 完全不知道 C50
python3 run.py design-induction --doses 0.1 0.3 1 --c50 0.48           # 用户直接报剂量
```

自动决策（不需要用户理解“试点/正式”）：

| 输入 | 行为 |
|---|---|
| 给了 `--doses` | 直接用；若同时知道 C50，会检查上下覆盖（不足会警告） |
| 没给剂量、有 `--c50`（或 `--strain` 读到的 `inducer.c50_mM`） | 在 C50 的 0.2–6 倍间对数布点，每点按预期 m 算管数 |
| 都没有 | 在 `--range` 内宽间距布点，每点固定管数（定位试点） |

**agent 义务**：用户只说“测 5 个浓度”时，先看 `strain_card` 里有没有 `inducer.c50_mM`——
有就带 `--strain 菌株卡.yaml`（不带 --strain 就带 `--c50 <值>`），没有就不带；
**绝不要把“C50 是多少”反问给用户**（那正是要帮他测的）。
⚠ `--ref-m` 只有在知道 C50 时才有意义：只给 `--ref-m` 不给 C50 会打印警告并忽略它。

## 9. 输出、边界与实测数字

- **部分涂板**只用 `--fraction-plated ζ`（二项稀释似然；Nt 填**整管总 CFU**），禁止手工 r/ζ——
  稀释会抬高零管比例，除以 ζ 补不回来：本机复跑 ζ=0.5/0.3/0.1 偏低 21%/39%/72%（winsor 时 13/30/68%）
  （二项稀释似然 0.99–1.04）。整管涂板是首选。
- **jackpot 管**（r 上千）先确认不是饱和/污染——引擎对 r>5000 会拒绝；
  部分涂板（`--fraction-plated`<1）叠加 jackpot（需要混合的 j=r/ζ >3000）也会明确报错，
  出路是整管涂板或 winsor 后重算。
- **r>150** 文献建议并为一类（Foster 2006 引 Asteris & Sarkar 1996，精度损失很小）；
  本引擎不强制（可精确算到 r=5000）。
- **MC 报告两个区间，别混用**：`只含参数不确定性`（P50 应≈点估计）与`再叠加克隆间离散`
  （`clone_dispersion_sigma_ln`，默认 1.0 → P5–P95 可差 ~20 倍）。后者是生物变异，
  **补测参数也不会收窄**——用前一栏判断“参数测准没”。
- **没有 CI 的参数**在 MC 里按 ±2× 对数正态假设（`inducer.*` 与 td/Nt_saturation 按常数），
  报告 §4 会点名列出来；要真区间就在卡片里写 `ci: [lo, hi]`（结构参数也支持）。
- **边界解**：最优解的剂量/T/N₀ 顶到网格端点时报告会标 ⚠ **边界解**——那不是内点最优，
  先放宽 `T_pulse_range_h`/`doses_available`/`N0_max` 再算一次。
- **复现**：报告附录 A 记录了引擎版本、剂量/T/N₀ 网格、MC 抽样数与 seed；同卡片同参数重跑
  可得到同一份报告。
- **克隆封顶 5000**：对外模拟器（`power --sim` / `compare`）与多项式系数法同域；内部性质测试 `simulate` 不封顶（仅测试用）。
