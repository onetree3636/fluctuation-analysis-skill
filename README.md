# fluctuation-analysis —— 波动分析 + mutator 表达窗口计算器

把 Luria–Delbrück 波动实验的菌落计数变成带 95%CI 的突变率（MSS-MLE），并在此基础上做
**相位分辨的突变预算优化**：在"筛到足够多抗性克隆"与"每克隆乘客突变负荷 ≤ 预算"之间，
算出最优的 mutator 诱导表达窗口（脉冲时长 × 诱导剂浓度 × 接种密度/涂板量）。

纯 numpy + pyyaml，无 R、无网络依赖；CLI 输出 Markdown 报告。

## 子命令一览

| 命令 | 回答的问题 |
|---|---|
| `rate` | 这批波动实验的突变率 μ 是多少（MSS-MLE + profile-likelihood CI） |
| `compare` | 两批速率有没有差（LRT + fold-change CI） |
| `power` | 要铺多少管（估速精度 / 两批比较功效） |
| `fit-fc` | F(c) 剂量响应 Hill 拟合（fmax/EC50/n，bootstrap CI） |
| `fit-induction` | 诱导表达开关拟合（Hill×Langmuir 机制模型） |
| `fit-tauoff` | 洗脱衰变尾 τ_off 拟合（关启动子 ≠ 关 mutator） |
| `design-induction` | 瞬变设计：诱导时机与时长 |
| `optimize` | 表达窗口三旋钮优化：min 乘客负荷 s.t. P(≥k 克隆) ≥ 0.95 |
| `multitarget` | 多靶点/膜活性药物的窗口——供给只认铺板实测 `--f-step` |
| `fstep` | 铺板记录 → 背景扣除后的 f_step 值（Wilson CI / rule of three） |

## 快速开始

```bash
# 依赖：python3 + numpy + pyyaml（或 bash install.sh 装成 windowcalc 短命令）
python -m windowcalc rate --input data/批次.csv --label on --nts 3e8
python -m windowcalc rate --input data/批次.csv --nt-mode per_tube     # 逐管 Nt（Ycart 2014）
python -m windowcalc compare --d1 data/wt.csv --d2 data/on.csv --labels WT,on
python -m windowcalc optimize --strain examples/strain_card.example.yaml \
    --protocol examples/protocol_card.example.yaml --out 报告.md
```

输入接受 CSV/Excel（`.xlsx/.xlsm`）；波动分析原始计数为 `r,Nt` 两列。卡片样例见 `examples/`。

## 核心思想：乘客负荷逐相位累积，静态公式会低估

设计规程的静态公式 `E[m] = n × M_on` 只算脉冲项。实际流程中负荷在五个相位累积——
预培养/菌落长出/传代都是 off（泄漏）态，洗脱尾在蛋白降解前持续突变：

| 相位 | 诱导状态 | 负荷贡献 |
|---|---|---|
| 0 预培养 | off（泄漏） | n₀ × M_off |
| 1 脉冲（含 τ_on 上升沿） | on | n_pulse × M_on |
| 2 洗脱衰变尾 | on→off（τ_off） | 蛋白不降解就持续 |
| 3 选择板菌落长出 | off（泄漏） | ~10–20 代 × M_off |
| 4 治愈/传代 | off（泄漏） | 2–4 代 × M_off |

全强度 mutator（F≈1000×）时脉冲项主导，静态公式够用；**弱 mutator（10–100×）＋多筛克隆**
时 off 相位与脉冲同量级——这正是本工具存在的理由。

三个方程（细节见 `docs/mutator-window.md`）：

```
开关活性 a(t)：0→1 上升（τ_on），1→e^(−t/τ_off) 衰减（τ_off）
M_on(c) = M_on_ref × [F(c)/F(ref_dose)]    锚定直接测量的 M_on_ref，剂量只给形状比
M(t) = M_off + (M_on(c) − M_off)·a(t)      每基因组每代乘客突变数
events = Σ μ_supply·[F(c_i)/F(ref)]·ā_i·ΔN_i    E[P] = Σ M(ā_i)·n_i    P = P(≥k)
```

两条标度律（优化器最值钱的输出）：

1. **供给 ∝ 窗口内总分裂次数**（≈ Nt − N₀），与脉冲长短基本无关；
2. **负荷 ∝ 每世系在窗口内经历的代数** n = log₂(Nt/N₀)——**高密度接种是最强降荷旋钮**。

推论：把 F 拉满不增加有效供给，只增加乘客。优化器据此在 (T_pulse, c, N₀) 里找
min E[P] s.t. P(≥k) ≥ 0.95，蒙特卡洛传播 CI 与克隆间离散（gkab583：Rif 频率差 15×）。

## 用它之前要测什么

| # | 参数 | 测法 | 用途 |
|---|---|---|---|
| 1 | M_base | 波动分析，Rif 板 | 负荷基线 / 泄漏红线（≤3×） |
| 2 | M_on @工作浓度 | 同上 | 脉冲项锚点 |
| 3 | M_off（泄漏） | 同上 | 相位 0/3/4 |
| 4 | td（±诱导剂、±抗生素） | 生长曲线 | 代数换算、表达负担 |
| 5 | τ_off 洗脱尾 | 洗脱后 0/2/4/8/16 h 小规模波动 | ★ 最可能推翻设计 |
| 6 | μ_supply @选定选择板 | 波动分析（5–10× MIC 板） | 供给主参数；不可按每碱基折算 |
| 7 | F(c) 剂量响应 | ≥3 浓度 × 波动分析 → Hill | 弱 mutator 策略的主要优化空间 |

完整规程（三梯队、记录表、典型问答）：`docs/experiment-design.md`。

## 实现要点与已知边界

- **MSS-MLE**：Ma–Sandri–Sarkar 递推 + Brent 最大化 + profile-likelihood CI，p₀ 法交叉核验；
  适用域 m ≥ 0.3（R&F 2000），实现支持到 ~300。
- **部分涂板**：二项稀释似然 P_obs(r)=Σ_j p_j·C(j,r)ζ^r(1−ζ)^(j−r)，与 rSalvador
  `newton.LD.plating` / flan `mutestim` 同一模型、数值等价。⚠️ 不要手工用 r/ζ 校正——
  稀释抬高 p₀，除以 ζ 补不回来（实测 ζ=0.5 → 0.78）。
- **诚实简化**：丝状菌"每代"按 per CFU per passage 口径；表达负担以 td_on−td_off 进代数；
  谱偏（PHP 缺陷产 indel，Rif 板测不到）不在模型内——用先导 WGS 对账校准 M_on。
- **多靶点路线**：供给只认实测 f_step（一次性吸收靶点个数×长度×有益比例）；
  μ_perbp 只进乘客支路（`docs/mu-genome-derivation.md`），不用于供给外推。

## 文档地图

| 想做什么 | 看哪里 |
|---|---|
| 设计波动实验 / 测参数 | `docs/experiment-design.md` |
| 理解窗口模型 | `docs/mutator-window.md` |
| 多靶点/膜活性药 | `docs/multitarget-model.md`（模型）· `docs/multitarget-extrapolation.md`（分类轴）· `docs/protocol-fstep.md`（f_step 铺板规程） |
| 单位点 μ → 全基因组 M 换算 | `docs/mu-genome-derivation.md` |
| 验证、golden 对拍、已知差距 | `docs/testing.md` |
| 术语与口径 | `docs/glossary.md` |
| agent 用法与领域护栏 | `SKILL.md` · `references/`（intake / guardrails / markers / commands） |

## 作为 Agent 技能分发

本工具打包为 agent 技能 `fluctuation-analysis`（dsh / opencode / Claude Code / codex 可导入）：
实验室成员说"帮我算这批波动分析 / 这两批有没有差 / 要铺多少管"，agent 加载 `SKILL.md` 后
自己装好环境（`run.py` 首次运行自举 `.venv`，默认源不通自动换镜像）再计算并解读。
`python3 run.py --doctor` 只诊断不安装。

## 测试与验证

`pytest tests/` 全量 + `selftest.py` 自检；golden 锚点：Rosche & Foster 2000 Table 3 复算、
Demerec 1945 / LD 1943 exp16 / 部分涂板 ζ 系列（与 rSalvador、flan 逐位对拍，
见 `docs/testing.md` §4/§8）。版本历史与完整变更记录见 [CHANGELOG.md](CHANGELOG.md)；
路线图：v1.0 先导脉冲 WGS 对账闭环（预测→实测→贝叶斯更新 M_on）。

## License

[MIT](LICENSE)
