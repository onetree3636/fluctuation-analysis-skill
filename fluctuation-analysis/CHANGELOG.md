# Changelog

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
- 文档去项目化：示例菌株/协议改为通用占位；M7/gkab583 等只作"文献示例"标注。
- `SKILL.md` 精简为路由 + 铁律 + 速查（细节移入 `references/`）。
- 新增 `README.md`。

### Migration
旧卡把 `M_WT:` 改成 `M_base:`、`mu_St:` 改成 `mu_supply:`；命令 `--mu-wt` 改成 `--mu-base`。
