# 技能测试方案（发布前必读）

> **级别：附录（维护者）** ｜ 什么时候读我：发版前跑自检/对拍，或查 2026-10-30 的历史核对 ｜ 不用读我：你只是使用者
> 入口与完整分级见 `INDEX.md`。

> 技能 = 文档 + 引擎 + 打包。**只测引擎不够**：文档会写错数字，打包会在别人的机器上崩，
> agent 会把两篇文档混着答。所以分五层，每层都有可执行的检查。

## 0. 一条命令

```bash
cd 技能目录
python selftest.py            # 28 项，约 40–100 s
python selftest.py --quick    # 只查环境与命令（约 2 s）
```

自检覆盖：环境 → **入口自举（`run.py --doctor`）** → **Agent Skills 规范** → 9 个子命令 → **Excel 与 CSV 结果一致** → 4 组**已知答案**（R&F 2000 文献值、部分涂板、LRT/θ、样本量）
→ 端到端（optimize 报告章节、fit 参数回收）→ 错误路径（坏数据/YAML 都不甩 traceback）→
文档（断链、命令是否存在、四个用途入口、**文档里的数字与代码是否一致**）→
**护栏**（§0 正向：每个子命令都有护栏；反向：每个 G 编号都挂在某一行，无孤儿规则）。

使用者日常只需 `python3 run.py --doctor`（2 s，查环境）；`selftest.py` 是深度验证，发布前 / 换机器后跑。

**在哪跑 selftest**：全量在**部署副本**（`~/.agents/skills/fluctuation-analysis/`）里跑——
入口文件（run.py/bootstrap.py）按设计只在副本（见 `sync_skill.sh`「不同步」清单）。
在项目源目录跑会显示一项 `[SKIP]`——入口项与 name=目录名检查需要入口文件存在才有意义，
不是失败（这里**不写具体项数**，免得每次加检查就漂一次）。

## 0.1 规范校验（Agent Skills 开放标准）

格式不合规会让技能**静默不加载**——这与机器环境无关，是纯格式问题，所以单独查一次：

```bash
pip install skills-ref          # 官方参考实现；命令名是 agentskills（不是 skills-ref）
agentskills validate /path/to/fluctuation-analysis
# → Valid skill: /path/to/fluctuation-analysis
```

⚠️ **`validate .` 会失败**（`Directory name '' must match skill name`）——它拿路径最后一段当目录名，
所以必须传**以技能目录结尾**的路径。

`selftest.py` 里有一条轻量等价检查（不依赖该包）：name 正则、name 必须等于目录名、
description ≤1024、compatibility ≤500、frontmatter 无未知键、SKILL.md <500 行。

各 agent 对 frontmatter 的差异（2026-10-29 核对官方文档）：

| 字段 | Agent Skills 规范 | Claude Code | opencode |
|---|---|---|---|
| `name` / `description` | 必填 | 必填 | 必填 |
| `license` | 可选 | 接受但不生效 | 可选 |
| `compatibility` | 可选（≤500，**写环境要求**） | 接受但不生效 | 可选 |
| `metadata` | 可选（string→string） | 支持 | 可选 |
| `allowed-tools` | 可选，**空格分隔字符串** | 空格/逗号字符串 或 YAML 列表 | 忽略 |
| 未知键 | 未规定 | **拒绝加载**（`Unexpected key(s)`） | 忽略 |

## 1. 五层测试

| 层 | 测什么 | 怎么跑 | 谁需要跑 |
|---|---|---|---|
| L1 引擎单测 | 算法细节（128 项，含不变量） | `cd 项目源/window-calculator && pytest tests/ -q` | 改代码的人 |
| L2 技能自检 | 安装 + 已知答案 + 文档 | `python selftest.py` | **每个用户** |
| L3 参考对对拍 | 与 rSalvador/flan 数值一致 | 见 §4（需 R 环境） | 发布前 |
| L4 agent 问答 | 回答是否正确 | 见 §5（15 题抽 5 题） | 发布前 |
| L5 跨平台 | Linux/macOS/Windows | 见 §6 | 发布前 |

## 2. L1 引擎单测

196 项，覆盖：MSS 递推 vs 复合泊松、R&F 2000 验收、profile CI、逐管 Nt、部分涂板
（与显式二项混合逐点对拍 + rSalvador golden 值）、LRT/θ、样本量解析 vs 模拟、拟合退化保护，
以及**跨模块不变量**（MC 中位数≈点估计、decay 相位计入负荷与供给、报告恒等式、坏卡片被拦）。

```bash
cd <项目源>/window-calculator          # 仅维护者本机有；使用者不必跑 L1
python3 -m pytest tests/ -q          # 196 passed，约 2 分钟
```

## 3. L2 技能自检（重点）

`selftest.py` 的已知答案全部来自**技能自带**的 `examples/acceptance/`：

| 数据 | 真值/文献值 | 判据 |
|---|---|---|
| `RF2000.xlsx`（同一批 60 管，Excel 版） | 同 `RF2000.csv` | 与 `.csv` 结果完全一致（m=1.2274±0.002） |
| `LD1943_exp16.csv`（20 管，ε=0.4） | Luria & Delbrück 1943 第 16 次实验；Zheng 2017 给 m=1.18636 [0.5803, 2.0908] | m=1.1864±0.002 |
| `Demerec1945.csv`（30 管） | Demerec 1945；rSalvador 原文示例：m=10.84383，95%CI 8.650538–13.194765 | m=10.8438±0.002；CI 逐位一致 |
| `RF2000.csv`（60 管，完整 Table 3） | R&F 2000 Table 4：MSS m=1.23（中位数法才是 1.11）；本文 profile-LR → 0.931–1.576 | m=1.2274±0.002 |
| `plating_z02.csv` | ζ=0.2 部分涂板，模拟真值 m=3.0（seed 2026） | m̂=3.174±0.006；stderr 必须提示 Nt 口径 |
| `groupA/B.csv` | μ 相差 10 倍（同一 Nt 波动序列） | p<1e-10，θ∈[8,14] |
| 解析样本量 | μ=2e-9、Nt=3e8、95% 半宽 25% | C=131（旧实现双重取整给 135，已修） |

> 这些值冻结在 `selftest.py` 里。若引擎改动导致失败，**先确认是有意改动**再更新阈值。

## 4. L3 参考实现对拍（rSalvador / flan）

在一台装有 R 4.4 + rSalvador + flan 的机器上（内部示例环境）：

```bash
# ① 把我们的 CSV 转成 R 能读的 r 向量
Rscript -e 'library(rsalvador); d<-scan("RF2000.csv", what="", sep=",")[-c(1,2)]; ...'
# ② 同一份数据两边算 m 与 95%CI，要求相对差 < 1e-6（求解器容差）
```

已有验收记录：ζ=1 的 7 组数据点估计相对差 ≤2.5e-7；ζ<1 时 `newton.LD.plating(e=ζ)`
与我们的 m 差 2.2e-7、CI 差 1e-7（数值记录原在 `source_code_audit.md` §4（该文件 2026-11-13 已清理），本表自洽；
**2026-11-13 更正**：原指向的审计文件已清理，结论以本表与 CHANGELOG 为准）。

**★ 自建 R 常驻参考环境（2026-11-13 实跑复验）**：`conda env fluctuation`
（R 4.4.3 + **rsalvador 1.9** + flan；CRAN 上已 404，rSalvador 现以
[GitHub eeeeeric/rSalvador](https://github.com/eeeeeric/rSalvador) 发行，包名小写）。
数据 = `examples/acceptance/`，估计器 = `newton.LD` / `newton.LD.plating(e=…)`：

| 数据 | golden | R 参考实跑 | 判定 |
|---|---|---|---|
| `Demerec1945.csv`（30 管，整板） | m=10.84383（rSalvador 原文示例） | **10.84383** | ✅ 逐位一致 |
| `LD1943_exp16.csv`（20 管，**ε=0.4**） | m=1.18636（Zheng 2017） | `newton.LD.plating(e=0.4)` → **1.18636** | ✅ 逐位一致 |
| `plating_z02.csv`（ζ=0.2） | 本文引擎 m̂=3.174 | **3.173983** | ✅ 六位一致 |

⚠️ 当场踩到的坑：`LD1943_exp16` 用整板估计器 `newton.LD` 会得 **0.694**——
**部分涂板数据必须用 plating 估计器**，所有 golden 值都是修正后的。
另：rsalvador 的 `newton.LD.plating(data, e)` 只吃**一个公共 e、没有逐管 Nt**
（原审计文件已清理，见 CHANGELOG）；其 `foldchange`/`mcmc`/`LRT.LD.plating` 系列
是两批比较（F 及其 CI）的现成参照，可供 `compare` 模块日后对拍。

## 5. L4 agent 问答抽测

**做法**：把技能目录交给一个**干净上下文**的 agent（不让它看到本对话），只给题目，
要求每题 2–5 句、给关键数字、末尾标出处；允许它读文件/跑命令，但不许改文件。

**题库（15 题，只给题不给答案）**：
1. 为什么报速率不报频率？ 2. m≈1、20% 精度要多少管？ 3. m 的推荐范围？为什么不能说"1–15"？
4. m 偏离目标怎么调、优先级？ 5. N₀ 怎么定、有哪些约束？ 6. 只能涂一半培养物怎么算？
7. 各管 Nt 差异大怎么办？ 8. WT 与 mutator 怎么比、看什么？ 9. 检出 2 倍差异要多少管？
10. 窗口优化的目标与约束？ 11. 为什么缩短脉冲能减乘客、提高浓度不一定？ 12. 参数先测哪个、为什么？（报告是按 tornado 排的，不看死记的优先级）
13. 报告里的红灯有哪些、各意味着什么？ 14. 有一管 r=3000 能直接用吗？ 15. 与 rSalvador 不一致先查什么？

**判分**：每题按三个要点打分（结论正确 / 数字正确 / 出处正确），满分 3；总分 ≥80% 视为通过。
**答案要点不放在技能目录里**——否则被测 agent 可以直接读到答案（2026-10-28 首轮 L4 就是这样失效的）。
答案要点存于工作区 `archive/l4_answer_key.md`（仅维护者可见）。

## 6. L5 跨平台

| 平台 | 怎么测 | 判据 |
|---|---|---|
| Linux（本机） | `python selftest.py` | **21/21** ✅ 2026-10-29 |
| macOS | 同上（`install.sh` 建 venv） | 17/17（旧项数；新增 3 项未在 macOS 重测） |
| **Windows** ✅ 2026-10-29 实测 | conda 环境里 `run.cmd --doctor` → `run.cmd rate --input examples\acceptance\RF2000.csv --nts 3e8`（用户实机） | **跑通**：`run.cmd` 找到 conda/python，结果与 Linux 一致 |
| Windows（历史，2026-10-28） | conda + Python 3.10.22；`python selftest.py`，分别在 `C:\windows-test\` 与含空格+中文用户名的路径 | 13/13 × 2（当时项数）；⚠️ 双击 `.cmd` 无效（Explorer 的 PATH 无 conda python）→ 用 Anaconda Prompt 手动跑 |
| **无依赖的新机器**（Linux 等价环境） | 用不含 numpy 的 `python3 -m venv` 解释器跑 `run.py rate --input examples/acceptance/RF2000.csv --nts 3e8` | ✅ 自动建 `.venv`、装依赖、算出 m=1.2274；默认 PyPI 不可达时自动换镜像 |

Windows 上首次运行要确认四件事：`PYTHONUTF8=1` 生效（输出重定向不报 UnicodeEncodeError）、
`run.cmd` 能从非 conda 终端找到解释器（它会跳过 `WindowsApps` 那个假 python）、路径含空格时 `.cmd` 仍能用、
Excel 另存的 GBK CSV 能读（`_read_text` 回退）。

## 7. 发布检查单

- [x] L2 自检 21/21（Linux ✅ 2026-10-29；Windows 13/13 旧项数 ✅ 2026-10-28，py3.10.22 + conda）
- [x] **L5 Windows 实机**：`run.cmd` 跑通（2026-10-29，conda 环境，用户实测）
- [x] 自动安装：无依赖的新机器零配置跑通（2026-10-29；default PyPI 不可达时自动换镜像）
- [x] L1 单测 128/128（2026-10-29）
- [x] L3 对拍至少 2 组（ζ=1 差 2.5e-7、ζ<1 差 2.2e-7）
- [x] L4 盲测 45/45（2026-10-28，两轮）
- [x] 文档里的数字与当时的审计文件一致（已清理，结论见 CHANGELOG 对应条目；并已把"文档数字"纳入 selftest 自动检查）
- [x] 技能目录 `diff` 与项目源 `windowcalc/` 一致（2026-10-29；用 `bash sync_skill.sh` 同步）
- [x] **Agent Skills 规范校验通过**（`agentskills validate`；2026-10-29，name/长度/frontmatter 键全合规）
- [x] 版本/日期更新（v0.7.5，2026-10-30）；`SKILL.md` 的"四个用途"链接有效

## 8. 已知不可测项（如实告知用户）

1. **真实生物学数据**：没有实测的 strain_card，窗口结论只是模型预测（报告会标"★未测"）；
2. **Windows 边缘配置**：2026-10-29 已在 conda 环境用 `run.cmd` 实机跑通；
   未覆盖的是「conda 装在非默认盘且不在 PATH」「完全离线」这类少见配置；
3. **R 对拍**：参考环境已常驻（自建 R 参考机，`conda env fluctuation`：rsalvador 1.9 + flan，见 §4）；普通用户机器仍无 R；
4. **agent 行为**：只能抽测，无法证明"任何问法都答对"——文档写得越自包含，抽测越可信；
5. **完全离线的机器**：自动安装要能访问 PyPI 或其镜像。离线且本机没有 numpy 时无法零配置启动，
   但会打印明确的中文提示与三条出路（设镜像 / 指定已有解释器 / 手动装）。

---

## 9. 附：历史核对记录：2026-10-30 的软件行为核对（原 software-audit）

> 2026-11-15 由 `software_audit_2026-10-30.md` **整篇并入**，内容未删；入口见 `INDEX.md`。

> 对象：0.7.x 的 CLI 行为与文档声明是否一致。方法：逐条构造输入跑 CLI / 查代码常量 / 跑测试；
> 结果分 `✅ 一致`、`🔧 已修`、`⚠️ 有前提`。本轮**不**覆盖 Windows/macOS 实机与网络镜像路径。

| # | 文档声明 | 核对方法 | 结果 |
|---|---|---|---|
| 1 | `r>5000` 在两种模式都直接报错 | `r=5001` CSV | ✅ 解析层报错，提示列顺序可能反了 |
| 2 | m 搜索域上限 300，顶到边界给警示、不拒绝 | 3 管 `r=5000` | ✅ m=299.9999，CI 顶到 300，stderr 有“CI 触及适用域边界” |
| 3 | `C<20` 给精度警示 | 10 管数据 | ✅ stderr 提示按 m 算 C |
| 4 | `--nt-mode per_tube` 配标量 Nt 报错 | r-only CSV + `--nts` | ✅ 明确报“需要 CSV 带逐管 Nt 列” |
| 5 | 部分涂板 + jackpot（j=r/ζ>3000）报错 | r=1000, ζ=0.1 | ✅ 报“超出精确二项混合适用域”，建议整管/winsor |
| 6 | 部分涂板下 p0 交叉核验 | LD16, ζ=0.4 | 🔧 原未按稀释模型校正；已实现 `P0=G(1−ζ)`，LD16 得 0.9787（rSalvador 0.97868） |
| 7 | `power` 默认 0.95 = 95% CI 半宽；20% 在 m=1 需 138 管 | `power --mu 1/3e8` | ✅ 138 管；m=1.5 → 103 管 |
| 8 | `--conf-level 0` = 1σ；m=1 需 36 管 | 同上 | ✅ 36 管 |
| 9 | `--field` 只认规范名，旧名报错 | `--field M_WT` | ✅ 报“未知字段”，列可选值 |
| 10 | `fit-fc` 只认 `--mu-base` | `--mu-wt` | ✅ 报缺少 `--mu-base` |
| 11 | `optimize --protocol` 可省略；`--fill-defaults` 只补 4 项约定量 | 最小卡 | ✅ 用内置协议卡；6 项实测物理量仍报错并点名 |
| 12 | 菌株卡缺必填项直接报错（不静默默认） | 空卡 | ✅ 列缺失项 |
| 13 | `M_* < 1e-5` 提示“像单位点速率” | M_base=2.1e-9 | ✅ stderr 给 10³–10⁶ 倍说明 |
| 14 | `M_off > 3×M_base` 红灯 | M_off=0.02 | ⚠️ 只在有可行方案时出现在报告 §6；无可行解时只给瓶颈诊断（不误报） |
| 15 | `compare` 的 p 值解释 | groupA/B | ✅ stderr 提示“p>0.05 ≠ 无差异；看 CI” |
| 16 | `design-induction --ref-m` 无 C50 时被忽略并警告 | `--n-doses 3 --ref-m 6` | ✅ 明确警告并忽略 |
| 17 | `rate` 无 Nt 时不算 μ、片段留空 | r-only CSV | ✅ stderr 说明 |
| 18 | `.xls` 不支持 | 空 `.xls` | ✅ 提示另存为 .xlsx/.csv |
| 19 | 克隆封顶 5000（“模拟器与多项式法同域”） | 查代码 | ⚠️ 对外模拟（`power --sim`/`compare`）封顶 5000；内部性质测试 `simulate` 不封顶（代码注释已写明）；文档措辞已改 |
| 20 | 保护上限 binding / 不可行 | 单测 `test_protection_cap...` | ✅ 报告写明哪条 binding、实际 events、缺口 |
| 21 | `--doctor` 正常时 rc=0 | 直接跑 | ✅ rc=0 |
| 22 | Excel / GBK 输入 | selftest | ✅ 21/21 覆盖 |
| 23 | m<0.3 不建议（文档声明） | 低计数数据 | 🔧 原无提示；已加 stderr 警告（R&F 建议、精度差、提高 Nt/加 C） |

### 未覆盖 / 仍需注意

- **ζ 偏差数字**已重跑：m=1、C=60、8 次重复 → 21%/39%/72%（跳过 j>3000）；winsor 13%/30%/68%；稀释似然 1.01–1.05。
- **`power --sim` 的 MC 功效**已重跑：C=40 68–72%、60 87%、75 93%、90 96.5%（n_exp 100–200，seed 0–2）。
- **Windows/macOS** 实机未测；`run.cmd` 路径逻辑只在 Linux 上审过。
- **p0 稀释模型**已实现（`P0=G(1−ζ)`，bisection）；LD16 与 rSalvador 0.97868 一致。

### 结论

文档与 CLI 行为总体一致；本轮发现并修了 1 个真问题（部分涂板下 p0 假“有效”）+ 1 处措辞过宽（克隆封顶），其余为“有前提”而非错误。
