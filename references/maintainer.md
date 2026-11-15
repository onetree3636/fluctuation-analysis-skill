# 维护者参考

> 由 `SKILL.md` 按需加载。使用者不需要这些文件；改引擎、发版、复核时才读。

## 0. 自检

```bash
python3 run.py --doctor    # 2 s，只查环境
python3 selftest.py        # Windows: python selftest.py；约 20–60 s
```

`selftest.py` 覆盖：环境、入口自举、规范校验、**全部子命令**能否启动、Excel 输入、
已知答案对拍、端到端报告、文档断链与数字一致性、错误路径、**通用性扫描**、**护栏双向覆盖**、
**文档↔代码不变量**。红了按打印的复现命令查。

> ⚠ **本文件不写具体项数**。项数会变（每次加检查就 +1），写死过一次就被漏更新——
> 现在的项数以 `python3 selftest.py` 的末行为准（`N 通过，M 失败，K 跳过`）。

其中「规范校验」按 agentskills 的硬约束检查 `SKILL.md`：
frontmatter 必填 `name`/`description`、name 必须与目录名一致且为小写连字符、
description 1–1024、compatibility ≤500、未知 frontmatter 键拒绝、正文 <500 行。

「文档数字一致性」会现场算三个数字再比对：子命令个数、selftest 项数、解析样本量。
因此 `SKILL.md` 必须保留两个字面量：**`N 项：环境`** 与 **`N 个子命令`**（正则抓这两处）——
改代码时**不要手改 SKILL.md 里的数字**，跑一次 selftest，它会告诉你该填几。

## 1. 同步流程（方向：项目源 → 技能副本）

引擎开发源在维护者本机（不随技能分发），技能目录是运行时副本。改完**任何**文件后：

```bash
bash sync_skill.sh      # 单向覆盖，并做 diff 校验
cd ~/.agents/skills/fluctuation-analysis && python3 selftest.py
```

**同步范围（2026-11-11 起扩展为全量）**：`windowcalc/`、`examples/`、`docs/`、`references/`
以及 `SKILL.md`、`CHANGELOG.md`、`selftest.py`。
→ 也就是说**项目源是唯一真源**，`docs/` 与 `SKILL.md` 也在项目里改，改完 sync。

**不同步**（有意为之）：`README.md`（项目与技能内容不同）、`LICENSE`、入口脚本（`run.py` 等，稳定）、
`tests/`、`data/`、`review_*/audit` 等开发记录。

⚠ `sync_skill.sh` 用 `rsync --delete`，但**显式排除 `.venv/`**——它是 `bootstrap.py` 建的运行环境，
误删会导致用户下次运行重装依赖。改动同步脚本时**不要去掉这条排除**。

## 2. 与参照工具的关系（独立复核，非运行时依赖）

- rSalvador / flan（R）：拿同一批 CSV 两边算，点估计应一致（容差级）。
- 源码级审计、文献核对、历次专业复核都留在**开发工作区**（不随技能分发）：
  - `archive/window-calculator.bak-<日期>/` —— 历史快照：`source_code_audit.md`、
    `literature_audit_2026-10-28.md`、`review_2026-10-27.md`、`review_2026-10-29_professional.md`、
    `review_2026-10-30_skill-architecture.md`、`acceptance_zrway.md`；
  - `audit/fluctuation-analysis/` —— 技能两臂对照（带技能 / 不带技能）的契约、评测与逐次运行记录。
  - 这些文件名带日期，**别在这里写死**；目录里 `ls` 一下就有。
- 版本号只看 `windowcalc/__init__.py` 的 `__version__`（报告页眉与附录 A 都从它读）。

## 3. 发布检查单

- [ ] `python3 selftest.py` 全绿（末行 `0 失败`；`[SKIP]` 是设计行为）。
- [ ] `SKILL.md` 里的 `N 项：环境` / `N 个子命令` 与自检输出一致（**跑 selftest 让它告诉你，别手数**）。
- [ ] 所有相对链接可达（`SKILL.md` → `docs/`、`references/`；`docs/*.md` 内部）。
- [ ] description ≤ 500 字符（DSH 目录截断线）且 ≤ 1024（规范）。
- [ ] 不新增规范外的 frontmatter 键（如 `whenToUse` 不在 agentskills 键集，会被 Claude Code 拒载）。
- [ ] 打 zip / 拷目录时带上 `references/`（否则 SKILL.md 的指路是断的）。
