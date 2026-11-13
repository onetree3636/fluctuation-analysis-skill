# 维护者参考

> 由 `SKILL.md` 按需加载。使用者不需要这些文件；改引擎、发版、复核时才读。

## 0. 自检

```bash
python3 run.py --doctor    # 2 s，只查环境
python3 selftest.py        # Windows: python selftest.py；约 20–60 s
```

`selftest.py` 共 **21 项**：环境、入口自举、规范校验、8 个子命令、Excel 输入、
4 组已知答案、端到端报告、文档数字一致性、错误路径、文档断链。红了按打印的复现命令查。

其中「规范校验」按 agentskills 的硬约束检查 `SKILL.md`：
frontmatter 必填 `name`/`description`、name 必须与目录名一致且为小写连字符、
description 1–1024、compatibility ≤500、未知 frontmatter 键拒绝、正文 <500 行。

「文档数字一致性」会现场算三个数字再比对：子命令个数、selftest 项数、解析样本量。
因此 `SKILL.md` 必须保留字面量 **`21 项：环境`** 与 **`8 个子命令`**（正则直接抓这两处）。

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
- 源码级审计、功能矩阵、已验证等价范围与历次复核都在开发源的 `*.md` 里
  （`dev/source_code_audit.md`、`acceptance_*.md`、`review_2026-10-2*.md`、
  `dev/literature_audit_2026-10-28.md`）。
- 最新一轮专业复核与修复见 `dev/review_2026-10-29_professional.md`（引擎 v0.6.0 起，现 v0.7.0）。
- SKILL.md 架构评审见 `dev/review_2026-10-30_skill-architecture.md`（渐进式披露重构的依据）。

## 3. 发布检查单

- [ ] `python3 selftest.py` 21 项全绿。
- [ ] `SKILL.md` 仍含 `21 项：环境` 与 `8 个子命令` 两个字面量。
- [ ] 所有相对链接可达（`SKILL.md` → `docs/`、`references/`；`docs/*.md` 内部）。
- [ ] description ≤ 500 字符（DSH 目录截断线）且 ≤ 1024（规范）。
- [ ] 不新增规范外的 frontmatter 键（如 `whenToUse` 不在 agentskills 键集，会被 Claude Code 拒载）。
- [ ] 打 zip / 拷目录时带上 `references/`（否则 SKILL.md 的指路是断的）。
