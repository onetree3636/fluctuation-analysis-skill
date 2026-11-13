# 安装、平台与故障救援

> 由 `SKILL.md` 按需加载。只在首次使用、换平台或环境报错时需要。

## 0. 运行方式

唯一入口：

```bash
python3 run.py <子命令> …                        # Linux / macOS
"C:\…\fluctuation-analysis\run.cmd" <子命令> …  # Windows（含 Miniconda / Anaconda）
```

第一次运行会检查依赖；缺 `numpy` / `pyyaml` / `openpyxl` 就在技能目录建 `.venv` 装好，
然后自动重跑那条命令。**用户不需要执行任何安装步骤，也不需要知道发生过什么。**

- 首次慢（20–60 s，在下载依赖）是正常的，继续等；
  只有看到 `[windowcalc] 无法准备运行环境` 才算失败。
- 失败时**原样转达**它打印的「下一步」三选一（设镜像 / 指定解释器 / 手动安装），
  不要自己另发明安装方案。最常用的救援：
  `WINDOWCALC_PYTHON=$(command -v python3) python3 run.py …`（改用机器上已装 numpy 的 Python）。
- 不要在给用户的回答里汇报「我装了依赖 / 我建了 venv」这类过程，除非失败需要他决策。

## 1. 排查（不触发安装）

```bash
python3 run.py --doctor
```

输出解释器版本、缺哪些依赖、候选 venv 与结论。
退出码 `0` = 可直接用、`2` = 首次会自动装、`1` = 需人工介入。

## 2. 沙箱不能执行 shell 时

告诉用户跑一次上面的 `--doctor`，把输出贴回来。**不要**凭猜测判断依赖齐不齐。

## 3. Windows + Miniconda / Anaconda

用 **`run.cmd`**，它自己找 Python：

```bat
"C:\path\to\fluctuation-analysis\run.cmd" rate --input 计数.csv --nts 3e8
```

- **不要用 `py -3 run.py`**：`py.exe` 是 python.org 安装器才带的，Miniconda 机器上一般没有。
- **不要假设 `python` 在 PATH 上**：agent 从非 conda 终端（PowerShell / cmd 子进程）执行时，
  conda 的 python 通常不在 PATH。`run.cmd` 依次试：PATH → `py` → `%USERPROFILE%\miniconda3`
  → `%USERPROFILE%\anaconda3` → `%LOCALAPPDATA%` → `%ProgramData%`。
- 装在非默认盘时，指定一次绝对路径即可：
  `set WINDOWCALC_PYTHON=D:\miniconda3\python.exe`，然后再跑 `run.cmd`。
- **不会污染你的 conda**：依赖只装进技能目录的 `.venv\`，不往 base 或任何现有环境里 pip install。
- Windows 有 260 字符路径上限：技能目录太深时 `.venv` 会自动改到
  `%LOCALAPPDATA%\fluctuation-analysis\venv`。
- conda 的 Python 若 ≥3.9 且已有 numpy，则**直接用**、不建 venv（零副作用）。

## 4. 把技能导入不同 agent

整个目录复制到所用 agent 的 skills 目录，目录名保持 `fluctuation-analysis`：

| agent | 通常位置 |
|---|---|
| dsh | `~/.agents/skills/fluctuation-analysis/` |
| opencode | `~/.config/opencode/skills/fluctuation-analysis/` |
| Claude Code | `~/.claude/skills/fluctuation-analysis/`（项目级 `.claude/skills/`） |
| 其他（codex 等） | 见其文档；本技能只用通用 frontmatter 字段，不依赖专有能力 |

## 5. venv、镜像与自包含清单

- 技能目录：`~/.agents/skills/fluctuation-analysis/`（Windows：
  `C:\Users\<你>\.agents\skills\fluctuation-analysis\`）——整个目录自包含，
  可直接拷到别的机器或 agent（`examples/` 自带验收数据）。
- venv：优先**技能目录内 `.venv/`**（跟着技能走，删目录即卸载）；技能目录只读时回退到
  `~/.local/share/fluctuation-analysis/venv`；`WINDOWCALC_VENV=/path/to/venv` 可指定别处。
- 依赖安装默认走官方源；不通时自动换清华/阿里镜像。
- **目录符合 Agent Skills 规范**（https://agentskills.io/specification）：
  `SKILL.md` + `scripts/` 等价物（`run.py`、`bootstrap.py`、`windowcalc/`）
  + `references/` 等价物（`references/`、`docs/`）+ `assets/` 等价物（`examples/`）。
  `agentskills validate` 已通过。
- **新手入口**：`docs/tutorial.md`（10 分钟：40 管计数 → 一个 μ）；术语 `docs/glossary.md`。
- 自包含内容：`docs/knowledge-base.md`（原理/估计器/文献）、`docs/experiment-design.md`（设计与数据）、
  `docs/mutator-window.md`（窗口与乘客预算）、`docs/testing.md`（测试方案与发布检查单）、
  `examples/`（计数/记录表模板 + `acceptance/` 已知答案数据）。
- **Windows 注意**：各入口都设了 `PYTHONUTF8=1`（Linux/macOS 本就是 UTF-8，统一设置）——
  不设的话，输出重定向到文件时 `⚠`/`★` 等字符会 `UnicodeEncodeError`；路径含空格要加引号；
  CSV 用 UTF-8（带 BOM 也可）。
- 引擎本身纯 Python（numpy + pyyaml），无 R、无平台相关代码。
- 人想在终端敲短命令（Linux/macOS）：跑一次 `bash install.sh`，装用户级 venv 与
  `~/.local/bin/windowcalc`。**agent 不要用这条**，一律 `python3 run.py`。
- Windows 想敲短命令：跑一次
  `powershell -ExecutionPolicy Bypass -File install.ps1`，之后可用 `windowcalc.cmd`；
  **agent 直接用 `run.cmd` 即可**。
