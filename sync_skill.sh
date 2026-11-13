#!/bin/bash
# 把项目源（本目录）同步到技能运行时副本。
#
# 方向是**单向**的：项目源 → 技能副本。以前是手工 cp，结果技能里改过一行 note、
# 项目源没改，两边分叉（见 review_2026-10-29 §P2-14）。改完就跑这个脚本，再跑 selftest：
#
#     bash sync_skill.sh && cd ~/.agents/skills/fluctuation-analysis && python3 selftest.py
#
# 同步范围（2026-11-11 扩展）：
#   windowcalc/  examples/  docs/  references/   +   SKILL.md  CHANGELOG.md  selftest.py
#
# **不同步**（有意为之，别"顺手补全"）：
#   README.md      项目与技能的 README 内容不同（一个是开发笔记，一个是 GitHub 门面）
#   LICENSE        技能侧独有
#   run.py 等入口  跨平台自举脚本，稳定；改动频率为零
#   tests/ data/   项目专属（tests 是开发测试，data 是验收数据）
#   review_*/audit  项目专属的开发评审记录
#
# **排除**（★ 最关键的一条）：
#   .venv/         bootstrap.py 自动建的运行环境。rsync --delete 一旦把它算进去，
#                  用户下次运行要重装依赖（20–60 s，且可能受网络影响）。
#   __pycache__/  字节码
#   .pytest_cache/
#   *.bak*         手工备份
#
set -euo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
SKILL="${FLUCTUATION_SKILL_DIR:-$HOME/.agents/skills/fluctuation-analysis}"

if [ ! -d "$SKILL" ]; then
    echo "找不到技能目录：$SKILL" >&2
    echo "用 FLUCTUATION_SKILL_DIR=<路径> 指定" >&2
    exit 2
fi

DIRS=(windowcalc examples docs references)
FILES=(SKILL.md CHANGELOG.md selftest.py)
EXCLUDES=(--exclude '__pycache__' --exclude '.pytest_cache' --exclude '.venv' --exclude '*.bak*')

for d in "${DIRS[@]}"; do
    [ -d "$SRC/$d" ] || { echo "缺目录：$SRC/$d" >&2; exit 2; }
done

mkdir -p "$SKILL/windowcalc" "$SKILL/examples" "$SKILL/docs" "$SKILL/references"
if command -v rsync >/dev/null 2>&1; then
    for d in "${DIRS[@]}"; do
        rsync -a --delete "${EXCLUDES[@]}" "$SRC/$d/" "$SKILL/$d/"
    done
    for f in "${FILES[@]}"; do
        rsync -a "$SRC/$f" "$SKILL/$f"
    done
else
    # 没有 rsync 就逐文件覆盖（不删多余文件；.venv 自然不受影响）
    for d in "${DIRS[@]}"; do
        cp -a "$SRC/$d/." "$SKILL/$d/"
    done
    for f in "${FILES[@]}"; do
        cp -a "$SRC/$f" "$SKILL/$f"
    done
fi

echo "已同步 ${DIRS[*]} + ${FILES[*]} → $SKILL"
if command -v diff >/dev/null 2>&1; then
    rc=0
    for d in "${DIRS[@]}"; do
        diff -r -q --exclude=__pycache__ --exclude=.pytest_cache --exclude=.venv --exclude='*.bak*' \
            "$SRC/$d" "$SKILL/$d" || rc=1
    done
    for f in "${FILES[@]}"; do
        diff -q "$SRC/$f" "$SKILL/$f" || rc=1
    done
    if [ "$rc" -eq 0 ]; then
        echo "副本一致 ✅（含 docs/references/SKILL.md/CHANGELOG/selftest）"
    else
        echo "⚠ 仍有差异，见上" >&2
    fi
fi

# 运行环境不能被同步删掉——显式检查一次
if [ -d "$SKILL/.venv" ]; then
    echo "运行环境保留：$SKILL/.venv ✅"
fi
echo
echo "下一步：cd \"$SKILL\" && python3 selftest.py"
