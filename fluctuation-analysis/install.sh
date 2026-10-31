#!/bin/bash
# fluctuation-analysis 技能安装（**给人用**）：venv + 依赖 + ~/.local/bin/windowcalc
#
# agent 不需要跑这个：agent 直接用 `python3 run.py <子命令>`，缺依赖会自动装。
# 只有当你想在终端里敲 `windowcalc` 时才需要本脚本。
#
#   bash install.sh
#   WINDOWCALC_VENV=/custom/venv bash install.sh
set -e
SKILL_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV="${WINDOWCALC_VENV:-$HOME/.local/share/fluctuation-analysis/venv}"
BIN="${WINDOWCALC_BIN:-$HOME/.local/bin/windowcalc}"

die() { echo "" >&2; echo "安装失败：$1" >&2; shift; for l in "$@"; do echo "  $l" >&2; done; exit 2; }

# ---- 1. 找一个 Python 3.9+ -------------------------------------------------
PY=""
for cand in python3 python; do
    if command -v "$cand" >/dev/null 2>&1; then PY="$cand"; break; fi
done
[ -n "$PY" ] || die "没找到 Python 3。" \
    "装一个 Python 3.9 或更新的版本再重试：" \
    "  Ubuntu/Debian: sudo apt install python3 python3-venv" \
    "  macOS:         brew install python" \
    "  Windows:       https://www.python.org/downloads/ （安装时勾 Add to PATH）"

"$PY" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)' \
    || die "Python 版本太旧（$( "$PY" -c 'import sys;print(sys.version.split()[0])' )），需要 3.9+。" \
           "如果机器上有多个 Python，可指定其一：WINDOWCALC_PYTHON=/path/to/python3 bash install.sh"

# ---- 2. venv 模块是否可用（Debian 最小安装常缺） ---------------------------
if ! "$PY" -c 'import venv, ensurepip' >/dev/null 2>&1; then
    die "当前 Python 缺少 venv / ensurepip 模块，无法创建独立环境。" \
        "  Ubuntu/Debian: sudo apt install python3-venv" \
        "  CentOS/RHEL:   sudo dnf install python3-libs" \
        "  或者：把已有 numpy/pyyaml 的解释器指定给它" \
        "    WINDOWCALC_PYTHON=/path/to/python3 bash install.sh"
fi

# ---- 3. 建 venv 并装依赖（与 run.py 共用同一套逻辑） ------------------------
echo "使用 Python: $("$PY" -c 'import sys;print(sys.executable, sys.version.split()[0])')"
echo "venv: $VENV"
WINDOWCALC_VENV="$VENV" "$PY" "$SKILL_DIR/bootstrap.py" --force \
    || die "依赖安装失败（见上方提示）。" \
           "可设置镜像后重试：export WINDOWCALC_PIP_INDEX=https://pypi.tuna.tsinghua.edu.cn/simple"

# ---- 4. 写全局命令 windowcalc ---------------------------------------------
if [ -w "$(dirname "$BIN")" ] || mkdir -p "$(dirname "$BIN")" 2>/dev/null; then
    cat > "$BIN" <<WRAPEOF
#!/bin/bash
# windowcalc —— fluctuation-analysis 技能入口（由 install.sh 生成）
export PYTHONUTF8=1
export PYTHONPATH="$SKILL_DIR"
exec "$VENV/bin/python" -m windowcalc "\$@"
WRAPEOF
    chmod +x "$BIN"
    echo "已生成命令: $BIN"
    case ":$PATH:" in
        *":$(dirname "$BIN"):"*) ;;
        *) echo "提示: $(dirname "$BIN") 不在 PATH 里，用绝对路径或在 shell 配置里加上它。" ;;
    esac
else
    echo "提示: 无法写入 $BIN（权限不足），跳过全局命令；仍可用 python3 $SKILL_DIR/run.py" >&2
fi

# ---- 5. 冒烟测试 -----------------------------------------------------------
"$VENV/bin/python" -m windowcalc --help >/dev/null 2>&1 \
    || die "安装完成但无法运行，请把下面的输出发给维护者：" \
           "  $VENV/bin/python $SKILL_DIR/run.py --doctor"

echo ""
echo "OK: windowcalc 已安装。试一下： windowcalc rate --input 计数.csv --nts 3e8"
echo "自检: python3 $SKILL_DIR/selftest.py"
