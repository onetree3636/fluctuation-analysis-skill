# -*- coding: utf-8 -*-
"""生成 docs/figures/propagation-topologies.png —— 板传代 vs 液传代（符号化对照）。

通用性政策：图内不出现任何数值参数；只画流程、判据与公式符号。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

CJK = ["Noto Sans CJK SC", "WenQuanYi Micro Hei", "WenQuanYi Zen Hei",
       "AR PL UKai CN", "AR PL KaitiM GB", "Droid Sans Fallback"]
avail = {f.name for f in fm.fontManager.ttflist}
for c in CJK:
    if c in avail:
        plt.rcParams["font.sans-serif"] = [c]
        break
plt.rcParams["axes.unicode_minus"] = False

BLUE, ORANGE = "#5b84b1", "#e0954a"
FIG_W, FIG_H = 13.0, 8.2

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[0.66, 0.34],
                      left=0.02, right=0.98, top=0.850, bottom=0.03,
                      wspace=0.10, hspace=0.14)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axT = fig.add_subplot(gs[1, :])
for ax in (axA, axB, axT):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def box(ax, x, y, w, h, text, fc, ec, fs=9.2, fw="normal", tc="#1a1a1a", ls="-"):
    ax.add_patch(mpatches.FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                         boxstyle="round,pad=0.008,rounding_size=0.02",
                                         facecolor=fc, edgecolor=ec, lw=1.1, ls=ls, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc,
            fontweight=fw, zorder=4, linespacing=1.35)


def arrow(ax, x1, y1, x2, y2, label="", color="#555555"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13,
                                 lw=1.3, color=color, zorder=2,
                                 shrinkA=1.5, shrinkB=1.5))
    if label:
        ax.text((x1 + x2) / 2 + 0.012, (y1 + y2) / 2, label, fontsize=8.4,
                color=color, ha="left", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.9))


# ================= A. 板传代 =================
axA.set_title("A · 板传代（板既是测量、也是传代）", fontsize=11.5, fontweight="bold",
              color=BLUE, pad=8)
SA, WA = 0.29, 0.52
yA = [0.905, 0.760, 0.615, 0.470, 0.325, 0.180]
textsA = ["脉冲收尾（窗口结束）",
          "涂 $c_{N+1}$ 板（克隆数要够）",
          "板内长 $n_{plate}$ 代（$\\log_2$ 每克隆细胞数）",
          "数克隆 → $f_{step}$ = 克隆数/CFU",
          "整板刮洗合并（不挑单克隆）",
          "悬液即下一档种子 → 开下一窗"]
for y, t in zip(yA, textsA):
    box(axA, SA, y, WA, 0.105, t, "#eaf1f7", BLUE, fs=8.8)
for i in range(len(yA) - 1):
    arrow(axA, SA, yA[i] - 0.053, SA, yA[i + 1] + 0.053)

# 旁注
box(axA, 0.775, 0.760, 0.37, 0.150,
    "瓶颈 = 板上那 $C$ 个克隆\n（每个都是独立谱系）", "#fff6e8", "#c9954a", fs=8.4)
box(axA, 0.775, 0.470, 0.37, 0.150,
    "偏差来源：刮不干净、\n只挑大克隆 → 系统偏移", "#fdeaea", "#c0392b", fs=8.4)
box(axA, 0.775, 0.180, 0.37, 0.150,
    "档间乘客 = $n_{plate}\\cdot M_{off}$\n（$n_{plate}$ 是结构下限）", "#eef4ea", "#6b8e5a", fs=8.4)
for yy in (0.760, 0.470, 0.180):
    axA.plot([SA + WA / 2 + 0.005, 0.59], [yy, yy], color="#999999", lw=0.9,
             ls=":", zorder=1)

# ================= B. 液传代 =================
axB.set_title("B · 液传代（板只做审计）", fontsize=11.5, fontweight="bold",
              color=ORANGE, pad=8)
SB, WB = 0.29, 0.52
yB = [0.905, 0.760, 0.615, 0.470, 0.325, 0.180]
textsB = ["脉冲收尾（窗口结束）",
          "离心洗一次（去诱导剂）+ 回收率",
          "分装：90–95% 液体 ｜ 5–10% 涂板",
          "接 $c_{N+1}$ 液体（主线开跑，不等板）",
          "存活者 $f\\cdot N_t$ 长到 $N_0$（$n_{grow}$ 代）",
          "长到 $N_0$ → 开下一窗"]
for y, t in zip(yB, textsB):
    box(axB, SB, y, WB, 0.105, t, "#fdf1e2", ORANGE, fs=8.8)
for i in range(len(yB) - 1):
    arrow(axB, SB, yB[i] - 0.053, SB, yB[i + 1] + 0.053)

box(axB, 0.775, 0.760, 0.37, 0.170,
    "无瓶颈：$f\\cdot N_t$ 个存活者\n全部保留", "#eef4ea", "#6b8e5a", fs=8.4)
box(axB, 0.775, 0.470, 0.37, 0.150,
    "风险：液相选择不严\n（步长 < MIC → WT 也长）", "#fdeaea", "#c0392b", fs=8.4)
box(axB, 0.775, 0.180, 0.37, 0.160,
    "档间乘客 = $n_{grow}\\cdot M_{off}$\n$n_{grow}=\\log_2(N_0/(f N_t))$",
    "#eef4ea", "#6b8e5a", fs=8.4)
for yy in (0.760, 0.470, 0.180):
    axB.plot([SB + WB / 2 + 0.005, 0.59], [yy, yy], color="#999999", lw=0.9,
             ls=":", zorder=1)

# 板臂（并行）画在 B 面板的“分装”旁
axB.add_patch(FancyArrowPatch((SB + WB / 2 + 0.004, 0.615), (0.643, 0.615),
                              arrowstyle="-|>", mutation_scale=11, lw=1.1,
                              color="#c9954a", ls="--", zorder=2))
axB.text(0.648, 0.615, "涂板（只数，不刮）\n→ 账本回填 $f_{step}$",
         fontsize=8.0, va="center", ha="left", color="#7a4a10", linespacing=1.3)

# ================= 对照表 =================
axT.set_xlim(0, 1); axT.set_ylim(0, 1)
axT.add_patch(mpatches.FancyBboxPatch((0.005, 0.03), 0.99, 0.94,
                                      boxstyle="round,pad=0.004,rounding_size=0.012",
                                      facecolor="#fbfaf6", edgecolor="#d8d8d2", lw=1.0))
rows = [
    ("选择严格性", "高：必须走完全部克隆世代", "中：能启动生长即可（可加严）"),
    ("群体瓶颈", "$C$ 个独立谱系", "无（$f\\cdot N_t$ 全保留）"),
    ("档间乘客", "$n_{plate}\\cdot M_{off}$", "$n_{grow}\\cdot M_{off}$"),
    ("什么时候必须用它", "杀菌药 / 要单克隆表型 / 液相选择漏 / 液体成团",
     "默认方案（文献五篇渐升协议都是液体群体传代）"),
]
xcol = [0.02, 0.30, 0.68]
axT.text(xcol[0], 0.92, "判据", fontsize=9.4, fontweight="bold", va="center")
axT.text(xcol[1], 0.92, "A · 板传代", fontsize=9.4, fontweight="bold", va="center", color=BLUE)
axT.text(xcol[2], 0.92, "B · 液传代", fontsize=9.4, fontweight="bold", va="center", color="#b8791f")
yy = 0.74
for name, a, b in rows:
    axT.text(xcol[0], yy, name, fontsize=8.9, va="center", fontweight="bold")
    axT.text(xcol[1], yy, a, fontsize=8.7, va="center")
    axT.text(xcol[2], yy, b, fontsize=8.7, va="center")
    yy -= 0.215
    if yy > 0.05:
        axT.plot([0.015, 0.985], [yy + 0.105, yy + 0.105], color="#e2e2dc", lw=0.8)

fig.suptitle("两种传代拓扑：板传代 vs 液传代", fontsize=14.5, fontweight="bold", y=0.972)
fig.text(0.5, 0.930,
         "符号化示意：图中不含任何数值；$n_{plate}=\\log_2$(每克隆细胞数)，$n_{grow}=\\log_2(N_0/(f N_t))$",
         ha="center", fontsize=9.8, color="#555555")

out = "propagation-topologies.png"
fig.savefig(out, dpi=300, facecolor="white")
print("saved:", out)
