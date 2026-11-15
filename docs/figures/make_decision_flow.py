# -*- coding: utf-8 -*-
"""生成 docs/figures/decide-plating-and-timing.png —— 两个决策：涂多少？何时能往下走？

通用性政策：图内不出现任何数值参数（只用 $C_{target}$ / $f_{step}$ / $N_t$ 这类符号）。
布局：每个面板 = 标题 + 副题 + 一句提示 + 一行流程 + 下方分支框。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Polygon
import matplotlib.pyplot as plt

CJK = ["Noto Sans CJK SC", "WenQuanYi Micro Hei", "WenQuanYi Zen Hei",
       "AR PL UKai CN", "AR PL KaitiM GB", "Droid Sans Fallback"]
avail = {f.name for f in fm.fontManager.ttflist}
for c in CJK:
    if c in avail:
        plt.rcParams["font.sans-serif"] = [c]
        break
plt.rcParams["axes.unicode_minus"] = False

BLUE, ORANGE, GREEN, RED = "#5b84b1", "#e0954a", "#6b8e5a", "#c0392b"

fig = plt.figure(figsize=(13.2, 8.0), dpi=300)
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0],
                      left=0.02, right=0.98, top=0.868, bottom=0.030, hspace=0.08)
axU, axD = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
for ax in (axU, axD):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def rbox(ax, x, y, w, h, text, fc, ec, fs=8.8):
    ax.add_patch(mpatches.FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                         boxstyle="round,pad=0.006,rounding_size=0.018",
                                         facecolor=fc, edgecolor=ec, lw=1.1, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color="#1a1a1a",
            zorder=4, linespacing=1.35)


def diamond(ax, x, y, w, h, text, fc, ec, fs=8.8):
    ax.add_patch(Polygon([[x, y + h / 2], [x + w / 2, y], [x, y - h / 2], [x - w / 2, y]],
                         closed=True, facecolor=fc, edgecolor=ec, lw=1.1, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color="#1a1a1a", zorder=4,
            linespacing=1.3)


def arr(ax, xy1, xy2, label="", color="#555555", ls="-", fs=8.3, dy=0.030):
    ax.add_patch(FancyArrowPatch(xy1, xy2, arrowstyle="-|>", mutation_scale=12, lw=1.25,
                                 color=color, linestyle=ls, zorder=2, shrinkA=0, shrinkB=0))
    if label:
        ax.text((xy1[0] + xy2[0]) / 2, max(xy1[1], xy2[1]) + dy, label,
                fontsize=fs, color=color, ha="center", va="bottom", zorder=5,
                bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.92))


# ======================= 上：涂多少 =======================
axU.text(0.012, 0.955, "问题一：这一档要涂多少？", fontsize=11.8, fontweight="bold", color=BLUE)
axU.text(0.012, 0.885, "答案永远是一个「算出来的比例」，不是「全涂」——精度只由数到的克隆数决定",
         fontsize=9.2, color="#555555")
axU.text(0.012, 0.815,
         "⚠ 常见误区：把「全涂」当保险。全涂既不提高精度，又会让板数随培养物总量爆炸。",
         fontsize=8.9, color=RED)

yU, hU = 0.545, 0.185
X1, W1 = 0.078, 0.140
X2, W2 = 0.268, 0.200
X3, W3 = 0.505, 0.200
X4, W4 = 0.855, 0.235

rbox(axU, X1, yU, W1, hU, "① 定 $C_{target}$\n要数到几个克隆\n（如 50–100）", "#eaf1f7", BLUE)
diamond(axU, X2, yU, W2, 0.30, "② 本档 $f_{step}$\n已知？", "#fff6e8", ORANGE)
rbox(axU, X3, yU, W3, hU + 0.02,
     "$N_{plated}=C_{target}/f_{step}$\n$V=N_{plated}/N_t$\n涂板比例 $=V/V_{tot}$", "#eaf1f7", BLUE)
rbox(axU, X4, yU, W4 + 0.05, hU + 0.02,
     "③ 铺 2 块 $c_{N+1}$ 板\n（= 本档 $f_{post}$ + 下一档 $f_{pre}$）\n另加 CFU 分母 3 块",
     "#eaf1f7", BLUE)

arr(axU, (X1 + W1 / 2 + 0.004, yU), (X2 - W2 / 2 - 0.004, yU))
arr(axU, (X2 + W2 / 2 + 0.004, yU), (X3 - W3 / 2 - 0.004, yU), "是")
arr(axU, (X3 + W3 / 2 + 0.004, yU), (X4 - (W4 + 0.05) / 2 - 0.004, yU))

rbox(axU, 0.268, 0.200, 0.320, 0.225,
     "否：三档夹逼 $10^5/10^6/10^7$ CFU 各 1–2 块\n"
     "取落在可数区（30–300 克隆）的那一块\n→ 反算 $f_{step}$；第二档起一次到位",
     "#fdf1e2", ORANGE, fs=8.4)
arr(axU, (X2, yU - 0.152), (X2, 0.315), "否", color=ORANGE, dy=0.012)

# ======================= 下：何时能往下走 =======================
axD.text(0.012, 0.955, "问题二：什么时候能接下一档？", fontsize=11.8, fontweight="bold",
         color=ORANGE)
axD.text(0.012, 0.885, "答案：不等板——阶梯事先登记，主线靠浊度闸判，板只做审计",
         fontsize=9.2, color="#555555")
axD.text(0.012, 0.815,
         "★ 只有四种情况才真的等板：首档（要定 $f$）、红灯（$f_{pre}\\geq0.5f_{post}$）、末档收尾、换过流程。",
         fontsize=8.9, color="#555555")

yD, hD = 0.545, 0.185
Y1, V1 = 0.078, 0.140
Y2, V2 = 0.258, 0.150
Y3, V3 = 0.452, 0.175
Y4, V4 = 0.658, 0.185
Y5, V5 = 0.858, 0.170

rbox(axD, Y1, yD, V1, hD, "① 脉冲收尾\n洗菌（回收率对照）", "#fdf1e2", ORANGE)
rbox(axD, Y2, yD, V2, hD, "② 同一管分装\n90–95% ｜ 5–10%", "#fdf1e2", ORANGE)
rbox(axD, Y3, yD, V3, hD, "③ 接下一档浓度液体\n立即开跑（主线）", "#fdf1e2", ORANGE)
diamond(axD, Y4, yD, V4, 0.30, "④ 24–48 h\n浑浊？", "#fff6e8", ORANGE)
rbox(axD, Y5, yD, V5, hD, "⑤ 推进下一窗\n（阶梯已登记）", "#eef4ea", GREEN)

arr(axD, (Y1 + V1 / 2 + 0.004, yD), (Y2 - V2 / 2 - 0.004, yD))
arr(axD, (Y2 + V2 / 2 + 0.004, yD), (Y3 - V3 / 2 - 0.004, yD))
arr(axD, (Y3 + V3 / 2 + 0.004, yD), (Y4 - V4 / 2 - 0.004, yD))
arr(axD, (Y4 + V4 / 2 + 0.004, yD), (Y5 - V5 / 2 - 0.004, yD), "是")

rbox(axD, 0.300, 0.205, 0.310, 0.140,
     "板臂（并行）：涂板 → 板孵育 → 双时间点读\n→ $f_{step}$ 回填账本（审计，不挡推进）",
     "#fbf3e3", "#c9954a", fs=8.3)
arr(axD, (Y2, yD - 0.102), (Y2, 0.281), color="#c9954a", ls="--")

rbox(axD, 0.700, 0.270, 0.300, 0.115,
     "浑浊 → 步长偏小（$f\\approx1$，等于没选）\n记录；下一档把步长加大", "#fdeaea", RED, fs=8.3)
rbox(axD, 0.700, 0.090, 0.300, 0.115,
     "48–72 h 仍不浑浊 → 步长偏大\n从上一档甘油/备份回滚，步长减半", "#fdeaea", RED, fs=8.3)
arr(axD, (Y4, yD - 0.152), (Y4, 0.330), "否", color=RED, dy=0.010)
arr(axD, (Y4, 0.212), (Y4, 0.150), color=RED, ls=":")

fig.suptitle("两个决策：涂多少 · 什么时候能接下一档", fontsize=14.5, fontweight="bold", y=0.968)
fig.text(0.5, 0.922,
         "符号化示意：图中不含任何数值；$C_{target}$=目标克隆数，$f_{step}$=实测升档频率，$N_t$=周期末密度",
         ha="center", fontsize=9.4, color="#555555")

out = "decide-plating-and-timing.png"
fig.savefig(out, dpi=300, facecolor="white")
print("saved:", out)
