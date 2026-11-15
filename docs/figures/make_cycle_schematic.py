# -*- coding: utf-8 -*-
"""生成 docs/figures/cycle-liquid-propagation.png —— 一个浓度周期的**符号化**示意图。

设计约束（docs/generic-policy.md）：
  · 图里**不出现任何数值参数**（不写 M_on=…、不写 f=…、不写代数）；
  · 纵轴只画形状，刻度标成符号；瀑布只标公式，不标数值；
  · 目的是"教结构"，数值交给用户自己的菌株卡。

用法：python make_cycle_schematic.py     （需 matplotlib + 中文字体）
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

CJK = ["Noto Sans CJK SC", "WenQuanYi Micro Hei", "WenQuanYi Zen Hei",
       "AR PL UKai CN", "AR PL KaitiM GB", "Droid Sans Fallback"]
avail = {f.name for f in fm.fontManager.ttflist}
for c in CJK:
    if c in avail:
        plt.rcParams["font.sans-serif"] = [c]
        break
plt.rcParams["axes.unicode_minus"] = False

# 形状参数（只是画图用的相对高度，不是任何生物量）
Y_BASE, Y_OFF, Y_ON = 0.0, 0.42, 1.0

fig = plt.figure(figsize=(12.8, 7.6), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.85],
                      left=0.092, right=0.975, top=0.845, bottom=0.075,
                      width_ratios=[1.25, 1.0], hspace=0.42, wspace=0.16)
ax1 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[1, 1])

# ---------------- 时间轴（示意） ----------------
seg = {"p0": (0.0, 0.75), "p1": (0.75, 3.05), "p2": (3.05, 3.75),
       "p4": (3.75, 3.95), "p3": (3.95, 6.55), "p5": (6.55, 7.55)}
for key, col in [("p0", "#e3ecf4"), ("p1", "#fde9d0"), ("p2", "#ece2f0"),
                 ("p4", "#fff2cc"), ("p3", "#e3ecf4"), ("p5", "#f2f2f2")]:
    ax1.axvspan(*seg[key], color=col, zorder=0)

xs, ys = [], []
def segadd(a, b, f, n=60):
    t = np.linspace(a, b, n); xs.append(t); ys.append(f(t))
segadd(*seg["p0"], lambda t: np.full_like(t, Y_OFF))
r = 0.35
segadd(seg["p1"][0], seg["p1"][0] + r, lambda t: Y_OFF + (Y_ON - Y_OFF) * (t - seg["p1"][0]) / r)
segadd(seg["p1"][0] + r, seg["p1"][1], lambda t: np.full_like(t, Y_ON))
t2 = np.linspace(*seg["p2"], 60)
xs.append(t2); ys.append(Y_OFF + (Y_ON - Y_OFF) * np.exp(-(t2 - seg["p2"][0]) / 1.8))
xs.append(np.array([seg["p4"][0], seg["p4"][1]]))
resid = Y_OFF + (Y_ON - Y_OFF) * 0.25
ys.append(np.array([resid, resid]))
t3 = np.linspace(seg["p3"][0], seg["p3"][0] + 0.9, 60)
xs.append(t3); ys.append(Y_OFF + (resid - Y_OFF) * np.exp(-(t3 - seg["p3"][0]) * 2.2))
segadd(seg["p3"][0] + 0.9, seg["p3"][1], lambda t: np.full_like(t, Y_OFF))
segadd(*seg["p5"], lambda t: np.full_like(t, Y_OFF))
ax1.plot(np.concatenate(xs), np.concatenate(ys), color="#2c3e50", lw=2.6, zorder=5)

ax1.axhline(Y_ON, color="#d35400", lw=1.0, alpha=0.85, zorder=2)
ax1.axhline(Y_OFF, color="#2e6da4", lw=1.0, ls="--", alpha=0.85, zorder=2)
ax1.axhline(Y_BASE, color="#7f8c8d", lw=1.0, ls=":", alpha=0.9, zorder=2)
ax1.set_ylim(-0.30, 1.46); ax1.set_yticks([Y_BASE, Y_OFF, Y_ON])
ax1.set_yticklabels(["WT 基线\n$M_{base}$", "泄漏\n$M_{off}$", "诱导\n$M_{on}$"], fontsize=9)
for lbl, cc in zip(ax1.get_yticklabels(), ["#7f8c8d", "#2e6da4", "#d35400"]):
    lbl.set_color(cc)
ax1.set_xlim(-0.05, 7.62); ax1.set_xticks([])
ax1.set_ylabel("全基因组突变负荷 M（示意，不标数值）", fontsize=10)

phases = [("①", "非诱导", 0.36), ("②", "诱导（+诱导剂）", 1.80), ("③", "洗脱尾", 3.28),
          ("④ 洗菌", None, 3.86), ("⑤", "液扩（下一档浓度）", 5.05), ("⑥", "门控", 7.05)]
for num, name, x in phases:
    ax1.text(x, 0.995, num, transform=ax1.get_xaxis_transform(), ha="center", va="top",
             fontsize=(10.0 if name is None else 12.5), fontweight="bold", color="#111111", zorder=7)
    if name is not None:
        ax1.text(x, 0.935, name, transform=ax1.get_xaxis_transform(), ha="center", va="top",
                 fontsize=9.2, color="#333333", zorder=7)

ax1.axvline(seg["p4"][0], ymin=0.0, ymax=0.10, color="#2e6da4", lw=1.1, ls="--", alpha=0.85, zorder=3)

tr = blended_transform_factory(ax1.transData, ax1.transAxes)
ax1.add_patch(mpatches.Rectangle((3.78, 0.035), 7.40 - 3.78, 0.075, transform=tr,
                                 facecolor="#fbe9ca", edgecolor="#c9954a", lw=1.0, ls="-.", zorder=4))
ax1.text(3.83, 0.0725, "⑤′ 板臂（并行审计，不挡下一档）：涂 $c_{N+1}$ 板 + 分母 → 板孵育 → 双时间点读 → $f_{step}$ 回填账本",
         transform=tr, fontsize=8.6, va="center", ha="left", color="#7a4a10", zorder=8,
         bbox=dict(boxstyle="round,pad=0.20", fc="#fbe9ca", ec="none"))

ax1.annotate("", xy=(seg["p5"][1] - 0.05, 0.20), xytext=(seg["p5"][0] + 0.12, 0.20),
             arrowprops=dict(arrowstyle="->", color="#555555", lw=1.6))
ax1.text(seg["p5"][1] - 0.09, 0.235, "下一档", fontsize=9, color="#555555", va="bottom", ha="right")

# ---------------- 符号化瀑布 ----------------
names = ["② 诱导窗口\n$n_{win}M_{on}$", "⑤ 液扩\n$n_{grow}M_{off}$",
         "⑤ $\\tau_{off}$ 尾\n$(M_{on}{-}M_{off})\\Sigma$", "周期后\n$n_{post}M_{off}$"]
frac = [0.62, 0.46, 0.72, 0.22]          # 只为排版：柱子相对高度
cols = ["#e8a04c", "#5b84b1", "#7fb3d5", "#95a5a6"]
b = 0.0
for i, (n, v, c) in enumerate(zip(names, frac, cols)):
    ax2.bar(i, v, bottom=b, width=0.62, color=c, edgecolor="#444444", lw=0.8)
    b += v
ax2.axhline(b, color="#a02c2c", lw=1.2, ls="--")
ax2.text(-0.52, b + 0.05, "合计 $E[P]$（= 四段之和）", color="#a02c2c", fontsize=9.6,
         ha="left", va="bottom", fontweight="bold")
ax2.set_xticks(range(4)); ax2.set_xticklabels(names, fontsize=8.6)
ax2.set_yticks([])
ax2.set_ylabel("乘客（仅示意）", fontsize=9.5)
ax2.set_ylim(0, b * 1.28); ax2.set_xlim(-0.6, 3.6)
ax2.set_title("每克隆乘客四段账（只标公式，不标数值）", fontsize=10.5)
ax2.spines[["top", "right", "left"]].set_visible(False)

# ---------------- 右侧：分装与红线 ----------------
ax3.axis("off")
txt = (
 "④ 那一管分五份（同一管、t=0）\n"
 "─────────────\n"
 "① 90–95% → 接下一档液体（传代主线）\n"
 "② 5–10%  → 涂板 + CFU 分母（审计）\n"
 "③ 3 管   → 甘油 −80 °C（推进前存）\n"
 "④ 一份   → 测序沉淀\n"
 "⑤ 1–2 mL → 备份活菌 4 °C（回滚）\n"
 "\n"
 "★ 四条红线\n"
 "─────────────\n"
 "· 不全涂；不从已长起来的液体取 $f$ 板\n"
 "· 洗前/洗后的分母不混用；刮下的克隆\n"
 "   不得在含诱导剂条件下扩增"
)
ax3.text(0.0, 1.0, txt, fontsize=8.9, va="top", linespacing=1.30,
         bbox=dict(boxstyle="round,pad=0.55", fc="#f7f7f5", ec="#cccccc"))

fig.suptitle("一个浓度周期长什么样（液传代 + 涂板只做审计）", fontsize=14, fontweight="bold", y=0.965)
fig.text(0.5, 0.912,
         "示意图：纵轴只表示形状，横轴不按时间等比；图中没有任何数值参数——数由你自己的菌株卡决定",
         ha="center", fontsize=9.4, color="#555555")

out = "cycle-liquid-propagation.png"
fig.savefig(out, dpi=300, facecolor="white")
print("saved:", out)
