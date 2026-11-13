"""fstep：铺板记录 → 可粘贴的升档供给频率 f_step（含 Wilson / Newcombe 区间）。

输入 CSV（examples/fstep_template.csv）：
    step,plate,when,cfu_plated,clones
    1,1,pre,3.0e8,0
    1,2,post,3.0e8,12
when ∈ {pre, post}：pre = 脉冲**前**铺的下一档浓度板（背景）；
                    post = 脉冲**后**铺的下一档浓度板（升档频率）。

方法（文献）：
  - 单比例 95%CI：Wilson（Wilson 1927, JASA 22:209；小样本下优于正态近似：
    Brown, Cai & DasGupta 2001, Statist. Sci. 16:101）。
    k=0 用 rule of three：95% 上界 = 3/n（Hanley & Lippman-Hand 1983, JAMA 249:1743）。
  - 两比例之差（f_post − f_pre）：Newcombe (1998) 法 10——由两个 Wilson 区间混合
    （Newcombe, Statist. Med. 17:873）。
红灯（docs/protocol-fstep.md §3）：f_pre ≥ 0.5 × f_post → 背景淹没，本档不可用。
"""

from __future__ import annotations

import math


def wilson_ci(k: int, n: int, z: float = 1.959963985) -> tuple[float, float, float]:
    """频率 k/n 的点估计与 95%CI（Wilson；k=0 用 rule of three 上界 3/n）。

    返回 (p, lo, hi)。
    """
    if n <= 0:
        raise ValueError(f"涂板 CFU 必须 > 0（当前 {n!r}）")
    if k < 0:
        raise ValueError(f"克隆数不能为负（当前 {k!r}）")
    if k > n:
        raise ValueError(f"克隆数 ({k}) > 涂板 CFU ({n})——数据有误")
    p = k / n
    if k == 0:
        return p, 0.0, 3.0 / n            # rule of three [Hanley & Lippman-Hand 1983]
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z / den * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, max(0.0, center - half), min(1.0, center + half)


def newcombe_diff(post: tuple[float, float, float],
                  pre: tuple[float, float, float]) -> tuple[float, float]:
    """两独立比例之差 f_post − f_pre 的 95%CI（Newcombe 1998 法 10）。

    post/pre 各是 (p, lo, hi)（Wilson 区间）。
    """
    p2, l2, u2 = post
    p1, l1, u1 = pre
    d = p2 - p1
    lo = d - math.sqrt((p2 - l2) ** 2 + (u1 - p1) ** 2)
    hi = d + math.sqrt((u2 - p2) ** 2 + (p1 - l1) ** 2)
    return d, max(-1.0, lo), min(1.0, hi)


def run_step(step: str, plates: list[tuple[str, float, int]]) -> dict:
    """一个 step 的全部板 → f_pre / f_post / f_step（扣背景）与判定。

    plates: [(when, cfu_plated, clones), ...]；when ∈ {"pre", "post"}。
    """
    def pool(w: str) -> tuple[int, int]:
        sel = [(c, k) for w_, c, k in plates if w_ == w]
        if not sel:
            raise ValueError(f"step {step} 缺 when={w!r} 的板（配对设计两边都要铺）")
        return sum(c for c, _ in sel), sum(k for _, k in sel)

    n1, k1 = pool("pre")
    n2, k2 = pool("post")
    f_pre = wilson_ci(k1, n1)
    f_post = wilson_ci(k2, n2)
    d, dlo, dhi = newcombe_diff(f_post, f_pre)
    f_step = max(0.0, d)
    red = (f_pre[0] >= 0.5 * f_post[0]) if f_post[0] > 0 else False
    return dict(step=step, k_pre=k1, n_pre=n1, f_pre=f_pre,
                k_post=k2, n_post=n2, f_post=f_post,
                f_step=f_step, ci=(dlo, dhi), red=red)


def render(results: list[dict]) -> str:
    """把各 step 的结果渲染成报告文本（含可粘贴的 --f-step 行）。"""
    L: list[str] = []
    for r in results:
        p1, lo1, hi1 = r["f_pre"]
        p2, lo2, hi2 = r["f_post"]
        L.append(f"== step {r['step']} ==")
        L.append(f"  f_pre  = {r['k_pre']} / {r['n_pre']:.3g} = {p1:.3g}"
                 f"   95%CI [{lo1:.2g}, {hi1:.2g}]")
        L.append(f"  f_post = {r['k_post']} / {r['n_post']:.3g} = {p2:.3g}"
                 f"   95%CI [{lo2:.2g}, {hi2:.2g}]（Wilson）")
        L.append(f"  f_step = f_post − f_pre = {r['f_step']:.3g}"
                 f"   Newcombe95 [{r['ci'][0]:.2g}, {r['ci'][1]:.2g}]")
        if r["red"]:
            L.append("  ⛔ 背景淹没：f_pre ≥ 0.5×f_post——本档不可用"
                     "（先解决假抗性/交叉污染，protocol-fstep.md §3）")
        elif r["f_step"] > 0:
            L.append(f"  → 粘贴：--f-step {r['f_step']:.3g}")
        else:
            L.append("  → f_step ≤ 0：本步没有可用的升档供给（见上）")
        L.append("")
    return "\n".join(L)
