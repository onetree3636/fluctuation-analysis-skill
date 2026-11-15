"""fstep：比例估计的两个数学底座（Wilson 区间 + Newcombe 差值）。

**本模块只做算术，不碰数据**——pre/post 的配对、多板池化、背景扣除、
红旗判定全部在 `ramp.run_ledger` 里（账本解释器，唯一实跑路径）。
（2026-11-15 删除了 `run_step()` / `render()`：它们是 fstep 子命令时代的
单步引擎，自 2026-11-14 并入 `multitarget` 后已无生产调用者，
其行为测试已迁到 `tests/test_ledger_formula_pin.py` 的 run_ledger 用例。）

方法（文献）：
  - 单比例 95%CI：Wilson（Wilson 1927, JASA 22:209；小样本下优于正态近似：
    Brown, Cai & DasGupta 2001, Statist. Sci. 16:101）。
    k=0 用 rule of three：95% 上界 = 3/n（Hanley & Lippman-Hand 1983, JAMA 249:1743）。
  - 两比例之差（f_post − f_pre）：Newcombe (1998) 法 10——由两个 Wilson 区间混合
    （Newcombe, Statist. Med. 17:873）。
红灯口径（docs/protocol-fstep.md §3）：f_pre ≥ 0.5 × f_post → 背景淹没，本档不可用。
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

