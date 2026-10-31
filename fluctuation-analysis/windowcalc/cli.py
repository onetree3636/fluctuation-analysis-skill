"""CLI：windowcalc rate（波动分析原始计数 → 速率）与 windowcalc optimize（窗口优化）。

用法：
    python -m windowcalc rate --input data/wt_rif.csv --label WT
    python -m windowcalc optimize --strain examples/strain_card.example.yaml \
        --protocol examples/protocol_card.example.yaml --out 报告.md
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import os
import sys
from pathlib import Path

from .mss_mle import analyze, analyze_per_tube
from .compare import foldchange, lrt_equal_rates
from .model import (InducerSpec, Phase, RateParam, StrainCard,
                    card_problems, phase_problems)
from .optimize import Decision, optimize
from .report import render_optimize_report


class CardError(ValueError):
    """卡片（YAML/字段）相关的输入错误——错误提示要指向卡片文档，而不是 CSV 教程。"""


# ----------------------------------------------------------------------------- rate

_EXCEL_EXTS = (".xlsx", ".xlsm", ".xltx", ".xltm")
_EXCEL_OLD_EXTS = (".xls", ".xlsb")


def _read_text(path: str) -> str:
    """读成文本：CSV 走编码探测，Excel（.xlsx/.xlsm）先转成 CSV 文本。

    中文 Windows 上 Excel「另存为 CSV」默认写本地代码页（GBK），直接按 UTF-8
    读会 UnicodeDecodeError——不能只报一句"编码错"就完事。

    Excel 文件**不要求用户另存**：直接读第一个工作表，转成 CSV 文本，交给下游
    同一套解析（表头跳过、列检查、行号提示全部复用）。
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in _EXCEL_OLD_EXTS:
        raise ValueError(
            f"{path}：{ext} 是 Excel 97-2003 / 二进制工作簿，读不了——"
            "请在 Excel 里「另存为」.xlsx 或 .csv 后重试")
    if ext in _EXCEL_EXTS:
        return _xlsx_to_csv_text(path)

    last = None
    for enc in ("utf-8-sig", "gb18030"):
        try:
            with open(path, encoding=enc) as f:
                return f.read()
        except UnicodeDecodeError as exc:
            last = exc
    raise ValueError(f"{path}: 无法解码（已试 UTF-8 与 GB18030）——请在 Excel/VSCode 里"
                     "另存为 UTF-8 编码的 CSV 后重试") from last


def _xlsx_to_csv_text(path: str) -> str:
    """Excel → CSV 文本（只读第一个工作表）。与 CSV 共用下游全部解析逻辑。"""
    rows = list(_iter_xlsx_rows(path))
    if not rows:
        raise ValueError(f"{path}: 工作表是空的（没读到任何行）")
    buf = io.StringIO()
    csv.writer(buf).writerows(rows)
    return buf.getvalue()


def _iter_xlsx_rows(path: str):
    """逐行产出单元格字符串（转义交给 csv.writer，含逗号/引号/换行的单元格都不会坏）。"""
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise ValueError(
            f"{path}: 读 Excel 需要 openpyxl，但当前环境没装。两条路：\n"
            "  1) 在 Excel 里「另存为」CSV 后重试（最快）\n"
            "  2) pip install openpyxl（用 run.py 启动会自动装）") from exc

    try:
        wb = load_workbook(path, data_only=True, read_only=True)
    except Exception as exc:                 # 加密 / 损坏 / 扩展名名不副实
        raise ValueError(
            f"{path}: 打不开这个 Excel 文件（{type(exc).__name__}: {exc}）。"
            "确认它没加密码、没损坏；或另存为 .xlsx 再试") from exc

    try:
        if not wb.sheetnames:
            raise ValueError(f"{path}: 工作簿里没有工作表")
        ws = wb[wb.sheetnames[0]]           # 多个 sheet 时默认第一个
        for row in ws.iter_rows(values_only=True):
            yield ["" if v is None else str(v) for v in row]
    finally:
        try:
            wb.close()
        except Exception:
            pass


def _looks_like_tube_index(raw: list) -> bool:
    """三列且「第一列小整数、第二列小整数、第三列像 Nt(CFU)」→ 多半多了一列管号。"""
    if len(raw) < 3:
        return False
    try:
        first = float(raw[0])
        second = float(raw[2])
    except (TypeError, ValueError):
        return False
    return (abs(first - round(first)) < 1e-9 and 0 <= first <= 999
            and second > 1e4)


def _row_hint(rv: float, nv: float | None) -> str:
    """根据读到的两列给一句"最可能的错因"。"""
    if rv > 5000:
        return ("第 1 列数值像 Nt（单板菌落数不可能 >5000）——列顺序可能反了，应为 r,Nt")
    if nv is None:
        return ""
    if nv <= 0:
        return ("第二列是 0/负数：如果它是**频率**（r/Nt）请换成原始 Nt（CFU/管）；"
                "如果它是**菌落数**，说明列顺序反了（应为 r,Nt）")
    if rv > 5000:
        return ("第 1 列数值像 Nt（单板菌落数不可能 >5000）——列顺序可能反了，应为 r,Nt")
    if 0 <= rv <= 50 and 0 < nv <= 50 and abs(nv - int(nv)) < 1e-9:
        return ("两列都是小数：第 1 列像管号、第 2 列像菌落数——请删掉管号列，只留 r 和 Nt")
    if nv < 1e4 or nv > 1e11:
        return ("Nt 不在常见区间 1e4–1e11（CFU/管）——确认是不是 CFU/mL 或频率")
    return ""


def _read_xy_csv(path: str, what: str) -> list:
    """两列数值 CSV（x,y）→ [(x, y), …]。

    自动跳过表头/注释行（任何"第一列不是数字"的首行，包括中文表头），
    报错时指出行号与期望格式——不要甩 "could not convert string to float"。
    """
    out = []
    for lineno, row in enumerate(csv.reader(io.StringIO(_read_text(path), newline="")), 1):
        cells = [c.strip() for c in row]
        if not any(cells) or cells[0].startswith("#"):
            continue
        try:
            x = float(cells[0])
        except ValueError:
            if not out:
                continue                     # 表头（中文/英文都算）
            raise ValueError(f"{path} 第 {lineno} 行：第一列 {cells[0]!r} 不是数字"
                             f"（{what} 需要两列数值，可带表头）")
        if len(cells) < 2 or not cells[1]:
            raise ValueError(f"{path} 第 {lineno} 行只有一列——{what} 需要两列：x,y")
        try:
            y = float(cells[1])
        except ValueError:
            raise ValueError(f"{path} 第 {lineno} 行：第二列 {cells[1]!r} 不是数字"
                             f"（{what} 需要两列数值）")
        out.append((x, y))
    if not out:
        raise ValueError(f"{path}: 没读到数据行——{what} 需要两列数值（可带表头）")
    return out


def _read_counts_csv(path: str):
    """CSV：r[,Nt[,volume_mL]] 每行一管。

    返回 (rs, nts, vols, notes)：
      - nts 为 **CFU/管**（若给了 volume_mL 列，则按 Nt 是 CFU/mL 自动换算）；
      - vols 为逐管体积（mL）或 None；
      - notes 为解析提示（单位可疑等），由调用方打印到 stderr。
    """
    rs, nts, vols, notes = [], [], [], []
    vol_seen_declared = None
    with io.StringIO(_read_text(path), newline="") as f:
        for lineno, row in enumerate(csv.reader(f), 1):
            if not row or not row[0].strip():
                continue
            raw = [x.strip() for x in row]
            try:
                rv = float(raw[0])
            except ValueError:
                if not rs:                     # 数据行之前的非数字行 = 表头
                    continue
                raise ValueError(f"{path} 第 {lineno} 行：r={raw[0]!r} 不是数字")
            if rv < 0 or rv != int(rv):
                raise ValueError(f"{path} 第 {lineno} 行读到 r={rv:g}：菌落数必须是非负整数"
                                 + (f"（提示：{_row_hint(rv, None)}）" if rv > 5000 else ""))
            if rv > 5000:
                raise ValueError(f"{path} 第 {lineno} 行读到 r={rv:g}：超出单板可数范围（>5000）"
                                 f"（提示：{_row_hint(rv, None)}）")
            rs.append(int(rv))

            nv = None
            if len(raw) > 1 and raw[1]:
                try:
                    nv = float(raw[1])
                except ValueError:
                    raise ValueError(f"{path} 第 {lineno} 行读到 (r={rv:g}, Nt={raw[1]!r})：Nt 不是数字")
                if nv <= 0:
                    hint = _row_hint(rv, nv)
                    if _looks_like_tube_index(raw):
                        hint = (f"看起来你多了一列**管号**：第一列是管号（{raw[0].strip()}），"
                                f"第二列才是菌落数、第三列才是 Nt。删掉第一列后应是 r,Nt"
                                f"（此行会读成 r={raw[1].strip()}, Nt={raw[2].strip()}）")
                    raise ValueError(f"{path} 第 {lineno} 行读到 (r={rv:g}, Nt={nv:g})：Nt 必须为正"
                                     f"（提示：{hint}）")
                nts.append(nv)

            if len(raw) > 2 and raw[2]:
                try:
                    vol = float(raw[2])
                except ValueError:
                    raise ValueError(f"{path} 第 {lineno} 行：体积={raw[2]!r} 不是数字")
                if vol <= 0:
                    raise ValueError(f"{path} 第 {lineno} 行：体积必须为正（mL）")
                vols.append(vol)
                vol_seen_declared = True
            elif vol_seen_declared:
                raise ValueError(f"{path} 第 {lineno} 行缺体积列：volume_mL 必须每行都给或整列不给")

    if not rs:
        raise ValueError(f"{path}: 没有读到任何计数行（应是每行一管的 r[,Nt[,volume_mL]]）")
    if nts and len(nts) != len(rs):
        raise ValueError(
            f"{path}: Nt 列只有 {len(nts)} 个值但计数有 {len(rs)} 行——Nt 必须逐管配对或整列留空")
    if vols and len(vols) != len(rs):
        raise ValueError(f"{path}: 体积列只有 {len(vols)} 个值但计数有 {len(rs)} 行")

    if vols:
        nts = [n * v for n, v in zip(nts, vols)]
        notes.append(f"Nt 按 CFU/mL × 体积换算为 CFU/管（体积均值 {sum(vols)/len(vols):.3g} mL）")
    elif nts:
        lo, hi = min(nts), max(nts)
        if hi < 1.0:
            notes.append(f"Nt={hi:.3g} 小于 1：这像**频率**（r/Nt）而不是每管细胞数——"
                         "请改成原始 Nt（CFU/管）")
        elif lo < 1e4 or hi > 1e11:
            notes.append(f"Nt 范围 {lo:.3g}–{hi:.3g} 不在常见区间 1e4–1e11（CFU/管）："
                         "若你测的是 CFU/mL，请加第三列 volume_mL 或改用 --nt-unit per_mL")
    return rs, (nts if nts else None), (vols if vols else None), notes


def cmd_rate(args) -> int:
    nt_unit = args.nt_unit or "per_tube"      # None = 用户没显式选 → 下面大声提示
    rs, nts, vols, notes = _read_counts_csv(args.input)
    for n in notes:
        print(f"注：{n}", file=sys.stderr)

    has_nt = bool(nts) or args.nts is not None
    kwargs = {}
    if vols and nt_unit == "per_mL":
        raise ValueError("CSV 已带 volume_mL 列（自动换算），不要再加 --nt-unit per_mL")
    if nts:
        vals = [float(x) for x in nts]
        src = "CSV 的 Nt 列"
    elif args.nts is not None:
        if args.nts <= 0:
            raise ValueError("--nts 必须为正")
        vals = [float(args.nts)]
        src = "--nts"
    else:
        vals = None
        src = None

    if vals is not None and nt_unit == "per_mL":
        if not args.volume_mL or args.volume_mL <= 0:
            raise ValueError("--nt-unit per_mL 需要同时给 --volume-mL（培养物体积，mL，需为正）")
        vals = [v * args.volume_mL for v in vals]
        print(f"注：{src} 按 CFU/mL × {args.volume_mL:g} mL 换算为 CFU/管", file=sys.stderr)

    if vals is None:
        print("警告：未提供 Nt，μ 无法计算（m 仍有效）；YAML 片段里 μ 留空", file=sys.stderr)
        kwargs["nts"] = 1.0
    elif len(vals) == 1:
        kwargs["nts"] = vals[0]
    else:
        kwargs["nts"] = vals

    # ---- 输入回显：让用户先核对"我读到了什么" ----
    if has_nt:
        nt_mean = (float(sum(kwargs["nts"]) / len(kwargs["nts"]))
                   if isinstance(kwargs["nts"], list) else float(kwargs["nts"]))
        nt_disp = f"{nt_mean:.3g} CFU/管"
    else:
        nt_disp = "未提供"
    print(f"输入：C={len(rs)} 管｜r 合计={sum(rs)}（最大 {max(rs)}）｜Nt={nt_disp}")
    if args.nts is not None and args.nt_unit is None:
        # 打到 stdout，紧跟在输入回显后面——新手是顺着往下读的
        print("      ⚠ --nts 按 **CFU/管** 解释。若你手里是 CFU/mL，"
              "请加 --nt-unit per_mL --volume-mL <培养体积 mL>（否则 μ 会差一个体积倍）")
    if not has_nt:
        print("      提示：没有 Nt 就无法算 μ；CFU/管 = CFU/mL × 培养体积(mL)")
    elif args.nt_mode != "per_tube" and isinstance(kwargs.get("nts"), list) \
            and len(kwargs["nts"]) > 1:
        mean_v = sum(kwargs["nts"]) / len(kwargs["nts"])
        cv = (sum((v - mean_v) ** 2 for v in kwargs["nts"]) / len(kwargs["nts"])) ** 0.5 / mean_v \
            if mean_v > 0 else 0.0
        if cv > 0.1:
            print(f"      提示：逐管 Nt 的 CV={cv:.0%}——换 --nt-mode per_tube 可消除 Nt 波动"
                  "带来的系统性低估（Ycart & Veziris 2014）", file=sys.stderr)

    if args.nt_mode == "per_tube":
        if not isinstance(kwargs["nts"], list) or len(kwargs["nts"]) != len(rs):
            print("错误：--nt-mode per_tube 需要 CSV 带逐管 Nt 列（与 r 一一配对）", file=sys.stderr)
            return 1
        res = analyze_per_tube(rs, kwargs["nts"], fraction_plated=args.fraction_plated,
                               label=args.label)
    else:
        res = analyze(rs, nts=kwargs["nts"], fraction_plated=args.fraction_plated,
                      label=args.label)

    print(f"== {args.label} ==")
    print(f"C={res.n_cultures}  mean_r={res.mean_r:.2f}  median_r={res.median_r:.0f}  p0={res.p0:.3f}")
    print(f"m (MSS-MLE)     = {res.m:.4f}   95%CI [{res.m_ci[0]:.4f}, {res.m_ci[1]:.4f}]")
    print(f"m (p0 交叉核验) = {res.m_p0:.4f}   {'有效区间内' if res.m_p0_valid else '无效区间（仅参考）'}")
    if has_nt:
        print(f"μ = m/mean(Nt)  = {res.mu:.3e} per CFU per passage"
              f"   95%CI [{res.mu_ci[0]:.3e}, {res.mu_ci[1]:.3e}]")
        print("     （参照量级：野生型 E. coli ~1e-9；典型 mutator 1e-8–1e-6——"
              "若差几个数量级，先查 Nt 单位）")
    for n in res.notes:
        print(f"注：{n}", file=sys.stderr)
    print("\n-- 粘到 strain_card --")
    print(res.as_yaml_snippet(args.label, suppress_mu=not has_nt,
                              key=getattr(args, "field", None)))
    if args.fraction_plated == 1.0:
        print("提示：本结果按**整管涂板**计算。若只涂了一部分（如 0.1/2 mL），"
              "请加 --fraction-plated ζ，否则 μ 会偏低；"
              "另外每行应是**一管**（多板 r 加总），不是一块板一行。", file=sys.stderr)
    return 0


# ----------------------------------------------------------------------------- optimize

_REQUIRED_STRAIN = {
    "M_base": "基线速率（每基因组每代）",
    "M_on_ref": "参考剂量下的诱导态速率",
    "M_off": "off 态泄漏速率",
    "mu_supply": "供给端（第二）标记的速率（每分裂，选定浓度实测）",
    "tau_on_h": "诱导上升沿（h）",
    "tau_off_h": "洗脱衰变尾（h）",
    "td_h_on": "诱导态倍增时间（h）",
    "td_h_off": "off 态倍增时间（h）",
    "Nt_saturation": "饱和密度（CFU/管）",
    "ref_dose_mM": "参考剂量（mM）",
}

_CARD_FIELDS = ("M_base", "M_on_ref", "M_off", "mu_supply", "background_false")


def _field_arg(s: str) -> str:
    """--field 的取值：必须是规范字段名，否则报错并列出可选项。"""
    if s not in _CARD_FIELDS:
        raise argparse.ArgumentTypeError(
            f"未知字段 {s!r}；可选：{' / '.join(_CARD_FIELDS)}")
    return s

_REQUIRED_DECISION = {
    "k_target_clones": "目标独立克隆数 k",
    "P_success_min": "成功率下限（如 0.95）",
    "T_pulse_range_h": "脉冲时长搜索区间 [lo, hi]（h）",
    "N0_min": "接种量下限（CFU）",
    "N0_max": "接种量上限（CFU）",
}


def _require(d: dict, required: dict, what: str) -> None:
    """缺必填参数就报错（不静默用默认值）。"""
    missing = [f"{k}：{v}" for k, v in required.items() if d.get(k) is None]
    if missing:
        raise CardError(
            f"{what}缺少必填参数，**不会用默认值代替**：\n  - "
            + "\n  - ".join(missing)
            + "\n请实测后填 `status: measured`；若只能先用假设值，请显式写 `status: assumed`，"
              "报告会把它们标成 ★未测。\n"
              "起步：技能目录 `examples/` 下有现成模板（strain_card*.yaml / protocol_card*.yaml），"
              "复制一份逐项替换即可——没测的字段保持 `assumed` 也能先跑。")


def _rate(d, key, default=None):
    v = d.get(key)
    if v is None:
        return default
    if isinstance(v, dict):
        ci = v.get("ci")
        return RateParam(float(v["value"]),
                         (float(ci[0]), float(ci[1])) if ci else None,
                         str(v.get("status", "assumed")), str(v.get("note", "")))
    return RateParam(float(v))


def _scalar(d, key, default=None):
    """结构参数（td / Nt_saturation / ref_dose）：裸数值或 `{value, ci, status}` 都吃。

    返回 float 或 RateParam——下游一律用 `model.val()` 取值，两种形式等价。
    这修掉了"文档说每个参数可带 ci，但 td_h_on 写成 dict 就 TypeError"。
    """
    v = d.get(key)
    if v is None:
        return default
    if isinstance(v, dict):
        if v.get("value") is None:
            raise ValueError(f"{key}: 写成字典时必须给 value（现在是 {v!r}）")
        return _rate(d, key)
    return float(v)


def _load_yaml(path: str) -> dict:
    """读 YAML 卡片：语法错误/空文件/顶层不是映射 → 一句中文错误（不甩 traceback）。"""
    import yaml
    try:
        d = yaml.safe_load(_read_text(path))
    except yaml.YAMLError as exc:
        raise CardError(f"{path}: YAML 语法错误——{exc}") from exc
    if d is None:
        raise CardError(f"{path}: 文件是空的（至少要有一行 key: value）")
    if not isinstance(d, dict):
        raise CardError(f"{path}: 顶层必须是一个映射（key: value），读到 {type(d).__name__}")
    return d


def _load_plasmids(d) -> list:
    """质粒清单：名称可省略（默认"质粒1/2…"）；单质粒可用简写 plasmid_kb/plasmid_copies。"""
    out = []
    for i, x in enumerate(d.get("plasmids") or []):
        if x.get("kb") is None:
            raise CardError(f"plasmids 第 {i + 1} 项缺 kb——写 {{name: 载体, kb: 16.0, copies: 20}}")
        out.append({"name": str(x.get("name") or f"质粒{i + 1}"),
                    "kb": float(x["kb"]), "copies": float(x.get("copies", 1))})
    if not out and d.get("plasmid_kb") is not None:
        out.append({"name": "质粒", "kb": float(d["plasmid_kb"]),
                    "copies": float(d.get("plasmid_copies", 1))})
    return out


def _load_protected(d, path: str = "") -> list:
    """保护区清单：每项必须显式给 kb 与 selection（不再静默默认 silent）。"""
    out = []
    for i, x in enumerate(d.get("protected") or []):
        nm = str(x.get("name") or f"保护区{i + 1}")
        if x.get("kb") is None:
            raise CardError(
                f"{path}: protected 第 {i + 1} 项（{nm}）缺 kb——"
                "写成 {name: 产物基因簇, kb: 22, copies: 1, selection: silent}")
        if "selection" not in x:
            raise CardError(
                f"{path}: protected 第 {i + 1} 项（{nm}）缺 selection——"
                "必须显式写 selection: silent（静默乘客，计入克隆质量）"
                "或 lethal（抗生素选择压下被淘汰，只影响产量）。两者结论含义不同，软件不替你猜。")
        out.append({"name": nm, "kb": float(x["kb"]),
                    "copies": float(x.get("copies", 1)),
                    "note": str(x.get("note", "")),
                    "selection": str(x["selection"]),
                    "role": str(x.get("role", "other"))})
    return out


#: --fill-defaults 用：菌株卡缺的必填项 → 默认假设值（全部 status: assumed，报告标 ★未测）
_DEFAULT_CARD = {
    "M_base": 0.004, "M_on_ref": 0.4, "M_off": 0.012, "mu_supply": 2.0e-8,
    "tau_on_h": 2.0, "tau_off_h": 4.0,
    "td_h_on": 2.5, "td_h_off": 2.0,
    "Nt_saturation": 3.0e8, "ref_dose_mM": 0.3,
}


def load_strain_card(path: str, fill_defaults: bool = False) -> StrainCard:
    d = _load_yaml(path)
    if fill_defaults:
        filled = []
        for key, value in _DEFAULT_CARD.items():
            if d.get(key) is None:
                d[key] = {"value": value, "status": "assumed",
                          "note": "默认值（--fill-defaults）"}
                filled.append(key)
        if filled:
            print("警告：菌株卡缺少 " + "、".join(filled)
                  + "——已用默认假设值填充（报告会标 ★未测）；有实测值请替换。",
                  file=sys.stderr)
    _require(d, _REQUIRED_STRAIN, f"{path}（菌株卡）")
    ind = d.get("inducer")
    inducer = None
    if ind:
        def _ind_val(key, default=None):
            """inducer 字段：裸数值或 `{value, ci}` 都吃 → (value, ci)。"""
            v = ind.get(key, default)
            if v is None:
                return None, None
            if isinstance(v, dict):
                if v.get("value") is None:
                    raise CardError(f"{path}: inducer.{key} 写成字典时必须给 value")
                ci = v.get("ci")
                return float(v["value"]), ((float(ci[0]), float(ci[1])) if ci else None)
            return float(v), None

        # 兼容两种写法：c50_mM（推荐，= 速率半饱和）与 ec50_mM（旧名）
        c50_key = "c50_mM" if "c50_mM" in ind else "ec50_mM"
        a_v, a_ci = _ind_val("amplitude")
        c_v, c_ci = _ind_val(c50_key)
        n_v, n_ci = _ind_val("n_hill")
        q_v, _ = _ind_val("q", 1.0)
        if a_v is None or c_v is None or n_v is None:
            raise CardError(f"{path}: inducer 块需要 amplitude / c50_mM / n_hill 三个数值")
        # 可选 ci 块：inducer: {..., ci: {amplitude: [lo,hi], c50: [...], n: [...]}}
        # 由 fit-induction 的 bootstrap 区间直接粘过来 → MC 才会传播剂量轴的不确定性
        block = ind.get("ci") or {}
        cis = {}
        for name, direct in (("amplitude", a_ci), ("c50", c_ci), ("n", n_ci)):
            if direct:
                cis[name] = direct
        for name, keys in (("amplitude", ("amplitude",)),
                           ("c50", ("c50", "c50_mM", "ec50_mM")),
                           ("n", ("n", "n_hill"))):
            if name in cis:
                continue
            for k in keys:
                if k in block:
                    lo, hi = block[k]
                    cis[name] = (float(lo), float(hi))
                    break
        inducer = InducerSpec(
            name=str(ind.get("name", "诱导剂")),
            amplitude=a_v,
            ec50_mM=c_v,
            n_hill=n_v,
            mixing=str(ind.get("mixing", "linear")),
            q=q_v if q_v is not None else 1.0,
            ci=cis,
        )

    return StrainCard(
        name=str(d.get("strain", "strain")),
        M_base=_rate(d, "M_base"), M_on_ref=_rate(d, "M_on_ref"), M_off=_rate(d, "M_off"),
        mu_supply=_rate(d, "mu_supply"), background_false=_rate(d, "background_false"),
        mu_supply_source=(str(d["mu_supply_source"]) if d.get("mu_supply_source") else None),
        selection_marker=(str(d["selection_marker"]) if d.get("selection_marker") else None),
        tau_on_h=_rate(d, "tau_on_h"), tau_off_h=_rate(d, "tau_off_h"),
        F_doses=[{"dose_mM": float(fd["dose_mM"]), "F": float(fd["F"]),
                  "status": str(fd.get("status", "assumed")),
                  "note": str(fd.get("note", ""))}
                 for fd in (d.get("F_doses", []) or [])],
        ref_dose_mM=_scalar(d, "ref_dose_mM"),
        td_h_off=_scalar(d, "td_h_off"), td_h_on=_scalar(d, "td_h_on"),
        Nt_saturation=_scalar(d, "Nt_saturation"),
        clone_dispersion_sigma_ln=float(d.get("clone_dispersion_sigma_ln", 1.0)),  # 见 notes
        passenger_scale=_scalar(d, "passenger_scale"),
        inducer=inducer,
        genome_mb=(float(d["genome_mb"]) if d.get("genome_mb") is not None else None),
        plasmids=_load_plasmids(d),
        protected=_load_protected(d, path),
    )


def _card_notes(d: dict, card: StrainCard | None = None) -> list:
    """列出"用了默认值/未配置"的可选项——不阻断计算，但会影响结论。"""
    notes = []
    src = d.get("mu_supply_source")
    tgt = d.get("selection_marker")
    if src and tgt and str(src).strip().lower() != str(tgt).strip().lower():
        notes.append(f"μ_supply 来源标记 = {src}，本次选择标记 = {tgt}——属于外推；"
                     "按 references/markers.md 判断：单拷贝蛋白可给先验（assumed），"
                     "rRNA 多拷贝/酶灭活/外排不可外推，必须实测。")
    if d.get("background_false") is None:
        notes.append("未提供 background_false：不做背景假抗性红灯检查")
    if d.get("clone_dispersion_sigma_ln") is None:
        notes.append("未提供 clone_dispersion_sigma_ln：按 1.0 处理（会明显放大 MC 区间）")
    if d.get("genome_mb") is None or not d.get("protected"):
        notes.append("未配置 genome_mb/protected：P_protected_max / P_mutator_max 不会生效"
                     "（保护上限无法计算）")
    elif d.get("plasmids") and not any(str(r.get("role", "other")) == "mutator"
                                       for r in (d.get("protected") or [])):
        notes.append("配置了质粒但 protected 里没有 role: mutator——mutator 自失活上限"
                     "（P_mutator_max）不会生效；若 mutator 在质粒上，请补一条 role: mutator 的保护区")
    if not d.get("F_doses") and d.get("inducer") is None:
        notes.append("未配置 dose 轴（F_doses 或 inducer）：按恒定表达处理，剂量旋钮失效")
    if card is not None:
        from .model import card_warnings, params_without_ci
        notes.extend(card_warnings(card))
        missing = params_without_ci(card)
        if missing:
            notes.append(
                "以下参数没有 CI：" + "、".join(missing)
                + "——MC 对 RateParam 类参数按 **±2× 对数正态假设**（不是测量值），"
                  "inducer.* 与结构参数（td/Nt_saturation）按**常数**处理。"
                  "要给区间就在卡片里写 ci: [lo, hi]。")
    return notes


def load_protocol(path: str) -> tuple[list, Decision]:
    d = _load_yaml(path)
    phases = []
    for i, p in enumerate(d.get("phases", [])):
        t = p.get("type", "off")
        if isinstance(t, bool):
            raise CardError(
                f"{path} 第{i}个相位的 type 被解析成了布尔值（YAML 1.1 把裸 off/on 当 False/True）"
                "——请给 type 加引号，如 type: \"off\"")
        # YAML 1.1 把 1.0e8（e 后无符号）解析成字符串——统一数值化
        def _num(key):
            v = p.get(key)
            return float(v) if v is not None else None
        phases.append(Phase(name=str(p.get("name", f"phase{i}")),
                            type=str(t),
                            duration_h=float(p.get("duration_h", 0.0)),
                            generations=_num("generations"),
                            N0=_num("N0"), Nt=_num("Nt"),
                            dose_mM=_num("dose_mM"), td_h=_num("td_h"),
                            a_bar_override=_num("a_bar_override"),
                            plasmid_loss_per_gen=_num("plasmid_loss_per_gen") or 0.0))
    if not phases:
        raise CardError(f"{path}: 没有 phases")
    dec_d = d.get("decision", {})
    _require(dec_d, _REQUIRED_DECISION, f"{path}（协议卡的 decision 段）")
    dec = Decision(
        k_target_clones=int(dec_d["k_target_clones"]),
        P_success_min=float(dec_d["P_success_min"]),
        T_pulse_range_h=tuple(float(x) for x in dec_d["T_pulse_range_h"]),
        N0_min=float(dec_d["N0_min"]), N0_max=float(dec_d["N0_max"]),
        doses_available=[float(x) for x in dec_d.get("doses_available", [])],
        P_protected_max=float(dec_d.get("P_protected_max", 0.10)),
        P_mutator_max=float(dec_d.get("P_mutator_max", 0.05)),
        mc_samples=int(dec_d.get("mc_samples", 3000)),
        seed=int(dec_d.get("seed", 0)),
    )
    return phases, dec


def cmd_optimize(args) -> int:
    card = load_strain_card(args.strain, fill_defaults=getattr(args, "fill_defaults", False))
    if args.protocol:
        phases, dec = load_protocol(args.protocol)
    else:
        default_protocol = Path(__file__).resolve().parent.parent / "examples" / "protocol_card.example.yaml"
        print(f"注：未给 --protocol，使用内置默认协议卡（{default_protocol.name}）", file=sys.stderr)
        phases, dec = load_protocol(str(default_protocol))
    problems = card_problems(card, dec) + phase_problems(phases)
    if problems:
        # 不带病计算：负的 td、写反的 M_on/M_off 必须在这里拦住，
        # 否则会一路算成 events=0，再以"μ_supply 太低"的假诊断出现在报告里。
        raise CardError("卡片不合法（**不会带病计算**，先修这些）：\n  - "
                         + "\n  - ".join(problems))
    for n in _card_notes(_load_yaml(args.strain), card):
        print(f"警告：{n}", file=sys.stderr)
    if args.mc_samples:
        dec.mc_samples = args.mc_samples
    cands, diag = optimize(card, phases, dec)
    _meta = getattr(cands, "meta", {}) or {}
    if _meta.get("protection_binding"):
        _bits = []
        if _meta.get("protected_binding"):
            _bits.append("BGC/其他保护区")
        if _meta.get("mutator_binding"):
            _bits.append("mutator 自失活")
        print(f"警告：保护上限 binding（{'/'.join(_bits) or '保护'}）——k 要求的窗口超过上限；"
              "已按上限截断，"
              f"实际 E[events]={_meta.get('actual_events', float('nan')):.1f}"
              f"（缺口 {_meta.get('shortfall', float('nan')):.1f}）。", file=sys.stderr)
    md = render_optimize_report(card, phases, cands, diag, dec)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(md)
        print(f"报告已写入 {args.out}")
    else:
        print(md)
    return 0


def cmd_compare(args) -> int:
    rs1, nts1, _, notes1 = _read_counts_csv(args.d1)
    rs2, nts2, _, notes2 = _read_counts_csv(args.d2)
    for n in notes1 + notes2:
        print(f"注：{n}", file=sys.stderr)
    const_nt = args.nt1 is not None or args.nt2 is not None
    if const_nt:
        if args.nt1 is None or args.nt2 is None:
            print("错误：--nt1 与 --nt2 必须同时给（两组各自的每管细胞数，CFU/管）",
                  file=sys.stderr)
            return 1
        if args.nt1 <= 0 or args.nt2 <= 0:
            print("错误：--nt1 / --nt2 必须为正（CFU/管）", file=sys.stderr)
            return 1
        nts1 = [float(args.nt1)] * len(rs1)
        nts2 = [float(args.nt2)] * len(rs2)
    else:
        for path, r_list, n_list in ((args.d1, rs1, nts1), (args.d2, rs2, nts2)):
            if not n_list or len(n_list) != len(r_list):
                print(f"错误：{path} 需要 r,Nt 两列且逐管配对；"
                      "只有均值 Nt 时用 --nt1 <WT 的 Nt> --nt2 <mutator 的 Nt>",
                      file=sys.stderr)
                return 1
    lab1, lab2 = (args.labels.split(",") + ["D1", "D2"])[:2]
    r1 = lrt_equal_rates(rs1, nts1, rs2, nts2, fraction_plated=args.fraction_plated)
    fc = foldchange(rs1, nts1, rs2, nts2, fraction_plated=args.fraction_plated)
    print(f"== 速率比较：{lab1} vs {lab2} ==")
    print(f"μ({lab1}) = {r1['mu1']:.3e}   μ({lab2}) = {r1['mu2']:.3e}  （每细胞每分裂）")
    print(f"LRT 速率齐性: LLR={r1['llr']:.2f}, p = {r1['p']:.4g}"
          f"  → {'两批速率有显著差异' if r1['p'] < 0.05 else '无显著差异（H0 不拒绝）'}")
    lo, hi = fc["ci"]
    print(f"fold-change θ = μ({lab2})/μ({lab1}) = {fc['theta']:.2f}   "
          f"95%CI [{lo:.2f}, {hi:.2f}]")
    if const_nt:
        print(f"（Nt 口径：常数 {args.nt1:g} / {args.nt2:g} CFU/管——各管 Nt 有波动时"
              "这会系统性低估速率，能逐管记录就换成 r,Nt 两列）")
    print("（小样本下 profile CI 偏乐观：C<50/批时结论需重复实验支持）")
    print("提示：p>0.05 不等于「没有差异」——看 fold-change 的 95%CI 是否跨 1；"
          "两批还必须是同标记、同选择浓度、同批操作。", file=sys.stderr)
    return 0


def _warn_nonmonotone(history: list, lower_is_better: bool) -> None:
    """模拟型样本量搜索的诚实提示：CV/功效估计自带 10–20% 的 MC 误差，
    历史上不单调是常态——结论要结合表格看，别只信那一行"建议 C"。"""
    vals = [v for _, v in history]
    if len(vals) < 3:
        return
    ordered = sorted(range(len(vals)), key=lambda i: history[i][0])
    seq = [vals[i] for i in ordered]
    mono = all(b <= a + 1e-12 for a, b in zip(seq, seq[1:])) if lower_is_better \
        else all(b >= a - 1e-12 for a, b in zip(seq, seq[1:]))
    if not mono:
        print("⚠ 模拟值随 C 非单调（Monte-Carlo 误差，非 bug）：请按上表取'刚达标且留余量'的 C；"
              "要确定性答案用不加 --sim 的解析式", file=sys.stderr)


def cmd_power(args) -> int:
    from .power import (power_lrt_sim, sample_size_power_sim,
                        sample_size_precision_analytic, sample_size_precision_sim)
    from .power import sigma_ln_from_cv
    sigma_ln = sigma_ln_from_cv(args.nt_cv)
    if args.mode == "precision":
        from statistics import NormalDist
        if args.conf_level > 0:
            z = NormalDist().inv_cdf(0.5 + args.conf_level / 2.0)
            level = f"{args.conf_level:.0%} CI 半宽"
        else:                                        # 0 = 1σ 口径（文档约定）
            z, level = 1.0, "1σ 相对标准误"
        # ★ 只取整一次：解析式内部按给定的置信水平算，不再"先 ceil 1σ 再乘 z²"
        C = sample_size_precision_analytic(args.mu, args.nt_mean, args.target_rse,
                                           conf_level=args.conf_level or None)
        C1 = sample_size_precision_analytic(args.mu, args.nt_mean, args.target_rse,
                                            conf_level=None)
        C95 = sample_size_precision_analytic(args.mu, args.nt_mean, args.target_rse,
                                             conf_level=0.95)
        print("== 精度型样本量（常数 Nt 解析式）==")
        print(f"目标：{level} ≤ {args.target_rse:.0%} × μ̂（μ={args.mu:.2e}，z={z:.3f}）")
        print(f"C = {C} 管")
        print(f"（同一精度数字的另外两个口径：1σ → {C1} 管；95% CI 半宽 → {C95} 管。"
              "design-induction 报的是 **1σ** 那一套，别混用）")
        if args.sim:
            target_1sigma = args.target_rse / z     # 与解析同口径：CI 半宽 → 1σ CV
            out = sample_size_precision_sim(args.mu, args.nt_mean, target_1sigma,
                                            sigma_ln_nt=sigma_ln, n_sim=args.n_sim,
                                            c_max=args.c_max, seed=args.seed)
            print("\n-- 模拟实测（含 Nt 波动/部分涂板，同一 CI 半宽口径）--")
            for c, v in out["history"]:
                print(f"  C={c:4d}  CV(μ̂)={v:.3f}（95% 半宽 {1.96*v:.1%}）")
            print(f"建议 C = {out['C']}   （{out['note']}）")
            _warn_nonmonotone([(c, v) for c, v in out["history"]], lower_is_better=True)
    else:
        mu2 = args.mu1 * args.fold
        print("== 检出型功效分析（LRT 两批比较）==")
        print(f"H1: μ2 = {args.fold:g}×μ1（{args.mu1:.2e} → {mu2:.2e}），α={args.alpha}")
        if args.C is not None:
            p = power_lrt_sim(args.mu1, mu2, args.C, args.nt_mean, sigma_ln_nt=sigma_ln,
                              alpha=args.alpha, n_exp=args.n_exp, seed=args.seed)
            print(f"C={args.C} 管/批 → 功效 = {p:.1%}"
                  f"  {'（足够）' if p >= args.target_power else f'（不足，目标 {args.target_power:.0%}）'}")
        else:
            out = sample_size_power_sim(args.mu1, args.fold, args.target_power,
                                        args.nt_mean, sigma_ln_nt=sigma_ln, alpha=args.alpha,
                                        n_exp=args.n_exp, c_max=args.c_max, seed=args.seed)
            for c, p in out["history"]:
                print(f"  C={c:4d}  功效={p:.1%}")
            print(f"建议 C = {out['C']}   （{out['note']}）" if out["C"] else out["note"])
            _warn_nonmonotone([(c, p) for c, p in out["history"]], lower_is_better=False)
    return 0


def cmd_fit_fc(args) -> int:
    from .fitting import fit_hill_fc
    pairs = _read_xy_csv(args.input, "fit-fc（剂量 → M_on）")
    doses = [x for x, _ in pairs]
    m_on = [y for _, y in pairs]
    fit = fit_hill_fc(doses, m_on, m_base=args.mu_base, n_boot=args.n_boot, seed=args.seed)
    print("== F(c) Hill 拟合 ==")
    print(f"F_max = {fit.fmax:.1f}\u00d7   EC50 = {fit.ec50:.3f} mM   n = {fit.n_exp:.2f}   R\u00b2 = {fit.r2:.3f}")
    if fit.ci:
        for k, (lo, hi) in fit.ci.items():
            print(f"  {k}: 95%CI [{lo:.3g}, {hi:.3g}]")
    for w in fit.warnings:
        print(f"\u26a0 {w}")
    print("\n-- 剂量查询（粘回 strain_card.F_doses）--")
    for q in args.query:
        print(f"  {{dose_mM: {q}, F: {fit.predict(q):.0f}, status: fitted}}")
    return 0


def cmd_fit_tau(args) -> int:
    from .fitting import fit_tau_off
    pairs = _read_xy_csv(args.input, "fit-tauoff（时间 h → μ）")
    hours = [x for x, _ in pairs]
    mus = [y for _, y in pairs]
    fit = fit_tau_off(hours, mus, mu_off=args.mu_off, td_h=args.td,
                      m_on=args.m_on, n_boot=args.n_boot, seed=args.seed)
    print("== \u03c4_off 洗脱衰变拟合 ==")
    print(f"\u03c4_off = {fit.tau_off_h:.2f} h （半衰期 {fit.halflife_h:.2f} h）")
    if fit.ci:
        print(f"95%CI [{fit.ci[0]:.2f}, {fit.ci[1]:.2f}] h")
    print(f"\u03bc_on(拟合) = {fit.mu_on_fit:.4g} /基因组/代；\u03bc_off(固定) = {fit.mu_off:.4g}")
    if fit.tail_passengers is not None:
        print(f"尾部负荷（每基因组）：洗脱后若继续培养至活性衰尽 (\u222b\u2080^\u221e)"
              f" \u2248 {fit.tail_passengers:.3g}；只再培养 1 个 \u03c4 \u2248 {fit.tail_1tau:.3g}")
    for w in fit.warnings:
        print(f"\u26a0 {w}")
    return 0


def cmd_fit_induction(args) -> int:
    """诱导剂剂量 → 速率：机制模型（表达×竞争占位）拟合 + 饱和分析。"""
    from .induction import adequacy, fit_induction, hill, psi, psi_max, saturation_table

    pairs = _read_xy_csv(args.input, "fit-induction（剂量 mM → M）")
    if len(pairs) < 3:
        raise ValueError(f"{args.input}: 至少 3 个剂量点（推荐 ≥5），当前读到 {len(pairs)} 行")
    doses = [x for x, _ in pairs]
    mus = [y for _, y in pairs]

    fit = fit_induction(doses, mus, args.m_off, mixing=args.mixing,
                        n_boot=args.n_boot, seed=args.seed)
    b = fit.best
    print("== 诱导剂剂量响应（机制模型：表达 × 竞争占位）==")
    print(f"M_off = {fit.m_off:.4g}（固定）｜点数 = {len(doses)}")
    print("\n模型比较（ΔAIC<2 视为分不出；按简约原则取更简单的）：")
    print(f"  {'模型':10s} {'A':>8s} {'C50(mM)':>9s} {'n':>6s} {'q':>5s} {'Δsat':>8s} {'AIC':>8s} {'ΔAIC':>7s}")
    for f in fit.models:
        qs = f"{f.q:.2f}" if f.mixing == "two_core" else "—"
        mark = " ← 选中" if f is b else ""
        print(f"  {f.mixing:10s} {f.amplitude:8.3f} {f.ec50_mM:9.3f} {f.n_hill:6.2f} "
              f"{qs:>5s} {f.delta_sat:8.3f} {f.aic:8.1f} {f.aic - b.aic:+7.1f}{mark}")

    print(f"\n最佳：{b.mixing}｜饱和幅度 Δsat = {b.delta_sat:.4g}"
          f"（= A·ψ(1)，不是占比 p）｜C50 = {b.ec50_mM:.3g} mM｜n = {b.n_hill:.2f}")
    if fit.ci:
        ci = fit.ci
        bits = []
        if "delta_sat" in ci:
            bits.append(f"Δsat 95%CI [{ci['delta_sat'][0]:.3g}, {ci['delta_sat'][1]:.3g}]")
        if "ec50_mM" in ci and ci["ec50_mM"][1] > ci["ec50_mM"][0] * (1 + 1e-9):
            bits.append(f"C50 95%CI [{ci['ec50_mM'][0]:.3g}, {ci['ec50_mM'][1]:.3g}]")
        if bits:
            print("  " + "；".join(bits))

    aq = adequacy(fit, doses)
    print(f"\n饱和分析：C50 = {aq['c50_mM']:.3g} mM，95% 饱和剂量 = {aq['c95_mM']:.3g} mM；"
          f"最高剂量 {aq['max_dose']:.3g} mM → "
          + ("已饱和" if aq["reaches_95"] else "未饱和（Δsat 仍属外推）"))

    print("\n剂量表（供给/负荷 = 活性 ÷ M，越大越划算）：")
    print(f"  {'c(mM)':>8s} {'h':>6s} {'活性':>7s} {'M 预测':>9s} {'供给/负荷':>9s}")
    for row in saturation_table(fit, doses):
        print(f"  {row['dose_mM']:8.3g} {row['h']:6.3f} {row['activity']:7.3f} "
              f"{row['M_pred']:9.4g} {row['supply_per_load']:9.3g}")

    if args.query:
        print("\n查询：")
        for q in args.query:
            h = float(hill(q, b.ec50_mM, b.n_hill))
            rel = float(psi(h, b.mixing, b.q)) / psi_max(b.mixing, b.q)
            print(f"  c={q:g} mM → M={fit.predict(q):.4g}"
                  f"（{rel:.1%} 饱和；再加浓度最多再涨 {(1 - rel) * b.delta_sat:.3g}）")

    for w in fit.warnings:
        print(f"⚠ {w}")
    print("\n-- 粘进 strain_card --")
    print("inducer:")
    print(f"  name: {args.name}")
    print(f"  mixing: {b.mixing}")
    print(f"  amplitude: {b.amplitude:.6g}")
    print(f"  c50_mM: {b.ec50_mM:.6g}   # 速率半饱和（旧名 ec50_mM 也认）")
    print(f"  n_hill: {b.n_hill:.6g}")
    print(f"  q: {b.q:.4g}")
    return 0


def cmd_design_induction(args) -> int:
    """诱导剂剂量轴方案：用户只给"几个浓度"，其余自动决定。

    决策顺序：
      1) 用户给了具体剂量 → 直接用，并检查覆盖；
      2) 没给剂量、但知道 C50 → 在 C50 上下按对数布点；
      3) 没给剂量、也不知道 C50 → 在给定范围内宽间距布点（定位试点）。
    每点管数：知道 C50 与参考 m 时按该点预期 m 计算；否则用固定管数。

    2026-10-29 修：`--ref-m` 不再在缺 C50 时被静默忽略；`--strain` 可以直接从菌株卡读
    C50/n（不必先把卡填满）；覆盖检查只看 C50（不再要求 ref-m）。
    """
    import numpy as _np

    from .induction import hill, psi
    from .power import sample_size_precision_analytic

    c50 = args.c50
    nh = args.n_hill
    card_note = ""
    if args.strain:
        c50_card, nh_card, amp_card = _peek_inducer(args.strain)
        if c50_card is not None:
            if c50 is None:
                c50 = c50_card
            if nh_card:            # 卡里有 n 就用卡的（--n-hill 只是"预期值"）
                nh = nh_card
            card_note = (f"（--strain 读到 C50={c50_card:g} mM"
                         + (f"、n={nh_card:g}" if nh_card else "")
                         + (f"、amplitude={amp_card:g}" if amp_card else "") + "）")
        else:
            print(f"注：{args.strain} 里没有 inducer.c50_mM——按 C50 未知处理", file=sys.stderr)
    if args.ref_m is not None and c50 is None:
        print("⚠ 给了 --ref-m 但没有 C50（既没用 --c50，也没从 --strain 读到）："
              "无法把 m 折算到各剂量，**--ref-m 被忽略**，管数只能给固定值。"
              "补 --c50 <半饱和剂量> 或 --strain <菌株卡.yaml> 即可按点预测。", file=sys.stderr)

    prec = args.precision if args.precision is not None else (0.30 if c50 is None else 0.20)
    nt = args.nt
    ref = args.ref_dose if args.ref_dose else c50

    # ---- 1) 决定剂量点 ----
    if args.doses:
        doses = sorted(set(float(x) for x in args.doses))
        basis = "按你给的剂量"
    elif c50:
        n = args.n_doses or 5
        doses = [round(float(x), 3) for x in _np.geomspace(c50 / 5.0, c50 * 6.0, n)]
        basis = f"已知 C50≈{c50:g} mM → 膝点上下布点"
    else:
        lo, hi = args.range
        if not (0 < lo < hi):
            raise ValueError("--range 需要 0 < 下限 < 上限（mM）")
        n = args.n_doses or 4
        doses = [round(float(x), 4) for x in _np.geomspace(lo, hi, n)]
        basis = f"C50 未知 → {lo:g}–{hi:g} mM 宽间距定位"

    # ---- 2) 每点管数 ----
    predict = c50 is not None and args.ref_m is not None
    h_ref = psi_ref = None
    if predict:
        if c50 <= 0 or nh <= 0 or ref <= 0:
            raise ValueError("--c50 / --n-hill / --ref-dose 必须为正")
        h_ref = float(hill(ref, c50, nh))
        psi_ref = float(psi(h_ref, "linear", 1.0))
        if psi_ref <= 0:
            raise ValueError("参考剂量处活性为 0：检查 --ref-dose 是否远低于 C50")

    print("== 诱导剂剂量轴方案 ==")
    print(f"依据：{basis}{card_note}")
    if predict:
        print(f"假设：C50={c50:g} mM、n={nh:g}、参考 {ref:g} mM 处 m={args.ref_m:g}、"
              f"Nt={nt:.2g}；精度目标 1σ={prec:.0%}")
    elif c50 is not None:
        print(f"每点固定 {max(args.c_min * 2, 24)} 管（有 C50 但缺 --ref-m：无法预测各点 m）")
    else:
        print(f"每点固定 {max(args.c_min * 2, 24)} 管（C50 未知，无法预测各点 m）")
    print()

    total = 0
    below = above = 0
    if predict:
        print(f"{'剂量(mM)':>9s} {'活性':>6s} {'预期 m':>7s} {'C(1σ)':>6s}  说明")
    else:
        print(f"{'剂量(mM)':>9s}   说明")

    for i, c in enumerate(doses):
        # 覆盖计数只看 C50：即使算不出管数，"有没有跨膝点"也是已知的
        if c50:
            if c < c50:
                below += 1
            elif c > c50:
                above += 1
        if not predict:
            note = "最低点：通常 m 偏小甚至 <0.3" if i == 0 else (
                   "最高点：可能已饱和，注意计数过多" if i == len(doses) - 1 else "")
            print(f"{c:9.4g}   {note}")
            total += max(args.c_min * 2, 24)
            continue

        h = float(hill(c, c50, nh))
        act = float(psi(h, "linear", 1.0)) / psi_ref
        m = args.ref_m * act
        if m < 0.3:
            print(f"{c:9.3g} {act:6.3f} {m:7.2f} {'—':>6s}  ⚠ m<0.3：任何估计器都不可靠，删掉或提高该点 Nt")
            continue
        c1 = sample_size_precision_analytic(m / nt, nt, prec, conf_level=None)
        n_use = max(args.c_min, c1)
        total += n_use
        note = f"取 {n_use} 管" if n_use != c1 else ""
        print(f"{c:9.3g} {act:6.3f} {m:7.2f} {c1:6d}  {note}")

    print()
    if c50:
        print(f"覆盖检查：C50 以下 {below} 点、以上 {above} 点（建议各 ≥2）")
        if below < 2 or above < 2:
            print("  ⚠ 覆盖不足：fit-induction 会因未跨膝点而给不出可信的 C50/饱和幅度")
    print(f"合计（不含对照）：约 {total} 管")
    print(f"对照建议：0 mM 未诱导 ×{args.c_min}（测 M_off）、WT ×{args.c_min}（确认诱导剂不影响 WT）")
    print()
    print(f"口径：上表 C 是 **1σ 相对标准误 ≤{prec:.0%}**；若你要的是 95% CI 半宽 ≤{prec:.0%}，"
          "把每个数字 ×3.84（`windowcalc power` 默认用后一种口径）。")
    print()
    if predict:
        print("下一步：测完 → windowcalc fit-induction --input 数据.csv --m-off <M_off>")
        print("        → windowcalc optimize --strain 菌株卡.yaml --protocol 协议卡.yaml（自动选浓度）")
    else:
        print("下一步：测完 → windowcalc fit-induction --input 数据.csv --m-off <M_off>")
        print("        得到粗 C50 后 → windowcalc design-induction --c50 <该值> --ref-m <参考点 m>")
    return 0


def _peek_inducer(path: str):
    """只读 inducer 块的 c50/n/amplitude（不要求菌株卡完整——还没实测时也能用）。"""
    d = _load_yaml(path)
    ind = d.get("inducer") or {}

    def _g(key):
        v = ind.get(key)
        if isinstance(v, dict):
            v = v.get("value")
        return float(v) if v is not None else None

    return (_g("c50_mM") if ind.get("c50_mM") is not None else _g("ec50_mM")), _g("n_hill"), _g("amplitude")


def build_parser() -> argparse.ArgumentParser:
    """构建 CLI 解析器（测试与文档数字自检要数子命令，所以单独暴露）。"""
    ap = argparse.ArgumentParser(
        prog="windowcalc",
        description="把「每管几个抗性菌落」变成「突变率」；mutator 线还能算表达窗口（开多久/多强）。",
        epilog="新手先跑：windowcalc rate --input 计数.csv --nts <每管 Nt>\n"
               "教程：docs/tutorial.md ｜ 术语：docs/glossary.md ｜ 自检：python3 selftest.py",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    # required=False：裸跑时由 main() 打一段中文引导，而不是 argparse 的英文报错
    sub = ap.add_subparsers(dest="cmd", required=False)

    ap_rate = sub.add_parser("rate", help="原始计数 → 突变率（最常用；有 r,Nt 两列就能算）")
    ap_rate.add_argument("--input", required=True, help="CSV：r[,Nt] 每行一管")
    ap_rate.add_argument("--label", default="group",
                         help="这组数据的名字（只用于打印/注释，不是卡片字段名）")
    ap_rate.add_argument("--field", default=None, type=_field_arg,
                         metavar="{M_base,M_on_ref,M_off,mu_supply,background_false}",
                         help="这段结果要填进菌株卡的哪个字段——给了就把片段键名写成它"
                              "（不给则用 --label，可能粘不进卡片）")
    ap_rate.add_argument("--fraction-plated", type=float, default=1.0,
                         help="每管涂板比例 ζ（整管=1.0）；ζ<1 用二项稀释似然，"
                              "此时 Nt 必须填整管总 CFU（不是涂板那部分的）")
    ap_rate.add_argument("--nts", type=float, default=None,
                         help="各管 Nt 均值（CSV 无 Nt 列时）；单位默认 CFU/管")
    ap_rate.add_argument("--nt-unit", choices=["per_tube", "per_mL"], default=None,
                         help="--nts 的单位：per_tube=CFU/管；per_mL=CFU/mL（需 --volume-mL）。"
                              "不给会按 CFU/管 解释并打印醒目提示（μ 差一个体积倍是最常见的错）")
    ap_rate.add_argument("--volume-mL", type=float, default=None,
                         help="培养物体积（mL），配合 --nt-unit per_mL 把 CFU/mL 换成 CFU/管")
    ap_rate.add_argument("--nt-mode", choices=["mean", "per_tube"], default="mean",
                         help="per_tube=Ycart 2014 逐管 Nt 校正（需 CSV 带 Nt 列）")
    ap_rate.set_defaults(func=cmd_rate)

    ap_opt = sub.add_parser("optimize", help="三旋钮窗口优化 → Markdown 报告")
    ap_opt.add_argument("--strain", required=True, help="strain_card YAML")
    ap_opt.add_argument("--protocol", default=None,
                        help="protocol_card YAML；缺省用技能自带的默认协议卡（只问 k/成功率/上限即可）")
    ap_opt.add_argument("--fill-defaults", action="store_true",
                        help="菌株卡缺的必填项用默认假设值填充（报告标 ★未测）；有实测值请替换")
    ap_opt.add_argument("--out", default=None, help="输出 Markdown 路径（缺省打印）")
    ap_opt.add_argument("--mc-samples", type=int, default=None)
    ap_opt.set_defaults(func=cmd_optimize)

    ap_cmp = sub.add_parser("compare", help="两批波动实验的速率比较（LRT + fold-change）")
    ap_cmp.add_argument("--d1", required=True, help="第一批 CSV（r,Nt）")
    ap_cmp.add_argument("--d2", required=True, help="第二批 CSV（r,Nt）")
    ap_cmp.add_argument("--labels", default="D1,D2", help="两批的名称，逗号分隔")
    ap_cmp.add_argument("--fraction-plated", type=float, default=1.0,
                        help="每管涂板比例 ζ（整管=1.0）；ζ<1 用二项稀释似然，Nt 填整管总 CFU")
    ap_cmp.add_argument("--nt1", type=float, default=None,
                        help="只有均值 Nt 时：第一批的每管细胞数（CFU/管）；与 --nt2 同时给")
    ap_cmp.add_argument("--nt2", type=float, default=None,
                        help="只有均值 Nt 时：第二批的每管细胞数（CFU/管）；与 --nt1 同时给")
    ap_cmp.set_defaults(func=cmd_compare)

    ap_pw = sub.add_parser("power", help="样本量/功效分析（解析 + 模拟）")
    ap_pw.add_argument("--mode", choices=["precision", "power"], default="precision",
                       help="precision=估速精度；power=两批比较的检出功效")
    ap_pw.add_argument("--mu", type=float, default=2e-9, help="precision 模式的 μ")
    ap_pw.add_argument("--mu1", type=float, default=2e-9, help="power 模式的基线 μ1")
    ap_pw.add_argument("--fold", type=float, default=2.0, help="power 模式要检出的倍数")
    ap_pw.add_argument("--nt-mean", type=float, default=3e8, help="终末细胞数均值")
    ap_pw.add_argument("--nt-cv", type=float, default=0.0,
                       help="Nt 的真 CV=std/mean（0=恒定）；内部换算 σ_ln=√ln(1+CV²)")
    ap_pw.add_argument("--target-rse", type=float, default=0.25, help="precision：目标 RSE")
    ap_pw.add_argument("--conf-level", type=float, default=0.95,
                       help="precision：置信水平 0–1；传 0 表示按 1σ 相对标准误")
    ap_pw.add_argument("--target-power", type=float, default=0.9, help="power：目标功效")
    ap_pw.add_argument("--alpha", type=float, default=0.05)
    ap_pw.add_argument("--C", type=int, default=None, help="power：指定 C 算功效（缺省搜索最小 C）")
    ap_pw.add_argument("--sim", action="store_true", help="precision：追加模拟实测")
    ap_pw.add_argument("--n-sim", type=int, default=150)
    ap_pw.add_argument("--n-exp", type=int, default=60)
    ap_pw.add_argument("--c-max", type=int, default=256)
    ap_pw.add_argument("--seed", type=int, default=0)
    ap_pw.set_defaults(func=cmd_power)
    ap_fc = sub.add_parser("fit-fc", help="F(c) 剂量响应 Hill 拟合（参数 7）")
    ap_fc.add_argument("--input", required=True, help="CSV：dose_mM,m_on（每行一个剂量）")
    ap_fc.add_argument("--mu-base", dest="mu_base", type=float, required=True,
                       help="M_base（F = m_on/M_base）")
    ap_fc.add_argument("--query", type=float, nargs="*", default=[],
                       help="查询这些剂量的 F 预测")
    ap_fc.add_argument("--n-boot", type=int, default=200)
    ap_fc.add_argument("--seed", type=int, default=0)
    ap_fc.set_defaults(func=cmd_fit_fc)

    ap_di = sub.add_parser("design-induction",
                           help="诱导剂剂量轴设计：剂量点 + 每点管数（按预期 C50/n）")
    ap_di.add_argument("--c50", type=float, default=None,
                       help="预期半饱和剂量（mM）；不知道就改用 --pilot")
    ap_di.add_argument("--strain", default=None,
                       help="从菌株卡读 inducer 的 C50/n（卡不必填满；也可显式给 --c50）")
    ap_di.add_argument("--pilot", action="store_true",
                       help="（兼容旧写法，等价于不给 --c50）")
    ap_di.add_argument("--n-doses", type=int, default=None,
                       help="剂量点个数（默认：已知 C50 用 5，未知用 4）")
    ap_di.add_argument("--range", type=float, nargs=2, default=[0.01, 10.0],
                       help="C50 未知时的剂量范围（mM），默认 0.01–10")
    ap_di.add_argument("--n-hill", type=float, default=1.5, help="预期 Hill 系数（默认 1.5）")
    ap_di.add_argument("--ref-dose", type=float, default=None,
                       help="参考剂量（默认取 C50）")
    ap_di.add_argument("--ref-m", type=float, default=None, help="参考剂量处的预期 m（正式模式必填）")
    ap_di.add_argument("--nt", type=float, default=3e8, help="终末细胞数（CFU/管）")
    ap_di.add_argument("--precision", type=float, default=None, help="目标 1σ 相对精度（默认 0.20）")
    ap_di.add_argument("--doses", type=float, nargs="*", default=None,
                       help="显式给剂量点（缺省按 C50 自动生成 5 点）")
    ap_di.add_argument("--c-min", type=int, default=12, help="每点最少管数（默认 12）")
    ap_di.set_defaults(func=cmd_design_induction)


    ap_ind = sub.add_parser("fit-induction",
                            help="诱导剂剂量→速率机制拟合（表达×竞争占位，含饱和分析）")
    ap_ind.add_argument("--input", required=True, help="CSV：dose_mM,M（每行一个剂量）")
    ap_ind.add_argument("--m-off", type=float, required=True, help="off 态速率（M_off，同单位）")
    ap_ind.add_argument("--mixing", choices=["auto", "linear", "two_core", "both_core"],
                        default="auto", help="混合规则；auto=AIC 三选一（默认）")
    ap_ind.add_argument("--name", default="诱导剂", help="诱导剂名称（写进 YAML 片段）")
    ap_ind.add_argument("--query", type=float, nargs="*", default=[],
                        help="查询这些剂量的 M 预测")
    ap_ind.add_argument("--n-boot", type=int, default=200)
    ap_ind.add_argument("--seed", type=int, default=0)
    ap_ind.set_defaults(func=cmd_fit_induction)

    ap_tau = sub.add_parser("fit-tauoff", help="τ_off 洗脱衰变拟合（参数 5）")
    ap_tau.add_argument("--input", required=True, help="CSV：h,mu（每行一个时间点）")
    ap_tau.add_argument("--mu-off", type=float, required=True, help="off 态速率（M_off）")
    ap_tau.add_argument("--td", type=float, default=2.0, help="倍增时间 h（尾部负荷折算）")
    ap_tau.add_argument("--m-on", type=float, default=None, help="诱导态速率（尾部负荷用）")
    ap_tau.add_argument("--n-boot", type=int, default=200)
    ap_tau.add_argument("--seed", type=int, default=0)
    ap_tau.set_defaults(func=cmd_fit_tau)

    return ap


NO_CMD_HELP = """你想做什么？

  ① 我有每管的抗性菌落数（r）和每管细胞数（Nt），要算突变率
       windowcalc rate --input 计数.csv --nts 3e8
  ② 我想知道这两批（WT / mutator）速率有没有差别
       windowcalc compare --d1 WT.csv --d2 M7.csv --nt1 3e8 --nt2 3e8
  ③ 我想知道要铺多少管
       windowcalc power --mode precision --mu 2e-9 --nt-mean 3e8
  ④ 我想设计诱导剂浓度 / 算 mutator 表达窗口
       windowcalc design-induction --n-doses 5     →   windowcalc optimize --strain 菌株卡.yaml --protocol 协议卡.yaml

第一次用先读：docs/tutorial.md（10 分钟）｜术语：docs/glossary.md｜自检：python3 selftest.py
每个子命令的帮助：windowcalc <子命令> --help"""


def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    if not getattr(args, "cmd", None):
        print(NO_CMD_HELP)
        return 0
    try:
        return args.func(args)
    except KeyboardInterrupt:
        return 130
    except (ValueError, KeyError, TypeError, ZeroDivisionError, OverflowError,
            FileNotFoundError, IsADirectoryError, PermissionError) as e:
        # 面向非程序员：数据/卡片写错时给一行中文错误，不要甩 traceback。
        # （此前只捕 ValueError/FileNotFoundError：坏 YAML、protected 缺 kb、
        #   td 写成 {value, ci} 都是原始 KeyError/TypeError 栈）
        if getattr(args, "debug", False):
            raise
        print(f"错误：{e}", file=sys.stderr)
        if isinstance(e, CardError):
            print("提示：卡片写法见 docs/mutator-window.md §4；加 --debug 可看完整 traceback。",
                  file=sys.stderr)
        else:
            print("提示：数据格式见 docs/tutorial.md（CSV 要两列 r,Nt）；加 --debug 可看完整 traceback。",
                  file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
