"""CLI 输入层测试：CSV 读取。

背景：旧实现用 csv.Sniffer().has_header()，对单列文件（只数菌落、没测 Nt——
最常见的输入）直接抛 `_csv.Error: Could not determine delimiter`。
"""

import pytest

from windowcalc.cli import _read_counts_csv


def test_single_column_without_header(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("0\n1\n0\n2\n1\n", encoding="utf-8")
    rs, nts, vols, notes = _read_counts_csv(str(p))
    assert rs == [0, 1, 0, 2, 1]
    assert nts is None


def test_header_with_paired_nt(tmp_path):
    p = tmp_path / "rnt.csv"
    p.write_text("r,Nt\n0,3e8\n1,2.9e8\n2,3.1e8\n", encoding="utf-8")
    rs, nts, vols, notes = _read_counts_csv(str(p))
    assert rs == [0, 1, 2]
    assert nts == pytest.approx([3e8, 2.9e8, 3.1e8])


def test_blank_lines_and_bom(tmp_path):
    p = tmp_path / "bom.csv"
    p.write_text("\ufeffr,Nt\n\n0,3e8\n\n1,3e8\n", encoding="utf-8")
    rs, nts, vols, notes = _read_counts_csv(str(p))
    assert rs == [0, 1] and len(nts) == 2


@pytest.mark.parametrize("body,msg", [
    ("r\n0\n-1\n", "非负整数"),
    ("r\n0\n1.5\n", "非负整数"),
    ("r\n0\nabc\n", "不是数字"),
    ("r,Nt\n0,\n1,3e8\n", "逐管配对"),
    ("r,Nt\n0,0\n", "必须为正"),
    ("r\n", "没有读到任何计数行"),
])
def test_rejects_bad_input(tmp_path, body, msg):
    p = tmp_path / "bad.csv"
    p.write_text(body, encoding="utf-8")
    with pytest.raises(ValueError, match=msg):
        _read_counts_csv(str(p))


def test_gbk_encoded_csv(tmp_path):
    """中文 Windows 上 Excel「另存为 CSV」是本地代码页（GBK）——必须能读。

    按 UTF-8 硬读会 UnicodeDecodeError；_read_text 回退 GB18030。
    """
    p = tmp_path / "gbk.csv"
    p.write_bytes("r,Nt\n0,3e8\n1,3e8\n0,2.9e8\n".encode("gb18030"))
    rs, nts, vols, notes = _read_counts_csv(str(p))
    assert rs == [0, 1, 0]
    assert nts == pytest.approx([3e8, 3e8, 2.9e8])


def test_undecodable_file_message(tmp_path):
    """两种编码都失败时给可执行的中文提示，而不是裸 UnicodeDecodeError。"""
    p = tmp_path / "bad.csv"
    p.write_bytes(b"r,Nt\n\xff\xfe\x00,3e8\n")
    with pytest.raises(ValueError, match="另存为 UTF-8"):
        _read_counts_csv(str(p))


# ---------------------------------------------------------------- Nt 单位（v0.4.1）
def test_volume_column_converts_cfu_per_ml(tmp_path):
    """第三列 volume_mL：Nt 按 CFU/mL 解释，自动乘体积。"""
    p = tmp_path / "vol.csv"
    p.write_text("r,Nt,volume_mL\n0,2.0e8,2\n1,2.0e8,2\n", encoding="utf-8")
    rs, nts, vols, notes = _read_counts_csv(str(p))
    assert rs == [0, 1]
    assert nts == pytest.approx([4.0e8, 4.0e8])
    assert vols == [2.0, 2.0]
    assert any("换算" in n for n in notes)


def test_volume_column_must_be_complete(tmp_path):
    p = tmp_path / "half.csv"
    p.write_text("r,Nt,volume_mL\n0,2e8,2\n1,2e8,\n", encoding="utf-8")
    with pytest.raises(ValueError, match="volume_mL"):
        _read_counts_csv(str(p))


def test_swapped_columns_hint(tmp_path):
    """列序反了（Nt,r）：r 会大到不可能，错误里应提示"像 Nt"。"""
    p = tmp_path / "swap.csv"
    p.write_text("r,Nt\n2.0e8,0\n2.0e8,1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="像 Nt"):
        _read_counts_csv(str(p))


def test_zero_nt_hint(tmp_path):
    """第二列是 0：提示可能是频率或列反了。"""
    p = tmp_path / "zero.csv"
    p.write_text("r,Nt\n0,0\n1,0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="频率|列顺序"):
        _read_counts_csv(str(p))


def test_frequency_like_nt_note(tmp_path):
    """Nt 写成频率（<1）时提示，但不阻断。"""
    p = tmp_path / "freq.csv"
    p.write_text("r,Nt\n0,1.0e-8\n1,2.0e-8\n", encoding="utf-8")
    _, _, _, notes = _read_counts_csv(str(p))
    assert any("频率" in n for n in notes)


def test_implausible_nt_range_note(tmp_path):
    """Nt 明显超出常见区间时提示（不阻断）。"""
    p = tmp_path / "big.csv"
    p.write_text("r,Nt\n0,5.0e-9\n1,2.0e13\n", encoding="utf-8")
    _, _, _, notes = _read_counts_csv(str(p))
    assert any(("CFU/mL" in n) or ("频率" in n) for n in notes)


# ---------------------------------------------------------------- 必填参数强制（v0.5.1）
def test_missing_required_strain_params_lists_all(tmp_path):
    """缺必填参数时报错并逐条列出，不静默用默认值。"""
    p = tmp_path / "card.yaml"
    p.write_text("strain: t\nM_base: {value: 0.004}\n", encoding="utf-8")
    from windowcalc.cli import load_strain_card
    with pytest.raises(ValueError) as ei:
        load_strain_card(str(p))
    msg = str(ei.value)
    for key in ("M_on_ref", "M_off", "mu_supply", "tau_on_h", "tau_off_h",
                "td_h_on", "td_h_off", "Nt_saturation", "ref_dose_mM"):
        assert key in msg
    assert "不会用默认值代替" in msg


def test_missing_required_decision_lists_all(tmp_path):
    p = tmp_path / "proto.yaml"
    p.write_text('protocol: x\nphases:\n  - {name: p, type: "pulse", duration_h: 6, N0: 1e7}\n'
                 'decision:\n  k_target_clones: 12\n', encoding="utf-8")
    from windowcalc.cli import load_protocol
    with pytest.raises(ValueError) as ei:
        load_protocol(str(p))
    msg = str(ei.value)
    for key in ("P_success_min", "T_pulse_range_h", "N0_min", "N0_max"):
        assert key in msg


def test_optional_defaults_are_announced(tmp_path):
    """可选项缺省时要给出警告清单，而不是静默。"""
    p = tmp_path / "card.yaml"
    p.write_text(
        "strain: t\nM_base: {value: 1e-3}\nM_on_ref: {value: 4e-3}\nM_off: {value: 1e-4}\n"
        "mu_supply: {value: 2e-8}\ntau_on_h: {value: 2}\ntau_off_h: {value: 4}\n"
        "td_h_on: 2.5\ntd_h_off: 2.0\nNt_saturation: 3e8\nref_dose_mM: 0.3\n",
        encoding="utf-8")
    import yaml
    from windowcalc.cli import _card_notes
    notes = _card_notes(yaml.safe_load(p.read_text(encoding="utf-8")))
    assert any("background_false" in n for n in notes)
    assert any("genome_mb" in n for n in notes)
    assert any("dose 轴" in n for n in notes)


# ---------------------------------------------------------------- Excel 输入
def test_xlsx_reads_same_as_csv(tmp_path):
    """新手直接给 .xlsx：必须与等价 CSV 得到同一组 (r, Nt)。"""
    pytest.importorskip("openpyxl")
    from openpyxl import Workbook

    wb = Workbook()
    ws = wb.active
    ws.title = "counts"
    ws.append(["r", "Nt"])                    # 中文/英文表头都要能跳过
    for r, nt in [(0, 3e8), (1, 3e8), (2, 2.9e8)]:
        ws.append([r, nt])
    p = tmp_path / "counts.xlsx"
    wb.save(p)

    rs, nts, vols, notes = _read_counts_csv(str(p))
    assert rs == [0, 1, 2]
    assert nts == pytest.approx([3e8, 3e8, 2.9e8])
    assert vols is None


def test_xlsx_with_text_header_and_blank_rows(tmp_path):
    """Excel 里常见：表头带文字、末尾空行。不能因此报错。"""
    pytest.importorskip("openpyxl")
    from openpyxl import Workbook

    wb = Workbook()
    ws = wb.active
    ws.append(["菌落数 r", "细胞数 Nt"])       # 中文表头
    ws.append([0, 3e8])
    ws.append([])                              # 空行
    ws.append([5, 3.1e8])
    p = tmp_path / "cn.xlsx"
    wb.save(p)

    rs, nts, _, _ = _read_counts_csv(str(p))
    assert rs == [0, 5]
    assert nts == pytest.approx([3e8, 3.1e8])


def test_old_xls_rejected_with_hint(tmp_path):
    """.xls（97-2003）读不了 → 必须点明「另存为」，而不是甩二进制解码错。"""
    p = tmp_path / "old.xls"
    p.write_bytes(b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1")   # OLE 头
    with pytest.raises(ValueError, match="另存为"):
        _read_counts_csv(str(p))


def test_acceptance_xlsx_matches_csv():
    """随包分发的验收 xlsx 必须与同名 CSV 一致（防打包漏文件、读错 sheet）。"""
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent / "examples" / "acceptance"
    a_rs, a_nts, _, _ = _read_counts_csv(str(root / "RF2000.csv"))
    b_rs, b_nts, _, _ = _read_counts_csv(str(root / "RF2000.xlsx"))
    assert a_rs == b_rs
    assert a_nts == pytest.approx(b_nts)



def test_fill_defaults_and_bundled_protocol(tmp_path):
    """无感路径：--fill-defaults 只补「与体系无关」的量（通用性政策规则二）。

    2026-11-15 修正：旧版会连 mu_supply / tau / td / ref_dose 一起编默认值——
    那正是政策禁止的"把别人的实验当通用推荐"。现在这六项缺就报错。
    """
    from pathlib import Path
    from windowcalc.cli import CardError, load_protocol, load_strain_card
    p = tmp_path / "min.yaml"
    p.write_text("strain: t\n", encoding="utf-8")
    with pytest.raises(CardError):
        load_strain_card(str(p))

    # 只补约定/文献量级的三项，测得的物理量不给编 → 仍然报错
    with pytest.raises(CardError) as ei:
        load_strain_card(str(p), fill_defaults=True)
    msg = str(ei.value)
    assert "mu_supply" in msg and "td_h_on" in msg, "应点名缺哪些实测项"

    # 把六项实测/物理量补齐后，--fill-defaults 才能跑
    q = tmp_path / "ok.yaml"
    q.write_text(
        "strain: t\n"
        "mu_supply: {value: 2e-8, status: assumed}\n"
        "tau_on_h: {value: 2.0, status: assumed}\n"
        "tau_off_h: {value: 4.0, status: assumed}\n"
        "td_h_on: {value: 2.5, status: assumed}\n"
        "td_h_off: {value: 2.0, status: assumed}\n"
        "ref_dose_mM: {value: 0.3, status: assumed}\n",
        encoding="utf-8")
    card = load_strain_card(str(q), fill_defaults=True)
    assert card.M_base.value == pytest.approx(0.004)
    assert card.M_base.status == "assumed"
    assert card.mu_supply.value == pytest.approx(2e-8)   # 来自卡片，不是内置默认
    bundled = Path(__file__).resolve().parent.parent / "examples" / "protocol_card.example.yaml"
    phases, dec = load_protocol(str(bundled))
    assert phases and dec.k_target_clones >= 1
    assert dec.P_protected_max > 0 and dec.P_mutator_max > 0


def test_marker_extrapolation_warning(tmp_path):
    """来源标记 ≠ 选择标记时，_card_notes 必须提示外推（不依赖用户主动问）。"""
    from windowcalc.cli import _card_notes, _load_yaml, load_strain_card
    p = tmp_path / "c.yaml"
    p.write_text(
        "strain: t\n"
        "M_base: 0.004\nM_on_ref: 0.4\nM_off: 0.012\n"
        "mu_supply: 2e-8\n"
        "mu_supply_source: tetracycline\nselection_marker: streptomycin\n"
        "tau_on_h: 2\ntau_off_h: 4\ntd_h_on: 2.5\ntd_h_off: 2\n"
        "Nt_saturation: 3e8\nref_dose_mM: 0.3\n", encoding="utf-8")
    card = load_strain_card(str(p))
    notes = _card_notes(_load_yaml(str(p)), card)
    assert any("外推" in n for n in notes)
    # 相同则不应提示
    p2 = tmp_path / "same.yaml"
    p2.write_text(p.read_text(encoding="utf-8").replace(
        "mu_supply_source: tetracycline\nselection_marker: streptomycin",
        "mu_supply_source: streptomycin\nselection_marker: streptomycin"), encoding="utf-8")
    notes2 = _card_notes(_load_yaml(str(p2)), load_strain_card(str(p2)))
    assert not any("外推" in n for n in notes2)


def test_plasmid_without_mutator_role_warns(tmp_path):
    """有质粒但 protected 没标 role: mutator → 自失活上限不会生效，必须提示。"""
    from windowcalc.cli import _card_notes, _load_yaml, load_strain_card
    p = tmp_path / "c.yaml"
    p.write_text(
        "strain: t\nM_base: 0.004\nM_on_ref: 0.4\nM_off: 0.012\nmu_supply: 2e-8\n"
        "tau_on_h: 2\ntau_off_h: 4\ntd_h_on: 2.5\ntd_h_off: 2\n"
        "Nt_saturation: 3e8\nref_dose_mM: 0.3\n"
        "genome_mb: 8.4\nplasmids: [{kb: 16, copies: 20}]\n"
        "protected: [{name: BGC, kb: 100, copies: 1, selection: silent, role: bgc}]\n",
        encoding="utf-8")
    notes = _card_notes(_load_yaml(str(p)), load_strain_card(str(p)))
    assert any("role: mutator" in n for n in notes)
