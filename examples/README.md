# examples/ 说明

| 文件 | 是什么 | 能用吗 |
|---|---|---|
| `raw_counts_template.csv` | 计数表模板（每行一管 r,Nt） | 空模板，照填 |
| `record-template.csv` | 实验记录表模板（含组名/批次注释） | 空模板 |
| `strain_card.example.yaml` | 菌株卡示例 | ⚠️ **数值是占位假设**（`status: assumed`），直接 optimize 只会得到"★未测"的占位结论；实测后逐项替换 |
| `protocol_card.example.yaml` | 协议卡示例（相位表 + 决策参数） | 相位结构可直接改；`decision` 里的 k/N₀ 范围按实验调整 |
| `protocol_card.dual_switch.example.yaml` | **双开关**构型（诱导启动子 + 温敏质粒） | 含切换/治愈相位；`plasmid_loss_per_gen` 按实测改 |
| `strain_card.inducer.example.yaml` | 机制模型版菌株卡（`inducer` 块） | 与 `M_on_ref` 自洽，可直接跑 |
| `acceptance/RF2000.csv` | Rosche & Foster 2000 Table 3（60 管，完整数据；MSS m=1.23；中位数法 1.11） | ✅ 可直接跑，用来验证安装 |
| `acceptance/gbk_rf2000.csv` | 同一份数据转存为 **GBK** 编码（测"Excel 另存的 CSV"这条路径） | ✅ 仅用于编码回归 |
| `acceptance/plating_z02.csv` | ζ=0.2 部分涂板模拟数据（真值 m=3.0） | ✅ `--nt-mode per_tube --fraction-plated 0.2` |
| `acceptance/groupA.csv` / `groupB.csv` | 两批数据（速率差 10 倍） | ✅ `compare` 用 |
| `acceptance/LD1943_exp16.csv` | Luria & Delbrück 1943 第 16 次实验（20 管，ε=0.4，含 jackpot） | ✅ 真实数据；`--fraction-plated 0.4` |
| `acceptance/Demerec1945.csv` | Demerec 1945（30 管，m≈10.84） | ✅ 真实数据 |
| `acceptance/markers_demo.csv` · `dose_induction.csv` | 合成演示（文件头写明生成参数） | ✅ 只演示列格式 |

**每个数据文件的文件头都写了出处**（`# 文献数据集：…` 或 `# 合成数据（非实测、非默认）：…`）。
这是 `docs/generic-policy.md` 规则二的落地要求：**不许让读者把别人的数当成通用推荐值**。
阅读器会跳过以 `#` 开头的注释行（CSV 与账本都支持）。

`acceptance/` 里的数据由 `selftest.py` 使用；改引擎后若这些验收失败，先确认是有意改动。
