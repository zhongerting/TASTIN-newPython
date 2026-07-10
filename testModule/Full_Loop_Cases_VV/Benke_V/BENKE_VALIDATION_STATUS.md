# Benke 热工水力验证状态报告

## 当前状态

Benke_V 已建立首版 `benke_thermal_network_v1`，并能对 Benke 典型保温工况进行热网络计算和文献范围校核。

本地 `TOPAZII_VV_public_figures` 中当前没有 Benke 轴向套筒温度图、水出口温度图或 regulated He gap 有效导热率图的裁图；已有资料只在 `TOPAZII_VV_public_experimental_data.md` 中整理了典型范围和验证对象。因此当前验证状态定义为：

```text
partial_missing_digitized_data
```

这不是模型失败，而是实验曲线数据尚未数字化。模型不得用 Venable 电输出结果反推 Benke 热边界。

## 已完成的验证

当前 runner 会检查：

- `P_az = 0.88 P_TISA` 是否落在典型 3003 W 附近。
- `regulated_he_effective_k_w_m_k` 是否位于 `0.073-0.087 W/(m K)`。
- `water_h_w_m2_k` 是否位于 `528-1012 W/(m2 K)`。
- 水侧能量闭合误差是否足够小。

## 待补充的实验数据

将 Benke 文献图像数字化后，放入 `experimental_data/`：

- `benke_sleeve_thermocouple_12pt_digitized.csv`
- `benke_water_balance_digitized.csv`

放入后，`run_benke_thermal_validation.py` 会自动计算套筒 12 点温度 MAE/RMSE，以及水出口温度或水温升误差。

## 下一步建议

1. 找到或补充 Benke 原始 PDF/图像页，优先数字化 12 点套筒轴向温度和水侧热平衡数据。
2. 用默认参数跑一次对比，记录 MAE/RMSE。
3. 只在 Benke_V 内做全局参数扫描：`water_h`、`water_mass_flow`、`regulated_he_effective_k`、`extra_resistance`。
4. 得到 Benke 热边界先验后，再传递给 Venable_V 做 I-V 和最大功率验证。

## 2026-07-02 文献范围参数包络

新增 `benke_parameter_scan.py`，在 Benke 文献整理范围内扫描：

- `regulated_he_effective_k_w_m_k = 0.073, 0.08, 0.087`
- `water_h_w_m2_k = 528, 800, 1012`

运行目录：`runs/20260702_benke_literature_range_envelope/`

结果摘要：

| 输出 | 最小值 | 最大值 |
| --- | ---:| ---:|
| water outlet K | 333.947 | 333.947 |
| water delta-T K | 23.947 | 23.947 |
| sleeve outer mean K | 885.426 | 1052.574 |
| sleeve outer max K | 958.453 | 1144.173 |
| collector inner mean K | 925.106 | 1092.254 |
| collector inner max K | 1002.543 | 1188.263 |

所有 range checks 通过，最大水侧能量闭合误差约 `1.53e-10 W`。该包络不是实验曲线验证，只是说明在 Benke 当前公开整理参数范围内，模型预测温度区间是多少；后续数字化 Benke 12 点套筒温度后，应判断实验曲线是否落入该包络并计算 MAE/RMSE。

## 2026-07-02 套筒温度参数反演工具

新增 `benke_calibration.py`，用于在 Benke 套筒 12 点热电偶数据数字化后，对 `regulated_he_effective_k_w_m_k` 和 `water_h_w_m2_k` 做网格反演。

当前工具只接受真实数字化实验数据：

- 输入文件：`experimental_data/benke_sleeve_thermocouple_12pt_digitized.csv`
- 默认扫描参数：`regulated_he_effective_k_w_m_k = 0.073, 0.08, 0.087`，`water_h_w_m2_k = 528, 800, 1012`
- 输出文件：`calibration_results/benke_sleeve_calibration_grid.csv` 和 `calibration_results/benke_sleeve_calibration_best.json`

如果缺少真实数字化 CSV，工具会显式报错；不得用模型生成温度、Venable 输出功率或 Venable I-V 曲线伪造 Benke 实验数据。该工具的目的，是在 Benke 热工边界先验确定后，再把可信的热边界传递给 Venable_V 的点输出特性验证。

## 2026-07-02 验证状态判据收紧

`evaluate_benke_validation()` 的状态判据已收紧：

- `partial_missing_digitized_data`：没有真实 Benke 数字化测点，只能做文献范围校核和能量闭合检查。
- `quantitative_partial_with_digitized_data`：已有套筒 12 点温度或水侧热平衡中的一类真实数字化数据，可以做部分量化对比，但不能称为完整 Benke 热工水力验证。
- `complete_with_digitized_data`：套筒 12 点温度和水侧热平衡两类真实数字化数据都已接入，才允许称为完整量化验证。

新增 `benke_report.py`，`run_benke_thermal_validation.py` 每次运行会生成 `validation_report.md`。该报告会汇总输入参数、主要输出、文献范围校核、套筒温度对比、水侧热平衡对比和当前结论。

## 2026-07-02 用户补充 Benke 实验数据后的状态

用户已补充 Benke Appendix B p.79、TISA input power 3412 W row 的实验数据：

- `experimental_data/benke_sleeve_thermocouple_12pt_digitized.csv`
- `experimental_data/benke_water_balance_digitized.csv`
- `experimental_data/README_benke_vv_data.md`

数据满足当前最低接收要求。T64 / thermocouple index 9 被 Benke 说明为 inoperative，CSV 中保留为 `NaN`；验证脚本已更新为跳过非有限测点并报告 `ignored_indices = [9]`。

匹配边界运行：`runs/20260702_benke_user_data_matched_boundary/`

关键结果：

| 指标 | 数值 |
| --- | ---:|
| validation status | `complete_with_digitized_data` |
| range checks | `passed` |
| active-zone power | 3002.56 W |
| sleeve points used | 11 / 12 |
| sleeve MAE | 232.379 K |
| sleeve RMSE | 249.852 K |
| sleeve max abs error | 365.553 K |
| water outlet abs error | 0.946 K |
| water delta-T abs error | 0.946 K |

结论：Benke_V 已具备完整 digitized-data 验证链路，并已完成首版量化对比。当前 v1 热网络能较好闭合水侧热平衡，但显著高估套筒温度；后续应围绕径向热阻、并联/旁路热流、热电偶轴向映射和 He gap 有效导热进行模型改进或参数校准。

## 2026-07-02 热阻/热流参数调整结果

本轮按可调参数原则，只在 Benke_V 内调整热阻和热流相关输入，不修改主程序 py/C++ 文件。

新增 `coolant_heat_fraction` 参数，默认 `1.0`，用于表示进入当前水冷径向热网络的 active-zone 热量比例。

参数调整记录见：`PARAMETER_ADJUSTMENT_ROUND_20260702.md`

主要结果：

| 工况 | k_He W/(m K) | water h W/(m2 K) | coolant heat fraction | sleeve RMSE K | water delta-T abs error K | range checks |
| --- | ---:| ---:| ---:| ---:| ---:| --- |
| baseline | 0.08 | 800 | 1.0 | 249.852 | 0.946 | passed |
| strict range best | 0.087 | 1012 | 0.90 | 143.034 | 0.704 | passed |
| expanded sensitivity best | 0.10 | 10000 | 0.94 | 94.311 | 0.044 | failed |

结论：误差已显著降低，但扩展敏感性最优需要超出当前 Benke 整理的水侧换热系数范围。后续若继续降低误差，应优先处理轴向热源形状、T56-T67 轴向映射、水流方向和并联/旁路热流，而不是继续单纯放大水侧 h。

## 2026-07-03 Benke 热电偶映射口径修正

根据用户补充的 Benke 试验口径，本轮修正了验证对比方式：

- 12 个热电偶轴向位置显式写入模型：`-205, -163, -108, -55, -55, 0, 0, 55, T64 failed, 108, 163, 205 mm`。
- 验证轴向域默认改为 `0.410 m`，与 `-205` 到 `+205 mm` 测点跨度一致；TISA heater 受热长度保持 `0.300 m`，对应 `-150` 到 `+150 mm`。
- 套筒温度采样点从原来的等效外表面改为 collector sleeve 内部测点，半径为 sleeve 内圆柱面外侧 `1.8 mm`。
- T64 / index 9 保留为失效测点，计算和实验对比均按 `NaN` 跳过。
- Benke 一维模型中的 average collector sleeve temperature 不再平均全部 11 个有效点，而是剔除 T56 和 T67 端部 helium chamber 影响点，并剔除失效 T64；采用 `T57-T63, T65-T66` 共 9 点，即 index `2,3,4,5,6,7,8,10,11`。

本轮重新运行了三组映射后结果：

| 工况 | 运行目录 | k_He W/(m K) | water h W/(m2 K) | coolant heat fraction | 11 有效点 RMSE K | Benke 9 点平均误差 K | water delta-T abs error K | range checks |
| --- | --- | ---:| ---:| ---:| ---:| ---:| ---:| --- |
| baseline | `runs/20260703_benke_mapped_baseline/` | 0.08 | 800 | 1.00 | 243.057 | +205.473 | 0.946 | passed |
| strict range best | `runs/20260703_benke_mapped_strict_range_best/` | 0.087 | 1012 | 0.90 | 173.803 | +83.780 | 0.704 | passed |
| expanded sensitivity best | `runs/20260703_benke_mapped_expanded_best/` | 0.10 | 10000 | 0.94 | 157.887 | -20.921 | 0.044 | failed |

解释：修正测点半径和 Benke 平均口径后，`expanded sensitivity best` 对 Benke 平均 collector sleeve temperature 的误差已经降到约 `-20.9 K`，水侧热平衡也很好；但它依赖 `water_h = 10000 W/(m2 K)` 和 `k_He = 0.10 W/(m K)`，超出当前整理的 Benke 文献范围，不能作为已确认边界。严格文献范围内最优仍高估 Benke 9 点平均温度约 `83.8 K`。

剩余 11 有效点 RMSE 仍较大，说明轴向形状尚未闭合。尤其端部和 heater 外侧点需要后续引入轴向导热、端部 helium chamber/支撑结构影响，或更真实的 TISA 轴向功率分布；不应继续只靠提高水侧 `h` 来拟合。
