# Benke 参数调整记录：热阻与热流

## 调整目的

根据用户要求，本轮只调整 Benke_V 内可变输入参数，重点放在热阻和热流分配，目标是降低与 Benke Appendix B p.79、TISA input power 3412 W row 的验证误差。

不修改主程序 Python/C++ 文件；不修改 Benke 原始 digitized 数据。

## 基线工况

基线运行：`runs/20260702_benke_user_data_matched_boundary/`

| 项目 | 数值 |
| --- | ---:|
| TISA input power | 3412 W |
| active-zone power | 3002.56 W |
| water inlet | 289.71 K |
| water mass flow | 0.043518 kg/s |
| regulated He effective k | 0.08 W/(m K) |
| water h | 800 W/(m2 K) |
| coolant heat fraction | 1.0 |
| sleeve MAE | 232.379 K |
| sleeve RMSE | 249.852 K |
| sleeve max abs error | 365.553 K |
| sleeve mean error | 205.383 K |
| water outlet abs error | 0.946 K |
| water delta-T abs error | 0.946 K |

误差形态：水侧总热量基本闭合，但套筒温度明显偏高，说明主要问题在套筒到水侧的等效热阻、轴向热流分布/映射，或未建模的并联/旁路热流。

## 新增可调参数

新增 `BenkeThermalNetworkConfig.coolant_heat_fraction`，默认值为 `1.0`。

含义：active-zone power 仍按 `P_az = 0.88 P_TISA` 记录，但进入当前水冷径向热网络的热量为：

```text
Q_coolant_network = P_az * coolant_heat_fraction
```

该参数用于表达未进入当前冷却水径向热网络的旁路散热、端部损失或试验台环境热损失。它不改变 Benke active-zone 修正本身。

## 扫描范围

运行目录：`calibration_results/20260702_adjust_heat_resistance_flow/`

扫描参数：

| 参数 | 扫描值 |
| --- | --- |
| regulated He effective k W/(m K) | 0.073, 0.08, 0.087, 0.10, 0.12, 0.15, 0.20 |
| water h W/(m2 K) | 528, 800, 1012, 1500, 2500, 5000, 10000, 20000 |
| coolant heat fraction | 0.90, 0.94, 0.943, 0.95, 1.0 |

## 严格文献范围内最优

运行目录：`runs/20260702_benke_adjusted_strict_range_best/`

| 项目 | 数值 |
| --- | ---:|
| regulated He effective k | 0.087 W/(m K) |
| water h | 1012 W/(m2 K) |
| coolant heat fraction | 0.90 |
| range checks | passed |
| sleeve MAE | 132.818 K |
| sleeve RMSE | 143.034 K |
| sleeve max abs error | 230.597 K |
| sleeve mean error | 85.303 K |
| water outlet abs error | 0.704 K |
| water delta-T abs error | 0.704 K |

相对基线：

| 指标 | 基线 | 调整后 | 降低 |
| --- | ---:| ---:| ---:|
| sleeve MAE | 232.379 K | 132.818 K | 99.562 K |
| sleeve RMSE | 249.852 K | 143.034 K | 106.818 K |
| water delta-T abs error | 0.946 K | 0.704 K | 0.242 K |

## 扩展敏感性范围最优

运行目录：`runs/20260702_benke_adjusted_expanded_best/`

| 项目 | 数值 |
| --- | ---:|
| regulated He effective k | 0.10 W/(m K) |
| water h | 10000 W/(m2 K) |
| coolant heat fraction | 0.94 |
| range checks | failed |
| sleeve MAE | 70.717 K |
| sleeve RMSE | 94.311 K |
| sleeve max abs error | 221.745 K |
| sleeve mean error | -18.881 K |
| water outlet abs error | 0.044 K |
| water delta-T abs error | 0.044 K |

相对基线：

| 指标 | 基线 | 调整后 | 降低 |
| --- | ---:| ---:| ---:|
| sleeve MAE | 232.379 K | 70.717 K | 161.662 K |
| sleeve RMSE | 249.852 K | 94.311 K | 155.541 K |
| water delta-T abs error | 0.946 K | 0.044 K | 0.902 K |

## 解释与约束

1. 严格文献范围内的最优参数仍明显高估套筒温度，说明仅靠 Benke 整理的 `k_He` 和 `h_water` 上限不足以完全解释实验套筒温度。
2. 扩展敏感性最优大幅降低误差，但 `water_h = 10000 W/(m2 K)` 明显超出当前整理的 `528-1012 W/(m2 K)`，不能直接作为已确认物理边界。
3. `coolant_heat_fraction = 0.94` 与水侧热平衡更一致，表示约 6% active-zone 热量没有进入当前水冷径向网络；这可能对应端部损失、环境散热或未建模旁路热流。
4. 剩余 RMSE 仍约 94 K，且 max abs error 约 222 K，说明还有轴向形状问题。下一步不应继续盲目加大 `h_water`，而应检查：
   - T56-T67 与轴向坐标的对应关系；
   - 水流方向与 thermocouple 编号方向；
   - TISA 轴向热源是否偏心或非平台分布；
   - sleeve thermocouple 是否测在槽内/局部位置而非当前 `sleeve_outer` 等效面；
   - 是否需要显式加入并联辐射/端部导热路径。

## 当前推荐

用于后续 Venable 热边界传递时，优先采用两套参数并列报告：

1. 保守参数：`k_He = 0.087 W/(m K)`, `h_water = 1012 W/(m2 K)`, `coolant_heat_fraction = 0.90`。
2. 误差最小敏感性参数：`k_He = 0.10 W/(m K)`, `h_water = 10000 W/(m2 K)`, `coolant_heat_fraction = 0.94`。

第二套只能作为敏感性最优，不应直接称为文献确认参数。
