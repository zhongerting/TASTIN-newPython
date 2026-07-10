# Benke 参数调整记录：热电偶映射口径修正

## 目的

本轮不是继续盲目调参，而是先修正 Benke sleeve thermocouple 的验证口径：轴向位置、测温半径、失效测点和 Benke 9 点平均规则。修正后重新评估 2026-07-02 的 baseline、严格文献范围最优和扩展敏感性最优三组参数。

不修改主程序 Python/C++ 文件；不修改 Benke 原始 digitized 数据。

## 口径修正

| 项目 | 修正后口径 |
| --- | --- |
| 验证轴向域 | `active_length_m = 0.410 m` |
| TISA heater 长度 | `tisa_heated_length_m = 0.300 m` |
| 热电偶轴向位置 | `-205, -163, -108, -55, -55, 0, 0, 55, failed, 108, 163, 205 mm` |
| 热电偶测温半径 | sleeve 内圆柱面外侧 `1.8 mm`，位于 stainless sleeve 内部 |
| 失效测点 | T64 / index 9 保留为 `NaN`，不参与统计 |
| Benke 平均温度 | index `2,3,4,5,6,7,8,10,11` 共 9 点 |

## 重新运行结果

| 工况 | 运行目录 | k_He W/(m K) | water h W/(m2 K) | coolant heat fraction | sleeve MAE K | sleeve RMSE K | sleeve mean error K | Benke 9 点平均误差 K | water delta-T abs error K | range checks |
| --- | --- | ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| --- |
| baseline | `runs/20260703_benke_mapped_baseline/` | 0.08 | 800 | 1.00 | 226.836 | 243.057 | +109.393 | +205.473 | 0.946 | passed |
| strict range best | `runs/20260703_benke_mapped_strict_range_best/` | 0.087 | 1012 | 0.90 | 144.196 | 173.803 | +8.234 | +83.780 | 0.704 | passed |
| expanded sensitivity best | `runs/20260703_benke_mapped_expanded_best/` | 0.10 | 10000 | 0.94 | 101.295 | 157.887 | -78.622 | -20.921 | 0.044 | failed |

## 对比解释

1. Benke 9 点平均温度是当前最接近文献一维热模型 average collector sleeve temperature 的指标。按这个口径，扩展敏感性最优已经从严重高估转为轻微低估，误差约 `20.9 K`。
2. 严格文献范围内最优仍高估 Benke 9 点平均温度约 `83.8 K`，说明仅靠当前范围内的 `k_He`、`h_water` 和 `coolant_heat_fraction` 不能完全闭合。
3. 11 有效点 RMSE 仍较大，特别反映轴向形状未闭合。这个误差不能只用径向热阻解释，更可能涉及 heater 外侧轴向导热、端部 helium chamber 影响、支撑结构热流和水流/测点方向映射。
4. 扩展敏感性最优的 `water_h = 10000 W/(m2 K)` 明显超出当前 Benke 整理的 `528-1012 W/(m2 K)`，只能作为误差下限敏感性，不应作为已确认物理参数。

## 后续可调方向

优先级从高到低：

1. 建立显式轴向导热或端部热泄漏模型，替代当前一维径向独立节点近似。
2. 重新审查 TISA heater 外侧功率形状：当前 v1 保留端部渐变热源作为工程近似；若轴向导热显式建模，应把 heater 外侧直接热源与轴向传热分开。
3. 在不超出文献范围的前提下继续扫描 `coolant_heat_fraction`、`extra_resistance_k_per_w` 和 regulated He 有效导热率。
4. 保留扩展敏感性参数作为 Venable 后续验证的不确定性包络，而非单一推荐值。
