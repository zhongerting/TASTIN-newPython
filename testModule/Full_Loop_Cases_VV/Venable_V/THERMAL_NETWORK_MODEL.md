# Venable thermal_network_v1 热网络说明

本文记录 `thermal_network_v1` 的建模口径，用于替代纯经验给定发射极/收集极温度的第一版物理闭合。

## 目标

`thermal_network_v1` 从 Table 7-1 的 `Q_az` 出发，计算轴向发射极温度 `T_emitter(z)` 和收集极温度 `T_collector(z)`，再交给现有 ThermoCalc 电性能模型进行最大输出功率扫描。

该模型仍是 Venable_V 算例局部前处理模型，不改变核心 `py` 或 C++ 求解器。

## 当前模型口径

- `Q_az` 直接使用 Table 7-1 active-zone power，不再乘 Benke 的 `0.88`。
- 工作段长度为 `0.375 m`，TISA 加热段为 `0.300 m`，居中布置。
- 轴向热源按 TISA 加热形状分配，并严格归一化到 `Q_az * thermal_network_heat_pickup_fraction`。
- 冷却水自下而上逐节点升温，满足 `m_dot cp (T_out - T_in) = Q_to_water`。
- 每个轴向节点使用径向热阻链：Cs gap、collector、Al2O3、unregulated He、sleeve、regulated He、water film、extra resistance。
- ThermoCalc 使用热网络计算出的 emitter Cs-gap 面温度和 collector 内表面温度。

## 默认参数

| 参数 | 默认值 | 说明 |
| --- | ---:| --- |
| `thermal_model_mode` | `empirical` | 默认仍保留既有经验闭合，显式指定后才使用热网络。 |
| `cooling_water_inlet_temperature_k` | `310 K` | 试验冷却水入口温度占位值。 |
| `cooling_water_mass_flow_kg_s` | `0.03 kg/s` | 试验流量未完全明确时的诊断默认值。 |
| `water_heat_transfer_coefficient_w_m2_k` | `800 W/m2/K` | 位于 Benke 整理范围 `528-1012 W/m2/K` 中部。 |
| `regulated_he_gap_effective_k_w_m_k` | `0.08 W/m/K` | 对应 Benke 推断的调节 He gap 有效导热量级。 |
| `unregulated_he_gap_effective_k_w_m_k` | `0.276 W/m/K` | 近似高压 He 连续介质导热率。 |
| `cs_gap_effective_k_w_m_k` | `0.12 W/m/K` | 首版等效热导占位值，后续需由温度/功率对照校核。 |
| `thermal_network_heat_pickup_fraction` | `1.0` | 首版假定 active-zone 热输入最终由水侧带走。 |
| `collector_extra_resistance_k_per_w` | `0` | 用于后续诊断未建模接触热阻或结构热阻。 |

## 使用方式

示例：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_VV\Venable_V\venable_validation_runner.py --stage single_case_smoke --thermal-model-mode thermal_network_v1 --run-id thermal_network_v1_smoke
```

常用调整入口：

```text
--thermal-model-mode thermal_network_v1
--water-inlet-temperature-k
--water-mass-flow-kg-s
--water-h-w-m2-k
--regulated-he-gap-effective-k-w-m-k
--unregulated-he-gap-effective-k-w-m-k
--thermal-network-heat-pickup-fraction
--cs-gap-effective-k-w-m-k
--extra-thermal-resistance-k-w
```

## 后续调整规则

1. 先固定 `Q_az`、Table 7-1 输出目标和 Cs 压力分段，不用输出功率反推输入。
2. 先检查热网络能量闭合、水温单调性和温度量级。
3. 若电功率偏差较大，优先调整有试验依据的冷却边界：入口温度、流量、水侧换热系数。
4. 其次调整 He gap 有效导热率和额外热阻，并记录是否仍在文献可解释范围内。
5. 最后才考虑 Cs gap 等效热导或表面发射参数；这些属于更强模型不确定性，必须单独标记。
6. 如果误差无法继续降低，应保留真实结果并分析原因，不进行逐工况任意拟合。
