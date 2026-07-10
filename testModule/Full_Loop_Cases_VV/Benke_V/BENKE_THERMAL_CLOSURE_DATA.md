# Benke 热模型闭合数据分类

本文用于把 Benke 单电池 TFE 试验台中可用于热工水力/热阻模型闭合的数据分成三类：可直接锁定、可反演确定、只能作为敏感性参数。后续 Venable 电输出验证只能使用这里形成的热边界先验，不能用 Venable 输出功率逐点反推热边界。

## 数据来源与适用范围

主要来源是 `TOPAZII_VV_public_experimental_data.md` 中 Benke 单电池 TFE 试验台整理内容：

- Benke, *Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand*, NPS thesis, 1994。
- Benke & Venable AIP 论文作为辅助来源。
- Venable 1995 与 Benke 使用相同或高度相近的 TOPAZ-II single-cell TFE 试验台，但 Venable 数据主要用于电输出验证，不用于反推 Benke 热边界。

Benke 试验台为非核电加热单电池 TFE：未装燃料，用 TISA 钨电阻加热器替代核释热。

## 可直接锁定的数据

这些数据来自结构、几何或明确实验口径。除非找到更高优先级原文，否则建模时应固定。

| 数据项 | 数值/说明 | 建模用途 |
| --- | ---:| --- |
| TISA 加热器类型 | 钨电阻加热器 | 内部热源模型 |
| TISA 最大电参数 | 约 29 VAC / 170 A，最大约 4500 W | 功率范围 sanity check |
| TISA 加热长度 | 300.0 mm | 轴向热源分布 |
| TISA 内导体 | W，直径约 6.5 mm | 结构说明，首版可不显式建模 |
| TISA 外导体 | W，外径约 7.0 mm，厚度约 0.4 mm | 结构说明，首版可不显式建模 |
| thermionic 工作段长度 | 375.0 mm | 轴向建模长度 |
| emitter 材料 | Mo-3%Nb 单晶，W 涂层 | 径向结构/材料说明 |
| emitter 几何 | 外径约 19.6 mm，厚度约 1.15 mm，W 涂层约 0.1 mm | 电极几何和导热路径 |
| Cs vapor gap | 约 0.5 mm | emitter-collector 间隙热阻 |
| collector 材料 | 多晶 Mo | collector 导热层 |
| collector 几何 | 内径约 20.6 mm，厚度约 1.4 mm | collector 导热层 |
| Al2O3 绝缘层 | 厚度约 0.15 mm | collector 外侧绝缘热阻 |
| unregulated He gap | 约 0.05 mm，通常 200-300 torr | Al2O3 到 sleeve 的热阻 |
| collector sleeve | 1X18H10T stainless，厚度约 3.0 mm，外径约 29.9 mm | sleeve 导热层 |
| sleeve 热电偶槽 | 12 个槽，深约 2.0 mm | 轴向温度验证点位置说明 |
| regulated He gap | 约 0.5 mm，正常运行 1-10 torr | 主要可控保温热阻 |
| water jacket | 1X18H10T，内壁 2.5 mm，外壁 1.0 mm | 冷却边界结构 |
| spiral coolant channel | pitch 约 35 mm，约 6.5 turns，水自下而上流动 | 轴向分段水温推进方向 |
| active-zone 修正 | `P_az = 0.88 P_TISA` | Benke 热模型中从 TISA 输入到 active-zone 热输入的换算 |
| 水侧热平衡公式 | `Q_w = m_dot cp (T_out - T_in)` | 能量闭合验证 |
| 圆筒径向导热公式 | `R = ln(r2/r1)/(2*pi*k*L)` | 径向热阻网络 |
| 圆筒辐射公式 | 同轴圆筒辐射换热 | 后续非线性/并联热阻扩展 |
| 总热阻链 | `R_Cs + R_collector + R_Al2O3 + R_He_unreg + R_sleeve + R_He_reg + R_water` | 首版热网络骨架 |

## 可反演确定的数据

这些数据不是单纯输入，而是可以由 Benke 的实验输出或典型工况反推，用来校准/验证热模型。

| 数据项 | 公开整理值/输出 | 反演或验证方式 |
| --- | ---:| --- |
| active-zone power | 典型约 3003 W | 由 `P_TISA` 和 0.88 修正得到；用于热网络输入 |
| regulated He gap effective k | 约 0.073-0.087 W/(m K) | 由温度分布/热流反推，与模型值比较 |
| water-side h | 约 528-1012 W/(m2 K) | 由水侧热平衡、几何和壁温反推 |
| cooling water Re | 约 1480 | 由流量和通道水力直径反推；用于判断层流/过渡流 |
| cooling water outlet temperature | Benke 输出量之一 | 由 `m_dot cp DeltaT` 与输入功率闭合 |
| sleeve 12 点轴向温度 | Benke 输出量之一 | 数字化后比较 MAE/RMSE 和峰值位置 |
| axial temperature distribution | 图像输出 | 数字化后比较趋势和热源段峰值位置 |
| active-zone input vs water heat uptake | Benke 输出量之一 | 能量闭合误差 |

## 只能作为敏感性参数的数据

这些量公开资料不足以唯一确定，不能逐工况任意拟合。后续只能做全局敏感性或有界校准。

| 参数 | 为什么不能直接锁定 | 建议用法 |
| --- | --- | --- |
| cooling water inlet temperature | 当前整理未给出逐工况入口温度表 | 固定为全局假设，做敏感性扫描 |
| cooling water mass flow | 当前整理只给 Re 约 1480，未给逐工况流量表 | 由 Re/水力直径反推或做全局扫描，不逐点拟合 |
| Cs gap effective k | Cs 间隙内导热/辐射/电子热输运复杂 | 首版给工程占位，后续由热/电联合验证约束 |
| unregulated He gap effective k | 压力较高但装配、稀薄效应和接触不确定 | 默认接近连续 He k，做有界敏感性 |
| regulated He gap rarefaction model | 公开整理只给有效 k 范围 | 先直接使用有效 k，后续再建立 pressure-to-k 关系 |
| 材料导热率随温度变化 | 资料整理未给完整温度函数 | 首版用常数，后续替换为材料库或分段值 |
| 辐射换热 | 公式已给，但表面发射率/温度非线性未充分给定 | 首版可关闭或作为并联敏感性 |
| 接触热阻/装配偏差 | 文献未完整给定 | 只能作为全局 extra resistance，严禁逐点拟合 |
| spiral channel 详细水力直径 | 仅有 pitch 和 turns，通道截面细节不足 | 首版用等效水侧 h，后续再补水力细节 |

## 对 Venable 验证的传递规则

1. Benke 用于确定热边界和热阻先验，Venable 用于验证电输出。
2. Venable Table 7-1 的输出功率不能反向作为 Benke 热模型输入。
3. 如果 Benke 热模型不能唯一确定某参数，该参数进入 Venable 时必须保持全局一致或来自明确工况记录。
4. 若 Benke 复现误差较大，应如实记录边界不足，而不是通过 Venable 逐点调参掩盖。

## 2026-07-03 热电偶测点与平均温度口径补充

用户补充的 Benke 试验口径应作为后续 Benke_V 热工验证的硬约束：

| 项目 | 口径 | 建模/验证处理 |
| --- | --- | --- |
| 12 点热电偶轴向位置 | `-205, -163, -108, -55, -55, 0, 0, 55, T64 failed, 108, 163, 205 mm` | 写入 `BENKE_THERMOCOUPLE_POSITIONS_MM`；T64 用 `None/NaN` 保留占位。 |
| TISA heater 受热长度 | `300 mm`，对应 `-150` 到 `+150 mm` | `tisa_heated_length_m = 0.300`。当前 v1 仍用端部渐变热源作为轴向导热/端部效应的工程近似，后续若建立显式轴向导热，应重新审查该形状。 |
| 验证轴向域 | 热电偶覆盖 `-205` 到 `+205 mm` | `active_length_m = 0.410`，用于轴向采样和分段水温推进。 |
| sleeve 测温半径 | 热电偶位于 stainless collector sleeve 内部，距 sleeve 内圆柱面 `1.8 mm` | 计算 `sleeve_thermocouple_temperature_k`，不再把 `sleeve_outer_temperature_k` 直接当实验测点。 |
| Benke 平均 collector sleeve temperature | 不平均全部 11 个有效点；T56/T67 受上下端 helium chamber 端部效应影响而剔除，T64 失效 | 使用 index `2,3,4,5,6,7,8,10,11` 共 9 点计算 `experimental_benke_average_k` 和 `calculated_benke_average_k`。 |

这些测点是 test stand 内部实际测得的唯一温度量。后续参数调整报告必须同时列出普通 11 有效点误差和 Benke 9 点平均误差，避免把端部形状误差与 Benke 平均温度验证目标混在一起。
