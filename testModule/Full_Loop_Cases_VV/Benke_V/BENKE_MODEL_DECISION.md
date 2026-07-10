# Benke 热模型形式决策

## 结论

Benke_V 首版应采用 **轴向分段水冷径向热阻网络**，而不是单一 lumped 总热阻模型。

理由：

- Benke 可验证输出包含 sleeve 12 点轴向温度和轴向温度分布，lumped 模型无法比较峰值位置和轴向趋势。
- TISA 加热长度 300 mm，而 thermionic 工作段 375 mm，热源本身具有轴向非均匀性。
- 冷却水沿 spiral channel 自下而上流动，水温边界应随轴向推进。
- regulated He gap 是核心保温对象，其有效导热率可作为径向热阻链中的明确参数。
- 首版仍可保持简单：每个轴向节点用局部径向热阻，水温一维推进；暂不引入完整二维固体导热和非线性辐射。

## 首版模型范围

`benke_thermal_network_v1` 计算以下量：

- TISA 输入功率 `P_TISA`。
- active-zone 热输入 `P_az = 0.88 P_TISA`。
- 轴向热源分布 `q_i`，TISA 300 mm 居中于 375 mm 工作段。
- 每个轴向节点的径向热阻：
  - Cs gap
  - collector
  - Al2O3
  - unregulated He gap
  - collector sleeve
  - regulated He gap
  - water film
  - optional extra resistance
- 冷却水 bulk 温度自下而上推进。
- collector inner surface、sleeve outer surface、water bulk 的轴向温度。
- 水侧能量闭合误差。

## 首版暂不纳入

以下内容先记录为后续扩展，不进入首版 smoke 验证：

- 真实 spiral coolant channel 的二维/三维几何水力求解。
- 轴向固体导热耦合。
- 同轴圆筒辐射换热的非线性迭代。
- regulated He gap pressure-to-effective-k 的稀薄气体模型。
- ThermoCalc 电输出反馈热源。

## 验证指标

首版 Benke 验证先采用可由整理资料直接约束的指标：

| 指标 | 目标/检查 |
| --- | --- |
| active-zone power | `P_az = 0.88 P_TISA`，典型工况约 3003 W |
| regulated He effective k | 默认落在 0.073-0.087 W/(m K) |
| water-side h | 默认落在 528-1012 W/(m2 K) |
| water energy balance | `m_dot cp DeltaT` 与 `P_az` 闭合 |
| water outlet temperature | 由模型输出，后续与 Benke 数字化数据比较 |
| sleeve axial temperature | 输出 12 点采样值，后续与 Benke 12 点热电偶数字化数据比较 |
| temperature trend | 中部加热段温度高，两端温度低；水流方向导致上游/下游存在合理偏置 |

## 与 Venable_V 的关系

Benke_V 的输出应作为 Venable_V 热边界先验：

- `collector_temperature_k(z)` 可传递给 Venable 电输出模型的 collector 温度。
- 若后续扩展 emitter 内部热源到 emitter surface 的导热，也可传递 `emitter_temperature_k(z)`。
- 进入 Venable 后，优先保持 Benke 已确定参数不变；只允许全局小范围敏感性调整。
