# V13 CaseA 流网与固体耦合拓扑说明

生成日期：2026-06-30

本文对应 `testModule/test_core_assemble_v13_caseA.py` 和默认运行入口 `testModule/run_v13_caseA_closed_loop.py`。命名以源码对象名为准；尺寸只采用当前 V13/V12/V8 源码中的构造参数和实例化对象字段。若水力对象只定义了 `flow_area` 和 `hydraulic_diam`，本文不推断其壁厚或外径，而标记为 `not modeled`。

![V13 flow and solid coupling topology](./v13_flow_solid_coupling_topology.svg)

## 1. 适用范围

V13 是闭式泵驱动的“堆芯 + TOPAZ-II 管翅式辐射器”回路。与 V11 的集流环/热管/翅片散热路径不同，V13 将堆芯出口直接接入 78 根管翅式辐射管。

```text
V13_PumpOutletNode
  -> Pipe11_CoreInletHeader
  -> V12_CoreInletDistribution
  -> V12_CoreInletBranch_1
     + V12_CoreInletBranch_2_3_Rep(multiplier=2)
  -> V12_CoreInletConnector
  -> Chan_Center / Chan_Ring1 / Chan_Ring2 / Chan_Ring3_TEC / Chan_Ring3_Open
  -> V12_CoreOutletConnector
  -> Pipe05_CoreOutletToRadiator
  -> V12_RadiatorInletSplit
  -> V12_RadiatorUpperHeader_01...78
  -> V12_RadiatorTubeFluid_01...78
  -> V12_RadiatorLowerHeader_01...78
  -> V12_RadiatorOutletMix
  -> Pipe06_RadiatorOutlet
  -> Pipe07_HeatExchangerHotSide
  -> Pipe08_ReturnInnerPipe
  -> Pipe09_ValveSegment
  -> V13_PumpMidNode
  -> V13_PumpOutletNode
```

压力约束：

- `V12_CoreInletConnector.is_pressure_reference=True`，默认 `target_P=207927.58 Pa`。
- V13 闭式回路中没有 `is_pressure_boundary=True` 的固定压力流动边界体。
- `J_Pipe09_to_V13_PumpA` 和 `J_V13_PumpA_to_PumpB` 是两台串联泵。默认 `pump_total_head_pa=7900 Pa`，每台 `3950 Pa`。

默认运行入口参数：

| 参数 | 默认值 | 说明 |
| --- | ---: | --- |
| `--total-inlet-flow-kg-s` | `1.3 kg/s` | 设计总质量流量，用于初始化与诊断 |
| `--target-flow-kg-s` | `1.3 kg/s` | 泵流量控制目标 |
| `--pump-total-head-pa` | `7900 Pa` | 两台泵合计压头 |
| `--inlet-temperature-k` | `727.0 K` | 构建/注入初始温度；闭式回路中不是固定温度边界 |
| `--tube-emissivity` | `0.80` | 辐射管裸管发射率 |
| `--fin-emissivity` | `0.80` | 翅片发射率 |
| `--t-space-k` | `3.0 K` | 管翅式辐射器默认深空背景温度 |
| `--thermo-update-interval` | `1.0 s` | TEC 计算调用间隔 |
| `--wire-resistance-scale` | `0.5` | TEC 连接导线电阻倍率 |
| `--solid-ode-method` | `RK45` | runner 默认把公共默认值切换为 `RK45` |
| `--fluid-solid-coupling-scheme` | `local_implicit` | runner 默认应用到具备固体热容的 `FluidSolidCouple` |

直接实例化默认 V13 的核对结果：

| 项目 | 数量 |
| --- | ---: |
| hydraulic volumes | 1036 |
| hydraulic junctions | 1122 |
| registered components | 79 (`TASTIN_Core_V8_CaseA` + 78 个辐射管组件) |
| registered solids | 114 |
| registered couplers | 113 |

## 2. 水力体与几何尺寸

表中 `ID` 对圆管取 `hydraulic_diam` 或实际内径。若对象只建模了流动面积和水力直径，没有壁厚或外径字段，则 `OD` 标为 `not modeled`。

### 2.1 主回路管段与连接体

| 对象名 | 类型/倍率 | 节点数 | 长度 L (m) | 流通面积 A (m2) | ID 或 Dh (m) | OD (m) | 固体换热 |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| `V13_PumpOutletNode` | `IncompressibleFluidVolume` | 1 | `0.02` 构造参数 | `3.800000e-4` | `0.014000` | not modeled | none |
| `Pipe11_CoreInletHeader` | 圆管水力通道 | 8 | `0.130000` | `3.800000e-4` | `0.014000` | not modeled | none |
| `V12_CoreInletDistribution` | 分配体 | 1 | `0.02` 构造参数 | `3.800000e-4` | `0.014000` | not modeled | none |
| `V12_CoreInletBranch_1` | 圆管，显式 1 根 | 8 | `1.890210` | `5.979816e-4` | `0.027600` | not modeled | none |
| `V12_CoreInletBranch_2_3_Rep` | 圆管，宏观代表 2 根 | 8 | `2.507050` | `5.979816e-4` | `0.027600` | not modeled | none |
| `V12_CoreInletConnector` | 压力参考连接体 | 1 | `0.02` 构造参数 | `2.050455e-3` | `0.001400` | not modeled | 只通过堆芯通道流动连接 |
| `V12_CoreOutletConnector` | 堆芯出口连接体 | 1 | `0.02` 构造参数 | `2.050455e-3` | `0.001400` | not modeled | 只通过堆芯通道流动连接 |
| `Pipe05_CoreOutletToRadiator` | 圆管水力通道 | 8 | `0.130000` | `3.800000e-4` | `0.014000` | not modeled | none |
| `V12_RadiatorInletSplit` | 辐射器入口分流体 | 1 | `0.02` 构造参数 | `3.800000e-4` | `0.014000` | not modeled | none |
| `V12_RadiatorOutletMix` | 辐射器出口汇流体 | 1 | `0.02` 构造参数 | `3.800000e-4` | `0.014000` | not modeled | none |
| `Pipe06_RadiatorOutlet` | 圆管水力通道 | 8 | `0.043408` | `3.800000e-4` | `0.014000` | not modeled | none |
| `Pipe07_HeatExchangerHotSide` | 圆管水力通道 | 8 | `0.005426` | `3.800000e-4` | `0.067900` | not modeled | none |
| `Pipe08_ReturnInnerPipe` | 圆管水力通道 | 8 | `0.130000` | `3.800000e-4` | `0.047000` | not modeled | none |
| `Pipe09_ValveSegment` | 圆管水力通道 | 8 | `0.130000` | `3.800000e-4` | `0.047000` | not modeled | none |
| `V13_PumpMidNode` | 泵间节点 | 1 | `0.02` 构造参数 | `3.800000e-4` | `0.014000` | not modeled | none |

`Pipe07_HeatExchangerHotSide` 当前只是水力网络保留段。V13 没有给它挂接换热器固体，也没有给它设置换热边界。

### 2.2 堆芯冷却剂通道

5 个代表性 TFE 通道使用相同轴向网格和冷却剂环隙几何。

| 对象名 | 代表倍率 | TEC 倍率 | 节点数 | 长度 L (m) | 流通面积 A (m2) | 冷却剂环隙 ID/OD (m) | Dh (m) | 固体耦合 |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |
| `Chan_Center` | 1 | 1 | 37 | `0.507000` | `5.541769e-5` | `0.02450/0.02590` | `0.001400` | `Center_iclad_coolant_couple`, `Center_oclad_coolant_couple` |
| `Chan_Ring1` | 6 | 6 | 37 | `0.507000` | `5.541769e-5` | `0.02450/0.02590` | `0.001400` | `Ring1_iclad_coolant_couple`, `Ring1_oclad_coolant_couple` |
| `Chan_Ring2` | 12 | 12 | 37 | `0.507000` | `5.541769e-5` | `0.02450/0.02590` | `0.001400` | `Ring2_iclad_coolant_couple`, `Ring2_oclad_coolant_couple` |
| `Chan_Ring3_TEC` | 15 | 15 | 37 | `0.507000` | `5.541769e-5` | `0.02450/0.02590` | `0.001400` | `Ring3_TEC_iclad_coolant_couple`, `Ring3_TEC_oclad_coolant_couple` |
| `Chan_Ring3_Open` | 3 | 0 | 37 | `0.507000` | `5.541769e-5` | `0.02450/0.02590` | `0.001400` | `Ring3_Open_iclad_coolant_couple`, `Ring3_Open_oclad_coolant_couple` |

TFE 轴向分段：

| 区段 | 长度 (m) | 节点数 | 单节点长度 (m) |
| --- | ---: | ---: | ---: |
| 下反射/包层段 | `0.065` | 6 | `0.010833333` |
| 活性段 | `0.377` | 25 | `0.015080000` |
| 上反射/包层段 | `0.065` | 6 | `0.010833333` |

### 2.3 管翅式辐射器水力几何

V13 辐射器由 78 根管、上集流环 78 个单节点段、下集流环 78 个单节点段组成。

| 对象模式 | 数量 | 每个节点数 | 每段长度 L (m) | 流通面积 A (m2) | ID/Dh (m) | OD (m) | 固体耦合 |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `V12_RadiatorUpperHeader_01...78` | 78 | 1 | `0.0331881070` | `3.141593e-4` | `0.020000` | not modeled | none |
| `V12_RadiatorTubeFluid_01...78` | 78 | 8 | `1.850000` | `3.848451e-5` | `0.007000` | `0.008000` through wall object | `V12_RadiatorTube_XX_FluidSolid` |
| `V12_RadiatorLowerHeader_01...78` | 78 | 1 | `0.0542126117` | `3.141593e-4` | `0.020000` | not modeled | none |

辐射管固体壁面：

| 固体对象模式 | 数量 | 轴向节点 | 径向壁节点 | 管长 (m) | 内径 (m) | 外径 (m) | 壁厚 (m) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `V12_RadiatorTube_01...78_Wall` | 78 | 8 | 1 | `1.850000` | `0.007000` | `0.008000` | `0.000500` |

管翅式辐射器缩阶翅片几何：

| 参数 | 值 |
| --- | ---: |
| `fin_thickness_m` | `0.000400` |
| `n_fin_width` | 12 |
| `fin_width_upper_m` | `0.033190` |
| `fin_width_lower_m` | `0.054210` |
| `fin_area_scale` | `0.35` |
| `fin_view_factor` | `1.0` |
| `fin_conductivity_w_m_k` | `348.9` |
| `fin_contact_resistance_m2k_w` | `0.0` |

集流环段当前是水力-only 对象；没有单独的集流环壁面导热固体，也没有集流环对深空辐射边界。

### 2.4 泵段

| 连接名 | from -> to | 类型 | 默认压头 | 默认流量目标 |
| --- | --- | --- | ---: | ---: |
| `J_Pipe09_to_V13_PumpA` | `Pipe09_ValveSegment_Vol_08` -> `V13_PumpMidNode` | `FlowControlledPumpJunction` | `3950 Pa` | `1.3 kg/s` |
| `J_V13_PumpA_to_PumpB` | `V13_PumpMidNode` -> `V13_PumpOutletNode` | `FlowControlledPumpJunction` | `3950 Pa` | `1.3 kg/s` |
| `J_V13_PumpOutlet_to_Pipe11` | `V13_PumpOutletNode` -> `Pipe11_CoreInletHeader_Vol_01` | `FlowJunction` | none | loop flow |

## 3. 外部主连接

### 3.1 堆芯入口、出口和冷端分配

| 连接名 | 类型 | from -> to | 面积 (m2) | K/倍率 |
| --- | --- | --- | ---: | --- |
| `J_V13_PumpOutlet_to_Pipe11` | `FlowJunction` | `V13_PumpOutletNode` -> `Pipe11_CoreInletHeader_Vol_01` | `3.800000e-4` | `K=0` |
| `J_Pipe11_to_CoreInletDistribution` | `FlowJunction` | `Pipe11_CoreInletHeader_Vol_08` -> `V12_CoreInletDistribution` | `3.800000e-4` | `K=0` |
| `J_CoreInletDistribution_to_CoreInletBranch_1` | `FlowJunction` | `V12_CoreInletDistribution` -> `V12_CoreInletBranch_1_Vol_01` | `3.800000e-4` | `K=0` |
| `J_CoreInletDistribution_to_CoreInletBranch_2_3_Rep` | `MacroFlowJunction` | `V12_CoreInletDistribution` -> `V12_CoreInletBranch_2_3_Rep_Vol_01` | `3.800000e-4` | `multiplier=2`, `K=0` |
| `J_V12_CoreInletBranch_1_to_CoreInletConnector` | `FlowJunction` | `V12_CoreInletBranch_1_Vol_08` -> `V12_CoreInletConnector` | `5.979816e-4` | branch outlet |
| `J_V12_CoreInletBranch_2_3_Rep_to_CoreInletConnector` | `MacroFlowJunction` | `V12_CoreInletBranch_2_3_Rep_Vol_08` -> `V12_CoreInletConnector` | `5.979816e-4` | `multiplier=2` |
| `J_PlenumIn_*` | macro inlet junctions | `V12_CoreInletConnector` -> `Chan_*_Vol_01` | channel area | representative multiplier |
| `J_PlenumOut_*` | macro outlet junctions | `Chan_*_Vol_37` -> `V12_CoreOutletConnector` | channel area | representative multiplier |
| `J_CoreOutletConnector_to_Pipe05` | `FlowJunction` | `V12_CoreOutletConnector` -> `Pipe05_CoreOutletToRadiator_Vol_01` | `3.800000e-4` | `K=0` |
| `J_Pipe05_to_RadiatorInletSplit` | `FlowJunction` | `Pipe05_CoreOutletToRadiator_Vol_08` -> `V12_RadiatorInletSplit` | `3.800000e-4` | `K=0` |

### 3.2 辐射器集流环和管束连接

| 连接模式 | 数量 | from -> to | K |
| --- | ---: | --- | ---: |
| `J_RadiatorInletSplit_to_UpperHeader_A` | 1 | `V12_RadiatorInletSplit` -> `V12_RadiatorUpperHeader_01_Vol_01` | `connector_k_loss=0` |
| `J_RadiatorInletSplit_to_UpperHeader_B` | 1 | `V12_RadiatorInletSplit` -> `V12_RadiatorUpperHeader_40_Vol_01` | `connector_k_loss=0` |
| `J_RadiatorUpperRing_XX_to_YY` | 78 | upper header `XX` -> next upper header segment | `radiator_header_k_loss=1` |
| `J_RadiatorUpper_to_Tube_XX` | 78 | upper header `XX` -> `V12_RadiatorTubeFluid_XX_Vol_01` | `radiator_tube_inlet_k_loss=100` |
| `J_RadiatorTube_XX_to_Lower` | 78 | `V12_RadiatorTubeFluid_XX_Vol_08` -> lower header `XX` | `radiator_tube_outlet_k_loss=100` |
| `J_RadiatorLowerRing_XX_to_YY` | 78 | lower header `XX` -> next lower header segment | `radiator_header_k_loss=1` |
| `J_LowerHeader_A_to_RadiatorOutletMix` | 1 | `V12_RadiatorLowerHeader_01_Vol_01` -> `V12_RadiatorOutletMix` | `connector_k_loss=0` |
| `J_LowerHeader_B_to_RadiatorOutletMix` | 1 | `V12_RadiatorLowerHeader_40_Vol_01` -> `V12_RadiatorOutletMix` | `connector_k_loss=0` |

`XX` 从 `01` 到 `78`；`YY` 为下一个编号并首尾闭合。

### 3.3 冷端回流连接

| 连接名 | from -> to | 说明 |
| --- | --- | --- |
| `J_V12_RadiatorOutletMix_to_Pipe06_RadiatorOutlet` | `V12_RadiatorOutletMix` -> `Pipe06_RadiatorOutlet_Vol_01` | 水力-only |
| `J_Pipe06_RadiatorOutlet_Vol_08_to_Pipe07_HeatExchangerHotSide` | `Pipe06_RadiatorOutlet_Vol_08` -> `Pipe07_HeatExchangerHotSide_Vol_01` | 水力-only |
| `J_Pipe07_HeatExchangerHotSide_Vol_08_to_Pipe08_ReturnInnerPipe` | `Pipe07_HeatExchangerHotSide_Vol_08` -> `Pipe08_ReturnInnerPipe_Vol_01` | 水力-only |
| `J_Pipe08_ReturnInnerPipe_Vol_08_to_Pipe09_ValveSegment` | `Pipe08_ReturnInnerPipe_Vol_08` -> `Pipe09_ValveSegment_Vol_01` | 水力-only |
| `J_Pipe09_to_V13_PumpA` | `Pipe09_ValveSegment_Vol_08` -> `V13_PumpMidNode` | 泵 A |
| `J_V13_PumpA_to_PumpB` | `V13_PumpMidNode` -> `V13_PumpOutletNode` | 泵 B |

## 4. 堆芯固体与边界条件

### 4.1 TFE 代表元件尺寸

每个代表 TFE 的轴向长度均为 `0.507 m`。径向几何如下，所有代表元件相同：

| 固体/间隙 | 对象或接口 | 内半径 (m) | 外半径 (m) | 厚度/间隙 (m) | ID/OD (m) | 径向节点数 |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| 燃料芯块 | `{name}_Pellet` | `0.00400` | `0.00850` | `0.00450` | `0.00800/0.01700` | 5 |
| 裂变气隙 | `pellet_emitter_gap` | `0.00850` | `0.00865` | `0.00015` | - | 简化间隙 |
| 发射极 | `{name}_Emitter` | `0.00865` | `0.00980` | `0.00115` | `0.01730/0.01960` | 1 |
| TEC 间隙 | `tec_couple` | `0.00980` | `0.01030` | `0.00050` | - | 简化间隙 |
| 接收极 | `{name}_Collector` | `0.01030` | `0.01185` | `0.00155` | `0.02060/0.02370` | 1 |
| He 间隙 | `collector_iclad_gap` | `0.01185` | `0.01190` | `0.00005` | - | 简化间隙 |
| 内套管 | `{name}_InnerClad` | `0.01190` | `0.01225` | `0.00035` | `0.02380/0.02450` | 1 |
| 冷却剂环隙 | `Chan_*` | `0.01225` | `0.01295` | `0.00070` | `0.02450/0.02590` | 37 轴向 |
| 外套管 | `{name}_OuterClad` | `0.01295` | `0.01330` | `0.00035` | `0.02590/0.02660` | 1 |
| CO2 间隙 | `oclad_mod_gap` | `0.01330` | `0.01352` | `0.00022` | - | 简化间隙 |
| 虚拟慢化剂 | `{name}_Moderator` | `0.01352` | `0.01627` | `0.00275` | `0.02704/0.03254` | 3 |

其中 `{name}` 为：

```text
Center
Ring1
Ring2
Ring3_TEC
Ring3_Open
```

### 4.2 TFE 固体边界条件

构造阶段每个 `HeatConduction2D` 边界都有默认 `ResistanceBC(R_ext=1e15 K/W)`，等效为占位绝热边界。下表只列物理有效边界。

| 固体边界 | 物理边界条件 | 耦合对象/源项 |
| --- | --- | --- |
| `{name}_Pellet.left` | 内半径绝热 | 默认占位绝热 |
| `{name}_Pellet.right` -> `{name}_Emitter.left` | 裂变气隙导热 + 辐射 | `GapCouple2D`，字典键 `pellet_emitter_gap`，gap `0.00015 m`，`h_eq=5678 W/m2/K` |
| `{name}_Pellet` 体源 | 核热源 | `ReactorCore.update_neutronic_power()` 写入燃料源项 |
| `{name}_Emitter.right` -> `{name}_Collector.left` | TEC 间隙导热/辐射 + 电子热流 | `TECCouple2D`，字典键 `tec_couple`，gap `0.00050 m`，默认 `h_eq=29 W/m2/K` |
| `{name}_Emitter` 体源 | 电极焦耳热 | 来自 ThermoCalc/C++ 结果映射 |
| `{name}_Collector.right` -> `{name}_InnerClad.left` | He 间隙导热 + 辐射 | `GapCouple2D`，字典键 `collector_iclad_gap`，gap `0.00005 m`，`h_eq=5678 W/m2/K` |
| `{name}_Collector` 体源 | 电极焦耳热 | 来自 ThermoCalc/C++ 结果映射 |
| `{name}_InnerClad.right` -> `Chan_{name}` | 流固换热 | `FluidSolidCouple`，名称 `{name}_iclad_coolant_couple`，湿周 `0.0769690200 m` |
| `{name}_OuterClad.left` -> `Chan_{name}` | 流固换热 | `FluidSolidCouple`，名称 `{name}_oclad_coolant_couple`，湿周 `0.0813672497 m` |
| `{name}_OuterClad.right` -> `{name}_Moderator.left` | CO2 间隙导热 + 辐射 | `GapCouple2D`，字典键 `oclad_mod_gap`，gap `0.00022 m`，`h_eq=53.6 W/m2/K` |
| `{name}_Moderator.right` | 到全局慢化剂环的软边界 | `mod_outer_bc`，`R_ext=1/(500*A)`，由 `ReactorCore.pre_step()` 更新并映射热流 |
| 各 TFE 固体 top/bottom | 轴向端面闭合 | 默认占位绝热 |

runner 选择 `local_implicit` 后，具备 `solid_node_capacitance` 的 `FluidSolidCouple` 会切换到局部隐式流固换热时间离散。它改变数值时间离散方式，不改变几何拓扑。

## 5. 全局慢化剂、筒体和反射层

全局慢化剂和外层结构均为圆柱坐标 2D 固体，轴向长度 `0.507 m`，37 个轴向节点。

| 固体对象名 | 材料 | 内半径 (m) | 外半径 (m) | ID/OD (m) | 径向节点数 | 边界条件 |
| --- | --- | ---: | ---: | --- | ---: | --- |
| `TASTIN_Core_V8_CaseA_ModRing_0` | `ZrH` | `0.01627` | `0.0272025` | `0.03254/0.054405` | 3 | 左侧由 TFE 虚拟慢化剂热流映射，右侧固固连接到 `ModRing_1` |
| `TASTIN_Core_V8_CaseA_ModRing_1` | `ZrH` | `0.0272025` | `0.0381350` | `0.054405/0.076270` | 3 | 左右与相邻 `ModRing` 固固连接 |
| `TASTIN_Core_V8_CaseA_ModRing_2` | `ZrH` | `0.0381350` | `0.0490675` | `0.076270/0.098135` | 3 | 左右与相邻 `ModRing` 固固连接 |
| `TASTIN_Core_V8_CaseA_ModRing_3` | `ZrH` | `0.0490675` | `0.0600000` | `0.098135/0.120000` | 3 | 右侧经 `0.005 m` 简化间隙到 `Barrel` |
| `TASTIN_Core_V8_CaseA_Barrel` | `StainlessSteel` | `0.06500` | `0.06800` | `0.13000/0.13600` | 3 | 左侧慢化剂-筒体间隙，右侧 `0.002 m` 间隙到反射层 |
| `TASTIN_Core_V8_CaseA_Reflector` | `BerylliumOxide` | `0.07000` | `0.10200` | `0.14000/0.20400` | 8 | 右边界对 `T_env=200 K` 动态辐射，`emissivity=0.2` |

耦合关系：

- `TASTIN_Core_V8_CaseA_ModRing_0` 到 `TASTIN_Core_V8_CaseA_ModRing_3` 之间使用 `SolidSolidCouple2D` 串联。
- `TASTIN_Core_V8_CaseA_ModRing_3` 到 `TASTIN_Core_V8_CaseA_Barrel` 使用 `GapCouple2D`，间隙 `0.005 m`，`h_eq=0`，辐射发射率 `0.8/0.8`。
- `TASTIN_Core_V8_CaseA_Barrel` 到 `TASTIN_Core_V8_CaseA_Reflector` 使用 `GapCouple2D`，间隙 `0.002 m`，`h_eq=0`，辐射发射率 `0.8/0.8`。
- `TASTIN_Core_V8_CaseA_Reflector.right` 挂 `DynamicRadiationResistanceBC`，对外环境 `200 K`。

## 6. 管翅式辐射器固体耦合

每根 `V12_RadiatorTube_XX` 组件包含：

- 一个水力通道 `V12_RadiatorTubeFluid_XX`。
- 一个管壁固体 `V12_RadiatorTube_XX_Wall`。
- 一个管内流体到管壁内侧的 `FluidSolidCouple`。
- 管壁外侧裸管辐射边界和准稳态翅片等效热阻边界。

| 固体边界 | 物理边界条件 |
| --- | --- |
| `V12_RadiatorTube_XX_Wall.left` | `FluidSolidCouple` 到 `V12_RadiatorTubeFluid_XX`，名称 `V12_RadiatorTube_XX_FluidSolid`，湿周 `pi*0.007 = 0.0219911486 m`，单轴向节点换热面积 `0.0050854531 m2` |
| `V12_RadiatorTube_XX_Wall.right` | 裸管动态辐射到背景温度，默认 `T_bg=3 K`，管发射率 `0.80`，单轴向节点裸管面积 `0.0058119464 m2` |
| `V12_RadiatorTube_XX_Wall.right` | 翅片根部准稳态等效热阻，由 `RadiatorPipeWithFin.pre_step()` 每步更新 |
| `V12_RadiatorTube_XX_Wall.top/bottom` | 无外部物理部件，保留求解器端面闭合条件 |

翅片不是独立 ODE 固体，不计入 `SystemManager.solid_components`，也不提供单独储能项。其热行为由 `RadiatorPipeWithFin` 内部缩阶模型在每个 `pre_step()` 更新为等效边界。

可选遮热罩：

```text
V13_RadiatorThermalShield
```

基础 V13 闭式 builder 默认不挂遮热罩。runner 或启动算例可调用 `attach_radiator_thermal_shield(...)` 添加它。遮热罩是辐射边界修正器，不是注册到全局 ODE 的固体；它通过 `RadiatorPipeWithFin.set_radiation_background_temperature(...)` 改变每根管看到的等效背景温度。

## 7. 固体对象与耦合器数量核对

按默认 V13 构建，`SystemManager.components` 包含：

```text
TASTIN_Core_V8_CaseA
V12_RadiatorTube_01
...
V12_RadiatorTube_78
```

固体对象共 114 个：

- 全局堆芯固体 6 个：`TASTIN_Core_V8_CaseA_ModRing_0` 到 `TASTIN_Core_V8_CaseA_ModRing_3`、`TASTIN_Core_V8_CaseA_Barrel`、`TASTIN_Core_V8_CaseA_Reflector`。
- TFE 代表元件固体 30 个：5 个代表元件，每个 `Pellet/Emitter/Collector/InnerClad/OuterClad/Moderator`。
- 辐射管壁固体 78 个：`V12_RadiatorTube_01_Wall` 到 `V12_RadiatorTube_78_Wall`。

耦合器共 113 个：

- `FluidSolidCouple` 88 个：10 个堆芯冷却剂-套管耦合 + 78 个辐射管内流体-管壁耦合。
- `GapCouple2D` 17 个：15 个 TFE 径向间隙 + 2 个全局慢化剂/筒体/反射层间隙。
- `TECCouple2D` 5 个：每个代表 TFE 的发射极-接收极 TEC 间隙。
- `SolidSolidCouple2D` 3 个：全局慢化剂环之间的径向固固导热。

当前 V13 中真正与固体换热的水力对象只有：

```text
Chan_Center
Chan_Ring1
Chan_Ring2
Chan_Ring3_TEC
Chan_Ring3_Open
V12_RadiatorTubeFluid_01...78
```

其余主回路管段、连接体、泵节点、辐射器上/下集流环段均为水力-only 或压力/混合节点，不挂固体导热对象。

## 8. 默认流量分配口径

默认设计总流量 `1.3 kg/s`：

| 路径 | 默认值 |
| --- | ---: |
| 闭式回路总流量 | `1.3 kg/s` |
| 每根辐射管设计流量 | `1.3 / 78 = 0.0166667 kg/s` |
| 每根真实 TFE 设计流量 | `1.3 / 37 = 0.0351351 kg/s` |
| `V12_CoreInletBranch_1` 宏观支路 | `1.3 / 3 = 0.433333 kg/s` |
| `V12_CoreInletBranch_2_3_Rep` 存储代表支路 | `0.433333 kg/s`，宏观倍率 2 |

## 9. 源码依据

主要依据：

- `testModule/test_core_assemble_v13_caseA.py`
- `testModule/run_v13_caseA_closed_loop.py`
- `testModule/test_core_assemble_v12_caseA.py`
- `testModule/test_core_assemble_v8_caseA.py`
- `testModule/test_core_assemble_v7_caseA.py`
- `Components/RadiatorPipeWithFin.py`
- `Components/RadiatorThermalShield.py`
- `Components/TFEUnit.py`
- `Components/ReactorCore.py`
- `Solvers/Couplers.py`
- `Solvers/Hydrodynamics/Components.py`

核对命令使用项目指定解释器实例化 V13 并读取实际对象字段：

```powershell
@'
import os, sys
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "testModule"))
from testModule.test_core_assemble_v13_caseA import build_v13_case_a_system
b = build_v13_case_a_system()
s = b["system"]
net = s.fluid_solver
print(len(net.volumes_obj), len(net.junctions_obj), len(s.components), len(s.solid_components), len(s.couplers))
print([v.name for v in net.volumes_obj if getattr(v, "is_pressure_reference", False)])
print([v.name for v in net.volumes_obj if getattr(v, "is_pressure_boundary", False)])
'@ | & "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -
```
