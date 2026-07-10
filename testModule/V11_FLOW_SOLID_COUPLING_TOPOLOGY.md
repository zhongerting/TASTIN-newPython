# V11 CaseA 流网与固体耦合拓扑说明

生成日期：2026-06-30

本文对应 `testModule/test_core_assemble_v11_caseA.py` 和默认运行入口 `testModule/run_v11_caseA_closed_loop.py`。命名以源码对象名为准；尺寸只采用 V11 源码、V8/V9 复用常量和 `CoolantLoop/model_collector_ring_6segment_v9_interface.py` 中实际传入对象的值。

![V11 flow and solid coupling topology](./v11_flow_solid_coupling_topology.svg)

## 1. 适用范围

V11 是闭式泵驱动回路：

```text
CoreInletConnector
  -> Chan_Center / Chan_Ring1 / Chan_Ring2 / Chan_Ring3_TEC / Chan_Ring3_Open
  -> CoreOutletConnector
  -> HotOutletBranch_1/2/3
  -> InletMix_I1/I2/I3
  -> A1_I1_to_O1_Channel ... A6_O3_to_I1_Channel
  -> OutletMix_O1/O2/O3
  -> Manifold_1/2/3
  -> V10_RadiatorManifoldMerge
  -> RadiatorInnerHeader_53
  -> RadiatorOuterHeader_52
  -> V11_PumpMidNode
  -> V11_PumpOutletDistributor_51
  -> ColdReturnBranch_1 / ColdReturnBranch_2_3_Rep
  -> V10_ColdReturnOutletMerge
  -> CoreInletConnector
```

压力约束：

- `CoreInletConnector.is_pressure_reference=True`，`target_P=166471.52 Pa`。
- V11 不含 `is_pressure_boundary=True` 的固定压力流动边界体。
- `J_RadiatorOuterHeader_52_to_PumpA` 和 `J_PumpA_to_PumpOutletDistributor_51` 是两台串联相同 `PumpJunction`，默认总压头 `6466.56 Pa`，每台 `3233.28 Pa`。

默认运行入口参数：

| 参数 | 默认值 | 说明 |
| --- | ---: | --- |
| `--total-inlet-flow-kg-s` | `1.3 kg/s` | 设计总质量流量，用于初始化与诊断 |
| `--target-flow-kg-s` | `1.3 kg/s` | 可选泵压头控制目标 |
| `--pump-total-head-pa` | `6466.56 Pa` | 两台泵合计压头 |
| `--inlet-temperature-k` | `727.0 K` | 仅用于构建/注入时的初始温度，不是闭式固定温度边界 |
| `--ring-emissivity` | `0.2` | 覆盖集流环壁面辐射率 |
| `--hp-emissivity` | `0.75` | 覆盖热管裸管辐射率 |
| `--fin-emissivity` | `0.75` | 覆盖翅片辐射率 |
| `--outer-header-emissivity` | `0.2` | `RadiatorOuterHeader_52_RadiationSink` 直接对流体扣热 |
| `--solid-ode-method` | `RK45` | 运行器默认会从公共默认值切换到 `RK45` |
| `--fluid-solid-coupling-scheme` | `local_implicit` | 运行器默认会从 `current` 切换到 `local_implicit` |

注意：裸调用 `build_v11_case_a_system()` 且不传辐射率时，会使用集流环接口算例旧默认值 `RING_EMISSIVITY=0.05`、`HP_EMISSIVITY=0.6`、`FIN_EMISSIVITY=0.6`，且 `outer_header_emissivity=0.0`。本文的运行工况说明以 `run_v11_caseA_closed_loop.py` 默认参数为准。

## 2. 水力体与几何尺寸

表中 `ID` 对圆管取 `hydraulic_diam`。若对象只建模了流动面积和水力直径，没有壁厚或外径，则 `OD` 标为 `not modeled`。这不是遗漏，而是 V11 水力对象没有该字段。矩形集流环同时列出矩形尺寸和等效圆柱固体尺寸。

| 对象名 | 类型/倍率 | 节点数 | 长度 L (m) | 流通面积 A (m2) | ID 或等效内径 (m) | OD 或等效外径 (m) | 水力直径 Dh (m) | 备注 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `CoreInletConnector` | `IncompressibleFluidVolume`, 压力参考 | 1 | `0.02` 构造参数 | `0.002050454693` | not modeled | not modeled | `0.0014` | 体积 `1.0e-5 m3`，对象不保存 length 字段 |
| `CoreOutletConnector` | `IncompressibleFluidVolume` | 1 | `0.02` 构造参数 | `0.002050454693` | not modeled | not modeled | `0.0014` | 体积 `1.0e-5 m3` |
| `Chan_Center` | TFE 冷却剂环隙，倍率 1 | 37 | `0.507` | `5.541769441e-5` | `0.02450` | `0.02590` | `0.0014` | 下/上段各 6 节点，活性段 25 节点 |
| `Chan_Ring1` | TFE 冷却剂环隙，倍率 6 | 37 | `0.507` | `5.541769441e-5` | `0.02450` | `0.02590` | `0.0014` | 同上 |
| `Chan_Ring2` | TFE 冷却剂环隙，倍率 12 | 37 | `0.507` | `5.541769441e-5` | `0.02450` | `0.02590` | `0.0014` | 同上 |
| `Chan_Ring3_TEC` | TFE 冷却剂环隙，倍率 15 | 37 | `0.507` | `5.541769441e-5` | `0.02450` | `0.02590` | `0.0014` | TEC 倍率 15 |
| `Chan_Ring3_Open` | TFE 冷却剂环隙，倍率 3 | 37 | `0.507` | `5.541769441e-5` | `0.02450` | `0.02590` | `0.0014` | 默认 TEC 倍率 0，主动 TEC 源清零 |
| `HotOutletBranch_1` | 圆管 | 8 | `2.19632` | `0.00059798` | `0.0276` | not modeled | `0.0276` | 从 `CoreOutletConnector` 到 `InletMix_I1` |
| `HotOutletBranch_2` | 圆管 | 8 | `2.19632` | `0.00059798` | `0.0276` | not modeled | `0.0276` | 到 `InletMix_I2` |
| `HotOutletBranch_3` | 圆管 | 8 | `2.19632` | `0.00059798` | `0.0276` | not modeled | `0.0276` | 到 `InletMix_I3` |
| `InletMix_I1` | 混合节点 | 1 | `0.0552` 推导 | `0.00059798` | `0.0276` | not modeled | `0.0276` | 体积 `3.3008496e-5 m3`，长度 `2*Dh` |
| `InletMix_I2` | 混合节点 | 1 | `0.0552` 推导 | `0.00059798` | `0.0276` | not modeled | `0.0276` | 同上 |
| `InletMix_I3` | 混合节点 | 1 | `0.0552` 推导 | `0.00059798` | `0.0276` | not modeled | `0.0276` | 同上 |
| `A1_I1_to_O1_Channel` | 矩形集流环通道 | 3 | `0.793` | `0.0044` | `0.110 x 0.040` 矩形 | 壁厚 `0.002` | `0.058666667` | 等效固体 ID/OD `0.095493/0.099493` |
| `A2_O1_to_I2_Channel` | 矩形集流环通道 | 3 | `0.793` | `0.0044` | `0.110 x 0.040` 矩形 | 壁厚 `0.002` | `0.058666667` | 同上 |
| `A3_I2_to_O2_Channel` | 矩形集流环通道 | 3 | `0.793` | `0.0044` | `0.110 x 0.040` 矩形 | 壁厚 `0.002` | `0.058666667` | 同上 |
| `A4_O2_to_I3_Channel` | 矩形集流环通道 | 3 | `0.793` | `0.0044` | `0.110 x 0.040` 矩形 | 壁厚 `0.002` | `0.058666667` | 同上 |
| `A5_I3_to_O3_Channel` | 矩形集流环通道 | 3 | `0.793` | `0.0044` | `0.110 x 0.040` 矩形 | 壁厚 `0.002` | `0.058666667` | 同上 |
| `A6_O3_to_I1_Channel` | 矩形集流环通道 | 3 | `0.793` | `0.0044` | `0.110 x 0.040` 矩形 | 壁厚 `0.002` | `0.058666667` | 同上 |
| `OutletMix_O1` | 混合节点 | 1 | `0.036` 推导 | `0.00025434` | `0.018` | not modeled | `0.018` | 体积 `9.15624e-6 m3`，长度 `2*Dh` |
| `OutletMix_O2` | 混合节点 | 1 | `0.036` 推导 | `0.00025434` | `0.018` | not modeled | `0.018` | 同上 |
| `OutletMix_O3` | 混合节点 | 1 | `0.036` 推导 | `0.00025434` | `0.018` | not modeled | `0.018` | 同上 |
| `Manifold_1` | 圆管 | 5 | `0.40911` | `0.00025434` | `0.018` | not modeled | `0.018` | 到 `V10_RadiatorManifoldMerge` |
| `Manifold_2` | 圆管 | 17 | `1.41912` | `0.00025434` | `0.018` | not modeled | `0.018` | 到 `V10_RadiatorManifoldMerge` |
| `Manifold_3` | 圆管 | 17 | `1.41912` | `0.00025434` | `0.018` | not modeled | `0.018` | 到 `V10_RadiatorManifoldMerge` |
| `V10_RadiatorManifoldMerge` | 汇合体 | 1 | `0.02` 构造参数 | `0.001734` | `0.047` | not modeled | `0.047` | 体积 `1.0e-5 m3` |
| `RadiatorInnerHeader_53` | 圆管 | 18 | `1.50969` | `0.001734` | `0.047` | not modeled | `0.047` | 节点长 `0.083871667 m` |
| `RadiatorOuterHeader_52` | 圆管 | 1 | `0.0915` | `0.001734` | `0.047` | not modeled | `0.047` | 可挂 `RadiatorOuterHeader_52_RadiationSink` |
| `V11_PumpMidNode` | 泵间节点 | 1 | `0.02` 构造参数 | `0.001734` | `0.047` | not modeled | `0.047` | 体积 `1.0e-5 m3` |
| `V11_PumpOutletDistributor_51` | 泵出口分流节点 | 1 | `0.02` 构造参数 | `0.001734` | `0.047` | not modeled | `0.047` | 体积 `1.0e-5 m3` |
| `ColdReturnBranch_1` | 圆管，显式 1 根 | 8 | `1.89021` | `0.0005979816` | `0.0276` | not modeled | `0.0276` | 接回 `V10_ColdReturnOutletMerge` |
| `ColdReturnBranch_2_3_Rep` | 圆管，宏观代表 2 根 | 8 | `2.50705` | `0.0005979816` | `0.0276` | not modeled | `0.0276` | 入口/出口均用 `MacroFlowJunction(multiplier=2)` |
| `V10_ColdReturnOutletMerge` | 汇合体 | 1 | `0.02` 构造参数 | `0.00179394` | not modeled | not modeled | `0.0276` | 体积 `1.0e-5 m3`，`A=3*AREA_CORE_BRANCH` |

内部节点命名规则：

- 每个 `IncompressibleFluidChannel` 内部体积命名为 `{ChannelName}_Vol_XX`。
- 每个通道内部连接命名为 `{ChannelName}_Junc_XX_YY`。
- 堆芯入口/出口旧连接保留 V8/V7 风格：`J_PlenumIn_{name}`、`J_PlenumOut_{name}`，其中 `{name}` 为 `Center`、`Ring1`、`Ring2`、`Ring3_TEC`、`Ring3_Open`。

## 3. 外部主连接

| 连接名 | 类型 | from -> to | 面积 (m2) | K | 倍率/压头 |
| --- | --- | --- | ---: | ---: | --- |
| `J_CoreOutletConnector_to_HotOutletBranch_1` | `FlowJunction` | `CoreOutletConnector` -> `HotOutletBranch_1_Vol_01` | `0.00059798` | `0.0` | - |
| `J_CoreOutletConnector_to_HotOutletBranch_2` | `FlowJunction` | `CoreOutletConnector` -> `HotOutletBranch_2_Vol_01` | `0.00059798` | `0.0` | - |
| `J_CoreOutletConnector_to_HotOutletBranch_3` | `FlowJunction` | `CoreOutletConnector` -> `HotOutletBranch_3_Vol_01` | `0.00059798` | `0.0` | - |
| `J_HotOutletBranch_1_to_InletMix_I1` | `MacroFlowJunction` | `HotOutletBranch_1_Vol_08` -> `InletMix_I1` | `0.00059798` | `1.0` | `multiplier=2` |
| `J_HotOutletBranch_2_to_InletMix_I2` | `MacroFlowJunction` | `HotOutletBranch_2_Vol_08` -> `InletMix_I2` | `0.00059798` | `1.0` | `multiplier=2` |
| `J_HotOutletBranch_3_to_InletMix_I3` | `MacroFlowJunction` | `HotOutletBranch_3_Vol_08` -> `InletMix_I3` | `0.00059798` | `1.0` | `multiplier=2` |
| `J_I1_to_A1_I1_to_O1` | `FlowJunction` | `InletMix_I1` -> `A1_I1_to_O1_Channel_Vol_01` | `0.0044` | `0.5` | - |
| `J_A1_I1_to_O1_to_O1` | `FlowJunction` | `A1_I1_to_O1_Channel_Vol_03` -> `OutletMix_O1` | `0.0044` | `0.7` | 含热管阵列动态阻力参数 |
| `J_O1_to_A2_O1_to_I2` | `FlowJunction` | `OutletMix_O1` -> `A2_O1_to_I2_Channel_Vol_01` | `0.0044` | `0.5` | - |
| `J_A2_O1_to_I2_to_I2` | `FlowJunction` | `A2_O1_to_I2_Channel_Vol_03` -> `InletMix_I2` | `0.0044` | `0.7` | 含热管阵列动态阻力参数 |
| `J_I2_to_A3_I2_to_O2` | `FlowJunction` | `InletMix_I2` -> `A3_I2_to_O2_Channel_Vol_01` | `0.0044` | `0.5` | - |
| `J_A3_I2_to_O2_to_O2` | `FlowJunction` | `A3_I2_to_O2_Channel_Vol_03` -> `OutletMix_O2` | `0.0044` | `0.7` | 含热管阵列动态阻力参数 |
| `J_O2_to_A4_O2_to_I3` | `FlowJunction` | `OutletMix_O2` -> `A4_O2_to_I3_Channel_Vol_01` | `0.0044` | `0.5` | - |
| `J_A4_O2_to_I3_to_I3` | `FlowJunction` | `A4_O2_to_I3_Channel_Vol_03` -> `InletMix_I3` | `0.0044` | `0.7` | 含热管阵列动态阻力参数 |
| `J_I3_to_A5_I3_to_O3` | `FlowJunction` | `InletMix_I3` -> `A5_I3_to_O3_Channel_Vol_01` | `0.0044` | `0.5` | - |
| `J_A5_I3_to_O3_to_O3` | `FlowJunction` | `A5_I3_to_O3_Channel_Vol_03` -> `OutletMix_O3` | `0.0044` | `0.7` | 含热管阵列动态阻力参数 |
| `J_O3_to_A6_O3_to_I1` | `FlowJunction` | `OutletMix_O3` -> `A6_O3_to_I1_Channel_Vol_01` | `0.0044` | `0.5` | - |
| `J_A6_O3_to_I1_to_I1` | `FlowJunction` | `A6_O3_to_I1_Channel_Vol_03` -> `InletMix_I1` | `0.0044` | `0.7` | 含热管阵列动态阻力参数 |
| `J_OutletMix_O1_Manifold_1` | `FlowJunction` | `OutletMix_O1` -> `Manifold_1_Vol_01` | `0.00025434` | `1.0` | - |
| `J_OutletMix_O2_Manifold_2` | `FlowJunction` | `OutletMix_O2` -> `Manifold_2_Vol_01` | `0.00025434` | `1.0` | - |
| `J_OutletMix_O3_Manifold_3` | `FlowJunction` | `OutletMix_O3` -> `Manifold_3_Vol_01` | `0.00025434` | `1.0` | - |
| `J_Manifold_1_to_RadiatorMerge` | `MacroFlowJunction` | `Manifold_1_Vol_05` -> `V10_RadiatorManifoldMerge` | `0.00025434` | `1.1` | `multiplier=2` |
| `J_Manifold_2_to_RadiatorMerge` | `MacroFlowJunction` | `Manifold_2_Vol_17` -> `V10_RadiatorManifoldMerge` | `0.00025434` | `1.1` | `multiplier=2` |
| `J_Manifold_3_to_RadiatorMerge` | `MacroFlowJunction` | `Manifold_3_Vol_17` -> `V10_RadiatorManifoldMerge` | `0.00025434` | `1.1` | `multiplier=2` |
| `J_RadiatorMerge_to_RadiatorInnerHeader_53` | `FlowJunction` | `V10_RadiatorManifoldMerge` -> `RadiatorInnerHeader_53_Vol_01` | `0.001734` | `0.0` | - |
| `J_RadiatorInnerHeader_53_to_RadiatorOuterHeader_52` | `FlowJunction` | `RadiatorInnerHeader_53_Vol_18` -> `RadiatorOuterHeader_52_Vol_01` | `0.001734` | `0.0` | - |
| `J_RadiatorOuterHeader_52_to_PumpA` | `PumpJunction` | `RadiatorOuterHeader_52_Vol_01` -> `V11_PumpMidNode` | `0.001734` | `0.0` | `delta_p=3233.28 Pa` |
| `J_PumpA_to_PumpOutletDistributor_51` | `PumpJunction` | `V11_PumpMidNode` -> `V11_PumpOutletDistributor_51` | `0.001734` | `0.0` | `delta_p=3233.28 Pa` |
| `J_PumpOutletDistributor_51_to_ColdReturnBranch_1` | `FlowJunction` | `V11_PumpOutletDistributor_51` -> `ColdReturnBranch_1_Vol_01` | `0.00059798` | `0.0` | - |
| `J_PumpOutletDistributor_51_to_ColdReturnBranch_2_3_Rep` | `MacroFlowJunction` | `V11_PumpOutletDistributor_51` -> `ColdReturnBranch_2_3_Rep_Vol_01` | `0.00059798` | `0.0` | `multiplier=2` |
| `J_ColdReturnBranch_1_to_ColdReturnOutletMerge` | `FlowJunction` | `ColdReturnBranch_1_Vol_08` -> `V10_ColdReturnOutletMerge` | `0.00059798` | `0.0` | - |
| `J_ColdReturnBranch_2_3_Rep_to_ColdReturnOutletMerge` | `MacroFlowJunction` | `ColdReturnBranch_2_3_Rep_Vol_08` -> `V10_ColdReturnOutletMerge` | `0.00059798` | `0.0` | `multiplier=2` |
| `J_ColdReturnOutletMerge_to_CoreInletConnector` | `FlowJunction` | `V10_ColdReturnOutletMerge` -> `CoreInletConnector` | `0.00179394` | `0.0` | - |

## 4. 堆芯固体与边界条件

### 4.1 TFE 代表元件尺寸

每个代表 TFE 的轴向长度均为 `0.507 m`：

- 下反射/包层段：`0.065 m`，6 节点，每节点 `0.010833333 m`。
- 活性段：`0.377 m`，25 节点，每节点 `0.01508 m`。
- 上反射/包层段：`0.065 m`，6 节点，每节点 `0.010833333 m`。

径向几何如下，所有代表元件相同：

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

构造阶段每个 `HeatConduction2D` 边界都有一个默认 `ResistanceBC(R_ext=1e15 K/W)`，等价于占位绝热边界。下表列的是物理有效边界；默认占位不再重复列出。

| 固体边界 | 物理边界条件 | 耦合对象/源项 |
| --- | --- | --- |
| `{name}_Pellet.left` | 内半径绝热 | 默认占位绝热 |
| `{name}_Pellet.right` -> `{name}_Emitter.left` | 裂变气隙导热 + 辐射 | `GapCouple2D`，字典键 `pellet_emitter_gap` |
| `{name}_Pellet.top/bottom` | 轴向绝热 | 默认占位绝热 |
| `{name}_Pellet` 体源 | 核热源 | `ReactorCore` 按代表倍率和功率因子写入 `Fuel.set_nuclear_power()` |
| `{name}_Emitter.left` | 裂变气隙侧 | `pellet_emitter_gap` |
| `{name}_Emitter.right` -> `{name}_Collector.left` | TEC 间隙导热/辐射 + 电子热流 | `TECCouple2D`，字典键 `tec_couple` |
| `{name}_Emitter` 体源 | 电极焦耳热 | `ReactorCore` 从 ThermoCalc/C++ 结果映射 |
| `{name}_Collector.left` | TEC 间隙侧 | `tec_couple` |
| `{name}_Collector.right` -> `{name}_InnerClad.left` | He 间隙导热 + 辐射 | `GapCouple2D`，字典键 `collector_iclad_gap` |
| `{name}_Collector` 体源 | 电极焦耳热 | `ReactorCore` 从 ThermoCalc/C++ 结果映射 |
| `{name}_InnerClad.left` | He 间隙侧 | `collector_iclad_gap` |
| `{name}_InnerClad.right` -> `Chan_{name}` | 流固换热 | `FluidSolidCouple`，名称 `{name}_iclad_coolant_couple`，湿周 `0.0769690200 m` |
| `{name}_OuterClad.left` -> `Chan_{name}` | 流固换热 | `FluidSolidCouple`，名称 `{name}_oclad_coolant_couple`，湿周 `0.0813672497 m` |
| `{name}_OuterClad.right` -> `{name}_Moderator.left` | CO2 间隙导热 + 辐射 | `GapCouple2D`，字典键 `oclad_mod_gap` |
| `{name}_Moderator.right` | 与全局慢化剂环的软边界 | `mod_outer_bc`，`R_soft=1/(500*A_surf)`，`T_ext` 每步由 `ReactorCore.pre_step()` 更新 |
| `{name}_Emitter.top/bottom`、`{name}_Collector.top/bottom`、`{name}_InnerClad.top/bottom`、`{name}_OuterClad.top/bottom`、`{name}_Moderator.top/bottom` | 轴向绝热 | 默认占位绝热 |

运行器选择 `local_implicit` 后，所有具有 `solid_node_capacitance` 的 `FluidSolidCouple` 会在对应固体边界额外挂接局部隐式热流边界。它改变时间离散方式，不改变几何拓扑。

## 5. 全局慢化剂、筒体和反射层

全局慢化剂和外层结构均为圆柱坐标 2D 固体，轴向长度 `0.507 m`，37 个轴向节点。

| 固体对象名 | 材料 | 内半径 (m) | 外半径 (m) | ID/OD (m) | 径向节点数 | 边界条件 |
| --- | --- | ---: | ---: | --- | ---: | --- |
| `TASTIN_Core_V8_CaseA_ModRing_0` | `ZrH` | `0.01627` | `0.0272025` | `0.03254/0.054405` | 3 | 左侧由 TFE 虚拟慢化剂热流映射，右侧固固连接到 `ModRing_1` |
| `TASTIN_Core_V8_CaseA_ModRing_1` | `ZrH` | `0.0272025` | `0.038135` | `0.054405/0.07627` | 3 | 左右与相邻 `ModRing` 固固连接 |
| `TASTIN_Core_V8_CaseA_ModRing_2` | `ZrH` | `0.038135` | `0.0490675` | `0.07627/0.098135` | 3 | 左右与相邻 `ModRing` 固固连接 |
| `TASTIN_Core_V8_CaseA_ModRing_3` | `ZrH` | `0.0490675` | `0.06000` | `0.098135/0.12000` | 3 | 右侧经 `5.0e-3 m` 简化间隙到 `Barrel` |
| `TASTIN_Core_V8_CaseA_Barrel` | `StainlessSteel` | `0.06500` | `0.06800` | `0.13000/0.13600` | 3 | 左侧经慢化剂-筒体间隙，右侧经 `2.0e-3 m` 间隙到反射层 |
| `TASTIN_Core_V8_CaseA_Reflector` | `BerylliumOxide` | `0.07000` | `0.10200` | `0.14000/0.20400` | 8 | 右边界对 `T_space=200 K` 动态辐射，`emissivity=0.2` |

耦合关系：

- `TASTIN_Core_V8_CaseA_ModRing_0` 到 `TASTIN_Core_V8_CaseA_ModRing_3` 之间使用 `SolidSolidCouple2D` 串联。
- `TASTIN_Core_V8_CaseA_ModRing_3` 到 `TASTIN_Core_V8_CaseA_Barrel` 使用 `GapCouple2D`，间隙宽度 `0.005 m`，简化气隙导热系数 `h_eq=0`，辐射发射率 `0.8/0.8`。
- `TASTIN_Core_V8_CaseA_Barrel` 到 `TASTIN_Core_V8_CaseA_Reflector` 使用 `GapCouple2D`，间隙宽度 `0.002 m`，简化气隙导热系数 `h_eq=0`，辐射发射率 `0.8/0.8`。
- `TASTIN_Core_V8_CaseA_Reflector.right` 挂 `DynamicRadiationResistanceBC`，对外环境 `200 K`。

## 6. 集流环、热管和翅片固体耦合

### 6.1 集流环段

V11 接入单套显式 6 段集流环，并在 `HotOutletBranch_*` 到 `InletMix_*`、`Manifold_*` 到 `V10_RadiatorManifoldMerge` 的接口用 `MacroFlowJunction(multiplier=2)` 表示第二套对称集流环。

每个集流环段都包含：

- 一个水力通道 `{segment}_Channel`。
- 一个壁面固体 `{segment}_Solid`。
- 一个组件 `{segment}_RingHP`，负责 `{segment}_Channel` 到 `{segment}_Solid` 的换热，以及每个通道节点上的代表热管换热。

| segment | 流向 | 热管倍率 per node | 壁面固体对象 | RingHP 组件 |
| --- | --- | --- | --- | --- |
| `A1_I1_to_O1` | `I1 -> O1` | `[5, 6, 6]` | `A1_I1_to_O1_Solid` | `A1_I1_to_O1_RingHP` |
| `A2_O1_to_I2` | `O1 -> I2` | `[5, 5, 6]` | `A2_O1_to_I2_Solid` | `A2_O1_to_I2_RingHP` |
| `A3_I2_to_O2` | `I2 -> O2` | `[5, 6, 6]` | `A3_I2_to_O2_Solid` | `A3_I2_to_O2_RingHP` |
| `A4_O2_to_I3` | `O2 -> I3` | `[5, 5, 6]` | `A4_O2_to_I3_Solid` | `A4_O2_to_I3_RingHP` |
| `A5_I3_to_O3` | `I3 -> O3` | `[5, 6, 6]` | `A5_I3_to_O3_Solid` | `A5_I3_to_O3_RingHP` |
| `A6_O3_to_I1` | `O3 -> I1` | `[5, 6, 6]` | `A6_O3_to_I1_Solid` | `A6_O3_to_I1_RingHP` |

集流环壁面固体几何：

| 项 | 值 |
| --- | ---: |
| 矩形流道尺寸 | `0.110 m x 0.040 m` |
| 流通面积 | `0.0044 m2` |
| 湿周 | `0.300 m` |
| 水力直径 | `0.058666667 m` |
| 等效固体内半径 | `0.0477464829 m` |
| 等效固体外半径 | `0.0497464829 m` |
| 等效固体 ID/OD | `0.095492966/0.099492966 m` |
| 壁厚 | `0.002 m` |
| 每段长度 | `0.793 m` |
| 每段固体网格 | 径向 1 节点，轴向 3 节点 |

集流环壁面边界条件：

| 固体边界 | 物理边界条件 |
| --- | --- |
| `{segment}_Solid.left` | `FluidSolidCouple` 到 `{segment}_Channel`，名称 `{segment}_RingHP_coupler_header`，湿周 `0.300 m` |
| `{segment}_Solid.right` | 对 `T_SPACE=200 K` 动态辐射；运行器默认 `emissivity=0.2` |
| `{segment}_Solid.top/bottom` | 轴向绝热占位 |

### 6.2 热管和翅片

每个 `{segment}_RingHP` 在三个集流环节点上各建立一个 `HPwithFin` 代表热管，固体对象名：

```text
{segment}_RingHP_HP_node0_HP_inner
{segment}_RingHP_HP_node1_HP_inner
{segment}_RingHP_HP_node2_HP_inner
```

热管几何：

| 项 | 值 |
| --- | ---: |
| 蒸汽腔半径 `R_VAPOR_HP` | `0.0075 m` |
| 吸液芯厚度 | `0.0005 m` |
| 管壁厚度 | `0.0010 m` |
| 管壁内半径 `R_IN_HP` | `0.0080 m` |
| 管壁外半径 `R_OUT_HP` | `0.0090 m` |
| 蒸汽腔 ID | `0.0150 m` |
| 热管 OD | `0.0180 m` |
| 蒸发段 `L_EVA` | `0.100 m` |
| 绝热段 `L_ABA` | `0.0 m` |
| 冷凝段 `L_CON` | `0.500 m` |
| 总长 | `0.600 m` |
| 轴向节点 | 蒸发段 1，冷凝段 12，总 13 |
| 径向节点 | 吸液芯 1，管壁 2，总 3 |

翅片等效几何：

| 项 | 值 |
| --- | ---: |
| `THIN_FIN` | `0.0004 m` |
| `FIN_HEIGHT` | `0.020 m` |
| `N_FIN_HEIGHT` | 15 |
| `fin_wrap_ratio` | `0.0141471` |
| 裸管运行器发射率 | `0.75` |
| 翅片运行器发射率 | `0.75` |
| `up_view_factor/down_view_factor` | `0.0 / 0.3` |
| 裸管等效辐射发射率 | `0.75*(1+0.3)/2 = 0.4875` |
| 翅片等效辐射发射率 | `0.75*(1+0.3)/2 = 0.4875` |

热管边界条件：

| 固体边界 | 物理边界条件 |
| --- | --- |
| `{segment}_RingHP_HP_node{i}_HP_inner.outer_eva` | `FluidSolidCouple` 到对应集流环流体节点，名称 `{segment}_RingHP_coupler_hp_{i}`，湿周 `2*pi*R_OUT_HP = 0.0565486678 m` |
| `{segment}_RingHP_HP_node{i}_HP_inner.outer_aba` | 长度为 0，边界数组为空；源码仍挂动态辐射条件但物理贡献为 0 |
| `{segment}_RingHP_HP_node{i}_HP_inner.outer_con` | 裸管动态辐射到 `T_SPACE=200 K`，并有翅片等效热阻边界，`HPwithFin.pre_step()` 每步更新 |
| `{segment}_RingHP_HP_node{i}_HP_inner.left` | 蒸汽腔内侧绝热占位 |
| `{segment}_RingHP_HP_node{i}_HP_inner.top/bottom` | 轴向端面绝热占位 |

`RadiatorOuterHeader_52_RadiationSink` 不是固体边界条件。它是 `FluidChannelRadiationSink`，运行器默认 `outer_header_emissivity=0.2` 时直接按 `RadiatorOuterHeader_52` 流体温度扣除辐射热，周长取 `pi*DH_HEADER`，环境温度默认 `T_SPACE=200 K`。

## 7. 固体对象与耦合器数量核对

按默认 V11 构建，`SystemManager.components` 包含：

```text
TASTIN_Core_V8_CaseA
A1_I1_to_O1_RingHP
A2_O1_to_I2_RingHP
A3_I2_to_O2_RingHP
A4_O2_to_I3_RingHP
A5_I3_to_O3_RingHP
A6_O3_to_I1_RingHP
```

固体对象共 60 个：

- 全局堆芯固体 6 个：`TASTIN_Core_V8_CaseA_ModRing_0` 到 `TASTIN_Core_V8_CaseA_ModRing_3`、`TASTIN_Core_V8_CaseA_Barrel`、`TASTIN_Core_V8_CaseA_Reflector`。
- TFE 代表元件固体 30 个：5 个代表元件，每个 `Pellet/Emitter/Collector/InnerClad/OuterClad/Moderator`。
- 集流环壁面固体 6 个：`A1_I1_to_O1_Solid` 到 `A6_O3_to_I1_Solid`。
- 热管固体 18 个：6 段 x 3 个 `HP_node{i}_HP_inner`。

耦合器共 59 个：

- 全局慢化剂环间和外层结构耦合：5 个。
- TFE 内部每个代表元件 6 个，共 30 个：裂变气隙、TEC 间隙、He 间隙、内套管-冷却剂、外套管-冷却剂、CO2 间隙。
- 集流环每段 4 个，共 24 个：1 个集流环壁面-流体换热 + 3 个热管蒸发段-流体换热。

## 8. 源码依据

主要依据：

- `testModule/test_core_assemble_v11_caseA.py`
- `testModule/run_v11_caseA_closed_loop.py`
- `testModule/test_core_assemble_v8_caseA.py`
- `testModule/test_core_assemble_v9_caseA.py`
- `testModule/test_core_assemble_v7_caseA.py`
- `CoolantLoop/model_collector_ring_6segment_v9_interface.py`
- `CoolantLoop/model_collector_ring_full_ringhp_geometry100hp_potassium_mixed.py`
- `Components/TFEUnit.py`
- `Components/ReactorCore.py`
- `Components/RingHP.py`
- `Components/HPwithFin.py`
- `Solvers/Couplers.py`
- `Solvers/Hydrodynamics/Components.py`

核对命令使用指定解释器实例化 V11 并读取实际对象字段：

```powershell
@'
import os, sys
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), "testModule"))
from testModule.test_core_assemble_v11_caseA import build_v11_case_a_system
b = build_v11_case_a_system(
    inlet_temperature_k=727.0,
    ring_emissivity=0.2,
    hp_emissivity=0.75,
    fin_emissivity=0.75,
    outer_header_emissivity=0.2,
)
'@ | & "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -
```
