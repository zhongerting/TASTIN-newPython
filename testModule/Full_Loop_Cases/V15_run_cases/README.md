# V15 运行算例说明

本目录保存 V15 管翅式辐射管辐射器 TOPAZ-II 系统的可运行 smoke 算例和已检查结果文件。

## 算例定位

V15 的目标是脱离 V13/V12 历史测试算例，按物理语义重新装配管翅式辐射器闭式回路。V15 只复用共用层提供的堆芯、堆芯进出口连接件、辐射器进出口接口和两台串联主泵；辐射器入口分配、78 根辐射管、上下集流环环向流动、辐射器内外总管、泵出口分配器和三条冷回流支路都在 V15 adapter 中本地声明。

## 拓扑边界

共用层提供：

- `CoreInletConnector`
- `CoreOutletConnector`
- `RadiatorInletHeader`
- `RadiatorOutletHeader`，在 V15 adapter 中重命名为 `RadiatorOuterHeader`
- `J_PumpA -> PumpMidNode -> J_PumpB -> PumpOutletNode`

V15 调用共用层时使用 `connect_pump_outlet_to_core=False`，因此不会生成 `CoreInletSegment`，也不会把 `PumpOutletNode` 直接连接到堆芯入口。

V15 本地构建：

- `RadiatorInletHeader -> RadiatorInletDistributor`
- `RadiatorUpperHeader_01...78`，包含完整环向连接
- `RadiatorTubeFluid_01...78`，每根管对应一个 `RadiatorPipeWithFin` 固体/翅片等效模型
- `RadiatorLowerHeader_01...78`，包含完整环向连接
- `RadiatorInnerHeader -> RadiatorOuterHeader -> PumpA -> PumpMidNode -> PumpB -> PumpOutletNode`
- `PumpOutletNode -> PumpOutletDistributor -> ColdReturnBranch_1/2/3 -> CoreInletConnector`

## 堆芯三进三出说明

V15 当前没有使用共用层的 `CoreInletSegment`。三条堆芯入口回流支路由 `V15PipeFinRadiatorConfig.cold_return_branch_*` 控制，默认尺寸为：

- `ColdReturnBranch_1/2/3`: `length=1.89021 m`, `area=pi*0.0138^2 m2`, `Dh=0.0276 m`, `n_nodes=1`

V15 的这三条冷回流支路与 V14 的 `HotOutletBranch_1/2/3` 截面积和水力直径相同，但长度和节点数不同。若后续确认三进三出连接段应作为 V14/V15 共同物理部件，需要把这些参数上提到共用层配置并统一拓扑边界。

## 命名约束

V15 不再使用 V13/V12 的历史遗留对象名和拆分方式，禁止重新引入：

- `V12_` 前缀
- `Pipe05`
- `Pipe06`
- `Pipe07`
- `Pipe08`
- `Pipe09`
- `Pipe11`

## 水力-only smoke

结果文件：`v15_flow_path_smoke_result.json`

生成命令：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases.V15_run_cases.test_v15_flow_path_smoke
```

测试范围：

- 构建 V15 共用堆芯/主泵层和本地管翅式辐射器。
- 将管壁和翅片发射率设为 `0.0`。
- 不调用 `SystemManager.initialize_system()`。
- 不执行热耦合器。
- 不推进任何固体导热模型。
- 只运行 `HydraulicNetwork.initialize_hydraulics()` 和一次 `step_hydraulic()`。

当前 smoke 参数：

- `hydraulic_init_dt_s = 0.005`
- `hydraulic_tol_kg_s = 1.0e-4`
- `hydraulic_max_iter = 1500`
- `hydraulic_step_dt_s = 1.0e-4`

当前结果摘要：

- `radiator_tube_count = 78`
- `cold_return_branch_count = 3`
- `n_volumes = 977`
- `n_junctions = 1064`
- 唯一压力参考为 `CoreInletConnector`
- 无固定压力边界
- 总泵压头 `6466.56 Pa`

该测试只验证闭式流网连通性和有限流量/压力，不代表热性能或稳态结果。