# V14 运行算例说明

本目录保存 V14 热管辐射器 TOPAZ-II 系统的可运行 smoke 算例和已检查结果文件。

## 算例定位

V14 的目标是脱离 V11 历史测试算例，重新在 `Full_Loop_Cases` 中装配热管辐射器闭式回路。V14 只复用共用层提供的堆芯、堆芯进出口连接件、辐射器进出口接口和两台串联主泵；热管辐射器、热端三出支路、集流环、出口混合节点和回流 manifold 都在 V14 adapter 中本地声明。

## 拓扑边界

共用层提供：

- `CoreInletConnector`
- `CoreOutletConnector`
- `RadiatorInletHeader`
- `RadiatorOutletHeader`
- `CoreInletSegment`
- `J_PumpA -> PumpMidNode -> J_PumpB -> PumpOutletNode`

V14 本地构建：

- `RadiatorInletHeader -> HotOutletBranch_1/2/3 -> InletMix_I1/I2/I3`
- 单个显式热管集流环：`I1 -> A1 -> O1 -> A2 -> I2 -> A3 -> O2 -> A4 -> I3 -> A5 -> O3 -> A6 -> I1`
- `OutletMix_O1/O2/O3 -> Manifold_1/2/3 -> RadiatorOutletHeader`

当前 V14 用 `MacroFlowJunction(multiplier=2)` 表示第二套对称集流环，因此热端三出到集流环、以及出口 manifold 回到 `RadiatorOutletHeader` 的宏观流量与显式单环流量不同。

## 堆芯三进三出说明

V14 当前没有把“三进三出连接段”整体放入共用层。共用层只保留 `CoreInletSegment` 这一条泵出口到堆芯入口的短连接段。V14 的三条热端出口支路由 `V14HeatPipeRadiatorConfig.hot_branch_*` 控制，默认尺寸为：

- `HotOutletBranch_1/2/3`: `length=2.19632 m`, `area=pi*0.0138^2 m2`, `Dh=0.0276 m`, `n_nodes=8`
- `CoreInletSegment`: `length=0.13 m`, `area=3.8e-4 m2`, `Dh=0.014 m`, `n_nodes=1`

这与 V15 的三条冷回流支路不是同一组对象，尺寸也不完全一致。

## 水力-only smoke

结果文件：`v14_flow_path_smoke_result.json`

生成命令：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases.V14_run_cases.test_v14_flow_path_smoke
```

测试范围：

- 构建 V14 共用堆芯/主泵层和本地热管辐射器。
- 将辐射器、热管和翅片发射率设为 `0.0`。
- 不调用 `SystemManager.initialize_system()`。
- 不执行热耦合器。
- 不推进任何固体导热模型。
- 只运行 `HydraulicNetwork.initialize_hydraulics()` 和一次 `step_hydraulic()`。

当前 smoke 参数：

- `hydraulic_init_dt_s = 0.01`
- `hydraulic_tol_kg_s = 1.0e-4`
- `hydraulic_max_iter = 1000`
- `hydraulic_step_dt_s = 1.0e-4`

该测试只验证闭式流网连通性和有限流量/压力，不代表热性能或稳态结果。