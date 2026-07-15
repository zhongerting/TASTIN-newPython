# AGENTS.md - TASTIN-Python 总体概括与 Codex 接管导航

> 面向新的 Codex 对话和开发者。  
> 目标：先通过本文建立全局认识，再按任务定向阅读子模块说明；不要默认重新遍历整个代码仓库。  
> 本文基于各模块说明整理，基线日期：2026-06-12。源码与文档冲突时，以当前源码为准，并在任务中同步修正文档。
>
> Codex 在本仓库内开始任务时，应将本文作为项目级首读指令，并按本文导航定向读取子模块说明。

## 0. Python 运行环境

在本仓库中运行 Python 脚本时，必须使用以下 Conda 环境中的解释器：

```text
E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe
```

PowerShell 示例：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" path\to\script.py
```

## 1. 新对话的最小阅读协议

处理本仓库任务时，按以下顺序建立上下文：

1. 先读本文，只判断任务涉及哪些模块。
2. 阅读“模块导航”中对应的首读文档；不要一开始扫描全部源码。
3. 依据子模块手册中的任务索引，只打开目标源码和直接依赖。
4. 遇到物理公式、单位、公开接口、构造参数、状态保存格式、执行顺序、网格映射或跨模块调用变化时，必须局部核验源码。
5. 完成代码修改后，同步更新受影响的模块手册；跨模块关系变化时还要更新本文。
6. 按 `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md` 的影响地图先跑窄测试，再决定是否进入 smoke、审计或长时运行。

本文是一级入口，不替代模块手册，更不替代源码。历史分析、运行产物、旧测试和优化重构目录只能作为参考，不能自动视为当前事实基准。

## 2. 项目定位

TASTIN-Python 是空间核电源热离子反应堆的多物理场瞬态仿真程序。Python 负责系统装配、瞬态调度、水力、导热、中子学、热管散热和测试运行；部分热离子能量转换（TEC）计算通过 pybind11 调用 C++ 后端。

当前没有统一生产入口：根目录 `main.py` 仍为空。实际运行通常从两类用例层进入：

| 用例层 | 用途 | 当前首读文档 |
| --- | --- | --- |
| `testModule/` | 堆芯全系统组装、v7/v8/v9 CaseA、断点续算、审计、结果提取、分层测试 | [`testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`](./testModule/AI_AGENT_TESTMODULE_ANALYSIS.md) |
| `CoolantLoop/` | 集流环冷却回路生产模型、包装器、诊断、性能分析 | [`CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md`](./CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md) |

## 3. 架构总览

```text
用例与运行层
  testModule/                    堆芯全系统、测试、审计、续算
  CoolantLoop/                   集流环冷却回路模型和运行包装器
        |
        v
宏观组件层 Components/
  ReactorCore / TFEUnit / Pipe / AnnularPipe
  RingHP / HPwithFin / TECCircuitManager
  职责：组装底层实体，向 SystemManager 暴露 solids、couplers 和生命周期钩子
        |
        v
求解器层 Solvers/
  SystemManager                  全局瞬态调度、Picard、回滚、断点
  Hydrodynamics/                 流体控制体、连接、水力网络
  HeatConduction/                1D/2D 固体有限体积导热和边界
  Couplers.py                    流固、固固、间隙、TEC 耦合
  Neutronics/PointReactor.py     点堆动力学和衰变热
        |
        v
基础物理与工具层
  Materials/                     流体、固体、气隙、热管和吸液芯物性
  Correlations/                  压降、空泡份额和换热经验关联式
  MathSolvers/                   ODE 与水力稀疏雅可比辅助
  profiler.py                    轻量性能剖析

TEC C++ 支路
  ReactorCore 或 TECCircuitManager
    -> ThermoCalc/ThermoCalcWrapper.py
    -> bindings.cpp
    -> circuitTECs
    -> singleThermionicEnergyConversion
    -> thermionicEmission
```

更长的依赖图和历史架构说明可参考 [`ARCHITECTURE.md`](./ARCHITECTURE.md)。定位辅助目录时参考 [`AI_AGENT_MISCELLANEOUS_ANALYSIS.md`](./AI_AGENT_MISCELLANEOUS_ANALYSIS.md)。

## 4. 模块导航

### 4.1 一级首读文档

| 模块 | 职责 | 进入该模块时先读 |
| --- | --- | --- |
| `Components/` | 宏观组件装配；堆芯、TFE、TEC、电热映射、热管、集流环、管道、外热源 | [`Components/COMPONENTS_DETAILED_INTRO.md`](./Components/COMPONENTS_DETAILED_INTRO.md) |
| `Solvers/` | 全局调度、水力网络、导热、耦合器、点堆、状态持久化 | [`Solvers/AI_AGENT_SOLVERS_ANALYSIS.md`](./Solvers/AI_AGENT_SOLVERS_ANALYSIS.md) |
| `Materials/` | 材料基类、Na/NaK/K 流体、堆芯固体、气隙、热管和吸液芯物性 | [`Materials/AI_AGENT_MATERIALS_ANALYSIS.md`](./Materials/AI_AGENT_MATERIALS_ANALYSIS.md) |
| `Correlations/` | 单相/两相压降、棒束模型、空泡份额和换热关联式 | [`Correlations/AI_AGENT_CORRELATIONS_ANALYSIS.md`](./Correlations/AI_AGENT_CORRELATIONS_ANALYSIS.md) |
| `ThermoCalc/` | TEC Python 包装、pybind11 绑定和 C++ 热离子电路后端 | [`ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`](./ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md) |
| `testModule/` | 当前堆芯组装主链、分层测试、restart、审计、长时运行 | [`testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`](./testModule/AI_AGENT_TESTMODULE_ANALYSIS.md) |
| `CoolantLoop/` | 集流环冷却回路模型、包装器、断点和诊断 | [`CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md`](./CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md) |

### 4.2 二级专题文档

进入一级文档后，再按任务继续读取：

| 任务范围 | 二级文档 |
| --- | --- |
| `Fuel`、电极、`TECPair`、`HeatPipe2D`、翅片 | [`Components/BASICCOMPONENTS_DETAILED_INTRO.md`](./Components/BASICCOMPONENTS_DETAILED_INTRO.md) |
| 基础组件公式和底层调用链 | [`Components/BASICCOMPONENTS_ANALYSIS.md`](./Components/BASICCOMPONENTS_ANALYSIS.md) |
| 轨道外热源、热流查表、外热边界 | [`Components/EXTERNALHEATSOURCES_DETAILED_INTRO.md`](./Components/EXTERNALHEATSOURCES_DETAILED_INTRO.md) |
| 水力网络、泵、稳压器、动态阻力 | [`Solvers/HYDRODYNAMICS_DETAILED_INTRO.md`](./Solvers/HYDRODYNAMICS_DETAILED_INTRO.md) |
| 固体网格、边界、辐射和热流单位 | [`Solvers/HEATCONDUCTION_DETAILED_INTRO.md`](./Solvers/HEATCONDUCTION_DETAILED_INTRO.md) |
| 流固、固固、间隙和 TEC 耦合 | [`Solvers/COUPLERS_DETAILED_INTRO.md`](./Solvers/COUPLERS_DETAILED_INTRO.md) |
| 点堆、衰变热、提交和恢复 | [`Solvers/NEUTRONICS_DETAILED_INTRO.md`](./Solvers/NEUTRONICS_DETAILED_INTRO.md) |
| Picard、源项、回滚、快照和自适应步长 | [`Solvers/SYSTEMMANAGER_DETAILED_INTRO.md`](./Solvers/SYSTEMMANAGER_DETAILED_INTRO.md) |
| Solvers 跨模块深度分析 | [`Solvers/SOLVERS_ANALYSIS.md`](./Solvers/SOLVERS_ANALYSIS.md) |

`ThermoCalc/THERMOCALC_ANALYSIS.md` 和 `ThermoCalc/NONUNIFORM_GRID_GUIDE.md` 是历史或专题材料，不是当前 ThermoCalc 源码事实基准。

涉及热离子查表数据库、runtime 表、dense runtime v2、扫描/优化/导出/加载流程时，先读 [`ThermoCalc/EMISSION_SCAN_GUIDE.md`](./ThermoCalc/EMISSION_SCAN_GUIDE.md) 的 `End-to-End Lookup Workflow`，不要在本文展开流程细节。

## 5. 两条主要系统路径

### 5.1 堆芯全系统路径

当前主要组装入口是 `testModule/test_core_assemble_v7_caseA.py`：

```text
test_core_assemble_v7_caseA.py
  -> build_v7_case_a_system()
  -> SystemManager
  -> ReactorCore
       -> 多个代表性 TFEUnit
       -> 全局慢化剂和结构区域
       -> PointReactor
       -> ThermoCalcModel
  -> HydraulicNetwork
  -> 瞬态运行、restart、审计和结果提取脚本
```

排查 v7 CaseA 时先读 `testModule` 手册，不要从旧版 `test_core_assemble_v1.py` 到 `v6.py` 顺序通读。旧版本用于演进参考。

V8 CaseA 位于 `testModule/test_core_assemble_v8_caseA.py`，用于把外圈拆为 `Ring3_TEC` 和 `Ring3_Open` 两个代表元件。它复用 V7 CaseA 几何和边界，但拥有独立水力拓扑与 restart；旧 V7 restart 必须通过 `testModule/migrate_v7_caseA_restart_to_v8.py` 迁移，不能直接加载。2026-06-10 起，V8 默认冷却剂已改为 `Materials.Fluids.SodiumPotassium78.SodiumPotassium78`，与 `CoolantLoop` 集流环算例保持一致；从旧 Sodium V8 restart 续算前应先运行 `testModule/migrate_v8_caseA_sodium_restart_to_nak.py`，该脚本会按当前 NaK78 材料由已保存的 `T/P` 重算 `h_vec` 和流体物性。

V9 CaseA 位于 `testModule/test_core_assemble_v9_caseA.py`，是基于 V8 堆芯和 TEC 配置的开式外管路骨架。它重建水力网络为“固定流量入口 -> 辐射器出口支路/总管 -> 冷回流支路 -> 堆芯 -> 三条热出口支路 -> 固定压力出口”，暂不包含集流环、热管、泵或局部阻力系数；运行入口为 `testModule/run_v9_caseA_open_loop.py`。V8 restart 不能直接加载到 V9，除非后续新增专门迁移器。

后续把 V9 与 `CoolantLoop` 的 6 段集流环连接时，不能把三条 V9 热出口支路直接接到单个显式集流环的 `I1/I2/I3`。当前集流环模型只显式建一套物理集流环，并用 `MacroFlowJunction(multiplier=2)` 代表第二套对称集流环；集成链条中热出口支路到 `I1/I2/I3` 需要宏观到单环的 `multiplier=2` 分流，`O1/O2/O3` 回到 V9 冷侧管路也需要匹配的单环到宏观汇流。

V10 CaseA 位于 `testModule/test_core_assemble_v10_caseA.py`，运行入口为 `testModule/run_v10_caseA_open_loop.py`。V10 已把 V9 堆芯侧与 `CoolantLoop/model_collector_ring_6segment_v9_interface.py` 的 6 段集流环连接为开式大链条，沿用固定入口总流量和出口定压边界，不添加泵或稳压器。V10 初始化应通过 runner 从 V9 带 TEC restart 和集流环接口 restart 注入状态；V8/V9 restart 不能直接作为 V10 `--restart-in` 使用。注入后必须保留集流环内部环段和 Manifold 的 restart 流量分布，只重置外部边界、堆芯通道、热出口和冷回流设计流量。

V11 CaseA 位于 `testModule/test_core_assemble_v11_caseA.py`，运行入口为 `testModule/run_v11_caseA_closed_loop.py`。V11 复用当前 V10 调参后的几何、NaK78 工质、辐射率、半导线电阻、RK45 固体求解和局部隐式流固换热设置；删除 V10 的开式入口/出口边界，在 `CoreInletConnector` 使用 `is_pressure_reference=True` 作为唯一被动压力参考点，并在 `RadiatorOuterHeader_52` 与 `V11_PumpOutletDistributor_51` 之间加入两台串联相同 `PumpJunction`。默认总泵压头为 `6466.56 Pa`，每台泵 `3233.28 Pa`。V11 应先从当前 V10 restart 注入状态，runner 会在注入时重建泵后冷回流段压力场；不要在注入前做冷态闭式水力初始化，也不要把实际流动节点设为 `is_pressure_boundary=True`，否则会冻结该节点焓/温度。

2026-06-25 起，V11/V13 默认仍只启用主串联 TEC 电路，`Ring3_Open` 代表的 3 根 TEC 默认断开；需要把这三根预留 TEC 接入独立并联电路时，在 runner 中使用 `--enable-reserved-parallel-tec` 以及 `--reserved-parallel-mode fixed_u|fixed_i|load_curve`。详细接口、诊断字段和验证命令见 `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md` 与 `ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`。

2026-06-25 起，V13 可选启用第一版准稳态遮热罩模型，默认关闭。该模型只通过 `RadiatorPipeWithFin` 的等效辐射背景温度改变管翅式辐射器对外散热，不引入全局 ODE，也不影响 V10/V11。入口参数和验证记录见 `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`。

2026-06-11 V9 已完成无 TEC 预热和带 TEC 续算验证：`testModule/v9_caseA_open_loop_tec_3000s/` 最终绝对时间约 `5010 s`，入口温度约 `743.000 K`，出口温度约 `838.944 K`，进出口温差约 `95.944 K`，冷却剂焓升约 `108745.775 W`，端电功率约 `4870.107 W`。该目录是本地运行产物，不应默认提交；需要复现实验时优先使用 V9 runner 和 V9 restart。

### 5.2 集流环冷却回路路径

`CoolantLoop/` 当前有两个生产模型：

| 模型 | 特点 | 优先用途 |
| --- | --- | --- |
| `model_collector_ring_6segment.py` | 6 个首尾拼接的 `1/6 RingHP` 扇区；内置 profiler | 分扇区定位、性能分析 |
| `model_collector_ring_full_ringhp.py` | 单一 360 度 `RingHP`，共 24 个环节点 | 完整环模型 |

两者都通过 `SystemManager.save_global_state()` / `load_global_state()` 使用 `.npz` 断点，并依赖 `Components/RingHP.py`、`HPwithFin.py`、水力网络、导热边界和材料物性。

## 6. 核心跨模块契约

修改代码前优先确认以下契约是否受影响：

| 契约 | 当前约定 |
| --- | --- |
| 宏观组件接入 | `BaseComponent` 通过 `get_solids()`、`get_couplers()`、`pre_step()`、`post_step()` 接入 `SystemManager`；支持迭代回滚和断点的组件还要维护状态接口。 |
| 全局推进 | `SystemManager.step()` 负责组件钩子、耦合器、流体、固体、可选点堆、Picard 回滚和提交。调整顺序时先读 Solvers 专题文档。 |
| 热源单位 | 离散节点和边界累计热功率通常使用 `W`；外部热流和 TEC 面热流使用 `W/m2`，由边界或耦合层在正确位置乘面积。 |
| TEC 焦耳热 | 电势或电场通过 `Components/tec_electric.py` 映射为每节点 `W`，非均匀网格必须按节点中心坐标处理。 |
| 流体物性 | 水力层依赖 `enthalpy()`、`temperature_from_enthalpy()`、`density()`、`viscosity()`、`heat_capacity()` 和密度导数；关联式还依赖液相/气相兼容接口。 |
| 状态恢复 | 当前统一使用 `save_global_state("*.npz")` 和 `load_global_state("*.npz")`。加载成功不代表重建系统的几何、倍率和求解器配置一定匹配。 |
| 网络拓扑 | `HydraulicNetwork` 初始化后视为拓扑固定；改变节点、连接或定压边界集合时应重建网络。 |
| 宏观倍率 | `ReactorCore` 和 `RingHP` 使用代表单元倍率恢复物理总量。修改倍率时必须核对功率、热源、反馈、流量和统计口径。 |

## 7. 当前高风险区域

这些问题在进一步修改或运行前需要主动核验：

1. `ThermoCalc` 已于 2026-06-01 闭合逐节点侧面积、`phiE/phiC/Vd` 结果读取和构建后 `Tcs` 热更新，并重新构建 `te_solver.cp312-win_amd64.pyd`。串联 `fixed_I` 已支持给定非负总电流；不可发电时回退零电流开路状态并标记未收敛。V11/V13 的并联固定总电流仍只用于可选预留并联电路，应通过 `ReactorCore.setup_reserved_parallel_tec_circuit(mode_str="fixed_i", ...)` 映射到 `parallel_fixed_i`。
2. 默认 `python` 仍为 3.9.13，而主要扩展为 `te_solver.cp312-win_amd64.pyd`。如果使用 Python 运行，请使用 `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`；运行 ThermoCalc 必须使用 ABI 匹配的 Python 3.12 环境；本轮验证解释器记录在 `ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`。
3. `Correlations.h_single_crossflow_pipe()` 当前调用参数不足，直接使用会触发运行时 `TypeError`。
4. `Materials/Fluids/NaK78.py` 仅部分实现；主用液态 NaK78 类是 `SodiumPotassium78`。不要仅按文件名推断模型使用的冷却剂。
5. `Components/basicComponents/Electord.py` 的文件名拼写已被现有导入路径依赖，不要直接重命名。
6. `ExternalHeatSources` 返回 `W/m2`，`ExternalHeatFluxBC` 负责乘面积。调用方不得重复乘面积。
   V14 的 `HPwithFin/RingHP` 外热吸收还要乘 `0.992 * epsilon_surface`，只作用于外侧受照面积；内侧只参与角系数等效辐射，不吸收外热，`hp_multipliers` 只能在汇总层使用。详见 `Components/EXTERNALHEATSOURCES_DETAILED_INTRO.md`。
7. `HPwithFin` 使用内部降维准稳态翅片模型；`basicComponents/FinConduction.py` 是独立二维翅片求解器，两者不要混用。
8. `CoolantLoop` 的历史测试、旧 bug report、CSV、PNG、日志和 profiler 报告不代表当前生产模型默认工况。

风险的证据、边界条件和处理原则应回到对应一级文档继续阅读。

## 8. 按任务选择阅读路径

| 任务 | 先读 | 再按需读取 |
| --- | --- | --- |
| 修改 v7 CaseA 系统组装、续算或审计 | `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md` | `test_core_assemble_v7_caseA.py`、目标 `run_*` / `audit_*` 脚本、`Components` 和 `SystemManager` 文档 |
| 修改全局 Picard、生命周期、回滚、restart | `Solvers/AI_AGENT_SOLVERS_ANALYSIS.md` | `SYSTEMMANAGER_DETAILED_INTRO.md`、`SystemManager.py`、`testModule` 影响地图 |
| 修改压力、流量、泵、稳压器或水力拓扑 | `Solvers/AI_AGENT_SOLVERS_ANALYSIS.md` | `HYDRODYNAMICS_DETAILED_INTRO.md`、水力源码、对应水力测试 |
| 修改导热、边界、辐射或热流单位 | `Solvers/AI_AGENT_SOLVERS_ANALYSIS.md` | `HEATCONDUCTION_DETAILED_INTRO.md`、目标组件文档、对应导热和能量测试 |
| 修改 TFE、ReactorCore、倍率或反馈 | `Components/COMPONENTS_DETAILED_INTRO.md` | `TFEUnit.py`、`ReactorCore.py`、基础组件文档、v7 CaseA 审计 |
| 修改 TEC 电势、焦耳热、面热流或 C++ 后端 | `ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md` 和 `Components/COMPONENTS_DETAILED_INTRO.md` | `tec_electric.py`、`TECCircuitManager.py`、`ReactorCore.py`、绑定和 C++ 源码、TEC 定向测试 |
| 修改热管、翅片、集流环或轨道外热 | `Components/COMPONENTS_DETAILED_INTRO.md` | 基础组件和外热专题文档、`CoolantLoop` 手册、热管能量测试 |
| 修改材料公式、物性边界或向量化 | `Materials/AI_AGENT_MATERIALS_ANALYSIS.md` | 目标材料源码、调用方、标量和数组验证 |
| 修改摩阻、压降或换热关联式 | `Correlations/AI_AGENT_CORRELATIONS_ANALYSIS.md` | `Correlations.py` 目标函数、材料契约、调用方测试 |
| 修改集流环冷却回路生产模型 | `CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md` | 目标 `model_*.py`、相关上级模块文档、最小诊断脚本 |
| 分析性能 | `AI_AGENT_MISCELLANEOUS_ANALYSIS.md` | `profiler.py`、`CoolantLoop` profiler 路径、各优化重构目录中的 benchmark |

## 9. 测试与验证原则

所有验证入口优先从 [`testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`](./testModule/AI_AGENT_TESTMODULE_ANALYSIS.md) 的“修改影响地图”选择。

基本顺序：

1. 先运行与改动直接相关的最小测试。
2. 涉及跨模块耦合时，再运行对应集成测试或审计脚本。
3. 涉及 v7/v8/v9 CaseA 时，优先使用低成本 smoke、短时续算、拓扑测试和审计脚本。
4. 只有任务确实需要长时行为时，才启动过夜或 `100000 s` 级续算。
5. 使用 restart 前，核对重建系统的几何、倍率、材料、TEC 配置和状态同步流程。

`CoolantLoop/` 的最小热管/header 验证、壁面辐射隔离诊断、完整集流环回归和 profiler 路径在其模块手册中单独维护。

## 10. 辅助目录和非首选入口

| 路径 | 定位 |
| --- | --- |
| [`AI_AGENT_MISCELLANEOUS_ANALYSIS.md`](./AI_AGENT_MISCELLANEOUS_ANALYSIS.md) | `MathSolvers/`、`profiler.py`、空 `main.py`、优化重构目录、遗留 `inputs/` 的速查手册 |
| [`ARCHITECTURE.md`](./ARCHITECTURE.md) | 较完整的架构和依赖关系参考；适合需要更长全局图时阅读 |
| `SystemManager优化重构/`、`HeatPipe优化重构/`、`HeatConduction优化重构/`、`HydraulicNetwork优化重构/` | 性能验证和重构过程脚本，不是核心运行入口 |
| `inputs/` | Fortran 时代输入文件；Python 版是否直接使用要按调用点核验 |
| `*.npz`、`*.json`、`*.csv`、`*.out`、`*.err`、`*.log`、`*.png` | restart、摘要、审计和运行产物，不是源码事实基准 |
| `testModule/v7_caseA_*`、`testModule/v8_caseA_*`、`testModule/v9_caseA_*` | 长时算例输出目录；除非用户明确要求，不应加入提交 |

## 11. 文档维护规则

后续任务中，若出现以下变化，应更新本文：

- 新增或删除一级模块、生产模型、主要运行入口或首读文档。
- 改变 `SystemManager`、`BaseComponent`、ThermoCalc 或 restart 的跨模块契约。
- 改变堆芯主链、CoolantLoop 生产主链或推荐测试入口。
- 修复或新增需要所有新会话知晓的高风险约束。

局部公式、类字段和详细 API 变化应优先更新对应模块手册，本文只保留接管时必须知道的全局信息和阅读分流。

## 12. 2026-06-02 TEC 焦耳热映射补充

TEC 生产焦耳热当前以 C++ `VcalcFVM()` 输出的逐轴向节点功率 `joulePowerE/C [W]` 为权威值。Python 层通过 `Components/tec_electric.py::distribute_axial_power_by_volume()` 将轴向功率按列内二维控制体体积比例分配到电极导热网格。

`UE / UC / rhoE / rhoC` 和节点中心梯度函数继续保留为诊断数据，但不得重新作为生产焦耳热源。2026-06-02 单 TFE TEC `1 s` 基线中，二维映射与 C++ 节点功率差为 `0 W`，TEC 转换闭合差约 `0.0250 W`，最终全局残差约 `0.0891 W`。

v7/V8 CaseA 多 TEC 串联路径仍会重复报告 `Failed to converge after 100000 iterations.`。2026-06-09 V8 `LSODA` smoke 中该信息出现 `4440` 次，全部位于进入长算主循环前的 `core.thermo_calc.calculate(verbose=False)` 阶段；后续 `system.step(0.01)` 未再输出。该问题属于 ThermoCalc C++ 电路收敛专项，不是固体导热 ODE 收敛失败，也不能视为焦耳热映射已经解决的范围。

2026-06-10 起，V8 CaseA 公共加载运行路径 `testModule/run_v8_caseA_common.py` 会向 ThermoCalc 施加导线电阻 `[0.00155199999999970, 0.00102400000000000, 0.000336000000000000, 0.000608000000000000] ohm`，并写入 `latest_state.json` 的 `wire_resistance_ohm`。该电阻必须在 restart 加载、TEC fixed voltage 设置和 `core.post_step(...)` 电极温度同步之后再重建 ThermoCalc；提前计算会显著拖慢首次 TEC 求解。V8 默认冷却剂同日改为 `SodiumPotassium78`；`testModule/v8_caseA_nak_wire_2000s/` 已从绝对时间 `41184 s` 续算到 `43184 s`，末段端电功率约 `4868.580 W`，冷却剂焓升约 `108560.205 W`，进出口温差约 `95.781 K`，thermal-model residual 约 `-0.132 W`，日志未再出现 `Failed to converge after 100000 iterations.`。

V9 复用该导线电阻和 NaK78 冷却剂约定。冷态直接启动 TEC 可能使首次 ThermoCalc 求解返回非有限焦耳热；冷态拓扑 smoke 应使用 `--disable-tec-coupled`，带 TEC 长算应从 V9 兼容的温热 restart 或专门迁移器输出继续。

## 13. 2026-06-02 全局慢化剂映射顺序补充

`ReactorCore.pre_step()` 必须先执行 `TFEUnit.pre_step()` 更新内部等效 moderator 外边界温度，再刷新该 moderator 的边界热流缓存，最后将外流按 `tfe_multipliers` 聚合到全局慢化剂环。不得直接复用旧时间层 `BoundaryRegion.current_flux`；该错误在 restart 后会放大为首步非物理源项脉冲。

## 14. 2026-06-28 ThermoCalc 低温零发射保护

`ThermoCalcModel.calculate()` 在进入 C++ 电路迭代前会进行低温零发射快速判定：当所有发射极温度低于默认 `THERMOCALC_ZERO_EMISSION_TE_MAX_K=1000 K` 时，直接返回开路零电流、零 TEC 热源结果，并在 `get_global_results()` 中标记 `zero_emission_skipped=True`。这不是禁用 TEC，而是避免启动低温下无发电能力的 fixed-voltage 电路迭代卡死；详细说明见 [`ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`](./ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md)。

补充：2026-06-28 起，除了 Python 层低温零发射预判，`ThermoCalc` C++ 电路层也加入了退化迭代快速跳出；因此启动/过渡工况可以默认保持 TEC 开启，无法形成有效电解的状态会返回有限的零电流或未收敛诊断，而不是长期卡在 `ThermoCalcModel.calculate()`。详细见 ThermoCalc 手册。


## 15. 2026-07-14 当前生产默认值

通用固体导热默认 `implicit_euler`，流固时间耦合默认 `local_implicit`，`SystemManager.step()` 默认外部耦合迭代次数为 `1`。标准生产运行应通过 `compute_adaptive_dt()` 选步长；显式固定步长的 benchmark、能量审计和历史复现脚本继续按其参数运行。点堆默认 `Radau`。

TEC 默认串联 `series`、定电压 `fixed_u`，发射计算优先自动发现本地 dense runtime v2 查表库，查询 miss 或本地无库时使用解析法。外热源默认关闭；显式启用时，翅片默认按单侧几何投影面积将外热加入各轴向段的准稳态翅片方程。遮热罩、TFE/堆芯间隙保持原配置，泵默认固定压头。
