# TASTIN-Python

TASTIN-Python 是空间核电源热离子反应堆的多物理场瞬态仿真程序。Python 负责系统装配、瞬态调度、水力、导热、中子学、热管散热和测试运行；部分热离子能量转换（TEC）计算通过 pybind11 调用 C++ 后端。

当前仓库没有统一生产入口：根目录 `main.py` 仍为空。实际运行通常从 `testModule/` 或 `CoolantLoop/` 的用例脚本进入。

## Python 运行环境

在本仓库中运行 Python 脚本时，请使用以下 Conda 环境中的解释器：

```text
E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe
```

PowerShell 示例：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" path\to\script.py
```

## 主要入口

| 用例层 | 用途 | 首读文档 |
| --- | --- | --- |
| `testModule/` | 堆芯全系统组装、v7/v8/v9 CaseA、断点续算、审计、结果提取、分层测试 | [`testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`](./testModule/AI_AGENT_TESTMODULE_ANALYSIS.md) |
| `CoolantLoop/` | 集流环冷却回路生产模型、包装器、诊断、性能分析 | [`CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md`](./CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md) |

## 架构概览

```text
用例与运行层
  testModule/                    堆芯全系统、测试、审计、续算
  CoolantLoop/                   集流环冷却回路模型和运行包装器
        |
        v
宏观组件层 Components/
  ReactorCore / TFEUnit / Pipe / AnnularPipe
  RingHP / HPwithFin / TECCircuitManager
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

更完整的依赖关系可参考 [`ARCHITECTURE.md`](./ARCHITECTURE.md)。

## 模块导航

| 模块 | 职责 | 说明文档 |
| --- | --- | --- |
| `Components/` | 宏观组件装配；堆芯、TFE、TEC、电热映射、热管、集流环、管道、外热源 | [`Components/COMPONENTS_DETAILED_INTRO.md`](./Components/COMPONENTS_DETAILED_INTRO.md) |
| `Solvers/` | 全局调度、水力网络、导热、耦合器、点堆、状态持久化 | [`Solvers/AI_AGENT_SOLVERS_ANALYSIS.md`](./Solvers/AI_AGENT_SOLVERS_ANALYSIS.md) |
| `Materials/` | 材料基类、Na/NaK/K 流体、堆芯固体、气隙、热管和吸液芯物性 | [`Materials/AI_AGENT_MATERIALS_ANALYSIS.md`](./Materials/AI_AGENT_MATERIALS_ANALYSIS.md) |
| `Correlations/` | 单相/两相压降、棒束模型、空泡份额和换热关联式 | [`Correlations/AI_AGENT_CORRELATIONS_ANALYSIS.md`](./Correlations/AI_AGENT_CORRELATIONS_ANALYSIS.md) |
| `ThermoCalc/` | TEC Python 包装、pybind11 绑定和 C++ 热离子电路后端 | [`ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`](./ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md) |
| `testModule/` | 堆芯组装主链、分层测试、restart、审计、长时运行 | [`testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`](./testModule/AI_AGENT_TESTMODULE_ANALYSIS.md) |
| `CoolantLoop/` | 集流环冷却回路模型、包装器、断点和诊断 | [`CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md`](./CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md) |

## 主要系统路径

### 堆芯全系统

当前主要组装入口是 [`testModule/test_core_assemble_v7_caseA.py`](./testModule/test_core_assemble_v7_caseA.py)：

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

V8 CaseA 位于 [`testModule/test_core_assemble_v8_caseA.py`](./testModule/test_core_assemble_v8_caseA.py)，用于把第四圈燃料元件拆成 `Ring3_TEC` 和 `Ring3_Open` 两类代表元件：接入电路的 34 根元件和不接入电路的 3 根元件。V8 复用 V7 CaseA 的主要几何和边界，但拥有独立水力拓扑与 restart；旧 V7 restart 需要通过迁移脚本转换后使用。当前 V8 默认冷却剂为 `SodiumPotassium78`，公共运行入口会施加标准导线电阻并使用 fixed-voltage TEC。

V9 CaseA 位于 [`testModule/test_core_assemble_v9_caseA.py`](./testModule/test_core_assemble_v9_caseA.py)，运行入口是 [`testModule/run_v9_caseA_open_loop.py`](./testModule/run_v9_caseA_open_loop.py)。V9 基于 V8 堆芯和 TEC 配置，重建外部开式水力骨架为：

```text
固定流量入口
  -> 辐射器出口支路/总管
  -> 冷回流支路
  -> 堆芯
  -> 三条热出口支路
  -> 固定压力出口
```

V9 暂不包含集流环、热管、泵和局部阻力系数；`CoreInletConnector` / `CoreOutletConnector` 只是数值连接节点，不代表真实进出口箱。V8 restart 不能直接加载到 V9，后续若需要迁移应使用专门迁移器。

V9 与集流环模型连接时需要保留 `CoolantLoop` 当前“一套显式集流环 + `MacroFlowJunction(multiplier=2)` 代表第二套对称集流环”的约定。三条 V9 热出口支路不能直接硬接到单个显式集流环的 `I1/I2/I3`，应先通过 `multiplier=2` 宏观到单环分流；`O1/O2/O3` 返回 V9 冷侧管路时也需要匹配汇流。

### 集流环冷却回路

`CoolantLoop/` 当前有两个生产模型：

| 模型 | 特点 | 优先用途 |
| --- | --- | --- |
| `model_collector_ring_6segment.py` | 6 个首尾拼接的 `1/6 RingHP` 扇区；内置 profiler | 分扇区定位、性能分析 |
| `model_collector_ring_full_ringhp.py` | 单一 360 度 `RingHP`，共 24 个环节点 | 完整环模型 |

两者都通过 `SystemManager.save_global_state()` / `load_global_state()` 使用 `.npz` 断点。

## 常用运行入口

V9 拓扑和轻量验证：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\test_core_assemble_v9_caseA.py testModule\run_v9_caseA_open_loop.py testModule\test_v9_caseA_topology.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_v9_caseA_topology
```

V9 开式骨架短算：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v9_caseA_open_loop.py --duration 0.1 --record-interval 0.1 --restart-interval 0.1 --max-dt 0.01 --disable-tec-coupled
```

V9 带 TEC 长算应从 V9 兼容的温热 restart 启动，不要直接使用 V8 restart。

## 运行环境

- 如果使用 Python 运行，请使用 `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`。
- Python：普通热工/水力脚本通常可用当前项目环境运行。
- ThermoCalc：主用扩展为 `te_solver.cp312-win_amd64.pyd`，运行 TEC 相关功能时需要 ABI 匹配的 Python 3.12 环境。
- 常用依赖：NumPy、SciPy；部分后处理或绘图脚本可能需要 Matplotlib。

注意：默认 `python` 未必是 Python 3.12。运行 ThermoCalc 前先确认解释器和 `.pyd` ABI 匹配。

## 测试与验证

验证入口优先参考 [`testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`](./testModule/AI_AGENT_TESTMODULE_ANALYSIS.md) 的影响地图。基本顺序：

1. 先运行与改动直接相关的最小测试。
2. 涉及跨模块耦合时，再运行对应集成测试或审计脚本。
3. 涉及 v7/v8/v9 CaseA 时，优先使用低成本 smoke、短时续算和能量审计脚本。
4. 只有任务确实需要长时行为时，才启动长时续算。
5. 使用 restart 前，核对重建系统的几何、倍率、材料、TEC 配置和状态同步流程。

## 当前注意事项

- `ThermoCalc` 已支持逐节点侧面积、`phiE/phiC/Vd` 结果读取和构建后 `Tcs` 热更新；`fixed_I` 模式仍不支持。
- TEC 焦耳热生产映射以 C++ `VcalcFVM()` 输出的逐轴向节点功率 `joulePowerE/C [W]` 为权威值，再由 `Components/tec_electric.py` 按二维控制体体积分配。
- `ExternalHeatSources` 返回 `W/m2`，`ExternalHeatFluxBC` 负责乘面积，调用方不要重复乘面积。
- `Components/basicComponents/Electord.py` 的文件名拼写已被现有导入路径依赖，不要直接重命名。
- `Materials/Fluids/NaK78.py` 仅部分实现；主用液态 NaK78 类是 `SodiumPotassium78`。
- V8/V9 CaseA 当前使用 `SodiumPotassium78` 冷却剂；V8 旧 Sodium restart 续算前需要先迁移。
- V8/V9 的 TEC fixed-voltage 路径会施加标准导线电阻 `[0.00155199999999970, 0.00102400000000000, 0.000336000000000000, 0.000608000000000000] ohm`。
- V7/V8 多 TEC 串联路径可能在 ThermoCalc 初始化阶段输出 `Failed to converge after 100000 iterations.`；该信息来自 C++ 电路求解，不是固体导热 ODE 收敛失败。
- `Correlations.h_single_crossflow_pipe()` 当前直接使用可能触发参数不足的 `TypeError`，调用前需要核验。
- 运行产物、restart、审计 CSV/JSON、PNG 和日志不是源码事实基准，不应作为代码提交内容。

## Git 与运行产物

仓库通过 `.gitignore` 排除常见运行产物：

```text
*.npz
*.csv
*.png
*.log
*.err
*.out
testModule/v7_caseA_*/
testModule/v8_caseA_*/
testModule/v9_caseA_*/
testModule/single_tfe_energy_conservation_v7*/
```

提交时应只记录源码、配置和说明文档。长时算例的 restart 点和审计输出保留在本地工作区即可。

## 开发导航

新开发任务建议按以下顺序建立上下文：

1. 先读 [`AGENTS.md`](./AGENTS.md)，判断任务涉及哪些模块。
2. 阅读模块导航中对应的首读文档。
3. 只打开目标源码和直接依赖，不默认遍历整个仓库。
4. 遇到物理公式、单位、公开接口、构造参数、状态保存格式、执行顺序、网格映射或跨模块调用变化时，必须局部核验源码。
5. 修改代码后，同步更新受影响模块的说明文档。

## 许可

本项目仅供研究使用。
