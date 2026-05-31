# testModule Codex 速查手册

> 本目录是 TASTIN 的测试、系统组装、长时运行、断点续算、诊断审计和结果提取入口，不是核心求解器源码目录。
>
> 当前顶层共有 **63 个 Python 脚本**、**23 个运行产物目录**。统计不包含 `__pycache__/`。

## 1. Codex 阅读顺序

1. 先读本手册，确认任务属于哪条链路。
2. 需要理解全项目架构时，再读父目录 [`ARCHITECTURE.md`](../ARCHITECTURE.md)。
3. 仅在修改功能、核对参数或排查异常时，定向打开对应脚本和它直接依赖的源码。
4. 本手册用于快速建立上下文；参数和行为的最终事实来源仍是当前源码。

## 2. 当前主链路

`test_core_assemble_v7_caseA.py` 是当前主要系统组装入口。它构建完整 `SystemManager`、水力网络、`ReactorCore`、TFE、慢化剂区域和 TEC 电热耦合，并提供加热瞬态运行入口。

```text
test_core_assemble_v7_caseA.py
  ├── build_v7_case_a_system()
  ├── run_test_v7_case_a_flow_only()
  ├── run_test_v7_case_a_heated()
  ├── _case_a_flow_diagnostics()
  ├── _case_a_electric_diagnostics()
  └── _case_a_reset_design_flows_after_restart()

run_v7_caseA_*.py
  ├── 短时续算、长时续算和过夜运行
  ├── restart 快照管理
  └── TEC、流量和能量状态记录

audit_v7_caseA_*.py / diagnose_v7_caseA_*.py / extract_v7_caseA_*.py
  ├── 基于 v7_caseA 构建系统
  ├── 加载 restart 状态
  └── 同步耦合状态后输出审计或提取结果
```

`run_v7_caseA_multipliers_short.py` 同时具有两种用途：

- 可执行短时检查：加载 restart，按自定义环乘子续算，写出 TEC 分布和摘要。
- 被其他脚本复用的工具模块：提供 `parse_multipliers()`、`json_default()`、`axial_node_volumes()`、`collect_tec_stats()` 和 `build_loaded_case()`。

## 3. v7_caseA 关键事实

### 3.1 基础参数

| 参数 | 当前值 | 源码位置 |
|---|---:|---|
| 加热运行默认总功率 | `115000.0 W` | `run_test_v7_case_a_heated()` |
| 进口温度 | `743.0 K` | `build_v7_case_a_system()` |
| 热工环乘子 | `(1, 6, 12, 18)`，合计 `37` | `CASE_A_RING_MULTIPLIERS` |
| TEC 环乘子 | `(1, 6, 12, 15)`，合计 `34` | `CASE_A_TEC_RING_MULTIPLIERS` |
| 设计总流量 | `1.30 kg/s` | `CASE_A_DESIGN_TOTAL_FLOW_KG_S` |
| 单热工通道设计流量 | `1.30 / 37 kg/s` | `CASE_A_DESIGN_CHANNEL_FLOW_KG_S` |
| TEC 固定电压 | `27.2 V` | `build_v7_case_a_system()` 内 `core.setup_tec_circuit()` |
| TEC 初始电流猜测 | `284.0 A` | `build_v7_case_a_system()` 内 `core.setup_tec_circuit()` |
| 默认管段节点数 | `8` | `build_v7_case_a_system()` |
| FastSteady 固体热容缩放 | `0.01`，作用域 `global_outer` | `test_core_assemble_v7_caseA_faststeady.py` |

注意：

- `test_core_assemble_v7.py` 是通用 v7 入口，使用 `SodiumPotassium78`；当前主用 `test_core_assemble_v7_caseA.py` 是独立 CaseA 构建入口，导入的是 `Sodium`。修改冷却剂时不要只看文件名推断。
- `build_v7_case_a_system()` 接收热工环乘子和 TEC 环乘子；`run_test_v7_case_a_heated()` 当前只把 `ring_multipliers` 透传到构建函数。需要自定义两组乘子时优先使用 `run_v7_caseA_multipliers_short.py` 及其复用链路。
- `run_v7_caseA_continue_100000s.py` 额外定义导线电阻 `[0.001552, 0.001024, 0.000336, 0.000608] Ω`。这是该续算脚本的专用处理，不是所有 restart 的通用步骤。

### 3.2 restart 与同步流程

`SystemManager` 的当前快照 API 是：

```python
system.save_global_state("restart.npz")
system.load_global_state("restart.npz")
```

`.npz` 中包含系统时间、流体求解器状态、固体状态、宏组件状态和点堆状态。`load_global_state()` 加载后还会同步固体时间、刷新边界缓存、准备流体耦合源项并运行耦合器。

v7_caseA 续算通常按以下路径工作：

```text
按目标参数重新 build_v7_case_a_system()
  -> system.load_global_state(restart_path)
  -> _case_a_reset_design_flows_after_restart(build)
  -> 按具体运行器需要刷新 TEC 或附加参数
  -> 继续推进并保存新的 restart / JSON / CSV
```

不同入口还有额外动作：

| 入口 | 加载后的附加动作 |
|---|---|
| `run_v7_caseA_multipliers_short.py` | 重置设计流量；如 TEC 存在则执行 `thermo_calc.calculate(verbose=False)` |
| `run_v7_caseA_continue_100000s.py` | 多次施加该脚本专用导线电阻，重置设计流量并刷新 TEC |
| `extract_v7_caseA_electric_node_power.py` | 重置设计流量并刷新 TEC，再提取节点级电功率 |
| `audit_v7_caseA_interface_energy.py` | 通过 `sync_for_audit()` 执行 `post_step -> TEC calculate -> fluid source -> boundary cache -> pre_step -> couplers` 同步 |
| `audit_v7_caseA_fluid_energy_network.py` | 复用 `build_loaded_case()` 和 `sync_for_audit()` |

## 4. 脚本索引

### 4.1 堆芯组装演进链

| 文件 | 定位 |
|---|---|
| `test_core_assemble_v1.py` | 早期 TFE、流体通道和 restart 原型 |
| `test_core_assemble_v2.py` | 早期多组件组装迭代 |
| `test_core_assemble_v3.py` | 环形布局迭代 |
| `test_core_assemble_v4.py` | 非均匀通道和宏观水力连接迭代 |
| `test_core_assemble_v5.py` | 提供轴向功率、环功率因子、慢化剂网格等公共构建工具 |
| `test_core_assemble_v6.py` | v6 统一构建和运行逻辑；提供 `get_time_dependent_dt_cap_v6()` |
| `test_core_assemble_v6_caseA.py` | v6 Case A 包装入口 |
| `test_core_assemble_v6_caseB.py` | v6 Case B：默认 `+5e-5` 反应性阶跃 |
| `test_core_assemble_v6_caseC.py` | v6 Case C：默认 `-5e-5` 反应性阶跃 |
| `test_core_assemble_v6_caseD.py` | v6 Case D：从 v5 长时间状态重启并启用 TEC 的包装入口 |
| `test_core_assemble_v7.py` | 通用 v7：NaK78、非均匀流体通道和 TEC 电路 |
| `test_core_assemble_v7_caseA.py` | 当前主要 CaseA 系统组装和加热瞬态入口 |
| `test_core_assemble_v7_caseA_faststeady.py` | 外围固体热容缩放和 FastSteady 能量审计 |

### 4.2 长时运行与续算

| 文件 | 定位 |
|---|---|
| `run_v7_caseA_heated_to_3600.py` | 加热到绝对时间 `3600 s`；自动寻找已有加热 restart |
| `run_v7_caseA_heated_to_6000.py` | 从 `3600 s` restart 续算到 `6000 s` |
| `run_v7_caseA_heated_to_7800.py` | 从 `6000 s` restart 续算到 `7800 s` |
| `run_v7_caseA_faststeady_continue.py` | FastSteady restart 续算；支持目标时间或增量时长 |
| `run_v7_caseA_electric_dt08_outerCp001_overnight.py` | 电气耦合过夜分段运行；默认 `max_dt=0.8 s`、外围 `Cp x 0.01` |
| `run_v7_caseA_multipliers_short.py` | 可执行短时乘子检查，也是多个脚本复用的工具模块 |
| `run_v7_caseA_newpyd_long200000.py` | restart 驱动的长时分段记录器；脚本名保留历史命名，当前默认增量时长为 `20000 s` |
| `run_v7_caseA_continue_100000s.py` | 带专用导线电阻的 `100000 s` 增量续算入口 |

### 4.3 审计、诊断和结果提取

| 文件 | 定位 |
|---|---|
| `audit_v7_caseA_interface_energy.py` | 流固、固固、TEC 和慢化剂界面能量审计 |
| `audit_v7_caseA_fluid_energy_network.py` | 流体源项、焓输运、控制面闭合、矩阵残差和回归审计 |
| `diagnose_v7_caseA_energy_continue.py` | 复用长时推进逻辑，记录续算过程中的能量诊断 |
| `extract_v7_caseA_electric_node_power.py` | 从 restart 提取 TEC 虚拟元件节点级电势、电场和功率数据 |
| `extract_center_channel_temperature_field_v5.py` | 提取 v5 中心通道温度场 |

### 4.4 基础、组件和物理级测试

| 分组 | 文件 |
|---|---|
| 导热历史用例 | `testHeadConduction1D.py`、`testHeadConduction2D.py`、`test_case_3_1.py`、`test_case_3_2.py`、`test_case_3_3.py` |
| 热管和集流环 | `testHPwithFin.py`、`test_HP_with_external_heat_source.py`、`test_ringHP_buzzin.py`、`test_single_hp_fin_energy_conservation.py`、`test_single_hp_header_energy_audit.py`、`test_ringhp_node_coupling_energy_conservation.py` |
| TEC、TFE 和 ThermoCalc | `test_TEC_with_heat.py`、`test_tecboundary.py`、`test_tec_joule_nonuniform.py`、`test_tfe_thermal_flow.py`、`test_thermionic_wrapper.py` |
| 管道组件 | `test_component_pipe.py`、`test_component_annular_pipe.py` |

### 4.5 水力、耦合和 SystemManager 测试

| 分组 | 文件 |
|---|---|
| 流动和换热 | `test_flow_heat.py`、`test_coupled_heating.py`、`test_pipe_heat_transfer.py`、`test_incompressible_heat_flow.py` |
| 流固和固固耦合 | `test_fluid_solid_couple.py`、`test_simple_solid_fluid_couple.py`、`test_verify_solid_couplers.py` |
| 通道和网络拓扑 | `test_open_channel.py`、`test_open_channel_transient.py`、`test_single_channel_transient.py`、`test_parallel_channels.py`、`test_phase_topology.py` |
| 泵和稳压器 | `test_fixed_pressure_pump_hydrodynamics.py`、`test_pump_junction_hydraulic_network.py`、`test_pressurizer_volume.py`、`test_pressurizer_pumped_closed_loop.py` |
| 系统生命周期 | `test_system_manager.py`、`test_system_manager_lifecycle.py` |
| 点堆动力学 | `test_PointReactor.py` |

## 5. 依赖关系

```text
test_core_assemble_v7_caseA.py
  ├── test_core_assemble_v5.py
  │     ├── build_axial_power_profile()
  │     ├── build_ring_power_factors()
  │     └── build_global_moderator_meshes()
  ├── test_core_assemble_v6.py
  │     └── get_time_dependent_dt_cap_v6()
  ├── Components/ReactorCore.py
  ├── Components/TFEUnit.py
  ├── Solvers/SystemManager.py
  ├── Solvers/Hydrodynamics/*
  └── Materials/*

run_v7_caseA_newpyd_long200000.py
  └── run_v7_caseA_multipliers_short.py

diagnose_v7_caseA_energy_continue.py
  └── run_v7_caseA_newpyd_long200000.py

audit_v7_caseA_fluid_energy_network.py
  ├── audit_v7_caseA_interface_energy.py
  └── run_v7_caseA_multipliers_short.py
```

排查主链路问题时，先从入口脚本沿上述关系定向阅读。不要默认扫描整个父目录。

## 6. 运行产物目录索引

当前 23 个顶层产物目录按用途维护。目录可为空，也可包含多轮运行的历史文件；本手册不逐个枚举 `.json`、`.csv`、`.npz`、`.out`、`.err` 和 `.log`。

### 6.1 smoke / debug 快速验证

| 目录 | 用途 |
|---|---|
| `v7_caseA_electric_dt08_outerCp0001_smoke/` | 电气耦合 `dt=0.8`、外围低热容 smoke |
| `v7_caseA_electric_dt08_outerCp001_debug20/` | 电气耦合短时 debug |
| `v7_caseA_electric_dt08_outerCp001_smoke/` | 电气耦合外围 `Cp x 0.01` smoke |
| `v7_caseA_electric_dt08_physicalCp_smoke/` | 物理热容 smoke |
| `v7_caseA_newpyd_long200000_debug/` | 长时记录器 debug |
| `v7_caseA_newpyd_long200000_ps_smoke/` | PowerShell 后台启动 smoke |
| `v7_caseA_newpyd_long200000_smoke/` | 长时记录器 smoke |
| `v7_caseA_newpyd_long20000_after_tec_heatfix_short1000/` | TEC heatfix 后短时验证 |
| `v7_caseA_newpyd_long20000_after_tec_heatfix_startcheck/` | TEC heatfix 后启动检查 |
| `v7_caseA_newpyd_long20000_after_tec_heatfix_startcheck_bg/` | TEC heatfix 后后台启动检查 |

### 6.2 overnight / long 长时运行

| 目录 | 用途 |
|---|---|
| `v7_caseA_100k_run/` | `100000 s` 续算产物 |
| `v7_caseA_electric_dt08_outerCp001_overnight/` | 外围 `Cp x 0.01` 电气耦合过夜日志 |
| `v7_caseA_electric_dt08_physicalCp_overnight/` | 物理热容电气耦合过夜运行产物 |
| `v7_caseA_newpyd_long200000/` | 长时后台运行日志 |
| `v7_caseA_newpyd_long20000_from_t23800/` | 从已有绝对时间状态继续推进的长时产物 |

### 6.3 geometry、area-sync 和 TEC heatfix 验证

| 目录 | 用途 |
|---|---|
| `v7_caseA_area_sync_continue500/` | 面积同步修复续算验证 |
| `v7_caseA_geometry_fix_300step_test/` | 几何修复短时验证 |
| `v7_caseA_geometry_fix_continue2000/` | 几何修复续算验证 |
| `tec_joule_fix_audit_tmp/` | TEC 非均匀焦耳热修复审计 |
| `v7_caseA_newpyd_long20000_after_tec_heatfix/` | TEC heatfix 后主验证产物及嵌套审计结果 |
| `v7_caseA_newpyd_long20000_after_tec_heatfix_bgtest/` | TEC heatfix 后后台运行检查 |

### 6.4 interface / fluid-energy 审计

| 目录 | 用途 |
|---|---|
| `v7_caseA_interface_audit/` | v7 CaseA 界面能量审计结果 |

### 6.5 multiplier 参数扫描

| 目录 | 用途 |
|---|---|
| `v7_caseA_multiplier_short/` | 环乘子短时扫描结果 |

### 6.6 产物读取优先级

| 文件模式 | 用途 | 建议 |
|---|---|---|
| `*.npz` | 完整 restart 状态 | 需要续算、重建内存状态或精确审计时优先读取 |
| `*_latest_state.json` | 最新运行摘要 | 先快速判断进度、时间和主要诊断量 |
| `*_summary.json` | 审计、提取或扫描摘要 | 优先阅读结论，再按需进入 CSV |
| `*.csv` | 时序、节点、界面和分项数据 | 用于定位具体通道、节点或时间段 |
| `*.out`、`*.err`、`*.log` | 后台运行诊断 | 进程异常、无输出或数值中断时优先查看 |

`__pycache__/` 是 Python 缓存目录，不纳入产物索引。

## 7. 修改影响地图

| 修改区域 | 优先运行的快速测试 | 再做的定向验证 |
|---|---|---|
| TEC 电势、电场、焦耳热、面积基准 | `test_tec_joule_nonuniform.py`、`test_tecboundary.py`、`test_TEC_with_heat.py` | `run_v7_caseA_multipliers_short.py`、`extract_v7_caseA_electric_node_power.py`、`audit_v7_caseA_interface_energy.py` |
| TFE 或 TEC 电热耦合 | `test_tfe_thermal_flow.py`、`test_verify_solid_couplers.py` | `run_v7_caseA_multipliers_short.py`、`audit_v7_caseA_interface_energy.py` |
| 热管、翅片、集流环节点耦合 | `testHPwithFin.py`、`test_single_hp_fin_energy_conservation.py`、`test_single_hp_header_energy_audit.py` | `test_ringHP_buzzin.py`、`test_ringhp_node_coupling_energy_conservation.py` |
| 水力网络、泵、稳压器、拓扑 | `test_phase_topology.py`、`test_fixed_pressure_pump_hydrodynamics.py`、`test_pump_junction_hydraulic_network.py`、`test_pressurizer_volume.py` | `test_pressurizer_pumped_closed_loop.py`、`audit_v7_caseA_fluid_energy_network.py` |
| 流固换热和耦合器 | `test_simple_solid_fluid_couple.py`、`test_fluid_solid_couple.py`、`test_verify_solid_couplers.py` | `test_pipe_heat_transfer.py`、`audit_v7_caseA_interface_energy.py` |
| restart 保存和恢复 | `test_system_manager_lifecycle.py` | 用已有 `.npz` 运行 `run_v7_caseA_multipliers_short.py`，确认加载、流量重置和 TEC 刷新 |
| `SystemManager` 生命周期、步进顺序或耦合刷新 | `test_system_manager_lifecycle.py`、`test_system_manager.py` | `test_simple_solid_fluid_couple.py`、`audit_v7_caseA_fluid_energy_network.py` |
| v7 CaseA 几何、环乘子或全系统组装 | 相关组件测试、`run_v7_caseA_multipliers_short.py` | `audit_v7_caseA_interface_energy.py`、`audit_v7_caseA_fluid_energy_network.py`，最后再启动长时运行 |

长时运行脚本成本较高。除非任务本身要求长时数据，不要把 `run_v7_caseA_continue_100000s.py` 或过夜运行当作第一层验证。

## 8. 常见任务速查

| 任务 | 优先入口 |
|---|---|
| 修改当前 v7 CaseA 系统组装 | `test_core_assemble_v7_caseA.py` |
| 做低成本 v7 CaseA restart 验证 | `run_v7_caseA_multipliers_short.py` |
| 运行 FastSteady | `test_core_assemble_v7_caseA_faststeady.py`、`run_v7_caseA_faststeady_continue.py` |
| 启动长时分段记录 | `run_v7_caseA_newpyd_long200000.py` |
| 续算 `100000 s` 并施加专用导线电阻 | `run_v7_caseA_continue_100000s.py` |
| 审计界面能量 | `audit_v7_caseA_interface_energy.py` |
| 审计流体能量网络 | `audit_v7_caseA_fluid_energy_network.py` |
| 提取 TEC 节点电功率 | `extract_v7_caseA_electric_node_power.py` |
| 检查 `SystemManager` restart 和生命周期 | `test_system_manager_lifecycle.py` |

## 9. 使用注意

1. 多个 CLI 默认路径写成 `testModule/...`，通常应从项目根目录执行；运行前先检查脚本的默认 restart 和输出目录。
2. restart 必须与重建系统所使用的几何、乘子和求解器配置匹配。旧快照加载成功不等于物理配置一定一致。
3. `.json` 是摘要或审计输出，不替代 `.npz` 完整 restart。
4. 审计脚本依赖同步后的内存状态。直接读取 restart 数组不能替代 `sync_for_audit()` 的耦合刷新流程。
5. 修改核心求解器时，先跑影响地图中的窄测试，再进入 v7 CaseA smoke、审计和长时运行。
