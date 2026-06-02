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
| `test_single_tfe_energy_conservation_v7.py` | 从 FastSteady 快照选择性提取中心通道，运行严格绝热单 TFE 能量守恒诊断；先读 `SINGLE_TFE_ENERGY_CONSERVATION_GUIDE.md` |

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
| TFE 或 TEC 电热耦合 | `test_tfe_thermal_flow.py`、`test_verify_solid_couplers.py` | `test_single_tfe_energy_conservation_v7.py`、`run_v7_caseA_multipliers_short.py`、`audit_v7_caseA_interface_energy.py` |
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
| 运行严格绝热单 TFE 能量守恒基线 | 先读 `SINGLE_TFE_ENERGY_CONSERVATION_GUIDE.md`，再运行 `test_single_tfe_energy_conservation_v7.py` |
| 检查 `SystemManager` restart 和生命周期 | `test_system_manager_lifecycle.py` |

## 9. 使用注意

1. 多个 CLI 默认路径写成 `testModule/...`，通常应从项目根目录执行；运行前先检查脚本的默认 restart 和输出目录。
2. restart 必须与重建系统所使用的几何、乘子和求解器配置匹配。旧快照加载成功不等于物理配置一定一致。
3. `.json` 是摘要或审计输出，不替代 `.npz` 完整 restart。
4. 审计脚本依赖同步后的内存状态。直接读取 restart 数组不能替代 `sync_for_audit()` 的耦合刷新流程。
5. 修改核心求解器时，先跑影响地图中的窄测试，再进入 v7 CaseA smoke、审计和长时运行。

## 10. 单 TFE 守恒基线状态

2026-06-01 已完成严格绝热中心单 TFE 热工闭环和首次 TEC 基线。细节见 `SINGLE_TFE_ENERGY_CONSERVATION_GUIDE.md`。

- `thermal-baseline --duration-s 10 --max-dt-s 0.1`：最后 `1 s` 平均相对残差约 `0.279%`。
- 外套管热损失逐节点严格为 `0 W`。
- 普通固固、流固映射和流体矩阵残差均低于 `1e-6 W` 门槛。
- `tec --duration-s 1 --max-dt-s 0.01`：首次记录最后 `1 s` 平均相对残差约 `0.089%`，暂不固化 TEC 硬阈值。
- TEC 模式必须使用 ABI 匹配的 Python 3.12 和重建后的 `ThermoCalc/te_solver.cp312-win_amd64.pyd`。

本轮定位到 `HeatConduction.step()` 的 `solve_ivp(y0=self.T)` 初值别名问题。该入口已改为 `y0=self.T.copy()`；后续若调整固体积分器，必须保留这一约束。

## 11. 2026-06-02 TEC FVM 焦耳热基线

TEC 生产路径已改为使用 ThermoCalc C++ 输出的逐轴向节点焦耳热 `joulePowerE/C [W]`。Python 层只执行轴向列内二维体积加权分配，旧节点电势梯度法保留为诊断对照。

验证入口：

```text
test_tec_joule_nonuniform.py
test_thermocalc_interface.py
test_single_tfe_energy_conservation_v7.py --mode tec --duration-s 1 --max-dt-s 0.01
```

正式单 TFE `1 s` 结果：二维映射与 C++ 节点功率差为 `0 W`，TEC 转换闭合差为 `0.025015 W`，最终全局残差为 `0.089127 W`。后续若修改 TEC 电势离散或外层电路迭代阈值，必须继续运行该审计。

v7 CaseA `1 s` smoke 和静态接口审计已运行。普通固固接口累计误差约 `5.24e-9 W`，慢化剂映射最大误差约 `1.82e-12 W`；但多 TEC 串联电路会重复输出 `Failed to converge after 100000 iterations.`，静态审计即时 TEC 闭合差约 `5.20 W`。该问题属于后续电路收敛专项，不得通过调整热源映射掩盖。

## 12. 2026-06-02 v7 CaseA 全局储能连续审计

`audit_v7_caseA_global_storage_continue.py` 用于在不中断推进的同一内存状态中记录全堆芯有限差分储能。它不修改生产推进顺序，并在每个内部 `SystemManager.step()` 上对外部功率执行梯形积分。

全局审计同时输出三种残差：

```text
R_thermal_model =
    Q_core
  - Q_fluid_external_net_out
  - Q_radiation
  - Q_TEC_applied_removed_scaled_by_thermal_multiplier
  - dE_solid/dt
  - dE_fluid/dt

R_TEC_count =
    Q_core
  - Q_fluid_external_net_out
  - Q_radiation
  - Q_TEC_applied_removed_scaled_by_TEC_multiplier
  - dE_solid/dt
  - dE_fluid/dt

R_terminal =
    Q_core
  - Q_fluid_external_net_out
  - Q_radiation
  - P_terminal
  - dE_solid/dt
  - dE_fluid/dt
```

其中：

- 固体有限差分储能使用 `sum(0.5 * (cap_old + cap_new) * (T_new - T_old)) / dt`。
- 流体有限差分储能使用有限且非定压控制体上的 `sum(0.5 * (m_old + m_new) * (h_new - h_old)) / dt`。
- 流体外部焓流使用固定压力边界控制面，并复用 `HydraulicNetwork` 的迎风供体与倍率行定义。
- `R_thermal_model` 使用实际下发到代表热工网格的 TEC 净移热，并按热工倍率还原。
- `R_TEC_count` 按电路中的 TEC 数量还原，`R_terminal` 使用端电功率。三者差值用于识别部分 TEC 配置、热源松弛和电路闭合误差。
- `moderator_mapping_source_minus_boundary_out_w_avg` 记录全局慢化剂环源项与代表 TFE 慢化剂外流之间的即时差值，用于定位跨层映射时序滞后。
- `solid_equation_residual_w`、`fluid_equation_residual_w` 和 `fluid_solid_mapping_residual_w` 将 `R_thermal_model` 拆为固体区间闭合、流体区间闭合和流固映射区间差。三者之和记录为 `decomposed_thermal_model_residual_w`。

注意：`ReactorCore.pre_step()` 当前使用 `electric_alpha = alpha_tec * tec_mult / thermal_mult`。该比例只改变热源缓存趋近目标值的速度，不改变最终目标幅值。对于热工倍率与 TEC 倍率不同的代表环，必须同时检查 `tec_thermal_model_minus_electric_count_w`，不能默认认为部分 TEC 配置已经在稳态幅值上做了平均。

首次连续审计结果：

- 从 `t=36940 s` restart 连续推进 `100 s`，每 `10 s` 记录一次。restart 后第一段存在一次性外层固体暖启动跳变，不应纳入稳态残差统计。
- 进一步从新 restart 连续推进 `20 s`，每 `2 s` 记录一次。去掉第一段暖启动后，`global_residual_using_thermal_model_tec_heat_w` 为 `3.08-5.23 kW`。
- `tec_thermal_model_minus_electric_count_w` 稳定约为 `525.5 W`，说明热工代表元件实际源项倍率与电路 TEC 数量倍率不同。
- `moderator_mapping_source_minus_boundary_out_w_avg` 从约 `-2.90 kW` 变化到 `-0.62 kW`，说明全局慢化剂跨层映射存在明显时序滞后。
- 当前不应把旧快速口径中的 `1%-2%` 直接视为导热泄漏，也不应直接启动过夜长算。后续应先修复部分 TEC 热源幅值平均和慢化剂映射时序，再重复本审计。

推荐先运行：

```powershell
python testModule\audit_v7_caseA_global_storage_continue.py `
  --restart-in <restart.npz> `
  --duration 100 `
  --record-interval 10 `
  --max-dt 0.8 `
  --inner-iter 1
```

## 13. 2026-06-02 v7 CaseA 慢化剂映射时间层修复

全局储能审计进一步确认，旧版 `ReactorCore.pre_step()` 在调用 `TFEUnit.pre_step()` 更新慢化剂外边界温度之前，就读取了内部等效 moderator 的 `BoundaryRegion.current_flux`。这会使用旧时间层缓存：

- 正常连续推进中，`dt=0.1 s` 探针仍可观察到约 `100 W` 的慢化剂映射时序差。
- 从 restart 重建后，首步曾把错误缓存 `+790.2 kW` 注入全局慢化剂环；真实边界在首步后立即翻转到约 `-12.0 kW`。
- 首区间的巨大外层固体储能跳变主要集中在四个全局慢化剂环，不是外边界辐射损失。

生产逻辑现已改为：先执行所有 `TFEUnit.pre_step()`，再刷新内部等效 moderator 的物性、热阻、边界状态和热流缓存，最后按 `tfe_multipliers` 聚合到全局慢化剂环。后续 restart 审计仍应单独检查第一段，不得默认忽略首段异常。

补充检查：

- 解析辐射功率与反射层实际 `BoundaryRegion.current_flux` 差约 `1.4e-11 W`。
- 堆芯重构热功率与芯块 `Q_source` 均为 `115000 W`。
- 修复后普通固固接口累计残差约 `5.27e-9 W`，流固边界放热与注入流体源项差约 `2.2e-10 W`，同步后的慢化剂映射最大误差约 `2.84e-14 W`。
- 水力倍率行无不匹配，有限非定压控制体矩阵残差维持在 `1e-9 W` 量级。`audit_v7_caseA_fluid_energy_network.py` 中先汇总绝对焓能再相减的字段会受浮点消差影响；判断水力闭合时应优先使用逐控制体矩阵残差。

修复后继续用 `max_dt=0.1 s` 运行分解审计，`R_thermal_model` 从约 `4.14 kW` 降到 `3.72 kW`。三项分解表明：

- `solid_equation_residual_w` 约为 `-0.26` 至 `-0.33 kW`。
- `fluid_equation_residual_w` 除首段外接近 `0 W`。
- `fluid_solid_mapping_residual_w` 约为 `4.03` 至 `4.32 kW`，是当前主要剩余项。

同一起点、单个 `0.1 s` 步长的 Picard 对照中，流固差随 `inner_iter=1 / 2 / 3 / 5 / 10` 依次约为 `3.84 / 3.23 / 2.65 / 1.80 / 0.68 kW`。`inner_iter=20` 时在第 `16` 次按默认温度阈值提前收敛，流固差约 `0.21 kW`。审计脚本提供 `--inner-iter` 用于受控对照；生产默认值仍保持 `1`，后续应评估共享区间换热量或正式 Picard 策略，不能仅靠提高默认迭代次数掩盖成本。

## 14. 2026-06-02 v7 CaseA 径向热流链审计

`audit_v7_caseA_global_storage_continue.py` 现已同时记录径向热流链。正号统一表示径向向外传热；如果字段为负值，则热流实际沿径向向内。普通固固接口同时记录左侧外流、右侧流入和两侧差值。TEC 接口不是普通守恒导热面，必须同时检查：

- 芯块到发射极边界热流。
- 发射极侧外流和接收极侧流入。
- 发射极电子冷却、接收极电子加热、发射极焦耳热、接收极焦耳热和端电功率。
- 接收极到内套管、内套管到流体、流体到外套管和外套管到虚拟慢化剂的边界热流。
- 虚拟慢化剂外流、真实慢化剂环源项及其映射残差。
- 真实慢化剂到套筒、套筒到反射层和反射层外边界辐射热流。

脚本还将固体有限差分储能拆为芯块、发射极、接收极、内套管、外套管、虚拟慢化剂、真实慢化剂、套筒和反射层，并输出各层 `chain_balance_*_residual_w`。各层口径为“流入 + 内热源 - 流出 - 有限差分储能”。外部功率在每个内部 `SystemManager.step()` 上采用梯形积分，而固体由 BDF 推进，因此逐层残差用于定位误差集中区域，不等同于严格的 BDF 离散方程残差。严格闭合仍应结合普通接口两侧差、同步后的映射误差和流体逐控制体矩阵残差判断。
