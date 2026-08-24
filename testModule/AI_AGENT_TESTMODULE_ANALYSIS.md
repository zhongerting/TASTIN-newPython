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

## 15. 2026-06-02 V8 CaseA 外圈代表元件拆分

V8 保持 V7 CaseA 的几何、水力边界、总热工元件数、四个真实慢化剂环和 `115000 W` 定功率口径，只把外圈代表元件拆分：

| 代表 TFE | 热工倍率 | TEC 倍率 | 真实慢化剂环 |
|---|---:|---:|---:|
| `Center` | 1 | 1 | 0 |
| `Ring1` | 6 | 6 | 1 |
| `Ring2` | 12 | 12 | 2 |
| `Ring3_TEC` | 15 | 15 | 3 |
| `Ring3_Open` | 3 | 0 | 3 |

`Ring3_Open` 保留铯隙被动导热和辐射，但电子冷却、电子加热、焦耳热、电流密度和电诊断缓存必须严格清零。组件层通过 `TFEUnit.clear_tec_sources()` 实现该契约，`ReactorCore.pre_step()` 每步对 `tec_multiplier == 0` 的代表件重新执行清零。

旧 V7 restart 不能直接加载到 V8，因为水力拓扑由 `176/179` 个控制体/连接扩展到 `213/217`。首次运行应执行：

```powershell
python testModule\migrate_v7_caseA_restart_to_v8.py
```

迁移器按对象语义复制公共状态，并将原 `Ring3` 热工和水力状态复制给两个新外圈代表件。迁移完成后必须先刷新固体边界缓存，再向 ThermoCalc 下发电极表面温度；否则 C++ 电路会收到占位温度并产生非有限焦耳热。

V8 正式长算入口：

```powershell
python testModule\run_v8_caseA_long_energy.py `
  --restart-in testModule\v8_caseA_migrated\v8_caseA_migrated_latest_restart.npz `
  --duration 20000 `
  --record-interval 200 `
  --restart-interval 200 `
  --max-dt 0.8 `
  --inner-iter 1 `
  --solid-ode-method LSODA
```

每 `200 s` 覆盖 latest restart，并向固定 CSV 续写完整径向热流链、分层储能、三种全局残差和 `Ring3_Open` 零源检查。restart 后再次启动时，运行器从 CSV 恢复统一记录起点，保持相对时间连续。V8 公共入口默认 `solid_ode_method = LSODA`，并在构建后和加载 restart 后对 `SystemManager.solid_components` 中的全部 36 个堆芯固体调用 `set_ode_method()`；`latest_state.json` 会记录 `solid_ode_method` 和逐固体 `solid_ode_methods`。

2026-06-10 补充：`run_v8_caseA_common.py` 现在会在 V8 加载运行中向每个 ThermoCalc virtual element 施加导线电阻 `[0.00155199999999970, 0.00102400000000000, 0.000336000000000000, 0.000608000000000000] ohm`，并在 `latest_state.json` 记录 `wire_resistance_ohm`。该电阻应在 `system.load_global_state()`、`core.setup_tec_circuit("fixed_u", 27.2, I_guess=150.0)` 和 `core.post_step(...)` 同步电极温度之后重建并计算 ThermoCalc；若在温度同步前施加非零导线电阻，首次 `thermo_calc.calculate()` 可能极慢。已从绝对时间 `39184 s` 继续计算 `2000 s` 到 `41184 s`，结果保存在 `testModule/v8_caseA_lsoda_wire_2000s/`：末段端电功率 `4873.557 W`，TEC applied removed heat `5760.694 W`，冷却剂焓升 `108524.931 W`，combined storage `-16.346 W`，thermal-model residual `-0.335 W`。terminal-power residual 约 `+886.802 W`，基本等于 TEC applied removed heat 与端电输出之差，说明带导线电阻后该口径包含当前尚未单独沉积的导线/端部损耗，需与 thermal-model residual 分开解读。

V8 定向审计入口：

```text
test_v8_caseA_topology.py
audit_v8_caseA_global_storage_continue.py
audit_v8_caseA_interface_energy.py
audit_v8_caseA_fluid_energy_network.py
```

首次短验证结果：五个代表 TFE、四个真实慢化剂环、`37/34` 倍率和 `Ring3_Open` 六项 TEC 零源检查均通过；`QTEC,thermal - QTEC,count = 0 W`。静态普通固固接口累计误差约 `5.29e-9 W`，同步慢化剂映射最大误差约 `7.11e-15 W`，流体逐控制体矩阵残差保持在 `1e-9 W` 附近。2026-06-09 的 `LSODA` smoke 中，`Failed to converge after 100000 iterations.` 共出现 `4440` 次，全部发生在进入 V8 长算主循环前的 `core.thermo_calc.calculate(verbose=False)` 阶段；`system.step(0.01)`、记录输出和结束阶段未再出现。该信息来自 ThermoCalc C++ 电路迭代，不是 `BaseHeatConduction` 的固体 ODE 收敛失败，该电路收敛专项尚未修复。

迁移后 `100 s` 受控延续结果：

- `QTEC,thermal - QTEC,count` 全程为 `0 W`，V7 中约 `527 W` 的倍率混用差已消除。
- 最后 `10 s` 的 `R_thermal_model = -1864.39 W`，其中固体区间残差 `-187.11 W`、流体方程残差 `-0.28 W`、流固时间层差 `-1677.01 W`。
- 最终静态普通固固接口累计误差约 `5.19e-9 W`，同步慢化剂映射最大误差约 `3.55e-15 W`。
- V8 不再被 TEC 倍率口径差掩盖；当前剩余项主要来自迁移后外圈热状态重新分化期间的流固时间层差。正式 `20000 s` 长算应继续记录其趋势，不得把某一时刻的全局抵消值直接视为严格闭合。

## 16. 2026-06-10 V8 CaseA NaK78 coolant migration

V8 CaseA now defaults to `SodiumPotassium78` coolant in `test_core_assemble_v8_caseA.py`, while V7 CaseA keeps its default `Sodium` behavior. The shared V7 builder accepts `coolant_material` and still stores the selected coolant object under the legacy `"Sodium"` material key so existing TFE and hydraulic construction code remains compatible.

Do not directly continue an old V8 Sodium restart with the NaK78 model. Use:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\migrate_v8_caseA_sodium_restart_to_nak.py
```

The migration builds the current V8 NaK78 model, loads `testModule/v8_caseA_lsoda_wire_2000s/v8_caseA_lsoda_wire_2000s_latest_restart.npz`, recomputes `HydraulicNetwork.h_vec` from saved `T_vec/P_vec` using the current NaK78 material, refreshes fluid properties, reapplies LSODA, fixed voltage, design flows, `115000 W` core power and the standard wire resistance, then writes:

- `testModule/v8_caseA_nak_migrated/v8_caseA_nak_migrated_latest_restart.npz`
- `testModule/v8_caseA_nak_migrated/v8_caseA_nak_migrated_migration_summary.json`

The NaK78 2000 s continuation was then run with:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v8_caseA_long_energy.py `
  --restart-in testModule\v8_caseA_nak_migrated\v8_caseA_nak_migrated_latest_restart.npz `
  --duration 2000 `
  --record-interval 200 `
  --restart-interval 200 `
  --max-dt 0.8 `
  --output-dir testModule\v8_caseA_nak_wire_2000s `
  --case-prefix v8_caseA_nak_wire_2000s `
  --solid-ode-method LSODA
```

Final state: absolute time `43184.00000000066 s`, relative time `2000 s`, terminal electric power `4868.580 W`, coolant enthalpy rise `108560.205 W`, coolant solid-to-fluid heat `108560.067 W`, outer-wall radiation `729.210 W`, combined storage `-11.771 W`, thermal-model residual `-0.132 W`, terminal-power residual `853.775 W`. Rebuilt final-state inlet/outlet plenum temperatures are `743.004650 K` and `838.785235 K`, so `DeltaT = 95.780585 K`; the NaK78 enthalpy rise is `83507.805 J/kg`, matching the `1.30 kg/s` loop heat pickup. `run.err` is empty and `run.out` contains no `Failed to converge`, `Traceback`, or `Error`.

`run_v8_caseA_long_energy.py` currently writes restart, JSON and CSV products only; it does not append `Diag/` fields into the `.npz` restart. The `Diag/` append behavior belongs to the separate v7 100000 s continuation path unless a future task explicitly extends it to V8.

## 17. 2026-06-11 V9 CaseA open external-piping skeleton

V9 CaseA is implemented in `testModule/test_core_assemble_v9_caseA.py` and `testModule/run_v9_caseA_open_loop.py`. It does not change the V8 default builder or V8 restart path. V9 first builds the V8 core and TFE/TEC stack, then replaces the pre-initialization hydraulic network with an open external-piping skeleton for the pre-collector-ring integration stage.

V9 keeps the V8 representative core (`Center/Ring1/Ring2/Ring3_TEC/Ring3_Open`), NaK78 coolant, fixed-voltage TEC mode, LSODA solid option, `Ring3_Open` zero-TEC-source contract, and standard wire resistance. The external hydraulic path is:

```text
fixed-flow inlet
  -> V9_InletConnector
  -> RadiatorOutletBranch_38
     + RadiatorOutletBranch_44_50_Rep(multiplier=2)
  -> RadiatorInnerHeader_53
  -> RadiatorOuterHeader_52
  -> ColdReturnBranch_1
     + ColdReturnBranch_2_3_Rep(multiplier=2)
  -> CoreInletConnector
  -> V8/V9 representative TFE channels
  -> CoreOutletConnector
  -> HotOutletBranch_1/2/3
  -> fixed-pressure outlet
```

Important constraints:

- V9 has no pump, no pressurizer, no collector-ring heat pipes, and no first-round local-loss K values; all new external junctions use `k_loss=0.0`.
- The inlet is a fixed mass-flow `InletJunction` with default `1.3 kg/s` and `743 K`; the only fixed pressure boundary is the outlet, default `160000 Pa`.
- #38 is explicit, while #44/#50 are represented by `RadiatorOutletBranch_44_50_Rep(multiplier=2)` because their lengths differ from #38.
- The cold return #2/#3 pair uses one representative branch with `multiplier=2`; the three hot outlet branches stay explicit for later connection to collector-ring `I1/I2/I3`.
- `CoreInletConnector` and `CoreOutletConnector` are numerical mixing/connection nodes, not physical inlet/outlet boxes.
- V8 restart files cannot be loaded directly into V9 because the hydraulic topology and fixed-pressure-boundary set are different. Use only V9 restart files with `run_v9_caseA_open_loop.py --restart-in ...`.

V9-to-collector-ring integration note: the current V9 hot outlet branches are macro-system branches. The current 6-segment collector-ring model explicitly builds one physical collector ring and uses `MacroFlowJunction(multiplier=2)` to represent the second symmetric ring. Therefore the future integrated open chain must not connect `HotOutletBranch_1/2/3` directly to the single-ring `I1/I2/I3` nodes. Each hot branch should connect through a macro-to-single-ring split, so a macro branch flow of about `0.433 kg/s` becomes about `0.2167 kg/s` into the explicitly modeled ring. The `O1/O2/O3` side needs the matching single-ring-to-macro merge before returning to V9 cold-side piping. Building two explicit collector rings is possible later, but the first integrated version should keep the one-ring-plus-`multiplier=2` convention used by `CoolantLoop`.

Recommended first checks:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\test_core_assemble_v9_caseA.py testModule\run_v9_caseA_open_loop.py testModule\test_v9_caseA_topology.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_v9_caseA_topology
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v9_caseA_open_loop.py --duration 0.1 --record-interval 0.1 --restart-interval 0.1 --max-dt 0.01 --disable-tec-coupled
```

Cold-start V9 hydraulic/topology smoke should use `--disable-tec-coupled`, because the uniform cold initial TFE state can make the first ThermoCalc solve return non-finite Joule heat before a migrated warm thermal state is available. TEC-coupled V9 runs should start from a V9-compatible warm restart or a future dedicated V8-to-V9 thermal migration path.

## 18. 2026-06-12 V10 CaseA integrated open loop

V10 CaseA is implemented in `testModule/test_core_assemble_v10_caseA.py` and `testModule/run_v10_caseA_open_loop.py`. It is the first integrated open-loop chain that combines the V9 core/external-piping skeleton with the current 6-segment one-ring collector model from `CoolantLoop/model_collector_ring_6segment_v9_interface.py`.

The hydraulic path is:

```text
fixed-flow inlet
  -> CoreInletConnector
  -> V8/V9 representative TFE channels
  -> CoreOutletConnector
  -> HotOutletBranch_1/2/3
  -> MacroFlowJunction(multiplier=2) into InletMix_I1/I2/I3
  -> A1..A6 one explicit collector ring with RingHP heat pipes
  -> OutletMix_O1/O2/O3
  -> Manifold_1/2/3
  -> MacroFlowJunction(multiplier=2) into V10_RadiatorManifoldMerge
  -> RadiatorInnerHeader_53
  -> RadiatorOuterHeader_52
  -> ColdReturnBranch_1 + ColdReturnBranch_2_3_Rep(multiplier=2)
  -> fixed-pressure outlet
```

Important constraints:

- V10 keeps V9's open-loop boundary convention: the inlet is fixed mass flow and the only fixed pressure reference is `V10_OutletBoundary_FixedPressure`; it does not add a pump or pressurizer.
- The V9 duplicate radiator-outlet-side branches before the collector ring are removed. The collector ring supplies `InletMix_I*`, `A1..A6`, `OutletMix_O*`, and `Manifold_*`.
- One explicit collector ring represents two symmetric physical rings. Hot macro branches of about `0.433 kg/s` enter the explicit ring through `MacroFlowJunction(multiplier=2)` at about `0.2167 kg/s`; the manifold side uses the matching single-ring-to-macro merge.
- `reset_v10_design_flows(..., preserve_ring_restart_flows=True)` must be used after V9/ring restart injection and after V10 restart loading. It resets the external open-chain boundary/core/cold-return flows while preserving collector-ring internal `A*`, `J_I*/J_A*`, `OutletMix -> Manifold`, and `Manifold_*_Junc_*` flow distribution copied from the ring restart.
- V10 restart files are topology-specific. Do not load V8 or V9 restart files directly through `--restart-in`; use the V10 runner's V9+ring injection path to create an initial V10 restart first.

Initial restart creation:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v10_caseA_open_loop.py `
  --create-init-only `
  --output-dir testModule\v10_caseA_open_loop_init `
  --case-prefix v10_caseA_open_loop_init `
  --max-dt 0.05
```

The 2026-06-12 initialization copied 36 core solids from `testModule/v9_caseA_open_loop_tec_3000s/v9_caseA_open_loop_tec_3000s_latest_restart.npz`, 24 collector-ring solids from `CoolantLoop/collector_ring_6segment_v9_interface_500s_from200s_restart.npz`, 211 V9 fluid volumes, and 63 ring fluid volumes. The initial V10 absolute time is `5010 s`.

Verified checks:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\test_core_assemble_v10_caseA.py testModule\run_v10_caseA_open_loop.py testModule\test_v10_caseA_topology.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_v10_caseA_topology
```

Both passed. A TEC-coupled `40 s` smoke from the injected V10 restart also completed with `max_dt=0.05 s`; the final record was approximately `Tin=743.000 K`, `Tcore_out=838.918 K`, `Tring_out=762.823 K`, and `Pel=4874.558 W`. A later `max_dt=0.1 s` continuation did not numerically fail, but was too slow for the 300 s command timeout and only recorded the first `20 s`. Current V10 long runs should therefore use conservative time steps and be treated as expensive until the integrated hydraulic/thermal startup is further optimized.

## 19. 2026-06-15 V10 optional local implicit fluid-solid exchange

`testModule/run_v10_caseA_open_loop.py` now accepts:

```powershell
--fluid-solid-coupling-scheme current
--fluid-solid-coupling-scheme local_implicit
```

The default is `current`, which preserves previous V10 behavior and remains compatible with old restart files. `local_implicit` switches all `FluidSolidCouple` objects that provide `solid_node_capacitance` to the local two-capacitance exchange scheme implemented in `Solvers/Couplers.py`. The runner records `fluid_solid_coupling_scheme` and `fluid_solid_coupler_count` in `latest_state.json`; history rows record the selected scheme.

Recommended first V10 local-implicit smoke:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v10_caseA_open_loop.py `
  --restart-in testModule\v10_caseA_open_loop_rk45_1000s\v10_caseA_open_loop_rk45_1000s_latest_restart.npz `
  --duration 100 `
  --record-interval 50 `
  --restart-interval 50 `
  --max-dt 0.1 `
  --solid-ode-method RK45 `
  --fluid-solid-coupling-scheme local_implicit
```

Treat `max_dt=0.5 s` as an exploratory setting, not an acceptance criterion, because the previous explicit-coupling V10 run developed hydraulic nonconvergence near absolute time `12237.6 s`.

## 20. 2026-06-15 V10 727 K radiator tuning preset

`testModule/run_v10_caseA_open_loop.py` now supports a 727 K radiator tuning preset:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v10_caseA_open_loop.py `
  --preset-727-radiator-tuning `
  --restart-in <v10_restart.npz> `
  --duration 500 `
  --record-interval 50 `
  --restart-interval 50
```

The preset sets `inlet_temperature_k=727.0`, `wire_resistance_scale=0.5`, `ring_emissivity=0.15`, `outer_header_emissivity=0.2`, `hp_emissivity=0.6`, `fin_emissivity=0.6`, `solid_ode_method=RK45`, `fluid_solid_coupling_scheme=local_implicit`, and `max_dt=0.1`.

Implementation notes:

- The inlet boundary temperature is re-applied after V10 restart loading or V9/ring restart injection, because restart files contain the previous boundary temperature.
- `RadiatorOuterHeader_52` remains a fluid channel; its radiation is an equivalent fluid-side sink `epsilon*sigma*A*(T^4-T_space^4)` applied in `pre_step()`.
- `ring_emissivity`, `hp_emissivity`, `fin_emissivity`, `outer_header_emissivity`, `outer_header_radiation_w`, `ring_wall_radiation_w`, and `hp_fin_radiation_w` are written to diagnostics/history.
- First tune with HP/fin emissivity fixed at `0.6`; adjust both together only after the 727 K case reaches a relative steady trend.

Tuning results as of 2026-06-15:

- `ring_emissivity=0.20`, `hp_emissivity=0.70`, `fin_emissivity=0.70` reached near steady state after `1000 s`; `RadiatorOuterHeader_52_out` was about `732.45 K`, still about `5.45 K` above the 727 K inlet target.
- `ring_emissivity=0.20`, `hp_emissivity=0.75`, `fin_emissivity=0.75` is the current preferred candidate. After the latest continuation to absolute time about `15530 s`, `RadiatorOuterHeader_52_out` was about `728.65 K`, while the core inlet boundary and connector stayed at about `727.000 K`.
- `ring_emissivity=0.20`, `hp_emissivity=0.80`, `fin_emissivity=0.80` over-cooled the outer header outlet to about `725.05 K`; use it only as an upper bracket.

The current `0.75/0.75` candidate restart is:

```text
testModule/v10_caseA_tune727_ring020_hpfin075_steady_plus1000s/v10_caseA_tune727_ring020_hpfin075_steady_plus1000s_latest_restart.npz
```

For future closed-loop construction, the pressure audit from the `0.75/0.75` near-steady case gives:

```text
CoreInletConnector pressure  ≈ 166471.52 Pa
CoreOutletConnector pressure ≈ 162614.31 Pa
Core pressure drop           ≈   3857.21 Pa
ColdReturnOutletMerge pressure ≈ 160004.97 Pa
Suggested initial pump head from ColdReturnOutletMerge to CoreInletConnector ≈ 6.5 kPa at 1.3 kg/s
```

The detailed pressure summary for the earlier `0.75/0.75` steady checkpoint is in the generated run directory as `pressure_summary.json`; run artifacts remain untracked and should not be committed unless explicitly requested.

## 21. 2026-06-15 V11 CaseA closed pumped loop

V11 CaseA is implemented in `testModule/test_core_assemble_v11_caseA.py` and `testModule/run_v11_caseA_closed_loop.py`. It reuses the V10 core, collector-ring, radiator-header, cold-return geometry, NaK78 coolant, tuned radiator emissivities, half wire resistance, RK45 solids, and local implicit fluid-solid coupling, but replaces the V10 open inlet/outlet boundaries with a closed pumped loop.

The V11 hydraulic path is:

```text
CoreInletConnector passive pressure reference
  -> V8/V9 representative TFE channels
  -> CoreOutletConnector
  -> HotOutletBranch_1/2/3
  -> MacroFlowJunction(multiplier=2) into InletMix_I1/I2/I3
  -> A1..A6 one explicit collector ring with RingHP heat pipes
  -> OutletMix_O1/O2/O3
  -> Manifold_1/2/3
  -> MacroFlowJunction(multiplier=2) into V10_RadiatorManifoldMerge
  -> RadiatorInnerHeader_53
  -> RadiatorOuterHeader_52
  -> PumpJunction A
  -> V11_PumpMidNode
  -> PumpJunction B
  -> V11_PumpOutletDistributor_51
  -> ColdReturnBranch_1 + ColdReturnBranch_2_3_Rep(multiplier=2)
  -> V10_ColdReturnOutletMerge
  -> CoreInletConnector
```

Important constraints:

- `V10_InletBoundary_FixedFlow`, `V10_OutletBoundary_FixedPressure`, `J_V10_InletBoundary_to_CoreInletConnector`, and `J_ColdReturnOutletMerge_to_OutletBoundary` are absent in V11.
- V11 has no fixed-pressure boundary volumes. `CoreInletConnector` is marked with `is_pressure_reference=True` and `target_P` copied from the V10 restart, so pressure is anchored without freezing the node enthalpy or temperature.
- The two pumps are series `PumpJunction`s between `RadiatorOuterHeader_52` and `V11_PumpOutletDistributor_51`. Default total pump head is `6466.56 Pa`, split equally as `3233.28 Pa` per pump.
- V11 must not run a cold hydraulic initialization before V10 state injection. The runner first maps the latest V10 state into V11, then rebuilds the pump-after cold-return pressure field so the pump outlet distributor is about one pump total head above `RadiatorOuterHeader_52`.
- `--enable-pump-head-control` optionally adjusts total pump head at `--pump-control-interval` using a bounded quadratic flow/head correction toward `--target-flow-kg-s` (default `1.3 kg/s`).

Default V10 injection restart:

```text
testModule/v10_caseA_tune727_ring020_hpfin075_steady_plus1000s/v10_caseA_tune727_ring020_hpfin075_steady_plus1000s_latest_restart.npz
```

Verified checks:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\test_core_assemble_v11_caseA.py testModule\run_v11_caseA_closed_loop.py testModule\test_v11_caseA_topology.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_v11_caseA_topology
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v11_caseA_closed_loop.py --create-init-only --output-dir testModule\v11_caseA_closed_loop_init2 --case-prefix v11_caseA_closed_loop_init2
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v11_caseA_closed_loop.py --restart-in testModule\v11_caseA_closed_loop_init2\v11_caseA_closed_loop_init2_latest_restart.npz --duration 10 --record-interval 5 --restart-interval 5 --max-dt 0.05 --output-dir testModule\v11_caseA_closed_loop_smoke2_10s --case-prefix v11_caseA_closed_loop_smoke2_10s
```

The corrected `10 s` smoke completed. The final record was approximately `Tin=727.000 K`, `Tcore_out=823.955 K`, `RadiatorOuterHeader_52_out=728.154 K`, `Wpump=1.297338 kg/s`, total pump head `6466.56 Pa`, and `Pel=5527.353 W`.

## 22. 2026-06-16 V12 CaseA open core + TOPAZ-II pipe-fin radiator

V12 CaseA is implemented in `testModule/test_core_assemble_v12_caseA.py` and `testModule/run_v12_caseA_open_loop.py`. It is an independent open-loop case definition and does not call the V10 system builder. The core side reuses the V8 representative-TFE assembly as the common core component source, but V12 defines its own external hydraulic path and replaces the V10 collector-ring/radiator-header chain with the TOPAZ-II 78-tube pipe-fin radiator model from `Components/RadiatorPipeWithFin.py`.

The V12 hydraulic path is:

```text
V12_InletBoundary_FixedFlow
  -> Pipe11_CoreInletHeader       (#11, flow-network core inlet header)
  -> V12_CoreInletDistribution
  -> V12_CoreInletBranch_1 + V12_CoreInletBranch_2_3_Rep(multiplier=2)
  -> V12_CoreInletConnector
  -> V8/V9 representative TFE channels, TEC disabled
  -> V12_CoreOutletConnector
  -> Pipe05_CoreOutletToRadiator  (#5, flow-network core outlet header)
  -> V12_RadiatorInletSplit
  -> 78 TOPAZ-II pipe-fin radiator tubes with upper/lower ring headers
  -> V12_RadiatorOutletMix
  -> Pipe06_RadiatorOutlet        (#6)
  -> Pipe07_HeatExchangerHotSide  (#7, simplified as one pipe segment)
  -> Pipe08_ReturnInnerPipe       (#8)
  -> Pipe09_ValveSegment          (#9, no pump in V12)
  -> V12_OutletBoundary_FixedPressure
```

Important constraints:

- V12 remains an open loop: the fixed-flow/fixed-temperature boundary is upstream of flow-network pipe #11, and the only fixed-pressure boundary is `V12_OutletBoundary_FixedPressure` after pipe #9.
- V12 defaults to no-TEC. Passing `--enable-tec-coupled` uses TEC multipliers `1,6,12,15,0`, fixed voltage `--target-voltage` default `27.2 V`, and `--thermo-update-interval` default `0.5 s`.
- The radiator defaults to 78 tubes, 8 axial fluid nodes per tube, `K_in=100`, `K_out=100`, `tube_emissivity=0.80`, `fin_emissivity=0.80`, and `fin_area_scale=0.35`.
- The runner re-applies `--inlet-temperature-k` after hydraulic initialization and after `--restart-in` loading, so the fixed-flow inlet boundary is not overwritten by restart temperature state.
- `--fluid-solid-coupling-scheme current` is the default. If `local_implicit` is requested, the runner initializes with `current` first and switches compatible fluid-solid couplers after hydraulic initialization, because `SystemManager.initialize_system()` does not pass a positive coupling `dt`.
- Flow-network pipe dimensions are stored in `latest_state.json` under `flow_network_pipe_specs`. The current values follow the supplied flow-network summary where usable; inconsistent diameter fields were kept conservative for the primary small pipes (`#5/#6/#11` use `dh=0.014 m`, `area=3.8e-4 m2`).

Verified checks:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\test_core_assemble_v12_caseA.py testModule\run_v12_caseA_open_loop.py testModule\test_v12_caseA_topology.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_v12_caseA_topology
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --create-init-only --output-dir testModule\v12_caseA_open_loop_no_tec_init --case-prefix v12_caseA_open_loop_no_tec_init --max-dt 0.05
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_no_tec_init\v12_caseA_open_loop_no_tec_init_latest_restart.npz --duration 1 --record-interval 0.5 --restart-interval 0.5 --max-dt 0.05 --output-dir testModule\v12_caseA_open_loop_no_tec_smoke_1s --case-prefix v12_caseA_open_loop_no_tec_smoke_1s
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_no_tec_smoke_1s\v12_caseA_open_loop_no_tec_smoke_1s_latest_restart.npz --duration 9 --record-interval 1 --restart-interval 3 --max-dt 0.05 --output-dir testModule\v12_caseA_open_loop_no_tec_smoke_10s --case-prefix v12_caseA_open_loop_no_tec_smoke_10s
```

All checks completed. The 10 s no-TEC smoke reached absolute time `10.0 s`. Final key values were approximately `Tin=753.327 K`, `Tcore_out=749.605 K`, `Tradiator_out_mix=762.519 K`, `radiator_tube_total_flow=1.294995 kg/s`, `tube_mean=0.016603 kg/s`, `tube_min=0.015862 kg/s`, `tube_max=0.018733 kg/s`, and `Qrad=99.783 kW`. The early no-TEC transient is still dominated by initial stored heat in the radiator and flow-network pipes, so this short run is a topology/numerics smoke, not a steady thermal-performance result.

2026-06-17 inlet-temperature restart handling was corrected and the 727 K no-TEC case was rerun from a fresh init:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --create-init-only --inlet-temperature-k 727 --output-dir testModule\v12_caseA_open_loop_no_tec_init_727K --case-prefix v12_caseA_open_loop_no_tec_init_727K --max-dt 0.05
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_no_tec_init_727K\v12_caseA_open_loop_no_tec_init_727K_latest_restart.npz --inlet-temperature-k 727 --duration 500 --record-interval 50 --restart-interval 50 --max-dt 0.05 --output-dir testModule\v12_caseA_open_loop_no_tec_727K_500s --case-prefix v12_caseA_open_loop_no_tec_727K_500s
```

The corrected 727 K run reached absolute time `500.0 s`. Final key values were approximately `Tin=727.000 K`, `Tcore_out=824.394 K`, `Tradiator_in=824.394 K`, `Tradiator_out_mix=727.094 K`, `core_delta_p=3857.04 Pa`, `coolant_enthalpy_rise=110.427 kW`, `Qrad=110.083 kW`, `radiator_tube_total_flow=1.300053 kg/s`, `tube_mean=0.0166673 kg/s`, `tube_min=0.0159407 kg/s`, and `tube_max=0.0187879 kg/s`.

TEC-coupled V12 was then enabled from the 727 K no-TEC `500 s` checkpoint. A `1 s` `max_dt=0.05 s` smoke and a `5 s` `max_dt=0.5 s` smoke both completed, then the run continued for `994 s` with `max_dt=0.5 s`, so TEC-coupled physical time totaled `1000 s` and the final absolute time was `1500 s`.

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_no_tec_727K_500s\v12_caseA_open_loop_no_tec_727K_500s_latest_restart.npz --inlet-temperature-k 727 --enable-tec-coupled --thermo-update-interval 0.5 --duration 1 --record-interval 0.5 --restart-interval 0.5 --max-dt 0.05 --output-dir testModule\v12_caseA_open_loop_tec_smoke_1s_from727K500s --case-prefix v12_caseA_open_loop_tec_smoke_1s_from727K500s
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_tec_smoke_1s_from727K500s\v12_caseA_open_loop_tec_smoke_1s_from727K500s_latest_restart.npz --inlet-temperature-k 727 --enable-tec-coupled --thermo-update-interval 0.5 --duration 5 --record-interval 1 --restart-interval 1 --max-dt 0.5 --output-dir testModule\v12_caseA_open_loop_tec_dt05_smoke_5s --case-prefix v12_caseA_open_loop_tec_dt05_smoke_5s
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_tec_dt05_smoke_5s\v12_caseA_open_loop_tec_dt05_smoke_5s_latest_restart.npz --inlet-temperature-k 727 --enable-tec-coupled --thermo-update-interval 0.5 --duration 994 --record-interval 50 --restart-interval 50 --max-dt 0.5 --output-dir testModule\v12_caseA_open_loop_tec_1000s_from727K500s --case-prefix v12_caseA_open_loop_tec_1000s_from727K500s
```

The final TEC-coupled record at absolute time `1500.0 s` was approximately `Tin=727.000 K`, `Tcore_out=821.148 K`, `Tradiator_in=821.157 K`, `Tradiator_out_mix=725.128 K`, `coolant_enthalpy_rise=106.749 kW`, `Qrad=108.618 kW`, `Pel=4.787 kW`, `I=176.004 A`, `core_delta_p=3856.29 Pa`, `radiator_tube_total_flow=1.299170 kg/s`, `tube_mean=0.0166560 kg/s`, `tube_min=0.0159289 kg/s`, and `tube_max=0.0187753 kg/s`.

To match V11, the TEC wire resistance scale should be `0.5`. Continuing from the V12 TEC `1500 s` checkpoint with `--wire-resistance-scale 0.5` for `200 s` completed:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v12_caseA_open_loop.py --restart-in testModule\v12_caseA_open_loop_tec_1000s_from727K500s\v12_caseA_open_loop_tec_1000s_from727K500s_latest_restart.npz --inlet-temperature-k 727 --enable-tec-coupled --thermo-update-interval 0.5 --wire-resistance-scale 0.5 --duration 200 --record-interval 25 --restart-interval 50 --max-dt 0.5 --output-dir testModule\v12_caseA_open_loop_tec_wire05_200s_from1500s --case-prefix v12_caseA_open_loop_tec_wire05_200s_from1500s
```

The final record at absolute time `1700.0 s` was approximately `Tin=727.000 K`, `Tcore_out=820.471 K`, `Tradiator_in=820.471 K`, `Tradiator_out_mix=724.995 K`, `coolant_enthalpy_rise=105.982 kW`, `Qrad=108.522 kW`, `Pel=5.416 kW`, `I=199.120 A`, `core_delta_p=3857.19 Pa`, `radiator_tube_total_flow=1.300180 kg/s`, `tube_mean=0.0166690 kg/s`, `tube_min=0.0159419 kg/s`, and `tube_max=0.0187900 kg/s`.

## 23. 2026-06-18 V13 CaseA closed core + TOPAZ-II pipe-fin radiator

V13 is implemented in `testModule/test_core_assemble_v13_caseA.py` and
`testModule/run_v13_caseA_closed_loop.py`. It converts the V12 open-loop
pipe-fin radiator model into a closed pumped loop by removing the fixed
temperature/pressure external boundaries and using `V12_CoreInletConnector`
as the passive pressure reference.

The default TEC update interval in the V13 runner is `1.0 s`. This reduces
the number of expensive `ThermoCalcModel.calculate()` calls compared with the
previous exploratory `0.5 s` setting. The TOPAZ-II upper/lower ring headers in
V13 are currently hydraulic volumes only; they do not have separate solid
wall radiation boundaries, so radiator thermal tuning is exposed through
`--tube-emissivity` and `--fin-emissivity` rather than a header emissivity.
`testModule/run_v13_caseA_eps088_wait_then_to10000.py` is a local helper that
waits for the active `eps=0.88/0.88` long run to finish and then continues the
same case to about `10000 s` absolute time.

## 24. 2026-06-25 V11/V13 optional reserved parallel TEC circuit

V11 and V13 now support one optional reserved parallel TEC circuit without changing the default case. The default remains:

```text
thermal ring multipliers:  Center=1, Ring1=6, Ring2=12, Ring3_TEC=15, Ring3_Open=3
main TEC multipliers:      Center=1, Ring1=6, Ring2=12, Ring3_TEC=15, Ring3_Open=0
```

So the main circuit is still the existing 34-TEC fixed-voltage series circuit. `Ring3_Open` stays disconnected and the runner keeps the zero-source check unless the reserved circuit is explicitly enabled.

The optional circuit uses only the three `Ring3_Open` TECs and is configured as a second independent `ThermoCalcModel`:

```powershell
# Fixed bus voltage for the reserved parallel circuit
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_caseA_closed_loop.py --duration 0.1 --record-interval 0.1 --restart-interval 0.1 --max-dt 0.1 --enable-reserved-parallel-tec --reserved-parallel-mode fixed_u --reserved-parallel-voltage 0.8

# Fixed total current for the reserved parallel circuit
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_caseA_closed_loop.py --duration 0.1 --record-interval 0.1 --restart-interval 0.1 --max-dt 0.1 --enable-reserved-parallel-tec --reserved-parallel-mode fixed_i --reserved-parallel-current 1000.0

# External load curve, CSV or NPZ, interpreted as U_load=f(I_total)
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_caseA_closed_loop.py --duration 0.1 --record-interval 0.1 --restart-interval 0.1 --max-dt 0.1 --enable-reserved-parallel-tec --reserved-parallel-mode load_curve --reserved-parallel-load-curve path\to\curve.csv
```

The same flags are available in `run_v11_caseA_closed_loop.py`. For `load_curve`, CSV may contain either `current_a,voltage_v` or `current,voltage` headers; `.npz` may contain `current_a/voltage_v` or `current/voltage`.

Diagnostics:

- `tec_total_voltage_v` and `tec_total_current_a` remain the main series-circuit terminal voltage/current for backward compatibility.
- `tec_total_electric_power_w` is the sum of main series power and reserved parallel power.
- `tec_main_*` records the main series circuit.
- `tec_reserved_parallel_*` records the optional `Ring3_Open` parallel circuit and is `null` when disabled.
- `reserved_parallel_tec_enabled` and `reserved_parallel_tec_mode` are written to `latest_state.json` and CSV history records.

Validation completed on 2026-06-25:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile Components\ReactorCore.py testModule\run_v8_caseA_common.py testModule\test_core_assemble_v7_caseA.py testModule\run_v11_caseA_closed_loop.py testModule\run_v13_caseA_closed_loop.py testModule\test_reactorcore_tec_topology.py ThermoCalc\ThermoCalcWrapper.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_reactorcore_tec_topology.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_parallel.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_interface.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_lookup.py
```

V13 `0.1 s` smoke results:

| Case | Result |
| --- | --- |
| Default, reserved circuit disabled | completed; `tec_reserved_parallel_* = null`, `Pel ~= 5403.650 W` |
| Reserved `fixed_u`, `0.8 V` | completed; reserved branch `I ~= 1021.511 A`, `P ~= 817.209 W`, total `Pel ~= 6220.859 W` |
| Reserved `fixed_i`, `1000 A` | completed; reserved branch `U ~= 0.812 V`, `I ~= 1000.148 A`, total `Pel ~= 6215.914 W` |
| Reserved `load_curve` with `U=0.0008 I` | completed; reserved branch `U ~= 0.807 V`, `I ~= 1009.013 A`, total `Pel ~= 6218.004 W` |

V11 `--create-init-only` completed both with the default disabled reserved circuit and with `--enable-reserved-parallel-tec --reserved-parallel-mode fixed_u --reserved-parallel-voltage 0.8`.

## 25. 2026-06-25 V13 first-version radiator thermal shield

V13 now has an optional first-version quasi-steady thermal shield for the TOPAZ-II pipe-fin radiator. It is default-off and only affects the V13 `RadiatorPipeWithFin` path. The model does not add a solid ODE and does not use the future `8 x 12` view-factor matrix; it updates each radiator unit's equivalent radiation background temperature before radiator `pre_step()`.

Runner flags:

```powershell
--enable-radiation-shield
--shield-active-until-s 10
--shield-inner-emissivity 0.8
--shield-outer-emissivity 0.8
--shield-conductivity-w-m-k 1.0
--shield-thickness-m 0.002
--shield-view-factor 0.8
--shield-solar-heat-flux-w-m2 0.0
--shield-background-temperature-k 3.0
--shield-relaxation 1.0
```

`--shield-active-until-s` is interpreted as relative time from the current run start; the runner stores the absolute cutoff in `radiation_shield_active_until_abs_s`.

New diagnostics are included in V13 JSON/CSV:

- `radiation_shield_enabled`
- `radiation_shield_active`
- `radiation_shield_effective_background_mean_k`
- `radiation_shield_inner_temperature_mean_k`
- `radiation_shield_outer_temperature_mean_k`
- `radiation_shield_q_from_radiator_w`
- `radiation_shield_q_to_space_w`
- `radiation_shield_view_factor`
- `radiation_shield_conductivity_w_m_k`
- `radiation_shield_thickness_m`

Validation completed:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_radiator_thermal_shield.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile Components\RadiatorThermalShield.py Components\RadiatorPipeWithFin.py testModule\test_core_assemble_v13_caseA.py testModule\run_v13_caseA_closed_loop.py testModule\test_radiator_thermal_shield.py
```

V13 `0.1 s` smoke from the current V12 warm restart:

| Case | Result |
| --- | --- |
| Shield off | completed; `q_radiator_total_w ~= 105380.709 W` |
| Shield on, `view_factor=0.8`, `thickness=0.002 m`, `k=1.0 W/m/K` | completed; `q_radiator_total_w ~= 93284.738 W`, effective background `~= 441.039 K`, outer shield mean `~= 424.709 K`, `q_to_space ~= 10641.871 W` |

This first version is intended to prove the boundary-coupling workflow. A later version should replace the equivalent background model with axial aggregation to 12 radiator segments and an 8-shield-segment view-factor matrix.

## 26. 2026-06-27 V13-start cold startup runner

`testModule/run_v13_start_case.py` is the first V13-start cold-start runner built from the V13 closed-loop pipe-fin radiator case. It is intentionally separate from `run_v13_caseA_closed_loop.py` and does not load V10/V11/V13 warm restarts by default.

Main startup assumptions:

- Engineering profile starts from `373 K`; `--startup-profile titam` switches the initial fluid/solid temperature to `300 K` and uses `1.5 kg/s` target flow.
- The runner explicitly cold-initializes all registered solids through `reset_solid_temperatures(...)`; V13 base builders otherwise keep several `ReactorCore/TFEUnit` default solid temperatures around `743 K`, `850 K`, `1200 K`, or `1550 K`, which is not valid for a cold-start calculation.
- The simplified startup controller follows the cold-start document milestones: safety drums withdraw over `8 s`, control drums advance at `1.4 deg/s`, first critical is represented near `125 deg`, the low-power hold is `5 kW`, then ramps to `35 kW` at `600 W/s` and to `110 kW` at `80 W/s`.
- TEC coupling is deferred. The runner does not build/apply ThermoCalc circuits during cold hydraulic initialization; it configures TEC, applies the V11/V13 half wire-resistance scale, and enables `core.enable_tec_coupled` only after both the post-critical time and emitter-temperature thresholds are satisfied. The default ThermoCalc update interval is `0.5 s`.
- The emitter-collector vapor/gas gap is simplified by changing the TEC gap coupler's equivalent gas conductance from a helium value to a cesium value with a smooth transition.

Radiator startup boundary handling:

- The V13 radiator thermal shield is attached by default with the `fortran_shield2` model and remains active until the startup controller latches shield jettison at `core_inlet_temperature_k >= 400 K`.
- The 78-tube external heat matrix in `Components/ExternalHeatSources` is applied to the `RadiatorPipeWithFin` wall outer boundary via one embedded matrix column per tube.
- The 6-partition shield external heat matrix is sampled each step and mapped to the first six entries of the shield `qsss_w_m2` vector; the upper/lower entries are currently set to zero. This is a first modeling assumption and should be revisited if the source table gains explicit upper/lower shield partitions.

Useful command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_start_case.py --duration 10 --record-interval 1 --restart-interval 5 --max-dt 0.1 --output-dir testModule\v13_start_smoke_10s_stepiter300_20260627 --case-prefix v13_start_smoke_10s_stepiter300_20260627
```

Validation completed on 2026-06-27:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile Solvers\SystemManager.py testModule\v13_startup_control.py testModule\test_v13_startup_control.py testModule\run_v13_start_case.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_v13_startup_control.py
```

Current 10 s smoke result:

- Completed from a `373 K` cold solid/fluid initialization; `cold_initialized_solid_count = 114`.
- Mean emitter temperature started near `373.2 K`, confirming that the previous hot-solid initialization problem is fixed.
- Core inlet temperature exceeded the `400 K` shield-jettison threshold within the first second when the full external heat matrices were applied, so the shield became inactive early. This is physically controlled by the current external-heat scale/area and low-temperature loop thermal inventory, not by a hard-coded time cutoff.
- Several fluid steps still reported small Picard residual plateaus even with `--step-hydraulic-max-iter 300`; later calibration should check external-heat magnitude/area mapping and time-step limits before treating this as a topology failure.

### 26.1 2026-06-28 V13-start W/m external heat and full cold initialization update

User clarified that the 78-tube startup external-heating table is a line load in `W/m`, not a heat-flux density in `W/m2`. `testModule/v13_startup_control.py` now provides `MatrixColumnLineHeatSource`, which converts one matrix column as:

```text
Q_node [W] = q_line(t) [W/m] * node_length [m]
q_equiv [W/m2] = Q_node / area_used_by_ExternalHeatFluxBC
```

This preserves the intended per-node heat power after `ExternalHeatFluxBC` multiplies by boundary area. `run_v13_start_case.py` defaults to `--tube-external-heat-units W/m`; `W/m2` remains available for compatibility.

The previous V13-start cold-smoke still showed fast inlet heating even with external heat disabled. The root cause was incomplete cold initialization: solids were reset, but the hydraulic control volumes and `HydraulicNetwork.T_vec/h_vec` could retain warm V12/V13 defaults. `reset_fluid_temperatures(...)` now resets every hydraulic control volume `T/h`, clears transient fluid heat sources, and resynchronizes the network vectors before `initialize_system()`.

Validation after this fix:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile Solvers\SystemManager.py testModule\v13_startup_control.py testModule\test_v13_startup_control.py testModule\run_v13_start_case.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_v13_startup_control.py
```

Short-run results:

- `external_heat_scale=0`, `10 s`: inlet stayed near `373 K`, shield stayed active, proving full fluid/solid cold initialization.
- default `W/m` external heat, `10 s`: inlet reached about `402 K` near `8 s`, shield jettisoned by the temperature trigger, mean emitter remained near `375 K`.
- default `W/m` external heat, `120 s`, `max_dt=0.5 s`: reached `REACTIVITY_PULLBACK`, prescribed thermal power about `3.79 kW`, inlet about `511.6 K`, emitter about `480.1 K`. This larger time step produced sustained hydraulic residual warnings after about `57 s`, so the pre-TFE long run uses `max_dt=0.1 s`.

Current pre-TFE long run launched in the background on 2026-06-28:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_start_case.py --duration 1500 --record-interval 10 --restart-interval 100 --max-dt 0.1 --step-hydraulic-max-iter 300 --output-dir testModule\v13_start_pre_tfe_1500s_20260628 --case-prefix v13_start_pre_tfe_1500s_20260628
```

This run intentionally stops before the default TEC enable time (`critical_time + 1500 s`, about `1597 s`) and is intended to cover cold startup through high-power hold before TFE ignition.

### 26.2 2026-06-28 pre-TFE h-gap iteration

The default `helium_gap_h_eq_w_m2_k=1200` full cold-to-1500 s run completed, but the pre-TFE state was cooler than the emitter ignition target:

```text
t = 1500 s, phase = CRITICAL_POWER_HOLD, Q = 110 kW
core inlet ~= 826.6 K, core outlet ~= 918.5 K, mean emitter ~= 984.2 K
TEC disabled, shield inactive
```

A 1500 s restart continuation to 1590 s with the same `h_eq=1200` cooled further to mean emitter about `970.7 K`. This indicates the case had already passed a thermal peak and would not approach the nominal `1050 K` emitter ignition window by simply holding the same settings.

Short restart-based sensitivity from the 1500 s state showed the expected trend: reducing the pre-ignition He-side emitter-collector gap conductance raises emitter temperature and lowers coolant temperature.

| He gap h_eq [W/m2/K] | 1590 s mean emitter [K] | 1590 s inlet [K] | Note |
| --- | ---: | ---: | --- |
| 800 | ~1002.1 | ~789.3 | Warmer emitter, still below ignition window |
| 600 | ~1031.0 | ~787.1 | Near pre-ignition target without crossing 1050 K |
| 400 | ~1081.6 | ~782.8 | Crosses ignition-window temperature |

Based on this, a full cold-start `h_eq=600` pre-TFE run was launched in the background:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_start_case.py --duration 1590 --record-interval 10 --restart-interval 100 --max-dt 0.1 --step-hydraulic-max-iter 300 --helium-gap-h-eq-w-m2-k 600 --tfe-start-after-critical-s 100000 --output-dir testModule\v13_start_pre_tfe_h600_1590s_20260628 --case-prefix v13_start_pre_tfe_h600_1590s_20260628
```

The high `--tfe-start-after-critical-s` intentionally keeps TEC disabled so this remains a pre-TFE ignition thermal-hydraulic run.

### 26.3 2026-06-28 V13-start restart continuation and active h600 run

`run_v13_start_case.py` now supports `--restart-in` for V13-start restart continuation. The parser uses `allow_abbrev=False` so `--restart-in` cannot be accidentally interpreted as `--restart-interval`. Restart continuation rebuilds the same V13-start topology, attaches the radiator shield and external-heat boundary modifiers, loads `SystemManager.load_global_state(...)`, restores V13 design flows and pump head, and then lets the absolute-time startup controller select the current phase.

The `h_eq=600` full cold-start pre-TFE run is active in the background:

```text
output: testModule/v13_start_pre_tfe_h600_1590s_20260628
pid:    67352
```

Early monitor result:

```text
t = 10 s: Tin ~= 412.6 K, mean emitter ~= 376.4 K, TEC disabled, gap h_eq = 600 W/m2/K
t = 30 s: Tin ~= 453.9 K, mean emitter ~= 402.9 K, TEC disabled
```

The run is progressing, not hung, but the early hydraulic residual plateau makes it slow: a 20 s wall-clock monitor showed about 19.6 s CPU progress and the log advanced from `t=31.7 s` to `t=34.2 s`. The current command intentionally keeps `--step-hydraulic-max-iter 300`; if this proves too slow, a later repeat can reduce the per-step limit after confirming that the residual plateau does not materially affect thermal trends.

### 26.4 2026-06-28 h600 run speed mitigation

The full cold `h_eq=600`, `step_hydraulic_max_iter=300` run progressed but was too slow in the early hydraulic residual plateau. A 60 s cold-start comparison with `--step-hydraulic-max-iter 50` produced the same 10/20/30/40 s temperatures as the interrupted 300-iteration run and completed the cold-start segment to 60 s:

```text
t = 60 s: Tin ~= 489.4 K, mean emitter ~= 443.4 K, TEC disabled, h_eq = 600 W/m2/K
```

That 60 s state is now used as the restart base for the remaining pre-TFE continuation to 1590 s:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_start_case.py --restart-in testModule\v13_start_h600_iter50_60s_20260628\v13_start_h600_iter50_60s_20260628_latest_restart.npz --duration 1530 --record-interval 10 --restart-interval 100 --max-dt 0.1 --step-hydraulic-max-iter 50 --helium-gap-h-eq-w-m2-k 600 --tfe-start-after-critical-s 100000 --output-dir testModule\v13_start_pre_tfe_h600_60to1590s_iter50_20260628 --case-prefix v13_start_pre_tfe_h600_60to1590s_iter50_20260628
```

Active background process:

```text
pid: 63684
```

Initial continuation monitor:

```text
absolute t ~= 70 s, relative continuation time ~= 10 s
Tin ~= 497.3 K, mean emitter ~= 454.1 K, h_eq = 600 W/m2/K, TEC disabled
```

### 26.5 2026-06-28 active h600 continuation monitor

The `v13_start_pre_tfe_h600_60to1590s_iter50_20260628` continuation remains active and is progressing normally. Monitor snapshot:

```text
pid: 63684
absolute t ~= 140 s, relative continuation time ~= 80 s
phase = FAST_POWER_RAMP
Q ~= 12.63 kW
Tin ~= 533.4 K
mean emitter ~= 515.8 K
TEC disabled
```

A 60 s wall-clock monitor showed about 57.9 s CPU progress and advanced the case from the pullback/low-power region into the fast power ramp. No `run.err` messages have appeared for this continuation so far.

### 26.6 2026-06-28 completed segmented h600 pre-TFE startup chain

The pure-cold to pre-TFE `h_eq=600 W/m2/K` startup calculation was completed as a segmented restart chain. This is a thermal-hydraulic startup demonstration driven by the engineering startup power schedule; TEC remains disabled by `--tfe-start-after-critical-s 100000`.

Segment chain:

| Segment | Output directory | Time span [s] | max_dt [s] | Hydraulic max iter | End state |
| --- | --- | ---: | ---: | ---: | --- |
| Cold initialization | `testModule/v13_start_h600_iter50_60s_20260628` | 0 -> 60 | 0.1 | 50 | `Tin ~= 489.4 K`, emitter `~= 443.4 K` |
| Early continuation | `testModule/v13_start_pre_tfe_h600_60to1590s_iter50_20260628` | 60 -> 260 | 0.1 | 50 | `Tin ~= 594.1 K`, emitter `~= 618.4 K` |
| Larger-step check | `testModule/v13_start_h600_260to460_dt05_test_20260628` | 260 -> 460 | 0.5 | 50 | `Tin ~= 698.1 K`, emitter `~= 814.3 K` |
| Final pre-TFE continuation | `testModule/v13_start_pre_tfe_h600_460to1590s_dt05_20260628` | 460 -> 1590 | 0.5 | 50 | `Tin ~= 807.5 K`, emitter `~= 1131.2 K` |

Final state at absolute `t ~= 1590 s`:

```text
phase = CRITICAL_POWER_HOLD
core thermal power = 110000 W
TEC coupled = false
shield jettisoned = true
core inlet = 807.517 K
core outlet = 913.060 K
core delta-T = 105.543 K
mean emitter = 1131.202 K
closed-loop flow = 1.29999996 kg/s
radiator tube external heat = 11568.035 W
radiator heat rejection total = 144715.930 W
78-tube flow mean/min/max = 0.0166656 / 0.0159366 / 0.0187932 kg/s
```

The final segment `run.err` contains one startup-step hydraulic residual warning at `t=460 s`:

```text
Hydraulic step NOT converged after 50 iterations. Max residual: 3.760874e-04
Fluid solver NOT converged at t=460.0000s
```

The warning did not stop the run and no later fatal error was recorded. Because the final emitter temperature is above the nominal `1050 K` ignition threshold while TEC is intentionally disabled, this chain is best treated as a successful pre-TFE thermal-hydraulic startup reachability run, not as the final calibrated ignition-boundary case. A later TFE ignition calculation should restart from this state or from a retuned endpoint, then enable TEC and cesium-gap behavior explicitly.

### 26.7 2026-06-28 zero-power shield-only cooling check and external-heat unit correction

The embedded `is58p5_w0_8p12_N78_sum` tube matrix should be treated as `W/m2`, not `W/m`. The earlier `W/m` pre-TFE runs are therefore diagnostic only and should not be used as a physical startup baseline. `run_v13_start_case.py` now defaults `--tube-external-heat-units W/m2`; the legacy `W/m` conversion remains available only for explicit sensitivity checks.

A zero-core-power cooling case was run to isolate the radiator/shield thermal behavior without startup power ramp:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_start_case.py --duration 500 --record-interval 10 --restart-interval 100 --max-dt 0.5 --step-hydraulic-max-iter 50 --fixed-startup-power-w 0 --source-power-w 0 --shield-jettison-temperature-k 1000000000 --tube-external-heat-units W/m2 --tube-external-heat-area-fraction 0 --output-dir testModule\v13_start_zero_power_shield_only_500s_20260628 --case-prefix v13_start_zero_power_shield_only_500s_20260628
```

This case keeps the shield active, applies the N6 shield external heat table, and disables direct N78 tube-wall heating by setting `--tube-external-heat-area-fraction 0`. Results show no reverse heating of the coolant:

| t [s] | Tin [K] | Tout [K] | Radiator out [K] | Mean emitter [K] | q_rad total [W] | Shield effective background [K] | Tube external heat [W] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 372.704 | 372.608 | 372.564 | 372.953 | 500.827 | 364.842 | 0.0 |
| 60 | 371.853 | 372.193 | 371.819 | 372.459 | 612.455 | 362.402 | 0.0 |
| 120 | 371.061 | 371.478 | 371.032 | 371.774 | 679.765 | 360.446 | 0.0 |
| 240 | 369.221 | 369.674 | 369.222 | 369.945 | 603.667 | 359.712 | 0.0 |
| 360 | 368.627 | 368.937 | 368.612 | 369.016 | 466.613 | 361.331 | 0.0 |
| 500 | 367.352 | 367.793 | 367.271 | 367.973 | 1013.908 | 350.567 | 0.0 |

Interpretation: with direct tube-wall external heat disabled and the shield retained, the coolant cools slowly from `373 K` to about `367.35 K` over `500 s`. The shield background stays around `350-365 K`, so it reduces net radiative loss but does not heat the loop above its initial temperature. This supports treating previous rapid warmup as a boundary-application error, not as a deep-space radiation-temperature error.

### 26.8 2026-06-28 corrected shield-startup run after external-heat fix

After correcting the startup boundary interpretation, the startup run was repeated with direct N78 tube-wall external heat disabled and the external heat table unit treated as `W/m2`. The shield still receives the N6 external heat table and is jettisoned by the inlet-temperature trigger.

Command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_start_case.py --duration 1590 --record-interval 20 --restart-interval 100 --max-dt 0.5 --step-hydraulic-max-iter 50 --helium-gap-h-eq-w-m2-k 600 --tfe-start-after-critical-s 100000 --tube-external-heat-units W/m2 --tube-external-heat-area-fraction 0 --shield-jettison-temperature-k 400 --output-dir testModule\v13_start_corrected_shield_startup_1590s_20260628 --case-prefix v13_start_corrected_shield_startup_1590s_20260628
```

Representative records:

| t [s] | Phase | Q_core [kW] | Tin [K] | Tout [K] | Mean emitter [K] | Shield active | Shield bg [K] | q_rad total [kW] | Tube external heat [W] |
| ---: | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 20 | CONTROL_DRUM_APPROACH | 0.001 | 372.314 | 372.552 | 372.896 | true | 364.439 | 0.513 | 0.0 |
| 100 | INITIAL_SUPERCRITICAL_RAMP | 0.453 | 371.246 | 371.661 | 372.103 | true | 359.871 | 0.727 | 0.0 |
| 140 | FAST_POWER_RAMP | 12.629 | 371.259 | 372.559 | 378.566 | true | 354.840 | 1.045 | 0.0 |
| 180 | SLOW_POWER_RAMP | 35.217 | 376.982 | 382.998 | 411.018 | true | 362.500 | 1.109 | 0.0 |
| 220 | SLOW_POWER_RAMP | 38.417 | 396.645 | 406.274 | 448.429 | true | 384.362 | 1.272 | 0.0 |
| 240 | SLOW_POWER_RAMP | 40.017 | 404.256 | 417.375 | 464.935 | false | 3.000 | 9.045 | 0.0 |
| 500 | SLOW_POWER_RAMP | 60.817 | 524.357 | 551.232 | 642.740 | false | 3.000 | 25.740 | 0.0 |
| 800 | SLOW_POWER_RAMP | 84.817 | 623.017 | 672.355 | 816.146 | false | 3.000 | 51.899 | 0.0 |
| 1000 | SLOW_POWER_RAMP | 100.817 | 671.252 | 736.292 | 912.425 | false | 3.000 | 70.371 | 0.0 |
| 1120 | CRITICAL_POWER_HOLD | 110.000 | 695.466 | 769.981 | 964.763 | false | 3.000 | 81.385 | 0.0 |
| 1300 | CRITICAL_POWER_HOLD | 110.000 | 719.217 | 802.593 | 1007.717 | false | 3.000 | 93.193 | 0.0 |
| 1590 | CRITICAL_POWER_HOLD | 110.000 | 733.807 | 823.482 | 1031.198 | false | 3.000 | 101.128 | 0.0 |

Final state at `1590 s`:

```text
phase = CRITICAL_POWER_HOLD
core power = 110 kW
TEC coupled = false
shield jettisoned = true
core inlet/outlet = 733.807 / 823.482 K
mean emitter = 1031.198 K
radiator heat rejection = 101.128 kW
coolant enthalpy rise = 101.661 kW
closed-loop flow = 1.300000 kg/s
78-tube flow mean/min/max = 0.0166667 / 0.0159400 / 0.0187875 kg/s
```

Compared with the earlier incorrect `W/m` + direct tube heating chain, the corrected run is much cooler and no longer shows artificial external warmup. At `1590 s`, the mean emitter is close to but still below the nominal `1050 K` TFE ignition threshold. This endpoint is a better pre-TFE restart candidate than the previous `1131 K` endpoint, but a real ignition run may need additional high-power hold time, a slightly different gap conductance, or an explicit ignition criterion sweep.

Caveat: with `--step-hydraulic-max-iter 50`, the run produced repeated hydraulic residual warnings in the mid-transient, mainly after about `294 s`. The run completed and the thermal trajectory is smooth, but a final production startup baseline should repeat a narrower segment with stricter hydraulic convergence or a smaller `max_dt`.

### 26.9 2026-06-28 ThermoCalc low-temperature zero-emission guard

A ThermoCalc regression was added for low-temperature fixed-voltage TEC startup states. The reproduced bad case was one TEC with `TE=800 K`, `TC=600 K`, `Tcs=520 K`, `fixed_u=0.8 V`, and 37 axial nodes; before the fix it entered `ThermoCalcModel.calculate()` and did not return within the test timeout.

`ThermoCalc/ThermoCalcWrapper.py` now skips C++ circuit iteration when all emitter temperatures are below the default zero-emission cutoff `1000 K`. The skipped result is an open-circuit zero-current state with zero TEC heat terms and diagnostic fields `zero_emission_skipped` / `zero_emission_reason` in `get_global_results()`.

Relevant verification:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_interface.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_parallel.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_lookup.py
```

### 26.10 2026-06-28 V13-start cesium conditioning before TEC electrical startup

The V13-start TEC ignition path now treats cesium conditioning as a prerequisite for electrical generation. This avoids forcing the external circuit solver into a low-temperature, low-emission operating point that is not a meaningful power-generation state.

Current startup sequence:

1. Helium-gap thermal startup with TEC electrical calculation disabled.
2. Cesium conditioning with TEC electrical calculation still disabled. The controller ramps `startup_cs_fraction` toward `1` and lowers the TEC gap equivalent conductance toward the V13 cesium-side value.
3. Electrical TEC startup only after cesium and emitter-temperature gates are satisfied. The main TEC circuit starts in fixed-resistance mode and switches to fixed total voltage only if the finite main voltage reaches the configured switch voltage.

Relevant runner controls in `testModule/run_v13_start_case.py`:

- `--tec-electrical-start-after-cesium-s`
- `--tec-electrical-start-cs-fraction`
- `--tec-electrical-start-emitter-temperature-k`
- `--startup-main-tec-initial-mode fixed_r|fixed_u`
- `--startup-main-tec-load-resistance-ohm`
- `--startup-main-tec-load-resistance-scope total|per_tec`
- `--startup-main-tec-switch-voltage-v`
- `--initial-cs-fraction`

The startup runner now records TEC solver diagnostics in history and `latest_state.json`: `tec_solver_mode`, `tec_solver_converged`, `tec_solver_iteration_count`, `tec_solver_zero_emission_skipped`, `tec_solver_zero_emission_reason`, and `tec_solver_output_finite`.

`Components/ReactorCore.py` also synchronizes TEC group temperatures immediately before `ThermoCalc.calculate()`. This is required for restart-based ignition runs: without the explicit sync, the wrapper could see stale/default low emitter temperatures after the thermal model had already reached the intended startup state.

Stable cesium-conditioned, pre-electrical restart:

```text
testModule/v13_start_cesium_conditioning_plus1000s_20260628/v13_start_cesium_conditioning_plus1000s_20260628_latest_restart.npz
```

Final state of that pre-electrical segment:

```text
t = 4090 s
startup_cs_fraction = 0.9999983895
startup_tec_gap_h_eq_w_m2_k = 250.0006
mean emitter ~= 1198.96 K
core inlet/outlet ~= 746.56 / 842.32 K
TEC electrical calculation disabled
```

The stable fixed-resistance electrical continuation used ThermoCalc lookup tables and interpreted `0.0044 ohm` as a per-TEC load for the 34-TEC main series chain:

```text
R_total = 34 * 0.0044 = 0.1496 ohm
```

Stable fixed-R restart:

```text
testModule/v13_start_tec_fixedr_totalR_lookup_500s_20260628/v13_start_tec_fixedr_totalR_lookup_500s_20260628_latest_restart.npz
```

Final state of that fixed-R segment:

```text
t = 4612 s
Rload = 0.1496 ohm
tec_main_voltage_v ~= 9.70466 V
tec_main_current_a ~= 64.8707 A
tec_main_electric_power_w ~= 629.55 W
mean emitter ~= 1188.26 K
core inlet/outlet ~= 745.85 / 841.21 K
fixed-voltage switch not triggered
```

A fixed-R load scan from this thermal state did not reach the target `27.2 V`; even high-load/open-circuit-direction cases remained near `20 V`. Direct fixed-voltage `27.2 V` is not yet a reliable continuation point with the current local backend because a lookup-enabled diagnostic returned non-finite voltage output with `iteration_count=100`.

Therefore the current continuation gate is:

- If the ThermoCalc backend is updated to make fixed-U `27.2 V` finite and converged, continue from the stable fixed-R restart above.
- Otherwise, keep TEC electrical startup gated until the thermal state or physical assumptions are adjusted enough for the fixed-R ramp to approach the target voltage.

Detailed run notes and the load-scan table are in `testModule/V13_STARTUP_TEC_IGNITION_STATUS_20260628.md`.

### 26.10 2026-06-28 ThermoCalc C++ no-hang iteration fallback

The ThermoCalc low-temperature regression now has two checks:

- default path: Python zero-emission guard skips a known non-generating low-temperature fixed-voltage case;
- diagnostic path: with `THERMOCALC_DISABLE_ZERO_EMISSION_GUARD=1`, the same case still returns through the rebuilt C++ circuit solver instead of hanging.

`testModule/test_thermocalc_interface.py` covers both cases. The second case verifies that users can keep TEC enabled while exploring uncertain startup states: if the state cannot produce a meaningful electrical solution, ThermoCalc returns finite diagnostics instead of blocking the whole transient run.

### 26.11 2026-06-28 updated ThermoCalc backend fixed-U smoke and TEC temperature gate

After the C++ iteration guard update, `ThermoCalc/te_solver.cp312-win_amd64.pyd` was rebuilt locally at `2026-06-28 23:42:33` (`396288` bytes). A `0.5 s` V13-start smoke from the stable fixed-R restart completed without hanging:

```text
output: testModule/v13_start_tec_fixedu27p2_smoke_0p5s_20260628
restart_in: testModule/v13_start_tec_fixedr_totalR_lookup_500s_20260628/v13_start_tec_fixedr_totalR_lookup_500s_20260628_latest_restart.npz
tec_solver_mode = fixed_u
tec_solver_output_finite = True
tec_solver_converged = False
tec_solver_iteration_count = 45
tec_main_voltage_v = 27.2 V
tec_main_current_a = 0.0 A
tec_main_electric_power_w = 0.0 W
mean_emitter_temperature_k ~= 1224.24 K
```

This proves the local backend now returns a finite non-converged state instead of entering the previous dead loop, but it does not prove that fixed-U `27.2 V` is physically usable at the current startup thermal state.

A no-time-advance ThermoCalc sensitivity from the same restart estimated the temperature margin to the fixed-R switch point. With the selected startup resistance `R_total = 0.1496 ohm`, the fixed-R voltage rises from about `9.84 V` at the current synced temperature field to about `27.66 V` only after adding a synthetic `+400 K` emitter offset. In mean-temperature terms this corresponds to roughly `1.58e3 K` for this diagnostic. Higher startup resistance reaches the voltage threshold earlier (`R=1 ohm` crosses near a `+250 K` offset; high-load/open-direction `R=100 ohm` crosses near `+100 K`), but that is a modeling/electrical-startup choice.

Therefore the current production continuation should not switch to fixed total voltage yet. Continue fixed-R only if the thermal model can plausibly raise the emitter field toward the switch condition, or explicitly revise the startup load strategy before launching a long steady-state run.

### 26.12 2026-06-29 R=100 ohm fixed-R startup continuation check

A higher startup resistance was tested after the ThermoCalc no-hang backend update. The `R_total=100 ohm` fixed-R path continued stably for `100 s` from the earlier 20 s R=100 ohm smoke and did not force fixed-U mode.

```text
output: testModule/v13_start_tec_fixedr_R100_20to120s_20260629
t = 4732 s
tec_solver_mode = fixed_r
tec_solver_converged = True
tec_solver_iteration_count = 2
tec_main_voltage_v ~= 23.18796 V
tec_main_current_a ~= 0.23188 A
tec_main_electric_power_w ~= 5.3768 W
mean_emitter_temperature_k ~= 1196.56 K
core inlet/outlet ~= 745.41 / 840.72 K
fixed-voltage switch not triggered
```

The voltage and emitter temperature both increased during this segment, but the endpoint is still below the `27.2 V` fixed-voltage switching threshold. This supports using the updated backend for continued fixed-R startup exploration, while still keeping the fixed-U transition gated by actual voltage rather than forcing it at the current thermal state.

### 26.13 2026-06-29 R=100 ohm fixed-R plateau below rated voltage

The `R_total=100 ohm` fixed-resistance startup path was continued for another `500 s` with the `27.2 V` automatic switch gate still enabled.

```text
output: testModule/v13_start_tec_fixedr_R100_120to620s_20260629
time span: 4732 s -> 5232 s
tec_solver_mode = fixed_r
tec_solver_converged = True
tec_solver_iteration_count = 2
tec_main_voltage_v ~= 23.32207 V
tec_main_current_a ~= 0.23322 A
tec_main_electric_power_w ~= 5.4392 W
mean_emitter_temperature_k ~= 1199.03 K
core inlet/outlet ~= 746.65 / 842.47 K
fixed-voltage switch not triggered
```

The voltage plateaued around `23.32 V`; the final 100 s slope was slightly negative despite a small ongoing mean-emitter temperature increase. Therefore, under the current `110 kW` startup thermal state and `R=100 ohm` startup load, continuing fixed-R alone is unlikely to reach the `27.2 V` switch gate. The remaining limitation is now physical/modeling, not the previous ThermoCalc dead loop.

### 26.14 2026-06-29 150 kW sensitivity reaches fixed-R voltage gate but fixed-U stalls

A diagnostic 150 kW sensitivity from the `R=100 ohm` plateau state showed that the rated-voltage gate is thermally reachable when the emitter field is hotter:

```text
output: testModule/v13_start_tec_fixedr_R100_power150k_200s_20260629
fixed-startup-power-w = 150000
Rload = 100 ohm
20 s: U ~= 26.656 V, mean emitter ~= 1230.20 K
40 s: U ~= 28.301 V, mean emitter ~= 1255.42 K
```

The `40 s` fixed-R record exceeded the `27.2 V` automatic switch threshold. The process then stopped writing records and consumed CPU until it was manually stopped, indicating trouble after the switch point.

To isolate the issue, the same case was repeated with the switch threshold set to `999 V`, producing a high-temperature fixed-R restart at `t=5272 s`. From this restart:

```text
fixed-R 0.5 s smoke: completed, U ~= 28.710 V, I ~= 0.2871 A, P ~= 8.24 W, converged=True, iteration_count=2
fixed-U 27.2 V 0.5 s smoke: did not return within 600 s and was manually stopped
```

Current implication: the startup fixed-R gate can be crossed by raising the thermal state, but the fixed-U circuit solve at that switch point is still not usable for long transient continuation. The ThermoCalc low-temperature dead-loop guard is not sufficient for this high-temperature fixed-U case; further fixed-U solver repair or a more gradual voltage-control transition is required before launching a final steady-state run.

### 26.15 2026-06-29 fixed-U switch-point stall root cause: lookup coverage gap

Systematic debugging narrowed the high-temperature fixed-U stall to ThermoCalc lookup coverage and repeated full series circuit solves. At the `150 kW / R=100 ohm` switch-point restart (`t=5272 s`), a fixed-R ThermoCalc solve completed but took about `23.6 s` even with lookup enabled:

```text
mode = fixed_r
Uout ~= 28.70997 V
Iout ~= 0.28710 A
converged = True
iteration_count = 2
lookup_found = 1221 / 1258 emission points
lookup_miss = 37 points
miss TE ~= 1769.97-1862.02 K
miss TC ~= 1100.49-1133.69 K
miss Vo ~= 0.88392-0.90588 V
```

The dense runtime v2 table currently leaves this region partly uncovered: `accident` reaches `TE=2400 K` but only `TC=1100 K`; `high_power` reaches higher `TC` only for `TE>=2160 K`; `core` only covers `TC<=900 K`. Thus the hottest collector nodes at the switch point fall back to the analytic thermionic solve. Fixed-R can still return because it converges in two outer iterations, but fixed-U repeatedly calls the full series `circuitCalc(I)` path and can spend many minutes in these fallback-heavy trial states.

Current implication: before a final automatic fixed-R-to-fixed-U steady run, regenerate or add lookup coverage for the switch-point band, at minimum approximately `TE=1700-1900 K`, `TC=1100-1150 K`, `Vo=0.8-1.0 V`, `Tcs=600 K`, or extend the `accident` region to `TC>=1150 K` over the relevant range. Then rerun the fixed-U `27.2 V` smoke from `testModule/v13_start_tec_fixedr_R100_power150k_60s_noswitch_20260629_latest_restart.npz`.

### 26.16 2026-06-29 augmented switch-point lookup and post-guard fixed-U result

After the ThermoCalc C++ no-hang update, a local augmented dense runtime manifest was created at:

```text
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr_augmented/runtime_dense_manifest.json
```

It references the existing `pcs_0p02_5torr` runtime `.tedb` files and adds a small `switch_point` region for the previously missed high-collector-temperature band. This is generated runtime data and should not be committed to git by default.

With environment variables:

```powershell
$env:THERMOCALC_ENABLE_LOOKUP='1'
$env:THERMOCALC_LOOKUP_DB='...\ThermoCalc\emission_runtime_db_v2\pcs_0p02_5torr_augmented'
$env:THERMOCALC_LOOKUP_REGIONS='switch_point,core,startup,high_power,accident'
```

the `150 kW / R=100 ohm` switch-point fixed-R diagnostic loaded 5 dense regions, converged in 2 outer iterations, and covered all emission points:

```text
Uout ~= 28.68281 V
Iout ~= 0.28683 A
ThermoCalc calculate ~= 0.61 s
lookup_found = 1258 / 1258
```

A fixed-U `27.2 V` `0.5 s` V13-start smoke from the same restart then completed rather than hanging, but remained non-converged and expensive:

```text
output = testModule/v13_start_tec_fixedu27p2_power150k_from5272_0p5s_auglookup_20260629
wall time ~= 430.6 s
tec_solver_mode = fixed_u
tec_solver_output_finite = True
tec_solver_converged = False
tec_solver_iteration_count = 100
reported Uout ~= 35.105 V
reported Iout ~= 172.932 A
reported electric power ~= 6070.83 W
```

The old dead loop is therefore fixed, and the switch-point lookup gap is closed, but the current fixed-voltage circuit branch is still not production-ready for automatic V13 startup continuation. It returns finite diagnostics at the iteration cap instead of solving the `27.2 V` operating point.

### 26.17 2026-06-29 fixed-U bracketing repair and hot switch-point lookup extension

The high-temperature `fixed_u` issue after the `R=100 ohm` startup voltage gate was traced to two independent causes:

1. The old `circuitTECs::uFixedCircuitCalc()` used an unbounded secant path and did not require voltage residual convergence before declaring or leaving output state.
2. As the V13 fixed-voltage segment warmed, the local switch-point lookup extension became too narrow. At `t=5314 s`, the miss band had moved to about `TE=1880-1925 K`, `TC=1164-1173 K`, `Vo~=0.918 V`, `Tcs=600 K`.

Source-side fixed-U repair:

```text
ThermoCalc/circuitTECs.cpp::uFixedCircuitCalc()
```

now samples non-negative current guesses, brackets `Utarget - circuitCalc(I)`, solves with a guarded secant/bisection hybrid, and only reports convergence when the voltage residual is satisfied. A small regression assertion was added to `testModule/test_thermocalc_interface.py` so the standard high-temperature fixed-U interface path must converge to its target voltage.

Lookup-side extension:

```text
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr_switch_point_hot/
  switch_point_hot.runtime.v2.npz
  runtime_dense_manifest.json
```

covers `TE=1700-2300 K`, `TC=1080-1300 K`, `Vo=0.70-1.20 V`, `Tcs=580/590/600/610 K` with `74152` safe points. The augmented manifest now loads:

```text
switch_point_hot, switch_point, accident, core, high_power, startup
```

Validation used the rebuilt test pyd via `THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release`; the root production pyd was not overwritten in this step.

Key results:

```text
ThermoCalc interface checks passed.
ThermoCalc lookup checks passed.

Initial high-temperature fixed-U isolated solve:
  Uout = 27.2 V
  Iout ~= 1.2148 A
  converged = True
  iteration_count = 10
  ThermoCalc wall ~= 0.736 s

Automatic fixed-R -> fixed-U 2 s smoke:
  switch at t=5272.5 s, fixed-R U ~= 28.7095 V
  subsequent fixed-U records converged at 27.2 V

Fixed-U 20 s continuation:
  all records converged, iteration_count 7-11
  final U = 27.2 V, I ~= 2.499 A, P ~= 67.98 W

At t=5314 s before hot lookup:
  fixed-U isolated solve still converged but took ~= 28.7 s
  lookup_found = 1254 / 1258

After switch_point_hot:
  lookup_found = 1258 / 1258
  fixed-U isolated solve took ~= 4.25 s

Fixed-U 40 s continuation from t=5314 s:
  all records converged
  final t = 5354 s
  U = 27.2 V, I ~= 4.026 A, P ~= 109.50 W
  mean emitter ~= 1310.85 K
  core inlet/outlet ~= 770.36 / 883.56 K
```

Current status: the fixed-R to fixed-U control chain is now numerically viable over the tested startup segment when using the rebuilt pyd and augmented local lookup data. The case is not yet steady; longer fixed-voltage continuation will likely need either the hot lookup band to remain adequate up to the eventual emitter/collector temperatures or a more complete dense runtime export covering this startup trajectory.

### 2026-06-29 correction: 150 kW runs are diagnostic only

The `power150k` V13 startup directories are not official V13 cold-start physics results. They were created only to force the fixed-R voltage gate above `27.2 V` so the fixed-voltage ThermoCalc branch and lookup coverage could be debugged.

The official V13 cold-start continuation remains startup_thermal_power_w = 110000 W; do not mix the power150k diagnostic trajectories into V13 thermal balance, startup-performance, or V11/V13 comparison conclusions.

Official V13 startup interpretation remains based on the startup-control `110 kW` hold unless the model assumptions are explicitly changed. The latest official `110 kW`, `R_total=100 ohm` fixed-resistance continuation is:

```text
output = testModule/v13_start_tec_fixedr_R100_120to620s_20260629
absolute_time = 5232 s
core/startup thermal power = 110000 W
mode = fixed_r
R_total = 100 ohm
U ~= 23.3221 V
I ~= 0.233221 A
P_e ~= 5.439 W
mean emitter ~= 1199.03 K
core inlet/outlet ~= 746.65 / 842.47 K
fixed-voltage switch = not triggered
last-100-s voltage slope ~= -1.28e-5 V/s
```

Therefore, under the current `110 kW` thermal state and `R_total=100 ohm` startup load, the formal fixed-R path has not reached the rated `27.2 V` switch condition. Any `150 kW` fixed-U continuation should be cited only as a numerical solver/lookup diagnostic, not as the V13 cold-start operating trajectory.

### 2026-06-29 official 110 kW load-ceiling probe

A no-time-advance electrical probe was run from the official `110 kW` restart:

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
absolute_time = 5232 s
startup power = 110000 W
mean emitter from rebuilt system ~= 1199.43 K
TEC gap h_eq ~= 250 W/m2/K
```

The root production pyd probe with a broad scan did not return and had to be stopped, so the two decisive points were repeated with the verified rebuilt test pyd via `THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release` and the augmented lookup regions. This did not overwrite the production pyd.

```text
R_total = 100 ohm:   U ~= 23.2822 V, I ~= 0.232822 A, P ~= 5.421 W, converged=True, iteration_count=2
R_total = 1e6 ohm:  U ~= 23.3653 V, I ~= 2.34e-5 A, P ~= 5.46e-4 W, converged=True, iteration_count=2
```

Interpretation: at the current official `110 kW` thermal state, even the near-open-circuit load voltage is only about `23.37 V`, below the `27.2 V` fixed-voltage switch threshold. Continuing the same fixed-R strategy cannot produce a legitimate switch to rated fixed-voltage generation unless the emitter/collector thermal state or the startup electrical/physical assumptions are changed. The earlier `150 kW` runs remain solver diagnostics only.

### 2026-06-29 official 110 kW R=100 continuation timeout check

A formal continuation was attempted from the official `5232 s` restart without changing the startup power:

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
output = testModule/v13_start_tec_fixedr_R100_620to720s_20260629
requested duration = 100 s
startup power = 110000 W
R_total = 100 ohm
switch voltage = 27.2 V
```

The process timed out after 30 min of wall time and was stopped manually. It had written history records through `t=5312 s` but the latest restart/state had only reached `t=5282 s` because the next restart interval had not completed.

Written history trend:

```text
t=5252 s: U ~= 23.120 V, mean emitter ~= 1198.55 K
t=5272 s: U ~= 23.103 V, mean emitter ~= 1198.20 K
t=5292 s: U ~= 22.981 V, mean emitter ~= 1198.46 K
t=5312 s: U ~= 22.941 V, mean emitter ~= 1198.78 K
```

The run remained in `fixed_r` and did not trigger the `27.2 V` fixed-voltage gate. This reinforces the official-path conclusion: under the current `110 kW`, `R_total=100 ohm` startup state, continuing the same fixed-resistance strategy is moving away from the voltage threshold, not toward it. The timeout also shows that the current production pyd/runtime path is not suitable for further long official startup continuation without first resolving the ThermoCalc performance/hang behavior or explicitly switching to a verified rebuilt backend.

### 2026-06-29 official 110 kW emitter-temperature margin probe

A no-time-advance ThermoCalc probe was run from the official `110 kW` restart to estimate how much hotter the emitter field must be before the startup voltage gate is physically reachable. The probe only changed the ThermoCalc input temperature matrix and did not modify the solid/restart state.

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
absolute_time = 5232 s
startup power = 110000 W
mean TEC emitter input = 1192.77 K
mean TEC collector input = 908.87 K
near-open load = R_total = 1e6 ohm
backend = verified rebuilt pyd via THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release
```

Near-open-circuit voltage versus uniform emitter offset:

```text
dTE =   0 K: U ~= 23.370 V
dTE =  25 K: U ~= 25.235 V
dTE =  50 K: U ~= 27.166 V
dTE =  75 K: U ~= 29.116 V
dTE = 100 K: U ~= 31.017 V
```

Linear interpolation between `25 K` and `50 K` gives a `27.2 V` crossing at approximately `dTE ~= 50.4 K`. Thus the official `110 kW` state is only about `50 K` short in emitter-side voltage margin for a near-open-circuit load. The formal transient, however, was trending slightly cooler/lower-voltage under `R_total=100 ohm`, so reaching the gate by simply continuing the current trajectory is still unlikely.

Interpretation: the limiting issue is primarily the startup thermal state seen by the TECs, not the fixed-voltage switch implementation. Plausible next modeling actions are to heat/hold the emitter field by revisiting startup power schedule, cesium-gap conductance trajectory, radiator/shield heat rejection during ignition, or the electrical load schedule; forcing fixed-U at the current unmodified `110 kW` restart remains unjustified.

### 2026-06-29 cesium-gap h_eq sensitivity for official 110 kW startup

To test whether the `27.2 V` gate is reachable without raising reactor power above the official `110 kW` hold, two thermal-only continuations were run from the official `5232 s` restart. TEC electrical calculation was deliberately disabled by delaying the electrical-start gate, but the cesium gap conductance was changed and applied to the thermal model.

Common setup:

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
startup power = 110000 W
TEC electrical calculation = disabled during thermal hold
initial Cs fraction ~= 1.0
```

`h_eq = 150 W/m2/K` for `200 s`:

```text
output = testModule/v13_start_h150_tec_off_200s_20260629
t=5432 s
console mean emitter ~= 1295.49 K
TEC input mean emitter ~= 1288.62 K
TEC input mean collector ~= 900.24 K
R_total=100 ohm probe: U ~= 37.514 V, I ~= 0.3751 A, P ~= 14.07 W
R_total=1e6 ohm probe: U ~= 38.264 V
```

`h_eq = 200 W/m2/K` for `200 s`:

```text
output = testModule/v13_start_h200_tec_off_200s_20260629
t=5432 s
console mean emitter ~= 1241.66 K
TEC input mean emitter ~= 1235.09 K
TEC input mean collector ~= 906.57 K
R_total=100 ohm probe: U ~= 29.542 V, I ~= 0.2954 A, P ~= 8.73 W
R_total=1e6 ohm probe: U ~= 29.881 V
```

Interpretation: lowering the Cs-filled gap conductance from `250` to `200 W/m2/K` during the post-cesium thermal hold is already enough to exceed the rated-voltage gate at official `110 kW`; `150 W/m2/K` overshoots substantially. This points to the Cs-gap thermal trajectory as a plausible startup-control lever, unlike the earlier diagnostic `150 kW` power increase.

### 2026-06-29 h_eq=200 automatic fixed-R to fixed-U smoke

A short automatic electrical-start smoke was then run from the `h_eq=200` thermal hold restart:

```text
restart = testModule/v13_start_h200_tec_off_200s_20260629/v13_start_h200_tec_off_200s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedr_to_fixedu_2s_20260629
startup power = 110000 W
initial mode = fixed_r
R_total = 100 ohm
switch voltage = 27.2 V
target fixed voltage = 27.2 V
backend = verified rebuilt pyd via THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release
```

The automatic switch triggered at `t=5432.5 s`:

```text
t=5432.5 s: fixed_r, U ~= 29.767 V, I ~= 0.2977 A, P ~= 8.86 W, converged=True, iter=2
t=5433.0 s: fixed_u, U = 27.2 V, I ~= 1.2977 A, P ~= 35.30 W, converged=True, iter=29
t=5433.5 s: fixed_u, U = 27.2 V, I ~= 1.2978 A, P ~= 35.30 W, converged=True, iter=42
t=5434.0 s: fixed_u, U ~= 27.178 V, I ~= 1.3272 A, P ~= 36.07 W, converged=False, iter=47
```

Conclusion: the official `110 kW` startup can reach and trigger the rated-voltage gate if the Cs-filled gap conductance trajectory is made more insulating, with `h_eq=200 W/m2/K` as a useful first candidate. However, the fixed-voltage branch is still too slow and marginal for long steady continuation: even with the rebuilt pyd and augmented lookup, a `2 s` smoke took several minutes and the final record was close to target but non-converged. Before launching a long fixed-U steady run, improve fixed-U iteration performance/robustness or add a staged voltage-control transition.

### 2026-06-29 fixed-U startup switch optimization after h_eq=200 smoke

The first `h_eq=200` automatic fixed-R to fixed-U smoke proved the thermal path can trigger the `27.2 V` gate, but fixed-U was too slow and marginal:

```text
baseline output = testModule/v13_start_h200_fixedr_to_fixedu_2s_20260629
wall ~= 390 s for 2 s physical time
fixed-U records: iter 29 / 42 / 47
last record: U ~= 27.178 V, converged=False
```

Root cause: after switching from `R_total=100 ohm`, the fixed-R current is only about `0.298 A`, while the fixed-U operating current is near `1 A`. The series fixed-U solver was using the switch current as its first guess and then spending many full `circuitCalc()` calls in broad bracket/secant search.

Source-side change in `ThermoCalc/circuitTECs.cpp::uFixedCircuitCalc()`:

- keep the bounded bracket/secant structure and low-temperature no-hang guards;
- use `0.05 V` as the fixed-voltage engineering residual tolerance;
- prioritize candidate currents around `I_guess + 1 A` before the wider fallback samples, which matches the V13 switch current jump;
- keep fixed-U public output semantics as `Uout=Utarget` on convergence.

Verification used only the rebuilt test pyd in `ThermoCalc/build_cp312/Release` through `THERMOCALC_PYD_DIR`; the root production pyd was not overwritten.

Regression checks:

```text
testModule/test_thermocalc_interface.py: passed
testModule/test_thermocalc_lookup.py: passed, lookup speedup ~= 37.6x
```

Optimized smoke:

```text
output = testModule/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629
wall ~= 88 s for 2 s physical time
switch at t=5432.5 s from fixed_r U ~= 29.767 V

fixed-U records:
t=5433.0 s: U=27.2 V, I ~= 0.9644 A, P ~= 26.23 W, converged=True, iter=15
t=5433.5 s: U=27.2 V, I ~= 1.1383 A, P ~= 30.96 W, converged=True, iter=3
t=5434.0 s: U=27.2 V, I ~= 1.3208 A, P ~= 35.93 W, converged=True, iter=4
```

Current implication: the official `110 kW`, `h_eq=200 W/m2/K` startup path is now numerically viable for a short fixed-R to fixed-U transition smoke with the rebuilt test pyd. It is still not ready to claim steady state: the root production pyd has not been replaced, the hydraulic solver still reports the known first-step residual warning, and a longer fixed-U continuation is needed to verify lookup coverage and stable energy balance.

### 2026-06-29 fixed-U short continuation after sample-order optimization

After the `I_guess + 1 A` prioritized sample order in `uFixedCircuitCalc()`, the `h_eq=200` fixed-R to fixed-U smoke was repeated:

```text
output = testModule/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629
wall ~= 88 s for 2 s physical time
switch at t=5432.5 s from fixed_r U ~= 29.767 V
fixed-U records all converged:
  t=5433.0 s: U=27.2 V, I ~= 0.9644 A, P ~= 26.23 W, iter=15
  t=5433.5 s: U=27.2 V, I ~= 1.1383 A, P ~= 30.96 W, iter=3
  t=5434.0 s: U=27.2 V, I ~= 1.3208 A, P ~= 35.93 W, iter=4
```

A fixed-voltage continuation from that restart also completed:

```text
restart = testModule/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus5s_20260629
mode = fixed_u
target voltage = 27.2 V
wall ~= 163 s for 5 s physical time

records:
t=5435 s: U=27.2 V, I ~= 1.4214 A, P ~= 38.66 W, iter=16, converged=True
t=5436 s: U=27.2 V, I ~= 1.8037 A, P ~= 49.06 W, iter=3, converged=True
t=5437 s: U=27.2 V, I ~= 2.2148 A, P ~= 60.24 W, iter=3, converged=True
t=5438 s: U=27.2 V, I ~= 1.8341 A, P ~= 49.89 W, iter=1, converged=True
t=5439 s: U=27.2 V, I ~= 1.7341 A, P ~= 47.17 W, iter=7, converged=True
```

Current status: with the rebuilt test pyd, the official `110 kW`, `h_eq=200 W/m2/K` startup route can pass the fixed-R to fixed-U transition and sustain a short fixed-U continuation. This is still not a steady result. The next calculation should extend fixed-U in moderate chunks while monitoring lookup coverage, hydraulic residuals, current/power oscillation, and energy balance before attempting an overnight run.

### 2026-06-29 fixed-U plus20s continuation status

The optimized fixed-U path was continued from `t=5439 s` for another `20 s`:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus5s_20260629/v13_start_h200_fixedu27p2_plus5s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus20s_20260629
mode = fixed_u
target voltage = 27.2 V
thermo update interval = 0.5 s
wall ~= 947 s for 20 s physical time
```

All 10 recorded fixed-U points converged. Endpoint:

```text
t = 5459 s
U = 27.2 V
I ~= 1.1384 A
P_e ~= 30.97 W
mean emitter ~= 1240.59 K
core inlet/outlet ~= 745.80 / 838.74 K
q_radiator_total ~= 107.46 kW
coolant enthalpy rise ~= 105.34 kW
core_heat - coolant_enthalpy - electric ~= 4.63 kW
```

Last-10-s trends:

```text
dI/dt ~= -2.69e-2 A/s
dP_e/dt ~= -0.732 W/s
dT_emitter/dt ~= +5.18e-2 K/s
dT_inlet/dt ~= -1.06e-1 K/s
dq_rad/dt ~= -72.5 W/s
```

Interpretation: the fixed-U continuation is numerically stable over this medium segment, but it is not near steady. Electrical output is still relaxing downward while the emitter temperature is slowly recovering. Continue in moderate chunks before attempting an overnight steady run.

### 2026-06-29 fixed-U plus40s/plus60s continuation status

The optimized `h_eq=200 W/m2/K`, `fixed_u=27.2 V` path was continued in two additional 20 s chunks, still with `thermo_update_interval=0.5 s` and the rebuilt test pyd.

`plus40s` segment:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus20s_20260629/v13_start_h200_fixedu27p2_plus20s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus40s_20260629
wall ~= 490 s for 20 s physical time
all records converged
endpoint t=5479 s:
  U = 27.2 V
  I ~= 1.0404 A
  P_e ~= 28.30 W
  mean emitter ~= 1240.80 K
  core inlet/outlet ~= 745.34 / 838.81 K
  q_radiator_total ~= 107.46 kW
  coolant enthalpy rise ~= 105.93 kW
  core_heat - coolant_enthalpy - electric ~= 4.04 kW
last-10-s trends:
  dI/dt ~= -1.93e-2 A/s
  dP_e/dt ~= -0.525 W/s
  dT_emitter/dt ~= +4.67e-2 K/s
```

`plus60s` segment:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus40s_20260629/v13_start_h200_fixedu27p2_plus40s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus60s_20260629
wall ~= 856 s for 20 s physical time
all records converged
endpoint t=5499 s:
  U = 27.2 V
  I ~= 1.0344 A
  P_e ~= 28.13 W
  mean emitter ~= 1240.94 K
  core inlet/outlet ~= 745.32 / 838.84 K
  q_radiator_total ~= 107.45 kW
  coolant enthalpy rise ~= 106.00 kW
  core_heat - coolant_enthalpy - electric ~= 3.97 kW
last-10-s trends:
  dI/dt ~= -1.71e-2 A/s
  dP_e/dt ~= -0.464 W/s
  dT_emitter/dt ~= +4.50e-2 K/s
```

Interpretation: fixed-U electrical convergence is now robust over the tested `60 s` continuation after switching, but the coupled thermal state is not steady. The residual thermal imbalance remains about `4 kW`, and the emitter is still rising slowly while electric output relaxes downward. The run is also still expensive at `0.5 s` TEC update frequency. Before an overnight run, either accept the cost and continue in longer chunks, or test a larger TEC update interval / further C++ inner-loop optimization.

### 2026-06-29 fixed-U TEC update interval 1.0 s check

A speed/accuracy check was run from the `t=5499 s` fixed-U restart by increasing the ThermoCalc update interval from `0.5 s` to `1.0 s` while keeping `max_dt=0.5 s`:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus60s_20260629/v13_start_h200_fixedu27p2_plus60s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus80s_tec1s_20260629
duration = 20 s
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
wall ~= 169 s for 20 s physical time
```

All records converged. Endpoint:

```text
t = 5519 s
U = 27.2 V
I ~= 1.0645 A
P_e ~= 28.95 W
mean emitter ~= 1241.06 K
core inlet/outlet ~= 745.35 / 838.87 K
q_radiator_total ~= 107.47 kW
coolant enthalpy rise ~= 106.00 kW
core_heat - coolant_enthalpy - electric ~= 3.97 kW
```

Last-10-s trends:

```text
dI/dt ~= -1.50e-2 A/s
dP_e/dt ~= -0.408 W/s
dT_emitter/dt ~= +4.41e-2 K/s
dT_inlet/dt ~= +5.06e-2 K/s
dq_rad/dt ~= -36.4 W/s
```

Comparison to the previous `0.5 s` segment indicates the physical trend is consistent, while wall time improved substantially. The case is still not steady, but `thermo_update_interval=1.0 s` is a reasonable setting for the next longer fixed-U continuation.

### 2026-06-29 fixed-U plus180s with 1.0 s TEC update

A longer fixed-U continuation was run after accepting `thermo_update_interval=1.0 s` as a speed-improving setting:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus80s_tec1s_20260629/v13_start_h200_fixedu27p2_plus80s_tec1s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus180s_tec1s_20260629
duration = 100 s
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
wall ~= 1267 s for 100 s physical time
```

All 10 recorded points converged. Endpoint:

```text
t = 5619 s
U = 27.2 V
I ~= 0.9685 A
P_e ~= 26.34 W
mean emitter ~= 1243.84 K
core inlet/outlet ~= 745.63 / 840.81 K
q_radiator_total ~= 107.98 kW
coolant enthalpy rise ~= 107.87 kW
core_heat - coolant_enthalpy - electric ~= 2.11 kW
```

Last-50-s trends:

```text
dI/dt ~= 0 A/s
dP_e/dt ~= 0 W/s
dT_emitter/dt ~= +2.92e-2 K/s
dT_inlet/dt ~= +1.53e-2 K/s
dq_rad/dt ~= +8.59 W/s
```

Interpretation: the fixed-voltage electrical solution is now stable over a `100 s` continuation and electric output has flattened near `26.3 W` for this low-power startup state. The thermal system is still not steady: emitter temperature is still rising and the residual heat balance remains about `2.1 kW`. Continue in longer but still bounded chunks, for example `300-500 s`, before declaring a near-steady startup fixed-U state or launching an overnight run.

### 2026-06-29 fixed-U plus480s near-steady continuation

A `300 s` continuation was run from the `t=5619 s` restart using the accepted `thermo_update_interval=1.0 s` setting:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus180s_tec1s_20260629/v13_start_h200_fixedu27p2_plus180s_tec1s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus480s_tec1s_20260629
duration = 300 s
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
wall ~= 5618 s for 300 s physical time
```

All 10 recorded points converged. Endpoint:

```text
t = 5919 s
U = 27.2 V
I ~= 0.9552 A
P_e ~= 25.98 W
mean emitter ~= 1248.11 K
core inlet/outlet ~= 747.86 / 843.78 K
q_radiator_total ~= 109.24 kW
coolant enthalpy rise ~= 108.72 kW
core_heat - coolant_enthalpy - electric ~= 1.26 kW
```

Last-150-s trends:

```text
dI/dt ~= +2.08e-4 A/s
dP_e/dt ~= +5.66e-3 W/s
dT_emitter/dt ~= +9.72e-3 K/s
dT_inlet/dt ~= -4.22e-3 K/s
dq_rad/dt ~= -2.06 W/s
d(coolant enthalpy rise)/dt ~= +5.95 W/s
```

Interpretation: this is close to a stable fixed-voltage startup state but not a strict steady state. Electrical output is nearly flat near `26 W`, while the thermal system still has about `1.26 kW` residual heat imbalance and a small positive emitter-temperature drift. Another several-hundred-second continuation should reduce the remaining drift, but the current run speed is about `18.7 wall-s / physical-s`, so an overnight run is appropriate only if this cost is acceptable.

### 2026-06-29 fixed-U plus1780s continuation

A `1000 s` official continuation from the `t=6219 s` restart completed under the formal `110 kW` startup power, `h_eq=200 W/m2/K`, `fixed_u=27.2 V`, and `thermo_update_interval=1.0 s` settings:

```text
output = testModule/v13_start_h200_fixedu27p2_plus1780s_tec1s_20260629
t = 7219 s
U = 27.2 V
I ~= 1.3873 A
P_e ~= 37.74 W
mean emitter ~= 1247.23 K
core inlet/outlet ~= 746.84 / 842.74 K
q_radiator_total ~= 108.69 kW
coolant enthalpy rise ~= 108.70 kW
core_heat - coolant_enthalpy - electric ~= 1.27 kW
```

All fixed-U records converged; the only stderr entry was the known first-step hydraulic residual warning at the restart boundary. Last-five-record slopes were small (`dP_e/dt ~= -1.86e-5 W/s`, `dT_emitter/dt ~= +1.21e-4 K/s`, `dq_radiator/dt ~= +4.89e-2 W/s`), so the state is close but still not a strict thermal steady state. A `5000 s` continuation is running from this restart in `testModule/v13_start_h200_fixedu27p2_plus6780s_tec1s_20260629` with PID `58580`.
### 2026-06-29 fixed-U plus6780s long continuation

A `5000 s` official continuation from the `t=7219 s` restart completed under the formal `110 kW` startup power, `h_eq=200 W/m2/K`, `fixed_u=27.2 V`, and `thermo_update_interval=1.0 s` settings:

```text
output = testModule/v13_start_h200_fixedu27p2_plus6780s_tec1s_20260629
t = 12219 s
U = 27.2 V
I ~= 1.3838 A
P_e ~= 37.64 W
mean emitter ~= 1247.53 K
core inlet/outlet ~= 747.12 / 843.17 K
q_radiator_total ~= 108.86 kW
coolant enthalpy rise ~= 108.86 kW
core_heat - coolant_enthalpy - electric ~= 1.10 kW
```

All history records converged in fixed-U mode with `iteration_count = 1` after the restart transient. The only stderr entry was the known first-step hydraulic residual warning at `t=7219 s`. Recent trends over the last five records were still positive but small:

```text
dP_e/dt ~= -1.86e-5 W/s
dT_emitter/dt ~= +5.86e-5 K/s
dT_inlet/dt ~= +5.41e-5 K/s
dq_radiator/dt ~= +3.27e-2 W/s
d(coolant enthalpy rise)/dt ~= +3.27e-2 W/s
```

Interpretation: the official fixed-U startup path remains stable and near steady, but the remaining storage term is still about `1.10 kW`. A longer `30000 s` steady approach run was launched from this restart:

```text
output = testModule/v13_start_h200_fixedu27p2_plus36780s_tec1s_20260629
pid = 19504
```
### 2026-06-29 fixed-U plus36780s first history checkpoint

The `30000 s` long steady-approach run has written its first formal history record, confirming that the restart-load `tec_solver_converged=False` flag was only an initial latest-state artifact. The first recorded point is converged:

```text
output = testModule/v13_start_h200_fixedu27p2_plus36780s_tec1s_20260629
pid = 19504
record = first history row
t = 14219 s
relative time = 2000 s
startup thermal power = 110000 W
mode = fixed_u
converged = True
iteration_count = 1
U = 27.2 V
I ~= 1.3826 A
P_e ~= 37.61 W
mean emitter ~= 1247.64 K
core inlet/outlet ~= 747.22 / 843.32 K
q_radiator_total ~= 108.92 kW
coolant enthalpy rise ~= 108.92 kW
core_heat - coolant_enthalpy - electric ~= 1.04 kW
```

Drift from the previous segment endpoint (`12219 -> 14219 s`):

```text
dP_e/dt ~= -1.60e-5 W/s
dT_emitter/dt ~= +5.46e-5 K/s
dT_inlet/dt ~= +5.04e-5 K/s
dq_radiator/dt ~= +3.05e-2 W/s
d(coolant enthalpy rise)/dt ~= +3.05e-2 W/s
```

Interpretation: the official `110 kW` fixed-U path is continuing normally and remains numerically stable. The thermal balance is still closing slowly, with the residual storage term reduced from about `1.10 kW` to about `1.04 kW` over this first `2000 s` checkpoint. Continue the long run before claiming strict steady state.
### 2026-06-29 correction: cesium TEC gap h_eq should be 29 W/m2/K

User review caught that the previous `h_eq=200 W/m2/K` fixed-U continuations used an artificial sensitivity value, not the physical cesium-vapor TEC gap setting. The V7 steady CaseA configuration uses `tec_gap_config h_eq=29.0 W/m2/K`, and the V13 cold-start cesium-filled gap should be consistent with that value unless explicitly running a sensitivity case.

Code correction:

```text
testModule/v13_startup_control.py: V13StartupControlConfig.cesium_gap_h_eq_w_m2_k default = 29.0
testModule/run_v13_start_case.py: --cesium-gap-h-eq-w-m2-k default = 29.0
testModule/test_v13_startup_control.py: added defaults tests for config and runner CLI
```

Verification:

```text
testModule/test_v13_startup_control.py: passed
python -m py_compile testModule/v13_startup_control.py testModule/run_v13_start_case.py testModule/test_v13_startup_control.py: passed
```

Therefore, the earlier `h_eq=200` results must be treated as numerical/control diagnostics only, not as official V13 physical startup results.

Corrected `h_eq=29` restart path:

```text
base restart = testModule/v13_start_cesium_conditioning_plus1000s_20260628/v13_start_cesium_conditioning_plus1000s_20260628_latest_restart.npz
base state: t=4090 s, Cs fraction ~= 1, TEC disabled, old h_eq ~= 250 W/m2/K
```

A corrected thermal hold with TEC disabled was run using `h_eq=29.0`:

```text
output = testModule/v13_start_h29_tec_off_200s_20260629
t = 4290 s
core power = 110000 W
TEC disabled
h_eq = 29.0 W/m2/K
core inlet/outlet ~= 722.28 / (history endpoint) K
mean emitter ~= 1513.22 K
```

Then a `2 s` fixed-R to fixed-U smoke was run:

```text
output = testModule/v13_start_h29_fixedr_to_fixedu_2s_20260629
fixed-R at t=4290.5 s: U ~= 65.04 V, I ~= 0.650 A, P_e ~= 42.3 W, converged=True
switch to fixed-U 27.2 V triggered immediately after the first fixed-R record
fixed-U at t=4292.0 s: U=27.2 V, I ~= 351.17 A, P_e ~= 9.55 kW, converged=True
```

Short fixed-U continuation:

```text
output = testModule/v13_start_h29_fixedu27p2_plus20s_20260629
t = 4312 s
U=27.2 V
I ~= 259.10 A
P_e ~= 7.05 kW
mean emitter ~= 1490.42 K
core inlet/outlet ~= 728.01 / 818.14 K
core_heat - coolant_enthalpy - electric ~= 0.76 kW
```

Medium fixed-U continuation:

```text
output = testModule/v13_start_h29_fixedu27p2_plus220s_20260629
t = 4512 s
U=27.2 V
I ~= 211.67 A
P_e ~= 5.76 kW
mean emitter ~= 1506.67 K
core inlet/outlet ~= 731.35 / 819.28 K
q_radiator_total ~= 99.65 kW
coolant enthalpy rise ~= 99.69 kW
core_heat - coolant_enthalpy - electric ~= 4.55 kW
```

A `1000 s` fixed-U continuation with `thermo_update_interval=1.0 s` completed:

```text
output = testModule/v13_start_h29_fixedu27p2_plus1220s_tec1s_20260629
t = 5512 s
U=27.2 V
I ~= 209.72 A
P_e ~= 5.70 kW
mean emitter ~= 1541.95 K
core inlet/outlet ~= 735.62 / 825.70 K
q_radiator_total ~= 102.08 kW
coolant enthalpy rise ~= 102.11 kW
core_heat - coolant_enthalpy - electric ~= 2.19 kW
```

Interpretation: with the corrected `h_eq=29.0 W/m2/K`, the V13 startup enters a meaningful TEC generation regime. The fixed-R voltage gate is crossed naturally, fixed-U solves converge, and electric power is now in the expected kilowatt range rather than the invalid tens-of-watts result from `h_eq=200`. The state is still not strict steady because the emitter and radiator heat rejection are drifting; a `5000 s` continuation is running:

```text
output = testModule/v13_start_h29_fixedu27p2_plus6220s_tec1s_20260629
pid = 53696
```
### 2026-06-29 correction: TFE ignition timing for cesium gap and fixed-R startup

User clarified the startup sequence: at `critical_time + 1500 s`, TFE ignition should immediately replace the emitter-collector gap equivalent heat-transfer coefficient with the cesium-vapor value `h_eq=29.0 W/m2/K`. The fixed-resistance external circuit should participate from this ignition point so voltage/current develop while the emitter warms; once the terminal voltage reaches `27.2 V`, the main circuit switches to fixed total voltage.

This supersedes the previous workflow that first ran a separate TEC-off thermal hold after cesium conditioning. Those TEC-off hold runs remain useful diagnostics, but are not the formal startup sequence.

Code updates:

```text
testModule/v13_startup_control.py
  - default cesium_gap_h_eq_w_m2_k = 29.0
  - TFE ignition latches by time after critical, not by emitter-temperature gate
  - once ignition latches, cs_fraction = 1.0 and h_eq immediately equals the cesium value
  - default electrical start gates are zero, so fixed-R TEC coupling starts at TFE ignition

testModule/run_v13_start_case.py
  - --cesium-gap-h-eq-w-m2-k default = 29.0
  - --tec-electrical-start-after-cesium-s default = 0.0
  - --tec-electrical-start-cs-fraction default = 0.0
  - --tec-electrical-start-emitter-temperature-k default = 0.0

testModule/test_v13_startup_control.py
  - added/updated tests for TFE ignition immediately setting h_eq=29 and enabling fixed-R TEC
```

Verification:

```text
testModule/test_v13_startup_control.py: passed
python -m py_compile testModule/v13_startup_control.py testModule/run_v13_start_case.py testModule/test_v13_startup_control.py: passed
```

Corrected sequence test from the `t=1590 s` pre-ignition restart:

```text
base restart = testModule/v13_start_corrected_shield_startup_1590s_20260628/v13_start_corrected_shield_startup_1590s_20260628_latest_restart.npz
base state: t=1590 s, time_after_critical ~= 1492.7 s, TEC off, h_eq=600 W/m2/K, mean emitter ~= 1031.2 K
```

With the corrected controller, TEC coupling enabled automatically at `t ~= 1597.342 s` and `h_eq=29.0` was applied from ignition.

Fixed-R load tests:

```text
R_total = 0.0044 ohm
output = testModule/v13_start_h29_tfe_ignition_fixedr_300s_20260629
stable fixed-R, but voltage only rose to ~= 2.23 V by t=1890 s; no switch.

R_total = 0.05 ohm
output = testModule/v13_start_h29_tfe_ignition_fixedr_R005_500s_20260629
stable fixed-R, voltage rose to ~= 16.40 V by t=2090 s; no switch.

R_total = 0.083 ohm from the warmed R=0.05 restart
output = testModule/v13_start_h29_fixedr_R0083_fromR005_200s_20260629
stable fixed-R, voltage rose to ~= 21.97 V by t=2290 s; no switch.

R_total = 0.105 ohm from the warmed R=0.083 restart
output = testModule/v13_start_h29_fixedr_R0105_fromR0083_100s_20260629
stable fixed-R, voltage rose to ~= 24.68 V by t=2390 s; no switch.

R_total = 0.12 ohm from the warmed R=0.105 restart
output = testModule/v13_start_h29_fixedr_R012_fromR0105_50s_20260629
stable fixed-R, voltage rose to ~= 26.20 V by t=2440 s; no switch.

R_total = 0.125 ohm from the warmed R=0.12 restart
output = testModule/v13_start_h29_fixedr_R0125_fromR012_20s_20260629
stable fixed-R, peak voltage ~= 26.72 V; no switch.

R_total = 0.131 ohm from the warmed R=0.125 restart
output = testModule/v13_start_h29_fixedr_R0131_fromR013_10s_20260629
fixed-R at t=2480.5 s: U ~= 27.246 V, I ~= 207.98 A, P_e ~= 5.67 kW
automatic switch to fixed-U 27.2 V succeeded; subsequent fixed-U records converged.
```

A direct `R_total=0.10 ohm` run from cold TFE ignition was attempted, but it consumed CPU without writing a first history record and was stopped. A `per_tec` interpretation of `0.0044 ohm` also failed by producing non-finite axial Joule heat after early records. Therefore, the currently stable route is staged fixed-R resistance from low value to approximately `0.131 ohm` as the emitter warms, followed by fixed-U.

Fixed-U continuation after successful switch:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus1000s_20260629
t = 3490 s
U = 27.2 V
I ~= 209.68 A
P_e ~= 5.70 kW
mean emitter ~= 1543.14 K
core inlet/outlet ~= 735.61 / 825.68 K
q_radiator_total ~= 102.07 kW
coolant enthalpy rise ~= 102.10 kW
core_heat - coolant_enthalpy - electric ~= 2.20 kW
```

The result is in the expected kilowatt range and comparable to the V11 electrical output scale, but it is not yet strict steady state. A corrected-sequence `5000 s` fixed-U continuation is now running:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus6000s_tec1s_20260629
pid = 77012
```
### 2026-06-29 corrected ignition fixed-U plus6000s checkpoint

The corrected-sequence fixed-U continuation from the successful `R_total=0.131 ohm` switch completed a further `5000 s` segment:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus6000s_tec1s_20260629
start = 3490 s
end = 8490 s
startup power = 110000 W
h_eq = 29.0 W/m2/K
mode = fixed_u
U = 27.2 V
```

Endpoint:

```text
t = 8490 s
I ~= 209.47 A
P_e ~= 5.70 kW
mean emitter ~= 1544.62 K
core inlet/outlet ~= 736.45 / 826.92 K
q_radiator_total ~= 102.55 kW
coolant enthalpy rise ~= 102.55 kW
core_heat - coolant_enthalpy - electric ~= 1.75 kW
```

All records converged with `tec_solver_iteration_count = 1` after restart. The only stderr entry was the known first-step hydraulic residual warning at the restart boundary. Recent slopes over the last five records:

```text
dP_e/dt ~= -8.16e-4 W/s
dT_emitter/dt ~= +3.10e-5 K/s
dT_inlet/dt ~= +5.07e-5 K/s
dq_radiator/dt ~= +2.93e-2 W/s
d(coolant enthalpy rise)/dt ~= +2.93e-2 W/s
```

Interpretation: the corrected TFE ignition path is stable and in the expected kilowatt electrical-output range, but it is still not a strict steady state because the residual storage term is about `1.75 kW`. A longer `30000 s` fixed-U continuation is now running:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629
pid = 79504
```
### 2026-06-29 V13 cold-start fixed-U handoff residual interpretation

The corrected V13 cold-start path now treats TFE ignition at `critical_time + 1500 s` as the point where the TEC gap is set directly to the V7 steady Cs-vapor value `h_eq = 29.0 W/m2/K` and the main TEC electrical solve is enabled. The stable demonstrated handoff is staged total fixed-R warmup to about `0.131 ohm`, followed by automatic switch to fixed-U `27.2 V`.

The latest long continuation was stopped by request instead of forcing strict steady state:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629
history final time = 32490 s
latest saved restart = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629_latest_restart.npz
```

At the final history row, `q_radiator_total - coolant_enthalpy_rise = -0.136 W`, while `core_heat - coolant_enthalpy_rise - TEC_electric = 1329.281 W`. Since the radiator/coolant balance is sub-watt and the TEC fixed-U solve converged in one iteration with finite outputs, this residual should be interpreted as transient storage in the still-warming core/TFE/structural solids, not as a hydraulic/radiator/ThermoCalc accounting failure.

For later long continuations, use the latest restart above and track `core_heat - coolant_dh - electric` as the storage residual. A practical near-steady acceptance gate is to reduce that residual below the task-specific tolerance, for example `<0.5 kW` or `<1%` core power, and confirm small final-window slopes in inlet/outlet coolant temperature, mean emitter temperature and radiator heat rejection.

Restart timestamp note: in the stopped `plus36000s` run, the CSV history reached `t = 32490 s`, but the latest restart file stores `System/global_time = 28490 s` and `System/last_dt = 0.5 s`. Use the restart for resumable continuation and the CSV final row for the latest residual diagnosis.

### 2026-06-29 residual diagnostic fields added to V13 startup runner

`testModule/run_v13_start_case.py` now writes derived energy residual diagnostics into every startup history row and `latest_state.json` `latest_record`:

```text
core_heat_minus_coolant_enthalpy_minus_electric_w
core_heat_minus_radiator_minus_electric_w
corrected_core_energy_residual_w = core_heat_power - coolant_enthalpy_rise - terminal_electric - tec_wire_joule_loss
corrected_loop_energy_residual_w = core_heat_power + radiator_tube_external_heat_w - q_radiator_total_w - terminal_electric - tec_wire_joule_loss
radiator_minus_coolant_enthalpy_w
tec_wire_joule_loss_w
core_energy_storage_residual_rel
radiator_coolant_balance_rel
```

`q_radiator_total_w - coolant_enthalpy_rise_w` remains the main loop closure check for radiator/coolant exchange. `radiator_tube_external_heat_w` is an explicit heat input term and is **not** treated as an error term. For transient runs, interpret `corrected_core_energy_residual_w` as the storage imbalance of core plus connected structure after excluding wire Joule loss; `corrected_loop_energy_residual_w` includes the explicit external heat input as a source. Fins in this runner are quasi-steady and do not add a separate storage term in the startup residual bookkeeping.

Verification:

```text
E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe testModule\test_v13_startup_control.py
E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m py_compile testModule\run_v13_start_case.py testModule\test_v13_startup_control.py
```

### 2026-06-30 local implicit fluid-solid dt limit

`SystemManager.compute_adaptive_dt()` now receives a strict fluid-solid coupling time-scale limit even when V11/V13 use `solid_ode_method=implicit_euler` and `fluid_solid_coupling_scheme=local_implicit`. The limiting coupler value is `safety_factor * min(C_eff / lambda)` after the first coupler execution. This prevents thin radiator/header walls from being advanced only by `max_dt` when local implicit exchange is enabled.

Verification added:

```text
E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m unittest testModule.test_local_implicit_heat_exchange testModule.test_system_manager_lifecycle
```

For V11/V13 timing comparisons, record `coupling_tau_min_s`, `coupling_dt_limit_s`, and `dt_over_coupling_tau_max` from coupler diagnostics when diagnosing unexpectedly small adaptive time steps or wall/fluid temperature reversals.

### 2026-07-10 V14_10kW powered debug baseline

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/` is the active 10 kW V14 heat-pipe-radiator package. It should be treated as a local branch of the full-loop work: keep new 10 kW modeling changes inside this directory unless the user explicitly asks to promote them back into `Full_Loop_Cases` or shared component code.

Primary docs:

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/README.md
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_debug/README.md
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_debug/TUNING_LOG.md
```

Primary powered runner:

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_debug/run_v14_210kw_debug.py
```

Run assumptions used by the current debug baseline:

```text
core power = 210000 W
loop flow target = 2.46 kg/s
solid solver = implicit_euler
space temperature = 4 K
external orbital heat flux = disabled
TEC lookup = enabled for powered tuning runs
main TEC target = about 206 A and 10.44 kW net electric power
```

Use direct runner invocations for `V14_210kW_debug`. Avoid `python -m unittest` for this staged debug workflow; it previously stalled and is not the intended entry. For syntax-level checks, use `py_compile` on the edited runner/builder files.

Current best restart:

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_debug/runs/final_eps07475_u50p65_wire0335_1200s_from7964/stage_01_restart.npz
```

Current best parameters:

```text
radiator_emissivity = 0.7475
tec_voltage = 50.65 V
wire_resistance_scale = 0.335
hp_up_view_factor = 0.0
upper_hp_down_view_factor = 0.3
lower_hp_down_view_factor = 0.4
```

Endpoint summary at `t=9163.85 s`:

```text
core inlet/outlet = 754.738 / 845.773 K
TEC current = 206.569 A
TEC net electric power = 10.463 kW
total radiator rejection = 195.212 kW
required pump head = 27.880 kPa
```

Interpretation: the result is close to the requested `754.45 K` inlet, `845.65 K` outlet, `206 A`, and `10.44 kW`, but remains slightly hot and slightly high in current/power. Treat it as a working continuation baseline. If tighter calibration is required, make small changes to radiator emissivity or wire-resistance scale and validate with at least a 1200 s window; 1 s electrical checks are not reliable final evidence for this case.

### 2026-07-19 V14 210 kW 全堆氦气瞬时完全失压事故

事故算例位于：

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/
```

首读 `README.md`，运行入口为 `run_v14_helium_depressurization.py`。独立初态是正常材料热容
长期计算后的约 `13864.2 s` restart，保存在该目录 `initial_state/`；不要从早期缩热容状态
或 0.1 s smoke 输出启动正式事故。

事故模型仅在算例层把 5 个代表性 TFE（倍率合计 58）的 `collector_iclad_gap.k_gas`
同步清零，即 `h_He: 5678 -> 0 W/(m2*K)`。间隙辐射、几何和发射率保持不变；这是一种
瞬时完全丧失气体导热的保守代理模型，不是显式氦气压力或泄漏流动模型。固定功率源关闭，
外界反应性为零，控制鼓关闭，点堆使用 ReactorCore 的相对温度反馈自动推进。

默认首轮为 `100 s`、固定步长 `0.05 s`、TEC 每 `0.05 s` 刷新、每 `0.1 s` 记录、
每 `10 s` 保存 checkpoint，
固体使用 `implicit_euler`，流固耦合沿用 `local_implicit`，目标总流量 `2.46 kg/s`，
外热关闭。每步检查壁温、芯块、接收极、慢化剂和反射层限值；触发后保存
`emergency_restart.npz` 和 `limit_trip.json`。

事故状态不保存在 `.npz` 的气隙参数中，必须与同目录 `run_config.json` 配套使用。事故续算
时 runner 根据 `helium_accident_active=true` 重新施加零气体导热，并保留原事故绝对时刻和
已保存点堆状态，不能只复制 `.npz` 后脱离配置运行。

2026-07-19 的 0.1 s 真实 smoke 与事故 restart 续算均通过。事故 runner 对 TEC 名义刷新
周期施加小量浮点调度容差，避免约 `13864.2 s` 绝对时间下 `0.05 s` 差值被表示为略小于
阈值而漏更；`run_config.json` 同时记录名义周期和内部调度阈值。正式计算在事故后
`1.75 s` 因 Ring1 接收极在 `z=0.29874 m` 达到 `1024.851 K`、超过 `1023 K` 限值而
安全终止。`dt=0.025 s` 的独立最终敏感性计算在 `1.725 s` 的同一位置以 `1023.120 K`
触发。两者 restart 中末次 TEC 更新时间分别只落后全局时间 `0.05 s` 和 `0.025 s`。
推荐结果目录为 `runs/accident_100s_final/`，敏感性目录为
`runs/sensitivity_dt0p025_final/`；其他同名前缀目录只保留作修正前对照。

事故态 restart 在推进前执行完整温限和数值预检；若已越限则以
`phase=restart_preflight` 原时刻停止。推进后若水力或 TEC 不收敛，或任一数值诊断出现
NaN/Inf，也会 fail-closed 保存紧急 restart 和 `limit_trip.json`。
检查同时直接覆盖水力 `T/P/h/rho/W` 原始数组及全部注册固体温度数组，避免汇总层
`nanmin/nanmax` 掩盖局部非有限值。事故 restart 预检前会先按当前温度显式刷新 TEC，
再要求 TEC 收敛。

## 2026-07-13 ThermoCalc 串联固定电流测试

`test_thermocalc_series_fixed_current.py` 验证串联 `fixed_i` 的四类状态：单根/多根零电流开路、与 `fixed_u` 交叉验证的可实现工作点、有限的不可发电目标回退到正开路电压，以及双重失败时安全返回有限零输出。测试扩展必须从独立构建目录加载，不得覆盖生产 `.pyd`。

### 2026-07-20 V14 210 kW external-heat continuation

The independent runner is `testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_reactivity_control_external_heat/run_v14_210kw_reactivity_control_external_heat.py`. It defaults to `V14_210kW_reactivity_control/反应性控制/checkpoint_t013864s.npz`, keeps the saved absolute system time, and uses that restart time as phase zero for the `5668.144369 s` N18 external-heat history. Upper and lower rings reuse N18 columns 0 through 17. The original reactivity-control and debug runners remain external-heat-off by default. A `0.05 s` smoke from `13864.2 s` completed with finite thermal/electrical outputs and converged fluid solve.
### 2026-07-20 V14 fixed-power two-orbit external-heat case

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_fixed_power_external_heat_2orbits/` keeps the core source fixed at `210000 W` and runs the N18 external-heat history for exactly two `5668.144369 s` periods (`11336.288738 s`). It loads the current `checkpoint_t013864s.npz`, uses its saved absolute time as external-heat phase zero, and preserves the adjacent run configuration for TEC lookup and calibrated radiator settings. The original debug runner remains unchanged by default; `case_prefix` now controls output case names as its dataclass field intended.

### 2026-08-13 V14 210 kW fast-shutdown baseline

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_fast_shutdown/`
starts from the completed two-orbit fixed-power `stage_01_restart.npz` at
absolute time `25200.488738 s`. It hands the saved `210000 W` total power to
`PointReactor`, calibrates temperature feedback at the loaded hot state, and
immediately applies a persistent `-2 dollar` external reactivity. The normal
`2.46 kg/s` flow, heat-pipe/radiator model, fixed-voltage TEC circuit, and
absolute orbital external-heat phase remain unchanged. TEC electrical coupling
opens only after current falls to or below `0.01 A`; passive gap heat transfer
is retained. The default duration is one `5668.144369 s` orbit, with five CSV
histories every `1 s` and restart checkpoints every `100 s`.

The `1 s` smoke handoff split `210000 W` into `197329.474 W` fission power and
`12670.526 W` decay heat. At `t+1 s`, total power was `70647.798 W`, both pump
flows remained approximately `2.46 kg/s`, TEC output was `10603.975 W`, and
the hydraulic solve remained converged.

The runner also accepts its own point-kinetics final restart for later-orbit
continuation. It restores the saved precursor/decay-heat and feedback-reference
state instead of reinitializing point kinetics, preserves the original shutdown
time for continuous `shutdown_elapsed_s`, and restores a previously zero-current
TEC as explicitly open. Continuations use `0.01 s` for the first `2 s` to settle
the rebuilt hydraulic state, then return to `0.05 s`. The first full shutdown
orbit completed at `t+5668.144369 s` with `3222.360 W` decay heat, negligible
fission power, coolant `574.475-581.680 K`, and zero TEC output.
### 2026-07-20 V14 fixed-power instantaneous LOCA-1

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_fixed_power_LOCA_1/` starts from the two-orbit
`checkpoint_t016864s.npz` and models instantaneous complete NaK loss. It removes all
`FluidSolidCouple` instances and their installed boundary conditions, adds five epsilon=0.8
vacuum-radiation `GapCouple2D` paths between each representative TFE inner/outer clad, zeros
fluid sources and flows, and replaces `HydraulicNetwork.step_Picard()` on the case instance
with a frozen successful no-op. The fixed `210000 W` core source, TEC, solid conduction, N18
external heat, heat-pipe stored state, and external radiation remain active. Post-accident fluid
`T/P/h` are recorded as absent (`NaN`) while the pre-accident liquid arrays remain in every
node snapshot as reference data. Summary CSV and compressed node snapshots are written every
`0.2 s`. The verified `runs/smoke_0p2s_final/` result completed four `0.05 s` solid steps with
finite solid/TEC fields, 58 removed fluid-solid couplers, five vacuum gaps, and no hydraulic solve.

### 2026-07-20 LOCA-1 split postprocessing histories

`V14_210kW_fixed_power_LOCA_1` now writes four long-form CSV files at every configured record time: `history_coolant.csv` (volume T/P/h and junction W/v, categorized as core/ordinary_pipe/collector_ring), `history_solids.csv` (all flattened solid temperatures, categorized as core_structure/pipe_wall/heat_pipe), `history_electrical.csv` (terminal metrics and per-TFE/per-axial-node current density, potentials, electron heat transfer, and emitter/collector Joule power), and `history_reactivity.csv` (fuel/electrode/moderator/reflector/total feedback). The existing `history.csv` remains the compact run and energy-audit summary; compressed snapshots remain the array-level reference.

LOCA-1 also supports temperature-limit termination after every solid step. Current thresholds are collector 1500 K, emitter 3000 K, existing coolant 1058 K, moderator 930 K, and reflector 1000 K. A triggered terminal state is written even when it falls between regular record times. Complete coolant loss makes the coolant-temperature criterion inactive and records `coolant_max_T_K=NaN`.

The ε=0.2 and ε=0.5 limit runs completed under `runs/LOCA_1_eps020_until_failure_record0p5s` and `runs/LOCA_1_eps050_until_failure_record0p5s`. Both stopped on `collector_temperature_limit`: ε=0.2 at accident elapsed 19.95 s and collector 1500.747 K; ε=0.5 at 24.25 s and collector 1500.512 K. Both retained 210 kW fixed power, 0.5 s regular records plus the exact terminal state, and empty stderr logs.

LOCA-1 can now hand a fixed-power restart to `PointReactor` with `--enable-reactivity-feedback`; fixed source reapplication is then disabled. `--scram-time 5 --scram-reactivity-dollars -2` applies a persistent external `rho=-2*beta_total=-0.0158642` from accident elapsed 5 s. `--staged-recording` selects 0.5 s records to 20 s, 2 s to 100 s, 5 s to 400 s, 10 s to 600 s, and 20 s afterward. The 5.1 s integration smoke reached 76.258 kW after the scram and wrote the exact endpoint; seven focused tests pass.

For scram runs, the runner now transitions TEC to a persistent open-circuit disabled state when main current is at or below the configurable threshold (default 0.01 A). The transition calls `TFEUnit.clear_tec_sources()` for every representative TFE, preserves passive TEC-gap heat transfer, and sets `core.enable_tec_coupled=False` so later steps do not enter ThermoCalc. The switch time and state are recorded. An integration smoke verified zero current/electric power and no subsequent TEC call; eight focused tests pass. Older scram outputs without the transition are superseded by `tecopen001A` reruns.
### 2026-07-28 V14 whole-core TEC open-circuit accidents

`V14_210kW_TEC_open_circuit_accident_fixed_power/` and
`V14_210kW_TEC_open_circuit_accident_reactive_feedback/` start from the two-orbit
`checkpoint_t019865s.npz`. Both permanently clear active TEC electrical/electron/Joule
sources while retaining the emitter, collector, passive TEC-gap heat transfer, NaK loop, and
orbital external heat with the saved period and time origin. The fixed-power case remains at
210 kW until a temperature trip, then initializes point kinetics at that exact state and adds
-2 dollars. The feedback case initializes point kinetics at accident start and adds -2 dollars
to the evolving feedback at a trip. Limits are channel wall 1058 K, fuel pellet 2700 K,
collector 1023 K, moderator 930 K, and reflector 1000 K; coolant temperature is diagnostic
only. Checkpoints are written every 50 s, and the 0.1/1/10 s history schedule resets after a
scram. Focused unit tests plus real restart normal/trip smoke runs verified both control paths,
zero TEC arrays, continued external heat, converged hydraulics, and persistent -2-dollar scram.

Both formal runs completed one 5668.144369 s orbit without a trip or numerical failure and
wrote 113 checkpoints each. Fixed-power full-run maxima were wall 893.362 K, fuel 2517.546 K,
collector 917.759 K, moderator 867.203 K, and reflector 808.874 K. Feedback-case maxima were
wall 875.492 K, fuel 2459.505 K, collector 898.084 K, moderator 850.731 K, and reflector
797.831 K; final power was 8.504 kW. The formal histories contain one extra 101 s sample due
to accumulated floating-point error at the 100 s schedule boundary; physical states and
restarts are unaffected. The runner now rounds phase elapsed time to microsecond resolution
before choosing the 0.1/1/10 s interval, with focused regression assertions at both boundaries.

### 2026-07-20 V14 orbital-state helium depressurization

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization_1/` starts from
`checkpoint_t019265s.npz` in the two-orbit fixed-power run. At accident time it sets all five
representative `collector_iclad_gap` gas conductivities to zero while retaining gap radiation,
NaK hydraulics, fluid-solid coupling, heat pipes, orbital external heat, TEC, and temperature
feedback. Fixed power is handed off to point kinetics. The two production runs use the LOCA
temperature limits and staged record schedule; the optional 5 s, -2 dollar scram also uses the
0.01 A persistent TEC open-circuit policy.

### 2026-07-22 V14 10 kW thermal-shield coupling

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/v14_heatpipe_radiator.py` can now attach the existing `RadiatorThermalShield` before the twelve `RingHP` components. The option remains disabled by default. With `thermal_shield_enabled=True`, the ordered 36 representative heat pipes map upper 18 -> shield sectors 0-5 and lower 18 -> sectors 6-11, with physical heat-pipe multipliers included in each sector's `T^4` average. The focused check is `python -m testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.test_v14_thermal_shield_coupling`.

### V14 heat-pipe transfer-failure accident runners

The three focused runners under `Full_Loop_Cases_10kW/` use the compatible
210 kW external-heat restart at `t=19864.7 s`. Their failure maps affect only
the evaporator fluid-solid coupling multiplier; nominal heat-pipe counts,
hydraulic loss parameters, flow areas, radiation areas, and orbital external
heat remain active. A multiplier of zero uses the exact zero-coupling path,
which removes the active heat-transfer source while retaining passive solid
conduction. TEC remains enabled and the thermal shield is disabled.

The cases are:

- `V14_210kW_heatpipe_partial_failure`: upper A5 local node 2 at 50% transfer.
- `V14_210kW_heatpipe_single_node_failure`: matching upper/lower A5 local node 2 at 0% transfer.
- `V14_210kW_heatpipe_sector_failure`: matching upper/lower A5 all three local nodes at 0% transfer.

Each runner holds fixed 210 kW until one solid temperature limit is reached,
then initializes point kinetics and applies -2 dollars. Coolant temperature
is diagnostic only; the channel-wall limit is separate. Default limits are
1058 K (channel wall), 2700 K (fuel), 1023 K (collector), 930 K (moderator),
and 1000 K (reflector). The default duration is one orbital period
(`5668.144369 s`), history tables are written every 1 s, periodic restart
files every 100 s, and accident-start/scram/final restart files are always
written. A scrammed case runs for at least another half orbit after scram.

Before long runs, use `testModule.test_fluid_solid_couple_multiplier`,
`testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.test_v14_thermal_shield_coupling`,
and the three runners with `--duration 2 --output-dir <temporary-directory>`.

With external heat enabled, the same V14 component now owns the atomic N6/N18 switch. N6 is the exact center-point subset of N18 columns `0,3,...,15`; both use `5668.14 s` and the same saved orbit origin. Shield-present steps apply `0.992*N6` only to the six side panels, with zero top/bottom heat and zero direct heat-pipe N18. `set_active(False)` reverses both sides of the switch in one pre-step. Global restart files save `active_override` and `orbit_time_origin_s`. The startup runner is `run_v14_shield_radiator_startup.py`; accepted Stage 0 results are under `Full_Loop_Cases_10kW/V14_210kW_start/phase_0_shielded_1800s_complete`. Stages 0-2 use a startup-only helium TEC gap with `h_eq=5678 W/m2/K`; other powered V14 cases retain the common cesium default `h_eq=29 W/m2/K`. Stage 0 runs zero power for `1800 s` from `300 K`, with TEC off, shield attached and total flow `0.615 kg/s`; the five CSV tables use a `10 s` interval and only the final NPZ is retained. Final coolant temperatures are `298.492-298.852 K`.

Stage 1 uses `Full_Loop_Cases_10kW/V14_210kW_start/run_v14_reactivity_startup.py`. It loads the Stage 0 restart, initializes point kinetics at `1 W`, holds `+0.50 $` without withdrawal, preserves the `0.615 kg/s` flow and TEC-off state, writes the five CSV tables every `1 s`, and stops at the first 10 kW crossing. The accepted helium-gap result in `stage_1_fixed_0p50_to_10kw` reached `10000.0009 W` after `141.45251 s` and contains only the final NPZ.

Stage 2 uses `Full_Loop_Cases_10kW/V14_210kW_start/run_v14_power_ramp.py`. It disables point kinetics, prescribes `600 W/s` from 10 to `70 kW`, and keeps TEC off. The five CSV tables use a `1 s` interval; restart files are written every `50 s` and at completion. The shield is removed when the minimum loop coolant reaches `373 K`, including an immediate start check, and the controlled total flow rises from `0.615` to `1.23 kg/s` when `CoreOutletConnector.T >= 500 K`. The accepted helium-gap `stage_2_power_ramp_10kw_to_70kw` run reached `70 kW` in `100.00000 s`; final minimum, core-outlet, and maximum coolant temperatures were `304.401 K`, `358.896 K`, and `379.636 K`, so neither transition fired.

### 2026-07-22 V14 10 kW prescribed-flow pump contract

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/common_flow_builder.py` 的两台串联泵仍各承担总设计压头的一半，但启用 `pump_flow_control` 时只有 `pump_a` 使用 `FlowControlledPumpJunction`；`pump_b` 必须保持普通 `PumpJunction`。闭式串联系统只施加一个总流量约束，不能在两台泵上重复施加，否则中间节点会出现幅值巨大、符号相反的泵压差。低电功率 fixed-I runner 只设置和验收实际具有 `set_flow_rate` 或 `target_W` 的受控泵。

### 2026-07-23 V14 20%电功率固定电流无外热工作点

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_20pct_electric_power_fixed_I/`
的当前工作点为堆芯热功率120.0 kW、总流量2.46 kg/s和TEC固定电流
213.4691467366893 A，轨道外热流关闭。正式调参采用零保持、零斜坡的直接设定，
不模拟降功率过程。

最终1500 s运行的末端300 s电功率为2039.144-2042.332 W，均值
2040.657 W，半极差0.0781%，线性漂移0.1557%；TEC和水力均收敛，未触发
温度安全限值。正式restart位于
`runs/noext_iter2_Q120000_1500s/final_restart.npz`，SHA256为
`21D39F835F4D62A42BE8437BF55147289610E1C5B4B3871F111BE74815BE825E`。
机器可读验收结果见`no_external_heat_summary.json`。

此前118.5 kW周期外热半周期验证未通过2.0 kW下限，保留在
`half_orbit_summary.json`中作为历史记录。当前工作点的验收范围仅为关闭外热流的
固定工况，不能据此声称已达到周期外热条件下的轨道周期稳态。

## 2026-08-03 V14 210 kW startup Stage 3 lookup and trial

`Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_start/run_v14_startup_stages.py` explicitly loads `core,startup,high_power,accident` from the default dense runtime v2 database. At 2802 s this changed lookup coverage from 0/2146 to 2146/2146 points. The 2800-2824 s verification run uses the series fixed-resistance no-generation fallback and records `tec_generating`; after a transient 0.606 A / 1.1 mW point at 2802 s, the tested interval returns finite zero electrical output without repeated failed outer iterations. The partial trial directory is retained under `startup_5000s_final/stage_3_stopped_24s_lookup_nogen_review`. Fixed-resistance TEC updates now return explicit zero current and voltage both for the early no-generation test and after exhausting 100 outer iterations; the wrapper also clears per-element electrical fields and Joule-power arrays before coupling, so failed trial state is not written into Stage 3 history or electrode heat sources. The focused regression is `test_thermocalc_interface.py::test_fixed_resistance_zero_output_clears_failed_node_state`.

## 2026-08-06 V14 210 kW startup fixed-current probes

`V14_210kW_start/run_v14_t2800_fixed_i_probe.py` evaluates a saved thermal state without advancing the coupled system. At the Stage 2 `t=2800 s` final restart, TEC surface temperatures were `TE=738.432-898.808 K` and `TC=738.238-863.219 K`; both the guarded and forced-C++ series `fixed_i=200 A` probes returned finite zero output with `converged=False`. At Stage 3 `checkpoint_t03700.000s.npz`, the same 58-element series probe used `TE=1063.459-2177.214 K` and `TC=749.900-876.671 K` and converged in one circuit iteration at `I=200 A`, `U=54.7510466 V`, and `Pe=10950.2093 W`. This demonstrates that the hot state can generate at 200 A; the continuing zero-output startup is caused by the `Rload=0.003 ohm` fixed-resistance strategy, not by insufficient high-state emitter temperature. A follow-up earliest-checkpoint probe found that `checkpoint_t02850.000s.npz` already converges at `I=200 A`, `U=2.8593556 V`, and `Pe=571.8711 W` with `TE=745.679-1723.969 K`; therefore the no-generation Stage 3 trajectory is already invalid before its first 50 s checkpoint and must be recomputed from the `t=2800 s` Stage 2 restart. A fixed-current U-I scan at the same temperature field also proves that the specified `Rload=0.003 ohm` has a physical load-line intersection near `I=218.7 A`: the TEC gives `U=0.6564785 V` while the load requires `0.6561 V`. The production fixed-resistance secant solver nevertheless failed to return within a 600 s diagnostic and the startup path fell back to zero output. The four configured wire resistances sum to `0.00352 ohm` per TEC and approximately `0.20416 ohm` over 58 series elements, 68 times the external load; the small external load places the root near the short-circuit knee and exposes the unbracketed fixed-resistance iteration weakness, but does not eliminate the physical operating point.
## 2026-08-06 V14 startup fixed-resistance solver acceptance

V14_210kW_start/run_v14_t2800_fixed_i_probe.py now also accepts
--mode fixed_r --resistance-ohm VALUE, so saved temperature fields can test
the production series load-line solver without advancing the coupled system.

With the rebuilt Python 3.12 extension, the Stage 2 t=2800 s restart returns
in about 1 ms through the low-temperature zero-emission guard:
TE=738.432-898.808 K, TC=738.238-863.219 K, finite zero output, and no C++
load-line iterations. At Stage 3 checkpoint_t02850.000s.npz,
Rload=0.003 ohm converges to I=218.703142 A, U=0.656109428 V, and
Pload=143.493193 W in 15 sampled current evaluations; all saved TEC arrays
are finite. The existing Stage 3 controller checks the result after each
coupled step and will switch immediately to strict fixed_i=216 A because
the accepted current exceeds the threshold. Full Stage 3 was not restarted
during this acceptance.
### 2026-08-06 V14 Stage 3 persistent fixed-R recovery

The first rebuilt Stage 3 run used max_dt=0.2 s and stopped at about 2804.9 s
after hydraulic NaN. A max_dt=0.05 s rerun reproduced the prior stable
hydraulic trajectory but remained at zero electrical output through 2877 s.
Its own checkpoint_t02850.000s.npz independently converged at Rload=0.003 ohm
to I=222.100115 A, U=0.666300345 V, and Pload=147.985383 W. This isolated the
problem to persistent circuit initialization after earlier zero-output
cleanup, not to the thermal state or load line. That run was stopped.

After adding fixed-R voltage-seed recovery and rebuilding the extension, the
formal restart from the original Stage 2 t=2800 s state is under
startup_5000s_fixed_r_recoveryfix_20260806 with max_dt=0.05 s. It retains
1 s five-CSV history and 50 s checkpoints. The old runs remain as diagnostics.
### 2026-08-06 V14 Stage 3 216 A event restart

The Stage 3 runner now saves checkpoint_tec_switch_tXXXXXXXXXs.npz immediately
after the fixed-R current first reaches 216 A and before rebuilding the circuit
as strict fixed_i=216 A. The event path is retained in the final summary as
tec_current_limit_switch_restart.

In startup_5000s_fixed_r_switchrestart_20260806, the persistent pre-generation
circuit was advanced to the regular checkpoint_t02850.000s.npz, then resumed
from that checkpoint to rebuild all TEC internal iteration state. The fixed-R
solve crossed the threshold at t=2850.030 s. The pre-switch event state is
checkpoint_tec_switch_t02850.030s.npz, SHA256
30450A32FE82BB9EF0411F85BD1A9F7227798CFE71898046131A99E40CA168F0.
The continuation then reported converged fixed_i=216 A; electrical power rose
from about 714.64 W at t=2851 s to 1557.98 W at t=2860 s.
### 2026-08-06 V14 Stage 3 completion at 5000 s

The supervised fixed-I continuation completed at absolute t=5000 s and wrote
final_restart.npz. Final values were I=216 A, U=41.659673 V, electrical power
8998.489 W, prescribed reactor power 210 kW, controlled flow 2.46 kg/s,
fluid temperature 718.864-844.184 K, and global solid temperature
582.484-2365.881 K. TEC convergence remained true. The final restart SHA256
is 1AF2A054BCECDAB604B4DB86B78FA3AD68A4B0CDF3E6833FEEE49039DA16749F.

The five history CSV files cover 2800-5000 s at 1 s summary intervals, and
50 s checkpoints are retained. This passes the requested Stage 3 control and
TEC workflow. It is not a strict hydraulic-convergence acceptance: the runner
uses fluid_max_iter=1 with fail_on_fluid_nonconvergence=False, and stderr
continues to contain finite-residual hydraulic non-convergence warnings.
### 2026-08-06 V14 Stage 3 fixed-I hold to 10000 s

V14_210kW_start/run_v14_stage3_hold_5000_to_10000.py continues from the
accepted 5000 s final restart while preserving 210 kW prescribed power,
2.46 kg/s flow, cesium gap heat transfer, enabled TEC, and strict fixed_i=216 A.
It writes a separate five-CSV history every 1 s and checkpoints every 50 s
under startup_10000s_fixed_i_continuation_20260806/. The initial continuation
check reached 5007 s with converged TEC and about 9.002 kW electrical power.
### 2026-08-07 V14 Stage 3 fixed-I hold completion at 10000 s

The 5000-10000 s fixed-I continuation completed and wrote final_restart.npz,
SHA256 CE6730645C695E030211C061CC92D14ACBD265A3E762A34E8834B16356D20340.
At 10000 s, TEC remained converged at 216 A and 42.443351 V, producing
9167.764 W. This is 169.274 W (1.88%) above the accepted 5000 s result.
Final coolant temperatures were 729.096-862.876 K and global solid
temperatures were 694.570-2368.575 K.

Over the final 100 s, coolant min/max still rose at about
0.00385/0.00255 K/s, global solid min/max at 0.00484/0.000353 K/s, and
electrical power at 0.0223 W/s. The state is closer to a plateau but is not a
strict thermal steady state. Hydraulic single-iteration warnings remain under
the same non-fail-closed runner policy.

### 2026-08-15 V14 20% electric-power endpoint calibration

The formal no-external-heat 20% endpoint restart was advanced for one full
external-heat period at fixed current 213.4691467366893 A and 2.46 kg/s. The
candidate thermal powers 120.5, 121.0, and 121.5 kW produced terminal electric
powers 2096.8718, 2153.1970, and 2210.5385 W, respectively, against the exact
20% target 2163.0637 W. All three runs completed with converged hydraulics and
TEC; linear interpolation gives the working endpoint setpoint
121086.033950196 W. The new complete trajectory uses
`V14_210kW_low_electric_power_fixed_I/runs/continuation_40pct_period10s_then_40to20_record1s_Q121086_20260815/`;
the older 120 kW trajectory is retained as a comparison only.

The resulting trajectory completed all three stages. The 40% hold used 10 s
history records, the 40% to 20% descent used 1 s records, and the following
20% cooling period used 10 s records; restart checkpoints were written every
60 s, 60 s, and 100 s respectively. After the final cooling period the
121086.033950196 W setpoint produced 2137.8100 W electrical output (19.7665%
of the 10815.3183 W full-power baseline), so the complete-trajectory endpoint
is 25.2536 W below the strict 20% target even though the isolated endpoint
calibration interpolated to that target.
