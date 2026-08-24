# ThermoCalc Codex 快速接管手册

> 更新时间：2026-06-25
> 适用范围：`ThermoCalc/` 模块及其必要的上层集成点。
> 使用方式：后续 Codex 首次接管 ThermoCalc 时先读本文件；只有在修改具体功能时，才按“修改场景索引”继续回查源码。

## 1. 先看结论

`ThermoCalc/` 是 TASTIN 中的热离子能量转换（TEC）求解模块。Python 上层向它提供多根 TFE 的轴向温度、几何、电路参数；C++ 核心求解串联或并联电路以及每根 TFE 内各轴向节点的发射、电流密度和电势分布。

当前源码的分层调用链是：

```text
Components/ReactorCore.py 或 Components/TECCircuitManager.py
    -> ThermoCalc/ThermoCalcWrapper.py
    -> ThermoCalc/bindings.cpp
    -> circuitTECs
    -> singleThermionicEnergyConversion
    -> thermionicEmission
```

2026-06-01 已闭合逐节点侧面积、`phiE/phiC/Vd` 结果读取和铯池温度运行时热更新，并使用 Python 3.12 重新构建 `te_solver.cp312-win_amd64.pyd`。串联公共模式支持 `fixed_u/fixed_r/fixed_i`；2026-06-25 起新增每根虚拟 TEC 一支路的并联模式 `parallel_fixed_u`、`parallel_fixed_i` 和 `parallel_load_curve`。`ReactorCore` 负责上层分组：V11/V13 默认主串联 34 根 TEC，可选把 `Ring3_Open` 代表的 3 根 TEC 单独接入预留并联电路。

2026-06-23 新增热离子发射查表加速实验路径。该路径保持原解析 `thermionicEmission::calc()` 可用，查表仅在显式启用时作为 `calc()` 的优先分支；表缺失或关闭时继续走原解析法。2026-06-25 并联验证通过后，根目录生产 `ThermoCalc/te_solver.cp312-win_amd64.pyd` 已替换为当前 `build_cp312/Release` 产物，旧版备份为 `ThermoCalc/te_solver.cp312-win_amd64.before_parallel_20260625.pyd`。

## 2. 首次接管阅读顺序

1. 先读本文件，建立当前源码基线。
2. 需要改 Python 公共接口或定位运行时故障时，读 [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py) 和 [`bindings.cpp`](./bindings.cpp)。
3. 需要改电路模式或全局迭代时，读 [`circuitTECs.h`](./circuitTECs.h) 和 [`circuitTECs.cpp`](./circuitTECs.cpp)。
4. 需要改单根 TFE 的轴向离散、电势、焦耳热相关输出时，读 [`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h) 和 [`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp)。
5. 需要改热离子发射物理时，读 [`thermionicEmission.h`](./thermionicEmission.h) 和 [`thermionicEmission.cpp`](./thermionicEmission.cpp)。
6. 需要判断 TASTIN 如何调用 ThermoCalc 时，再读 [`../Components/TECCircuitManager.py`](../Components/TECCircuitManager.py) 和 [`../Components/ReactorCore.py`](../Components/ReactorCore.py)。

[`THERMOCALC_ANALYSIS.md`](./THERMOCALC_ANALYSIS.md) 和 [`NONUNIFORM_GRID_GUIDE.md`](./NONUNIFORM_GRID_GUIDE.md) 保留为历史或专题资料。它们不是当前源码事实基准。

## 3. 文件清单

| 文件 | 职责 |
|---|---|
| [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py) | Python 公共包装层：准备输入、构建电路、调用计算、提取结果、提供更新方法 |
| [`bindings.cpp`](./bindings.cpp) | pybind11 边界：定义 `CalculationMode`、`InputData`、工厂函数和 Python 可见属性 |
| [`emissionLookup.h`](./emissionLookup.h)、[`emissionLookup.cpp`](./emissionLookup.cpp) | 可选热离子发射查表后端：管理分块表、四维插值、安全标志和全局启停 |
| [`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 多根 TFE 串联/并联电路：串联定电压/定电阻/定电流，并联定电压/定总电流/外部 U-I 负载曲线 |
| [`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h)、[`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp) | 单根 TFE：轴向节点、导线与电极电势、电流密度和电阻率 |
| [`thermionicEmission.h`](./thermionicEmission.h)、[`thermionicEmission.cpp`](./thermionicEmission.cpp) | 单个轴向节点的热离子发射 J-V 模型：阻塞、过渡、饱和分支 |
| [`NonLinerSolver.h`](./NonLinerSolver.h)、[`NonLinerSolver.cpp`](./NonLinerSolver.cpp) | C++ 辅助求解器；被构建并由头文件引用 |
| [`CMakeLists.txt`](./CMakeLists.txt) | C++17 + pybind11 构建配置，目标模块名为 `te_solver`，当前测试版构建包含 `emissionLookup.cpp` |
| [`tools/emission_database.py`](./tools/emission_database.py) | 离线热离子数据库计划、分块生成、统计、风险标志和优化表生成工具 |
| [`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md) | 热离子发射相图扫描、数据库和查表验证的专题说明 |
| [`test_real_case.py`](./test_real_case.py)、[`test_real_case_v2.py`](./test_real_case_v2.py)、[`test_real_case_v3.py`](./test_real_case_v3.py)、[`test_real_case_v4.py`](./test_real_case_v4.py) | ThermoCalc 独立脚本；运行前先确认扩展与解释器版本匹配 |

目录中当前存在 `te_solver.cp312-win_amd64.pyd` 及多个历史 `.pyd` 变体。不要仅凭文件存在就认为当前解释器可以加载它们。

## 4. C++ 分层职责

### 4.1 `circuitTECs`

`circuitTECs` 保存 `vector<singleThermionicEnergyConversion*> TECs`，负责多根 TFE 的串联和并联电路求解。

当前绑定可触发的顶层入口是 `circuitTECsCalc()`，Python 名为 `calc()`：

```text
isFixedU == true -> uFixedCircuitCalc()
isFixedR == true -> resistanceFixedCircuitCalc()
isParallelFixedU == true -> parallelUFixedCircuitCalc()
isParallelFixedI == true -> parallelIFixedCircuitCalc()
isParallelLoadCurve == true -> parallelLoadCurveCircuitCalc()
```

串联两个模式都会通过 `circuitCalc(I)` 在给定总电流下求串联电路，再用弦割式迭代修正电流。并联第一版采用“每个虚拟 TEC 独立接到同一母线”的拓扑，`parallelCircuitCalc(Ubus)` 在给定母线电压下逐支路求解并汇总 `Iout=sum(I_branch)`；并联定总电流和外部负载曲线模式再对 `Ubus` 做一维求根。

`parallel_load_curve` 的外部负载格式是 `U_load=f(I_total)`：Python 通过 `ThermoCalcModel.set_load_curve(current_a, voltage_v)` 传入严格递增的总电流轴和对应电压轴，C++ 线性插值求负载电压。若仅调用 `setup_circuit_mode("parallel_load_curve", R)` 而未显式提供曲线，包装层会生成一条线性欧姆曲线 `U=R*I`，便于与定电阻工况对比。

### 4.1.1 2026-06-25 parallel-circuit integration note

The current production `ThermoCalc/te_solver.cp312-win_amd64.pyd` was rebuilt from the C++ sources in this directory after adding the parallel circuit APIs. The previous root pyd was kept locally as `ThermoCalc/te_solver.cp312-win_amd64.before_parallel_20260625.pyd`.

Validated public modes:

| Python mode | C++ path | Intended upper-level use |
| --- | --- | --- |
| `fixed_u` | series fixed voltage | Default main V11/V13 circuit |
| `fixed_r` | series fixed resistance | Legacy series support |
| `fixed_i` | series fixed current with open-circuit fallback | Prescribed-current generation |
| `parallel_fixed_u` | each virtual TEC connected to a common bus voltage | Reserved Ring3_Open parallel circuit |
| `parallel_fixed_i` | solve common bus voltage for target total current | Reserved Ring3_Open fixed-current circuit |
| `parallel_load_curve` | solve common bus voltage against `U_load=f(I_total)` | Reserved Ring3_Open external load curve |

Series `fixed_i` reuses `circuitCalc(Itarget)`. A finite positive generated voltage is accepted; otherwise the solver reinitializes and evaluates `circuitCalc(0)`, returns zero current and the open-circuit voltage, and reports `converged=false`. If both solves fail, it returns finite zero current/voltage and emits a warning. `ReactorCore` should not switch all V11/V13 TECs to global parallel. It creates one main series `ThermoCalcModel` for `Center/Ring1/Ring2/Ring3_TEC` and, only when requested, one separate reserved parallel `ThermoCalcModel` for `Ring3_Open`.

2026-07-13 series fixed-current verification used an isolated extension under `ThermoCalc/build_series_fixed_i_test/Release`; the root production pyd was not replaced. The focused check is `testModule/test_thermocalc_series_fixed_current.py`, including single/multiple-TEC zero-current open circuit, a fixed-U cross-check operating point, a finite excessive-current fallback to positive open-circuit voltage, and the finite zero-output guard when both solves fail.

Verification commands used for this change:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_series_fixed_current.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_parallel.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_interface.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_lookup.py
```

### 4.2 `singleThermionicEnergyConversion`

每个实例表示一根轴向离散为 `n_node` 个节点的 TFE。构造函数接收 `input[13]`：

| 槽位 | 内容 | 当前 C++ 维度 |
|---|---|---|
| `[0]` | `Temitter` 发射极温度 `[K]` | `n_node` |
| `[1]` | `Tcollector` 接收极温度 `[K]` | `n_node` |
| `[2]` | `dlE` 发射极节点长度 `[m]` | `n_node` |
| `[3]` | `dlC` 接收极节点长度 `[m]` | `n_node` |
| `[4]` | `{crossAreaE, crossAreaC}` 电极横截面积 `[m^2]` | `2` |
| `[5]` | `sideAreaE` 发射极侧面积 `[m^2]` | `n_node` |
| `[6]` | `sideAreaC` 接收极侧面积 `[m^2]` | `n_node` |
| `[7]` | `resistanceWire` 导线电阻 `[Ohm]` | `4` |
| `[8]` | `{U_init, d_gap}` 初始总电压 `[V]`、极间距 `[mm]` | `2` |
| `[9]` | `Tcs` 铯池温度 `[K]` | `n_node` |
| `[10]` | `V_init` 极板电势差初值 `[V]` | `n_node` |
| `[11]` | `{Itarget}` 目标电流 `[A]` | `1` |
| `[12]` | `wireU` 导线电压初值 `[V]` | `4` |

`initial()` 初始化电阻率、电阻和 `thermionicUnits`。`Icalc()` 反复执行 `Jcalc()`、`UwireCalc()` 和 `VcalcFVM()`，更新节点电流密度、电势与内部截面电流。`sideAreaE` 和 `sideAreaC` 当前都是长度为 `n_node` 的向量，所有编译范围内的电流积分和电势公式均按节点索引使用面积。

### 4.3 `thermionicEmission`

每个实例表示一个轴向节点。输入包括发射极温度 `TE`、接收极温度 `TC`、铯池温度 `Tcs`、极间距 `d` 和极板电势差 `Vo`。`calc()` 先计算阻塞区结果，再按势垒条件进入过渡区和饱和区逻辑，输出节点电流密度 `J [A/cm^2]`。功函数 `phiE`、`phiC` 和电弧降 `Vd` 存在于这一最底层对象中。

2026-06-23 起，`calc()` 顶部会先检查全局查表开关。若查表启用且 `queryEmissionLookup(TE, TC, Vo, Tcs, d)` 命中安全插值点，则直接写回 `J/Vd/delta_V/phiE/phiC` 并返回；否则继续执行原解析逻辑。`calcDiagnostics()` 是诊断入口，用于离线生成表和记录 `regime/converged/iteration_count` 等元数据，不是生产电路主路径。

### 4.4 `emissionLookup`

`emissionLookup` 是进程内单例式查表仓库。Python 通过 pybind11 调用 `add_emission_lookup_block()` 加载多个分块；查询时按 `priority` 和块覆盖范围选择表，并对 `TE/TC/Vo/Tcs` 做四维线性插值。所有参与插值的角点必须通过 `lookup_safe` 标志，否则该块不命中。轴长度为 1 的退化维度允许存在，用于最后一段切片或固定维度块。

当前暴露给 Python 的辅助函数包括：

```text
clear_emission_lookup()
set_emission_lookup_enabled(enabled)
is_emission_lookup_enabled()
emission_lookup_block_count()
add_emission_lookup_block(...)
lookup_emission_point(TE, TC, Vo, Tcs, d_gap=0.5)
lookup_emission_points(TE_array, TC_array, Vo_array, Tcs_array, d_gap=0.5)
calc_emission_point(...)
calc_emission_point_production(...)
```

`calc_emission_point()` 调用诊断解析路径；`calc_emission_point_production()` 走生产 `thermionicEmission::calc()`，可用于检查查表是否已经接入生产分支。

## 5. Python 公共入口

`ThermoCalcModel` 位于 [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)。

| 入口 | 行为 | 当前注意事项 |
|---|---|---|
| `ThermoCalcModel(n_elements, n_nodes, lookup_db=None, enable_lookup=None, lookup_regions=None)` | 创建输入并可显式控制 lookup | 默认优先发现本地 `emission_runtime_db_v2/pcs_0p02_5torr`；miss 或无库时使用解析法 |
| `set_temperatures(T_em, T_co)` | 校验 `(N_elem, n_node)`，保存副本；电路已构建时写入每根 `SingleTEC` | 当前绑定暴露了 `Temitter`、`Tcollector` |
| `setup_circuit_mode(mode_str, target_value, I_guess=150.0)` | 接受 `fixed_R`、`fixed_U`、`fixed_I`、`parallel_fixed_u`、`parallel_fixed_i`、`parallel_load_curve`；`fixed_I` 要求有限非负目标电流 | `parallel_load_curve` 推荐配合 `set_load_curve()` 使用 |
| `set_load_curve(current_a, voltage_v)` | 设置并联外部负载 `U_load=f(I_total)` 曲线 | 电流轴必须一维、有限且严格递增 |
| `build()` | 把温度写入 `InputData`，调用 `te_solver.create_circuit()` | 绑定层在 `unchecked<>` 前执行完整形状校验 |
| `calculate(verbose=False)` | 必要时自动 `build()`，再调用 `_circuit.calc()`；返回耗时 `[ms]` | 不是物理时间步长度 |
| `load_emission_lookup_database(db_dir, enable=True, force=False)` | 从 `manifest.json`、`chunk_plan.json` 和 chunk `.npz` 加载查表数据库到 C++ 单例 | 优先加载同名 `.optimized.npz`；2026-06-25 起根目录生产 pyd 已暴露该 API |
| `get_global_results()` | 返回 `Iout`、`Uout`、`Rload`，并附带 `mode/converged/iteration_count/branch_currents/branch_voltages/effective_rload` | 电路未构建时返回 `None` |
| `get_tec_results(idx)` | 返回指定 TFE 的详细字典 | 包装层读取了绑定未暴露的字段，见风险清单 |
| `set_tcs(tcs_val)` | 接受标量或 `(N_elem, n_node)`，更新 `_input_data.Tcs` | 构建后调用电路级 `set_tcs()` 逐 TEC 热更新 |
| `set_rload(rload_val)` | 尝试更新输入结构；构建后尝试更新 C++ 电路负载 | 构建后的 `CircuitTECs.Rload` 属性已绑定；构建前输入结构没有对应可写绑定 |

## 6. `InputData` 字段、维度与默认值

`bindings.cpp` 定义 `InputData`，`ThermoCalcWrapper.py` 负责填充。`create_circuit()` 在进入 `unchecked<1>()` 和 `unchecked<2>()` 前对所有数组执行显式形状校验。

包装层支持通过环境变量选择测试扩展和查表数据库：

```powershell
$env:THERMOCALC_PYD_DIR = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\build_cp312\Release"
$env:THERMOCALC_ENABLE_LOOKUP = "1"
$env:THERMOCALC_LOOKUP_DB = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\emission_database"
```

`THERMOCALC_PYD_DIR` 会被插入到 `sys.path` 的最高优先级。查表控制优先级为显式参数、环境变量、本地推荐库；`enable_lookup=False` 或 `THERMOCALC_ENABLE_LOOKUP=0` 可强制解析法。

| 字段 | 绑定层期望维度 | 单位 | 包装层默认值 |
|---|---|---|---|
| `N_elements` | 标量 | - | `n_elements` |
| `n_axi` | 标量 | - | `n_nodes`；当前主要作为记录值 |
| `Temitter` | `(N_elem, n_node)` | `K` | `build()` 时写入内部 `600.0` 初值或最新温度 |
| `Tcollector` | `(N_elem, n_node)` | `K` | `build()` 时写入内部 `600.0` 初值或最新温度 |
| `dlE`、`dlC` | `(N_elem, n_node)` | `m` | 每节点 `0.507 / n_node` |
| `Tcs` | `(N_elem, n_node)` | `K` | `600.0` |
| `V_init` | `(N_elem, n_node)` | `V` | `0.2` |
| `crossAreaE` | `(N_elem,)` | `m^2` | `6.667e-5` |
| `crossAreaC` | `(N_elem,)` | `m^2` | `1.0786e-4` |
| `sideAreaE` | `(N_elem, n_node)` | `m^2` | 每格 `0.00092855424159680002 * 25 / n_node` |
| `sideAreaC` | `(N_elem, n_node)` | `m^2` | 每格 `0.00097592945800480005 * 25 / n_node` |
| `U_init` | `(N_elem,)` | `V` | `1.6` |
| `d_gap` | `(N_elem,)` | `mm` | `0.5` |
| `Itarget` | `(N_elem,)` | `A` | `200.0` |
| `resistanceWire` | `(N_elem, 4)` | `Ohm` | 全零 |
| `wireU` | `(N_elem, 4)` | `V` | 每根 `[0.8, 0.8, 0.0, 0.0]` |
| `mode` | 枚举 | - | `CalculationMode.FixedVoltage` |
| `target_val` | 标量 | 模式相关 | `0.89 * 34 = 30.26`，默认作为目标电压 |
| `I_total_init` | 标量 | `A` | `284.0` |
| `loadCurveCurrent` | `(n_curve,)` | `A` | 默认 `[0, 1000]`，仅 `parallel_load_curve` 使用 |
| `loadCurveVoltage` | `(n_curve,)` | `V` | 默认 `[0, 100]`，仅 `parallel_load_curve` 使用 |
| `R_load_init` | 标量 | `Ohm` | 仅在 C++ 结构体中声明；当前未绑定、未由包装层赋值、未参与工厂函数逻辑 |

## 7. 结果提取

### 7.1 已由当前绑定暴露

`CircuitTECs` 暴露：

| 字段 | 含义 |
|---|---|
| `TECs` | 单根 TFE 对象列表 |
| `Utarget`、`Rload`、`Itarget` | 目标电压、负载电阻、目标总电流 |
| `Iout`、`Uout` | 总电流、总电压 |
| `isFixedU`、`isFixedR`、`isParallelFixedU`、`isParallelFixedI`、`isParallelLoadCurve` | 模式标志 |
| `converged`、`iterationCount` | 顶层电路迭代状态 |
| `branchCurrents`、`branchVoltages` | 并联模式的各支路电流和电压 |
| `set_load_curve(current, voltage)` | 运行时更新并联外部 U-I 负载曲线 |
| `calc()` | 顶层计算入口 |

`SingleTEC` 暴露：

| 字段 | 维度 | 含义 |
|---|---|---|
| `I`、`U` | 标量 | 单根 TFE 电流 `[A]`、电压 `[V]` |
| `J` | `n_node` | 电流密度 `[A/cm^2]` |
| `V` | `n_node` | 极板电势差 `[V]` |
| `UE`、`UC` | `n_node` | 发射极、接收极电势 `[V]` |
| `rhoE`、`rhoC` | `n_node` | 发射极、接收极电阻率 `[Ohm*m]` |
| `IEsecSingle`、`ICsecSingle` | `n_node` | 发射极、接收极内部截面电流 `[A]` |
| `Temitter`、`Tcollector`、`Tcs` | `n_node` | 温度输入，可由 Python 读写 |
| `phiE`、`phiC`、`Vd` | `n_node` | 节点功函数和电弧降结果 |

### 7.2 包装层声明的 `get_tec_results()`

包装层当前尝试返回：

```text
I, U, J, V, UE, UC, rhoE, rhoC,
IEsecSingle, ICsecSingle,
phiE, phiC, Vd,
TE, TC
```

其中 `TE` 和 `TC` 分别从绑定后的 `Temitter`、`Tcollector` 读取。`Jcalc()` 会把节点级 `thermionicEmission` 的 `phiE`、`phiC`、`Vd` 同步到 `SingleTEC` 向量，并由绑定层向 Python 暴露。

## 8. 两条上层集成路径

### 8.1 `TECCircuitManager`：显式电热耦合路径

[`../Components/TECCircuitManager.py`](../Components/TECCircuitManager.py) 的 `pre_step()` 调用 `sync_thermo_electric()`：

```text
可选宏观参数更新 Tcs / R_load
    -> 从各 TECPair 收集间隙表面温度
    -> ThermoCalcModel.set_temperatures()
    -> ThermoCalcModel.calculate()
    -> 逐根 get_tec_results()
    -> 由 UE / UC / rhoE / rhoC 计算并下发焦耳热相关源项
    -> 由 J / phiE / TE / (UE - UC) 计算并下发等离子体热流
```

这是纯显式算子分裂路径。它依赖 `get_tec_results()` 中的 `phiE`，因此也受当前绑定缺口影响。

### 8.2 `ReactorCore`：直接耦合路径

[`../Components/ReactorCore.py`](../Components/ReactorCore.py) 直接持有 `ThermoCalcModel`：

```text
_build_thermo_calc()
    -> 按 tec_multipliers 创建虚拟 TEC 元件
    -> _configure_thermo_calc_geometry()

setup_tec_circuit()
    -> 主 TEC 电路默认沿用 fixed_u/fixed_r 串联模式

setup_reserved_parallel_tec_circuit()
    -> Ring3_Open 预留 TEC 映射到 parallel_fixed_u / parallel_fixed_i / parallel_load_curve

pre_step()
    -> 按 thermo_update_interval 调用 calculate()
    -> 从代表性虚拟元件读取结果
    -> 计算电场、焦耳热相关源项和等离子体热流
    -> 按 tec_mult / thermal_mult 缩放并下发到物理 TFE

post_step()
    -> 收集收敛后的发射极、接收极表面温度
    -> set_temperatures() 写回虚拟元件矩阵
```

`_configure_thermo_calc_geometry()` 会从热网格复制逐节点 `dlE/dlC` 和逐节点 `sideAreaE/sideAreaC`。当前 C++ 绑定和单 TEC 后端已支持这一输入。

## 9. 事实等级

### A. 已由当前源码确认

- 模块调用链为 `ReactorCore / TECCircuitManager -> ThermoCalcWrapper -> bindings.cpp -> circuitTECs -> singleThermionicEnergyConversion -> thermionicEmission`。
- 当前绑定枚举包含串联 `FixedVoltage/FixedResistance`，以及并联 `ParallelFixedVoltage/ParallelFixedCurrent/ParallelLoadCurve`。
- 当前绑定层把 `sideAreaE/sideAreaC` 作为 `(N_elem, n_axi)` 二维数组逐行送入单根 TFE。
- 当前单根 TFE 的 `sideAreaE/sideAreaC` 是逐节点 `vector<double>`。
- `dlE/dlC` 在单根 TFE 内是向量，当前 `VcalcFVM()` 会读取节点长度。
- 当前 `SingleTEC` 绑定暴露 `phiE`、`phiC`、`Vd`。
- 当前 `CircuitTECs` 绑定提供 `set_tcs()`，并逐根更新 `SingleTEC.Tcs`。
- 当前 `CircuitTECs` 绑定提供 `set_load_curve()`、`branchCurrents/branchVoltages` 和电路收敛状态；并联第一版拓扑为每根虚拟 TEC 一个独立支路。
- 上层两条路径都依赖 `ThermoCalcModel`，`ReactorCore` 还会写入二维逐节点侧面积。
- 测试版 `te_solver` 暴露热离子查表 API；`ThermoCalcWrapper.py` 可通过环境变量加载 `ThermoCalc/emission_database`。
- 查表数据库的全量计划、优化表和验证摘要见 [`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md)。

### B. 需要继续验证

- 更长时 TEC 瞬态下，电子边界功率差与端功率、焦耳热之间的稳定离散误差门槛。
- `set_rload()` 在目标运行环境中的构建前、构建后行为。
- 定电压、定电阻模式在目标案例中的数值结果和收敛性。
- 查表插值与原解析法的全工况误差门槛；V13 1 s 短测中 TEC 端电功率约有 `1.14%` 差异，热工量差异约 `1e-5` 相对量级。

### C. 历史文档中的旧结论

- [`NONUNIFORM_GRID_GUIDE.md`](./NONUNIFORM_GRID_GUIDE.md) 声称逐节点侧面积和任意非均匀网格支持已经完成。
- [`THERMOCALC_ANALYSIS.md`](./THERMOCALC_ANALYSIS.md) 以及本手册旧版本曾把 `fixed_I`、`phiE/phiC/Vd` 提取和若干热更新路径描述为可用能力。

这些旧结论不能覆盖当前源码。修改前必须重新核对对应实现，并使用匹配解释器的扩展做运行验证。

## 10. 已知风险

| 风险 | 当前源码证据 | 处理原则 |
|---|---|---|
| 默认解释器与扩展 ABI 不匹配 | 默认 `python --version` 仍为 `3.9.13`，主扩展为 `te_solver.cp312-win_amd64.pyd` | TEC 测试和运行使用 Python 3.12 |
| 串联 `fixed_i` 不允许外加强迫耗电 | 目标电流求解仅接受有限正端电压；否则回退到零电流开路状态并标记未收敛 | 仍需区分串联 `fixed_i` 与并联定总电流 `parallel_fixed_i` |
| 并联拓扑第一版较简单 | `parallel_*` 采用每根虚拟 TEC 一支路，暂不支持“支路内部多根串联再并联”的 branch groups | 复杂接线需要新增分组输入和上层映射，不能直接复用当前语义 |
| 非均匀网格完整数学验证仍有限 | `dlE/dlC` 和侧面积已逐节点化，但 `VcalcFVM()` 的界面距离仍沿用现有离散公式 | 修改电势离散时补专项守恒验证 |
| `set_rload()` 构建前行为不完整 | 包装层尝试写 `_input_data.Rload/R_load`，当前 `InputData` 未绑定这些字段 | 构建后可写 `CircuitTECs.Rload`；构建前优先使用 `setup_circuit_mode('fixed_R', ...)` |
| 原始指针生命周期风险 | `create_circuit()` 为每根 TFE `new` 对象；相关析构函数当前为空 | 若处理长时运行内存问题，专项检查所有权与释放逻辑 |
| 查表默认优先、miss 回退解析法 | 根目录生产 pyd 已包含当前查表接口；数据库可由 `ThermoCalcModel` 显式参数或环境变量加载 | 算例层优先使用显式 `tec_lookup_*` 配置；旧 runner 仍可设置 `THERMOCALC_ENABLE_LOOKUP=1` 和 `THERMOCALC_LOOKUP_DB` |
| 查表表外点会回退解析 | `thermionicEmission::calc()` 查表 miss 后继续原解析法 | 若希望完全避免解析失败输出，应扩大/修正表覆盖或改电路层策略 |
| setup 阶段首次 TEC 仍可能慢且打印失败 | V13 `apply_wire_resistance()` 会重建电路并立即 `calculate()`；30 s 查表计时中 setup 首算约 `7.81 s` 且打印失败信息 | 正式推进 warm-start 后 TEC 单次约 `1.8 s`；setup 首算需单独优化 |

## 11. 修改场景索引

| 修改场景 | 必读文件 | 重点检查 |
|---|---|---|
| 改 Python 公共 API、结果字段、热更新 | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp) | Python 名称、绑定字段、形状、构建前后行为必须成对闭合 |
| 增加或修复 `fixed_I` | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp)、[`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 枚举、`build()` 分支、顶层分发、目标电流语义 |
| 改并联电路模式或负载曲线 | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp)、[`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | `parallelCircuitCalc(Ubus)`、支路状态、`branchCurrents/Voltages`、`U_load=f(I_total)` 插值和求根收敛 |
| 暴露 `phiE/phiC/Vd` | [`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h)、[`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp)、[`bindings.cpp`](./bindings.cpp)、[`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py) | 节点级对象如何汇总成 `SingleTEC` 向量 |
| 修复逐节点侧面积或非均匀网格 | [`bindings.cpp`](./bindings.cpp)、[`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h)、[`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp)、[`../Components/ReactorCore.py`](../Components/ReactorCore.py) | 所有积分、电势方程、截面电流、界面距离和测试输入 |
| 改全局串联电路迭代 | [`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 定电压/定电阻分支、弦割迭代、`Iout/Uout/Rload/Utarget` |
| 改热离子发射物理 | [`thermionicEmission.h`](./thermionicEmission.h)、[`thermionicEmission.cpp`](./thermionicEmission.cpp) | 阻塞、过渡、饱和分支以及单位 |
| 改查表数据库或插值逻辑 | [`emissionLookup.h`](./emissionLookup.h)、[`emissionLookup.cpp`](./emissionLookup.cpp)、[`tools/emission_database.py`](./tools/emission_database.py)、[`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md) | 表覆盖范围、`lookup_safe`、退化轴、`.optimized.npz` 优先级、解析回退行为 |
| 改上层显式耦合 | [`../Components/TECCircuitManager.py`](../Components/TECCircuitManager.py) | gather、calculate、scatter、热流符号、电场计算 |
| 改 ReactorCore 直接耦合 | [`../Components/ReactorCore.py`](../Components/ReactorCore.py) | 虚拟元件映射、几何同步、更新频率、乘数缩放、post-step 温度回写 |
| 改构建配置或重编 `.pyd` | [`CMakeLists.txt`](./CMakeLists.txt)、[`bindings.cpp`](./bindings.cpp) | Python ABI、模块名 `te_solver`、pybind11 获取方式、C++17、MSVC 选项 |

## 12. 验证要求

修改 ThermoCalc 后至少执行以下检查：

1. 使用与目标 `.pyd` ABI 匹配的 Python 解释器导入 `te_solver`。
2. 核对 `CalculationMode`、`InputData`、`CircuitTECs`、`SingleTEC` 的运行时属性。
3. 对构建前后分别验证温度、`Tcs` 和负载更新。
4. 运行至少一个均匀网格独立案例。
5. 若修改面积或网格，增加逐节点面积、非均匀长度、电荷守恒和结果回归检查。
6. 若修改上层耦合，分别检查 `TECCircuitManager` 和 `ReactorCore` 路径，避免只修复其中一条。
7. 若修改查表路径，至少运行 `testModule/test_thermocalc_lookup.py`，再用 V13 restart 做 `1 s` smoke，确认 `tec_coupled_enabled=True` 且未出现 `disabling TEC coupling`。
8. 若修改并联电路，运行 `testModule/test_thermocalc_parallel.py`，确认并联定电压、定总电流、U-I 负载曲线和构建后温度/Tcs 更新均通过。

2026-06-01 已使用：

```text
C:\Users\HC Zhao\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe
Python 3.12.13
```

重新构建并验证 `te_solver.cp312-win_amd64.pyd`。`testModule/test_thermocalc_interface.py` 覆盖形状拒绝、均匀与非均匀节点面积、`phiE/phiC/Vd`、构建前后温度和 `Tcs` 更新、`fixed_I` 显式拒绝。单 TFE TEC `1 s` 基线也已运行。

2026-06-25 已使用仓库 Conda Python 3.12 环境重建 `ThermoCalc/build_cp312/Release/te_solver.cp312-win_amd64.pyd`，新增并联三模式绑定和 wrapper 接口；验证命令包括 `testModule/test_thermocalc_parallel.py`、`testModule/test_thermocalc_interface.py`、`testModule/test_thermocalc_lookup.py`，均通过。根目录生产 `.pyd` 已替换为该构建产物，旧版备份为 `ThermoCalc/te_solver.cp312-win_amd64.before_parallel_20260625.pyd`。

## 13. 2026-06-02 FVM 一致焦耳热输出

`singleThermionicEnergyConversion::VcalcFVM()` 当前直接输出：

```text
joulePowerE[n_axi] [W]
joulePowerC[n_axi] [W]
```

每个内部电阻面的耗散功率按 `0.5 / 0.5` 分配给相邻轴向热节点；两端半单元电阻的耗散功率全部分配给端点节点。绑定和 `ThermoCalcWrapper.get_tec_results()` 同时暴露四个端点电势，供审计重构面电导耗散。

生产耦合层必须使用 `joulePowerE/C` 下发焦耳热。`UE / UC / rhoE / rhoC` 和 Python 节点梯度函数继续保留为诊断数据，不得再次作为生产焦耳热权威值。

2026-06-02 已使用 Python `3.12.13` 重建本地扩展。单 TFE TEC `1 s` 审计中，二维映射与 C++ 节点功率总差为 `0 W`；TEC 转换闭合差由旧梯度口径的约 `0.280 W` 降为 `0.0250 W`。剩余量主要受当前外层电路电流停止条件影响，本轮未修改该阈值。

v7 CaseA `1 s` smoke 和静态接口审计可以完成，但底层 TEC 会重复报告 `Failed to converge after 100000 iterations.`。静态接口审计的即时 TEC 闭合差约为 `5.20 W`。2026-06-09 的 V8 CaseA `LSODA` smoke 中，该信息共出现 `4440` 次，分段探针确认全部发生在进入长算主循环前的 `core.thermo_calc.calculate(verbose=False)` 调用；随后 `system.step(0.01)`、记录输出和结束阶段未再出现。这不是二维热源映射误差，也不是 `BaseHeatConduction` 固体 ODE 收敛失败；后续应单独检查多 TEC 串联电路的收敛条件、失败状态传播和告警限流。

## 14. 2026-06-23 查表加速阶段总结

本阶段目标是把局部热离子发射函数

```text
f(TE, TC, Vo, phiE, phiC, d_gap, Tcs) -> J / Vd / delta_V / phiE / phiC
```

的解析迭代替换为可选查表路径，同时保留原解析法。已完成的源码改动包括：

- `thermionicEmission::calcDiagnostics()`：新增诊断计算入口，返回 `J/Vd/delta_V/phiE/phiC/regime/converged/iteration_count` 等元数据。
- `thermionicEmission::calc()`：新增查表优先分支；命中后直接返回，未命中时执行原解析法。
- `emissionLookup.*`：新增 C++ 查表仓库、四维线性插值、安全角点检查和全局启停。
- `bindings.cpp`：暴露查表 API、批量查表 API、诊断单点 API 和生产单点 API。
- `ThermoCalcWrapper.py`：支持 `THERMOCALC_PYD_DIR` 选择测试版 pyd，支持 `THERMOCALC_ENABLE_LOOKUP` / `THERMOCALC_LOOKUP_DB` 自动加载数据库。
- `CMakeLists.txt`：测试版构建纳入 `emissionLookup.cpp`。
- `testModule/test_thermocalc_lookup.py`：覆盖单块加载、优化块加载、生产 `calc()` 查表分支和查表/解析速度对比。

旧数据库位于 `ThermoCalc/emission_database/`，使用的是旧铯压范围。旧全量计划含 `18,737,388` 点、`78` 个 chunk；已完成 `startup` 和 `accident` 风险点优化表，`.optimized.npz` 会被包装层优先加载。优化后 `startup+accident` 中原始无效点 `55,506` 个，其中 `43,104` 个按零发射处理，`12,402` 个由邻域插补，未解决点为 `0`。该旧库已不再作为当前压力范围基准。

新的 `0.02-5.0 torr` 全量数据库已经生成在 `ThermoCalc/emission_database/pcs_0p02_5torr/`。它保留原 Pcs 点数和 log spacing 类型，唯一物理网格仍为 `18,737,388` 点，计划 chunk 数为 `76` 个。原始 chunk 文件内部包含相邻块共用的 TE 右边界平面，因此 `summarize --scan-chunks` 的 chunk 计数为 `25,756,400` 点；这不是物理网格变大。

当前新库已完成 `core/startup/high_power/accident` 全 region 优化：原始无效点 `655,530` 个，其中 `107,272` 个按安全零电流处理，`548,258` 个由邻域插补，未解决点为 `0`，优化后安全点覆盖率为 `1.0`。dense runtime v2 已导出到 `ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/`，包含 `.npz` 和 C++ 直接加载的 `.tedb`。

验证结果：

- `testModule/test_thermocalc_lookup.py` 通过；单点批量查表约 `1.53e6 points/s`，解析法约 `5.84e4 points/s`，局部函数约 `26x` 加速。
- V13 `1 s` 查表 smoke 可完成，`tec_coupled_enabled=True`；热工量与解析基准接近，TEC 端电功率高约 `62.2 W`，约 `1.14%`。
- V13 查表分段长算已从 `21000 s` 推进到 `22000 s`，最终 `tec_coupled_enabled=True`，未检出 `disabling TEC`、`Traceback` 或运行阶段解析收敛失败文本。
- V13 真实推进 `30 s` 计时：推进总耗时 `422.02 s`，TEC 计算 `51.24 s` 占 `12.16%`，导热 `252.70 s` 占 `59.98%`，流动 `3.41 s` 占 `0.81%`，其他系统层开销 `27.04%`。
- 查表 warm-start 后 TEC 单次更新约 `1.8 s`；setup 阶段导线电阻重建后的首次 TEC 约 `7.81 s`。

当前判断：查表路径已经把 TEC 从主要瓶颈之一降为次要瓶颈；V13 长算主耗时转移到导热求解和系统层调度。后续若继续优化速度，优先检查导热 solid 数量、辐射器管壁求解、coupler/组件调度开销，以及 `circuitTECs` 外层迭代次数，而不是继续只优化局部 `thermionicEmission` 单点。

## 15. 2026-06-23 runtime 查表压缩与索引

本轮在不删除原解析法的前提下，新增了运行时专用查表格式和 C++ 查询索引：

- `tools/emission_database.py export-runtime-dense` 从 `ThermoCalc/emission_database/pcs_0p02_5torr/` 导出 `ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/`，优先读取 `.optimized.npz`，输出 `runtime_dense_manifest.json`、按 region 分开的 `*.runtime.v2.npz` 和 C++ 直接加载的 `*.runtime.v2.tedb`。
- runtime 表只保留 `TE_axis/TC_axis/Vo_axis/Tcs_axis`、`J/Vd/delta_V/phiE/phiC`、`lookup_safe/zero_mask`；默认字段精度为 `float32`，但 `phiE/phiC` 保留，供边界条件继续调用。
- `zero_mask` 标记安全零电流区；启用 `--zero-compress` 后这些点的 `J` 在 runtime 表中直接写为 `0`，但电压和功函数字段仍参与插值。
- `ThermoCalcWrapper.load_emission_lookup_database(..., regions=...)` 同时支持旧全量库、legacy runtime 库和 dense runtime v2；默认只加载 `core`，可用环境变量 `THERMOCALC_LOOKUP_REGIONS=core,startup,high_power,accident` 扩展覆盖。
- `emissionLookup.*` 内部按 `region_id/priority` 建索引，对每块维护 bbox，用 TE chunk 直接定位候选块，并缓存上一次命中的块，减少全表线性扫描。

当前定向验证：

```text
cmake --build ThermoCalc\build_cp312 --config Release
python -m py_compile ThermoCalc\ThermoCalcWrapper.py ThermoCalc\tools\emission_database.py testModule\test_thermocalc_lookup.py
testModule/test_thermocalc_lookup.py
  passed
  runtime export/load path covered
  lookup batch: about 1.12e6 points/s in this run
  analytic local solver: about 5.60e4 points/s
  local speedup: about 20x
```

`ThermoCalc/emission_database/`、`ThermoCalc/emission_runtime_db/` 和 `ThermoCalc/emission_runtime_db_v2/` 都是生成数据，不应提交到 git。当前推荐运行路径是 `ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/`；使用全 region 覆盖时，设置 `THERMOCALC_LOOKUP_DB` 指向该目录，并显式设置 `THERMOCALC_LOOKUP_REGIONS=core,startup,high_power,accident`。

2026-06-23 后续修复：

- `chunk_te_ranges()` 现在为每个 TE chunk 保留右侧边界平面，避免新生成数据库出现 `1300-1310 K` 后直接跳到 `1320-1330 K` 的插值空隙。
- `export-runtime` 对旧数据库自动拼接下一 chunk 的第一个 TE 平面；新数据库在原始生成阶段已经保留右边界平面，dense runtime v2 导出时再去除重复拼接平面，唯一物理点数保持 `18,737,388`。
- 修复 `lookup_emission_points()` 输出数组 stride 为 `0` 的绑定问题；批量 API 现在与单点 `lookup_emission_point()` 一致，可以重新作为 benchmark 使用。
- 重新导出的 core runtime 表为 `43` 个 chunk、约 `15,276,928` runtime 点、`129.49 MB`；连续 core 随机采样 `200000/200000` 命中，批量查表约 `1.05e6 points/s`。
- `testModule/test_thermocalc_lookup.py` 已覆盖 runtime 右边界拼接、TE 空隙点命中、批量数组 stride 和批量/单点一致性；本轮局部 benchmark 为查表约 `3.72e6 points/s`、解析约 `9.44e4 points/s`、约 `39x`。

## 16. 热离子查表完整流程索引

热离子查表当前分为离线数据生成和运行时调用两条链，详细步骤维护在 [`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md) 的 `End-to-End Lookup Workflow`。

最小流程如下：

```text
离线生成:
  emission_database.py plan
  -> emission_database.py worker
  -> ThermoCalc/emission_database/chunks/*.npz
  -> summarize / verify / optimize-table
  -> export-runtime-dense
  -> ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/*.runtime.v2.tedb

运行调用:
  ThermoCalcModel.__init__()
  -> load_emission_lookup_database()
  -> te_solver.load_emission_dense_file()
  -> emissionLookup.cpp 内存索引
  -> thermionicEmission::calc()
  -> queryEmissionLookup()
  -> 命中则返回 J/Vd/delta_V/phiE/phiC
  -> 未命中则回退原解析 calc()
```

关键边界：

- `ThermoCalc/emission_database/pcs_0p02_5torr/` 是当前 `0.02-5.0 torr` 原始/审计库，保留诊断字段和 `.optimized.npz` sidecar，不提交 git。
- `ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/` 是当前推荐运行库，只保留 `J/Vd/delta_V/phiE/phiC/lookup_safe/zero_mask` 和轴，并提供 `.tedb` 给 C++ 直接加载，不提交 git。
- 默认优先自动发现 `ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/`；也可通过环境变量或算例层 `tec_lookup_*` 覆盖。默认只加载 `core`，其他 region miss 时回退解析法。
- 根目录生产 `.pyd` 已包含查表和并联接口；查表命中时直接返回，miss 时沿用 C++ 解析计算。

## 17. 2026-06-23 dense runtime v2 补充

当前推荐的运行时查表格式是 `export-runtime-dense` 生成的 dense runtime v2：

```text
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/
  runtime_dense_manifest.json
  core.runtime.v2.npz / core.runtime.v2.tedb
  startup.runtime.v2.npz / startup.runtime.v2.tedb
  high_power.runtime.v2.npz / high_power.runtime.v2.tedb
  accident.runtime.v2.npz / accident.runtime.v2.tedb
```

该格式按 region 存储一个连续四维张量，字段为 `J/Vd/delta_V/phiE/phiC`，并把 `lookup_safe` 和 `zero_mask` 压缩为 bit-packed mask。`.npz` 是可移植格式，`.tedb` 是 C++ 直接加载格式；包装层发现 `runtime_dense_manifest.json` 后会优先加载 `.tedb`，否则回退到 `.npz`。

2026-06-23 压力范围修正后的全量库将 `core/startup/high_power/accident` 的铯压范围改为 `0.02-5.0 torr`，但保留原 Pcs 点数和 log spacing 类型：`core/high_power=41`、`startup=21`、`accident=31`。这里 `Pcs` 单位明确为 torr，不是 Pa；换算沿用 C++ 生产模型公式 `Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)`；该范围对应 `Tcs ≈ 441.44-614.62 K`，覆盖当前 `Tcs=600 K` 算例。

当前 corrected dense runtime v2 汇总：

```text
total_points: 18,737,388
total_size_bytes: 537,914,570
zero_compress: true
zero_j_threshold: 1e-3

core        shape 86 x 41 x 71 x 41, points 10,264,186, NPZ 87,074,650 bytes, TEDB 207,851,764 bytes
startup     shape 31 x 31 x 36 x 21, points    726,516, NPZ  5,368,120 bytes, TEDB  14,712,989 bytes
high_power  shape 25 x 26 x 71 x 41, points  1,892,150, NPZ 18,317,488 bytes, TEDB  38,317,432 bytes
accident    shape 86 x 61 x 36 x 31, points  5,854,536, NPZ 47,715,973 bytes, TEDB 118,556,154 bytes
```

加载 smoke 已通过：`load_emission_lookup_database("ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr", regions=["core","startup","high_power","accident"], force=True)` 返回 `4`，`emission_lookup_dense_region_count()` 返回 `4`，`lookup_emission_point(1800, 800, 1.0, 600, 0.5)` 命中 `core`，`lookup_emission_point(1000, 650, 0.5, 600, 0.5)` 命中 `startup`。

`ThermoCalc/emission_runtime_db_v2/` 是生成数据，不提交 git。完整复现命令、字段说明和 v1/v2 对比维护在 [`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md)。

## 18. 2026-06-28 low-temperature zero-emission guard

`ThermoCalcModel.calculate()` now has a Python-layer fast guard before entering the C++ circuit iteration. If all emitter temperatures are below the zero-emission cutoff, the wrapper treats the TEC set as a non-generating open-circuit state and does not call `CircuitTECs.calc()`.

Default behavior:

```text
THERMOCALC_ZERO_EMISSION_TE_MAX_K = 1000.0
THERMOCALC_DISABLE_ZERO_EMISSION_GUARD unset/false
```

When the guard triggers:

- `Iout = 0.0`
- fixed-voltage modes keep `Uout = target_val`
- `J`, `UE`, `UC`, `V`, `phiE`, `phiC`, `Vd`, `joulePowerE`, and `joulePowerC` are written as zero arrays
- `converged = True`, `iteration_count = 0`
- `get_global_results()` reports `zero_emission_skipped=True` and a human-readable `zero_emission_reason`

The intent is to avoid meaningless low-temperature fixed-voltage TEC solves where the physical emission is effectively zero but the series circuit secant iteration attempts to force a current solution. This protection is deliberately conservative and can be disabled for diagnostics with:

```powershell
$env:THERMOCALC_DISABLE_ZERO_EMISSION_GUARD = "1"
```

Regression coverage is in `testModule/test_thermocalc_interface.py::test_low_temperature_fixed_voltage_auto_skips_zero_emission_case`.

### 18.1 2026-06-28 C++ iteration guard for uncertain TEC computability

The zero-emission guard is now backed by C++-level iteration protection in `circuitTECs.cpp`. This covers cases where the user does not know ahead of time whether a TEC state is computable and therefore wants to keep TEC enabled by default.

Protected paths:

- `initialSingleTECU()` and `singleTECU()` now return early when the fixed-current secant slope is non-finite or effectively zero.
- The per-TEC secant iteration cap was reduced from `1000` to `100` for these voltage-search loops.
- `uFixedCircuitCalc()` and `resistanceFixedCircuitCalc()` now detect non-finite or zero voltage-response slopes and return a finite non-converged zero-current state instead of propagating NaN or spending a long time in nested iterations.
- `circuitCalc()` no longer evaluates the series-voltage convergence ratio when `Uout` is zero or non-finite.

This means production runs can leave TEC coupling enabled during startup/transition exploration: clearly non-generating low-temperature states are skipped by the Python guard, while harder-to-predict pathological electrical states are stopped inside C++ and surfaced through `converged=False` / finite global results rather than hanging the process.

The production `ThermoCalc/te_solver.cp312-win_amd64.pyd` was rebuilt from the current C++ sources after this change. Verification included:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_interface.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_series_fixed_current.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_parallel.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_lookup.py
```

### 18.2 2026-06-29 series fixed-voltage bracketing update

A high-temperature V13 startup switch-point exposed a different failure mode from the low-temperature dead loop: the old series `fixed_u` outer iteration used two unbounded secant guesses (`Iout` and `Iout + 10 A`) and treated small current-step change as convergence without checking the voltage residual. At the V13 switch point this could run to the iteration cap, leave the last `circuitCalc(I)` state in `Uout/Iout`, and report finite but non-converged values far from the target voltage.

`circuitTECs::uFixedCircuitCalc()` has been changed in source to:

- sample a conservative set of non-negative current guesses around the requested initial current;
- seek a sign-change bracket for `Utarget - circuitCalc(I)`;
- use a guarded secant/bisection hybrid once a bracket is found;
- require voltage residual convergence before reporting `converged=True`;
- return the best finite sampled state with `converged=False` when no bracket exists.

Validation was performed with the rebuilt test extension in `ThermoCalc/build_cp312/Release` via `THERMOCALC_PYD_DIR`; the root production `.pyd` was not overwritten by this verification step.

Relevant checks:

```powershell
$env:THERMOCALC_PYD_DIR='...\ThermoCalc\build_cp312\Release'
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_interface.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_lookup.py
```

The V13 high-temperature fixed-U isolated switch-point solve changed from a non-converged `~430 s` full-step behavior to an isolated `~0.736 s` ThermoCalc solve with `Uout=27.2 V`, `Iout~=1.215 A`, `converged=True`, and `iteration_count=10`, when paired with the augmented local lookup table.

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
### 2026-06-29 V13 cold-start residual is not a ThermoCalc fixed-U residual

In the corrected V13 cold-start continuation using `h_eq = 29.0 W/m2/K`, fixed-R warmup successfully switched to fixed-U `27.2 V` and the subsequent fixed-U continuation reported finite outputs with `tec_solver_converged = true` and `tec_solver_iteration_count = 1` at the recorded checkpoints.

For the stopped long continuation:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629
history final time = 32490 s
TEC electric power ~= 5.693 kW
radiator heat rejection - coolant enthalpy rise ~= -0.136 W
core heat - coolant enthalpy rise - TEC electric power ~= 1.329 kW
```

Because the TEC solve converged and the radiator/coolant energy balance is sub-watt, the remaining kilowatt-scale term is a system transient storage term in the still-warming core/TFE/structural solids. Do not diagnose this specific residual as a ThermoCalc C++ fixed-U convergence failure unless future runs show non-finite TEC outputs, failed convergence flags, large circuit iteration counts, or discontinuous jumps in electrical power.

## 2026-07-07 UE/UC persistence in upper-level TFE state

`ThermoCalcModel.get_tec_results(idx)` already exposes `UE`, `UC`, `terminalPointUE1`, `terminalPointUE2`, `terminalPointUC1`, and `terminalPointUC2`. The upper-level `ReactorCore` path now stores those values into each `TFEUnit` restart state as raw electrode-potential diagnostics:

```text
.../electric/emitter_potential
.../electric/collector_potential
.../electric/emitter_collector_voltage_drop
.../electric/terminal_point_ue1
.../electric/terminal_point_ue2
.../electric/terminal_point_uc1
.../electric/terminal_point_uc2
```

This is a Python-side persistence change only. It does not alter the C++ TEC solve, Joule heat authority, plasma heat-flux formula, or circuit convergence logic.

## 2026-07-20 emission lookup TC extension

The production accident lookup envelope has been extended from TC 1100 K to
TC 1500 K while retaining 10 K spacing. The preferred local runtime directory
is now ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr_tc1500; the previous
pcs_0p02_5torr directory remains a compatibility fallback.

The full database contains 22,576,428 unique points. Optimization leaves zero
unresolved points, and the dense runtime v2 NPZ plus TEDB artifacts occupy
644,446,789 bytes. Exact-grid checks with the root production pyd passed at
TC 1100, 1110, 1300, and 1500 K; TC 1500.1 K correctly returns a lookup miss.

Within the 3,839,040 newly added TC 1110-1500 K points, 13,399 raw analytic
points were invalid: 4,447 were safely zero-filled and 8,952 were
neighbor-imputed. The failures are dominated by TE >= 1900 K, Vo < 1.0 V,
and Pcs 0.02-0.05 torr; only 74 have TE <= TC. No unresolved point remains
in the runtime table.

For commands, per-region shapes, provenance notes, and reproduction details,
read EMISSION_SCAN_GUIDE.md section 2026-07-20 TC 1500 K Production Extension.

## 2026-07-20 emission lookup TE extension

The current production accident lookup envelope is TE 700-3000 K at 20 K
spacing and TC 500-1500 K at 10 K spacing. The preferred runtime directory is
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr_tc1500_te3000, with the tc1500
and legacy directories retained as fallbacks.

The full table contains 25,957,908 points. The added TE 2420-3000 K range has
3,381,480 points, including 302,197 raw analytic failures: 9,796 were safely
zero-filled and 292,401 were neighbor-imputed. The final table has zero
unresolved points. Treat TE > 2400 K as an accident lookup extension rather
than validation of the empirical high-temperature model. Detailed counts,
artifact sizes, and boundary checks are in EMISSION_SCAN_GUIDE.md.

## 2026-07-09 case-level lookup control

`ReactorCore` now accepts `tec_lookup_enabled`, `tec_lookup_db`, and `tec_lookup_regions` and passes them to every `ThermoCalcModel` it creates. `testModule/Full_Loop_Cases/Full_Loop_Cases_10kW.FullLoopCoreConfig` exposes the same fields so V14_10kW/V15-style full-loop cases can enable or disable lookup without relying on process-wide environment variables.

Precedence:

```text
explicit enable_lookup True/False > THERMOCALC_ENABLE_LOOKUP
explicit lookup_db              > THERMOCALC_LOOKUP_DB
explicit lookup_regions         > THERMOCALC_LOOKUP_REGIONS
```

Verification:

```text
E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m unittest testModule.test_thermocalc_lookup_config testModule.test_full_loop_10kw_core_geometry.FullLoop10kWCoreGeometryTests.test_v14_10kw_tec_lookup_config_reaches_thermocalc
```

## 2026-08-06 Series fixed-resistance bracketed load-line solve

circuitTECs::resistanceFixedCircuitCalc() now solves
U_TEC(I) - I*Rload = 0 with a bracketed secant/bisection search centered on
the previous valid current or configured current guess. Every sampled
circuitCalc(I) must report a finite voltage and converged=true; acceptance
requires both a 1e-3 V load-line residual and a 1e-2 A bracket. Invalid loads or
a missing valid bracket return finite Iout=Uout=0 with converged=false
instead of reporting a failed solve as converged.

ThermoCalcModel.calculate() clears all per-element electrical and Joule-power
fields whenever series fixed resistance returns zero output, while preserving
the C++ convergence flag and iteration count. The normal low-temperature
zero-emission guard remains a distinct successful skipped state.

The Python 3.12 production extension was rebuilt as
ThermoCalc/te_solver.cp312-win_amd64.pyd, SHA256
AC513C79BDB9A08E69CC1381A3A5763CCD6A4FCA329A3519EC697AC8CC96355D.
Focused interface and series-circuit tests pass. At the V14 t=2850 s
restart with 58 TECs and Rload=0.003 ohm, the solver converges in 15 sampled
current evaluations to I=218.703142 A, U=0.656109428 V, and
Pload=143.493193 W; all exported node fields are finite.
### 2026-08-06 fixed-R recovery after zero-output cleanup

A persistent startup circuit can legitimately return zero output before the
emitter field is hot enough. The Python cleanup then sets Uout to zero. The
next fixed-R solve now restores its first circuit voltage seed from the load
line before entering circuitCalc(); an invalid sampled point also forces the
following point to reinitialize. This allows the same circuit object to
recover when later temperature updates create a valid operating point.

The focused series regression explicitly zeros Iout/Uout on a hot fixed-R
model and verifies recovery to the fixed-U reference load line. The rebuilt
production extension SHA256 is
AC513C79BDB9A08E69CC1381A3A5763CCD6A4FCA329A3519EC697AC8CC96355D.
