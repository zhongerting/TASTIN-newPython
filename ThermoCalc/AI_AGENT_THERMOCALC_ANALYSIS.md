# ThermoCalc Codex 快速接管手册

> 更新时间：2026-06-23  
> 适用范围：`ThermoCalc/` 模块及其必要的上层集成点。  
> 使用方式：后续 Codex 首次接管 ThermoCalc 时先读本文件；只有在修改具体功能时，才按“修改场景索引”继续回查源码。

## 1. 先看结论

`ThermoCalc/` 是 TASTIN 中的热离子能量转换（TEC）求解模块。Python 上层向它提供多根 TFE 的轴向温度、几何、电路参数；C++ 核心求解串联电路以及每根 TFE 内各轴向节点的发射、电流密度和电势分布。

当前源码的分层调用链是：

```text
Components/ReactorCore.py 或 Components/TECCircuitManager.py
    -> ThermoCalc/ThermoCalcWrapper.py
    -> ThermoCalc/bindings.cpp
    -> circuitTECs
    -> singleThermionicEnergyConversion
    -> thermionicEmission
```

2026-06-01 已闭合逐节点侧面积、`phiE/phiC/Vd` 结果读取和铯池温度运行时热更新，并使用 Python 3.12 重新构建 `te_solver.cp312-win_amd64.pyd`。`fixed_I` 不在本轮实现范围内，包装层会明确拒绝该模式。

2026-06-23 新增热离子发射查表加速实验路径。该路径保持原解析 `thermionicEmission::calc()` 可用，查表仅在显式启用时作为 `calc()` 的优先分支；表缺失或关闭时继续走原解析法。当前测试版扩展位于 `ThermoCalc/build_cp312/Release/te_solver.cp312-win_amd64.pyd`，根目录生产 `ThermoCalc/te_solver.cp312-win_amd64.pyd` 未被替换。

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
| [`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 多根 TFE 串联电路：定电压和定电阻模式的全局迭代 |
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

`circuitTECs` 保存 `vector<singleThermionicEnergyConversion*> TECs`，负责多根 TFE 的串联电路求解。

当前绑定可触发的顶层入口是 `circuitTECsCalc()`，Python 名为 `calc()`：

```text
isFixedU == true -> uFixedCircuitCalc()
isFixedR == true -> resistanceFixedCircuitCalc()
```

两个模式都会通过 `circuitCalc(I)` 在给定总电流下求串联电路，再用弦割式迭代修正电流。C++ 内部确实有“给定电流计算”能力 `circuitCalc(double I)`，但当前 pybind11 枚举和 Python 构建流程没有暴露完整的定电流公共模式。

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
| `ThermoCalcModel(n_elements, n_nodes)` | 创建 `te_solver.InputData()` 并填默认值 | 若 `te_solver` 未成功导入，实例化仍会失败 |
| `set_temperatures(T_em, T_co)` | 校验 `(N_elem, n_node)`，保存副本；电路已构建时写入每根 `SingleTEC` | 当前绑定暴露了 `Temitter`、`Tcollector` |
| `setup_circuit_mode(mode_str, target_value, I_guess=150.0)` | 接受 `fixed_R`、`fixed_U`；对 `fixed_I` 明确抛出 `ValueError` | `fixed_I` 未暴露 |
| `build()` | 把温度写入 `InputData`，调用 `te_solver.create_circuit()` | 绑定层在 `unchecked<>` 前执行完整形状校验 |
| `calculate(verbose=False)` | 必要时自动 `build()`，再调用 `_circuit.calc()`；返回耗时 `[ms]` | 不是物理时间步长度 |
| `load_emission_lookup_database(db_dir, enable=True, force=False)` | 从 `manifest.json`、`chunk_plan.json` 和 chunk `.npz` 加载查表数据库到 C++ 单例 | 优先加载同名 `.optimized.npz`；仅测试版 pyd 暴露该 API |
| `get_global_results()` | 返回 `Iout`、`Uout`、`Rload` | 电路未构建时返回 `None` |
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

`THERMOCALC_PYD_DIR` 会被插入到 `sys.path` 的最高优先级，用于在不覆盖根目录 `.pyd` 的情况下测试新扩展。只有同时设置 `THERMOCALC_ENABLE_LOOKUP=1` 和 `THERMOCALC_LOOKUP_DB` 时，`ThermoCalcModel.__init__()` 才会自动加载查表数据库。

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
| `R_load_init` | 标量 | `Ohm` | 仅在 C++ 结构体中声明；当前未绑定、未由包装层赋值、未参与工厂函数逻辑 |

## 7. 结果提取

### 7.1 已由当前绑定暴露

`CircuitTECs` 暴露：

| 字段 | 含义 |
|---|---|
| `TECs` | 单根 TFE 对象列表 |
| `Utarget`、`Rload` | 目标电压、负载电阻 |
| `Iout`、`Uout` | 总电流、总电压 |
| `isFixedU`、`isFixedR` | 模式标志 |
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
    -> 按 tec_multipliers 创建虚拟串联元件
    -> _configure_thermo_calc_geometry()

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
- 当前绑定枚举只有 `FixedVoltage` 和 `FixedResistance`。
- 当前绑定层把 `sideAreaE/sideAreaC` 作为 `(N_elem, n_axi)` 二维数组逐行送入单根 TFE。
- 当前单根 TFE 的 `sideAreaE/sideAreaC` 是逐节点 `vector<double>`。
- `dlE/dlC` 在单根 TFE 内是向量，当前 `VcalcFVM()` 会读取节点长度。
- 当前 `SingleTEC` 绑定暴露 `phiE`、`phiC`、`Vd`。
- 当前 `CircuitTECs` 绑定提供 `set_tcs()`，并逐根更新 `SingleTEC.Tcs`。
- 上层两条路径都依赖 `ThermoCalcModel`，`ReactorCore` 还会写入二维逐节点侧面积。
- 测试版 `te_solver` 暴露热离子查表 API；`ThermoCalcWrapper.py` 可通过环境变量加载 `ThermoCalc/emission_database`。
- 查表数据库的全量计划、优化表和验证摘要见 [`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md)。

### B. 需要 Python 3.12 编译产物进一步验证

- 根目录生产 `te_solver.cp312-win_amd64.pyd` 是否切换到当前查表源码构建；截至 2026-06-23 仍保持旧生产 pyd 不动。
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
| `fixed_I` 公共模式未实现 | 包装层明确拒绝 `fixed_I`；绑定枚举只有 `FixedVoltage`、`FixedResistance` | 不要宣称 Python 支持定电流模式 |
| 非均匀网格完整数学验证仍有限 | `dlE/dlC` 和侧面积已逐节点化，但 `VcalcFVM()` 的界面距离仍沿用现有离散公式 | 修改电势离散时补专项守恒验证 |
| `set_rload()` 构建前行为不完整 | 包装层尝试写 `_input_data.Rload/R_load`，当前 `InputData` 未绑定这些字段 | 构建后可写 `CircuitTECs.Rload`；构建前优先使用 `setup_circuit_mode('fixed_R', ...)` |
| 原始指针生命周期风险 | `create_circuit()` 为每根 TFE `new` 对象；相关析构函数当前为空 | 若处理长时运行内存问题，专项检查所有权与释放逻辑 |
| 查表仅在测试版 pyd 中验证 | 新扩展位于 `ThermoCalc/build_cp312/Release`，根目录生产 pyd 未替换 | 运行查表路径必须显式设置 `THERMOCALC_PYD_DIR` |
| 查表表外点会回退解析 | `thermionicEmission::calc()` 查表 miss 后继续原解析法 | 若希望完全避免解析失败输出，应扩大/修正表覆盖或改电路层策略 |
| setup 阶段首次 TEC 仍可能慢且打印失败 | V13 `apply_wire_resistance()` 会重建电路并立即 `calculate()`；30 s 查表计时中 setup 首算约 `7.81 s` 且打印失败信息 | 正式推进 warm-start 后 TEC 单次约 `1.8 s`；setup 首算需单独优化 |

## 11. 修改场景索引

| 修改场景 | 必读文件 | 重点检查 |
|---|---|---|
| 改 Python 公共 API、结果字段、热更新 | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp) | Python 名称、绑定字段、形状、构建前后行为必须成对闭合 |
| 增加或修复 `fixed_I` | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp)、[`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 枚举、`build()` 分支、顶层分发、目标电流语义 |
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

2026-06-01 已使用：

```text
C:\Users\HC Zhao\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe
Python 3.12.13
```

重新构建并验证 `te_solver.cp312-win_amd64.pyd`。`testModule/test_thermocalc_interface.py` 覆盖形状拒绝、均匀与非均匀节点面积、`phiE/phiC/Vd`、构建前后温度和 `Tcs` 更新、`fixed_I` 显式拒绝。单 TFE TEC `1 s` 基线也已运行。

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

旧数据库位于 `ThermoCalc/emission_database/`，使用的是旧铯压范围。旧全量计划含 `18,737,388` 点、`78` 个 chunk；已完成 `startup` 和 `accident` 风险点优化表，`.optimized.npz` 会被包装层优先加载。优化后 `startup+accident` 中原始无效点 `55,506` 个，其中 `43,104` 个按零发射处理，`12,402` 个由邻域插补，未解决点为 `0`。新的 `0.02-3.0 torr` 全量计划为 `31,716,828` 点、`108` 个 chunk，需要重新生成数据库和 dense runtime v2 表。

验证结果：

- `testModule/test_thermocalc_lookup.py` 通过；单点批量查表约 `1.53e6 points/s`，解析法约 `5.84e4 points/s`，局部函数约 `26x` 加速。
- V13 `1 s` 查表 smoke 可完成，`tec_coupled_enabled=True`；热工量与解析基准接近，TEC 端电功率高约 `62.2 W`，约 `1.14%`。
- V13 查表分段长算已从 `21000 s` 推进到 `22000 s`，最终 `tec_coupled_enabled=True`，未检出 `disabling TEC`、`Traceback` 或运行阶段解析收敛失败文本。
- V13 真实推进 `30 s` 计时：推进总耗时 `422.02 s`，TEC 计算 `51.24 s` 占 `12.16%`，导热 `252.70 s` 占 `59.98%`，流动 `3.41 s` 占 `0.81%`，其他系统层开销 `27.04%`。
- 查表 warm-start 后 TEC 单次更新约 `1.8 s`；setup 阶段导线电阻重建后的首次 TEC 约 `7.81 s`。

当前判断：查表路径已经把 TEC 从主要瓶颈之一降为次要瓶颈；V13 长算主耗时转移到导热求解和系统层调度。后续若继续优化速度，优先检查导热 solid 数量、辐射器管壁求解、coupler/组件调度开销，以及 `circuitTECs` 外层迭代次数，而不是继续只优化局部 `thermionicEmission` 单点。

## 15. 2026-06-23 runtime 查表压缩与索引

本轮在不删除原解析法、不替换根目录生产 `.pyd` 的前提下，新增了运行时专用查表格式和 C++ 查询索引：

- `tools/emission_database.py export-runtime` 从 `ThermoCalc/emission_database/` 导出 `ThermoCalc/emission_runtime_db/`，优先读取 `.optimized.npz`，输出 `runtime_manifest.json` 和按 region 分开的 `*.runtime.npz`。
- runtime 表只保留 `TE_axis/TC_axis/Vo_axis/Tcs_axis`、`J/Vd/delta_V/phiE/phiC`、`lookup_safe/zero_mask`；默认字段精度为 `float32`，但 `phiE/phiC` 保留，供边界条件继续调用。
- `zero_mask` 标记安全零电流区；启用 `--zero-compress` 后这些点的 `J` 在 runtime 表中直接写为 `0`，但电压和功函数字段仍参与插值。
- `ThermoCalcWrapper.load_emission_lookup_database(..., regions=...)` 同时支持旧全量库和 runtime 库；默认只加载 `core`，可用环境变量 `THERMOCALC_LOOKUP_REGIONS=core,startup,high_power,accident` 扩展覆盖。
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

`ThermoCalc/emission_runtime_db/` 和 `ThermoCalc/emission_database/` 都是生成数据，不应提交到 git。若后续需要完整工况覆盖，先导出全 region runtime 表，再设置 `THERMOCALC_LOOKUP_DB` 指向 runtime 目录，并显式设置 `THERMOCALC_LOOKUP_REGIONS`。

2026-06-23 后续修复：

- `chunk_te_ranges()` 现在为每个 TE chunk 保留右侧边界平面，避免新生成数据库出现 `1300-1310 K` 后直接跳到 `1320-1330 K` 的插值空隙。
- `export-runtime` 对旧数据库自动拼接下一 chunk 的第一个 TE 平面；现有 core runtime 表由 `1300-1310 K`、`1320-1330 K` 等旧块导出为 `1300-1320 K`、`1320-1340 K` 等连续块，不需要重算原始 1873 万点。
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
  -> export-runtime
  -> ThermoCalc/emission_runtime_db/*.runtime.npz

运行调用:
  ThermoCalcModel.__init__()
  -> load_emission_lookup_database()
  -> te_solver.add_emission_runtime_block()
  -> emissionLookup.cpp 内存索引
  -> thermionicEmission::calc()
  -> queryEmissionLookup()
  -> 命中则返回 J/Vd/delta_V/phiE/phiC
  -> 未命中则回退原解析 calc()
```

关键边界：

- `ThermoCalc/emission_database/` 是原始/审计库，保留诊断字段和 `.optimized.npz` sidecar，不提交 git。
- `ThermoCalc/emission_runtime_db/` 是运行库，只保留 `J/Vd/delta_V/phiE/phiC/lookup_safe/zero_mask` 和轴，不提交 git。
- 自动加载需要同时设置 `THERMOCALC_ENABLE_LOOKUP=1` 和 `THERMOCALC_LOOKUP_DB`；默认只加载 `core`，更广覆盖由 `THERMOCALC_LOOKUP_REGIONS` 控制。
- 当前查表仅在 `ThermoCalc/build_cp312/Release` 测试版 `.pyd` 中验证，根目录生产 `.pyd` 仍未替换。

## 17. 2026-06-23 dense runtime v2 补充

当前推荐的运行时查表格式是 `export-runtime-dense` 生成的 dense runtime v2：

```text
ThermoCalc/emission_runtime_db_v2/
  runtime_dense_manifest.json
  core.runtime.v2.npz
  core.runtime.v2.tedb
```

该格式按 region 存储一个连续四维张量，字段为 `J/Vd/delta_V/phiE/phiC`，并把 `lookup_safe` 和 `zero_mask` 压缩为 bit-packed mask。`.npz` 是可移植格式，`.tedb` 是 C++ 直接加载格式；包装层发现 `runtime_dense_manifest.json` 后会优先加载 `.tedb`，否则回退到 `.npz`。

旧 core dense v2 表从旧压力范围本地全量库导出，形状为 `86 x 41 x 71 x 41`，共 `10,264,186` 点；`NPZ` 约 `86.87 MiB`，`TEDB` 约 `198.22 MiB`，`TEDB` 加载约 `0.167 s`。连续 core 随机 `200000` 点批量查表约 `1.49e6 points/s`；聚焦回归 `testModule/test_thermocalc_lookup.py` 中查表约 `3.55e6 points/s`、解析法约 `9.87e4 points/s`、约 `36x`。

2026-06-23 压力范围修正：新的全量 plan 将 `core/startup/high_power/accident` 的铯压轴统一为 `0.02-3.0 torr`、`61` 个 log-spaced 点。此前已经生成的 `ThermoCalc/emission_database/` 和 `ThermoCalc/emission_runtime_db_v2/` 属于旧压力范围产物，后续正式查表使用前需要重新执行 `plan -> worker -> summarize/verify -> optimize-table -> export-runtime-dense`。

`ThermoCalc/emission_runtime_db_v2/` 是生成数据，不提交 git。完整复现命令、字段说明和 v1/v2 对比维护在 [`EMISSION_SCAN_GUIDE.md`](./EMISSION_SCAN_GUIDE.md)。
