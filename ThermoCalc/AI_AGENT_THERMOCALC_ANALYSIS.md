# ThermoCalc Codex 快速接管手册

> 更新时间：2026-05-31  
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

当前工作树中存在接口演进未闭合的问题。尤其不要把“逐节点侧面积”“`fixed_I` Python 模式”“功函数结果读取”“铯池温度运行时热更新”写成已经可用的能力。详见“事实等级”和“已知风险”。

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
| [`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 多根 TFE 串联电路：定电压和定电阻模式的全局迭代 |
| [`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h)、[`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp) | 单根 TFE：轴向节点、导线与电极电势、电流密度和电阻率 |
| [`thermionicEmission.h`](./thermionicEmission.h)、[`thermionicEmission.cpp`](./thermionicEmission.cpp) | 单个轴向节点的热离子发射 J-V 模型：阻塞、过渡、饱和分支 |
| [`NonLinerSolver.h`](./NonLinerSolver.h)、[`NonLinerSolver.cpp`](./NonLinerSolver.cpp) | C++ 辅助求解器；被构建并由头文件引用 |
| [`CMakeLists.txt`](./CMakeLists.txt) | C++17 + pybind11 构建配置，目标模块名为 `te_solver` |
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
| `[5]` | `{sideAreaE}` 发射极侧面积 `[m^2]` | `1`，标量 |
| `[6]` | `{sideAreaC}` 接收极侧面积 `[m^2]` | `1`，标量 |
| `[7]` | `resistanceWire` 导线电阻 `[Ohm]` | `4` |
| `[8]` | `{U_init, d_gap}` 初始总电压 `[V]`、极间距 `[mm]` | `2` |
| `[9]` | `Tcs` 铯池温度 `[K]` | `n_node` |
| `[10]` | `V_init` 极板电势差初值 `[V]` | `n_node` |
| `[11]` | `{Itarget}` 目标电流 `[A]` | `1` |
| `[12]` | `wireU` 导线电压初值 `[V]` | `4` |

`initial()` 初始化电阻率、电阻和 `thermionicUnits`。`Icalc()` 反复执行 `Jcalc()`、`UwireCalc()` 和 `VcalcFVM()`，更新节点电流密度、电势与内部截面电流。当前头文件中 `sideAreaE` 和 `sideAreaC` 仍是 `double`。

### 4.3 `thermionicEmission`

每个实例表示一个轴向节点。输入包括发射极温度 `TE`、接收极温度 `TC`、铯池温度 `Tcs`、极间距 `d` 和极板电势差 `Vo`。`calc()` 先计算阻塞区结果，再按势垒条件进入过渡区和饱和区逻辑，输出节点电流密度 `J [A/cm^2]`。功函数 `phiE`、`phiC` 和电弧降 `Vd` 存在于这一最底层对象中。

## 5. Python 公共入口

`ThermoCalcModel` 位于 [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)。

| 入口 | 行为 | 当前注意事项 |
|---|---|---|
| `ThermoCalcModel(n_elements, n_nodes)` | 创建 `te_solver.InputData()` 并填默认值 | 若 `te_solver` 未成功导入，实例化仍会失败 |
| `set_temperatures(T_em, T_co)` | 校验 `(N_elem, n_node)`，保存副本；电路已构建时写入每根 `SingleTEC` | 当前绑定暴露了 `Temitter`、`Tcollector` |
| `setup_circuit_mode(mode_str, target_value, I_guess=150.0)` | 接受 `fixed_R`、`fixed_U`、`fixed_I` 字符串 | `fixed_I` 与当前绑定不一致，见风险清单 |
| `build()` | 把温度写入 `InputData`，调用 `te_solver.create_circuit()` | 当前侧面积维度不一致可能在这里触发故障 |
| `calculate(verbose=False)` | 必要时自动 `build()`，再调用 `_circuit.calc()`；返回耗时 `[ms]` | 不是物理时间步长度 |
| `get_global_results()` | 返回 `Iout`、`Uout`、`Rload` | 电路未构建时返回 `None` |
| `get_tec_results(idx)` | 返回指定 TFE 的详细字典 | 包装层读取了绑定未暴露的字段，见风险清单 |
| `set_tcs(tcs_val)` | 接受标量或 `(N_elem, n_node)`，更新 `_input_data.Tcs` | 构建后的热更新路径与当前绑定不一致 |
| `set_rload(rload_val)` | 尝试更新输入结构；构建后尝试更新 C++ 电路负载 | 构建后的 `CircuitTECs.Rload` 属性已绑定；构建前输入结构没有对应可写绑定 |

## 6. `InputData` 字段、维度与默认值

`bindings.cpp` 定义 `InputData`，`ThermoCalcWrapper.py` 负责填充。`create_circuit()` 使用 `unchecked<1>()` 和 `unchecked<2>()` 提取 NumPy 数据，没有完整的显式形状校验；调用方必须自行保证维度正确。

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
| `sideAreaE` | 当前绑定读取 `(N_elem,)` | `m^2` | 包装层当前却创建 `(N_elem, n_node)`，每格 `0.00092855424159680002 * 25 / n_node` |
| `sideAreaC` | 当前绑定读取 `(N_elem,)` | `m^2` | 包装层当前却创建 `(N_elem, n_node)`，每格 `0.00097592945800480005 * 25 / n_node` |
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

### 7.2 包装层声明的 `get_tec_results()`

包装层当前尝试返回：

```text
I, U, J, V, UE, UC, rhoE, rhoC,
IEsecSingle, ICsecSingle,
phiE, phiC, Vd,
TE, TC
```

其中 `TE` 和 `TC` 分别从绑定后的 `Temitter`、`Tcollector` 读取。`phiE`、`phiC`、`Vd` 虽存在于节点级 `thermionicEmission`，但当前 `bindings.cpp` 没有把它们汇总到 `SingleTEC`，也没有向 Python 暴露。调用 `get_tec_results()` 时需要重点验证这三个字段。

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

`_configure_thermo_calc_geometry()` 会从热网格复制逐节点 `dlE/dlC` 和逐节点 `sideAreaE/sideAreaC`。这描述的是 Python 上层意图，不代表当前 C++ 绑定已支持逐节点侧面积。

## 9. 事实等级

### A. 已由当前源码确认

- 模块调用链为 `ReactorCore / TECCircuitManager -> ThermoCalcWrapper -> bindings.cpp -> circuitTECs -> singleThermionicEnergyConversion -> thermionicEmission`。
- 当前绑定枚举只有 `FixedVoltage` 和 `FixedResistance`。
- 当前绑定层把 `sideAreaE/sideAreaC` 送入单根 TFE 时按每根元件一个标量读取。
- 当前单根 TFE 的 `sideAreaE/sideAreaC` 是标量 `double`。
- `dlE/dlC` 在单根 TFE 内是向量，当前 `VcalcFVM()` 会读取节点长度。
- 当前 `SingleTEC` 绑定没有 `phiE`、`phiC`、`Vd`。
- 当前 `CircuitTECs` 绑定没有 `set_tcs()` 方法或 `Tcs` 属性，但 `SingleTEC` 有 `Tcs` 属性。
- 上层两条路径都依赖 `ThermoCalcModel`，`ReactorCore` 还会写入二维逐节点侧面积。

### B. 需要 Python 3.12 编译产物进一步验证

- `te_solver.cp312-win_amd64.pyd` 是否确实由当前工作树源码构建。
- 当前二维 `sideAreaE/sideAreaC` 传入 `unchecked<1>()` 后的实际失败形式。
- `get_tec_results()` 读取 `phiE/phiC/Vd` 时的实际异常位置。
- `set_tcs()` 在构建后的运行时热更新是否必然进入包装层的 `AttributeError` 分支。
- `set_rload()` 在目标运行环境中的构建前、构建后行为。
- 定电压、定电阻模式在目标案例中的数值结果和收敛性。

### C. 历史文档中的旧结论

- [`NONUNIFORM_GRID_GUIDE.md`](./NONUNIFORM_GRID_GUIDE.md) 声称逐节点侧面积和任意非均匀网格支持已经完成。
- [`THERMOCALC_ANALYSIS.md`](./THERMOCALC_ANALYSIS.md) 以及本手册旧版本曾把 `fixed_I`、`phiE/phiC/Vd` 提取和若干热更新路径描述为可用能力。

这些旧结论不能覆盖当前源码。修改前必须重新核对对应实现，并使用匹配解释器的扩展做运行验证。

## 10. 已知风险

| 风险 | 当前源码证据 | 处理原则 |
|---|---|---|
| 当前解释器与扩展 ABI 不匹配 | 当前 `python --version` 为 `3.9.13`，目录中的主扩展为 `te_solver.cp312-win_amd64.pyd` | 使用 Python 3.12 环境重新核验运行时接口 |
| `fixed_I` 公共模式未闭合 | 包装层引用 `CalculationMode.FixedCurrent`，绑定枚举只有 `FixedVoltage`、`FixedResistance` | 不要宣称 Python 支持定电流模式 |
| 侧面积维度冲突 | 包装层和 `ReactorCore` 写二维数组；绑定层 `get_scalar()` 按一维读取；单根 TFE 保存标量 | 不要宣称逐节点侧面积已经实现 |
| 非均匀网格正确性未证明 | `dlE/dlC` 是向量，但逐节点侧面积未闭合，`VcalcFVM()` 的界面距离仍直接使用 `dl[i]` | 把指南视为历史资料；修改时补专项验证 |
| 详细结果字段缺口 | 包装层读取 `tec.phiE/phiC/Vd`，当前 `SingleTEC` 未暴露 | 上层等离子体热流路径需要目标环境验证或后续修复 |
| `set_tcs()` 热更新路径不一致 | 包装层查找 circuit 级 `set_tcs` 或 `Tcs`；当前 `CircuitTECs` 均未绑定 | 构建前更新输入可记录，构建后热更新不要当作已实现 |
| `set_rload()` 构建前行为不完整 | 包装层尝试写 `_input_data.Rload/R_load`，当前 `InputData` 未绑定这些字段 | 构建后可写 `CircuitTECs.Rload`；构建前优先使用 `setup_circuit_mode('fixed_R', ...)` |
| 原始指针生命周期风险 | `create_circuit()` 为每根 TFE `new` 对象；相关析构函数当前为空 | 若处理长时运行内存问题，专项检查所有权与释放逻辑 |

## 11. 修改场景索引

| 修改场景 | 必读文件 | 重点检查 |
|---|---|---|
| 改 Python 公共 API、结果字段、热更新 | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp) | Python 名称、绑定字段、形状、构建前后行为必须成对闭合 |
| 增加或修复 `fixed_I` | [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)、[`bindings.cpp`](./bindings.cpp)、[`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 枚举、`build()` 分支、顶层分发、目标电流语义 |
| 暴露 `phiE/phiC/Vd` | [`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h)、[`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp)、[`bindings.cpp`](./bindings.cpp)、[`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py) | 节点级对象如何汇总成 `SingleTEC` 向量 |
| 修复逐节点侧面积或非均匀网格 | [`bindings.cpp`](./bindings.cpp)、[`singleThermionicEnergyConversion.h`](./singleThermionicEnergyConversion.h)、[`singleThermionicEnergyConversion.cpp`](./singleThermionicEnergyConversion.cpp)、[`../Components/ReactorCore.py`](../Components/ReactorCore.py) | 所有积分、电势方程、截面电流、界面距离和测试输入 |
| 改全局串联电路迭代 | [`circuitTECs.h`](./circuitTECs.h)、[`circuitTECs.cpp`](./circuitTECs.cpp) | 定电压/定电阻分支、弦割迭代、`Iout/Uout/Rload/Utarget` |
| 改热离子发射物理 | [`thermionicEmission.h`](./thermionicEmission.h)、[`thermionicEmission.cpp`](./thermionicEmission.cpp) | 阻塞、过渡、饱和分支以及单位 |
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

当前环境只能完成静态源码核对：`python` 为 3.9.13，无法直接导入 `te_solver.cp312-win_amd64.pyd` 做 Python 3.12 运行时验证。
