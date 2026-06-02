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

2026-06-01 已闭合逐节点侧面积、`phiE/phiC/Vd` 结果读取和铯池温度运行时热更新，并使用 Python 3.12 重新构建 `te_solver.cp312-win_amd64.pyd`。`fixed_I` 不在本轮实现范围内，包装层会明确拒绝该模式。

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

## 5. Python 公共入口

`ThermoCalcModel` 位于 [`ThermoCalcWrapper.py`](./ThermoCalcWrapper.py)。

| 入口 | 行为 | 当前注意事项 |
|---|---|---|
| `ThermoCalcModel(n_elements, n_nodes)` | 创建 `te_solver.InputData()` 并填默认值 | 若 `te_solver` 未成功导入，实例化仍会失败 |
| `set_temperatures(T_em, T_co)` | 校验 `(N_elem, n_node)`，保存副本；电路已构建时写入每根 `SingleTEC` | 当前绑定暴露了 `Temitter`、`Tcollector` |
| `setup_circuit_mode(mode_str, target_value, I_guess=150.0)` | 接受 `fixed_R`、`fixed_U`；对 `fixed_I` 明确抛出 `ValueError` | `fixed_I` 未暴露 |
| `build()` | 把温度写入 `InputData`，调用 `te_solver.create_circuit()` | 绑定层在 `unchecked<>` 前执行完整形状校验 |
| `calculate(verbose=False)` | 必要时自动 `build()`，再调用 `_circuit.calc()`；返回耗时 `[ms]` | 不是物理时间步长度 |
| `get_global_results()` | 返回 `Iout`、`Uout`、`Rload` | 电路未构建时返回 `None` |
| `get_tec_results(idx)` | 返回指定 TFE 的详细字典 | 包装层读取了绑定未暴露的字段，见风险清单 |
| `set_tcs(tcs_val)` | 接受标量或 `(N_elem, n_node)`，更新 `_input_data.Tcs` | 构建后调用电路级 `set_tcs()` 逐 TEC 热更新 |
| `set_rload(rload_val)` | 尝试更新输入结构；构建后尝试更新 C++ 电路负载 | 构建后的 `CircuitTECs.Rload` 属性已绑定；构建前输入结构没有对应可写绑定 |

## 6. `InputData` 字段、维度与默认值

`bindings.cpp` 定义 `InputData`，`ThermoCalcWrapper.py` 负责填充。`create_circuit()` 在进入 `unchecked<1>()` 和 `unchecked<2>()` 前对所有数组执行显式形状校验。

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

### B. 需要 Python 3.12 编译产物进一步验证

- `te_solver.cp312-win_amd64.pyd` 是否确实由当前工作树源码构建。
- 更长时 TEC 瞬态下，电子边界功率差与端功率、焦耳热之间的稳定离散误差门槛。
- `set_rload()` 在目标运行环境中的构建前、构建后行为。
- 定电压、定电阻模式在目标案例中的数值结果和收敛性。

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

2026-06-01 已使用：

```text
C:\Users\HC Zhao\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe
Python 3.12.13
```

重新构建并验证 `te_solver.cp312-win_amd64.pyd`。`testModule/test_thermocalc_interface.py` 覆盖形状拒绝、均匀与非均匀节点面积、`phiE/phiC/Vd`、构建前后温度和 `Tcs` 更新、`fixed_I` 显式拒绝。单 TFE TEC `1 s` 基线也已运行。
