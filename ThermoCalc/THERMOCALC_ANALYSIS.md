# ThermoCalc 模块分析文档

## 1. 模块概述

`ThermoCalc/` 是 TASTIN 的热离子能量转换 (Thermionic Energy Conversion, TEC) 计算模块。它封装了一个 **C++ 原生计算库**，通过 **pybind11** 暴露给 Python 端使用。

**核心职责：**
1. 接受各 TFE 元件沿轴向的**发射极/接收极温度分布**及几何电学参数
2. 求解**全局串联电路**（定电压/定电阻/定电流模式）
3. 返回每根元件的**电流密度分布 J(z)**、**电压分布 V(z)**、**电极电势 UE/UC**、**电阻率分布**等物理场

**架构层次：**

```
Python 层 (TASTIN)
    │
    ├── TECCircuitManager.pre_step()    ← 每时间步调用
    │     └── ThermoCalcModel.calculate()
    │           └── t0..circuit.calc()   ← 触发 C++ 求解
    │
    ├── ThermoCalcModel (Python 封装)
    │     ├── InputData 准备
    │     ├── te_solver.create_circuit() → C++
    │     └── 结果提取 (get_tec_results)
    │
    └── te_solver.pyd (编译产物, 不在仓库中)
          └── bindings.cpp (pybind11 绑定)
                └── circuitTECs / singleThermionicEnergyConversion (C++)
```

---

## 2. C++ 层数据结构

### 2.1 `thermionicEmission` — 热离子发射单元

**文件：** `thermionicEmission.h/cpp`

这是最底层的物理模型，描述**单段轴向切片**上发射极-接收极间的物理过程。

**关键成员变量：**

| 变量 | 类型 | 含义 |
|---|---|---|
| `TE` | double | 发射极表面温度 [K] |
| `TC` | double | 接收极表面温度 [K] |
| `Tcs` | double | 铯池温度 [K] |
| `d` | double | 电极间距 [mm] |
| `phiE` | double | 发射极功函数 [eV] |
| `phiC` | double | 接收极功函数 [eV] |
| `Vd` | double | 电弧降 (arc voltage drop) [V] |
| `Vo` | double | 电极电势差 [V] |
| `J` | double | 电流密度 [A/cm²] |
| `P` | double | 铯压力 |
| `JE` / `JC` | double | 发射极/接收极电子电流密度 |
| `TeE` / `Te` / `TeC` | double | 发射极表面/平均/接收极表面电子温度 |

**核心方法：**
- `calc()` → 根据放电模式（阻塞/过渡/饱和）计算 J
- `obstructedCalc()` — 阻塞模式 (低电流密度)
- `transitionCalc()` — 过渡模式
- `saturationCalc()` — 饱和模式 (高电流密度)

**关键物理常数：**
- Richardson 常数 `A = 120 A/(cm²·K²)`
- Boltzmann 常数倒数 `k = 1/11605 eV/K`

---

### 2.2 `singleThermionicEnergyConversion` — 单根 TFE 元件

**文件：** `singleThermionicEnergyConversion.h/cpp`

代表**一根完整的热离子燃料元件**，沿轴向离散为 `n_node` 个 `thermionicEmission` 单元串。

**构造函数参数 `input[13]`：**

| 索引 | 含义 | 维度 |
|---|---|---|
| `[0]` | Temitter — 发射极温度分布 | n_node |
| `[1]` | Tcollector — 接收极温度分布 | n_node |
| `[2]` | dlE — 发射极单元长度分布 | n_node |
| `[3]` | dlC — 接收极单元长度分布 | n_node |
| `[4]` | {crossAreaE, crossAreaC} — 发射极/接收极截面积 | 2 |
| `[5]` | {sideAreaE} — 发射极侧面积 (用于计算电子冷却) | 1 |
| `[6]` | {sideAreaC} — 接收极侧面积 | 1 |
| `[7]` | resistanceWire — 导线电阻 | 4 |
| `[8]` | {U_init, d_gap} — 电压初值 + 电极间距 | 2 |
| `[9]` | Tcs — 铯池温度分布 | n_node |
| `[10]` | V_init — 极板电势差初值分布 | n_node |
| `[11]` | {Itarget} — 目标电流 | 1 |
| `[12]` | wireU — 导线电压 | 4 |

**关键输出成员变量（pybind11 暴露给 Python）：**

| 变量 | 维度 | 含义 |
|---|---|---|
| `J` | n_node | 电流密度分布 [A/cm²] |
| `V` | n_node | 电极电势差分布 V(z) [V] |
| `UE` | n_node | 发射极电势分布 [V] |
| `UC` | n_node | 接收极电势分布 [V] |
| `rhoE` | n_node | 发射极电阻率 [Ω·m] |
| `rhoC` | n_node | 接收极电阻率 [Ω·m] |
| `IEsecSingle` | n_node | 发射极内部截面电流 [A] |
| `ICsecSingle` | n_node | 接收极内部截面电流 [A] |
| `phiE` | — | 发射极功函数 [eV] (由 thermionicEmission 提供) |
| `phiC` | — | 接收极功函数 [eV] |
| `Vd` | — | 电弧降 [V] |
| `I` | 标量 | 元件总电流 [A] |
| `U` | 标量 | 元件总电压 [V] |
| `P` | 标量 | 元件总功率 [W] |
| `Temitter` | n_node | 发射极温度 (输入，也可读回) |
| `Tcollector` | n_node | 接收极温度 (输入，也可读回) |
| `Tcs` | n_node | 铯池温度 (输入，也可读回) |

**核心方法：**
- `initial()` — 初始化 thermionicUnits 池，分配内存
- `Icalc()` — **电流计算**，根据当前电极电势分布确定电流 (串联一致性)
- `Jcalc()` — 调用每个 thermionicEmission 单元的 `calc()`，更新 J 分布
- `Vcalc()` / `VcalcFVM()` — 电压分布计算 (有限体积法)
- `ICIEcalc()` — 内部截面电流分布计算

---

### 2.3 `circuitTECs` — 全局串联电路

**文件：** `circuitTECs.h/cpp`

管理 **N_elem 根 TFE 元件的串联电路**。

**关键成员变量：**

| 变量 | 含义 |
|---|---|
| `TECs` | `vector<singleThermionicEnergyConversion*>` — 所有元件指针 |
| `nTECs` | 元件数量 |
| `Iout` | 电路总电流 [A] |
| `Uout` | 电路总电压 [V] |
| `Rload` | 外部负载电阻 [Ω] |
| `Utarget` | 目标电压 [V] (定电压模式) |
| `isFixedU` | 是否定电压模式 |
| `isFixedR` | 是否定电阻模式 |
| `IE` | 各元件电流分布 |
| `deltaU1` / `deltaU2` | 导线压降 |

**核心方法：**
- `circuitTECsCalc()` — **顶层入口**，根据模式调用对应求解函数
- `circuitCalc(I)` — 固定电流电路计算
- `uFixedCircuitCalc()` — 固定电压电路计算 (牛顿迭代求 I)
- `resistanceFixedCircuitCalc()` — 固定电阻计算 (迭代求公共电流)
- `singleTECU(deltaV, n)` — 给第 n 根元件施加电压 deltaV，计算其 I

**求解模式：**

| 模式 | 说明 | 迭代变量 |
|---|---|---|
| `FixedResistance` | 给定 R_load，求稳态工作点 | 迭代 I_out 使 U_out = I_out × R_load |
| `FixedVoltage` | 给定 U_target，求电流 | 牛顿迭代求 I 使 U_out = U_target |
| `FixedCurrent` | 给定 I_target | 直接计算各元件 U |

---

## 3. Python 层接口 (ThermoCalcWrapper.py)

### 3.1 `ThermoCalcModel` 类

这是 TASTIN 中直接使用的 Python 封装类。

#### 3.1.1 初始化与配置

```python
model = ThermoCalcModel(n_elements=6, n_nodes=50)
```

初始化时自动填充**默认几何与电学参数**：

| 参数 | 默认值 | 含义 |
|---|---|---|
| `Tcs` | 600.0 K | 铯池温度 (全局均匀) |
| `V_init` | 0.2 V | 极板电势差初值 |
| `dlE / dlC` | 0.507/n_node m | 发射极/接收极单元长度 |
| `crossAreaE` | 6.667e-5 m² | 发射极截面积 |
| `crossAreaC` | 1.0786e-4 m² | 接收极截面积 |
| `sideAreaE` | 0.00092855×25/n_node m² | 发射极侧面积 |
| `sideAreaC` | 0.00097593×25/n_node m² | 接收极侧面积 |
| `U_init` | 1.6 V | 元件总电压初值 |
| `d_gap` | 0.5 mm | 电极间距 |
| `Itarget` | 200.0 A | 目标电流初值 |
| `resistanceWire` | zeros(N, 4) | 导线电阻 |
| `wireU` | [0.8, 0.8, 0.0, 0.0] | 导线电压边界 |
| `mode` | FixedVoltage | 默认定电压模式 |
| `target_val` | 0.89×34 ≈ 30.26 V | 目标电压 |
| `I_total_init` | 284.0 A | 初始电流猜测 |

#### 3.1.2 核心方法调用流程

```
ThermoCalcModel 生命周期:
  │
  ├── __init__()               创建 InputData，填充默认参数
  │
  ├── set_temperatures(T_em, T_co)   更新各元件轴向温度分布
  │     └── 如果 circuit 已构建 → 直接写 C++ 对象的 Temitter/Tcollector
  │
  ├── setup_circuit_mode(mode, value)    设置电路模式
  │     ├── 'fixed_R'  → 定电阻 (Rload = value)
  │     ├── 'fixed_U'  → 定电压 (Utarget = value)
  │     └── 'fixed_I'  → 定电流 (Itarget = value)
  │
  ├── set_tcs(val)             更新铯池温度 (标量或二维数组)
  │     └── 写 InputData.Tcs，若 circuit 已构建则热更新
  │
  ├── set_rload(val)           更新负载电阻
  │     └── 写 InputData，若 circuit 已构建则热更新
  │
  ├── build()                  首次构建 C++ circuit 对象
  │     ├── 将 InputData 中的温度场推入
  │     ├── 调用 te_solver.create_circuit(input_data)
  │     └── 根据模式设置 isFixedU/isFixedR 及 Uout/Iout/Rload
  │
  ├── calculate(verbose)       触发 C++ 计算 (瞬态每步调用)
  │     ├── 若 circuit 未构建 → 自动 build()
  │     ├── 调用 circuit.calc()  → C++ 迭代求解
  │     └── 返回计算耗时 [ms]
  │
  ├── get_global_results()     获取系统级结果
  │     └── {"Iout": I, "Uout": U, "Rload": R}
  │
  └── get_tec_results(idx)     获取第 idx 根元件的详细结果
        └── {
              "I": float,          总电流 [A]
              "U": float,          总电压 [V]
              "J": ndarray[n_node],   电流密度 [A/cm²]
              "V": ndarray[n_node],   电极电势差 [V]
              "UE": ndarray[n_node],  发射极电势 [V]
              "UC": ndarray[n_node],  接收极电势 [V]
              "rhoE": ndarray[n_node],  发射极电阻率 [Ω·m]
              "rhoC": ndarray[n_node],  接收极电阻率 [Ω·m]
              "IEsecSingle": ndarray[n_node],  发射极截面电流 [A]
              "ICsecSingle": ndarray[n_node],  接收极截面电流 [A]
              "phiE": ndarray,       发射极功函数 [eV]
              "phiC": ndarray,       接收极功函数 [eV]
              "Vd": ndarray,         电弧降 [V]
              "TE": ndarray[n_node], 发射极温度 [K]
              "TC": ndarray[n_node]  接收极温度 [K]
            }
```

---

## 4. Pybind11 绑定层 (bindings.cpp)

### 4.1 数据结构

**`InputData` (struct → py::class_)：** 纯数据容器，Python 端创建填充后传入 `create_circuit()`

**暴露的模式：**
- `CalculationMode::FixedVoltage` — 定电压
- `CalculationMode::FixedResistance` — 定电阻

### 4.2 核心工厂函数 `create_circuit(InputData)`

```
create_circuit(data)
  ├── 创建 circuitTECs 实例
  ├── for i in [0, N_elements):
  │     ├── 从 InputData 中提取第 i 行的 13 组输入数据
  │     ├── new singleThermionicEnergyConversion(input)
  │     ├── tec->initial()  ← 分配 thermionicUnits 池
  │     └── circuit->TECs.push_back(tec)
  ├── 设置 circuit->nTECs, Iout
  ├── 根据 data.mode 设置 isFixedU / isFixedR / Utarget / Rload
  └── 返回 unique_ptr<circuitTECs> (转移所有权给 Python)
```

### 4.3 暴露的 C++ 类 (py::class_)

**`SingleTEC`：**
- 暴露了所有输出物理场为可读写属性
- 暴露了 `Icalc()` 方法

**`CircuitTECs`：**
- 暴露了 `TECs` 列表（Python 可以直接索引到每个 SingleTEC）
- 暴露了全局控制参数 (`Utarget`, `Rload`, `Iout`, `Uout`, `isFixedU`, `isFixedR`)
- 暴露了 `calc()` 方法 → 调用 `circuitTECsCalc()`

---

## 5. 在 TASTIN 中的调用链

### 5.1 TECCircuitManager 调用流程

```python
# TECCircuitManager.pre_step(dt, current_time):
#   每个物理时间步的"显式电热计算"

# 1. 收集所有 TFE 的发射极/接收极温度
for tec in self.tfe_list:
    self._T_emit_matrix[i, :] = tec.emitter.T_surface
    self._T_coll_matrix[i, :] = tec.collector.T_surface

# 2. 设置温度到 ThermoCalcModel
self.circuit.set_temperatures(self._T_emit_matrix, self._T_coll_matrix)

# 3. 触发 C++ 电路求解
self.circuit.calculate()

# 4. 提取结果并转换为热流
for i, tec in enumerate(self.tfe_list):
    results = self.circuit.get_tec_results(i)
    
    # 电流密度 J [A/cm²] → 电子冷却/加热热流 [W/m²]
    q_cooling = compute_electron_cooling(results['J'], results['phiE'], ...)
    q_heating  = compute_electron_heating(results['J'], results['phiC'], ...)
    
    # 焦耳热 [W/m³] (电极内部的欧姆热)
    q_joule_E = joule_power_from_electric_field(results['IEsecSingle'], 
                                                  results['rhoE'], tec.crossAreaE)
    q_joule_C = joule_power_from_electric_field(results['ICsecSingle'],
                                                  results['rhoC'], tec.crossAreaC)
    
    # 下发热流给 TFEUnit → 固体导热求解器
    tec.apply_plasma_flux(q_cooling, q_heating)
    tec.apply_joule_heat(q_joule_E, q_joule_C)
```

### 5.2 数据流图

```
TFEUnit (每根燃料元件)
  │  .emitter (HeatConduction2D)  →  T_surface (发射极表面温度)
  │  .collector (HeatConduction2D) →  T_surface (接收极表面温度)
  │
  ▼
TECCircuitManager.pre_step()
  │  收集温度 → T_emit_matrix[N_elem, n_node], T_coll_matrix[N_elem, n_node]
  │
  ▼
ThermoCalcModel
  │  .set_temperatures(T_emit, T_coll)
  │  .calculate()
  │     └── C++ circuit.calc()
  │           ├── 定电压: 牛顿迭代求 I_out 使 U_out = U_target
  │           ├── 定电阻: 迭代求 I_out 使 U_out = I_out × R_load
  │           └── 各元件同步: Icalc() → Jcalc() → Vcalc()
  │
  ▼
ThermoCalcModel.get_tec_results(idx)
  │  → J(z), V(z), UE(z), UC(z), rhoE(z), rhoC(z), IEsec(z), ICsec(z), phiE, phiC, Vd
  │
  ▼
TECCircuitManager._compute_plasma_fluxes() + _compute_joule_heat()
  │
  ├──→ 电子冷却热流 [W/m²] → TFEUnit.emitter (BoundaryRegion, FluxBC)
  ├──→ 电子加热热流 [W/m²] → TFEUnit.collector (BoundaryRegion, FluxBC)
  ├──→ 发射极焦耳热 [W/m³] → TFEUnit.emitter (Q_source)
  └──→ 接收极焦耳热 [W/m³] → TFEUnit.collector (Q_source)
```

---

## 6. 文件清单

| 文件 | 类型 | 职责 |
|---|---|---|
| `thermionicEmission.h/cpp` | C++ 核心 | 热离子发射物理模型 (阻塞/过渡/饱和) |
| `singleThermionicEnergyConversion.h/cpp` | C++ 核心 | 单根 TFE 元件 (沿轴向离散 + 电势求解) |
| `circuitTECs.h/cpp` | C++ 核心 | 串联电路管理器 (全局迭代求解) |
| `NonLinerSolver.h/cpp` | C++ 辅助 | 非线性方程求解器 (牛顿法) |
| `bindings.cpp` | Pybind11 绑定 | InputData ↔ struct, C++ ↔ Python |
| `ThermoCalcWrapper.py` | Python 封装 | 参数管理 + 构建 + 计算 + 结果提取 |
| `__init__.py` | Python 包 | 包初始化 (空) |
| `te_solver.pyd` | 编译产物 | CMake 构建的 Python 扩展模块 (不在仓库中) |
| `CMakeLists.txt` | 构建配置 | C++ 编译配置 |
| `test_real_case*.py` | 测试脚本 | 独立测试热离子电路求解 |
| `NONUNIFORM_GRID_GUIDE.md` | 文档 | 非均匀网格配置指南 |

---

## 7. 关键约束与注意事项

1. **维度一致性：** 所有二维数组必须是 `(N_elem, n_node)`，一维数组必须是 `(N_elem,)`，否则 C++ 端数组越界。

2. **温度单位：** 统一使用 [K]，铯池温度也是 K。

3. **电流密度单位：** J 的单位是 [A/cm²]，在计算热流时需要转换为 [A/m²]。

4. **首次构建 vs 热更新：**
   - `build()` 仅在首次或几何/连接改变时调用
   - `calculate()` 可在每个瞬态时间步直接调用 (如果 circuit 已构建)
   - `set_temperatures()` 支持热更新已构建 circuit 的 C++ 对象属性

5. **定电压模式 (默认)：** `target_val = 0.89 * 34 ≈ 30.26 V`，迭代求解电流使总电压 = 目标值。

6. **编译产物：** `te_solver.pyd` 需要在本地通过 CMake 编译生成，不在 Git 仓库中。`ThermoCalcWrapper.py` 中有 `HAS_TE_SOLVER` 标志用于宽容导入。