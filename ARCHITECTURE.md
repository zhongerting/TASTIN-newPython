# TASTIN-Python 程序架构与调用依赖关系文档

## 项目概述

TASTIN (Thermionic Space Nuclear Power System) 是一款**空间核电源热离子反应堆**的全系统瞬态仿真程序，采用 Python 实现，部分热离子电路计算通过 pybind11 调用 C++ 后端。

---

## 一、三层架构总览

```
Layer 3: 宏观组件层 (Components/)
  ReactorCore / TFEUnit / Pipe / AnnularPipe / RingHP / HPwithFin / TECCircuitManager
  职责：组装底层物理实体 + 耦合器，暴露 get_solids()/get_couplers()

Layer 2: 求解器层 (Solvers/)
  SystemManager (顶层调度器)、HydraulicNetwork、HeatConduction、
  PointReactor (中子学)、Couplers (流固/固固/间隙耦合)
  职责：统一调度多物理场瞬态推进

Layer 1: 基础物理层
  Materials/ (物性)、Correlations/ (经验关联式)、MathSolvers/ (ODE/雅可比)、profiler.py (性能剖析)
  职责：提供材料物性、传热/流动关系式、数学求解工具、性能诊断
```

---

## 二、顶层入口与调度

### 2.1 入口文件

- **`main.py`**：当前为空文件，预留作为主程序入口。
- **`profiler.py`**：`TEASAProfiler` 轻量级性能剖析器（装饰器模式），被 Couplers.py 中的耦合器方法使用。
- 实际仿真启动通过 `testModule/` 和 `CoolantLoop/` 下的脚本完成，典型模式：

```python
from Components.ReactorCore import ReactorCore
from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import FluidVolume, FlowJunction
```

### 2.2 核心调度器：`Solvers/SystemManager.py`

`SystemManager` 是整个仿真的**全局协调中心**。

#### 内部持有对象

| 成员 | 类型 | 来源文件 |
|---|---|---|
| `fluid_solver` | `HydraulicNetwork` | `Solvers/Hydrodynamics/HydraulicNetwork.py` |
| `solid_components` | `Dict[str, BaseHeatConduction]` | `Solvers/HeatConduction/HeatConduction.py` |
| `couplers` | `List[Coupler]` | `Solvers/Couplers.py` |
| `point_reactor` | `PointReactor` (可选) | `Solvers/Neutronics/PointReactor.py` |
| `components` | `List[BaseComponent]` | `Components/BaseComponent.py` |

#### 注册流程 (Registration API)

```
SystemManager.add_component(component)
  ├── component.get_solids()     → add_solid_component()   注册固体导热求解器
  ├── component.get_couplers()   → add_coupler()           注册耦合器
  └── add_point_reactor()        → 注册点堆中子动力学模块
```

#### 时间步推进流程

```
SystemManager.step(dt)
  ├── pre_step():  调用所有 Component.pre_step(dt, t)    (更新功率、TEC等)
  ├── Picard 内迭代:
  │     ├── _prepare_fluid_sources_for_coupling()        清零流体源项 + 恢复持久源
  │     ├── _run_couplers()         同步所有耦合器边界
  │     ├── _solve_solids(dt)       推进固体导热 ODE
  │     ├── _solve_fluid(dt)        推进流体网络 (压力 + 焓)
  │     └── 收敛判断
  ├── 更新全局时间
  └── post_step(): 调用所有 Component.post_step(dt, t)
```

---

## 三、Solvers 层详细依赖

### 3.1 `Solvers/Hydrodynamics/` — 液力网络求解器

```
HydraulicNetwork.py
  ├── 依赖 Components.py       → FluidVolume, FlowJunction
  ├── 依赖 BoundaryVolume.py   → IncompressibleBoundaryVolume
  ├── 依赖 scipy.sparse.linalg.spsolve
  └── 被 SystemManager 调用
```

```
Components.py
  ├── FluidVolume         控制体类 (P, h, T, rho)
  ├── FlowJunction        连接件 (动量方程)
  ├── FluidChannel        流体通道 (向量化封装)
  ├── 依赖 Correlations.py → friction_single_phase
  └── 依赖 Materials/Base.py → FluidMaterial (类型提示)
```

```
BoundaryVolume.py → 压力边界 / 稳压器建模
```

**调用链**: `SystemManager` → `HydraulicNetwork` → `FluidVolume.get_volume_derivatives()` → `FluidMaterial`

### 3.2 `Solvers/HeatConduction/` — 固体导热求解器

```
HeatConduction.py
  ├── BaseHeatConduction       抽象基类 (1D/2D 共用)
  ├── HeatConduction1D         一维导热 (柱坐标)
  ├── HeatConduction2D         二维导热 (r-z 轴对称)
  ├── 依赖 Mesh.py             Mesh1D, Mesh2D
  ├── 依赖 Boundary.py         BoundaryRegion, ResistanceBC, FluxBC
  ├── 依赖 Materials/Base.py   SolidMaterial
  └── 使用 scipy.integrate.solve_ivp (BDF/Radau)
```

**调用链**: `SystemManager._solve_solids()` → `BaseHeatConduction.solve_step()` → `scipy.integrate.solve_ivp`

### 3.3 `Solvers/Couplers.py` — 多物理场耦合器

| 耦合器类 | 作用 | 连接双方 | 核心方法 |
|---|---|---|---|
| `FluidSolidCouple` | 流固对流换热 | `FluidChannel` ↔ `BoundaryRegion` | `sync()` |
| `SolidSolidCouple1D` | 一维固固导热 | `HeatConduction1D` ↔ `HeatConduction1D` | `sync()` |
| `SolidSolidCouple2D` | 二维固固导热 | `HeatConduction2D` ↔ `HeatConduction2D` | `sync()` |
| `GapCouple2D` | 气隙传热（导热+辐射） | `BoundaryRegion` ↔ `BoundaryRegion` | `sync()` |
| `TECCouple2D` | 热离子电极间隙 | 发射极 ↔ 接收极 | `sync()` |

**工作原理**: `sync()` 方法在每个 Picard 迭代中双向更新边界条件——将对方表面温度和内部热阻传递为 `ResistanceBC` 的外部参数。

### 3.4 `Solvers/Neutronics/PointReactor.py` — 点堆中子动力学

```
PointReactor
  ├── Numba JIT: _prke_rhs()  裂变功率 + 6组缓发中子 + 4组衰变热 = 11维ODE
  ├── Numba JIT: _prke_jac()  解析雅可比矩阵 (加速 BDF/Radau)
  ├── 依赖 scipy.integrate.solve_ivp
  └── 被 SystemManager 调用，通过 ReactorCore 提供反应性反馈
```

---

## 四、Components 层详细依赖

### 4.1 组件继承体系

```
BaseComponent (ABC, Components/BaseComponent.py)
  ├── 抽象方法: get_solids() → list
  ├── 抽象方法: get_couplers() → list
  ├── 生命周期钩子: pre_step(dt, t)
  └── 生命周期钩子: post_step(dt, t)
```

### 4.2 `Components/Pipe.py` — 单管道组件

```
Pipe(BaseComponent)
  ├── 持有: solid_wall (HeatConduction2D)
  ├── 持有: fluid_channel (FluidChannel)
  ├── 自动创建: FluidSolidCouple × 1
  └── get_solids() → [solid_wall]
      get_couplers() → [FluidSolidCouple]
```

### 4.3 `Components/AnnularPipe.py` — 环形管道组件

```
AnnularPipe(BaseComponent)
  ├── 持有: solid_inner (HeatConduction2D, 内管壁)
  ├── 持有: solid_outer (HeatConduction2D, 外管壁)
  ├── 持有: fluid_channel (FluidChannel)
  ├── 自动创建: FluidSolidCouple × 2
  └── get_solids() → [solid_inner, solid_outer]
      get_couplers() → [coupler_inner, coupler_outer]
```

### 4.4 `Components/HPwithFin.py` — 热管+降维翅片

```
HPwithFin(BaseComponent)
  ├── 内部创建: Mesh2D (轴向三段 + 径向吸液芯/管壁)
  ├── 持有: hp (HeatPipe2D, 热管本体)
  ├── 降维翅片模型: 沿翅高方向的一维准稳态导热
  ├── 等效外边界: ResistanceBC (翅片热阻) + DynamicRadiationBC (辐射)
  ├── 支持: ExternalHeatFluxBC (翅片直接受照轨道热流)
  └── 依赖: Components/basicComponents/HeatPipe2D.py
```

### 4.5 `Components/RingHP.py` — 集流环+代表性热管阵列

```
RingHP(BaseComponent)
  ├── 持有: fluid_channel (集流环流体通道)
  ├── 持有: solid_header (集流环壁, HeatConduction2D)
  ├── 自动创建: FluidSolidCouple (header 流固换热) × 1
  ├── 为每个流体节点创建: HPwithFin × N (代表性热管)
  ├── 自动创建: FluidSolidCouple (蒸发段换热) × N
  ├── 支持: ExternalHeatSource 配置 (轨道外热流)
  ├── 依赖: Components/HPwithFin.py
  ├── 依赖: Components/ExternalHeatSources/ (轨道热源)
  └── 依赖: Solvers/Couplers.py → FluidSolidCouple
```

### 4.6 `Components/TFEUnit.py` — 热离子燃料元件（最复杂组件）

```
TFEUnit(BaseComponent)
  ├── 内部创建多层次网格:
  │     ├── Fuel (UO2 芯块, HeatConduction2D)
  │     ├── Emitter (Mo 发射极, HeatConduction2D)
  │     ├── Collector (Nb 接收极, HeatConduction2D)
  │     ├── InnerClad (内套管, HeatConduction2D)
  │     ├── OuterClad (外套管, HeatConduction2D)
  │     └── Moderator (ZrH 慢化剂, HeatConduction2D)
  ├── 内部创建多层间隙耦合:
  │     ├── GapCouple2D: 裂变气隙 (Pellet-Emitter)
  │     ├── TECCouple2D: 铯气隙/极间隙 (Emitter-Collector)
  │     ├── GapCouple2D: 氦气隙 (Collector-InnerClad)
  │     └── GapCouple2D: CO2气隙 (OuterClad-Moderator)
  ├── 自动创建: SolidSolidCouple2D (固固直接接触面)
  ├── 自动创建: FluidSolidCouple (冷却剂通道换热)
  ├── 电场/热离子计算接口:
  │     ├── ElectricFieldData: 发射极/接收极电压、电流密度、焦耳热
  │     └── PlasmaCouplingData: 功函数、阻塞电压、电子冷却/加热热流
  ├── 依赖: Components/basicComponents/Fuel.py
  ├── 依赖: Components/basicComponents/Electord.py (Emitter, Collector)
  ├── 依赖: Solvers/Couplers.py → SolidSolidCouple2D, GapCouple2D, TECCouple2D, FluidSolidCouple
  ├── 依赖: Solvers/HeatConduction/HeatConduction.py → HeatConduction2D
  ├── 依赖: Correlations/Correlations.py → nu_ringpipe
  └── 依赖: Materials/Solids/CompositePelletMaterial.py
```

### 4.7 `Components/TECCircuitManager.py` — 全局电路管理器

```
TECCircuitManager(BaseComponent)
  ├── 持有: TECPair × N (电极对列表)
  ├── 持有: ThermoCalcModel (C++ 电路求解器封装)
  ├── pre_step(): 收集电极温度 → 求解全局电路 → 下发等离子体热流与焦耳热
  ├── get_solids() → 所有 Emitter + Collector
  ├── get_couplers() → 所有 TECCouple2D
  ├── 依赖: Components/basicComponents/TECPair.py
  ├── 依赖: ThermoCalc/ThermoCalcWrapper.py → ThermoCalcModel
  └── 依赖: Components/tec_electric.py → electric_field_from_node_potential
```

### 4.8 `Components/ReactorCore.py` — 堆芯宏观容器（顶层组装）

```
ReactorCore(BaseComponent)
  ├── 持有: TFEUnit × N (热离子燃料元件)
  ├── 持有: PointReactor (点堆中子动力学)
  ├── 持有: 全局结构:
  │     ├── 慢化剂环 (HeatConduction2D × N_ring)
  │     ├── 筒体 Barrel (HeatConduction2D)
  │     ├── 反射层 Reflector (HeatConduction2D)
  │     ├── 慢化剂-筒体间隙 (GapCouple2D)
  │     └── 筒体-反射层间隙 (GapCouple2D)
  ├── 功率分配: 将总裂变功率按 power_factors 分配给各 TFE
  ├── 反应性反馈: 计算燃料/电极/慢化剂/反射层平均温度 → 反馈反应性
  ├── TEC 更新: 协调 TECCircuitManager 的电热耦合
  ├── get_solids() → 所有 TFE 固体 + 全局结构固体
  ├── get_couplers() → 所有 TFE 耦合器 + 全局结构耦合器
  ├── 依赖: Components/TFEUnit.py
  ├── 依赖: Components/TECCircuitManager.py
  ├── 依赖: Solvers/Neutronics/PointReactor.py
  └── 依赖: ThermoCalc/ThermoCalcWrapper.py
```

### 4.9 `Components/basicComponents/` — 基础物理元件

| 文件 | 类 | 作用 | 依赖 |
|---|---|---|---|
| `Fuel.py` | `Fuel` | UO2 燃料芯块导热 | `HeatConduction2D`, `SolidMaterial` |
| `Electord.py` | `Emitter`, `Collector` | 发射极/接收极电极 | `HeatConduction2D`, `SolidMaterial` |
| `HeatPipe2D.py` | `HeatPipe2D` | 二维热管求解器 | `HeatConduction2D`, `FluidMaterial`, `SolidMaterial` |
| `FinConduction.py` | `FinConduction` | 翅片导热 | `HeatConduction1D` |
| `TECPair.py` | `TECPair` | TEC 电极对封装 | `Emitter`, `Collector`, `TECCouple2D` |

### 4.10 `Components/ExternalHeatSources/` — 轨道外热流

```
ExternalHeatSources/
  ├── embedded_flux_tables.py      预编译轨道热流查表数据 (FORTRAN 格式)
  └── __init__.py
        ├── BaseExternalHeatSource     抽象基类 → get_heat_flux(time) → [W/m²]
        ├── OrbitalHeatSource          太阳直射热流
        ├── AlbedoHeatSource           地球反照热流
        ├── EarthIRHeatSource          地球红外热流
        ├── CompositeHeatSource        组合多种热源
        ├── OrbitalTableHeatSource     查表模式热源
        └── ExternalHeatFluxBC         边界条件封装 → 挂载到 BoundaryRegion
```

**调用链**: `RingHP` / `HPwithFin` → `ExternalHeatFluxBC` → `BoundaryRegion.add_flux_condition()`

---

## 五、Materials 层 — 材料物性库

### 5.1 基类

```
Materials/Base.py
  ├── FluidMaterial  流体物性基类 (enthalpy, density, viscosity, heat_capacity, ...)
  └── SolidMaterial  固体物性基类 (conductivity, density, heat_capacity, ...)
```

### 5.2 流体工质

| 文件 | 材料 | 被调用者 |
|---|---|---|
| `Fluids/NaK78.py` | 钠钾合金 NaK-78 | `HydraulicNetwork` → `FluidVolume` |
| `Fluids/Sodium.py` | 液态钠 | `HydraulicNetwork` → `FluidVolume` |
| `Fluids/Potassium.py` | 液态钾 | `HydraulicNetwork` → `FluidVolume` |
| `Fluids/SodiumPotassium78.py` | 钠钾合金（别名） | 同上 |

### 5.3 固体材料

| 文件 | 材料 | 被调用者 |
|---|---|---|
| `Solids/UO2.py` | 二氧化铀燃料 | `Fuel` |
| `Solids/Molybdenum.py` | 钼 | `Emitter` |
| `Solids/MoNb.py` | 钼铌合金 | `Collector` |
| `Solids/ZrH.py` | 氢化锆慢化剂 | `TFEUnit`, `ReactorCore` |
| `Solids/B4C.py` | 碳化硼 | 吸收体 |
| `Solids/Beryllium.py` | 铍 | 反射层 |
| `Solids/BerylliumOxide.py` | 氧化铍 | 反射层/慢化剂 |
| `Solids/StainlessSteel.py` | 不锈钢 | 结构材料 |
| `Solids/WallMaterial.py` | 热管管壁材料 | `HPwithFin`, `HeatPipe2D` |
| `Solids/WickMaterial.py` | 吸液芯材料 | `HPwithFin`, `HeatPipe2D` |
| `Solids/WickStructure.py` | 吸液芯结构 | `HPwithFin`, `HeatPipe2D` |
| `Solids/KHP.py` | 钾热管工质 | `HeatPipe2D` |
| `Solids/NaHP.py` | 钠热管工质 | `HeatPipe2D` |
| `Solids/GasGaps.py` | 气隙等效物性 | `TFEUnit` (间隙建模) |
| `Solids/CompositePelletMaterial.py` | 复合芯块材料 | `TFEUnit` |

---

## 六、辅助模块

### 6.1 `Correlations/Correlations.py` — 经验关联式

| 函数 | 用途 | 被调用者 |
|---|---|---|
| `friction_single_phase()` | 单相摩擦系数 | `FluidVolume`, `FlowJunction` |
| `nu_ringpipe()` | 环形管道 Nu 数 | `TFEUnit` (FluidSolidCouple) |

### 6.2 `MathSolvers/` — 数学求解工具

| 文件 | 类 | 作用 | 被调用者 |
|---|---|---|---|
| `solver_module.py` | `NuclearODESolver` | 封装 `scipy.solve_ivp` 用于刚性 ODE | `SystemManager`, `PointReactor` |
| `optimization_utils.py` | `FluidJacobianBlockLayout` | 构建水力网络块状雅可比稀疏矩阵 | `SystemManager` (Jacobian 求解) |

### 6.3 `ThermoCalc/` — C++ 热离子电路求解器

```
ThermoCalc/
  ├── circuitTECs.cpp/h              C++ 电路求解核心 (串联/并联 TEC)
  ├── thermionicEmission.cpp/h       热离子发射物理模型
  ├── singleThermionicEnergyConversion.cpp/h  单管能量转换
  ├── NonLinerSolver.cpp/h           非线性方程求解器
  ├── bindings.cpp                   pybind11 Python 绑定
  ├── ThermoCalcWrapper.py           Python 封装层
  ├── CMakeLists.txt                 C++ 构建配置
  └── __init__.py                    包初始化

  调用链: TECCircuitManager → ThermoCalcWrapper → bindings.cpp → C++ solver
```

### 6.4 `profiler.py` — 性能分析

`TEASAProfiler` 装饰器，用于关键函数的计时分析。

| 被装饰的模块 / 类 |
|---|
| `SystemManager` |
| `BaseHeatConduction` |
| `HydraulicNetwork` |
| `TECCircuitManager` |
| `Couplers` |

### 6.5 `inputs/` — 输入参数文件

| 文件 | 内容 |
|---|---|
| `CoreInput.txt` | 堆芯参数 |
| `GloInput.txt` | 全局参数 |
| `HPInput.txt` | 热管参数 |
| `PipeInput.txt` | 管道参数 |
| `PipeNetInput.txt` | 管网参数 |
| `PumpInput.txt` | 泵参数 |
| `RadiatorInput.txt` | 散热器参数 |
| `steady.txt` | 稳态工况 |
| `transient.txt` | 瞬态工况 |
| `trip.txt` | 触发条件 |
| `density1-3.txt` | 密度分布 |
| `sheildtemp.txt` | 屏蔽层温度 |

---

## 七、完整依赖关系图

```
                        SystemManager (顶层调度)
                       /    |     |     \      \
                      /     |     |      \      \
            Hydraulic   HeatConduction  Couplers  PointReactor
            Network      (1D/2D)       (流固/固固/间隙)
           /   |   \        |             |           |
    FluidVolume  Junction   |             |    (Numba JIT)
    FlowJunction            |             |
         |                  |             |
    FluidMaterial    SolidMaterial   BoundaryRegion
    (NaK78/Sodium/   (UO2/Mo/ZrH/   (ResistanceBC/
     Potassium)       SS/Beryllium)   FluxBC/RadiationBC)
         |                  |             |
    [Correlations.py]  [Correlations.py]  [ExternalHeatSources]
    (friction/nu)      (conductivity)     (轨道外热流)


    ReactorCore ──→ TFEUnit × N ──→ Fuel / Emitter / Collector
        │               │                (basicComponents)
        │               ├── GapCouple2D × 4 (裂变气隙/极间隙/氦气隙/CO2气隙)
        │               ├── TECCouple2D × 1
        │               └── FluidSolidCouple (冷却剂换热)
        │
        ├── TECCircuitManager ──→ TECPair × N ──→ ThermoCalc (C++/pybind11)
        │
        ├── 全局慢化剂环 (HeatConduction2D × N_ring)
        ├── 筒体 Barrel (HeatConduction2D)
        ├── 反射层 Reflector (HeatConduction2D)
        ├── 层间间隙 (GapCouple2D)
        └── PointReactor (中子学反馈: 温度 → 反应性)


    Pipe / AnnularPipe ──→ FluidSolidCouple + HeatConduction2D
    RingHP ──→ HPwithFin × N + FluidSolidCouple + ExternalHeatSources
    HPwithFin ──→ HeatPipe2D + 翅片准稳态模型 + ExternalHeatSources
```

---

## 八、关键调用关系速查表

| 被调用模块 | 调用者模块 |
|---|---|
| `Solvers/SystemManager.py` | `testModule/test_core_assemble_v*.py`, `CoolantLoop/run_*.py` |
| `Solvers/Hydrodynamics/HydraulicNetwork.py` | `SystemManager` |
| `Solvers/Hydrodynamics/Components.py` | `HydraulicNetwork`, `Pipe`, `AnnularPipe`, `RingHP` |
| `Solvers/Hydrodynamics/BoundaryVolume.py` | `HydraulicNetwork` |
| `Solvers/HeatConduction/HeatConduction.py` | `Pipe`, `AnnularPipe`, `TFEUnit`, `ReactorCore`, `HPwithFin`, `RingHP` |
| `Solvers/HeatConduction/Mesh.py` | `HeatConduction`, `TFEUnit`, `HPwithFin`, `ReactorCore` |
| `Solvers/HeatConduction/Boundary.py` | `HeatConduction`, `Couplers`, `HPwithFin`, `RingHP` |
| `Solvers/Couplers.py` | `Pipe`, `AnnularPipe`, `TFEUnit`, `ReactorCore`, `RingHP` |
| `Solvers/Neutronics/PointReactor.py` | `SystemManager` (通过 `ReactorCore`) |
| `Components/BaseComponent.py` | 所有组件类 (Pipe, TFEUnit, ReactorCore 等) |
| `Components/ReactorCore.py` | `test_core_assemble_v*.py` |
| `Components/TFEUnit.py` | `ReactorCore` |
| `Components/TECCircuitManager.py` | `ReactorCore` |
| `Components/RingHP.py` | `CoolantLoop/model_collector_ring_*.py` |
| `Components/HPwithFin.py` | `RingHP` |
| `Components/Pipe.py` | 各种管道测试/组装脚本 |
| `Components/AnnularPipe.py` | 冷却回路模型 |
| `Components/basicComponents/Fuel.py` | `TFEUnit` |
| `Components/basicComponents/Electord.py` | `TFEUnit`, `TECPair` |
| `Components/basicComponents/HeatPipe2D.py` | `HPwithFin` |
| `Components/basicComponents/TECPair.py` | `TECCircuitManager` |
| `Components/basicComponents/FinConduction.py` | `HPwithFin` |
| `Components/ExternalHeatSources/` | `RingHP`, `HPwithFin` |
| `Components/tec_electric.py` | `TECCircuitManager`, `TFEUnit` |
| `ThermoCalc/ThermoCalcWrapper.py` | `TECCircuitManager`, `ReactorCore` |
| `Materials/Base.py` | 所有 `Materials/Fluids/*.py`, `Materials/Solids/*.py`, `HeatConduction`, `FluidVolume` |
| `Materials/Fluids/NaK78.py` | `HydraulicNetwork` → `FluidVolume` |
| `Materials/Solids/UO2.py` | `Fuel`, `TFEUnit` |
| `Materials/Solids/Molybdenum.py` | `Emitter`, `TFEUnit` |
| `Materials/Solids/MoNb.py` | `Collector`, `TFEUnit` |
| `Materials/Solids/ZrH.py` | `TFEUnit`, `ReactorCore` |
| `Materials/Solids/WallMaterial.py` | `HPwithFin`, `HeatPipe2D` |
| `Materials/Solids/WickMaterial.py` | `HPwithFin`, `HeatPipe2D` |
| `Materials/Solids/WickStructure.py` | `HPwithFin`, `HeatPipe2D` |
| `Materials/Solids/KHP.py` | `HeatPipe2D` |
| `Materials/Solids/NaHP.py` | `HeatPipe2D` |
| `Materials/Solids/GasGaps.py` | `TFEUnit` |
| `Materials/Solids/CompositePelletMaterial.py` | `TFEUnit` |
| `Correlations/Correlations.py` | `FluidVolume`, `FluidSolidCouple` |
| `MathSolvers/solver_module.py` | `SystemManager`, `PointReactor` |
| `MathSolvers/optimization_utils.py` | `SystemManager` |
| `profiler.py` | `SystemManager`, `HeatConduction`, `HydraulicNetwork`, `TECCircuitManager`, `Couplers` |

---

## 九、测试与运行脚本

### 9.1 `testModule/` — 单元测试与集成测试

| 测试脚本 | 测试对象 |
|---|---|
| `test_core_assemble_v1~v7*.py` | 堆芯组装 (逐步迭代版本 v1→v7) |
| `test_component_pipe.py` | Pipe 组件 |
| `test_component_annular_pipe.py` | AnnularPipe 组件 |
| `test_HP_with_external_heat_source.py` | 热管外热流 |
| `testHPwithFin.py` | 带翅片热管 |
| `test_ringHP_buzzin.py` | 环形热管 |
| `test_TEC_with_heat.py` | TEC 电热耦合 |
| `test_tec_joule_nonuniform.py` | TEC 焦耳热非均匀分布 |
| `test_tecboundary.py` | TEC 边界条件 |
| `test_thermionic_wrapper.py` | ThermoCalc 封装 |
| `test_tfe_thermal_flow.py` | TFE 热工水力 |
| `test_system_manager.py` | SystemManager |
| `test_system_manager_lifecycle.py` | SystemManager 生命周期 |
| `test_flow_heat.py` | 流动换热 |
| `test_fluid_solid_couple.py` | 流固耦合 |
| `test_simple_solid_fluid_couple.py` | 简化流固耦合 |
| `test_open_channel.py` | 开式通道 |
| `test_open_channel_transient.py` | 开式通道瞬态 |
| `test_single_channel_transient.py` | 单通道瞬态 |
| `test_parallel_channels.py` | 并联通道 |
| `test_pressurizer_volume.py` | 稳压器 |
| `test_pressurizer_pumped_closed_loop.py` | 泵驱闭式回路 |
| `test_pump_junction_hydraulic_network.py` | 泵节点水力网络 |
| `test_PointReactor.py` | 点堆动力学 |
| `testHeadConduction1D.py` | 一维导热 |
| `testHeadConduction2D.py` | 二维导热 |
| `test_verify_solid_couplers.py` | 固固耦合器验证 |

### 9.2 `CoolantLoop/` — 冷却回路模型

| 脚本 | 说明 |
|---|---|
| `model_collector_ring_6segment.py` | 6 段集流环模型定义 |
| `model_collector_ring_full_ringhp.py` | 含完整 RingHP 的集流环模型 |
| `run_collector_ring_6segment_160s.py` | 6 段模型运行 160s |
| `run_collector_ring_full_ringhp_200s_resume.py` | 完整模型续算 200s |
| `run_collector_ring_full_ringhp_500s_resume.py` | 完整模型续算 500s |
| `test_coolant_loop_v4*.py` | 冷却回路 v4 版本测试 |
| `test_coolant_loop_v5*.py` | 冷却回路 v5 版本测试 |
| `test_full_collector_ring.py` | 完整集流环测试 |
| `test_single_header_cell_one_hp.py` | 单集流环-单热管测试 |
| `verify_dudt_long_tail_case.py` | 长尾衰减验证 |

### 9.3 优化重构目录

| 目录 | 内容 |
|---|---|
| `HeatConduction优化重构/` | `benchmark_systemmanager_heatconduction.py`：导热求解器性能基准 |
| `HeatPipe优化重构/` | `benchmark_ringhp_v5.py`：热管 v5 性能基准 |
| `HydraulicNetwork优化重构/` | `benchmark_open_loop_complex_300cv.py`：水力网络性能基准 |
| `SystemManager优化重构/` | `validate_inner_iter_picard.py`：Picard 内迭代收敛验证 |