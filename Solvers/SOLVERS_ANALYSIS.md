# TASTIN Solvers 层深度解析

> 本文件是跨模块深度参考，不是 Codex 首次进入 `Solvers/` 时的默认必读入口。默认先读 `AI_AGENT_SOLVERS_ANALYSIS.md`，再按任务选择模块文档。

## 目录

1. [概述](#1-概述)
2. [HydraulicNetwork 水力网络求解器](#2-hydraulicnetwork-水力网络求解器)
   - [2.1 FluidVolume / FlowJunction / FluidChannel (Components.py)](#21-fluidvolume--flowjunction--fluidchannel-componentspy)
   - [2.2 IncompressibleFluidVolume / PressurizerVolume](#22-incompressiblefluidvolume--pressurizervolume)
   - [2.3 FixedPressurePumpVolume / PumpJunction](#23-fixedpressurepumpvolume--pumpjunction)
   - [2.4 MacroFlowJunction (宏观看管)](#24-macroflowjunction-宏观看管)
   - [2.5 BoundaryVolumes (BoundaryVolume.py)](#25-boundaryvolumes-boundaryvolumepy)
   - [2.6 HydraulicNetwork 总体逻辑](#26-hydraulicnetwork-总体逻辑)
3. [HeatConduction 固体导热求解器](#3-heatconduction-固体导热求解器)
   - [3.1 Mesh1D / Mesh2D (Mesh.py)](#31-mesh1d--mesh2d-meshpy)
   - [3.2 BoundaryRegion / BoundaryConditions (Boundary.py)](#32-boundaryregion--boundaryconditions-boundarypy)
   - [3.3 BaseHeatConduction / HeatConduction1D / HeatConduction2D](#33-baseheatconduction--heatconduction1d--heatconduction2d)
4. [Couplers 多物理场耦合器](#4-couplers-多物理场耦合器)
   - [4.1 FluidSolidCouple (流固耦合)](#41-fluidsolidcouple-流固耦合)
   - [4.2 SolidSolidCouple1D / SolidSolidCouple2D (固固耦合)](#42-solidsolidcouple1d--solidsolidcouple2d-固固耦合)
5. [PointReactor 点堆中子动力学](#5-pointreactor-点堆中子动力学)
6. [SystemManager 全局调度器](#6-systemmanager-全局调度器)
   - [6.1 step() 方法完整流程](#61-step-方法完整流程)
7. [关键数据流图](#7-关键数据流图)
8. [源码校验信息](#8-源码校验信息)

***

## 1. 概述

TASTIN 的 `Solvers/` 层包含所有物理求解器，按功能分为四个子系统：

| 子系统   | 目录                | 核心文件                  | 职责                              |
| ----- | ----------------- | --------------------- | ------------------------------- |
| 水力网络  | `Hydrodynamics/`  | `HydraulicNetwork.py` | 求解冷却剂流体的压力-流量-温度瞬态分布            |
| 固体导热  | `HeatConduction/` | `HeatConduction.py`   | 求解固体域 (燃料/电极/管壁等) 的 1D/2D 瞬态温度场 |
| 中子动力学 | `Neutronics/`     | `PointReactor.py`     | 点堆动力学 (裂变功率 + 缓发中子 + 衰变热)       |
| 全局调度  | `Solvers/`        | `SystemManager.py`    | Picard 内迭代调度所有子求解器 + 耦合器        |

***

## 2. HydraulicNetwork 水力网络求解器

### 2.1 FluidVolume / FlowJunction / FluidChannel (Components.py)

#### FluidVolume — 控制体类 (标量容器)

**物理职责：存储 P (压力)、h (焓)、T (温度)、rho (密度)，求解质量/能量守恒。**

**核心方法：`get_volume_derivatives(material) -> (dPdt, dhdt)`**

```python
# 方程组 (2×2 解耦):
# Eq1 (质量守恒): drho_dP * dPdt + drho_dh * dhdt = (ΣW_in - ΣW_out) / V
# Eq2 (能量守恒): rho * dhdt = (Σ(W*h_net) + Q_total) / V
```

**关键几何属性：**

- `vol` — 节点体积 \[m³]
- `len` — 节点长度 \[m]
- `area` — 流通面积 \[m²]
- `d_h` — 水力直径 \[m]

**关键源项接口：**

- `Q_wall` — 壁面热流 \[W] (正值加热流体)
- `Q_vol` — 体积热源 \[W]
- `implicit_coeff` — 隐式换热系数 λ \[W/K]，用于半隐式处理对流换热

**拓扑连接：**

- `inlet_junctions: List[FlowJunction]` — 流入本节点的连接
- `outlet_junctions: List[FlowJunction]` — 流出本节点的连接

**能量通量采用施主元胞 (Donor Cell) 策略：**

- 入口连接：`h_donor = from_vol.h if W >= 0 else self.h`
- 出口连接：`h_donor = self.h if W >= 0 else to_vol.h`

**人工可压缩性：**

- 设置虚拟声速 `c_art = 200.0 m/s`，`drho_dP_artificial = 1/c²`
- 取 `min(drho_dP_artificial, denom_P)` 作为最终压缩性，确保数值稳定

***

#### FlowJunction — 连接件 (矢量通道)

**物理职责：连接两个控制体，求解动量守恒，计算 W (质量流量)。**

**核心方法：`get_momentum_derivative(material) -> dWdt`**

**动量方程：**

```
(L/A) * dW/dt = (P_in - P_out) + Gravity + Pump - Friction - FormLoss
```

**阻力计算分为两类：**

1. **摩擦阻力** (`calculate_friction_pressure_drop`): Darcy-Weisbach 公式，调用 `Correlations.friction_single_phase(Re)`
2. **形阻** (`calculate_form_loss_pressure_drop`): K 系数法；配置 `dynamic_loss_params` 后可叠加动态阻力模型

**动态阻力模型 (`_compute_dynamic_k_loss`)：**

- 当 `dynamic_loss_params['model'] != 'none'` 时激活：
  - `'single_crossflow_pipe'` — 集流环内横掠管束模型
    - 计算最窄截面流速 `V_max = V_nominal * A_flow / A_min`
    - 根据 Re 数分段计算阻力系数 C\_D
    - 换算为等效 K\_eq = C\_D \* (A\_proj / A\_flow) \* (A\_flow / A\_min)²
- 动态模型启用时，负的基础 `k_loss` 先按 `0` 处理，再叠加动态阻力。负值本身不是模型开关。

**几何属性：**

- `area` — 连接截面积 (默认取两端较小值)
- `length` — 惯性长度 (默认 = 0.5\*(from\_vol.len + to\_vol.len))

**拓扑自动注册：** 构造函数中自动将 `self` 添加到 `from_vol.outlet_junctions` 和 `to_vol.inlet_junctions`

***

#### FluidChannel — 流体通道 (复合组件)

**物理职责：自动化创建并连接一系列 FluidVolume 和 FlowJunction，形成一维离散通道。**

**构造函数流程：**

```
FluidChannel.__init__(name, n_nodes, total_length, flow_area, hydraulic_diam, ...)
  ├── 计算 node_length = total_length / n_nodes
  ├── 计算 node_volume = flow_area * node_length
  ├── for i in range(n_nodes):
  │     └── 创建 FluidVolume (均匀网格)
  └── for i in range(n_nodes-1):
        └── 创建 FlowJunction (连接相邻节点)
```

**向量化访问接口 (供 Coupler 使用)：**

- `temperature_vector` — np.array(\[vol.T for vol in volumes])
- `pressure_vector` — np.array(\[vol.P for vol in volumes])
- `density_vector` — np.array(\[vol.rho for vol in volumes])
- `velocity_vector` — 节点入口/出口平均流速

**批量源项操作：**

- `add_coupling_source_distribution(explicit_arr, implicit_arr)` — 逐节点累加耦合源项
- `add_heat_source_distribution(heat_dist)` — 逐节点累加热源
- `clear_sources()` — 每步开始时清零所有源项

***

### 2.2 IncompressibleFluidVolume / PressurizerVolume

#### IncompressibleFluidVolume

**设计目的：** 实现压力-流量解耦的代数压力法。

**关键差异 (相比 FluidVolume)：**

- `get_volume_derivatives` 仅返回 `dydt_h` (float)，不再返回 `dPdt`
- 压力 P 由外部通过 `update_pressure_algebraic(new_P)` 代数更新
- 物性 (T, rho, mu) 在 `update_pressure_algebraic` 中同步刷新

**能量方程采用非守恒输运形式：**

```
rho * V * dh/dt = Σ(W_in * (h_donor - h_self)) + Q_total
```

此形式消除了质量不平衡带来的误差，比原始形式更适合不可压缩流动。

***

#### PressurizerVolume (被动压力参考点)

**作用：** 封闭回路的绝对压力锚点，但不作为质量源/汇。

**关键特性：**

- `is_pressure_reference = True` — 标记为压力参考点
- `is_pressure_boundary = False` — 不是固定压力边界
- `set_pressure(pressure)` — 设置固定目标压力
- `set_pressure_table(times, pressures)` — 设置时间-压力插值表

***

### 2.3 FixedPressurePumpVolume / PumpJunction

#### FixedPressurePumpVolume

**作用：** 提供固定压升 ΔP 的泵控制体。

- `delta_p` — 泵压升 \[Pa]
- 在 `update_pressure_distribution_downstream` 中，`pressure_rise` 被叠加到递推公式

#### PumpJunction

**作用：** 在 FlowJunction 的动量方程中叠加泵扬程。

- `compute_pump_head(time)` — 返回当前泵扬程，支持时间和压力表格插值
- `is_pump_junction = True` — HydraulicNetwork 据此识别并在动量系数计算中叠加 dP\_pump

***

### 2.4 MacroFlowJunction (宏观看管)

**物理职责：** 处理"代表管"模型中的**跨尺度缩放**。

**核心原理：** 系统中可能用一根管子代表 N 根相同的并联管。在动量方程中保持单管物理 (1×)，但在质量/能量守恒中对宏观端施加乘子 M。

**关键属性：**

- `macro_vol` — 被放大的宏观端 (如联箱)
- `multiplier` — 放大倍数 N
- `multiplier_from` / `multiplier_to` — 根据流向自动判断

**核心方法：`get_mass_flow_for(target_vol) -> float`**

```python
if target_vol == self.from_vol:
    return self.W * self.multiplier_from  # 宏观端看到 N×W
elif target_vol == self.to_vol:
    return self.W * self.multiplier_to
```

***

### 2.5 BoundaryVolumes (BoundaryVolume.py)

#### BoundaryVolume (可压缩版本)

**作用：** 固定压力/温度边界的瞬态松弛驱动。

**核心修改：** 不再直接强制修改 P 和 T，而是返回导数驱动：

- `dP/dt = K * (P_target - P_current)`，默认 `K = 50.0`
- `dh/dt = K * (h_target - h_current)`

#### IncompressibleBoundaryVolume

**作用：** 不可压缩流体的压力参考点 + 温度边界。

**两种模式：**

1. **无限大热池** (`mixing_enabled=False`): 温度恒定，dh/dt = 0
2. **有限体积混合** (`mixing_enabled=True`): 基于流入流出能量平衡计算 dh/dt

**接口：**

- `set_boundary_state(P, T)` — 动态修改边界条件
- `forcing_dydt_h` — 外部强制变温速率

#### InletJunction

**作用：** 质量流量边界，强制执行目标流量。

**工作原理：** 松弛法驱动

```python
dW/dt = K * (W_target - W_current)， K = 0.1
```

***

### 2.6 HydraulicNetwork 总体逻辑

**时间传播：** `HydraulicNetwork.set_time(time)` 保存网络时间，并向所有提供 `set_time()` 的控制体和连接传播。当前主要用于刷新 `PressurizerVolume` 的压力参考目标和 `PumpJunction` 的时间表泵压升。

#### Phase 1: 拓扑构建 (`_build_topology()`)

```
1. 建立 vol_to_idx 映射表
2. 提取并缓存所有节点的几何常数 (V_vec, L_node_vec, A_node_vec, Dh_node_vec, z_vec)
3. 识别定压边界 (fixed_pressure_indices) 和压力参考点 (pressure_reference_idx)
4. 解析所有连接关系，缓存连接几何常数
5. 提取 MacroFlowJunction 的乘子 (M_from_vec, M_to_vec)
```

#### 预构建矩阵结构缓存

**`_build_pressure_system_cache()`：** 预计算压力方程稀疏矩阵 CSR 结构

- 对角线：每个节点与自身的关系
- 非对角线：每个 Junction 产生的 from↔to 连接
- 定压边界行被特殊标记，在求解时固定为狄利克雷条件

**`_build_energy_system_cache()`：** 预计算能量方程稀疏矩阵结构

- 类似压力矩阵，用于半隐式焓方程求解

#### Phase 2: 物性更新 (`_update_fluid_properties()`)

**向量化批量计算：**

```python
cp_vec = material.heat_capacity(T_vec, P_vec)       # 热容
rho_vec = material.density(T_vec, P_vec)             # 密度
mu_vec = material.viscosity(T_vec, P_vec)             # 粘度
drho_dp_vec = material.liquid_density_derivative_P(P_vec)  # 压缩性
drho_dt_vec = material.liquid_density_derivative_T(T_vec)  # 热膨胀
```

所有计算使用单次 NumPy 向量化调用，避免 Python 循环。

#### Phase 3: 动量系数计算 (`_calc_momentum_coeffs_fast()`)

**物理公式：**

```
W^{n+1} = a_j * (P_in - P_out) + b_j

a_j = 1 / (I_term + K_linear)           # 导纳
b_j = a_j * (I_term*W_old + Grav + Acc + Pump)  # 源项

I_term = L_inertial / (A_flow * dt)     # 惯性项
K_linear = (f*L/D + K_loss)*|W|/(2*rho*A²)  # 线性化阻力
```

**摩擦系数** **`_calc_friction_factor_static_vec(Re)`：**

- Re ≤ 1000: 层流 f = 64/Re
- 2300 < Re < 1e5: Blasius f = 0.3164/Re^0.25
- Re ≥ 1e5: Karman-Prandtl 隐式迭代

**处理入口边界：**

```python
W_new = W_old + dt*K*(W_target - W_old)
```

#### Phase 4: 压力求解 (`step_Picard()`)

通过组装并求解稀疏线性方程组 A\*P = B，得到当前步压力分布。

#### Phase 5: 能量求解

通过半隐式方法求解焓方程，得到 T\_vec 和 h\_vec，完成水力网络的时间步推进。

***

## 3. HeatConduction 固体导热求解器

### 3.1 Mesh1D / Mesh2D (Mesh.py)

#### Mesh1D

**支持的几何类型：**

- `'cartesian'` — 笛卡尔平板
- `'cylindrical'` — 圆柱坐标 (燃料棒/管道)

**构造函数方法：**

1. **标准均匀网格：** `Mesh1D(total_dim, n_volumes, geometry_type, inner_radius, height)`
2. **自定义非均匀网格：** `Mesh1D.from_custom_faces(face_locations, geometry_type, height)`

**计算的几何数据 (存储于 self 属性)：**

- `node_centers` — 节点中心坐标 \[N]
- `face_locations` — 界面坐标 \[N+1]
- `volumes` — 控制体体积 \[N] (圆柱: π(r\_outer² - r\_inner²)\*H)
- `face_areas` — 界面面积 \[N+1] (圆柱: 2π*r*H)
- `dr_node_to_node` — 相邻节点距离 \[N-1]
- `dr_node_to_face` — 节点到自身上下界面的距离 \[N, 2]

***

#### Mesh2D

**支持的几何类型：** `'cartesian'` 或 `'cylindrical'` (r-z 轴对称)

**构造函数方法：**

1. **标准均匀网格：** `Mesh2D(x_dim, n_x, y_dim, n_y, geometry_type, inner_radius)`
2. **自定义非均匀网格：** `Mesh2D(..., x_faces=custom, y_faces=custom)`

**关键设计：** 所有数据扁平化为 1D 数组存储，但保留 2D 视图

- `shape_nodes = (n_x, n_y)` — 节点形状
- `shape_faces_x = (n_x+1, n_y)` — X 方向界面形状
- `shape_faces_y = (n_x, n_y+1)` — Y 方向界面形状

**GeometricData2D 容器：** 包含 volumes、node\_centers\_x/y、area\_x/y、距离数组等，全部为扁平数组。

***

### 3.2 BoundaryRegion / BoundaryConditions (Boundary.py)

#### BaseBoundaryCondition (抽象基类)

**核心方法：** `compute_flux_from_node(T_node, R_internal) -> flux`

#### ResistanceBC (热阻型 Robin 边界)

**公式：** `Flux = (T_ext - T_node) / (R_internal + R_ext + R_add)`

**三个参数：**

- `T_ext` — 外部温度 \[K]
- `R_ext` — 外部热阻 \[K/W]
- `R_add` — 附加热阻 (污垢/接触) \[K/W]

#### FluxBC (固定热流 Neumann 边界)

**公式：** `Flux = q_flux` (直接返回设定热流)

#### DynamicRadiationResistanceBC (动态非线性辐射边界)

**核心设计：** 在每个残差评估步动态计算 T⁴ 辐射热流

**线性化策略：**

```
h_rad = ε * σ * (T_surf + T_env) * (T_surf² + T_env²)
G_rad = h_rad * area
R_ext = 1 / G_rad
Flux = (T_env - T_node) / (R_internal + R_ext)
```

通过将非线性的 T⁴ 项局部线性化为 `h_rad * (T_env - T)`，可以在隐式 ODE 求解器中使用热阻网络法。

**`update_state(T_node, T_surface)`** **方法：** 每次残差评估时被 `BoundaryRegion.update_internal_state` 调用，重新计算 `h_rad`、`G_rad`、`R_ext`、`T_ext`。零长度边界数组直接返回，不调用 `min()`。

***

#### BoundaryRegion (边界区域管理器)

**核心原理：诺顿等效 (Parallel Conductance Superposition)**

多个边界条件可以叠加在同一个边界面上：

```
G_total = Σ(1/R_i)           # 总热导 = 各通道热导之和
J_total = Σ(T_i/R_i) + ΣQ_j  # 总诺顿源 = 温度驱动源 + 热流源

R_eff = 1 / G_total          # 等效热阻
T_eff = J_total / G_total    # 等效戴维南温度
```

**核心方法：**

1. **`update_internal_state(T_node, R_int, current_time)`** — 从固体求解器推送内部节点温度和内部热阻
2. **`compute_net_flux_for_solver()`** — 计算叠加后的净热流和表面温度
   ```
   for each bc in conditions:
     if resistance: 累加 G_sum += 1/R_ext, J_sum += T_ext/R_ext
     if flux:       累加 Q_sum_flux += q_flux

   R_eff = 1/G_sum
   T_eff = J_sum/G_sum + Q_sum_flux * R_eff  # 戴维南等效
   Flux = (T_eff - T_node) / (R_eff + R_internal)
   T_surface = T_node + Flux * R_internal
   ```
3. **`add_resistance_condition(T_ext, R_ext, R_add)`** — 添加热阻型边界条件，返回 `ResistanceBC` 对象
4. **`add_dynamic_radiation_condition(emissivity, bare_area_array, T_env)`** — 添加动态辐射条件
5. **`get_coupling_snapshot()`** — 返回 `(T_adj_node, R_internal)`，供 Coupler 获取固体侧状态
6. **`get_coupling_surface_snapshot()`** — 返回 `(T_surface, R_internal)`

***

### 3.3 BaseHeatConduction / HeatConduction1D / HeatConduction2D

#### BaseHeatConduction (抽象基类)

**核心状态：**

- `T` — 节点温度数组 \[N]
- `dTdt` — 温度导数 \[N]
- `Q_source` — 内热源 \[W]
- `boundaries` — 边界管理器字典

**核心接口：**

1. **`get_derivatives(t, T_current) -> dTdt`** — ODE 求解器接口，返回 RHS
   ```
   T = T_current
   _update_properties()         → k, rho, cp, thermal_capacitance
   _compute_internal_resistance()   → 内部界面热阻/热导
   _update_boundaries_state(t)     → 推送状态到边界
   Q_net = _compute_fluxes(t)      → 净导热热流
   _update_sources(t)              → 更新内热源
   dTdt = (Q_net + Q_source) / thermal_capacitance
   ```
2. **`step(dt, method) -> bool`** — 执行 ODE 积分一个时间步
   ```
   solve_ivp(fun=get_derivatives, t_span=[t, t+dt], y0=T, method=solid.ode_method, ...)
   ```
3. **`save_state_dict(prefix)`** **/** **`load_state_dict(data, prefix)`** — 持久化接口
4. **`link_source_buffer(external_buffer)`** — 内存绑定外部热源数组 (高性能耦合)
5. **`set_source_term(source_func)`** — 设置内热源回调函数

***

#### HeatConduction1D (一维导热)

**边界初始化：**

- `'inner'` — 左/内边界 (index 0)
- `'outer'` — 右/外边界 (index N-1)

**`_compute_internal_resistance()`** — 计算几何热阻

圆柱坐标使用精确对数热阻：

```
R_inner = ln(r_next / r_i) / (2 * π * k_interface * H)
R_boundary_inner = |ln(r_node / r_face)| / (2 * π * k * H)
```

笛卡尔坐标使用线性热阻：

```
R = dx / (k * area)
```

**`_compute_fluxes(t)`** — 计算各界面热流

```
内部界面: Flux_i = (T_i - T_i+1) / R_geom_i
内边界:   Flux_0 = Boundary.compute_net_flux_for_solver()
外边界:   Flux_N = -Boundary.compute_net_flux_for_solver()
净热流:   Q_net_i = Flux_i - Flux_i+1       (节点能量平衡)
```

***

#### HeatConduction2D (二维导热)

**边界初始化 (四个方向)：**

- `'left'` — X=0 边界，shape=(ny,)
- `'right'` — X=nx 边界，shape=(ny,)
- `'bottom'` — Y=0 边界，shape=(nx,)
- `'top'` — Y=ny 边界，shape=(nx,)

**`_compute_internal_resistance()`** — 直接计算热导 G = 1/R

X 方向 (径向)：

```
k_interface = 2*k_left*k_right / (k_left + k_right)  # 调和平均
G_x = k_interface * area_x / dx
```

Y 方向 (轴向)：

```
G_y = k_interface * area_y / dy
```

采用热导而非热阻，避免 k=0 时的除零错误。

**`get_jac_sparsity()`** — 返回 2D 五点差分格式的稀疏雅可比矩阵模式，用于 BDF/Radau 隐式求解器。

***

## 4. Couplers 多物理场耦合器

### 4.1 FluidSolidCouple (流固耦合)

**物理功能：**

1. 计算对流换热系数：Nu → h → λ = h\*A
2. 更新固体边界 (ResistanceBC)：`T_ext = T_fluid`, `R_ext = 1/λ`
3. 更新流体源项 (半隐式)：`Q = λ*T_wall - λ*T_fluid`

**`execute(interface_relaxation)`** **方法流程：**

```
1. 获取向量化状态: T_f, P_f, rho, vel, mu, k_f, Cp_f
2. 计算 Re = rho*|vel|*Dh / mu
3. 计算 Pr (调用 material.prandtl_number)
4. 计算 Nu = correlation_func(Re, Pr, P_D_ratio)
5. 计算 h = Nu * k_f / Dh
6. 计算 λ = h * A_node
7. 更新固体 BC: R_convection = 1/λ
8. 获取固体壁温 T_wall
9. 更新流体源项:
     explicit = λ * T_wall
     implicit = λ
     fluid.add_coupling_source_distribution(explicit, implicit)
10. 可选：界面松弛 (interface_relaxation < 1.0)
```

**自适应步长：** `get_max_stable_dt(safety_factor, max_limit)`

```
dt_max = min(C_solid / λ) * safety_factor
```

确保显式流固耦合不破坏 CFL 条件。

***

### 4.2 SolidSolidCouple1D / SolidSolidCouple2D (固固耦合)

**工作原理：** `sync()` 方法在每次 Picard 迭代中交叉更新边界条件

**1D 耦合器：**

```python
def sync():
    # 获取对方快照
    T1_surf, R1_int = bound1.get_coupling_snapshot()  # shape=(1,)
    T2_surf, R2_int = bound2.get_coupling_snapshot()

    # 交叉更新
    bc1.update_params(T_ext=T2_surf, R_ext=R2_int)
    bc2.update_params(T_ext=T1_surf, R_ext=R1_int)
```

**2D 耦合器：**

- 额外校验边界维度匹配
- 支持向量化整条边界的同步更新
- 方向映射表：`{'right': ('right','left'), 'left': ('left','right'), 'top': ('top','bottom'), 'bottom': ('bottom','top')}`

***

### 4.3 GapCouple2D (二维间隙耦合器)

**继承关系：** `GapCouple2D` → `SolidSolidCouple2D`

**物理功能：** 在 SolidSolidCouple2D 基础上，增加**间隙热阻**，同时考虑两种并联传热机制：

1. **表面对表面辐射** (Surface-to-Surface Radiation)
2. **间隙气体导热** (Gap Gas Conduction)

**构造函数参数：**

| 参数                 | 类型               | 含义                         |
| ------------------ | ---------------- | -------------------------- |
| `obj1, obj2`       | HeatConduction2D | 两侧固体                       |
| `direction`        | str              | 方向 (right/left/top/bottom) |
| `gap_width`        | float            | 间隙宽度 δ \[m]                |
| `gas_conductivity` | float            | 气体导热系数 k\_gas \[W/(m·K)]   |
| `emissivity1`      | float            | obj1 侧发射率 ε₁ (默认 0.8)      |
| `emissivity2`      | float            | obj2 侧发射率 ε₂ (默认 0.8)      |

**`sync()`** **方法流程：**

```
1. 获取表面温度 T1_surf, T2_surf (get_coupling_surface_snapshot)
2. 几何判定：自动识别内外表面 (面积小的为内表面 A_in)
3. 辐射电导计算：
     h_rad* = σ × (T₁² + T₂²) × (T₁ + T₂)
     分母 = 1/ε_in + (1/ε_out - 1) × (A_in / A_out)
     G_rad = (A_in × h_rad*) / 分母 [W/K]
4. 气体导热电导：
     G_cond = k_gas × A_in / δ [W/K]
5. 总热导 = G_rad + G_cond (并联)
6. 总间隙热阻 = 1 / G_total
7. 交叉更新边界：
     bc1: T_ext = T2_node, R_ext = R_gap_total + R2_int
     bc2: T_ext = T1_node, R_ext = R_gap_total + R1_int
```

**关键成员变量：**

- `self.gap` — 间隙宽度 \[m]
- `self.k_gas` — 气体导热系数 \[W/(m·K)]
- `self.R_gap_total` — 当前步计算的总间隙热阻 \[K/W] (向量)
- `self.eps1 / self.eps2` — 两侧发射率
- `self.sigma` — Stefan-Boltzmann 常数 5.670374419×10⁻⁸

**物理本质：** 辐射与气体导热并联 → 等效热阻串联在两侧固体之间。每个 sync() 调用时根据当前温度重新计算辐射热阻 (非线性 T⁴ 项在每个时间步重新线性化)。

***

### 4.4 ActiveGapCouple2D (带源项的间隙耦合器)

**继承关系：** `ActiveGapCouple2D` → `GapCouple2D` → `SolidSolidCouple2D`

**物理功能：** 在 GapCouple2D 基础上，额外增加**随时间变化的固定热流源** (如裂变气体电子冷却/加热)。利用**戴维南等效 (Thevenin Equivalent)** 将并联热流源转化为串联温差源，避免温差为 0 时的数值奇异点。

**新增接口：**

```python
set_active_heat_source(Q_array: Union[float, np.ndarray])
```

**参数约定：**

| 参数            | 含义                                     |
| ------------- | -------------------------------------- |
| `Q_array > 0` | 热量从 Obj1 流出，流入 Obj2 (Obj1 冷却, Obj2 加热) |
| `Q_array < 0` | 热量从 Obj2 流出，流入 Obj1                    |

**`sync()`** **方法流程：**

```
1. 调用父类 GapCouple2D.sync() → 计算基础热阻 (辐射∥气体导热)
2. 判断 Q_source 是否为零，若全零则跳过
3. 戴维南温差修正:
     ΔT = Q_source × R_gap_total
4. 应用修正:
     bc1.T_ext_new = T2_node - ΔT   (Obj1 侧：对方看起来更冷 → 多流出热)
     bc2.T_ext_new = T1_node + ΔT   (Obj2 侧：对方看起来更热 → 多接收热)
5. 保持热阻不变 (R_ext 沿用父类计算结果)
```

**戴维南等效原理：**

```
原始问题: 间隙中并联一个固定热流源 Q_source
等效变换: Q_source → 串联温差 ΔT = Q_source × R_gap
结果: 间隙两端等效温差被修正，而热阻网络结构不变
```

**关键成员变量：**

- `self.shape` — 边界节点形状 (如 (30,))
- `self._current_Q_source` — 当前热流源 \[W] (向量)

***

### 4.5 TECCouple2D (热电转换耦合器)

**继承关系：** `TECCouple2D` → `GapCouple2D` → `SolidSolidCouple2D`

**物理功能：** 专为**热离子燃料元件 (TFE) 电极间隙**设计。继承 GapCouple2D 的辐射+气体导热，同时允许在间隙**两侧分别施加独立且非守恒的热流源**：

- **发射极侧 (Side 1 / Obj1):** 电子冷却 → 热流流出 (负值)
- **接收极侧 (Side 2 / Obj2):** 电子加热 → 热流流入 (正值)

**核心接口：**

```python
set_tec_sources(Q_emitter: np.ndarray, Q_collector: np.ndarray)
```

**参数约定：**

| 参数            | 类型      | 含义               | 典型符号        |
| ------------- | ------- | ---------------- | ----------- |
| `Q_emitter`   | ndarray | 作用于发射极表面的热流 \[W] | 负值 (电子带走能量) |
| `Q_collector` | ndarray | 作用于接收极表面的热流 \[W] | 正值 (电子带来能量) |

**`sync()`** **方法流程：**

```
1. 调用父类 GapCouple2D.sync() → 计算间隙热阻网络和基础边界条件
2. 判断两侧 Q_source 是否全零，若全零则跳过
3. 获取当前的内部/外部热阻:
     R_ext_1, R_int_1 = Obj1 侧的外部/内部热阻
     R_ext_2, R_int_2 = Obj2 侧的外部/内部热阻
4. 获取当前的节点温度 (父类 sync 已设置 T_ext = T_neighbor_node):
     T_node_1 → bc2.T_ext
     T_node_2 → bc1.T_ext
5. 非对称戴维南修正:
     T_ext_1_new = T_node_2 + (Q_emitter × R_ext_1) + (Q_collector × R_int_2)
     T_ext_2_new = T_node_1 + (Q_collector × R_ext_2) + (Q_emitter × R_int_1)
6. 更新边界条件 (热阻保持不变)
```

**与 ActiveGapCouple2D 的区别：**

| 特性    | ActiveGapCouple2D | TECCouple2D                      |
| ----- | ----------------- | -------------------------------- |
| 热流源数量 | 1 个 (对称)          | 2 个 (非对称独立)                      |
| 能量守恒  | 守恒 (Q₁ = -Q₂)     | 非守恒 (Q\_emitter ≠ -Q\_collector) |
| 应用场景  | 裂变气体等效热源          | TFE 电极间隙 (电子冷却/加热不相等)            |
| 物理含义  | 间隙内部有热源           | 两侧表面分别有独立热流注入                    |

**关键成员变量：**

- `self.shape` — 边界节点形状
- `self.Q_source_1` — 发射极热流源 \[W] (向量)
- `self.Q_source_2` — 接收极热流源 \[W] (向量)

**典型调用链 (在 TECCircuitManager 中)：**

```python
# 从 C++ ThermoCalc 获取轴向分布结果
res = circuit.get_tec_results(idx)  # J[A/cm²], phiE[eV], phiC[eV]
# 计算电子冷却/加热热流 [W/m²] → [W]
q_emit = -f(J, phiE, TE) * area_per_node   # 负值: 电子冷却
q_coll = +f(J, phiC, TC) * area_per_node   # 正值: 电子加热
# 设置到耦合器
gap_coupler.set_tec_sources(q_emit, q_coll)
# sync() 在 SystemManager._run_couplers() 中被自动调用
```

**耦合器继承层次总结：**

```
SolidSolidCouple2D (基础固固耦合)
  ├── GapCouple2D (+ 间隙辐射 + 气体导热)
  │     ├── ActiveGapCouple2D (+ 单一热流源)
  │     └── TECCouple2D (+ 非对称双热流源)
```

***

## 5. PointReactor 点堆中子动力学

**物理模型：** 11 维 ODE 系统

- `y[0]` — 裂变功率 P\_fiss
- `y[1:7]` — 6 组缓发中子先驱核 C\_i
- `y[7:11]` — 4 组衰变热 W\_j

**方程：**

```
dP/dt = ((ρ - β) / Λ) * P + Σ(λ_i * C_i)
dC_i/dt = (β_i / Λ) * P - λ_i * C_i
dW_j/dt = γ_j * P - λ_dj * W_j
```

**关键参数：**

- `Λ = 0.2e-4` — 中子代时间
- `β_total = 0.0079321` — 总缓发中子份额
- `γ_total = 0.06421` — 总衰变热份额

**性能优化：** 使用 Numba JIT 编译 `_prke_rhs()` 和 `_prke_jac()`，解析雅可比矩阵加速 BDF/Radau 隐式求解。

**双缓冲状态管理：**

- `_y_committed` — 上一个收敛时间步的真实状态 (回滚起点)
- `_y_trial` — 当前 Picard 迭代中的试探状态

**关键方法：**

- `step(dt, reactivity_control, reactivity_feedback)` — 从 committed 重新积分
- `commit()` — 收敛后固化试探状态
- `initialize_steady_state(total_power)` — 解析计算绝对稳态

***

## 6. SystemManager 全局调度器

### 6.1 step() 方法完整流程

```
SystemManager.step(dt, inner_iter, convergence_tol, reactivity_control,
                   fail_on_fluid_nonconvergence, interface_relaxation,
                   interface_convergence_tol)

  ├── [保存入口状态] _save_system_state() → 流体+固体+组件+点堆
  │
  ├── [准备] 清除所有耦合器松弛状态
  │
  ├── [组件预处理] component.pre_step(dt, t_start)
  │     │  ← ReactorCore.pre_step: 更新功率分配、TEC 电热计算
  │     │  ← TECCircuitManager.pre_step: 求解全局电路、下发热流
  │     └  ← RingHP.pre_step: 更新翅片辐射热阻、外热流
  │
  ├── 保存预处理后的流体源项基准 (base_fluid_sources)
  │
  ├── [Picard 内迭代] for k in range(inner_iter):
  │     │
  │     ├── 恢复流体源项到基准状态
  │     │
  │     ├── _run_couplers(interface_relaxation)
  │     │     ├── _refresh_solid_boundary_cache()  → 刷新固体物性+热阻
  │     │     ├── for coupler in couplers: sync()  → 固固间隙耦合器双向同步
  │     │     ├── _refresh_solid_boundary_cache()  → 再次刷新
  │     │     └── for coupler in couplers: execute() → 流固耦合器计算 h,更新 BC,下发热源
  │     │
  │     ├── 中子学推进 (_advance_neutronics_for_iteration)
  │     │     ├── 组件自己处理 (如 ReactorCore.advance_neutronics)
  │     │     └── 否则 SystemManager 直接调用 point_reactor.step()
  │     │
  │     ├── 回滚+应用核功率
  │     │     ├── if k > 0: _rollback_system_state()  (只回滚流体+固体)
  │     │     └── _apply_pending_nuclear_power(fission, decay, total)
  │     │
  │     ├── fluid_solver.set_time(coupling_time)  → 向时间相关控制体和连接传播
  │     │
  │     ├── fluid_solver.step_Picard(dt)  → 水力网络压力+焓求解
  │     │
  │     ├── for solid in solid_components:
  │     │     └── solid.step(dt)  → 按 solid.ode_method 推进固体温度场，默认 BDF
  │     │
  │     ├── [收敛检查] if inner_iter > 1:
  │     │     ├── err_f = max|T_f_new - T_f_old|
  │     │     ├── err_s = max|T_s_new - T_s_old|
  │     │     ├── if err_f < tol and err_s < tol: break
  │     │     └── 可选: interface_residual 检查
  │     │
  │     └── T_f_prev, T_s_prev = T_f_curr, T_s_curr
  │
  ├── _commit_neutronics()  → 固化中子学试探状态
  │
  ├── global_time += dt
  │
  ├── _sync_solid_times_to_global()  → 同步所有固体时间
  │
  ├── component.post_step(dt, global_time)  → 组件后处理
  │
  └── [异常处理] 任何步骤失败 → _rollback_system_state() + 恢复流体源项
```

**自适应步长控制 (`compute_adaptive_dt`)：**

```
dt_target = min(dt_fluid, dt_coupler)
  ├── dt_fluid = fluid_solver.get_max_stable_dt() * safety_factor
  └── dt_coupler = min(coupler.get_max_stable_dt() across all couplers)

收敛控制:
  ├── 流体不收敛 → dt *= 0.5
  ├── Picard 迭代次数用尽 → dt *= 0.8
  └── dt_growth_limit → min(dt_target, 1.2 * dt_last)
```

***

## 7. 关键数据流图

```
SystemManager.step(dt)
     │
     ├── Component.pre_step()  ←── 功率分配、TEC求解
     │
     ├── Couplers:
     │     ├── SolidSolidCouple.sync()     ←→ BoundaryRegion ↔ BoundaryRegion
     │     └── FluidSolidCouple.execute()  ←→ FluidChannel  ↔ BoundaryRegion
     │                                            │
     │                                   ┌────────┴────────┐
     │                                   │ 计算 Nu → h → λ  │
     │                                   │ 更新固体 R_ext   │
     │                                   │ 更新流体源项     │
     │                                   └─────────────────┘
     │
     ├── PointReactor.step(dt, ρ)  ←── 反应性控制 + 温度反馈
     │     │                              (ReactorCore提供反馈ρ)
     │     └── 裂变功率 → 衰变热 → 总功率
     │
     ├── HydraulicNetwork.step_Picard(dt)
     │     ├── _update_fluid_properties()  ←── FluidMaterial (T, P → rho, mu, cp)
     │     ├── _calc_momentum_coeffs()     ←── 摩擦/形阻 → a_j, b_j
     │     ├── 组装并求解 A*P = B           ←── scipy.sparse.linalg.spsolve
     │     └── 半隐式焓方程求解            → P, T, h, W 更新
     │
     ├── HeatConduction.step(dt)   ←── 每个固体
     │     └── scipy.integrate.solve_ivp (solid.ode_method, default BDF)
     │           ├── get_derivatives(t, T):
     │           │     ├── _update_properties()    ←── SolidMaterial (T → k, rho, cp)
     │           │     ├── _compute_internal_resistance()  → R_geom 或 G
     │           │     ├── _update_boundaries_state(t)     → BoundaryRegion
     │           │     │     └── DynamicRadiationResistanceBC.update_state()  ← 动态T⁴辐射
     │           │     ├── _compute_fluxes(t)  → Q_net
     │           │     └── dTdt = (Q_net + Q_source) / C
     │           └── 返回 T(t+dt)
     │
     └── 收敛判断 → 循环或退出
```

***

## 8. 源码校验信息

最后按源码校验日期：2026-05-31。

涉及文件：

- `Hydrodynamics/Components.py`
- `Hydrodynamics/BoundaryVolume.py`
- `Hydrodynamics/HydraulicNetwork.py`
- `HeatConduction/Mesh.py`
- `HeatConduction/Boundary.py`
- `HeatConduction/HeatConduction.py`
- `Couplers.py`
- `Neutronics/PointReactor.py`
- `SystemManager.py`
