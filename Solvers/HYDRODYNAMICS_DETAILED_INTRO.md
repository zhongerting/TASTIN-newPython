# Hydrodynamics 水力网络模块详细介绍

## 1. 模块定位

`Hydrodynamics/` 是 TASTIN `Solvers` 层中的流体水力与热工输运求解子系统，主要负责冷却剂回路中的压力、质量流量、温度、焓、密度和粘度等状态量计算。

该模块面向一维流体网络建模，使用控制体和连接件构造任意拓扑的管路、联箱、泵、边界和代表通道。它与 `HeatConduction`、`Couplers.py`、`Neutronics` 和 `SystemManager` 共同组成多物理场瞬态求解框架。

目录结构如下：

```text
Hydrodynamics/
  Components.py        # 流体控制体、连接件、泵、通道等基础组件
  BoundaryVolume.py    # 压力/温度/流量边界组件
  HydraulicNetwork.py  # 半隐式水力网络求解器
```

模块的核心职责包括：

- 构建流体控制体和连接件的拓扑网络；
- 求解压力场和质量流量场；
- 计算摩擦压降、局部压降、泵压升、重力项和加速项；
- 支持可压缩和不可压缩流体建模方式；
- 求解焓/温度输运方程；
- 接收来自固体导热耦合器的换热源项；
- 提供自适应时间步长、状态保存、回滚和断点续算接口。

## 2. 建模思想

Hydrodynamics 使用两类基本对象描述流体网络：

| 对象 | 物理含义 | 典型状态 |
| --- | --- | --- |
| `FluidVolume` | 流体控制体，表示一个节点或一段管道单元 | `P`, `h`, `T`, `rho`, `mu` |
| `FlowJunction` | 两个控制体之间的流动连接 | `W`, `vel` |

控制体负责质量和能量守恒，连接件负责动量方程。复杂通道通过 `FluidChannel` 自动生成一串 `FluidVolume` 和 `FlowJunction`。

求解器 `HydraulicNetwork` 会把这些对象转译为向量和稀疏矩阵：

```text
volumes_obj[]   -> P_vec, T_vec, h_vec, rho_vec, mu_vec
junctions_obj[] -> W_vec, A_coeffs, B_coeffs
topology        -> CSR 压力矩阵和能量矩阵
```

因此，组件层保持面向对象和物理可读性，求解器层使用 NumPy/SciPy 向量化和稀疏线性代数提高效率。

## 3. `Components.py`：基础水力组件

### 3.1 `FluidVolume`

`FluidVolume` 是可压缩流体控制体，存储单个节点的流体状态。

核心属性：

| 属性 | 含义 |
| --- | --- |
| `name` | 控制体名称 |
| `vol` | 控制体体积，单位 m3 |
| `len` | 控制体长度，单位 m |
| `area` | 流通面积，单位 m2 |
| `d_h` | 水力直径，单位 m |
| `P` | 压力，单位 Pa |
| `h` | 比焓，单位 J/kg |
| `T` | 温度，单位 K |
| `rho` | 密度，单位 kg/m3 |
| `mu` | 动力粘度，单位 Pa.s |
| `Q_wall` | 壁面热源，单位 W |
| `Q_vol` | 体积热源，单位 W |
| `implicit_coeff` | 半隐式热源系数，单位 W/K |

拓扑连接通过两个列表维护：

```python
self.inlet_junctions
self.outlet_junctions
```

`FluidVolume.get_volume_derivatives()` 求解质量和能量守恒：

```text
质量:
  d(rho)/dt = (sum(W_in) - sum(W_out)) / V

能量:
  rho * dh/dt = (sum(W * h_donor) + Q_total) / V

源项:
  Q_total = Q_wall + Q_vol - implicit_coeff * T
```

焓输运采用施主元胞策略：

- 流入当前节点时，使用上游节点焓；
- 流出当前节点时，使用当前节点焓；
- 发生倒流时，自动切换供体节点。

压力导数由密度偏导数给出：

```text
drho/dP * dP/dt + drho/dh * dh/dt = mass_rhs
```

代码中使用人工可压缩性兜底，虚拟声速 `c_art = 200 m/s`，防止不可压缩液体导致压力方程过硬或除零。

### 3.2 `IncompressibleFluidVolume`

`IncompressibleFluidVolume` 继承自 `FluidVolume`，用于不可压缩液体的压力-流量代数解耦。

关键区别：

- 压力 `P` 不作为 ODE 状态演化；
- `get_volume_derivatives()` 只返回 `dh/dt`；
- 压力由外部网络通过 `update_pressure_algebraic(new_P)` 更新；
- 更新压力后立即刷新 `T`、`rho`、`mu`。

不可压缩能量方程采用非守恒输运形式：

```text
rho * V * dh/dt =
  sum(W_in * (h_donor - h_self)) + Q_total
```

这种写法避免质量不平衡误差直接污染能量方程，适合液态金属等低可压缩性冷却剂。

### 3.3 `PressurizerVolume`

`PressurizerVolume` 是不可压缩闭式回路中的被动压力参考点。

关键标记：

```python
is_pressure_reference = True
is_pressure_boundary = False
```

它不是质量源或固定压力边界，而是用于锚定闭式网络的绝对压力水平。`HydraulicNetwork` 求解完压力分布后，会通过压力参考目标做整体平移，使参考节点压力等于目标压力。

支持接口：

```python
set_pressure(pressure)
set_pressure_table(times, pressures)
clear_pressure_table()
compute_target_pressure(time)
```

时间表通过 `np.interp` 插值，超出范围时钳位到端点。

### 3.4 `FixedPressurePumpVolume`

`FixedPressurePumpVolume` 是带固定压升的不可压缩控制体，适用于代数压力递推通道。

核心接口：

```python
set_delta_p(delta_p)
get_pressure_rise()
```

它的压升不直接改变能量方程，而是在通道压力分布递推时作为局部压力增量使用。

### 3.5 `FlowJunction`

`FlowJunction` 连接两个控制体，存储质量流量和速度，并自动注册到两端控制体：

```python
from_vol.outlet_junctions.append(self)
to_vol.inlet_junctions.append(self)
```

核心属性：

| 属性 | 含义 |
| --- | --- |
| `from_vol` | 正方向上游控制体 |
| `to_vol` | 正方向下游控制体 |
| `area` | 连接流通面积 |
| `length` | 惯性长度 |
| `k_loss` | 局部阻力系数 |
| `W` | 质量流量，正方向为 `from_vol -> to_vol` |
| `vel` | 按施主密度计算的流速 |
| `dynamic_loss_params` | 动态阻力模型参数 |

动量方程形式为：

```text
(L / A) * dW/dt =
  (P_from - P_to) + gravity + pump - friction - form_loss
```

摩擦压降：

```text
dP_fric = f * (L / D_h) * 0.5 * rho * v^2
```

其中 `f` 来自 `Correlations.Correlations.friction_single_phase(Re)`。

局部压降：

```text
dP_form = K * 0.5 * rho * v^2
```

当 `dynamic_loss_params` 启用且模型不是 `none` 时，会叠加动态局部阻力。

### 3.6 动态局部阻力模型

当前内置的动态模型包括 `single_crossflow_pipe` 和 `inline_tube_bank_euler`。动态阻力由 `dynamic_loss_params["model"]` 启用，并在水力 Picard 迭代中按当前流量刷新。

#### 3.6.1 `single_crossflow_pipe`

`single_crossflow_pipe` 用于描述单根横掠管对流道的阻塞和绕流阻力。

主要计算过程：

```text
A_proj = D_out * L_pipe
A_min = A_flow - A_proj
v_max = v_nominal * A_flow / A_min
Re_D = rho * v_max * D_out / mu
```

根据 `Re_D` 分段计算圆管绕流阻力系数 `C_D`：

| Re 范围 | C_D |
| --- | --- |
| `Re < 1` | `24 / Re` |
| `1 <= Re < 1000` | `24/Re + 3/sqrt(Re) + 0.34` |
| `1000 <= Re < 2e5` | `1.2` |
| `Re >= 2e5` | `0.3` |

再换算为网络连接件可用的等效阻力：

```text
K_eq = C_D * (A_proj / A_flow) * (A_flow / A_min)^2
```

若投影面积接近堵塞流道，返回极大的阻力系数 `1e5`。

#### 3.6.2 `inline_tube_bank_euler`

`inline_tube_bank_euler` 用于顺排横掠管束压降。当前 `RingHP` 用它描述集流环内热管蒸发段造成的附加局部阻力，公式为：

```text
A_proj = D_out * L_pipe
A_min = A_flow - A_proj
S_T = A_flow / L_pipe
v_max = |W| / (rho * A_min)
Re_D = rho * v_max * D_out / mu
Eu = 0.67 * (S_T / D_out - 1)^(-0.5) * Re_D^(-0.15)
K_eq = Eu * N_rows * (A_flow / A_min)^2
Delta p = K_eq * 0.5 * rho * v_nominal^2
```

`dynamic_loss_params` 可显式传入 `pitch_ratio = S_T/D_out`；否则按 `A_flow / L_pipe / D_out` 计算。`N_rows` 表示该连接件所在控制体内沿流向串联的管排数。若 `A_min <= 0`、投影面积接近堵塞流道或 `S_T/D_out <= 1`，模型返回极大的阻力系数 `1e5`。

### 3.7 `PumpJunction`

`PumpJunction` 继承自 `FlowJunction`，在动量方程中额外提供泵压升。

关键标记：

```python
is_pump_junction = True
```

接口：

```python
set_delta_p(delta_p)
set_pressure_table(times, pressures)
compute_pump_head(time)
```

正的 `delta_p` 驱动 `from_vol -> to_vol` 方向流动。若设置时间表，泵压升随 `current_time` 插值变化。

### 3.8 `MacroFlowJunction`

`MacroFlowJunction` 用于“代表管”或“单通道代表多根并联管”的跨尺度建模。

设计原则：

- 动量方程仍按单根管计算；
- 质量和能量守恒中，宏观端看到 `multiplier * W`；
- 微观端仍看到单管流量 `W`。

核心接口：

```python
get_mass_flow_for(target_vol)
```

如果 `target_vol` 是宏观端，则返回放大后的流量；否则返回真实单管流量。`HydraulicNetwork` 会把这个乘子缓存到 `M_from_vec` 和 `M_to_vec`，用于压力方程和能量方程组装。

### 3.9 `FluidChannel`

`FluidChannel` 是复合组件，用于自动创建均匀一维通道。

构造时完成：

```text
1. 根据 total_length / n_nodes 得到 node_length
2. 创建 n_nodes 个 FluidVolume
3. 创建 n_nodes - 1 个内部 FlowJunction
4. 为每个节点写入 mesh_index 和 z_coordinate
```

常用向量化访问接口：

```python
temperature_vector
pressure_vector
density_vector
velocity_vector
```

耦合源项接口：

```python
add_coupling_source_distribution(explicit_arr, implicit_arr)
add_heat_source_distribution(heat_dist)
clear_sources()
```

其中流固耦合器通常写入：

```text
Q_fluid = explicit_part - implicit_factor * T_fluid
explicit_part = h * A * T_wall
implicit_factor = h * A
```

### 3.10 `IncompressibleFluidChannel`

`IncompressibleFluidChannel` 自动创建 `IncompressibleFluidVolume` 节点，并提供代数压力递推方法：

```python
update_pressure_distribution_downstream(P_inlet)
update_pressure_distribution_upstream(P_outlet)
```

递推中考虑：

- 摩擦压降；
- 局部压降；
- 可选压力泵体压升；
- 当前实现中空间工况下重力项设为 `g = 0`。

### 3.11 `NonUniformIncompressibleFluidChannel`

`NonUniformIncompressibleFluidChannel` 支持非均匀轴向网格。它接收 `node_lengths` 数组而不是总长和节点数。

特点：

- `node_length` 是数组；
- `node_volume = area * node_length`；
- 节点中心 `z_coordinate` 由累积长度计算；
- 内部连接惯性长度为相邻两个节点半长之和；
- 压力递推方法继承自 `IncompressibleFluidChannel`。

该类适合轴向局部加密通道，例如堆芯功率峰值区、换热入口区或热电耦合强梯度区。

## 4. `BoundaryVolume.py`：边界组件

### 4.1 `BoundaryVolume`

`BoundaryVolume` 继承自 `FluidVolume`，用于可压缩求解中的压力/温度边界。

它不直接强制覆盖状态，而是通过松弛导数驱动当前状态追随目标：

```text
dP/dt = K * (P_target - P)
dh/dt = K * (h_target - h)
```

当前松弛增益：

```text
K = 50 1/s
```

接口：

```python
set_state(P=None, T=None)
```

该接口更新目标压力和目标温度，实际状态由积分过程逐步逼近。

### 4.2 `IncompressibleBoundaryVolume`

`IncompressibleBoundaryVolume` 用于不可压缩回路的边界或集管。

它支持两种热工模式：

| 模式 | `mixing_enabled` | 行为 |
| --- | --- | --- |
| 无限大热池 | `False` | 忽略流入流出对自身温度影响，`dh/dt = forcing_dydt_h` |
| 有限体积混合 | `True` | 调用父类能量方程，根据流入流出焓差计算混合温度 |

关键接口：

```python
set_boundary_state(P=None, T=None)
update_pressure_algebraic(new_P)
```

`forcing_dydt_h` 可由外部设置，用于实现边界温度随时间变化。

### 4.3 `InletJunction`

`InletJunction` 是质量流量边界连接，继承自 `FlowJunction`，但不求解普通动量方程。

它通过松弛法驱动流量追随目标：

```text
dW/dt = K * (W_target - W)
```

当前源码中组件级 `gain_k = 0.1`，而 `HydraulicNetwork` 的向量化快速路径对入口连接使用 `inlet_relaxation_gain = 50.0` 生成离散更新：

```text
W_new = W_old + dt * gain * (W_target - W_old)
```

实际网络求解中主要由 `HydraulicNetwork` 的快速路径处理入口流量。

## 5. `HydraulicNetwork.py`：半隐式水力网络求解器

`HydraulicNetwork` 是 Hydrodynamics 的核心求解器。它接收 `volumes` 和 `junctions` 列表，在初始化阶段构建固定拓扑缓存，时间推进时只更新矩阵数据和右端项。

### 5.1 核心状态向量

| 向量 | 长度 | 含义 |
| --- | --- | --- |
| `P_vec` | `n_vol` | 节点压力 |
| `T_vec` | `n_vol` | 节点温度 |
| `h_vec` | `n_vol` | 节点比焓 |
| `rho_vec` | `n_vol` | 节点密度 |
| `mu_vec` | `n_vol` | 节点粘度 |
| `W_vec` | `n_junc` | 连接件质量流量 |
| `A_coeffs` | `n_junc` | 动量线性化导纳 |
| `B_coeffs` | `n_junc` | 动量线性化源项 |

几何缓存：

```text
V_vec, L_node_vec, A_node_vec, Dh_node_vec, z_vec
A_junc_vec, L_junc_vec, K_loss_vec
idx_from_vec, idx_to_vec
M_from_vec, M_to_vec
```

矩阵缓存：

```text
pressure_matrix     # CSR 压力方程结构
energy_matrix       # CSR 焓方程结构
*_ptrs              # CSR data 数组中的固定写入位置
```

### 5.2 拓扑构建

`_build_topology()` 完成三件事：

1. 为每个 `FluidVolume` 建立节点索引；
2. 提取节点和连接件几何常数；
3. 识别定压边界、被动压力参考点、入口连接、泵连接和动态阻力连接。

压力相关节点分为：

| 类型 | 标记 | 作用 |
| --- | --- | --- |
| 定压边界 | `is_pressure_boundary=True` | 压力矩阵中对应行为 Dirichlet 条件 |
| 被动压力参考 | `is_pressure_reference=True` | 闭式回路求解后整体平移压力 |
| 普通节点 | 默认 | 参与质量守恒压力方程 |

源码明确禁止同一节点同时作为定压边界和压力参考点。

初始化完成后，`fixed_pressure_indices` 被锁定，不能再修改。这是因为压力矩阵的 CSR 非零结构已经按固定边界集合预构建。

### 5.3 压力矩阵缓存

`_build_pressure_system_cache()` 在固定拓扑下预构建压力矩阵结构。

普通节点的质量守恒离散形式可理解为：

```text
(V / dt) * drho/dP * P_new + sum(W_new) = RHS
```

连接件流量通过动量线性化表示：

```text
W_j = a_j * (P_from - P_to) + b_j
```

代入质量守恒后，每个连接在压力矩阵中贡献对角项和非对角项：

```text
from 行: +M_from*a_j*P_from - M_from*a_j*P_to
to 行:   +M_to*a_j*P_to   - M_to*a_j*P_from
```

`M_from` 和 `M_to` 来自 `MacroFlowJunction`，普通连接默认为 1。

定压边界行被设置为：

```text
P_i = P_target
```

### 5.4 能量矩阵缓存

`_build_energy_system_cache()` 预构建全隐式焓方程矩阵结构。矩阵包含：

- 每个节点自身的惯性项；
- 根据流向激活的迎风输运对角/非对角项；
- 宏观乘子对质量流和能量流的缩放。

实际每个时间步只更新 CSR `data` 和右端项。

### 5.5 物性更新

`_update_fluid_properties()` 使用主物性对象对整条网络向量化计算：

```python
cp_vec  = material.heat_capacity(T_vec, P_vec)
rho_vec = material.density(T_vec, P_vec)
mu_vec  = material.viscosity(T_vec, P_vec)
```

若物性对象提供密度偏导，则同步更新：

```python
drho_dp_vec = material.liquid_density_derivative_P(P_vec)
drho_dt_vec = material.liquid_density_derivative_T(T_vec)
```

数值保护：

- `rho >= 1e-1`；
- `mu >= 1e-10`；
- `drho_dp >= 1e-11`。

同时从底层 `Volume` 对象读取热源：

```text
Q_expl_vec = Q_wall + Q_vol
lam_imp_vec = implicit_coeff
```

### 5.6 摩擦因子与动态阻力

`HydraulicNetwork` 内部提供向量化摩擦因子计算：

```text
Re <= 2300:          f = 64 / Re
2300 < Re < 1e5:    f = 0.3164 / Re^0.25
Re >= 1e5:          Karman-Prandtl 迭代
```

动态阻力通过 `_refresh_effective_k_loss()` 刷新到 `effective_K_loss_vec`。对于 `single_crossflow_pipe`，求解器实现了与 `FlowJunction._compute_dynamic_k_loss()` 等价的向量路径。

### 5.7 动量方程线性化

`_calc_momentum_coeffs_fast()` 把每个连接的动量方程离散为：

```text
W_new = a_j * (P_from - P_to) + b_j
```

其中：

```text
I_term = L_inertial / (A_flow * dt)
K_linear = (f*L/D + K_loss) * |W_iter| / (2*rho*A_flow^2)

a_j = 1 / (I_term + K_linear)
b_j = a_j * (I_term*W_old + dP_grav + dP_acc + dP_pump)
```

代码还考虑：

- 迎风密度和粘度的平滑切换；
- 上下游半节点摩擦长度；
- 重力压差 `rho * g * dz`；
- 动压加速项；
- 泵连接的 `compute_pump_head(current_time)`；
- 入口连接的目标流量松弛。

### 5.8 压力求解与流量回代

`_assemble_pressure_system()` 更新压力矩阵和右端项：

```text
1. 清零 CSR data
2. 写入连接件贡献 a_j
3. 写入普通节点压缩性项
4. 写入定压边界行
5. 写入热膨胀源项和 b_j 源项
```

`_solve_linear_system()` 使用：

```python
scipy.sparse.linalg.spsolve(A_sparse, B)
```

得到 `P_vec` 后，流量回代：

```text
W_vec = A_coeffs * (P_from - P_to) + B_coeffs
```

若存在被动压力参考点，求解器会在压力求解后做整体平移，使参考节点压力等于目标压力。

### 5.9 热膨胀源项

压力方程可包含热膨胀源项：

```text
S_thermal = -V * d(rho)/dt
          = -V * (d rho / d h) * dh/dt
```

其中：

```text
d rho / d h = (d rho / d T) / Cp
```

`dh/dt` 由当前流量和热源显式预测得到。这个源项让加热导致密度降低、体积膨胀并驱动压力/流量变化。

### 5.10 全隐式焓方程

当前 `step_Picard()` 使用 `_step_energy_implicit()` 推进能量方程。

核心形式：

```text
M * (h_new - h_old) / dt
  + upwind_advection(h_new)
  + lambda_h * h_new
  = M/dt * h_old + Q_modified
```

其中：

```text
M = rho * V
lambda_h = implicit_coeff / Cp
Q_modified = Q_expl - lambda*T_old + lambda_h*h_old
```

这种做法使用焓作为未知量，并通过稀疏矩阵求解：

```python
h_new_vec = spsolve(energy_matrix, B_enth)
```

随后优先调用物性对象的向量化接口反算温度：

```python
T_vec = material.temperature_from_enthalpy(h_new_vec, P_vec)
```

如果反算失败，则退化为 `Cp` 近似更新温度。

### 5.11 `step_Picard()` 主流程

`step_Picard(dt, max_iter, tol)` 是当前全局调度中最重要的接口。

流程：

```text
1. 更新物性
2. 计算并冻结热膨胀源项
3. 保存时间步入口流量 W_old
4. Picard 迭代:
   4.1 备份当前 W
   4.2 用当前 W 计算摩擦因子、动态阻力和动量系数
   4.3 组装压力矩阵
   4.4 求解压力
   4.5 回代更新流量
   4.6 检查 max(abs(W_new - W_old_iter)) < tol
5. 流场收敛或达到最大迭代次数后，推进全隐式焓方程
6. 将压力、流量、温度、焓同步回对象层
```

即使未在最大迭代次数内收敛，函数也会推进能量方程并返回 `False`，由上层 `SystemManager` 决定是否回滚、缩短步长或继续。

### 5.12 液力初始化

`initialize_hydraulics()` 用于冻结温度场，只迭代压力和流量，以消除初始流场非物理震荡。

特点：

- 不推进能量方程；
- 使用虚拟时间步作为松弛尺度；
- 每轮压力求解后对流量做亚松弛：

```text
W = omega * W_new + (1 - omega) * W_old
```

适合在瞬态计算开始前建立与当前阻力和边界匹配的初始水力状态。

### 5.13 自适应步长接口

`get_max_stable_dt(max_limit)` 基于 CFL 条件估算流体输运时间步：

```text
dt < L_node / v_max
```

对每个非虚拟节点，取其相邻连接中的最大流速，返回全网最小限制。该值通常由 `SystemManager.compute_adaptive_dt()` 进一步乘安全因子并与耦合器稳定步长共同取最小。

### 5.14 状态同步与持久化

求解器内部以向量为主，但对象层仍需要与之同步：

```python
_sync_vectors_to_objects(
    sync_pressure=True,
    sync_flow=True,
    sync_energy=True,
    sync_properties=True,
)
```

状态保存接口：

```python
save_state()
load_state()
```

用于运行时回滚。

断点续算接口：

```python
get_state_dict(prefix)
load_state_dict(data, prefix)
```

持久化内容包括：

- 拓扑指纹 `[n_vol, n_junc]`；
- 定压边界索引集合；
- `P_vec`、`T_vec`、`h_vec`、`W_vec`；
- 边界目标压力；
- 入口目标流量。

恢复时会检查拓扑和固定压力边界集合是否匹配，因为 CSR 矩阵结构依赖这些信息。

### 5.15 网络时间传播

`HydraulicNetwork.set_time(time)` 保存 `current_time`，并向所有提供 `set_time()` 的控制体和连接传播时间：

```python
network.set_time(time)
```

当前直接受益的对象包括：

- `PressurizerVolume`：刷新时间表插值得到的绝对压力参考目标；
- `PumpJunction`：刷新时间表插值使用的本地时间。

网络快速路径计算泵压升时还会显式调用：

```python
junction.compute_pump_head(self.current_time)
```

`SystemManager.step()` 在流体求解前同步耦合时刻。首轮全局 Picard 使用步首时间，后续轮次使用步末时间。

## 6. 与 Couplers 和 SystemManager 的关系

### 6.1 流固耦合

`FluidSolidCouple` 通过 `FluidChannel` 的向量接口读取流体状态：

```python
temperature_vector
pressure_vector
density_vector
velocity_vector
```

然后计算：

```text
Re = rho * |v| * D_h / mu
Pr = material.prandtl_number(T, P)
Nu = correlation_func(Re, Pr, ...)
h = Nu * k_fluid / D_h
lambda = h * A
```

写回流体源项：

```text
Q = lambda * T_wall - lambda * T_fluid
```

对应到 `FluidVolume`：

```text
Q_wall += lambda * T_wall
implicit_coeff += lambda
```

这样能量方程中出现：

```text
Q_wall - implicit_coeff * T_fluid
```

属于半隐式换热处理。

### 6.2 与全局 Picard 迭代

在 `SystemManager.step()` 中，Hydrodynamics 通常按如下顺序参与：

```text
1. 组件 pre_step 更新外部热源、泵、TEC 等
2. 运行耦合器，向流体写入换热源项
3. 推进中子学并分配功率
4. 调用 fluid_solver.step_Picard(dt)
5. 推进所有固体导热对象
6. 检查流体和固体温度收敛
```

因此，Hydrodynamics 的热源通常不是静态输入，而是在每个 Picard 内迭代中由耦合器重新写入。

## 7. 典型使用流程

### 7.1 创建流体通道

```python
from Solvers.Hydrodynamics.Components import FluidChannel

channel = FluidChannel(
    name="coolant_channel",
    n_nodes=20,
    total_length=1.0,
    flow_area=1.0e-4,
    hydraulic_diam=0.01,
    initial_P=1.0e5,
    initial_T=700.0,
    material=fluid_material,
)
```

### 7.2 创建边界和入口流量

```python
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction

inlet = IncompressibleBoundaryVolume(
    name="inlet_plenum",
    material=fluid_material,
    P=1.0e5,
    T=650.0,
)

inlet_junc = InletJunction(
    name="inlet_flow",
    from_vol=inlet,
    to_vol=channel.volumes[0],
    W_initial=0.1,
)
inlet_junc.set_flow_rate(0.1)
```

### 7.3 创建网络并推进

```python
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork

volumes = [inlet] + channel.volumes
junctions = [inlet_junc] + channel.internal_junctions

network = HydraulicNetwork(volumes, junctions)
network.initialize_hydraulics(dt=0.1, tol=1e-5)

ok = network.step_Picard(dt=0.01, max_iter=20, tol=1e-4)
```

### 7.4 添加流固耦合源项

```python
explicit = hA * wall_temperature
implicit = hA
channel.add_coupling_source_distribution(explicit, implicit)
```

下一次 `network.step_Picard()` 中，能量方程将使用：

```text
Q = explicit - implicit * T_fluid
```

## 8. 重要实现细节与注意事项

### 8.1 拓扑在初始化后固定

`HydraulicNetwork` 初始化时会构建压力矩阵和能量矩阵的 CSR 结构。初始化后不能再改变固定压力边界集合，也不应随意增删节点和连接。

如果需要改变拓扑，应重新创建 `HydraulicNetwork`。

### 8.2 热源单位

`Q_wall`、`Q_vol`、`explicit_arr` 都是节点总热功率，单位 W。若外部给定热流密度 W/m2，需要先乘以对应节点换热面积。

### 8.3 流量符号

`FlowJunction.W > 0` 表示从 `from_vol` 流向 `to_vol`。所有施主元胞、阻力方向、泵压升和宏观乘子都基于这个方向定义。

### 8.4 定压边界与压力参考点不同

- 定压边界会把压力矩阵对应行变成 Dirichlet 条件；
- 被动压力参考点只用于闭式回路压力整体平移；
- 同一节点不能同时是二者。

### 8.5 入口连接的处理路径

`InletJunction` 在组件层提供 `get_momentum_derivative()`，但在当前 `HydraulicNetwork` 快速路径中，入口连接通过 `target_W_vec` 和 `inlet_relaxation_gain` 直接更新 `B_coeffs`。调试入口流量响应时应优先查看网络层参数。

### 8.6 物性对象应支持向量化

`HydraulicNetwork` 会向材料对象传入整个数组：

```python
material.heat_capacity(T_vec, P_vec)
material.density(T_vec, P_vec)
material.viscosity(T_vec, P_vec)
```

为了性能和正确性，流体材料类应支持 NumPy 数组输入。

### 8.7 `step()` 与 `step_Picard()`

源码中保留了 `step()`，但当前全局调度和文档分析都以 `step_Picard()` 为主。`step_Picard()` 会对非线性阻力进行外部迭代，并使用全隐式焓方程推进温度，更适合强耦合瞬态。

### 8.8 动态阻力的启用条件

动态阻力由 `dynamic_loss_params` 中非 `none` 的 `model` 启用，不由负的 `k_loss` 单独启用。启用动态模型后，若基础 `k_loss < 0`，组件层和网络快速路径都会先把基础值按 `0` 处理，再叠加动态阻力。

## 9. 模块优点

Hydrodynamics 当前实现具有以下特点：

- 控制体/连接件建模直观，适合搭建复杂回路；
- 求解器层将对象网络转为向量和稀疏矩阵，性能更好；
- 压力矩阵和能量矩阵结构预缓存，避免每步重复构建拓扑；
- 支持入口流量、泵压升、定压边界、压力参考点和宏观乘子；
- 支持动态局部阻力，适合热管阵列、横掠障碍等局部结构；
- 全隐式焓方程打破显式对流 CFL 对能量方程的严格限制；
- 能与固体导热边界通过半隐式换热源项稳定耦合；
- 提供初始化、回滚、断点续算和自适应步长接口。

## 10. 推荐阅读顺序

维护或扩展 Hydrodynamics 时，建议按以下顺序阅读：

1. `Hydrodynamics/Components.py`
   - 先理解 `FluidVolume`、`FlowJunction`、`FluidChannel` 的物理职责；
2. `Hydrodynamics/BoundaryVolume.py`
   - 理解压力、温度和流量边界如何表示；
3. `Hydrodynamics/HydraulicNetwork.py`
   - 重点阅读 `_build_topology()`、`_calc_momentum_coeffs_fast()`、`_assemble_pressure_system()`、`_step_energy_implicit()` 和 `step_Picard()`；
4. `Couplers.py`
   - 查看流固耦合器如何向 `FluidVolume` 写入半隐式热源；
5. `SystemManager.py`
   - 理解水力网络在全局 Picard 内迭代中的执行顺序。

## 11. 源码校验信息

最后按源码校验日期：2026-05-31。

涉及文件：

- `Hydrodynamics/Components.py`
- `Hydrodynamics/BoundaryVolume.py`
- `Hydrodynamics/HydraulicNetwork.py`
- `Couplers.py`
- `SystemManager.py`
