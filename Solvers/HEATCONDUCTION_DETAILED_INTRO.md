# HeatConduction 固体导热模块详细介绍

## 1. 模块定位

`HeatConduction/` 是 TASTIN `Solvers` 层中的固体瞬态导热求解子系统，主要用于求解燃料、包壳、电极、管壁、热管翅片等固体区域的温度场。它与水力网络、中子动力学以及多物理场耦合器共同工作，在 `SystemManager` 的 Picard 内迭代中被逐步推进。

该模块的核心目标是：

- 用有限体积法离散固体导热方程；
- 支持 1D 与 2D 固体域；
- 支持笛卡尔平板和圆柱/轴对称几何；
- 用统一的边界热阻网络表达对流、定热流、辐射、固固接触和耦合边界；
- 通过 `scipy.integrate.solve_ivp` 的隐式积分方法推进刚性瞬态热问题；
- 为 `Couplers.py` 提供稳定、向量化的边界状态接口。

目录结构如下：

```text
HeatConduction/
  Mesh.py             # 1D/2D 网格与几何量生成
  Boundary.py         # 热边界条件与边界区域管理
  HeatConduction.py   # 1D/2D 固体导热求解器
```

## 2. 控制方程与数值形式

固体导热模块求解的基本能量方程可以写成：

```text
rho * cp * V * dT/dt = Q_conduction + Q_source
```

其中：

- `T` 是控制体中心温度；
- `rho`、`cp`、`k` 来自 `SolidMaterial`，可随温度变化；
- `V` 是控制体体积；
- `Q_conduction` 是内部界面和边界贡献的净热流；
- `Q_source` 是体内热源，单位为 W，已按控制体积分配。

求解器内部统一将温度场扁平化为一维数组：

```python
self.T = np.full(self.N, initial_temp)
self.dTdt = np.zeros(self.N)
self.Q_source = np.zeros(self.N)
```

即使是 2D 网格，ODE 求解器接收和返回的仍然是一维状态向量。二维结构通过 `reshape(self.shape_nodes)` 得到视图。

每次 RHS 评估的流程是：

```text
get_derivatives(t, T_current)
  1. 接收试探温度 T_current
  2. 更新物性 k, rho, cp
  3. 更新控制体热容 C = rho * cp * V
  4. 计算内部热阻或热导
  5. 将边界节点温度和内部热阻推送给 BoundaryRegion
  6. 计算内部导热与边界热流形成的 Q_net
  7. 更新内热源 Q_source
  8. 返回 dTdt = (Q_net + Q_source) / C
```

## 3. `Mesh.py`：网格与几何数据

### 3.1 `Mesh1D`

`Mesh1D` 负责生成一维有限体积网格，支持两类几何：

| 几何类型 | 说明 | 典型用途 |
| --- | --- | --- |
| `cartesian` | 平板/直角坐标一维导热 | 平板、简化壁面 |
| `cylindrical` | 圆柱径向导热 | 燃料棒、管壁、圆柱电极 |

构造方式有两种：

```python
Mesh1D(
    total_dim=0.01,
    n_volumes=20,
    geometry_type="cylindrical",
    inner_radius=0.0,
    height=0.5,
)
```

或使用自定义界面坐标创建非均匀网格：

```python
Mesh1D.from_custom_faces(
    face_locations=np.array([0.0, 0.001, 0.003, 0.006, 0.01]),
    geometry_type="cartesian",
    height=1.0,
)
```

核心几何属性：

| 属性 | 形状 | 含义 |
| --- | --- | --- |
| `face_locations` | `N+1` | 控制体界面位置 |
| `node_centers` | `N` | 控制体中心位置 |
| `volumes` | `N` | 控制体体积 |
| `face_areas` | `N+1` | 每个界面的传热面积 |
| `dr_node_to_node` | `N-1` | 相邻控制体中心距离 |
| `dr_node_to_face` | `N, 2` | 节点中心到左右界面的距离 |

笛卡尔几何中，体积和面积按单位宽度/深度处理：

```text
V_i = dx_i * H
A_face = H
```

圆柱几何中，体积和界面面积为：

```text
V_i = pi * (r_outer^2 - r_inner^2) * H
A_face = 2 * pi * r_face * H
```

### 3.2 `Mesh2D`

`Mesh2D` 用于二维有限体积网格，支持：

- `cartesian`：二维平面导热；
- `cylindrical`：轴对称 `r-z` 导热。

构造函数支持均匀网格，也支持 `x_faces` 和 `y_faces` 非均匀网格：

```python
Mesh2D(
    x_dim=0.01,
    n_x=20,
    y_dim=1.0,
    n_y=50,
    geometry_type="cylindrical",
    inner_radius=0.0,
)
```

二维网格的关键设计是：核心数据扁平化存储，但保留矩阵视图的访问接口。

| 属性 | 形状 | 含义 |
| --- | --- | --- |
| `shape_nodes` | `(n_x, n_y)` | 节点矩阵形状 |
| `shape_faces_x` | `(n_x+1, n_y)` | X/R 方向界面形状 |
| `shape_faces_y` | `(n_x, n_y+1)` | Y/Z 方向界面形状 |
| `geom_data.volumes` | `N` | 扁平化控制体体积 |
| `geom_data.area_x` | `(n_x+1)*n_y` | 垂直于 X/R 方向的界面面积 |
| `geom_data.area_y` | `n_x*(n_y+1)` | 垂直于 Y/Z 方向的界面面积 |
| `geom_data.dx_node_to_node` | `(n_x+1)*n_y` | 跨 X/R 界面的节点距离 |
| `geom_data.dy_node_to_node` | `n_x*(n_y+1)` | 跨 Y/Z 界面的节点距离 |

扁平化索引规则是：

```text
k = i * n_y + j
```

其中 `i` 是 X/R 方向索引，`j` 是 Y/Z 方向索引。

辅助接口：

```python
mesh.flatten_index(i, j)
mesh.unravel_index(k)
mesh.volumes_matrix
mesh.area_x_matrix
mesh.area_y_matrix
mesh.dx_matrix
mesh.dy_matrix
```

## 4. `Boundary.py`：边界条件与热阻网络

`Boundary.py` 将边界条件分成两层：

1. 物理边界条件层：`ResistanceBC`、`FluxBC`、`DynamicRadiationResistanceBC`；
2. 边界区域管理层：`BoundaryRegion`，负责在同一个边界面上叠加多个条件。

### 4.1 边界热流符号约定

`BoundaryRegion.compute_net_flux_for_solver()` 返回的是“流入固体节点”的热流，单位为 W：

```text
Flux > 0  表示外部向固体节点输入热量
Flux < 0  表示固体节点向外部散热
```

求解器再根据边界方向将该热流加入对应节点的能量平衡。

### 4.2 `ResistanceBC`

`ResistanceBC` 是通用热阻型边界，可表达对流、指定环境温度、固固接触、等效 Dirichlet 边界等。

单个热阻边界的形式为：

```text
Flux = (T_ext - T_node) / (R_internal + R_ext + R_add)
```

其中：

- `T_ext`：外部等效温度；
- `R_ext`：外部热阻；
- `R_add`：附加热阻，如污垢热阻、接触热阻；
- `R_internal`：固体节点中心到边界面的内部导热热阻，由求解器计算并传入。

对流边界通常通过 `BoundaryRegion.add_convection_condition(T_fluid, h_coeff)` 添加，内部会转换为：

```text
R_conv = 1 / (h * A)
```

### 4.3 `FluxBC`

`FluxBC` 表示固定热流边界：

```text
Flux = q_flux
```

这里的 `q_flux` 是该边界离散面上的总热流，单位 W，而不是 W/m2 的热流密度。

### 4.4 `DynamicRadiationResistanceBC`

`DynamicRadiationResistanceBC` 用热阻网络表达非线性辐射边界。辐射热流原本具有四次方形式：

```text
Q_rad = epsilon * sigma * A * (T_env^4 - T_surface^4)
```

代码中将其局部线性化为等效热导：

```text
h_rad = epsilon * sigma * (T_ref + T_env) * (T_ref^2 + T_env^2)
G_rad = h_rad * A
R_ext = 1 / G_rad
```

随后仍然通过 `ResistanceBC` 的形式进入边界热阻网络。每次 `BoundaryRegion.update_internal_state()` 被调用时，动态辐射边界会根据当前节点温度或表面温度更新 `h_rad`、`G_rad`、`R_ext` 和 `q_flux` 诊断值；它在叠加时仍按热阻型边界处理。

该实现还包含低温保护：

- 试探温度低于 0 K 时记录诊断标记；
- 用 `max(T, 1e-3)` 避免四次方线性化中的非法温度；
- 对 `nan`、`inf` 做数值清理。

### 4.5 `BoundaryRegion`

`BoundaryRegion` 是边界面管理器。一个 `BoundaryRegion` 可以包含多个并联的边界条件，例如：

- 对流换热；
- 外表面辐射；
- 固定热流；
- 接触热阻；
- 耦合器动态更新的等效边界。

其核心思想是诺顿等效/戴维南等效的热网络叠加：

```text
G_total = sum(1 / R_i)
J_total = sum(T_i / R_i) + sum(Q_flux_j)

R_eff = 1 / G_total
T_eff = J_total / G_total + Q_flux * R_eff
Flux = (T_eff - T_node) / (R_eff + R_internal)
T_surface = T_node + Flux * R_internal
```

关键状态：

| 属性 | 含义 |
| --- | --- |
| `T_adj_node` | 与边界相邻的固体节点温度 |
| `R_internal` | 节点中心到边界面的内部热阻 |
| `T_surface` | 当前估算的边界表面温度 |
| `current_flux` | 当前边界净热流，流入固体为正 |
| `conditions` | 已挂载的边界条件列表 |
| `G_sum`、`J_sum` | 并联叠加缓存 |

常用接口：

```python
region.add_resistance_condition(T_ext, R_ext, R_add=0.0)
region.add_convection_condition(T_fluid, h_coeff)
region.add_flux_condition(q_flux)
region.add_dynamic_radiation_condition(emissivity, bare_area_array, T_env)

region.update_internal_state(T_node, R_int, current_time)
region.compute_net_flux_for_solver()

region.get_coupling_snapshot()
region.get_coupling_surface_snapshot()
```

`BoundaryRegion` 构造时默认添加一个极大热阻的弱边界，相当于近似绝热，避免初始状态为空。

## 5. `HeatConduction.py`：导热求解器

### 5.1 `BaseHeatConduction`

`BaseHeatConduction` 抽取了 1D 与 2D 的共同逻辑，包括：

- 状态向量 `T`、`dTdt`；
- 物性缓存 `k_node`、`rho_node`、`cp_node`；
- 热容缓存 `thermal_capacitance`；
- 内热源 `Q_source`；
- 边界容器 `boundaries`；
- ODE RHS 接口 `get_derivatives()`；
- 时间推进接口 `step()`；
- 状态持久化接口。

物性更新逻辑：

```python
self.k_node[:] = self.material.conductivity(self.T)
self.rho_node[:] = self.material.density(self.T)
self.cp_node[:] = self.material.heat_capacity(self.T)
self.thermal_capacitance[:] = self.rho_node * self.cp_node * volumes
```

热源支持两种模式：

1. 回调函数模式：

```python
solid.set_source_term(lambda t, T: source_array)
```

2. 外部数组绑定模式：

```python
solid.link_source_buffer(external_buffer)
```

第二种方式适合与反应堆功率分布、TEC 热源等高频耦合场景共享内存。

时间推进由 `solve_ivp` 完成。代码默认使用 BDF 隐式方法，并在 2D 情况下可注入稀疏雅可比结构。

`step()` 必须把温度初值副本传给积分器：

```python
y0=self.T.copy()
```

不能直接传 `self.T` 本体。`get_derivatives()` 会把积分器试探状态原地写回 `self.T`；若 `y0` 与 `self.T` 共享底层数组，试探点评估可能污染积分初值并制造伪储能。2026-06-01 的单 TFE 严格绝热守恒诊断曾由此出现 kW 到 MW 级残差。

### 5.2 `HeatConduction1D`

`HeatConduction1D` 适用于一维固体导热。默认创建两个边界：

| 边界名 | 位置 | 节点 |
| --- | --- | --- |
| `inner` | 左边界/内边界 | `T[0]` |
| `outer` | 右边界/外边界 | `T[-1]` |

内部缓存：

```python
self.R_geom = np.zeros(self.N + 1)
self.interface_flux = np.zeros(self.N + 1)
```

`R_geom` 包含全部界面热阻：

```text
R_geom[0]      内边界到第一个节点的内部热阻
R_geom[1:-1]   相邻节点之间的内部界面热阻
R_geom[-1]     最后一个节点到外边界的内部热阻
```

笛卡尔坐标下：

```text
R = dx / (k * A)
```

圆柱坐标下，内部与边界热阻使用对数形式：

```text
R = ln(r_2 / r_1) / (2 * pi * k * H)
```

相邻节点界面的导热系数采用调和平均：

```text
k_interface = 2 * k_left * k_right / (k_left + k_right)
```

热流计算逻辑：

```text
内部界面:
  Flux_i = (T_i - T_{i+1}) / R_i

内边界:
  interface_flux[0] = q_inner_in

外边界:
  interface_flux[-1] = -q_outer_in

节点净热流:
  Q_net_i = interface_flux[i] - interface_flux[i+1]
```

其中 `interface_flux` 的正方向定义为从左到右。

### 5.3 `HeatConduction2D`

`HeatConduction2D` 适用于二维平面或轴对称 `r-z` 导热。默认创建四个边界：

| 边界名 | 位置 | shape |
| --- | --- | --- |
| `left` | X/R 最小侧 | `(n_y,)` |
| `right` | X/R 最大侧 | `(n_y,)` |
| `bottom` | Y/Z 最小侧 | `(n_x,)` |
| `top` | Y/Z 最大侧 | `(n_x,)` |

二维求解器不显式保存内部热阻，而是保存热导：

```python
self.G_x_inner  # shape = (n_x-1, n_y)
self.G_y_inner  # shape = (n_x, n_y-1)
```

这样可以避免反复除法，并降低 `k=0` 或距离极小导致的数值问题。

内部热导计算：

```text
k_face = 2 * k_left * k_right / (k_left + k_right)
G_x = k_face * area_x / dx
G_y = k_face * area_y / dy
```

净热流组装采用五点模板：

```text
X方向:
  q_x = (T_left - T_right) * G_x
  左节点减 q_x，右节点加 q_x

Y方向:
  q_y = (T_bottom - T_top) * G_y
  下节点减 q_y，上节点加 q_y

边界:
  对应边界返回的 q_in 直接加到边界节点
```

注意：按照当前代码的符号实现，内部导热项中 `q_x = (T_left - T_right) * G_x` 后，对左节点执行减法、对右节点执行加法。这相当于把从高温节点流向低温节点的热流从出流侧扣除、在入流侧增加，最终形成每个节点的净热流。

### 5.4 稀疏雅可比结构

`HeatConduction2D.get_jac_sparsity()` 返回隐式 ODE 求解器使用的稀疏雅可比非零结构。对于二维五点格式，一个节点最多与以下节点耦合：

```text
自身:      offset = 0
Y方向邻居: offset = -1, +1
X方向邻居: offset = -n_y, +n_y
```

当 `n_x == 1` 或 `n_y == 1` 时，函数会自动退化为较低维的稀疏结构。

该优化对 BDF/Radau 等隐式方法很重要，可以减少有限差分扰动次数和稠密矩阵处理成本。

## 6. 与 Couplers 和 SystemManager 的关系

HeatConduction 模块通常不单独运行，而是被 `SystemManager` 调度，并通过 `Couplers.py` 与其他物理场交换信息。

### 6.1 流固耦合

`FluidSolidCouple` 负责连接水力通道和固体边界，其主要工作是：

```text
1. 从流体通道读取 T_f, P_f, rho, velocity, mu, k_f, cp
2. 计算 Re, Pr, Nu
3. 得到 h 和 lambda = h * A
4. 更新固体边界:
   T_ext = T_fluid
   R_ext = 1 / lambda
5. 读取固体壁面温度 T_wall
6. 向流体侧写入半隐式热源:
   Q = lambda * T_wall - lambda * T_fluid
```

固体侧通过 `BoundaryRegion` 暴露：

```python
region.get_coupling_snapshot()
region.get_coupling_surface_snapshot()
```

### 6.2 固固耦合与间隙耦合

固固耦合器通过两个固体边界互相交换等效温度和热阻：

```text
solid1.boundary.T_ext = solid2.surface_or_node_temperature
solid1.boundary.R_ext = solid2.R_internal + optional_gap_resistance

solid2.boundary.T_ext = solid1.surface_or_node_temperature
solid2.boundary.R_ext = solid1.R_internal + optional_gap_resistance
```

`GapCouple2D` 在此基础上增加：

- 间隙气体导热；
- 表面对表面辐射；
- 温度相关的等效间隙热阻。

`ActiveGapCouple2D` 和 `TECCouple2D` 进一步将固定热流源转化为等效戴维南温差，适用于带电子冷却/加热的热离子电极间隙。

### 6.3 `SystemManager` 调度流程中的位置

在 `SystemManager.step(dt)` 中，固体导热通常位于 Picard 内迭代中：

```text
1. 组件 pre_step 更新功率、TEC、热管等外部源项
2. 运行 couplers，刷新固体边界条件和流体源项
3. 推进中子学
4. 推进水力网络
5. 对每个 solid 调用 solid.step(dt)
6. 检查流体和固体温度收敛
```

因此，HeatConduction 的边界条件不是静态输入，而是会在每个耦合迭代中被动态刷新。

## 7. 典型使用流程

### 7.1 创建 1D 圆柱导热对象

```python
from Solvers.HeatConduction.Mesh import Mesh1D
from Solvers.HeatConduction.HeatConduction import HeatConduction1D

mesh = Mesh1D(
    total_dim=0.005,
    n_volumes=30,
    geometry_type="cylindrical",
    inner_radius=0.0,
    height=1.0,
)

solid = HeatConduction1D(
    mesh=mesh,
    material=solid_material,
    initial_temp=900.0,
    name="fuel_pin",
)
```

### 7.2 添加外边界对流

```python
outer = solid.boundaries["outer"]
outer.clear_conditions()
outer.add_convection_condition(T_fluid=700.0, h_coeff=5000.0)
```

### 7.3 添加外表面辐射

```python
outer.add_dynamic_radiation_condition(
    emissivity=0.8,
    bare_area_array=outer.area,
    T_env=300.0,
)
```

### 7.4 设置体热源并推进

```python
def source_func(t, T):
    return np.full_like(T, 100.0)  # 每个控制体 100 W

solid.set_source_term(source_func)
ok = solid.step(dt=0.01, rtol=1e-5, atol=1e-7)
```

### 7.5 创建 2D 轴对称导热对象

```python
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D

mesh = Mesh2D(
    x_dim=0.004,
    n_x=20,
    y_dim=0.5,
    n_y=60,
    geometry_type="cylindrical",
    inner_radius=0.0,
)

solid = HeatConduction2D(
    mesh=mesh,
    material=solid_material,
    initial_temp=800.0,
    name="rz_solid",
)
```

## 8. 重要实现细节与注意事项

### 8.1 热流单位

边界条件中的 `q_flux`、`current_flux`、`Q_source` 都按离散控制体或离散边界面的总热功率处理，单位是 W。若输入来自热流密度 W/m2，需要先乘以对应边界面积。

### 8.2 边界默认近似绝热

`BoundaryRegion` 初始化时会添加一个 `T_ext=300 K`、`R_ext=1e15 K/W` 的默认热阻边界。实际使用时，如果希望边界条件完全由用户定义，建议先调用：

```python
solid.boundaries["outer"].clear_conditions()
```

然后再添加目标边界条件。

### 8.2.1 边界叠加缓存

`BoundaryRegion._accumulate_bc()` 仍在使用：`add_resistance_condition()`、`add_convection_condition()`、`add_flux_condition()` 和 `add_dynamic_radiation_condition()` 添加条件时会调用它，为 `G_sum/J_sum/Q_sum_flux` 提供初始化缓存。正式求解路径中的 `compute_net_flux_for_solver()` 会在每次计算边界热流时重新清零并重建这些缓存；隐式导热矩阵装配随后读取 `G_sum/R_eff/T_eff/Q_sum_flux`。

热阻到电导的转换遵循以下约定：

```text
R_ext = inf  -> G = 0       # 绝热
R_ext = 0    -> G = 1e20    # 定温/Dirichlet 极限
R_ext > 0    -> G = 1/R_ext
```

因此 `R_ext=0` 是有效的定温边界表达，不应被当作绝热。固定热流 `FluxBC` 进入 `Q_sum_flux`，不进入 `J_sum`。

### 8.3 物性函数需要支持向量输入

`BaseHeatConduction._update_properties()` 直接向材料对象传入数组：

```python
material.conductivity(self.T)
material.density(self.T)
material.heat_capacity(self.T)
```

因此 `SolidMaterial` 的实现应支持 NumPy 向量化输入。

### 8.4 2D 数据布局

2D 温度数组的扁平化规则是：

```text
k = i * n_y + j
```

这意味着 `T.reshape((n_x, n_y))` 后：

- 第一维是 X/R 方向；
- 第二维是 Y/Z 方向；
- 相邻 Y/Z 节点在扁平数组中相邻；
- 相邻 X/R 节点在扁平数组中的间隔是 `n_y`。

### 8.5 时间积分方法

`BaseHeatConduction` now supports a per-solid default ODE method through `solid.ode_method` and `solid.set_ode_method(method)`. The default is still `BDF`, so existing `SystemManager` calls to `solid.step(dt)` keep their previous behavior unless a component sets a different method on that solid.

`BaseHeatConduction.step(dt, method=None, **kwargs)` uses `self.ode_method` when `method is None`; passing `method` explicitly overrides only that single call. Supported methods are the SciPy `solve_ivp` methods `RK45`, `RK23`, `DOP853`, `Radau`, `BDF`, and `LSODA`, plus `implicit_euler` for generic `HeatConduction1D/2D`. `BDF` and `Radau` continue to receive available sparse Jacobian structure.

`implicit_euler` performs a backward-Euler sparse algebraic solve once per global `SystemManager.step()` time step. In this mode, resistance/convection/dynamic-radiation boundary terms are assembled into the implicit matrix, while pure `FluxBC` terms stay on the explicit RHS.

Implementation note for the optimized path: generic `HeatConduction1D/2D` cache the CSC sparse matrix pattern after the first implicit step for a fixed mesh shape. Later steps update only numeric `data` entries, including conductance and boundary-linearization diagonal terms. This removes Python edge-loop assembly from the hot path while keeping material and radiation coefficients free to change each step. There is no cross-solid assembly and no LU/factorization cache in this round.

If implicit solve fails, `HeatConduction.step()` warns and falls back to `solve_ivp`.

Because `implicit_euler` is first-order and intentionally dissipative, compare it against `solve_ivp` with the same global `max_dt` when validating system cases. In the 2026-06-24 V13 no-TEC 1 s smoke, `max_dt=0.1 s` kept the core outlet difference to about `0.11 K`; `max_dt=0.5 s` showed visibly larger backward-Euler damping.

For local performance regression checks, use:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\benchmark_heatconduction_implicit.py --steps-1d 40 --steps-2d 20 --n1d 4000 --n2d-x 90 --n2d-y 90
```

2026-06-24 post-optimization reference on this machine: 1D 4000-node implicit stepping `0.001869 s/step`; 2D 90x90 implicit stepping `0.026584 s/step`.

V8 CaseA 的公共入口 `testModule/run_v8_caseA_common.py` 会在构建后和加载 restart 后统一调用 `solid.set_ode_method()`。默认 `--solid-ode-method` 为 `LSODA`，并覆盖 `SystemManager.solid_components` 中注册的全部堆芯固体，包括 TFE 内部固体、全局慢化剂环、筒体、反射层和网格化间隙固体。显式传入 `BDF`、`Radau` 等合法方法可恢复或切换该算例的固体积分器。

### 8.6 初始化同步

`HeatConduction1D` 和 `HeatConduction2D` 构造结束时都会调用 `initialize_state()`，用于：

1. 根据初始温度更新物性；
2. 计算内部热阻/热导；
3. 将边界相邻节点温度和内部热阻推送给 `BoundaryRegion`；
4. 计算一次边界热流和表面温度。

这避免了耦合器在第一个时间步读取到未初始化的边界缓存。

### 8.7 动态辐射边界允许空数组

`DynamicRadiationResistanceBC.update_state()` 会先检查参考温度数组大小。对于零长度边界数组，函数直接返回，不调用 `np.min()`：

```python
if t_ref_raw.size == 0:
    return
```

因此，空边界数组是允许的边缘情况，不应被当作动态辐射异常。

## 9. 模块优点

HeatConduction 当前实现具有以下特点：

- 1D 和 2D 共享统一的 ODE、物性和热源接口；
- 网格类单独负责几何量，求解器只消费面积、体积和距离；
- 边界条件通过热阻网络统一表达，便于叠加对流、辐射、接触和固定热流；
- 2D 求解器使用热导缓存和预分配数组，减少高频 RHS 调用中的内存分配；
- 支持非均匀网格，适合边界层、间隙、电极等局部加密场景；
- 支持稀疏雅可比结构，适合隐式求解较大的固体导热问题；
- 与 `Couplers.py` 和 `SystemManager.py` 的接口清晰，适合多物理场 Picard 内迭代。

## 10. 推荐阅读顺序

若需要进一步维护或扩展该模块，建议按以下顺序阅读源码：

1. `HeatConduction/Mesh.py`
   - 先理解控制体体积、界面面积和距离的生成方式；
2. `HeatConduction/Boundary.py`
   - 理解 `BoundaryRegion` 如何叠加多个边界条件；
3. `HeatConduction/HeatConduction.py`
   - 阅读 `BaseHeatConduction.get_derivatives()`，再分别阅读 `HeatConduction1D` 和 `HeatConduction2D`；
4. `Couplers.py`
   - 查看流固、固固、间隙和 TEC 耦合如何动态更新边界；
5. `SystemManager.py`
   - 理解 HeatConduction 在全局 Picard 内迭代中的执行时机。

## 11. 源码校验信息

最后按源码校验日期：2026-05-31。

涉及文件：

- `HeatConduction/Mesh.py`
- `HeatConduction/Boundary.py`
- `HeatConduction/HeatConduction.py`
- `Couplers.py`
- `SystemManager.py`

## 12. 2026-06-08 HeatPipe2D heat-pipe-only implicit boundary linearization

`DynamicRadiationResistanceBC` keeps its generic `BoundaryRegion` behavior: `compute_net_flux_for_solver()` still evaluates the thermal-resistance network and returns the heat flow into the adjacent solid node. `Components/basicComponents/HeatPipe2D.py` adds a heat-pipe-only option for the sparse `implicit_euler` and `theta_implicit` paths:

```python
hp.set_implicit_boundary_linearization(True)
```

The option is disabled by default and does not affect BDF, generic `HeatConduction1D/2D`, reactor-core radiation boundaries, or shared `BoundaryRegion` semantics. When enabled, heat-pipe resistance-type boundary terms are assembled into the heat-pipe matrix as `Q = G * (T_eff - T_node)`. This includes bare-wall `DynamicRadiationResistanceBC` and the `HPwithFin` equivalent fin-branch `ResistanceBC`; flux-only terms remain explicit RHS contributions.
