# Couplers 多物理场耦合器详细介绍

## 1. 模块定位

`Couplers.py` 是 TASTIN `Solvers` 层中的多物理场耦合桥梁。它不直接求解某一个物理场的 ODE 或线性方程，而是负责在不同求解器之间同步边界条件、热源项和等效热阻。

它连接的主要对象包括：

- `HeatConduction` 固体导热求解器；
- `Hydrodynamics` 水力网络和流体通道；
- 热间隙、辐射、气体导热、热电转换等界面物理；
- `SystemManager` 的 Picard 内迭代调度。

当前 `Couplers.py` 中包含以下耦合器：

| 类 | 作用 |
| --- | --- |
| `SolidSolidCouple1D` | 一维固体-固体边界耦合 |
| `SolidSolidCouple2D` | 二维固体-固体边界耦合 |
| `FluidSolidCouple` | 流体-固体对流换热耦合 |
| `GapCouple2D` | 二维间隙耦合，辐射与气体导热并联 |
| `ActiveGapCouple2D` | 带固定热流源的有源间隙耦合 |
| `TECCouple2D` | 热电转换耦合，支持两侧非对称热流源 |

## 2. 总体设计思想

Couplers 的核心思想是：不直接修改求解器内部离散方程，而是通过各求解器公开的边界和源项接口实现物理耦合。

典型数据流如下：

```text
HeatConduction.BoundaryRegion  <->  SolidSolidCouple  <->  HeatConduction.BoundaryRegion

Hydrodynamics.FluidChannel     <->  FluidSolidCouple  <->  HeatConduction.BoundaryRegion

HeatConduction.BoundaryRegion  <->  Gap/TEC Coupler   <->  HeatConduction.BoundaryRegion
```

所有固体侧耦合最终都落到 `BoundaryRegion.add_resistance_condition()` 创建的 `ResistanceBC` 上，通过动态更新：

```text
T_ext
R_ext
R_add
```

来表达对方节点温度、对方内部热阻、接触热阻、间隙热阻或等效热源。

流体侧耦合则通过 `FluidChannel.add_coupling_source_distribution()` 写入半隐式热源：

```text
Q_fluid = explicit_part - implicit_factor * T_fluid
```

## 3. 与 SystemManager 的关系

根据 `SOLVERS_ANALYSIS.md` 和 `AI_AGENT_SOLVERS_ANALYSIS.md` 的说明，耦合器通常在 `SystemManager.step()` 的 Picard 内迭代中执行。

典型顺序：

```text
SystemManager.step(dt)
  1. 保存入口状态
  2. component.pre_step(dt, t_start)
  3. Picard 内迭代:
     3.1 恢复流体源项基准
     3.2 _run_couplers(interface_relaxation)
         - 刷新固体边界缓存
         - SolidSolidCouple.sync()
         - GapCouple2D.sync()
         - 再次刷新固体边界缓存
         - FluidSolidCouple.execute()
     3.3 推进中子学
     3.4 应用核功率/衰变热源
     3.5 HydraulicNetwork.step_Picard(dt)
     3.6 HeatConduction.step(dt)
     3.7 检查流体/固体收敛
  4. 收敛后提交状态
```

因此，Couplers 的输出通常只在当前 Picard 迭代内有效，下一次迭代会根据新的温度、流量和边界状态重新计算。

## 4. `SolidSolidCouple1D`

`SolidSolidCouple1D` 用于连接两个 `HeatConduction1D` 对象的边界。

### 4.1 方向映射

构造参数：

```python
SolidSolidCouple1D(obj1, obj2, direction, contact_resistance=0.0)
```

方向含义：

| `direction` | 连接关系 |
| --- | --- |
| `right` / `outer` | `obj1.outer` 连接 `obj2.inner` |
| `left` / `inner` | `obj1.inner` 连接 `obj2.outer` |

### 4.2 初始化逻辑

构造时会取出两个固体对象的 `BoundaryRegion`：

```python
self.bound1 = obj1.boundaries[self.loc1]
self.bound2 = obj2.boundaries[self.loc2]
```

然后分别挂载一个热阻型边界条件：

```python
self.bc1 = self.bound1.add_resistance_condition(T_ext=300.0, R_ext=0.0, R_add=R_contact)
self.bc2 = self.bound2.add_resistance_condition(T_ext=300.0, R_ext=0.0, R_add=R_contact)
```

### 4.3 同步逻辑

`sync()` 会交叉更新两侧边界：

```text
obj1 边界看到:
  T_ext = obj2 边界相邻节点/表面温度
  R_ext = obj2 内部热阻

obj2 边界看到:
  T_ext = obj1 边界相邻节点/表面温度
  R_ext = obj1 内部热阻
```

源码中使用：

```python
T1_surf, R1_int = bound1.get_coupling_snapshot()
T2_surf, R2_int = bound2.get_coupling_snapshot()
bc1.update_params(T_ext=T2_surf, R_ext=R2_int)
bc2.update_params(T_ext=T1_surf, R_ext=R1_int)
```

其中 `contact_resistance` 通过 `R_add` 进入边界条件。

## 5. `SolidSolidCouple2D`

`SolidSolidCouple2D` 用于连接两个 `HeatConduction2D` 对象的整条边界。

### 5.1 方向映射

构造参数：

```python
SolidSolidCouple2D(obj1, obj2, direction, contact_resistance=0.0)
```

方向映射：

| `direction` | 连接关系 |
| --- | --- |
| `right` | `obj1.right` 连接 `obj2.left` |
| `left` | `obj1.left` 连接 `obj2.right` |
| `top` | `obj1.top` 连接 `obj2.bottom` |
| `bottom` | `obj1.bottom` 连接 `obj2.top` |

### 5.2 维度校验

二维耦合必须保证两侧边界节点数相同：

```python
if self.bound1.shape != self.bound2.shape:
    raise ValueError(...)
```

这保证了后续数组化同步可以一一对应。

### 5.3 同步逻辑

`sync()` 与 1D 类似，但直接对整条边界数组操作：

```python
T1_surf, R1_int = self.bound1.get_coupling_snapshot()
T2_surf, R2_int = self.bound2.get_coupling_snapshot()

self.bc1.update_params(T_ext=T2_surf, R_ext=R2_int)
self.bc2.update_params(T_ext=T1_surf, R_ext=R1_int)
```

## 6. `FluidSolidCouple`

`FluidSolidCouple` 是流固热耦合器，用于把流体通道和固体边界通过对流换热连接起来。

### 6.1 构造参数

```python
FluidSolidCouple(
    name,
    fluid,
    solid_boundary_region,
    heated_perimeter,
    correlation_func,
    solid_node_capacitance=None,
)
```

参数说明：

| 参数 | 含义 |
| --- | --- |
| `name` | 耦合器名称 |
| `fluid` | `FluidChannel` 或兼容对象 |
| `solid_boundary_region` | 固体边界 `BoundaryRegion` |
| `heated_perimeter` | 湿周，单位 m |
| `correlation_func` | Nu 数关联式函数 |
| `solid_node_capacitance` | 固体边界节点热容，用于稳定步长估计 |

`correlation_func` 的接口约定为：

```python
Nu = correlation_func(Re, Pr, P_D_ratio)
```

### 6.2 网格一致性

流体节点数必须与固体边界节点数一致：

```python
solid_bound.shape[0] == fluid.n_nodes
```

如果提供 `solid_node_capacitance`，其长度也必须等于 `fluid.n_nodes`。

### 6.3 换热面积

每个流体控制体的换热面积为：

```text
A_node = heated_perimeter * fluid.node_length
```

对于非均匀通道，`fluid.node_length` 可以是数组，因此 `node_areas` 也自然成为数组。

### 6.4 固体侧边界条件

构造时会在固体边界上创建一个 `ResistanceBC`：

```python
self.solid_bc = solid_bound.add_resistance_condition(
    T_ext=T_fluid_init,
    R_ext=R_init,
)
```

后续每次 `execute()` 都会更新：

```text
T_ext = T_fluid
R_ext = 1 / (h * A)
```

对应物理对流：

```text
Q_solid = h * A * (T_fluid - T_wall)
```

### 6.5 执行流程

`execute(interface_relaxation=1.0)` 的核心流程：

```text
1. 读取流体温度、压力、密度、速度
2. 从流体材料计算 mu、k、Cp、Pr
3. 计算 Re
4. 调用 correlation_func 得到 Nu
5. 计算 h = Nu * k / D_h
6. 计算 lambda = h * A
7. 固体侧:
   R_ext = 1 / lambda
   T_ext = T_fluid
8. 读取固体壁面温度 T_wall
9. 流体侧:
   explicit = lambda * T_wall
   implicit = lambda
   fluid.add_coupling_source_distribution(explicit, implicit)
10. 记录耦合诊断量和缓存 lambda
```

### 6.6 半隐式源项

流体方程中写入的源项为：

```text
Q_fluid = explicit - implicit * T_fluid
        = lambda * T_wall - lambda * T_fluid
        = h * A * (T_wall - T_fluid)
```

这与固体侧对流边界方向相反，保证界面热量交换方向一致。

### 6.7 界面松弛

`execute()` 支持界面松弛：

```python
execute(interface_relaxation=0.5)
```

当 `interface_relaxation < 1.0` 且存在上一轮缓存状态时，会对以下界面状态做松弛：

```text
lambda
T_f
T_wall
```

松弛形式：

```text
x_relaxed = alpha * x_current + (1 - alpha) * x_previous
```

缓存状态保存在：

```python
self._interface_relaxation_previous
```

可通过：

```python
reset_interface_relaxation()
```

清除上一轮松弛记忆。

### 6.8 耦合诊断量

`FluidSolidCouple` 会记录最近一次耦合诊断：

```python
get_coupling_diagnostics()
```

包含：

- `interface_residual`；
- `source_residual`；
- `delta_lambda`；
- `delta_T_f`；
- `delta_T_wall`；
- `delta_explicit`；
- `delta_implicit`；
- `max_lambda`；
- `max_explicit`；
- 是否实际执行了松弛。

这些信息可用于 `SystemManager` 评估界面收敛。

### 6.9 稳定步长估计

`get_max_stable_dt()` 基于流固界面热容和耦合热导估计稳定步长：

```text
dt < C_eff / lambda
```

其中当前实现使用固体热容和流体热容的调和组合：

```text
C_eff = C_solid * C_fluid / (C_solid + C_fluid)
```

流体热容：

```text
C_fluid = rho * V_fluid * Cp
V_fluid = flow_area * node_length
```

最终返回：

```text
min(min(C_eff / lambda) * safety_factor, max_limit)
```

## 7. `GapCouple2D`

`GapCouple2D` 继承自 `SolidSolidCouple2D`，用于二维固体之间存在间隙时的耦合。

它同时考虑：

1. 表面对表面辐射；
2. 间隙气体导热。

二者为并联关系。

### 7.1 构造参数

```python
GapCouple2D(
    obj1,
    obj2,
    direction,
    gap_width,
    gas_conductivity,
    emissivity1=0.8,
    emissivity2=0.8,
)
```

参数说明：

| 参数 | 含义 |
| --- | --- |
| `gap_width` | 间隙宽度，单位 m |
| `gas_conductivity` | 间隙气体导热系数，单位 W/(m K) |
| `emissivity1` | obj1 表面发射率 |
| `emissivity2` | obj2 表面发射率 |

### 7.2 辐射电导

辐射热流被线性化为：

```text
h_rad* = sigma * (T1^2 + T2^2) * (T1 + T2)
```

考虑内外面积和两侧发射率后：

```text
G_rad = A_in * h_rad* / denom_rad
```

其中：

```text
denom_rad = 1/eps_in + (1/eps_out - 1) * A_in/A_out
```

代码通过面积大小判断哪一侧是内表面：

```python
is_1_inner = A1 < A2
```

### 7.3 气体导热电导

间隙气体导热采用薄间隙近似：

```text
G_cond = k_gas * A_in / gap_width
```

### 7.4 总间隙热阻

辐射与气体导热并联：

```text
G_total = G_rad + G_cond
R_gap_total = 1 / G_total
```

该结果缓存于：

```python
self.R_gap_total
```

### 7.5 边界更新

`GapCouple2D.sync()` 获取两侧节点状态后，更新两侧边界：

```text
obj1 看到:
  T_ext = T2_node
  R_ext = R_gap_total + R2_int

obj2 看到:
  T_ext = T1_node
  R_ext = R_gap_total + R1_int
```

也就是说，间隙热阻与对方固体的内部热阻串联。

## 8. `ActiveGapCouple2D`

`ActiveGapCouple2D` 继承自 `GapCouple2D`，在辐射和气体导热基础上增加一个分布式固定热流源。

### 8.1 热源接口

```python
set_active_heat_source(Q_array)
```

符号约定：

| `Q_array` | 含义 |
| --- | --- |
| `Q > 0` | 热量从 Obj1 流出，流入 Obj2 |
| `Q < 0` | 热量从 Obj2 流出，流入 Obj1 |

`Q_array` 可以是标量，也可以是与边界 shape 一致的数组，单位为 W。

### 8.2 戴维南等效

`sync()` 先调用父类 `GapCouple2D.sync()` 得到基础热阻，然后将固定热流源转化为等效温差：

```text
delta_T = Q_source * R_gap
```

更新方式：

```text
Obj1:
  T_ext = T2_node - delta_T

Obj2:
  T_ext = T1_node + delta_T
```

热阻保持父类计算结果不变。

这种处理避免在温差接近 0 时直接用热流反算等效热阻导致数值奇异。

## 9. `TECCouple2D`

`TECCouple2D` 继承自 `GapCouple2D`，用于热电转换或热离子电极间隙，支持两侧独立且非对称的表面热流源。

### 9.1 热源接口

```python
set_tec_sources(Q_emitter, Q_collector)
```

参数含义：

| 参数 | 作用侧 | 典型符号 |
| --- | --- | --- |
| `Q_emitter` | Obj1 / Emitter 侧 | 电子冷却通常为负 |
| `Q_collector` | Obj2 / Collector 侧 | 电子加热通常为正 |

源码中的符号约定：

```text
Q > 0: 热量流入固体
Q < 0: 热量流出固体
```

输入可以是标量或与边界 shape 一致的数组。

### 9.2 同步逻辑

`sync()` 先调用 `GapCouple2D.sync()`，获得基础间隙热阻、辐射和气体导热边界。若两侧热源全为 0，则直接返回。

非对称热源通过等效外部温度修正：

```text
T_ext_1_new =
  T_node_2
  + Q_source_1 * R_ext_1
  + Q_source_2 * R_int_2

T_ext_2_new =
  T_node_1
  + Q_source_2 * R_ext_2
  + Q_source_1 * R_int_1
```

其中：

- `R_ext_1` 是 Obj1 侧看到的外部总热阻；
- `R_ext_2` 是 Obj2 侧看到的外部总热阻；
- `R_int_1`、`R_int_2` 是两侧内部热阻。

该模型允许 `Q_emitter` 与 `Q_collector` 不相等，因此可以表示非守恒的电子冷却/加热过程。

## 10. 耦合器继承关系

当前类层次如下：

```text
SolidSolidCouple2D
  └── GapCouple2D
        ├── ActiveGapCouple2D
        └── TECCouple2D
```

含义：

- `SolidSolidCouple2D` 只做两侧固体热阻边界交叉同步；
- `GapCouple2D` 增加间隙辐射和气体导热；
- `ActiveGapCouple2D` 增加单个守恒方向的固定热流源；
- `TECCouple2D` 增加两侧独立的非对称热流源。

## 11. 典型使用流程

### 11.1 固固耦合

```python
couple = SolidSolidCouple2D(
    obj1=solid_a,
    obj2=solid_b,
    direction="right",
    contact_resistance=1.0e-4,
)

couple.sync()
```

### 11.2 流固耦合

```python
couple = FluidSolidCouple(
    name="coolant_to_wall",
    fluid=channel,
    solid_boundary_region=solid.boundaries["left"],
    heated_perimeter=perimeter,
    correlation_func=nu_correlation,
    solid_node_capacitance=solid.get_boundary_node_capacitance("left"),
)

couple.execute(interface_relaxation=0.7)
```

### 11.3 间隙耦合

```python
gap = GapCouple2D(
    obj1=emitter,
    obj2=collector,
    direction="right",
    gap_width=1.0e-3,
    gas_conductivity=0.02,
    emissivity1=0.8,
    emissivity2=0.8,
)

gap.sync()
```

### 11.4 TEC 热电耦合

```python
tec = TECCouple2D(
    obj1=emitter,
    obj2=collector,
    direction="right",
    gap_width=1.0e-3,
    gas_conductivity=0.02,
)

tec.set_tec_sources(Q_emitter=q_emit, Q_collector=q_coll)
tec.sync()
```

## 12. 重要实现细节与注意事项

### 12.1 热流单位

`Q_array`、`Q_emitter`、`Q_collector`、`explicit_arr` 都是离散节点上的总热功率，单位 W。若外部结果是热流密度 W/m2，应先乘以对应面积。

### 12.2 热阻单位

`ResistanceBC` 中的 `R_ext` 和 `R_add` 按绝对热阻处理，单位 K/W。若输入的是单位面积热阻 m2 K/W，需要除以或换算到每个边界节点的面积。

### 12.3 `sync()` 与 `execute()` 的调用时机

- 固固、间隙、TEC 耦合器使用 `sync()`；
- 流固耦合器使用 `execute()`；
- 它们通常应在固体和流体求解器推进之前调用；
- 在 Picard 内迭代中，每轮都应重新调用。

### 12.4 流固耦合前应清理流体源项

`FluidSolidCouple.execute()` 会向流体通道累加源项。因此上层调度器通常需要在每轮 Picard 迭代开始时恢复或清空流体源项基准，避免重复累加。

### 12.5 动态边界依赖固体边界缓存

`FluidSolidCouple` 会调用 `solid_bound.compute_net_flux_for_solver()` 来刷新壁面温度。上层调度中通常还会在耦合前后刷新固体边界缓存，保证 `T_surface`、`R_internal` 和 `T_adj_node` 与当前温度场一致。

### 12.6 间隙宽度不能为零

`GapCouple2D` 中：

```text
G_cond = k_gas * A_in / gap_width
```

因此 `gap_width` 应为正数。若需要接触导热，应使用 `SolidSolidCouple2D` 的 `contact_resistance` 或显式接触热阻模型。

### 12.7 发射率为零的处理

辐射计算中若发射率为 0，会产生无穷辐射热阻，代码通过 `np.nan_to_num` 将对应辐射电导清理为 0。

### 12.8 界面松弛只在形状匹配且有历史状态时生效

`FluidSolidCouple` 的界面松弛要求上一次缓存状态与当前状态键和数组形状一致。否则即使传入 `interface_relaxation < 1.0`，也会使用当前状态。

## 13. 模块优点

`Couplers.py` 的实现具有以下特点：

- 使用边界热阻网络统一表达固固、间隙、辐射、气体导热和热电源项；
- 流固耦合采用半隐式源项，降低强换热导致的数值刚性；
- 支持向量化整条边界耦合，适合轴向分布模型；
- 支持 Picard 内迭代界面松弛和耦合残差诊断；
- 支持固体/流体热容共同约束的稳定步长估计；
- 支持代表热电转换的非对称非守恒热源；
- 与 `HeatConduction.BoundaryRegion` 和 `Hydrodynamics.FluidChannel` 的接口边界清晰。

## 14. 推荐阅读顺序

维护或扩展耦合器时，建议按以下顺序阅读：

1. `HeatConduction/Boundary.py`
   - 理解 `BoundaryRegion`、`ResistanceBC`、`get_coupling_snapshot()`；
2. `Couplers.py`
   - 先看 `SolidSolidCouple2D`，再看 `GapCouple2D`、`ActiveGapCouple2D` 和 `TECCouple2D`；
3. `Hydrodynamics/Components.py`
   - 理解 `FluidChannel.add_coupling_source_distribution()`；
4. `Hydrodynamics/HydraulicNetwork.py`
   - 理解流体半隐式源项如何进入能量方程；
5. `SystemManager.py`
   - 理解耦合器在 Picard 内迭代中的执行顺序。

## 15. 源码校验信息

最后按源码校验日期：2026-05-31。

涉及文件：

- `Couplers.py`
- `HeatConduction/Boundary.py`
- `HeatConduction/HeatConduction.py`
- `Hydrodynamics/Components.py`
- `Hydrodynamics/HydraulicNetwork.py`
- `SystemManager.py`
