# SystemManager 全局调度器详细介绍

## 1. 模块定位

`SystemManager.py` 是 TASTIN `Solvers` 层的全局瞬态调度器。它本身不直接求解某一个物理场方程，而是负责把水力网络、固体导热、中子动力学、耦合器和宏观组件组织到同一个时间步推进流程中。

它协调的主要对象包括：

- `HydraulicNetwork`：流体压力、流量、焓和温度求解；
- `BaseHeatConduction` 派生类：固体导热求解；
- `Couplers.py` 中的耦合器：流固、固固、间隙、TEC 等界面耦合；
- `PointReactor`：点堆中子动力学；
- `BaseComponent` 派生宏观组件：如堆芯、热电电路、热管等组件级逻辑。

`SystemManager` 的核心职责是：

- 管理全局时间 `global_time`；
- 注册固体、耦合器、宏观组件和点堆对象；
- 初始化水力场和初始耦合状态；
- 在每个时间步执行 Picard 内迭代；
- 管理流体源项清零、恢复和持久源项；
- 调度耦合器同步边界条件；
- 调度中子学功率计算和功率写入；
- 推进流体和固体求解器；
- 检查温度和界面收敛；
- 失败时回滚系统状态；
- 提供自适应时间步长和全局断点续算。

## 2. 类结构概览

当前文件只定义一个核心类：

```python
class SystemManager:
    ...
```

构造函数：

```python
SystemManager(fluid_network: HydraulicNetwork, start_time: float = 0.0)
```

内部主要成员：

| 成员 | 含义 |
| --- | --- |
| `global_time` | 当前全局时间 |
| `fluid_solver` | 水力网络求解器 |
| `solid_components` | 固体求解器字典，键为固体名称 |
| `couplers` | 耦合器列表 |
| `point_reactor` | 可选点堆求解器 |
| `components` | 宏观组件列表 |
| `_persistent_fluid_sources` | 每轮清源后重新施加的持久流体源 |
| `last_step_diagnostics` | 最近一次时间步诊断信息 |
| `_component_step_states` | 组件回滚状态缓存 |
| `_point_reactor_step_state` | 点堆回滚状态缓存 |

构造时若水力求解器提供 `set_time()`，会同步到 `start_time`。

## 3. 注册接口

### 3.1 注册点堆求解器

```python
add_point_reactor(reactor)
```

将 `PointReactor` 或兼容对象挂到 `SystemManager`。后续如果没有组件自行处理中子学，`SystemManager` 会直接调用：

```python
point_reactor.step(dt, reactivity_control, total_feedback)
point_reactor.commit()
```

### 3.2 注册持久流体源

```python
add_persistent_fluid_source(source_func)
```

`source_func` 必须是可调用对象，签名通常为：

```python
def source_func(system_manager):
    ...
```

持久流体源用于在每轮流体源项清零后重新施加基准源项。例如外部固定加热、背景体积热源或不应被耦合器清除的源项。

### 3.3 注册固体组件

```python
add_solid_component(component: BaseHeatConduction)
```

要求对象必须继承自 `BaseHeatConduction`。注册时会：

- 读取 `component.name`；
- 若名称不存在则自动分配；
- 检查名称唯一性；
- 写入 `solid_components` 字典。

### 3.4 注册耦合器

```python
add_coupler(coupler)
```

有效耦合器至少需要提供：

```text
execute()
```

或：

```text
sync()
```

流固耦合器通常使用 `execute()`，固固/间隙/TEC 耦合器通常使用 `sync()`。

### 3.5 注册宏观组件

```python
add_component(component: BaseComponent)
```

宏观组件注册时，`SystemManager` 会从组件中提取：

```python
component.get_solids()
component.get_couplers()
```

并分别注册到固体列表和耦合器列表中。

组件还可选择实现以下接口：

```text
pre_step()
post_step()
advance_neutronics()
commit_neutronics()
save_step_state()
load_step_state()
get_state_dict()
load_state_dict()
```

这些接口会在全局时间步、回滚和断点续算中被调用。

## 4. 初始化流程

```python
initialize_system(dt_init=0.1, tol=1e-5, max_iter=500)
```

初始化包含两步：

```text
1. fluid_solver.initialize_hydraulics(dt, tol, max_iter)
2. 清理/施加流体源项并运行耦合器同步
```

如果水力初始化失败，会抛出：

```python
RuntimeError("System initialization failed at hydraulic stage.")
```

这一步的目的通常是先建立合理的初始压力/流量分布，再同步固体边界、流体源项和界面状态。

## 5. 流体源项生命周期

流体源项由 `SystemManager` 统一管理，避免 Picard 内迭代中耦合源项重复累加。

### 5.1 遍历流体节点

```python
_iter_fluid_volumes()
```

优先从：

```python
fluid_solver.volumes_obj
```

读取节点；若不存在，则尝试：

```python
fluid_solver.volumes
```

### 5.2 清空流体源项

```python
_clear_fluid_sources()
```

对每个流体节点清零：

```text
Q_wall = 0
Q_vol = 0
implicit_coeff = 0
```

### 5.3 施加持久源项

```python
_apply_persistent_fluid_sources()
```

逐个调用通过 `add_persistent_fluid_source()` 注册的函数。

### 5.4 准备耦合源项

```python
_prepare_fluid_sources_for_coupling()
```

等价于：

```text
清空流体源项 -> 重新施加持久源项
```

### 5.5 捕获和恢复源项

```python
_capture_fluid_sources()
_restore_fluid_sources(snapshot)
```

快照内容为每个节点的：

```text
(Q_wall, Q_vol, implicit_coeff)
```

在每轮 Picard 迭代开始时，系统会恢复到基准源项，再运行耦合器重新累加当轮界面源项。

## 6. 固体边界缓存与耦合器调度

### 6.1 刷新固体边界缓存

```python
_refresh_solid_boundary_cache(update_flux=False, current_time=None)
```

对每个固体依次调用：

```text
_update_properties()
_compute_internal_resistance()
_update_boundaries_state(current_time)
```

若 `update_flux=True`，还会调用：

```text
_compute_fluxes(current_time)
```

或逐个边界调用：

```text
boundary.compute_net_flux_for_solver()
```

该缓存刷新保证耦合器读取到最新的：

- 固体节点温度；
- 固体物性；
- 内部热阻；
- 边界表面温度；
- 边界热流。

### 6.2 运行耦合器

```python
_run_couplers(interface_relaxation=1.0, current_time=None)
```

执行顺序是：

```text
1. 刷新固体边界缓存
2. 对所有带 sync() 的耦合器调用 sync()
   - 固固
   - 间隙
   - TEC
3. 再次刷新固体边界缓存
4. 对 FluidSolidCouple 调用 execute(interface_relaxation)
5. 对其他带 execute() 的耦合器调用 execute()
```

这种顺序保证固固/间隙耦合先更新固体边界，随后流固耦合基于最新壁面状态计算对流换热。

## 7. `step()` 主时间步流程

`step()` 是 `SystemManager` 的核心接口：

```python
step(
    dt,
    inner_iter=1,
    convergence_tol=1e-3,
    reactivity_control=0.0,
    fail_on_fluid_nonconvergence=False,
    interface_relaxation=1.0,
    interface_convergence_tol=None,
)
```

参数含义：

| 参数 | 含义 |
| --- | --- |
| `dt` | 当前时间步长 |
| `inner_iter` | Picard 最大内迭代次数 |
| `convergence_tol` | 流体/固体温度收敛阈值 |
| `reactivity_control` | 外部控制反应性 |
| `fail_on_fluid_nonconvergence` | 流体不收敛时是否直接失败 |
| `interface_relaxation` | 界面松弛因子，范围 `(0, 1]` |
| `interface_convergence_tol` | 可选界面残差收敛阈值 |

### 7.1 参数校验

`step()` 会检查：

```text
inner_iter >= 1
0 < interface_relaxation <= 1
interface_convergence_tol >= 0 或 None
若启用 interface_convergence_tol，则 inner_iter >= 2
```

### 7.2 入口状态保存

时间步开始时：

```text
t_start = global_time
entry_fluid_sources = 当前流体源项快照
保存流体、固体、组件、点堆状态
创建 diagnostics
```

状态保存用于：

- Picard 内迭代回滚；
- 时间步异常时恢复入口状态；
- 保证中子学双缓冲状态不会被失败步污染。

### 7.3 时间步预处理

在 Picard 迭代前：

```text
1. 重置所有耦合器的界面松弛历史
2. 清理并重新施加持久流体源
3. 调用每个组件的 pre_step(dt, t_start)
4. 捕获 base_fluid_sources
```

`pre_step()` 通常用于组件级更新，例如：

- 功率分布；
- TEC 电热计算；
- 热管边界更新；
- 外部工况设置。

### 7.4 Picard 内迭代

每轮 `k` 执行：

```text
1. 恢复流体源项到 base_fluid_sources
2. 计算 coupling_time:
   k == 0 -> t_start
   k > 0  -> t_start + dt
3. _run_couplers(interface_relaxation, coupling_time)
4. 收集耦合器诊断
5. 推进中子学或让组件处理中子学
6. 如果 k > 0，回滚流体和固体到时间步入口状态
7. 应用待处理核功率到固体
8. 设置流体求解器时间
9. fluid_solver.step_Picard(dt)
10. 逐个 solid.step(dt)
11. 若 inner_iter > 1，检查温度和界面收敛
```

关键点：从第二轮 Picard 开始，系统会先完成耦合和中子学试探，然后回滚流体/固体求解状态，再应用新的功率源并重新推进。这样每轮迭代都从同一个时间步入口状态出发，只改变耦合条件和源项。

在流体求解前，若流体求解器提供 `set_time()`，还会执行：

```python
fluid_solver.set_time(coupling_time)
```

对于 `HydraulicNetwork`，这会继续向支持时间更新的控制体和连接传播，使 `PressurizerVolume` 和 `PumpJunction` 的时间表与当前耦合时刻一致。

### 7.5 流体求解

流体求解调用：

```python
fluid_solver.step_Picard(
    dt,
    max_iter=20 if inner_iter > 1 else 100,
)
```

若返回 `False`：

- 诊断中记录 warning；
- 若 `fail_on_fluid_nonconvergence=True`，抛出异常并回滚；
- 否则继续推进固体，但诊断会标记流体未收敛。

### 7.6 固体求解

对所有已注册固体：

```python
success = solid.step(dt)
```

若任一固体返回失败，直接抛出异常并回滚整个时间步。

### 7.7 温度收敛检查

当 `inner_iter > 1` 时，从第二轮开始检查：

```text
err_fluid = max(abs(T_f_curr - T_f_prev))
err_solid = 所有固体 max(abs(T_s_curr - T_s_prev)) 的最大值
```

温度收敛条件：

```text
err_fluid < convergence_tol
and
err_solid < convergence_tol
```

### 7.8 界面收敛检查

若传入 `interface_convergence_tol`，系统会从耦合器诊断中提取：

```text
interface_residual
```

取所有耦合器最大值作为界面残差：

```text
interface_residual = max(coupler.interface_residual)
```

界面收敛条件：

```text
interface_residual < interface_convergence_tol
```

如果启用了界面收敛检查，但没有耦合器提供 `interface_residual`，会抛出异常。

### 7.9 时间步完成

若没有异常：

```text
1. 提交中子学状态
2. global_time = t_start + dt
3. 同步固体 current_time 到 global_time
4. 刷新固体边界缓存并更新热流
5. 调用组件 post_step(dt, global_time)
6. diagnostics.status = "completed"
```

### 7.10 异常回滚

若时间步中任何步骤抛出异常：

```text
1. 回滚流体、固体、组件和点堆状态
2. 恢复入口流体源项
3. global_time 恢复为 t_start
4. 同步固体时间
5. 刷新固体边界缓存
6. diagnostics.status = "failed"
7. 记录异常类型和消息
8. 重新抛出异常
```

## 8. 中子学与功率处理

### 8.1 推进中子学

```python
_advance_neutronics_for_iteration(dt, reactivity_control, iteration_index)
```

优先让宏观组件处理中子学：

```python
comp.advance_neutronics(...)
```

若任一组件返回已处理，则 `SystemManager` 不再直接调用 `point_reactor`。

如果没有组件处理，且存在 `point_reactor`：

```text
1. 汇总所有 solid.get_reactivity_feedback()
2. point_reactor.step(dt, reactivity_control, total_feedback)
3. 返回 (fission_power, decay_power, total_power)
```

### 8.2 应用核功率

```python
_apply_pending_nuclear_power(fallback_power)
```

如果 `fallback_power` 不为 `None`，则对每个固体调用：

```python
solid.set_nuclear_power(p_fiss, p_decay, p_total)
```

这要求固体或组件自行实现功率分配和源项映射。

### 8.3 提交中子学

```python
_commit_neutronics(component_neutronics_handled)
```

优先调用组件的：

```python
comp.commit_neutronics()
```

如果没有组件处理和提交，且存在 `point_reactor`，则调用：

```python
point_reactor.commit()
```

## 9. 诊断信息

每个时间步会创建 `last_step_diagnostics`，主要字段包括：

| 字段 | 含义 |
| --- | --- |
| `status` | `running`、`completed` 或 `failed` |
| `t_start` / `t_end` | 时间步起止时间 |
| `dt` | 时间步长 |
| `inner_iter_limit` | Picard 最大迭代次数 |
| `iterations` | 实际执行迭代数 |
| `converged` | 整体收敛标志 |
| `temperature_converged` | 温度收敛标志 |
| `interface_converged` | 界面残差收敛标志 |
| `err_fluid_temperature` | 流体温度最大变化 |
| `err_solid_temperature` | 固体温度最大变化 |
| `interface_residual` | 最大界面残差 |
| `fluid_converged_by_iteration` | 每轮流体求解是否收敛 |
| `warnings` | 警告列表 |
| `couplers` | 最近一轮耦合器诊断 |
| `coupler_diagnostics_by_iteration` | 每轮耦合器诊断 |
| `component_neutronics_handled` | 中子学是否由组件处理 |
| `exception` | 失败时的异常信息 |

`FluidSolidCouple` 可提供界面残差诊断，`SystemManager` 通过：

```python
_collect_coupler_diagnostics()
```

收集这些信息。

## 10. 自适应时间步长

```python
compute_adaptive_dt(
    min_dt=1e-4,
    max_dt=0.5,
    safety_factor=0.8,
    use_last_step_diagnostics=True,
)
```

计算过程：

```text
1. 从 fluid_solver.get_max_stable_dt(max_limit=max_dt) 获取流体 CFL 限制
2. 乘以 safety_factor
3. 遍历 couplers，若提供 get_max_stable_dt()，取耦合稳定步长
4. dt_target = min(dt_fluid, dt_coupler)
5. 限制不超过 max_dt
6. 根据上一时间步诊断调整:
   - 流体未收敛或 Picard 未收敛 -> dt *= 0.5
   - 迭代次数达到上限 -> dt *= 0.8
7. 若已有上一时间步长，则增长率限制为 1.2 倍
8. 保存 `_last_dt` 并返回
```

注意：当物理目标步长低于 `min_dt` 时，函数会记录警告，但返回物理目标值，而不是强行钳到 `min_dt`。

## 11. 状态保存、回滚与断点续算

### 11.1 时间步内保存

```python
_save_system_state(include_components=False, include_point_reactor=False)
```

会调用：

```text
fluid_solver.save_state()
solid.save_state()
component.save_step_state()
point_reactor.save_step_state()
```

组件和点堆是否保存由参数控制。

### 11.2 时间步内回滚

```python
_rollback_system_state(include_components=False, include_point_reactor=False)
```

会调用：

```text
fluid_solver.load_state()
solid.load_state()
component.load_step_state(state)
point_reactor.load_step_state(state)
```

在 Picard 内迭代中，通常只回滚流体和固体；在时间步异常时，会连组件和点堆一起回滚。

### 11.3 同步固体时间

```python
_sync_solid_times_to_global()
```

将每个固体的 `current_time` 设置为 `global_time`。

### 11.4 保存全局状态

```python
save_global_state(filepath)
```

保存为压缩 `.npz` 文件，包含：

```text
System/global_time
System/last_dt
Fluid/*
Solid_{name}/*
Macro_{component.name}/*
PointReactor/*
```

具体保存内容由各对象的 `get_state_dict(prefix)` 决定。

### 11.5 加载全局状态

```python
load_global_state(filepath)
```

加载后依次恢复：

```text
1. global_time 和 _last_dt
2. fluid_solver.load_state_dict()
3. solid.load_state_dict()
4. component.load_state_dict()
5. point_reactor.load_state_dict()
6. 同步固体时间
7. 刷新固体边界缓存
8. 准备流体源项
9. 运行耦合器
```

加载过程要求各子系统的拓扑和状态字典兼容，否则由子系统抛出错误。

当前实现加载完成后会同步固体时间，但不会立即调用：

```python
fluid_solver.set_time(global_time)
```

下一次 `step()` 会在流体求解前自动同步。如果调用方需要在恢复后、下一次推进前读取时间表泵压升或稳压器目标，应显式执行：

```python
manager.fluid_solver.set_time(manager.global_time)
```

## 12. 典型使用流程

### 12.1 创建并注册对象

```python
manager = SystemManager(fluid_network=hydraulic_network, start_time=0.0)

manager.add_solid_component(fuel_solid)
manager.add_solid_component(wall_solid)
manager.add_coupler(fluid_solid_coupler)
manager.add_point_reactor(point_reactor)
```

若使用宏观组件：

```python
manager.add_component(reactor_core_component)
```

### 12.2 初始化

```python
manager.initialize_system(dt_init=0.1, tol=1e-5, max_iter=500)
```

### 12.3 推进一个时间步

```python
manager.step(
    dt=0.01,
    inner_iter=5,
    convergence_tol=1e-3,
    reactivity_control=0.0,
    interface_relaxation=0.7,
    interface_convergence_tol=1e-2,
)
```

### 12.4 自适应步长运行

```python
dt = manager.compute_adaptive_dt(
    min_dt=1e-4,
    max_dt=0.5,
    safety_factor=0.8,
)

manager.step(dt, inner_iter=5, convergence_tol=1e-3)
```

### 12.5 断点续算

```python
manager.save_global_state("restart_state.npz")
manager.load_global_state("restart_state.npz")
```

## 13. 重要实现细节与注意事项

### 13.1 耦合源项会累加，必须管理生命周期

`FluidSolidCouple.execute()` 会向流体节点累加 `Q_wall` 和 `implicit_coeff`。因此 `SystemManager` 在每轮 Picard 内迭代开始时会恢复基准流体源项，避免重复累加。

### 13.2 Picard 内迭代从同一入口状态重新推进

从第二轮迭代开始，`SystemManager` 会回滚流体和固体状态，再用新的耦合源项和功率源重新推进。这保证每轮迭代比较的是同一时间步入口条件下的不同耦合试探结果。

### 13.3 组件可接管中子学

如果组件实现并返回 `advance_neutronics()` 已处理，`SystemManager` 不会直接调用 `point_reactor.step()`。这允许复杂堆芯组件自行管理功率分配和反馈。

### 13.4 `interface_convergence_tol` 依赖耦合器诊断

只有提供 `get_coupling_diagnostics()` 且返回 `interface_residual` 的耦合器才能参与界面收敛判断。当前 `FluidSolidCouple` 支持该诊断。

### 13.5 异常时会恢复入口流体源项

如果时间步失败，系统会恢复到进入 `step()` 前的流体源项，而不是 Picard 基准源项。这保证失败步不污染外部源项状态。

### 13.6 全局状态文件依赖各子系统实现

`save_global_state()` 和 `load_global_state()` 是统一调度接口，实际字段由流体、固体、组件和点堆各自的 `get_state_dict()` / `load_state_dict()` 决定。

## 14. 模块优点

`SystemManager` 当前实现具有以下特点：

- 把流体、固体、中子学、耦合器和宏观组件统一到一个时间步流程；
- 使用 Picard 内迭代处理强耦合热工问题；
- 支持界面松弛和界面残差诊断；
- 明确管理流体源项生命周期，避免耦合源项重复累加；
- 支持时间步失败自动回滚；
- 支持组件级扩展点，不把所有物理逻辑硬编码在调度器中；
- 支持自适应时间步和增长率限制；
- 支持压缩 `.npz` 全局断点续算。

## 15. 推荐阅读顺序

维护或扩展 `SystemManager.py` 时，建议按以下顺序阅读：

1. `SystemManager.step()`
   - 先理解完整全局时间步和 Picard 内迭代流程；
2. `_run_couplers()` 和 `_refresh_solid_boundary_cache()`
   - 理解耦合器调用顺序和边界缓存生命周期；
3. `_advance_neutronics_for_iteration()` 和 `_commit_neutronics()`
   - 理解点堆和组件中子学的优先级；
4. `_save_system_state()` / `_rollback_system_state()`
   - 理解异常恢复和内迭代回滚；
5. `compute_adaptive_dt()`
   - 理解时间步控制逻辑；
6. `save_global_state()` / `load_global_state()`
   - 理解断点续算边界。

## 16. 源码校验信息

最后按源码校验日期：2026-05-31。

涉及文件：

- `SystemManager.py`
- `Hydrodynamics/HydraulicNetwork.py`
- `Hydrodynamics/Components.py`
- `HeatConduction/HeatConduction.py`
- `HeatConduction/Boundary.py`
- `Couplers.py`
- `Neutronics/PointReactor.py`
