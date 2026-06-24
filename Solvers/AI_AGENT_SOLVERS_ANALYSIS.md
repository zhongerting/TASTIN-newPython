# Solvers - Codex 必读入口

> 适用范围：`Solvers/` 目录。  
> 最后按源码校验日期：2026-05-31。

## 1. Codex 使用说明

进入 `Solvers/` 后，默认先读本文件，再根据任务读取对应模块文档。不要默认重复扫描全部源码，也不要把 `SOLVERS_ANALYSIS.md` 当作首次进入目录时的必读文件。

推荐路径：

| 任务 | 先读文档 | 必要时回看源码 |
| --- | --- | --- |
| 压力、流量、泵、稳压器、动态阻力 | `HYDRODYNAMICS_DETAILED_INTRO.md` | `Hydrodynamics/*.py` |
| 固体网格、边界、辐射、热流单位 | `HEATCONDUCTION_DETAILED_INTRO.md` | `HeatConduction/*.py` |
| 流固、固固、间隙、TEC 耦合 | `COUPLERS_DETAILED_INTRO.md` | `Couplers.py` |
| 点堆、衰变热、提交与恢复 | `NEUTRONICS_DETAILED_INTRO.md` | `Neutronics/PointReactor.py` |
| 全局 Picard、源项、回滚、快照 | `SYSTEMMANAGER_DETAILED_INTRO.md` | `SystemManager.py` |
| 需要跨模块公式和实现细节 | `SOLVERS_ANALYSIS.md` | 对应源码 |

文档与源码冲突时，以源码为准，并在当前任务中同步修正文档。

## 2. 核心文件映射

`Solvers/` 有 9 个核心 Python 文件：

| 文件 | 职责 | 公开类 |
| --- | --- | --- |
| `Hydrodynamics/Components.py` | 流体控制体、连接、泵、通道 | `FluidVolume`, `IncompressibleFluidVolume`, `PressurizerVolume`, `FixedPressurePumpVolume`, `FlowJunction`, `PumpJunction`, `MacroFlowJunction`, `FluidChannel`, `IncompressibleFluidChannel`, `NonUniformIncompressibleFluidChannel` |
| `Hydrodynamics/BoundaryVolume.py` | 压力、温度和流量边界 | `BoundaryVolume`, `IncompressibleBoundaryVolume`, `InletJunction` |
| `Hydrodynamics/HydraulicNetwork.py` | 水力网络向量化与稀疏求解 | `HydraulicNetwork` |
| `HeatConduction/Mesh.py` | 1D/2D 有限体积网格 | `GeometricData`, `Mesh1D`, `GeometricData2D`, `Mesh2D` |
| `HeatConduction/Boundary.py` | 热边界与边界热阻网络 | `BaseBoundaryCondition`, `ResistanceBC`, `FluxBC`, `DynamicRadiationResistanceBC`, `BoundaryRegion` |
| `HeatConduction/HeatConduction.py` | 1D/2D 固体瞬态导热 | `BaseHeatConduction`, `HeatConduction1D`, `HeatConduction2D` |
| `Couplers.py` | 流固、固固、间隙和 TEC 耦合 | `SolidSolidCouple1D`, `SolidSolidCouple2D`, `FluidSolidCouple`, `GapCouple2D`, `ActiveGapCouple2D`, `TECCouple2D` |
| `Neutronics/PointReactor.py` | 点堆动力学与衰变热 | `PointReactor` |
| `SystemManager.py` | 全局瞬态调度 | `SystemManager` |

## 3. 跨模块数据流

主调用方向：

```text
SystemManager.step()
  -> Couplers
       -> HydraulicNetwork / FluidChannel
       -> HeatConduction / BoundaryRegion
  -> PointReactor 或由宏观组件接管中子学
  -> HydraulicNetwork.step_Picard()
  -> HeatConduction.step()
  -> 收敛检查、提交或回滚
```

流固换热：

```text
HydraulicNetwork.FluidChannel
  -> FluidSolidCouple: Re -> Nu -> h -> lambda = h * A
  -> HeatConduction.BoundaryRegion: T_ext = T_fluid, R_ext = 1 / lambda
  -> FluidChannel: Q_fluid = lambda * T_wall - lambda * T_fluid
```

固固与间隙：

```text
BoundaryRegion
  <-> SolidSolidCouple / GapCouple2D / ActiveGapCouple2D / TECCouple2D
  <-> BoundaryRegion
```

## 4. 当前实现的重要约束

### 4.1 水力网络

- `HydraulicNetwork.set_time(time)` 会保存网络时间，并向所有提供 `set_time()` 的控制体和连接传播。
- `SystemManager.step()` 在调用流体求解前执行 `fluid_solver.set_time(coupling_time)`。首轮 Picard 使用步首时间，后续轮次使用步末时间。
- `PressurizerVolume` 是封闭不可压缩网络的被动绝对压力参考点，不是定压质量边界。网络求解压力差后，通过统一压力平移锚定绝对压力。
- 每个 `HydraulicNetwork` 最多支持一个被动压力参考点；同一节点不能同时是定压边界和被动参考点。
- `PumpJunction.compute_pump_head(time)` 支持固定压升和时间表插值，正压升驱动 `from_vol -> to_vol`。
- `FixedPressurePumpVolume` 只用于不可压缩通道的代数压力分布递推，不替代网络动量方程中的 `PumpJunction`。
- 动态局部阻力由 `dynamic_loss_params` 启用。网络快速路径在每轮水力 Picard 中按当前迭代流量刷新 `effective_K_loss_vec`。
- 动态阻力启用时，负的基础 `k_loss` 会按 `0` 处理；不要把负值本身理解为动态模型开关。
- 网络拓扑、定压边界集合和 CSR 结构在初始化后视为固定。改变节点、连接或定压边界时应重建网络。

### 4.2 固体导热与边界

- `BoundaryRegion` 和 `FluidChannel` 中的热源、热流接口按离散面或离散节点总功率处理，单位是 W，不是 W/m2。
- `DynamicRadiationResistanceBC.update_state()` 允许空数组。空边界会直接返回，不应触发 `np.min()` 异常。
- 动态辐射边界依赖边界缓存刷新；修改边界更新顺序时回看 `BoundaryRegion.update_internal_state()` 和 `SystemManager._run_couplers()`。

### 4.3 耦合与全局推进

- `FluidSolidCouple.execute()` 会累加流体半隐式源项；每轮 Picard 必须先恢复基准流体源项。
- `SystemManager._run_couplers()` 的顺序是：刷新固体边界缓存、运行 `sync()` 耦合器、再次刷新缓存、运行 `execute()` 耦合器。
- 从第二轮全局 Picard 起，流体和固体会回滚到时间步入口状态，再用新的耦合条件重新推进。
- `PointReactor.step()` 只写试探状态；全局时间步完成后才由 `commit()` 固化。
- `SystemManager.load_global_state()` 当前不会立即调用 `fluid_solver.set_time(global_time)`。恢复后若要在下一次 `step()` 前读取时间表泵压升或稳压器目标，先显式调用 `fluid_solver.set_time(manager.global_time)`。

## 5. 典型调用片段

时间相关泵和稳压器：

```python
pressurizer = PressurizerVolume(
    name="pressurizer",
    volume=1.0e-3,
    length=0.1,
    flow_area=1.0e-3,
    hydraulic_diam=0.03,
    material=fluid_material,
)
pressurizer.set_pressure_table([0.0, 10.0], [1.0e5, 1.2e5])

pump = PumpJunction(
    name="main_pump",
    from_vol=vol_a,
    to_vol=vol_b,
    delta_p=2.0e4,
)
pump.set_pressure_table([0.0, 5.0], [0.0, 2.0e4])

network = HydraulicNetwork(volumes, junctions)
network.set_time(2.5)
```

全局调度：

```python
manager = SystemManager(fluid_network=network, start_time=0.0)
manager.add_solid_component(solid)
manager.add_coupler(coupler)
manager.step(dt=0.1, inner_iter=5, convergence_tol=1.0e-3)
```

## 6. 其他文档定位

| 文档 | 定位 |
| --- | --- |
| `HYDRODYNAMICS_DETAILED_INTRO.md` | 水力模块详细说明 |
| `HEATCONDUCTION_DETAILED_INTRO.md` | 导热、网格和热边界详细说明 |
| `COUPLERS_DETAILED_INTRO.md` | 耦合器详细说明 |
| `NEUTRONICS_DETAILED_INTRO.md` | 点堆详细说明 |
| `SYSTEMMANAGER_DETAILED_INTRO.md` | 全局调度详细说明 |
| `SOLVERS_ANALYSIS.md` | 跨模块深度参考，不是默认必读入口 |

## 7. 修改后维护清单

- 修改核心类、公开方法、执行顺序、状态持久化格式或物理模型后，必须同步更新对应模块文档。
- 修改跨模块调用关系时，必须同步更新本文件的数据流和阅读路径。
- Codex 后续处理局部任务时，默认先读本文件，再读对应模块文档；只有文档不足或代码已变更时才重新扫描源码。
- 文档与源码冲突时，以源码为准，并在当前任务中同步修正文档。

## 8. 本次校验涉及文件

- `Hydrodynamics/Components.py`
- `Hydrodynamics/BoundaryVolume.py`
- `Hydrodynamics/HydraulicNetwork.py`
- `HeatConduction/Mesh.py`
- `HeatConduction/Boundary.py`
- `HeatConduction/HeatConduction.py`
- `Couplers.py`
- `Neutronics/PointReactor.py`
- `SystemManager.py`

## 9. 2026-06-15 local implicit fluid-solid coupling

`FluidSolidCouple` has an opt-in `coupling_time_scheme="local_implicit"` mode. The default remains `"current"`, so existing builders, tests, and restart workflows keep their previous Robin-boundary plus fluid semi-implicit-source behavior.

When local implicit mode is enabled, `SystemManager.step(dt, ...)` passes `dt` into `FluidSolidCouple.execute(...)`. The coupler solves the local two-capacitance heat exchange analytically for each interface node, applies `q_to_solid = -q_to_fluid` through a solid-side `FluxBC`, and adds `q_to_fluid` to the fluid channel with zero implicit coefficient. This removes the explicit fluid-solid exchange time constant from `FluidSolidCouple.get_max_stable_dt()`, but does not override fluid hydraulic stability, `max_dt`, or convergence-based time-step reductions.

Do not enable local implicit mode on a `FluidSolidCouple` without `solid_node_capacitance`; the coupler intentionally raises an error instead of silently falling back.

## 10. 2026-06-24 generic HeatConduction `implicit_euler` path

`HeatConduction1D` / `HeatConduction2D` now support a per-solid ODE mode `implicit_euler` via `set_ode_method("implicit_euler")`.

`implicit_euler` performs a backward-Euler sparse algebraic solve once per global `SystemManager.step()` time step. Existing `solve_ivp` methods remain supported (`RK45`, `RK23`, `DOP853`, `Radau`, `BDF`, `LSODA`) and continue to be the default behavior for generic solids when not overridden.

Boundary mapping for the sparse path:
- resistance/convection/dynamic-radiation boundary terms are assembled into the implicit linear system;
- pure flux (`FluxBC`) terms stay on the explicit RHS.

On failure of the implicit solve, `HeatConduction.step()` emits a warning and falls back to a `solve_ivp` step.

No `SystemManager` scheduling changes are required: it still calls `solid.step(dt)` exactly as before.

Initial V13 no-TEC smoke validation used `testModule/run_v13_caseA_closed_loop.py --disable-tec-coupled --solid-ode-method implicit_euler`. With `max_dt=0.1 s`, the 1 s result stayed close to the RK45 baseline: core outlet `813.245 K` vs `813.355 K`, core connector delta-T `85.707 K` vs `85.812 K`, and coolant enthalpy rise `97184.680 W` vs `97303.408 W`. With `max_dt=0.5 s`, backward-Euler numerical dissipation was visibly larger, so V13 accuracy comparisons should use the same global `max_dt` and should not assume large-step `implicit_euler` matches RK45 within `0.5 K`.

## 11. 2026-06-24 implicit matrix assembly optimization

This round keeps the same generic `implicit_euler` numerical contract and does not introduce `SolidAssembly`.

Implementation changes:
- `HeatConduction1D._assemble_theta_implicit_matrix()` and `HeatConduction2D._assemble_theta_implicit_matrix()` now assemble internal conduction edges with NumPy vector operations instead of Python double loops.
- The sparse CSC matrix pattern is cached per solid shape. Each time step updates only the cached numeric `data` arrays and boundary/source-dependent diagonal values.
- Sparse storage remains CSC for direct `scipy.sparse.linalg.spsolve`; there is no LU/factorization cache in this round, because material and radiation coefficients can change every step.
- `testModule/benchmark_heatconduction_implicit.py` is the lightweight reproducible benchmark for this path.

Recorded stepwise micro-benchmark on the compact validation case (`1D n=96, 80 steps`; `2D 12x96, 20 steps`):

| stage | 1D per step | 2D per step | note |
| --- | ---: | ---: | --- |
| before vectorization | `0.000459 s` | `0.006884 s` | Python loop assembly |
| vectorized assembly | `0.000353 s` | `0.003043 s` | main 2D gain |
| CSC pattern/data cache | `0.000262 s` | `0.002947 s` | main 1D gain, 2D solve now dominates |

Fresh post-optimization checks:
- `python -m py_compile Solvers\HeatConduction\HeatConduction.py Solvers\SystemManager.py Components\basicComponents\HeatPipe2D.py testModule\test_heatconduction_implicit.py testModule\benchmark_heatconduction_implicit.py` passed.
- `python -m unittest testModule.test_heatconduction_implicit testModule.test_system_manager_lifecycle` ran 25 tests and passed.
- `testModule/benchmark_heatconduction_implicit.py --steps-1d 40 --steps-2d 20 --n1d 4000 --n2d-x 90 --n2d-y 90` measured `0.001869 s/step` for 1D 4000 nodes and `0.026584 s/step` for 2D 8100 nodes on this machine.
- Supplemental same-script solve-method comparison measured BDF at about `1.971 s/step` for 1D 4000 nodes and `0.0574 s/step` for 2D 8100 nodes; RK45 measured about `0.0101 s/step` for 1D and `0.00282 s/step` for 2D, but RK45 is a different explicit method and should not be treated as an accuracy-equivalent replacement for the fixed-step implicit path.
- A direct `HeatPipe2D` `implicit_euler` smoke step with real `HPwithFin`/`SodiumHP` materials passed. The script-level `SystemManager.step()` path returns `None` on success for that dummy network, so smoke checks should not require a strict `True` return value.
- V13 no-TEC `implicit_euler`, `max_dt=0.1 s`, `duration=1 s` completed in `4.961685 s` wall time. Compared with the pre-optimization `implicit_euler` result, key deltas were numerical roundoff level: core outlet `-1.79e-08 K`, coolant enthalpy rise `-1.54e-05 W`, radiator total heat `-8.84e-07 W`.
