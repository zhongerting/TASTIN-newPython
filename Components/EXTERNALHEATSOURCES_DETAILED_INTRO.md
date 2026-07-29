# ExternalHeatSources Detailed Intro

本文档总结 `Components/ExternalHeatSources` 目录中的外部热流模块，面向程序使用和二次开发。

## 模块定位

`ExternalHeatSources` 用于把轨道环境、太阳辐照、地球反照率、地球红外、查表热流和用户自定义热流统一封装成边界热流源。

模块的核心约定是：

- 热源模型返回热流密度，单位为 W/m2。
- 正值表示热量流入固体边界。
- `ExternalHeatFluxBC` 负责把热流密度乘以边界面积，转换为求解器使用的热流功率，单位为 W。
- 热流源形状 `shape` 必须与目标 `BoundaryRegion.shape` 一致，或能安全广播到该形状。

## 文件概览

| 文件 | 内容 |
| --- | --- |
| `__init__.py` | 外热流模型、边界条件封装、组合热源和预留积分模型 |
| `embedded_flux_tables.py` | 内嵌 Fortran 轨道热流查表数据和查表库 |
| `README.md` | 现有模块定位、推荐路线和限制说明 |
| `PARAMETERS.md` | Fortran 参数映射、轨道角度和表面法向角说明 |

## 基础接口

### `BaseExternalHeatSource`

所有外热流模型的抽象基类。

主要接口：

- `get_heat_flux(time)`：返回当前时间的热流密度数组。
- `update_params(**kwargs)`：运行中更新轨道、姿态、热流或缩放参数。
- `_broadcast_flux(flux)`：把标量或数组转换成与边界一致的数组形状。

使用要求：

- 子类必须实现 `get_heat_flux()` 和 `update_params()`。
- 返回数组单位为 W/m2。
- 返回数组形状应为初始化时传入的 `shape`。

### `ExternalHeatFluxBC`

外热流边界条件封装类，用于挂到 `BoundaryRegion.conditions`。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `heat_source` | `BaseExternalHeatSource` 子类实例 |
| `area_array` | 边界节点面积数组，单位 m2 |

关键行为：

- `update_state(T_node)` 从 `BoundaryRegion` 接收当前节点状态。
- 内部读取 `current_time`，调用 `heat_source.get_heat_flux(current_time)`。
- 将 W/m2 转换为 W：

```text
q_flux = heat_flux_density * area_array
```

- `compute_flux_from_node()` 返回当前热流功率数组。
- `update_params()` 会把参数更新转发给内部热源。

典型挂接方式：

```python
from Components.ExternalHeatSources import OrbitalHeatSource, ExternalHeatFluxBC

boundary = hp.hp.boundaries["outer_con"]
source = OrbitalHeatSource(
    shape=boundary.shape,
    surface_normal_angles=(1.57079632679, 0.0),
)

bc = ExternalHeatFluxBC(
    heat_source=source,
    area_array=boundary.area,
)
boundary.conditions.append(bc)
```

## 已实现热源

### `OrbitalHeatSource`

解析式轨道太阳辐射热源。

用途：

- 根据轨道周期、轨道高度、轨道倾角和表面法向角计算太阳入射因子。
- 适合做单根热管或简化辐射器受照分析。

重要参数：

| 参数 | 说明 |
| --- | --- |
| `solar_constant` | 太阳常数，默认 1361 W/m2 |
| `orbit_height` | 轨道高度，默认 800 km |
| `orbit_period` | 轨道周期，默认 7644 s |
| `orbit_inclination` | 轨道倾角，单位 deg |
| `surface_normal_angles` | 表面法向角 `(w0, l0)`，单位 rad |

热流计算：

```text
q = solar_constant * X
```

其中 `X` 为太阳方向与表面法向的入射余弦因子，小于 0 的部分被截断为 0。

### `OrbitalTableHeatSource`

Fortran 风格的轨道热流查表源。

用途：

- 复现旧 Fortran 输入卡中辐射器分区的热流表。
- 不再依赖运行时读取外部 `RadiatorInput.txt`，数据已内嵌在 Python 文件中。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `shape` | 边界形状 |
| `table_ids` | 表编号，可以是标量、列表或数组 |
| `table_library` | 查表库，默认 `FORTRAN_ORBITAL_HEAT_TABLE_LIBRARY` |
| `scale_factor` | 热流缩放系数 |
| `offset` | 热流偏置 |
| `periodic` | 是否按表周期循环采样 |

使用方式：

```python
from Components.ExternalHeatSources import OrbitalTableHeatSource

source = OrbitalTableHeatSource(
    shape=boundary.shape,
    table_ids=1,
    scale_factor=1.0,
)
```

如果 `table_ids` 是数组，则每个边界节点可使用不同热流表。

### `AlbedoHeatSource`

简化地球反照率热源。

当前模型：

```text
q_albedo = albedo_factor * solar_constant
```

默认 `albedo_factor = 0.3`，`solar_constant = 1361 W/m2`。

注意：这是工程简化模型，不是严格地球反照率积分模型。

### `EarthIRHeatSource`

简化地球红外热源。

默认热流：

```text
earth_ir_flux = 237 W/m2
```

`emissivity` 参数被保存，但当前 `get_heat_flux()` 直接返回 `q_ir`，并未额外乘以发射率。

### `TimeVaryingHeatSource`

通用时间函数热源。

用途：

- 包装用户自定义函数。
- 实现脉冲热流、控制信号热流、测试曲线等。

示例：

```python
source = TimeVaryingHeatSource(
    shape=boundary.shape,
    time_func=lambda t: 200.0 if t < 100.0 else 0.0,
    amplitude=1.0,
)
```

### `FixedPowerHeatSource`

恒定热流密度源。

示例：

```python
source = FixedPowerHeatSource(
    shape=boundary.shape,
    q_flux=50.0,
)
```

这里 `q_flux` 单位仍是 W/m2，不是总功率。

### `CompositeHeatSource`

组合热源，用于把多个热流源叠加为一个边界源。

示例：

```python
from Components.ExternalHeatSources import (
    CompositeHeatSource,
    OrbitalHeatSource,
    AlbedoHeatSource,
    EarthIRHeatSource,
)

source = CompositeHeatSource(shape=boundary.shape)
source.add_source(OrbitalHeatSource(shape=boundary.shape))
source.add_source(AlbedoHeatSource(shape=boundary.shape))
source.add_source(EarthIRHeatSource(shape=boundary.shape))
```

注意：

- 所有子热源的 `shape` 必须一致。
- `update_params(**kwargs)` 会广播给所有子热源。

## 预留但尚未实现的严格积分模型

### `OrbitalIntegralHeatSource`

严格轨道积分热流模型的预留基类。

### `IntegralAlbedoHeatSource`

严格地球反照率积分模型的预留类，当前调用 `get_heat_flux()` 会抛出 `NotImplementedError`。

### `IntegralEarthIRHeatSource`

严格地球红外积分模型的预留类，当前调用 `get_heat_flux()` 会抛出 `NotImplementedError`。

这些类用于保留未来扩展接口，不应直接用于生产算例。

## 查表数据结构

### `EmbeddedFluxTable`

不可变数据类，字段包括：

| 字段 | 说明 |
| --- | --- |
| `table_id` | 表编号 |
| `name` | 表名称 |
| `time` | 时间采样点 |
| `values` | 热流密度值 |
| `periodic` | 是否周期循环 |

### `EmbeddedFluxTableLibrary`

查表库容器。

主要方法：

- `get_table(table_id)`：返回指定表。
- `available_ids()`：返回可用表编号。
- `has_table(table_id)`：检查表是否存在。

当前内嵌库：

- `FORTRAN_ORBITAL_HEAT_TABLES`
- `FORTRAN_ORBITAL_HEAT_TABLE_LIBRARY`

其中普通辐射器分区表从 `1` 开始编号，`1001` 为顶部参考面热流表。

## Fortran 参数映射

`PARAMETERS.md` 中区分了两类角度：

| 角度 | Python 参数 | 物理含义 |
| --- | --- | --- |
| 热管倾斜角 | `hp_tilt_angle` | 热管轴线与水平面的夹角，主要用于重力压降 |
| 表面法向极角 | `surface_normal_angles[0]` 或 `w0` | 表面法向与 z 轴夹角 |
| 表面法向方位角 | `surface_normal_angles[1]` 或 `l0` | 表面法向在 x-y 平面投影方位 |

当前外热流模块真正使用的是 `(w0, l0)` 表面法向角；`hp_tilt_angle` 属于热管水力或安装姿态参数，不应与受照面法向角混用。

## 与 `HPwithFin` 的关系

`HPwithFin` 支持两种外热流接入方式：

1. 把外热流作为冷凝段裸壁边界热流挂到 `outer_con`。
2. 使用 `set_fin_external_heat_source()`，让外热流直接进入降维翅片准稳态方程。

推荐在需要区分“热管壁吸热”和“翅片直接受照”时使用第二种方式，并用 `configure_external_heat_accounting()` 做后处理核算。

## 2026-07-15 V14 热管外热吸收口径

V14 的 N18 外热表对应 18 个代表热管区，分别挂到六个 `RingHP` 的三个流体节点；表格本身仍只返回 `W/m2`，不包含面积和热管数量倍率。V14 启用 `distributed_fin_absorption` 时：

- 管壁只通过外侧受照面积接收外热，默认 `wall_illumination_factor=0.5`；
- 翅片只通过单侧投影/受照面积接收外热，内侧翅片面不接收外热；
- 管壁有效吸收面积为 `A_outer_wall * illumination * 0.992 * epsilon_hp`；
- 翅片有效吸收面积为 `A_fin_one_side * illumination * 0.992 * epsilon_fin`；
- 同一个代表热管的源项只挂一次：管壁进入 `ExternalHeatFluxBC`，翅片进入 `HPwithFin` 准稳态翅片方程；
- `hp_multipliers` 和 V14 的 `symmetric_ring_multiplier` 只在 RingHP/整机汇总层使用，不在单根热管吸收面积中重复使用。

`0.992` 仅表示外热吸收效率，`alpha=epsilon` 仅用于入射外热；它不修改热管或翅片向空间辐射的有效发射率。外侧和内侧辐射角系数属于 `HPwithFin` 的辐射排热模型，不能拿来替代外热受照面积。

## 使用建议

优先路线：

1. 工程复现实验或旧程序对齐：使用 `OrbitalTableHeatSource`。
2. 简化太阳受照分析：使用 `OrbitalHeatSource`。
3. 组合环境热流：使用 `CompositeHeatSource` 叠加太阳、反照率和地球红外。
4. 临时工况或测试：使用 `FixedPowerHeatSource` 或 `TimeVaryingHeatSource`。

单位检查：

- 热源返回 W/m2。
- 边界面积为 m2。
- `ExternalHeatFluxBC.q_flux` 为 W。
- 若使用面积缩放，缩放应施加到 `area_array` 或热源 `scale_factor`，不要重复缩放。

## 当前限制

- `IntegralAlbedoHeatSource` 和 `IntegralEarthIRHeatSource` 未实现。
- `AlbedoHeatSource` 与 `EarthIRHeatSource` 是简化常值模型。
- 查表法输出总热流，不分离太阳、反照率和地球红外分量。
- `OrbitalHeatSource` 的几何模型适合简化分析，复杂遮挡、姿态变化和真实可见因子需要后续扩展。



### `OrbitalMatrixHeatSource`

Matrix lookup external heat source for the fixed case `w0=8.12deg, i_s=58.5deg`. Only the pre-summed `sum` heat flux is packaged, so the runtime does not add solar, albedo, and Earth IR components separately.

Packaged matrices are embedded in `Components/ExternalHeatSources/w0_8p12_sum_data.py`:

| `matrix_key` | Columns |
| --- | --- |
| `is58p5_w0_8p12_N6_sum` | 6 |
| `is58p5_w0_8p12_N18_sum` | 18 |
| `is58p5_w0_8p12_N78_sum` | 78 |
| `is58p5_w0_8p12_N200_sum` | 200 |

`OrbitalMatrixHeatSource` returns `W/m2`. `ExternalHeatFluxBC` remains responsible for multiplying by the boundary or fin illuminated area and converting the load to `W`. The matrix column count must match `np.prod(shape)`; this class performs time interpolation and periodic wrapping only, not spatial resampling.

```python
source = OrbitalMatrixHeatSource(
    shape=boundary.shape,
    matrix_key="is58p5_w0_8p12_N6_sum",
)
```

### CSV circumferential tables

load_csv_flux_table_library(path, orbit_period_s=...) reads CSV files formatted as a sample index followed by theta columns and exposes each circumferential column as a scalar OrbitalTableHeatSource table. The file must contain at least two rows, finite values, and a strictly increasing sample column. When orbit_period_s is provided, the sample axis is scaled to that physical period. The packaged legacy matrix keeps its 6552 s axis; the V14 startup CSV path explicitly rescales all 360 samples to `5668.14 s`. Values remain in W/m2. V14 uses this path for 18 representative heat-pipe zones and V15 for 78 radiator tubes.

### External-heat loading and no-double-counting contract

The runtime path is:

```text
sum table [W/m2]
  -> OrbitalTableHeatSource / OrbitalMatrixHeatSource
  -> ExternalHeatFluxBC multiplies each axial segment area once
  -> BoundaryRegion adds the resulting segment powers [W] once
```

The `*_sum` tables are complete external-heat inputs. Do not add separate solar, albedo, or Earth-IR sources on top of them. `RingHP._build_external_heat_source()` treats matrix/table modes as exclusive and returns before the analytic component path.

Repeating one circumferential table value over several axial nodes does not multiply total power by the node count. Each node receives the same heat-flux density but uses only its own segment area, so the boundary total remains `q_density * sum(segment_area)`.

`ExternalHeatFluxBC.update_state()` overwrites `q_flux`; it does not increment it. `BoundaryRegion.compute_net_flux_for_solver()` clears and rebuilds `Q_sum_flux` on every call, so Picard iterations and ODE trial evaluations do not accumulate external heat across evaluations. `HPwithFin.configure_external_heat_accounting()` is diagnostic only and never adds another boundary source.

For representative V14 heat pipes, `hp_multipliers` scale the heat exchanged with the collector-ring fluid and scaled diagnostics; they do not create additional external-heat boundary conditions. V15 maps one table column to one physical radiator tube. The topology tests require exactly one `ExternalHeatFluxBC` per V14 representative heat pipe and per V15 radiator tube.

V15 的 RadiatorPipeWithFin 额外支持可选的 distributed_fin 模式。该模式仍只保留一个管壁 ExternalHeatFluxBC，但该边界只使用管壁外侧受照面积；翅片部分通过同一个热流源进入准稳态翅片方程。管壁与翅片的有效吸热面积分别为 `A_outer * illumination * 0.992 * epsilon_tube` 和 `fin_strip_width * fin_height * illumination * 0.992 * epsilon_fin`，即采用 `alpha=epsilon`。它们不乘仅用于双面辐射几何修正的 fin_area_scale、fin_view_factor 或内外面辐射角系数。总吸热诊断必须调用 get_external_heat_absorption_distribution()，不能只累计管壁边界的 q_flux。

### 2026-07-20 V14 10 kW two-ring N18 case

`Full_Loop_Cases_10kW/V14_210kW_reactivity_control_external_heat` explicitly loads the N18 CSV with a `5668.144369 s` period. Upper and lower rings reuse the same columns: sector `n` maps its three representative heat pipes to columns `3n`, `3n+1`, and `3n+2`. `OrbitalTableHeatSource.time_origin_s` is subtracted from absolute system time, so the loaded restart remains at its saved `System/global_time` while that instant is external-heat phase zero. The default remains `time_origin_s=0` and external heat remains disabled in existing 10 kW V14 cases.

The V14 shield startup path uses `5668.14 s` for both histories. Its N6 values are the exact center-point subset of N18 columns `0,3,6,9,12,15`. Shield-present steps apply `qsss=(alpha/epsilon)q_incident=0.992*N6` to the six side panels only; they do not multiply outer emissivity a second time and do not apply N18 to the heat pipes.
