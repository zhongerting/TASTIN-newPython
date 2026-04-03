# TASTIN ExternalHeatSources 说明

## 1. 模块定位

`Components/ExternalHeatSources` 用来描述轨道环境下的外热流，并将这些外热流通过 `ExternalHeatFluxBC` 挂接到现有 `BoundaryRegion` 边界体系中。

本模块当前主要服务两类需求：

1. 工程计算需求  
   使用查表法，直接按时间读取轨道热流表，快速稳定。
2. 重构研究需求  
   保留解析/积分法接口，逐步把原有理论模型补齐。

---

## 2. 当前实现状态

这是目前最重要的结论。

### 2.1 已实现

1. `ExternalHeatFluxBC`  
   已经可以与 `BoundaryRegion` 正常耦合，并支持时间相关热流自动更新。

2. `OrbitalHeatSource`  
   已实现太阳直射热流的解析计算。  
   该实现参考了 Fortran 版本中的太阳热流计算思路，包含轨道周期、轨道高度、面元朝向和地影影响。

3. `OrbitalTableHeatSource`  
   已实现 Fortran 风格的查表法。  
   当前做法是：
   - 直接在 Python 程序内部嵌入热流表
   - 运行时按表号选择
   - 按当前时间做周期折返与线性插值

4. `AlbedoHeatSource`  
   已实现简化版地球反照热流模型。  
   当前不是严格几何积分，而是一个简化常值模型。

5. `EarthIRHeatSource`  
   已实现简化版地球红外热流模型。  
   当前也不是严格几何积分，而是一个简化常值模型。

6. `CompositeHeatSource`  
   已实现多热源叠加。

### 2.2 已预留但尚未实现

1. `IntegralAlbedoHeatSource`
2. `IntegralEarthIRHeatSource`

这两个类代表“严格解析/积分法”路线的预留接口。  
当前如果直接调用，会抛出 `NotImplementedError`。

也就是说，目前模块状态应理解为：

- 太阳直射：解析法已实现
- 地球反照：当前为简化模型，严格积分法未实现
- 地球红外：当前为简化模型，严格积分法未实现
- 查表法：已实现并可直接用于工程计算

---

## 3. 两条主要路线

### 3.1 查表法

查表法对应 `OrbitalTableHeatSource`。

它的特点是：

1. 不在运行时现场求解轨道几何公式
2. 不区分太阳、反照、红外各自的中间物理量
3. 直接把“某个位置在某个时刻的总轨道热流”作为已知输入

当前实现流程：

1. 通过表号选择内置热流表
2. 根据当前时间做周期折返
3. 通过 `np.interp` 进行线性插值
4. 返回热流密度 `W/m^2`

适用场景：

1. 对齐旧 Fortran 工程工作流
2. 快速做系统级瞬态计算
3. 当外热流表已经由上游工具或历史模型给定时

### 3.2 解析/半解析路线

当前这条路线并不是“三种轨道热流全都严格解析求解”。

目前实际状态是：

1. 太阳直射  
   使用 `OrbitalHeatSource`，属于已实现的解析模型。

2. 地球反照  
   如果使用 `AlbedoHeatSource`，当前走的是简化常值模型。

3. 地球红外  
   如果使用 `EarthIRHeatSource`，当前走的也是简化常值模型。

因此，更准确的说法应该是：

“当前模块支持太阳热流解析计算，以及反照/红外的简化模型；严格积分版反照和红外尚未实现。”

---

## 4. 类结构

### 4.1 热流模型层

- `BaseExternalHeatSource`  
  所有外热流模型的抽象基类，返回热流密度 `W/m^2`

- `OrbitalHeatSource`  
  太阳直射解析模型

- `AlbedoHeatSource`  
  地球反照简化模型

- `EarthIRHeatSource`  
  地球红外简化模型

- `IntegralAlbedoHeatSource`  
  严格地球反照积分模型预留接口

- `IntegralEarthIRHeatSource`  
  严格地球红外积分模型预留接口

- `OrbitalTableHeatSource`  
  Fortran 风格查表热流模型

- `CompositeHeatSource`  
  热流叠加器

### 4.2 边界封装层

- `ExternalHeatFluxBC`  
  将热流密度 `W/m^2` 乘以边界面积，转换成求解器真正使用的热流 `W`

---

## 5. 与 BoundaryRegion 的关系

`ExternalHeatFluxBC` 的设计目标是保持与当前导热边界体系兼容。

工作方式如下：

1. `BaseHeatConduction.get_derivatives()` 会更新边界内部状态
2. `BoundaryRegion.update_internal_state(...)` 会遍历边界条件
3. 如果某个边界条件带有 `update_state()`，则会在每次边界刷新时自动调用
4. `ExternalHeatFluxBC.update_state()` 内部再调用 `heat_source.get_heat_flux(current_time)`

这样就实现了：

1. 时间相关轨道热流自动更新
2. 用户不需要在外部手工调用热流更新函数

---

## 6. 当前推荐用法

### 6.1 用查表法

```python
from Components.ExternalHeatSources import OrbitalTableHeatSource, ExternalHeatFluxBC

shape = boundary.shape

table_source = OrbitalTableHeatSource(
    shape=shape,
    table_ids=4,
    scale_factor=1.0,
    offset=0.0,
    periodic=True,
)

external_bc = ExternalHeatFluxBC(
    heat_source=table_source,
    area_array=boundary.area,
)

boundary.conditions.append(external_bc)
```

### 6.2 用太阳解析法

```python
from Components.ExternalHeatSources import OrbitalHeatSource, ExternalHeatFluxBC

shape = boundary.shape

solar_source = OrbitalHeatSource(
    shape=shape,
    solar_constant=1361.0,
    orbit_height=800.0,
    orbit_period=7644.0,
    orbit_inclination=0.0,
    surface_normal_angles=(0.0, 0.0),
)

external_bc = ExternalHeatFluxBC(
    heat_source=solar_source,
    area_array=boundary.area,
)

boundary.conditions.append(external_bc)
```

### 6.3 用组合热流

注意：这里的反照和红外当前仍是简化模型。

```python
from Components.ExternalHeatSources import (
    OrbitalHeatSource,
    AlbedoHeatSource,
    EarthIRHeatSource,
    CompositeHeatSource,
    ExternalHeatFluxBC,
)

shape = boundary.shape
source = CompositeHeatSource(shape)
source.add_source(OrbitalHeatSource(shape=shape))
source.add_source(AlbedoHeatSource(shape=shape, albedo_factor=0.3))
source.add_source(EarthIRHeatSource(shape=shape, earth_ir_flux=237.0))

external_bc = ExternalHeatFluxBC(
    heat_source=source,
    area_array=boundary.area,
)

boundary.conditions.append(external_bc)
```

---

## 7. 关于 HPwithFin 的配合关系

在 `HPwithFin` 中，轨道外热流的来源和翅片如何吸收这部分热流是两层独立选择：

1. 热流来源  
   - 查表法
   - 太阳解析法
   - 太阳解析 + 反照简化 + 红外简化

2. 翅片加载方式  
   - `lumped_root_area`：把翅片受照折算回冷凝段根部边界
   - `distributed_fin_absorption`：翅片受照直接进入翅片降维方程

因此可以组合成：

1. 查表法 + 简化翅片加载
2. 查表法 + 翅片直接受照
3. 太阳解析法 + 简化翅片加载
4. 太阳解析法 + 翅片直接受照

---

## 8. 当前限制

1. `IntegralAlbedoHeatSource` 尚未实现
2. `IntegralEarthIRHeatSource` 尚未实现
3. `AlbedoHeatSource` 与 `EarthIRHeatSource` 目前只是简化模型，不能等同于论文或教材中的严格积分公式实现
4. 查表法不提供三种热流分量的分离信息，只负责按时间输出总热流

---

## 9. 后续建议

如果后续继续完善解析路线，推荐按下面顺序推进：

1. 先实现 `IntegralAlbedoHeatSource`
2. 再实现 `IntegralEarthIRHeatSource`
3. 最后再考虑把姿态、遮挡、受照投影修正进一步细化

这样可以保持：

1. 当前工程算例继续可跑
2. 研究版解析模型逐步替换简化模型
3. 现有 `BoundaryRegion` 和 `HPwithFin` 耦合接口不需要再重构
