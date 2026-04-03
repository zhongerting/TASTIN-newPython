# TASTIN 外热流模块参数说明

## 一、问题背景

TASTIN（热离子空间核电源核热电耦合计算分析软件）在进行轨道空间核电源系统仿真时，需要考虑来自太空环境的外部热流边界条件。在 Fortran 原版程序中，热管辐射器与太阳辐射之间通过角度参数建立几何关系。

## 二、Fortran 原版角度体系分析

### 2.1 两种角度概念

根据对 Fortran 源代码的分析，原版程序中存在**两种不同的角度概念**：

| 角度 | 变量名 | 物理含义 | 用途 |
|------|--------|----------|------|
| **热管倾斜角** | `HPAngle` / `GLoHPAngle` | 热管轴线与水平面的夹角 | 重力压降计算 |
| **表面法向角** | `(w0, l0)` | 热管外表面的法向方向 | 太阳辐射入射角计算 |

### 2.2 HPAngle - 热管倾斜角

**定义位置**：[Fortran源文件/HPInput.for:20](file:///e:/项目任务/五院-电源/source_code/Fortran源文件/HPInput.for#L20)

```fortran
real(8) :: GLoHPAngle  ! 倾斜角 [rad]
```

**使用位置**：[Fortran源文件/HPSSInitial.for:110](file:///e:/项目任务/五院-电源/source_code/Fortran源文件/HPSSInitial.for#L110)

```fortran
dPGra = dPGra + Dens*Gravity*temp2*sin(HPAngle*Pi/180.)
```

**物理意义**：
- 热管轴线与水平面的夹角
- 用于计算**重力压降**分量
- `sin(HPAngle)` 表示重力在工质流动方向上的分量

### 2.3 (w0, l0) - 表面法向角

**定义位置**：[Fortran源文件/sun_heat.f90:13](file:///e:/项目任务/五院-电源/source_code/Fortran源文件/sun_heat.f90#L13)

```fortran
real(8)::l0,w0  ! 微元法向
```

**使用公式**：[Fortran源文件/sun_heat.f90:80](file:///e:/项目任务/五院-电源/source_code/Fortran源文件/sun_heat.f90#L80)

```fortran
X = cos(iz_s)*sin(w0) + sin(iz_s)*cos(w0)*cos(l0-l_s)
```

**物理意义**：

| 参数 | 名称 | 定义 |
|------|------|------|
| **w0** | 极角 (Polar angle) | 表面法向与 z 轴的夹角 [rad] |
| **l0** | 方位角 (Azimuthal angle) | 表面法向在 x-y 平面投影与 x 轴的夹角 [rad] |

**太阳入射因子 X**：
- X = 1：太阳直射（法向完全朝向太阳）
- X = 0：完全阴影（法向背向太阳）
- 0 < X < 1：部分入射

### 2.4 Fortran 原版太阳热流计算特点

从 `sun_heat.f90` 来看，原版程序中：

1. 太阳热流计算是**独立模块**（标量函数 `qqqqq`）
2. 需要外部确定每根热管的 `(w0, l0)`
3. 或使用**简化假设**（所有热管使用同一代表性角度）

## 三、Python 重构版参数对照

### 3.1 参数映射表

| Fortran 参数 | Python 参数 | 数据类型 | 说明 |
|-------------|------------|----------|------|
| `GLoHPAngle` | `hp_tilt_angle` | float | 热管倾斜角 [rad] |
| `w0` | `surface_normal_angles[0]` | float/array | 极角 [rad] |
| `l0` | `surface_normal_angles[1]` | float/array | 方位角 [rad] |
| `S` | `solar_constant` | float | 太阳常数 [W/m²]，默认 1361 |
| `h` | `orbit_height` | float | 轨道高度 [km]，默认 800 |
| `TT` | `orbit_period` | float | 轨道周期 [s]，默认 7644 |
| `i_s` | `orbit_inclination` | float | 轨道倾角 [deg]，默认 0 |

### 3.2 坐标系说明

```
        z (轨道平面法向)
        ↑
        |   / 表面法向
        |  /
        | /
        |/_________→ x (轨道方向)
       /\
      /  \
     /    \
    /      \
    ↓       → y (指向地球)
   (热管)
```

### 3.3 典型参数组合

| 场景 | w0 | l0 | 说明 |
|------|-----|-----|------|
| 热管垂直于轨道平面 | π/2 | 0 | 法向在 xy 平面指向 +x |
| 热管指向地球方向 | π/2 | π | 法向在 xy 平面指向 -y |
| 热管倾斜向太阳 | π/2 - 10° | 0 | 略偏离 +x 方向 |

## 四、Python 参数设置方案

### 4.1 方案1：单根代表性热管（简化计算）

适用于初步验证或单根热管测试：

```python
from Components.ExternalHeatSources import OrbitalHeatSource

solar_source = OrbitalHeatSource(
    shape=(n_con,),
    orbit_height=800.0,           # 800 km LEO
    orbit_period=7644.0,          # 约127分钟
    orbit_inclination=0.0,        # 太阳同步轨道
    surface_normal_angles=(np.pi/2, 0.0)  # w0=90°, l0=0° → 法向指向 +x
)
```

### 4.2 方案2：多根热管扇形排列（精确计算）

适用于完整辐射器阵列模拟：

```python
import numpy as np

n_pipes = 8          # 8根热管
n_con = 12           # 每根管12个节点
angle_step = np.pi / 4  # 45度间隔

for i in range(n_pipes):
    l0 = i * angle_step          # 0°, 45°, 90°, 135°, ...
    w0 = np.pi / 2               # 都在 xy 平面内

    solar_source_i = OrbitalHeatSource(
        shape=(n_con,),
        orbit_height=800.0,
        orbit_period=7644.0,
        orbit_inclination=0.0,
        surface_normal_angles=(w0, l0)
    )

    # 创建边界条件并集成到对应热管
    external_bc_i = ExternalHeatFluxBC(
        heat_source=solar_source_i,
        area_array=hp_radiator_i.hp.boundaries['outer_con'].area * 0.5
    )
    hp_radiator_i.hp.boundaries['outer_con'].conditions.append(external_bc_i)
```

### 4.3 方案3：考虑热管安装倾角

如果热管安装时有倾斜（如向地球方向倾斜10度）：

```python
tilt_angle = 10 * np.pi / 180  # 10度倾角

for i in range(n_pipes):
    l0 = i * np.pi / 4              # 基础方位角
    w0 = np.pi / 2 - tilt_angle     # 考虑倾斜后的极角

    solar_source_i = OrbitalHeatSource(
        shape=(n_con,),
        orbit_height=800.0,
        orbit_period=7644.0,
        orbit_inclination=0.0,
        surface_normal_angles=(w0, l0)
    )
```

## 五、与其他模块的集成

### 5.1 与 HPwithFin 的集成

当前 `HPwithFin` 类已支持通过 `ExternalHeatFluxBC` 添加外热流边界：

```python
# 在 HPwithFin 创建后添加外热流
outer_con_boundary = hp_radiator.hp.boundaries['outer_con']

solar_source = OrbitalHeatSource(
    shape=outer_con_boundary.shape,
    orbit_height=800.0,
    orbit_period=7644.0,
    orbit_inclination=0.0,
    surface_normal_angles=(np.pi/2, 0.0)
)

solar_bc = ExternalHeatFluxBC(
    heat_source=solar_source,
    area_array=outer_con_boundary.area * 0.5  # 0.5系数
)

outer_con_boundary.conditions.append(solar_bc)
```

### 5.2 与 HPAngle 的集成建议

未来可在 `HPwithFin` 中添加 `hp_tilt_angle` 参数用于重力压降计算：

```python
class HPwithFin:
    def __init__(self,
                 ...,
                 hp_tilt_angle: float = 0.0,  # 新增
                 ...):
        self.hp_tilt_angle = hp_tilt_angle
```

## 六、版本信息

- **创建日期**：2026-04-02
- **参考版本**：TASTIN Fortran 原版
- **文档版本**：v1.0
