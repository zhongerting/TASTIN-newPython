# Materials 子系统维护手册

> 本文是 Codex 进入 `Materials/` 后的首选阅读入口。内容按当前源码核验，覆盖材料类、数值边界、上级调用关系和维护风险。修改公式、边界条件、构造参数或方法签名时，必须回看源码并同步更新本文。

## 导航

- [Codex 启动阅读指南](#1-codex-启动阅读指南)
- [目录和继承体系](#2-目录和继承体系)
- [基类接口矩阵](#3-基类接口矩阵)
- [流体材料](#4-流体材料)
- [固体和等效材料](#5-固体和等效材料)
- [上级调用关系](#6-上级调用关系)
- [向量化约定](#7-向量化约定)
- [已知风险与维护规则](#8-已知风险与维护规则)
- [最小验证清单](#9-最小验证清单)

## 1. Codex 启动阅读指南

先读本文，再按任务只打开必要源码：

| 任务 | 首选源码 |
| --- | --- |
| 修改基类接口、兼容回退或抽象约束 | `Base.py` |
| 修改冷却剂、水力网络所需物性 | `Fluids/Sodium.py`、`Fluids/Potassium.py`、`Fluids/SodiumPotassium78.py` |
| 修改热管冻结/融化或等效热导 | `Solids/NaHP.py`、`Solids/KHP.py`、`Solids/WickMaterial.py`、`Components/basicComponents/HeatPipe2D.py` |
| 修改堆芯固体或反射层 | 对应 `Solids/*.py`；复合芯块还需看 `Solids/CompositePelletMaterial.py` 和 `Components/TFEUnit.py` |
| 修改壁材或管壁导热 | `Solids/WallMaterial.py` |
| 修改流体调用约定 | `Solvers/Hydrodynamics/Components.py`、`Solvers/Hydrodynamics/HydraulicNetwork.py`、`Correlations/Correlations.py` |
| 修改固体调用约定 | `Solvers/HeatConduction/HeatConduction.py` |

以下场景不能只依赖本文：修改经验公式、温度分段、截断规则、查表范围、单位换算、向量化行为、缓存失效规则、公开方法签名或上级调用路径。

## 2. 目录和继承体系

```text
Materials/
|-- Base.py
|-- Fluids/
|   |-- Sodium.py
|   |-- Potassium.py
|   |-- SodiumPotassium78.py
|   `-- NaK78.py
`-- Solids/
    |-- B4C.py
    |-- Beryllium.py
    |-- BerylliumOxide.py
    |-- CompositePelletMaterial.py
    |-- GasGaps.py
    |-- KHP.py
    |-- MoNb.py
    |-- Molybdenum.py
    |-- NaHP.py
    |-- StainlessSteel.py
    |-- UO2.py
    |-- WallMaterial.py
    |-- WickMaterial.py
    |-- WickStructure.py
    `-- ZrH.py
```

```text
Material(ABC)
|-- SolidMaterial(ABC)
|   `-- Solids/ 下的全部材料
`-- FluidMaterial(ABC)
    `-- Sodium、Potassium、SodiumPotassium78、未完成的 NaK78

TwoPhaseFluidMaterial = FluidMaterial
```

`TwoPhaseFluidMaterial` 只是兼容别名，不是额外的一层类型。

## 3. 基类接口矩阵

### 3.1 `SolidMaterial`

| 方法 | 类型 | 约定 |
| --- | --- | --- |
| `conductivity(T)` | 抽象接口 | `k [W/(m*K)]` |
| `heat_capacity(T)` | 抽象接口 | `Cp [J/(kg*K)]` |
| `density(T)` | 抽象接口 | `rho [kg/m^3]` |
| `emissivity(T)` | 默认实现 | `epsilon = 0.8`；数组输入返回同形数组 |
| `diffusivity(T)` | 派生属性 | `alpha = k / (rho * Cp) [m^2/s]` |
| `saturation_pressure(T)`、`vapor_viscosity(T)`、`latent_heat(T)` | 可选占位 | 基类为 `pass`；供热管工质扩展 |

### 3.2 `FluidMaterial`

必须由子类实现：

| 方法 | 单位或用途 |
| --- | --- |
| `density(T, P)` | `rho [kg/m^3]` |
| `viscosity(T, P=0)` | `mu [Pa*s]` |
| `heat_capacity(T, P=0)` | `Cp [J/(kg*K)]` |
| `conductivity(T, P=0)` | `k [W/(m*K)]` |
| `enthalpy(T, P=0)` | `h [J/kg]` |
| `temperature_from_enthalpy(h, P)` | `T [K]` |
| `saturation_pressure(T)`、`saturation_temperature(P)`、`latent_heat(T)` | 两相接口 |
| `liquid_density_derivative_P(P)`、`liquid_density_derivative_T(T)` | 水力网络所需密度导数 |

基类提供：

| 方法 | 类型 | 行为 |
| --- | --- | --- |
| `prandtl_number(T, P=0)` | 派生属性 | `Pr = mu * Cp / k` |
| `liquid_density(T, P=0)` | 兼容层 | 优先 `1 / specific_volume_liquid_sat(T)`，否则回退 `density(T, P)` |
| `vapor_density(T, P=0)` | 兼容层 | 需要 `specific_volume_vapor_sat(T)`，否则抛出 `NotImplementedError` |
| `liquid_viscosity(T, P=0)` | 兼容层 | 优先 `viscosity_liquid(T)`，否则回退 `viscosity(T, P)` |
| `vapor_viscosity(T, P=0)` | 兼容层 | 优先 `viscosity_vapor(T)`，否则抛出 `NotImplementedError` |
| `liquid_enthalpy_sat(T)` | 兼容层 | 优先 `enthalpy_saturated_liquid(T)`，否则回退 `enthalpy(T, saturation_pressure(T))` |
| `vapor_enthalpy_sat(T)` | 兼容层 | 优先 `enthalpy_saturated_vapor(T)`，否则回退 `liquid_enthalpy_sat(T) + latent_heat(T)` |

## 4. 流体材料

### 4.1 总览

| 类 | 构造 | 状态与用途 | 关键注意事项 |
| --- | --- | --- | --- |
| `Sodium` | `Sodium()` | 完整度最高的钠流体模型；用于可压缩/不可压缩水力示例 | `P <= 0` 时 `saturation_temperature()` 返回 `371 K` |
| `Potassium` | `Potassium()` | 钾流体模型；含固、液、气分段和若干扩展性质 | 部分接口仍是返回 `0.0` 的占位实现 |
| `SodiumPotassium78` | `SodiumPotassium78()` | 当前 NaK78 液态冷却剂主用类 | 只应按液态冷却剂使用；继承的两相兼容层并不完整 |
| `NaK78` | 不应实例化 | 未完成草稿 | 只有 `density()`，且引用未定义的 `specific_volume()`；仍缺多个抽象接口 |

### 4.2 `Sodium`

核心边界和公式：

| 性质 | 实现 |
| --- | --- |
| `T_crit`、`P_crit` | `2503.7 K`、`25.64e6 Pa` |
| 液态 `Cp` | `1097.73 - 0.556577*T + 3.43167e-4*T^2` |
| 饱和压 | `exp(11.9463 - 12633.73/T - 0.4672*ln(T)) * 1e6 Pa`；`T <= 0` 返回 `0` |
| 汽化潜热 | `[393.37*x + 4398.6*x^0.29302] * 1000`，`x = max(1 - T/T_crit, 0)` |
| 饱和液比体积 | `1 / [219 + 275.32*x + 511.58*sqrt(x)]` |
| 液态导热 | `abs(13.9723 + 0.0331088*T - 2.22398e-5*T^2)` |
| 液态粘度 | 三次多项式绝对值 |
| 表面张力 | `240.5 * max(1 - T/T_crit, 0)^1.126 / 1000` |

`specific_volume(T, P)` 在 `T < T_sat(P)` 时使用液态饱和值，否则使用 20 次迭代的过热蒸汽维里模型。`temperature_from_enthalpy(h, P)` 按过冷液、两相区、过热气体分段：两相区直接返回 `T_sat`，液相拟合结果下限为 `300 K`。

### 4.3 `Potassium`

核心边界：

| 性质 | 实现 |
| --- | --- |
| 熔点 | `T_melt = 336.35 K` |
| 饱和压 | `4.0168e11 * 10^(-4625.3/T) / T^0.7 * 1e6 Pa`；`T <= 0` 返回 `0` |
| 饱和温度 | 牛顿迭代反求；`P <= 1e-5 Pa` 返回熔点 |
| 密度 | `t_c < 63.2` 固相；`63.2 <= t_c < 756.5` 液相；`756.5 <= t_c < 2100` 气相；更高温返回 `0` |
| 液态粘度 | `336.35 < T < 653.15 K` 与其他温区使用不同参数 |
| 通用导热 | 液相公式 `43.8 - 2.22e-2*T + 3.95e3/max(T,1)` |
| 通用 `Cp` | 液相公式 |

扩展接口包括固态导热、蒸汽粘度、蒸汽导热、压缩系数、声速、表面张力和过热蒸汽压力。`specific_volume_gas_working()`、`temperature_from_pressure_enthalpy()`、`entropy()`、`conductivity_gas()` 仍返回 `0.0`，不能视为已实现。

### 4.4 `SodiumPotassium78`

该类是 NaK78 冷却剂主用实现，按液态条件使用：

| 性质 | 实现 |
| --- | --- |
| 组成 | 摩尔分数 `N_Na=0.324`、`N_K=0.676`；质量分数 `m_Na=0.22`、`m_K=0.78` |
| 比体积 | `1.003 * (N_Na/rho_Na + N_K/rho_K)` |
| 密度 | `1 / specific_volume(T, P)` |
| 粘度 | 以 `673.15 K` 为中心，用 `tanh(0.1*(T-673.15))` 平滑混合两套经验式 |
| 导热 | `21.4 + 2.07e-2*t_c - 2.2e-5*t_c^2` |
| `Cp` | `970.688 - 0.3690288*t_c + 3.43088e-4*t_c^2` |
| 焓 | `145.607 + 970.688*t_c - 0.1845144*t_c^2 + 1.1436267e-4*t_c^3`，量纲按实现用途为 `J/kg` |
| 焓反求温度 | 对 `h/1e5` 使用六阶拟合 |
| 声速 | Na、K 声速按质量分数加权 |
| 绝热压缩系数 | `1 / (rho*c^2)` |

风险：

- 类注释明确只考虑液态；不要把继承来的饱和压、潜热等钠式近似误解为完整 NaK78 两相模型。
- 构造函数当前仍设置 `name="Sodium"`、`formula="Na"`，这是现状，不代表材料组成。
- `surface_tension()` 注释称占位拟合，代码实际为 `0.119 - 8e-5*T`。
- 没有 `specific_volume_vapor_sat()`；调用基类 `vapor_density()` 会失败。

## 5. 固体和等效材料

### 5.1 堆芯与结构材料

| 类 | 构造 | `rho [kg/m^3]` | `k [W/(m*K)]` 与边界 | `Cp [J/(kg*K)]` 与边界 |
| --- | --- | --- | --- | --- |
| `UO2` | `UO2()` | `10400` | 公式输入钳位到 `T >= 1500 K` | 公式输入钳位到 `T >= 700 K` |
| `MoNb` | `MoNb()` | `11506` | `120 - 0.014*T` | `65.05 + 0.0133*T` |
| `Molybdenum` | `Molybdenum()` | `10200` | `(110 - 0.015*T)*0.9` | `101.8 + 0.025*T` |
| `Beryllium` | `Beryllium()` | `1830` | 多项式输入钳位到 `[198.1, 1556] K` | `T <= 223.1 K` 返回 `142.3`；`T >= 1556 K` 返回 `357.7`；中间为多项式 |
| `BerylliumOxide` | `BerylliumOxide()` | `2800` | `4e7/T^2`；`T <= 0` 返回 `4e7` | `1343 + 1.8469*T` |
| `ZirconiumHydride` | `ZirconiumHydride()` | `5615` | `20` | `187 + 0.745*T` |
| `BoronCarbide` | `BoronCarbide()` | `2300` | `400..1000 K` 查表，越界钳到端点 | `300..2500 K` 查表，越界钳到端点 |
| `AusteniticStainlessSteel` | `AusteniticStainlessSteel()` | `8084 - 0.4209*T - 3.894e-5*T^2` | `9.2 + 0.0175*T - 2e-6*T^2` | `T <= 0` 返回 `472`，否则 `472 + 0.136*T - 2.82e6/T^2` |

### 5.2 壁材

`WallMaterial.py` 先将温度钳到 `T >= 1e-3 K`：

| 类 | 构造 | 关键行为 |
| --- | --- | --- |
| `SS321` | `SS321(name="SS321")` | `rho=8090`，`k=14.5+0.015*(T-273)`，`Cp=500` |
| `SS316` | `SS316(name="SS316")` | `rho=7900`，`k=max(15.61342+0.01324*(T-273), 1e-6)`，`Cp=509.77978+0.14008*(T-273)` |
| `SS316H` | `SS316H(name="SS316H")` | `rho=7900`，导热最小值 `1e-6`，`Cp=472+0.136*T-2820000/T^2` |
| `Haynes` | `Haynes(name="Haynes")` | `rho=8180`；`k`、`Cp` 在 `773.15 K` 分段 |

### 5.3 气隙等效材料

`GasGaps.py` 中的气体按固体等效材料使用，不是水力流体：

| 类 | `rho` | `k` | `Cp` |
| --- | --- | --- | --- |
| `CarbonDioxide()` | `1.1` | `0.0064 + 1.44e-4*T` | `880` |
| `Helium()` | `4/22.4` | `3.1e-4*T + 0.059` | `5190` |
| `Xenon()` | `1.0` | `0.3` | `158` |
| `Cesium()` | `1.0` | `0.15` | `100` |

### 5.4 热管工质

`SodiumHP` 和 `PotassiumHP` 继承 `SolidMaterial`，用于热管导热模型，不是水力网络流体。两者都在熔点附近使用表观热容法：

```text
f = clip((T - T_melt + dt_thaw) / (2*dt_thaw), 0, 1)
Cp_mushy = Cp_solid + f*(Cp_liquid - Cp_solid) + H_sf/(2*dt_thaw)
```

| 类 | 默认构造 | `T_melt` | `H_sf [J/kg]` | `dt_thaw` | 糊状区 |
| --- | --- | --- | --- | --- | --- |
| `SodiumHP` | `SodiumHP(name="Sodium_Wick_Fluid")` | `371.0 K` | `114.7e3` | `1.0 K` | `[370, 372] K` |
| `PotassiumHP` | `PotassiumHP(name="Potassium_Wick_Fluid")` | `336.35 K` | `553.8 cal/mol = 5.926e4 J/kg` | `1.0 K` | `[335.35, 337.35] K` |

两者还提供 `saturation_pressure(T)`、`vapor_viscosity(T)`、`latent_heat(T)` 和 `molar_mass`，供 `WickMaterial` 计算赝热导；`PotassiumHP` 另提供 `saturated_liquid_density(T)`、`vapor_density(T)`、`specific_volume_liquid_sat(T)`、`specific_volume_vapor_sat(T)`、`vapor_heat_capacity(T)`、`vapor_gas_constant()`、`vapor_heat_capacity_ratio()`、`vapor_ideal_cp()`、`vapor_ideal_cv()`、`vapor_sound_speed(T)`、`liquid_viscosity(T)`、`viscosity(T)` 和 `surface_tension(T)`。

`PotassiumHP` 的基础常数及下列物性已按 2026-06-04/2026-06-05 用户提供数据逐项更新；每个公式的适用范围和参考文献已同步写入源码注释。2026-06-05 起，`molar_mass` 使用 NIST 分子量 `3.90983e-2 kg/mol`，替代此前 `39.1e-3 kg/mol` 的圆整值。

`PotassiumHP.density()` 的固态钾密度已按 2026-06-04 用户提供公式更新：`rho_s=857.6/[1+2.39e-4*(T-300)] kg/m3`，来源建议适用范围为 `270-320 K`；液态/饱和液态钾密度已按用户最新提供公式更新：`rho_l_sat=890.29-2.113e-1*T kg/m3`，来源建议适用范围为 `623-1123 K`，高于常压沸点时应理解为相应饱和压力下的饱和液状态。`PotassiumHP.saturated_liquid_density()` 返回同一液态公式，`specific_volume_liquid_sat()` 返回其倒数。
`PotassiumHP.conductivity()` 的固态钾导热系数已按 2026-06-04 用户提供公式更新：`lambda_s=102.5-1.04e-1*(T-298.2) W/(m*K)`，来源建议适用范围为 `270-336.8 K`；液态钾导热系数已按用户提供公式更新：`lambda_l=66.09-3.579e-2*T W/(m*K)`，来源建议适用范围为 `336.8-1000 K`。
`PotassiumHP.heat_capacity()` 的固态钾定压比热已按 2026-06-04 用户提供 Shomate 公式更新：`tau=T/1000`，`cp_s=(-63.47410-3226.340*tau+14644.60*tau^2-16229.50*tau^3+16.29410*tau^-2)/3.90983e-2 J/(kg*K)`，来源建议适用范围为 `298-336.35 K`；液态钾定压比热也已按用户提供 Shomate 公式更新：`cp_l=(40.27113-30.54542*tau+26.49505*tau^2-5.727854*tau^3-0.063477*tau^-2)/3.90983e-2 J/(kg*K)`，来源建议适用范围为 `336.35-1039.54 K`。
`PotassiumHP.saturation_pressure()` 已按 2026-06-05 用户提供宽温区公式更新为单一关系式：`ln(p_sat)=25.109-10488/T-0.448*ln(T)`，`p_sat` 单位为 `Pa`，来源建议适用范围为 `350-1000 K` 的钾气液共存线。该式替代此前低温式、NIST Antoine 式和平滑拼接方案；`900 K` 以上给出总饱和压力，若需单原子/双原子分压需额外引入蒸汽化学平衡。
`PotassiumHP.vapor_density()` 已按 2026-06-04 用户最新提供饱和汽密度公式更新：`rho_v_sat=2.398e3*exp(-8.698e3/T) kg/m3`，来源建议适用范围为 `623-1123 K` 的饱和钾蒸汽；高于常压沸点时应理解为相应饱和压力下的饱和汽状态，固态区域、临界区、高精度状态方程或明显二聚体/多聚体效应工况不建议使用。`PotassiumHP.specific_volume_vapor_sat()` 返回该密度的倒数。
`PotassiumHP.vapor_viscosity()` 已按 2026-06-05 用户提供饱和蒸汽粘度公式更新：`mu_v=5.450e-6+2.830e-8*T-5.600e-12*T^2 Pa*s`，来源建议适用范围为 `350-1000 K`，用于低压钾蒸汽气相流动阻力和动量输运；该式严禁直接用于液态钾。`900 K` 以上钾蒸汽中 `K2` 缔合体比例增加，本多项式为平均非理想性平滑拟合；强非平衡等离子体、电离区或超高压稠密区需额外碰撞积分修正。
`PotassiumHP.vapor_heat_capacity()` 已按 2026-06-04 用户提供气相 Shomate 公式新增：`tau=T/1000`，`cp_v=(20.66122+0.391869*tau-0.417344*tau^2+0.145582*tau^3+0.003764*tau^-2)/3.90983e-2 J/(kg*K)`，来源建议适用范围为 `1039.54-1800 K`；低压或中低密度钾蒸汽可近似为单原子理想气体 `cp_v≈531.6 J/(kg*K)`，临界区、高压稠密蒸汽或明显二聚体/多聚体效应工况不建议使用。
`PotassiumHP` 的气相校核接口已按 2026-06-05 用户提供公式新增：`vapor_gas_constant()=R_u/M_K=212.65 J/(kg*K)`，`vapor_heat_capacity_ratio()=5/3`，`vapor_ideal_cp()=5R_K/2≈531.6 J/(kg*K)`，`vapor_ideal_cv()=3R_K/2≈319.0 J/(kg*K)`，`vapor_sound_speed(T)=sqrt(gamma_K*R_K*T)≈18.83*sqrt(T) m/s`。这些接口仅适用于 `700-1800 K` 低压或中低密度单原子钾蒸汽工程校核；液态、两相区、临界区、高压稠密蒸汽不得直接使用，饱和两相流应采用两相声速模型。
`PotassiumHP.liquid_viscosity()` 已按 2026-06-04 用户提供公式新增：`mu_l=4.293e-5*exp(1017.72/T) Pa*s`，来源建议适用范围为 `623-1123 K`；高于常压沸点时需保证系统压力维持液态或饱和液态。`PotassiumHP.viscosity()` 作为液态动力黏度别名返回同一结果。
`PotassiumHP.surface_tension()` 已按 2026-06-04 用户提供公式新增：`sigma_l=1.3794e-1-6.927e-5*T N/m`，来源建议适用范围为 `623-1123 K`；高于常压沸点时需保证系统压力维持液态或饱和液态，固态、临界区或高精度界面稳定性计算不建议使用。
`PotassiumHP.latent_heat()` 已按 2026-06-05 用户提供宽温区多项式更新为单一关系式：`h_fg=2.487169e6-396.5976*T-0.102412*T^2 J/kg`，来源建议适用范围为 `350-900 K` 的饱和线数据拟合，并以 Watson 关联式 `h_fg,b=1.966838e6 J/kg`、`T_b=1032 K`、`T_c=2223 K` 作为校核。该式替代此前低温线性式、高温推导式和平滑拼接方案；不适用于过热蒸汽、过冷液体、固态升华、临界区或高压稠密超临界钾。

### 5.5 `WickMaterial`

构造：

```python
WickMaterial(name, solid_mat, fluid_mat, porosity, r_vapor, r_in_wall)
```

其中 `fluid_mat` 通常为 `SodiumHP` 或 `PotassiumHP`。基础物性：

```text
rho_eff = rho_fluid
Cp_eff = [phi*rho_f*Cp_f + (1-phi)*rho_s*Cp_s] / rho_f
k_axial = min(k_structural + k_pseudothermal, conductivity_cap)
k_radial = k_structural
```

结构等效导热：

```text
k_structural = k_f * [k_s+k_f-(1-phi)*(k_f-k_s)]
                     / [k_f+k_s+(1-phi)*(k_f-k_s)]
```

赝热导：

```text
pse1 = [r_vapor^4 / (r_in_wall^2-r_vapor^2)]
       * (h_fg*M_g*P_sat)^2 / [4*mu_v*R_gas^2*T^3]
v_ave = sqrt(8*R_gas*T / (pi*M_g))
k_pseudothermal = pse1 / [1 + (8/3)*mu_v*v_ave/(r_vapor*P_sat)]
```

数值规则：

- 默认 `conductivity_cap = 5e7 W/(m*K)`；`set_conductivity_cap(None)` 可取消上限。
- 默认启用惰性构建的 `8192` 点查表，范围为 `[250, min(fluid.T_crit, 2500)] K`。
- 总导热和结构导热按对数空间插值；查表范围外的温度钳到端点。
- `set_conductivity_cap()` 会调用 `invalidate_lookup_table()`。若直接修改流体、骨架、几何、孔隙率或查表参数，必须手动调用 `invalidate_lookup_table()`。
- `is_high_nonlinearity_temperature(T)` 使用 `|d(log(k))/dT| >= 0.05 * max_gradient` 的查表区间，供 `HeatPipe2D` 决定是否强制刷新物性。
- `r_in_wall^2 == r_vapor^2` 会使赝热导公式奇异；构造时没有显式几何校验。

`WickStructure()` 是另一个简化等效材料：`rho=1.0`、`k=0.15`、`Cp=100.0`。其常数方法返回标量并依赖 NumPy 广播。

### 5.6 `CompositePelletMaterial`

```python
CompositePelletMaterial(
    fuel_mat, reflector_mat, shape_nodes,
    n_lower, n_active, n_upper
)
```

该类按二维网格的轴向切片分发 `conductivity()`、`heat_capacity()`、`density()`：

```text
[:, :idx_1]       -> reflector_mat
[:, idx_1:idx_2]  -> fuel_mat
[:, idx_2:]       -> reflector_mat
```

构造时要求 `n_lower + n_active + n_upper == shape_nodes[1]`，否则抛出 `ValueError`。输出始终扁平化，以适配导热求解器缓存。

## 6. 上级调用关系

实际源码调用路径：

| 上级模块 | Materials 契约 |
| --- | --- |
| `Solvers/HeatConduction/HeatConduction.py` | 对 `SolidMaterial` 调用 `conductivity(T)`、`density(T)`、`heat_capacity(T)`；输入是温度数组 |
| `Solvers/Hydrodynamics/Components.py` | 对流体调用 `enthalpy(T,P)`、`temperature_from_enthalpy(h,P)`、`density(T,P)`、`viscosity(T,P)`、`heat_capacity(T,P)` 和密度导数 |
| `Solvers/Hydrodynamics/HydraulicNetwork.py` | 批量调用 `Cp/rho/mu`、密度导数和焓反求温度 |
| `Correlations/Correlations.py` | 调用兼容层 `liquid_density()`、`vapor_density()`、`liquid_viscosity()`、`vapor_viscosity()`、`liquid_enthalpy_sat()`、`vapor_enthalpy_sat()` |
| `Components/basicComponents/HeatPipe2D.py` | 内部构造 `WickMaterial`；调用轴向/径向导热、密度、比热和高非线性区判断 |
| `Components/TFEUnit.py` | 可构造 `CompositePelletMaterial`；材料字典使用 `UO2`、`BerylliumOxide`/`BeO`、`MoNb`、`Molybdenum` |
| `Components/basicComponents/Fuel.py`、`Electord.py` | 缺省实例分别为 `UO2()`、`MoNb()`、`Molybdenum()` |

冷却回路和测试脚本中，旧场景仍常见 `Sodium()`；较新的 NaK78 场景直接实例化 `SodiumPotassium78()`。判断模型用途时应以实际构造点为准。

## 7. 向量化约定

- 主流材料方法普遍接受标量和 `np.ndarray`，但返回类型不完全统一：有些标量输入返回 NumPy 0 维数组，有些显式转回 `float`，常量材料可能直接返回 Python 标量。
- 导热求解器依赖数组赋值和 NumPy 广播；不要假设每个材料方法都返回与输入严格同形的数组。
- `CompositePelletMaterial` 会将数组重塑为构造时的 `shape_nodes`，形状不匹配会失败。
- `Sodium.temperature_from_enthalpy()` 会展平输入并恢复 `h` 的原始形状；压力数组应与焓数组兼容。
- 修改向量化实现后，至少测试标量、0 维数组、一维数组和求解器实际网格数组。

## 8. 已知风险与维护规则

1. `Fluids/NaK78.py` 仅部分实现，不能替代 `SodiumPotassium78`。
2. `SodiumPotassium78` 是液态 NaK78 主用类，但仍带有钠命名和钠式两相近似；气相关联式不得直接复用。
3. `FluidMaterial` 兼容层会静默回退到通用方法；新增流体时要逐项确认回退是否符合物理语义。
4. `Potassium` 中存在明确返回 `0.0` 的占位接口。
5. `SodiumHP`、`PotassiumHP` 的 `Cp` 在 2 K 糊状区内陡增，会直接影响刚性求解器。
6. `WickMaterial` 的导热上限、查表端点钳位、高非线性刷新区和缓存失效规则都会影响热管数值稳定性。
7. `GasGaps.py` 的材料是等效固体，不应传入水力网络。
8. 多个经验式在边界外继续外推，或通过 `np.clip`、`np.interp` 钳位；修改边界前必须检查求解器工作温区。
9. Materials 源码变更后必须同步更新本文，并用仓库级搜索重新核对调用路径。

## 9. 最小验证清单

修改 Materials 后至少检查：

1. 构造函数、抽象方法和兼容层是否仍可调用。
2. 单位是否保持为 SI：`K`、`Pa`、`kg/m^3`、`Pa*s`、`W/(m*K)`、`J/(kg*K)`、`J/kg`。
3. 分段边界两侧、熔点糊状区两侧、临界温度附近和非正输入保护是否符合预期。
4. 标量与数组输入是否都能工作。
5. `WickMaterial` 参数变化后查表是否失效并重建。
6. `Solvers/`、`Components/`、`Correlations/` 的实际调用点是否与本文一致。
