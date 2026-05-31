# Correlations - AI Agent 速查手册

> 面向后续 AI Agent 的模块级索引。优先阅读本文件，再按需定位 `Correlations.py` 中的具体函数；通常无需重新通读整份源码。

## 0. 阅读指引

### 0.1 模块职责

`Correlations/` 是 TASTIN 的经验关联式库，负责：

- 单相摩阻系数和两相摩阻倍增因子。
- 绕丝棒束摩擦压降：6 种模型及统一分派器。
- 空泡份额：6 种模型。
- 圆管单相/两相摩擦压降。
- 环管、棒束、Aoki、IHX 和横掠单管换热关联式。

### 0.2 源码规模与建议阅读顺序

本目录只有 1 个 Python 源文件：

| 文件 | 规模 | 用途 |
|------|------|------|
| `Correlations.py` | 1143 行 | 全部经验关联式、`BundleGeometry` 和内部辅助函数 |

建议阅读顺序：

1. 先读本文 §1 API 总表和 §6 真实调用关系。
2. 需要接入压降模型时，再读 §2、§3、§5。
3. 需要接入换热模型时，再读 §4。
4. 仅在修改公式、排查数值异常或扩展材料对象时定位源码实现。

### 0.3 适用范围

本文只描述当前 `Correlations` 模块，不是整个 TASTIN-python 项目总览。公式说明以当前 Python 实现为准；已知但未修复的问题单列在 §7。

---

## 1. API 总表

### 1.1 数据结构和压降 API

| # | 接口 | 类别 | Fortran 对应 | 输入摘要 | 输出 | 备注 |
|---|------|------|-------------|----------|------|------|
| 1 | `BundleGeometry` | 数据结构 | 分散的几何参数 | 棒束整体和 3 类子通道几何 | dataclass | 见 §2 |
| 2 | `friction_single_phase(Re)` | 单相摩阻 | `FRICTION(RE)` | `Re` | 达西摩阻因子 `f` | 圆管路径和外部水力元件使用 |
| 3 | `failo2(method_idx, P_sys, x, G, material)` | 两相倍增 | `FAILO2(N, PP, XX, GG)` | 方法 1-4、压力、干度、质量流速、材料 | `phi^2_LO` | 棒束和圆管两相路径硬编码调用方法 4 |
| 4 | `_get_bundle_state(P_fluid, H_fluid, G, material)` | 内部辅助 | - | 压力、焓、质量流速、材料 | `T_fluid, x, rho, mu` | 仅供棒束压降模型使用 |
| 5 | `friction_novendstern(P_fluid, H_fluid, G, length, geom, material)` | 棒束压降 | `NOVENDSTERN` | 基础状态、长度、几何、材料 | `delta_P [Pa]` | 模型 1 |
| 6 | `friction_rehme(P_fluid, H_fluid, G, length, geom, nr, material)` | 棒束压降 | `REHME` | 同上，另含绕丝圈数 `nr` | `delta_P [Pa]` | 模型 2 |
| 7 | `friction_cheng(P_fluid, H_fluid, G, length, geom, material)` | 棒束压降 | `CHENG` | 基础状态、长度、几何、材料 | `delta_P [Pa]` | 模型 3 |
| 8 | `friction_engel(P_fluid, H_fluid, G, length, geom, material)` | 棒束压降 | `ENGEL` | 基础状态、长度、几何、材料 | `delta_P [Pa]` | 模型 4 |
| 9 | `friction_sobolev(P_fluid, H_fluid, G, length, geom, material)` | 棒束压降 | `SOBOLEV` | 基础状态、长度、几何、材料 | `delta_P [Pa]` | 模型 5 |
| 10 | `friction_crt(P_fluid, H_fluid, G, length, geom, material)` | 棒束压降 | `CRT` | 基础状态、长度、几何、材料 | `delta_P [Pa]` | 模型 6，含 Python 版修正 |
| 11 | `pressure_drop_bundle_dispatcher(method_flag, P_fluid, H_fluid, G, length, geom, nr, material)` | 分派器 | `PRESSUREDPFBUNDLE` | `method_flag=1..6` 和棒束参数 | `delta_P [Pa]` | 非法 flag 返回 `0.0` |
| 12 | `void_fraction(method_idx, P_fluid, x, G, D_hyd, material)` | 空泡份额 | `VOID(N, PP, X, G, DE)` | `method_idx=1..6`、压力、干度、质量流速、水力直径、材料 | `alpha` | 边界行为见 §3.4 和 §7 |
| 13 | `pressure_drop_tube(P_fluid, H_fluid, length, D_hyd, G, material)` | 圆管压降 | `PRESSUREDPFTUBE` | 压力、焓、长度、水力直径、质量流速、材料 | `delta_P [Pa]` | 单相和两相统一入口 |

### 1.2 换热 API

| # | 接口 | 场景 | Fortran 对应 | 输入摘要 | 输出 | 备注 |
|---|------|------|-------------|----------|------|------|
| 14 | `nu_ringpipe(R_out, R_in, Re, Pr)` | 环形通道 | `NU_RINGPIPE` | 内外半径或直径、`Re`、`Pr` | `Nu` | 主程序当前使用 |
| 15 | `nu_fftf(P, D, Pe)` | FFTF 棒束 | `NU_FFTF` | 棒间距、棒径、`Pe` | `Nu` | 推荐范围：`P/D=1.15..1.3`，`Pe=10..1500` |
| 16 | `nu_calamai(P, D, Pe)` | Calamai 棒束 | `NU_CALAMAI` | 棒间距、棒径、`Pe` | `Nu` | 当前公式与 FFTF 相同 |
| 17 | `nu_fuel_bundle_dispatcher(method_idx, P, D, Pe)` | 棒束换热分派器 | `NU_FUELBUNDLE` | `method_idx=1..2` 和棒束参数 | `Nu` | 非法 method 回退到 FFTF |
| 18 | `nu_aoki(Re, Pr)` | Aoki 液态金属管内换热 | `H_AOKI` 核心 | 标量或 NumPy 数组 | `Nu` | 支持广播和向量化 |
| 19 | `h_aoki(Re, Pr, D_hyd, T_fluid, material)` | Aoki 换热系数 | `H_AOKI` | `Re`、`Pr`、水力直径、温度、材料 | `h [W/(m^2*K)]` | `h=Nu*k/D` |
| 20 | `nu_ihx_tube(Pe)` | IHX 管束 | `NU_IHXTube` | `Pe` | `Nu` | 中间热交换器 |
| 21 | `nu_single_crossflow_pipe(Re_channel, Pr, D_channel, D_out, A_flow, A_min)` | 横掠单根热管 | Python 新增 | 名义 `Re`、流道和热管尺寸、流通面积 | `Nu` | 内部换算最小截面处 `Re_max` |
| 22 | `h_single_crossflow_pipe(Re_max, Pr, D_out, T_fluid, material)` | 横掠单根热管换热系数 | Python 新增 | `Re_max`、`Pr`、热管外径、温度、材料 | `h [W/(m^2*K)]` | 当前实现存在运行时错误，见 §7 |

`Numeric = Union[float, np.ndarray]` 仅用于 `nu_aoki()` 的向量化类型标注。

---

## 2. BundleGeometry

```python
@dataclass
class BundleGeometry:
    P: float       # 棒间距 Pitch [m]
    D: float       # 棒直径 Diameter [m]
    PH: float      # 绕丝螺距 Wire Pitch / Lead [m]
    DW: float      # 绕丝直径 Wire Diameter [m]
    PWT: float     # 总湿周 Wetted Perimeter [m]

    N1: int        # 内部子通道数量
    N2: int        # 边缘子通道数量
    N3: int        # 角部子通道数量
    A1: float      # 内部子通道面积 [m^2]
    A2: float      # 边缘子通道面积 [m^2]
    A3: float      # 角部子通道面积 [m^2]
    DE1: float     # 内部子通道水力直径 [m]
    DE2: float     # 边缘子通道水力直径 [m]
    DE3: float     # 角部子通道水力直径 [m]
```

源码注释将三类子通道记为 `1: Inner, 2: Edge, 3: Corner`。各棒束模型使用同一数据结构，但几何组合方式不同。

---

## 3. 压降和空泡份额

### 3.1 单相摩阻 `friction_single_phase(Re)`

```text
Re <= 1000                 -> f = 64 / Re
2300 < Re < 100000         -> f = 0.3164 / Re^0.25
1000 < Re <= 2300          -> Karman-Prandtl 隐式迭代
Re >= 100000               -> Karman-Prandtl 隐式迭代
```

函数先执行 `Re_calc = max(Re, 1e-5)`，再进行分段计算。它只被圆管压降和外部水力元件调用；6 个棒束模型使用各自摩阻公式，不直接调用它。

### 3.2 两相倍增因子 `failo2(...)`

| method | 名称 | 当前实现要点 |
|--------|------|--------------|
| 1 | Chisholm | `G <= 2000` 和 `G > 2000` 两个分支；`x <= 0` 或 `x >= 1` 时直接返回 `1.0` |
| 2 | USSR Standard | `1 + x * (v_g / v_l - 1)` |
| 3 | Martinelli-Nelson | 使用滑移比和 `X_tt`；边界干度存在风险 |
| 4 | Homo-Flow | `phi^2 = (1 + x*(v_g-v_l)/v_l) * (1 + x*(mu_l-mu_g)/mu_g)^(-0.25)` |

所有棒束模型和 `pressure_drop_tube()` 的两相分支都硬编码调用 `failo2(4, ...)`。

### 3.3 六种棒束压降模型

| flag | 模型 | 特征尺寸 | 摩阻公式特征 | 两相处理 |
|------|------|----------|--------------|----------|
| 1 | Novendstern | `DE1` | 中心通道流量分配，`M * 0.3164 / Re1^0.25` | `0 <= x <= 1` 时乘 `failo2(4, ...)` |
| 2 | Rehme | 总体 `DE = 4*A_total/PWT` | 绕丝几何因子 `F`，另乘 `nr*pi*(D+DW)/PWT` | `0 <= x <= 1` 时乘 `failo2(4, ...)` |
| 3 | Cheng & Todreas | 总体 `DE` | 层流、过渡和湍流三段式 | `0 <= x <= 1` 时乘 `failo2(4, ...)` |
| 4 | Engel | 总体 `DE` | `Re < 400`、`400 <= Re < 5000`、`Re >= 5000` 三段式 | `0.001 <= x <= 1` 时乘 `failo2(4, ...)` |
| 5 | Sobolev | 总体 `DE` | `f = term1 * term2`，含 `D/PH` 和 `P/D` | `0 <= x < 1` 时乘 `failo2(4, ...)` |
| 6 | CRT | `DE1` 分配流量，总体 `DE` 计算压降 | Chiu-Rohsenow-Todreas 子通道模型 | `0 <= x <= 1` 时乘 `failo2(4, ...)` |

统一入口：

```python
pressure_drop_bundle_dispatcher(
    method_flag, P_fluid, H_fluid, G, length, geom, nr, material
)
```

### 3.4 空泡份额 `void_fraction(...)`

| method | 名称 | 当前实现要点 |
|--------|------|--------------|
| 1 | Osimaqkih | Froude 类参数、体积干度 `beta` 和滑移修正 |
| 2 | Bankoff | 迭代求解，`K = 0.71 + 0.0145*P_MPa` |
| 3 | Bankoff-Jones | 迭代求解，`K_BJ = K_B + (1-K_B)*alpha^R` |
| 4 | Smith | 滑移比显式计算 |
| 5 | Chisholm | 滑移比显式计算 |
| 6 | Homogeneous | 均相流，公式等价于滑移比为 `1` |

边界行为必须准确区分：

- `x <= 0.0`：进入模型计算前直接返回 `0.0`。
- `x >= 1.0`：进入模型计算前直接返回 `1.0`。
- `0.0 < x < 1.0`：模型计算完成后截断到 `[0.0, 0.99]`。
- 非法 `method_idx`：返回 `0.0`。

### 3.5 圆管压降 `pressure_drop_tube(...)`

```text
先计算 x = (H - h_f) / (h_g - h_f)

x < 0       -> 液相物性 -> friction_single_phase(Re) -> delta_P
0 <= x <= 1 -> 饱和液物性 -> friction_single_phase(Re_lo)
               -> delta_P_lo * failo2(4, ...)
x > 1       -> 气相物性 -> friction_single_phase(Re) -> delta_P
```

---

## 4. 换热关联式

| 接口 | 场景 | 当前公式摘要 |
|------|------|--------------|
| `nu_ringpipe()` | 环形通道 | `Nu = AA + BB*(PF*Re*Pr)^gamma`，几何比为 `R_out/R_in` |
| `nu_fftf()` | FFTF 棒束 | `4 + 0.16*(P/D)^5 + 0.33*(P/D)^3.8*(Pe/100)^0.86` |
| `nu_calamai()` | Calamai 棒束 | 当前与 FFTF 完全相同，保留独立入口 |
| `nu_aoki()` | Aoki 液态金属管内换热 | `Re < 3000` 返回 `4.36`；否则返回湍流公式 |
| `h_aoki()` | Aoki 换热系数 | `nu_aoki(Re, Pr) * material.conductivity(T_fluid) / D_hyd` |
| `nu_ihx_tube()` | IHX 管束 | `4.5 + 0.014*Pe^0.8` |
| `nu_single_crossflow_pipe()` | 横掠单根热管 | `Pr < 0.1` 使用液态金属分支；否则使用 Churchill-Bernstein |
| `h_single_crossflow_pipe()` | 横掠单根热管换热系数 | 设计目标为 `Nu*k/D_out`，但当前实现不可正常调用，见 §7 |

---

## 5. material 对象契约

### 5.1 汇总

按调用路径组合，`material` 对象可能需要提供以下方法：

```python
# 热力状态
material.saturation_temperature(P)          # -> T_sat [K]
material.temperature_from_enthalpy(H, P)    # -> T [K]

# 饱和焓：源码中存在两套命名
material.liquid_enthalpy_sat(T_sat)         # -> h_f [J/kg]
material.vapor_enthalpy_sat(T_sat)          # -> h_g [J/kg]
material.enthalpy_saturated_liquid(T_sat)   # -> h_f [J/kg]
material.enthalpy_saturated_vapor(T_sat)    # -> h_g [J/kg]

# 密度
material.density(T, P)                      # -> rho [kg/m^3]
material.liquid_density(T, P)               # -> rho_l [kg/m^3]
material.vapor_density(T, P)                # -> rho_g [kg/m^3]

# 动力黏度
material.liquid_viscosity(T, P)             # -> mu_l [Pa*s]
material.vapor_viscosity(T, P)              # -> mu_g [Pa*s]

# 导热系数
material.conductivity(T)                    # -> k [W/(m*K)]
```

### 5.2 两套饱和焓接口不可混淆

| 调用路径 | 实际调用的方法 |
|----------|----------------|
| `_get_bundle_state()`，即所有棒束压降模型 | `liquid_enthalpy_sat()` / `vapor_enthalpy_sat()` |
| `pressure_drop_tube()` | `enthalpy_saturated_liquid()` / `enthalpy_saturated_vapor()` |

`Materials/Base.py` 中的 `FluidMaterial` 已提供前一组兼容包装，并可转调后一组名称。但独立实现或自定义材料对象如果不继承该基类，仍必须兼容实际使用路径需要的接口。

### 5.3 按功能拆分

| 功能 | 依赖的方法 |
|------|------------|
| `failo2()` | `saturation_temperature()`、液/气相密度、液/气相黏度 |
| `_get_bundle_state()` | `temperature_from_enthalpy()`、`saturation_temperature()`、`liquid_enthalpy_sat()`、`vapor_enthalpy_sat()`、`density()`、按相态选择的黏度 |
| `pressure_drop_tube()` | `saturation_temperature()`、`enthalpy_saturated_liquid()`、`enthalpy_saturated_vapor()`、`temperature_from_enthalpy()`、按相态选择的密度和黏度 |
| `h_aoki()`、`h_single_crossflow_pipe()` | `conductivity()` |

压力入口统一按 `Pa` 传入；`void_fraction()` 内部会将压力转为 `MPa` 供经验公式使用。

---

## 6. 真实调用关系

### 6.1 主程序当前使用

对上级目录 Python 调用点检索后，主程序实际使用两个入口：

| 接口 | 调用位置 | 用途 |
|------|----------|------|
| `friction_single_phase()` | `Solvers/Hydrodynamics/Components.py` | 水力元件 Darcy-Weisbach 摩擦压降 |
| `nu_ringpipe()` | `Components/TFEUnit.py` | TFE 环形流道换热，经过局部 adapter 传入耦合组件 |

`Solvers/Hydrodynamics/BoundaryVolume.py` 也导入了 `friction_single_phase()`，但当前未发现实际调用。

### 6.2 测试中使用或导入

| 接口 | 测试位置 | 状态 |
|------|----------|------|
| `nu_ringpipe()` | `testModule/test_component_annular_pipe.py` | adapter 中实际调用 |
| `nu_aoki()` | `test_flow_heat.py`、`test_fluid_solid_couple.py`、`test_pipe_heat_transfer.py` 等 | 多处导入，部分测试中实际调用 |
| `nu_fftf()` | `testModule/test_component_pipe.py` | 当前导入，未发现函数调用 |

### 6.3 当前未发现外部调用，但保留为库能力

以下接口在当前上级目录 Python 代码中未发现外部调用，仍应作为关联式库能力保留：

- `failo2()`、`pressure_drop_bundle_dispatcher()` 和 6 个棒束压降模型。
- `void_fraction()`、`pressure_drop_tube()`。
- `nu_calamai()`、`nu_fuel_bundle_dispatcher()`。
- `h_aoki()`、`nu_ihx_tube()`。
- `nu_single_crossflow_pipe()`、`h_single_crossflow_pipe()`。

### 6.4 内部关系图

```text
pressure_drop_bundle_dispatcher(flag=1..6, ...)
|-- friction_novendstern()
|-- friction_rehme()
|-- friction_cheng()
|-- friction_engel()
|-- friction_sobolev()
`-- friction_crt()
    共同路径：_get_bundle_state()；两相时 failo2(4, ...)

pressure_drop_tube(...)
|-- friction_single_phase(Re)        # 单相和 liquid-only 基础压降
`-- failo2(4, ...)                   # 两相倍增

nu_fuel_bundle_dispatcher(method=1..2, ...)
|-- nu_fftf()
`-- nu_calamai()

h_aoki() -> nu_aoki() -> h = Nu*k/D

h_single_crossflow_pipe() -> nu_single_crossflow_pipe()
# 当前调用参数不足，见 §7
```

---

## 7. 已知风险，不在本次修改范围内

| # | 风险 | 影响 |
|---|------|------|
| 1 | `h_single_crossflow_pipe()` 只向 `nu_single_crossflow_pipe()` 传入 `Re_max, Pr` 两个参数，但后者要求 6 个参数。 | 调用 `h_single_crossflow_pipe()` 会触发运行时 `TypeError`。 |
| 2 | 棒束和圆管压降依赖的饱和焓接口命名不同。 | 材料实现必须兼容 §5.2 中的实际路径，否则部分压降函数运行失败。 |
| 3 | `void_fraction()` 对 `x >= 1.0` 在模型计算前直接返回 `1.0`；只有内部计算分支结果才会截断至 `0.99`。 | 不要将函数整体输出范围误写成 `[0.0, 0.99]`。 |
| 4 | `failo2()` 的部分非默认公开分支存在边界输入风险，例如方法 3 在边界干度附近含除法。 | 当前棒束和圆管两相压降硬编码使用 `failo2(4, ...)`；直接调用其他方法时需要额外验证输入。 |
| 5 | `friction_crt()` 对原 Fortran 的子通道流量分配公式做了 `DE2 -> DE3` 修正。 | 这是有意保留的物理逻辑修正，后续迁移或对照 Fortran 时不要误删。 |

---

## 8. 快速查找表

| 需求 | 首选接口 |
|------|----------|
| 光管单相摩阻因子 | `friction_single_phase(Re)` |
| 两相摩阻倍增因子 | `failo2(method_idx, P_sys, x, G, material)` |
| 棒束压降，6 种模型任选 | `pressure_drop_bundle_dispatcher(method_flag, P_fluid, H_fluid, G, length, geom, nr, material)` |
| 圆管单相/两相压降 | `pressure_drop_tube(P_fluid, H_fluid, length, D_hyd, G, material)` |
| 空泡份额 | `void_fraction(method_idx, P_fluid, x, G, D_hyd, material)` |
| 环管换热 `Nu` | `nu_ringpipe(R_out, R_in, Re, Pr)` |
| 棒束换热 `Nu` | `nu_fuel_bundle_dispatcher(method_idx, P, D, Pe)` |
| Aoki 液态金属换热 `Nu` / `h` | `nu_aoki(Re, Pr)` / `h_aoki(Re, Pr, D_hyd, T_fluid, material)` |
| IHX 管束换热 `Nu` | `nu_ihx_tube(Pe)` |
| 横掠单根热管换热 `Nu` | `nu_single_crossflow_pipe(Re_channel, Pr, D_channel, D_out, A_flow, A_min)` |
| 横掠单根热管换热 `h` | `h_single_crossflow_pipe(...)`，使用前先修复 §7 风险 1 |

