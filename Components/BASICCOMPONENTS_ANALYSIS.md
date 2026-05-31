# Components/basicComponents 基础物理元件分析

## 1. 概述

`Components/basicComponents/` 是 TASTIN 组件层中最底层的物理元件集合。它们直接继承自 `Solvers/` 层的求解器基类，封装了单一物理实体（燃料芯块、电极、热管、翅片、电极对）的几何建模、物性关联和数值求解方法。

## 文件清单

| 文件 | 类 | 继承自 | 职责 |
|---|---|---|---|
| `Fuel.py` | `Fuel` | `HeatConduction2D` | 燃料芯块导热 + 功率分配 + 反应性反馈 |
| `Electord.py` | `Electrode` / `Emitter` / `Collector` | `HeatConduction2D` | 电极导热 + 焦耳热注入 + 反应性反馈 |
| `HeatPipe2D.py` | `HeatPipe2D` | `HeatConduction2D` | 热管 (吸液芯/管壁) 二维导热 + 分段外边界 |
| `FinConduction.py` | `FinConduction` | `HeatConduction2D` | 翅片准稳态导热 (稳态稀疏求解替代 ODE) |
| `TECPair.py` | `TECPair` | 独立类 | 发射极-接收极配对 + 间隙耦合 + 等离子体热流映射 |

---

## 2. Fuel 燃料芯块

**文件：** `Fuel.py`
**继承：** `HeatConduction2D`

### 2.1 物理定位

燃料芯块是整个系统中裂变能量的**最终热源接收者**。它继承了 2D 导热求解器的所有能力，并在其基础上扩展了：

1. **功率空间分配**：将点堆产生的标量总功率按轴向和径向分布分配到每个网格节点
2. **反应性反馈输出**：提供体积加权平均温度 → 反馈反应性多项式
3. **轴向接触热阻**：在指定的轴向交界面串联额外的热阻

### 2.2 构造函数参数

```python
Fuel(name, mesh, material=None, initial_temp=600.0,
     power_fraction=1.0, axial_power_profile=None,
     contact_resistance_interfaces=None, axial_contact_resistance=0.0)
```

| 参数 | 默认值 | 含义 |
|---|---|---|
| `name` | — | 组件唯一名称 |
| `mesh` | — | Mesh2D 圆柱网格 (r-z) |
| `material` | `UO2()` | 燃料材料 (默认二氧化铀) |
| `initial_temp` | 600.0 | 初始均匀温度 [K] |
| `power_fraction` | 1.0 | 该燃料元件占总堆功率的份额 |
| `axial_power_profile` | 均匀分布 | 轴向功率分布数组 [ny] |
| `contact_resistance_interfaces` | None | 存在接触热阻的轴向交界面索引列表 |
| `axial_contact_resistance` | 0.0 | 接触界面总热阻 [K/W] |

### 2.3 核心方法

#### `set_axial_power_profile(profile_array)`

预计算功率分配权重矩阵，流程：
```
1. 入参归一化使 Σprofile = 1.0
2. 对每个轴向层 j:
     weights_2d[:, j] = profile[j] * (vols_2d[:, j] / vols_layer[j])
     即：轴向分配 × 径向按体积分
3. 展平为 self.power_allocation_weights [N]
```

#### `set_nuclear_power(p_fiss, p_decay, p_total)`

```
Q_source = p_total × power_fraction × power_allocation_weights
Q_vol = Q_source / vols_flat
```

关键设计：使用预计算的权重矩阵，避免每次调用时循环计算，实现 O(1) 向量化分配。

#### `get_reactivity_feedback() -> float`

```
T_avg = Σ(T_i × vol_i) / total_vol   # 体积加权平均温度
ρ_fb = 0.001360811 - 6.47927757e-6×T + 2.321231e-9×T² - 3.52e-13×T³
```

#### 拦截方法：`_compute_internal_resistance()`

在父类计算的界面热导基础上，对指定轴向交界面串联接触热阻：
```python
R_contact_ring = r_c / A_ring          # 在截面上分摊
G_y_inner[:, j] = 1 / (R_orig + R_contact_ring)
```

### 2.4 调用链

```
SystemManager → ReactorCore
  → Fuel.set_nuclear_power(fission, decay, total)
  → Fuel.step(dt)  (继承自 HeatConduction2D)
  → Fuel.get_reactivity_feedback()  → 返回反馈反应性
```

---

## 3. 电极组件 (Electrode / Emitter / Collector)

**文件：** `Electord.py`
**继承：** `Electrode` → `HeatConduction2D`，`Emitter`/`Collector` → `Electrode`

### 3.1 物理定位

电极是热离子燃料元件的核心部分：
- **Emitter（发射极）**：靠近燃料侧，高温发射电子
- **Collector（接收极）**：靠近冷却剂侧，接收电子

两者都承受**焦耳热**（电流通过产生），并向点堆提供温度反应性反馈。

### 3.2 继承体系

```
HeatConduction2D
    └── Electrode (基类)
          ├── Emitter  (材料默认: MoNb 钼铌合金)
          └── Collector (材料默认: Molybdenum 钼)
```

### 3.3 Electrode 基类

**构造函数参数：** `Electrode(name, mesh, material, initial_temp)`

**核心方法：**

#### `set_joule_heating(q_joule_array: np.ndarray)`

TEC 计算完成后，接收展平的焦耳热功率分布 [W/node]：

```python
Q_source[:] = q_joule_array                       # 原地赋值
use_external_source_buffer = True                  # 激活外部热源标记
Q_vol = Q_source / vols_flat                      # 体积热源 [W/m³]
```

关键设计：设置 `use_external_source_buffer = True` 告诉基类 `BaseHeatConduction`：**这个热源由外部显式维护，不要在 `_update_sources` 中将其清零！**

#### `get_reactivity_feedback()` — 在基类中抛 NotImplementedError，强制子类实现

### 3.4 Emitter 与 Collector

两者结构完全对称，仅材料默认值和反应性反馈多项式不同。

**反应性反馈多项式：**

```
T_avg = Σ(T_i × vol_i) / total_vol
ρ_fb = 1e-4 × (3.46455e-6 × T² - 0.03232167 × T + 0.74202216)
```

> ⚠️ 注意：Collector 当前使用的反馈多项式与 Emitter 相同，代码中有 `TODO` 标记提醒修改。

### 3.5 调用链

```
TECCircuitManager.pre_step()
  → TECPair.set_joule_heating(...)
      → emitter.set_joule_heating(q_e_array)
      → collector.set_joule_heating(q_c_array)

ReactorCore
  → emitter.get_reactivity_feedback()
  → collector.get_reactivity_feedback()
```

---

## 4. HeatPipe2D 热管二维求解器

**文件：** `HeatPipe2D.py`
**继承：** `HeatConduction2D`

### 4.1 物理定位

热管是空间核电源**废热排放**的核心部件。TASTIN 的热管模型基于**纯导热假设**（相变换热等效为极高轴向热导率），采用 2D 圆柱网格 (r-z)。

### 4.2 网格结构

```
径向分层 (r direction):
  ┌─────────────────────────────┐
  │  vapor core (虚边界, 绝热)   │  r = 0 .. n_wick
  │  吸液芯区 (wick)              │  n_wick 个径向节点
  │  管壁区 (wall)                │  n_x - n_wick 个径向节点
  └─────────────────────────────┘

轴向分段 (z direction):
  [ 蒸发段 (eva) ][ 绝热段 (aba) ][ 冷凝段 (con) ]
   n_eva 节点       n_aba 节点       n_con 节点
```

### 4.3 构造函数参数

```python
HeatPipe2D(mesh, solid1, solid2, solid3, n_wick, porosity,
           n_eva, n_aba, n_con, name, emissivity,
           up_view_factor, down_view_factor, initial_temp)
```

| 参数 | 含义 |
|---|---|
| `solid1` | 管壁材料 |
| `solid2` | 热管工质 (液态金属) |
| `solid3` | 吸液芯骨架材料 |
| `n_wick` | 吸液芯区径向节点数 |
| `porosity` | 吸液芯孔隙率 |
| `n_eva / n_aba / n_con` | 蒸发段/绝热段/冷凝段轴向节点数 |
| `emissivity` | 外壁面发射率 |
| `up_view_factor / down_view_factor` | 上下表面角系数 |

### 4.4 关键设计特点

#### 吸液芯复合物性 (`WickMaterial`)

吸液芯不是单一材料，而是通过 `WickMaterial` 类将**骨架固体 + 液态工质**按孔隙率复合：
```
k_wick = f(k_solid, k_fluid, porosity)
ρ_wick = ρ_solid×(1-porosity) + ρ_fluid×porosity
cp_wick = (cp_solid×ρ_solid×(1-porosity) + cp_fluid×ρ_fluid×porosity) / ρ_wick
```

支持各向异性导热 (`use_anisotropic_wick_conductivity`)：轴向和径向使用不同的导热系数模型。

#### 分段外边界

热管外壁不再使用统一的 `right` 边界，而是按轴向段拆分为三个独立边界：
```python
self.boundaries = {
    'left':    ...,   # 内壁 (蒸汽腔, 绝热)
    'bottom':  ...,   # 轴向底
    'top':     ...,   # 轴向顶
    'outer_eva': BoundaryRegion(shape=(n_eva,)),  # 蒸发段外壁
    'outer_aba': BoundaryRegion(shape=(n_aba,)),  # 绝热段外壁
    'outer_con': BoundaryRegion(shape=(n_con,)),  # 冷凝段外壁
}
```

原始 `right` 边界在 `_setup_virtual_boundaries()` 中被删除。

#### 惰性物性更新 (`_update_properties`)

物性计算采用**增量更新策略**：
- 记录上次更新时的温度缓存 (`_wick_temperature_cache`, `_wall_temperature_cache`)
- 仅对温度变化超过 `_property_update_tol` (5%) 或进入非线性区的节点重新计算物性
- 大幅降低高频调用时的计算开销

#### 界面热导计算模式 (`face_conductance_mode`)

支持四种模式：
| 模式 | 说明 |
|---|---|
| `legacy_harmonic` | 调和平均 (默认，与 HeatConduction2D 一致) |
| `resistance_split_axial` | 轴向电阻分裂法 |
| `resistance_split_xy` | 径向和轴向均电阻分裂 |
| `resistance_split_full` | 完全电阻分裂 (圆柱坐标使用对数热阻) |

#### 冻结物性校正 (`step` 方法)

重写了 `step()`，支持 **frozen-property correction** 模式：
```
1. 用当前温度更新物性、热阻、边界
2. 冻结物性，用 BDF 推进 dt
3. 用推进后的温度重新更新物性，检查温度变化
4. 若变化 > outer_property_tol → 重复 correction 循环
5. 最多执行 max_outer_property_corrections 次
```

### 4.5 调用链

```
HPwithFin.__init__()
  → HeatPipe2D(name, mesh, wall_mat, fluid_mat, wick_mat, ...)

HPwithFin.pre_step()
  → hp.step(dt)  → BDF 推进温度场
```

---

## 5. FinConduction 翅片准稳态求解器

**文件：** `FinConduction.py`
**继承：** `HeatConduction2D`

### 5.1 物理定位

翅片是热管冷凝段外挂的扩展散热面。由于翅片极薄（热容≈0），采用准稳态假设：每个全局时间步内翅片瞬时达到热平衡 (dT/dt=0)，无需 ODE 时间推进。

### 5.2 降维策略

1. **忽略轴向 (Y) 导热**：将 2D 网格退化为 Ny 个独立的 1D 传热条
2. **辐射体热源化**：将表面辐射等效为体积热源/热漏

### 5.3 构造函数参数

```python
FinConduction(mesh, material, fin_thickness, emissivity=0.8,
              up_view_factor=1.0, down_view_factor=1.0,
              T_env=3.0, initial_temp=298.15)
```

| 参数 | 含义 |
|---|---|
| `fin_thickness` | 翅片厚度 [m] |
| `emissivity` | 表面发射率 |
| `up_view_factor` | 上表面角系数 |
| `down_view_factor` | 下表面角系数 |
| `T_env` | 环境辐射温度 [K] (太空 ≈ 3K) |

### 5.4 核心方法：`step(dt)` — 重写 ODE 求解

**流程：**

```
1. 更新物性 + 热阻 + 边界状态

2. 组装稀疏线性系统 A*T = b:

   内部X方向导热:
     k1 ← 左节点, k2 ← 右节点
     D[k1] += G_x, D[k2] += G_x    (主对角, 流出为正)
     A[k1,k2] = -G_x, A[k2,k1] = -G_x  (非对角, 流入为负)

   [不添加 Y 方向导热 → 自动解耦为 Ny 个独立块！]

   表面辐射散热:
     h_rad = ε×σ×(T²+T_env²)×(T+T_env)  ← Picard 线性化
     A_eff = (F_up + F_down) × vol/thickness
     D += h_rad × A_eff
     b += h_rad × A_eff × T_env

   四条边界条件:
     left (根部 → 接热管):  戴维南等效 → G_bound, S_bound
     right (远端边缘):      同上
     bottom / top (轴向端): 同上

3. 求解 spsolve(A, b)

4. 松弛更新: T_new = α×T_solved + (1-α)×T_old

5. 迭代至收敛 (Picard 迭代, 处理辐射非线性)
```

**关键设计点：**
- 直接使用 `scipy.sparse.linalg.spsolve` 求解稳态，而非 ODE 积分
- 使用 `np.add.at` 进行稀疏矩阵组装（避免 Python 循环）
- 辐射非线性通过 Picard 迭代 + 松弛因子处理
- 边界条件通过 `_get_boundary_linear_terms()` 从 `BoundaryRegion` 的多个叠加条件中提取戴维南等效 `(G_bound, S_bound)`

### 5.5 调用链

```
HPwithFin.pre_step()
  → FinConduction.step(dt)  → 直接求解稳态稀疏系统
```

---

## 6. TECPair 热离子电极对

**文件：** `TECPair.py`
**继承：** 独立类

### 6.1 物理定位

`TECPair` 封装了一对发射极-接收极，是 `TFEUnit` 和 `TECCircuitManager` 之间**承上启下**的关键桥梁。

### 6.2 内部对象

```
TECPair
  ├── emitter (Emitter)         ← 继承 Electrode → HeatConduction2D
  ├── collector (Collector)     ← 继承 Electrode → HeatConduction2D
  ├── tec_gap (TECCouple2D)    ← 极间隙耦合器 (导热+辐射)
  ├── inner_boundary           ← emitter.boundaries['left']
  └── outer_boundary           ← collector.boundaries['right']
```

### 6.3 构造函数参数

```python
TECPair(name, L_node, n_node, R_e_in, delta_e, delta_gap, delta_c,
        n_rad_e, n_rad_c, mat_emitter, mat_collector,
        T_init_e=1600, T_init_c=800,
        k_gap_gas=0.0, emissivity_e=0.8, emissivity_c=0.8)
```

**自动推算的几何尺寸：**
```
R_e_out = R_e_in + delta_e              # 发射极外径
R_c_in = R_e_out + delta_gap            # 接收极内径 (间隙外边界)
R_c_out = R_c_in + delta_c              # 接收极外径
```

### 6.4 核心方法

#### `get_gap_surface_temperatures() -> (T_emit_surf, T_coll_surf)`

提取间隙两侧表面温度，供 `TECCircuitManager` 传入 C++ 电路求解器。

#### `set_joule_heating(dU_emit, rho_emit, dU_coll, rho_coll, alpha)` / `set_joule_heating_fields(E_emit, rho_emit, E_coll, rho_coll, alpha)`

接收 C++ 电路求解器返回的电场数据，转换为焦耳热并下发。

```
Q_vol = E² / ρ             # 体积焦耳热率 [W/m³]
E = dU / L_node            # 电场强度 [V/m]

q_watts_flat = E² / ρ × vols_flat  # 每个节点的焦耳热功率 [W]

支持亚松弛: q_new = α×q_calc + (1-α)×q_old
```

调用 `joule_power_from_electric_field()` 辅助函数 (`Components/tec_electric.py`)。

#### `update_plasma_heat_flux(q_e_flux, q_c_flux, alpha)`

接收 C++ 电路求解器返回的等离子体放电热流密度 [W/m²]，转换为功率 [W] 并施加到间隙耦合器。

```
Q_e_watts = q_e_new × A_emit_surf      # 发射极电子冷却功率
Q_c_watts = q_c_new × A_emit_surf      # 接收极电子加热功率 (基于发射极面积)

tec_gap.set_tec_sources(Q_emitter=Q_e_watts, Q_collector=Q_c_watts)
```

**面积换算说明：** 电流密度 J [A/cm²] 基于发射极外表面积计算，因此等离子体热流也统一基于发射极面积折算，保证能量守恒。

#### `sync_and_step(dt)`

```python
tec_gap.sync()            # 间隙耦合器双向同步
emitter.step(dt)          # BDF 推进发射极温度
collector.step(dt)        # BDF 推进接收极温度
```

### 6.5 调用链

```
TFEUnit.__init__()
  → TECPair(name, L_node, n_node, R_e_in, ...)

TECCircuitManager.pre_step()
  → TECPair.get_gap_surface_temperatures()  → 收集温度
  → (C++ 电路求解)
  → TECPair.update_plasma_heat_flux(q_e, q_c)  → 下发热流
  → TECPair.set_joule_heating(...)             → 下发焦耳热

TFEUnit
  → TECPair.sync_and_step(dt)  → 推进固体温度
```

---

## 7. 依赖关系图

```
HeatConduction2D (Solvers/HeatConduction/HeatConduction.py)
    ├── Fuel
    │     ├── 重写 _compute_internal_resistance() (接触热阻)
    │     ├── set_nuclear_power()  → self.Q_source
    │     └── get_reactivity_feedback()  → 反馈多项式
    │
    ├── Electrode (ABC)
    │     ├── set_joule_heating()  → self.Q_source (外部缓冲)
    │     ├── get_reactivity_feedback()  → 纯虚 (子类实现)
    │     │
    │     ├── Emitter (默认: MoNb)
    │     │     └── get_reactivity_feedback()
    │     │
    │     └── Collector (默认: Molybdenum)
    │           └── get_reactivity_feedback()  [TODO: 修正多项式]
    │
    ├── HeatPipe2D
    │     ├── 吸液芯复合物性 (WickMaterial)
    │     ├── 分段外边界 (outer_eva/aba/con)
    │     ├── 惰性物性更新
    │     ├── 冻结物性校正 (frozen-property correction)
    │     └── 多模态界面热导 (legacy/resistance_split)
    │
    └── FinConduction
          ├── 重写 step()  → spsolve(A, b) 准稳态求解
          ├── 忽略轴向导热 (1D 条降维)
          └── Picard 迭代处理辐射非线性

TECPair (独立类)
    ├── Emitter + Collector + TECCouple2D
    ├── get_gap_surface_temperatures()  → C++ 电路输入
    ├── update_plasma_heat_flux()  → C++ 电路输出  → tec_gap
    ├── set_joule_heating()        → C++ 电路输出  → electrode.Q_source
    └── sync_and_step(dt)  → 推进固体温度

辅助函数:
    Components/tec_electric.py::joule_power_from_electric_field()
        → 电场 + 电阻率 + 体积 → 焦耳热功率
```

---

## 8. 每个组件在 TASTIN 中的角色总览

| 组件 | 输入 | 输出 | 核心物理过程 |
|---|---|---|---|
| `Fuel` | 裂变功率 (标量) | 温度场 + 反馈反应性 | 导热 + 内热源 + 反应性反馈 |
| `Emitter` | 焦耳热 + 电子冷却热流 | 温度场 + 反馈反应性 | 导热 + 焦耳热 |
| `Collector` | 焦耳热 + 电子加热热流 | 温度场 + 反馈反应性 | 导热 + 焦耳热 |
| `HeatPipe2D` | 蒸发段壁面热流 | 温度场 (蒸发→冷凝) | 相变等效导热 (极高的轴向热导率) |
| `FinConduction` | 冷凝段根部温度 | 翅片温度场 | 准稳态导热 + 表面辐射散热 |
| `TECPair` | C++ 电路输出 (J/dU/rho) | 间隙表面温度 + 热流下发 | 电极对封装 + 数据中转 |
