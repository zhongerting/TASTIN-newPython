# basicComponents Detailed Intro

本文档根据 `basicComponents` 目录下源码整理，说明各基础组件的物理职责、类结构、主要接口和耦合关系。

## 目录概览

`basicComponents` 提供 TASTIN 组件层中的基础热物理部件，主要围绕二维导热、热离子电极、燃料芯块、热管和翅片散热展开。

| 文件 | 主要类 | 职责 |
| --- | --- | --- |
| `Electord.py` | `Electrode`, `Emitter`, `Collector` | 电极导热基类及发射极、接收极实现 |
| `TECPair.py` | `TECPair` | 组合发射极和接收极，管理极间隙耦合及热离子源项 |
| `Fuel.py` | `Fuel` | 燃料芯块导热、核功率分配和反应性反馈 |
| `HeatPipe2D.py` | `HeatPipe2D` | 二维柱坐标热管导热求解，区分吸液芯和管壁 |
| `FinConduction.py` | `FinConduction` | 准稳态二维翅片导热和辐射散热求解 |

> 注：文件名 `Electord.py` 可能是 `Electrode.py` 的拼写变体，当前代码中其他模块按 `Components.basicComponents.Electord` 导入，应保持现有名称以避免破坏导入路径。

## 公共设计

这些组件大多继承或组合底层导热求解器：

- `HeatConduction2D`：二维导热瞬态求解基类，维护温度场、热容、内部导热、边界条件和体热源。
- `Mesh2D`：二维网格对象，可用于柱坐标几何。
- `SolidMaterial`：固体材料物性接口，提供导热率、密度、比热等温度相关属性。
- `BoundaryRegion`：边界区域对象，承载热流、热阻、耦合器等边界条件。

组件层的主要作用是把具体物理对象封装成可耦合单元，包括：

- 根据几何生成或接收网格。
- 设置默认材料和初始温度。
- 提供多物理场输入接口，例如核功率、焦耳热、等离子体热流。
- 提供多物理场输出接口，例如表面温度和反应性反馈。
- 将外部物理模型输出转换成导热求解器可接受的节点热源或边界热流。

## `Electord.py`

### `Electrode`

`Electrode` 是发射极和接收极的公共基类，继承自 `HeatConduction2D`。

物理意义：

- 表示带有内部焦耳热源的二维柱坐标电极实体。
- 电极本身只负责热传导和热源映射，不定义具体材料默认值和反应性反馈多项式。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `name` | 组件名称 |
| `mesh` | 二维网格对象，通常应为柱坐标网格 |
| `material` | 电极材料对象 |
| `initial_temp` | 初始均匀温度，单位 K |

关键属性：

- `vols_flat`：展平的一维单元体积数组。
- `total_vol`：电极总体积。
- `nx`, `ny`：二维网格节点数。
- `Q_vol`：体积发热率，单位 W/m3，用于后处理。

主要方法：

#### `set_joule_heating(q_joule_array)`

接收外部计算得到的焦耳热功率数组，并写入导热求解器的 `Q_source`。

要求：

- `q_joule_array` 必须是一维数组。
- 数组长度必须等于电极网格节点数 `self.N`。
- 数组元素单位为 W，表示每个节点分配到的焦耳热功率。

实现要点：

- 使用原地赋值更新 `self.Q_source[:]`，避免改变数组底层内存。
- 设置 `use_external_source_buffer = True`，防止基类在更新热源时清空外部维护的源项。
- 同步计算 `Q_vol = Q_source / volume`。

#### `get_reactivity_feedback()`

电极基类不实现具体反馈关系，直接抛出 `NotImplementedError`。该方法必须由 `Emitter` 或 `Collector` 实现。

### `Emitter`

`Emitter` 表示热离子发射极。

特点：

- 继承自 `Electrode`。
- 默认材料为 `MoNb`。
- 实现发射极温度反应性反馈多项式。

主要方法：

#### `get_reactivity_feedback()`

计算体积加权平均温度：

```text
T_avg = sum(T * volume) / total_volume
```

然后使用多项式计算反应性反馈：

```text
rho_fb = 1e-4 * (3.46455e-6 * T^2 - 0.03232167 * T + 0.74202216)
```

### `Collector`

`Collector` 表示热离子接收极。

特点：

- 继承自 `Electrode`。
- 默认材料为 `Molybdenum`。
- 当前实现使用与 `Emitter` 相同形式的温度反应性反馈多项式。

注意事项：

- 代码中保留了 `TODO 修改接收极反应性反馈函数`，说明接收极反馈关系后续可能需要替换为专属关联式。

## `TECPair.py`

### `TECPair`

`TECPair` 封装一个热离子电极对，包括发射极、极间隙和接收极。

物理职责：

- 构建发射极和接收极的二维柱坐标导热模型。
- 维护极间隙传热、辐射和热离子源项耦合。
- 提供与底层热离子/放电模型交互的标准接口。
- 处理基于面积守恒的热流缩放。
- 将轴向电场或电压降转换为电极内部焦耳热源。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `name` | 电极对名称 |
| `L_node` | 单个轴向节点长度 |
| `n_node` | 轴向节点数 |
| `R_e_in` | 发射极内半径 |
| `delta_e` | 发射极厚度 |
| `delta_gap` | 极间隙厚度 |
| `delta_c` | 接收极厚度 |
| `n_rad_e` | 发射极径向节点数 |
| `n_rad_c` | 接收极径向节点数 |
| `mat_emitter` | 发射极材料，可选 |
| `mat_collector` | 接收极材料，可选 |
| `T_init_e` | 发射极初始温度 |
| `T_init_c` | 接收极初始温度 |
| `k_gap_gas` | 间隙气体导热系数 |
| `emissivity_e` | 发射极表面发射率 |
| `emissivity_c` | 接收极表面发射率 |

几何关系：

```text
R_e_out = R_e_in + delta_e
R_c_in  = R_e_out + delta_gap
R_c_out = R_c_in + delta_c
```

内部对象：

- `self.emitter`：`Emitter` 实例。
- `self.collector`：`Collector` 实例。
- `self.tec_gap`：`TECCouple2D` 实例，连接发射极右边界与接收极左边界。
- `self.inner_boundary`：发射极左边界，通常用于燃料发热耦合。
- `self.outer_boundary`：接收极右边界，通常用于冷却剂或外部散热耦合。

主要方法：

#### `get_gap_surface_temperatures()`

返回靠近极间隙两侧的表面温度分布：

- 发射极右边界温度。
- 接收极左边界温度。

该方法用于向外部热离子或放电模型传递当前电极表面温度。

#### `set_joule_heating(dU_emit, rho_emit, dU_coll, rho_coll, alpha=1.0)`

根据轴向电压降和电阻率分布设置焦耳热。

计算逻辑：

```text
E = dU / L_node
```

随后调用 `set_joule_heating_fields()`，由电场和电阻率计算每个网格的焦耳热功率。

#### `set_joule_heating_fields(E_emit, rho_emit, E_coll, rho_coll, alpha=1.0)`

根据轴向电场直接施加焦耳热。

实现流程：

1. 调用 `joule_power_from_electric_field()` 计算发射极和接收极每个网格单元的焦耳热功率。
2. 对焦耳热源进行亚松弛：

   ```text
   q_new = alpha * q_current + (1 - alpha) * q_old
   ```

3. 调用 `Emitter.set_joule_heating()` 和 `Collector.set_joule_heating()` 写入源项。

#### `update_plasma_heat_flux(q_e_flux, q_c_flux, alpha=1.0)`

施加等离子体放电造成的表面热流。

重要约定：

- `q_e_flux` 和 `q_c_flux` 是面热流密度，单位 W/m2。
- 二者均基于发射极外表面积定义。
- 为保持功率守恒，接收极侧热流也乘以发射极面积转换为总功率。

实现流程：

1. 读取发射极右边界面积和接收极左边界面积。
2. 对传入热流密度进行亚松弛。
3. 将面热流密度转换为功率：

   ```text
   Q_e = q_e * A_emit_surface
   Q_c = q_c * A_emit_surface
   ```

4. 通过 `tec_gap.set_tec_sources()` 设置极间隙热离子源项。
5. 保存 `plasma_area_diagnostics`，用于检查面积基准和功率差异。

#### `sync_and_step(dt)`

执行一个时间步：

1. `tec_gap.sync()` 同步极间隙耦合。
2. `emitter.step(dt)` 推进发射极导热。
3. `collector.step(dt)` 推进接收极导热。

## `Fuel.py`

### `Fuel`

`Fuel` 表示燃料芯块，继承自 `HeatConduction2D`。

物理意义：

- 带内热源的二维柱坐标导热体。
- 与点堆模型进行双向耦合。
- 输入核功率，输出燃料温度反应性反馈。
- 可在指定轴向交界面加入接触热阻。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `name` | 燃料组件名称 |
| `mesh` | 二维网格对象 |
| `material` | 燃料材料，缺省为 `UO2` |
| `initial_temp` | 初始温度，单位 K |
| `power_fraction` | 该燃料组件分配到的总功率比例 |
| `axial_power_profile` | 轴向功率分布数组 |
| `contact_resistance_interfaces` | 存在接触热阻的轴向交界面索引列表 |
| `axial_contact_resistance` | 交界面的总接触热阻，单位 K/W |

关键属性：

- `vols_flat`：展平体积数组。
- `vols_2d`：二维体积矩阵。
- `total_vol`：总体积。
- `vols_layer`：每个轴向层的体积总和。
- `Q_vol`：体积发热率，单位 W/m3。
- `power_allocation_weights`：每个节点的功率分配权重。

主要方法：

#### `_compute_internal_resistance()`

重载父类内部热阻计算，在指定轴向交界面串联接触热阻。

处理逻辑：

1. 先调用父类计算原始内部热阻。
2. 如果未启用接触热阻则直接返回。
3. 对每个 `contact_resistance_interfaces` 中的轴向交界面，计算环向节点对应的接触热阻。
4. 将原始热阻和接触热阻串联，再更新轴向导热导纳 `G_y_inner`。

#### `set_axial_power_profile(profile_array=None)`

设置或更新轴向功率分布，并预计算节点功率分配权重。

规则：

- 如果未传入功率分布，默认各轴向层均匀分配。
- 如果传入数组，长度必须等于轴向节点数 `ny`。
- 分布会自动归一化，防止总功率被错误缩放。
- 每个轴向层内部再按该层各径向单元体积比例分配。

最终满足：

```text
sum(power_allocation_weights) = 1
```

#### `set_nuclear_power(p_fiss, p_decay, p_total)`

接收点堆模型输出的核功率，并写入导热热源。

当前实现使用：

```text
component_power = p_total * power_fraction
Q_source = component_power * power_allocation_weights
Q_vol = Q_source / volume
```

其中 `p_fiss` 和 `p_decay` 当前未单独参与分配，接口保留了这两个参数以适配多物理场输入。

#### `get_reactivity_feedback()`

计算燃料体积加权平均温度，并代入燃料温度反馈多项式：

```text
rho_fb = 0.001360811
       - 6.47927757e-6 * T
       + 2.321231e-9 * T^2
       - 3.52e-13 * T^3
```

返回值为该燃料组件引入的反应性反馈。

## `HeatPipe2D.py`

### `HeatPipe2D`

`HeatPipe2D` 是二维柱坐标热管导热求解器，继承自 `HeatConduction2D`。

物理特点：

- 径向上区分吸液芯区域和管壁区域。
- 轴向外壁被划分为蒸发段、绝热段和冷凝段。
- 支持吸液芯等效材料 `WickMaterial`。
- 支持不同面导热导纳计算模式。
- 支持冻结物性外迭代修正，以提升非线性物性问题的稳定性。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `mesh` | 二维网格对象 |
| `solid1` | 管壁材料 |
| `solid2` | 吸液芯内工质材料 |
| `solid3` | 吸液芯固体骨架材料 |
| `n_wick` | 吸液芯径向节点数 |
| `porosity` | 吸液芯孔隙率 |
| `n_eva` | 蒸发段轴向节点数 |
| `n_aba` | 绝热段轴向节点数 |
| `n_con` | 冷凝段轴向节点数 |
| `name` | 组件名称 |
| `emissivity` | 表面发射率 |
| `up_view_factor` | 上表面角系数 |
| `down_view_factor` | 下表面角系数 |
| `initial_temp` | 初始温度 |

输入约束：

- `n_wick` 必须大于 0 且小于总径向节点数 `mesh.n_x`。
- `n_eva + n_aba + n_con` 必须等于轴向节点数 `mesh.n_y`。

边界划分：

`HeatPipe2D` 删除默认外侧 `right` 边界，并创建三个外壁虚拟边界：

- `outer_eva`：蒸发段外壁。
- `outer_aba`：绝热段外壁。
- `outer_con`：冷凝段外壁。

其他默认边界如 `left`、`bottom`、`top` 仍按父类逻辑存在。

关键属性和模式：

- `use_anisotropic_wick_conductivity`：是否使用吸液芯各向异性导热率。
- `face_conductance_mode`：面导热导纳计算模式。
- `enable_frozen_property_correction`：是否启用冻结物性外迭代。
- `minimum_physical_temperature`、`maximum_physical_temperature`：物理温度上下限。

### 面导热导纳模式

通过 `set_face_conductance_mode(mode)` 设置：

| 模式 | 说明 |
| --- | --- |
| `legacy_harmonic` | 径向和轴向均使用传统调和平均 |
| `resistance_split_axial` | 径向使用传统调和平均，轴向使用半单元热阻串联 |
| `resistance_split_xy` | 径向和轴向均使用热阻拆分 |
| `resistance_split_full` | 内部面和边界半单元均使用更完整的热阻拆分 |

### 主要方法

#### `_setup_virtual_boundaries()`

根据外壁面积和轴向分段建立 `outer_eva`、`outer_aba`、`outer_con` 三个边界区域。

#### `set_wick_conductivity_mode(anisotropic)`

设置吸液芯是否使用各向异性导热率：

- `False`：使用单一等效导热率。
- `True`：分别调用 `conductivity_axial()` 和 `conductivity_radial()`。

设置后会使物性缓存失效。

#### `_update_properties()`

根据当前温度更新材料物性。

实现特点：

- 将一维温度和物性数组重塑为二维视图，提高访问效率。
- 吸液芯区域使用 `WickMaterial`。
- 管壁区域使用 `wall_mat`。
- 通过温度变化阈值 `_property_update_tol` 减少不必要的物性更新。
- 对高非线性温区强制更新吸液芯物性。
- 更新热容 `rho * cp * volume`。

#### `_compute_internal_resistance()`

根据 `face_conductance_mode` 选择内部导热导纳计算策略。

#### `_update_boundaries_state(current_time=None)`

更新各边界的内部状态，包括边界节点温度和内部半单元热阻。

对于外壁分段边界：

- `outer_eva` 使用外壁温度的蒸发段切片。
- `outer_aba` 使用外壁温度的绝热段切片。
- `outer_con` 使用外壁温度的冷凝段切片。

#### `_compute_fluxes(t)`

计算内部导热和边界条件贡献的净热流，返回展平的一维数组。

包含：

- 径向内部导热。
- 轴向内部导热。
- `left`、`bottom`、`top` 边界热流。
- `outer_eva`、`outer_aba`、`outer_con` 外壁分段边界热流。

#### `get_derivatives(t, T_current)`

计算温度导数。

当未冻结物性时，直接调用父类逻辑；当冻结物性时，使用当前缓存物性和导热导纳计算导数。

#### `step(dt, method='BDF', **kwargs)`

推进一个时间步。

两种路径：

- 未启用冻结物性修正：调用父类 `step()`，并检查求解成功和温度物理性。
- 启用冻结物性修正：执行外层物性修正迭代，在冻结物性条件下调用 `_solve_ivp_step()`，直到温度变化满足 `outer_property_tol`。

失败处理：

- 若求解器失败或温度非物理，将回滚到上一步状态。
- 失败信息会包含试算温度和边界条件诊断。

#### `get_boundary_node_capacitance(location)`

返回指定边界节点对应热容。

对 `outer_eva`、`outer_aba`、`outer_con` 做了特殊切片处理；其他边界调用父类方法。

## `FinConduction.py`

### `FinConduction`

`FinConduction` 是准稳态二维翅片导热求解器，继承自 `HeatConduction2D`，但重写了 `step()`。

物理假设：

- 翅片热容很小，在每个全局时间步内瞬时达到稳态。
- 忽略沿热管轴向，即 Y 方向的导热。
- 将正反两面辐射散热等效为节点热漏。
- 边缘散热仍通过标准边界条件处理。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `mesh` | 二维网格 |
| `material` | 翅片材料 |
| `fin_thickness` | 翅片厚度 |
| `emissivity` | 发射率 |
| `up_view_factor` | 上表面角系数 |
| `down_view_factor` | 下表面角系数 |
| `T_env` | 环境温度，可为标量或数组 |
| `initial_temp` | 初始温度 |

主要方法：

#### `_get_boundary_linear_terms(boundary_key)`

将复杂边界条件转换为稀疏稳态方程所需的线性项。

返回：

- `G_bound`：等效边界导纳。
- `S_bound`：等效源项。

处理的边界条件类型：

- `resistance`：热阻边界，转换为戴维南等效形式。
- `flux`：纯热流边界，直接加入源项。

等效形式：

```text
Flux = G_bound * T_eff - G_bound * T_node + Q_flux_only
```

#### `step(dt, max_iter=50, tol=1e-4, relaxation=0.8, **kwargs)`

重写基类瞬态 ODE 求解，改为非线性稳态 Picard 迭代。

求解流程：

1. 保存旧温度 `T_old`。
2. 更新材料物性、内部热阻和边界状态。
3. 组装稀疏线性系统 `A * T_new = b`。
4. 只加入 X 方向内部导热，忽略 Y 方向导热。
5. 对正反面辐射进行局部线性化：

   ```text
   h_rad = emissivity * sigma * (T_old^2 + T_env^2) * (T_old + T_env)
   ```

6. 根据上下角系数和翅片厚度计算有效辐射面积：

   ```text
   A_rad_eff = (up_view_factor + down_view_factor) * volume / fin_thickness
   ```

7. 加入四条边界条件的线性贡献。
8. 加入外部体热源 `Q_source`。
9. 使用 `scipy.sparse.linalg.spsolve()` 求解。
10. 使用松弛更新温度并检查收敛。

如果达到最大迭代次数仍未满足收敛准则，会打印警告，但仍更新时间并返回 `True`。

## 组件耦合关系

### 燃料到电极

`Fuel` 可作为内部热源区域，电极对通过 `TECPair.inner_boundary` 暴露发射极内边界，用于燃料和发射极之间的热耦合。

### 电极对内部

`TECPair` 管理：

- 发射极导热。
- 接收极导热。
- 极间隙导热、辐射和热离子源项。

典型数据流：

```text
外部热离子模型读取表面温度
    <- TECPair.get_gap_surface_temperatures()

外部热离子模型返回热流和电学结果
    -> TECPair.update_plasma_heat_flux()
    -> TECPair.set_joule_heating() 或 set_joule_heating_fields()

导热推进
    -> TECPair.sync_and_step(dt)
```

### 热管和翅片

`HeatPipe2D` 将外壁分为蒸发段、绝热段和冷凝段，冷凝段可与 `FinConduction` 或其他散热组件通过边界耦合器连接。

`FinConduction` 将翅片作为准稳态散热部件，适合在全局瞬态系统中快速响应热管冷凝段传入的热量。

## 典型使用顺序

一个多物理场时间步通常遵循以下顺序：

1. 从导热组件读取当前边界或表面温度。
2. 调用外部物理模型计算核功率、热离子热流、电场、电阻率或外部换热条件。
3. 将外部模型输出写回组件：
   - `Fuel.set_nuclear_power()`
   - `TECPair.update_plasma_heat_flux()`
   - `TECPair.set_joule_heating()`
   - 边界条件或耦合器接口
4. 调用各组件 `step(dt)` 或 `sync_and_step(dt)` 推进热状态。
5. 收集温度、热流、体热源和反应性反馈等结果。

## 注意事项

- 所有温度默认单位为 K。
- `Q_source` 通常表示节点总功率，单位 W；`Q_vol` 表示体积发热率，单位 W/m3。
- 电极等离子体热流输入按发射极面积基准定义，代码中通过面积乘法保证总功率守恒。
- `Fuel.set_nuclear_power()` 当前使用 `p_total` 分配功率，`p_fiss` 和 `p_decay` 作为接口兼容参数保留。
- `HeatPipe2D` 的冻结物性修正可提高强非线性物性问题的鲁棒性，但会增加每步计算量。
- `FinConduction.step()` 是稳态迭代求解，不使用父类瞬态 ODE 积分。
- 修改 `Electord.py` 文件名或类名时，需要同步修改所有导入路径。
