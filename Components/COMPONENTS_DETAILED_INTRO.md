# Components Detailed Intro

本文档是 Codex 助手和开发者进入 `Components` 的唯一入口文档。常规定位、小范围修改和调用关系确认应优先只读本文；涉及物理公式、底层实现或接口变更时，再按本文链接读取专题文档和对应源码。

## 阅读导航

### 目录树

```text
Components/
├── BaseComponent.py
├── Pipe.py
├── AnnularPipe.py
├── HPwithFin.py
├── RingHP.py
├── TFEUnit.py
├── TECCircuitManager.py
├── tec_electric.py
├── COMPONENTS_DETAILED_INTRO.md
├── BASICCOMPONENTS_DETAILED_INTRO.md
├── BASICCOMPONENTS_ANALYSIS.md
├── EXTERNALHEATSOURCES_DETAILED_INTRO.md
├── basicComponents/
│   ├── Electord.py
│   ├── TECPair.py
│   ├── Fuel.py
│   ├── HeatPipe2D.py
│   └── FinConduction.py
└── ExternalHeatSources/
    ├── __init__.py
    ├── embedded_flux_tables.py
    ├── README.md
    └── PARAMETERS.md
```

### 文档职责

| 文档 | 用途 | 何时阅读 |
| --- | --- | --- |
| 本文 | 目录结构、组件关系、接口、配置、数据流、单位和风险总览 | 每次进入 `Components` 时优先阅读 |
| [`BASICCOMPONENTS_DETAILED_INTRO.md`](./BASICCOMPONENTS_DETAILED_INTRO.md) | `Fuel`、电极、`TECPair`、`HeatPipe2D`、`FinConduction` 的接口说明 | 修改基础热工元件或 TEC 细节时 |
| [`EXTERNALHEATSOURCES_DETAILED_INTRO.md`](./EXTERNALHEATSOURCES_DETAILED_INTRO.md) | 轨道外热源、查表热流和边界封装 | 修改外热模型或辐射器外热加载时 |
| [`BASICCOMPONENTS_ANALYSIS.md`](./BASICCOMPONENTS_ANALYSIS.md) | 基础元件物理实现和调用链分析 | 修改公式、求解方法或排查物理结果时 |

外热源的补充资料位于 [`ExternalHeatSources/README.md`](./ExternalHeatSources/README.md) 和 [`ExternalHeatSources/PARAMETERS.md`](./ExternalHeatSources/PARAMETERS.md)。

### 何时必须重新阅读源码

常规定位和小范围修改无需重复遍历整个目录。出现以下情况时，必须打开专题文档和目标源码局部核验：

- 修改物理公式、经验关联式、热阻、边界符号或单位换算。
- 修改构造函数、返回对象、状态保存键或 `SystemManager` 生命周期接口。
- 修改 `ReactorCore` 倍率、功率分配、全局慢化剂环映射或点堆推进顺序。
- 修改非均匀网格、电场梯度、焦耳热映射或 TEC 面热流。
- 修改 `RingHP.external_heat_config`、查表热流、翅片受照面积或净散热统计。
- 新增 Python 文件、组件类型或外部依赖。

## 模块定位

`Components` 根目录提供更高层的宏观组件和系统组件。它们通常不直接实现单一导热方程，而是把多个基础固体、流体通道、气隙、耦合器和外部物理模型组装成可交给 `SystemManager` 管理的对象。

统一接口来自 `BaseComponent`：

- `get_solids()`：返回组件内部所有固体导热求解器。
- `get_couplers()`：返回组件内部所有耦合器。
- `pre_step(dt, current_time)`：全局时间步求解前执行。
- `post_step(dt, current_time)`：全局时间步求解后执行。

对于支持内迭代回滚和断点续算的组件，还应提供：

- `save_step_state()` / `load_step_state(state)`：保存和恢复当前时间步内迭代状态。
- `get_state_dict(prefix)` / `load_state_dict(data, prefix)`：生成和恢复可序列化的全局断点状态。

## 架构总览

```text
ReactorCore
  -> TFEUnit
    -> Fuel / Emitter / Collector / Clads / Moderator
    -> GapCouple2D / TECCouple2D / FluidSolidCouple
  -> ThermoCalcModel
  -> PointReactor

RingHP
  -> HPwithFin
    -> HeatPipe2D
    -> ExternalHeatSources
  -> FluidSolidCouple
```

两条主路径：

1. 堆芯路径：`ReactorCore` 组织代表性 `TFEUnit`，根据倍率还原全堆统计，并在需要时接入 `ThermoCalcModel` 和 `PointReactor`。
2. 散热器路径：`RingHP` 用一个流体节点对应一根代表性 `HPwithFin`，通过 `hp_multipliers` 放大到真实热管数量。

`BaseComponent` 与 `SystemManager` 的契约是“宏观组件负责暴露底层实体，系统管理器负责统一推进”。新增宏观组件时，必须完整返回内部固体和耦合器，并保证生命周期钩子可重复调用。

## 文件概览

| 文件 | 主要内容 |
| --- | --- |
| `BaseComponent.py` | 宏观组件基类 |
| `Pipe.py` | 单壁管道流固耦合组件 |
| `AnnularPipe.py` | 环形管道双壁流固耦合组件 |
| `HPwithFin.py` | 带降维翅片的热管散热器 |
| `RingHP.py` | 集流环与代表性热管阵列组件 |
| `TFEUnit.py` | 热离子燃料元件装配体 |
| `TECCircuitManager.py` | 多个 `TECPair` 的热离子电路管理器 |
| `ReactorCore.py` | 堆芯级容器、全局结构、点堆和 TEC 耦合 |
| `tec_electric.py` | 电场和焦耳热映射工具函数 |

## `BaseComponent.py`

### `BaseComponent`

所有宏观组件的抽象基类。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `name` | 组件名称 |

默认行为：

- `get_solids()` 返回空列表。
- `get_couplers()` 返回空列表。
- `pre_step()` 和 `post_step()` 为空钩子。

子类应按自身组合关系重写这些接口，使 `SystemManager` 能自动收集底层求解对象。

## `Pipe.py`

### `Pipe`

单壁管道组件，由一个固体管壁和一个流体通道组成。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `name` | 组件名称 |
| `solid_wall` | 固体管壁对象，通常为 `HeatConduction2D` |
| `fluid_channel` | 流体通道对象 |
| `heated_perimeter` | 换热湿周，单位 m |
| `correlation_func` | 对流换热关联式 |
| `coupled_boundary_name` | 固体与流体耦合的边界名，默认 `left` |

构造流程：

1. 保存固体和流体对象。
2. 调用 `solid.initialize_state()`。
3. 获取耦合边界节点热容。
4. 构造 `FluidSolidCouple`。

对外暴露：

- `get_solids()` 返回 `[solid]`。
- `get_couplers()` 返回 `[coupler]`。

## `AnnularPipe.py`

### `AnnularPipe`

环形管道组件，由内管壁、外管壁和两者之间的环形流体通道组成。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `solid_inner` | 内侧固体壁面 |
| `solid_outer` | 外侧固体壁面 |
| `fluid_channel` | 环形流体通道 |
| `heated_perimeter_inner` | 内壁湿周 |
| `heated_perimeter_outer` | 外壁湿周 |
| `correlation_func` | 对流换热关联式 |
| `boundary_inner_solid` | 内壁耦合边界，默认 `right` |
| `boundary_outer_solid` | 外壁耦合边界，默认 `left` |

构造结果：

- `coupler_inner`：流体与内壁的换热耦合器。
- `coupler_outer`：流体与外壁的换热耦合器。

对外暴露：

- `get_solids()` 返回内外两个固体。
- `get_couplers()` 返回两个流固耦合器。

## `tec_electric.py`

电场和焦耳热计算工具模块，供 `TECPair`、`TFEUnit`、`TECCircuitManager` 和 `ReactorCore` 使用。

### `node_centers_from_lengths(node_lengths)`

由轴向单元长度计算单元中心坐标。

要求：

- `node_lengths` 必须是一维数组。
- 所有长度必须大于 0。

### `node_centers_from_faces(y_faces)`

由轴向面坐标计算单元中心坐标。

要求：

- `y_faces` 必须是一维数组。
- 至少包含两个面坐标。
- 坐标必须严格递增。

### `electric_field_from_node_potential(potential, y_faces=None, node_lengths=None)`

根据轴向节点电势计算中心处电场梯度。

特点：

- 支持非均匀轴向网格。
- 优先使用 `y_faces`，也可使用 `node_lengths`。
- 使用 `np.gradient()` 计算 `dU/dy`。
- 单节点时返回 0 数组。

### `joule_power_from_electric_field(electric_field, resistivity, volumes_flat, shape_nodes)`

把一维轴向电场映射为二维网格节点焦耳热功率。

计算关系：

```text
q_vol_1d = E^2 / rho
q_watts_flat = broadcast(q_vol_1d) * volume
```

返回：

- `q_watts_flat`：展平后的每节点焦耳热功率，单位 W。
- `q_vol_1d`：轴向体积发热率，单位 W/m3。

### `joule_power_from_node_potential(...)`

组合接口：先由电势求电场，再由电场求焦耳热。

返回：

- 节点焦耳热功率。
- 轴向体积发热率。
- 轴向电场。

## `TECCircuitManager.py`

### `TECCircuitManager`

热离子电热耦合电路管理器，用于管理多个 `TECPair`。

工作范式：

- 在 `pre_step()` 中读取各电极对表面温度。
- 调用 `ThermoCalcModel` 计算全局 TEC 电路。
- 将电子冷却/加热热流和焦耳热下发到各个 `TECPair`。

初始化参数：

| 参数 | 说明 |
| --- | --- |
| `name` | 管理器名称 |
| `tfe_list` | `TECPair` 列表 |
| `Tcs_init` | 铯池初始温度 |
| `V_init` | 初始电压估值 |

使用约束：

- `tfe_list` 不能为空。
- 所有 `TECPair.n_node` 必须一致。

主要接口：

- `get_solids()`：返回所有发射极和接收极固体。
- `get_couplers()`：返回所有极间隙 `TECCouple2D`。
- `pre_step()`：调用 `sync_thermo_electric()`。
- `sync_thermo_electric(alpha, macro_params)`：执行 gather-calculate-scatter。

热流符号约定：

- 发射极电子冷却热流为负。
- 接收极电子加热热流为正。
- 热流密度基准为发射极外表面积。

## `TFEUnit.py`

### 数据类

#### `TFEGeometry`

热离子燃料元件几何参数，从内到外描述燃料芯块、裂变气隙、发射极、接收极、套管、冷却剂流道、外套管和慢化剂。

主要字段：

- `r_pellet_inner`, `r_pellet_outer`
- `r_fission_gas_outer`
- `r_emitter_outer`
- `r_collector_inner`, `r_collector_outer`
- `r_inner_clad_inner`, `r_inner_clad_outer`
- `r_coolant_inner`, `r_coolant_outer`
- `r_outer_clad_outer`
- `r_moderator_inner`, `r_moderator_outer`
- `height`

#### `TFEMeshParams`

TFE 径向和轴向网格划分参数。

字段包括：

- `n_axial`
- `n_r_pellet`
- `n_r_emitter`
- `n_r_collector`
- `n_r_inner_clad`
- `n_r_outer_clad`
- `n_r_moderator`

#### `GapConfig`

单个气隙配置。

模式：

- `simplified`：使用等效换热系数和辐射参数，通过气隙耦合器表示。
- `meshed`：把气隙建成真实网格固体，需要提供 `material` 和 `n_radial_nodes`。

字段：

- `mode`
- `h_eq`
- `material`
- `emissivity_inner`
- `emissivity_outer`
- `n_radial_nodes`

### `TFEUnit`

热离子燃料元件装配体，集成：

- 燃料芯块 `Fuel`
- 发射极 `Emitter`
- 接收极 `Collector`
- 内套管、外套管、慢化剂等固体
- 裂变气隙、极间隙、氦气隙、CO2 气隙
- 冷却剂流固耦合
- 核功率、焦耳热和等离子体热流输入接口

主要初始化参数：

| 参数 | 说明 |
| --- | --- |
| `geometry` | `TFEGeometry` |
| `mesh_params` | `TFEMeshParams` |
| `materials` | 材料字典 |
| `coolant_channel` | 冷却剂流体通道 |
| `fission_gas_config` | 裂变气隙配置 |
| `tec_gap_config` | 极间隙配置 |
| `he_gap_config` | 氦气隙配置 |
| `co2_gap_config` | CO2 气隙配置 |
| `power_fraction` | 功率份额 |
| `axial_power_profile` | 轴向功率分布 |
| `axial_length_allocation` | 轴向非均匀长度分配 |
| `axial_node_allocation` | 轴向非均匀节点分配 |
| `axial_contact_resistance` | 轴向接触热阻 |

主要接口：

- `update_neutronic_power()`：把外部总功率映射到燃料芯块。
- `update_electric_fields()`：由电压降和电阻率设置电极焦耳热。
- `update_electric_field_sources()`：由电场和电阻率设置电极焦耳热。
- `update_plasma_flux()`：设置极间隙等离子体表面热流。
- `get_solids()`：返回所有内部固体。
- `get_couplers()`：返回所有内部耦合器。
- `save_step_state()` / `load_step_state()`：保存和恢复时间步状态。
- `get_state_dict()` / `load_state_dict()`：断点续算序列化接口。

## `HPwithFin.py`

### `HPwithFin`

带降维翅片的热管散热器。

建模思想：

- 热管本体由 `HeatPipe2D` 显式求解。
- 翅片不作为二维或三维网格，而是降维为每个冷凝段轴向切片上的一维准稳态导热问题。
- 翅片解算结果转换为挂在热管冷凝段外壁上的等效热阻支路。

主要初始化参数：

| 参数 | 说明 |
| --- | --- |
| `r_out_wall`, `r_in_wall`, `r_vapor` | 热管径向几何 |
| `L_eva`, `L_aba`, `L_con` | 蒸发段、绝热段、冷凝段长度 |
| `n_eva`, `n_aba`, `n_con` | 轴向分段节点数 |
| `n_wick`, `n_wall` | 吸液芯和管壁径向节点数 |
| `wall_mat`, `fluid_mat`, `wick_struct_mat` | 管壁、工质、吸液芯骨架材料 |
| `porosity` | 吸液芯孔隙率 |
| `fin_thickness`, `fin_height`, `n_fin_height` | 翅片厚度、高度和降维节点数 |
| `fin_wrap_ratio` | 冷凝段周向覆盖比例 |
| `emissivity` | 发射率 |
| `up_view_factor`, `down_view_factor` | 内侧上下角系数 |
| `T_env` | 环境温度 |

外边界：

- `outer_aba`：绝热段裸壁辐射。
- `outer_con`：冷凝段裸壁辐射。
- 冷凝段还挂有动态更新的翅片等效热阻边界。

主要接口：

- `set_fin_external_heat_source()`：给翅片直接受照模型挂外热流。
- `configure_external_heat_accounting()`：配置外热流吸收后处理核算。
- `get_fin_illuminated_area_array()`：返回翅片受照面积数组。
- `get_external_heat_absorption_distribution()`：返回当前外热流吸收分布。
- `pre_step()`：每个时间步前求解翅片准稳态问题并更新等效热阻。
- `get_heat_rejection_distribution()`：返回绝热段与冷凝段散热分布。
- `get_heat_exchange_breakdown()`：返回裸壁、翅片辐射、外热吸收和净散热拆分。
- `get_temperature_distribution()`：返回热管二维温度场。

## `RingHP.py`

### `SingleVolumeProxy`

单节点流体代理类。

用途：

- 把一个流体控制体伪装成单节点流体通道。
- 让 `FluidSolidCouple` 能复用在“一个控制体对应一根代表热管”的场景。
- 将单根代表热管耦合源项按 `N_hp` 放大回真实热管数量。

### `RingHP`

集流环加代表性热管阵列组件。

内部包含：

- 集流环流体与集流环壁面的流固换热。
- 每个流体控制体位置的代表性 `HPwithFin`。
- 每根代表热管蒸发段与对应流体控制体的流固耦合。
- 可选轨道外热流配置。

主要初始化参数：

| 参数 | 说明 |
| --- | --- |
| `fluid_channel` | 集流环流体通道 |
| `solid_header` | 集流环壁面固体 |
| `hp_multipliers` | 每个流体节点代表的热管数量 |
| `header_flow_area`, `header_dh`, `header_heated_perimeter` | 集流环水力几何 |
| `hp_*` | 热管几何、网格和材料参数 |
| `fin_*` | 翅片几何参数 |
| `header_correlation_func` | 集流环换热关联式 |
| `hp_crossflow_base_func` | 热管蒸发段横掠换热关联式 |
| `external_heat_config` | 外热流配置 |

主要接口：

- `get_solids()`：返回集流环壁面和所有存在的热管固体。
- `get_couplers()`：返回集流环耦合器和热管蒸发段耦合器。
- `pre_step()` / `post_step()`：转发给内部热管。
- `get_total_heat_rejection()`：返回未乘真实数量的代表热管总散热。
- `get_total_heat_rejection_scaled()`：按 `hp_multipliers` 放大后的总散热。
- `get_total_external_heat_absorption_scaled()`：按真实数量统计外热流吸收。
- `get_total_net_heat_rejection_scaled()`：真实尺度净排热。
- `get_hp_status_summary()`：返回热管存在性、散热、吸热、压降等状态摘要。

## `ReactorCore.py`

### 数据类

#### `GlobalAnnulusStructureConfig`

堆芯外部全局环形固体层配置，用于筒体层或反射层。

支持直接传入 `Mesh2D`，也支持用半径、厚度和径向节点数自动生成柱坐标环形网格。

#### `GlobalGapStructureConfig`

堆芯外部环形间隙配置。

模式：

- `simplified`：使用 `GapCouple2D`。
- `meshed`：把间隙本身建成 `HeatConduction2D`。

#### `ReactivityFeedbackTemperatureSummary`

保存燃料、发射极、接收极、慢化剂、反射层和筒体的平均温度。

#### `ReactivityFeedbackResult`

保存各项反应性反馈和总反馈。

### `ReactorCore`

堆芯级宏观容器。

主要职责：

- 组织多个代表性 `TFEUnit`。
- 按真实数量倍率统计功率、温度和反应性反馈。
- 构建全局慢化剂环、筒体、反射层及层间间隙。
- 管理 TEC 电路耦合和点堆耦合。
- 提供断点续算状态保存。

主要初始化参数：

| 参数 | 说明 |
| --- | --- |
| `tfe_dict` | 代表性 TFE 字典 |
| `tfe_multipliers` | 每个代表 TFE 对应的真实数量 |
| `tec_multipliers` | TEC 电路数量倍率 |
| `tfe_power_factors` | 功率分配因子 |
| `mod_meshes`, `mod_material`, `ring_mapping` | 全局慢化剂环建模 |
| `barrel_config` | 筒体配置 |
| `reflector_config` | 反射层配置 |
| `moderator_barrel_gap_config` | 慢化剂-筒体间隙 |
| `barrel_reflector_gap_config` | 筒体-反射层间隙 |
| `T_space` | 外部空间温度 |
| `alpha_tec` | TEC 热源亚松弛系数 |
| `enable_tec_coupled` | 是否启用 TEC 耦合 |

主要接口：

- `collect_reactivity_feedback_temperatures()`：收集反馈用结构平均温度。
- `compute_reactivity_feedback()`：计算燃料、电极、慢化剂、反射层等反馈。
- `get_reactivity_feedback()`：返回总反馈。
- `calibrate_reactivity_feedback_reference()`：建立反馈参考态。
- `get_effective_reactivity_feedback()`：返回扣除参考态后的有效反馈。
- `setup_tec_circuit()`：配置 ThermoCalc 电路模式。
- `attach_point_reactor()` / `initialize_point_reactor()`：挂接并初始化点堆模型。
- `update_neutronic_power()`：将堆功率分配给代表性 TFE。
- `advance_neutronics()` / `commit_neutronics()`：推进并提交点堆状态。
- `get_solids()` / `get_couplers()`：向 `SystemManager` 暴露所有底层对象。
- `pre_step()`：执行 TEC 计算、慢化剂源项转移和 TFE 预处理。
- `post_step()`：同步 TEC 表面温度和全局慢化剂温度。
- `get_state_dict()` / `load_state_dict()`：堆芯级断点续算。

## TFE 配置手册

### 数据域

`TFEUnit.py` 在构造体旁定义了四个运行期数据域：

| 数据类 | 作用 | 关键内容 |
| --- | --- | --- |
| `NeutronicData` | 中子学热源缓存 | 当前总功率和上一次功率，用于亚松弛 |
| `ElectricFieldData` | 电场和焦耳热缓存 | 发射极/接收极节点电势、电阻率、电流密度和体积焦耳热 |
| `PlasmaCouplingData` | 极间等离子体接口 | 功函数、阻塞压降、发射极温度、发射极失热和接收极得热 |
| `BoundaryConditionData` | 外部边界输入 | 等效慢化剂外部温度分布 |

### 材料键

`TFEUnit.materials` 使用固定字符串键。基础模型必须提供：

| 键 | 用途 |
| --- | --- |
| `UO2` | 燃料芯块 |
| `MoNb` | 发射极 |
| `Molybdenum` | 接收极 |
| `StainlessSteel` | 内外套管 |
| `ZrH` | TFE 内部等效慢化剂 |

当 `axial_length_allocation` 和 `axial_node_allocation` 引入上下轴向反射段时，还必须提供 `BerylliumOxide` 或别名 `BeO`。该材料用于构造 `CompositePelletMaterial`。

### 四类气隙

| 参数 | 位置 | `simplified` 模式 | `meshed` 模式 |
| --- | --- | --- | --- |
| `fission_gas_config` | 芯块到发射极 | `GapCouple2D` | 网格化气隙固体、固固导热和并联辐射 |
| `tec_gap_config` | 发射极到接收极 | `TECCouple2D` 及气隙处理 | 网格化气隙固体及对应耦合 |
| `he_gap_config` | 接收极到内套管 | `GapCouple2D` | 网格化气隙固体、固固导热和并联辐射 |
| `co2_gap_config` | 外套管到慢化剂 | `GapCouple2D` | 网格化气隙固体、固固导热和并联辐射 |

`GapConfig(mode='simplified', h_eq=...)` 使用等效换热系数。`GapConfig(mode='meshed', material=..., n_radial_nodes=...)` 将气隙本身建成固体网格；此模式必须同时提供材料和正整数径向节点数。

### 严格绝热单 TFE 模式

`TFEUnit(..., strict_adiabatic_single_tfe=True)` 仅用于单根 TFE 能量守恒诊断。默认值为 `False`，现有堆芯行为不变。

启用后：

- 不创建内部等效 moderator 固体。
- 不创建 CO2 间隙网格，也不注册外套管到 moderator 的耦合器。
- 清空外套管右边界条件，使该边界逐节点严格为零热流。

标准堆芯装配不得开启该选项。专题入口见 [`../testModule/SINGLE_TFE_ENERGY_CONSERVATION_GUIDE.md`](../testModule/SINGLE_TFE_ENERGY_CONSERVATION_GUIDE.md)。

### 最小 `TFEUnit` 模板

下面模板可直接复制。`coolant_channel` 需要先按水力网络接口构造，通常为 `IncompressibleFluidChannel`。

```python
from Components.TFEUnit import TFEGeometry, TFEMeshParams, GapConfig, TFEUnit
from Materials.Solids.UO2 import UO2
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.ZrH import ZirconiumHydride

geometry = TFEGeometry(
    r_pellet_inner=4.00e-3,
    r_pellet_outer=8.50e-3,
    r_fission_gas_outer=8.65e-3,
    r_emitter_outer=9.80e-3,
    r_collector_inner=10.30e-3,
    r_collector_outer=11.85e-3,
    r_inner_clad_inner=11.90e-3,
    r_inner_clad_outer=12.25e-3,
    r_coolant_inner=12.25e-3,
    r_coolant_outer=12.95e-3,
    r_outer_clad_outer=13.30e-3,
    r_moderator_inner=13.52e-3,
    r_moderator_outer=14.52e-3,
    height=0.375,
)
mesh = TFEMeshParams(
    n_axial=30,
    n_r_pellet=5,
    n_r_emitter=1,
    n_r_collector=1,
    n_r_inner_clad=2,
    n_r_outer_clad=2,
    n_r_moderator=3,
)
materials = {
    "UO2": UO2(),
    "MoNb": MoNb(),
    "Molybdenum": Molybdenum(),
    "StainlessSteel": AusteniticStainlessSteel(),
    "ZrH": ZirconiumHydride(),
}
tfe = TFEUnit(
    name="TFE_Center",
    geometry=geometry,
    mesh_params=mesh,
    materials=materials,
    coolant_channel=coolant_channel,
    fission_gas_config=GapConfig(mode="simplified", h_eq=5678.0),
    tec_gap_config=GapConfig(mode="simplified", h_eq=168.0),
    he_gap_config=GapConfig(mode="simplified", h_eq=168.0),
    co2_gap_config=GapConfig(mode="simplified", h_eq=1.0e-20),
)
```

## 堆芯配置手册

### 四类映射参数

| 参数 | 含义 | 约束和影响 |
| --- | --- | --- |
| `tfe_multipliers` | 每个代表性 TFE 对应的真实热工/水力 TFE 数量 | 必须为正整数；参与功率守恒回推、慢化剂统计和反应性反馈统计 |
| `tec_multipliers` | 每个代表性 TFE 对应的真实 TEC 电路数量 | 默认复制 `tfe_multipliers`；允许为 0，但不能大于对应 `tfe_multipliers` |
| `tfe_power_factors` | 单根代表 TFE 获得的全堆功率份额 | 必须满足 `sum(tfe_power_factor * tfe_multiplier) == 1` |
| `ring_mapping` | 代表性 TFE 到全局慢化剂环编号的映射 | 用于将同一环内多个代表件按 `tfe_multipliers` 加权聚合 |

`ReactorCore.update_neutronic_power()` 下发的是“单根代表件功率”，不是整组功率。全堆功率通过 `代表件功率 * tfe_multipliers` 回推。

### 全局结构关系

```text
TFE 内部等效 moderator
  -> 按 ring_mapping 聚合到全局 moderator rings
  -> moderator rings 之间径向 SolidSolidCouple2D
  -> 可选 moderator-barrel gap
  -> 可选 barrel
  -> 可选 barrel-reflector gap
  -> 可选 reflector
  -> 最外层空间辐射边界
```

全局慢化剂环由 `mod_meshes`、`mod_material` 和 `ring_mapping` 共同启用。筒体和反射层分别使用 `GlobalAnnulusStructureConfig`；层间间隙使用 `GlobalGapStructureConfig`。外层间隙同样支持 `simplified` 和 `meshed` 两种模式。

### 最小堆芯装配

```python
from Components.ReactorCore import ReactorCore

core = ReactorCore(
    name="Core",
    tfe_dict={"TFE_Center": tfe},
    tfe_multipliers={"TFE_Center": 6},
    tec_multipliers={"TFE_Center": 6},
    tfe_power_factors={"TFE_Center": 1.0 / 6.0},
    enable_tec_coupled=False,
)
system.add_component(core)
```

启用 TEC 时，在加入系统后或初始化前配置电路模式：

```python
core.setup_tec_circuit(mode_str="fixed_u", target_value=V_target, I_guess=I_guess)
```

### 点堆推进和断点续算顺序

正常时间步由 `SystemManager` 统一调度：

```text
component.save_step_state()
  -> component.pre_step(dt, t_start)
  -> ReactorCore.advance_neutronics(dt, ...)
  -> 固体、流体和耦合器求解
  -> ReactorCore.commit_neutronics()
  -> component.post_step(dt, t_end)
```

首次启用点堆时调用 `core.initialize_point_reactor(total_power_initial=...)`。它会建立反馈参考态，并立即将初始功率下发到各个 TFE。长时间瞬态使用 `SystemManager.save_global_state()` 和 `SystemManager.load_global_state()`；底层会调用组件的 `get_state_dict()` 和 `load_state_dict()`。

## 热管散热器配置手册

### 散热和吸热口径

`HPwithFin` 的翅片是每个冷凝段轴向切片上的一维准稳态降维模型，不是独立固体网格，也不会出现在 `get_solids()` 中。

| 口径 | 含义 |
| --- | --- |
| 裸壁辐射 | 热管绝热段和冷凝段裸露外壁向环境辐射 |
| 翅片辐射 | 降维翅片模型计算出的辐射散热 |
| 外热吸收 | 太阳、查表热流、反照率或地球红外进入壁面和翅片的功率 |
| 毛散热 | 裸壁辐射与翅片辐射之和 |
| 净散热 | 毛散热减去外热吸收 |

`RingHP` 的 `get_total_*_scaled()` 系列接口会按 `hp_multipliers` 放大到真实热管数量。

### `external_heat_config` 字段表

| 字段 | 默认值 | 作用 |
| --- | --- | --- |
| `use_embedded_table` | `False` | 使用内嵌 Fortran 查表热流 |
| `table_ids` | `1` | 查表编号，可为标量或可广播到边界形状的数组 |
| `table_scale_factor` | `1.0` | 查表值缩放系数 |
| `table_offset` | `0.0` | 查表值偏置，单位 W/m2 |
| `table_periodic` | `True` | 是否按查表周期循环 |
| `add_solar` | `False` | 增加解析太阳热流 |
| `solar_constant` | `1361.0` | 太阳常数，单位 W/m2 |
| `orbit_height` | `800.0` | 轨道高度，单位 km |
| `orbit_period` | `7644.0` | 轨道周期，单位 s |
| `orbit_inclination` | `0.0` | 轨道倾角，单位 deg |
| `surface_normal_angles` | `(0.0, 0.0)` | 表面法向角 |
| `add_albedo` | `False` | 增加简化地球反照率热流 |
| `albedo_factor` | `0.3` | 简化反照率系数 |
| `add_earth_ir` | `False` | 增加简化地球红外热流 |
| `earth_ir_flux` | `237.0` | 地球红外热流密度，单位 W/m2 |
| `wall_illumination_factor` | `0.5` | 冷凝段裸壁受照面积系数 |
| `fin_illuminated_area_scale` | `1.0` | 翅片受照面积缩放系数 |
| `fin_loading_mode` | `"lumped_root_area"` | 翅片外热加载模式 |
| `node_configs` | 无 | 按代表热管节点覆盖统一配置，可用字典或与流体节点等长的列表 |
| `*_by_node` | 无 | 对单个字段逐节点覆盖；支持下表字段 |

`*_by_node` 支持：`table_ids`、`surface_normal_angles`、`solar_constant`、`orbit_height`、`orbit_period`、`orbit_inclination`、`albedo_factor`、`earth_ir_flux`、`wall_illumination_factor`、`fin_illuminated_area_scale`、`fin_loading_mode`、`table_scale_factor`、`table_offset`、`table_periodic`、`add_solar`、`add_albedo`、`add_earth_ir`、`use_embedded_table`。

当 `use_embedded_table=True` 时，当前实现优先使用查表热源并直接返回，不再叠加 `add_solar`、`add_albedo` 和 `add_earth_ir`。

### 翅片加载模式

| 模式 | 行为 | 适用场景 |
| --- | --- | --- |
| `lumped_root_area` | 将壁面受照面积和翅片受照面积合并后挂到冷凝段根部边界 | 默认兼容模式 |
| `distributed_fin_absorption` | 壁面外热仍挂在冷凝段边界；翅片外热进入降维翅片求解器 | 需要让翅片直接吸收外热时 |

### 查表热流和节点覆盖示例

统一使用 Fortran 查表热流：

```python
external_heat_config = {
    "use_embedded_table": True,
    "table_ids": 1,
    "table_scale_factor": 1.0,
    "table_offset": 0.0,
    "table_periodic": True,
    "wall_illumination_factor": 0.5,
    "fin_illuminated_area_scale": 1.0,
    "fin_loading_mode": "lumped_root_area",
}
ring_hp = RingHP(..., external_heat_config=external_heat_config)
```

切换为翅片直接吸热，并让不同代表热管使用不同查表编号：

```python
external_heat_config = {
    "use_embedded_table": True,
    "table_ids_by_node": [1, 2, 3, 4],
    "fin_loading_mode": "distributed_fin_absorption",
    "fin_illuminated_area_scale_by_node": [1.0, 1.0, 0.8, 0.8],
    "node_configs": {
        0: {"wall_illumination_factor": 0.6},
        3: {"table_scale_factor": 0.9},
    },
}
ring_hp = RingHP(..., external_heat_config=external_heat_config)
```

列表长度必须与集流环流体节点数一致。也可以使用以节点编号为键的字典，仅覆盖指定节点。

## 全局数据流

```text
PointReactor -> ReactorCore -> TFEUnit.update_neutronic_power()
ThermoCalcModel -> electric field -> Joule heat -> Electrode.Q_source
ThermoCalcModel -> plasma flux -> TECCouple2D
ExternalHeatSources -> ExternalHeatFluxBC -> BoundaryRegion
SystemManager -> pre_step -> solve -> post_step
```

TEC 热流和焦耳热进入导热方程的路径不同：

- 面热流：`ThermoCalcModel` 生成发射极失热和接收极得热的 W/m2 数组，`TECCouple2D` 将其作为极间边界热流。
- 焦耳热：电势或电场经 [`tec_electric.py`](./tec_electric.py) 转为每节点 W，随后 `Electrode.set_joule_heating()` 写入 `Q_source`。

## 单位和符号约定

| 数据 | 单位 | 约定 |
| --- | --- | --- |
| 温度 | K | 所有热工组件统一使用 |
| 节点热源 `Q_source` | W | 每个节点的总功率 |
| 体积热源 `Q_vol` | W/m3 | 用于物理计算或后处理 |
| 外热流密度 | W/m2 | 正值表示流入边界 |
| 等离子体发射极热流 | W/m2 | 负值表示发射极失热 |
| 等离子体接收极热流 | W/m2 | 正值表示接收极得热 |
| TEC 热流面积基准 | m2 | 统一按发射极外表面积 |
| 电场 | V/m | 非均匀网格必须按节点中心坐标求导 |
| 电阻率 | Ohm*m | 焦耳热计算使用 |

## 修改代码前必读

- 不要直接修正 [`basicComponents/Electord.py`](./basicComponents/Electord.py) 的文件名，否则会破坏现有导入路径。
- TEC 焦耳热优先复用 [`tec_electric.py`](./tec_electric.py)，避免重新实现非均匀网格梯度逻辑。
- 外热源模型返回 W/m2；`ExternalHeatFluxBC` 负责乘面积转换为 W。调用方不要再次乘面积。
- `ReactorCore` 倍率会影响功率、慢化剂热源、TEC 热源和反应性反馈统计。修改倍率前先核验四类映射参数。
- [`basicComponents/Electord.py`](./basicComponents/Electord.py) 中的 `Collector.get_reactivity_feedback()` 当前仍复用发射极反馈多项式，并保留待完善项。`ReactorCore` 自身的堆芯级反馈公式已经分开处理两类电极；修改反馈时必须同时核对两处。
- 严格反照率和地球红外积分模型尚未实现，不应作为生产模型调用。
- [`basicComponents/FinConduction.py`](./basicComponents/FinConduction.py) 是独立二维翅片求解器；`HPwithFin` 使用的是内部降维准稳态翅片模型。两者用途不同，修改时不得混淆。

## Codex 快速定位

| 修改目标 | 优先阅读 |
| --- | --- |
| 宏观组件契约、生命周期 | [`BaseComponent.py`](./BaseComponent.py) |
| 单壁或环形流道组装 | [`Pipe.py`](./Pipe.py)、[`AnnularPipe.py`](./AnnularPipe.py) |
| 燃料功率分配、反馈 | [`basicComponents/Fuel.py`](./basicComponents/Fuel.py) |
| 电极焦耳热、TEC 面热流 | [`basicComponents/Electord.py`](./basicComponents/Electord.py)、[`basicComponents/TECPair.py`](./basicComponents/TECPair.py)、[`tec_electric.py`](./tec_electric.py) |
| 多个独立 `TECPair` 的 TEC 电路管理 | [`TECCircuitManager.py`](./TECCircuitManager.py) |
| TFE 内部层级和气隙 | [`TFEUnit.py`](./TFEUnit.py) |
| 堆芯倍率、慢化剂、点堆 | [`ReactorCore.py`](./ReactorCore.py) |
| 热管导热和物性冻结 | [`basicComponents/HeatPipe2D.py`](./basicComponents/HeatPipe2D.py) |
| 独立二维翅片求解 | [`basicComponents/FinConduction.py`](./basicComponents/FinConduction.py) |
| 热管内部降维翅片散热 | [`HPwithFin.py`](./HPwithFin.py) |
| 集流环和代表热管阵列 | [`RingHP.py`](./RingHP.py) |
| 轨道外热流 | [`ExternalHeatSources/__init__.py`](./ExternalHeatSources/__init__.py) |
| Fortran 热流查表 | [`ExternalHeatSources/embedded_flux_tables.py`](./ExternalHeatSources/embedded_flux_tables.py) |

## 典型集成顺序

1. 构造基础材料、网格、流体通道和换热关联式。
2. 构造 `TFEUnit`、`Pipe`、`AnnularPipe`、`HPwithFin` 或 `RingHP`。
3. 将多个 TFE 组织进 `ReactorCore`，配置真实数量倍率和全局结构。
4. 将宏观组件交给 `SystemManager`，由其调用 `get_solids()` 和 `get_couplers()`。
5. 每个时间步：
   - `SystemManager` 调用组件 `pre_step()`。
   - 固体、流体和耦合器推进。
   - `SystemManager` 调用组件 `post_step()`。
6. 按需读取散热、外热吸收、反应性反馈、电热耦合结果或断点状态。

## 使用注意事项

- 根目录组件是装配层，底层导热和边界细节主要来自 `basicComponents` 与 `Solvers`。
- `get_solids()` 和 `get_couplers()` 是系统集成的关键接口，新增组件时必须保持返回对象完整。
- 外热流应使用 `ExternalHeatSources` 中的热源模型，并注意 W/m2 与 W 的转换边界。
- TEC 焦耳热映射应优先使用 `tec_electric.py` 中支持非均匀网格的函数。
- `ReactorCore` 中的倍率会影响功率分配、慢化剂热源转移和总反馈统计，设置时应与物理代表单元一致。
- `HPwithFin` 的翅片是准稳态降维模型，不会作为独立固体出现在 `get_solids()` 中。
- `TFEUnit` 和 `ReactorCore` 都提供状态保存接口，长时间瞬态计算应优先使用这些接口实现断点续算。

## 2026-06-02 TEC 焦耳热映射更新

TEC 生产路径不再用节点电势的 `np.gradient()` 结果重新推导焦耳热。当前权威输入是 ThermoCalc C++ `VcalcFVM()` 输出的逐轴向节点功率：

```text
joulePowerE / joulePowerC [W]
    -> distribute_axial_power_by_volume()
    -> Electrode.set_joule_heating()
    -> Electrode.Q_source [W]
```

`TFEUnit.update_joule_power_sources()` 和 `TECPair.set_joule_heating_axial_power()` 按每个轴向列内的二维控制体体积比例分配功率，并保持每列总功率不变。旧电场接口继续保留，用于兼容和诊断，不再作为 TEC 生产热源的权威路径。

## 2026-06-02 全局慢化剂映射时间层修复

`ReactorCore.pre_step()` 中的全局慢化剂源项转移必须按以下顺序执行：

```text
TFEUnit.pre_step()
  -> 更新内部等效 moderator 的外边界温度
  -> 刷新内部 moderator 的物性、热阻、边界状态和热流缓存
  -> 读取 moderator 外流
  -> 按 tfe_multipliers 聚合到全局 moderator rings
```

不得在 `TFEUnit.pre_step()` 之前读取内部 moderator 的 `BoundaryRegion.current_flux`。旧顺序会在正常推进中引入一步滞后，并在 restart 重建后把旧边界缓存作为首步源项注入全局慢化剂环。
