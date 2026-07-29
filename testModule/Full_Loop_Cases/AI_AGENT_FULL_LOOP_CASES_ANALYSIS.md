# Full_Loop_Cases Codex Notes

`testModule/Full_Loop_Cases/` contains the new shared TOPAZ-II full-loop builder layer.

Rules:

- Do not import legacy case builders from `test_core_assemble_v7_caseA.py`, `test_core_assemble_v8_caseA.py`, `test_core_assemble_v11_caseA.py`, `test_core_assemble_v12_caseA.py`, or `test_core_assemble_v13_caseA.py`.
- Shared builders may import only lower-level project modules such as `Components`, `Materials`, and `Solvers`.
- The shared layer builds the common closed-loop skeleton: reactor core, representative TFE coolant channels, core inlet/outlet connectors, radiator inlet/outlet header interfaces, and two series pumps.
- V14 and V15 radiator implementations should attach through the common radiator interface rather than modifying the common core builder.
- `close_with_placeholder_bridge=True` is only for topology tests and closes `RadiatorInletHeader` to `RadiatorOutletHeader` with a simple hydraulic junction.
- `connect_pump_outlet_to_core=False` keeps `PumpOutletNode` detached from the common `CoreInletSegment`. Use this for cases such as V15 where the pump outlet must pass through case-specific distributor and cold return branches before reaching `CoreInletConnector`.

Default representative TFEs:

```text
Center, Ring1, Ring2, Ring3_TEC, Ring3_Open
```

Default multipliers:

```text
ring_multipliers = 1, 6, 12, 15, 3
tec_ring_multipliers = 1, 6, 12, 15, 0
```

`tec_ring_multipliers` is configurable, including `Ring3_Open=3`, so later V14/V15 runners can attach the three reserved outer-ring TECs through a separate parallel circuit.

## V14 heat-pipe radiator case

V14 是热管辐射器 TOPAZ-II 闭式回路算例。它不导入 V11 或 CoolantLoop 的历史算例 builder，只在 `Full_Loop_Cases` 中复用共用层对象并本地声明热管辐射器部件。

Implemented entry points:

- `v14_heatpipe_radiator.py`: 本地声明 V14 热管辐射器 adapter，构建三条热端出口支路、三个入口混合节点、六段显式单环 `RingHP` 集流环、三个出口混合节点和三条 manifold。
- `v14_case.py`: 暴露 `build_v14_case_a_system(...)` 和 `v14_basic_diagnostics(...)`。
- `V14_run_cases/run_v14_flow_path_smoke.py`: 水力-only smoke，关闭辐射发射率，只初始化闭式流网并推进一次水力步。

Common-layer boundary:

- 共用层提供 `CoreInletConnector`、`CoreOutletConnector`、`RadiatorInletHeader`、`RadiatorOutletHeader`、`CoreInletSegment` 和两台串联主泵。
- V14 的 `HotOutletBranch_1/2/3`、入口混合节点、热管集流环、出口混合节点和 `Manifold_1/2/3` 均由 V14 adapter 本地构建。
- 当前 V14 没有把三进三出连接段整体放入共用层；只有 `CoreInletSegment` 是共用层里的泵出口到堆芯入口短连接段。

Topology notes:

- `RadiatorInletHeader -> HotOutletBranch_1/2/3 -> InletMix_I1/I2/I3`。
- 集流环顺序为 `I1 -> A1 -> O1 -> A2 -> I2 -> A3 -> O2 -> A4 -> I3 -> A5 -> O3 -> A6 -> I1`。
- `OutletMix_O1/O2/O3 -> Manifold_1/2/3 -> RadiatorOutletHeader`。
- 热端支路到集流环、以及 manifold 到出口总管使用 `MacroFlowJunction(multiplier=2)`，因为 V14 只显式建一套物理集流环，并用倍率表示第二套对称集流环。

Current validation:

- `testModule/test_v14_caseA_topology.py` 检查 V14 拓扑、macro multiplier、泵出口到堆芯入口分离、压力参考策略、诊断字段和禁止历史依赖。
- `testModule/Full_Loop_Cases/V14_run_cases/test_v14_flow_path_smoke.py` 检查水力-only smoke 结果。

## 2026-07-15 V14 热管辐射器外热与角系数更新

V14 默认仍关闭直接辐射器外热；显式设置 `V14HeatPipeRadiatorConfig(external_heat_enabled=True)` 后，N18 表加载到 18 个代表热管区。当前 V14 的热管管壁和翅片共用以下面角系数：外侧上/下 `1.0/1.0`，内侧上/下 `0.0/0.3`。因此 `F_outer=1.0`、`F_inner=0.3`、双面等效 `F_face=0.65`，默认 `epsilon=0.75` 时有效辐射发射率为 `0.4875`。

外热吸收与辐射排热分开：外侧管壁和翅片分别使用外侧受照面积，并乘 `external_heat_absorption_efficiency=0.992` 与各自表面发射率；内侧不吸收外热。默认 `distributed_fin_absorption` 中，管壁有效面积进入一个边界热流条件，翅片单侧有效面积进入降维翅片方程。代表热管的 `hp_multipliers` 以及第二套对称环倍率只在汇总层使用，不能再次乘到单根吸收面积。

本轮没有修改 `Components/basicComponents/HeatPipe2D.py`；它继续作为热管内部二维求解器接口。相关单元覆盖位于 `testModule/test_v14_caseA_topology.py`，包括默认 V14 角系数、有效发射率、外热吸收面积、效率和 RingHP 倍率口径。

## V15 pipe-fin radiator case

V15 是管翅式辐射管辐射器 TOPAZ-II 闭式回路算例。它不照搬 V13/V12 的 `Pipe05/06/07/08/09/Pipe11` 命名和拆分方式，而是按物理语义重建辐射器、总管、泵出口分配和堆芯入口回流支路。

Implemented entry points:

- `v15_pipefin_radiator.py`: 本地声明 V15 管翅式辐射器 adapter，构建入口分配器、78 段上集流环、78 根 `RadiatorTubeFluid` 与 `RadiatorPipeWithFin`、78 段下集流环、辐射器内外总管、泵出口分配器和三条冷回流支路。
- `v15_case.py`: 暴露 `build_v15_case_a_system(...)` 和 `v15_basic_diagnostics(...)`。
- `V15_run_cases/run_v15_flow_path_smoke.py`: 水力-only smoke，关闭管壁和翅片发射率，只初始化闭式流网并推进一次水力步。

Common-layer boundary:

- V15 使用 `build_full_loop_common_base(..., connect_pump_outlet_to_core=False)`，因此共用层主泵段仍然存在，但不会生成 `CoreInletSegment`，也不会把 `PumpOutletNode` 直接接到 `CoreInletConnector`。
- V15 adapter 负责连接 `PumpOutletNode -> PumpOutletDistributor -> ColdReturnBranch_1/2/3 -> CoreInletConnector`。
- V15 adapter 会把共用层的 `RadiatorOutletHeader` 重命名为 `RadiatorOuterHeader`，以匹配物理语义。

Topology notes:

- `CoreOutletConnector -> RadiatorInletHeader -> RadiatorInletDistributor`。
- `RadiatorUpperHeader_01...78` 和 `RadiatorLowerHeader_01...78` 都包含完整环向流动连接。
- `RadiatorUpperHeader_XX -> RadiatorTubeFluid_XX -> RadiatorLowerHeader_XX`，每根管挂接一个 `RadiatorPipeWithFin` 固体/翅片等效模型。
- `RadiatorInnerHeader -> RadiatorOuterHeader -> PumpA -> PumpMidNode -> PumpB -> PumpOutletNode -> PumpOutletDistributor`。
- `PumpOutletDistributor -> ColdReturnBranch_1/2/3 -> CoreInletConnector`。

Dimension notes:

- V14 的 `HotOutletBranch_1/2/3` 默认 `length=2.19632 m`、`area=pi*0.0138^2 m2`、`Dh=0.0276 m`、`n_nodes=8`。
- V15 的 `ColdReturnBranch_1/2/3` 默认 `length=1.89021 m`、`area=pi*0.0138^2 m2`、`Dh=0.0276 m`、`n_nodes=1`。
- 两者截面积和水力直径相同，但长度和节点数不同；后续若确认三进三出连接段应共用，应把相关参数上提到共用层。

Naming constraints:

- V15 不得重新引入 `V12_` 前缀。
- V15 不得重新引入 `Pipe05`、`Pipe06`、`Pipe07`、`Pipe08`、`Pipe09`、`Pipe11`。

Current validation:

- `testModule/test_v15_caseA_topology.py` 检查 V15 对象数量、完整 78 管/上下集流环拓扑、泵出口分离、压力参考策略、诊断字段和禁止历史导入/命名。
- `testModule/Full_Loop_Cases/V15_run_cases/test_v15_flow_path_smoke.py` 检查水力-only smoke 结果，包括收敛、有限压力/流量、唯一被动压力参考、无固定压力边界和 V15 计数。

## 2026-07-13 V14/V15 orbital external heat

V14 and V15 default to no external heat. Set `external_heat_enabled=True` explicitly to load the periodic w0=8.12 deg, i_s=58.5 deg tables.

- V14 reads Components/ExternalHeatSources/is58p5_w0_8p12_N18_sum.csv. Its six ring sectors each contain three representative heat-pipe zones mapped in order to theta_000...theta_017. The second symmetric ring represented by symmetric_ring_multiplier=2 uses the same circumferential heat history.
- V15 reads Components/ExternalHeatSources/is58p5_w0_8p12_N78_sum.csv. theta_000...theta_077 map one-to-one to RadiatorTube_01...78; each tube applies its value uniformly over eight axial nodes.
- V14 defaults to `distributed_fin_absorption`: wall heat uses 0.5 of bare-wall area and fin heat enters the reduced fin equation through one illuminated face.
- V15 默认 `external_heat_fin_loading_mode=distributed_fin`；翅片使用单侧几何投影面积，旧 `lumped_root_area` 仍可显式选择。外热吸收统一乘 `external_heat_absorption_efficiency=0.992`，并按 Kirchhoff 关系分别取 `alpha_tube=epsilon_tube`、`alpha_fin=epsilon_fin`。
- `distributed_fin` 只把裸管外侧有效吸收面积挂到管壁边界，并将同一个 N78 热流密度作为源项加入翅片准稳态方程。每根管各轴向段使用相同热流密度，各段只乘本段面积；有效吸收面积为 `A_projected * illumination * 0.992 * epsilon_surface`。
- 分布模式的总外热输入应通过 RadiatorPipeWithFin.get_external_heat_absorption_distribution() 统计；external_heat_bc.q_flux 只代表管壁部分。
- V15 全局能量审计的向空间排热项使用 gross_rejection；外热输入在审计式中单独加入，不能再使用已扣除翅片吸热的 net_rejection。
- V15 管壁和翅片共用显式内外面角系数：外面上/下 1/1，内面上/下 0/0.17；一维等效为 F_outer=1、F_inner=0.17、双面平均 F=0.585，所以默认表面发射率 0.8 对应有效辐射发射率 0.468。
- 0.992 只作用于入射外热吸收，不作用于表面向外自发辐射。


The sum files already contain direct solar, Earth IR, and Earth-reflected heat. Do not add those components again. V15 applies one surface absorptivity (`alpha=epsilon`) to the combined table because its three components are not separated at runtime.

## V15 V71 center-0.30 m uniform-heating variant

`V15_run_cases_V71` is a V15 pipe-fin radiator case copy whose loop topology and hydraulic-only smoke settings remain the same as `V15_run_cases`. Its core-profile difference is the axial power profile: `build_v15_v71_case_a_system(...)` wraps `build_v15_case_a_system(...)` and replaces every representative TFE fuel profile with `center_0p30m_uniform` before system initialization. V71 overrides the case-local default total pump head to `20200 Pa` (`10100 Pa` per series pump) and sets both radiator-tube inlet/outlet loss coefficients to `K=38` for an actual total flow near `1.18 kg/s`; the `1.3 kg/s` flow setting remains the design/initial-flow guess.

The V71 profile is generated by `build_center_uniform_axial_power_profile(6, 25, 6, 0.30)`: within the full TFE axial length `0.065 + 0.377 + 0.065 = 0.507 m`, it centers a `0.30 m` heater interval, allocates power by each axial cell's overlap length with that interval, and normalizes the profile to sum to `1.0`. Ring power shares, representative multipliers, TEC multipliers, V15 radiator topology, and cold-return branches are unchanged. V71 uses a case-local default total pump head of `20200 Pa` (`10100 Pa` per series pump) and radiator-tube inlet/outlet loss coefficients `K=38`, calibrated against the warm coupled restart for a total loop flow of approximately `1.18 kg/s`.

Validation entry points:

- `testModule/test_v15_v71_core_profile.py`: checks V71 profile construction and application to all representative TFEs.
- `testModule/Full_Loop_Cases/V15_run_cases_V71/test_v15_v71_flow_path_smoke.py`: V71 hydraulic-only smoke and profile metadata check.

## 2026-07-10 Full_Loop_Cases_10kW branch

`testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/` is a copied V14 heat-pipe-radiator full-loop package for the 10 kW / five-ring core work. It is intentionally local to that directory: avoid changing the shared `Full_Loop_Cases` builders unless the same behavior is required by both the original V14/V15 and the 10 kW branch.

Key differences from the original V14 package:

- Core representative rings are `Center, Ring1, Ring2, Ring3, Ring4` with physical counts `1, 6, 9, 18, 24`.
- Main TEC series voltage is in the 10 kW tuning range around `50.5 V`, not the old 37-TFE `27.2 V` setting.
- The V14 heat-pipe radiator explicitly builds both upper and lower collector rings. It does not use the old single-ring plus `multiplier=2` shortcut.
- Heat-pipe counts are `154` upper and `186` lower. The design hydraulic split follows `154/340` and `186/340`, but the transient flow is still solved by the hydraulic network.
- Heat-pipe geometry is `L_eva=0.0605 m`, `L_aba=0.0415 m`, `L_con=0.47 m`; fins attach only to `L_con`.
- `V14HeatPipeRadiatorConfig` supports separate lower-side view factors for upper/lower rings: `upper_hp_down_view_factor` and `lower_hp_down_view_factor`. The current tuning path keeps `hp_up_view_factor=0.0`.

Powered debug runner:

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_debug/run_v14_210kw_debug.py
```

This runner uses fixed core power (`210000 W`), fixed target loop flow (`2.46 kg/s`), implicit-Euler solid conduction, disabled orbital external heat flux, and optional TEC lookup. Use direct runner invocations for this workflow; do not use `python -m unittest` as the debug entry.

Current best powered baseline:

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/V14_210kW_debug/runs/final_eps07475_u50p65_wire0335_1200s_from7964
restart = stage_01_restart.npz
radiator_emissivity = 0.7475
tec_voltage = 50.65 V
wire_resistance_scale = 0.335
upper/lower hp down view factor = 0.3 / 0.4
```

Endpoint at `t=9163.85 s`:

```text
core inlet/outlet = 754.738 / 845.773 K
TEC current = 206.569 A
TEC net electric power = 10.463 kW
total radiator rejection = 195.212 kW
required pump head = 27.880 kPa
```

This is close to the current targets (`754.45 K`, `845.65 K`, `206 A`, `10.44 kW`) but should be treated as the working baseline, not a final exact steady calibration. Local details are documented in `Full_Loop_Cases_10kW/README.md` and `Full_Loop_Cases_10kW/V14_210kW_debug/TUNING_LOG.md`.
The case-local `V15_run_cases_V71/run_v15_v71_cold_start.py` supports a `723 K`, `106 kW`, `4 K`-space cold start with TEC disabled and `implicit_euler` solid conduction, followed by a separate one-shot `27.2 V` lookup-based TEC solve from the saved `1000 s` restart. The thermal and electrical stages are separate processes so a stalled ThermoCalc call cannot destroy the thermal restart. Use the continuous thermal history, not the immediately reloaded radiator diagnostic, for final pipe-fin rejection because the reduced-order fin warm-start arrays are not restart state.

## 2026-07-13 V14/V15 pre-start cooldown

`Prestart_run_cases/run_prestart_cooldown.py` runs original V14 or V15 from a uniform `350 K` state with zero core power, disabled TEC/point kinetics, fixed `0.26 kg/s` loop flow, and a `4 K` space background. It stops when the minimum coolant temperature reaches `260 K` or after `10000 s`.

- Direct radiator external heat is disabled. V14 N18 and V15 N78 values are averaged into six circumferential shield sectors and applied only to `RadiatorThermalShield`.
- The shield uses `segment_balance`, inner/outer emissivity `0.8/0.1`, conductivity `0.0008 W/(m K)`, thickness `0.01 m`, and view factor `0.8`.
- V14 uses its 18 representative `HPwithFin` units and scales shield area by heat-pipe count and the two-ring symmetry multiplier.
- Complete fluid arrays and every registered solid temperature field are written every `1 s` in compressed 60-s NPZ chunks. CSV/JSON diagnostics use the same interval; global restart files are saved every `60 s`.
- Since 2026-07-14 this runner switches all capacitance-backed `FluidSolidCouple` objects to `local_implicit` after initialization/restart loading. It calls `compute_adaptive_dt(..., respect_fluid_cfl=False)` because `step_Picard()` advances enthalpy with the fully implicit sparse solve; other runners retain the default CFL behavior.
- Restart continuations passed staged `max_dt=0.02 s`, `0.05 s`, and `0.10 s` checks for both V14 and V15. The active long-run convention is `max_dt=0.10 s`, with the existing 1.2x growth limiter smoothing each restart transition.
- Fixed-flow operation constrains PumpA only. PumpB remains a normal pressure-head junction; constraining both serial pumps leaves `PumpMidNode` pressure underdetermined.

The same runner also supports direct radiator heating with `--external-heat-target radiator`. In this mode no `RadiatorThermalShield` is created: V14 maps the N18 table to its 18 representative heat-pipe zones, while V15 maps N78 one-to-one to the 78 pipe-fin units. The 360 CSV samples span the physical `6552 s` orbit; use `--duration 13104 --stop-temperature 260.55` for two complete external-heat periods terminated early only at the NaK-78 freezing threshold. Diagnostics include the six aggregated circumferential fluxes and `radiator_external_heat_W` in addition to the full 1-s state history.

V14 radiator power diagnostics use the complete two-ring system basis. `radiator_rejection_W` is the sum of double-ring heat-pipe/fin rejection and double-ring collector-wall radiation; the two contributions are also recorded as `radiator_heatpipe_rejection_W` and `radiator_ring_wall_rejection_W`. Temperatures and local heat fluxes remain representative-unit values and are not multiplied.

`run_prestart_cooldown.py` accepts `--initial-temperature` for uniform fluid and solid initialization. Its w0=8.12 external-heat tables treat the 1...360 CSV coordinate as orbital degrees and scale it to a `6552 s` physical period. Direct V14/V15 radiator heating uses table scale `1.0`; each radiator model then applies its own `0.992 * 0.8` absorption multiplier. Shield heating applies only the shield outer-surface absorptivity `alpha=epsilon=0.1` and does not apply `0.992`. The 1-s state chunks include representative fin temperature, radiation, absorption, root heat flow, and radiation-background arrays; these quasi-steady fin fields are intentionally not added to restart files. `run_temperature_batch.py` runs the requested 15-case 300-340 K matrix: V14/V15 direct-radiator external heat plus V15 shield-only external heat, all with a `12000 s` limit, `260.55 K` coolant stop threshold, 1-s full-state records, and 60-s restart files.
