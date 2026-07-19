# V13 CaseA Flow Network and Solid Coupling Map

> Scope: current `testModule/test_core_assemble_v13_caseA.py` default V13 CaseA closed-loop model and the V13 runners that use it.
> The names below intentionally match the model object names in the current code. Dimensions are taken from the current source defaults. Do not infer missing outer diameters for hydraulic-only channels.

Related figure: [`V13_FLOW_NETWORK_AND_SOLID_COUPLING.svg`](./V13_FLOW_NETWORK_AND_SOLID_COUPLING.svg)

## 1. Model Identity

V13 CaseA is built by `build_v13_case_a_system(...)`.

Default case string:

```text
v13_closed_core_pipefin_radiator_pumped_loop
```

Default closed-loop counts from a direct build:

| Item | Count |
| --- | ---: |
| Hydraulic volumes | 1036 |
| Hydraulic junctions | 1122 |
| Registered components | 79 (`TASTIN_Core_V8_CaseA` + 78 radiator tube components) |
| Registered solids | 114 |
| Registered couplers | 113 |

Pressure boundary/reference status in the closed-loop build:

| Object | Status |
| --- | --- |
| `V12_CoreInletConnector` | `is_pressure_reference=True` |
| Fixed pressure boundary volumes | none |
| Fixed flow / temperature boundary volumes | none in V13 closed loop; V12 open-loop boundary volumes are removed |

## 2. Hydraulic Topology

The V13 closed-loop hydraulic path is:

```text
V13_PumpOutletNode
  -> Pipe11_CoreInletHeader
  -> V12_CoreInletDistribution
  -> V12_CoreInletBranch_1
     + V12_CoreInletBranch_2_3_Rep(multiplier=2)
  -> V12_CoreInletConnector
  -> Chan_Center / Chan_Ring1 / Chan_Ring2 / Chan_Ring3_TEC / Chan_Ring3_Open
  -> V12_CoreOutletConnector
  -> Pipe05_CoreOutletToRadiator
  -> V12_RadiatorInletSplit
  -> V12_RadiatorUpperHeader_01...78
  -> V12_RadiatorTubeFluid_01...78
  -> V12_RadiatorLowerHeader_01...78
  -> V12_RadiatorOutletMix
  -> Pipe06_RadiatorOutlet
  -> Pipe07_HeatExchangerHotSide
  -> Pipe08_ReturnInnerPipe
  -> Pipe09_ValveSegment
  -> V13_PumpMidNode
  -> V13_PumpOutletNode
```

### 2.1 Main Loop Pipes and Connector Volumes

Hydraulic-only channels do not have solid wall objects in V13. For those channels, the code defines `flow_area` and `hydraulic_diam`; an outer diameter is not part of the model.

| Object | Nodes | Length (m) | Flow Area (m2) | Hydraulic Diameter / Inner Diameter (m) | Outer Diameter in V13 | Solid Coupling |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `Pipe11_CoreInletHeader` | 8 | 0.130000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `V12_CoreInletDistribution` | volume | connector input length 0.020000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `V12_CoreInletBranch_1` | 8 | 1.890210 | 5.979816e-4 | `D_h=0.027600` | not modeled | none |
| `V12_CoreInletBranch_2_3_Rep` | 8 | 2.507050 | 5.979816e-4 | `D_h=0.027600` | not modeled | none; macro multiplier = 2 |
| `V12_CoreInletConnector` | volume | connector input length 0.020000 | 2.050455e-3 | `D_h=0.001400` | not modeled | connected to core channel junctions |
| `V12_CoreOutletConnector` | volume | connector input length 0.020000 | 2.050455e-3 | `D_h=0.001400` | not modeled | connected to core channel junctions |
| `Pipe05_CoreOutletToRadiator` | 8 | 0.130000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `V12_RadiatorInletSplit` | volume | connector input length 0.020000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `V12_RadiatorOutletMix` | volume | connector input length 0.020000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `Pipe06_RadiatorOutlet` | 8 | 0.043408 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `Pipe07_HeatExchangerHotSide` | 8 | 0.005426 | 3.800000e-4 | `D_h=0.067900` | not modeled | none |
| `Pipe08_ReturnInnerPipe` | 8 | 0.130000 | 3.800000e-4 | `D_h=0.047000` | not modeled | none |
| `Pipe09_ValveSegment` | 8 | 0.130000 | 3.800000e-4 | `D_h=0.047000` | not modeled | none |
| `V13_PumpMidNode` | volume | connector input length 0.020000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |
| `V13_PumpOutletNode` | volume | connector input length 0.020000 | 3.800000e-4 | `D_h=0.014000` | not modeled | none |

`Pipe07_HeatExchangerHotSide` is retained as a hydraulic flow-network segment. The current V13 model does not attach a heat-exchanger solid or heat-transfer boundary to it.

### 2.2 Core Coolant Channels

All core coolant channels use the same non-uniform axial mesh and annular coolant geometry.

| Object | Representative Multiplier | TEC Multiplier | Nodes | Length (m) | Flow Area (m2) | Hydraulic Diameter (m) | Solid Coupling |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `Chan_Center` | 1 | 1 | 37 | 0.507000 | 5.541769e-5 | 0.001400 | `Center_iclad_coolant_couple`, `Center_oclad_coolant_couple` |
| `Chan_Ring1` | 6 | 6 | 37 | 0.507000 | 5.541769e-5 | 0.001400 | `Ring1_iclad_coolant_couple`, `Ring1_oclad_coolant_couple` |
| `Chan_Ring2` | 12 | 12 | 37 | 0.507000 | 5.541769e-5 | 0.001400 | `Ring2_iclad_coolant_couple`, `Ring2_oclad_coolant_couple` |
| `Chan_Ring3_TEC` | 15 | 15 | 37 | 0.507000 | 5.541769e-5 | 0.001400 | `Ring3_TEC_iclad_coolant_couple`, `Ring3_TEC_oclad_coolant_couple` |
| `Chan_Ring3_Open` | 3 | 0 | 37 | 0.507000 | 5.541769e-5 | 0.001400 | `Ring3_Open_iclad_coolant_couple`, `Ring3_Open_oclad_coolant_couple` |

The TFE axial layout is:

| Region | Length (m) | Nodes |
| --- | ---: | ---: |
| lower | 0.065 | 6 |
| active | 0.377 | 25 |
| upper | 0.065 | 6 |
| total | 0.507 | 37 |

### 2.3 Pipe-Fin Radiator Hydraulic Geometry

The radiator contains 78 identical tube fluid channels and two one-node-per-segment header rings.

| Object Pattern | Count | Nodes Each | Length Each (m) | Flow Area (m2) | Diameter in Code | Solid Coupling |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| `V12_RadiatorUpperHeader_01...78` | 78 | 1 | 0.0331881070 | 3.141593e-4 | header inner diameter = 0.020000 m | none |
| `V12_RadiatorTubeFluid_01...78` | 78 | 8 | 1.850000 | 3.848451e-5 | tube inner diameter = 0.007000 m | `V12_RadiatorTube_XX_FluidSolid` |
| `V12_RadiatorLowerHeader_01...78` | 78 | 1 | 0.0542126117 | 3.141593e-4 | header inner diameter = 0.020000 m | none |

Radiator tube solid wall geometry:

| Object Pattern | Count | Axial Nodes | Radial Wall Nodes | Tube Length (m) | Tube Inner Diameter (m) | Tube Outer Diameter (m) | Wall Thickness (m) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `V12_RadiatorTube_01...78_Wall` | 78 | 8 | 1 | 1.850000 | 0.007000 | 0.008000 | 0.000500 |

Radiator fin reduced-order geometry:

| Field | Value |
| --- | ---: |
| `fin_thickness_m` | 0.000400 |
| `n_fin_width` | 12 |
| `fin_width_upper_m` | 0.033190 |
| `fin_width_lower_m` | 0.054210 |
| `fin_area_scale` | 0.35 |
| `fin_view_factor` | 1.0 |
| `fin_conductivity_w_m_k` | 348.9 |
| `fin_contact_resistance_m2k_w` | 0.0 |

The header segment outer diameter is not modeled. Header segments are hydraulic-only volumes/junctions; there is no header solid conduction object and no header-to-space radiation boundary in the V13 pipe-fin radiator path.

### 2.4 Pump Segment

The V13 pump is represented by two serial pump junctions with one intermediate fluid volume.

| Junction | From | To | Type | Default Head (Pa) | Default Flow Target |
| --- | --- | --- | --- | ---: | ---: |
| `J_Pipe09_to_V13_PumpA` | `Pipe09_ValveSegment_Vol_08` | `V13_PumpMidNode` | `FlowControlledPumpJunction` by default | 3950 | 1.3 kg/s |
| `J_V13_PumpA_to_PumpB` | `V13_PumpMidNode` | `V13_PumpOutletNode` | `FlowControlledPumpJunction` by default | 3950 | 1.3 kg/s |
| `J_V13_PumpOutlet_to_Pipe11` | `V13_PumpOutletNode` | `Pipe11_CoreInletHeader_Vol_01` | `FlowJunction` | none | loop flow |

Default total pump head is `7900 Pa`. Each pump receives half of the total head.

## 3. Main Junction Groups

### 3.1 Core Inlet and Outlet Junctions

| Junction | From | To | Notes |
| --- | --- | --- | --- |
| `J_V13_PumpOutlet_to_Pipe11` | `V13_PumpOutletNode` | `Pipe11_CoreInletHeader_Vol_01` | closes pumped return into cold feed |
| `J_Pipe11_to_CoreInletDistribution` | `Pipe11_CoreInletHeader_Vol_08` | `V12_CoreInletDistribution` | connector length 0.02 m |
| `J_CoreInletDistribution_to_CoreInletBranch_1` | `V12_CoreInletDistribution` | `V12_CoreInletBranch_1_Vol_01` | one physical branch |
| `J_CoreInletDistribution_to_CoreInletBranch_2_3_Rep` | `V12_CoreInletDistribution` | `V12_CoreInletBranch_2_3_Rep_Vol_01` | macro branch, multiplier = 2 |
| `J_V12_CoreInletBranch_1_to_CoreInletConnector` | `V12_CoreInletBranch_1_Vol_08` | `V12_CoreInletConnector` | inherited V8/V12 branch outlet |
| `J_V12_CoreInletBranch_2_3_Rep_to_CoreInletConnector` | `V12_CoreInletBranch_2_3_Rep_Vol_08` | `V12_CoreInletConnector` | macro branch outlet, multiplier = 2 |
| `J_PlenumIn_*` | `V12_CoreInletConnector` | `Chan_*_Vol_01` | macro TFE inlet by representative multiplier |
| `J_PlenumOut_*` | `Chan_*_Vol_37` | `V12_CoreOutletConnector` | macro TFE outlet by representative multiplier |
| `J_CoreOutletConnector_to_Pipe05` | `V12_CoreOutletConnector` | `Pipe05_CoreOutletToRadiator_Vol_01` | hot leg start |
| `J_Pipe05_to_RadiatorInletSplit` | `Pipe05_CoreOutletToRadiator_Vol_08` | `V12_RadiatorInletSplit` | radiator inlet split node |

### 3.2 Radiator Header and Tube Junctions

| Junction Pattern | Count | From | To | K Loss |
| --- | ---: | --- | --- | ---: |
| `J_RadiatorInletSplit_to_UpperHeader_A` | 1 | `V12_RadiatorInletSplit` | `V12_RadiatorUpperHeader_01_Vol_01` | `connector_k_loss`, default 0 |
| `J_RadiatorInletSplit_to_UpperHeader_B` | 1 | `V12_RadiatorInletSplit` | `V12_RadiatorUpperHeader_40_Vol_01` | `connector_k_loss`, default 0 |
| `J_RadiatorUpperRing_XX_to_YY` | 78 | upper header segment `XX` | next upper header segment | `radiator_header_k_loss`, default 1 |
| `J_RadiatorUpper_to_Tube_XX` | 78 | upper header `XX` | `V12_RadiatorTubeFluid_XX_Vol_01` | `radiator_tube_inlet_k_loss`, default 100 |
| `J_RadiatorTube_XX_to_Lower` | 78 | `V12_RadiatorTubeFluid_XX_Vol_08` | lower header `XX` | `radiator_tube_outlet_k_loss`, default 100 |
| `J_RadiatorLowerRing_XX_to_YY` | 78 | lower header segment `XX` | next lower header segment | `radiator_header_k_loss`, default 1 |
| `J_LowerHeader_A_to_RadiatorOutletMix` | 1 | `V12_RadiatorLowerHeader_01_Vol_01` | `V12_RadiatorOutletMix` | `connector_k_loss`, default 0 |
| `J_LowerHeader_B_to_RadiatorOutletMix` | 1 | `V12_RadiatorLowerHeader_40_Vol_01` | `V12_RadiatorOutletMix` | `connector_k_loss`, default 0 |

`XX` runs from `01` to `78`; `YY` is the next segment with wrap-around.

### 3.3 Cold Return Junctions

| Junction | From | To | Notes |
| --- | --- | --- | --- |
| `J_V12_RadiatorOutletMix_to_Pipe06_RadiatorOutlet` | `V12_RadiatorOutletMix` | `Pipe06_RadiatorOutlet_Vol_01` | generated as `J_{previous.name}_to_{channel.name}` |
| `J_Pipe06_RadiatorOutlet_Vol_08_to_Pipe07_HeatExchangerHotSide` | `Pipe06_RadiatorOutlet_Vol_08` | `Pipe07_HeatExchangerHotSide_Vol_01` | hydraulic-only |
| `J_Pipe07_HeatExchangerHotSide_Vol_08_to_Pipe08_ReturnInnerPipe` | `Pipe07_HeatExchangerHotSide_Vol_08` | `Pipe08_ReturnInnerPipe_Vol_01` | hydraulic-only |
| `J_Pipe08_ReturnInnerPipe_Vol_08_to_Pipe09_ValveSegment` | `Pipe08_ReturnInnerPipe_Vol_08` | `Pipe09_ValveSegment_Vol_01` | hydraulic-only |
| `J_Pipe09_to_V13_PumpA` | `Pipe09_ValveSegment_Vol_08` | `V13_PumpMidNode` | pump A |
| `J_V13_PumpA_to_PumpB` | `V13_PumpMidNode` | `V13_PumpOutletNode` | pump B |

## 4. Solid Coupling and Boundary Conditions

### 4.1 Registered Solid/Coupler Summary

| Coupler Type | Count | Role |
| --- | ---: | --- |
| `FluidSolidCouple` | 88 | 10 core coolant/clad couplers + 78 radiator tube fluid/wall couplers |
| `GapCouple2D` | 17 | TFE radial gaps + global moderator/barrel/reflector radiation gaps |
| `TECCouple2D` | 5 | emitter-collector TEC gap for each representative TFE |
| `SolidSolidCouple2D` | 3 | global moderator ring-to-ring radial conduction |

The V13 builder creates radiator `FluidSolidCouple`s with `local_implicit` and TFE coolant `FluidSolidCouple`s with the class default `local_implicit`. The V13 closed-loop runner default is `--fluid-solid-coupling-scheme local_implicit`, and it applies that scheme to all `FluidSolidCouple`s with available solid capacitance.

### 4.2 TFE Radial Solid Stack

Every representative TFE has the same radial stack:

```text
pellet -> fission gas gap -> emitter -> TEC gap -> collector
  -> helium gap -> inner_clad -> coolant annulus -> outer_clad
  -> CO2 gap -> virtual moderator -> global moderator rings
```

TFE geometry:

| Radius / Region | Value (m) |
| --- | ---: |
| `r_pellet_inner` | 0.00400 |
| `r_pellet_outer` | 0.00850 |
| `r_fission_gas_outer` | 0.00865 |
| `r_emitter_outer` | 0.00980 |
| `r_collector_inner` | 0.01030 |
| `r_collector_outer` | 0.01185 |
| `r_inner_clad_inner` | 0.01190 |
| `r_inner_clad_outer` | 0.01225 |
| `r_coolant_inner` | 0.01225 |
| `r_coolant_outer` | 0.01295 |
| `r_outer_clad_outer` | 0.01330 |
| `r_moderator_inner` | 0.01352 |
| `r_moderator_outer` | 0.01627 |

TFE physical couplers and boundary conditions:

| Coupler Key / Object | Boundary | Condition Type | Parameters |
| --- | --- | --- | --- |
| `pellet_emitter_gap` | `pellet.right` <-> `emitter.left` | `GapCouple2D` | gap = 0.00015 m; `h_eq=5678 W/m2/K`; effective `k_gas=0.8517 W/m/K`; emissivity 0.15 / 0.15 |
| `tec_couple` | `emitter.right` <-> `collector.left` | `TECCouple2D` | gap = 0.00050 m; default `h_eq=29 W/m2/K`; effective `k_gas=0.0145 W/m/K`; emissivity 0.15 / 0.60; TEC electron heat and Joule heat sources are added when TEC is enabled |
| `collector_iclad_gap` | `collector.right` <-> `inner_clad.left` | `GapCouple2D` | gap = 0.00005 m; `h_eq=5678 W/m2/K`; effective `k_gas=0.2839 W/m/K`; emissivity 0.60 / 0.80 |
| `iclad_coolant_fsc` | `inner_clad.right` <-> `Chan_*` | `FluidSolidCouple` -> Robin resistance BC on solid + fluid source | heated perimeter = 0.0769690200 m; Nu correlation = annular `nu_ringpipe` |
| `oclad_coolant_fsc` | `outer_clad.left` <-> `Chan_*` | `FluidSolidCouple` -> Robin resistance BC on solid + fluid source | heated perimeter = 0.0813672497 m; Nu correlation = annular `nu_ringpipe` |
| `oclad_mod_gap` | `outer_clad.right` <-> `moderator.left` | `GapCouple2D` | gap = 0.00022 m; `h_eq=53.6 W/m2/K`; effective `k_gas=0.011792 W/m/K`; emissivity 0.80 / 0.80 |
| `mod_outer_bc` | `moderator.right` | `ResistanceBC` to `boundary_data.moderator_temperature` | `R_ext = 1 / (500 A)`; `ReactorCore.pre_step()` updates and maps virtual moderator outer flux to global moderator rings |

The axial `top` and `bottom` boundaries of TFE solids are not connected to external solid/fluid objects in V13; they retain solver boundary conditions and should be treated as axial end closure conditions, not as separate physical heat exchangers.

### 4.3 Global Moderator, Barrel, and Reflector

Global outer core structures are part of `TASTIN_Core_V8_CaseA`.

| Solid | Geometry | Boundary / Coupling |
| --- | --- | --- |
| `TASTIN_Core_V8_CaseA_ModRing_0...3` | 4 radial moderator rings, inner radius starts at TFE moderator outer radius 0.01627 m, outer radius 0.06000 m | adjacent rings connected by `SolidSolidCouple2D` |
| `TASTIN_Core_V8_CaseA_Barrel` | inner radius 0.06500 m; thickness 0.00300 m; radial nodes = 3 | connected to outer moderator by simplified gap |
| `TASTIN_Core_V8_CaseA_Reflector` | outer radius 0.10200 m; radial nodes = 8 | connected to barrel by simplified gap; outer surface radiates to environment |

Global gap/radiation boundary conditions:

| Interface | Condition Type | Parameters |
| --- | --- | --- |
| outer moderator -> barrel | `GapCouple2D` | gap = 0.00500 m; `h_eq=0`; radiation only; emissivity 0.80 / 0.80 |
| barrel -> reflector | `GapCouple2D` | gap = 0.00200 m; `h_eq=0`; radiation only; emissivity 0.80 / 0.80 |
| reflector outer surface | `DynamicRadiationResistanceBC` | `T_env=200 K`; emissivity = 0.20 |

### 4.4 Radiator Pipe-Fin Solid Boundaries

Each `V12_RadiatorTube_XX` component owns one solid wall:

```text
V12_RadiatorTube_XX_Wall
```

Boundary conditions:

| Wall Boundary | Coupled Object / Boundary | Condition Type | Parameters |
| --- | --- | --- | --- |
| `left` / inner wall | `V12_RadiatorTubeFluid_XX` | `FluidSolidCouple` -> Robin resistance BC + fluid source | heated perimeter = pi * 0.007 = 0.0219911486 m; node heat-transfer area = 0.0050854531 m2 for each axial node; Nu correlation = `nak_internal_nu`; runner default scheme = `local_implicit` |
| `right` / outer bare tube | radiation background | `DynamicRadiationResistanceBC` | tube emissivity = 0.80; default background = 3 K; bare tube area per axial node = 0.0058119464 m2 |
| `right` / fin root | quasi-steady fin equivalent branch | `ResistanceBC` updated in `RadiatorPipeWithFin.pre_step()` | fin radiation solved by reduced-order tridiagonal model; no separate fin ODE solid |
| `top`, `bottom` | no external physical component | solver boundary conditions | no modeled axial end heat exchanger |

The reduced-order fins radiate to the same `radiation_background_temperature` as the bare tube. If a `RadiatorThermalShield` is attached, it changes this equivalent background temperature through `RadiatorPipeWithFin.set_radiation_background_temperature(...)`.

### 4.5 Optional Radiation Shield

The base V13 closed-loop builder does not add a shield by itself. `attach_radiator_thermal_shield(...)` can attach:

```text
V13_RadiatorThermalShield
```

The shield is a boundary modifier, not a registered solid ODE. It reads radiator surface temperatures and sets each radiator unit's equivalent radiation background temperature. It can use either `segment_balance` or `fortran_shield2` depending on runner options.

## 5. Heat Transfer Locations

The only hydraulic objects in V13 that exchange heat with registered solids are:

| Hydraulic Object Pattern | Solid Object Pattern | Coupler |
| --- | --- | --- |
| `Chan_Center`, `Chan_Ring1`, `Chan_Ring2`, `Chan_Ring3_TEC`, `Chan_Ring3_Open` | representative TFE `inner_clad` and `outer_clad` solids | `*_iclad_coolant_couple`, `*_oclad_coolant_couple` |
| `V12_RadiatorTubeFluid_01...78` | `V12_RadiatorTube_01...78_Wall` | `V12_RadiatorTube_01...78_FluidSolid` |

The following hydraulic segments are adiabatic/hydraulic-only in current V13:

```text
Pipe11_CoreInletHeader
V12_CoreInletDistribution
V12_CoreInletBranch_1
V12_CoreInletBranch_2_3_Rep
Pipe05_CoreOutletToRadiator
V12_RadiatorInletSplit
V12_RadiatorUpperHeader_01...78
V12_RadiatorLowerHeader_01...78
V12_RadiatorOutletMix
Pipe06_RadiatorOutlet
Pipe07_HeatExchangerHotSide
Pipe08_ReturnInnerPipe
Pipe09_ValveSegment
V13_PumpMidNode
V13_PumpOutletNode
```

## 6. Flow Distribution Defaults

For the default design flow `total_inlet_flow_kg_s = 1.3 kg/s`:

| Flow Path | Default Value |
| --- | ---: |
| total loop flow | 1.3 kg/s |
| each radiator tube target design flow | 1.3 / 78 = 0.0166667 kg/s |
| each core representative physical TFE flow | 1.3 / 37 = 0.0351351 kg/s per real TFE |
| `V12_CoreInletBranch_1` macro branch | 1.3 / 3 = 0.433333 kg/s |
| `V12_CoreInletBranch_2_3_Rep` stored representative branch flow | 0.433333 kg/s, macro multiplier = 2 |

## 7. Important Modeling Caveats

1. Many V13 channels are hydraulic-only; the model does not know their wall thickness or outer diameter. The diagram therefore labels those dimensions as `D_h` and "OD not modeled" rather than inventing an outer diameter.
2. Radiator header rings are hydraulic-only; only the 78 vertical tube walls are solid heat-conduction objects.
3. `Pipe07_HeatExchangerHotSide` is a named flow-network segment only. V13 does not attach a heat exchanger solid or heat-transfer boundary to it.
4. Fin thermal behavior in `RadiatorPipeWithFin` is quasi-steady. Fin energy storage is not a separate state variable.
5. The startup runner can add external heat fluxes and a shield; those are boundary modifiers on the radiator path and are not part of the bare V13 closed-loop builder.
