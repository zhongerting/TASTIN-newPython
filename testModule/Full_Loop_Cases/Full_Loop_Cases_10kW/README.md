# Full_Loop_Cases_10kW

This directory is the 10 kW / five-ring V14 heat-pipe-radiator case package copied from `testModule/Full_Loop_Cases`. Changes are intended to stay inside this directory and the directly related 10 kW tests.

## Current Scope

The core model has been changed to five representative rings. The V14 heat-pipe radiator topology has also been rebuilt for the 10 kW package. The current powered debug path uses fixed total core power, fixed target loop flow, implicit-Euler solid conduction, optional TEC lookup calculation, and disabled orbital external heat flux.

## Five-Ring Core

Representative TFE groups:

```text
Center, Ring1, Ring2, Ring3, Ring4
1,      6,     9,     18,    24
```

Total physical TFE count is `58`. TEC multipliers follow the same physical counts. The main series TEC voltage target is changed from the old 37-TFE `27.2 V` setting to the 10 kW tuning range around `50.5 V`; the current powered debug baseline uses `50.65 V`.

Radial power shares:

```python
[0.019568969, 0.120310302, 0.180465534, 0.319655034, 0.360000072]
```

Axial shape coefficients:

```python
[5.0392372538e-02, -3.2174418071e-02, 6.847842042e-03, -4.513204066e-03, 2.890683804e-03]
```

## Moderator, Barrel, Reflector

TFE center radii from outer to inner:

```text
137e-3, 104e-3, 69.5e-3, 37e-3, 5e-3
```

Moderator ring boundaries derived by midpoint method, inner to outer:

```text
0, 21.0e-3, 53.25e-3, 86.75e-3, 120.5e-3, 164.0e-3
```

Barrel radii:

```text
inner = 164.0e-3
outer = 166.0e-3
```

Reflector radii:

```text
inner = 166.0e-3
outer = 261.0e-3
```

The moderator/barrel/reflector interfaces no longer model explicit gaps. They use `GlobalGapStructureConfig(mode="simplified", width=0.0, h_eq=5678.0)` as the equivalent interface resistance.

## TEC Lookup Control

`FullLoopCoreConfig` controls ThermoCalc lookup loading for this 10 kW case:

```python
FullLoopCoreConfig(
    tec_lookup_enabled=True,
    tec_lookup_db=r"ThermoCalc\emission_runtime_db_v2\pcs_0p02_5torr",
    tec_lookup_regions=("core", "startup", "high_power", "accident"),
)
```

`tec_lookup_enabled=None` keeps the older environment-variable behavior. Set it to `False` to force analytic ThermoCalc even when `THERMOCALC_ENABLE_LOOKUP=1` is present in the environment.

The V14 powered debug runner exposes `--enable-tec-lookup` and applies the lookup choice through the core config rather than a global-only switch.

## V14 Heat-Pipe Radiator

V14_10kW no longer uses the old single explicit ring plus `multiplier=2` symmetric ring simplification. It explicitly builds one upper heat-pipe collector ring and one lower heat-pipe collector ring:

- `Upper_A1...Upper_A6`
- `Lower_A1...Lower_A6`

The external hot branches and manifolds remain shared:

- `HotOutletBranch_1/2/3`
- `Manifold_1/2/3`

Each `HotOutletBranch_i` splits into two links entering `Upper_InletMix_Ii` and `Lower_InletMix_Ii`. At the outlet side, the corresponding upper and lower ring outlets merge into one shared `OutletMix_Oi`, then enter `Manifold_i`.

Heat-pipe lengths:

```text
L_eva = 0.0605 m
L_aba = 0.0415 m
L_con = 0.47 m
```

Fins are attached only to the condenser section `L_con`.

Heat-pipe counts:

```text
Upper: 5 * [8, 9, 9] + 1 * [8, 8, 8] = 154
Lower: 6 * [10, 10, 11] = 186
Total: 340
```

The initial/design hydraulic split is proportional to heat-pipe count:

```text
Upper fraction = 154 / 340
Lower fraction = 186 / 340
```

This is only an initial/design flow setting, not a fixed flow boundary. Later flow distribution is still governed by total pump head and resistance coefficients.

The heat-pipe radiation view-factor interface supports separate upper/lower ring lower-side factors:

```python
V14HeatPipeRadiatorConfig(
    hp_up_view_factor=0.0,
    upper_hp_down_view_factor=0.3,
    lower_hp_down_view_factor=0.4,
)
```

The current tuning path adjusts heat-pipe and fin emissivity together and keeps the upper-side view factor at zero.

The V14 thermal shield is now coupled to all 36 representative `HPwithFin` units. It is opt-in so existing V14 restarts and baselines are unchanged:

```python
V14HeatPipeRadiatorConfig(
    thermal_shield_enabled=True,
    thermal_shield_active_until_s=None,  # active until a withdrawal time is assigned
)
```

The `fortran_shield2` mapping uses upper-ring units 0-17 for shield sectors 0-5 and lower-ring units 18-35 for sectors 6-11. Each three-unit radiator sector is averaged in `T^4` using its physical heat-pipe multiplier; the 36 representatives therefore retain the full 154 + 186 = 340 heat-pipe radiation area. The shield component runs before every `RingHP` and writes the equivalent background temperature into the condenser tube and reduced-order fin radiation boundaries. `thermal_shield_active_until_s` is an absolute system time; once exceeded, the background returns to `t_space_k`.

When `thermal_shield_enabled` and `external_heat_enabled` are both true, V14 uses the six center samples `(0, 3, 6, 9, 12, 15)` from N18 as the numerically identical N6 history. Both histories use `5668.14 s` and one `time_origin_s`. With the shield present, only `qsss[:6] = 0.992 * N6` is active and direct N18 is zero; after `shield.set_active(False)`, `qsss` becomes zero and N18 is restored in the same pre-step. The latch and orbit origin are included in global restart files.

Use `run_v14_shield_radiator_startup.py` for Stage 0. The startup-only TEC gap contains helium with `h_eq=5678 W/m2/K`; the common powered-case default remains cesium with `h_eq=29 W/m2/K`. The accepted run starts the zero-power, TEC-off, shield-attached loop near `300 K` at `0.615 kg/s` (`25%` rated flow), advances `1800 s`, writes the five `history*.csv` tables every `10 s`, and keeps only the final restart. Its directory is `V14_210kW_start/phase_0_shielded_1800s_complete`; final coolant temperatures are `298.492-298.852 K`.

Stage 1 uses `V14_210kW_start/run_v14_reactivity_startup.py`. It starts point kinetics from `1 W`, holds `+0.50 $` without withdrawal, keeps TEC off and flow at `0.615 kg/s`, writes all five CSV tables every `1 s`, and stops at the first `10 kW` crossing. The accepted helium-gap run reached `10000.0009 W` in `141.45251 s` and saved only `V14_210kW_start/stage_1_fixed_0p50_to_10kw/final_restart.npz`.

Stage 2 uses `V14_210kW_start/run_v14_power_ramp.py`. It removes point kinetics, prescribes a `600 W/s` ramp from 10 to `70 kW`, keeps TEC off, writes the five CSV tables every `1 s`, and saves restart state every `50 s` plus the final state. The shield is jettisoned when the minimum loop coolant temperature reaches `373 K`; the controlled flow rises from `0.615` to `1.23 kg/s` when `CoreOutletConnector.T >= 500 K`. The accepted helium-gap ramp ended after `100.00000 s` at `70 kW`; neither threshold was reached (minimum coolant `304.401 K`, core outlet `358.896 K`, maximum coolant `379.636 K`), so the shield remained attached and flow remained at 25%.

When both `thermal_shield_enabled` and `external_heat_enabled` are true, V14 uses the six center samples `(0, 3, 6, 9, 12, 15)` from the N18 CSV as the numerically identical N6 shield history. N6 and N18 share the `5668.14 s` period and `time_origin_s`. While the shield is present, only `qsss[:6] = 0.992 * N6` is applied and all direct N18 sources are zeroed; after `shield.set_active(False)`, `qsss` is zero and N18 is restored in the same pre-step. The active latch and orbit origin are included in the global restart.

The formal Stage 0 entry point is:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.run_v14_shield_radiator_startup --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_start\phase_0_shielded_1800s_complete --duration 1800 --max-dt 0.2 --record-interval 10 --initial-temperature 300 --target-flow 0.615
```

Pass `--restart-in <previous/final_restart.npz>` for the next S-stage and add `--withdraw-shield` only at the stage that removes the shield.

## V14 210 kW Powered Debug Folder

`V14_210kW_debug/` contains the fixed-power, fixed-flow staged runner used for the 10 kW electrical tuning work. It records `history.csv`, checkpoint restarts, per-stage summaries, and the final `run_summary.json`.

Current best long-window baseline:

```text
V14_210kW_debug/runs/final_eps07475_u50p65_wire0335_1200s_from7964
```

Final state of that run:

```text
core inlet = 754.738 K
core outlet = 845.773 K
TEC current = 206.569 A
TEC net electric power = 10.463 kW
total radiator rejection = 195.212 kW
required pump head = 27.880 kPa
```

It is close to, but not exactly at, the present targets of `754.45 K`, `845.65 K`, `206 A`, and `10.44 kW`.

## V14 210 kW Reactivity-Control Baseline

`V14_210kW_reactivity_control/` provides the handoff from the fixed-power
steady restart to point kinetics. Its first baseline mode deliberately fixes
external reactivity at zero, requires the control drum to remain disabled, and
uses only temperature-feedback changes relative to the loaded 210 kW steady
state. A restart and its sibling `run_config.json` are mandatory.

The existing `V14_210kW_debug/` behavior remains fixed-power by default. The
new runner calls its builder with fixed-power application disabled, initializes
and calibrates point kinetics only for the first handoff, and preserves the
saved point-reactor and feedback-reference state on continuation.

## Contents

- Common layer: `common_config.py`, `common_core_builder.py`, `common_flow_builder.py`, `common_builder.py`, `common_diagnostics.py`
- V14 case: `v14_case.py`, `v14_heatpipe_radiator.py`
- V14 hydraulic smoke case: `V14_run_cases/`
- V14 powered debug/tuning case: `V14_210kW_debug/`
- V14 zero-input temperature-feedback baseline: `V14_210kW_reactivity_control/`
- Package entry: `__init__.py`

## Verification

Core geometry test:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_full_loop_10kw_core_geometry
```

V14 hydraulic-only smoke:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_run_cases.test_v14_flow_path_smoke
```

For `V14_210kW_debug`, prefer direct runner invocations instead of `python -m unittest`; the debug workflow is a staged case runner, not a unit-test entry.

## V14 Heat-Pipe Transfer-Failure Accidents

The following three accident folders start from
`V14_210kW_fixed_power_external_heat_2orbits/runs/two_orbits_from13864_20260720/checkpoint_t019865s.npz`:

| Folder | Failure signature |
| --- | --- |
| `V14_210kW_heatpipe_partial_failure/` | Upper A5, local node 2, effective heat transfer 50% |
| `V14_210kW_heatpipe_single_node_failure/` | Matching upper/lower A5 local node 2, effective heat transfer 0% |
| `V14_210kW_heatpipe_sector_failure/` | Matching upper/lower A5, all three local nodes, effective heat transfer 0% |

The failure changes only fluid-to-evaporator heat-transfer coupling. Nominal heat-pipe counts, hydraulic loss maps, flow areas, radiation areas, and orbital external heat remain unchanged. TEC calculation remains enabled and the thermal shield remains disabled. Coolant temperature is recorded but is not a trip criterion. The initial controller holds 210 kW; the first solid-temperature limit crossing switches to point kinetics and applies -2 dollars. Limits are channel wall 1058 K, fuel pellet 2700 K, collector 1023 K, moderator 930 K, and reflector 1000 K.

Each runner defaults to one orbital period (`5668.144369 s`), records the five history tables every 1 s, writes periodic restart files every 100 s, and writes explicit accident-start, scram-event, and final restart files. If scram occurs, the run continues until at least half an orbital period after scram (and never shorter than the configured minimum duration).
## Half-Radiator Transfer-Failure Accident

`V14_210kW_heatpipe_half_radiator_failure/` disables sectors A1-A3 in both
upper and lower rings. This disables 6 of 12 radiator sectors and 171 of the
nominal 340 heat pipes (about 50.3 percent), while preserving the hydraulic
paths and orbital external-heat boundary.
