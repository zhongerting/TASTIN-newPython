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

## Contents

- Common layer: `common_config.py`, `common_core_builder.py`, `common_flow_builder.py`, `common_builder.py`, `common_diagnostics.py`
- V14 case: `v14_case.py`, `v14_heatpipe_radiator.py`
- V14 hydraulic smoke case: `V14_run_cases/`
- V14 powered debug/tuning case: `V14_210kW_debug/`
- Package entry: `__init__.py`

## Verification

Core geometry test:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_full_loop_10kw_core_geometry
```

V14 hydraulic-only smoke:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_run_cases.test_v14_flow_path_smoke
```

For `V14_210kW_debug`, prefer direct runner invocations instead of `python -m unittest`; the debug workflow is a staged case runner, not a unit-test entry.
