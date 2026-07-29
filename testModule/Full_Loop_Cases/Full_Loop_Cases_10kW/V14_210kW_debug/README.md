# V14_10kW 210 kW Powered Debug Run

This folder contains the staged debug runner for the V14_10kW heat-pipe radiator full-loop case.

## Run Assumptions

- Core total power is fixed at `210000 W`.
- Total loop flow is fixed at `2.46 kg/s` through pump target-flow control.
- TEC electrical calculation is disabled by default with `main_tec_enabled=False`; use `--enable-tec` for the powered electrical tuning path.
- Point kinetics is not attached.
- All fluid and solid temperatures are initialized to `754.15 K` unless a restart is supplied.
- Solid heat conduction uses `implicit_euler`.
- Space background temperature is `4 K`.
- External orbital heat flux is disabled in the V14 heat-pipe radiator builder.

When a restart is loaded, the saved positive last-step `dt` is supplied to the
initial `local_implicit` fluid-solid coupling refresh before normal system
initialization.

Default run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_debug\run_v14_210kw_debug.py
```

Hydraulic/thermal staged run while debugging:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_debug\run_v14_210kw_debug.py --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_debug\runs\trial_001 --stage-durations 3600 --dt 0.05 --record-interval 10 --checkpoint-interval 20 --min-fluid-temperature-stop 500
```

Powered TEC tuning run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_debug\run_v14_210kw_debug.py --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_debug\runs\trial_tec --restart-in testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_debug\runs\previous\stage_01_restart.npz --stage-durations 1200 --dt 0.05 --record-interval 10 --checkpoint-interval 20 --min-fluid-temperature-stop 500 --enable-tec --tec-voltage 50.65 --enable-tec-lookup --wire-resistance-scale 0.335 --radiator-emissivity 0.7475 --hp-up-view-factor 0.0 --upper-hp-down-view-factor 0.3 --lower-hp-down-view-factor 0.4 --tec-current-guess 206
```

## Tuning Knobs

- `--radiator-emissivity`: applied to both heat-pipe condenser and fin equivalent radiation.
- `--upper-hp-down-view-factor` / `--lower-hp-down-view-factor`: lower-side view factors for upper/lower heat-pipe rings. The upper-side view factor is kept at `0.0` in the current tuning path.
- `--wire-resistance-scale`: scales the four preserved wire-resistance ratios `[0.001552, 0.001024, 0.000336, 0.000608] ohm`.
- `--enable-tec-lookup`: enables ThermoCalc lookup tables through `FullLoopCoreConfig`; omit it to force analytic ThermoCalc.

## Outputs

Each stage writes:

- `checkpoint_tXXXXXXs.npz` every checkpoint interval
- `stage_XX_restart.npz`
- `stage_XX_summary.json`
- `history.csv`
- `latest_state.json`
- `run_summary.json`

Tracked diagnostics include core inlet/outlet temperature, heat-pipe rejection, collector-ring wall rejection, total external heat rejection, estimated net power, required pump pressure rise, pump flow, and fluid/solid temperature extrema. The runner also stops and writes `emergency_low_temp_restart.npz` if `min_fluid_T_K` drops below `--min-fluid-temperature-stop` (default `500 K`; use `0` to disable).

## Current Best Powered Baseline

The closest long-window result in this folder is:

```text
runs/final_eps07475_u50p65_wire0335_1200s_from7964
```

Parameters:

```text
radiator_emissivity = 0.7475
tec_voltage = 50.65 V
wire_resistance_scale = 0.335
hp_up_view_factor = 0.0
upper_hp_down_view_factor = 0.3
lower_hp_down_view_factor = 0.4
```

Final metrics at `t = 9163.85 s`:

```text
core_inlet_T = 754.738 K
core_outlet_T = 845.773 K
delta_T = 91.034 K
TEC current = 206.569 A
TEC net electric power = 10.463 kW
total radiator rejection = 195.212 kW
heat-pipe rejection = 188.936 kW
ring-wall rejection = 6.275 kW
required pump head = 27.880 kPa
```

Relative to the present tuning targets (`754.45 K`, `845.65 K`, `206 A`, `10.44 kW`), this run is still slightly hot and slightly high in current/power. Treat it as the current working baseline, not a final converged calibration.
