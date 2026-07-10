# V13 Startup TEC Ignition Status - 2026-06-28

## Current startup sequence

The V13 startup runner now separates cesium conditioning from electrical TEC startup:

1. Helium-gap thermal startup with TEC electrical calculation disabled.
2. Cesium conditioning with TEC electrical calculation still disabled; `startup_cs_fraction` ramps toward 1 and the gap equivalent heat-transfer coefficient approaches the V13 cesium value.
3. Electrical startup after cesium and emitter-temperature gates are satisfied. Main TEC starts in fixed resistance mode and switches to fixed total voltage only if the finite main voltage reaches the configured switch voltage.

Official V13 cold-start continuation in this note uses `startup_thermal_power_w = 110000 W`. Any `power150k` directory or `150 kW` paragraph below is a solver/lookup diagnostic only; it must not be cited as the V13 physical startup trajectory.

## Code fixes made for this workflow

- `testModule/v13_startup_control.py`
  - Split cesium conditioning and electrical TEC latches.
  - Added `seed_cesium_conditioning()` for restart continuation from an already cesium-conditioned state.
- `testModule/run_v13_start_case.py`
  - Added CLI gates for electrical startup after cesium conditioning.
  - Added `--initial-cs-fraction` for restart continuation.
  - Added `--startup-main-tec-load-resistance-scope total|per_tec`; `per_tec` multiplies the input resistance by the main series TEC count.
  - Records `tec_solver_mode`, `tec_solver_converged`, `tec_solver_iteration_count`, `tec_solver_zero_emission_skipped`, `tec_solver_zero_emission_reason`, and `tec_solver_output_finite` in startup history/latest state.
  - Configures startup wire resistance without forcing an unsynchronized first TEC calculation.
- `Components/ReactorCore.py`
  - `pre_step()` now calls `_sync_tec_group_temperatures(group)` before `ThermoCalc.calculate()`.

## Verification

Commands that passed:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile Components\ReactorCore.py testModule\run_v13_start_case.py testModule\test_v13_startup_control.py testModule\test_reactorcore_control_drum.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_v13_startup_control.py
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_reactorcore_control_drum.py
```

## Stable cesium-conditioned thermal state

Stable pre-electrical/cesium-conditioned restart:

```text
testModule/v13_start_cesium_conditioning_plus1000s_20260628/v13_start_cesium_conditioning_plus1000s_20260628_latest_restart.npz
```

Final state at `t = 4090 s`:

- `startup_cs_fraction = 0.9999983895`
- `startup_tec_gap_h_eq_w_m2_k = 250.0006`
- `mean_emitter_temperature_k = 1198.96 K`
- `core_inlet_connector_t_k = 746.56 K`
- `core_outlet_connector_t_k = 842.32 K`
- TEC electrical calculation disabled in that run.

## Fixed-R electrical startup tests

Lookup environment used for stable tests:

```powershell
$env:THERMOCALC_ENABLE_LOOKUP='1'
$env:THERMOCALC_LOOKUP_DB='E:\????\??-??\source_code\TASTIN-python\ThermoCalc\emission_runtime_db_v2\pcs_0p02_5torr'
$env:THERMOCALC_LOOKUP_REGIONS='core,startup,high_power,accident'
```

### Total-load interpretation of 0.0044 ohm

`0.0044 ohm` used directly as total series load was numerically fragile and produced non-finite TEC heat sources in a 100 s lookup run.

Interpreting `0.0044 ohm` as per-TEC load for a 34-TEC main series chain gives. This can now be selected directly with `--startup-main-tec-load-resistance-scope per_tec`:

```text
R_total = 34 * 0.0044 = 0.1496 ohm
```

This path was stable.

Stable fixed-R transition output:

```text
testModule/v13_start_tec_fixedr_totalR_lookup_500s_20260628/v13_start_tec_fixedr_totalR_lookup_500s_20260628_latest_restart.npz
```

Final state at `t = 4612 s`:

- `Rload = 0.1496 ohm`
- `tec_main_voltage_v = 9.70466 V`
- `tec_main_current_a = 64.8707 A`
- `tec_main_electric_power_w = 629.55 W`
- `mean_emitter_temperature_k = 1188.26 K`
- `core_inlet_connector_t_k = 745.85 K`
- `core_outlet_connector_t_k = 841.21 K`
- Fixed-voltage switch was not triggered.

## Load scan at the 4612 s thermal state

Fixed-R load scan with lookup enabled:

| Rload ohm | Voltage V | Current A | Power W |
| ---: | ---: | ---: | ---: |
| 0.05 | 6.5945 | 131.891 | 869.76 |
| 0.10 | 8.7083 | 87.083 | 758.34 |
| 0.121 | 9.2422 | 76.382 | 705.94 |
| 0.1496 | 9.8413 | 65.784 | 647.40 |
| 0.2 | 10.6402 | 53.201 | 566.07 |
| 0.3 | 11.5695 | 38.565 | 446.18 |
| 0.5 | 12.6449 | 25.290 | 319.79 |
| 0.8 | 13.6678 | 17.085 | 233.51 |
| 1.2 | 14.3326 | 11.944 | 171.19 |
| 1.5 | 14.6237 | 9.749 | 142.57 |
| 2.0 | 14.9269 | 7.463 | 111.41 |
| 3.0 | 15.2429 | 5.081 | 77.45 |
| 5.0 | 15.9583 | 3.192 | 50.93 |
| 8.0 | 16.9090 | 2.114 | 35.74 |
| 12.0 | 17.5564 | 1.463 | 25.69 |
| 20.0 | 18.1388 | 0.907 | 16.45 |
| 50.0 | 18.7914 | 0.376 | 7.06 |
| 100.0 | 19.9564 | 0.200 | 3.98 |

## Current gate

At the current V13 startup thermal state and ThermoCalc backend, the fixed-R load line does not reach `27.2 V`; even high load/open-circuit direction remains around `20 V`. Direct fixed-U `27.2 V` is not yet trustworthy: lookup-enabled fixed-U returned `Uout=nan` with `iteration_count=100` in a diagnostic.

The next safe continuation gate is one of:

1. Update/fix ThermoCalc so fixed-U `27.2 V` returns finite, converged results; then continue from the 500 s fixed-R restart above.
2. Revisit the physical target: power level, cesium pressure/table, TEC load definition, or expected 27.2 V condition, because the current V13 startup thermal state cannot reach that voltage by fixed-R ramp alone.

## 2026-06-28 late update: rebuilt ThermoCalc C++ iteration guard

The local production backend was updated after the earlier diagnostics:

```text
ThermoCalc/te_solver.cp312-win_amd64.pyd
LastWriteTime = 2026-06-28 23:42:33
Length = 396288 bytes
```

The new backend no longer hangs on the tested startup fixed-voltage case. A clean `0.5 s` V13-start smoke from the stable fixed-R restart completed with fixed total voltage `27.2 V`:

```text
output = testModule/v13_start_tec_fixedu27p2_smoke_0p5s_20260628
restart_in = testModule/v13_start_tec_fixedr_totalR_lookup_500s_20260628/v13_start_tec_fixedr_totalR_lookup_500s_20260628_latest_restart.npz
t = 4612.5 s
tec_solver_mode = fixed_u
tec_solver_output_finite = True
tec_solver_converged = False
tec_solver_iteration_count = 45
tec_main_voltage_v = 27.2 V
tec_main_current_a = 0.0 A
tec_main_electric_power_w = 0.0 W
mean_emitter_temperature_k = 1224.24 K
```

Interpretation: the C++ dead-loop/NaN behavior is improved, but this is still not a usable `27.2 V` power-generation point. The backend returns a finite non-converged zero-current state instead of hanging.

A no-time-advance ThermoCalc temperature sensitivity was then run from the same restart. The table below adds a uniform emitter-temperature offset while keeping collector temperatures unchanged. This is a diagnostic gate estimate, not a thermal transient result.

| Emitter offset [K] | Rload [ohm] | Uout [V] | Iout [A] | Pout [W] | Converged |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0.1496 | 9.841 | 65.784 | 647.4 | True |
| 100 | 0.1496 | 14.496 | 96.896 | 1404.6 | True |
| 200 | 0.1496 | 18.957 | 126.720 | 2402.3 | True |
| 300 | 0.1496 | 23.300 | 155.746 | 3628.8 | True |
| 400 | 0.1496 | 27.665 | 184.923 | 5115.8 | True |
| 100 | 100.0 | 27.037 | 0.270 | 7.31 | True |
| 250 | 1.0 | 28.780 | 28.780 | 828.3 | True |

Current implication: with the selected startup fixed resistance `R_total = 0.1496 ohm`, the automatic fixed-R-to-fixed-U switch will not occur until the emitter field is roughly `+400 K` hotter than the current fixed-R restart state. In mean-temperature terms this is about `1.58e3 K` for this synthetic offset. Higher startup resistance can reach the voltage threshold at lower emitter temperature, but that changes the electrical startup path and should be treated as a modeling choice rather than a numerical fix.

## 2026-06-29 R=100 ohm fixed-resistance continuation check

After confirming the ThermoCalc no-hang backend, a higher startup resistance was tested as an alternative fixed-R ignition strategy. The case continued from the stable `R_total=0.1496 ohm` restart but used `R_total=100 ohm` and did not force fixed-voltage mode.

```text
output = testModule/v13_start_tec_fixedr_R100_20to120s_20260629
restart_in = testModule/v13_start_tec_fixedr_R100_20s_20260629/v13_start_tec_fixedr_R100_20s_20260629_latest_restart.npz
time span = 4632 s -> 4732 s
```

Final state:

```text
t = 4732 s
tec_solver_mode = fixed_r
tec_solver_converged = True
tec_solver_iteration_count = 2
Rload = 100 ohm
tec_main_voltage_v = 23.18796 V
tec_main_current_a = 0.23188 A
tec_main_electric_power_w = 5.3768 W
mean_emitter_temperature_k = 1196.56 K
core inlet/outlet ~= 745.41 / 840.72 K
fixed-voltage switch not triggered
```

Trend over this 100 s continuation: voltage increased from about `22.14 V` to `23.19 V`, and mean emitter temperature increased from about `1193.93 K` to `1196.56 K`. The path is numerically stable and no longer blocked by the ThermoCalc dead-loop issue, but it still remains below the `27.2 V` switching threshold.

## 2026-06-29 R=100 ohm 500 s continuation plateau

The `R_total=100 ohm` startup path was continued for another `500 s` from the previous R=100 ohm restart while keeping the automatic `27.2 V` fixed-U switch gate enabled.

```text
output = testModule/v13_start_tec_fixedr_R100_120to620s_20260629
restart_in = testModule/v13_start_tec_fixedr_R100_20to120s_20260629/v13_start_tec_fixedr_R100_20to120s_20260629_latest_restart.npz
time span = 4732 s -> 5232 s
```

Final state:

```text
t = 5232 s
tec_solver_mode = fixed_r
tec_solver_converged = True
tec_solver_iteration_count = 2
Rload = 100 ohm
tec_main_voltage_v = 23.32207 V
tec_main_current_a = 0.23322 A
tec_main_electric_power_w = 5.4392 W
mean_emitter_temperature_k = 1199.03 K
core inlet/outlet ~= 746.65 / 842.47 K
q_radiator_total_w ~= 108.58 kW
coolant_enthalpy_rise_w ~= 108.59 kW
fixed-voltage switch not triggered
```

Voltage trend over the segment:

| relative time [s] | U [V] | mean emitter [K] |
| ---: | ---: | ---: |
| 50 | 23.31444 | 1197.44 |
| 100 | 23.32500 | 1197.92 |
| 200 | 23.32683 | 1198.47 |
| 300 | 23.32503 | 1198.76 |
| 400 | 23.32336 | 1198.92 |
| 500 | 23.32207 | 1199.03 |

The voltage is effectively plateaued around `23.32 V`, well below the `27.2 V` switching threshold. The last-100-s voltage slope is slightly negative (`-1.28e-5 V/s`), even though mean emitter temperature is still creeping upward. This indicates that simply continuing the `R=100 ohm` fixed-resistance path is unlikely to reach the rated-voltage gate under the current 110 kW startup thermal state.

Current gate after this run: the ThermoCalc dead-loop issue is no longer the limiting factor. The limiting factor is the selected thermal/electrical startup state. To proceed to rated fixed-voltage generation, the model needs either a hotter emitter field, a revised startup load strategy, or revised assumptions for the cesium-gap/electrical target condition.

## 2026-06-29 150 kW sensitivity: fixed-R can cross voltage gate, fixed-U still stalls

A controlled higher-power sensitivity was run from the `R=100 ohm` plateau restart to determine whether the rated-voltage gate is reachable by increasing the emitter thermal state. This is a diagnostic sensitivity, not a calibrated startup baseline.

First attempt used the normal `27.2 V` automatic switch gate:

```text
output = testModule/v13_start_tec_fixedr_R100_power150k_200s_20260629
restart_in = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
fixed-startup-power-w = 150000
Rload = 100 ohm
```

The fixed-R records before switching were:

| relative time [s] | U [V] | I [A] | P [W] | mean emitter [K] |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 26.65647 | 0.26656 | 7.1057 | 1230.20 |
| 40 | 28.30129 | 0.28301 | 8.0096 | 1255.42 |

The `40 s` record exceeded the `27.2 V` switch threshold. After that point the process continued consuming CPU without writing further records and was stopped manually. To isolate the cause, the same 150 kW case was repeated with the switch threshold set to `999 V`, producing a high-temperature fixed-R restart at `t=5272 s`. From that restart:

- A `0.5 s` fixed-R step completed in about `28 s` and returned finite converged results:

```text
output = testModule/v13_start_tec_fixedr_R100_power150k_from5272_0p5s_20260629
t = 5272.5 s
mode = fixed_r
U = 28.70997 V
I = 0.28710 A
P = 8.2426 W
mean emitter = 1255.92 K
converged = True
iteration_count = 2
```

- A `0.5 s` direct fixed-U `27.2 V` step from the same restart did not return within the 600 s command timeout and was stopped manually:

```text
output = testModule/v13_start_tec_fixedu27p2_power150k_from5272_0p5s_20260629
entry state: t = 5272.0 s, mean emitter ~= 1255.91 K
mode = fixed_u
U target = 27.2 V
I guess = 0.28 A
```

Interpretation: the thermal/electrical fixed-R gate can be crossed at higher power, so the startup strategy is meaningful in principle. The current blocker is the fixed-U circuit solve at the switch point: even after the no-hang C++ guard update, this high-temperature fixed-U case can still run for more than 10 minutes without completing. Do not launch a long automatic fixed-R-to-fixed-U run until fixed-U convergence at the switch-point restart is repaired or a different voltage-control transition method is implemented.

## 2026-06-29 root-cause note: fixed-U stall is dominated by lookup coverage gap at hot collector nodes

Systematic debugging narrowed the high-temperature fixed-U stall to the ThermoCalc layer. At the `150 kW / R=100 ohm` switch-point restart (`t=5272 s`), a fixed-R ThermoCalc solve completes, but it is already expensive:

```text
fixed-R calculate time ~= 23608 ms
mode = fixed_r
Uout = 28.70997 V
Iout = 0.28710 A
converged = True
iteration_count = 2
```

The solved TEC state has 1258 axial emission points. Lookup coverage at the fixed-R point was:

```text
lookup_found = 1221 / 1258
missed points = 37
missed TE range = 1769.97 - 1862.02 K
missed TC range = 1100.49 - 1133.69 K
missed Tcs = 600 K
missed Vo range = 0.88392 - 0.90588 V
```

The current dense runtime v2 region axes explain the miss:

```text
startup:   TE 700-1280.6 K, TC 500-800 K
core:      TE 1300-2150 K,  TC 700-900 K
high_power:TE 2160-2400 K,  TC 750-1000 K
accident:  TE 700-2400 K,   TC 500-1100 K
```

The missed points sit just above the current `accident` collector-temperature ceiling (`TC > 1100 K`) while remaining below the `high_power` emitter-temperature floor (`TE < 2160 K`). They therefore fall back to the original analytic emission solve. A single fixed-R top-level solve is still tolerable because it converges in two outer iterations, but fixed-U repeatedly calls the same full series `circuitCalc(I)` path and can spend many minutes in fallback-heavy trial states.

Current root-cause conclusion: the remaining fixed-U switch-point stall is not the old low-temperature zero-emission dead loop. It is a high-temperature lookup-coverage/performance problem combined with the expensive series fixed-voltage secant iteration. A practical next fix is to extend or add a runtime lookup region covering at least:

```text
TE ~= 1700-1900 K
TC ~= 1100-1150 K
Vo ~= 0.8-1.0 V
Tcs ~= 600 K
```

or more conservatively extend `accident` to `TC >= 1150 K` for the relevant `TE/Vo/Tcs` range. After regenerating/loading that table, repeat the `fixed_u 27.2 V` smoke from `testModule/v13_start_tec_fixedr_R100_power150k_60s_noswitch_20260629_latest_restart.npz` before launching any automatic fixed-R-to-fixed-U long run.

### 2026-06-29 augmented lookup smoke after circuit no-hang update

After confirming the production `ThermoCalc/te_solver.cp312-win_amd64.pyd` timestamp (`2026-06-28 23:49:48`, `396288` bytes), a narrow local dense runtime lookup band was attached to the existing `pcs_0p02_5torr` runtime database through an augmented manifest:

```text
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr_augmented/runtime_dense_manifest.json
regions = switch_point, accident, core, high_power, startup
```

This manifest references the existing base `.tedb` files and the small generated switch-point `.npz`; it does not copy the full runtime database.

At the `150 kW / R=100 ohm` high-temperature restart (`t=5272 s`), a direct fixed-R ThermoCalc diagnostic with the augmented lookup produced:

```text
loaded dense regions = 5
mode = fixed_r
Uout ~= 28.68281 V
Iout ~= 0.28683 A
converged = True
iteration_count = 2
ThermoCalc calculate wall ~= 0.61 s
lookup_found = 1258 / 1258
```

So the previous lookup miss band is covered; the fixed-R switch-point solve is no longer dominated by analytic fallback.

A fixed-U `27.2 V` V13-start `0.5 s` smoke from the same restart completed without hanging:

```text
output = testModule/v13_start_tec_fixedu27p2_power150k_from5272_0p5s_auglookup_20260629
wall time ~= 430.6 s
absolute time = 5272.5 s
mean emitter ~= 1255.92 K
core inlet/outlet ~= 751.54 / 834.47 K
tec_solver_mode = fixed_u
tec_solver_output_finite = True
tec_solver_converged = False
tec_solver_iteration_count = 100
reported Uout ~= 35.105 V
reported Iout ~= 172.932 A
reported electric power ~= 6070.83 W
```

Current conclusion: the updated backend now prevents the old fixed-U dead loop and returns finite diagnostics, but the fixed-U branch is still not a usable production continuation at this switch point. It reaches the iteration cap and does not satisfy the target-voltage solution. The immediate blocker has moved from lookup coverage/dead-loop behavior to the fixed-voltage circuit solve strategy or the physical load-control transition.

### 2026-06-29 fixed-U bracketing repair and hot lookup continuation

The high-temperature fixed-U problem was separated into two parts.

First, the series fixed-voltage solver itself was repaired in source. The old `uFixedCircuitCalc()` used an unbounded secant update from `Iout` and `Iout+10 A` and did not check voltage residual convergence. The updated source brackets `Utarget - circuitCalc(I)` over non-negative current samples and then uses a guarded secant/bisection solve. Validation was done with the rebuilt test pyd in `ThermoCalc/build_cp312/Release` through `THERMOCALC_PYD_DIR`; the root production `.pyd` was not overwritten during this verification.

Second, the local startup lookup table was extended again:

```text
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr_switch_point_hot/
```

The new `switch_point_hot` region covers `TE=1700-2300 K`, `TC=1080-1300 K`, `Vo=0.70-1.20 V`, and `Tcs=580-610 K`. The augmented manifest now loads six regions: `switch_point_hot`, `switch_point`, `accident`, `core`, `high_power`, and `startup`.

Validated sequence with the rebuilt pyd and augmented lookup:

```text
fixed-U isolated at t=5272 s:
  U = 27.2 V, I ~= 1.2148 A, converged=True, iteration_count=10, wall ~= 0.736 s

fixed-R -> fixed-U automatic 2 s smoke:
  switch at t=5272.5 s when fixed-R U ~= 28.7095 V
  fixed-U records at 1.0/1.5/2.0 s all converged

fixed-U 20 s continuation:
  final t = 5294 s
  U = 27.2 V, I ~= 2.499 A, P ~= 67.98 W
  mean emitter ~= 1274.66 K

fixed-U isolated at t=5314 s before hot lookup:
  converged=True but wall ~= 28.7 s
  lookup_found = 1254 / 1258
  miss band TE ~= 1880-1925 K, TC ~= 1164-1173 K, Vo ~= 0.918 V

fixed-U isolated at t=5314 s after hot lookup:
  converged=True, iteration_count=11, wall ~= 4.25 s
  lookup_found = 1258 / 1258

fixed-U 40 s continuation from t=5314 s:
  final t = 5354 s
  U = 27.2 V, I ~= 4.026 A, P ~= 109.50 W
  mean emitter ~= 1310.85 K
  core inlet/outlet ~= 770.36 / 883.56 K
  all records converged
```

The startup fixed-R to fixed-U chain is now numerically viable for the tested segment, but this is not yet a steady-state result. Continuing toward steady will require monitoring lookup coverage as collector temperature rises; otherwise the run will fall back to the analytic emission solver and slow down again.

### 2026-06-29 correction: 150 kW runs are diagnostic only

The `power150k` V13 startup directories are not official V13 cold-start physics results. They were created only to force the fixed-R voltage gate above `27.2 V` so the fixed-voltage ThermoCalc branch and lookup coverage could be debugged.

Official V13 startup interpretation remains based on the startup-control `110 kW` hold unless the model assumptions are explicitly changed. The latest official `110 kW`, `R_total=100 ohm` fixed-resistance continuation is:

```text
output = testModule/v13_start_tec_fixedr_R100_120to620s_20260629
absolute_time = 5232 s
core/startup thermal power = 110000 W
mode = fixed_r
R_total = 100 ohm
U ~= 23.3221 V
I ~= 0.233221 A
P_e ~= 5.439 W
mean emitter ~= 1199.03 K
core inlet/outlet ~= 746.65 / 842.47 K
fixed-voltage switch = not triggered
last-100-s voltage slope ~= -1.28e-5 V/s
```

Therefore, under the current `110 kW` thermal state and `R_total=100 ohm` startup load, the formal fixed-R path has not reached the rated `27.2 V` switch condition. Any `150 kW` fixed-U continuation should be cited only as a numerical solver/lookup diagnostic, not as the V13 cold-start operating trajectory.

### 2026-06-29 official 110 kW load-ceiling probe

A no-time-advance electrical probe was run from the official `110 kW` restart:

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
absolute_time = 5232 s
startup power = 110000 W
mean emitter from rebuilt system ~= 1199.43 K
TEC gap h_eq ~= 250 W/m2/K
```

The root production pyd probe with a broad scan did not return and had to be stopped, so the two decisive points were repeated with the verified rebuilt test pyd via `THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release` and the augmented lookup regions. This did not overwrite the production pyd.

```text
R_total = 100 ohm:   U ~= 23.2822 V, I ~= 0.232822 A, P ~= 5.421 W, converged=True, iteration_count=2
R_total = 1e6 ohm:  U ~= 23.3653 V, I ~= 2.34e-5 A, P ~= 5.46e-4 W, converged=True, iteration_count=2
```

Interpretation: at the current official `110 kW` thermal state, even the near-open-circuit load voltage is only about `23.37 V`, below the `27.2 V` fixed-voltage switch threshold. Continuing the same fixed-R strategy cannot produce a legitimate switch to rated fixed-voltage generation unless the emitter/collector thermal state or the startup electrical/physical assumptions are changed. The earlier `150 kW` runs remain solver diagnostics only.

### 2026-06-29 official 110 kW R=100 continuation timeout check

A formal continuation was attempted from the official `5232 s` restart without changing the startup power:

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
output = testModule/v13_start_tec_fixedr_R100_620to720s_20260629
requested duration = 100 s
startup power = 110000 W
R_total = 100 ohm
switch voltage = 27.2 V
```

The process timed out after 30 min of wall time and was stopped manually. It had written history records through `t=5312 s` but the latest restart/state had only reached `t=5282 s` because the next restart interval had not completed.

Written history trend:

```text
t=5252 s: U ~= 23.120 V, mean emitter ~= 1198.55 K
t=5272 s: U ~= 23.103 V, mean emitter ~= 1198.20 K
t=5292 s: U ~= 22.981 V, mean emitter ~= 1198.46 K
t=5312 s: U ~= 22.941 V, mean emitter ~= 1198.78 K
```

The run remained in `fixed_r` and did not trigger the `27.2 V` fixed-voltage gate. This reinforces the official-path conclusion: under the current `110 kW`, `R_total=100 ohm` startup state, continuing the same fixed-resistance strategy is moving away from the voltage threshold, not toward it. The timeout also shows that the current production pyd/runtime path is not suitable for further long official startup continuation without first resolving the ThermoCalc performance/hang behavior or explicitly switching to a verified rebuilt backend.

### 2026-06-29 official 110 kW emitter-temperature margin probe

A no-time-advance ThermoCalc probe was run from the official `110 kW` restart to estimate how much hotter the emitter field must be before the startup voltage gate is physically reachable. The probe only changed the ThermoCalc input temperature matrix and did not modify the solid/restart state.

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
absolute_time = 5232 s
startup power = 110000 W
mean TEC emitter input = 1192.77 K
mean TEC collector input = 908.87 K
near-open load = R_total = 1e6 ohm
backend = verified rebuilt pyd via THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release
```

Near-open-circuit voltage versus uniform emitter offset:

```text
dTE =   0 K: U ~= 23.370 V
dTE =  25 K: U ~= 25.235 V
dTE =  50 K: U ~= 27.166 V
dTE =  75 K: U ~= 29.116 V
dTE = 100 K: U ~= 31.017 V
```

Linear interpolation between `25 K` and `50 K` gives a `27.2 V` crossing at approximately `dTE ~= 50.4 K`. Thus the official `110 kW` state is only about `50 K` short in emitter-side voltage margin for a near-open-circuit load. The formal transient, however, was trending slightly cooler/lower-voltage under `R_total=100 ohm`, so reaching the gate by simply continuing the current trajectory is still unlikely.

Interpretation: the limiting issue is primarily the startup thermal state seen by the TECs, not the fixed-voltage switch implementation. Plausible next modeling actions are to heat/hold the emitter field by revisiting startup power schedule, cesium-gap conductance trajectory, radiator/shield heat rejection during ignition, or the electrical load schedule; forcing fixed-U at the current unmodified `110 kW` restart remains unjustified.

### 2026-06-29 cesium-gap h_eq sensitivity for official 110 kW startup

To test whether the `27.2 V` gate is reachable without raising reactor power above the official `110 kW` hold, two thermal-only continuations were run from the official `5232 s` restart. TEC electrical calculation was deliberately disabled by delaying the electrical-start gate, but the cesium gap conductance was changed and applied to the thermal model.

Common setup:

```text
restart = testModule/v13_start_tec_fixedr_R100_120to620s_20260629/v13_start_tec_fixedr_R100_120to620s_20260629_latest_restart.npz
startup power = 110000 W
TEC electrical calculation = disabled during thermal hold
initial Cs fraction ~= 1.0
```

`h_eq = 150 W/m2/K` for `200 s`:

```text
output = testModule/v13_start_h150_tec_off_200s_20260629
t=5432 s
console mean emitter ~= 1295.49 K
TEC input mean emitter ~= 1288.62 K
TEC input mean collector ~= 900.24 K
R_total=100 ohm probe: U ~= 37.514 V, I ~= 0.3751 A, P ~= 14.07 W
R_total=1e6 ohm probe: U ~= 38.264 V
```

`h_eq = 200 W/m2/K` for `200 s`:

```text
output = testModule/v13_start_h200_tec_off_200s_20260629
t=5432 s
console mean emitter ~= 1241.66 K
TEC input mean emitter ~= 1235.09 K
TEC input mean collector ~= 906.57 K
R_total=100 ohm probe: U ~= 29.542 V, I ~= 0.2954 A, P ~= 8.73 W
R_total=1e6 ohm probe: U ~= 29.881 V
```

Interpretation: lowering the Cs-filled gap conductance from `250` to `200 W/m2/K` during the post-cesium thermal hold is already enough to exceed the rated-voltage gate at official `110 kW`; `150 W/m2/K` overshoots substantially. This points to the Cs-gap thermal trajectory as a plausible startup-control lever, unlike the earlier diagnostic `150 kW` power increase.

### 2026-06-29 h_eq=200 automatic fixed-R to fixed-U smoke

A short automatic electrical-start smoke was then run from the `h_eq=200` thermal hold restart:

```text
restart = testModule/v13_start_h200_tec_off_200s_20260629/v13_start_h200_tec_off_200s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedr_to_fixedu_2s_20260629
startup power = 110000 W
initial mode = fixed_r
R_total = 100 ohm
switch voltage = 27.2 V
target fixed voltage = 27.2 V
backend = verified rebuilt pyd via THERMOCALC_PYD_DIR=ThermoCalc/build_cp312/Release
```

The automatic switch triggered at `t=5432.5 s`:

```text
t=5432.5 s: fixed_r, U ~= 29.767 V, I ~= 0.2977 A, P ~= 8.86 W, converged=True, iter=2
t=5433.0 s: fixed_u, U = 27.2 V, I ~= 1.2977 A, P ~= 35.30 W, converged=True, iter=29
t=5433.5 s: fixed_u, U = 27.2 V, I ~= 1.2978 A, P ~= 35.30 W, converged=True, iter=42
t=5434.0 s: fixed_u, U ~= 27.178 V, I ~= 1.3272 A, P ~= 36.07 W, converged=False, iter=47
```

Conclusion: the official `110 kW` startup can reach and trigger the rated-voltage gate if the Cs-filled gap conductance trajectory is made more insulating, with `h_eq=200 W/m2/K` as a useful first candidate. However, the fixed-voltage branch is still too slow and marginal for long steady continuation: even with the rebuilt pyd and augmented lookup, a `2 s` smoke took several minutes and the final record was close to target but non-converged. Before launching a long fixed-U steady run, improve fixed-U iteration performance/robustness or add a staged voltage-control transition.

### 2026-06-29 fixed-U startup switch optimization after h_eq=200 smoke

The first `h_eq=200` automatic fixed-R to fixed-U smoke proved the thermal path can trigger the `27.2 V` gate, but fixed-U was too slow and marginal:

```text
baseline output = testModule/v13_start_h200_fixedr_to_fixedu_2s_20260629
wall ~= 390 s for 2 s physical time
fixed-U records: iter 29 / 42 / 47
last record: U ~= 27.178 V, converged=False
```

Root cause: after switching from `R_total=100 ohm`, the fixed-R current is only about `0.298 A`, while the fixed-U operating current is near `1 A`. The series fixed-U solver was using the switch current as its first guess and then spending many full `circuitCalc()` calls in broad bracket/secant search.

Source-side change in `ThermoCalc/circuitTECs.cpp::uFixedCircuitCalc()`:

- keep the bounded bracket/secant structure and low-temperature no-hang guards;
- use `0.05 V` as the fixed-voltage engineering residual tolerance;
- prioritize candidate currents around `I_guess + 1 A` before the wider fallback samples, which matches the V13 switch current jump;
- keep fixed-U public output semantics as `Uout=Utarget` on convergence.

Verification used only the rebuilt test pyd in `ThermoCalc/build_cp312/Release` through `THERMOCALC_PYD_DIR`; the root production pyd was not overwritten.

Regression checks:

```text
testModule/test_thermocalc_interface.py: passed
testModule/test_thermocalc_lookup.py: passed, lookup speedup ~= 37.6x
```

Optimized smoke:

```text
output = testModule/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629
wall ~= 88 s for 2 s physical time
switch at t=5432.5 s from fixed_r U ~= 29.767 V

fixed-U records:
t=5433.0 s: U=27.2 V, I ~= 0.9644 A, P ~= 26.23 W, converged=True, iter=15
t=5433.5 s: U=27.2 V, I ~= 1.1383 A, P ~= 30.96 W, converged=True, iter=3
t=5434.0 s: U=27.2 V, I ~= 1.3208 A, P ~= 35.93 W, converged=True, iter=4
```

Current implication: the official `110 kW`, `h_eq=200 W/m2/K` startup path is now numerically viable for a short fixed-R to fixed-U transition smoke with the rebuilt test pyd. It is still not ready to claim steady state: the root production pyd has not been replaced, the hydraulic solver still reports the known first-step residual warning, and a longer fixed-U continuation is needed to verify lookup coverage and stable energy balance.

### 2026-06-29 fixed-U short continuation after sample-order optimization

After the `I_guess + 1 A` prioritized sample order in `uFixedCircuitCalc()`, the `h_eq=200` fixed-R to fixed-U smoke was repeated:

```text
output = testModule/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629
wall ~= 88 s for 2 s physical time
switch at t=5432.5 s from fixed_r U ~= 29.767 V
fixed-U records all converged:
  t=5433.0 s: U=27.2 V, I ~= 0.9644 A, P ~= 26.23 W, iter=15
  t=5433.5 s: U=27.2 V, I ~= 1.1383 A, P ~= 30.96 W, iter=3
  t=5434.0 s: U=27.2 V, I ~= 1.3208 A, P ~= 35.93 W, iter=4
```

A fixed-voltage continuation from that restart also completed:

```text
restart = testModule/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629/v13_start_h200_fixedr_to_fixedu_2s_plus1first_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus5s_20260629
mode = fixed_u
target voltage = 27.2 V
wall ~= 163 s for 5 s physical time

records:
t=5435 s: U=27.2 V, I ~= 1.4214 A, P ~= 38.66 W, iter=16, converged=True
t=5436 s: U=27.2 V, I ~= 1.8037 A, P ~= 49.06 W, iter=3, converged=True
t=5437 s: U=27.2 V, I ~= 2.2148 A, P ~= 60.24 W, iter=3, converged=True
t=5438 s: U=27.2 V, I ~= 1.8341 A, P ~= 49.89 W, iter=1, converged=True
t=5439 s: U=27.2 V, I ~= 1.7341 A, P ~= 47.17 W, iter=7, converged=True
```

Current status: with the rebuilt test pyd, the official `110 kW`, `h_eq=200 W/m2/K` startup route can pass the fixed-R to fixed-U transition and sustain a short fixed-U continuation. This is still not a steady result. The next calculation should extend fixed-U in moderate chunks while monitoring lookup coverage, hydraulic residuals, current/power oscillation, and energy balance before attempting an overnight run.

### 2026-06-29 fixed-U plus20s continuation status

The optimized fixed-U path was continued from `t=5439 s` for another `20 s`:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus5s_20260629/v13_start_h200_fixedu27p2_plus5s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus20s_20260629
mode = fixed_u
target voltage = 27.2 V
thermo update interval = 0.5 s
wall ~= 947 s for 20 s physical time
```

All 10 recorded fixed-U points converged. Endpoint:

```text
t = 5459 s
U = 27.2 V
I ~= 1.1384 A
P_e ~= 30.97 W
mean emitter ~= 1240.59 K
core inlet/outlet ~= 745.80 / 838.74 K
q_radiator_total ~= 107.46 kW
coolant enthalpy rise ~= 105.34 kW
core_heat - coolant_enthalpy - electric ~= 4.63 kW
```

Last-10-s trends:

```text
dI/dt ~= -2.69e-2 A/s
dP_e/dt ~= -0.732 W/s
dT_emitter/dt ~= +5.18e-2 K/s
dT_inlet/dt ~= -1.06e-1 K/s
dq_rad/dt ~= -72.5 W/s
```

Interpretation: the fixed-U continuation is numerically stable over this medium segment, but it is not near steady. Electrical output is still relaxing downward while the emitter temperature is slowly recovering. Continue in moderate chunks before attempting an overnight steady run.

### 2026-06-29 fixed-U plus40s/plus60s continuation status

The optimized `h_eq=200 W/m2/K`, `fixed_u=27.2 V` path was continued in two additional 20 s chunks, still with `thermo_update_interval=0.5 s` and the rebuilt test pyd.

`plus40s` segment:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus20s_20260629/v13_start_h200_fixedu27p2_plus20s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus40s_20260629
wall ~= 490 s for 20 s physical time
all records converged
endpoint t=5479 s:
  U = 27.2 V
  I ~= 1.0404 A
  P_e ~= 28.30 W
  mean emitter ~= 1240.80 K
  core inlet/outlet ~= 745.34 / 838.81 K
  q_radiator_total ~= 107.46 kW
  coolant enthalpy rise ~= 105.93 kW
  core_heat - coolant_enthalpy - electric ~= 4.04 kW
last-10-s trends:
  dI/dt ~= -1.93e-2 A/s
  dP_e/dt ~= -0.525 W/s
  dT_emitter/dt ~= +4.67e-2 K/s
```

`plus60s` segment:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus40s_20260629/v13_start_h200_fixedu27p2_plus40s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus60s_20260629
wall ~= 856 s for 20 s physical time
all records converged
endpoint t=5499 s:
  U = 27.2 V
  I ~= 1.0344 A
  P_e ~= 28.13 W
  mean emitter ~= 1240.94 K
  core inlet/outlet ~= 745.32 / 838.84 K
  q_radiator_total ~= 107.45 kW
  coolant enthalpy rise ~= 106.00 kW
  core_heat - coolant_enthalpy - electric ~= 3.97 kW
last-10-s trends:
  dI/dt ~= -1.71e-2 A/s
  dP_e/dt ~= -0.464 W/s
  dT_emitter/dt ~= +4.50e-2 K/s
```

Interpretation: fixed-U electrical convergence is now robust over the tested `60 s` continuation after switching, but the coupled thermal state is not steady. The residual thermal imbalance remains about `4 kW`, and the emitter is still rising slowly while electric output relaxes downward. The run is also still expensive at `0.5 s` TEC update frequency. Before an overnight run, either accept the cost and continue in longer chunks, or test a larger TEC update interval / further C++ inner-loop optimization.

### 2026-06-29 fixed-U TEC update interval 1.0 s check

A speed/accuracy check was run from the `t=5499 s` fixed-U restart by increasing the ThermoCalc update interval from `0.5 s` to `1.0 s` while keeping `max_dt=0.5 s`:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus60s_20260629/v13_start_h200_fixedu27p2_plus60s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus80s_tec1s_20260629
duration = 20 s
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
wall ~= 169 s for 20 s physical time
```

All records converged. Endpoint:

```text
t = 5519 s
U = 27.2 V
I ~= 1.0645 A
P_e ~= 28.95 W
mean emitter ~= 1241.06 K
core inlet/outlet ~= 745.35 / 838.87 K
q_radiator_total ~= 107.47 kW
coolant enthalpy rise ~= 106.00 kW
core_heat - coolant_enthalpy - electric ~= 3.97 kW
```

Last-10-s trends:

```text
dI/dt ~= -1.50e-2 A/s
dP_e/dt ~= -0.408 W/s
dT_emitter/dt ~= +4.41e-2 K/s
dT_inlet/dt ~= +5.06e-2 K/s
dq_rad/dt ~= -36.4 W/s
```

Comparison to the previous `0.5 s` segment indicates the physical trend is consistent, while wall time improved substantially. The case is still not steady, but `thermo_update_interval=1.0 s` is a reasonable setting for the next longer fixed-U continuation.

### 2026-06-29 fixed-U plus180s with 1.0 s TEC update

A longer fixed-U continuation was run after accepting `thermo_update_interval=1.0 s` as a speed-improving setting:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus80s_tec1s_20260629/v13_start_h200_fixedu27p2_plus80s_tec1s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus180s_tec1s_20260629
duration = 100 s
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
wall ~= 1267 s for 100 s physical time
```

All 10 recorded points converged. Endpoint:

```text
t = 5619 s
U = 27.2 V
I ~= 0.9685 A
P_e ~= 26.34 W
mean emitter ~= 1243.84 K
core inlet/outlet ~= 745.63 / 840.81 K
q_radiator_total ~= 107.98 kW
coolant enthalpy rise ~= 107.87 kW
core_heat - coolant_enthalpy - electric ~= 2.11 kW
```

Last-50-s trends:

```text
dI/dt ~= 0 A/s
dP_e/dt ~= 0 W/s
dT_emitter/dt ~= +2.92e-2 K/s
dT_inlet/dt ~= +1.53e-2 K/s
dq_rad/dt ~= +8.59 W/s
```

Interpretation: the fixed-voltage electrical solution is now stable over a `100 s` continuation and electric output has flattened near `26.3 W` for this low-power startup state. The thermal system is still not steady: emitter temperature is still rising and the residual heat balance remains about `2.1 kW`. Continue in longer but still bounded chunks, for example `300-500 s`, before declaring a near-steady startup fixed-U state or launching an overnight run.

### 2026-06-29 fixed-U plus480s near-steady continuation

A `300 s` continuation was run from the `t=5619 s` restart using the accepted `thermo_update_interval=1.0 s` setting:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus180s_tec1s_20260629/v13_start_h200_fixedu27p2_plus180s_tec1s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus480s_tec1s_20260629
duration = 300 s
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
wall ~= 5618 s for 300 s physical time
```

All 10 recorded points converged. Endpoint:

```text
t = 5919 s
U = 27.2 V
I ~= 0.9552 A
P_e ~= 25.98 W
mean emitter ~= 1248.11 K
core inlet/outlet ~= 747.86 / 843.78 K
q_radiator_total ~= 109.24 kW
coolant enthalpy rise ~= 108.72 kW
core_heat - coolant_enthalpy - electric ~= 1.26 kW
```

Last-150-s trends:

```text
dI/dt ~= +2.08e-4 A/s
dP_e/dt ~= +5.66e-3 W/s
dT_emitter/dt ~= +9.72e-3 K/s
dT_inlet/dt ~= -4.22e-3 K/s
dq_rad/dt ~= -2.06 W/s
d(coolant enthalpy rise)/dt ~= +5.95 W/s
```

Interpretation: this is close to a stable fixed-voltage startup state but not a strict steady state. Electrical output is nearly flat near `26 W`, while the thermal system still has about `1.26 kW` residual heat imbalance and a small positive emitter-temperature drift. Another several-hundred-second continuation should reduce the remaining drift, but the current run speed is about `18.7 wall-s / physical-s`, so an overnight run is appropriate only if this cost is acceptable.

### 2026-06-29 fixed-U plus1780s continuation

A further `1000 s` official continuation was run from the `t=6219 s` restart:

```text
restart = testModule/v13_start_h200_fixedu27p2_plus780s_tec1s_20260629/v13_start_h200_fixedu27p2_plus780s_tec1s_20260629_latest_restart.npz
output = testModule/v13_start_h200_fixedu27p2_plus1780s_tec1s_20260629
duration = 1000 s
startup thermal power = 110000 W
mode = fixed_u
target voltage = 27.2 V
thermo_update_interval = 1.0 s
```

Endpoint:

```text
t = 7219 s
U = 27.2 V
I ~= 1.3873 A
P_e ~= 37.74 W
mean emitter ~= 1247.23 K
core inlet/outlet ~= 746.84 / 842.74 K
q_radiator_total ~= 108.69 kW
coolant enthalpy rise ~= 108.70 kW
core_heat - coolant_enthalpy - electric ~= 1.27 kW
```

All recorded fixed-U points converged. The only stderr message was the known first-step hydraulic residual warning at the restart boundary. Over the last five recorded points (`6819-7219 s`), the drift was small but not zero:

```text
dP_e/dt ~= -1.86e-5 W/s
dT_emitter/dt ~= +1.21e-4 K/s
dT_inlet/dt ~= +8.15e-5 K/s
dq_radiator/dt ~= +4.89e-2 W/s
d(coolant enthalpy rise)/dt ~= +4.72e-2 W/s
```

Interpretation: the official `110 kW`, `h_eq=200 W/m2/K`, `fixed_u=27.2 V` startup state is numerically stable and close to steady in temperatures and electrical output, but the residual storage term is still about `1.27 kW`. A longer `5000 s` continuation was launched from this restart to reduce the remaining thermal drift:

```text
output = testModule/v13_start_h200_fixedu27p2_plus6780s_tec1s_20260629
pid = 58580
```
### 2026-06-29 fixed-U plus6780s long continuation

A `5000 s` official continuation from the `t=7219 s` restart completed under the formal `110 kW` startup power, `h_eq=200 W/m2/K`, `fixed_u=27.2 V`, and `thermo_update_interval=1.0 s` settings:

```text
output = testModule/v13_start_h200_fixedu27p2_plus6780s_tec1s_20260629
t = 12219 s
U = 27.2 V
I ~= 1.3838 A
P_e ~= 37.64 W
mean emitter ~= 1247.53 K
core inlet/outlet ~= 747.12 / 843.17 K
q_radiator_total ~= 108.86 kW
coolant enthalpy rise ~= 108.86 kW
core_heat - coolant_enthalpy - electric ~= 1.10 kW
```

All history records converged in fixed-U mode with `iteration_count = 1` after the restart transient. The only stderr entry was the known first-step hydraulic residual warning at `t=7219 s`. Recent trends over the last five records were still positive but small:

```text
dP_e/dt ~= -1.86e-5 W/s
dT_emitter/dt ~= +5.86e-5 K/s
dT_inlet/dt ~= +5.41e-5 K/s
dq_radiator/dt ~= +3.27e-2 W/s
d(coolant enthalpy rise)/dt ~= +3.27e-2 W/s
```

Interpretation: the official fixed-U startup path remains stable and near steady, but the remaining storage term is still about `1.10 kW`. A longer `30000 s` steady approach run was launched from this restart:

```text
output = testModule/v13_start_h200_fixedu27p2_plus36780s_tec1s_20260629
pid = 19504
```
### 2026-06-29 fixed-U plus36780s first history checkpoint

The `30000 s` long steady-approach run has written its first formal history record, confirming that the restart-load `tec_solver_converged=False` flag was only an initial latest-state artifact. The first recorded point is converged:

```text
output = testModule/v13_start_h200_fixedu27p2_plus36780s_tec1s_20260629
pid = 19504
record = first history row
t = 14219 s
relative time = 2000 s
startup thermal power = 110000 W
mode = fixed_u
converged = True
iteration_count = 1
U = 27.2 V
I ~= 1.3826 A
P_e ~= 37.61 W
mean emitter ~= 1247.64 K
core inlet/outlet ~= 747.22 / 843.32 K
q_radiator_total ~= 108.92 kW
coolant enthalpy rise ~= 108.92 kW
core_heat - coolant_enthalpy - electric ~= 1.04 kW
```

Drift from the previous segment endpoint (`12219 -> 14219 s`):

```text
dP_e/dt ~= -1.60e-5 W/s
dT_emitter/dt ~= +5.46e-5 K/s
dT_inlet/dt ~= +5.04e-5 K/s
dq_radiator/dt ~= +3.05e-2 W/s
d(coolant enthalpy rise)/dt ~= +3.05e-2 W/s
```

Interpretation: the official `110 kW` fixed-U path is continuing normally and remains numerically stable. The thermal balance is still closing slowly, with the residual storage term reduced from about `1.10 kW` to about `1.04 kW` over this first `2000 s` checkpoint. Continue the long run before claiming strict steady state.
### 2026-06-29 correction: cesium TEC gap h_eq should be 29 W/m2/K

User review caught that the previous `h_eq=200 W/m2/K` fixed-U continuations used an artificial sensitivity value, not the physical cesium-vapor TEC gap setting. The V7 steady CaseA configuration uses `tec_gap_config h_eq=29.0 W/m2/K`, and the V13 cold-start cesium-filled gap should be consistent with that value unless explicitly running a sensitivity case.

Code correction:

```text
testModule/v13_startup_control.py: V13StartupControlConfig.cesium_gap_h_eq_w_m2_k default = 29.0
testModule/run_v13_start_case.py: --cesium-gap-h-eq-w-m2-k default = 29.0
testModule/test_v13_startup_control.py: added defaults tests for config and runner CLI
```

Verification:

```text
testModule/test_v13_startup_control.py: passed
python -m py_compile testModule/v13_startup_control.py testModule/run_v13_start_case.py testModule/test_v13_startup_control.py: passed
```

Therefore, the earlier `h_eq=200` results must be treated as numerical/control diagnostics only, not as official V13 physical startup results.

Corrected `h_eq=29` restart path:

```text
base restart = testModule/v13_start_cesium_conditioning_plus1000s_20260628/v13_start_cesium_conditioning_plus1000s_20260628_latest_restart.npz
base state: t=4090 s, Cs fraction ~= 1, TEC disabled, old h_eq ~= 250 W/m2/K
```

A corrected thermal hold with TEC disabled was run using `h_eq=29.0`:

```text
output = testModule/v13_start_h29_tec_off_200s_20260629
t = 4290 s
core power = 110000 W
TEC disabled
h_eq = 29.0 W/m2/K
core inlet/outlet ~= 722.28 / (history endpoint) K
mean emitter ~= 1513.22 K
```

Then a `2 s` fixed-R to fixed-U smoke was run:

```text
output = testModule/v13_start_h29_fixedr_to_fixedu_2s_20260629
fixed-R at t=4290.5 s: U ~= 65.04 V, I ~= 0.650 A, P_e ~= 42.3 W, converged=True
switch to fixed-U 27.2 V triggered immediately after the first fixed-R record
fixed-U at t=4292.0 s: U=27.2 V, I ~= 351.17 A, P_e ~= 9.55 kW, converged=True
```

Short fixed-U continuation:

```text
output = testModule/v13_start_h29_fixedu27p2_plus20s_20260629
t = 4312 s
U=27.2 V
I ~= 259.10 A
P_e ~= 7.05 kW
mean emitter ~= 1490.42 K
core inlet/outlet ~= 728.01 / 818.14 K
core_heat - coolant_enthalpy - electric ~= 0.76 kW
```

Medium fixed-U continuation:

```text
output = testModule/v13_start_h29_fixedu27p2_plus220s_20260629
t = 4512 s
U=27.2 V
I ~= 211.67 A
P_e ~= 5.76 kW
mean emitter ~= 1506.67 K
core inlet/outlet ~= 731.35 / 819.28 K
q_radiator_total ~= 99.65 kW
coolant enthalpy rise ~= 99.69 kW
core_heat - coolant_enthalpy - electric ~= 4.55 kW
```

A `1000 s` fixed-U continuation with `thermo_update_interval=1.0 s` completed:

```text
output = testModule/v13_start_h29_fixedu27p2_plus1220s_tec1s_20260629
t = 5512 s
U=27.2 V
I ~= 209.72 A
P_e ~= 5.70 kW
mean emitter ~= 1541.95 K
core inlet/outlet ~= 735.62 / 825.70 K
q_radiator_total ~= 102.08 kW
coolant enthalpy rise ~= 102.11 kW
core_heat - coolant_enthalpy - electric ~= 2.19 kW
```

Interpretation: with the corrected `h_eq=29.0 W/m2/K`, the V13 startup enters a meaningful TEC generation regime. The fixed-R voltage gate is crossed naturally, fixed-U solves converge, and electric power is now in the expected kilowatt range rather than the invalid tens-of-watts result from `h_eq=200`. The state is still not strict steady because the emitter and radiator heat rejection are drifting; a `5000 s` continuation is running:

```text
output = testModule/v13_start_h29_fixedu27p2_plus6220s_tec1s_20260629
pid = 53696
```
### 2026-06-29 correction: TFE ignition timing for cesium gap and fixed-R startup

User clarified the startup sequence: at `critical_time + 1500 s`, TFE ignition should immediately replace the emitter-collector gap equivalent heat-transfer coefficient with the cesium-vapor value `h_eq=29.0 W/m2/K`. The fixed-resistance external circuit should participate from this ignition point so voltage/current develop while the emitter warms; once the terminal voltage reaches `27.2 V`, the main circuit switches to fixed total voltage.

This supersedes the previous workflow that first ran a separate TEC-off thermal hold after cesium conditioning. Those TEC-off hold runs remain useful diagnostics, but are not the formal startup sequence.

Code updates:

```text
testModule/v13_startup_control.py
  - default cesium_gap_h_eq_w_m2_k = 29.0
  - TFE ignition latches by time after critical, not by emitter-temperature gate
  - once ignition latches, cs_fraction = 1.0 and h_eq immediately equals the cesium value
  - default electrical start gates are zero, so fixed-R TEC coupling starts at TFE ignition

testModule/run_v13_start_case.py
  - --cesium-gap-h-eq-w-m2-k default = 29.0
  - --tec-electrical-start-after-cesium-s default = 0.0
  - --tec-electrical-start-cs-fraction default = 0.0
  - --tec-electrical-start-emitter-temperature-k default = 0.0

testModule/test_v13_startup_control.py
  - added/updated tests for TFE ignition immediately setting h_eq=29 and enabling fixed-R TEC
```

Verification:

```text
testModule/test_v13_startup_control.py: passed
python -m py_compile testModule/v13_startup_control.py testModule/run_v13_start_case.py testModule/test_v13_startup_control.py: passed
```

Corrected sequence test from the `t=1590 s` pre-ignition restart:

```text
base restart = testModule/v13_start_corrected_shield_startup_1590s_20260628/v13_start_corrected_shield_startup_1590s_20260628_latest_restart.npz
base state: t=1590 s, time_after_critical ~= 1492.7 s, TEC off, h_eq=600 W/m2/K, mean emitter ~= 1031.2 K
```

With the corrected controller, TEC coupling enabled automatically at `t ~= 1597.342 s` and `h_eq=29.0` was applied from ignition.

Fixed-R load tests:

```text
R_total = 0.0044 ohm
output = testModule/v13_start_h29_tfe_ignition_fixedr_300s_20260629
stable fixed-R, but voltage only rose to ~= 2.23 V by t=1890 s; no switch.

R_total = 0.05 ohm
output = testModule/v13_start_h29_tfe_ignition_fixedr_R005_500s_20260629
stable fixed-R, voltage rose to ~= 16.40 V by t=2090 s; no switch.

R_total = 0.083 ohm from the warmed R=0.05 restart
output = testModule/v13_start_h29_fixedr_R0083_fromR005_200s_20260629
stable fixed-R, voltage rose to ~= 21.97 V by t=2290 s; no switch.

R_total = 0.105 ohm from the warmed R=0.083 restart
output = testModule/v13_start_h29_fixedr_R0105_fromR0083_100s_20260629
stable fixed-R, voltage rose to ~= 24.68 V by t=2390 s; no switch.

R_total = 0.12 ohm from the warmed R=0.105 restart
output = testModule/v13_start_h29_fixedr_R012_fromR0105_50s_20260629
stable fixed-R, voltage rose to ~= 26.20 V by t=2440 s; no switch.

R_total = 0.125 ohm from the warmed R=0.12 restart
output = testModule/v13_start_h29_fixedr_R0125_fromR012_20s_20260629
stable fixed-R, peak voltage ~= 26.72 V; no switch.

R_total = 0.131 ohm from the warmed R=0.125 restart
output = testModule/v13_start_h29_fixedr_R0131_fromR013_10s_20260629
fixed-R at t=2480.5 s: U ~= 27.246 V, I ~= 207.98 A, P_e ~= 5.67 kW
automatic switch to fixed-U 27.2 V succeeded; subsequent fixed-U records converged.
```

A direct `R_total=0.10 ohm` run from cold TFE ignition was attempted, but it consumed CPU without writing a first history record and was stopped. A `per_tec` interpretation of `0.0044 ohm` also failed by producing non-finite axial Joule heat after early records. Therefore, the currently stable route is staged fixed-R resistance from low value to approximately `0.131 ohm` as the emitter warms, followed by fixed-U.

Fixed-U continuation after successful switch:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus1000s_20260629
t = 3490 s
U = 27.2 V
I ~= 209.68 A
P_e ~= 5.70 kW
mean emitter ~= 1543.14 K
core inlet/outlet ~= 735.61 / 825.68 K
q_radiator_total ~= 102.07 kW
coolant enthalpy rise ~= 102.10 kW
core_heat - coolant_enthalpy - electric ~= 2.20 kW
```

The result is in the expected kilowatt range and comparable to the V11 electrical output scale, but it is not yet strict steady state. A corrected-sequence `5000 s` fixed-U continuation is now running:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus6000s_tec1s_20260629
pid = 77012
```
### 2026-06-29 corrected ignition fixed-U plus6000s checkpoint

The corrected-sequence fixed-U continuation from the successful `R_total=0.131 ohm` switch completed a further `5000 s` segment:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus6000s_tec1s_20260629
start = 3490 s
end = 8490 s
startup power = 110000 W
h_eq = 29.0 W/m2/K
mode = fixed_u
U = 27.2 V
```

Endpoint:

```text
t = 8490 s
I ~= 209.47 A
P_e ~= 5.70 kW
mean emitter ~= 1544.62 K
core inlet/outlet ~= 736.45 / 826.92 K
q_radiator_total ~= 102.55 kW
coolant enthalpy rise ~= 102.55 kW
core_heat - coolant_enthalpy - electric ~= 1.75 kW
```

All records converged with `tec_solver_iteration_count = 1` after restart. The only stderr entry was the known first-step hydraulic residual warning at the restart boundary. Recent slopes over the last five records:

```text
dP_e/dt ~= -8.16e-4 W/s
dT_emitter/dt ~= +3.10e-5 K/s
dT_inlet/dt ~= +5.07e-5 K/s
dq_radiator/dt ~= +2.93e-2 W/s
d(coolant enthalpy rise)/dt ~= +2.93e-2 W/s
```

Interpretation: the corrected TFE ignition path is stable and in the expected kilowatt electrical-output range, but it is still not a strict steady state because the residual storage term is about `1.75 kW`. A longer `30000 s` fixed-U continuation is now running:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629
pid = 79504
```
### 2026-06-29 corrected ignition long continuation stopped and residual diagnosis

The longer fixed-U continuation was stopped after the user decided that the remaining
long steady approach can be run manually:

```text
output = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629
stopped pid = 79504
history final time = 32490 s
history final relative time = 24000 s
latest saved restart = testModule/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629/v13_start_h29_fixedu_fromR0131_plus36000s_tec1s_20260629_latest_restart.npz
```

Final history row:

```text
core heat power = 110000.000 W
coolant enthalpy rise = 102978.091 W
TEC electric power = 5692.628 W
radiator heat rejection = 102977.955 W
  tube = 53330.259 W
  fin  = 49647.695 W
core inlet/outlet = 737.185 / 828.029 K
mean emitter = 1545.021 K
TEC mode = fixed_u, converged = true, iteration_count = 1
```

Energy diagnostics:

```text
radiator heat rejection - coolant enthalpy rise = -0.136 W
core heat - coolant enthalpy rise - TEC electric power = 1329.281 W
core heat - radiator heat rejection - TEC electric power = 1329.417 W
```

Interpretation:

- The radiator/coolant-loop balance is closed to sub-watt level, so the kilowatt-scale
  residual is not caused by the radiator heat rejection model or coolant enthalpy
  accounting.
- The TEC solver converged in one iteration and all reported TEC outputs are finite,
  so the residual is not a ThermoCalc fixed-U non-convergence symptom.
- The remaining about `1.33 kW` is the transient storage term in the still-warming
  core/TFE/structural solids. The last records still show slowly increasing inlet,
  outlet, emitter and radiator heat rejection values, so strict steady state has not
  been reached.

Recommended handling:

- Treat this run as a verified startup-to-fixed-U handoff, not as a final steady
  benchmark.
- Use the latest restart above for the user's own long continuation.
- Report `core_heat - coolant_dh - electric` as the transient solid-storage residual
  until its magnitude is small enough for the target use, for example below about
  `0.5 kW` or below `1%` of core power, together with small final-window slopes in
  `T_in`, `T_out`, emitter temperature and radiator heat rejection.

Restart timestamp note: the interrupted run wrote history records through `t = 32490 s`, but the latest saved restart file stores `System/global_time = 28490 s` with `System/last_dt = 0.5 s`. Continue from the restart path above when a restartable state is needed; use the CSV final row only as the latest diagnostic evidence from the interrupted run.

### 2026-06-29 residual diagnostic fields in V13 startup history

The V13 startup runner now records explicit derived residual fields in history CSV and latest-state output:

```text
core_heat_minus_coolant_enthalpy_minus_electric_w
core_heat_minus_radiator_minus_electric_w
radiator_minus_coolant_enthalpy_w
core_energy_storage_residual_rel
radiator_coolant_balance_rel
```

Use these fields for future long steady-approach checks instead of recomputing the residual manually. A large `core_heat_minus_coolant_enthalpy_minus_electric_w` with small `radiator_minus_coolant_enthalpy_w` and converged finite TEC output means the residual is dominated by transient solid/core storage, not by ThermoCalc fixed-U failure.
