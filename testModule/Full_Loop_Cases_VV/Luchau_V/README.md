# Luchau_V Single TFE Model

This case builds one standalone full-length `TFEUnit` and one single-TFE `ThermoCalcModel` fixed-voltage electrical circuit.

## Geometry And Mesh

The model copies the TFE geometry and axial allocation used by `testModule/Full_Loop_Cases/common_core_builder.py`:

- axial lengths: `0.065 m` lower reflector, `0.377 m` active zone, `0.065 m` upper reflector
- axial nodes: `6 + 25 + 6 = 37`
- total length: `0.507 m`
- radial TFE dimensions: `r_pellet_inner=4.0e-3 m` through `r_moderator_outer=16.27e-3 m`

The standalone component is created with `strict_adiabatic_single_tfe=True`, so it models the TFE body and coolant channel without wrapping it in a `ReactorCore` or adding global moderator/barrel/reflector rings.

## Runtime Parameters

These values are intentionally required at runtime:

- `thermal_power_w`: total central heater power, distributed uniformly over the center `0.30 m`
- `target_voltage_v`: single-TFE fixed terminal voltage for ThermoCalc

Fixed case inputs:

- coolant: `SodiumPotassium78`
- coolant inlet temperature: `727 K`
- coolant mass flow: `1.3 / 37 kg/s`
- Cs pressure: `0.4 torr`, converted to `Tcs` with `Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)`
- TEC gap equivalent heat transfer coefficient: `29 W/m2/K`

## Entry Points

- `luchau_single_tfe_model.py`: model builder and ThermoCalc setup helper
- `run_luchau_single_tfe.py`: CLI runner; requires `--thermal-power-w` and `--target-voltage-v`
- `test_luchau_single_tfe_model.py`: focused unit tests


The runner can also perform coupled transient advancement. When `--steady` is supplied, `--duration-s` is treated as the maximum simulated time and the run stops early only after `--steady-window-steps` consecutive steps satisfy `max(abs(dT))/dt <= --steady-dtemp-k-s`.

Example coupled transient with a smaller time step:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_VV\Luchau_V\run_luchau_single_tfe.py --thermal-power-w 3260 --target-voltage-v 0.30 --duration-s 50 --dt-s 0.02 --steady --steady-dtemp-k-s 1e-3 --steady-window-steps 20
```

Current verification note: the `3260 W`, `0.30 V`, `dt=0.02 s`, `50 s` coupled run remains transient (`steady_reached=false`) and ThermoCalc reports repeated C++ convergence failures after the emitter exceeds the low-temperature zero-emission cutoff. Treat that run as a diagnostic, not a valid steady-state operating point.
Example smoke without running the C++ ThermoCalc solve:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_VV\Luchau_V\run_luchau_single_tfe.py --thermal-power-w 3000 --target-voltage-v 1.2 --skip-thermocalc-calc --skip-transient
```
