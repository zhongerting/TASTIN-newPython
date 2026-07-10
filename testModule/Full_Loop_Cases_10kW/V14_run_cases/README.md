# V14 10kW Run Cases

This directory contains the hydraulic-only V14 heat-pipe-radiator smoke case and its result file for `Full_Loop_Cases_10kW`.

## Current Scope

- Uses the five-ring, 58-TFE 10 kW core configuration.
- Applies no heating power.
- Disables TEC electrical calculation through `main_tec_enabled=False`.
- Does not advance solid heat conduction or thermal couplers.
- Only checks that the hydraulic network initializes and advances one hydraulic step.

## Radiator Topology

The smoke case uses explicit upper and lower heat-pipe collector rings. `HotOutletBranch_1/2/3` remains one shared set, and each hot branch splits into upper/lower inlet mix nodes. Upper and lower outlets merge into the shared `OutletMix_O1/O2/O3` nodes before entering `Manifold_1/2/3`. The old `multiplier=2` symmetric collector ring is not used.

- Upper heat-pipe count: `154`
- Lower heat-pipe count: `186`
- Heat-pipe lengths: `L_eva=0.0605 m`, `L_aba=0.0415 m`, `L_con=0.47 m`

## Command

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_run_cases.test_v14_flow_path_smoke
```

Result file: `v14_flow_path_smoke_result.json`

The result only represents the current 10 kW copy package under no-heating and no-electrical-calculation hydraulic connectivity conditions. It is not a thermal steady-state or powered operation result.
## Hydraulic Stability Run

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_10kW\V14_run_cases\run_v14_hydraulic_stability.py
```

Result file: `v14_hydraulic_stability_result.json`

The stability criterion is the maximum step-to-step `J_PumpA` flow change in the final window. The current threshold is `1.0e-6 kg/s`.
