# V15 Operating-Steady Calibration

This directory is isolated from the existing V15/V71 runs. All generated logs, restart files,
fin histories and summaries stay under `runs/` here.

The runner uses the original V15 geometry and power profile with:

- uniform initial fluid/solid temperature `727 K`;
- fixed core power `115 kW`, point kinetics disabled;
- `implicit_euler` solids and `local_implicit` fluid-solid coupling;
- initial flow-controlled PumpA at `1.3 kg/s`, followed by two equal fixed-head pumps;
- main 34-TFE series circuit at `27.2 V`;
- reserved three-TFE parallel circuit at `0.35 V`;
- dense runtime v2 lookup and one common multiplier on the four wire resistances;
- one common tube/fin surface emissivity, with `epsilon_effective=0.585*epsilon_surface`;
- final direct N78 external heat, no shield.

Run the autonomous sequence with the required Python 3.12 environment:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" `
  testModule\Full_Loop_Cases\OperatingSteady_run_cases\V15\run_v15_operating_steady.py `
  auto --output-dir testModule\Full_Loop_Cases\OperatingSteady_run_cases\V15\runs\auto_727K
```

Individual restartable phases are `thermal-flow`, `fixed-head`, `coupled`, and `orbit`.
Each phase writes `history.csv`, `state_history/*.npz`, periodic checkpoints,
`final_restart.npz`, `run_config.json`, and `summary.json`. Quasi-steady fin temperatures,
radiation, absorption, root heat flow and radiation background are history only and are not
added to the global restart format.

Only the cleaned final confirmation run is retained under `runs/`; see
`CURRENT_RESULTS.md`.
