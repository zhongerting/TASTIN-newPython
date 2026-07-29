# V14 operating steady calibration

This directory is isolated from earlier run outputs. `run_v14_operating_steady.py`
supports four restartable stages:

1. `thermal-flow`: 727 K, 115 kW, fixed 1.3 kg/s, TEC/external heat off.
2. `thermal-head`: rebuild with two equal fixed-head pumps and calibrate the total
   head back to 1.3 kg/s.
3. `coupled`: enable dense-v2 lookup TEC, main series fixed at 27.2 V and the
   three reserved parallel TECs fixed at 0.35 V.
4. `external`: keep the selected parameters and enable direct N18 orbital heat
   for at least `2 * 6552 = 13104 s`.

All solids use `implicit_euler`; fluid-solid couplers use `local_implicit`.
Heat-pipe and fin emissivity move together. Collector-ring wall emissivity is
fixed at 0.2. Wire resistance is always
`wire_scale * [0.001552, 0.001024, 0.000336, 0.000608] ohm` for both circuits.

Each stage writes `history.csv`, `latest_state.json`, `latest_restart.npz`,
periodic restart files, and `summary.json` below its own output directory.

Only the cleaned final confirmation run is retained under `runs/`; see `RESULTS.md`.
