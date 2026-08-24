# V14 210 kW half-radiator heat-pipe failure

This case uses the same accident controls as the three existing heat-pipe
failure cases. It starts from the specified 210 kW external-heat restart:

`V14_210kW_fixed_power_external_heat_2orbits/runs/two_orbits_from13864_20260720/checkpoint_t019865s.npz`

Failure definition:

- Upper and lower collector rings: sectors A1, A2, and A3 are fully failed.
- Six of the twelve physical radiator sectors are disabled, corresponding to
  171 of the nominal 340 heat pipes (about 50.3 percent).
- Affected fluid-to-evaporator heat-transfer coupling is zero.
- Nominal hydraulic resistance, flow areas, heat-pipe solid conduction,
  radiation areas, orbital external heat, and TEC calculation are retained.
- The thermal shield is disabled.

Controls and limits:

- Fixed core power is 210 kW until a solid-temperature limit is reached.
- A limit crossing switches to point kinetics and applies -2 dollars.
- Channel wall: 1058 K; fuel pellet: 2700 K; collector: 1023 K;
  moderator: 930 K; reflector: 1000 K.
- Coolant temperature is recorded but is not a trip criterion.
- History tables are written every 1 s; periodic restarts are written every
  100 s. Accident-start, scram-event, and final restart files are explicit.

Run with the required environment:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" `
  testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_heatpipe_half_radiator_failure\run_v14_heatpipe_half_radiator_failure.py
```

The default run is one orbital period (`5668.144369 s`). Use `--restart-in`
and `--output-dir` for a continuation without overwriting this run.
