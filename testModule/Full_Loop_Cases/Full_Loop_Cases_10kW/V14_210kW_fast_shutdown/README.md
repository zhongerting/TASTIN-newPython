# V14 210 kW Fast Shutdown

This case starts from the completed two-orbit fixed-power V14 state. At the
initial time it hands 210 kW to point kinetics and applies a persistent
`-2 dollar` external reactivity. The normal 2.46 kg/s coolant circulation,
heat pipes, radiator, orbital external heat, and fixed-voltage TEC circuit are
otherwise unchanged. If the TEC current falls to 0.01 A or below, electrical
coupling is opened while passive cesium-gap heat transfer remains active.

The default calculation lasts one orbital period. It writes five CSV history
tables, a shutdown-event restart, periodic restarts, and a final restart.
The same runner can continue from its own final restart; point-kinetics and
decay-heat state are restored instead of being initialized again, and
`shutdown_elapsed_s` remains continuous from the original scram.
Continuation uses `0.01 s` steps for its first `2 s` by default to settle the
rebuilt hydraulic state, then returns to the normal `0.05 s` step.

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_fast_shutdown\run_v14_210kw_fast_shutdown.py
```
