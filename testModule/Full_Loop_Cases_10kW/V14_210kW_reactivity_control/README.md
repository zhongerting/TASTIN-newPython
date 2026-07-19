# V14 210 kW Temperature-Feedback Reactivity Control

This runner converts an existing 210 kW fixed-power steady restart into a
point-kinetics run. The control drum is disabled and external reactivity is
fixed at zero. Therefore, the only reactivity change during this baseline run
is the temperature-feedback increment relative to the loaded steady state.

## Handoff definition

For a restart made by `V14_210kW_debug`:

1. Rebuild the model from the sibling `run_config.json` and load the restart.
2. Do not reapply the fixed 210 kW source.
3. Initialize the point reactor at the loaded total core power.
4. Calibrate the current fuel, electrode, moderator, and reflector feedback as
   the zero-reactivity reference.
5. Advance every step with `reactivity_control=0.0`.

For a restart previously written by this runner, the saved point-reactor state
and feedback reference are preserved. They are not initialized again.

The runner rejects a restart if its point-kinetics state disagrees with the
sibling configuration or if the control drum is enabled.

## Run

`--restart-in` is required. A typical baseline command is:

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' testModule\Full_Loop_Cases_10kW\V14_210kW_reactivity_control\run_v14_210kw_reactivity_control.py --restart-in testModule\Full_Loop_Cases_10kW\V14_210kW_debug\runs\final_eps07475_u50p65_wire0335_1200s_from7964\stage_01_restart.npz --output-dir testModule\Full_Loop_Cases_10kW\V14_210kW_reactivity_control\runs\baseline_10s --duration 10 --dt 0.05 --record-interval 1 --checkpoint-interval 10
```

To continue from its output, use the generated `stage_01_restart.npz` as the
next `--restart-in` and choose a new output directory.

## Outputs and diagnostics

The output directory contains `history.csv`, `run_config.json`,
`run_summary.json`, `latest_state.json`, checkpoints, and
`stage_01_restart.npz`.

Reactivity diagnostics distinguish:

- absolute temperature-feedback components;
- the saved steady-state feedback reference;
- effective temperature feedback, defined as current minus reference;
- external reactivity, fixed at zero in this version;
- control-drum reactivity, required to remain zero;
- total reactivity represented by the current thermal state. After a completed
  step this is the feedback available to the next coupling advance, not a claim
  about the previous step's already-integrated input.

This first version intentionally has no accident-reactivity CLI option. Add an
external reactivity history only after the zero-input handoff shows negligible
initial feedback and acceptable 210 kW power stability.
