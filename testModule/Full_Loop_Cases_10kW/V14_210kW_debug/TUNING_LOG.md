# V14_10kW 210 kW Powered Tuning Log

This file summarizes the first powered V14_10kW debug/tuning pass.

## Fixed Conditions

```text
core power = 210000 W
loop flow = 2.46 kg/s
solid solver = implicit_euler
initial/restart temperature family = warm V14_10kW restart near 754 K
space temperature = 4 K
external orbital heat flux = disabled
TEC lookup = enabled for powered runs
main TEC target current = 206 A
main TEC target net power = 10.44 kW
core inlet target = 754.45 K
core outlet target = 845.65 K
```

## New Tuning Interfaces Used

```text
radiator_emissivity              heat-pipe condenser and fin emissivity together
hp_up_view_factor                upper-side heat-pipe radiation view factor
upper_hp_down_view_factor         lower-side view factor for upper collector ring
lower_hp_down_view_factor         lower-side view factor for lower collector ring
wire_resistance_scale             scale on the four preserved wire resistance ratios
tec_voltage                       main series TEC voltage
```

The tuned view-factor baseline is:

```text
hp_up_view_factor = 0.0
upper_hp_down_view_factor = 0.3
lower_hp_down_view_factor = 0.4
```

## Main Runs

| Run | Key parameters | Final result | Decision |
| --- | --- | --- | --- |
| `tune_eps075_u50p7_wire020_600s` | eps 0.75, U 50.7 V, wire 0.20 | Tin 753.417 K, Tout 844.057 K, I 214.898 A, P 10.895 kW | Too cold and too much current |
| `tune_eps070_u50p7_wire033_300s_from2564` | eps 0.70, U 50.7 V, wire 0.33 | Tin 761.630 K, Tout 851.073 K, I 208.231 A, P 10.557 kW | Too hot |
| `tune_eps074_u50p7_wire033_300s_from2564` | eps 0.74, U 50.7 V, wire 0.33 | Tin 754.946 K, Tout 845.373 K, I 206.564 A, P 10.473 kW | Close, but still evolving |
| `tune_eps07425_u50p7_wire0335_600s_from2564` | eps 0.7425, U 50.7 V, wire 0.335 | Tin 754.756 K, Tout 845.286 K, I 206.209 A, P 10.455 kW | Good intermediate |
| `final_eps07425_u50p7_wire0335_3600s_from3164` | eps 0.7425, U 50.7 V, wire 0.335 | Tin 755.530 K, Tout 846.397 K, I 206.456 A, P 10.467 kW | Long-window run drifted hot |
| `correct_eps075_u50p7_wire0335_1200s_from6764` | eps 0.75, U 50.7 V, wire 0.335 | Tin 754.185 K, Tout 845.228 K, I 206.112 A, P 10.450 kW | Stable but too cold |
| `final_eps07475_u50p65_wire0415_1200s_from7964` | eps 0.7475, U 50.65 V, wire 0.415 | Tin 754.850 K, Tout 845.936 K, I 201.644 A, P 10.213 kW | Wire scale too large |
| `final_eps07475_u50p65_wire0335_1200s_from7964` | eps 0.7475, U 50.65 V, wire 0.335 | Tin 754.738 K, Tout 845.773 K, I 206.569 A, P 10.463 kW | Current best baseline |

Short 1 s checks were not used as final evidence because the TEC update interval and electrical/thermal relaxation can make the first refreshed electrical state misleading. The final decision is based on 1200 s or longer windows.

## Current Best Baseline

Path:

```text
testModule/Full_Loop_Cases_10kW/V14_210kW_debug/runs/final_eps07475_u50p65_wire0335_1200s_from7964
```

Restart:

```text
testModule/Full_Loop_Cases_10kW/V14_210kW_debug/runs/final_eps07475_u50p65_wire0335_1200s_from7964/stage_01_restart.npz
```

Final metrics:

```text
time = 9163.85 s
core_inlet_T = 754.738 K
core_outlet_T = 845.773 K
core_delta_T = 91.034 K
TEC current = 206.569 A
TEC net electric power = 10462.741 W
upper heat-pipe rejection = 85317.474 W
lower heat-pipe rejection = 103618.741 W
heat-pipe rejection total = 188936.215 W
ring-wall rejection = 6275.383 W
total external rejection = 195211.598 W
pump total required head = 27880.038 Pa
minimum fluid temperature = 721.290 K
fluid converged = true
TEC converged = true
```

Residual against targets:

```text
core inlet: +0.288 K
core outlet: +0.123 K
current: +0.569 A
net electric power: +22.7 W
```

## Suggested Next Adjustment

If a tighter calibration is needed, the next conservative move is a small increase in radiator emissivity from `0.7475` or a small increase in wire scale from `0.335`, followed by another 1200 s or longer validation. A 1 s electrical check is not sufficient for final acceptance.
