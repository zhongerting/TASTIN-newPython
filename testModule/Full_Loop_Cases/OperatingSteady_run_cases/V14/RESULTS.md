# V14 final operating result

## Configuration

- Core power: `115 kW`
- Pump total head: `8516.4489 Pa`
- Heat-pipe/fin emissivity: `0.84`
- Ring-wall emissivity: `0.2`
- Wire-resistance scale: `1.5`
- Main/reserved TEC voltage: `27.2/0.35 V`
- Direct N18 external heat, no shield

## Final period

The retained final period is `33367-39919 s`.

| Metric | Period mean |
| --- | ---: |
| Core inlet temperature | `730.720 K` |
| Core outlet temperature | `825.875 K` |
| Minimum fluid temperature | `717.916 K` |
| Radiator rejection | `110.099 kW` |
| External heat | `2.289 kW` |

The coldest-solid boundary drift decreased to `-0.513 K/cycle`. Fluid, hydraulic,
radiator, and TEC results are periodic-near-steady; the local cold solid is still
converging slowly.

Retained artifacts:

```text
runs/external_N18_additional_3periods/final_restart.npz
runs/external_N18_additional_3periods/last_cycle_history.csv
runs/external_N18_additional_3periods/summary.json
runs/external_N18_additional_3periods/run.out
runs/external_N18_additional_3periods/run.err
```
