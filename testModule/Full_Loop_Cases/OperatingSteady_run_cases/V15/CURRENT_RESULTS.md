# V15 final operating result

## Configuration

- Core power: `115 kW`
- Pump total head: `40083.2884 Pa`
- Tube/fin surface emissivity: `0.815`
- Effective emissivity: `0.476775`
- Wire-resistance scale: `1.5`
- Main/reserved TEC voltage: `27.2/0.35 V`
- Direct N78 external heat, no shield

## Final period

The retained final period is `32709-39261 s`.

| Metric | Period mean |
| --- | ---: |
| Core inlet temperature | `728.696 K` |
| Core outlet temperature | `823.586 K` |
| Minimum fluid temperature | `724.124 K` |
| Radiator rejection | `110.051 kW` |
| External heat absorption | `2.251 kW` |
| Main-series TEC power | `4.741 kW` |

Mean-solid boundary drift decreased to `-0.060 K/cycle`, establishing a good
periodic near-steady state.

Retained artifacts:

```text
runs/orbit_N78_additional_3periods/final_restart.npz
runs/orbit_N78_additional_3periods/last_cycle_history.csv
runs/orbit_N78_additional_3periods/summary.json
runs/orbit_N78_additional_3periods/run_config.json
runs/orbit_N78_additional_3periods/run.out
runs/orbit_N78_additional_3periods/run.err
```

The final restart retains six unused non-finite point-kinetics feedback-reference
placeholders. Point kinetics is disabled in this fixed-power run.
