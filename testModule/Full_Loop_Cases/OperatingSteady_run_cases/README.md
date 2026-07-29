# V14/V15 Operating Steady-State Tuning

This directory isolates the long-running V14 and V15 operating-point tuning work.

## Fixed requirements

- Uniform initial temperature: `727 K`
- Fixed core power: `115 kW`; point kinetics disabled
- Initial flow control: `1.3 kg/s`
- Solid conduction: `implicit_euler`
- Main TEC: 34 series units, lookup enabled, `fixed_u=27.2 V`
- Reserved TEC: 3 parallel units, lookup enabled, `parallel_fixed_u=0.35 V`
- Wire resistance base vector: `[0.001552, 0.001024, 0.000336, 0.000608] ohm`
- Wire resistance tuning uses one scale factor and preserves the four-value ratio
- Initial tuning excludes orbital heat and the thermal shield; radiator background is `4 K`

## Acceptance targets

- Core inlet: `727 +/- 8 K`
- Core outlet: `823 +/- 8 K`
- Main-series terminal power: `5.0 kW +/- 8%`
- Fixed-head verification flow: `1.3 kg/s +/- 1%`
- Reserved parallel circuit is diagnostic-only apart from its `0.35 V` target

V14 and V15 have independent emissivity, wire-scale, and pump-head results. Final orbital runs use direct N18 and N78 external heat respectively, without a shield.

## Final retained results

V14 uses emissivity `0.84`, total pump head `8516.4489 Pa`, and direct N18 heat.
V15 uses surface emissivity `0.815`, total pump head `40083.2884 Pa`, and direct
N78 heat. Both use wire scale `1.5`.

Only the final three-period confirmation directories remain under each `runs/`
directory. Each retains one final restart, the final-period key-parameter CSV,
summary/configuration, and debug logs. V15 mean-solid drift fell to `0.060 K/cycle`;
V14's local cold-solid drift fell to `0.513 K/cycle`.
