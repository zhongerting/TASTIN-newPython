# V14 20% Fixed-I Low-Power Design

## Goal

Find an independent V14 low-power quasi-steady point with orbital external heat enabled,
main TEC current fixed at the source-checkpoint value, terminal electric power between
2.0 kW and 2.2 kW, core thermal power strictly below 138.6 kW, and final-window
variation within +/-3%.

## Selected Approach

Use the validated parameterized fixed-I runner as the numerical engine, but create a
new 20% case directory, candidate definitions, histories, restarts, and acceptance
analysis. Every candidate starts from
`V14_210kW_fixed_power_external_heat_2orbits/runs/two_orbits_from13864_20260720/checkpoint_t019265s.npz`;
no 40% candidate restart is allowed.

Rejected alternatives:

- Continuing a 40% restart is faster but does not produce an independent 20% case.
- Copying the full fixed-I runner creates a large duplicate with no physical benefit.

## Controls

- Main TEC current: source-checkpoint refreshed value, expected
  `213.4691467366893 A`.
- Primary control: final core thermal power, always `<138600 W`.
- Auxiliary control: total coolant flow; first batch keeps `2.46 kg/s`.
- Transition: 20 s initial hold followed by a smooth quintic thermal-power ramp.
- External heat: enabled with the source checkpoint's phase origin and period.

## Candidate Workflow

1. Run three thermal-power candidates from the same immutable checkpoint.
2. Hold each endpoint until the final 300 s window is sufficiently settled; no full
   orbit is required during screening.
3. Interpolate the next thermal-power bracket from the measured electric response.
4. Adjust flow only if the electric target is bracketed but thermal margins or
   radiator temperature require correction.
5. Run the selected point longer and qualify it over a final 600 s window; extend to a
   full external-heat period only after the endpoint is selected.

## Acceptance

- Electric power: every final-window record is `2000..2200 W`.
- Fixed current: TEC converged, no fallback, current equals the fixed target within the
  existing runner tolerance.
- Quasi-steady variation: for electric power and each positive tracked temperature,
  `(max-min)/(2*mean) <= 0.03` in the final window.
- Drift guard: absolute fitted end-to-end electric-power drift over the final window is
  at most 1% of the mean.
- Existing hard limits remain active: fuel 2700 K, collector 1500 K, emitter 3000 K,
  coolant 1058 K, moderator 930 K, reflector 1000 K.
- Hydraulic and TEC solves must converge and all recorded values must be finite.

## Outputs

The independent case directory contains the candidate JSON, a small acceptance
analyzer, tests, Chinese README, run histories, restarts, and a final parameter
summary. Existing 40% files and processes are not modified.

## No-External-Heat Retuning And Periodic Qualification

The half-orbit check rejected the 118.5 kW candidate because electric power reached
1.9118 kW. Retuning starts from that run's final_restart.npz. External heat is disabled
during each 1500 s screening iteration; total flow remains 2.46 kg/s and TEC current
remains 213.4691467366893 A.

Start at 122 kW. At the end of each iteration, increase thermal power when terminal
electric power is below 2.0 kW and decrease it when above 2.2 kW. Reuse the preceding
iteration's final restart so the thermal state continues to relax. Stop screening when
the final 300 s is within 2.0-2.2 kW, has electric half-range at most 3%, electric
drift at most 1%, finite fields, converged TEC/hydraulics, and no thermal trip.

Then enable periodic external heat without changing thermal power, flow, or current,
and run two complete 5668.144369 s periods. The new periodic phase starts at zero at
the handoff restart. The final decision compares the two periods at matching phases and
reports the power band, phase-matched differences, convergence, and thermal limits.
The periodic run is not adjusted mid-orbit.
