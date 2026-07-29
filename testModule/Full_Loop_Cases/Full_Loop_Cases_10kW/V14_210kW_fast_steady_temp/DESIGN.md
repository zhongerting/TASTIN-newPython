# V14 210 kW Fast-Steady Temporary Case Design

## Scope

Create an isolated temporary case under this directory. Do not modify the
existing powered-debug case, reactivity-control case, common builders, shared
materials, or solvers.

The case reuses the existing V14 10 kW assembly and the current best 210 kW
restart. Orbital external heat remains disabled. Core power remains fixed at
210 kW while finding the thermal steady state.

## Accelerated model

After assembling and loading the existing case, wrap materials only on this
in-memory system so solid heat capacities are multiplied by 0.01. Apply the
scale to core solids, moderator, barrel, reflector, collector-ring solids,
heat-pipe walls, and the solid skeleton of heat-pipe wicks.

Keep NaK coolant and potassium working-fluid heat capacities unchanged. Keep
conductivity, density, emissivity, geometry, flow, TEC settings, radiation,
and fixed core power unchanged.

## Calculation flow

1. Load the current best physical 210 kW restart and sibling run configuration.
2. Apply the in-memory solid heat-capacity scale and refresh property caches.
3. Run fixed-power accelerated stages and save normal SystemManager restarts.
4. Judge convergence by comparing all fluid and solid temperature arrays over
   the final 10 accelerated seconds; require maximum absolute change below
   0.05 K.
5. Rebuild the same system with physical heat capacities, load the accelerated
   temperature/state restart, and run a 100 s fixed-power confirmation.
6. Accept the state only if physical-capacity temperatures remain stable and a
   subsequent 10 s zero-external-reactivity handoff has less than 1% power
   drift.

## Safety and outputs

Reject nonpositive heat-capacity scales and restarts without their sibling
configuration. Preserve the existing low-fluid-temperature stop. Write all
history, summaries, and restarts below this temporary directory's `runs/`
folder. Record the heat-capacity scale and scaled-solid count in run metadata.

## Verification

Use one small test to prove scaling changes solid capacitance by 0.01 without
changing conductivity or fluid heat capacity. Run syntax checks, a short
accelerated smoke, the accelerated steady calculation, physical-capacity
confirmation, and finally the existing zero-reactivity handoff.
