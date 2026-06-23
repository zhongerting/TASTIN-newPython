# Thermionic Emission Diagnostic Scan

Date: 2026-06-22

This note documents the first-stage local thermionic-emission phase-map scan.
It does not change the default ThermoCalc production calculation path.

## Diagnostic Entry Point

The pybind module exposes:

```text
te_solver.calc_emission_point(TE, TC, Vo, Tcs, d_gap=0.5)
```

It evaluates one local `thermionicEmission` point through
`thermionicEmission::calcDiagnostics()` and returns:

```text
J
Vd
delta_V
phiE
phiC
regime
converged
finite_flag
iteration_count
obstructed_iterations
transition_iterations
saturation_iterations
obstructed_residual
transition_residual
saturation_residual
```

`regime` values:

```text
0  obstructed branch only
1  transition/saturation branch
2  blend zone near abs(delta_V) < 0.05
-1 failed or non-finite result
```

The existing production `thermionicEmission::calc()` path remains available and
unchanged for `singleThermionicEnergyConversion`.

## First-Stage Coarse Scan

Tools:

```text
ThermoCalc/tools/scan_emission_map.py
ThermoCalc/tools/analyze_emission_map.py
```

Default grid:

```text
TE: 800-2400 K, 25 points
TC: 600-1000 K, 15 points
Vo: 0-2.5 V, 25 points
Tcs: 500-800 K, 20 points
d_gap: 0.5
```

Total default points: 187,500.

Outputs are written under:

```text
ThermoCalc/emission_scan_outputs/
```

This directory is ignored by git because full scan outputs can become large.

## Intended Use

The scan is for phase-map and failure-boundary discovery. It is not a production
lookup table yet. Future lookup-table acceleration must remain optional and
fall back to the analytic solver for table misses, failed neighborhoods, or
non-finite interpolation results.

## Emission Database Generator

The larger optional lookup-table dataset is managed by:

```text
ThermoCalc/tools/emission_database.py
```

It keeps the analytic C++ solver as the source of truth and stores chunked
diagnostic outputs under:

```text
ThermoCalc/emission_database/
```

The default full plan covers four separately stored regions:

```text
core        TE 1300-2150 K, TC 700-900 K,  Vo 0-3.5 V, Pcs 0.5-2 torr
startup     TE 700-1300 K,  TC 500-800 K,  Vo 0-3.5 V, Pcs 0.5-2 torr
high_power  TE 2150-2400 K, TC 750-1000 K, Vo 0-3.5 V, Pcs 0.5-2 torr
accident    TE 700-2400 K,  TC 500-1100 K, Vo 0-3.5 V, Pcs 0.1-4 torr
```

Generate the full manifest and chunk plan:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py plan --db-dir ThermoCalc\emission_database --preset full
```

The current full plan contains 18,737,388 points in 78 chunks. Use multiple
workers by assigning stable worker indices:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py worker --db-dir ThermoCalc\emission_database --worker-id w0 --worker-index 0 --worker-count 6
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py worker --db-dir ThermoCalc\emission_database --worker-id w1 --worker-index 1 --worker-count 6
```

After chunks are generated, summarize and verify:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py summarize --db-dir ThermoCalc\emission_database --scan-chunks
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py verify --db-dir ThermoCalc\emission_database --samples 200
```

For quick regression of the generator itself, use the smoke preset:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py plan --db-dir ThermoCalc\emission_database\smoke --preset smoke
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py worker --db-dir ThermoCalc\emission_database\smoke --worker-id smoke --force
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py summarize --db-dir ThermoCalc\emission_database\smoke --scan-chunks
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py verify --db-dir ThermoCalc\emission_database\smoke --samples 7 --seed 22
```

## Current Database Status

As of 2026-06-23, the full local database under
`ThermoCalc/emission_database/` has been generated once:

```text
total points: 18,737,388
chunks: 78
done points: 18,737,388
converged points: 18,681,882
failed or non-finite points: 55,506
finite points: 18,737,388
zero-emission points: 9,465,664
J max: 58.189945 A/cm2
J mean: 0.697151 A/cm2
```

The generated chunks are intentionally ignored by git. They are runtime data,
not source. Keep `ThermoCalc/.gitignore` ignoring:

```text
emission_scan_outputs/
emission_scan_outputs_pressure_0p01_4torr/
emission_database/
```

## Optimized Lookup Sidecars

The original analytic scan contains non-converged points in startup and accident
regions. Direct analytic fallback is not useful there because the fallback often
fails for the same reason. The current table therefore supports optimized
sidecars:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py optimize-table --db-dir ThermoCalc\emission_database --zero-j-threshold 1e-3 --region startup --region accident
```

Current optimized result:

```text
regions optimized: startup, accident
raw invalid points: 55,506
zero-filled points: 43,104
neighbor-imputed points: 12,402
unresolved points: 0
```

The optimizer writes `.optimized.npz` files next to raw chunks. The Python
loader prefers these sidecars when present and otherwise loads the raw chunk.

## Runtime Lookup Export

The full analytic database keeps diagnostic fields needed for audit and
post-processing. It is not the preferred runtime artifact. Export a compact
runtime-only table after the raw chunks and optimized sidecars are ready:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py export-runtime --db-dir ThermoCalc\emission_database --out-dir ThermoCalc\emission_runtime_db --dtype float32 --zero-compress
```

Use `--region core`, `--region startup`, `--region high_power`, or
`--region accident` to export only selected scenarios. The runtime export keeps
only:

```text
TE_axis, TC_axis, Vo_axis, Tcs_axis
J, Vd, delta_V, phiE, phiC
lookup_safe, zero_mask
```

`J/Vd/delta_V/phiE/phiC` are stored as `float32` by default. `phiE/phiC` are
kept because upper-level boundary handling still uses them. `zero_mask` marks
safe near-zero current points; with `--zero-compress`, those `J` entries are
stored as exactly zero while the voltage and work-function fields remain
available for interpolation.

The runtime directory contains `runtime_manifest.json` and per-region
`*.runtime.npz` chunks. It is intentionally ignored by git:

```text
emission_runtime_db/
```

Runtime export stitches the first TE plane from the next source chunk onto the
current chunk when an old non-overlapping source plan is used. This turns legacy
source chunks such as `1300-1310 K` and `1320-1330 K` into runtime chunks such
as `1300-1320 K` and `1320-1340 K`, so interpolation no longer misses the
`1310-1320 K` interval. New chunk plans also include the right TE boundary
plane directly.

## Production Lookup Path

The optional C++ lookup path is implemented in:

```text
ThermoCalc/emissionLookup.h
ThermoCalc/emissionLookup.cpp
ThermoCalc/thermionicEmission.cpp
ThermoCalc/bindings.cpp
ThermoCalc/ThermoCalcWrapper.py
```

It is enabled only when both environment variables are set:

```powershell
$env:THERMOCALC_ENABLE_LOOKUP = "1"
$env:THERMOCALC_LOOKUP_DB = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\emission_database"
```

Use `THERMOCALC_PYD_DIR` to test the new extension without replacing the root
production `.pyd`:

```powershell
$env:THERMOCALC_PYD_DIR = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\build_cp312\Release"
```

`thermionicEmission::calc()` checks the lookup table first. A safe table hit
sets `J/Vd/delta_V/phiE/phiC` and returns immediately. A table miss continues
through the original analytic calculation.

`ThermoCalcWrapper.load_emission_lookup_database()` supports both the legacy
full database (`manifest.json` + `chunk_plan.json`) and the compact runtime
database (`runtime_manifest.json`). For production-style testing, point
`THERMOCALC_LOOKUP_DB` at `ThermoCalc/emission_runtime_db` rather than the full
analytic scan directory.

The default loaded scenario is `core` to reduce memory and startup cost. Set
`THERMOCALC_LOOKUP_REGIONS` for broader coverage, for example:

```powershell
$env:THERMOCALC_LOOKUP_REGIONS = "core,startup,high_power,accident"
```

The C++ lookup store uses bounding-box filtering, region indexes, direct TE
chunk location, and a last-block cache before running the four-dimensional
interpolator.

## Validation Snapshot

Validated commands and results:

```text
testModule/test_thermocalc_lookup.py
  passed
  lookup batch: about 3.72e6 points/s
  analytic local solver: about 9.44e4 points/s
  local speedup: about 39x
  runtime core continuous random sample: 200000 / 200000 found
```

V13 restart smoke with lookup enabled:

```text
1 s smoke completed
tec_coupled_enabled: True
no "disabling TEC coupling" message
coolant enthalpy and core delta-T differed from analytic baseline by about 1e-5 relative
TEC terminal electric power was about 62.2 W higher than analytic baseline, about 1.14%
```

V13 long-running validation:

```text
lookup run advanced from 21000 s to 22000 s in two segments
final tec_coupled_enabled: True
no runtime Traceback or disabling-TEC message was detected
```

30 s V13 timing after lookup:

```text
total wall time: 422.02 s
TEC calculate: 51.24 s, 12.16%
solid heat conduction: 252.70 s, 59.98%
fluid step_Picard: 3.41 s, 0.81%
other system overhead: 27.04%
warm-start TEC calculate: about 1.8 s per update
setup first TEC calculate after wire-resistance rebuild: about 7.81 s
```

The current performance bottleneck for V13 is no longer local thermionic
emission. It is dominated by solid heat conduction, especially the many
radiator tube wall solids in the pipe-fin radiator model.
