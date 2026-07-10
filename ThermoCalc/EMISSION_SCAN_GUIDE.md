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

## End-to-End Lookup Workflow

The thermionic lookup workflow has two separate phases: an offline data
generation phase and a runtime calculation phase. The original analytic solver
remains the source of truth in the offline phase and remains the fallback in the
runtime phase.

```text
offline:
  thermionicEmission::calcDiagnostics()
    -> emission_database.py plan / worker
    -> ThermoCalc/emission_database/
    -> summarize / verify / optimize-table
    -> emission_database.py export-runtime or export-runtime-dense
    -> ThermoCalc/emission_runtime_db/ or emission_runtime_db_v2/

runtime:
  ThermoCalcWrapper.py
    -> load_emission_lookup_database()
    -> te_solver.add_emission_runtime_block() or load_emission_dense_file()
    -> emissionLookup.cpp in-memory store
    -> thermionicEmission::calc()
    -> queryEmissionLookup()
    -> hit: return lookup J/Vd/delta_V/phiE/phiC
    -> miss: execute original analytic calc()
```

### 1. Define The Scan Grid

`emission_database.py plan` builds the scan axes and chunk plan. The axes are
stored under `ThermoCalc/emission_database/axes/`, and the chunk layout is
stored in `chunk_plan.json`.

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py plan --db-dir ThermoCalc\emission_database --preset full
```

The scan variables are:

```text
TE   emitter temperature
TC   collector temperature
Vo   output/electrode voltage parameter
Pcs  cesium pressure
Tcs  cesium temperature converted from Pcs
```

The table is split into regions: `core`, `startup`, `high_power`, and
`accident`. Region priority is used at runtime when regions overlap; lower
priority numbers are preferred.

### 2. Generate Raw Diagnostic Chunks

`emission_database.py worker` reads `chunk_plan.json`, calculates assigned
chunks, and writes raw `.npz` files under `ThermoCalc/emission_database/chunks/`.
Each point calls the diagnostic C++ path:

```text
te_solver.calc_emission_point()
  -> thermionicEmission::calcDiagnostics()
```

The raw chunk files contain both runtime fields and diagnostic fields:

```text
runtime fields:
  J, Vd, delta_V, phiE, phiC

diagnostic fields:
  converged, finite_flag, done
  regime, iteration_count, error_code
  valid_for_interpolation, near_failed_region
  zero_emission_flag, high_risk_flag, source_region_id
```

Typical parallel generation uses multiple workers:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py worker --db-dir ThermoCalc\emission_database --worker-id w0 --worker-index 0 --worker-count 6
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py worker --db-dir ThermoCalc\emission_database --worker-id w1 --worker-index 1 --worker-count 6
```

### 3. Summarize, Verify, And Optimize

After raw chunks are generated, `summarize` writes audit outputs under
`ThermoCalc/emission_database/summaries/`, and `verify` compares sampled table
points against the diagnostic analytic solver.

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py summarize --db-dir ThermoCalc\emission_database --scan-chunks
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py verify --db-dir ThermoCalc\emission_database --samples 200
```

`optimize-table` is used for regions where the raw analytic scan has known
non-converged points. It writes `.optimized.npz` sidecars next to raw chunks.
Those sidecars are preferred by the runtime exporter and loader.

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py optimize-table --db-dir ThermoCalc\emission_database --zero-j-threshold 1e-3 --region startup --region accident
```

The optimization policy is:

```text
safe zero-emission points:
  set J = 0 and keep voltage/work-function fields

remaining isolated invalid points:
  use neighbor imputation where available

unresolved points:
  stay unsafe and will not be used for interpolation
```

### 4. Export The Runtime Database

The raw database is audit-friendly but too large and too detailed for normal
runtime use. `export-runtime` creates the compact runtime database under
`ThermoCalc/emission_runtime_db/`.

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py export-runtime --db-dir ThermoCalc\emission_database --out-dir ThermoCalc\emission_runtime_db --dtype float32 --zero-compress
```

The runtime export keeps only:

```text
TE_axis, TC_axis, Vo_axis, Tcs_axis
J, Vd, delta_V, phiE, phiC
lookup_safe, zero_mask
```

It also stitches the first TE plane from the next source chunk onto the current
runtime chunk when needed. This is what makes legacy chunks such as
`1300-1310 K` and `1320-1330 K` become continuous runtime chunks such as
`1300-1320 K` and `1320-1340 K`.

Dense runtime v2 is the preferred compact local runtime artifact. It stores one
dense tensor per region instead of many stitched TE chunks:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py export-runtime-dense --db-dir ThermoCalc\emission_database --out-dir ThermoCalc\emission_runtime_db_v2 --format both --region core --dtype float32 --zero-compress
```

Dense runtime v2 writes:

```text
runtime_dense_manifest.json
<region>.runtime.v2.npz
<region>.runtime.v2.tedb
```

The `.npz` file is the portable Python artifact. The `.tedb` file is a simple
binary artifact that C++ can load directly through `load_emission_dense_file()`.
When both are present, `ThermoCalcWrapper.load_emission_lookup_database()`
prefers the `.tedb` path because it avoids constructing large NumPy objects in
Python.

Dense runtime v2 keeps only:

```text
TE_axis, TC_axis, Vo_axis, Tcs_axis
J, Vd, delta_V, phiE, phiC
lookup_safe_bits, zero_mask_bits
```

`lookup_safe` and `zero_mask` are bit-packed using little-endian bit order.
The dense exporter also avoids duplicate stitched TE boundary planes, so the
core dense v2 artifact is smaller than the previous stitched runtime export
while preserving the same continuous interpolation range. `phiE/phiC` are kept
because upper-level boundary handling still uses them.

### 5. Load The Runtime Database

At runtime, use the test extension and enable lookup explicitly:

```powershell
$env:THERMOCALC_PYD_DIR = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\build_cp312\Release"
$env:THERMOCALC_ENABLE_LOOKUP = "1"
$env:THERMOCALC_LOOKUP_DB = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\emission_runtime_db_v2"
$env:THERMOCALC_LOOKUP_REGIONS = "core"
```

`ThermoCalcModel.__init__()` calls `load_emission_lookup_database()` when
`THERMOCALC_ENABLE_LOOKUP=1` and `THERMOCALC_LOOKUP_DB` is set. The loader
detects `runtime_manifest.json`, reads selected region chunks, and calls:

```text
te_solver.clear_emission_lookup()
te_solver.add_emission_runtime_block(...)
te_solver.set_emission_lookup_enabled(True)
```

The selected regions default to `core`. Set
`THERMOCALC_LOOKUP_REGIONS=core,startup,high_power,accident` when the system
needs full startup, high-power, and accident coverage.

### 6. Query During TEC Calculation

The production path enters C++ through `CircuitTECs.calc()` and eventually calls
`thermionicEmission::calc()` for each local emission solve. `calc()` first asks
the in-memory lookup store:

```text
queryEmissionLookup(TE, TC, Vo, Tcs, d_gap)
```

The C++ store checks:

```text
d_gap support
last-block cache
region priority
TE chunk index
chunk bounding box
axis location
lookup_safe on all interpolation corners
finite interpolated fields
```

A safe hit directly fills:

```text
J, Vd, delta_V, phiE, phiC
```

and returns without running the analytic iteration. A miss falls through to the
original analytic `thermionicEmission::calc()` logic. This fallback is deliberate
so the lookup path remains optional and does not remove the existing solver.

### 7. Verify The Full Chain

Use the focused lookup regression after changing table generation, runtime
export, bindings, or C++ interpolation:

```powershell
$env:THERMOCALC_PYD_DIR = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\build_cp312\Release"
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\test_thermocalc_lookup.py
```

This test covers:

```text
raw chunk loading
optimized sidecar loading
runtime export/load
TE boundary stitching
single-point lookup
batch lookup output layout
production calc() lookup branch
lookup-vs-analytic speed comparison
```

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

The default full plan covers four separately stored regions. All production
regions use `Pcs = 0.02-5.0 torr` with log spacing, while preserving the
original per-region pressure point counts:

```text
core/high_power: 41 Pcs points
startup:         21 Pcs points
accident:        31 Pcs points
```

`Pcs` is in torr, not Pa. The conversion to `Tcs` uses the same formula as the
C++ production model:

```text
Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)
```

```text
core        TE 1300-2150 K, TC 700-900 K,  Vo 0-3.5 V, Pcs 0.02-5.0 torr, 41 Pcs points
startup     TE 700-1300 K,  TC 500-800 K,  Vo 0-3.5 V, Pcs 0.02-5.0 torr, 21 Pcs points
high_power  TE 2150-2400 K, TC 750-1000 K, Vo 0-3.5 V, Pcs 0.02-5.0 torr, 41 Pcs points
accident    TE 700-2400 K,  TC 500-1100 K, Vo 0-3.5 V, Pcs 0.02-5.0 torr, 31 Pcs points
```

Generate the full manifest and chunk plan:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py plan --db-dir ThermoCalc\emission_database --preset full
```

The corrected full plan contains 18,737,388 points in 76 chunks with the
current right-boundary-preserving chunk layout. Use multiple
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

As of 2026-06-23, the corrected `0.02-5.0 torr` local database has been
generated under:

```text
ThermoCalc/emission_database/pcs_0p02_5torr/
```

Generation status:

```text
planned unique grid points: 18,737,388
planned chunks: 76
raw chunk files: 76 / 76
raw summary files: 76 / 76
missing chunks: 0
```

The raw `summarize --scan-chunks` result counts duplicated TE boundary planes
inside chunk files, so its `done_points` is larger than the unique grid count:

```text
raw chunk-counted done points: 25,756,400
raw chunk-counted converged points: 24,857,196
raw chunk-counted failed or non-finite points: 899,204
finite rate over done points: 1.0
```

The unique physical grid count remains `18,737,388`, from `manifest.json`.

The previously generated full local database under `ThermoCalc/emission_database/`
used the older pressure ranges `0.5-2.0 torr` for `core/startup/high_power` and
`0.1-4.0 torr` for `accident`. It is now an obsolete pressure-range artifact.

The old generated database summary was:

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

The analytic scan can contain non-converged points. Direct analytic fallback is
not useful there because the fallback often fails for the same reason. The table
therefore supports optimized sidecars:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py optimize-table --db-dir ThermoCalc\emission_database\pcs_0p02_5torr --zero-j-threshold 1e-3 --region core --region startup --region high_power --region accident
```

Current `pcs_0p02_5torr` optimized result:

```text
regions optimized: core, startup, high_power, accident
total points: 18,737,388
raw invalid points: 655,530
zero-filled points: 107,272
neighbor-imputed points: 548,258
unresolved points: 0
optimized safe rate: 1.0
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

Dense runtime v2 is now the recommended local runtime format:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" ThermoCalc\tools\emission_database.py export-runtime-dense --db-dir ThermoCalc\emission_database\pcs_0p02_5torr --out-dir ThermoCalc\emission_runtime_db_v2\pcs_0p02_5torr --format both --dtype float32 --zero-compress
```

It writes `runtime_dense_manifest.json`, one `*.runtime.v2.npz` per region,
and optionally one `*.runtime.v2.tedb` per region. The `.npz` file is portable;
the `.tedb` file is loaded directly by C++ and is preferred by the wrapper when
available. Dense v2 stores `J/Vd/delta_V/phiE/phiC` as `float32` by default and
stores `lookup_safe` / `zero_mask` as bit-packed masks:

```text
lookup_safe_bits
zero_mask_bits
```

This removes boolean-array overhead and avoids duplicate stitched TE planes.
`phiE/phiC` remain stored because they are still needed by upper-level boundary
handling.

The current corrected dense runtime v2 output is:

```text
ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/
  runtime_dense_manifest.json
  core.runtime.v2.npz / core.runtime.v2.tedb
  startup.runtime.v2.npz / startup.runtime.v2.tedb
  high_power.runtime.v2.npz / high_power.runtime.v2.tedb
  accident.runtime.v2.npz / accident.runtime.v2.tedb
```

Current output summary:

```text
region_count: 4
total_points: 18,737,388
total_size_bytes: 537,914,570
zero_compress: true
zero_j_threshold: 1e-3

core        shape 86 x 41 x 71 x 41, points 10,264,186, NPZ 87,074,650 bytes, TEDB 207,851,764 bytes
startup     shape 31 x 31 x 36 x 21, points    726,516, NPZ  5,368,120 bytes, TEDB  14,712,989 bytes
high_power  shape 25 x 26 x 71 x 41, points  1,892,150, NPZ 18,317,488 bytes, TEDB  38,317,432 bytes
accident    shape 86 x 61 x 36 x 31, points  5,854,536, NPZ 47,715,973 bytes, TEDB 118,556,154 bytes
```

`ThermoCalc/emission_runtime_db_v2/` is generated runtime data and is ignored by
git.

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
$env:THERMOCALC_LOOKUP_DB = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\emission_runtime_db_v2\pcs_0p02_5torr"
```

Use `THERMOCALC_PYD_DIR` to test the new extension without replacing the root
production `.pyd`:

```powershell
$env:THERMOCALC_PYD_DIR = "E:\项目任务\五院-电源\source_code\TASTIN-python\ThermoCalc\build_cp312\Release"
```

`thermionicEmission::calc()` checks the lookup table first. A safe table hit
sets `J/Vd/delta_V/phiE/phiC` and returns immediately. A table miss continues
through the original analytic calculation.

`ThermoCalcWrapper.load_emission_lookup_database()` supports the legacy full
database (`manifest.json` + `chunk_plan.json`), the legacy compact runtime
database (`runtime_manifest.json`), and dense runtime v2
(`runtime_dense_manifest.json`). For production-style testing, point
`THERMOCALC_LOOKUP_DB` at `ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr`,
`ThermoCalc/emission_runtime_db_v2`, or
`ThermoCalc/emission_runtime_db` rather than the full analytic scan directory.

The default loaded scenario is `core` to reduce memory and startup cost. Set
`THERMOCALC_LOOKUP_REGIONS` for broader coverage, for example:

```powershell
$env:THERMOCALC_LOOKUP_REGIONS = "core,startup,high_power,accident"
```

The C++ lookup store uses bounding-box filtering, region indexes, direct TE
chunk location for legacy runtime chunks, dense-region bounding boxes, and a
last-hit cache before running the four-dimensional interpolator.

## Validation Snapshot

Current corrected database validation:

```text
raw database: ThermoCalc/emission_database/pcs_0p02_5torr/
runtime database: ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr/
planned unique grid points: 18,737,388
planned chunks: 76
raw chunk files: 76 / 76
raw summary files: 76 / 76
optimized safe rate: 1.0
unresolved optimized points: 0
dense runtime v2 total size: 537,914,570 bytes
```

Dense runtime v2 loading smoke:

```text
load_emission_lookup_database(..., regions=["core","startup","high_power","accident"], force=True)
  loaded regions: 4
  dense region count: 4
  legacy block count: 0
  core sample lookup: found=True, source=core
  startup sample lookup: found=True, source=startup
```

Earlier focused regressions measured dense lookup at about `3.55e6 points/s`
inside `testModule/test_thermocalc_lookup.py`, versus analytic local solver
speed of about `9.87e4 points/s`. Re-run that benchmark after replacing or
selecting the desired `.pyd` if exact machine-local speed is needed.

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
