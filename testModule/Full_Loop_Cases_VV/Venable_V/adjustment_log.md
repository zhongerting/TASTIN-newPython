# Venable_V Adjustment Log

## 2026-07-02 - baseline and thermal-closure adjustment round

### Context

The initial 14-case baseline used the placeholder thermal closure in `venable_single_tfe_model.py`:

- emitter mean: 1450 K to 1900 K
- collector mean: 760 K to 870 K
- axial shape amplitude: 0.04
- voltage scan: 0.02 V to 0.55 V, 15 points

The baseline significantly overpredicted Table 7-1 maximum output power. This indicated that the placeholder emitter/collector temperature closure was too hot for this validation target.

### Adjustment 1

Parameter: thermal closure, emitter/collector mean temperatures
Previous value: emitter 1450-1900 K, collector 760-870 K
New value: emitter 1300-1650 K, collector 700-780 K
Scope: global
Source status: fitted sensitivity, diagnostic only
Rationale: reduce systematic overprediction from the default placeholder thermal closure.
Affected files: runner CLI inputs only; no subject Python/C++ code changed.
Baseline run: `20260702_000400_baseline_14_cases_v002_055_n15`
New run: `20260702_001200_baseline_14_cases_cooler_1300_1650`
Baseline error summary: mean abs relative power error 2.7255, max abs relative power error 4.7356
New error summary: mean abs relative power error 0.2941, max abs relative power error 0.5160
Decision: keep as useful diagnostic; not final because low/mid power points remain underpredicted.
Next action: adjust global slope of emitter temperature closure.

### Adjustment 2

Parameter: thermal closure, emitter mean range
Previous value: emitter 1300-1650 K, collector 700-780 K
New value: emitter 1360-1625 K, collector 700-780 K
Scope: global
Source status: fitted sensitivity, diagnostic only
Rationale: raise low-power emitter temperature while lowering high-power emitter temperature.
Affected files: runner CLI inputs only; no subject Python/C++ code changed.
Baseline run: `20260702_001200_baseline_14_cases_cooler_1300_1650`
New run: `20260702_001700_baseline_14_cases_closure_1360_1625`
Baseline error summary: mean abs relative power error 0.2941, max abs relative power error 0.5160
New error summary: mean abs relative power error 0.2151, max abs relative power error 0.6482
Decision: reject as final candidate because the lowest-power point becomes overpredicted by 64.8%.
Next action: try an intermediate closure.

### Adjustment 3

Parameter: thermal closure, emitter mean range
Previous value: emitter 1300-1650 K / 1360-1625 K sensitivity bounds
New value: emitter 1320-1635 K, collector 700-780 K
Scope: global
Source status: fitted sensitivity, diagnostic only
Rationale: reduce both the lowest-point overprediction and high-power underprediction.
Affected files: runner CLI inputs only; no subject Python/C++ code changed.
New run: `20260702_002200_baseline_14_cases_closure_1320_1635`
New error summary: mean abs relative power error 0.2480, max abs relative power error 0.4604
Decision: reject as best candidate because mean error is worse than the next intermediate run.
Next action: try a midpoint between adjustment 2 and adjustment 3.

### Adjustment 4

Parameter: thermal closure, emitter mean range
Previous value: emitter 1320-1635 K / 1360-1625 K sensitivity bounds
New value: emitter 1340-1630 K, collector 700-780 K
Scope: global
Source status: fitted sensitivity, diagnostic only
Rationale: minimize mean and max error with one global linear thermal closure.
Affected files: runner CLI inputs only; no subject Python/C++ code changed.
New run: `20260702_002700_baseline_14_cases_closure_1340_1630`
New error summary: mean abs relative power error 0.2090, max abs relative power error 0.3812
Decision: best current diagnostic candidate, but not a final validated model. Several low/mid-power points remain underpredicted and the closure is fitted rather than source-derived.
Next action: replace fitted placeholder thermal closure with a Benke/Venable thermal-resistance or energy-balance model if source details are available.
### Adjustment 5

Parameter: Cs pressure to Cs reservoir temperature mapping (`Pcs -> Tcs`)
Previous value: placeholder table, 0.4/0.5/0.8/1.0 torr -> 560/570/590/600 K
New value: production ThermoCalc pressure formula, `Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)`
Scope: global
Source status: directly sourced from `ThermoCalc/EMISSION_SCAN_GUIDE.md` and `ThermoCalc/tools/scan_emission_map.py`
Rationale: replace an explicit placeholder with the same Cs vapor-pressure closure used by the current ThermoCalc emission lookup workflow.
Affected files: `venable_single_tfe_model.py`, `venable_validation_runner.py`, Venable_V tests and output documentation only; no subject Python/C++ code changed.
Baseline run: `20260702_002700_baseline_14_cases_closure_1340_1630`
New run: `20260702_003200_baseline_14_cases_closure_1340_1630_tcs_formula`
Baseline error summary: mean abs relative power error 0.2090, max abs relative power error 0.3812
New error summary: mean abs relative power error 0.1772, max abs relative power error 0.3125
Decision: keep as the current diagnostic input mapping because it is source-derived and improves the overall error. It is still not final V&V success because several points remain above 20% relative error.
Next action: reconstruct the emitter/collector thermal closure from Venable/Benke test-stand heat-transfer evidence rather than continue unconstrained fitting.
### Adjustment 6

Parameter: voltage scan resolution
Previous value: 15 linear voltage points from 0.02 V to 0.55 V
New value: 41 linear voltage points from 0.02 V to 0.55 V
Scope: global numerical setting
Source status: numerical convergence check
Rationale: verify whether the remaining error was caused by coarse maximum-power extraction.
Affected files: run outputs and logs only; no subject Python/C++ code changed.
Baseline run: `20260702_003200_baseline_14_cases_closure_1340_1630_tcs_formula`
New run: `20260702_004500_baseline_14_cases_closure_1340_1630_tcs_formula_scan41`
Baseline error summary: mean abs relative power error 0.1772, max abs relative power error 0.3125
New error summary: mean abs relative power error 0.1770, max abs relative power error 0.3118
Decision: scan resolution is not the dominant error source. Keep 41 points for final diagnostic confirmation runs, but do not treat scan refinement as a physical improvement.
Next action: inspect source-derived geometry and thermal closure.

### Adjustment 7

Parameter: TFE geometry area inputs
Previous value: inherited ThermoCalc-style defaults, emitter side area 0.0232138560 m2, collector side area 0.0243982365 m2, collector cross area 1.0786e-4 m2
New value: source-derived Benke/Venable single-cell dimensions: active length 0.375 m, emitter OD 19.6 mm, emitter thickness 1.15 mm, collector ID 20.6 mm, collector thickness 1.4 mm, Cs gap 0.5 mm
Scope: global source-derived geometry
Source status: derived from `TOPAZII_VV_public_experimental_data.md` single-cell structure table
Rationale: replace inherited geometry defaults with dimensions traceable to the reference test stand.
Affected files: `venable_single_tfe_model.py`, `test_venable_model_setup.py`, `cases_table71.csv`, `model_config_summary.json`; no subject Python/C++ code changed.
Baseline run: `20260702_004500_baseline_14_cases_closure_1340_1630_tcs_formula_scan41`
New run: `20260702_005500_baseline_14_cases_source_geometry_tcs_formula_scan41`
Baseline error summary: mean abs relative power error 0.1770, max abs relative power error 0.3118
New error summary: mean abs relative power error 0.1763, max abs relative power error 0.3171
Decision: keep because it is source-derived, but it does not materially solve the validation error.
Next action: evaluate whether the linear emitter-temperature closure is the dominant limitation.

### Adjustment 8

Parameter: global nonlinear emitter-temperature closure
Previous value: linear emitter mean 1340-1630 K, no quadratic term
New value: emitter mean 1325-1570 K plus 130 K quadratic peak, collector 700-780 K
Scope: global diagnostic sensitivity
Source status: fitted sensitivity, diagnostic only
Rationale: test whether raising the mid-power emitter temperature while lowering both endpoints can fix the underprediction/overprediction pattern.
Affected files: `venable_single_tfe_model.py`, `venable_validation_runner.py`, tests and run outputs; no subject Python/C++ code changed.
Baseline run: `20260702_005500_baseline_14_cases_source_geometry_tcs_formula_scan41`
New run: `20260702_010500_baseline_14_cases_source_geometry_tcs_formula_quadT_diag`
Baseline error summary: mean abs relative power error 0.1763, max abs relative power error 0.3171
New error summary: mean abs relative power error 0.3478, max abs relative power error 0.6709
Decision: reject. The peak is too large and overpredicts mid/high power cases.
Next action: try bounded smaller peak terms only.

### Adjustment 9

Parameter: global nonlinear emitter-temperature closure
Previous value: linear source-geometry closure and rejected 130 K peak sensitivity
New value: emitter mean 1325-1580 K plus 50 K quadratic peak, collector 700-780 K
Scope: global diagnostic sensitivity
Source status: fitted sensitivity, diagnostic only
Rationale: bounded smaller nonlinear correction to test the linear-closure limitation without using case-wise parameters.
Affected files: runner CLI inputs only after implementation; no subject Python/C++ code changed.
New run: `20260702_012000_baseline_14_cases_source_geometry_tcs_formula_quadT_mild_scan41`
New error summary: mean abs relative power error 0.0794, max abs relative power error 0.2030
Decision: useful diagnostic; substantially better than the linear closure but still has a 20% error point.
Next action: test one intermediate global peak value.

### Adjustment 10

Parameter: global nonlinear emitter-temperature closure
Previous value: emitter mean 1325-1580 K plus 50 K quadratic peak
New value: emitter mean 1325-1580 K plus 60 K quadratic peak, collector 700-780 K
Scope: global diagnostic sensitivity
Source status: fitted sensitivity, diagnostic only
Rationale: reduce the remaining 1405-1580 W underprediction while checking that high-power cases are not overfit.
Affected files: runner CLI inputs only; no subject Python/C++ code changed.
Baseline run: `20260702_012000_baseline_14_cases_source_geometry_tcs_formula_quadT_mild_scan41`
New run: `20260702_014500_baseline_14_cases_source_geometry_tcs_formula_quadT_peak60_scan41`
Baseline error summary: mean abs relative power error 0.0794, max abs relative power error 0.2030
New error summary: mean abs relative power error 0.0727, max abs relative power error 0.1377
Decision: current best diagnostic candidate. It is still not a fully source-derived V&V result because the nonlinear thermal closure is fitted; do not continue to case-wise fitting.
Next action: stage summary and identify missing experimental thermal-boundary information.
### Adjustment 11

Parameter: axial heat-release / temperature profile
Previous value: cosine profile over the full 375 mm active length
New value: centered TISA profile with 300 mm heated length in the 375 mm thermionic working section
Scope: global source-informed axial shape
Source status: directly sourced geometry feature from `TOPAZII_VV_public_experimental_data.md`; temperature amplitude remains diagnostic
Rationale: Venable/Benke test stand uses a 300 mm tungsten TISA heater inside a 375 mm thermionic working section, so the axial temperature shape should not assume full-length uniform heating.
Affected files: `venable_single_tfe_model.py`, `venable_validation_runner.py`, tests and output documentation; no subject Python/C++ code changed.
Baseline run: `20260702_014500_baseline_14_cases_source_geometry_tcs_formula_quadT_peak60_scan41`
New run: `20260702_033000_tisa300_profile_meanlift_scan41`
Baseline error summary: mean abs relative power error 0.0727, max abs relative power error 0.1377
New error summary: mean abs relative power error 0.0685, max abs relative power error 0.1538
Decision: keep as a source-informed axial profile. It improves mean error but increases the maximum relative error, so it remains diagnostic rather than final source-derived validation.
Next action: tune the global emitter peak under the TISA profile.

### Adjustment 12

Parameter: global emitter peak under TISA axial profile
Previous value: emitter mean 1360-1620 K plus 62 K quadratic peak
New value: emitter mean 1360-1620 K plus 55 K quadratic peak
Scope: global diagnostic sensitivity
Source status: fitted sensitivity, diagnostic only
Rationale: reduce high-power overprediction introduced by the TISA profile while preserving the lower mean error.
Affected files: runner CLI inputs only; no subject Python/C++ code changed.
Baseline run: `20260702_033000_tisa300_profile_meanlift_scan41`
New run: `20260702_035000_tisa300_profile_peak55_scan41`
Baseline error summary: mean abs relative power error 0.0685, max abs relative power error 0.1538
New error summary: mean abs relative power error 0.0664, max abs relative power error 0.1566
Decision: current lowest mean-error diagnostic candidate. It sacrifices maximum relative error compared with the previous stage but improves the overall 14-point fit.
Next action: evaluate Benke water-jacket and He-gap closures.

### Adjustment 13

Parameter: Benke water-jacket collector boundary
Previous value: linear collector mean 710-790 K
New value: `collector_boundary_mode=benke_water_jacket`, using water inlet temperature, water mass flow, water-side h, regulated He gap width/conductivity, and heat-pickup fraction
Scope: global source-structured diagnostic model
Source status: structure directly sourced from Benke summary; numerical boundary values remain diagnostic because Venable Table 7-1 does not provide per-run cooling-water or sleeve-temperature data
Rationale: test whether the collector temperature can be generated from the Benke water-jacket/He-gap heat-resistance path instead of a fitted linear collector temperature.
Affected files: `venable_single_tfe_model.py`, `venable_validation_runner.py`, tests and output documentation; no subject Python/C++ code changed.
Representative runs: `20260702_022000_tisa300_benke_water_default_diag`, `20260702_024000_tisa300_benke_water_equiv650_diag`, `20260702_025000_tisa300_benke_water_frac010_tin650_diag`, `20260702_030000_tisa300_benke_water_frac015_tin650_diag`, `20260702_031000_tisa300_benke_water_frac013_em1365_1610_diag`, `20260702_032000_tisa300_benke_water_frac013_em1365_1630_peak50_diag`
Best water/He diagnostic summary: mean abs relative power error 0.0756 for `20260702_031000_tisa300_benke_water_frac013_em1365_1610_diag`; lowest maximum relative error among tested water/He runs was about 0.156.
Decision: do not keep as current best. The simplified Benke water-jacket closure is physically structured but underconstrained; without measured sleeve/water boundary data it did not beat the empirical collector closure.
Next action: keep the implemented water/He parameters for future source-derived thermal-boundary work, but use the TISA-profile empirical collector closure as the current best electrical-validation candidate.

## Current best diagnostic candidate

Run ID: `20260702_035000_tisa300_profile_peak55_scan41`

Thermal, Cs, and geometry closure:

- source-derived single-cell geometry from Benke/Venable dimensions
- axial profile: centered 300 mm TISA heated length inside 375 mm active length
- emitter mean min: 1360 K
- emitter mean max: 1620 K
- emitter quadratic peak: 55 K
- collector mean min: 710 K
- collector mean max: 790 K
- collector boundary mode: `linear` retained as best empirical electrical-validation candidate
- axial shape amplitude: 0.04
- `Tcs` mode: `pressure_formula`
- voltage scan: 0.02-0.55 V, 41 points

Error summary:

- case count: 14
- mean absolute relative power error: 6.64%
- maximum absolute relative power error: 15.66%
- largest remaining relative error: 1405 W case, -15.66%
- largest remaining absolute error: 2474 W case, +14.11 W

Status: best current mean-error diagnostic candidate. It is more source-informed axially than the previous best, but it is still not a fully source-derived V&V result because emitter/collector mean temperatures remain fitted.
