# Venable_V Validation Result Summary

## Status

Current status: Venable Table 7-1 validation calculation has been extended with a source-informed axial TISA heating profile and a Benke-style water-jacket/He-gap collector-boundary option. The best current electrical-validation candidate improves the mean 14-point error, but it is still a diagnostic fit rather than a fully source-derived predictive V&V model.

No subject Python/C++ source files were modified. Changes were confined to `Venable_V` case scripts, tests, outputs, and documentation.

## Current Best Run

Run ID: `20260702_035000_tisa300_profile_peak55_scan41`

Input summary:

| Item | Value |
| --- | --- |
| Geometry | active length 0.375 m, emitter OD 19.6 mm, emitter wall 1.15 mm, collector ID 20.6 mm, collector wall 1.4 mm, Cs gap 0.5 mm |
| Axial profile | `tisa_300mm`, centered 300 mm heated length inside 375 mm thermionic working section |
| `Pcs -> Tcs` | `Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)` |
| Emitter mean closure | 1360-1620 K linear with +55 K global quadratic peak |
| Collector mean closure | 710-790 K linear empirical electrical-validation closure |
| Voltage scan | 0.02-0.55 V, 41 points |

Error summary:

| Metric | Value |
| --- | ---: |
| Case count | 14 |
| Finite cases | 14 |
| Mean absolute relative power error | 6.64% |
| Maximum absolute relative power error | 15.66% |
| Mean signed power error | +0.44 W |
| Maximum absolute power error | 14.11 W |
| Finite scan points | all finite |
| All scan points converged and target-matched | false |

`converged_all_cases=false` is retained because some scan points report non-convergence or target-voltage mismatch, but maximum-power extraction uses only finite target-matched scan points.

## Detailed Results

| Q_az W | Pcs torr | P_exp W | P_calc W | Relative error | eta_exp % | eta_calc % |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 892 | 0.4 | 10.23 | 9.96 | -2.6% | 1.15 | 1.12 |
| 1062 | 0.4 | 17.80 | 16.83 | -5.5% | 1.68 | 1.58 |
| 1237 | 0.4 | 30.13 | 26.41 | -12.4% | 2.44 | 2.13 |
| 1405 | 0.4 | 45.00 | 37.95 | -15.7% | 3.20 | 2.70 |
| 1580 | 0.5 | 63.25 | 54.56 | -13.7% | 4.01 | 3.45 |
| 1755 | 0.5 | 77.28 | 70.95 | -8.2% | 4.40 | 4.04 |
| 1933 | 0.5 | 86.26 | 88.03 | +2.1% | 4.46 | 4.55 |
| 2112 | 0.5 | 103.97 | 103.02 | -0.9% | 4.92 | 4.88 |
| 2281 | 0.8 | 115.44 | 125.95 | +9.1% | 5.06 | 5.52 |
| 2474 | 0.8 | 129.87 | 143.98 | +10.9% | 5.25 | 5.82 |
| 2637 | 0.8 | 146.75 | 156.85 | +6.9% | 5.57 | 5.95 |
| 2813 | 1.0 | 167.06 | 168.98 | +1.1% | 5.94 | 6.01 |
| 2999 | 1.0 | 178.16 | 179.82 | +0.9% | 5.94 | 6.00 |
| 3162 | 1.0 | 192.46 | 186.57 | -3.1% | 6.09 | 5.90 |

## Run Comparison

| Run ID | Main change | Mean abs rel error | Max abs rel error | Decision |
| --- | --- | ---: | ---: | --- |
| `20260702_014500_baseline_14_cases_source_geometry_tcs_formula_quadT_peak60_scan41` | previous best, cosine axial profile | 7.27% | 13.77% | lower max error, higher mean error |
| `20260702_020000_tisa300_profile_peak60_scan41` | direct TISA 300 mm axial profile, no mean retune | 35.00% | 100.00% | rejected; axial profile lowered effective emission too much |
| `20260702_033000_tisa300_profile_meanlift_scan41` | TISA profile with global mean lift and 62 K peak | 6.85% | 15.38% | useful, superseded by 55 K peak |
| `20260702_035000_tisa300_profile_peak55_scan41` | TISA profile with 55 K peak | 6.64% | 15.66% | current best mean-error diagnostic candidate |
| `20260702_031000_tisa300_benke_water_frac013_em1365_1610_diag` | best tested Benke water/He boundary candidate | 7.56% | 15.61% | not retained; physical structure useful but underconstrained |

## Benke Water-Jacket / He-Gap Assessment

Implemented collector-boundary mode: `benke_water_jacket`.

The model includes:

- cooling-water inlet/effective boundary temperature;
- cooling-water mass flow and heat capacity;
- coolant heat-pickup fraction;
- water-side heat-transfer coefficient;
- regulated He gap width;
- regulated He gap effective conductivity;
- optional extra collector resistance.

The structure follows the Benke summary: water-jacket heat balance, water-side heat transfer, regulated He gap, and radial heat-resistance path. However, the available local Venable Table 7-1 data do not provide the per-run sleeve temperature, true water inlet/outlet temperature, water flow, or regulated He gap condition. As a result, the explicit water/He closure was underconstrained and did not outperform the empirical collector temperature closure.

Important diagnostic result: using room-temperature water-like conditions made collector temperature far too low, causing severe underprediction and even zero output in the lowest-power point. Raising the effective boundary temperature or reducing He conductivity can recover low/mid-power points but tends to overpredict high-power points unless additional source data constrain the heat path.

## Interpretation

The TISA 300 mm axial profile is a real source-informed improvement over the earlier full-length cosine assumption. It lowers the mean absolute relative error from 7.27% to 6.64%. The tradeoff is that the maximum relative error increases from 13.77% to 15.66%, mainly around 1405 W and 2474 W.

The Benke water-jacket and regulated-He-gap model is now available for future work, but the current public summary is insufficient to turn it into a fully predictive collector boundary for the 14 Venable electrical runs. It should not be forced into the final model by using an unrealistically high water inlet temperature. If future data provide sleeve thermocouple profiles, water inlet/outlet temperatures, or regulated He pressure for each electrical run, this mode should replace the empirical collector closure.

## Conclusion

This round made concrete progress toward the requested direction:

- axial heat-release distribution was changed to a source-informed 300 mm TISA profile;
- water-jacket and He-gap collector-boundary parameters were implemented and tested;
- water flow/temperature and He-gap effective resistance were explored;
- the best retained result is `20260702_035000_tisa300_profile_peak55_scan41`, with 6.64% mean absolute relative error and 15.66% maximum relative error.

The next useful modification is not more unconstrained fitting. It is to obtain or digitize Benke/Venable thermal-boundary data: sleeve axial temperature, regulated He gap pressure/state, and cooling-water inlet/outlet conditions for the electrical runs.
