# Benke user-provided experimental data assessment

## Data admission result

The newly provided files under `experimental_data/` satisfy the minimum admission requirements for Benke thermal-hydraulic validation:

- `benke_sleeve_thermocouple_12pt_digitized.csv`
  - Contains thermocouple indices 1..12.
  - Contains 11 finite sleeve temperature measurements in K.
  - Index 9 / T64 is explicitly marked `NaN` because Benke reported it inoperative.
  - Source and notes identify Benke 1994 Appendix B, p.79, TISA input power 3412 W row, DTIC ADA293595 PDF.
- `benke_water_balance_digitized.csv`
  - Contains finite `water_outlet_k` and `water_delta_t_k`.
  - Notes identify inlet temperature, two outlet thermocouples, flow rate, and the same Appendix B row.

These files are therefore acceptable as traceable experimental data. They are not model outputs.

## Required validation handling

Because T64 is inoperative, validation must ignore non-finite sleeve values instead of including them in MAE/RMSE. `benke_validation.py` has been updated accordingly:

- `point_count = 11`
- `expected_point_count = 12`
- `ignored_indices = [9]`

## Matched-boundary validation run

Run directory:

```text
runs/20260702_benke_user_data_matched_boundary/
```

Boundary inputs used to match the provided Benke Appendix B row:

| Input | Value |
| --- | ---:|
| TISA input power | 3412 W |
| active-zone power | 3002.56 W |
| regulated He pressure | 10 torr |
| water inlet temperature | 289.71 K |
| water mass flow | 0.043518 kg/s |
| water-side h | 800 W/(m2 K) |
| regulated He effective k | 0.08 W/(m K) |

Validation status:

```text
complete_with_digitized_data
```

Range checks:

```text
passed
```

## Quantitative comparison

| Metric | Value |
| --- | ---:|
| Sleeve point count used | 11 |
| Ignored sleeve index | 9 |
| Sleeve MAE | 232.379 K |
| Sleeve RMSE | 249.852 K |
| Sleeve max abs error | 365.553 K |
| Sleeve mean error | 205.383 K |
| Water outlet abs error | 0.946 K |
| Water delta-T abs error | 0.946 K |
| Energy balance error | -5.32e-11 W |

## Interpretation

The user-provided data is sufficient to run a complete digitized-data validation. The current v1 thermal-network model does not yet match the sleeve thermocouple temperatures: it overpredicts sleeve temperatures by about 205 K on average for the matched-boundary run.

The water-side heat balance is much closer once the Benke row inlet temperature and flow rate are used. The remaining water delta-T difference is about 0.95 K with the current constant water heat capacity and active-zone power assumption.

The next model-improvement step should focus on the radial/axial heat path before using Benke boundaries in Venable electrical validation. Candidate adjustment targets include:

- regulated He effective conductivity within or near the Benke inferred range;
- water-side h within the Benke inferred range;
- active axial heat-source distribution / thermocouple axial mapping;
- extra parallel or bypass heat paths not represented in the current v1 series-only radial network;
- exact sleeve thermocouple physical positions and whether the reported T56-T67 ordering matches model axial index ordering.
