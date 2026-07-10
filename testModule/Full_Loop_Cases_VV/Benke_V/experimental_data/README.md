# Benke digitized experimental data

This directory stores digitized Benke thermal-validation data. Keep raw digitized values separate from model outputs.

Expected files:

## benke_sleeve_thermocouple_12pt_digitized.csv

Required columns:

```csv
thermocouple_index,sleeve_outer_k,source,notes
1,NaN,Benke figure/page TBD,replace NaN after digitization
```

- `thermocouple_index` must be 1..12.
- `sleeve_outer_k` is the digitized sleeve thermocouple temperature in K.
- Keep source and notes for traceability.

## benke_water_balance_digitized.csv

Accepted columns:

```csv
water_outlet_k,water_delta_t_k,source,notes
NaN,NaN,Benke figure/page TBD,replace one or both NaN values after digitization
```

At least one of `water_outlet_k` or `water_delta_t_k` must be populated before this file is used for comparison.

Current status: no digitized Benke sleeve/water data is available in the local figure package. Validation currently performs literature range checks only and reports `partial_missing_digitized_data`.

Status interpretation:

- `partial_missing_digitized_data`: neither sleeve nor water digitized data is available.
- `quantitative_partial_with_digitized_data`: only sleeve or water digitized data is available.
- `complete_with_digitized_data`: both sleeve and water digitized data are available and compared.

Minimum data quality requirements:

1. The sleeve CSV should contain measurements for the same Benke thermal condition being modeled, preferably the typical approximately 3003 W active-zone, 10 torr regulated-He case.
2. If the source reports Celsius, convert to K before entry and note the conversion in `notes`.
3. If the source is digitized from a plot, record the figure/page, axis calibration, and digitization tool in `source` or `notes`.
4. Water data must correspond to the same condition as the sleeve data, or the mismatch must be explicitly documented before claiming complete validation.
5. Do not copy values from `runs/*/results/*.csv`; those are model outputs, not experimental data.

After real data is added, run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_VV\Benke_V\run_benke_thermal_validation.py --run-id benke_real_digitized_validation
```

The run can be considered a complete Benke thermal-hydraulic validation only if `run_summary.json` reports:

```json
"status": "complete_with_digitized_data"
```
