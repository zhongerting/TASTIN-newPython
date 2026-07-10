# Benke_V goal completion audit

This audit maps the active goal requirements to current evidence in `Benke_V`. It is intentionally strict: the goal is not considered complete until the Benke experiment is quantitatively compared against traceable Benke data.

## Requirement status

| Requirement | Current evidence | Status |
| --- | --- | --- |
| Create `Full_Loop_Cases_VV/Benke_V` | Directory exists with model, validation, report, tests, data templates, run outputs, and process log. | Complete |
| Classify Benke thermal-closure data into direct / inferable / sensitivity-only groups | `BENKE_THERMAL_CLOSURE_DATA.md` contains the three-class table and Venable transfer rules. | Complete |
| Decide between simple lumped resistance and axial segmented water-cooled model | `BENKE_MODEL_DECISION.md` chooses axial segmented water-cooled radial resistance network and records why. | Complete |
| Store modeling data in `.md` files under `Benke_V` | `BENKE_THERMAL_CLOSURE_DATA.md`, `BENKE_MODEL_DECISION.md`, `BENKE_VALIDATION_STATUS.md`, `BENKE_SOURCE_ACQUISITION_LOG.md`, and this audit file are under `Benke_V`. | Complete |
| Build Benke thermal-hydraulic model | `benke_thermal_network.py` implements the axial segmented water-cooled thermal network with active-zone power correction, radial resistance chain, water temperature marching, and 12 sleeve thermocouple sampling. | Complete for v1 model |
| Implement validation workflow | `run_benke_thermal_validation.py`, `benke_validation.py`, `benke_report.py`, `benke_parameter_scan.py`, and `benke_calibration.py` implement range checks, CSV comparison, reporting, parameter envelope scanning, and sleeve-data calibration. | Complete as data-ready workflow |
| Quantitatively validate against Benke experimental sleeve/water data | Current local figure package and public searches have not found Benke sleeve thermocouple or water-balance digitized data. `experimental_data/` contains templates only. | Incomplete |

## Current validation level

Current status is `partial_missing_digitized_data`:

- Literature range checks can be performed.
- Water-side energy balance is internally closed.
- The code can compare against real Benke digitized CSVs when they are supplied.
- The code must not claim completed curve-level validation until both `benke_sleeve_thermocouple_12pt_digitized.csv` and `benke_water_balance_digitized.csv` are populated from traceable Benke sources.

## Evidence from latest smoke run

Latest smoke run:

```text
runs/20260702_benke_report_status_smoke/
```

Key outputs:

| Output | Value |
| --- | ---:|
| active-zone power | 3003 W |
| water outlet temperature | 333.947 K |
| water delta-T | 23.947 K |
| energy balance error | -1.53e-10 W |
| collector inner mean/max | 988.877 / 1073.399 K |
| sleeve outer mean/max | 949.196 / 1029.309 K |
| validation status | `partial_missing_digitized_data` |
| range checks | `passed` |

## Data required to close the goal

To mark the goal complete, add traceable Benke experimental data under `experimental_data/`:

1. `benke_sleeve_thermocouple_12pt_digitized.csv`
   - Required columns: `thermocouple_index,sleeve_outer_k`
   - Required content: 12 traceable sleeve thermocouple temperatures from Benke experiment/figure/table.
2. `benke_water_balance_digitized.csv`
   - Required columns: at least one of `water_outlet_k` or `water_delta_t_k`
   - Required content: traceable Benke water outlet temperature or water temperature rise for the same or explicitly matched condition.
3. Source traceability
   - Record figure/page/table identifier, digitization method, and any unit conversions in the CSV `source`/`notes` columns or an adjacent note.

After these files are added, run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_VV\Benke_V\run_benke_thermal_validation.py --run-id <new_run_id>
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest discover -s testModule\Full_Loop_Cases_VV\Benke_V -p "test_*.py"
```

Completion evidence will be:

- `run_summary.json` reports `validation.status = complete_with_digitized_data`.
- `validation_report.md` includes sleeve MAE/RMSE and water-side error metrics.
- Local Benke_V tests pass.

## Non-admission rule

Do not populate `experimental_data/*_digitized.csv` from:

- model-generated output;
- Venable electrical output or maximum-power data;
- inferred values chosen to make the model pass;
- any source that cannot be traced to Benke or a clearly matched Benke/Venable single-cell test-stand thermal measurement.
