from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable

from venable_single_tfe_model import (
    build_all_case_models,
    model_config_summary,
    model_summary_row,
)


CSV_FIELDS = (
    "case_id",
    "q_az_w",
    "p_out_exp_w",
    "eta_exp_percent",
    "pcs_torr",
    "tcs_k",
    "tcs_mode",
    "n_nodes",
    "active_length_m",
    "gap_m",
    "temitter_min_k",
    "temitter_max_k",
    "tcollector_min_k",
    "tcollector_max_k",
    "runs_thermocalc_calculation",
)


def _csv_value(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.6f}"
    return value


def _write_csv(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row[field]) for field in CSV_FIELDS})


def write_model_setup_outputs(output_dir: Path) -> Dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    models = build_all_case_models()
    case_rows = [model_summary_row(model) for model in models]

    cases_csv = output_dir / "cases_table71.csv"
    summary_json = output_dir / "model_config_summary.json"

    _write_csv(cases_csv, case_rows)
    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(model_config_summary(models), handle, indent=2, sort_keys=True)
        handle.write("\n")

    return {
        "cases_csv": cases_csv,
        "summary_json": summary_json,
    }


def main() -> int:
    default_output = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Build Venable Table 7-1 single-TFE model setup files only."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output,
        help="Directory for cases_table71.csv and model_config_summary.json.",
    )
    args = parser.parse_args()

    result = write_model_setup_outputs(args.output_dir)
    print(f"Wrote {result['cases_csv']}")
    print(f"Wrote {result['summary_json']}")
    print("No ThermoCalc calculation was run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
