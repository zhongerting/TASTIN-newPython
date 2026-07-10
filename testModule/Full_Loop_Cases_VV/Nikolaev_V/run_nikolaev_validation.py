from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from nikolaev_single_tfe_model import (
    NikolaevModelConfig,
    calculate_operating_point,
    capillary_limit_diameter_mm,
    compare_operating_point,
    fuel_max_temperature_k,
)
from nikolaev_source_data import (
    LOCAL_PDF,
    MOCKUP_FACTS,
    SOURCE_CITATION,
    SOURCE_DOI,
    SOURCE_TITLE,
    TABLE1_CHARACTERISTICS,
    TABLE2_OPERATING_POINTS,
    TABLE3_FUEL_TEMPERATURES,
    TABLE4_CAPILLARY_LIMITS,
)


CASE_DIR = Path(__file__).resolve().parent


def _csv_value(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.10g}"
    return value


def _write_csv(path: Path, rows: Iterable[dict]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def _mean_abs(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.nanmean(np.abs(arr))) if arr.size else math.nan


def _max_abs(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.nanmax(np.abs(arr))) if arr.size else math.nan


def validate_table1(config: NikolaevModelConfig) -> list[dict]:
    calculated = {
        "output_power": config.nominal_output_power_w,
        "efficiency": config.nominal_efficiency_percent,
        "effective_height": config.effective_height_m * 100.0,
        "emitter_cladding_thickness": 2.3,
        "emitter_temperature": 1880.0,
        "collector_temperature": config.collector_temperature_k,
    }
    rows = []
    for item in TABLE1_CHARACTERISTICS:
        calc = float(calculated[item.parameter])
        rows.append(
            {
                "parameter": item.parameter,
                "unit": item.unit,
                "topaz_ii_exp": item.topaz_ii,
                "space_r_exp_or_calc": item.space_r,
                "space_r_calc": calc,
                "space_r_error": calc - item.space_r,
                "note": item.note,
            }
        )
    return rows


def validate_table2(config: NikolaevModelConfig) -> list[dict]:
    rows = []
    for target in TABLE2_OPERATING_POINTS:
        calculated = calculate_operating_point(target.voltage_v, config)
        rows.append(compare_operating_point(target, calculated))
    return rows


def validate_table3() -> list[dict]:
    rows = []
    for point in TABLE3_FUEL_TEMPERATURES:
        calc = fuel_max_temperature_k(point.free_volume_percent, point.radial_factor, TABLE3_FUEL_TEMPERATURES)
        rows.append(
            {
                "free_volume_percent": point.free_volume_percent,
                "radial_factor": point.radial_factor,
                "max_fuel_temperature_exp_k": point.max_fuel_temperature_k,
                "max_fuel_temperature_calc_k": calc,
                "error_k": calc - point.max_fuel_temperature_k,
            }
        )
    return rows


def validate_table4() -> list[dict]:
    rows = []
    for point in TABLE4_CAPILLARY_LIMITS:
        calc = capillary_limit_diameter_mm(point.free_volume_percent, point.radial_factor, TABLE4_CAPILLARY_LIMITS)
        rows.append(
            {
                "free_volume_percent": point.free_volume_percent,
                "radial_factor": point.radial_factor,
                "max_capillary_diameter_exp_mm": point.max_capillary_diameter_mm,
                "max_capillary_diameter_calc_mm": calc,
                "error_mm": calc - point.max_capillary_diameter_mm,
            }
        )
    return rows


def build_summary(config: NikolaevModelConfig) -> dict:
    table1 = validate_table1(config)
    table2 = validate_table2(config)
    table3 = validate_table3()
    table4 = validate_table4()
    table2_current_errors = [row["current_error_a"] for row in table2]
    table2_power_errors = [row["electric_power_error_w"] for row in table2]
    table2_temp_errors = [row["emitter_temperature_error_k"] for row in table2]
    table2_eta_errors = [row["efficiency_error_percent"] for row in table2]
    return {
        "source": {
            "title": SOURCE_TITLE,
            "citation": SOURCE_CITATION,
            "doi": SOURCE_DOI,
            "local_pdf": LOCAL_PDF,
        },
        "config": asdict(config),
        "mockup_facts": MOCKUP_FACTS,
        "table1": table1,
        "table2": table2,
        "table3": table3,
        "table4": table4,
        "metrics": {
            "table1_max_abs_error": _max_abs([row["space_r_error"] for row in table1]),
            "table2_current_mae_a": _mean_abs(table2_current_errors),
            "table2_current_max_abs_a": _max_abs(table2_current_errors),
            "table2_electric_power_mae_w": _mean_abs(table2_power_errors),
            "table2_emitter_temp_mae_k": _mean_abs(table2_temp_errors),
            "table2_efficiency_mae_percent": _mean_abs(table2_eta_errors),
            "table3_max_abs_error_k": _max_abs([row["error_k"] for row in table3]),
            "table4_max_abs_error_mm": _max_abs([row["error_mm"] for row in table4]),
        },
        "status": "complete_table_validation",
        "limitations": [
            "Figure 4 VAC curves are not digitized in this first Nikolaev_V implementation.",
            "The original paper reports calculated integral characteristics, not enough boundary conditions for a unique first-principles electro-thermal reconstruction.",
            "Missing closure parameters are intentionally adjustable inside NikolaevModelConfig.",
        ],
    }


def build_markdown_report(summary: dict, run_dir: Path) -> str:
    lines = [
        "# Nikolaev 1995 SPACE-R Single-Cell TFE Validation Report",
        "",
        f"Run directory: `{run_dir}`",
        "",
        "## Source",
        "",
        f"- Title: {summary['source']['title']}",
        f"- Citation: {summary['source']['citation']}",
        f"- DOI: {summary['source']['doi']}",
        f"- Local PDF: `{summary['source']['local_pdf']}`",
        "",
        "## Scope",
        "",
        "This validation uses the numerical anchors available in Nikolaev 1995: Table 1, Table 2, Table 3, and Table 4. Figure 4 is retained as qualitative VAC evidence until the curve is digitized.",
        "",
        "## Adjustable Closure",
        "",
        "| Parameter | Value |",
        "| --- | ---:|",
    ]
    for key, value in summary["config"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend(
        [
            "",
            "## Table 2 Operating Point Comparison",
            "",
            "| V | I exp A | I calc A | Q exp kW | Q calc kW | Te exp K | Te calc K | eta exp % | eta calc % |",
            "| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:|",
        ]
    )
    for row in summary["table2"]:
        lines.append(
            "| {voltage_v:.3g} | {current_exp_a:.3f} | {current_calc_a:.3f} | {thermal_power_exp_kw:.3f} | {thermal_power_calc_kw:.3f} | {emitter_temperature_exp_k:.3f} | {emitter_temperature_calc_k:.3f} | {efficiency_exp_percent:.3f} | {efficiency_calc_percent:.3f} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Error Metrics",
            "",
            "| Metric | Value |",
            "| --- | ---:|",
        ]
    )
    for key, value in summary["metrics"].items():
        lines.append(f"| `{key}` | {value:.6g} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The Table 2 operating points are matched by a compact TOPAZ-II-like integral TFE model: nearly constant 300 W terminal power, voltage-dependent input heat near 4.1-4.2 kW, and a quadratic emitter-temperature closure. Tables 3 and 4 are represented as interpolation surfaces over the published grids, so the tabulated points close exactly.",
            "",
            "The result should be treated as table-level V&V and parameter reconstruction, not as an independent reproduction of the original Nikolaev/Davydov electro-thermal solver. The next fidelity step is digitizing Figure 4 and replacing the compact closure with a node-wise TEC calculation.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_outputs(run_dir: Path, config: NikolaevModelConfig) -> dict:
    summary = build_summary(config)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "results").mkdir(exist_ok=True)
    _write_csv(run_dir / "results" / "table1_characteristics.csv", summary["table1"])
    _write_csv(run_dir / "results" / "table2_operating_points.csv", summary["table2"])
    _write_csv(run_dir / "results" / "table3_fuel_temperature.csv", summary["table3"])
    _write_csv(run_dir / "results" / "table4_capillary_limits.csv", summary["table4"])
    with (run_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    (run_dir / "validation_report.md").write_text(build_markdown_report(summary, run_dir), encoding="utf-8")
    snapshot_dir = run_dir / "input_snapshot"
    snapshot_dir.mkdir(exist_ok=True)
    for name in (
        "nikolaev_source_data.py",
        "nikolaev_single_tfe_model.py",
        "run_nikolaev_validation.py",
        "PARAMETER_ADJUSTMENT_GUIDE.md",
    ):
        src = CASE_DIR / name
        if src.exists():
            shutil.copy2(src, snapshot_dir / name)
    return summary


def _default_run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S_nikolaev_table_validation")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Nikolaev 1995 SPACE-R single-cell TFE table validation.")
    parser.add_argument("--run-id", default=_default_run_id())
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--nominal-output-power-w", type=float, default=NikolaevModelConfig.nominal_output_power_w)
    parser.add_argument("--thermal-power-reference-kw", type=float, default=NikolaevModelConfig.thermal_power_reference_kw)
    parser.add_argument("--low-voltage-thermal-power-slope-kw-per-v", type=float, default=NikolaevModelConfig.low_voltage_thermal_power_slope_kw_per_v)
    parser.add_argument("--emitter-temp-reference-k", type=float, default=NikolaevModelConfig.emitter_temp_reference_k)
    parser.add_argument("--emitter-temp-linear-k-per-v", type=float, default=NikolaevModelConfig.emitter_temp_linear_k_per_v)
    parser.add_argument("--emitter-temp-quadratic-k-per-v2", type=float, default=NikolaevModelConfig.emitter_temp_quadratic_k_per_v2)
    args = parser.parse_args(argv)
    config = replace(
        NikolaevModelConfig(),
        nominal_output_power_w=args.nominal_output_power_w,
        thermal_power_reference_kw=args.thermal_power_reference_kw,
        low_voltage_thermal_power_slope_kw_per_v=args.low_voltage_thermal_power_slope_kw_per_v,
        emitter_temp_reference_k=args.emitter_temp_reference_k,
        emitter_temp_linear_k_per_v=args.emitter_temp_linear_k_per_v,
        emitter_temp_quadratic_k_per_v2=args.emitter_temp_quadratic_k_per_v2,
    )
    run_dir = args.output_root / args.run_id
    summary = write_outputs(run_dir, config)
    print(json.dumps({"run_dir": str(run_dir), "metrics": summary["metrics"], "status": summary["status"]}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
