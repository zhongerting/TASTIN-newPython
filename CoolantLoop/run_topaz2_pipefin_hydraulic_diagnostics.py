import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from model_topaz2_tube_fin_radiator import (
    build_model,
    collect_diagnostics,
    header_flow_diagnostics,
    initialize_flow_guess,
    make_default_args,
    pressure_budget_diagnostics,
    tube_flow_distribution,
)


def _json_default(obj):
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _parse_float_list(text: str) -> List[float]:
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def _write_csv(path: Path, rows: Iterable[Dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _run_diagnostic_case(label: str, base_args: argparse.Namespace, overrides: Dict[str, Any]) -> Dict[str, Any]:
    case_args = argparse.Namespace(**vars(base_args))
    for key, value in overrides.items():
        setattr(case_args, key, value)

    build = build_model(case_args)
    initialize_flow_guess(build)
    system = build["system"]
    system.initialize_system(
        dt_init=float(case_args.init_dt),
        tol=float(case_args.hydraulic_tol),
        max_iter=int(case_args.hydraulic_max_iter),
    )

    thermal_duration = float(getattr(case_args, "run_thermal_duration", 0.0))
    while system.global_time < thermal_duration - 1.0e-12:
        dt = min(float(case_args.max_dt), thermal_duration - system.global_time)
        system.step(
            dt=dt,
            inner_iter=int(case_args.inner_iter),
            convergence_tol=float(case_args.convergence_tol),
        )

    tube_diag = tube_flow_distribution(build)
    header_diag = header_flow_diagnostics(build)
    pressure_diag = pressure_budget_diagnostics(build)
    thermal_diag = collect_diagnostics(build)

    return {
        "label": label,
        "overrides": overrides,
        "thermal_duration_s": thermal_duration,
        "tube_inlet_k_loss": float(case_args.tube_inlet_k_loss),
        "tube_outlet_k_loss": float(case_args.tube_outlet_k_loss),
        "header_inner_diameter_m": float(case_args.header_inner_diameter_m),
        "hydraulic_calibrated": bool(getattr(case_args, "hydraulic_calibrated", False)),
        "tube_distribution": tube_diag,
        "header_diagnostics": header_diag,
        "pressure_budget": pressure_diag,
        "thermal_diagnostics": thermal_diag,
    }


def _summary_row(result: Dict[str, Any]) -> Dict[str, Any]:
    tube = result["tube_distribution"]
    header = result["header_diagnostics"]
    return {
        "label": result["label"],
        "thermal_duration_s": result["thermal_duration_s"],
        "hydraulic_calibrated": result["hydraulic_calibrated"],
        "tube_inlet_k_loss": result["tube_inlet_k_loss"],
        "tube_outlet_k_loss": result["tube_outlet_k_loss"],
        "header_inner_diameter_m": result["header_inner_diameter_m"],
        "mean_tube_flow_kg_s": tube["mean_tube_flow_kg_s"],
        "min_tube_flow_kg_s": tube["min_tube_flow_kg_s"],
        "min_tube_index": tube["min_tube_index"],
        "max_tube_flow_kg_s": tube["max_tube_flow_kg_s"],
        "max_tube_index": tube["max_tube_index"],
        "flow_spread_over_mean": tube["flow_spread_over_mean"],
        "max_over_min": tube["max_over_min"],
        "symmetry_error_kg_s": tube["symmetry_error_kg_s"],
        "upper_header_max_abs_flow_kg_s": header["upper_header"]["max_abs_flow_kg_s"],
        "lower_header_max_abs_flow_kg_s": header["lower_header"]["max_abs_flow_kg_s"],
        "upper_header_reversal_count": header["upper_header"]["direction_reversal_count"],
        "lower_header_reversal_count": header["lower_header"]["direction_reversal_count"],
    }


def _plot_selected_cases(output_dir: Path, case_prefix: str, results: List[Dict[str, Any]]) -> List[Path]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []

    selected = [
        result
        for result in results
        if result["label"] in {"baseline", "hydraulic_calibrated"}
    ]
    paths = []
    for result in selected:
        flows = np.array(
            [row["mass_flow_kg_s"] for row in result["tube_distribution"]["tube_flows"]],
            dtype=float,
        )
        tubes = np.arange(1, len(flows) + 1)
        mean = float(np.mean(flows))
        min_idx = int(np.argmin(flows))
        max_idx = int(np.argmax(flows))

        fig, ax = plt.subplots(figsize=(12, 5.5), dpi=160)
        ax.bar(tubes, flows * 1000.0, color="#4C78A8", alpha=0.72, width=0.82)
        ax.plot(tubes, flows * 1000.0, color="#1F4E79", linewidth=1.4)
        ax.axhline(mean * 1000.0, color="#F58518", linestyle="--", linewidth=1.5)
        ax.scatter([tubes[max_idx]], [flows[max_idx] * 1000.0], color="#C62828", s=40, zorder=5)
        ax.scatter([tubes[min_idx]], [flows[min_idx] * 1000.0], color="#2E7D32", s=40, zorder=5)
        ax.set_title(f"TOPAZ-II tube flow distribution: {result['label']}")
        ax.set_xlabel("Radiator tube index")
        ax.set_ylabel("Mass flow (g/s)")
        ax.set_xlim(0.3, len(flows) + 0.7)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()

        path = output_dir / f"{case_prefix}_{result['label']}_tube_flow.png"
        fig.savefig(path)
        plt.close(fig)
        paths.append(path)
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Diagnose TOPAZ-II pipe-fin radiator hydraulic distribution.")
    parser.add_argument("--output-dir", default="CoolantLoop/topaz2_pipefin_hydraulic_diagnostics")
    parser.add_argument("--case-prefix", default="topaz2_pipefin_hydraulic_diagnostics")
    parser.add_argument("--run-thermal-duration", type=float, default=0.0)
    parser.add_argument("--max-dt", type=float, default=0.2)
    parser.add_argument("--scan-tube-k", default="2,5,10,20,50,100,200")
    parser.add_argument("--scan-header-id-mm", default="20,22,24,26,30,40")
    parser.add_argument("--no-scans", action="store_true")
    parser.add_argument("--tube-emissivity", type=float, default=0.80)
    parser.add_argument("--fin-emissivity", type=float, default=0.80)
    parser.add_argument("--n-axial", type=int, default=8)
    parser.add_argument("--n-fin-width", type=int, default=12)
    parser.add_argument("--fin-area-scale", type=float, default=0.35)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_args = make_default_args(
        output_dir=args.output_dir,
        case_prefix=args.case_prefix,
        max_dt=args.max_dt,
        tube_emissivity=args.tube_emissivity,
        fin_emissivity=args.fin_emissivity,
        n_axial=args.n_axial,
        n_fin_width=args.n_fin_width,
        fin_area_scale=args.fin_area_scale,
    )
    base_args.run_thermal_duration = float(args.run_thermal_duration)

    cases = [
        ("baseline", {}),
        ("hydraulic_calibrated", {"hydraulic_calibrated": True}),
    ]
    if not args.no_scans:
        for k_loss in _parse_float_list(args.scan_tube_k):
            cases.append(
                (
                    f"tube_K_each_{k_loss:g}",
                    {"tube_inlet_k_loss": k_loss, "tube_outlet_k_loss": k_loss},
                )
            )
        for header_id_mm in _parse_float_list(args.scan_header_id_mm):
            cases.append(
                (
                    f"header_id_{header_id_mm:g}mm",
                    {"header_inner_diameter_m": header_id_mm / 1000.0},
                )
            )

    results = [_run_diagnostic_case(label, base_args, overrides) for label, overrides in cases]
    summary_rows = [_summary_row(result) for result in results]

    summary_path = output_dir / f"{args.case_prefix}_summary.json"
    summary_csv_path = output_dir / f"{args.case_prefix}_summary.csv"
    tube_csv_path = output_dir / f"{args.case_prefix}_tube_flows.csv"
    pressure_csv_path = output_dir / f"{args.case_prefix}_pressure_budget.csv"

    with summary_path.open("w", encoding="utf-8") as f:
        json.dump({"cases": results}, f, indent=2, ensure_ascii=False, default=_json_default)
    _write_csv(summary_csv_path, summary_rows)

    tube_rows = []
    pressure_rows = []
    for result in results:
        for row in result["tube_distribution"]["tube_flows"]:
            tube_rows.append({"label": result["label"], **row})
        for row in result["pressure_budget"]:
            pressure_rows.append({"label": result["label"], **row})
    _write_csv(tube_csv_path, tube_rows)
    _write_csv(pressure_csv_path, pressure_rows)

    plot_paths = _plot_selected_cases(output_dir, args.case_prefix, results)

    printed = {
        "summary_json": summary_path,
        "summary_csv": summary_csv_path,
        "tube_flows_csv": tube_csv_path,
        "pressure_budget_csv": pressure_csv_path,
        "plots": plot_paths,
        "summary": summary_rows,
    }
    print(json.dumps(printed, indent=2, ensure_ascii=False, default=_json_default))


if __name__ == "__main__":
    main()
