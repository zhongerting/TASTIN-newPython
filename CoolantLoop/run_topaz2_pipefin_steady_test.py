import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from model_topaz2_tube_fin_radiator import make_default_args, run_case


def _json_default(obj):
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _read_history(path: Path) -> List[Dict[str, float]]:
    rows = []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def _slope_over_tail(rows: List[Dict[str, float]], field: str, window_s: float) -> float:
    if len(rows) < 2:
        return float("nan")
    t_end = rows[-1]["time_s"]
    tail = [row for row in rows if row["time_s"] >= t_end - window_s]
    if len(tail) < 2:
        tail = rows[-2:]
    dt = tail[-1]["time_s"] - tail[0]["time_s"]
    if dt <= 0.0:
        return float("nan")
    return (tail[-1][field] - tail[0][field]) / dt


def summarize_case(latest: Dict[str, Any], tail_window_s: float, wall_time_s: float) -> Dict[str, Any]:
    history_path = Path(latest["history_path"])
    rows = _read_history(history_path)
    record = latest["latest_record"]
    q_total = float(record["q_total_radiation_w"])
    energy_residual = float(record["energy_residual_w"])
    area = float(record["total_effective_area_m2"])
    q_tube = float(record["q_tube_radiation_w"])
    q_fin = float(record["q_fin_radiation_w"])

    summary = {
        "case_prefix": latest["case_prefix"],
        "history_path": history_path,
        "time_s": float(record["time_s"]),
        "wall_time_s": float(wall_time_s),
        "speed_ratio_sim_s_per_wall_s": float(record["time_s"]) / wall_time_s if wall_time_s > 0.0 else float("nan"),
        "mixed_outlet_temperature_k": float(record["mixed_outlet_temperature_k"]),
        "outlet_delta_from_target_k": float(record["outlet_delta_from_target_k"]),
        "outlet_temperature_slope_k_per_s": _slope_over_tail(
            rows,
            "mixed_outlet_temperature_k",
            tail_window_s,
        ),
        "tail_window_s": float(tail_window_s),
        "q_total_radiation_w": q_total,
        "q_tube_radiation_w": q_tube,
        "q_fin_radiation_w": q_fin,
        "fin_radiation_fraction": q_fin / q_total if q_total else float("nan"),
        "energy_residual_w": energy_residual,
        "energy_residual_fraction": energy_residual / q_total if q_total else float("nan"),
        "total_effective_area_m2": area,
        "area_within_reference_band": 7.15 <= area <= 7.30,
        "tube_flow_rel_spread": float(record["tube_flow_rel_spread"]),
        "fin_root_to_tip_delta_mean_k": float(record["fin_root_to_tip_delta_mean_k"]),
        "fin_iteration_mean": float(record.get("fin_iteration_mean", float("nan"))),
        "fin_iteration_max": float(record.get("fin_iteration_max", float("nan"))),
        "fin_max_delta_mean_k": float(record.get("fin_max_delta_mean_k", float("nan"))),
        "fin_max_delta_max_k": float(record.get("fin_max_delta_max_k", float("nan"))),
        "fin_warm_start_fraction": float(record.get("fin_warm_start_fraction", float("nan"))),
        "n_fin_width": int(record["n_fin_width"]),
        "finite_numeric_outputs": all(
            math.isfinite(float(record[key]))
            for key in (
                "mixed_outlet_temperature_k",
                "q_total_radiation_w",
                "energy_residual_w",
                "fin_root_to_tip_delta_mean_k",
            )
        ),
    }
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the TOPAZ-II pipe-fin radiator steady test case.")
    parser.add_argument("--output-dir", default="CoolantLoop/topaz2_pipefin_steady_eps080_500s")
    parser.add_argument("--case-prefix", default="topaz2_pipefin_steady_eps080_500s")
    parser.add_argument("--duration", type=float, default=500.0)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--max-dt", type=float, default=0.2)
    parser.add_argument("--tail-window-s", type=float, default=100.0)
    parser.add_argument("--tube-emissivity", type=float, default=0.80)
    parser.add_argument("--fin-emissivity", type=float, default=0.80)
    parser.add_argument("--solid-ode-method", default="RK45", choices=("RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"))
    parser.add_argument("--n-fin-width", type=int, default=12)
    parser.add_argument("--n-axial", type=int, default=8)
    parser.add_argument("--fin-area-scale", type=float, default=0.35)
    parser.add_argument("--hydraulic-calibrated", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case_args = make_default_args(
        output_dir=args.output_dir,
        case_prefix=args.case_prefix,
        duration=args.duration,
        record_interval=args.record_interval,
        max_dt=args.max_dt,
        tube_emissivity=args.tube_emissivity,
        fin_emissivity=args.fin_emissivity,
        solid_ode_method=args.solid_ode_method,
        n_fin_width=args.n_fin_width,
        n_axial=args.n_axial,
        fin_area_scale=args.fin_area_scale,
        hydraulic_calibrated=args.hydraulic_calibrated,
    )
    wall_start = time.perf_counter()
    latest = run_case(case_args)
    wall_time_s = time.perf_counter() - wall_start
    summary = summarize_case(latest, tail_window_s=args.tail_window_s, wall_time_s=wall_time_s)

    summary_path = Path(args.output_dir) / f"{args.case_prefix}_steady_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=_json_default)

    print(json.dumps(summary, indent=2, ensure_ascii=False, default=_json_default))


if __name__ == "__main__":
    main()
