import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from run_v7_caseA_multipliers_short import (
    build_loaded_case,
    collect_tec_stats,
    json_default,
    parse_multipliers,
)
from test_core_assemble_v7_caseA import (
    _case_a_electric_diagnostics,
    _case_a_flow_diagnostics,
)
from test_core_assemble_v7_caseA_faststeady import compute_faststeady_energy_audit


DEFAULT_RESTART = (
    "testModule/v7_caseA_newpyd_long20000_from_t23800/"
    "v7_caseA_newpyd_long20000_from_t23800_latest_restart.npz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Long V7 CaseA run with new non-uniform TEC side-area interface."
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v7_caseA_newpyd_long20000_after_tec_heatfix")
    parser.add_argument("--case-prefix", default="v7_caseA_newpyd_long20000_after_tec_heatfix")
    parser.add_argument("--ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,18"))
    parser.add_argument("--tec-ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,15"))
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--duration", type=float, default=20000.0)
    parser.add_argument("--target-time", type=float, default=None)
    parser.add_argument("--record-start-time", type=float, default=None)
    parser.add_argument("--record-interval", type=float, default=200.0)
    parser.add_argument("--snapshot-interval", type=float, default=2000.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    return parser.parse_args()


def write_json(path: Path, data: Dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=json_default)


def time_label(value: float) -> str:
    if abs(value - round(value)) < 1.0e-6:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def ring_temp_columns(tec_stats: Dict[str, Any]) -> Dict[str, float]:
    row: Dict[str, float] = {}
    for name, stats in tec_stats.get("by_ring", {}).items():
        row[f"{name}_emitter_mean_k"] = stats["emitter_temperature_mean_k"]
        row[f"{name}_emitter_max_k"] = stats["emitter_temperature_max_k"]
        row[f"{name}_collector_mean_k"] = stats["collector_temperature_mean_k"]
        row[f"{name}_collector_max_k"] = stats["collector_temperature_max_k"]
        row[f"{name}_j_mean_nonzero_a_cm2"] = stats["j_mean_nonzero_a_cm2"]
        row[f"{name}_j_max_a_cm2"] = stats["j_max_a_cm2"]
    return row


def fuel_pellet_max_columns(build: Dict[str, Any]) -> Dict[str, float]:
    row: Dict[str, float] = {}
    for name, tfe in build["tfes"].items():
        pellet = tfe.solids["pellet"]
        row[f"{name}_fuel_pellet_max_k"] = float(np.max(pellet.T))
    return row


def collect_record(build: Dict[str, Any], start_time: float, restart_path: Path) -> Dict[str, Any]:
    system = build["system"]
    core = build["core"]
    electric = _case_a_electric_diagnostics(core)
    flow = _case_a_flow_diagnostics(build)
    energy = compute_faststeady_energy_audit({
        "build": build,
        "system": system,
        "core": core,
    })
    tec_stats = collect_tec_stats(build)
    tec_totals = tec_stats.get("totals", {})
    terminal_power = float(electric["tec_total_electric_power_w"])
    electron_boundary = float(tec_totals.get("electron_boundary_heat_difference_w") or 0.0)
    total_joule = float(tec_totals.get("total_joule_heat_w") or 0.0)
    terminal_plus_joule = terminal_power + total_joule

    record: Dict[str, Any] = {
        "absolute_time_s": float(system.global_time),
        "relative_time_s": float(system.global_time) - start_time,
        "last_dt_s": float(getattr(system, "last_dt", np.nan)),
        "restart_file": str(restart_path),
        "voltage_v": electric["tec_total_voltage_v"],
        "current_a": electric["tec_total_current_a"],
        "electric_power_w": terminal_power,
        "terminal_electric_power_w": terminal_power,
        "mass_flow_kg_s": flow["inlet_total_macro_flow_kg_s"],
        "tfe_total_macro_flow_kg_s": flow["tfe_total_macro_flow_kg_s"],
        "inlet_plenum_temperature_k": energy["inlet_plenum_temperature_k"],
        "outlet_plenum_temperature_k": energy["outlet_plenum_temperature_k"],
        "core_heat_power_w": energy["core_heat_power_w"],
        "coolant_heat_pickup_w": energy["coolant_heat_pickup_w"],
        "outer_wall_radiation_w": energy["outer_wall_radiation_w"],
        "balance_residual_w": energy["balance_residual_w"],
        "balance_residual_using_terminal_power_w": energy["balance_residual_w"],
        "balance_residual_percent": energy["balance_residual_percent"],
        "emitter_electron_heat_source_signed_w": tec_totals.get("emitter_electron_heat_source_signed_w"),
        "emitter_electron_heat_removed_positive_w": tec_totals.get("emitter_electron_heat_removed_positive_w"),
        "collector_electron_heat_source_signed_w": tec_totals.get("collector_electron_heat_source_signed_w"),
        "electron_boundary_heat_difference_w": electron_boundary,
        "emitter_joule_heat_w": tec_totals.get("emitter_joule_heat_w"),
        "collector_joule_heat_w": tec_totals.get("collector_joule_heat_w"),
        "total_joule_heat_w": total_joule,
        "terminal_plus_joule_w": terminal_plus_joule,
        "tec_energy_gap_w": terminal_plus_joule - electron_boundary,
    }
    record.update(fuel_pellet_max_columns(build))
    record.update(ring_temp_columns(tec_stats))
    return record


def append_tec_distribution_records(
    build: Dict[str, Any],
    start_time: float,
    distribution_csv: Path,
    write_header: bool,
) -> None:
    core = build["core"]
    thermo_calc = core.thermo_calc
    if thermo_calc is None:
        return

    reference_tfe = next(iter(build["tfes"].values()))
    y_faces = np.asarray(reference_tfe.common_y_faces, dtype=float)
    y_centers = 0.5 * (y_faces[:-1] + y_faces[1:])
    n_lower = int(getattr(reference_tfe, "n_lower", 0))
    n_active = int(getattr(reference_tfe, "n_active", len(y_centers)))
    active_start = n_lower
    active_stop = n_lower + n_active

    input_data = thermo_calc._input_data
    fieldnames = [
        "absolute_time_s",
        "relative_time_s",
        "representative_tfe",
        "thermal_multiplier",
        "tec_multiplier",
        "axial_node",
        "y_center_m",
        "is_active_node",
        "T_emitter_k",
        "T_collector_k",
        "J_A_cm2",
        "node_current_A",
        "UE_minus_UC_v",
        "Vd_v",
    ]

    absolute_time = float(build["system"].global_time)
    relative_time = absolute_time - start_time
    ring_multipliers = build["ring_multipliers"]
    tec_ring_multipliers = build["tec_ring_multipliers"]

    idx = 0
    with distribution_csv.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        for name, tec_multiplier in tec_ring_multipliers.items():
            tec_multiplier = int(tec_multiplier)
            if tec_multiplier <= 0:
                continue
            res = thermo_calc.get_tec_results(idx)
            area = np.asarray(input_data.sideAreaE[idx], dtype=float)
            j = np.asarray(res["J"], dtype=float)
            te = np.asarray(res["TE"], dtype=float)
            tc = np.asarray(res["TC"], dtype=float)
            ue = np.asarray(res["UE"], dtype=float)
            uc = np.asarray(res["UC"], dtype=float)
            vd = np.asarray(res["Vd"], dtype=float)
            for node in range(len(j)):
                writer.writerow({
                    "absolute_time_s": absolute_time,
                    "relative_time_s": relative_time,
                    "representative_tfe": name,
                    "thermal_multiplier": int(ring_multipliers[name]),
                    "tec_multiplier": tec_multiplier,
                    "axial_node": node,
                    "y_center_m": float(y_centers[node]),
                    "is_active_node": active_start <= node < active_stop,
                    "T_emitter_k": float(te[node]),
                    "T_collector_k": float(tc[node]),
                    "J_A_cm2": float(j[node]),
                    "node_current_A": float(j[node] * 1.0e4 * area[node]),
                    "UE_minus_UC_v": float(ue[node] - uc[node]),
                    "Vd_v": float(vd[node]),
                })
            idx += tec_multiplier


def advance_to(build: Dict[str, Any], stop_time: float, args: argparse.Namespace) -> None:
    system = build["system"]
    core = build["core"]
    while float(system.global_time) < stop_time - 1.0e-10:
        core.update_neutronic_power(
            p_total=115000.0,
            p_fiss=115000.0,
            p_decay=0.0,
            alpha=1.0,
        )
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=float(args.max_dt),
            safety_factor=float(args.safety_factor),
        )
        dt = min(float(dt), float(args.max_dt), stop_time - float(system.global_time))
        system.step(dt)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    build = build_loaded_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    record_start_time = (
        start_time
        if args.record_start_time is None
        else float(args.record_start_time)
    )
    target_time = (
        start_time + float(args.duration)
        if args.target_time is None
        else float(args.target_time)
    )
    if target_time < start_time - 1.0e-10:
        raise ValueError(
            f"target_time={target_time} is earlier than restart time {start_time}."
        )
    record_interval = float(args.record_interval)

    latest_restart = output_dir / f"{args.case_prefix}_latest_restart.npz"
    latest_state = output_dir / f"{args.case_prefix}_latest_state.json"
    history_csv = output_dir / f"{args.case_prefix}_history.csv"
    distribution_csv = output_dir / f"{args.case_prefix}_tec_distribution_history.csv"

    elapsed_for_records = max(0.0, start_time - record_start_time)
    next_record_index = int(np.floor(elapsed_for_records / record_interval + 1.0e-10)) + 1
    next_record_time = record_start_time + float(next_record_index) * record_interval
    if next_record_time <= start_time + 1.0e-10:
        next_record_time += record_interval

    fieldnames = None
    if history_csv.exists():
        with history_csv.open("r", newline="") as f:
            reader = csv.reader(f)
            fieldnames = next(reader, None)

    print("=== V7 CaseA new-pyd long run ===", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"record_start_time={record_start_time:.6f}", flush=True)
    print(f"snapshot_interval={args.snapshot_interval:.6f}", flush=True)
    print(f"history_csv={history_csv}", flush=True)
    print(f"distribution_csv={distribution_csv}", flush=True)
    print(f"latest_restart={latest_restart}", flush=True)

    while next_record_time <= target_time + 1.0e-10:
        advance_to(build, next_record_time, args)
        system.save_global_state(str(latest_restart))
        restart_for_record = latest_restart
        rel_time = float(system.global_time) - record_start_time
        snapshot_interval = float(args.snapshot_interval)
        if snapshot_interval > 0.0:
            snapshot_index = round(rel_time / snapshot_interval)
            if abs(rel_time - snapshot_index * snapshot_interval) < 1.0e-6:
                snapshot_path = output_dir / (
                    f"{args.case_prefix}_restart_rel{time_label(rel_time)}"
                    f"_abs{time_label(float(system.global_time))}.npz"
                )
                system.save_global_state(str(snapshot_path))
                restart_for_record = snapshot_path
        record = collect_record(build, record_start_time, restart_for_record)

        if fieldnames is None:
            fieldnames = list(record.keys())
            with history_csv.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
        elif set(record.keys()) != set(fieldnames):
            missing = sorted(set(fieldnames) - set(record.keys()))
            extra = sorted(set(record.keys()) - set(fieldnames))
            raise ValueError(
                f"Existing history CSV columns do not match current record. "
                f"missing={missing}, extra={extra}"
            )

        with history_csv.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writerow(record)
        append_tec_distribution_records(
            build,
            record_start_time,
            distribution_csv,
            write_header=not distribution_csv.exists(),
        )

        write_json(latest_state, {
            "case_prefix": args.case_prefix,
            "restart_in": args.restart_in,
            "latest_restart": str(latest_restart),
            "history_csv": str(history_csv),
            "distribution_csv": str(distribution_csv),
            "ring_multipliers": args.ring_multipliers,
            "tec_ring_multipliers": args.tec_ring_multipliers,
            "target_voltage_v": args.target_voltage,
            "duration_target_s": args.duration,
            "absolute_target_time_s": target_time,
            "record_start_time_s": record_start_time,
            "record_interval_s": args.record_interval,
            "snapshot_interval_s": args.snapshot_interval,
            "record": record,
        })

        print(
            f"t_rel={record['relative_time_s']:.1f}s "
            f"I={record['current_a']:.3f}A "
            f"P={record['electric_power_w']:.3f}W "
            f"Qcool={record['coolant_heat_pickup_w']:.3f}W "
            f"Qe_emit_removed={record['emitter_electron_heat_removed_positive_w']:.3f}W "
            f"Qe_coll={record['collector_electron_heat_source_signed_w']:.3f}W "
            f"TECgap={record['tec_energy_gap_w']:.3f}W "
            f"Rterm={record['balance_residual_using_terminal_power_w']:.3f}W",
            flush=True,
        )
        next_record_time += record_interval

    print("Long run completed.", flush=True)


if __name__ == "__main__":
    main()
