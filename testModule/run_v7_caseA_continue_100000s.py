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

from test_core_assemble_v7_caseA import (
    build_v7_case_a_system,
    _case_a_electric_diagnostics,
    _case_a_flow_diagnostics,
    _case_a_reset_design_flows_after_restart,
)
from run_v7_caseA_multipliers_short import (
    collect_tec_stats,
    json_default,
    parse_multipliers,
)
from test_core_assemble_v7_caseA_faststeady import compute_faststeady_energy_audit

TOTAL_POWER_W = 115000.0
WIRE_RESISTANCE_OHM = [
    0.00155199999999970,
    0.00102400000000000,
    0.000336000000000000,
    0.000608000000000000,
]
DEFAULT_RESTART = (
    "testModule/v7_caseA_newpyd_long20000_after_tec_heatfix/"
    "v7_caseA_newpyd_long20000_after_tec_heatfix_latest_restart.npz"
)


def write_json(path: Path, data: Dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=json_default)


def apply_wire_resistance(core: Any) -> None:
    if core.thermo_calc is None:
        return
    n_elem = core.thermo_calc.N_elem
    wire_res = np.asarray(WIRE_RESISTANCE_OHM, dtype=float)
    core.thermo_calc._input_data.resistanceWire = np.tile(wire_res, (n_elem, 1))
    core.thermo_calc.build()
    core.thermo_calc.calculate(verbose=False)

def build_loaded_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v7_case_a_system(
        pipe_n_nodes=args.pipe_n_nodes,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=260.0)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )

    system.initialize_system()
    print(f"Applying wire resistance: {WIRE_RESISTANCE_OHM}", flush=True)
    apply_wire_resistance(core)

    if not os.path.exists(args.restart_in):
        raise FileNotFoundError(f"Restart file not found: {args.restart_in}")
    
    print(f"Loading restart state from: {args.restart_in}", flush=True)
    system.load_global_state(args.restart_in)
    
    # Re-apply setups that might be overwritten or need reset
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=260.0)
    apply_wire_resistance(core)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)
    core.post_step(0.0, float(system.global_time))
    apply_wire_resistance(core)
        
    return build

def collect_record(build: Dict[str, Any], start_time: float) -> Dict[str, Any]:
    system = build["system"]
    core = build["core"]
    
    # Refresh solid boundaries to get accurate current_flux
    system._refresh_solid_boundary_cache(update_flux=True, current_time=float(system.global_time))
    
    electric = _case_a_electric_diagnostics(core)
    energy = compute_faststeady_energy_audit({
        "build": build,
        "system": system,
        "core": core,
    })
    tec_stats = collect_tec_stats(build)
    
    terminal_power = float(electric["tec_total_electric_power_w"])
    current = float(electric["tec_total_current_a"])
    core_heat = float(energy["core_heat_power_w"])
    coolant_heat = float(energy["coolant_heat_pickup_w"])
    outer_wall_radiation = float(energy["outer_wall_radiation_w"])
    balance_residual = core_heat - terminal_power - coolant_heat - outer_wall_radiation
    
    record = {
        "absolute_time_s": float(system.global_time),
        "relative_time_s": float(system.global_time) - start_time,
        "last_dt_s": float(getattr(system, "_last_dt", np.nan)),
        "wire_resistance_ohm": list(WIRE_RESISTANCE_OHM),
        "electric_power_w": terminal_power,
        "current_a": current,
        "outer_wall_radiation_w": outer_wall_radiation,
        "core_heat_power_w": core_heat,
        "coolant_heat_pickup_w": coolant_heat,
        "balance_residual_w": balance_residual,
    }
    
    # - 不同通道的燃料最高温度，燃料平均温度
    for tfe_name, tfe in build["tfes"].items():
        pellet = tfe.solids["pellet"]
        record[f"{tfe_name}_fuel_pellet_max_k"] = float(np.max(pellet.T))
        record[f"{tfe_name}_fuel_pellet_mean_k"] = float(np.mean(pellet.T))
        
    # - 发射极最高温度，发射极平均温度
    # - 接收极最高温度，接收极平均温度
    for ring_name, stats in tec_stats.get("by_ring", {}).items():
        record[f"{ring_name}_emitter_max_k"] = stats["emitter_temperature_max_k"]
        record[f"{ring_name}_emitter_mean_k"] = stats["emitter_temperature_mean_k"]
        record[f"{ring_name}_collector_max_k"] = stats["collector_temperature_max_k"]
        record[f"{ring_name}_collector_mean_k"] = stats["collector_temperature_mean_k"]
        
    # - 发射极外表面总热流，发射极外表面电子导热热流
    # - 接收极内表面总热流，接收极内表面电子导热热流
    for tfe_name, tfe in build["tfes"].items():
        mult = float(build["ring_multipliers"][tfe_name])
        # Emitter right boundary (heat leaving is positive for our reporting)
        q_em_out = -np.asarray(tfe.solids["emitter"].boundaries["right"].current_flux, dtype=float)
        # Collector left boundary (heat entering is positive for our reporting)
        q_col_in = np.asarray(tfe.solids["collector"].boundaries["left"].current_flux, dtype=float)
        ring_energy = (
            tec_stats.get("by_ring", {})
            .get(tfe_name, {})
            .get("tec_group_energy", {})
        )
        record[f"{tfe_name}_emitter_outer_total_heat_w"] = float(np.sum(q_em_out)) * mult
        record[f"{tfe_name}_emitter_outer_electron_heat_w"] = float(
            ring_energy.get("emitter_electron_heat_removed_positive_w", 0.0)
        )
        record[f"{tfe_name}_collector_inner_total_heat_w"] = float(np.sum(q_col_in)) * mult
        record[f"{tfe_name}_collector_inner_electron_heat_w"] = float(
            ring_energy.get(
                "collector_electron_heat_source_if_sideAreaC_w",
                ring_energy.get("collector_electron_heat_source_signed_w", 0.0),
            )
        )
    
    return record


def _diag_array(value: Any) -> np.ndarray:
    if isinstance(value, (list, tuple)):
        return np.asarray(value)
    return np.asarray([value])


def append_diag_to_restart(restart_path: Path, record: Dict[str, Any]) -> None:
    tmp_path = restart_path.with_suffix(".tmp.npz")
    with np.load(restart_path, allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    for key, value in record.items():
        arrays[f"Diag/{key}"] = _diag_array(value)
    np.savez_compressed(tmp_path, **arrays)
    os.replace(tmp_path, restart_path)


def advance_to(build: Dict[str, Any], stop_time: float, args: argparse.Namespace) -> None:
    system = build["system"]
    core = build["core"]
    while float(system.global_time) < stop_time - 1.0e-10:
        core.update_neutronic_power(
            p_total=TOTAL_POWER_W,
            p_fiss=TOTAL_POWER_W,
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

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="100,000s Continuation Run with modified wire resistance.")
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v7_caseA_100k_run")
    parser.add_argument("--case-prefix", default="v7_caseA_100k_run")
    parser.add_argument("--ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,18"))
    parser.add_argument("--tec-ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,15"))
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--duration", type=float, default=100000.0)
    parser.add_argument("--record-interval", type=float, default=200.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    build = build_loaded_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)
    record_interval = float(args.record_interval)

    # Only one restart file that gets overwritten
    latest_restart = output_dir / f"{args.case_prefix}_latest_restart.npz"
    latest_state = output_dir / f"{args.case_prefix}_latest_state.json"
    history_csv = output_dir / f"{args.case_prefix}_history.csv"

    next_record_time = start_time + record_interval
    fieldnames = None
    recorded_times = set()
    if history_csv.exists():
        with history_csv.open("r", newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            for row in reader:
                try:
                    recorded_times.add(round(float(row["absolute_time_s"]), 9))
                except (KeyError, TypeError, ValueError):
                    pass

    print(f"=== Starting 100k run ===", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"latest_restart={latest_restart}", flush=True)
    print(f"latest_state={latest_state}", flush=True)
    print(f"history_csv={history_csv}", flush=True)

    while next_record_time <= target_time + 1.0e-10:
        advance_to(build, next_record_time, args)
        
        # Save and overwrite the exact same file
        system.save_global_state(str(latest_restart))
        
        record = collect_record(build, start_time)
        append_diag_to_restart(latest_restart, record)
        
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
                    
        record_time_key = round(float(record["absolute_time_s"]), 9)
        if record_time_key not in recorded_times:
            with history_csv.open("a", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writerow(record)
            recorded_times.add(record_time_key)
        else:
            print(
                f"history_csv already has t_abs={record['absolute_time_s']:.1f}s; "
                "skipping duplicate row",
                flush=True,
            )

        write_json(latest_state, {
            "case_prefix": args.case_prefix,
            "restart_in": args.restart_in,
            "latest_restart": str(latest_restart),
            "history_csv": str(history_csv),
            "ring_multipliers": args.ring_multipliers,
            "tec_ring_multipliers": args.tec_ring_multipliers,
            "target_voltage_v": args.target_voltage,
            "duration_target_s": args.duration,
            "absolute_start_time_s": start_time,
            "absolute_target_time_s": target_time,
            "record_interval_s": args.record_interval,
            "max_dt_s": args.max_dt,
            "safety_factor": args.safety_factor,
            "wire_resistance_ohm": WIRE_RESISTANCE_OHM,
            "record": record,
        })
            
        print(
            f"t_abs={record['absolute_time_s']:.1f}s "
            f"P={record['electric_power_w']:.3f}W "
            f"I={record['current_a']:.3f}A "
            f"Rterm={record['balance_residual_w']:.3f}W",
            flush=True,
        )
        next_record_time += record_interval

    print("Long run completed.", flush=True)

if __name__ == "__main__":
    main()
