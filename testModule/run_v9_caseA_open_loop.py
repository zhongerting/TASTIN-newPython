import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from run_v8_caseA_common import (
    DEFAULT_COOLANT_MATERIAL,
    DEFAULT_SOLID_ODE_METHOD,
    TOTAL_POWER_W,
    apply_solid_ode_method,
    apply_wire_resistance,
    get_solid_ode_methods,
    get_wire_resistance,
    json_default,
    parse_solid_ode_method,
    parse_v8_multipliers,
    passive_tec_source_totals,
)
from test_core_assemble_v9_caseA import (
    V9_CASE_VERSION,
    build_v9_case_a_system,
    reset_v9_design_flows,
    v9_basic_diagnostics,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V9 CaseA open external-piping model.")
    parser.add_argument("--restart-in", default=None)
    parser.add_argument("--output-dir", default="testModule/v9_caseA_open_loop")
    parser.add_argument("--case-prefix", default="v9_caseA_open_loop")
    parser.add_argument("--duration", type=float, default=200.0)
    parser.add_argument("--record-interval", type=float, default=20.0)
    parser.add_argument("--restart-interval", type=float, default=20.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--external-pipe-n-nodes", type=int, default=5)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--thermo-update-interval", type=float, default=0.8)
    parser.add_argument(
        "--disable-tec-coupled",
        action="store_true",
        help="Disable TEC coupling for cold-start hydraulic/topology smoke runs.",
    )
    parser.add_argument("--inlet-temperature-k", type=float, default=743.0)
    parser.add_argument("--total-inlet-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--outlet-pressure-pa", type=float, default=160000.0)
    parser.add_argument("--connector-volume-m3", type=float, default=1.0e-5)
    parser.add_argument("--connector-length-m", type=float, default=0.02)
    parser.add_argument("--coolant-material", default=DEFAULT_COOLANT_MATERIAL)
    parser.add_argument(
        "--solid-ode-method",
        type=parse_solid_ode_method,
        default=DEFAULT_SOLID_ODE_METHOD,
    )
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=False),
        default=parse_v8_multipliers("1,6,12,15,3"),
    )
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=True),
        default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True),
    )
    return parser.parse_args()


def flatten_for_csv(record: Dict[str, Any]) -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    for key, value in record.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                flat[f"{key}_{sub_key}"] = sub_value
        elif isinstance(value, (list, tuple)):
            for idx, sub_value in enumerate(value, start=1):
                flat[f"{key}_{idx}"] = sub_value
        else:
            flat[key] = value
    return flat


def append_row(path: Path, fieldnames: List[str], row: Dict[str, Any], *, write_header: bool) -> None:
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def build_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v9_case_a_system(
        inlet_temperature_k=args.inlet_temperature_k,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        outlet_pressure_pa=args.outlet_pressure_pa,
        pipe_n_nodes=args.pipe_n_nodes,
        external_pipe_n_nodes=args.external_pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = not bool(args.disable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    apply_solid_ode_method(build, args.solid_ode_method)
    system.initialize_system()
    if args.restart_in:
        system.load_global_state(args.restart_in)
        apply_solid_ode_method(build, args.solid_ode_method)
    core.point_reactor = None
    core.enable_tec_coupled = not bool(args.disable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    if core.enable_tec_coupled:
        core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=150.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    reset_v9_design_flows(build)
    core.post_step(0.0, float(system.global_time))
    if core.enable_tec_coupled and core.thermo_calc is not None:
        apply_wire_resistance(core)
        core._last_thermo_update_time = float(system.global_time)
    core.pre_step(0.0, float(system.global_time))
    build["solid_ode_method"] = args.solid_ode_method
    build["solid_ode_methods"] = get_solid_ode_methods(build)
    build["wire_resistance_ohm"] = get_wire_resistance(core)
    build["tec_coupled_enabled"] = bool(core.enable_tec_coupled)
    return build


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0:
        raise ValueError("duration must be positive.")
    if args.record_interval <= 0.0:
        raise ValueError("record-interval must be positive.")
    if args.restart_interval <= 0.0:
        raise ValueError("restart-interval must be positive.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_history.csv"
    latest_state_path = output_dir / f"{args.case_prefix}_latest_state.json"
    latest_restart_path = output_dir / f"{args.case_prefix}_latest_restart.npz"

    build = build_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)
    next_record_time = min(start_time + float(args.record_interval), target_time)
    next_restart_time = min(start_time + float(args.restart_interval), target_time)
    fieldnames: Optional[List[str]] = None

    print("=== V9 CaseA open external piping ===", flush=True)
    print(f"case_version={V9_CASE_VERSION}", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)
    print(f"latest_restart={latest_restart_path}", flush=True)
    print(f"coolant_material={build['coolant_material']}", flush=True)
    print(f"wire_resistance_ohm={build['wire_resistance_ohm']}", flush=True)

    last_record: Dict[str, Any] = {}
    while float(system.global_time) < target_time - 1.0e-10:
        build["core"].update_neutronic_power(
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
        dt = min(float(dt), float(args.max_dt), next_record_time - float(system.global_time), target_time - float(system.global_time))
        system.step(dt, inner_iter=int(args.inner_iter))

        if float(system.global_time) >= next_record_time - 1.0e-10:
            passive = passive_tec_source_totals(build)
            if any(value != 0.0 for value in passive.values()):
                raise RuntimeError(f"Ring3_Open TEC sources are not zero: {passive}")
            record = {
                **v9_basic_diagnostics(build),
                "relative_time_s": float(system.global_time) - start_time,
                "max_dt_s": float(args.max_dt),
                "inner_iter": int(args.inner_iter),
                "solid_ode_method": build["solid_ode_method"],
                "wire_resistance_ohm": build["wire_resistance_ohm"],
                "tec_coupled_enabled": build["tec_coupled_enabled"],
            }
            flat = flatten_for_csv(record)
            if fieldnames is None:
                fieldnames = list(flat.keys())
            append_row(history_path, fieldnames, flat, write_header=not history_path.exists())
            last_record = record
            print(
                f"t_rel={record['relative_time_s']:.1f}s "
                f"Tin={record['core_inlet_connector_t_k']:.3f}K "
                f"Tout={record['core_outlet_connector_t_k']:.3f}K "
                f"Qcool={record['coolant_enthalpy_rise_w']:.3f}W "
                f"Pel={float(record['tec_total_electric_power_w'] or 0.0):.3f}W",
                flush=True,
            )
            next_record_time = min(next_record_time + float(args.record_interval), target_time)

        if float(system.global_time) >= next_restart_time - 1.0e-10:
            system.save_global_state(str(latest_restart_path))
            latest = {
                "case_version": V9_CASE_VERSION,
                "restart_in": args.restart_in,
                "restart_out": str(latest_restart_path),
                "history_csv": str(history_path),
                "start_time_s": start_time,
                "end_time_s": float(system.global_time),
                "target_time_s": target_time,
                "record_interval_s": float(args.record_interval),
                "restart_interval_s": float(args.restart_interval),
                "max_dt_s": float(args.max_dt),
                "inner_iter": int(args.inner_iter),
                "coolant_material": build["coolant_material"],
                "solid_ode_method": build["solid_ode_method"],
                "solid_ode_methods": get_solid_ode_methods(build),
                "wire_resistance_ohm": build["wire_resistance_ohm"],
                "tec_coupled_enabled": build["tec_coupled_enabled"],
                "ring_multipliers": build["ring_multipliers"],
                "tec_ring_multipliers": build["tec_ring_multipliers"],
                "latest": last_record or v9_basic_diagnostics(build),
            }
            with latest_state_path.open("w", encoding="utf-8") as f:
                json.dump(latest, f, indent=2, ensure_ascii=False, default=json_default)
            next_restart_time = min(next_restart_time + float(args.restart_interval), target_time)

    if not latest_restart_path.exists():
        system.save_global_state(str(latest_restart_path))
    print("V9 open-loop run completed.", flush=True)


if __name__ == "__main__":
    main()
