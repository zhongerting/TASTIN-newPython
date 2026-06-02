import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from audit_v7_caseA_fluid_energy_network import build_volume_scale_map, classify_volume
from run_v7_caseA_multipliers_short import build_loaded_case, json_default, parse_multipliers
from test_core_assemble_v7_caseA import _case_a_electric_diagnostics


DEFAULT_RESTART = (
    "testModule/v7_caseA_odecopy_fvm_continue100/"
    "v7_caseA_odecopy_fvm_continue100_latest_restart.npz"
)
TOTAL_POWER_W = 115000.0
SIGMA_W_M2_K4 = 5.670374419e-8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Continue V7 CaseA without interruption and audit global finite-difference "
            "solid/fluid storage against external energy flows."
        )
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v7_caseA_global_storage_continue")
    parser.add_argument("--case-prefix", default="v7_caseA_global_storage_continue")
    parser.add_argument("--duration", type=float, default=100.0)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,18"))
    parser.add_argument("--tec-ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,15"))
    return parser.parse_args()


def solid_multiplier(build: Dict[str, Any], solid_name: str) -> float:
    for ring_name, multiplier in build["ring_multipliers"].items():
        if solid_name.startswith(f"{ring_name}_"):
            return float(multiplier)
    return 1.0


def solid_category(build: Dict[str, Any], solid_name: str) -> str:
    for ring_name in build["ring_multipliers"]:
        if solid_name.startswith(f"{ring_name}_"):
            return "tfe"
    return "global_outer"


def capture_storage_state(build: Dict[str, Any]) -> Dict[str, Any]:
    solids = {}
    for name, solid in build["system"].solid_components.items():
        if hasattr(solid, "_update_properties"):
            solid._update_properties()
        solids[name] = {
            "T": np.asarray(solid.T, dtype=float).copy(),
            "cap": np.asarray(solid.thermal_capacitance, dtype=float).copy(),
        }

    fluids = {}
    scale_map = build_volume_scale_map(build)
    for idx, vol in enumerate(build["system"].fluid_solver.volumes_obj):
        name = getattr(vol, "name", f"vol_{idx}")
        fluids[name] = {
            "h": float(getattr(vol, "h", np.nan)),
            "mass": float(getattr(vol, "rho", 0.0) * getattr(vol, "vol", 0.0)),
            "scale": float(scale_map.get(id(vol), 1.0)),
            "is_fixed": bool(getattr(vol, "is_pressure_boundary", False)),
            "type": classify_volume(vol),
        }
    return {"solids": solids, "fluids": fluids}


def finite_difference_storage_rates(
    build: Dict[str, Any],
    before: Dict[str, Any],
    after: Dict[str, Any],
    dt: float,
) -> Dict[str, float]:
    rates = {
        "solid_storage_w": 0.0,
        "solid_tfe_storage_w": 0.0,
        "solid_global_outer_storage_w": 0.0,
        "fluid_storage_w": 0.0,
        "fluid_storage_exact_delta_mh_w": 0.0,
        "fluid_plenum_storage_w": 0.0,
        "fluid_pipe_storage_w": 0.0,
        "fluid_tfe_channel_storage_w": 0.0,
    }

    for name, old in before["solids"].items():
        new = after["solids"][name]
        multiplier = solid_multiplier(build, name)
        storage = float(
            np.sum(0.5 * (old["cap"] + new["cap"]) * (new["T"] - old["T"]))
            / dt
            * multiplier
        )
        rates["solid_storage_w"] += storage
        rates[f"solid_{solid_category(build, name)}_storage_w"] += storage

    for name, old in before["fluids"].items():
        new = after["fluids"][name]
        if old["is_fixed"]:
            continue
        scale = float(old["scale"])
        storage = (
            0.5
            * (float(old["mass"]) + float(new["mass"]))
            * (float(new["h"]) - float(old["h"]))
            / dt
            * scale
        )
        exact_delta_mh = (
            (float(new["mass"]) * float(new["h"]) - float(old["mass"]) * float(old["h"]))
            / dt
            * scale
        )
        rates["fluid_storage_w"] += storage
        rates["fluid_storage_exact_delta_mh_w"] += exact_delta_mh
        vol_type = str(old["type"])
        if vol_type == "plenum":
            rates["fluid_plenum_storage_w"] += storage
        elif vol_type == "pipe":
            rates["fluid_pipe_storage_w"] += storage
        elif vol_type == "tfe_channel":
            rates["fluid_tfe_channel_storage_w"] += storage

    rates["combined_storage_w"] = rates["solid_storage_w"] + rates["fluid_storage_w"]
    return rates


def external_fluid_enthalpy_net_out_w(build: Dict[str, Any]) -> float:
    net = build["system"].fluid_solver
    scale_map = build_volume_scale_map(build)
    net_out = 0.0
    for j_idx, (_, idx_from, idx_to) in enumerate(net.junction_descriptors):
        from_vol = net.volumes_obj[idx_from]
        to_vol = net.volumes_obj[idx_to]
        from_fixed = bool(getattr(from_vol, "is_pressure_boundary", False))
        to_fixed = bool(getattr(to_vol, "is_pressure_boundary", False))
        if from_fixed == to_fixed:
            continue

        w = float(net.W_vec[j_idx])
        donor_idx = idx_from if w >= 0.0 else idx_to
        flux_from_to = w * float(getattr(net.volumes_obj[donor_idx], "h", np.nan))
        finite_vol = to_vol if from_fixed else from_vol
        scale = float(scale_map.get(id(finite_vol), 1.0))
        finite_row_advection = (
            float(net.M_to_vec[j_idx]) * flux_from_to
            if from_fixed
            else -float(net.M_from_vec[j_idx]) * flux_from_to
        )
        net_out -= scale * finite_row_advection
    return net_out


def fluid_effective_source_w(build: Dict[str, Any]) -> float:
    scale_map = build_volume_scale_map(build)
    total = 0.0
    for vol in build["system"].fluid_solver.volumes_obj:
        scale = float(scale_map.get(id(vol), 1.0))
        total += scale * (
            float(getattr(vol, "Q_wall", 0.0))
            + float(getattr(vol, "Q_vol", 0.0))
            - float(getattr(vol, "implicit_coeff", 0.0)) * float(getattr(vol, "T", 0.0))
        )
    return total


def coolant_solid_to_fluid_w(build: Dict[str, Any]) -> float:
    total = 0.0
    for name, multiplier in build["ring_multipliers"].items():
        tfe = build["tfes"][name]
        for coupler_name in ("iclad_coolant_fsc", "oclad_coolant_fsc"):
            boundary = tfe.couplers[coupler_name].solid_bound
            total += -float(np.sum(np.asarray(boundary.current_flux, dtype=float))) * float(multiplier)
    return total


def outer_wall_radiation_w(build: Dict[str, Any]) -> float:
    boundary = build["core"].reflector.boundaries["right"]
    surface_temperature = np.asarray(boundary.T_surface, dtype=float)
    area = np.asarray(boundary.area, dtype=float)
    return float(np.sum(0.2 * SIGMA_W_M2_K4 * area * (surface_temperature**4 - 200.0**4)))


def moderator_mapping_source_minus_boundary_out_w(build: Dict[str, Any]) -> float:
    core = build["core"]
    if not getattr(core, "has_global_moderator", False):
        return 0.0
    actual_ring_source = float(sum(np.sum(np.asarray(ring.Q_source, dtype=float)) for ring in core.mod_rings))
    expected_from_tfe_boundary = 0.0
    for name, tfe in build["tfes"].items():
        boundary = tfe.solids["moderator"].boundaries["right"]
        expected_from_tfe_boundary += (
            -float(np.sum(np.asarray(boundary.current_flux, dtype=float)))
            * float(build["ring_multipliers"][name])
        )
    return actual_ring_source - expected_from_tfe_boundary


def applied_tec_heat_removed_w(build: Dict[str, Any], multiplier_key: str) -> float:
    total = 0.0
    for name, multiplier in build[multiplier_key].items():
        if int(build["tec_ring_multipliers"].get(name, 0)) <= 0:
            continue
        tfe = build["tfes"][name]
        area = np.asarray(tfe.solids["emitter"].boundaries["right"].area, dtype=float)
        emitter_removed = -float(np.sum(np.asarray(tfe.plasma_data.electron_cooling_flux) * area))
        collector_source = float(np.sum(np.asarray(tfe.plasma_data.electron_heating_flux) * area))
        joule = float(np.sum(tfe.electric_data.emitter_joule_heat)) + float(
            np.sum(tfe.electric_data.collector_joule_heat)
        )
        total += (emitter_removed - collector_source - joule) * float(multiplier)
    return total


def external_power_snapshot(build: Dict[str, Any]) -> Dict[str, float]:
    core = build["core"]
    electric = _case_a_electric_diagnostics(core)
    q_core = sum(
        float(tfe.neutronic_data.total_power) * float(core.tfe_multipliers[name])
        for name, tfe in build["tfes"].items()
    )
    return {
        "core_heat_w": q_core,
        "fluid_external_net_out_w": external_fluid_enthalpy_net_out_w(build),
        "fluid_effective_source_w": fluid_effective_source_w(build),
        "coolant_solid_to_fluid_w": coolant_solid_to_fluid_w(build),
        "outer_wall_radiation_w": outer_wall_radiation_w(build),
        "moderator_mapping_source_minus_boundary_out_w": moderator_mapping_source_minus_boundary_out_w(build),
        "terminal_electric_power_w": float(electric["tec_total_electric_power_w"] or 0.0),
        "applied_tec_heat_removed_electric_count_w": applied_tec_heat_removed_w(
            build,
            "tec_ring_multipliers",
        ),
        "applied_tec_heat_removed_thermal_model_w": applied_tec_heat_removed_w(
            build,
            "ring_multipliers",
        ),
    }


def add_trapezoid_integral(
    integral: Dict[str, float],
    before: Dict[str, float],
    after: Dict[str, float],
    dt: float,
) -> None:
    for key in before:
        integral[key] += 0.5 * (float(before[key]) + float(after[key])) * dt


def advance_interval(
    build: Dict[str, Any],
    stop_time: float,
    args: argparse.Namespace,
    initial_power: Dict[str, float],
) -> Dict[str, Any]:
    system = build["system"]
    core = build["core"]
    power = dict(initial_power)
    integral = {key: 0.0 for key in power}
    step_count = 0
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
        system.step(dt, inner_iter=int(args.inner_iter))
        next_power = external_power_snapshot(build)
        add_trapezoid_integral(integral, power, next_power, dt)
        power = next_power
        step_count += 1
    return {"integral": integral, "power": power, "step_count": step_count}


def build_record(
    build: Dict[str, Any],
    start_time: float,
    interval_start: float,
    interval_end: float,
    integral: Dict[str, float],
    storage: Dict[str, float],
    step_count: int,
) -> Dict[str, float]:
    dt = interval_end - interval_start
    avg = {f"{key}_avg": float(value) / dt for key, value in integral.items()}
    q_core = avg["core_heat_w_avg"]
    q_fluid = avg["fluid_external_net_out_w_avg"]
    q_fluid_source = avg["fluid_effective_source_w_avg"]
    q_solid_to_fluid = avg["coolant_solid_to_fluid_w_avg"]
    q_rad = avg["outer_wall_radiation_w_avg"]
    p_terminal = avg["terminal_electric_power_w_avg"]
    q_tec_electric_count = avg["applied_tec_heat_removed_electric_count_w_avg"]
    q_tec_thermal_model = avg["applied_tec_heat_removed_thermal_model_w_avg"]
    combined_storage = storage["combined_storage_w"]
    model_residual = q_core - q_fluid - q_rad - q_tec_thermal_model - combined_storage
    physical_tec_residual = q_core - q_fluid - q_rad - q_tec_electric_count - combined_storage
    terminal_residual = q_core - q_fluid - q_rad - p_terminal - combined_storage
    solid_equation_residual = q_core - q_rad - q_tec_thermal_model - q_solid_to_fluid - storage["solid_storage_w"]
    fluid_equation_residual = q_fluid_source - q_fluid - storage["fluid_storage_w"]
    fluid_solid_mapping_residual = q_solid_to_fluid - q_fluid_source
    return {
        "absolute_time_s": interval_end,
        "relative_time_s": interval_end - start_time,
        "interval_start_s": interval_start,
        "interval_dt_s": dt,
        "internal_step_count": step_count,
        **avg,
        **storage,
        "tec_thermal_model_minus_electric_count_w": q_tec_thermal_model - q_tec_electric_count,
        "tec_terminal_minus_electric_count_heat_w": p_terminal - q_tec_electric_count,
        "solid_equation_residual_w": solid_equation_residual,
        "fluid_equation_residual_w": fluid_equation_residual,
        "fluid_solid_mapping_residual_w": fluid_solid_mapping_residual,
        "decomposed_thermal_model_residual_w": (
            solid_equation_residual
            + fluid_equation_residual
            + fluid_solid_mapping_residual
        ),
        "fast_residual_without_storage_using_terminal_w": q_core - q_fluid - q_rad - p_terminal,
        "global_residual_using_thermal_model_tec_heat_w": model_residual,
        "global_residual_using_electric_count_tec_heat_w": physical_tec_residual,
        "global_residual_using_terminal_power_w": terminal_residual,
        "global_relative_residual_using_thermal_model_tec_heat": model_residual / q_core,
        "global_relative_residual_using_electric_count_tec_heat": physical_tec_residual / q_core,
        "global_relative_residual_using_terminal_power": terminal_residual / q_core,
    }


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0:
        raise ValueError("duration must be positive.")
    if args.record_interval <= 0.0:
        raise ValueError("record-interval must be positive.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_global_storage_history.csv"
    summary_path = output_dir / f"{args.case_prefix}_global_storage_summary.json"
    restart_path = output_dir / f"{args.case_prefix}_latest_restart.npz"

    build = build_loaded_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)
    interval_start = start_time
    storage_before = capture_storage_state(build)
    power = external_power_snapshot(build)
    rows: List[Dict[str, float]] = []

    print("=== V7 CaseA global storage continuation audit ===", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)

    while interval_start < target_time - 1.0e-10:
        interval_end = min(interval_start + float(args.record_interval), target_time)
        advanced = advance_interval(build, interval_end, args, power)
        power = advanced["power"]
        storage_after = capture_storage_state(build)
        storage = finite_difference_storage_rates(
            build,
            storage_before,
            storage_after,
            interval_end - interval_start,
        )
        row = build_record(
            build,
            start_time,
            interval_start,
            interval_end,
            advanced["integral"],
            storage,
            advanced["step_count"],
        )
        rows.append(row)
        print(
            f"t_rel={row['relative_time_s']:.1f}s "
            f"Rfast={row['fast_residual_without_storage_using_terminal_w']:.3f}W "
            f"Ssolid={row['solid_storage_w']:.3f}W "
            f"Sfluid={row['fluid_storage_w']:.3f}W "
            f"Rmodel={row['global_residual_using_thermal_model_tec_heat_w']:.3f}W "
            f"Rtec={row['global_residual_using_electric_count_tec_heat_w']:.3f}W "
            f"Rterminal={row['global_residual_using_terminal_power_w']:.3f}W",
            flush=True,
        )
        interval_start = interval_end
        storage_before = storage_after

    with history_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    system.save_global_state(str(restart_path))
    summary = {
        "restart_in": args.restart_in,
        "restart_out": str(restart_path),
        "history_csv": str(history_path),
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "record_interval_s": float(args.record_interval),
        "max_dt_s": float(args.max_dt),
        "inner_iter": int(args.inner_iter),
        "ring_multipliers": build["ring_multipliers"],
        "tec_ring_multipliers": build["tec_ring_multipliers"],
        "latest": rows[-1],
        "notes": {
            "solid_storage_definition": "sum(0.5 * (cap_old + cap_new) * (T_new - T_old)) / interval_dt",
            "fluid_storage_definition": "sum(0.5 * (mass_old + mass_new) * (h_new - h_old)) / interval_dt over finite non-fixed volumes",
            "fluid_boundary_definition": "net advective enthalpy outflow across fixed-pressure boundaries using HydraulicNetwork donor and multiplier rows",
            "fluid_effective_source_definition": "sum(scale * (Q_wall + Q_vol - implicit_coeff * T)) over hydraulic volumes",
            "coolant_solid_to_fluid_definition": "negative inner/outer clad coolant BoundaryRegion.current_flux scaled by thermal multipliers",
            "external_power_integration": "trapezoidal integration over every internal SystemManager step",
            "moderator_mapping_diagnostic": "instantaneous global moderator ring source minus thermal-multiplier-scaled representative TFE moderator outer-boundary outflow",
            "thermal_model_residual_definition": "Qcore - Qfluid_external_net_out - Qradiation - Qtec_applied_removed_scaled_by_thermal_multiplier - dEsolid_dt - dEfluid_dt",
            "thermal_model_residual_decomposition": "solid_equation_residual + fluid_equation_residual + fluid_solid_mapping_residual",
            "electric_count_tec_residual_definition": "Qcore - Qfluid_external_net_out - Qradiation - Qtec_applied_removed_scaled_by_tec_multiplier - dEsolid_dt - dEfluid_dt",
            "terminal_residual_definition": "Qcore - Qfluid_external_net_out - Qradiation - Pterminal - dEsolid_dt - dEfluid_dt",
        },
    }
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)

    print(f"CSV: {history_path}", flush=True)
    print(f"Summary: {summary_path}", flush=True)
    print(f"Restart: {restart_path}", flush=True)


if __name__ == "__main__":
    main()
