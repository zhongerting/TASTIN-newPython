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

from test_core_assemble_v7_caseA import (
    CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    build_v7_case_a_system,
    _case_a_electric_diagnostics,
    _case_a_flow_diagnostics,
    _case_a_reset_design_flows_after_restart,
)
from test_core_assemble_v7_caseA_faststeady import compute_faststeady_energy_audit


DEFAULT_RESTART = (
    "testModule/v7_caseA_geometry_fix_continue2000/"
    "test_core_assemble_v7_caseA_geometry_fix_continue2000_restart_rel2000_abs36840.npz"
)
TOTAL_POWER_W = 115000.0


def parse_multipliers(text: str) -> List[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if len(values) != 4:
        raise argparse.ArgumentTypeError("Use four comma-separated integers, e.g. 1,6,12,15.")
    if any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("All multipliers must be positive.")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Short V7 CaseA run with custom ring multipliers.")
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v7_caseA_multiplier_short")
    parser.add_argument("--case-prefix", default="test_core_assemble_v7_caseA_multiplier_short")
    parser.add_argument("--ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,18"))
    parser.add_argument("--tec-ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,15"))
    parser.add_argument("--duration", type=float, default=40.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    return parser.parse_args()


def json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return str(value)


def axial_node_volumes(solid, n_nodes: int) -> np.ndarray:
    vols = np.asarray(solid.vols_flat, dtype=float).reshape(solid.nx, n_nodes)
    return np.sum(vols, axis=0)


def collect_tec_stats(build: Dict[str, Any]) -> Dict[str, Any]:
    core = build["core"]
    ring_multipliers = build["ring_multipliers"]
    tec_ring_multipliers = build["tec_ring_multipliers"]
    thermo_calc = core.thermo_calc
    if thermo_calc is None:
        return {}

    input_data = thermo_calc._input_data
    idx = 0
    per_ring = {}
    totals = {
        "node_current_sum_a": 0.0,
        "emitter_electron_heat_source_signed_w": 0.0,
        "emitter_electron_heat_removed_positive_w": 0.0,
        "collector_electron_heat_source_signed_w": 0.0,
        "collector_electron_heat_source_if_sideAreaC_w": 0.0,
        "collector_sideAreaC_minus_emitter_basis_w": 0.0,
        "electron_boundary_heat_difference_w": 0.0,
        "emitter_joule_heat_w": 0.0,
        "collector_joule_heat_w": 0.0,
        "total_joule_heat_w": 0.0,
    }
    for name, tec_multiplier in tec_ring_multipliers.items():
        tec_multiplier = int(tec_multiplier)
        if tec_multiplier <= 0:
            per_ring[name] = {
                "thermal_multiplier": int(ring_multipliers[name]),
                "tec_multiplier": tec_multiplier,
                "node_current_sum_a": 0.0,
                "j_nonzero_count": 0,
            }
            continue
        res = thermo_calc.get_tec_results(idx)
        j = np.asarray(res["J"], dtype=float)
        j_a_m2 = j * 1.0e4
        nonzero = np.abs(j) > 1.0e-12
        active_j = j[nonzero]
        if active_j.size == 0:
            active_j = j
        area = np.asarray(input_data.sideAreaE[idx], dtype=float)
        area_c = np.asarray(input_data.sideAreaC[idx], dtype=float)
        node_current = j_a_m2 * area
        te = np.asarray(res["TE"], dtype=float)
        phi_e = np.asarray(res["phiE"], dtype=float)
        ue = np.asarray(res["UE"], dtype=float)
        uc = np.asarray(res["UC"], dtype=float)

        q_e_flux = -j_a_m2 * (phi_e + 2.0 * 8.617e-5 * te)
        q_c_flux = j_a_m2 * (phi_e + 2.0 * 8.617e-5 * te - (ue - uc))
        emitter_electron_source = q_e_flux * area
        emitter_electron_removed = -emitter_electron_source
        collector_electron_source = q_c_flux * area
        collector_electron_source_if_sideAreaC = q_c_flux * area_c
        boundary_difference = emitter_electron_removed - collector_electron_source

        emitter_joule = np.asarray(res["joulePowerE"], dtype=float)
        collector_joule = np.asarray(res["joulePowerC"], dtype=float)
        total_joule = emitter_joule + collector_joule

        single_energy = {
            "node_current_sum_a": float(np.sum(node_current)),
            "emitter_electron_heat_source_signed_w": float(np.sum(emitter_electron_source)),
            "emitter_electron_heat_removed_positive_w": float(np.sum(emitter_electron_removed)),
            "collector_electron_heat_source_signed_w": float(np.sum(collector_electron_source)),
            "collector_electron_heat_source_if_sideAreaC_w": float(np.sum(collector_electron_source_if_sideAreaC)),
            "collector_sideAreaC_minus_emitter_basis_w": float(
                np.sum(collector_electron_source_if_sideAreaC - collector_electron_source)
            ),
            "electron_boundary_heat_difference_w": float(np.sum(boundary_difference)),
            "emitter_joule_heat_w": float(np.sum(emitter_joule)),
            "collector_joule_heat_w": float(np.sum(collector_joule)),
            "total_joule_heat_w": float(np.sum(total_joule)),
        }
        group_energy = {
            key: value * float(tec_multiplier)
            for key, value in single_energy.items()
        }
        for key, value in group_energy.items():
            if key not in totals:
                totals[key] = 0.0
            totals[key] += value

        per_ring[name] = {
            "thermal_multiplier": int(ring_multipliers[name]),
            "tec_multiplier": tec_multiplier,
            "sideAreaE_min_m2": float(np.min(area)),
            "sideAreaE_mean_m2": float(np.mean(area)),
            "sideAreaE_max_m2": float(np.max(area)),
            "sideAreaE_sum_m2": float(np.sum(area)),
            "sideAreaC_min_m2": float(np.min(area_c)),
            "sideAreaC_mean_m2": float(np.mean(area_c)),
            "sideAreaC_max_m2": float(np.max(area_c)),
            "sideAreaC_sum_m2": float(np.sum(area_c)),
            "node_current_sum_a": float(np.sum(node_current)),
            "j_nonzero_count": int(np.count_nonzero(nonzero)),
            "j_mean_nonzero_a_cm2": float(np.mean(active_j)),
            "j_max_a_cm2": float(np.max(j)),
            "j_min_a_cm2": float(np.min(j)),
            "emitter_temperature_mean_k": float(np.mean(res["TE"])),
            "emitter_temperature_max_k": float(np.max(res["TE"])),
            "collector_temperature_mean_k": float(np.mean(res["TC"])),
            "collector_temperature_max_k": float(np.max(res["TC"])),
            "single_generator_energy": single_energy,
            "tec_group_energy": group_energy,
        }
        idx += tec_multiplier
    return {
        "by_ring": per_ring,
        "totals": totals,
    }


def write_tec_distribution_csv(build: Dict[str, Any], csv_path: Path) -> None:
    core = build["core"]
    thermo_calc = core.thermo_calc
    if thermo_calc is None:
        return

    ring_multipliers = build["ring_multipliers"]
    tec_ring_multipliers = build["tec_ring_multipliers"]
    reference_tfe = next(iter(build["tfes"].values()))
    y_faces = np.asarray(reference_tfe.common_y_faces, dtype=float)
    y_centers = 0.5 * (y_faces[:-1] + y_faces[1:])
    n_lower = int(getattr(reference_tfe, "n_lower", 0))
    n_active = int(getattr(reference_tfe, "n_active", len(y_centers)))
    active_start = n_lower
    active_stop = n_lower + n_active

    input_data = thermo_calc._input_data
    fieldnames = [
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

    idx = 0
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
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
                node_current = float(j[node] * 1.0e4 * area[node])
                writer.writerow({
                    "representative_tfe": name,
                    "thermal_multiplier": int(ring_multipliers[name]),
                    "tec_multiplier": tec_multiplier,
                    "axial_node": node,
                    "y_center_m": float(y_centers[node]),
                    "is_active_node": active_start <= node < active_stop,
                    "T_emitter_k": float(te[node]),
                    "T_collector_k": float(tc[node]),
                    "J_A_cm2": float(j[node]),
                    "node_current_A": node_current,
                    "UE_minus_UC_v": float(ue[node] - uc[node]),
                    "Vd_v": float(vd[node]),
                })
            idx += tec_multiplier


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
    if not os.path.exists(args.restart_in):
        raise FileNotFoundError(f"Restart file not found: {args.restart_in}")
    system.load_global_state(args.restart_in)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=260.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)
    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is not None:
        core.thermo_calc.calculate(verbose=False)
    return build


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    build = build_loaded_case(args)
    system = build["system"]
    core = build["core"]
    start_time = float(system.global_time)
    stop_time = start_time + float(args.duration)

    history = []
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
        dt = min(float(dt), stop_time - float(system.global_time), float(args.max_dt))
        system.step(dt)

        electric = _case_a_electric_diagnostics(core)
        history.append({
            "time_s": float(system.global_time),
            "relative_time_s": float(system.global_time) - start_time,
            "dt_s": dt,
            "current_a": electric["tec_total_current_a"],
            "voltage_v": electric["tec_total_voltage_v"],
            "power_w": electric["tec_total_electric_power_w"],
        })

    electric = _case_a_electric_diagnostics(core)
    flow = _case_a_flow_diagnostics(build)
    energy = compute_faststeady_energy_audit({
        "build": build,
        "system": system,
        "core": core,
    })
    tec_stats = collect_tec_stats(build)

    summary = {
        "restart_in": args.restart_in,
        "ring_multipliers": args.ring_multipliers,
        "tec_ring_multipliers": args.tec_ring_multipliers,
        "total_thermal_elements": int(sum(args.ring_multipliers)),
        "total_virtual_elements": int(sum(args.tec_ring_multipliers)),
        "total_design_flow_kg_s": CASE_A_DESIGN_TOTAL_FLOW_KG_S,
        "single_tfe_design_flow_kg_s": CASE_A_DESIGN_TOTAL_FLOW_KG_S / float(sum(args.ring_multipliers)),
        "target_voltage_v": args.target_voltage,
        "start_time_s": start_time,
        "final_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "electric": electric,
        "flow": flow,
        "energy": energy,
        "tec_stats": tec_stats,
        "history_last": history[-1] if history else None,
    }

    history_path = output_dir / f"{args.case_prefix}_history.csv"
    with history_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()) if history else ["time_s"])
        writer.writeheader()
        writer.writerows(history)

    summary_path = output_dir / f"{args.case_prefix}_summary.json"
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=json_default)
    distribution_path = output_dir / f"{args.case_prefix}_tec_distribution.csv"
    write_tec_distribution_csv(build, distribution_path)

    print(f"History: {history_path}")
    print(f"Summary: {summary_path}")
    print(f"TEC distribution: {distribution_path}")
    print(json.dumps({
        "ring_multipliers": args.ring_multipliers,
        "tec_ring_multipliers": args.tec_ring_multipliers,
        "total_thermal_elements": int(sum(args.ring_multipliers)),
        "total_virtual_elements": int(sum(args.tec_ring_multipliers)),
        "current_a": electric["tec_total_current_a"],
        "voltage_v": electric["tec_total_voltage_v"],
        "power_w": electric["tec_total_electric_power_w"],
        "single_tfe_design_flow_kg_s": summary["single_tfe_design_flow_kg_s"],
        "tec_stats": tec_stats,
    }, indent=2, default=json_default))


if __name__ == "__main__":
    main()
