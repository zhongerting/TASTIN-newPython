import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

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
    _case_a_reset_design_flows_after_restart,
)
from Components.tec_electric import electric_field_from_node_potential


DEFAULT_RESTART = (
    "testModule/v7_caseA_electric_dt08_physicalCp_overnight/"
    "test_core_assemble_v7_caseA_electric_dt08_physicalCp_restart_rel26800_abs34600.npz"
)
DEFAULT_OUTPUT_DIR = "testModule/v7_caseA_electric_dt08_physicalCp_overnight"
TOTAL_POWER_W = 115000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract per-virtual-element per-node TEC voltage/current-density data."
    )
    parser.add_argument("--restart", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--prefix", default=None)
    return parser.parse_args()


def json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return str(value)


def build_and_load(restart: str) -> Dict[str, Any]:
    build = build_v7_case_a_system(
        pipe_n_nodes=8,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
    )
    system = build["system"]
    core = build["core"]
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
    system.load_global_state(restart)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)

    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is None:
        raise RuntimeError("TEC calculation is not available.")
    core.thermo_calc.calculate(verbose=False)
    return build


def virtual_element_mapping(core) -> List[Tuple[int, str, int]]:
    mapping: List[Tuple[int, str, int]] = []
    idx = 0
    for tfe_name, multiplier in getattr(core, "tec_multipliers", core.tfe_multipliers).items():
        for local_index in range(int(multiplier)):
            mapping.append((idx, tfe_name, local_index))
            idx += 1
    return mapping


def axial_node_volumes(solid, n_nodes: int) -> np.ndarray:
    vols = np.asarray(solid.vols_flat, dtype=float).reshape(solid.nx, n_nodes)
    return np.sum(vols, axis=0)


def extract_node_rows(build: Dict[str, Any]) -> Tuple[List[Dict[str, float]], Dict[str, Any]]:
    core = build["core"]
    rows: List[Dict[str, float]] = []
    totals = {
        "node_power_gap_voltage_w": 0.0,
        "node_power_ue_minus_uc_w": 0.0,
        "node_power_barrier_vd_w": 0.0,
        "node_current_sum_a": 0.0,
        "node_electron_net_source_w": 0.0,
        "emitter_electron_heat_source_signed_w": 0.0,
        "emitter_electron_heat_removed_positive_w": 0.0,
        "collector_electron_heat_source_signed_w": 0.0,
        "electron_boundary_heat_difference_w": 0.0,
        "emitter_joule_heat_w": 0.0,
        "collector_joule_heat_w": 0.0,
        "total_joule_heat_w": 0.0,
        "applied_emitter_electron_heat_source_signed_w": 0.0,
        "applied_emitter_electron_heat_removed_positive_w": 0.0,
        "applied_collector_electron_heat_source_signed_w": 0.0,
        "applied_electron_boundary_heat_difference_w": 0.0,
        "applied_emitter_joule_heat_w": 0.0,
        "applied_collector_joule_heat_w": 0.0,
        "applied_total_joule_heat_w": 0.0,
    }
    by_representative: Dict[str, Dict[str, float]] = {}

    for virtual_idx, tfe_name, local_index in virtual_element_mapping(core):
        tfe = build["tfes"][tfe_name]
        result = core.thermo_calc.get_tec_results(virtual_idx)
        if result is None:
            continue

        emitter = tfe.solids["emitter"]
        collector = tfe.solids["collector"]
        area = np.asarray(emitter.boundaries["right"].area, dtype=float)
        emitter_node_volumes = axial_node_volumes(emitter, len(area))
        collector_node_volumes = axial_node_volumes(collector, len(area))

        j_a_cm2 = np.asarray(result["J"], dtype=float)
        j_a_m2 = j_a_cm2 * 1.0e4
        v_gap = np.asarray(result["V"], dtype=float)
        ue = np.asarray(result["UE"], dtype=float)
        uc = np.asarray(result["UC"], dtype=float)
        vd = np.asarray(result["Vd"], dtype=float)
        phi_e = np.asarray(result["phiE"], dtype=float)
        phi_c = np.asarray(result["phiC"], dtype=float)
        rho_e = np.asarray(result["rhoE"], dtype=float)
        rho_c = np.asarray(result["rhoC"], dtype=float)
        te = np.asarray(result["TE"], dtype=float)

        y_faces = np.asarray(getattr(tfe, "common_y_faces", emitter.mesh.y_faces), dtype=float)
        e_field = electric_field_from_node_potential(ue, y_faces=y_faces)
        c_field = electric_field_from_node_potential(uc, y_faces=y_faces)
        cpp_joule_e = np.asarray(result["joulePowerE"], dtype=float)
        cpp_joule_c = np.asarray(result["joulePowerC"], dtype=float)

        q_electron_e_flux = -1.0 * j_a_m2 * (phi_e + 2.0 * 8.617e-5 * te)
        q_electron_c_flux = 1.0 * j_a_m2 * (
            phi_e + 2.0 * 8.617e-5 * te - (ue - uc)
        )

        applied_q_e_flux = np.asarray(tfe.plasma_data.electron_cooling_flux, dtype=float)
        applied_q_c_flux = np.asarray(tfe.plasma_data.electron_heating_flux, dtype=float)
        applied_joule_e_flat = np.asarray(tfe.electric_data.emitter_joule_heat, dtype=float)
        applied_joule_c_flat = np.asarray(tfe.electric_data.collector_joule_heat, dtype=float)
        applied_joule_e_by_node = applied_joule_e_flat.reshape(emitter.nx, len(area)).sum(axis=0)
        applied_joule_c_by_node = applied_joule_c_flat.reshape(collector.nx, len(area)).sum(axis=0)

        n = min(
            len(area),
            len(emitter_node_volumes),
            len(collector_node_volumes),
            len(j_a_cm2),
            len(v_gap),
            len(ue),
            len(uc),
            len(vd),
            len(phi_e),
            len(phi_c),
            len(rho_e),
            len(rho_c),
            len(te),
            len(e_field),
            len(c_field),
            len(cpp_joule_e),
            len(cpp_joule_c),
            len(q_electron_e_flux),
            len(q_electron_c_flux),
            len(applied_q_e_flux),
            len(applied_q_c_flux),
            len(applied_joule_e_by_node),
            len(applied_joule_c_by_node),
        )
        rep_totals = by_representative.setdefault(
            tfe_name,
            {
                "virtual_element_count": 0.0,
                "node_power_gap_voltage_w": 0.0,
                "node_power_ue_minus_uc_w": 0.0,
                "node_power_barrier_vd_w": 0.0,
                "node_current_sum_a": 0.0,
                "node_electron_net_source_w": 0.0,
                "emitter_electron_heat_source_signed_w": 0.0,
                "emitter_electron_heat_removed_positive_w": 0.0,
                "collector_electron_heat_source_signed_w": 0.0,
                "electron_boundary_heat_difference_w": 0.0,
                "emitter_joule_heat_w": 0.0,
                "collector_joule_heat_w": 0.0,
                "total_joule_heat_w": 0.0,
                "applied_emitter_electron_heat_source_signed_w": 0.0,
                "applied_emitter_electron_heat_removed_positive_w": 0.0,
                "applied_collector_electron_heat_source_signed_w": 0.0,
                "applied_electron_boundary_heat_difference_w": 0.0,
                "applied_emitter_joule_heat_w": 0.0,
                "applied_collector_joule_heat_w": 0.0,
                "applied_total_joule_heat_w": 0.0,
            },
        )
        rep_totals["virtual_element_count"] += 1.0

        for node in range(n):
            node_current = float(j_a_m2[node] * area[node])
            delta_ue_uc = float(ue[node] - uc[node])
            p_gap = float(v_gap[node] * node_current)
            p_ue_uc = float(delta_ue_uc * node_current)
            p_vd = float(vd[node] * node_current)
            electron_net = float((ue[node] - uc[node]) * node_current)
            emitter_electron_source = float(q_electron_e_flux[node] * area[node])
            emitter_electron_removed = float(-emitter_electron_source)
            collector_electron_source = float(q_electron_c_flux[node] * area[node])
            boundary_heat_difference = emitter_electron_removed - collector_electron_source
            emitter_joule = float(cpp_joule_e[node])
            collector_joule = float(cpp_joule_c[node])
            total_joule = emitter_joule + collector_joule
            applied_emitter_electron_source = float(applied_q_e_flux[node] * area[node])
            applied_emitter_electron_removed = float(-applied_emitter_electron_source)
            applied_collector_electron_source = float(applied_q_c_flux[node] * area[node])
            applied_boundary_heat_difference = (
                applied_emitter_electron_removed - applied_collector_electron_source
            )
            applied_emitter_joule = float(applied_joule_e_by_node[node])
            applied_collector_joule = float(applied_joule_c_by_node[node])
            applied_total_joule = applied_emitter_joule + applied_collector_joule

            row = {
                "virtual_element_index": virtual_idx,
                "representative_tfe": tfe_name,
                "representative_local_index": local_index,
                "axial_node_index": node,
                "node_area_m2": float(area[node]),
                "J_A_cm2": float(j_a_cm2[node]),
                "J_A_m2": float(j_a_m2[node]),
                "node_current_A": node_current,
                "V_gap_v": float(v_gap[node]),
                "UE_v": float(ue[node]),
                "UC_v": float(uc[node]),
                "UE_minus_UC_v": delta_ue_uc,
                "Vd_v": float(vd[node]),
                "phiE_v": float(phi_e[node]),
                "phiC_v": float(phi_c[node]),
                "emitter_electron_heat_flux_w_m2_signed": float(q_electron_e_flux[node]),
                "emitter_electron_heat_source_w_signed": emitter_electron_source,
                "emitter_electron_heat_removed_w_positive": emitter_electron_removed,
                "collector_electron_heat_flux_w_m2_signed": float(q_electron_c_flux[node]),
                "collector_electron_heat_source_w_signed": collector_electron_source,
                "electron_boundary_heat_difference_w": boundary_heat_difference,
                "rhoE_ohm_m": float(rho_e[node]),
                "rhoC_ohm_m": float(rho_c[node]),
                "E_emitter_axial_v_m": float(e_field[node]),
                "E_collector_axial_v_m": float(c_field[node]),
                "emitter_node_volume_m3": float(emitter_node_volumes[node]),
                "collector_node_volume_m3": float(collector_node_volumes[node]),
                "emitter_joule_heat_w": emitter_joule,
                "collector_joule_heat_w": collector_joule,
                "total_joule_heat_w": total_joule,
                "applied_emitter_electron_heat_flux_w_m2_signed": float(applied_q_e_flux[node]),
                "applied_emitter_electron_heat_source_w_signed": applied_emitter_electron_source,
                "applied_emitter_electron_heat_removed_w_positive": applied_emitter_electron_removed,
                "applied_collector_electron_heat_flux_w_m2_signed": float(applied_q_c_flux[node]),
                "applied_collector_electron_heat_source_w_signed": applied_collector_electron_source,
                "applied_electron_boundary_heat_difference_w": applied_boundary_heat_difference,
                "applied_emitter_joule_heat_w": applied_emitter_joule,
                "applied_collector_joule_heat_w": applied_collector_joule,
                "applied_total_joule_heat_w": applied_total_joule,
                "node_power_gap_voltage_w": p_gap,
                "node_power_ue_minus_uc_w": p_ue_uc,
                "node_power_barrier_vd_w": p_vd,
                "node_electron_net_source_w": electron_net,
            }
            rows.append(row)
            for key, value in (
                ("node_power_gap_voltage_w", p_gap),
                ("node_power_ue_minus_uc_w", p_ue_uc),
                ("node_power_barrier_vd_w", p_vd),
                ("node_current_sum_a", node_current),
                ("node_electron_net_source_w", electron_net),
                ("emitter_electron_heat_source_signed_w", emitter_electron_source),
                ("emitter_electron_heat_removed_positive_w", emitter_electron_removed),
                ("collector_electron_heat_source_signed_w", collector_electron_source),
                ("electron_boundary_heat_difference_w", boundary_heat_difference),
                ("emitter_joule_heat_w", emitter_joule),
                ("collector_joule_heat_w", collector_joule),
                ("total_joule_heat_w", total_joule),
                ("applied_emitter_electron_heat_source_signed_w", applied_emitter_electron_source),
                ("applied_emitter_electron_heat_removed_positive_w", applied_emitter_electron_removed),
                ("applied_collector_electron_heat_source_signed_w", applied_collector_electron_source),
                ("applied_electron_boundary_heat_difference_w", applied_boundary_heat_difference),
                ("applied_emitter_joule_heat_w", applied_emitter_joule),
                ("applied_collector_joule_heat_w", applied_collector_joule),
                ("applied_total_joule_heat_w", applied_total_joule),
            ):
                totals[key] += value
                rep_totals[key] += value

    electric_global = _case_a_electric_diagnostics(core)
    terminal_power = float(electric_global.get("tec_total_electric_power_w") or 0.0)
    totals["terminal_power_plus_total_joule_heat_w"] = (
        terminal_power + totals["total_joule_heat_w"]
    )
    totals["boundary_difference_minus_terminal_power_w"] = (
        totals["electron_boundary_heat_difference_w"] - terminal_power
    )
    totals["boundary_difference_minus_terminal_power_plus_joule_w"] = (
        totals["electron_boundary_heat_difference_w"]
        - totals["terminal_power_plus_total_joule_heat_w"]
    )
    totals["terminal_power_plus_applied_total_joule_heat_w"] = (
        terminal_power + totals["applied_total_joule_heat_w"]
    )
    totals["applied_boundary_difference_minus_terminal_power_w"] = (
        totals["applied_electron_boundary_heat_difference_w"] - terminal_power
    )
    totals["applied_boundary_difference_minus_terminal_power_plus_joule_w"] = (
        totals["applied_electron_boundary_heat_difference_w"]
        - totals["terminal_power_plus_applied_total_joule_heat_w"]
    )

    summary = {
        "time_s": float(build["system"].global_time),
        "electric_global": electric_global,
        "node_integrated_totals": totals,
        "by_representative_tfe": by_representative,
        "notes": [
            "J from ThermoCalc is stored in A/cm2; J_A_m2 = J_A_cm2 * 1e4.",
            "node_current_A = J_A_m2 * emitter_outer_surface_node_area_m2.",
            "node_power_gap_voltage_w = V_gap_v * node_current_A.",
            "node_power_ue_minus_uc_w = (UE_v - UC_v) * node_current_A.",
            "emitter_electron_heat_source_w_signed = -J_A_m2 * (phiE + 2*kB_eV*TE) * emitter_area; negative means emitter cooling.",
            "collector_electron_heat_source_w_signed = J_A_m2 * (phiE + 2*kB_eV*TE - (UE-UC)) * emitter_area; positive means collector heating.",
            "electron_boundary_heat_difference_w = -emitter_electron_heat_source_w_signed - collector_electron_heat_source_w_signed.",
            "Joule heat follows the non-uniform axial-grid helper: q_vol = (dU/dy)^2 / rho on node-center coordinates, integrated over electrode cell volumes.",
            "Fields prefixed with applied_ are the cached heat sources stored in the loaded TFE state and actually applied to thermal boundaries; non-prefixed heat flux/Joule fields are recomputed from the instantaneous TEC solution.",
            "These node-integrated powers are local plasma/electrode voltage-current sums and are not identical to terminal Uout*Iout.",
        ],
    }
    return rows, summary


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    restart_stem = Path(args.restart).stem
    prefix = args.prefix or f"{restart_stem}_electric_nodes"

    build = build_and_load(args.restart)
    rows, summary = extract_node_rows(build)

    csv_path = output_dir / f"{prefix}.csv"
    json_path = output_dir / f"{prefix}_summary.json"
    if rows:
        with csv_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    json_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=json_default),
        encoding="utf-8",
    )

    print(f"CSV: {csv_path}")
    print(f"Summary: {json_path}")
    print(json.dumps(summary["electric_global"], ensure_ascii=False, indent=2))
    print(json.dumps(summary["node_integrated_totals"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
