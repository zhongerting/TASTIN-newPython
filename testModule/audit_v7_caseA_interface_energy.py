import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from Solvers.Couplers import FluidSolidCouple, GapCouple2D, TECCouple2D
from run_v7_caseA_multipliers_short import (
    build_loaded_case,
    collect_tec_stats,
    json_default,
)
from test_core_assemble_v7_caseA import _case_a_electric_diagnostics, _case_a_flow_diagnostics


DEFAULT_RESTART = (
    "testModule/v7_caseA_newpyd_long20000_from_t23800/"
    "v7_caseA_newpyd_long20000_from_t23800_latest_restart.npz"
)


def parse_ring_multipliers(text: str, *, allow_zero: bool = False) -> List[int]:
    values = [int(part.strip()) for part in str(text).split(",") if part.strip()]
    if len(values) != 4:
        raise argparse.ArgumentTypeError("Use four comma-separated integers, e.g. 1,6,12,18.")
    if allow_zero:
        bad = [value for value in values if value < 0]
        message = "Multipliers must be non-negative."
    else:
        bad = [value for value in values if value <= 0]
        message = "Multipliers must be positive."
    if bad:
        raise argparse.ArgumentTypeError(message)
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit V7 CaseA interface energy closure at one restart state."
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v7_caseA_interface_audit")
    parser.add_argument("--case-prefix", default="v7_caseA_interface_audit")
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_ring_multipliers(text, allow_zero=False),
        default=parse_ring_multipliers("1,6,12,18"),
    )
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_ring_multipliers(text, allow_zero=True),
        default=parse_ring_multipliers("1,6,12,15", allow_zero=True),
    )
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    return parser.parse_args()


def finite_or_none(value: Any) -> Optional[float]:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(value):
        return None
    return value


def sum_or_zero(values: Iterable[float]) -> float:
    return float(sum(float(value) for value in values))


def axial_centers(tfe) -> np.ndarray:
    y_faces = np.asarray(tfe.common_y_faces, dtype=float)
    return 0.5 * (y_faces[:-1] + y_faces[1:])


def axial_sum_from_flat(solid, flat_values: np.ndarray, n_nodes: int) -> np.ndarray:
    arr = np.asarray(flat_values, dtype=float)
    if arr.size == n_nodes:
        return arr.copy()
    return arr.reshape(solid.nx, n_nodes).sum(axis=0)


def snapshot_moderator_mapping(build: Dict[str, Any], label: str) -> List[Dict[str, Any]]:
    core = build["core"]
    snapshots: List[Dict[str, Any]] = []
    if not getattr(core, "has_global_moderator", False):
        return snapshots

    for ring_idx, ring in enumerate(core.mod_rings):
        member_names = core.get_ring_member_names(ring_idx)
        if not member_names:
            continue
        q_expected = None
        tfe_entries = []
        for tfe_name in member_names:
            tfe = build["tfes"][tfe_name]
            mult = float(build["ring_multipliers"][tfe_name])
            boundary = tfe.solids["moderator"].boundaries["right"]
            q_out = -np.asarray(boundary.current_flux, dtype=float).copy()
            if q_expected is None:
                q_expected = np.zeros_like(q_out, dtype=float)
            q_expected += q_out * mult
            tfe_entries.append({
                "tfe_name": tfe_name,
                "thermal_multiplier": mult,
                "q_out_single": q_out,
                "surface_temperature": np.asarray(boundary.T_surface, dtype=float).copy(),
                "boundary_temperature": np.asarray(tfe.boundary_data.moderator_temperature, dtype=float).copy(),
            })

        snapshots.append({
            "label": label,
            "ring_idx": ring_idx,
            "expected_total": np.zeros(ring.shape_nodes[1], dtype=float)
            if q_expected is None
            else q_expected,
            "actual_ring_source": ring_q_source_axial(ring).copy(),
            "global_average_temperature": ring_average_temperature_axial(ring).copy(),
            "tfe_entries": tfe_entries,
        })
    return snapshots


def sync_for_audit(build: Dict[str, Any]) -> Dict[str, Any]:
    """Refresh one restart state into the same in-memory coupling口径 used by the solver."""
    system = build["system"]
    core = build["core"]
    current_time = float(system.global_time)

    core.post_step(0.0, current_time)
    if core.thermo_calc is not None:
        core.thermo_calc.calculate(verbose=False)

    system._prepare_fluid_sources_for_coupling()
    system._refresh_solid_boundary_cache(update_flux=True, current_time=current_time)
    core.pre_step(0.0, current_time)
    moderator_pre_step = snapshot_moderator_mapping(build, "pre_step_write")
    system._run_couplers(interface_relaxation=1.0, current_time=current_time)
    system._refresh_solid_boundary_cache(update_flux=True, current_time=current_time)
    moderator_post_coupler = snapshot_moderator_mapping(build, "post_coupler_flux")
    return {
        "moderator_pre_step": moderator_pre_step,
        "moderator_post_coupler": moderator_post_coupler,
    }


def gap_passive_terms(coupler: GapCouple2D) -> Dict[str, np.ndarray]:
    T1_surf, _ = coupler.bound1.get_coupling_surface_snapshot()
    T2_surf, _ = coupler.bound2.get_coupling_surface_snapshot()
    T1_surf = np.asarray(T1_surf, dtype=float)
    T2_surf = np.asarray(T2_surf, dtype=float)

    A1 = np.maximum(np.asarray(coupler.bound1.area, dtype=float), 1.0e-12)
    A2 = np.maximum(np.asarray(coupler.bound2.area, dtype=float), 1.0e-12)
    is_1_inner = A1 < A2
    A_in = np.where(is_1_inner, A1, A2)
    A_out = np.where(is_1_inner, A2, A1)
    eps_in = np.where(is_1_inner, coupler.eps1, coupler.eps2)
    eps_out = np.where(is_1_inner, coupler.eps2, coupler.eps1)

    h_rad_star = coupler.sigma * (T1_surf**2 + T2_surf**2) * (T1_surf + T2_surf)
    h_rad_star = np.maximum(h_rad_star, 1.0e-20)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom_rad = (1.0 / eps_in) + ((1.0 / eps_out) - 1.0) * (A_in / A_out)
        G_rad = (A_in * h_rad_star) / denom_rad
        G_rad = np.nan_to_num(G_rad, nan=0.0, posinf=0.0, neginf=0.0)
    G_cond = (coupler.k_gas * A_in) / coupler.gap
    dT = T1_surf - T2_surf
    return {
        "T1_surf": T1_surf,
        "T2_surf": T2_surf,
        "G_cond": G_cond,
        "G_rad": G_rad,
        "q_cond_1_to_2": G_cond * dT,
        "q_rad_1_to_2": G_rad * dT,
        "q_total_1_to_2": (G_cond + G_rad) * dT,
    }


def tec_result_by_tfe(build: Dict[str, Any]) -> Dict[str, Tuple[int, Optional[Dict[str, Any]]]]:
    core = build["core"]
    thermo_calc = core.thermo_calc
    idx = 0
    result = {}
    for name, mult in build["tec_ring_multipliers"].items():
        tec_mult = int(mult)
        if thermo_calc is None or tec_mult <= 0:
            result[name] = (idx, None)
            continue
        result[name] = (idx, thermo_calc.get_tec_results(idx))
        idx += tec_mult
    return result


def tec_node_electric_terms(build: Dict[str, Any], name: str, res: Dict[str, Any]) -> Dict[str, np.ndarray]:
    n_nodes = len(np.asarray(res["J"], dtype=float))
    input_data = build["core"].thermo_calc._input_data
    idx = tec_result_by_tfe(build)[name][0]
    area = np.asarray(input_data.sideAreaE[idx], dtype=float)
    j_a_m2 = np.asarray(res["J"], dtype=float) * 1.0e4
    node_current = j_a_m2 * area

    te = np.asarray(res["TE"], dtype=float)
    phi_e = np.asarray(res["phiE"], dtype=float)
    ue = np.asarray(res["UE"], dtype=float)
    uc = np.asarray(res["UC"], dtype=float)

    q_e_flux = -j_a_m2 * (phi_e + 2.0 * 8.617e-5 * te)
    q_c_flux = j_a_m2 * (phi_e + 2.0 * 8.617e-5 * te - (ue - uc))
    emitter_electron = q_e_flux * area
    collector_electron = q_c_flux * area
    electron_diff = -emitter_electron - collector_electron

    emitter_joule = np.asarray(res["joulePowerE"], dtype=float)
    collector_joule = np.asarray(res["joulePowerC"], dtype=float)
    joule_total = emitter_joule + collector_joule

    terminal_single = float(res.get("U", 0.0)) * float(res.get("I", 0.0))
    weights = np.abs(node_current)
    if float(np.sum(weights)) <= 0.0:
        weights = np.ones(n_nodes, dtype=float)
    terminal_alloc = terminal_single * weights / float(np.sum(weights))

    return {
        "node_current_a": node_current,
        "emitter_electron_source_w": emitter_electron,
        "collector_electron_source_w": collector_electron,
        "electron_boundary_diff_w": electron_diff,
        "emitter_joule_w": emitter_joule,
        "collector_joule_w": collector_joule,
        "total_joule_w": joule_total,
        "terminal_power_alloc_w": terminal_alloc,
    }


def collect_interface_rows(build: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    tec_results = tec_result_by_tfe(build)
    electric = _case_a_electric_diagnostics(build["core"])
    global_terminal_power_w = float(electric["tec_total_electric_power_w"] or 0.0)
    summary = {
        "ordinary_scaled_residual_w": 0.0,
        "ordinary_scaled_abs_residual_w": 0.0,
        "ordinary_max_abs_single_node_residual_w": 0.0,
        "tec": {
            "passive_closure_scaled_w": 0.0,
            "instant_electron_boundary_diff_w": 0.0,
            "instant_joule_heat_w": 0.0,
            "instant_terminal_power_alloc_w": 0.0,
            "instant_tec_energy_gap_w": 0.0,
            "applied_electron_boundary_diff_w": 0.0,
            "applied_joule_heat_w": 0.0,
            "applied_terminal_power_alloc_w": 0.0,
            "applied_tec_energy_gap_w": 0.0,
        },
    }

    for name, tfe in build["tfes"].items():
        thermal_mult = float(build["ring_multipliers"][name])
        tec_mult = float(build["tec_ring_multipliers"][name])
        centers = axial_centers(tfe)
        for coupler_name, coupler in tfe.couplers.items():
            if isinstance(coupler, FluidSolidCouple):
                continue
            if not hasattr(coupler, "bound1") or not hasattr(coupler, "bound2"):
                continue

            b1 = coupler.bound1
            b2 = coupler.bound2
            q1_in = np.asarray(b1.current_flux, dtype=float)
            q2_in = np.asarray(b2.current_flux, dtype=float)
            n_nodes = min(q1_in.size, q2_in.size, centers.size)
            passive = gap_passive_terms(coupler) if isinstance(coupler, GapCouple2D) else None

            applied_q1 = np.zeros(n_nodes, dtype=float)
            applied_q2 = np.zeros(n_nodes, dtype=float)
            applied_joule = np.zeros(n_nodes, dtype=float)
            instant = None
            if isinstance(coupler, TECCouple2D):
                applied_q1 = np.asarray(coupler.Q_source_1, dtype=float)[:n_nodes]
                applied_q2 = np.asarray(coupler.Q_source_2, dtype=float)[:n_nodes]
                applied_joule = (
                    axial_sum_from_flat(
                        tfe.solids["emitter"],
                        tfe.electric_data.emitter_joule_heat,
                        n_nodes,
                    )
                    + axial_sum_from_flat(
                        tfe.solids["collector"],
                        tfe.electric_data.collector_joule_heat,
                        n_nodes,
                    )
                )
                _, res = tec_results[name]
                if res is not None:
                    instant = tec_node_electric_terms(build, name, res)

            for node in range(n_nodes):
                row_multiplier = tec_mult if isinstance(coupler, TECCouple2D) else thermal_mult
                residual_single = float(q1_in[node] + q2_in[node])
                row = {
                    "interface_type": type(coupler).__name__,
                    "representative_tfe": name,
                    "coupler_name": coupler_name,
                    "axial_node": node,
                    "y_center_m": float(centers[node]),
                    "thermal_multiplier": thermal_mult,
                    "tec_multiplier": tec_mult,
                    "row_multiplier": row_multiplier,
                    "side1_solid": getattr(coupler.obj1, "name", ""),
                    "side2_solid": getattr(coupler.obj2, "name", ""),
                    "side1_boundary_inflow_w_single": float(q1_in[node]),
                    "side2_boundary_inflow_w_single": float(q2_in[node]),
                    "boundary_inflow_sum_w_single": residual_single,
                    "boundary_inflow_sum_w_scaled": residual_single * row_multiplier,
                    "passive_cond_1_to_2_w_single": None,
                    "passive_rad_1_to_2_w_single": None,
                    "passive_total_1_to_2_w_single": None,
                    "passive_side1_inflow_w_single": None,
                    "passive_side2_inflow_w_single": None,
                    "passive_closure_w_single": None,
                    "emitter_electron_source_w_single": None,
                    "collector_electron_source_w_single": None,
                    "electron_boundary_diff_w_single": None,
                    "terminal_power_alloc_w_single": None,
                    "node_current_a_single": None,
                    "joule_heat_w_single": None,
                    "tec_energy_gap_w_single": None,
                    "applied_emitter_electron_source_w_single": None,
                    "applied_collector_electron_source_w_single": None,
                    "applied_electron_boundary_diff_w_single": None,
                    "applied_joule_heat_w_single": None,
                    "applied_tec_energy_gap_w_single": None,
                }

                if passive is not None:
                    q_passive = float(passive["q_total_1_to_2"][node])
                    row.update({
                        "passive_cond_1_to_2_w_single": float(passive["q_cond_1_to_2"][node]),
                        "passive_rad_1_to_2_w_single": float(passive["q_rad_1_to_2"][node]),
                        "passive_total_1_to_2_w_single": q_passive,
                        "passive_side1_inflow_w_single": -q_passive,
                        "passive_side2_inflow_w_single": q_passive,
                        "passive_closure_w_single": 0.0,
                    })

                if isinstance(coupler, TECCouple2D):
                    applied_diff = float(-applied_q1[node] - applied_q2[node])
                    row.update({
                        "applied_emitter_electron_source_w_single": float(applied_q1[node]),
                        "applied_collector_electron_source_w_single": float(applied_q2[node]),
                        "applied_electron_boundary_diff_w_single": applied_diff,
                        "applied_joule_heat_w_single": float(applied_joule[node]),
                    })
                    summary["tec"]["applied_electron_boundary_diff_w"] += applied_diff * tec_mult
                    summary["tec"]["applied_joule_heat_w"] += float(applied_joule[node]) * tec_mult

                    if instant is not None:
                        diff = float(instant["electron_boundary_diff_w"][node])
                        joule = float(instant["total_joule_w"][node])
                        row.update({
                            "node_current_a_single": float(instant["node_current_a"][node]),
                            "emitter_electron_source_w_single": float(
                                instant["emitter_electron_source_w"][node]
                            ),
                            "collector_electron_source_w_single": float(
                                instant["collector_electron_source_w"][node]
                            ),
                            "electron_boundary_diff_w_single": diff,
                            "joule_heat_w_single": joule,
                        })
                        summary["tec"]["instant_electron_boundary_diff_w"] += diff * tec_mult
                        summary["tec"]["instant_joule_heat_w"] += joule * tec_mult
                else:
                    scaled = residual_single * thermal_mult
                    summary["ordinary_scaled_residual_w"] += scaled
                    summary["ordinary_scaled_abs_residual_w"] += abs(scaled)
                    summary["ordinary_max_abs_single_node_residual_w"] = max(
                        summary["ordinary_max_abs_single_node_residual_w"],
                        abs(residual_single),
                    )

                if isinstance(coupler, TECCouple2D) and passive is not None:
                    closure = float(row["passive_closure_w_single"] or 0.0)
                    summary["tec"]["passive_closure_scaled_w"] += closure * tec_mult

                rows.append(row)

    terminal_weight = 0.0
    for row in rows:
        if row.get("interface_type") != "TECCouple2D":
            continue
        current = finite_or_none(row.get("node_current_a_single"))
        if current is None:
            continue
        terminal_weight += abs(current) * float(row["tec_multiplier"])

    if terminal_weight > 0.0:
        for row in rows:
            if row.get("interface_type") != "TECCouple2D":
                continue
            current = finite_or_none(row.get("node_current_a_single"))
            if current is None:
                continue
            tec_mult = float(row["tec_multiplier"])
            term_single = global_terminal_power_w * abs(current) / terminal_weight
            row["terminal_power_alloc_w_single"] = term_single

            diff = finite_or_none(row.get("electron_boundary_diff_w_single"))
            joule = finite_or_none(row.get("joule_heat_w_single"))
            if diff is not None and joule is not None:
                gap = term_single + joule - diff
                row["tec_energy_gap_w_single"] = gap
                summary["tec"]["instant_terminal_power_alloc_w"] += term_single * tec_mult
                summary["tec"]["instant_tec_energy_gap_w"] += gap * tec_mult

            applied_diff = finite_or_none(row.get("applied_electron_boundary_diff_w_single"))
            applied_joule = finite_or_none(row.get("applied_joule_heat_w_single"))
            if applied_diff is not None and applied_joule is not None:
                applied_gap = term_single + applied_joule - applied_diff
                row["applied_tec_energy_gap_w_single"] = applied_gap
                summary["tec"]["applied_terminal_power_alloc_w"] += term_single * tec_mult
                summary["tec"]["applied_tec_energy_gap_w"] += applied_gap * tec_mult

    return rows, summary


def fsc_source_terms(coupler: FluidSolidCouple) -> Dict[str, np.ndarray]:
    if coupler._last_lambda is not None:
        lam = np.asarray(coupler._last_lambda, dtype=float)
    else:
        T_f = np.asarray(coupler.fluid.temperature_vector, dtype=float)
        P_f = coupler.fluid.pressure_vector
        rho = coupler.fluid.density_vector
        vel = coupler.fluid.velocity_vector
        mu = np.maximum(coupler.fluid.material.viscosity(T_f, P_f), 1.0e-10)
        k_f = coupler.fluid.material.conductivity(T_f, P_f)
        Pr = coupler.fluid.material.prandtl_number(T_f, P_f)
        Re = (rho * np.abs(vel) * coupler.fluid.d_h) / mu
        Nu = coupler.correlation_func(Re, Pr, 1.1)
        lam = np.maximum(Nu * k_f / coupler.fluid.d_h * coupler.node_areas, 1.0e-8)

    T_f = np.asarray(coupler.fluid.temperature_vector, dtype=float)
    T_wall = np.asarray(coupler.solid_bound.T_surface, dtype=float)
    source = lam * (T_wall - T_f)
    return {
        "lambda_w_per_k": lam,
        "T_fluid_k": T_f,
        "T_wall_k": T_wall,
        "source_w": source,
    }


def collect_coolant_rows(build: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    summary = {
        "solid_to_fluid_w": 0.0,
        "recomputed_fsc_source_w": 0.0,
        "merged_fluid_source_w": 0.0,
        "channel_enthalpy_pickup_w": 0.0,
        "solid_minus_recomputed_source_w": 0.0,
        "source_minus_enthalpy_w": 0.0,
        "plenum_enthalpy_pickup_w": None,
    }

    for name, tfe in build["tfes"].items():
        channel = build["fluid_channels"][name]
        mult = float(build["ring_multipliers"][name])
        inner = tfe.couplers["iclad_coolant_fsc"]
        outer = tfe.couplers["oclad_coolant_fsc"]
        inner_terms = fsc_source_terms(inner)
        outer_terms = fsc_source_terms(outer)
        inner_solid_to_fluid = -np.asarray(inner.solid_bound.current_flux, dtype=float)
        outer_solid_to_fluid = -np.asarray(outer.solid_bound.current_flux, dtype=float)
        merged = np.array([
            vol.Q_wall + vol.Q_vol - vol.implicit_coeff * vol.T
            for vol in channel.volumes
        ], dtype=float)

        n_nodes = channel.n_nodes
        centers = axial_centers(tfe)
        for node in range(n_nodes):
            source_node = float(inner_terms["source_w"][node] + outer_terms["source_w"][node])
            solid_node = float(inner_solid_to_fluid[node] + outer_solid_to_fluid[node])
            rows.append({
                "row_type": "node",
                "representative_tfe": name,
                "thermal_multiplier": mult,
                "axial_node": node,
                "y_center_m": float(centers[node]) if node < centers.size else None,
                "inner_wall_temperature_k": float(inner_terms["T_wall_k"][node]),
                "outer_wall_temperature_k": float(outer_terms["T_wall_k"][node]),
                "fluid_temperature_k": float(inner_terms["T_fluid_k"][node]),
                "inner_lambda_w_per_k": float(inner_terms["lambda_w_per_k"][node]),
                "outer_lambda_w_per_k": float(outer_terms["lambda_w_per_k"][node]),
                "inner_solid_to_fluid_w_single": float(inner_solid_to_fluid[node]),
                "outer_solid_to_fluid_w_single": float(outer_solid_to_fluid[node]),
                "solid_to_fluid_w_single": solid_node,
                "inner_recomputed_fluid_source_w_single": float(inner_terms["source_w"][node]),
                "outer_recomputed_fluid_source_w_single": float(outer_terms["source_w"][node]),
                "recomputed_fluid_source_w_single": source_node,
                "merged_fluid_source_w_single": float(merged[node]),
                "solid_minus_recomputed_w_single": solid_node - source_node,
                "merged_minus_recomputed_w_single": float(merged[node]) - source_node,
                "channel_single_flow_kg_s": float(channel.internal_junctions[0].W),
                "channel_enthalpy_pickup_w_scaled": None,
                "source_minus_enthalpy_w_scaled": None,
            })

        w_single = float(channel.internal_junctions[0].W)
        h_in = float(channel.volumes[0].h)
        h_out = float(channel.volumes[-1].h)
        enthalpy = w_single * mult * (h_out - h_in)
        solid_sum = float(np.sum(inner_solid_to_fluid + outer_solid_to_fluid)) * mult
        source_sum = float(np.sum(inner_terms["source_w"] + outer_terms["source_w"])) * mult
        merged_sum = float(np.sum(merged)) * mult
        rows.append({
            "row_type": "channel_total",
            "representative_tfe": name,
            "thermal_multiplier": mult,
            "axial_node": -1,
            "y_center_m": None,
            "inner_wall_temperature_k": None,
            "outer_wall_temperature_k": None,
            "fluid_temperature_k": None,
            "inner_lambda_w_per_k": None,
            "outer_lambda_w_per_k": None,
            "inner_solid_to_fluid_w_single": float(np.sum(inner_solid_to_fluid)),
            "outer_solid_to_fluid_w_single": float(np.sum(outer_solid_to_fluid)),
            "solid_to_fluid_w_single": float(np.sum(inner_solid_to_fluid + outer_solid_to_fluid)),
            "inner_recomputed_fluid_source_w_single": float(np.sum(inner_terms["source_w"])),
            "outer_recomputed_fluid_source_w_single": float(np.sum(outer_terms["source_w"])),
            "recomputed_fluid_source_w_single": float(np.sum(inner_terms["source_w"] + outer_terms["source_w"])),
            "merged_fluid_source_w_single": float(np.sum(merged)),
            "solid_minus_recomputed_w_single": float(
                np.sum(inner_solid_to_fluid + outer_solid_to_fluid)
                - np.sum(inner_terms["source_w"] + outer_terms["source_w"])
            ),
            "merged_minus_recomputed_w_single": float(np.sum(merged) - np.sum(inner_terms["source_w"] + outer_terms["source_w"])),
            "channel_single_flow_kg_s": w_single,
            "channel_inlet_temperature_k": float(channel.volumes[0].T),
            "channel_outlet_temperature_k": float(channel.volumes[-1].T),
            "channel_inlet_enthalpy_j_kg": h_in,
            "channel_outlet_enthalpy_j_kg": h_out,
            "channel_enthalpy_pickup_w_scaled": enthalpy,
            "source_minus_enthalpy_w_scaled": source_sum - enthalpy,
        })

        summary["solid_to_fluid_w"] += solid_sum
        summary["recomputed_fsc_source_w"] += source_sum
        summary["merged_fluid_source_w"] += merged_sum
        summary["channel_enthalpy_pickup_w"] += enthalpy
        summary["solid_minus_recomputed_source_w"] += solid_sum - source_sum
        summary["source_minus_enthalpy_w"] += source_sum - enthalpy

    flow = _case_a_flow_diagnostics(build)
    mass_flow = float(flow["inlet_total_macro_flow_kg_s"])
    summary["plenum_enthalpy_pickup_w"] = mass_flow * (
        float(build["outlet_plenum"].h) - float(build["inlet_plenum"].h)
    )
    summary["channel_sum_minus_plenum_enthalpy_w"] = (
        summary["channel_enthalpy_pickup_w"] - summary["plenum_enthalpy_pickup_w"]
    )
    return rows, summary


def compute_global_energy_snapshot(build: Dict[str, Any], electric: Dict[str, Any]) -> Dict[str, float]:
    core = build["core"]
    system = build["system"]
    sigma = 5.670374419e-8

    flow = _case_a_flow_diagnostics(build)
    mass_flow = float(flow["inlet_total_macro_flow_kg_s"])
    q_fluid = mass_flow * (float(build["outlet_plenum"].h) - float(build["inlet_plenum"].h))
    q_electric = float(electric["tec_total_electric_power_w"] or 0.0)

    reflector_outer = core.reflector.boundaries["right"]
    surface_temperature = np.asarray(reflector_outer.T_surface, dtype=float)
    radiation_area = np.asarray(reflector_outer.area, dtype=float)
    q_rad_nodes = (
        0.2
        * sigma
        * radiation_area
        * (surface_temperature**4 - 200.0**4)
    )
    q_radiation = float(np.sum(q_rad_nodes))

    q_core = sum(
        float(tfe.neutronic_data.total_power) * float(core.tfe_multipliers[name])
        for name, tfe in build["tfes"].items()
    )
    residual = q_core - q_fluid - q_electric - q_radiation

    return {
        "time_s": float(system.global_time),
        "solid_heat_capacity_scale": float(build["solid_heat_capacity_scale"]),
        "core_heat_power_w": q_core,
        "coolant_heat_pickup_w": q_fluid,
        "electric_power_w": q_electric,
        "outer_wall_radiation_w": q_radiation,
        "balance_residual_w": residual,
        "balance_residual_percent": 100.0 * residual / q_core,
        "mass_flow_kg_s": mass_flow,
        "inlet_plenum_temperature_k": float(build["inlet_plenum"].T),
        "outlet_plenum_temperature_k": float(build["outlet_plenum"].T),
        "outer_wall_area_m2": float(np.sum(radiation_area)),
        "outer_wall_radiation_area_avg_flux_w_m2": q_radiation / float(np.sum(radiation_area)),
        "tec_total_voltage_v": electric["tec_total_voltage_v"],
        "tec_total_current_a": electric["tec_total_current_a"],
    }


def ring_q_source_axial(ring) -> np.ndarray:
    nx, ny = ring.shape_nodes
    return np.asarray(ring.Q_source, dtype=float).reshape(nx, ny).sum(axis=0)


def ring_average_temperature_axial(ring) -> np.ndarray:
    nx, ny = ring.shape_nodes
    return np.asarray(ring.T, dtype=float).reshape(nx, ny).mean(axis=0)


def collect_moderator_rows(
    build: Dict[str, Any],
    snapshots: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    core = build["core"]
    rows: List[Dict[str, Any]] = []
    summary = {
        "max_abs_mapping_residual_w": 0.0,
        "sum_mapping_residual_w": 0.0,
        "sum_expected_q_w": 0.0,
        "sum_actual_ring_q_source_w": 0.0,
        "max_abs_temperature_backfill_delta_k": 0.0,
    }

    if not getattr(core, "has_global_moderator", False):
        return rows, summary

    if snapshots is None:
        snapshots = snapshot_moderator_mapping(build, "current")

    for snapshot in snapshots:
        label = snapshot["label"]
        ring_idx = int(snapshot["ring_idx"])
        actual = np.asarray(snapshot["actual_ring_source"], dtype=float)
        expected_total = np.asarray(snapshot["expected_total"], dtype=float)
        global_t = np.asarray(snapshot["global_average_temperature"], dtype=float)

        for tfe_entry in snapshot["tfe_entries"]:
            tfe_name = tfe_entry["tfe_name"]
            mult = float(tfe_entry["thermal_multiplier"])
            q_out = np.asarray(tfe_entry["q_out_single"], dtype=float)
            surface_t = np.asarray(tfe_entry["surface_temperature"], dtype=float)
            boundary_t = np.asarray(tfe_entry["boundary_temperature"], dtype=float)
            for node in range(actual.size):
                residual = float(actual[node] - expected_total[node])
                backfill_delta = float(boundary_t[node] - global_t[node])
                rows.append({
                    "snapshot_label": label,
                    "ring_idx": ring_idx,
                    "ring_name": tfe_name if len(snapshot["tfe_entries"]) == 1 else f"Ring{ring_idx}",
                    "representative_tfe": tfe_name,
                    "thermal_multiplier": mult,
                    "axial_node": node,
                    "q_out_of_tfe_moderator_w_single": float(q_out[node]),
                    "q_expected_from_this_tfe_w": float(q_out[node] * mult),
                    "q_expected_ring_total_w": float(expected_total[node]),
                    "q_actual_global_ring_source_w": float(actual[node]),
                    "mapping_residual_actual_minus_expected_w": residual,
                    "tfe_moderator_outer_surface_temperature_k": float(surface_t[node]),
                    "global_moderator_axial_average_temperature_k": float(global_t[node]),
                    "tfe_boundary_data_moderator_temperature_k": float(boundary_t[node]),
                    "boundary_minus_global_average_k": backfill_delta,
                    "surface_minus_boundary_k": float(surface_t[node] - boundary_t[node]),
                })
                if label == "pre_step_write":
                    summary["max_abs_mapping_residual_w"] = max(
                        summary["max_abs_mapping_residual_w"],
                        abs(residual),
                    )
                    summary["max_abs_temperature_backfill_delta_k"] = max(
                        summary["max_abs_temperature_backfill_delta_k"],
                        abs(backfill_delta),
                    )
        if label == "pre_step_write":
            summary["sum_mapping_residual_w"] += float(np.sum(actual - expected_total))
            summary["sum_expected_q_w"] += float(np.sum(expected_total))
            summary["sum_actual_ring_q_source_w"] += float(np.sum(actual))

    return rows, summary


def infer_sign_convention(interface_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    for row in interface_rows:
        if row["interface_type"] == "TECCouple2D":
            continue
        passive = finite_or_none(row.get("passive_total_1_to_2_w_single"))
        side1 = finite_or_none(row.get("side1_boundary_inflow_w_single"))
        side2 = finite_or_none(row.get("side2_boundary_inflow_w_single"))
        if passive is None or side1 is None or side2 is None:
            continue
        if abs(passive) < 1.0e-12:
            continue
        return {
            "BoundaryRegion.current_flux": "positive means heat inflow to the owning solid node",
            "solid_to_neighbor_for_outward_boundary": "-current_flux",
            "validation_interface": {
                "representative_tfe": row["representative_tfe"],
                "coupler_name": row["coupler_name"],
                "axial_node": row["axial_node"],
                "passive_1_to_2_w": passive,
                "side1_current_flux_w": side1,
                "side2_current_flux_w": side2,
                "observed_relation": (
                    "side1 current_flux has opposite sign to passive 1_to_2; "
                    "side2 current_flux has same sign"
                ),
            },
        }
    return {
        "BoundaryRegion.current_flux": "positive means heat inflow to the owning solid node",
        "validation_interface": None,
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    build = build_loaded_case(args)
    sync_snapshots = sync_for_audit(build)

    interface_rows, interface_summary = collect_interface_rows(build)
    coolant_rows, coolant_summary = collect_coolant_rows(build)
    moderator_rows, moderator_summary = collect_moderator_rows(
        build,
        sync_snapshots["moderator_pre_step"] + sync_snapshots["moderator_post_coupler"],
    )

    interface_nodes_path = output_dir / f"{args.case_prefix}_interface_nodes.csv"
    coolant_path = output_dir / f"{args.case_prefix}_coolant_channel_audit.csv"
    moderator_path = output_dir / f"{args.case_prefix}_moderator_mapping_audit.csv"
    summary_path = output_dir / f"{args.case_prefix}_interface_summary.json"

    write_csv(interface_nodes_path, interface_rows)
    write_csv(coolant_path, coolant_rows)
    write_csv(moderator_path, moderator_rows)

    system = build["system"]
    core = build["core"]
    electric = _case_a_electric_diagnostics(core)
    flow = _case_a_flow_diagnostics(build)
    global_energy = compute_global_energy_snapshot(build, electric)

    summary = {
        "restart_in": args.restart_in,
        "time_s": float(system.global_time),
        "target_voltage_v": float(args.target_voltage),
        "ring_multipliers": build["ring_multipliers"],
        "tec_ring_multipliers": build["tec_ring_multipliers"],
        "sign_convention": infer_sign_convention(interface_rows),
        "outputs": {
            "interface_nodes_csv": str(interface_nodes_path),
            "coolant_channel_audit_csv": str(coolant_path),
            "moderator_mapping_audit_csv": str(moderator_path),
            "interface_summary_json": str(summary_path),
        },
        "electric": electric,
        "flow": flow,
        "global_energy": global_energy,
        "interface_summary": interface_summary,
        "coolant_summary": coolant_summary,
        "moderator_summary": moderator_summary,
        "tec_stats_reference": collect_tec_stats(build),
        "notes": {
            "sync_sequence": (
                "load restart -> core.post_step(0,t) -> thermo_calc.calculate(False) -> "
                "clear fluid sources -> refresh boundary fluxes -> core.pre_step(0,t) "
                "and snapshot moderator mapping -> run couplers once -> refresh boundary fluxes"
            ),
            "tec_energy_gap_definition": (
                "terminal_power_alloc_w + joule_heat_w - electron_boundary_diff_w; "
                "terminal power is allocated over axial nodes by absolute node current"
            ),
            "coolant_source_definition": (
                "fluid_source_w = Q_wall + Q_vol - implicit_coeff*T; per-coupler source is "
                "recomputed as lambda*(T_wall - T_fluid)"
            ),
        },
    }

    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)

    print(f"Summary: {summary_path}")
    print(f"Interface nodes: {interface_nodes_path}")
    print(f"Coolant audit: {coolant_path}")
    print(f"Moderator audit: {moderator_path}")
    print(json.dumps({
        "time_s": float(system.global_time),
        "ordinary_solid_interface_abs_residual_w": interface_summary["ordinary_scaled_abs_residual_w"],
        "tec_instant_energy_gap_w": interface_summary["tec"]["instant_tec_energy_gap_w"],
        "tec_applied_energy_gap_w": interface_summary["tec"]["applied_tec_energy_gap_w"],
        "coolant_source_minus_enthalpy_w": coolant_summary["source_minus_enthalpy_w"],
        "moderator_mapping_max_abs_residual_w": moderator_summary["max_abs_mapping_residual_w"],
    }, indent=2))


if __name__ == "__main__":
    main()
