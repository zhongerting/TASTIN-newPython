import argparse
import csv
import json
import math
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

from Materials.Fluids.Sodium import Sodium
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
    MacroFlowJunction,
)
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from audit_v7_caseA_interface_energy import sync_for_audit
from run_v7_caseA_multipliers_short import build_loaded_case, json_default, parse_multipliers
from test_core_assemble_v7_caseA import _case_a_electric_diagnostics, _case_a_flow_diagnostics


DEFAULT_RESTART = (
    "testModule/v7_caseA_newpyd_long20000_after_tec_heatfix/"
    "v7_caseA_newpyd_long20000_after_tec_heatfix_latest_restart.npz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read-only fluid source, enthalpy, and macro multiplier audit for V7 CaseA."
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--previous-restart-in", default=None)
    parser.add_argument("--output-dir", default="testModule/v7_caseA_newpyd_long20000_after_tec_heatfix/fluid_energy_audit")
    parser.add_argument("--case-prefix", default="after_tec_heatfix_latest")
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,18"))
    parser.add_argument("--tec-ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,15"))
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--dt-probe", type=float, default=None)
    parser.add_argument("--dt-probe-sweep", default="1e-4,1e-2,0.8,10.0")
    parser.add_argument("--run-regression-tests", action="store_true")
    return parser.parse_args()


def parse_float_list(value: str) -> List[float]:
    out = []
    for item in str(value or "").split(","):
        item = item.strip()
        if item:
            out.append(float(item))
    return out


def finite_float(value: Any) -> Optional[float]:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_load_args(args: argparse.Namespace, restart_path: str) -> argparse.Namespace:
    return argparse.Namespace(
        restart_in=restart_path,
        target_voltage=args.target_voltage,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        pipe_n_nodes=args.pipe_n_nodes,
    )


def load_synced_build(args: argparse.Namespace, restart_path: str) -> Dict[str, Any]:
    build = build_loaded_case(make_load_args(args, restart_path))
    sync_for_audit(build)
    return build


def classify_volume(vol: Any) -> str:
    name = getattr(vol, "name", "")
    if bool(getattr(vol, "is_pressure_boundary", False)):
        return "fixed_pressure_boundary"
    if "Plenum" in name:
        return "plenum"
    if name.startswith("InletPipe") or name.startswith("OutletPipe"):
        return "pipe"
    if name.startswith("Chan_"):
        return "tfe_channel"
    return type(vol).__name__


def build_volume_scale_map(build: Dict[str, Any]) -> Dict[int, float]:
    scale: Dict[int, float] = {}
    for vol in build["system"].fluid_solver.volumes_obj:
        scale[id(vol)] = 1.0

    for vol in build["inlet_pipe_23"].volumes:
        scale[id(vol)] = 2.0
    for vol in build["outlet_pipe"].volumes:
        scale[id(vol)] = 3.0
    for name, channel in build["fluid_channels"].items():
        mult = float(build["ring_multipliers"][name])
        for vol in channel.volumes:
            scale[id(vol)] = mult
    return scale


def volume_state_by_name(build: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
    state: Dict[str, Dict[str, float]] = {}
    for vol in build["system"].fluid_solver.volumes_obj:
        state[getattr(vol, "name", "")] = {
            "T": float(getattr(vol, "T", np.nan)),
            "h": float(getattr(vol, "h", np.nan)),
            "rho": float(getattr(vol, "rho", np.nan)),
            "volume": float(getattr(vol, "vol", np.nan)),
            "mass": float(getattr(vol, "rho", 0.0) * getattr(vol, "vol", 0.0)),
        }
    return state


def collect_volume_sources(build: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, float]]:
    net = build["system"].fluid_solver
    volume_scale = build_volume_scale_map(build)
    rows: List[Dict[str, Any]] = []
    total_source_single = 0.0
    total_source_macro = 0.0
    unexpected_rows = []
    for idx, vol in enumerate(net.volumes_obj):
        q_wall = float(getattr(vol, "Q_wall", 0.0))
        q_vol = float(getattr(vol, "Q_vol", 0.0))
        lam = float(getattr(vol, "implicit_coeff", 0.0))
        temp = float(getattr(vol, "T", np.nan))
        effective = q_wall + q_vol - lam * temp
        scale = float(volume_scale.get(id(vol), 1.0))
        total_source_single += effective
        total_source_macro += effective * scale
        row = {
            "idx": idx,
            "name": getattr(vol, "name", f"vol_{idx}"),
            "type": classify_volume(vol),
            "is_fixed_pressure": bool(getattr(vol, "is_pressure_boundary", False)),
            "representative_scale": scale,
            "Q_wall_w": q_wall,
            "Q_vol_w": q_vol,
            "implicit_coeff_w_per_k": lam,
            "T_k": temp,
            "effective_source_w_single": effective,
            "effective_source_w_macro_scaled": effective * scale,
            "volume_m3": float(getattr(vol, "vol", np.nan)),
            "rho_kg_m3": float(getattr(vol, "rho", np.nan)),
            "h_j_kg": float(getattr(vol, "h", np.nan)),
            "mass_kg_single": float(getattr(vol, "rho", 0.0) * getattr(vol, "vol", 0.0)),
            "mass_kg_macro_scaled": float(getattr(vol, "rho", 0.0) * getattr(vol, "vol", 0.0) * scale),
        }
        rows.append(row)
        if abs(effective) > 1.0e-7 and row["type"] not in {"tfe_channel"}:
            unexpected_rows.append(row["name"])
    summary = {
        "sum_effective_source_w_single": total_source_single,
        "sum_effective_source_w_macro_scaled": total_source_macro,
        "unexpected_nonzero_source_count": float(len(unexpected_rows)),
    }
    return rows, {**summary, "unexpected_nonzero_source_names": unexpected_rows}


def energy_matrix_multiplier_maps(net: HydraulicNetwork) -> Tuple[Dict[int, float], Dict[int, float]]:
    from_map = {
        int(j_idx): float(mult)
        for j_idx, mult in zip(np.asarray(net.energy_from_junc_idx, dtype=int), np.asarray(net.energy_from_multiplier, dtype=float))
    }
    to_map = {
        int(j_idx): float(mult)
        for j_idx, mult in zip(np.asarray(net.energy_to_junc_idx, dtype=int), np.asarray(net.energy_to_multiplier, dtype=float))
    }
    return from_map, to_map


def enthalpy_junction_rows(build: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, float], Dict[str, float]]:
    net = build["system"].fluid_solver
    volume_scale = build_volume_scale_map(build)
    adv_single = {getattr(vol, "name", f"vol_{idx}"): 0.0 for idx, vol in enumerate(net.volumes_obj)}
    adv_macro_scaled = dict(adv_single)
    rows: List[Dict[str, Any]] = []
    for j_idx, (junc, idx_from, idx_to) in enumerate(net.junction_descriptors):
        from_vol = net.volumes_obj[idx_from]
        to_vol = net.volumes_obj[idx_to]
        w = float(junc.W)
        donor = from_vol if w >= 0.0 else to_vol
        h_donor = float(getattr(donor, "h", np.nan))
        mult_from = float(getattr(junc, "multiplier_from", 1.0))
        mult_to = float(getattr(junc, "multiplier_to", 1.0))
        from_adv = -mult_from * w * h_donor
        to_adv = mult_to * w * h_donor
        from_name = getattr(from_vol, "name", f"vol_{idx_from}")
        to_name = getattr(to_vol, "name", f"vol_{idx_to}")
        adv_single[from_name] += from_adv
        adv_single[to_name] += to_adv
        adv_macro_scaled[from_name] += from_adv * float(volume_scale.get(id(from_vol), 1.0))
        adv_macro_scaled[to_name] += to_adv * float(volume_scale.get(id(to_vol), 1.0))
        rows.append({
            "junction_idx": j_idx,
            "junction_name": getattr(junc, "name", f"junc_{j_idx}"),
            "junction_type": type(junc).__name__,
            "from_vol": from_name,
            "to_vol": to_name,
            "donor_vol": getattr(donor, "name", ""),
            "W_single_kg_s": w,
            "multiplier_from": mult_from,
            "multiplier_to": mult_to,
            "from_volume_scale": float(volume_scale.get(id(from_vol), 1.0)),
            "to_volume_scale": float(volume_scale.get(id(to_vol), 1.0)),
            "h_donor_j_kg": h_donor,
            "single_enthalpy_flux_w": w * h_donor,
            "from_row_advective_w_single": from_adv,
            "to_row_advective_w_single": to_adv,
            "from_row_advective_w_macro_scaled": from_adv * float(volume_scale.get(id(from_vol), 1.0)),
            "to_row_advective_w_macro_scaled": to_adv * float(volume_scale.get(id(to_vol), 1.0)),
        })
    return rows, adv_single, adv_macro_scaled


def find_junction_by_name(net: HydraulicNetwork, name: str) -> Optional[Any]:
    for junc in net.junctions_obj:
        if getattr(junc, "name", "") == name:
            return junc
    return None


def collect_control_surface_enthalpy(build: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    inlet_plenum = build["inlet_plenum"]
    outlet_plenum = build["outlet_plenum"]
    net = build["system"].fluid_solver
    rows: List[Dict[str, Any]] = []
    total = {
        "channel_vol01_to_vol37_pickup_w": 0.0,
        "inlet_plenum_to_channel_last_pickup_w": 0.0,
        "channel_first_to_outlet_plenum_pickup_w": 0.0,
        "inlet_plenum_to_outlet_plenum_pickup_w": 0.0,
        "inlet_plenum_to_channel_first_gap_w": 0.0,
        "channel_last_to_outlet_plenum_gap_w": 0.0,
        "merged_fluid_source_w": 0.0,
        "net_mass_flow_kg_s": 0.0,
    }

    for name, channel in build["fluid_channels"].items():
        mult = float(build["ring_multipliers"][name])
        j_in = find_junction_by_name(net, f"J_PlenumIn_{name}")
        j_out = find_junction_by_name(net, f"J_PlenumOut_{name}")
        if j_in is not None:
            w_single = float(j_in.W)
        elif channel.internal_junctions:
            w_single = float(channel.internal_junctions[0].W)
        elif j_out is not None:
            w_single = float(j_out.W)
        else:
            w_single = 0.0
        w_scaled = w_single * mult

        h_inlet_plenum = float(inlet_plenum.h)
        h_channel_first = float(channel.volumes[0].h)
        h_channel_last = float(channel.volumes[-1].h)
        h_outlet_plenum = float(outlet_plenum.h)
        source = float(sum(
            float(getattr(vol, "Q_wall", 0.0))
            + float(getattr(vol, "Q_vol", 0.0))
            - float(getattr(vol, "implicit_coeff", 0.0)) * float(getattr(vol, "T", 0.0))
            for vol in channel.volumes
        ) * mult)

        row = {
            "row_type": "ring",
            "representative_tfe": name,
            "thermal_multiplier": mult,
            "channel_first_volume": getattr(channel.volumes[0], "name", ""),
            "channel_last_volume": getattr(channel.volumes[-1], "name", ""),
            "single_flow_kg_s": w_single,
            "scaled_flow_kg_s": w_scaled,
            "inlet_plenum_h_j_kg": h_inlet_plenum,
            "channel_first_h_j_kg": h_channel_first,
            "channel_last_h_j_kg": h_channel_last,
            "outlet_plenum_h_j_kg": h_outlet_plenum,
            "channel_vol01_to_vol37_pickup_w": w_scaled * (h_channel_last - h_channel_first),
            "inlet_plenum_to_channel_last_pickup_w": w_scaled * (h_channel_last - h_inlet_plenum),
            "channel_first_to_outlet_plenum_pickup_w": w_scaled * (h_outlet_plenum - h_channel_first),
            "inlet_plenum_to_outlet_plenum_pickup_w": w_scaled * (h_outlet_plenum - h_inlet_plenum),
            "inlet_plenum_to_channel_first_gap_w": w_scaled * (h_channel_first - h_inlet_plenum),
            "channel_last_to_outlet_plenum_gap_w": w_scaled * (h_outlet_plenum - h_channel_last),
            "merged_fluid_source_w": source,
        }
        row["channel_vol01_to_vol37_minus_source_w"] = row["channel_vol01_to_vol37_pickup_w"] - source
        row["inlet_plenum_to_outlet_plenum_minus_source_w"] = row["inlet_plenum_to_outlet_plenum_pickup_w"] - source
        rows.append(row)

        for key in total:
            if key == "net_mass_flow_kg_s":
                total[key] += w_scaled
            else:
                total[key] += float(row[key])

    total_row = {"row_type": "total", "representative_tfe": "all", **total}
    total_row["channel_vol01_to_vol37_minus_source_w"] = (
        total["channel_vol01_to_vol37_pickup_w"] - total["merged_fluid_source_w"]
    )
    total_row["inlet_plenum_to_outlet_plenum_minus_source_w"] = (
        total["inlet_plenum_to_outlet_plenum_pickup_w"] - total["merged_fluid_source_w"]
    )
    rows.append(total_row)

    summary = {
        **total,
        "channel_vol01_to_vol37_minus_source_w": total_row["channel_vol01_to_vol37_minus_source_w"],
        "inlet_plenum_to_outlet_plenum_minus_source_w": total_row["inlet_plenum_to_outlet_plenum_minus_source_w"],
        "gap_sum_w": (
            total["inlet_plenum_to_channel_first_gap_w"]
            + total["channel_last_to_outlet_plenum_gap_w"]
        ),
        "control_surface_identity_w": (
            total["inlet_plenum_to_outlet_plenum_pickup_w"]
            - total["channel_vol01_to_vol37_pickup_w"]
            - total["inlet_plenum_to_channel_first_gap_w"]
            - total["channel_last_to_outlet_plenum_gap_w"]
        ),
    }
    return rows, summary


def advective_enthalpy_by_volume(
    net: HydraulicNetwork,
    h_vec: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    adv = np.zeros(net.n_vol)
    flux = np.zeros(net.n_junc)
    donor_idx = np.zeros(net.n_junc, dtype=np.int32)
    for j_idx, (_, idx_from, idx_to) in enumerate(net.junction_descriptors):
        w = float(net.W_vec[j_idx])
        donor = idx_from if w >= 0.0 else idx_to
        donor_idx[j_idx] = donor
        flux[j_idx] = w * float(h_vec[donor])
        adv[idx_from] -= float(net.M_from_vec[j_idx]) * flux[j_idx]
        adv[idx_to] += float(net.M_to_vec[j_idx]) * flux[j_idx]
    return adv, flux, donor_idx


def collect_energy_matrix_probe(
    build: Dict[str, Any],
    dt: float,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    net = build["system"].fluid_solver
    scale_map = build_volume_scale_map(build)
    sources = capture_sources(net)
    net.save_state()
    try:
        net._update_fluid_properties()
        h_old = net.h_vec.copy()
        T_old = net.T_vec.copy()
        mass_old = (net.rho_vec * net.V_vec).copy()
        energy_before_single = total_fluid_enthalpy_energy(build, macro_scaled=False)
        energy_before_macro = total_fluid_enthalpy_energy(build, macro_scaled=True)
        adv_old, _, _ = advective_enthalpy_by_volume(net, h_old)

        net._step_energy_implicit(float(dt))
        h_new = net.h_vec.copy()
        T_new = net.T_vec.copy()
        A = net.energy_matrix.copy()
        B = net.energy_rhs_buffer.copy()
        residual = np.asarray(A.dot(h_new) - B, dtype=float)
        adv_new, _, _ = advective_enthalpy_by_volume(net, h_new)
        energy_after_single = float(np.sum(mass_old * h_new))
        energy_after_macro = 0.0

        rows: List[Dict[str, Any]] = []
        by_type: Dict[str, Dict[str, float]] = {}
        finite_by_type: Dict[str, Dict[str, float]] = {}
        max_abs_residual = 0.0
        max_abs_residual_macro = 0.0
        finite_sum = 0.0
        finite_sum_macro = 0.0
        finite_storage_sum = 0.0
        finite_storage_sum_macro = 0.0
        finite_source_old_sum = 0.0
        finite_source_old_sum_macro = 0.0
        finite_adv_old_sum = 0.0
        finite_adv_old_sum_macro = 0.0
        finite_adv_new_sum = 0.0
        finite_adv_new_sum_macro = 0.0

        for idx, vol in enumerate(net.volumes_obj):
            name = getattr(vol, "name", f"vol_{idx}")
            vol_type = classify_volume(vol)
            is_fixed = bool(getattr(vol, "is_pressure_boundary", False))
            scale = float(scale_map.get(id(vol), 1.0))
            q_eff_old = (
                float(getattr(vol, "Q_wall", 0.0))
                + float(getattr(vol, "Q_vol", 0.0))
                - float(getattr(vol, "implicit_coeff", 0.0)) * float(T_old[idx])
            )
            storage = float(mass_old[idx] * (h_new[idx] - h_old[idx]) / float(dt))
            row_residual = float(residual[idx])
            row_residual_macro = row_residual * scale
            energy_after_macro += float(mass_old[idx] * h_new[idx] * scale)
            max_abs_residual = max(max_abs_residual, abs(row_residual))
            max_abs_residual_macro = max(max_abs_residual_macro, abs(row_residual_macro))
            rows.append({
                "dt_s": float(dt),
                "idx": idx,
                "name": name,
                "type": vol_type,
                "is_fixed_pressure": is_fixed,
                "representative_scale": scale,
                "mass_kg_single": float(mass_old[idx]),
                "h_old_j_kg": float(h_old[idx]),
                "h_new_j_kg": float(h_new[idx]),
                "T_old_k": float(T_old[idx]),
                "T_new_k": float(T_new[idx]),
                "storage_w_single": storage,
                "old_enthalpy_advection_w_single": float(adv_old[idx]),
                "new_enthalpy_advection_w_single": float(adv_new[idx]),
                "effective_source_old_w_single": q_eff_old,
                "matrix_residual_w": row_residual,
                "matrix_residual_w_macro_scaled": row_residual_macro,
            })
            target = by_type.setdefault(vol_type, {
                "count": 0.0,
                "sum_matrix_residual_w": 0.0,
                "sum_abs_matrix_residual_w": 0.0,
                "max_abs_matrix_residual_w": 0.0,
                "sum_matrix_residual_w_macro_scaled": 0.0,
            })
            target["count"] += 1.0
            target["sum_matrix_residual_w"] += row_residual
            target["sum_abs_matrix_residual_w"] += abs(row_residual)
            target["max_abs_matrix_residual_w"] = max(target["max_abs_matrix_residual_w"], abs(row_residual))
            target["sum_matrix_residual_w_macro_scaled"] += row_residual_macro
            if not is_fixed:
                finite_sum += row_residual
                finite_sum_macro += row_residual_macro
                finite_storage_sum += storage
                finite_storage_sum_macro += storage * scale
                finite_source_old_sum += q_eff_old
                finite_source_old_sum_macro += q_eff_old * scale
                finite_adv_old_sum += float(adv_old[idx])
                finite_adv_old_sum_macro += float(adv_old[idx]) * scale
                finite_adv_new_sum += float(adv_new[idx])
                finite_adv_new_sum_macro += float(adv_new[idx]) * scale
                finite_target = finite_by_type.setdefault(vol_type, {
                    "count": 0.0,
                    "sum_matrix_residual_w": 0.0,
                    "sum_abs_matrix_residual_w": 0.0,
                    "max_abs_matrix_residual_w": 0.0,
                    "sum_matrix_residual_w_macro_scaled": 0.0,
                })
                finite_target["count"] += 1.0
                finite_target["sum_matrix_residual_w"] += row_residual
                finite_target["sum_abs_matrix_residual_w"] += abs(row_residual)
                finite_target["max_abs_matrix_residual_w"] = max(
                    finite_target["max_abs_matrix_residual_w"], abs(row_residual)
                )
                finite_target["sum_matrix_residual_w_macro_scaled"] += row_residual_macro

        summary = {
            "dt_s": float(dt),
            "max_abs_matrix_residual_w": max_abs_residual,
            "max_abs_matrix_residual_w_macro_scaled": max_abs_residual_macro,
            "sum_matrix_residual_w": float(np.sum(residual)),
            "sum_finite_nonfixed_matrix_residual_w": finite_sum,
            "sum_finite_nonfixed_matrix_residual_w_macro_scaled": finite_sum_macro,
            "sum_finite_nonfixed_storage_w_single": finite_storage_sum,
            "sum_finite_nonfixed_storage_w_macro_scaled": finite_storage_sum_macro,
            "sum_finite_nonfixed_source_old_w_single": finite_source_old_sum,
            "sum_finite_nonfixed_source_old_w_macro_scaled": finite_source_old_sum_macro,
            "sum_finite_nonfixed_old_advection_w_single": finite_adv_old_sum,
            "sum_finite_nonfixed_old_advection_w_macro_scaled": finite_adv_old_sum_macro,
            "sum_finite_nonfixed_new_advection_w_single": finite_adv_new_sum,
            "sum_finite_nonfixed_new_advection_w_macro_scaled": finite_adv_new_sum_macro,
            "by_type": by_type,
            "finite_nonfixed_by_type": finite_by_type,
            "d_fluid_enthalpy_energy_w_single": (energy_after_single - energy_before_single) / float(dt),
            "d_fluid_enthalpy_energy_w_macro_scaled": (energy_after_macro - energy_before_macro) / float(dt),
            "max_abs_temperature_change_k": float(np.max(np.abs(T_new - T_old))),
            "max_abs_enthalpy_change_j_kg": float(np.max(np.abs(h_new - h_old))),
        }
        return rows, summary
    finally:
        net.load_state()
        restore_sources(net, sources)


def collect_dt_probe_sweep(build: Dict[str, Any], dts: Iterable[float]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    all_rows: List[Dict[str, Any]] = []
    sweep_rows: List[Dict[str, Any]] = []
    for dt in dts:
        rows, summary = collect_energy_matrix_probe(build, float(dt))
        all_rows.extend(rows)
        sweep_rows.append({
            "dt_s": float(dt),
            "d_fluid_enthalpy_energy_w_single": summary["d_fluid_enthalpy_energy_w_single"],
            "d_fluid_enthalpy_energy_w_macro_scaled": summary["d_fluid_enthalpy_energy_w_macro_scaled"],
            "max_abs_temperature_change_k": summary["max_abs_temperature_change_k"],
            "max_abs_enthalpy_change_j_kg": summary["max_abs_enthalpy_change_j_kg"],
            "max_abs_matrix_residual_w": summary["max_abs_matrix_residual_w"],
            "max_abs_matrix_residual_w_macro_scaled": summary["max_abs_matrix_residual_w_macro_scaled"],
            "sum_finite_nonfixed_matrix_residual_w": summary["sum_finite_nonfixed_matrix_residual_w"],
            "sum_finite_nonfixed_matrix_residual_w_macro_scaled": summary[
                "sum_finite_nonfixed_matrix_residual_w_macro_scaled"
            ],
            "sum_finite_nonfixed_storage_w_single": summary["sum_finite_nonfixed_storage_w_single"],
            "sum_finite_nonfixed_storage_w_macro_scaled": summary["sum_finite_nonfixed_storage_w_macro_scaled"],
            "sum_finite_nonfixed_source_old_w_single": summary["sum_finite_nonfixed_source_old_w_single"],
            "sum_finite_nonfixed_source_old_w_macro_scaled": summary["sum_finite_nonfixed_source_old_w_macro_scaled"],
            "sum_finite_nonfixed_old_advection_w_single": summary["sum_finite_nonfixed_old_advection_w_single"],
            "sum_finite_nonfixed_old_advection_w_macro_scaled": summary[
                "sum_finite_nonfixed_old_advection_w_macro_scaled"
            ],
            "sum_finite_nonfixed_new_advection_w_single": summary["sum_finite_nonfixed_new_advection_w_single"],
            "sum_finite_nonfixed_new_advection_w_macro_scaled": summary[
                "sum_finite_nonfixed_new_advection_w_macro_scaled"
            ],
        })
    return all_rows, sweep_rows


def collect_volume_energy_balance(
    build: Dict[str, Any],
    previous_build: Optional[Dict[str, Any]],
    adv_single: Dict[str, float],
    adv_macro_scaled: Dict[str, float],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    current_time = float(build["system"].global_time)
    previous_time = None
    previous_state = None
    dt_storage = None
    if previous_build is not None:
        previous_time = float(previous_build["system"].global_time)
        dt_storage = current_time - previous_time
        previous_state = volume_state_by_name(previous_build)
        if dt_storage <= 0.0:
            dt_storage = None

    rows: List[Dict[str, Any]] = []
    volume_scale = build_volume_scale_map(build)
    focus_names = {
        "Global_Inlet_Plenum",
        "Global_Outlet_Plenum",
        "Global_Inlet_Boundary",
        "Global_Outlet_Boundary",
    }
    totals = {
        "storage_w_single": 0.0,
        "net_advective_enthalpy_w_single": 0.0,
        "source_w_single": 0.0,
        "residual_w_single": 0.0,
        "storage_w_macro_scaled": 0.0,
        "net_advective_enthalpy_w_macro_scaled": 0.0,
        "source_w_macro_scaled": 0.0,
        "residual_w_macro_scaled": 0.0,
    }
    finite_totals = {key: 0.0 for key in totals}
    suspicious = []

    for idx, vol in enumerate(build["system"].fluid_solver.volumes_obj):
        name = getattr(vol, "name", f"vol_{idx}")
        scale = float(volume_scale.get(id(vol), 1.0))
        q_eff = (
            float(getattr(vol, "Q_wall", 0.0))
            + float(getattr(vol, "Q_vol", 0.0))
            - float(getattr(vol, "implicit_coeff", 0.0)) * float(getattr(vol, "T", 0.0))
        )
        mass = float(getattr(vol, "rho", 0.0) * getattr(vol, "vol", 0.0))
        is_fixed = bool(getattr(vol, "is_pressure_boundary", False))
        storage = None
        if previous_state is not None and dt_storage is not None and name in previous_state:
            storage = mass * (float(getattr(vol, "h", np.nan)) - previous_state[name]["h"]) / dt_storage
        adv = float(adv_single.get(name, 0.0))
        residual = None if storage is None or is_fixed else storage - adv - q_eff
        storage_macro = None if storage is None else storage * scale
        adv_macro = float(adv_macro_scaled.get(name, 0.0))
        q_macro = q_eff * scale
        residual_macro = None if storage_macro is None or is_fixed else storage_macro - adv_macro - q_macro

        row = {
            "idx": idx,
            "name": name,
            "type": classify_volume(vol),
            "is_focus_volume": bool(
                name in focus_names
                or name.startswith("InletPipe1")
                or name.startswith("InletPipe23Rep")
                or name.startswith("OutletPipeRep")
                or name.startswith("Chan_")
            ),
            "is_fixed_pressure": is_fixed,
            "representative_scale": scale,
            "T_k": float(getattr(vol, "T", np.nan)),
            "h_j_kg": float(getattr(vol, "h", np.nan)),
            "mass_kg_single": mass,
            "source_w_single": q_eff,
            "net_advective_enthalpy_w_single": adv,
            "storage_w_single": storage,
            "residual_w_single": residual,
            "source_w_macro_scaled": q_macro,
            "net_advective_enthalpy_w_macro_scaled": adv_macro,
            "storage_w_macro_scaled": storage_macro,
            "residual_w_macro_scaled": residual_macro,
        }
        rows.append(row)

        for key, value in (
            ("source_w_single", q_eff),
            ("net_advective_enthalpy_w_single", adv),
            ("source_w_macro_scaled", q_macro),
            ("net_advective_enthalpy_w_macro_scaled", adv_macro),
        ):
            totals[key] += value
            if not is_fixed:
                finite_totals[key] += value
        if storage is not None:
            totals["storage_w_single"] += storage
            totals["storage_w_macro_scaled"] += float(storage_macro)
            if not is_fixed:
                finite_totals["storage_w_single"] += storage
                finite_totals["storage_w_macro_scaled"] += float(storage_macro)
                finite_totals["residual_w_single"] += float(residual)
                finite_totals["residual_w_macro_scaled"] += float(residual_macro)
        if (
            storage is not None
            and abs(q_eff) < 1.0e-7
            and residual is not None
            and residual > 1.0
            and not bool(getattr(vol, "is_pressure_boundary", False))
        ):
            suspicious.append({"name": name, "residual_w_single": residual})

    summary = {
        "current_time_s": current_time,
        "previous_time_s": previous_time,
        "storage_dt_s": dt_storage,
        "totals_including_fixed_boundaries": totals,
        "totals_finite_nonfixed_volumes": finite_totals,
        "suspicious_positive_no_source_residuals": suspicious,
    }
    return rows, summary


def collect_macro_multiplier_audit(build: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    net = build["system"].fluid_solver
    from_matrix, to_matrix = energy_matrix_multiplier_maps(net)
    rows: List[Dict[str, Any]] = []
    mismatches: List[Dict[str, Any]] = []
    for j_idx, (junc, idx_from, idx_to) in enumerate(net.junction_descriptors):
        if not isinstance(junc, MacroFlowJunction):
            continue
        expected_from = float(getattr(junc, "multiplier_from", 1.0))
        expected_to = float(getattr(junc, "multiplier_to", 1.0))
        matrix_from = from_matrix.get(j_idx)
        matrix_to = to_matrix.get(j_idx)
        row = {
            "junction_idx": j_idx,
            "junction_name": getattr(junc, "name", f"junc_{j_idx}"),
            "from_vol": getattr(junc.from_vol, "name", ""),
            "to_vol": getattr(junc.to_vol, "name", ""),
            "macro_vol": getattr(junc.macro_vol, "name", ""),
            "W_single_kg_s": float(junc.W),
            "multiplier_from": expected_from,
            "multiplier_to": expected_to,
            "energy_matrix_multiplier_from": matrix_from,
            "energy_matrix_multiplier_to": matrix_to,
            "from_row_is_fixed_pressure": bool(getattr(junc.from_vol, "is_pressure_boundary", False)),
            "to_row_is_fixed_pressure": bool(getattr(junc.to_vol, "is_pressure_boundary", False)),
            "matches_energy_matrix_from": matrix_from is None or math.isclose(matrix_from, expected_from, rel_tol=0.0, abs_tol=1.0e-12),
            "matches_energy_matrix_to": matrix_to is None or math.isclose(matrix_to, expected_to, rel_tol=0.0, abs_tol=1.0e-12),
        }
        rows.append(row)
        if not row["matches_energy_matrix_from"] or not row["matches_energy_matrix_to"]:
            mismatches.append(row)
    return rows, {"macro_junction_count": len(rows), "matrix_multiplier_mismatches": mismatches}


def capture_sources(net: HydraulicNetwork) -> List[Tuple[float, float, float]]:
    return [
        (
            float(getattr(vol, "Q_wall", 0.0)),
            float(getattr(vol, "Q_vol", 0.0)),
            float(getattr(vol, "implicit_coeff", 0.0)),
        )
        for vol in net.volumes_obj
    ]


def restore_sources(net: HydraulicNetwork, sources: List[Tuple[float, float, float]]) -> None:
    for vol, (q_wall, q_vol, lam) in zip(net.volumes_obj, sources):
        vol.Q_wall = q_wall
        vol.Q_vol = q_vol
        vol.implicit_coeff = lam


def total_fluid_enthalpy_energy(build: Dict[str, Any], macro_scaled: bool) -> float:
    scale_map = build_volume_scale_map(build)
    total = 0.0
    for vol in build["system"].fluid_solver.volumes_obj:
        scale = float(scale_map.get(id(vol), 1.0)) if macro_scaled else 1.0
        total += float(getattr(vol, "rho", 0.0) * getattr(vol, "vol", 0.0) * getattr(vol, "h", 0.0) * scale)
    return total


def run_dt_probe(build: Dict[str, Any], dt: float) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    sources = capture_sources(net)
    energy_before_single = total_fluid_enthalpy_energy(build, macro_scaled=False)
    energy_before_macro = total_fluid_enthalpy_energy(build, macro_scaled=True)
    net.save_state()
    try:
        net._update_fluid_properties()
        net._step_energy_implicit(float(dt))
        net._sync_vectors_to_objects(sync_pressure=False, sync_flow=False, sync_energy=True, sync_properties=False)
        energy_after_single = total_fluid_enthalpy_energy(build, macro_scaled=False)
        energy_after_macro = total_fluid_enthalpy_energy(build, macro_scaled=True)
        return {
            "dt_s": float(dt),
            "d_fluid_enthalpy_energy_w_single": (energy_after_single - energy_before_single) / float(dt),
            "d_fluid_enthalpy_energy_w_macro_scaled": (energy_after_macro - energy_before_macro) / float(dt),
            "max_abs_temperature_change_k": float(np.max(np.abs(net.T_vec - net.T_backup))),
            "max_abs_enthalpy_change_j_kg": float(np.max(np.abs(net.h_vec - net.h_backup))),
        }
    finally:
        net.load_state()
        restore_sources(net, sources)


def set_network_flow(net: HydraulicNetwork, flow: float) -> None:
    for idx, junc in enumerate(net.junctions_obj):
        junc.W = float(flow)
        net.W_vec[idx] = float(flow)
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()


def make_pipe_network(
    internal_temperatures: Iterable[float],
    inlet_temperature: float,
    outlet_temperature: float,
    flow: float,
    macro_multiplier: Optional[int] = None,
) -> Tuple[HydraulicNetwork, List[Any]]:
    sodium = Sodium()
    area = 1.0e-4
    dh = 0.01
    inlet = IncompressibleBoundaryVolume("Test_Inlet", sodium, P=1.2e5, T=inlet_temperature, flow_area=area, hydraulic_diam=dh)
    outlet = IncompressibleBoundaryVolume("Test_Outlet", sodium, P=1.0e5, T=outlet_temperature, flow_area=area, hydraulic_diam=dh)
    inlet.is_pressure_boundary = True
    outlet.is_pressure_boundary = True
    vols = [inlet]
    internal = []
    for idx, temp in enumerate(internal_temperatures):
        vol = IncompressibleFluidVolume(
            f"Test_Vol_{idx + 1}",
            volume=1.0e-4,
            length=0.1,
            flow_area=area,
            hydraulic_diam=dh,
            initial_P=1.1e5,
            initial_T=float(temp),
            material=sodium,
        )
        internal.append(vol)
        vols.append(vol)
    vols.append(outlet)
    juncs = []
    for idx in range(len(vols) - 1):
        if macro_multiplier is not None and idx == len(vols) - 2:
            junc = MacroFlowJunction(
                f"Test_J_{idx}",
                vols[idx],
                vols[idx + 1],
                macro_vol=vols[idx],
                multiplier=macro_multiplier,
                flow_area=area,
            )
        else:
            junc = FlowJunction(f"Test_J_{idx}", vols[idx], vols[idx + 1], flow_area=area)
        junc.W = float(flow)
        juncs.append(junc)
    net = HydraulicNetwork(vols, juncs, gravity_vector=0.0)
    set_network_flow(net, flow)
    return net, internal


def regression_total_energy_rate(net: HydraulicNetwork, before_h: np.ndarray, dt: float) -> float:
    mass = np.asarray(net.rho_vec, dtype=float) * np.asarray(net.V_vec, dtype=float)
    return float(np.sum(mass * (net.h_vec - before_h)) / float(dt))


def run_regression_tests() -> Dict[str, Any]:
    results: Dict[str, Any] = {}

    dt = 0.05
    net, internal = make_pipe_network([743.0, 743.0, 743.0], 743.0, 743.0, flow=0.05)
    h0 = net.h_vec.copy()
    t0 = net.T_vec.copy()
    net._step_energy_implicit(dt)
    results["A_same_temperature_no_source"] = {
        "max_abs_delta_T_k": float(np.max(np.abs(net.T_vec - t0))),
        "total_energy_rate_w": regression_total_energy_rate(net, h0, dt),
    }

    net, internal = make_pipe_network([743.0, 743.0, 743.0], 800.0, 700.0, flow=0.05)
    h0 = net.h_vec.copy()
    net._step_energy_implicit(dt)
    results["B_hot_inlet_no_source"] = {
        "total_energy_rate_w": regression_total_energy_rate(net, h0, dt),
        "internal_mean_temperature_k": float(np.mean([vol.T for vol in internal])),
    }

    sodium = Sodium()
    area = 1.0e-4
    dh = 0.01
    inlet = IncompressibleBoundaryVolume("Macro_Inlet", sodium, P=1.2e5, T=743.0, flow_area=area, hydraulic_diam=dh)
    outlet = IncompressibleBoundaryVolume("Macro_Outlet", sodium, P=1.0e5, T=743.0, flow_area=area, hydraulic_diam=dh)
    inlet.is_pressure_boundary = True
    outlet.is_pressure_boundary = True
    plenum = IncompressibleFluidVolume("Macro_Plenum", 1.0e-4, 0.1, area, dh, 1.1e5, 743.0, sodium)
    pipe = IncompressibleFluidVolume("Macro_RepPipe", 1.0e-4, 0.1, area, dh, 1.1e5, 743.0, sodium)
    j_in = MacroFlowJunction("Macro_Inlet_to_Plenum", inlet, plenum, macro_vol=plenum, multiplier=3, flow_area=area)
    j_out = MacroFlowJunction("Macro_Plenum_to_RepPipe", plenum, pipe, macro_vol=plenum, multiplier=3, flow_area=area)
    j_bound = FlowJunction("Macro_RepPipe_to_Outlet", pipe, outlet, flow_area=area)
    net = HydraulicNetwork([inlet, plenum, pipe, outlet], [j_in, j_out, j_bound], gravity_vector=0.0)
    set_network_flow(net, 0.05)
    h0 = net.h_vec.copy()
    t0 = net.T_vec.copy()
    net._step_energy_implicit(dt)
    results["C_macro_same_temperature_no_source"] = {
        "max_abs_delta_T_k": float(np.max(np.abs(net.T_vec - t0))),
        "total_energy_rate_w": regression_total_energy_rate(net, h0, dt),
    }

    vol = IncompressibleFluidVolume(
        "Heated_Box",
        volume=1.0e-4,
        length=0.1,
        flow_area=1.0e-4,
        hydraulic_diam=0.01,
        initial_P=1.0e5,
        initial_T=743.0,
        material=sodium,
    )
    net = HydraulicNetwork([vol], [], gravity_vector=0.0)
    vol.Q_wall = 100.0
    net.Q_expl_vec[:] = 100.0
    h0 = net.h_vec.copy()
    net._step_energy_implicit(dt)
    energy_rate = regression_total_energy_rate(net, h0, dt)
    results["D_single_volume_fixed_source"] = {
        "expected_source_w": 100.0,
        "total_energy_rate_w": energy_rate,
        "error_w": energy_rate - 100.0,
    }

    inlet_plenum = IncompressibleFluidVolume("Mini_Inlet_Plenum", 2.0e-4, 0.1, area, dh, 1.1e5, 730.0, sodium)
    channel_first = IncompressibleFluidVolume("Mini_Channel_Vol01", 1.0e-4, 0.1, area, dh, 1.1e5, 743.0, sodium)
    channel_last = IncompressibleFluidVolume("Mini_Channel_Vol02", 1.0e-4, 0.1, area, dh, 1.1e5, 750.0, sodium)
    outlet_plenum = IncompressibleFluidVolume("Mini_Outlet_Plenum", 2.0e-4, 0.1, area, dh, 1.1e5, 760.0, sodium)
    j0 = FlowJunction("Mini_PlenumIn", inlet_plenum, channel_first, flow_area=area)
    j1 = FlowJunction("Mini_ChannelInternal", channel_first, channel_last, flow_area=area)
    j2 = FlowJunction("Mini_PlenumOut", channel_last, outlet_plenum, flow_area=area)
    net = HydraulicNetwork([inlet_plenum, channel_first, channel_last, outlet_plenum], [j0, j1, j2], gravity_vector=0.0)
    set_network_flow(net, 0.05)
    h0 = net.h_vec.copy()
    net._step_energy_implicit(dt)
    w = 0.05
    channel_pickup = w * (float(channel_last.h) - float(channel_first.h))
    plenum_to_plenum_pickup = w * (float(outlet_plenum.h) - float(inlet_plenum.h))
    results["E_plenum_channel_plenum_control_surface_offset"] = {
        "initial_channel_pickup_w": channel_pickup,
        "initial_plenum_to_plenum_pickup_w": plenum_to_plenum_pickup,
        "initial_control_surface_gap_w": plenum_to_plenum_pickup - channel_pickup,
        "implicit_step_total_energy_rate_w": regression_total_energy_rate(net, h0, dt),
    }
    return results


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    build = load_synced_build(args, args.restart_in)
    previous_build = load_synced_build(args, args.previous_restart_in) if args.previous_restart_in else None

    volume_source_rows, volume_source_summary = collect_volume_sources(build)
    junction_rows, adv_single, adv_macro_scaled = enthalpy_junction_rows(build)
    balance_rows, balance_summary = collect_volume_energy_balance(build, previous_build, adv_single, adv_macro_scaled)
    macro_rows, macro_summary = collect_macro_multiplier_audit(build)
    control_rows, control_summary = collect_control_surface_enthalpy(build)

    dt_sweep = parse_float_list(args.dt_probe_sweep)
    if args.dt_probe is not None and not any(math.isclose(args.dt_probe, item, rel_tol=0.0, abs_tol=1.0e-15) for item in dt_sweep):
        dt_sweep.insert(0, float(args.dt_probe))
    matrix_rows, dt_probe_sweep_rows = collect_dt_probe_sweep(build, dt_sweep)
    matrix_residual_summary = {
        str(row["dt_s"]): row
        for row in dt_probe_sweep_rows
    }

    source_path = output_dir / f"{args.case_prefix}_fluid_volume_sources.csv"
    junction_path = output_dir / f"{args.case_prefix}_fluid_junction_enthalpy.csv"
    balance_path = output_dir / f"{args.case_prefix}_fluid_volume_energy_balance.csv"
    macro_path = output_dir / f"{args.case_prefix}_macro_multiplier_audit.csv"
    control_path = output_dir / f"{args.case_prefix}_fluid_control_surface_enthalpy.csv"
    matrix_path = output_dir / f"{args.case_prefix}_energy_matrix_residual.csv"
    dt_sweep_path = output_dir / f"{args.case_prefix}_dt_probe_sweep.csv"
    summary_path = output_dir / f"{args.case_prefix}_fluid_energy_summary.json"

    write_csv(source_path, volume_source_rows)
    write_csv(junction_path, junction_rows)
    write_csv(balance_path, balance_rows)
    write_csv(macro_path, macro_rows)
    write_csv(control_path, control_rows)
    write_csv(matrix_path, matrix_rows)
    write_csv(dt_sweep_path, dt_probe_sweep_rows)

    regression_summary = run_regression_tests() if args.run_regression_tests else None

    summary = {
        "restart_in": args.restart_in,
        "previous_restart_in": args.previous_restart_in,
        "time_s": float(build["system"].global_time),
        "target_voltage_v": float(args.target_voltage),
        "ring_multipliers": build["ring_multipliers"],
        "tec_ring_multipliers": build["tec_ring_multipliers"],
        "outputs": {
            "fluid_volume_sources_csv": str(source_path),
            "fluid_junction_enthalpy_csv": str(junction_path),
            "fluid_volume_energy_balance_csv": str(balance_path),
            "macro_multiplier_audit_csv": str(macro_path),
            "fluid_control_surface_enthalpy_csv": str(control_path),
            "energy_matrix_residual_csv": str(matrix_path),
            "dt_probe_sweep_csv": str(dt_sweep_path),
            "fluid_energy_summary_json": str(summary_path),
        },
        "flow": _case_a_flow_diagnostics(build),
        "electric": _case_a_electric_diagnostics(build["core"]),
        "volume_source_summary": volume_source_summary,
        "volume_energy_balance_summary": balance_summary,
        "macro_multiplier_summary": macro_summary,
        "control_surface_closure": control_summary,
        "matrix_residual_summary": matrix_residual_summary,
        "dt_probe_sweep_summary": dt_probe_sweep_rows,
        "regression_tests": regression_summary,
        "notes": {
            "source_definition": "effective_source_w = Q_wall + Q_vol - implicit_coeff*T",
            "net_advective_enthalpy_definition": (
                "Uses the same upwind donor and multiplier_from/to convention as "
                "HydraulicNetwork._step_energy_implicit matrix rows."
            ),
            "storage_definition": "mass_current * (h_current - h_previous) / (t_current - t_previous)",
            "macro_scaled_definition": (
                "Representative pipe/TFE volume rows are multiplied by their physical representative scale; "
                "matrix-side macro multipliers are still reported separately per junction."
            ),
        },
    }
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)

    print(f"Summary: {summary_path}")
    print(f"Volume sources: {source_path}")
    print(f"Junction enthalpy: {junction_path}")
    print(f"Volume balance: {balance_path}")
    print(f"Macro multipliers: {macro_path}")
    print(f"Control surfaces: {control_path}")
    print(f"Energy matrix residual: {matrix_path}")
    print(f"dt probe sweep: {dt_sweep_path}")
    print(json.dumps({
        "time_s": summary["time_s"],
        "source_single_w": volume_source_summary["sum_effective_source_w_single"],
        "source_macro_scaled_w": volume_source_summary["sum_effective_source_w_macro_scaled"],
        "unexpected_nonzero_source_count": len(volume_source_summary["unexpected_nonzero_source_names"]),
        "storage_dt_s": balance_summary["storage_dt_s"],
        "macro_multiplier_mismatch_count": len(macro_summary["matrix_multiplier_mismatches"]),
        "control_surface_closure": control_summary,
        "dt_probe_sweep_summary": dt_probe_sweep_rows,
    }, indent=2, ensure_ascii=False, default=json_default))


if __name__ == "__main__":
    main()
