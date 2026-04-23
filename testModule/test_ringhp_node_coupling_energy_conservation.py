import argparse
import csv
import os
import sys
import time
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Components.RingHP import RingHP
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.BoundaryVolume import (
    IncompressibleBoundaryVolume,
    InletJunction,
)
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidChannel
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager


@dataclass(frozen=True)
class AuditConfig:
    n_nodes: int = 4
    hp_multipliers: Tuple[float, ...] = (1.0, 2.0, 0.0, 3.0)

    inlet_temperature: float = 863.0
    initial_fluid_temperature: float = 863.0
    initial_header_temperature: float = 863.0
    initial_hp_temperature: float = 800.0
    inlet_mass_flow: float = 0.1815
    outlet_pressure: float = 161000.0

    header_flow_area: float = 0.0016065
    header_dh: float = 0.04167
    header_wall_thickness: float = 0.002

    hp_r_out: float = 0.0085
    hp_r_in: float = 0.0081
    hp_r_vapor: float = 0.0075
    hp_l_eva: float = 0.06
    hp_l_con: float = 0.482
    hp_n_eva: int = 1
    hp_n_con: int = 12
    hp_n_wick: int = 1
    hp_n_wall: int = 1
    hp_porosity: float = 0.5

    fin_thickness: float = 0.0003
    fin_height: float = 22.65e-3
    n_fin_height: int = 15
    emissivity: float = 0.93
    up_view_factor: float = 0.0
    down_view_factor: float = 0.3
    space_temperature: float = 3.0

    t_end: float = 500.0
    min_dt: float = 1.0e-3
    max_dt: float = 0.2
    safety_factor: float = 0.5
    inner_iter: int = 2
    print_every: int = 100
    rolling_window_steps: int = 100
    csv_path: Optional[str] = os.path.join(
        current_dir, "ringhp_node_coupling_energy_audit.csv"
    )
    restart_path: Optional[str] = os.path.join(
        current_dir, "ringhp_node_coupling_restart_500s.npz"
    )
    monitor_csv_path: Optional[str] = os.path.join(
        current_dir, "ringhp_node0_failure_monitor.csv"
    )
    monitor_snapshot_path: Optional[str] = os.path.join(
        current_dir, "ringhp_node0_failure_snapshot.npz"
    )


DEFAULT_CONFIG = AuditConfig()


def lyon_martinelli(Re, Pr, p_d_ratio=1.0):
    _ = p_d_ratio
    pe = np.maximum(np.asarray(Re, dtype=float) * np.asarray(Pr, dtype=float), 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


def write_csv(path: str, rows: List[Dict[str, float]]):
    if not path or not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def save_restart(
    path: str,
    config: AuditConfig,
    case,
    rows: List[Dict[str, float]],
    trends: List[Dict[str, float]],
    criteria: Dict[str, float],
):
    if not path or not rows:
        return

    os.makedirs(os.path.dirname(path), exist_ok=True)
    system = case["system"]
    ring_hp = case["ring_hp"]
    channel = case["channel"]
    header_solid = case["header_solid"]
    inlet = case["inlet"]
    outlet = case["outlet"]
    j_in = case["j_in"]
    j_out = case["j_out"]
    final = rows[-1]

    hp_T = np.stack([np.array(hp.hp.T, copy=True) for hp in ring_hp.hp_units])
    hp_current_time = np.array([float(hp.hp.current_time) for hp in ring_hp.hp_units])
    fin_temperature = np.stack(
        [np.array(hp.last_fin_temperature, copy=True) for hp in ring_hp.hp_units]
    )
    fin_radiation = np.stack(
        [
            np.array(hp.last_fin_radiation_distribution, copy=True)
            for hp in ring_hp.hp_units
        ]
    )
    fin_absorption = np.stack(
        [
            np.array(hp.last_fin_absorption_distribution, copy=True)
            for hp in ring_hp.hp_units
        ]
    )
    fin_net_from_root = np.stack(
        [
            np.array(hp.last_fin_net_from_root_distribution, copy=True)
            for hp in ring_hp.hp_units
        ]
    )

    trend_mean_abs = np.array(
        [item["mean_abs_residual_hp_scaled"] for item in trends], dtype=float
    )
    trend_max_rel = np.array(
        [item["max_relative_hp_scaled"] for item in trends], dtype=float
    )

    np.savez_compressed(
        path,
        global_time=np.array([float(system.global_time)]),
        last_dt=np.array([float(getattr(system, "_last_dt", 0.0))]),
        config_n_nodes=np.array([config.n_nodes]),
        config_hp_multipliers=np.asarray(config.hp_multipliers, dtype=float),
        channel_T=np.array([vol.T for vol in channel.volumes], dtype=float),
        channel_P=np.array([vol.P for vol in channel.volumes], dtype=float),
        channel_h=np.array([vol.h for vol in channel.volumes], dtype=float),
        channel_rho=np.array([vol.rho for vol in channel.volumes], dtype=float),
        inlet_T=np.array([float(inlet.T)]),
        inlet_P=np.array([float(inlet.P)]),
        inlet_h=np.array([float(inlet.h)]),
        outlet_T=np.array([float(outlet.T)]),
        outlet_P=np.array([float(outlet.P)]),
        outlet_h=np.array([float(outlet.h)]),
        W_in=np.array([float(j_in.W)]),
        W_out=np.array([float(j_out.W)]),
        internal_W=np.array([float(j.W) for j in channel.internal_junctions]),
        header_T=np.array(header_solid.T, copy=True),
        hp_T=hp_T,
        hp_current_time=hp_current_time,
        fin_temperature=fin_temperature,
        fin_radiation_distribution=fin_radiation,
        fin_absorption_distribution=fin_absorption,
        fin_net_from_root_distribution=fin_net_from_root,
        final_Q_fluid_enthalpy_drop=np.array(
            [float(final["Q_fluid_enthalpy_drop"])]
        ),
        final_Q_hp_eva_scaled_total=np.array(
            [float(final["Q_hp_eva_scaled_total"])]
        ),
        final_Q_rej_net_scaled=np.array([float(final["Q_rej_net_scaled"])]),
        final_dU_hp_scaled_dt=np.array([float(final["dU_hp_scaled_dt"])]),
        final_residual_hp_scaled=np.array([float(final["residual_hp_scaled"])]),
        final_relative_hp_scaled=np.array([float(final["relative_hp_scaled"])]),
        final_T_hp_min=np.array([float(final["T_hp_min"])]),
        trend_mean_abs_residual_hp_scaled=trend_mean_abs,
        trend_max_relative_hp_scaled=trend_max_rel,
        criteria_temperature_valid=np.array(
            [1.0 if criteria.get("temperature_valid", False) else 0.0]
        ),
        criteria_coupling_valid=np.array(
            [1.0 if criteria.get("coupling_valid", False) else 0.0]
        ),
        criteria_final_window_valid=np.array(
            [1.0 if criteria.get("final_window_valid", False) else 0.0]
        ),
        criteria_final_dU_valid=np.array(
            [1.0 if criteria.get("final_dU_valid", False) else 0.0]
        ),
    )
    print(f"Restart written to: {path}")


def save_monitor_snapshot(path: str, case, previous_hp0_T, current_hp0_T, rows):
    if not path or not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    system = case["system"]
    ring_hp = case["ring_hp"]
    channel = case["channel"]
    final = rows[-1]

    np.savez_compressed(
        path,
        global_time=np.array([float(system.global_time)]),
        previous_hp0_T=np.array(previous_hp0_T, copy=True),
        current_hp0_T=np.array(current_hp0_T, copy=True),
        hp_T=np.stack([np.array(hp.hp.T, copy=True) for hp in ring_hp.hp_units]),
        channel_T=np.array([vol.T for vol in channel.volumes], dtype=float),
        channel_h=np.array([vol.h for vol in channel.volumes], dtype=float),
        channel_P=np.array([vol.P for vol in channel.volumes], dtype=float),
        final_step=np.array([int(final["step"])]),
        final_t=np.array([float(final["t"])]),
        final_dt=np.array([float(final["dt"])]),
        final_stop_reason=np.array([str(final["stop_reason"])]),
        final_hp0_T_min=np.array([float(final["hp0_T_min"])]),
        final_hp0_T_max=np.array([float(final["hp0_T_max"])]),
        final_hp0_max_abs_delta_T=np.array(
            [float(final["hp0_max_abs_delta_T"])]
        ),
        final_Q_hp0_eva=np.array([float(final["hp0_Q_eva"])]),
        final_Q_hp0_rej_net=np.array([float(final["hp0_Q_rej_net"])]),
        final_dU_hp0_dt=np.array([float(final["hp0_dU_dt"])]),
        final_residual_hp0=np.array([float(final["hp0_residual"])]),
    )
    print(f"Monitor snapshot written to: {path}")


def solid_energy(solid) -> float:
    if hasattr(solid, "_update_properties"):
        solid._update_properties()
    return float(np.sum(solid.thermal_capacitance * solid.T))


def hp_scaled_energy(ring_hp: RingHP) -> float:
    total = 0.0
    for node_index, hp_pos, multiplier in iter_present_hps(ring_hp):
        _ = node_index
        total += multiplier * solid_energy(ring_hp.hp_units[hp_pos].hp)
    return float(total)


def iter_present_hps(ring_hp: RingHP):
    summary = ring_hp.get_hp_status_summary()
    mask = np.asarray(summary["hp_presence_mask"], dtype=bool)
    multipliers = np.asarray(ring_hp.audit_hp_multipliers, dtype=float)
    hp_pos = 0
    for node_index, present in enumerate(mask):
        if present:
            yield node_index, hp_pos, float(multipliers[node_index])
            hp_pos += 1


def build_ringhp_node_coupling_case(config: AuditConfig):
    if len(config.hp_multipliers) != config.n_nodes:
        raise ValueError("hp_multipliers length must match n_nodes.")

    mat_fluid = Sodium()
    mat_wall = SS316(name="Test5_SS316_Wall")
    mat_hp_fluid = SodiumHP(name="Test5_HP_Fluid_Na")
    mat_wick = WickMaterial(
        name="Test5_HP_Wick_Composite",
        solid_mat=SS316(),
        fluid_mat=mat_hp_fluid,
        porosity=config.hp_porosity,
        r_vapor=config.hp_r_vapor,
        r_in_wall=config.hp_r_in,
    )

    inlet = IncompressibleBoundaryVolume(
        name="Test5_Inlet",
        material=mat_fluid,
        P=config.outlet_pressure + 1000.0,
        T=config.inlet_temperature,
        flow_area=config.header_flow_area,
        hydraulic_diam=config.header_dh,
    )
    outlet = IncompressibleBoundaryVolume(
        name="Test5_Outlet",
        material=mat_fluid,
        P=config.outlet_pressure,
        T=config.initial_fluid_temperature,
        flow_area=config.header_flow_area,
        hydraulic_diam=config.header_dh,
    )
    outlet.is_pressure_boundary = True

    header_length = config.n_nodes * config.hp_l_eva
    channel = IncompressibleFluidChannel(
        name="Test5_Header_Channel",
        n_nodes=config.n_nodes,
        total_length=header_length,
        flow_area=config.header_flow_area,
        hydraulic_diam=config.header_dh,
        initial_P=config.outlet_pressure + 500.0,
        initial_T=config.initial_fluid_temperature,
        material=mat_fluid,
    )

    j_in = InletJunction(
        name="Test5_J_Inlet",
        from_vol=inlet,
        to_vol=channel.volumes[0],
        W_initial=config.inlet_mass_flow,
    )
    j_out = FlowJunction(
        name="Test5_J_Outlet",
        from_vol=channel.volumes[-1],
        to_vol=outlet,
        flow_area=config.header_flow_area,
        k_loss=0.0,
    )

    r_in_header = config.header_dh / 2.0
    header_perimeter = 2.0 * np.pi * r_in_header
    header_mesh = Mesh2D(
        x_dim=config.header_wall_thickness,
        n_x=1,
        y_dim=header_length,
        n_y=config.n_nodes,
        geometry_type="cylindrical",
        inner_radius=r_in_header,
    )
    header_solid = HeatConduction2D(
        mesh=header_mesh,
        material=mat_wall,
        initial_temp=config.initial_header_temperature,
        name="Test5_Header_Wall",
    )
    header_solid.boundaries["right"].add_flux_condition(q_flux=0.0)
    header_solid.boundaries["top"].add_flux_condition(q_flux=0.0)
    header_solid.boundaries["bottom"].add_flux_condition(q_flux=0.0)

    fin_wrap_ratio = (2.0 * config.fin_thickness) / (
        2.0 * np.pi * config.hp_r_out
    )

    ring_hp = RingHP(
        name="Test5_RingHP",
        fluid_channel=channel,
        solid_header=header_solid,
        hp_multipliers=list(config.hp_multipliers),
        header_flow_area=config.header_flow_area,
        header_dh=config.header_dh,
        header_heated_perimeter=header_perimeter,
        hp_r_out=config.hp_r_out,
        hp_r_in=config.hp_r_in,
        hp_r_vapor=config.hp_r_vapor,
        hp_L_eva=config.hp_l_eva,
        hp_L_con=config.hp_l_con,
        hp_n_eva=config.hp_n_eva,
        hp_n_con=config.hp_n_con,
        hp_n_wick=config.hp_n_wick,
        hp_n_wall=config.hp_n_wall,
        porosity_hp=config.hp_porosity,
        HP_initial_temp=config.initial_hp_temperature,
        fin_thickness=config.fin_thickness,
        fin_height=config.fin_height,
        n_fin_height=config.n_fin_height,
        fin_wrap_ratio=fin_wrap_ratio,
        emissivity=config.emissivity,
        up_view_factor=config.up_view_factor,
        down_view_factor=config.down_view_factor,
        T_space=config.space_temperature,
        hp_wall_mat=mat_wall,
        hp_fluid_mat=mat_hp_fluid,
        hp_wick_mat=mat_wick,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=None,
    )
    ring_hp.audit_hp_multipliers = np.asarray(config.hp_multipliers, dtype=float)
    j_out.k_loss = ring_hp.outlet_k_loss

    volumes = [inlet, outlet] + channel.volumes
    junctions = [j_in, j_out] + channel.internal_junctions
    network = HydraulicNetwork(volumes, junctions, gravity_vector=0.0)
    system = SystemManager(fluid_network=network)
    system.add_component(ring_hp)

    return {
        "system": system,
        "ring_hp": ring_hp,
        "channel": channel,
        "header_solid": header_solid,
        "inlet": inlet,
        "outlet": outlet,
        "j_in": j_in,
        "j_out": j_out,
    }


def check_mapping(case, expected_mask=(True, True, False, True)):
    ring_hp = case["ring_hp"]
    summary = ring_hp.get_hp_status_summary()
    mask = tuple(bool(x) for x in summary["hp_presence_mask"])
    if mask != tuple(expected_mask):
        raise AssertionError(f"hp_presence_mask mismatch: {mask} != {expected_mask}")
    if len(ring_hp.hp_units) != 3:
        raise AssertionError(f"Expected 3 hp_units, got {len(ring_hp.hp_units)}")
    if len(ring_hp.coupler_hps) != 3:
        raise AssertionError(
            f"Expected 3 hp couplers, got {len(ring_hp.coupler_hps)}"
        )

    present_indices = [idx for idx, present in enumerate(mask) if present]
    if present_indices != [0, 1, 3]:
        raise AssertionError(f"Unexpected present node indices: {present_indices}")


def collect_node_coupling_audit(case) -> Dict[str, float]:
    ring_hp = case["ring_hp"]
    channel = case["channel"]

    q_hp_eva_scaled_total = 0.0
    q_fluid_to_hp_scaled_total = 0.0
    q_rej_net_scaled = 0.0
    max_relative_couple = 0.0

    result: Dict[str, float] = {}
    for node_index, hp_pos, multiplier in iter_present_hps(ring_hp):
        hp_unit = ring_hp.hp_units[hp_pos]
        coupler = ring_hp.coupler_hps[hp_pos]

        q_hp_eva_single = float(
            np.sum(hp_unit.hp.boundaries["outer_eva"].current_flux)
        )
        if coupler._last_lambda is None:
            lambda_i = 0.0
        else:
            lambda_i = float(np.asarray(coupler._last_lambda).reshape(-1)[0])
        t_f_current = float(channel.volumes[node_index].T)
        t_f_bc = float(np.asarray(coupler.solid_bc.T_ext).reshape(-1)[0])
        t_wall = float(np.asarray(coupler.solid_bound.T_surface).reshape(-1)[0])
        q_fluid_source_single = lambda_i * (t_wall - t_f_bc)
        q_fluid_to_hp_single = -q_fluid_source_single

        q_hp_eva_scaled = multiplier * q_hp_eva_single
        q_fluid_to_hp_scaled = multiplier * q_fluid_to_hp_single
        residual_couple = q_hp_eva_scaled - q_fluid_to_hp_scaled
        relative_couple = abs(residual_couple) / max(
            abs(q_hp_eva_scaled), abs(q_fluid_to_hp_scaled), 1.0
        )

        breakdown = hp_unit.get_heat_exchange_breakdown()
        q_rej_net_single = float(np.sum(breakdown["net_rejection"]))
        q_rej_net_scaled_i = multiplier * q_rej_net_single

        prefix = f"node{node_index}"
        result[f"{prefix}_multiplier"] = multiplier
        result[f"{prefix}_lambda"] = lambda_i
        result[f"{prefix}_T_fluid_current"] = t_f_current
        result[f"{prefix}_T_fluid_bc"] = t_f_bc
        result[f"{prefix}_T_wall"] = t_wall
        result[f"{prefix}_Q_hp_eva_single"] = q_hp_eva_single
        result[f"{prefix}_Q_fluid_to_hp_single"] = q_fluid_to_hp_single
        result[f"{prefix}_Q_hp_eva_scaled"] = q_hp_eva_scaled
        result[f"{prefix}_Q_fluid_to_hp_scaled"] = q_fluid_to_hp_scaled
        result[f"{prefix}_Q_rej_net_scaled"] = q_rej_net_scaled_i
        result[f"{prefix}_residual_couple"] = residual_couple
        result[f"{prefix}_relative_couple"] = relative_couple

        q_hp_eva_scaled_total += q_hp_eva_scaled
        q_fluid_to_hp_scaled_total += q_fluid_to_hp_scaled
        q_rej_net_scaled += q_rej_net_scaled_i
        max_relative_couple = max(max_relative_couple, relative_couple)

    residual_couple_total = q_hp_eva_scaled_total - q_fluid_to_hp_scaled_total
    relative_couple_total = abs(residual_couple_total) / max(
        abs(q_hp_eva_scaled_total), abs(q_fluid_to_hp_scaled_total), 1.0
    )

    result.update(
        {
            "Q_hp_eva_scaled_total": q_hp_eva_scaled_total,
            "Q_fluid_to_hp_scaled_total": q_fluid_to_hp_scaled_total,
            "Q_rej_net_scaled": q_rej_net_scaled,
            "residual_couple_total": residual_couple_total,
            "relative_couple_total": relative_couple_total,
            "max_relative_couple": max_relative_couple,
        }
    )
    return result


def collect_fluid_enthalpy_audit(case) -> Dict[str, float]:
    inlet = case["inlet"]
    outlet = case["outlet"]
    channel = case["channel"]
    j_in = case["j_in"]
    j_out = case["j_out"]

    w_in = float(j_in.W)
    w_out = float(j_out.W)
    h_in = float(inlet.h if w_in >= 0.0 else channel.volumes[0].h)
    h_out = float(channel.volumes[-1].h if w_out >= 0.0 else outlet.h)
    q_enthalpy_drop = w_in * h_in - w_out * h_out

    return {
        "W_in": w_in,
        "W_out": w_out,
        "T_fluid_in": float(inlet.T),
        "T_fluid_out": float(channel.volumes[-1].T),
        "h_in": h_in,
        "h_out": h_out,
        "Q_fluid_enthalpy_drop": q_enthalpy_drop,
    }


def hp_temperature_stats(ring_hp: RingHP) -> Tuple[float, float, float]:
    arrays = [np.asarray(hp.hp.T, dtype=float).reshape(-1) for hp in ring_hp.hp_units]
    if not arrays:
        return np.nan, np.nan, np.nan
    all_t = np.concatenate(arrays)
    return float(np.min(all_t)), float(np.mean(all_t)), float(np.max(all_t))


def hp_unit_net_rejection(hp_unit) -> float:
    breakdown = hp_unit.get_heat_exchange_breakdown()
    return float(np.sum(breakdown["net_rejection"]))


def run_hp0_failure_monitor(
    config: AuditConfig,
    stop_temp_min: float = 300.0,
    warn_temp_min: float = 750.0,
    stop_temp_max: float = 2000.0,
    stop_delta_t: float = 500.0,
):
    case = build_ringhp_node_coupling_case(config)
    check_mapping(case)
    system = case["system"]
    ring_hp = case["ring_hp"]
    hp0 = ring_hp.hp_units[0]

    print("=" * 86)
    print("TEST5 monitor: RingHP node0 / first heat-pipe temperature failure")
    print("=" * 86)
    print(f"t_end / max_dt           : {config.t_end:.3f} s / {config.max_dt:.5f} s")
    print(f"inner_iter               : {config.inner_iter}")
    print(f"warn_temp_min            : {warn_temp_min:.3f} K")
    print(f"stop_temp_min/max        : {stop_temp_min:.3f} K / {stop_temp_max:.3f} K")
    print(f"stop_delta_t             : {stop_delta_t:.3f} K per step")
    print("-" * 86)

    system.initialize_system()
    previous_hp0_T = np.array(hp0.hp.T, copy=True)
    previous_hp0_U = solid_energy(hp0.hp)
    previous_time = float(system.global_time)
    rows: List[Dict[str, float]] = []
    step_count = 0
    stop_reason = ""
    snapshot_before = np.array(previous_hp0_T, copy=True)
    snapshot_after = np.array(previous_hp0_T, copy=True)
    cpu_start = time.time()

    while system.global_time < config.t_end - 1.0e-12:
        remaining = config.t_end - system.global_time
        dt_target = system.compute_adaptive_dt(
            min_dt=config.min_dt,
            max_dt=config.max_dt,
            safety_factor=config.safety_factor,
        )
        dt = min(dt_target, remaining)
        if dt <= max(1.0e-12, config.min_dt * 1.0e-3):
            stop_reason = "dt_too_small"
            break

        T_before = np.array(hp0.hp.T, copy=True)
        U_before = previous_hp0_U
        system.step(dt=dt, inner_iter=config.inner_iter)
        current_time = float(system.global_time)
        actual_dt = max(current_time - previous_time, 1.0e-30)

        T_after = np.array(hp0.hp.T, copy=True)
        snapshot_before = np.array(T_before, copy=True)
        snapshot_after = np.array(T_after, copy=True)
        U_after = solid_energy(hp0.hp)
        dU_dt = (U_after - U_before) / actual_dt
        delta = T_after - T_before
        abs_delta = np.abs(delta)
        idx_min = np.unravel_index(np.argmin(T_after), T_after.shape)
        idx_max = np.unravel_index(np.argmax(T_after), T_after.shape)
        idx_dmax = np.unravel_index(np.argmax(abs_delta), abs_delta.shape)

        q_eva = float(np.sum(hp0.hp.boundaries["outer_eva"].current_flux))
        q_rej_net = hp_unit_net_rejection(hp0)
        residual = q_eva - q_rej_net - dU_dt
        relative = abs(residual) / max(abs(q_eva), abs(q_rej_net), 1.0)

        row = {
            "step": step_count + 1,
            "t": current_time,
            "dt": actual_dt,
            "hp0_T_min": float(T_after[idx_min]),
            "hp0_T_max": float(T_after[idx_max]),
            "hp0_T_mean": float(np.mean(T_after)),
            "hp0_T_min_index": str(tuple(int(x) for x in idx_min)),
            "hp0_T_max_index": str(tuple(int(x) for x in idx_max)),
            "hp0_max_abs_delta_T": float(abs_delta[idx_dmax]),
            "hp0_max_abs_delta_T_index": str(tuple(int(x) for x in idx_dmax)),
            "hp0_T_before_at_dmax": float(T_before[idx_dmax]),
            "hp0_T_after_at_dmax": float(T_after[idx_dmax]),
            "hp0_Q_eva": q_eva,
            "hp0_Q_rej_net": q_rej_net,
            "hp0_dU_dt": dU_dt,
            "hp0_residual": residual,
            "hp0_relative_residual": relative,
            "node0_T_fluid": float(case["channel"].volumes[0].T),
            "node0_T_wall": float(
                np.asarray(ring_hp.coupler_hps[0].solid_bound.T_surface).reshape(-1)[0]
            ),
            "node0_lambda": float(
                np.asarray(ring_hp.coupler_hps[0]._last_lambda).reshape(-1)[0]
            )
            if ring_hp.coupler_hps[0]._last_lambda is not None
            else 0.0,
            "stop_reason": "",
        }

        reasons = []
        if row["hp0_T_min"] < warn_temp_min:
            reasons.append("below_warn_temp")
        if row["hp0_T_min"] < stop_temp_min:
            reasons.append("below_stop_temp")
        if row["hp0_T_min"] < 0.0:
            reasons.append("negative_temperature")
        if row["hp0_T_max"] > stop_temp_max:
            reasons.append("above_stop_temp")
        if row["hp0_max_abs_delta_T"] > stop_delta_t:
            reasons.append("large_temperature_jump")
        if not np.all(np.isfinite(T_after)):
            reasons.append("nonfinite_temperature")

        row["stop_reason"] = "|".join(reasons)
        rows.append(row)

        print(
            f"step={row['step']:4d} t={row['t']:9.6f}s dt={row['dt']:.6f}s "
            f"Tmin={row['hp0_T_min']:11.3f}K@{row['hp0_T_min_index']} "
            f"Tmax={row['hp0_T_max']:10.3f}K@{row['hp0_T_max_index']} "
            f"dTmax={row['hp0_max_abs_delta_T']:10.3f}K@{row['hp0_max_abs_delta_T_index']} "
            f"Qeva={row['hp0_Q_eva']:10.3f}W "
            f"Qrej={row['hp0_Q_rej_net']:10.3f}W "
            f"dU={row['hp0_dU_dt']:10.3f}W "
            f"res={row['hp0_residual']:10.3e} "
            f"{row['stop_reason']}"
        )

        previous_hp0_T = T_after
        previous_hp0_U = U_after
        previous_time = current_time
        step_count += 1

        if row["stop_reason"] and any(
            reason in row["stop_reason"]
            for reason in (
                "below_stop_temp",
                "negative_temperature",
                "above_stop_temp",
                "large_temperature_jump",
                "nonfinite_temperature",
            )
        ):
            stop_reason = row["stop_reason"]
            break

    elapsed = time.time() - cpu_start
    print("-" * 86)
    print(f"Monitor completed {step_count} steps in {elapsed:.2f} s")
    print(f"Stop reason: {stop_reason or 'completed'}")

    if config.monitor_csv_path:
        write_csv(config.monitor_csv_path, rows)
        print(f"Monitor CSV written to: {config.monitor_csv_path}")
    if config.monitor_snapshot_path and rows:
        save_monitor_snapshot(
            config.monitor_snapshot_path,
            case,
            snapshot_before,
            snapshot_after,
            rows,
        )

    return rows


def add_rolling_metrics(
    row: Dict[str, float],
    rows: List[Dict[str, float]],
    window_steps: int,
):
    window = rows[-window_steps:] if window_steps > 0 else rows
    if not window:
        row["rolling_max_relative_hp_scaled"] = row["relative_hp_scaled"]
        row["rolling_mean_abs_residual_hp_scaled"] = row["abs_residual_hp_scaled"]
        return
    row["rolling_max_relative_hp_scaled"] = float(
        max(r["relative_hp_scaled"] for r in window)
    )
    row["rolling_mean_abs_residual_hp_scaled"] = float(
        np.mean([r["abs_residual_hp_scaled"] for r in window])
    )


def compute_trend(rows: List[Dict[str, float]], skip_steps: int = 5, n_windows: int = 5):
    valid = rows[skip_steps:]
    if len(valid) < n_windows:
        return []
    chunks = np.array_split(np.arange(len(valid)), n_windows)
    trends = []
    for idx, chunk in enumerate(chunks):
        chunk_rows = [valid[int(i)] for i in chunk]
        trends.append(
            {
                "window": idx + 1,
                "t_start": float(chunk_rows[0]["t"]),
                "t_end": float(chunk_rows[-1]["t"]),
                "mean_abs_residual_hp_scaled": float(
                    np.mean([r["abs_residual_hp_scaled"] for r in chunk_rows])
                ),
                "max_relative_hp_scaled": float(
                    max(r["relative_hp_scaled"] for r in chunk_rows)
                ),
                "mean_relative_hp_scaled": float(
                    np.mean([r["relative_hp_scaled"] for r in chunk_rows])
                ),
            }
        )
    return trends


def summarize_criteria(rows: List[Dict[str, float]], trends: List[Dict[str, float]]):
    if not rows:
        return {}
    skipped = rows[5:] if len(rows) > 5 else rows
    final_window_start = max(0, int(0.8 * len(skipped)))
    final_window = skipped[final_window_start:] if skipped else rows
    first_trend = trends[0] if trends else None
    last_trend = trends[-1] if trends else None

    final = rows[-1]
    min_t_hp = min(r["T_hp_min"] for r in rows)
    max_rel_couple = max(r["max_relative_couple"] for r in rows)
    max_rel_couple_total = max(r["relative_couple_total"] for r in rows)
    final_window_mean_relative_hp = float(
        np.mean([r["relative_hp_scaled"] for r in final_window])
    )
    final_dU_ratio = abs(final["dU_hp_scaled_dt"]) / max(
        abs(final["Q_hp_eva_scaled_total"]), 1.0
    )
    trend_abs_decreases = (
        bool(last_trend["mean_abs_residual_hp_scaled"] < first_trend["mean_abs_residual_hp_scaled"])
        if first_trend and last_trend
        else False
    )
    trend_rel_decreases = (
        bool(last_trend["max_relative_hp_scaled"] < first_trend["max_relative_hp_scaled"])
        if first_trend and last_trend
        else False
    )

    return {
        "temperature_valid": min_t_hp >= 750.0,
        "min_T_hp": float(min_t_hp),
        "coupling_valid": max_rel_couple < 1.0e-6
        and max_rel_couple_total < 1.0e-6,
        "max_relative_couple": float(max_rel_couple),
        "max_relative_couple_total": float(max_rel_couple_total),
        "trend_abs_decreases": trend_abs_decreases,
        "trend_rel_decreases": trend_rel_decreases,
        "final_window_mean_relative_hp_scaled": final_window_mean_relative_hp,
        "final_window_valid": final_window_mean_relative_hp < 1.0e-3,
        "final_dU_ratio": float(final_dU_ratio),
        "final_dU_valid": final_dU_ratio < 1.0e-4,
    }


def run_audit(config: AuditConfig):
    case = build_ringhp_node_coupling_case(config)
    system = case["system"]
    ring_hp = case["ring_hp"]

    check_mapping(case)

    print("=" * 86)
    print("TEST5: RingHP node coupling energy conservation")
    print("=" * 86)
    print(f"n_nodes / hp_multipliers : {config.n_nodes} / {list(config.hp_multipliers)}")
    print(f"external_heat_config     : None")
    print(f"t_end / max_dt           : {config.t_end:.3f} s / {config.max_dt:.5f} s")
    print(f"inner_iter               : {config.inner_iter}")
    print("-" * 86)

    system.initialize_system()
    previous_u_hp_scaled = hp_scaled_energy(ring_hp)
    previous_time = float(system.global_time)
    rows: List[Dict[str, float]] = []
    step_count = 0
    cpu_start = time.time()

    while system.global_time < config.t_end - 1.0e-12:
        remaining = config.t_end - system.global_time
        dt_target = system.compute_adaptive_dt(
            min_dt=config.min_dt,
            max_dt=config.max_dt,
            safety_factor=config.safety_factor,
        )
        dt = min(dt_target, remaining)
        if dt <= max(1.0e-12, config.min_dt * 1.0e-3):
            break

        system.step(dt=dt, inner_iter=config.inner_iter)
        current_time = float(system.global_time)
        actual_dt = max(current_time - previous_time, 1.0e-30)

        current_u_hp_scaled = hp_scaled_energy(ring_hp)
        dU_hp_scaled_dt = (current_u_hp_scaled - previous_u_hp_scaled) / actual_dt

        row = {
            "step": step_count + 1,
            "t": current_time,
            "dt": actual_dt,
        }
        row.update(collect_fluid_enthalpy_audit(case))
        row.update(collect_node_coupling_audit(case))
        row["Q_fluid_enthalpy_drop_minus_Q_hp_eva"] = (
            row["Q_fluid_enthalpy_drop"] - row["Q_hp_eva_scaled_total"]
        )
        row["Q_fluid_enthalpy_drop_minus_Q_rej_net"] = (
            row["Q_fluid_enthalpy_drop"] - row["Q_rej_net_scaled"]
        )
        row["U_hp_scaled"] = current_u_hp_scaled
        row["dU_hp_scaled_dt"] = dU_hp_scaled_dt
        row["residual_hp_scaled"] = (
            row["Q_hp_eva_scaled_total"]
            - row["Q_rej_net_scaled"]
            - row["dU_hp_scaled_dt"]
        )
        row["abs_residual_hp_scaled"] = abs(row["residual_hp_scaled"])
        row["relative_hp_scaled"] = row["abs_residual_hp_scaled"] / max(
            abs(row["Q_hp_eva_scaled_total"]), abs(row["Q_rej_net_scaled"]), 1.0
        )
        t_min, t_mean, t_max = hp_temperature_stats(ring_hp)
        row["T_hp_min"] = t_min
        row["T_hp_mean"] = t_mean
        row["T_hp_max"] = t_max

        rows.append(row)
        add_rolling_metrics(row, rows, config.rolling_window_steps)

        previous_u_hp_scaled = current_u_hp_scaled
        previous_time = current_time
        step_count += 1

        if step_count == 1 or step_count % max(config.print_every, 1) == 0:
            print(
                f"t={row['t']:9.3f}s "
                f"Qdh={row['Q_fluid_enthalpy_drop']:12.4f} W "
                f"Qeva={row['Q_hp_eva_scaled_total']:12.4f} W "
                f"Qrej={row['Q_rej_net_scaled']:12.4f} W "
                f"dU={row['dU_hp_scaled_dt']:12.4f} W "
                f"res_hp={row['residual_hp_scaled']:12.4e} W "
                f"rel_hp={row['relative_hp_scaled']:.3e} "
                f"roll_abs={row['rolling_mean_abs_residual_hp_scaled']:.4e} W "
                f"Tmin={row['T_hp_min']:.3f} K"
            )

    elapsed = time.time() - cpu_start
    trends = compute_trend(rows)
    criteria = summarize_criteria(rows, trends)
    final = rows[-1] if rows else None

    print("-" * 86)
    print(f"Completed {step_count} steps in {elapsed:.2f} s")
    if final:
        print("Final audit:")
        print(f"  Q_fluid_enthalpy_drop      = {final['Q_fluid_enthalpy_drop']: .6e} W")
        print(f"  Q_hp_eva_scaled_total       = {final['Q_hp_eva_scaled_total']: .6e} W")
        print(f"  Q_dh - Q_hp_eva             = {final['Q_fluid_enthalpy_drop_minus_Q_hp_eva']: .6e} W")
        print(f"  Q_fluid_to_hp_scaled_total  = {final['Q_fluid_to_hp_scaled_total']: .6e} W")
        print(f"  Q_rej_net_scaled            = {final['Q_rej_net_scaled']: .6e} W")
        print(f"  Q_dh - Q_rej_net            = {final['Q_fluid_enthalpy_drop_minus_Q_rej_net']: .6e} W")
        print(f"  dU_hp_scaled_dt             = {final['dU_hp_scaled_dt']: .6e} W")
        print(f"  residual_couple_total       = {final['residual_couple_total']: .6e} W")
        print(f"  relative_couple_total       = {final['relative_couple_total']: .6e}")
        print(f"  max_relative_couple         = {final['max_relative_couple']: .6e}")
        print(f"  residual_hp_scaled          = {final['residual_hp_scaled']: .6e} W")
        print(f"  relative_hp_scaled          = {final['relative_hp_scaled']: .6e}")
        print(f"  T_hp_min                    = {final['T_hp_min']: .6f} K")

    if trends:
        print("Trend windows after skipping first 5 steps:")
        for item in trends:
            print(
                f"  W{item['window']}: "
                f"t={item['t_start']:.3f}-{item['t_end']:.3f}s "
                f"mean_abs={item['mean_abs_residual_hp_scaled']:.6e} W "
                f"max_rel={item['max_relative_hp_scaled']:.6e}"
            )

    if criteria:
        print("Criteria summary:")
        print(f"  temperature_valid          = {criteria['temperature_valid']}")
        print(f"  min_T_hp                   = {criteria['min_T_hp']:.6f} K")
        print(f"  coupling_valid             = {criteria['coupling_valid']}")
        print(f"  max_relative_couple        = {criteria['max_relative_couple']:.6e}")
        print(f"  max_relative_couple_total  = {criteria['max_relative_couple_total']:.6e}")
        print(f"  trend_abs_decreases        = {criteria['trend_abs_decreases']}")
        print(f"  trend_rel_decreases        = {criteria['trend_rel_decreases']}")
        print(
            "  final_window_mean_relative = "
            f"{criteria['final_window_mean_relative_hp_scaled']:.6e}"
        )
        print(f"  final_window_valid         = {criteria['final_window_valid']}")
        print(f"  final_dU_ratio             = {criteria['final_dU_ratio']:.6e}")
        print(f"  final_dU_valid             = {criteria['final_dU_valid']}")

    if config.csv_path:
        write_csv(config.csv_path, rows)
        print(f"CSV written to: {config.csv_path}")

    if config.restart_path:
        save_restart(config.restart_path, config, case, rows, trends, criteria)

    return rows, trends


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run TEST5 RingHP node coupling energy conservation audit."
    )
    parser.add_argument("--t-end", type=float, default=None)
    parser.add_argument("--min-dt", type=float, default=None)
    parser.add_argument("--max-dt", type=float, default=None)
    parser.add_argument("--safety-factor", type=float, default=None)
    parser.add_argument("--inner-iter", type=int, default=None)
    parser.add_argument("--print-every", type=int, default=None)
    parser.add_argument("--csv-path", default=None)
    parser.add_argument("--no-csv", action="store_true")
    parser.add_argument("--restart-path", default=None)
    parser.add_argument("--no-restart", action="store_true")
    parser.add_argument("--monitor-hp0", action="store_true")
    parser.add_argument("--monitor-csv-path", default=None)
    parser.add_argument("--monitor-snapshot-path", default=None)
    parser.add_argument("--monitor-stop-temp-min", type=float, default=300.0)
    parser.add_argument("--monitor-warn-temp-min", type=float, default=750.0)
    parser.add_argument("--monitor-stop-temp-max", type=float, default=2000.0)
    parser.add_argument("--monitor-stop-delta-t", type=float, default=500.0)
    return parser.parse_args()


def config_from_args(args) -> AuditConfig:
    replacements = {}
    for arg_name, field_name in (
        ("t_end", "t_end"),
        ("min_dt", "min_dt"),
        ("max_dt", "max_dt"),
        ("safety_factor", "safety_factor"),
        ("inner_iter", "inner_iter"),
        ("print_every", "print_every"),
    ):
        value = getattr(args, arg_name)
        if value is not None:
            replacements[field_name] = value

    if args.csv_path is not None:
        replacements["csv_path"] = args.csv_path
    if args.no_csv:
        replacements["csv_path"] = None
    if args.restart_path is not None:
        replacements["restart_path"] = args.restart_path
    if args.no_restart:
        replacements["restart_path"] = None
    if args.monitor_csv_path is not None:
        replacements["monitor_csv_path"] = args.monitor_csv_path
    if args.monitor_snapshot_path is not None:
        replacements["monitor_snapshot_path"] = args.monitor_snapshot_path

    return replace(DEFAULT_CONFIG, **replacements) if replacements else DEFAULT_CONFIG


if __name__ == "__main__":
    parsed_args = parse_args()
    parsed_config = config_from_args(parsed_args)
    if parsed_args.monitor_hp0:
        run_hp0_failure_monitor(
            parsed_config,
            stop_temp_min=parsed_args.monitor_stop_temp_min,
            warn_temp_min=parsed_args.monitor_warn_temp_min,
            stop_temp_max=parsed_args.monitor_stop_temp_max,
            stop_delta_t=parsed_args.monitor_stop_delta_t,
        )
    else:
        run_audit(parsed_config)
