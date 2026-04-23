import argparse
import csv
import os
import sys
import time
from dataclasses import dataclass, replace
from typing import Dict, List, Optional

import numpy as np

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
    # Operating point. The default mass flow is the v4_2 single-branch flow
    # divided by 78 radiator tubes.
    inlet_temperature: float = 968.0
    inlet_mass_flow: float = 0.1815 / 78.0
    outlet_pressure: float = 161000.0
    initial_fluid_temperature: float = 863.0
    initial_hp_temperature: float = 800.0

    # Single-header-cell geometry.
    header_length: float = 0.06
    header_flow_area: float = 0.0016065
    header_dh: float = 0.04167
    header_wall_thickness: float = 0.002

    # Heat pipe geometry from v4_2.
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

    # Fin/radiation model from v4_2.
    fin_thickness: float = 0.0003
    fin_height: float = 22.65e-3
    n_fin_height: int = 15
    emissivity: float = 0.93
    up_view_factor: float = 0.0
    down_view_factor: float = 0.3
    space_temperature: float = 3.0

    # Time integration controls.
    t_end: float = 50.0
    min_dt: float = 1.0e-3
    max_dt: float = 0.2
    safety_factor: float = 0.5
    inner_iter: int = 2
    print_every: int = 20

    # Output controls.
    csv_path: Optional[str] = None


DEFAULT_CONFIG = AuditConfig()


def lyon_martinelli(Re, Pr, p_d_ratio=1.0):
    """Liquid-metal forced-convection correlation used by v4_2."""
    _ = p_d_ratio
    pe = np.maximum(np.asarray(Re, dtype=float) * np.asarray(Pr, dtype=float), 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


def _solid_energy(solid) -> float:
    if hasattr(solid, "_update_properties"):
        solid._update_properties()
    return float(np.sum(solid.thermal_capacitance * solid.T))


def _fluid_energy(volumes) -> float:
    total = 0.0
    for vol in volumes:
        total += float(vol.rho * vol.vol * vol.h)
    return total


def _hp_breakdown(hp_unit) -> Dict[str, float]:
    breakdown = hp_unit.get_heat_exchange_breakdown()
    return {key: float(np.sum(value)) for key, value in breakdown.items()}


def build_single_hp_header_case(config: AuditConfig):
    mat_fluid = Sodium()
    mat_wall = SS316(name="Audit_SS316_Wall")
    mat_hp_fluid = SodiumHP(name="Audit_HP_Fluid_Na")
    mat_wick = WickMaterial(
        name="Audit_HP_Wick_Composite",
        solid_mat=SS316(),
        fluid_mat=mat_hp_fluid,
        porosity=config.hp_porosity,
        r_vapor=config.hp_r_vapor,
        r_in_wall=config.hp_r_in,
    )

    inlet = IncompressibleBoundaryVolume(
        name="Audit_Inlet",
        material=mat_fluid,
        P=config.outlet_pressure + 1000.0,
        T=config.inlet_temperature,
        flow_area=config.header_flow_area,
        hydraulic_diam=config.header_dh,
    )
    outlet = IncompressibleBoundaryVolume(
        name="Audit_Outlet",
        material=mat_fluid,
        P=config.outlet_pressure,
        T=config.initial_fluid_temperature,
        flow_area=config.header_flow_area,
        hydraulic_diam=config.header_dh,
    )
    outlet.is_pressure_boundary = True

    channel = IncompressibleFluidChannel(
        name="Audit_Header_Channel",
        n_nodes=1,
        total_length=config.header_length,
        flow_area=config.header_flow_area,
        hydraulic_diam=config.header_dh,
        initial_P=config.outlet_pressure + 500.0,
        initial_T=config.initial_fluid_temperature,
        material=mat_fluid,
    )

    j_in = InletJunction(
        name="Audit_J_Inlet",
        from_vol=inlet,
        to_vol=channel.volumes[0],
        W_initial=config.inlet_mass_flow,
    )
    j_out = FlowJunction(
        name="Audit_J_Outlet",
        from_vol=channel.volumes[0],
        to_vol=outlet,
        flow_area=config.header_flow_area,
        k_loss=0.0,
    )

    r_in_header = config.header_dh / 2.0
    header_perimeter = 2.0 * np.pi * r_in_header
    header_mesh = Mesh2D(
        x_dim=config.header_wall_thickness,
        n_x=1,
        y_dim=config.header_length,
        n_y=1,
        geometry_type="cylindrical",
        inner_radius=r_in_header,
    )
    header_solid = HeatConduction2D(
        mesh=header_mesh,
        material=mat_wall,
        initial_temp=config.initial_fluid_temperature,
        name="Audit_Header_Wall",
    )
    header_solid.boundaries["right"].add_flux_condition(q_flux=0.0)
    header_solid.boundaries["top"].add_flux_condition(q_flux=0.0)
    header_solid.boundaries["bottom"].add_flux_condition(q_flux=0.0)

    fin_wrap_ratio = (2.0 * config.fin_thickness) / (2.0 * np.pi * config.hp_r_out)

    ring_hp = RingHP(
        name="Audit_Single_RingHP",
        fluid_channel=channel,
        solid_header=header_solid,
        hp_multipliers=[1.0],
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
    )
    j_out.k_loss = ring_hp.outlet_k_loss

    volumes = [inlet, outlet] + channel.volumes
    junctions = [j_in, j_out]
    network = HydraulicNetwork(volumes, junctions, gravity_vector=0.0)
    system = SystemManager(fluid_network=network)
    system.add_component(ring_hp)

    return {
        "system": system,
        "ring_hp": ring_hp,
        "hp_unit": ring_hp.hp_units[0],
        "channel": channel,
        "header_solid": header_solid,
        "inlet": inlet,
        "outlet": outlet,
        "j_in": j_in,
        "j_out": j_out,
        "finite_fluid_volumes": channel.volumes,
        "fluid_material": mat_fluid,
    }


def collect_energy_state(case) -> Dict[str, float]:
    fluid_u = _fluid_energy(case["finite_fluid_volumes"])
    header_u = _solid_energy(case["header_solid"])
    hp_u = _solid_energy(case["hp_unit"].hp)
    return {
        "fluid": fluid_u,
        "header": header_u,
        "hp": hp_u,
        "total": fluid_u + header_u + hp_u,
    }


def audit_step(case, before_energy: Dict[str, float], after_energy: Dict[str, float], dt: float):
    hp_unit = case["hp_unit"]
    channel_vol = case["channel"].volumes[0]
    inlet = case["inlet"]
    j_in = case["j_in"]
    j_out = case["j_out"]

    w_in = float(j_in.W)
    w_out = float(j_out.W)
    h_in = float(inlet.h)
    h_out = float(channel_vol.h)

    q_fluid_loss = w_in * h_in - w_out * h_out

    eva_flux_raw = float(np.sum(hp_unit.hp.boundaries["outer_eva"].current_flux))
    con_flux_raw = float(np.sum(hp_unit.hp.boundaries["outer_con"].current_flux))
    q_hp_eva = eva_flux_raw
    q_hp_rej_solver = -con_flux_raw

    breakdown = _hp_breakdown(hp_unit)
    dU_fluid_dt = (after_energy["fluid"] - before_energy["fluid"]) / dt
    dU_header_dt = (after_energy["header"] - before_energy["header"]) / dt
    dU_hp_dt = (after_energy["hp"] - before_energy["hp"]) / dt

    residual_local = (
        q_fluid_loss
        - q_hp_rej_solver
        - dU_fluid_dt
        - dU_header_dt
        - dU_hp_dt
    )
    residual_interface = q_fluid_loss - q_hp_eva - dU_fluid_dt - dU_header_dt
    residual_hp = q_hp_eva - q_hp_rej_solver - dU_hp_dt

    rel_residual = abs(residual_local) / max(abs(q_fluid_loss), 1.0)

    return {
        "t": float(case["system"].global_time),
        "dt": float(dt),
        "W_in": w_in,
        "W_out": w_out,
        "T_in": float(inlet.T),
        "T_fluid": float(channel_vol.T),
        "h_in": h_in,
        "h_out": h_out,
        "Q_fluid_loss": q_fluid_loss,
        "Q_hp_eva": q_hp_eva,
        "Q_hp_eva_boundary_raw": eva_flux_raw,
        "Q_hp_rej_solver": q_hp_rej_solver,
        "Q_hp_con_boundary_raw": con_flux_raw,
        "Q_hp_bare_rad": breakdown["bare_radiation"],
        "Q_hp_fin_rad": breakdown["fin_radiation"],
        "Q_hp_fin_net_from_root": breakdown["fin_net_from_root"],
        "Q_hp_gross_rejection": breakdown["gross_rejection"],
        "Q_hp_net_rejection": breakdown["net_rejection"],
        "dU_fluid_dt": dU_fluid_dt,
        "dU_header_dt": dU_header_dt,
        "dU_hp_dt": dU_hp_dt,
        "dU_total_dt": (
            after_energy["total"] - before_energy["total"]
        ) / dt,
        "U_fluid": after_energy["fluid"],
        "U_header": after_energy["header"],
        "U_hp": after_energy["hp"],
        "residual_local": residual_local,
        "residual_interface": residual_interface,
        "residual_hp": residual_hp,
        "relative_residual_local": rel_residual,
    }


def _write_csv(path: str, rows: List[Dict[str, float]]):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_audit(config: AuditConfig):
    case = build_single_hp_header_case(config)
    system = case["system"]

    print("=" * 78)
    print("Single HP Header Energy Audit")
    print("=" * 78)
    print(f"Inlet T             : {config.inlet_temperature:.3f} K")
    print(f"Inlet W             : {config.inlet_mass_flow:.8f} kg/s")
    print(f"Header length       : {config.header_length:.5f} m")
    print(f"HP multiplier       : 1.0")
    print(f"Emissivity          : {config.emissivity:.5f}")
    print(f"t_end / max_dt      : {config.t_end:.3f} s / {config.max_dt:.5f} s")
    print("-" * 78)

    system.initialize_system()
    previous_energy = collect_energy_state(case)
    previous_time = float(system.global_time)
    rows: List[Dict[str, float]] = []
    step_count = 0
    cpu_start = time.time()

    while system.global_time < config.t_end - 1.0e-12:
        remaining = config.t_end - system.global_time
        tiny_step = max(1.0e-12, config.min_dt * 1.0e-3)
        if remaining <= tiny_step:
            break

        dt_target = system.compute_adaptive_dt(
            min_dt=config.min_dt,
            max_dt=config.max_dt,
            safety_factor=config.safety_factor,
        )
        dt = min(dt_target, remaining)
        if dt <= tiny_step:
            break

        before = previous_energy
        system.step(dt=dt, inner_iter=config.inner_iter)
        after = collect_energy_state(case)

        actual_dt = max(float(system.global_time) - previous_time, 1.0e-30)
        row = audit_step(case, before, after, actual_dt)
        rows.append(row)

        previous_energy = after
        previous_time = float(system.global_time)
        step_count += 1

        if step_count == 1 or step_count % max(config.print_every, 1) == 0:
            print(
                f"t={row['t']:8.4f}s dt={row['dt']:7.4f}s "
                f"W={row['W_in']:.7f}/{row['W_out']:.7f} kg/s "
                f"T_f={row['T_fluid']:8.3f}K "
                f"Qloop={row['Q_fluid_loss']:10.3f} W "
                f"Qeva={row['Q_hp_eva']:10.3f} W "
                f"Qrej={row['Q_hp_rej_solver']:10.3f} W "
                f"dU={row['dU_total_dt']:10.3f} W "
                f"res={row['residual_local']:10.3f} W "
                f"rel={row['relative_residual_local']:.3e}"
            )

    elapsed = time.time() - cpu_start
    final = rows[-1] if rows else None
    print("-" * 78)
    print(f"Completed {step_count} steps in {elapsed:.2f} s")
    if final is not None:
        print("Final audit:")
        print(f"  Q_fluid_loss          = {final['Q_fluid_loss']: .6e} W")
        print(f"  Q_hp_eva              = {final['Q_hp_eva']: .6e} W")
        print(f"  Q_hp_rej_solver       = {final['Q_hp_rej_solver']: .6e} W")
        print(f"  Q_hp_bare_rad         = {final['Q_hp_bare_rad']: .6e} W")
        print(f"  Q_hp_fin_rad          = {final['Q_hp_fin_rad']: .6e} W")
        print(f"  Q_hp_fin_net_from_root= {final['Q_hp_fin_net_from_root']: .6e} W")
        print(f"  dU_fluid/dt           = {final['dU_fluid_dt']: .6e} W")
        print(f"  dU_header/dt          = {final['dU_header_dt']: .6e} W")
        print(f"  dU_hp/dt              = {final['dU_hp_dt']: .6e} W")
        print(f"  residual_local        = {final['residual_local']: .6e} W")
        print(f"  residual_interface    = {final['residual_interface']: .6e} W")
        print(f"  residual_hp           = {final['residual_hp']: .6e} W")
        print(f"  relative residual     = {final['relative_residual_local']: .6e}")

    if config.csv_path:
        _write_csv(config.csv_path, rows)
        print(f"CSV written to: {config.csv_path}")

    return rows


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a single heat-pipe/header-cell energy audit."
    )
    parser.add_argument(
        "--case",
        choices=("normal", "weak-radiation", "steady"),
        default="normal",
        help="Predefined scenario.",
    )
    parser.add_argument("--t-end", type=float, default=None)
    parser.add_argument("--min-dt", type=float, default=None)
    parser.add_argument("--max-dt", type=float, default=None)
    parser.add_argument("--safety-factor", type=float, default=None)
    parser.add_argument("--inner-iter", type=int, default=None)
    parser.add_argument("--print-every", type=int, default=None)
    parser.add_argument("--inlet-temperature", type=float, default=None)
    parser.add_argument(
        "--initial-temperature",
        type=float,
        default=None,
        help="Set both initial fluid/header and heat-pipe temperatures.",
    )
    parser.add_argument("--initial-fluid-temperature", type=float, default=None)
    parser.add_argument("--initial-hp-temperature", type=float, default=None)
    parser.add_argument("--mass-flow", type=float, default=None)
    parser.add_argument("--emissivity", type=float, default=None)
    parser.add_argument("--space-temperature", type=float, default=None)
    parser.add_argument("--csv-path", default=None)
    parser.add_argument("--no-csv", action="store_true")
    return parser.parse_args()


def config_from_args(args) -> AuditConfig:
    config = DEFAULT_CONFIG

    if args.case == "weak-radiation":
        config = replace(
            config,
            emissivity=0.01,
            t_end=10.0,
            max_dt=0.1,
            csv_path=os.path.join(current_dir, "single_hp_header_audit_weak.csv"),
        )
    elif args.case == "steady":
        config = replace(
            config,
            t_end=200.0,
            max_dt=0.1,
            csv_path=os.path.join(current_dir, "single_hp_header_audit_steady.csv"),
        )
    else:
        config = replace(
            config,
            csv_path=os.path.join(current_dir, "single_hp_header_audit_normal.csv"),
        )

    replacements = {}
    for arg_name, field_name in (
        ("t_end", "t_end"),
        ("min_dt", "min_dt"),
        ("max_dt", "max_dt"),
        ("safety_factor", "safety_factor"),
        ("inner_iter", "inner_iter"),
        ("print_every", "print_every"),
        ("inlet_temperature", "inlet_temperature"),
        ("mass_flow", "inlet_mass_flow"),
        ("emissivity", "emissivity"),
        ("space_temperature", "space_temperature"),
    ):
        value = getattr(args, arg_name)
        if value is not None:
            replacements[field_name] = value

    if args.csv_path is not None:
        replacements["csv_path"] = args.csv_path
    if args.no_csv:
        replacements["csv_path"] = None

    if args.initial_temperature is not None:
        replacements["initial_fluid_temperature"] = args.initial_temperature
        replacements["initial_hp_temperature"] = args.initial_temperature
    if args.initial_fluid_temperature is not None:
        replacements["initial_fluid_temperature"] = args.initial_fluid_temperature
    if args.initial_hp_temperature is not None:
        replacements["initial_hp_temperature"] = args.initial_hp_temperature

    if replacements:
        config = replace(config, **replacements)
    return config


if __name__ == "__main__":
    run_audit(config_from_args(parse_args()))
