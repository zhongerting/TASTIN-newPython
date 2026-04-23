"""
Single header-cell + one heat-pipe verification case.

This is a local, deliberately small case derived from coolant-loop v4_2:

    fixed-T/W inlet -> one header control volume -> pressure outlet
                               |
                               +-> one HPwithFin unit -> radiation sink

The goal is detailed observation of one inserted heat pipe without the full
loop topology.  It records fluid/header/HP temperatures, heat rates, and a
finite-domain energy audit to CSV and PNG.
"""

import csv
from datetime import datetime
import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
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
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
)
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager


logging.basicConfig(level=logging.WARNING)


def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    pe = np.maximum(Re * Pr, 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


def solid_energy(solid):
    if hasattr(solid, "_update_properties"):
        solid._update_properties()
    return float(np.sum(solid.thermal_capacitance * solid.T))


def fluid_energy(vol):
    return float(vol.rho * vol.vol * vol.h)


def hp_temperatures(hp_unit):
    field = hp_unit.hp.T.reshape(hp_unit.hp.shape_nodes)
    eva = field[:, : hp_unit.hp.n_eva]
    con = field[:, hp_unit.hp.n_eva + hp_unit.hp.n_aba :]
    t_con_surf, _ = hp_unit.hp.boundaries["outer_con"].get_coupling_surface_snapshot()
    t_eva_surf, _ = hp_unit.hp.boundaries["outer_eva"].get_coupling_surface_snapshot()
    min_idx = np.unravel_index(np.argmin(field), field.shape)
    return {
        "T_hp_min": float(np.min(field)),
        "T_hp_min_radial_idx": int(min_idx[0]),
        "T_hp_min_axial_idx": int(min_idx[1]),
        "T_hp_eva_avg": float(np.mean(eva)),
        "T_hp_con_avg": float(np.mean(con)),
        "T_hp_eva_surface": float(np.mean(t_eva_surf)),
        "T_hp_con_surface": float(np.mean(t_con_surf)),
        "T_fin_tip_avg": float(np.mean(hp_unit.last_fin_temperature[:, -1])),
    }


def hp_heat_rates(hp_unit):
    breakdown = hp_unit.get_heat_exchange_breakdown()
    q_eva_in = float(np.sum(hp_unit.hp.boundaries["outer_eva"].current_flux))
    q_aba_solver = float(-np.sum(hp_unit.hp.boundaries["outer_aba"].current_flux))
    q_con_solver = float(-np.sum(hp_unit.hp.boundaries["outer_con"].current_flux))
    return {
        "Q_hp_eva_in": q_eva_in,
        "Q_hp_aba_solver_rej": q_aba_solver,
        "Q_hp_con_solver_rej": q_con_solver,
        "Q_hp_solver_rej": q_aba_solver + q_con_solver,
        "Q_hp_bare_rad": float(np.sum(breakdown["bare_radiation"])),
        "Q_hp_fin_rad": float(np.sum(breakdown["fin_radiation"])),
        "Q_hp_fin_net_from_root": float(np.sum(breakdown["fin_net_from_root"])),
        "Q_hp_gross_rej": float(np.sum(breakdown["gross_rejection"])),
    }


def available_output_path(path):
    """Use the requested path, or add a timestamp if the file is locked/open."""
    try:
        with open(path, "a", encoding="utf-8"):
            pass
        return path
    except PermissionError:
        stem, ext = os.path.splitext(path)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"{stem}_{stamp}{ext}"


def main():
    print("=" * 72)
    print("Single Header Cell + One Heat Pipe Test")
    print("=" * 72)

    # Boundary and initial conditions.
    t_inlet = float(os.environ.get("SINGLE_HP_T_INLET", "968.0"))
    p_outlet = float(os.environ.get("SINGLE_HP_P_OUTLET", "161000.0"))
    t_init = float(os.environ.get("SINGLE_HP_T_INIT", "863.0"))
    hp_init_temp = float(os.environ.get("SINGLE_HP_HP_INIT", "800.0"))
    t_end = float(os.environ.get("SINGLE_HP_T_END", "500.0"))
    max_dt = float(os.environ.get("SINGLE_HP_MAX_DT", "0.2"))

    # v4_2 geometry, reduced to one sixth of one ring and one inserted HP.
    l_ring = 0.793
    n_ring_v42 = 6
    l_cell = l_ring / n_ring_v42
    n_hp_per_ring = 78.0
    w_total_design = 0.5445
    w_ring = w_total_design / 2.0
    w_single_hp = float(os.environ.get("SINGLE_HP_W", str(w_ring / n_hp_per_ring)))

    a_ring = 0.0016065
    dh_ring = 0.04167
    wall_thickness_ring = 0.002
    r_in_ring = dh_ring / 2.0
    perim_ring = 2.0 * np.pi * r_in_ring

    hp_r_out, hp_r_in, hp_r_vapor = 0.0085, 0.0081, 0.0075
    hp_l_eva, hp_l_con = 0.06, 0.482
    hp_n_eva, hp_n_con = 1, 12
    hp_n_wick, hp_n_wall = 1, 1
    hp_porosity = 0.5

    fin_thickness = 0.0003
    fin_height = 22.65e-3
    n_fin_height = 15
    fin_wrap_ratio = (2.0 * fin_thickness) / (2.0 * np.pi * hp_r_out)
    # Local diagnostic default: weaken radiation so the HP cold end stays away
    # from the Na melting/mushy range during detailed bookkeeping checks.
    emissivity = float(os.environ.get("SINGLE_HP_EMISSIVITY", "0.20"))
    up_vf = float(os.environ.get("SINGLE_HP_UP_VF", "1.0"))
    down_vf = float(os.environ.get("SINGLE_HP_DOWN_VF", "0.675"))
    t_space = 3.0

    nak = Sodium()
    hp_fluid = SodiumHP(name="HP_Fluid_Na")
    mat_wall = SS316(name="SS316_Wall")
    mat_wick = WickMaterial(
        name="HP_Wick_Composite",
        solid_mat=SS316(),
        fluid_mat=hp_fluid,
        porosity=hp_porosity,
        r_vapor=hp_r_vapor,
        r_in_wall=hp_r_in,
    )

    inlet = IncompressibleBoundaryVolume(
        name="SingleCell_Inlet",
        material=nak,
        P=p_outlet + 5000.0,
        T=t_inlet,
    )
    header_channel = IncompressibleFluidChannel(
        name="SingleHeaderCell_Channel",
        n_nodes=1,
        total_length=l_cell,
        flow_area=a_ring,
        hydraulic_diam=dh_ring,
        initial_P=p_outlet,
        initial_T=t_init,
        material=nak,
    )
    outlet = IncompressibleBoundaryVolume(
        name="SingleCell_Outlet",
        material=nak,
        P=p_outlet,
        T=t_init,
    )
    outlet.is_pressure_boundary = True

    mesh = Mesh2D(
        x_dim=wall_thickness_ring,
        n_x=1,
        y_dim=l_cell,
        n_y=1,
        geometry_type="cylindrical",
        inner_radius=r_in_ring,
    )
    header_wall = HeatConduction2D(mesh=mesh, material=mat_wall, initial_temp=t_init)
    header_wall.name = "SingleHeaderCell_Wall"
    header_wall.boundaries["right"].add_flux_condition(q_flux=0.0)
    header_wall.boundaries["top"].add_flux_condition(q_flux=0.0)
    header_wall.boundaries["bottom"].add_flux_condition(q_flux=0.0)

    ring_hp = RingHP(
        name="SingleHeaderCell_OneHP",
        fluid_channel=header_channel,
        solid_header=header_wall,
        hp_multipliers=[1],
        header_flow_area=a_ring,
        header_dh=dh_ring,
        header_heated_perimeter=perim_ring,
        hp_r_out=hp_r_out,
        hp_r_in=hp_r_in,
        hp_r_vapor=hp_r_vapor,
        hp_L_eva=hp_l_eva,
        hp_L_con=hp_l_con,
        hp_n_eva=hp_n_eva,
        hp_n_con=hp_n_con,
        hp_n_wick=hp_n_wick,
        hp_n_wall=hp_n_wall,
        porosity_hp=hp_porosity,
        HP_initial_temp=hp_init_temp,
        fin_thickness=fin_thickness,
        fin_height=fin_height,
        n_fin_height=n_fin_height,
        fin_wrap_ratio=fin_wrap_ratio,
        emissivity=emissivity,
        up_view_factor=up_vf,
        down_view_factor=down_vf,
        T_space=t_space,
        hp_wall_mat=mat_wall,
        hp_fluid_mat=hp_fluid,
        hp_wick_mat=mat_wick,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
    )
    hp_unit = ring_hp.hp_units[0]

    j_inlet = InletJunction(
        name="J_SingleCell_Inlet",
        from_vol=inlet,
        to_vol=header_channel.volumes[0],
        W_initial=w_single_hp,
    )
    j_outlet = FlowJunction(
        name="J_SingleCell_Outlet",
        from_vol=header_channel.volumes[0],
        to_vol=outlet,
        flow_area=a_ring,
        k_loss=0.5 + ring_hp.outlet_k_loss,
    )

    volumes = [inlet] + header_channel.volumes + [outlet]
    junctions = [j_inlet, j_outlet]
    net = HydraulicNetwork(volumes=volumes, junctions=junctions, gravity_vector=0.0)
    sys_mgr = SystemManager(fluid_network=net)
    sys_mgr.add_component(ring_hp)
    sys_mgr.initialize_system()

    print(f"Cell length:      {l_cell:.6f} m")
    print(f"Header area:      {a_ring:.7f} m2")
    print(f"One-HP flow:      {w_single_hp:.8f} kg/s")
    print(f"Inlet T:          {t_inlet:.2f} K")
    print(f"Radiation eps/F:  {emissivity:.3f} / {up_vf:.3f}, {down_vf:.3f}")
    print(f"Initial fluid/wall/HP T: {t_init:.2f} / {t_init:.2f} / {hp_init_temp:.2f} K")
    print(f"End time:         {t_end:g} s")

    history = []
    t = 0.0
    prev_time = 0.0
    prev_u = None
    step_count = 0

    while t < t_end:
        dt = sys_mgr.compute_adaptive_dt(min_dt=1.0e-4, max_dt=max_dt, safety_factor=1.0)
        sys_mgr.step(dt=dt, inner_iter=2)
        t = sys_mgr.global_time
        step_count += 1

        vol = header_channel.volumes[0]
        h_in_set = nak.enthalpy_saturated_liquid(t_inlet)
        h_out_set = nak.enthalpy_saturated_liquid(vol.T)
        q_loop_flux_set = j_inlet.W * h_in_set - j_outlet.W * h_out_set
        # Use the actual solver state for strict bookkeeping.  Boundary pressure
        # updates can shift inlet.h/T slightly from the nominal input values.
        q_loop_flux = j_inlet.W * inlet.h - j_outlet.W * vol.h

        u_fluid = fluid_energy(vol)
        u_header = solid_energy(header_wall)
        u_hp = solid_energy(hp_unit.hp)
        u_total = u_fluid + u_header + u_hp
        if prev_u is None:
            dudt = 0.0
        else:
            dudt = (u_total - prev_u) / max(t - prev_time, 1.0e-30)
        prev_u = u_total
        prev_time = t

        q_header_in = float(np.sum(header_wall.boundaries["left"].current_flux))
        hp_t = hp_temperatures(hp_unit)
        hp_q = hp_heat_rates(hp_unit)
        residual = q_loop_flux - hp_q["Q_hp_solver_rej"] - dudt

        row = {
            "t": t,
            "dt": dt,
            "W_in": j_inlet.W,
            "W_out": j_outlet.W,
            "T_inlet_actual": inlet.T,
            "h_inlet_actual": inlet.h,
            "h_outlet_actual": vol.h,
            "T_fluid": vol.T,
            "P_fluid": vol.P,
            "T_header_wall": float(np.mean(header_wall.T)),
            "Q_loop_flux": q_loop_flux,
            "Q_loop_flux_setpoint": q_loop_flux_set,
            "Q_header_in": q_header_in,
            "U_fluid": u_fluid,
            "U_header": u_header,
            "U_hp": u_hp,
            "U_total": u_total,
            "dUdt": dudt,
            "balance_residual": residual,
        }
        row.update(hp_t)
        row.update(hp_q)
        history.append(row)

        if step_count % 50 == 0 or t >= t_end:
            print(
                f"t={t:8.2f}s | T_fluid={vol.T:8.2f}K "
                f"| T_hp_con={hp_t['T_hp_con_avg']:8.2f}K "
                f"| Q_loop={q_loop_flux:8.2f}W "
                f"| Q_eva={hp_q['Q_hp_eva_in']:8.2f}W "
                f"| Q_rej={hp_q['Q_hp_solver_rej']:8.2f}W "
                f"| dUdt={dudt:8.2f}W "
                f"| residual={residual:8.2f}W"
            )

    csv_path = available_output_path(os.path.join(current_dir, "single_header_cell_one_hp.csv"))
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)

    t_arr = np.array([r["t"] for r in history])
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)

    axes[0].plot(t_arr, [r["T_fluid"] for r in history], label="fluid")
    axes[0].plot(t_arr, [r["T_header_wall"] for r in history], label="header wall")
    axes[0].plot(t_arr, [r["T_hp_eva_avg"] for r in history], label="HP evaporator")
    axes[0].plot(t_arr, [r["T_hp_con_avg"] for r in history], label="HP condenser")
    axes[0].plot(t_arr, [r["T_hp_min"] for r in history], label="HP min")
    axes[0].plot(t_arr, [r["T_fin_tip_avg"] for r in history], label="fin tip")
    axes[0].axhspan(370.0, 372.0, color="tab:red", alpha=0.12, label="Na mushy zone")
    axes[0].set_ylabel("Temperature [K]")
    axes[0].legend(loc="best", fontsize=9)
    axes[0].grid(True, linestyle="--", alpha=0.45)

    axes[1].plot(t_arr, [r["Q_loop_flux"] for r in history], label="Q_loop_flux")
    axes[1].plot(t_arr, [r["Q_hp_eva_in"] for r in history], label="Q_HP_eva_in")
    axes[1].plot(t_arr, [r["Q_hp_solver_rej"] for r in history], label="Q_HP_solver_rej")
    axes[1].plot(t_arr, [r["Q_hp_fin_rad"] for r in history], label="Q_fin_rad")
    axes[1].set_ylabel("Heat rate [W]")
    axes[1].legend(loc="best", fontsize=9)
    axes[1].grid(True, linestyle="--", alpha=0.45)

    axes[2].plot(t_arr, [r["U_fluid"] / 1000.0 for r in history], label="U_fluid")
    axes[2].plot(t_arr, [r["U_header"] / 1000.0 for r in history], label="U_header")
    axes[2].plot(t_arr, [r["U_hp"] / 1000.0 for r in history], label="U_HP")
    axes[2].set_ylabel("Stored energy [kJ]")
    axes[2].legend(loc="best", fontsize=9)
    axes[2].grid(True, linestyle="--", alpha=0.45)

    axes[3].plot(t_arr, [r["dUdt"] for r in history], label="dUdt")
    axes[3].plot(t_arr, [r["balance_residual"] for r in history], label="residual")
    axes[3].set_xlabel("Time [s]")
    axes[3].set_ylabel("Power [W]")
    axes[3].legend(loc="best", fontsize=9)
    axes[3].grid(True, linestyle="--", alpha=0.45)

    plt.tight_layout()
    plot_path = available_output_path(os.path.join(current_dir, "single_header_cell_one_hp.png"))
    plt.savefig(plot_path, dpi=130)
    plt.close()

    final = history[-1]
    print("\n[Final]")
    for key in [
        "T_fluid",
        "T_header_wall",
        "T_hp_min",
        "T_hp_min_radial_idx",
        "T_hp_min_axial_idx",
        "T_hp_eva_avg",
        "T_hp_con_avg",
        "T_fin_tip_avg",
        "Q_loop_flux",
        "Q_hp_eva_in",
        "Q_hp_solver_rej",
        "dUdt",
        "balance_residual",
    ]:
        print(f"  {key:22s}: {final[key]:.6g}")
    print(f"\nCSV saved:  {csv_path}")
    print(f"Plot saved: {plot_path}")
    print("=" * 72)


if __name__ == "__main__":
    main()
