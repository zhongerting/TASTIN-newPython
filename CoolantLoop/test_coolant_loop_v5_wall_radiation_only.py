import csv
import logging
import os
import sys
from datetime import datetime

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
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
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


# ============================================================================
# Geometry and boundary conditions, matched to test_coolant_loop_v5.py
# ============================================================================
T_SPACE = 3.0
T_INLET = 900.0
W_TOTAL = 2.2
W_BRANCH = W_TOTAL / 3.0
W_RING = W_BRANCH / 2.0 / 3.0
P_OUTLET = 160000.0
T_INIT = 900.0

L_RING = 0.793
N_RING_NODES = 10
R_IN_RING = 0.020835
R_OUT_RING = 0.022835
DIAM_RING = 0.04167
AREA_RING = 0.0016065
R_IN_HEADER = DIAM_RING / 2.0
PERIM_HEADER = 2.0 * np.pi * R_IN_HEADER
WALL_THICKNESS_RING = R_OUT_RING - R_IN_RING

# Dummy heat-pipe parameters required by the RingHP constructor. With
# HP_MULTIPLIERS = 0 for every node, no HPwithFin units are created.
R_OUT_HP = 0.0085
R_IN_HP = 0.0081
R_VAPOR_HP = 0.0075
L_EVA = 0.0605
L_ABA = 0.0415
L_CON = 0.47
POROSITY = 0.966
THIN_FIN = 0.0003
FIN_HEIGHT = 22.65e-3
N_FIN_HEIGHT = 15

EMISSIVITY_RING_WALL = 0.9
BARE_AREA_FACTOR = 0.5
HP_MULTIPLIERS_ZERO = [0.0] * N_RING_NODES
DEFAULT_OUTLET_BUFFER_NODES = 2
DEFAULT_OUTLET_BUFFER_LENGTH = L_RING / N_RING_NODES * DEFAULT_OUTLET_BUFFER_NODES


def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    pe = np.maximum(Re * Pr, 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


def boundary_outward_heat_by_node(boundary_region):
    """Return outward heat loss by boundary node [W], positive to environment."""
    q_node = np.zeros(boundary_region.shape, dtype=float)
    for condition in boundary_region.conditions:
        q_flux = getattr(condition, "last_q_flux", getattr(condition, "q_flux", None))
        if q_flux is not None:
            # Boundary convention is positive into the solid. Radiation loss is
            # therefore negative in q_flux, so outward heat is -q_flux.
            q_node += -np.array(q_flux, dtype=float)
    return q_node.reshape(-1)


def fluid_domain_energy(channel):
    return float(sum(vol.rho * vol.vol * vol.h for vol in channel.volumes))


def solid_domain_energy(solid):
    if hasattr(solid, "_update_properties"):
        solid._update_properties()
    return float(np.sum(solid.thermal_capacitance * solid.T))


def write_history_csv(csv_path, history):
    if not history:
        print("[WARN] No history to write.")
        return
    try:
        f = open(csv_path, "w", newline="", encoding="utf-8-sig")
    except PermissionError:
        stem, ext = os.path.splitext(csv_path)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_path = f"{stem}_{stamp}{ext}"
        f = open(csv_path, "w", newline="", encoding="utf-8-sig")

    with f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    print(f"CSV saved: {csv_path}")


def build_case(
    with_outlet_buffer=False,
    outlet_buffer_nodes=DEFAULT_OUTLET_BUFFER_NODES,
    outlet_buffer_length=DEFAULT_OUTLET_BUFFER_LENGTH,
):
    nak = SodiumPotassium78()
    mat_wall = SS316(name="SS316_wall")
    mat_hp_fluid = SodiumHP(name="HP_Fluid_Na")
    mat_wick = WickMaterial(
        name="WickMaterial",
        solid_mat=mat_wall,
        fluid_mat=mat_hp_fluid,
        porosity=POROSITY,
        r_vapor=R_VAPOR_HP,
        r_in_wall=R_IN_HP,
    )

    inlet_boundary = IncompressibleBoundaryVolume(
        name="WallRadOnly_Inlet_FlowBoundary",
        material=nak,
        P=P_OUTLET + 5000.0,
        T=T_INLET,
    )
    outlet_boundary = IncompressibleBoundaryVolume(
        name="WallRadOnly_Outlet_PressureBoundary",
        material=nak,
        P=P_OUTLET,
        T=T_INIT,
    )
    outlet_boundary.is_pressure_boundary = True

    chan_ring = IncompressibleFluidChannel(
        name="WallRadOnly_Ring_Channel",
        n_nodes=N_RING_NODES,
        total_length=L_RING,
        flow_area=AREA_RING,
        hydraulic_diam=DIAM_RING,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=nak,
    )

    mesh_ring = Mesh2D(
        x_dim=WALL_THICKNESS_RING,
        n_x=1,
        y_dim=L_RING,
        n_y=N_RING_NODES,
        geometry_type="cylindrical",
        inner_radius=R_IN_RING,
    )
    solid_ring = HeatConduction2D(
        mesh=mesh_ring,
        material=mat_wall,
        name="WallRadOnly_Solid_Ring",
        initial_temp=T_INIT,
    )

    # Only the outer wall radiates to space; other boundaries are adiabatic
    # except the left boundary, which is controlled by the header fluid coupler.
    bare_area_array = solid_ring.boundaries["right"].area * BARE_AREA_FACTOR
    solid_ring.boundaries["right"].add_dynamic_radiation_condition(
        emissivity=EMISSIVITY_RING_WALL,
        bare_area_array=bare_area_array,
        T_env=T_SPACE,
    )
    solid_ring.boundaries["top"].add_flux_condition(q_flux=0.0)
    solid_ring.boundaries["bottom"].add_flux_condition(q_flux=0.0)

    # RingHP is used here only as a header fluid-solid coupling wrapper.
    # HP_MULTIPLIERS_ZERO disables all heat pipes.
    ring_header_only = RingHP(
        name="WallRadOnly_HeaderCoupler",
        fluid_channel=chan_ring,
        solid_header=solid_ring,
        hp_multipliers=HP_MULTIPLIERS_ZERO,
        header_flow_area=AREA_RING,
        header_dh=DIAM_RING,
        header_heated_perimeter=PERIM_HEADER,
        hp_r_out=R_OUT_HP,
        hp_r_in=R_IN_HP,
        hp_r_vapor=R_VAPOR_HP,
        hp_L_eva=L_EVA,
        hp_L_con=L_CON,
        hp_L_aba=L_ABA,
        hp_n_eva=1,
        hp_n_con=1,
        hp_n_aba=12,
        hp_n_wick=1,
        hp_n_wall=2,
        porosity_hp=POROSITY,
        HP_initial_temp=800.0,
        fin_thickness=THIN_FIN,
        fin_height=FIN_HEIGHT,
        n_fin_height=N_FIN_HEIGHT,
        fin_wrap_ratio=(2.0 * THIN_FIN) / (2.0 * np.pi * R_OUT_HP),
        emissivity=0.9,
        up_view_factor=0.0,
        down_view_factor=0.3,
        T_space=T_SPACE,
        hp_wall_mat=mat_wall,
        hp_fluid_mat=mat_hp_fluid,
        hp_wick_mat=mat_wick,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=None,
    )

    j_in = InletJunction(
        name="J_WallRadOnly_Inlet_Ring",
        from_vol=inlet_boundary,
        to_vol=chan_ring.volumes[0],
        W_initial=W_RING,
    )
    outlet_buffer = None
    j_ring_to_buffer = None
    if with_outlet_buffer:
        outlet_buffer = IncompressibleFluidChannel(
            name="WallRadOnly_Adiabatic_Outlet_Buffer",
            n_nodes=outlet_buffer_nodes,
            total_length=outlet_buffer_length,
            flow_area=AREA_RING,
            hydraulic_diam=DIAM_RING,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        j_ring_to_buffer = FlowJunction(
            name="J_WallRadOnly_Ring_Buffer",
            from_vol=chan_ring.volumes[-1],
            to_vol=outlet_buffer.volumes[0],
            flow_area=AREA_RING,
            k_loss=0.0,
        )
        j_out_from_vol = outlet_buffer.volumes[-1]
    else:
        j_out_from_vol = chan_ring.volumes[-1]

    j_out = FlowJunction(
        name="J_WallRadOnly_Outlet_Boundary",
        from_vol=j_out_from_vol,
        to_vol=outlet_boundary,
        flow_area=AREA_RING,
        k_loss=0.0,
    )

    all_vols = [inlet_boundary] + chan_ring.volumes
    all_juncs = [j_in] + chan_ring.internal_junctions
    if with_outlet_buffer:
        all_vols += outlet_buffer.volumes
        all_juncs += [j_ring_to_buffer] + outlet_buffer.internal_junctions
    all_vols += [outlet_boundary]
    all_juncs += [j_out]
    network = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)

    sys_mgr = SystemManager(fluid_network=network)
    sys_mgr.add_component(ring_header_only)

    return {
        "sys_mgr": sys_mgr,
        "channel": chan_ring,
        "solid_ring": solid_ring,
        "j_in": j_in,
        "j_out": j_out,
        "j_ring_to_buffer": j_ring_to_buffer,
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "outlet_buffer": outlet_buffer,
        "ring_header_only": ring_header_only,
    }


def run_case(
    case_name="wall_radiation_only_10node",
    t_end=1000.0,
    min_dt=1e-3,
    max_dt=0.5,
    safety_factor=0.5,
    inner_iter=2,
    print_every=100,
    csv_path=None,
    with_outlet_buffer=False,
    outlet_buffer_nodes=DEFAULT_OUTLET_BUFFER_NODES,
    outlet_buffer_length=DEFAULT_OUTLET_BUFFER_LENGTH,
):
    case = build_case(
        with_outlet_buffer=with_outlet_buffer,
        outlet_buffer_nodes=outlet_buffer_nodes,
        outlet_buffer_length=outlet_buffer_length,
    )
    sys_mgr = case["sys_mgr"]
    chan_ring = case["channel"]
    solid_ring = case["solid_ring"]
    j_in = case["j_in"]
    j_out = case["j_out"]
    outlet_buffer = case["outlet_buffer"]

    print("=" * 72)
    print(f"Running V5 wall-radiation-only test: {case_name}")
    print("=" * 72)
    print(f"Ring nodes: {N_RING_NODES}")
    print(f"Heat-pipe multipliers: {HP_MULTIPLIERS_ZERO}")
    print(f"Inlet T: {T_INLET:.3f} K")
    print(f"Inlet W: {W_RING:.9f} kg/s")
    print(f"Outlet P: {P_OUTLET:.3f} Pa")
    print(f"Ring length: {L_RING:.6f} m")
    print(f"Bare outer wall area factor: {BARE_AREA_FACTOR:.3f}")
    if outlet_buffer is not None:
        print(
            f"Outlet buffer: {outlet_buffer_nodes} adiabatic nodes, "
            f"length={outlet_buffer_length:.6f} m"
        )

    sys_mgr.initialize_system()

    history = []
    prev_u = None
    prev_t = sys_mgr.global_time
    step_count = 0

    while sys_mgr.global_time < t_end:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=min_dt,
            max_dt=max_dt,
            safety_factor=safety_factor,
        )
        dt = min(dt, t_end - sys_mgr.global_time)
        sys_mgr.step(dt=dt, inner_iter=inner_iter)

        step_count += 1
        t = sys_mgr.global_time

        q_rad_nodes = boundary_outward_heat_by_node(solid_ring.boundaries["right"])
        q_rad_total = float(np.sum(q_rad_nodes))

        q_fluid_enthalpy = j_in.W * j_in.from_vol.h - j_out.W * j_out.from_vol.h
        delta_h = j_in.from_vol.h - j_out.from_vol.h

        u_fluid = fluid_domain_energy(chan_ring)
        u_solid = solid_domain_energy(solid_ring)
        u_total = u_fluid + u_solid
        if prev_u is None:
            dUdt = 0.0
        else:
            dUdt = (u_total - prev_u) / max(t - prev_t, 1.0e-30)
        prev_u = u_total
        prev_t = t

        residual = q_fluid_enthalpy - q_rad_total - dUdt
        residual_rel = residual / max(abs(q_fluid_enthalpy), abs(q_rad_total), 1.0)

        row = {
            "time": t,
            "dt": dt,
            "W_in": j_in.W,
            "W_out": j_out.W,
            "T_inlet_boundary": j_in.from_vol.T,
            "T_outlet_fluid": j_out.from_vol.T,
            "T_ring_exit_fluid": chan_ring.volumes[-1].T,
            "h_inlet_boundary": j_in.from_vol.h,
            "h_outlet_fluid": j_out.from_vol.h,
            "h_ring_exit_fluid": chan_ring.volumes[-1].h,
            "delta_h": delta_h,
            "Q_fluid_enthalpy": q_fluid_enthalpy,
            "Q_wall_rad_total": q_rad_total,
            "U_fluid": u_fluid,
            "U_solid_ring": u_solid,
            "U_total": u_total,
            "dUdt_total": dUdt,
            "energy_residual": residual,
            "energy_residual_rel": residual_rel,
            "T_fluid_avg": float(np.mean(chan_ring.temperature_vector)),
            "T_wall_avg": float(np.mean(solid_ring.T)),
        }

        for i, q in enumerate(q_rad_nodes):
            row[f"Q_wall_rad_node_{i + 1:02d}"] = float(q)
        for i, vol in enumerate(chan_ring.volumes):
            row[f"T_fluid_node_{i + 1:02d}"] = float(vol.T)
            row[f"P_fluid_node_{i + 1:02d}"] = float(vol.P)
        for i, temp in enumerate(solid_ring.T.reshape(-1)):
            row[f"T_wall_node_{i + 1:02d}"] = float(temp)
        if outlet_buffer is not None:
            for i, vol in enumerate(outlet_buffer.volumes):
                row[f"T_buffer_node_{i + 1:02d}"] = float(vol.T)
                row[f"P_buffer_node_{i + 1:02d}"] = float(vol.P)
                row[f"h_buffer_node_{i + 1:02d}"] = float(vol.h)

        history.append(row)

        if step_count % print_every == 0 or t >= t_end:
            print(
                f"t={t:8.3f}s | T_out={row['T_outlet_fluid']:.3f} K | "
                f"Q_h={q_fluid_enthalpy:.3f} W | Q_rad={q_rad_total:.3f} W | "
                f"dUdt={dUdt:.3f} W | residual={residual:.3f} W "
                f"({residual_rel:.3e})"
            )

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}.csv")
    write_history_csv(csv_path, history)

    final = history[-1]
    print("-" * 72)
    print("Final energy check")
    print(f"Q_fluid_enthalpy: {final['Q_fluid_enthalpy']:.6f} W")
    print(f"Q_wall_rad_total: {final['Q_wall_rad_total']:.6f} W")
    print(f"dUdt_total:       {final['dUdt_total']:.6f} W")
    print(f"residual:         {final['energy_residual']:.6f} W")
    print(f"relative residual:{final['energy_residual_rel']:.6e}")
    print(f"CSV: {csv_path}")
    print("=" * 72)

    return history


def run_four_dt_iter_cases(t_end=1000.0):
    cases = [
        ("dt0p5_iter2", 0.5, 2),
        ("dt0p1_iter2", 0.1, 2),
        ("dt0p1_iter5", 0.1, 5),
        ("dt0p05_iter5", 0.05, 5),
    ]

    results = {}
    for suffix, max_dt, inner_iter in cases:
        case_name = f"V5_wall_radiation_only_10node_{suffix}"
        csv_path = os.path.join(current_dir, f"{case_name}.csv")
        results[suffix] = run_case(
            case_name=case_name,
            t_end=t_end,
            max_dt=max_dt,
            inner_iter=inner_iter,
            print_every=1000,
            csv_path=csv_path,
        )
    return results


if __name__ == "__main__":
    run_four_dt_iter_cases(t_end=1000.0)
