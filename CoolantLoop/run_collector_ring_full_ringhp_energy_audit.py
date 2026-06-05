import argparse
import csv
import os
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_full_ringhp as model


DEFAULT_RESTART_CANDIDATES = [
    "collector_ring_full_ringhp_steady_debug_503s_restart.npz",
    "collector_ring_full_ringhp_steady_debug_502s_restart.npz",
    "collector_ring_full_ringhp_steady_debug_501s_restart.npz",
    "collector_ring_full_ringhp_buffered_half_ringflow_500s_resume_from200s_restart_t0500s.npz",
    "collector_ring_full_ringhp_buffered_half_ringflow_500s_resume_from200s_restart.npz",
]


def existing_default_restart():
    for name in DEFAULT_RESTART_CANDIDATES:
        path = os.path.join(current_dir, name)
        if os.path.exists(path):
            return path
    return os.path.join(current_dir, DEFAULT_RESTART_CANDIDATES[-2])


def restart_time(restart_path):
    with np.load(restart_path, allow_pickle=True) as data:
        if "System/global_time" not in data:
            raise KeyError(f"Restart has no System/global_time: {restart_path}")
        return float(data["System/global_time"][0])


def fluid_volume_energy(vol):
    if all(hasattr(vol, attr) for attr in ("rho", "vol", "h")):
        return float(vol.rho * vol.vol * vol.h)
    return 0.0


def fluid_channel_energy(channel):
    return float(sum(fluid_volume_energy(vol) for vol in channel.volumes))


def solid_energy(solid):
    if hasattr(solid, "_update_properties"):
        solid._update_properties()
    if hasattr(solid, "thermal_capacitance") and hasattr(solid, "T"):
        return float(np.sum(solid.thermal_capacitance * solid.T))
    return 0.0


def boundary_outward_power(boundary_region):
    # BoundaryRegion.current_flux is positive into the solid. Outward rejection is negative.
    flux = getattr(boundary_region, "current_flux", None)
    if flux is None:
        return 0.0
    return float(-np.sum(flux))


def iter_hp_with_multiplier(ring_hp):
    if hasattr(ring_hp, "_iter_present_hp_units_with_multiplier"):
        yield from ring_hp._iter_present_hp_units_with_multiplier()
        return

    multipliers = np.asarray(getattr(ring_hp, "_hp_multipliers"), dtype=float)
    hp_pos = 0
    for node_index, multiplier in enumerate(multipliers):
        if multiplier <= 0:
            continue
        yield node_index, ring_hp.hp_units[hp_pos], float(multiplier)
        hp_pos += 1


def hp_storage_and_rejection(ring_hp):
    hp_storage_representative = 0.0
    hp_storage_scaled = 0.0
    hp_aba_rejection_scaled = 0.0
    hp_con_rejection_scaled = 0.0
    hp_gross_rejection_scaled = 0.0

    for _, hp_unit, multiplier in iter_hp_with_multiplier(ring_hp):
        hp_u = solid_energy(hp_unit.hp)
        hp_storage_representative += hp_u
        hp_storage_scaled += multiplier * hp_u
        q_aba = boundary_outward_power(hp_unit.hp.boundaries["outer_aba"])
        q_con = boundary_outward_power(hp_unit.hp.boundaries["outer_con"])
        hp_aba_rejection_scaled += multiplier * q_aba
        hp_con_rejection_scaled += multiplier * q_con
        hp_gross_rejection_scaled += multiplier * (q_aba + q_con)

    return {
        "U_hp_representative": hp_storage_representative,
        "U_hp_scaled": hp_storage_scaled,
        "Q_hp_aba_rejection": hp_aba_rejection_scaled,
        "Q_hp_con_rejection": hp_con_rejection_scaled,
        "Q_hp_rejection": hp_gross_rejection_scaled,
    }


def collect_energy_state(m):
    inlet_boundary = m["inlet_boundary"]
    outlet_buffer = m["outlet_buffer_channel"]
    inlet_junction = m["inlet_junction"]
    outlet_junction = m["outlet_junction"]

    U_inlet_buffer = fluid_channel_energy(m["inlet_buffer_channel"])
    U_hot_legs = float(sum(fluid_channel_energy(channel) for channel in m["hot_legs"]))
    U_ring_channel = fluid_channel_energy(m["ring_channel"])
    U_manifolds = float(sum(fluid_channel_energy(channel) for channel in m["manifolds"]))
    U_outlet_buffer = fluid_channel_energy(m["outlet_buffer_channel"])

    U_fluid_solver = (
        U_inlet_buffer
        + U_hot_legs
        + U_ring_channel
        + U_manifolds
        + U_outlet_buffer
    )
    U_fluid_macro_scaled = (
        U_inlet_buffer
        + 2.0 * U_hot_legs
        + U_ring_channel
        + 2.0 * U_manifolds
        + U_outlet_buffer
    )

    U_ring_solid = solid_energy(m["ring_solid"])
    hp_terms = hp_storage_and_rejection(m["ring_hp"])

    q_ring_wall = boundary_outward_power(m["ring_solid"].boundaries["right"])
    q_rejection_total = q_ring_wall + hp_terms["Q_hp_rejection"]

    h_in = float(inlet_boundary.h)
    h_out = float(outlet_buffer.volumes[-1].h)
    w_in = float(inlet_junction.W)
    w_out = float(outlet_junction.W)
    q_loop_flux = w_in * h_in - w_out * h_out

    q_ring_enthalpy_flux = 0.0
    for junc in m["hot_leg_to_ring"]:
        if junc.W >= 0.0:
            q_ring_enthalpy_flux += float(junc.W * junc.from_vol.h)
        else:
            q_ring_enthalpy_flux += float(junc.W * junc.to_vol.h)
    for junc in m["ring_to_manifold"]:
        if junc.W >= 0.0:
            q_ring_enthalpy_flux -= float(junc.W * junc.from_vol.h)
        else:
            q_ring_enthalpy_flux -= float(junc.W * junc.to_vol.h)

    U_total_solver_fluid = U_fluid_solver + U_ring_solid + hp_terms["U_hp_scaled"]
    U_total_macro_fluid = U_fluid_macro_scaled + U_ring_solid + hp_terms["U_hp_scaled"]
    U_ring_hp_domain = U_ring_channel + U_ring_solid + hp_terms["U_hp_scaled"]

    return {
        "time": float(m["sys_mgr"].global_time),
        "W_in_total": w_in,
        "W_out_total": w_out,
        "T_inlet_boundary": float(inlet_boundary.T),
        "T_outlet_buffer_out": float(outlet_buffer.volumes[-1].T),
        "h_in": h_in,
        "h_out": h_out,
        "Q_loop_flux": q_loop_flux,
        "Q_loop_flux_half_boundary": 0.5 * q_loop_flux,
        "Q_ring_enthalpy_flux": q_ring_enthalpy_flux,
        "Q_ring_wall_rejection": q_ring_wall,
        **hp_terms,
        "Q_rejection_total": q_rejection_total,
        "U_fluid_solver": U_fluid_solver,
        "U_fluid_macro_scaled": U_fluid_macro_scaled,
        "U_ring_solid": U_ring_solid,
        "U_total_solver_fluid": U_total_solver_fluid,
        "U_total_macro_fluid": U_total_macro_fluid,
        "U_ring_hp_domain": U_ring_hp_domain,
    }


def append_rates(current, previous, dt):
    dt_safe = max(float(dt), 1.0e-30)
    for key in ("U_fluid_solver", "U_fluid_macro_scaled", "U_ring_solid", "U_hp_scaled"):
        current[f"d{key}_dt"] = (current[key] - previous[key]) / dt_safe
    current["dU_total_solver_fluid_dt"] = (
        current["U_total_solver_fluid"] - previous["U_total_solver_fluid"]
    ) / dt_safe
    current["dU_total_macro_fluid_dt"] = (
        current["U_total_macro_fluid"] - previous["U_total_macro_fluid"]
    ) / dt_safe
    current["dU_ring_hp_domain_dt"] = (
        current["U_ring_hp_domain"] - previous["U_ring_hp_domain"]
    ) / dt_safe
    current["Q_balance_residual_solver_fluid"] = (
        current["Q_loop_flux"]
        - current["Q_rejection_total"]
        - current["dU_total_solver_fluid_dt"]
    )
    current["Q_balance_residual_macro_fluid"] = (
        current["Q_loop_flux"]
        - current["Q_rejection_total"]
        - current["dU_total_macro_fluid_dt"]
    )
    current["Q_balance_residual_half_boundary"] = (
        current["Q_loop_flux_half_boundary"]
        - current["Q_rejection_total"]
        - current["dU_total_solver_fluid_dt"]
    )
    current["Q_balance_residual_ring_domain"] = (
        current["Q_ring_enthalpy_flux"]
        - current["Q_rejection_total"]
        - current["dU_ring_hp_domain_dt"]
    )
    current["relative_residual_solver_fluid"] = current["Q_balance_residual_solver_fluid"] / max(
        abs(current["Q_loop_flux"]), 1.0
    )
    current["relative_residual_macro_fluid"] = current["Q_balance_residual_macro_fluid"] / max(
        abs(current["Q_loop_flux"]), 1.0
    )
    current["relative_residual_half_boundary"] = current["Q_balance_residual_half_boundary"] / max(
        abs(current["Q_loop_flux_half_boundary"]), 1.0
    )
    current["relative_residual_ring_domain"] = current["Q_balance_residual_ring_domain"] / max(
        abs(current["Q_ring_enthalpy_flux"]), 1.0
    )


def write_csv(csv_path, rows):
    if not rows:
        return
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Energy audit CSV saved: {csv_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Energy-conservation audit for the full-ring collector-ring + RingHP case."
    )
    parser.add_argument("--restart-from", default=existing_default_restart())
    parser.add_argument("--duration", type=float, default=1.0)
    parser.add_argument("--target-time", type=float, default=None)
    parser.add_argument("--case-name", default="collector_ring_full_ringhp_energy_audit")
    parser.add_argument("--csv-path", default=None)
    parser.add_argument("--min-dt", type=float, default=1.0e-3)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--safety-factor", type=float, default=1.0)
    parser.add_argument("--inner-iter", type=int, default=2)
    parser.add_argument("--print-every", type=float, default=1.0)
    return parser.parse_args()


def main():
    args = parse_args()
    restart_path = os.path.abspath(args.restart_from)
    if not os.path.exists(restart_path):
        raise FileNotFoundError(restart_path)

    start_time = restart_time(restart_path)
    target_time = args.target_time if args.target_time is not None else start_time + args.duration
    if target_time <= start_time:
        raise ValueError(f"Target time must be greater than restart time: {target_time} <= {start_time}")

    m = model.build_model()
    sys_mgr = m["sys_mgr"]
    sys_mgr.load_global_state(restart_path)
    m["inlet_boundary"].set_boundary_state(P=model.P_OUTLET + 5000.0, T=model.T_INLET)
    model.sync_boundary_to_network(m["network"], m["inlet_boundary"])
    m["outlet_boundary"].set_boundary_state(P=model.P_OUTLET)
    model.sync_boundary_to_network(m["network"], m["outlet_boundary"])

    print("=" * 78)
    print(f"Energy audit case: {args.case_name}")
    print(f"Restart loaded: {restart_path}")
    print(f"Restart time: {sys_mgr.global_time:.6f} s")
    print(f"Target time : {target_time:.6f} s")
    print("=" * 78)

    previous = collect_energy_state(m)
    previous_time = previous["time"]
    rows = []
    next_print = model.next_event_time(sys_mgr.global_time, args.print_every)

    while sys_mgr.global_time < target_time:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=args.min_dt,
            max_dt=args.max_dt,
            safety_factor=args.safety_factor,
        )
        if next_print is not None and sys_mgr.global_time < next_print < sys_mgr.global_time + dt:
            dt = next_print - sys_mgr.global_time
        dt = min(dt, target_time - sys_mgr.global_time)

        sys_mgr.step(dt=dt, inner_iter=args.inner_iter)

        current = collect_energy_state(m)
        dt_energy = current["time"] - previous_time
        append_rates(current, previous, dt_energy)
        rows.append(current)
        previous = current
        previous_time = current["time"]

        if next_print is not None and current["time"] >= next_print - 1.0e-12:
            print(
                f"t = {current['time']:9.6f} s | "
                f"Q_loop = {current['Q_loop_flux'] / 1000.0:9.3f} kW | "
                f"Q_rej = {current['Q_rejection_total'] / 1000.0:9.3f} kW | "
                f"dUdt = {current['dU_total_solver_fluid_dt'] / 1000.0:9.3f} kW | "
                f"res = {current['Q_balance_residual_solver_fluid'] / 1000.0:9.3f} kW"
            )
            while next_print is not None and current["time"] >= next_print - 1.0e-12:
                next_print += args.print_every

    csv_path = args.csv_path
    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{args.case_name}.csv")
    write_csv(csv_path, rows)

    final = rows[-1]
    print("-" * 78)
    print("Final Energy Balance")
    print(f"  time                         : {final['time']:.9f} s")
    print(f"  T_outlet_buffer_out          : {final['T_outlet_buffer_out']:.9f} K")
    print(f"  W_in / W_out                 : {final['W_in_total']:.9f} / {final['W_out_total']:.9f} kg/s")
    print(f"  Q_loop_flux                  : {final['Q_loop_flux'] / 1000.0:.9f} kW")
    print(f"  Q_loop_flux_half_boundary    : {final['Q_loop_flux_half_boundary'] / 1000.0:.9f} kW")
    print(f"  Q_ring_enthalpy_flux         : {final['Q_ring_enthalpy_flux'] / 1000.0:.9f} kW")
    print(f"  Q_ring_wall_rejection        : {final['Q_ring_wall_rejection'] / 1000.0:.9f} kW")
    print(f"  Q_hp_rejection               : {final['Q_hp_rejection'] / 1000.0:.9f} kW")
    print(f"  Q_rejection_total            : {final['Q_rejection_total'] / 1000.0:.9f} kW")
    print(f"  dU/dt solver-fluid           : {final['dU_total_solver_fluid_dt'] / 1000.0:.9f} kW")
    print(f"  residual solver-fluid        : {final['Q_balance_residual_solver_fluid'] / 1000.0:.9f} kW")
    print(f"  relative residual solver     : {final['relative_residual_solver_fluid']:.9e}")
    print(f"  dU/dt macro-fluid            : {final['dU_total_macro_fluid_dt'] / 1000.0:.9f} kW")
    print(f"  residual macro-fluid         : {final['Q_balance_residual_macro_fluid'] / 1000.0:.9f} kW")
    print(f"  relative residual macro      : {final['relative_residual_macro_fluid']:.9e}")
    print(f"  dU/dt ring+HP domain         : {final['dU_ring_hp_domain_dt'] / 1000.0:.9f} kW")
    print(f"  residual half-boundary       : {final['Q_balance_residual_half_boundary'] / 1000.0:.9f} kW")
    print(f"  relative residual half       : {final['relative_residual_half_boundary']:.9e}")
    print(f"  residual ring+HP domain      : {final['Q_balance_residual_ring_domain'] / 1000.0:.9f} kW")
    print(f"  relative residual ring+HP    : {final['relative_residual_ring_domain']:.9e}")
    print("-" * 78)


if __name__ == "__main__":
    main()
