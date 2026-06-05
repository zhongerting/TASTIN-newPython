import argparse
import os
import sys
from types import SimpleNamespace

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
from CoolantLoop.run_collector_ring_full_ringhp_geometry100hp_test import (
    apply_geometry_overrides,
    restore_overrides,
)


def parse_combo(text):
    parts = [part.strip() for part in text.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("combo must be formatted as HP_EMISSIVITY,FIN_EMISSIVITY")
    hp_eps, fin_eps = map(float, parts)
    if hp_eps < 0.0 or fin_eps < 0.0:
        raise argparse.ArgumentTypeError("emissivity must be non-negative")
    return hp_eps, fin_eps


def heat_capacity_at_reference():
    for value in (model.T_INLET, 800.0):
        try:
            return float(model.nak.heat_capacity(value))
        except TypeError:
            continue
    return float(model.nak.heat_capacity())


def build_args(base_args, hp_eps, fin_eps):
    return SimpleNamespace(
        total_flow=base_args.total_flow,
        init_temp=base_args.init_temp,
        hp_init_temp=base_args.hp_init_temp,
        hp_emissivity=hp_eps,
        fin_emissivity=fin_eps,
    )


def evaluate_combo(base_args, hp_eps, fin_eps):
    old_values = apply_geometry_overrides(build_args(base_args, hp_eps, fin_eps))
    try:
        case_model = model.build_model()
        sys_mgr = case_model["sys_mgr"]
        network = case_model["network"]
        inlet_boundary = case_model["inlet_boundary"]
        outlet_boundary = case_model["outlet_boundary"]
        ring_hp = case_model["ring_hp"]

        sys_mgr.load_global_state(base_args.restart_from)
        inlet_boundary.set_boundary_state(P=model.P_OUTLET + 5000.0, T=model.T_INLET)
        model.sync_boundary_to_network(network, inlet_boundary)
        outlet_boundary.set_boundary_state(P=model.P_OUTLET)
        model.sync_boundary_to_network(network, outlet_boundary)

        ring_hp.pre_step(0.0, sys_mgr.global_time)

        q_aba = 0.0
        q_con_total = 0.0
        q_con_bare = 0.0
        q_fin = 0.0
        hp_iter = iter(ring_hp.hp_units)
        for multiplier in ring_hp._hp_multipliers:
            if multiplier <= 0.0:
                continue
            hp = next(hp_iter)
            aba_dist, con_dist = hp.get_heat_rejection_distribution()
            breakdown = hp.get_heat_exchange_breakdown()
            scale = float(multiplier)
            q_aba += scale * float(np.sum(aba_dist))
            q_con_total += scale * float(np.sum(con_dist))
            q_con_bare += scale * float(np.sum(breakdown["bare_radiation"]))
            q_fin += scale * float(np.sum(breakdown["fin_radiation"]))

        manifolds = case_model["manifolds"]
        t_out = float(np.mean([channel.volumes[-1].T for channel in manifolds]))
        cp = heat_capacity_at_reference()
        q_total = q_aba + q_con_total
        delta_t_equiv = q_total / (float(base_args.total_flow) * cp)

        return {
            "time": float(sys_mgr.global_time),
            "hp_eps": hp_eps,
            "fin_eps": fin_eps,
            "t_out": t_out,
            "delta_t_now": float(model.T_INLET - t_out),
            "q_total_kw": q_total / 1000.0,
            "q_bare_kw": (q_aba + q_con_bare) / 1000.0,
            "q_fin_kw": q_fin / 1000.0,
            "delta_t_equiv": delta_t_equiv,
        }
    finally:
        restore_overrides(old_values)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Instantaneous emissivity sensitivity scan for the geometry100hp collector-ring case."
    )
    parser.add_argument(
        "--restart-from",
        required=True,
        help="Restart file whose current temperature field is used for the instantaneous radiation scan.",
    )
    parser.add_argument("--total-flow", type=float, default=1.3)
    parser.add_argument("--init-temp", type=float, default=model.T_INLET)
    parser.add_argument("--hp-init-temp", type=float, default=800.0)
    parser.add_argument(
        "--combo",
        action="append",
        type=parse_combo,
        required=True,
        help="HP and fin emissivity pair, for example: --combo 0.7,0.9",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    rows = [evaluate_combo(args, hp_eps, fin_eps) for hp_eps, fin_eps in args.combo]
    print(
        "time_s,hp_eps,fin_eps,T_out_K,DeltaT_now_K,"
        "Q_total_kW,Q_bare_kW,Q_fin_kW,DeltaT_equiv_K"
    )
    for row in rows:
        print(
            f"{row['time']:.6f},{row['hp_eps']:.3f},{row['fin_eps']:.3f},"
            f"{row['t_out']:.3f},{row['delta_t_now']:.3f},"
            f"{row['q_total_kw']:.3f},{row['q_bare_kw']:.3f},"
            f"{row['q_fin_kw']:.3f},{row['delta_t_equiv']:.3f}"
        )


if __name__ == "__main__":
    main()
