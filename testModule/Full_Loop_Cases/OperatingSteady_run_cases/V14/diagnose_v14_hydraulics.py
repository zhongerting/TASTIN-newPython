"""Inspect V14 hydraulic state without changing shared model code."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from Solvers.Couplers import FluidSolidCouple  # noqa: E402
from run_v14_operating_steady import (  # noqa: E402
    _apply_core_power,
    _configure_numerics,
    _set_pump_head,
    build_case,
)


def _snapshot(build: dict, dt_s: float) -> dict:
    system = build["system"]
    net = system.fluid_solver
    net._update_fluid_properties()
    net._refresh_effective_k_loss(net.W_vec)

    mass_balance = np.zeros(net.n_vol, dtype=float)
    np.add.at(mass_balance, net.idx_from_vec, -net.M_from_vec * net.W_vec)
    np.add.at(mass_balance, net.idx_to_vec, net.M_to_vec * net.W_vec)

    pressure_order = np.argsort(net.P_vec)
    imbalance_order = np.argsort(np.abs(mass_balance))[::-1]
    selected = []
    for idx, junction in enumerate(net.junctions_obj):
        name = str(getattr(junction, "name", f"junction_{idx}"))
        if (
            bool(getattr(junction, "is_pump_junction", False))
            or isinstance(junction, type(build["hot_outlet_to_ring_junctions"][0]))
            or "RingHP" in name
            or "InletMix" in name
            or "OutletMix" in name
            or "Manifold" in name
        ):
            selected.append(
                {
                    "name": name,
                    "flow_kg_s": float(net.W_vec[idx]),
                    "effective_k_loss": float(net.effective_K_loss_vec[idx]),
                    "multiplier_from": float(net.M_from_vec[idx]),
                    "multiplier_to": float(net.M_to_vec[idx]),
                    "from": str(getattr(junction.from_vol, "name", net.idx_from_vec[idx])),
                    "to": str(getattr(junction.to_vol, "name", net.idx_to_vec[idx])),
                    "target_flow_kg_s": (
                        float(getattr(junction, "target_W"))
                        if hasattr(junction, "target_W")
                        else None
                    ),
                }
            )

    def pressure_item(index: int) -> dict:
        return {
            "name": str(getattr(net.volumes_obj[index], "name", index)),
            "pressure_Pa": float(net.P_vec[index]),
        }

    return {
        "time_s": float(system.global_time),
        "dt_s": float(dt_s),
        "pressure_range_Pa": float(np.ptp(net.P_vec)),
        "pressure_low": [pressure_item(int(i)) for i in pressure_order[:5]],
        "pressure_high": [pressure_item(int(i)) for i in pressure_order[-5:][::-1]],
        "max_abs_mass_balance_kg_s": float(np.max(np.abs(mass_balance))),
        "mass_balance_top": [
            {
                "name": str(getattr(net.volumes_obj[int(i)], "name", int(i))),
                "balance_kg_s": float(mass_balance[int(i)]),
            }
            for i in imbalance_order[:15]
        ],
        "selected_junctions": selected,
        "min_flow_kg_s": float(np.min(net.W_vec)),
        "max_flow_kg_s": float(np.max(net.W_vec)),
        "nonfinite": bool(
            not np.all(np.isfinite(net.P_vec))
            or not np.all(np.isfinite(net.W_vec))
            or not np.all(np.isfinite(net.effective_K_loss_vec))
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("flow", "head"), default="flow")
    parser.add_argument("--pump-head", type=float, default=7900.0)
    parser.add_argument("--duration", type=float, default=0.0)
    parser.add_argument("--max-dt", type=float, default=0.001)
    parser.add_argument("--sample-interval", type=float, default=0.1)
    args = parser.parse_args()

    build = build_case(
        emissivity=0.75,
        pump_head_pa=float(args.pump_head),
        flow_control=args.mode == "flow",
        external_heat=False,
    )
    system = build["system"]
    for coupler in system.couplers:
        if isinstance(coupler, FluidSolidCouple):
            coupler.set_coupling_time_scheme("current")
    system.load_global_state(str(args.restart))
    _configure_numerics(build)
    _apply_core_power(build)
    if args.mode == "head":
        _set_pump_head(build, float(args.pump_head))

    snapshots = [_snapshot(build, 0.0)]
    end_time = float(system.global_time) + float(args.duration)
    next_sample = float(system.global_time) + float(args.sample_interval)
    failure = None
    while system.global_time < end_time - 1.0e-12:
        dt = min(float(args.max_dt), end_time - float(system.global_time))
        try:
            system.step(
                dt,
                inner_iter=1,
                fail_on_fluid_nonconvergence=True,
                fluid_max_iter=300,
            )
        except RuntimeError as exc:
            failure = str(exc)
            break
        _apply_core_power(build)
        if system.global_time >= next_sample - 1.0e-12:
            snapshots.append(_snapshot(build, dt))
            next_sample += float(args.sample_interval)

    if snapshots[-1]["time_s"] != float(system.global_time):
        snapshots.append(_snapshot(build, 0.0))
    payload = {
        "mode": args.mode,
        "restart": str(args.restart),
        "requested_duration_s": float(args.duration),
        "failure": failure,
        "snapshots": snapshots,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failure is not None or any(item["nonfinite"] for item in snapshots):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
