import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from run_v8_caseA_common import (
    DEFAULT_COOLANT_MATERIAL,
    DEFAULT_SOLID_ODE_METHOD,
    TOTAL_POWER_W,
    apply_solid_ode_method,
    apply_wire_resistance,
    get_solid_ode_methods,
    get_wire_resistance,
    json_default,
    parse_solid_ode_method,
    parse_v8_multipliers,
)
from test_core_assemble_v8_caseA import (
    _case_a_reset_design_flows_after_restart,
    build_v8_case_a_system,
)


DEFAULT_RESTART_IN = (
    "testModule/v8_caseA_lsoda_wire_2000s/"
    "v8_caseA_lsoda_wire_2000s_latest_restart.npz"
)
DEFAULT_OUTPUT_DIR = "testModule/v8_caseA_nak_migrated"
DEFAULT_CASE_PREFIX = "v8_caseA_nak_migrated"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Migrate a V8 CaseA Sodium restart to the NaK78 coolant definition."
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART_IN)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--case-prefix", default=DEFAULT_CASE_PREFIX)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument(
        "--solid-ode-method",
        type=parse_solid_ode_method,
        default=DEFAULT_SOLID_ODE_METHOD,
    )
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=False),
        default=parse_v8_multipliers("1,6,12,15,3"),
    )
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=True),
        default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True),
    )
    return parser.parse_args()


def recompute_fluid_enthalpy_for_current_material(system: Any) -> Dict[str, float]:
    network = system.fluid_solver
    h_old = np.asarray(network.h_vec, dtype=float).copy()
    t_vec = np.asarray(network.T_vec, dtype=float)
    p_vec = np.asarray(network.P_vec, dtype=float)
    h_new = np.empty_like(h_old)
    for idx, vol in enumerate(network.volumes_obj):
        h_new[idx] = float(vol.material.enthalpy(float(t_vec[idx]), float(p_vec[idx])))
    network.h_vec[:] = h_new
    network._sync_vectors_to_objects(sync_pressure=True, sync_flow=True, sync_energy=True, sync_properties=False)
    network._update_fluid_properties()
    network._sync_vectors_to_objects()
    return {
        "h_old_min_j_kg": float(np.min(h_old)),
        "h_old_max_j_kg": float(np.max(h_old)),
        "h_new_min_j_kg": float(np.min(h_new)),
        "h_new_max_j_kg": float(np.max(h_new)),
        "h_delta_min_j_kg": float(np.min(h_new - h_old)),
        "h_delta_max_j_kg": float(np.max(h_new - h_old)),
    }


def build_migrated_case(args: argparse.Namespace) -> Dict[str, Any]:
    solid_ode_method = parse_solid_ode_method(args.solid_ode_method)
    build = build_v8_case_a_system(
        pipe_n_nodes=args.pipe_n_nodes,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=DEFAULT_COOLANT_MATERIAL,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    apply_solid_ode_method(build, solid_ode_method)
    system.initialize_system()
    system.load_global_state(args.restart_in)

    enthalpy_summary = recompute_fluid_enthalpy_for_current_material(system)
    apply_solid_ode_method(build, solid_ode_method)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=150.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)
    core.post_step(0.0, float(system.global_time))
    apply_wire_resistance(core)
    core.pre_step(0.0, float(system.global_time))

    build["solid_ode_method"] = solid_ode_method
    build["solid_ode_methods"] = get_solid_ode_methods(build)
    build["wire_resistance_ohm"] = get_wire_resistance(core)
    build["coolant_material"] = build.get("coolant_material", DEFAULT_COOLANT_MATERIAL)
    build["migration_enthalpy_summary"] = enthalpy_summary
    return build


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    restart_out = output_dir / f"{args.case_prefix}_latest_restart.npz"
    summary_out = output_dir / f"{args.case_prefix}_migration_summary.json"

    build = build_migrated_case(args)
    system = build["system"]
    system.save_global_state(str(restart_out))

    latest = {
        "restart_in": args.restart_in,
        "restart_out": str(restart_out),
        "summary_out": str(summary_out),
        "absolute_time_s": float(system.global_time),
        "coolant_material": build["coolant_material"],
        "solid_ode_method": build["solid_ode_method"],
        "solid_ode_methods": build["solid_ode_methods"],
        "wire_resistance_ohm": build["wire_resistance_ohm"],
        "ring_multipliers": build["ring_multipliers"],
        "tec_ring_multipliers": build["tec_ring_multipliers"],
        "fluid_enthalpy_recomputed_from_temperature": True,
        "migration_enthalpy_summary": build["migration_enthalpy_summary"],
    }
    with summary_out.open("w", encoding="utf-8") as f:
        json.dump(latest, f, indent=2, ensure_ascii=False, default=json_default)

    print("=== V8 CaseA Sodium -> NaK78 migration completed ===", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"restart_out={restart_out}", flush=True)
    print(f"absolute_time_s={system.global_time:.6f}", flush=True)
    print(f"coolant_material={build['coolant_material']}", flush=True)
    print(f"wire_resistance_ohm={build['wire_resistance_ohm']}", flush=True)


if __name__ == "__main__":
    main()
