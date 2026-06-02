import argparse
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from run_v7_caseA_multipliers_short import build_loaded_case as build_loaded_v7
from run_v8_caseA_common import json_default, passive_tec_source_totals
from test_core_assemble_v8_caseA import (
    V8_RING_MULTIPLIERS,
    V8_TEC_RING_MULTIPLIERS,
    _case_a_reset_design_flows_after_restart,
    build_v8_case_a_system,
)


DEFAULT_V7_RESTART = (
    "testModule/v7_caseA_heat_chain_continue100/"
    "v7_caseA_heat_chain_continue100_latest_restart.npz"
)
DEFAULT_OUTPUT = "testModule/v8_caseA_migrated/v8_caseA_migrated_latest_restart.npz"
TOTAL_POWER_W = 115000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Migrate a V7 CaseA restart into native V8 topology.")
    parser.add_argument("--restart-in", default=DEFAULT_V7_RESTART)
    parser.add_argument("--restart-out", default=DEFAULT_OUTPUT)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    return parser.parse_args()


def source_name_for_v8(name: str) -> str:
    return (
        str(name)
        .replace("Ring3_TEC", "Ring3")
        .replace("Ring3_Open", "Ring3")
        .replace("TASTIN_Core_V8_CaseA", "TASTIN_Core_V7_CaseA")
    )


def copy_solid_states(v7: Dict[str, Any], v8: Dict[str, Any]) -> None:
    source_solids = v7["system"].solid_components
    for target_name, target in v8["system"].solid_components.items():
        source_name = source_name_for_v8(target_name)
        if source_name not in source_solids:
            raise KeyError(f"V7 solid state not found for V8 target '{target_name}' via '{source_name}'.")
        source = source_solids[source_name]
        if target.T.shape != source.T.shape:
            raise ValueError(f"Solid shape mismatch for '{target_name}': {target.T.shape} vs {source.T.shape}.")
        target.T[:] = source.T
        target.dTdt[:] = source.dTdt


def copy_tfe_macro_states(v7: Dict[str, Any], v8: Dict[str, Any]) -> None:
    for target_name, target in v8["tfes"].items():
        source = v7["tfes"][source_name_for_v8(target_name)]
        target.neutronic_data.total_power = float(source.neutronic_data.total_power)
        target.neutronic_data._total_power_old = float(source.neutronic_data._total_power_old)
        for attr in vars(target.electric_data):
            getattr(target.electric_data, attr)[...] = getattr(source.electric_data, attr)
        for attr in vars(target.plasma_data):
            getattr(target.plasma_data, attr)[...] = getattr(source.plasma_data, attr)
        target.boundary_data.moderator_temperature[:] = source.boundary_data.moderator_temperature


def copy_fluid_states(v7: Dict[str, Any], v8: Dict[str, Any]) -> None:
    source_net = v7["system"].fluid_solver
    target_net = v8["system"].fluid_solver
    source_volumes = {vol.name: vol for vol in source_net.volumes_obj}
    source_junctions = {junc.name: junc for junc in source_net.junctions_obj}

    for target in target_net.volumes_obj:
        source_name = source_name_for_v8(target.name)
        if source_name not in source_volumes:
            raise KeyError(f"V7 fluid volume not found for V8 target '{target.name}' via '{source_name}'.")
        source = source_volumes[source_name]
        target.P = float(source.P)
        target.T = float(source.T)
        target.h = float(source.h)
        if hasattr(target, "target_P") and hasattr(source, "target_P"):
            target.target_P = float(source.target_P)

    for target in target_net.junctions_obj:
        source_name = source_name_for_v8(target.name)
        if source_name not in source_junctions:
            raise KeyError(f"V7 junction not found for V8 target '{target.name}' via '{source_name}'.")
        source = source_junctions[source_name]
        target.W = float(source.W)
        if hasattr(target, "target_W") and hasattr(source, "target_W"):
            target.target_W = float(source.target_W)

    target_net._initialize_state_from_objects()
    target_net._update_fluid_properties()
    target_net._sync_vectors_to_objects()
    if hasattr(target_net, "W_old"):
        target_net.W_old[:] = target_net.W_vec
    if hasattr(target_net, "W_iterate"):
        target_net.W_iterate[:] = target_net.W_vec
    target_net._refresh_cached_boundary_targets()


def migrate(args: argparse.Namespace) -> Dict[str, Any]:
    v7_args = SimpleNamespace(
        restart_in=args.restart_in,
        pipe_n_nodes=args.pipe_n_nodes,
        target_voltage=args.target_voltage,
        ring_multipliers=[1, 6, 12, 18],
        tec_ring_multipliers=[1, 6, 12, 15],
    )
    v7 = build_loaded_v7(v7_args)
    v8 = build_v8_case_a_system(
        pipe_n_nodes=args.pipe_n_nodes,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        ring_multipliers=V8_RING_MULTIPLIERS,
        tec_ring_multipliers=V8_TEC_RING_MULTIPLIERS,
    )
    system = v8["system"]
    core = v8["core"]
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    system.initialize_system()

    copy_solid_states(v7, v8)
    copy_tfe_macro_states(v7, v8)
    copy_fluid_states(v7, v8)
    system.global_time = float(v7["system"].global_time)
    system._sync_solid_times_to_global()

    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=260.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(v8)
    system._refresh_solid_boundary_cache(update_flux=True, current_time=float(system.global_time))
    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is not None:
        core.thermo_calc.calculate(verbose=False)
    core.pre_step(0.0, float(system.global_time))
    system._prepare_fluid_sources_for_coupling()
    system._refresh_solid_boundary_cache(update_flux=True, current_time=float(system.global_time))
    system._run_couplers(interface_relaxation=1.0, current_time=float(system.global_time))
    system._refresh_solid_boundary_cache(update_flux=True, current_time=float(system.global_time))

    passive_totals = passive_tec_source_totals(v8)
    if any(value != 0.0 for value in passive_totals.values()):
        raise RuntimeError(f"Ring3_Open TEC sources are not zero after migration: {passive_totals}")

    restart_out = Path(args.restart_out)
    restart_out.parent.mkdir(parents=True, exist_ok=True)
    system.save_global_state(str(restart_out))
    summary = {
        "restart_in": args.restart_in,
        "restart_out": str(restart_out),
        "global_time_s": float(system.global_time),
        "v7_fluid_shape": [v7["system"].fluid_solver.n_vol, v7["system"].fluid_solver.n_junc],
        "v8_fluid_shape": [system.fluid_solver.n_vol, system.fluid_solver.n_junc],
        "v8_ring_multipliers": v8["ring_multipliers"],
        "v8_tec_ring_multipliers": v8["tec_ring_multipliers"],
        "v8_physical_ring_count": v8["physical_ring_count"],
        "v8_ring3_members": core.get_ring_member_names(3),
        "ring3_open_passive_tec_source_totals": passive_totals,
    }
    summary_path = restart_out.with_name(f"{restart_out.stem}_migration_summary.json")
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return summary


def main() -> None:
    summary = migrate(parse_args())
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))


if __name__ == "__main__":
    main()
