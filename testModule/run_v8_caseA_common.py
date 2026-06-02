import argparse
import os
import sys
from typing import Any, Dict, List

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from run_v7_caseA_multipliers_short import json_default
from test_core_assemble_v8_caseA import (
    _case_a_reset_design_flows_after_restart,
    build_v8_case_a_system,
)


TOTAL_POWER_W = 115000.0


def parse_v8_multipliers(text: str, *, allow_zero: bool = False) -> List[int]:
    values = [int(part.strip()) for part in str(text).split(",") if part.strip()]
    if len(values) != 5:
        raise argparse.ArgumentTypeError("Use five comma-separated integers, e.g. 1,6,12,15,3.")
    if allow_zero:
        invalid = any(value < 0 for value in values)
        message = "Multipliers must be non-negative."
    else:
        invalid = any(value <= 0 for value in values)
        message = "Multipliers must be positive."
    if invalid:
        raise argparse.ArgumentTypeError(message)
    return values


def passive_tec_source_totals(build: Dict[str, Any]) -> Dict[str, float]:
    totals = {
        "electron_cooling_flux_w_m2_sum": 0.0,
        "electron_heating_flux_w_m2_sum": 0.0,
        "emitter_joule_heat_w": 0.0,
        "collector_joule_heat_w": 0.0,
        "coupler_emitter_source_w": 0.0,
        "coupler_collector_source_w": 0.0,
    }
    for name in build.get("passive_tfe_names", []):
        tfe = build["tfes"][name]
        totals["electron_cooling_flux_w_m2_sum"] += float(np.sum(np.abs(tfe.plasma_data.electron_cooling_flux)))
        totals["electron_heating_flux_w_m2_sum"] += float(np.sum(np.abs(tfe.plasma_data.electron_heating_flux)))
        totals["emitter_joule_heat_w"] += float(np.sum(np.abs(tfe.electric_data.emitter_joule_heat)))
        totals["collector_joule_heat_w"] += float(np.sum(np.abs(tfe.electric_data.collector_joule_heat)))
        totals["coupler_emitter_source_w"] += float(np.sum(np.abs(tfe.couplers["tec_couple"].Q_source_1)))
        totals["coupler_collector_source_w"] += float(np.sum(np.abs(tfe.couplers["tec_couple"].Q_source_2)))
    return totals


def build_loaded_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v8_case_a_system(
        pipe_n_nodes=args.pipe_n_nodes,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    system.initialize_system()
    system.load_global_state(args.restart_in)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=260.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)
    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is not None:
        core.thermo_calc.calculate(verbose=False)
    core.pre_step(0.0, float(system.global_time))
    return build


__all__ = [
    "TOTAL_POWER_W",
    "build_loaded_case",
    "json_default",
    "parse_v8_multipliers",
    "passive_tec_source_totals",
]
