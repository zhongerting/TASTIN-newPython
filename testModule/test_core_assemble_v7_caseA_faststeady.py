import os
import sys
from typing import Any, Dict

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from test_core_assemble_v7_caseA import (
    build_v7_case_a_system,
    run_test_v7_case_a_heated,
    _case_a_electric_diagnostics,
    _case_a_flow_diagnostics,
)


FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE = 0.01
FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE_SCOPE = "global_outer"
FASTSTEADY_RESTART_IN = "test_core_assemble_v7_caseA_heated_restart_t7800.npz"
FASTSTEADY_RESTART_OUT = "test_core_assemble_v7_caseA_faststeady_restart_t8800.npz"
FASTSTEADY_RUN_DURATION_S = 1000.0
FASTSTEADY_TOTAL_POWER_W = 115000.0


def compute_faststeady_energy_audit(result: Dict[str, Any]) -> Dict[str, float]:
    build = result["build"]
    system = result["system"]
    core = result["core"]
    sigma = 5.670374419e-8

    flow = _case_a_flow_diagnostics(build)
    mass_flow = float(flow["inlet_total_macro_flow_kg_s"])
    h_in = float(build["inlet_plenum"].h)
    h_out = float(build["outlet_plenum"].h)
    q_fluid = mass_flow * (h_out - h_in)

    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is not None:
        core.thermo_calc.calculate(verbose=False)
    electric = _case_a_electric_diagnostics(core)
    q_electric = float(electric["tec_total_electric_power_w"] or 0.0)

    reflector_outer = core.reflector.boundaries["right"]
    surface_temperature = np.asarray(reflector_outer.T_surface, dtype=float)
    radiation_area = np.asarray(reflector_outer.area, dtype=float)
    q_rad_nodes = (
        0.2
        * sigma
        * radiation_area
        * (surface_temperature**4 - 200.0**4)
    )
    q_radiation = float(np.sum(q_rad_nodes))

    q_core = sum(
        float(tfe.neutronic_data.total_power) * float(core.tfe_multipliers[name])
        for name, tfe in build["tfes"].items()
    )
    q_balance_residual = q_core - q_fluid - q_electric - q_radiation

    return {
        "time_s": float(system.global_time),
        "solid_heat_capacity_scale": float(build["solid_heat_capacity_scale"]),
        "core_heat_power_w": q_core,
        "coolant_heat_pickup_w": q_fluid,
        "electric_power_w": q_electric,
        "outer_wall_radiation_w": q_radiation,
        "balance_residual_w": q_balance_residual,
        "balance_residual_percent": 100.0 * q_balance_residual / q_core,
        "mass_flow_kg_s": mass_flow,
        "inlet_plenum_temperature_k": float(build["inlet_plenum"].T),
        "outlet_plenum_temperature_k": float(build["outlet_plenum"].T),
        "outer_wall_area_m2": float(np.sum(radiation_area)),
        "outer_wall_radiation_area_avg_flux_w_m2": q_radiation / float(np.sum(radiation_area)),
        "tec_total_voltage_v": electric["tec_total_voltage_v"],
        "tec_total_current_a": electric["tec_total_current_a"],
    }


def run_v7_case_a_faststeady() -> Dict[str, Any]:
    result = run_test_v7_case_a_heated(
        run_duration_s=FASTSTEADY_RUN_DURATION_S,
        total_power_w=FASTSTEADY_TOTAL_POWER_W,
        max_dt=5.0,
        safety_factor=20.0,
        restart_file=FASTSTEADY_RESTART_IN,
        save_interval=0.0,
        final_restart_file=FASTSTEADY_RESTART_OUT,
        keep_only_latest_restart=False,
        reset_design_flow_after_restart=True,
        solid_heat_capacity_scale=FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE,
        solid_heat_capacity_scale_scope=FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE_SCOPE,
    )
    audit = compute_faststeady_energy_audit(result)
    result["energy_audit"] = audit
    print("FASTSTEADY_ENERGY_AUDIT")
    print(audit)
    return result


if __name__ == "__main__":
    run_v7_case_a_faststeady()
