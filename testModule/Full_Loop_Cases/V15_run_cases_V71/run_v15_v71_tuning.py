"""Minimal V71 tuning runner with explicit no-external-heat and flow calibration."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases import (  # noqa: E402
    FullLoopCoreConfig,
    FullLoopPumpConfig,
    V15PipeFinRadiatorConfig,
    V15_V71_RADIATOR_TUBE_K_LOSS,
    build_v15_v71_case_a_system,
)
from testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_cold_start import (  # noqa: E402
    _apply_core_power,
    collect_metrics,
)
from testModule.run_v8_caseA_common import WIRE_RESISTANCE_OHM  # noqa: E402


LOOKUP_DB = REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr"


def _set_ode_method(system: Any) -> None:
    for solid in system.solid_components.values():
        solid.set_ode_method("implicit_euler")


def _set_pump_head(build: Dict[str, Any], total_head_pa: float) -> None:
    single_head = 0.5 * float(total_head_pa)
    build["pump_a"].delta_p = single_head
    build["pump_b"].delta_p = single_head


def _flow_after_hydraulic_solve(build: Dict[str, Any], total_head_pa: float) -> float:
    _set_pump_head(build, total_head_pa)
    net = build["system"].fluid_solver
    net.initialize_hydraulics(dt=0.1, tol=1.0e-8, max_iter=1500, omega=0.5)
    return 0.5 * (float(build["pump_a"].W) + float(build["pump_b"].W))


def calibrate_pump_head(
    build: Dict[str, Any], target_flow_kg_s: float, nominal_head_pa: float
) -> Dict[str, float]:
    target = float(target_flow_kg_s)
    low = max(1000.0, float(nominal_head_pa) - 4000.0)
    high = float(nominal_head_pa) + 4000.0
    flow_low = _flow_after_hydraulic_solve(build, low)
    flow_high = _flow_after_hydraulic_solve(build, high)
    while flow_low > target:
        low *= 0.8
        flow_low = _flow_after_hydraulic_solve(build, low)
    while flow_high < target:
        high *= 1.2
        flow_high = _flow_after_hydraulic_solve(build, high)
    for _ in range(32):
        mid = 0.5 * (low + high)
        flow_mid = _flow_after_hydraulic_solve(build, mid)
        if flow_mid < target:
            low = mid
        else:
            high = mid
    head = 0.5 * (low + high)
    flow = _flow_after_hydraulic_solve(build, head)
    return {
        "pump_total_head_pa": float(head),
        "pump_single_head_pa": 0.5 * float(head),
        "total_flow_kg_s": float(flow),
        "flow_error_kg_s": float(flow - target),
        "hydraulic_max_dW_kg_s": float(np.max(build["system"].fluid_solver.W_residual)),
    }


def configure_tec(build: Dict[str, Any], voltage_v: float, wire_scale: float) -> None:
    core = build["core"]
    core.tec_lookup_enabled = True
    core.tec_lookup_db = str(LOOKUP_DB)
    core.tec_lookup_regions = ("core", "startup", "high_power", "accident")
    core.enable_tec_coupled = True
    core._build_thermo_calc()
    core.setup_tec_circuit("fixed_u", float(voltage_v), I_guess=150.0, topology="series")
    wire = np.asarray(WIRE_RESISTANCE_OHM, dtype=float) * float(wire_scale)
    core.thermo_calc._input_data.resistanceWire = np.tile(
        wire, (core.thermo_calc.N_elem, 1)
    )
    core.thermo_update_interval = 0.8


def _refresh_radiator_cache(build: Dict[str, Any]) -> None:
    system = build["system"]
    for unit in build["radiator_units"]:
        unit.pre_step(0.0, float(system.global_time))
    for solid in system.solid_components.values():
        solid._update_boundaries_state(current_time=float(system.global_time))


def _tube_flow_metrics(build: Dict[str, Any]) -> Dict[str, float]:
    junctions = {j.name: j for j in build["system"].fluid_solver.junctions_obj}
    flows = np.array(
        [junctions[f"J_RadiatorUpper_to_Tube_{i:02d}"].W for i in range(1, 79)],
        dtype=float,
    )
    return {
        "min_kg_s": float(np.min(flows)),
        "max_kg_s": float(np.max(flows)),
        "mean_kg_s": float(np.mean(flows)),
        "cv_percent": float(100.0 * np.std(flows) / np.mean(flows)),
        "range_percent_of_mean": float(
            100.0 * (np.max(flows) - np.min(flows)) / np.mean(flows)
        ),
    }


def _external_heat_input(build: Dict[str, Any]) -> float:
    return float(
        sum(
            np.sum(
                unit.get_external_heat_absorption_distribution(
                    build["system"].global_time
                )[2]
            )
            for unit in build["radiator_units"]
            if getattr(unit, "external_heat_bc", None) is not None
        )
    )


def final_metrics(build: Dict[str, Any], flow_calibration: Dict[str, float]) -> Dict[str, Any]:
    _refresh_radiator_cache(build)
    metrics = collect_metrics(build, 0.0)
    results = build["core"].get_tec_circuit_global_results().get("main") or {}
    flow = float(flow_calibration["total_flow_kg_s"])
    q_ext = _external_heat_input(build)
    q_rad = float(metrics["radiator_heat_rejection_W"])
    tec_power = float(results.get("Iout", 0.0)) * float(results.get("Uout", 0.0))
    h_in = float(build["core_inlet_connector"].h)
    h_out = float(build["core_outlet_connector"].h)
    return {
        **metrics,
        "pump": flow_calibration,
        "tube_flow": _tube_flow_metrics(build),
        "tec_voltage_V": float(results.get("Uout", np.nan)),
        "tec_current_A": float(results.get("Iout", np.nan)),
        "tec_power_W": tec_power,
        "tec_converged": bool(results.get("converged", False)),
        "tec_iteration_count": int(results.get("iteration_count", 0)),
        "external_heat_enabled": bool(build.get("external_heat_enabled", False)),
        "external_heat_input_W": q_ext,
        "fluid_core_enthalpy_gain_W": flow * (h_out - h_in),
        "energy_residual_W": float(
            build["core"].last_total_core_power + q_ext - tec_power - q_rad
        ),
    }


def run(args: argparse.Namespace) -> Dict[str, Any]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    build = build_v15_v71_case_a_system(
        core_config=FullLoopCoreConfig(
            inlet_temperature_k=723.0,
            main_tec_enabled=False,
        ),
        pump_config=FullLoopPumpConfig(pump_total_head_pa=float(args.pump_head)),
        radiator_config=V15PipeFinRadiatorConfig(
            t_space_k=4.0,
            tube_emissivity=float(args.emissivity),
            fin_emissivity=float(args.emissivity),
            external_heat_enabled=False,
            solid_ode_method="implicit_euler",
            radiator_tube_inlet_k_loss=V15_V71_RADIATOR_TUBE_K_LOSS,
            radiator_tube_outlet_k_loss=V15_V71_RADIATOR_TUBE_K_LOSS,
        ),
    )
    system = build["system"]
    system.load_global_state(str(args.restart))
    _set_ode_method(system)
    _apply_core_power(build, float(args.power))
    configure_tec(build, float(args.voltage), float(args.wire_scale))
    build["core"].post_step(0.0, float(system.global_time))
    build["core"].thermo_calc.calculate(verbose=False)

    flow_calibration = calibrate_pump_head(
        build, float(args.target_flow), float(args.pump_head)
    )
    history = []
    record_interval = max(float(args.record_interval), float(args.dt))
    next_record = float(system.global_time)
    end_time = float(system.global_time) + float(args.duration)
    while system.global_time < end_time - 1.0e-12:
        dt = min(float(args.dt), end_time - float(system.global_time))
        _apply_core_power(build, float(args.power))
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=100)
        _apply_core_power(build, float(args.power))
        if system.global_time + 1.0e-12 >= next_record or system.global_time >= end_time - 1.0e-12:
            row = final_metrics(build, flow_calibration)
            row["time_s"] = float(system.global_time)
            history.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)
            next_record = system.global_time + record_interval

    flow_calibration = calibrate_pump_head(
        build, float(args.target_flow), float(flow_calibration["pump_total_head_pa"])
    )
    result = {
        "case": "V15_V71_tuning_no_external_heat",
        "restart": str(args.restart),
        "absolute_start_time_s": float(system.global_time - float(args.duration)),
        "absolute_end_time_s": float(system.global_time),
        "duration_s": float(args.duration),
        "dt_s": float(args.dt),
        "core_power_W": float(args.power),
        "radiator_emissivity": float(args.emissivity),
        "tec_voltage_target_V": float(args.voltage),
        "wire_resistance_scale": float(args.wire_scale),
        "external_heat_enabled": False,
        "final": final_metrics(build, flow_calibration),
        "history": history,
    }
    (output_dir / "tuning_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=0.0)
    parser.add_argument("--dt", type=float, default=0.1)
    parser.add_argument("--record-interval", type=float, default=20.0)
    parser.add_argument("--power", type=float, default=106000.0)
    parser.add_argument("--emissivity", type=float, default=0.8)
    parser.add_argument("--voltage", type=float, default=28.0)
    parser.add_argument("--wire-scale", type=float, default=0.0)
    parser.add_argument("--pump-head", type=float, default=20200.0)
    parser.add_argument("--target-flow", type=float, default=1.18)
    return parser


if __name__ == "__main__":
    started = time.perf_counter()
    result = run(_parser().parse_args())
    print(json.dumps({"final": result["final"], "elapsed_s": time.perf_counter() - started}, indent=2, sort_keys=True))
