"""V71 tuning runner with independent tube and fin emissivities."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import numpy as np

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopPumpConfig,
    V15PipeFinRadiatorConfig,
    V15_V71_RADIATOR_TUBE_K_LOSS,
    build_v15_v71_case_a_system,
)
from testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_tuning import (
    _apply_core_power,
    _external_heat_input,
    _refresh_radiator_cache,
    _set_ode_method,
    _tube_flow_metrics,
    calibrate_pump_head,
    collect_metrics,
    configure_tec,
)


def final_metrics(build: Dict[str, Any], flow_calibration: Dict[str, float]) -> Dict[str, Any]:
    _refresh_radiator_cache(build)
    metrics = collect_metrics(build, 0.0)
    results = build["core"].get_tec_circuit_global_results().get("main") or {}
    q_ext = _external_heat_input(build)
    q_rad = float(metrics["radiator_heat_rejection_W"])
    tec_power = float(results.get("Iout", 0.0)) * float(results.get("Uout", 0.0))
    flow = float(flow_calibration["total_flow_kg_s"])
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
        "external_heat_enabled": False,
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
        core_config=FullLoopCoreConfig(inlet_temperature_k=723.0, main_tec_enabled=False),
        pump_config=FullLoopPumpConfig(pump_total_head_pa=float(args.pump_head)),
        radiator_config=V15PipeFinRadiatorConfig(
            t_space_k=4.0,
            tube_emissivity=float(args.tube_emissivity),
            fin_emissivity=float(args.fin_emissivity),
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

    flow_calibration = calibrate_pump_head(build, float(args.target_flow), float(args.pump_head))
    history = []
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
            next_record = system.global_time + float(args.record_interval)

    flow_calibration = calibrate_pump_head(
        build, float(args.target_flow), float(flow_calibration["pump_total_head_pa"])
    )
    result = {
        "case": "V15_V71_split_emissivity_no_external_heat",
        "restart": str(args.restart),
        "absolute_start_time_s": float(system.global_time - float(args.duration)),
        "absolute_end_time_s": float(system.global_time),
        "duration_s": float(args.duration),
        "dt_s": float(args.dt),
        "core_power_W": float(args.power),
        "tube_emissivity": float(args.tube_emissivity),
        "fin_emissivity": float(args.fin_emissivity),
        "tec_voltage_target_V": float(args.voltage),
        "wire_resistance_scale": float(args.wire_scale),
        "external_heat_enabled": False,
        "final": final_metrics(build, flow_calibration),
        "history": history,
    }
    (output_dir / "split_tuning_summary.json").write_text(
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
    parser.add_argument("--power", type=float, default=98000.0)
    parser.add_argument("--tube-emissivity", type=float, default=0.85)
    parser.add_argument("--fin-emissivity", type=float, default=0.80)
    parser.add_argument("--voltage", type=float, default=28.0)
    parser.add_argument("--wire-scale", type=float, default=0.0)
    parser.add_argument("--pump-head", type=float, default=20200.0)
    parser.add_argument("--target-flow", type=float, default=1.18)
    return parser


if __name__ == "__main__":
    result = run(_parser().parse_args())
    print(json.dumps({"final": result["final"]}, indent=2, sort_keys=True))
