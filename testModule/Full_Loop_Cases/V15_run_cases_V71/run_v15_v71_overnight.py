"""Long V71 continuation with one rolling restart file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopPumpConfig,
    V15PipeFinRadiatorConfig,
    V15_V71_RADIATOR_TUBE_K_LOSS,
)
from testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_tuning import (
    _apply_core_power,
    _set_ode_method,
    build_v15_v71_case_a_system,
    calibrate_pump_head,
    configure_tec,
    final_metrics,
)


def _save_latest_restart(system, output_dir: Path) -> Path:
    restart_path = output_dir / "v15_v71_final_restart.npz"
    temp_path = output_dir / "v15_v71_final_restart.tmp.npz"
    temp_path.unlink(missing_ok=True)
    system.save_global_state(str(temp_path))
    temp_path.replace(restart_path)
    return restart_path


def run(args: argparse.Namespace) -> dict:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_interval = float(args.checkpoint_interval)
    if checkpoint_interval <= 0.0:
        raise ValueError("--checkpoint-interval must be positive")

    build = build_v15_v71_case_a_system(
        core_config=FullLoopCoreConfig(inlet_temperature_k=723.0, main_tec_enabled=False),
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

    start_time = float(system.global_time)
    flow_calibration = calibrate_pump_head(build, float(args.target_flow), float(args.pump_head))
    restart_out = _save_latest_restart(system, output_dir)
    history = []
    next_record = start_time
    next_checkpoint = start_time + checkpoint_interval
    end_time = start_time + float(args.duration)

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

        if system.global_time + 1.0e-12 >= next_checkpoint:
            restart_out = _save_latest_restart(system, output_dir)
            while next_checkpoint <= system.global_time + 1.0e-12:
                next_checkpoint += checkpoint_interval

    flow_calibration = calibrate_pump_head(
        build, float(args.target_flow), float(flow_calibration["pump_total_head_pa"])
    )
    final = final_metrics(build, flow_calibration)
    final["time_s"] = float(system.global_time)
    restart_out = _save_latest_restart(system, output_dir)
    result = {
        "case": "V15_V71_near_steady_no_external_heat",
        "restart_in": str(args.restart),
        "restart_out": str(restart_out),
        "absolute_start_time_s": start_time,
        "absolute_end_time_s": float(system.global_time),
        "duration_s": float(args.duration),
        "dt_s": float(args.dt),
        "checkpoint_interval_s": checkpoint_interval,
        "core_power_W": float(args.power),
        "radiator_emissivity": float(args.emissivity),
        "tec_voltage_target_V": float(args.voltage),
        "wire_resistance_scale": float(args.wire_scale),
        "external_heat_enabled": False,
        "final": final,
        "history": history,
    }
    (output_dir / "steady_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=1000.0)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--record-interval", type=float, default=100.0)
    parser.add_argument("--checkpoint-interval", type=float, default=10.0)
    parser.add_argument("--power", type=float, default=98000.0)
    parser.add_argument("--emissivity", type=float, default=0.83)
    parser.add_argument("--voltage", type=float, default=28.0)
    parser.add_argument("--wire-scale", type=float, default=0.02)
    parser.add_argument("--pump-head", type=float, default=20200.0)
    parser.add_argument("--target-flow", type=float, default=1.18)
    return parser


if __name__ == "__main__":
    result = run(_parser().parse_args())
    print(json.dumps({"final": result["final"]}, indent=2, sort_keys=True))
