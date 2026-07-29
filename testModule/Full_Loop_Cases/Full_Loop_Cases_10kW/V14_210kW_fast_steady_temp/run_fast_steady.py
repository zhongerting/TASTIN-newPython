"""Run an isolated 210 kW case with optionally reduced solid heat capacity."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Sequence

import numpy as np

REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Components.basicComponents.HeatPipe2D import HeatPipe2D  # noqa: E402
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (  # noqa: E402
    HISTORY_FIELDS,
    _apply_fixed_core_power,
    build_debug_case,
    collect_metrics,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (  # noqa: E402
    ReactivityControlRunConfig,
    load_baseline_debug_config,
)

CASE_NAME = "V14_10kW_210kW_fast_steady_temp"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "runs" / "default"
FAST_STEADY_HISTORY_FIELDS = HISTORY_FIELDS + ["solid_heat_capacity_scale"]


class ScaledHeatCapacityMaterial:
    """Delegate material properties while scaling only heat capacity."""

    def __init__(self, material: Any, scale: float):
        scale = float(scale)
        if not math.isfinite(scale) or scale <= 0.0:
            raise ValueError("heat-capacity scale must be finite and positive")
        self._material = material
        self.scale = scale

    def heat_capacity(self, temperature: Any) -> Any:
        return self.scale * self._material.heat_capacity(temperature)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._material, name)


@dataclass(frozen=True)
class FastSteadyRunConfig:
    restart_in: Path
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = 100.0
    dt_s: float = 0.05
    record_interval_s: float = 1.0
    checkpoint_interval_s: float = 10.0
    min_fluid_temperature_stop_k: Optional[float] = 500.0
    heat_capacity_scale: float = 0.01
    steady_window_s: float = 10.0
    steady_tolerance_k: float = 0.05
    inner_iter: Optional[int] = None
    interface_relaxation: float = 1.0


def scale_system_solids(system: Any, scale: float) -> list[str]:
    """Scale solid storage, leaving NaK and heat-pipe working-fluid Cp physical."""
    scale = float(scale)
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("heat-capacity scale must be finite and positive")
    if math.isclose(scale, 1.0, rel_tol=0.0, abs_tol=1.0e-15):
        return []

    scaled_names = []
    for name, solid in system.solid_components.items():
        if isinstance(solid, HeatPipe2D):
            solid.wall_mat = ScaledHeatCapacityMaterial(solid.wall_mat, scale)
            solid.material = solid.wall_mat
            solid.wick_mat.solid = ScaledHeatCapacityMaterial(
                solid.wick_mat.solid, scale,
            )
            solid._property_cache_initialized = False
        else:
            solid.material = ScaledHeatCapacityMaterial(solid.material, scale)
        solid._update_properties()
        scaled_names.append(str(name))
    return scaled_names


def _validate(config: FastSteadyRunConfig) -> None:
    for name in (
        "duration_s", "dt_s", "record_interval_s", "heat_capacity_scale",
        "steady_window_s", "steady_tolerance_k",
    ):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(float(config.checkpoint_interval_s)):
        raise ValueError("checkpoint_interval_s must be finite")
    if config.inner_iter is not None and int(config.inner_iter) < 1:
        raise ValueError("inner_iter must be positive")
    relaxation = float(config.interface_relaxation)
    if not math.isfinite(relaxation) or not 0.0 < relaxation <= 1.0:
        raise ValueError("interface_relaxation must be in (0, 1]")


def _temperatures(system: Any) -> Dict[str, np.ndarray]:
    values = {"fluid": np.asarray(system.fluid_solver.T_vec, dtype=float).copy()}
    values.update({
        f"solid:{name}": np.asarray(solid.T, dtype=float).copy()
        for name, solid in system.solid_components.items()
    })
    return values


def _maximum_change(
        start: Dict[str, np.ndarray], end: Dict[str, np.ndarray]) -> float:
    if start.keys() != end.keys():
        raise ValueError("temperature snapshots contain different components")
    changes = []
    for name, before in start.items():
        after = end[name]
        if before.shape != after.shape:
            raise ValueError(f"temperature shape changed for {name}")
        changes.append(float(np.max(np.abs(after - before))))
    return max(changes, default=0.0)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )


def _append_history(path: Path, row: Dict[str, Any]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=FAST_STEADY_HISTORY_FIELDS, extrasaction="ignore",
        )
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _print_progress(row: Dict[str, Any]) -> None:
    print(
        ("[t={time_s:.3f}s] Tin={core_inlet_T_K:.3f}K "
         "Tout={core_outlet_T_K:.3f}K Tsolid=[{min_solid_T_K:.3f}, "
         "{max_solid_T_K:.3f}]K Pnet={net_power_estimate_W:.2f}W "
         "TEC_ok={tec_main_converged}").format(**row),
        flush=True,
    )


def run_fast_steady(config: FastSteadyRunConfig) -> Dict[str, Any]:
    _validate(config)
    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in), output_dir=Path(config.output_dir),
        duration_s=float(config.duration_s), dt_s=float(config.dt_s),
        record_interval_s=float(config.record_interval_s),
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=config.min_fluid_temperature_stop_k,
    )
    debug, source = load_baseline_debug_config(runtime)
    if bool(source["point_kinetics_enabled"]):
        raise ValueError("fast-steady input must be a fixed-power restart")
    if not math.isclose(float(debug.power_w), 210000.0, rel_tol=0.0, abs_tol=1.0e-6):
        raise ValueError("fast-steady case requires exactly 210000 W")
    if bool(source.get("external_heat_enabled", False)):
        raise ValueError("fast-steady case requires external heat to be disabled")

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / "history.csv"
    if history_path.exists():
        raise FileExistsError(f"output history already exists: {history_path}")

    build = build_debug_case(debug, apply_fixed_power=True)
    system = build["system"]
    inner_iter = (
        int(debug.inner_iter)
        if config.inner_iter is None
        else int(config.inner_iter)
    )
    scaled_names = scale_system_solids(system, config.heat_capacity_scale)
    if scaled_names:
        system._prepare_fluid_sources_for_coupling()
        system._run_couplers(
            current_time=float(system.global_time), dt=float(config.dt_s),
        )
    _apply_fixed_core_power(build, debug.power_w)

    run_config = dict(source)
    run_config.update({
        "case": CASE_NAME,
        "duration_s": float(config.duration_s),
        "stage_durations_s": [float(config.duration_s)],
        "dt_s": float(config.dt_s),
        "record_interval_s": float(config.record_interval_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "min_fluid_temperature_stop_k": config.min_fluid_temperature_stop_k,
        "restart_in": str(config.restart_in),
        "source_run_config": str(Path(config.restart_in).parent / "run_config.json"),
        "solid_heat_capacity_scale": float(config.heat_capacity_scale),
        "scaled_solid_count": len(scaled_names),
        "scaled_solid_names": scaled_names,
        "steady_window_s": float(config.steady_window_s),
        "steady_tolerance_k": float(config.steady_tolerance_k),
        "inner_iter": inner_iter,
        "interface_relaxation": float(config.interface_relaxation),
        "external_heat_enabled": False,
        "point_kinetics_enabled": False,
    })
    _write_json(out_dir / "run_config.json", run_config)

    start_time = float(system.global_time)
    end_time = start_time + float(config.duration_s)
    window_target = end_time - min(
        float(config.steady_window_s), float(config.duration_s),
    )
    window_start_time = start_time
    window_start = None
    last_record_time = start_time
    last_checkpoint_time = start_time
    stop_reason = "completed"
    tec_nonconverged_steps = 0

    latest = collect_metrics(build, stage_index=1, dt_s=0.0)
    latest["solid_heat_capacity_scale"] = float(config.heat_capacity_scale)
    _append_history(history_path, latest)
    _print_progress(latest)

    while float(system.global_time) < end_time - 1.0e-9:
        current_time = float(system.global_time)
        if window_start is None and current_time >= window_target - 1.0e-9:
            window_start = _temperatures(system)
            window_start_time = current_time

        dt = min(float(config.dt_s), end_time - current_time)
        _apply_fixed_core_power(build, debug.power_w)
        system.step(
            dt, inner_iter=inner_iter,
            fail_on_fluid_nonconvergence=False,
            fluid_max_iter=int(debug.fluid_max_iter),
            interface_relaxation=float(config.interface_relaxation),
        )
        _apply_fixed_core_power(build, debug.power_w)

        latest = collect_metrics(build, stage_index=1, dt_s=dt)
        latest["solid_heat_capacity_scale"] = float(config.heat_capacity_scale)
        if not bool(latest["tec_main_converged"]):
            tec_nonconverged_steps += 1

        snapshot = _temperatures(system)
        if not all(np.all(np.isfinite(value)) for value in snapshot.values()):
            stop_reason = "nonfinite_temperature"
        else:
            threshold = config.min_fluid_temperature_stop_k
            if threshold is not None and float(threshold) > 0.0:
                if float(latest["min_fluid_T_K"]) < float(threshold):
                    stop_reason = "low_fluid_temperature"

        current_time = float(system.global_time)
        if (current_time - last_record_time
                >= float(config.record_interval_s) - 1.0e-9
                or current_time >= end_time - 1.0e-9
                or stop_reason != "completed"):
            _append_history(history_path, latest)
            _print_progress(latest)
            last_record_time = current_time

        if stop_reason != "completed":
            system.save_global_state(str(out_dir / "emergency_restart.npz"))
            break

        if (float(config.checkpoint_interval_s) > 0.0
                and current_time - last_checkpoint_time
                >= float(config.checkpoint_interval_s) - 1.0e-9):
            system.save_global_state(
                str(out_dir / f"checkpoint_t{current_time:.3f}s.npz")
            )
            last_checkpoint_time = current_time

    if window_start is None:
        window_start = _temperatures(system)
        window_start_time = float(system.global_time)
    window_end = _temperatures(system)
    maximum_change_k = _maximum_change(window_start, window_end)
    finite_final_state = all(
        np.all(np.isfinite(value)) for value in window_end.values()
    )
    steady_converged = (
        stop_reason == "completed" and finite_final_state
        and maximum_change_k <= float(config.steady_tolerance_k)
    )

    restart_path = out_dir / "stage_01_restart.npz"
    system.save_global_state(str(restart_path))
    result = {
        "case": CASE_NAME,
        "output_dir": str(out_dir),
        "history_path": str(history_path),
        "restart_path": str(restart_path),
        "source_restart_path": str(config.restart_in),
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "stop_reason": stop_reason,
        "solid_heat_capacity_scale": float(config.heat_capacity_scale),
        "scaled_solid_count": len(scaled_names),
        "steady_window_start_time_s": window_start_time,
        "steady_window_duration_s": float(system.global_time) - window_start_time,
        "max_temperature_change_over_window_k": maximum_change_k,
        "steady_tolerance_k": float(config.steady_tolerance_k),
        "steady_converged": steady_converged,
        "finite_final_state": finite_final_state,
        "tec_nonconverged_steps": tec_nonconverged_steps,
        "latest_metrics": latest,
    }
    summary_path = out_dir / "run_summary.json"
    _write_json(summary_path, result)
    _write_json(out_dir / "latest_state.json", {
        "case": CASE_NAME,
        "latest_restart_path": str(restart_path),
        "latest_summary_path": str(summary_path),
        "latest_metrics": latest,
    })
    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--duration", type=float, default=100.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--record-interval", type=float, default=1.0)
    parser.add_argument("--checkpoint-interval", type=float, default=10.0)
    parser.add_argument("--heat-capacity-scale", type=float, default=0.01)
    parser.add_argument("--steady-window", type=float, default=10.0)
    parser.add_argument("--steady-tolerance", type=float, default=0.05)
    parser.add_argument("--inner-iter", type=int, default=None)
    parser.add_argument("--interface-relaxation", type=float, default=1.0)
    parser.add_argument("--min-fluid-temperature-stop", type=float, default=500.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    threshold = (None if float(args.min_fluid_temperature_stop) <= 0.0
                 else float(args.min_fluid_temperature_stop))
    result = run_fast_steady(FastSteadyRunConfig(
        restart_in=args.restart_in,
        output_dir=args.output_dir,
        duration_s=float(args.duration),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        min_fluid_temperature_stop_k=threshold,
        heat_capacity_scale=float(args.heat_capacity_scale),
        steady_window_s=float(args.steady_window),
        steady_tolerance_k=float(args.steady_tolerance),
        inner_iter=args.inner_iter,
        interface_relaxation=float(args.interface_relaxation),
    ))
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
