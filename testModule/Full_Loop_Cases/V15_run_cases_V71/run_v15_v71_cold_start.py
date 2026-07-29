"""Cold-start V15 V71 thermal run and one-shot TEC lookup solve."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases import (  # noqa: E402
    FullLoopCoreConfig,
    FullLoopPumpConfig,
    V15_V71_PUMP_TOTAL_HEAD_PA,
    V15_V71_RADIATOR_TUBE_K_LOSS,
    V15PipeFinRadiatorConfig,
    build_v15_v71_case_a_system,
)

DEFAULT_LOOKUP_DB = REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "runs" / "cold_start_723k_1000s"


@dataclass(frozen=True)
class ColdStartConfig:
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = 1000.0
    dt_s: float = 0.5
    record_interval_s: float = 10.0
    checkpoint_interval_s: float = 100.0
    initial_temperature_k: float = 723.0
    space_temperature_k: float = 4.0
    core_power_w: float = 106000.0
    tec_voltage_v: float = 27.2
    tec_current_guess_a: float = 150.0
    lookup_db: Path = DEFAULT_LOOKUP_DB
    lookup_regions: Sequence[str] = ("core", "startup", "high_power", "accident")


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def _sync_fluid_temperature(net: Any, temperature_k: float) -> None:
    for volume in net.volumes_obj:
        volume.T = float(temperature_k)
        volume.h = float(volume.material.enthalpy(volume.T, volume.P))
        volume.update_properties(volume.material)
    net._initialize_state_from_objects()
    net._update_fluid_properties()
    net._sync_vectors_to_objects()


def _sync_solid_temperature(system: Any, temperature_k: float) -> None:
    for solid in system.solid_components.values():
        solid.T[:] = float(temperature_k)
        solid.current_time = float(system.global_time)
        solid._update_properties()
        solid._update_boundaries_state(current_time=float(system.global_time))


def _set_implicit_euler(system: Any) -> None:
    for solid in system.solid_components.values():
        solid.set_ode_method("implicit_euler")


def _apply_core_power(build: Dict[str, Any], power_w: float) -> None:
    build["core"].update_neutronic_power(
        p_total=float(power_w), p_fiss=float(power_w), p_decay=0.0, alpha=1.0
    )


def build_cold_start_case(
    config: ColdStartConfig,
    *,
    initialize_hydraulics: bool = True,
) -> Dict[str, Any]:
    build = build_v15_v71_case_a_system(
        pump_config=FullLoopPumpConfig(pump_total_head_pa=V15_V71_PUMP_TOTAL_HEAD_PA),
        core_config=FullLoopCoreConfig(
            inlet_temperature_k=float(config.initial_temperature_k),
            main_tec_enabled=False,
        ),
        radiator_config=V15PipeFinRadiatorConfig(
            t_space_k=float(config.space_temperature_k),
            solid_ode_method="implicit_euler",
            radiator_tube_inlet_k_loss=V15_V71_RADIATOR_TUBE_K_LOSS,
            radiator_tube_outlet_k_loss=V15_V71_RADIATOR_TUBE_K_LOSS,
        ),
    )
    system = build["system"]
    _sync_fluid_temperature(system.fluid_solver, config.initial_temperature_k)
    _sync_solid_temperature(system, config.initial_temperature_k)
    _set_implicit_euler(system)
    _apply_core_power(build, config.core_power_w)
    if initialize_hydraulics:
        system.initialize_system(dt_init=0.01, tol=1.0e-4, max_iter=1000)
        _apply_core_power(build, config.core_power_w)
    return build


def collect_metrics(build: Dict[str, Any], dt_s: float) -> Dict[str, Any]:
    system = build["system"]
    net = system.fluid_solver
    solids = np.concatenate([
        np.asarray(solid.T, dtype=float).reshape(-1)
        for solid in system.solid_components.values()
    ])
    radiator_q = sum(
        float(np.sum(unit.get_heat_exchange_breakdown()["gross_rejection"]))
        for unit in build["radiator_units"]
    )
    diag = system.last_step_diagnostics or {}
    fluid_flags = diag.get("fluid_converged_by_iteration", [])
    return {
        "time_s": float(system.global_time),
        "dt_s": float(dt_s),
        "core_power_W": float(build["core"].last_total_core_power),
        "core_inlet_T_K": float(build["core_inlet_connector"].T),
        "core_outlet_T_K": float(build["core_outlet_connector"].T),
        "core_delta_T_K": float(build["core_outlet_connector"].T - build["core_inlet_connector"].T),
        "radiator_heat_rejection_W": radiator_q,
        "min_fluid_T_K": float(np.min(net.T_vec)),
        "max_fluid_T_K": float(np.max(net.T_vec)),
        "min_solid_T_K": float(np.min(solids)),
        "max_solid_T_K": float(np.max(solids)),
        "mean_emitter_T_K": float(np.mean([
            np.mean(tfe.solids["emitter"].T) for tfe in build["tfes"].values()
        ])),
        "mean_collector_T_K": float(np.mean([
            np.mean(tfe.solids["collector"].T) for tfe in build["tfes"].values()
        ])),
        "fluid_converged": bool(all(fluid_flags)) if fluid_flags else True,
    }


def run_thermal(config: ColdStartConfig) -> Dict[str, Any]:
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    build = build_cold_start_case(config)
    system = build["system"]
    history_path = out_dir / "thermal_history.csv"
    fields = list(collect_metrics(build, 0.0))
    with history_path.open("w", newline="", encoding="utf-8") as handle:
        csv.DictWriter(handle, fieldnames=fields).writeheader()

    next_record = 0.0
    next_checkpoint = float(config.checkpoint_interval_s)
    final_time = float(config.duration_s)
    latest: Dict[str, Any] = {}
    while system.global_time < final_time - 1.0e-12:
        dt = min(float(config.dt_s), final_time - float(system.global_time))
        _apply_core_power(build, config.core_power_w)
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=100)
        _apply_core_power(build, config.core_power_w)

        if system.global_time + 1.0e-12 >= next_record or system.global_time >= final_time - 1.0e-12:
            latest = collect_metrics(build, dt)
            if not all(math.isfinite(float(value)) for key, value in latest.items() if key != "fluid_converged"):
                raise FloatingPointError(f"Non-finite state at t={system.global_time}: {latest}")
            with history_path.open("a", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=fields).writerow(latest)
            print(json.dumps(latest, sort_keys=True), flush=True)
            next_record = system.global_time + float(config.record_interval_s)

        if config.checkpoint_interval_s > 0.0 and system.global_time + 1.0e-12 >= next_checkpoint:
            system.save_global_state(str(out_dir / f"checkpoint_t{int(round(system.global_time)):04d}s.npz"))
            next_checkpoint += float(config.checkpoint_interval_s)

    restart_path = out_dir / "thermal_1000s_restart.npz"
    system.save_global_state(str(restart_path))
    summary = {
        "case": "V15_V71_cold_start_723K",
        "config": {**asdict(config), "output_dir": str(config.output_dir), "lookup_db": str(config.lookup_db)},
        "restart_path": str(restart_path),
        "latest_metrics": latest,
    }
    _write_json(out_dir / "thermal_summary.json", summary)
    return summary


def run_tec_once(config: ColdStartConfig, restart_path: Path) -> Dict[str, Any]:
    build = build_cold_start_case(config)
    system = build["system"]
    system.load_global_state(str(restart_path))
    core = build["core"]
    core.tec_lookup_enabled = True
    core.tec_lookup_db = str(config.lookup_db)
    core.tec_lookup_regions = tuple(config.lookup_regions)
    core.enable_tec_coupled = True
    core._build_thermo_calc()
    core.setup_tec_circuit("fixed_u", config.tec_voltage_v, I_guess=config.tec_current_guess_a, topology="series")
    core.post_step(0.0, float(system.global_time))
    elapsed_ms = float(core.thermo_calc.calculate(verbose=True))
    results = core.get_tec_circuit_global_results()["main"] or {}
    current = float(results.get("Iout", float("nan")))
    voltage = float(results.get("Uout", float("nan")))
    output = {
        "time_s": float(system.global_time),
        "lookup_enabled": bool(core.thermo_calc.lookup_enabled),
        "lookup_db": str(config.lookup_db),
        "mode": results.get("mode"),
        "target_voltage_V": float(config.tec_voltage_v),
        "output_voltage_V": voltage,
        "output_current_A": current,
        "electric_power_W": voltage * current,
        "converged": bool(results.get("converged", False)),
        "iteration_count": int(results.get("iteration_count", 0)),
        "zero_emission_skipped": bool(results.get("zero_emission_skipped", False)),
        "zero_emission_reason": results.get("zero_emission_reason"),
        "elapsed_ms": elapsed_ms,
        "thermal_metrics": collect_metrics(build, 0.0),
    }
    _write_json(Path(config.output_dir) / "tec_once_27p2V.json", output)
    print(json.dumps(output, indent=2, sort_keys=True), flush=True)
    return output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("thermal", "tec-once"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--duration", type=float, default=1000.0)
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--checkpoint-interval", type=float, default=100.0)
    parser.add_argument("--restart", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    config = ColdStartConfig(
        output_dir=args.output_dir,
        duration_s=float(args.duration),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
    )
    if args.mode == "thermal":
        run_thermal(config)
    else:
        if args.restart is None:
            raise SystemExit("--restart is required for tec-once")
        run_tec_once(config, args.restart)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
