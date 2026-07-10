from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path
from typing import Any
import sys

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
DEFAULT_LOOKUP_DB = REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr"
DEFAULT_LOOKUP_REGIONS = ("core", "startup", "high_power", "accident")
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases_VV.Luchau_V.luchau_single_tfe_model import (
    LuchauSingleTFEConfig,
    build_luchau_single_tfe,
    configure_luchau_thermocalc,
    tcs_from_cesium_pressure,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import (
    _advance_system_until_steady,
    json_default,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Preheat Luchau_V single TFE thermally, then run ThermoCalc lookup-backed fixed_u voltage sweep.")
    parser.add_argument("--thermal-power-w", type=float, required=True)
    parser.add_argument("--heater-length-m", type=float, default=0.30)
    parser.add_argument("--voltage-start-v", type=float, default=0.30)
    parser.add_argument("--voltage-end-v", type=float, default=0.40)
    parser.add_argument("--voltage-step-v", type=float, default=0.02)
    parser.add_argument("--duration-s", type=float, default=200.0, help="Maximum thermal preheat duration.")
    parser.add_argument("--dt-s", type=float, default=0.02)
    parser.add_argument("--steady-dtemp-k-s", type=float, default=1.0e-3)
    parser.add_argument("--steady-window-steps", type=int, default=20)
    parser.add_argument("--lookup-db", type=Path, default=DEFAULT_LOOKUP_DB)
    parser.add_argument("--lookup-regions", type=str, default=",".join(DEFAULT_LOOKUP_REGIONS))
    parser.add_argument("--output-label", type=str, default=None)
    parser.add_argument("--restart-in", type=Path, default=None)
    parser.add_argument("--restart-out", type=Path, default=None)
    return parser


def _voltage_points(start_v: float, end_v: float, step_v: float) -> list[float]:
    if step_v <= 0.0:
        raise ValueError("step_v must be positive.")
    if end_v < start_v:
        raise ValueError("end_v must be greater than or equal to start_v.")
    values: list[float] = []
    index = 0
    while True:
        value = round(float(start_v) + index * float(step_v), 10)
        if value > float(end_v) + 1.0e-10:
            break
        values.append(round(value, 2))
        index += 1
    if not values or abs(values[-1] - float(end_v)) > 1.0e-9:
        values.append(round(float(end_v), 2))
    return values


def _parse_lookup_regions(value: str | tuple[str, ...] | list[str]) -> tuple[str, ...]:
    if isinstance(value, str):
        regions = tuple(part.strip() for part in value.split(",") if part.strip())
    else:
        regions = tuple(str(part).strip() for part in value if str(part).strip())
    return regions or DEFAULT_LOOKUP_REGIONS


def _ensure_lookup_loaded(lookup_db: Path, lookup_regions: str | tuple[str, ...] | list[str]) -> dict[str, Any]:
    from ThermoCalc import ThermoCalcWrapper as tcw

    db_path = lookup_db.resolve()
    regions = _parse_lookup_regions(lookup_regions)
    loaded_count = tcw.load_emission_lookup_database(
        str(db_path),
        enable=True,
        force=False,
        regions=regions,
    )
    enabled = bool(tcw.te_solver.is_emission_lookup_enabled()) if hasattr(tcw.te_solver, "is_emission_lookup_enabled") else None
    dense_count = int(tcw.te_solver.emission_lookup_dense_region_count()) if hasattr(tcw.te_solver, "emission_lookup_dense_region_count") else None
    block_count = int(tcw.te_solver.emission_lookup_block_count()) if hasattr(tcw.te_solver, "emission_lookup_block_count") else None
    return {
        "enabled": enabled,
        "db_dir": str(db_path),
        "regions": list(regions),
        "loaded_count": int(loaded_count),
        "dense_region_count": dense_count,
        "legacy_block_count": block_count,
    }


def run_case(args: argparse.Namespace) -> dict[str, Any]:
    voltages = _voltage_points(float(args.voltage_start_v), float(args.voltage_end_v), float(args.voltage_step_v))
    config = LuchauSingleTFEConfig(
        thermal_power_w=float(args.thermal_power_w),
        target_voltage_v=float(voltages[0]),
        heater_length_m=float(args.heater_length_m),
    )
    build = build_luchau_single_tfe(config)
    if args.restart_in is not None:
        build["system"].load_global_state(str(args.restart_in))
    preheat_result = _advance_system_until_steady(
        build["system"],
        build,
        float(args.duration_s),
        float(args.dt_s),
        float(args.steady_dtemp_k_s),
        int(args.steady_window_steps),
    )
    if args.restart_out is not None:
        args.restart_out.parent.mkdir(parents=True, exist_ok=True)
        build["system"].save_global_state(str(args.restart_out))
    lookup_info = _ensure_lookup_loaded(Path(args.lookup_db), args.lookup_regions)
    sweep_results = _run_sweep(build, config, voltages)
    summary = _build_summary(config, build, voltages, preheat_result, sweep_results, lookup_info)
    summary["restart_in"] = None if args.restart_in is None else str(args.restart_in)
    summary["restart_out"] = None if args.restart_out is None else str(args.restart_out)
    output_dir = _write_summary(summary, args.output_label)
    summary["output_dir"] = str(output_dir)
    return summary


def _run_sweep(build: dict[str, Any], base_config: LuchauSingleTFEConfig, voltages: list[float]) -> list[dict[str, Any]]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    results: list[dict[str, Any]] = []
    for voltage in voltages:
        config = replace(base_config, target_voltage_v=float(voltage))
        thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
        configure_luchau_thermocalc(thermo_model, build, config)
        thermo_model.calculate(verbose=False)
        global_results = thermo_model.get_global_results()
        tec_results = thermo_model.get_tec_results(0)
        current = _finite_float(global_results.get("Iout")) if global_results else None
        voltage_out = _finite_float(global_results.get("Uout")) if global_results else None
        results.append({
            "target_voltage_v": float(voltage),
            "Iout_a": current,
            "Uout_v": voltage_out,
            "electric_power_w": None if current is None or voltage_out is None else float(current * voltage_out),
            "converged": None if not global_results else bool(global_results.get("converged")),
            "iteration_count": None if not global_results else global_results.get("iteration_count"),
            "zero_emission_skipped": None if not global_results else bool(global_results.get("zero_emission_skipped", False)),
            "zero_emission_reason": None if not global_results else global_results.get("zero_emission_reason"),
            "thermocalc_global_results": global_results,
            "tec_results_present": tec_results is not None,
        })
    return results


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _build_summary(
    config: LuchauSingleTFEConfig,
    build: dict[str, Any],
    voltages: list[float],
    preheat_result: dict[str, Any],
    sweep_results: list[dict[str, Any]],
    lookup_info: dict[str, Any],
) -> dict[str, Any]:
    tfe = build["tfe"]
    emitter_t = np.asarray(tfe.solids["emitter"].T, dtype=float)
    collector_t = np.asarray(tfe.solids["collector"].T, dtype=float)
    coolant_t = np.array([float(vol.T) for vol in build["channel"].volumes], dtype=float)
    return {
        "case": "Luchau_V_single_TFE_preheat_then_lookup_fixed_u_sweep",
        "config": asdict(config),
        "voltage_points_v": voltages,
        "preheat_mode": "thermal_only_no_thermocalc",
        "preheat_result": preheat_result,
        "lookup": lookup_info,
        "full_length_m": float(tfe.geom.height),
        "axial_node_count": int(tfe.mesh.n_axial),
        "pellet_power_w": float(np.sum(tfe.solids["pellet"].Q_source)),
        "coolant_inlet_temperature_k": float(config.inlet_temperature_k),
        "coolant_outlet_temperature_k": float(coolant_t[-1]),
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(config.cesium_pressure_torr)),
        "emitter_temperature_mean_k": float(np.mean(emitter_t)),
        "emitter_temperature_max_k": float(np.max(emitter_t)),
        "collector_temperature_mean_k": float(np.mean(collector_t)),
        "collector_temperature_max_k": float(np.max(collector_t)),
        "sweep_results": sweep_results,
        "final_time_s": float(build["system"].global_time),
        "finite": bool(np.isfinite([np.mean(emitter_t), np.max(emitter_t), np.mean(collector_t), coolant_t[-1]]).all()),
    }


def _write_summary(summary: dict[str, Any], label: str | None) -> Path:
    if label is None:
        label = datetime.now().strftime("%Y%m%d_%H%M%S_luchau_preheat_then_sweep")
    output_dir = CASE_DIR / "runs" / label
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return output_dir


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = run_case(args)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
