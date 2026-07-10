from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Any
import sys

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases_VV.Luchau_V.luchau_single_tfe_model import (
    LuchauSingleTFEConfig,
    LuchauSingleTFEThermoCalcCoupler,
    build_luchau_single_tfe,
    configure_luchau_thermocalc,
    tcs_from_cesium_pressure,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_preheat_then_sweep import (
    DEFAULT_LOOKUP_DB,
    DEFAULT_LOOKUP_REGIONS,
    _ensure_lookup_loaded,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import json_default


TARGET_CSV = CASE_DIR / "reference_image_all_power_iv_fit_range_250_410A.csv"
DEFAULT_3P26_VOLTAGE_LIST = "1.49793,1.48345,1.47103,1.45862,1.44414,1.43172,1.41517,1.39655,1.37586,1.35517,1.33448,1.31379,1.29517,1.27448,1.25793,1.23931,1.21862,1.19793,1.17724,1.15655,1.14207,1.12345,1.10483,1.08621,1.06759,1.0469,1.02828,1.00759,0.98897,0.96828,0.94966,0.93103,0.91034,0.88966,0.87103,0.85034,0.82966"
DEFAULT_3P26_CURRENT_LIST = "200.23256,208.60465,216.97674,225.34884,232.7907,241.16279,250.46512,257.90698,265.34884,272.7907,281.16279,288.60465,296.04651,303.48837,309.06977,315.5814,323.02326,330.46512,337.90698,345.34884,350.93023,357.44186,364.88372,371.39535,378.83721,387.2093,394.65116,402.09302,409.53488,417.90698,424.4186,431.86047,440.23256,447.67442,455.11628,464.4186,472.7907"
DEFAULT_PROGRESS_INTERVAL_S = 10.0

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fit one thermal power against manually provided voltage points, or fall back to a target CSV.")
    # 单次计算的热输入功率，单位 kWt；默认 3.26 表示总热功率 3.26 kW。
    parser.add_argument("--thermal-power-kwt", type=float, default=5.6)
    # 输出电压点列表，单位 V；多个值用英文逗号分隔。默认是 3.26 kWt 曲线在 250-410 A 范围内的电压点。
    parser.add_argument("--voltage-list", type=str, default=DEFAULT_3P26_VOLTAGE_LIST)
    # 可选目标电流列表，单位 A；用于计算误差，数量必须与 voltage-list 相同。没有目标电流时可设为空字符串。
    parser.add_argument("--target-current-list", type=str, default=DEFAULT_3P26_CURRENT_LIST)
    # 铯蒸汽压力扫描列表，单位 torr；多个值用英文逗号分隔，例如 "1.2,1.35,1.5"。
    parser.add_argument("--pressure-list", type=str, default="1.5,2.1,2.7")
    # 导线电阻扫描列表，单位 ohm；当前脚本会把每个标量电阻同时用于 ThermoCalc 的 4 段导线电阻。
    parser.add_argument("--wire-list", type=str, default="0.00008")
    # 备用参考 I-V 数据文件；只有当 voltage-list 为空字符串时才从 CSV 读取多个工况点。
    parser.add_argument("--target-csv", type=Path, default=TARGET_CSV)
    # 每个工况点的瞬态计算物理时长，单位 s；例如 200 表示每个电压点从预热 restart 继续算 200 s。
    parser.add_argument("--duration-s", type=float, default=200.0)
    # 单步时间步长，单位 s；值越小越稳但越慢。
    parser.add_argument("--dt-s", type=float, default=0.1)
    # 单进程详细进度输出的物理时间间隔，单位 s；仅 workers=1 时会输出每步进度。
    parser.add_argument("--progress-interval-s", type=float, default=DEFAULT_PROGRESS_INTERVAL_S)
    # 并行进程数；会同时计算不同电压点/参数组合。直接断点调试建议设为 1，批量扫描可设为 6 或 8。
    parser.add_argument("--workers", type=int, default=16)
    # 本次计算结果输出目录名，会写到 Luchau_V/runs/<output-label>/。
    parser.add_argument("--output-label", type=str, default="manual_debug_3p26_best_200s")
    return parser


def _parse_float_list(value: Any, *, positive: bool = False, nonnegative: bool = False, allow_empty: bool = False) -> list[float]:
    if value is None:
        values: list[float] = []
    elif isinstance(value, (int, float, np.floating)):
        values = [float(value)]
    elif isinstance(value, (list, tuple)):
        values = [float(item) for item in value]
    else:
        values = [float(part.strip()) for part in str(value).split(",") if part.strip()]
    if not values:
        if allow_empty:
            return []
        raise ValueError("list argument must contain at least one value.")
    for item in values:
        if positive and item <= 0.0:
            raise ValueError("values must be positive.")
        if nonnegative and item < 0.0:
            raise ValueError("values must be non-negative.")
    return values


def _load_targets(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for index, row in enumerate(reader, start=1):
            rows.append({
                "row_index": index,
                "power_kwt": float(row["power_kwt"]),
                "point_index": int(row["point_index"]),
                "voltage_v": float(row["voltage_v"]),
                "target_current_a": float(row["current_a"]),
            })
    if not rows:
        raise ValueError(f"no target rows found in {path}")
    return rows


def _manual_targets_from_args(
    *,
    power_kwt: float,
    voltage_list: list[float],
    target_current_list: list[float] | None = None,
) -> list[dict[str, Any]]:
    if not voltage_list:
        raise ValueError("voltage_list must contain at least one voltage point.")
    if target_current_list is not None and len(target_current_list) != len(voltage_list):
        raise ValueError("target_current_list length must match voltage_list length.")
    currents = target_current_list if target_current_list is not None else [0.0] * len(voltage_list)
    rows: list[dict[str, Any]] = []
    for index, (voltage, current) in enumerate(zip(voltage_list, currents), start=1):
        rows.append({
            "row_index": index,
            "power_kwt": float(power_kwt),
            "point_index": index,
            "voltage_v": float(voltage),
            "target_current_a": float(current),
        })
    return rows


def _progress_interval_steps(duration_s: float, dt_s: float, progress_interval_s: float) -> int:
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    total_steps = max(1, int(round(float(duration_s) / float(dt_s))))
    if progress_interval_s <= 0.0:
        return total_steps
    return max(1, min(total_steps, int(round(float(progress_interval_s) / float(dt_s)))))


def _format_optional_float(value: Any, fmt: str = ".3g") -> str:
    number = _finite_float(value)
    if number is None:
        return "nan"
    return format(number, fmt)

def _power_label(power_kw: float) -> str:
    return str(power_kw).replace(".", "p")


def _restart_token(value: float) -> str:
    text = f"{float(value):.8g}"
    return text.replace("-", "m").replace("+", "").replace(".", "p")


def _case_restart_path(output_dir: Path, task: dict[str, Any]) -> Path:
    return output_dir / "case_restarts" / (
        f"P{_restart_token(float(task['power_kwt']))}kWt_"
        f"U{_restart_token(float(task['voltage_v']))}V_"
        f"Cs{_restart_token(float(task['pressure_torr']))}torr_"
        f"Rw{_restart_token(float(task['wire_resistance_ohm']))}ohm.npz"
    )


def _thermal_restart_path(power_kw: float) -> Path:
    return CASE_DIR / "runs" / f"thermal_center0p30_uniform_{_power_label(power_kw)}kw_1000s_dt0p05" / "preheat_restart.npz"


def _temperature_vector(build: dict[str, Any]) -> np.ndarray:
    return np.concatenate(
        [np.asarray(solid.T, dtype=float).ravel() for solid in build["tfe"].solids.values() if hasattr(solid, "T")]
    )


def _temperature_stats(build: dict[str, Any]) -> dict[str, float]:
    tfe = build["tfe"]
    emitter = np.asarray(tfe.solids["emitter"].T, dtype=float)
    collector = np.asarray(tfe.solids["collector"].T, dtype=float)
    coolant = np.array([float(vol.T) for vol in build["channel"].volumes], dtype=float)
    return {
        "emitter_temperature_mean_k": float(np.mean(emitter)),
        "emitter_temperature_max_k": float(np.max(emitter)),
        "collector_temperature_mean_k": float(np.mean(collector)),
        "collector_temperature_max_k": float(np.max(collector)),
        "coolant_outlet_temperature_k": float(coolant[-1]),
    }


def _electric_summary(thermo_model: Any) -> dict[str, Any]:
    global_results = thermo_model.get_global_results() or {}
    res = thermo_model.get_tec_results(0) or {}
    j = np.asarray(res.get("J", []), dtype=float)
    current = _finite_float(global_results.get("Iout"))
    voltage = _finite_float(global_results.get("Uout"))
    return {
        "Iout_a": current,
        "Uout_v": voltage,
        "electric_power_w": None if current is None or voltage is None else float(current * voltage),
        "converged": bool(global_results.get("converged", False)),
        "iteration_count": global_results.get("iteration_count"),
        "zero_emission_skipped": bool(global_results.get("zero_emission_skipped", False)),
        "J_mean_A_cm2": None if j.size == 0 else float(np.mean(j)),
        "J_max_A_cm2": None if j.size == 0 else float(np.max(j)),
    }


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _run_one_case(task: dict[str, Any]) -> dict[str, Any]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    power_kw = float(task["power_kwt"])
    pressure_torr = float(task["pressure_torr"])
    wire_resistance_ohm = float(task["wire_resistance_ohm"])
    voltage_v = float(task["voltage_v"])
    duration_s = float(task["duration_s"])
    dt_s = float(task["dt_s"])
    restart_path = Path(task["restart_path"])
    case_restart_path = Path(task["case_restart_path"])
    progress_interval_s = float(task.get("progress_interval_s", DEFAULT_PROGRESS_INTERVAL_S))
    show_progress = bool(task.get("show_progress", False))

    config = LuchauSingleTFEConfig(
        thermal_power_w=power_kw * 1000.0,
        target_voltage_v=voltage_v,
        heater_length_m=0.30,
        cesium_pressure_torr=pressure_torr,
        i_guess_a=max(50.0, float(task["target_current_a"])),
        wire_resistance_ohm=wire_resistance_ohm,
    )
    build = build_luchau_single_tfe(config)
    load_restart_path = case_restart_path if case_restart_path.exists() else restart_path
    loaded_case_restart = case_restart_path.exists()
    build["system"].load_global_state(str(load_restart_path))
    build["tfe"].clear_tec_sources()
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)

    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name=f"Luchau_fit_{power_kw:g}kW_Pcs{pressure_torr:g}_Rw{wire_resistance_ohm:g}_U{voltage_v:g}",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    steps = max(1, int(round(duration_s / dt_s)))
    progress_steps = _progress_interval_steps(duration_s, dt_s, progress_interval_s)
    previous = _temperature_vector(build)
    last_rate = math.nan
    wall_start = time.perf_counter()
    if show_progress:
        print(
            f"[start] P={power_kw:.2f}kW U={voltage_v:.5g}V target={float(task['target_current_a']):.2f}A "
            f"Pcs={pressure_torr:.4g}torr Rw={wire_resistance_ohm:.4g}ohm steps={steps} dt={dt_s:g}s",
            flush=True,
        )
    for step_index in range(1, steps + 1):
        build["system"].step(dt_s)
        current = _temperature_vector(build)
        last_rate = float(np.max(np.abs(current - previous)) / dt_s)
        previous = current
        if show_progress and (step_index == 1 or step_index % progress_steps == 0 or step_index == steps):
            elapsed = time.perf_counter() - wall_start
            sim_time = step_index * dt_s
            percent = 100.0 * step_index / steps
            electric_now = _electric_summary(thermo_model)
            print(
                f"[progress] P={power_kw:.2f}kW U={voltage_v:.5g}V "
                f"sim={sim_time:.2f}/{duration_s:.2f}s ({percent:5.1f}%) "
                f"wall={elapsed:.1f}s I={_format_optional_float(electric_now['Iout_a'], '.2f')}A "
                f"dTdt_max={last_rate:.3e}K/s",
                flush=True,
            )
    wall_time_s = time.perf_counter() - wall_start
    case_restart_path.parent.mkdir(parents=True, exist_ok=True)
    build["system"].save_global_state(str(case_restart_path))

    electric = _electric_summary(thermo_model)
    iout = electric["Iout_a"]
    target = float(task["target_current_a"])
    error = None if iout is None else float(iout - target)
    return {
        "power_kwt": power_kw,
        "pressure_torr": pressure_torr,
        "wire_resistance_ohm": wire_resistance_ohm,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(pressure_torr),
        "row_index": int(task["row_index"]),
        "point_index": int(task["point_index"]),
        "target_voltage_v": voltage_v,
        "target_current_a": target,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "loaded_restart_path": str(load_restart_path),
        "loaded_case_restart": bool(loaded_case_restart),
        "saved_case_restart_path": str(case_restart_path),
        "wall_time_s": wall_time_s,
        "max_temperature_rate_k_s": last_rate,
        **_temperature_stats(build),
        **electric,
        "current_error_a": error,
        "current_abs_error_a": None if error is None else abs(error),
    }


def _metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    valid = [row for row in rows if row.get("Iout_a") is not None and bool(row.get("converged"))]
    errors = np.array([float(row["current_error_a"]) for row in valid], dtype=float)
    if errors.size == 0:
        return {"n_points": len(rows), "n_converged": 0, "rmse_a": None, "mae_a": None, "max_abs_error_a": None, "mean_abs_rel_error": None, "signed_mean_error_a": None}
    abs_errors = np.abs(errors)
    targets = np.array([float(row["target_current_a"]) for row in valid], dtype=float)
    return {
        "n_points": len(rows),
        "n_converged": len(valid),
        "rmse_a": float(np.sqrt(np.mean(errors**2))),
        "mae_a": float(np.mean(abs_errors)),
        "max_abs_error_a": float(np.max(abs_errors)),
        "mean_abs_rel_error": float(np.mean(abs_errors / np.maximum(np.abs(targets), 1.0))),
        "signed_mean_error_a": float(np.mean(errors)),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write empty CSV")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    pressures = _parse_float_list(args.pressure_list, positive=True)
    wires = _parse_float_list(args.wire_list, nonnegative=True)
    voltages = _parse_float_list(args.voltage_list, positive=True, allow_empty=True)
    if voltages:
        target_currents = _parse_float_list(args.target_current_list, allow_empty=True)
        targets = _manual_targets_from_args(
            power_kwt=float(args.thermal_power_kwt),
            voltage_list=voltages,
            target_current_list=target_currents or None,
        )
    else:
        targets = _load_targets(Path(args.target_csv))
    output_dir = CASE_DIR / "runs" / str(args.output_label)
    output_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    missing_restarts = []
    for target in targets:
        restart = _thermal_restart_path(float(target["power_kwt"]))
        if not restart.exists():
            missing_restarts.append(str(restart))
        for pressure in pressures:
            for wire in wires:
                tasks.append({
                    **target,
                    "pressure_torr": float(pressure),
                    "wire_resistance_ohm": float(wire),
                    "duration_s": float(args.duration_s),
                    "dt_s": float(args.dt_s),
                    "progress_interval_s": float(args.progress_interval_s),
                    "show_progress": False,
                    "restart_path": str(restart),
                    "case_restart_path": str(_case_restart_path(output_dir, {**target, "pressure_torr": float(pressure), "wire_resistance_ohm": float(wire)})),
                })
    if missing_restarts:
        raise FileNotFoundError("Missing thermal restarts: " + "; ".join(sorted(set(missing_restarts))))

    wall_start = time.perf_counter()
    rows: list[dict[str, Any]] = []
    if int(args.workers) == 1:
        for task in tasks:
            rows.append(_run_one_case(task))
    else:
        with ProcessPoolExecutor(max_workers=int(args.workers)) as executor:
            futures = [executor.submit(_run_one_case, task) for task in tasks]
            for future in as_completed(futures):
                rows.append(future.result())
    rows.sort(key=lambda r: (float(r["power_kwt"]), float(r["pressure_torr"]), float(r["wire_resistance_ohm"]), int(r["row_index"])))

    for row in sorted(rows, key=lambda r: float(r["target_voltage_v"])):
        print(f"{float(row['target_current_a']):.5f} {float(row['Iout_a']):.5f}", flush=True)


    metric_rows = []
    for power in sorted(set(float(row["power_kwt"]) for row in rows)):
        for pressure in pressures:
            for wire in wires:
                param_rows = [
                    row for row in rows
                    if abs(float(row["power_kwt"]) - power) < 1.0e-12
                    and abs(float(row["pressure_torr"]) - pressure) < 1.0e-12
                    and abs(float(row["wire_resistance_ohm"]) - wire) < 1.0e-15
                ]
                metric_rows.append({
                    "power_kwt": power,
                    "pressure_torr": float(pressure),
                    "wire_resistance_ohm": float(wire),
                    **_metrics(param_rows),
                })
    metric_rows.sort(key=lambda r: (float(r["power_kwt"]), float("inf") if r["rmse_a"] is None else float(r["rmse_a"])))
    best_by_power = []
    for power in sorted(set(float(row["power_kwt"]) for row in metric_rows)):
        best_by_power.append(next(row for row in metric_rows if abs(float(row["power_kwt"]) - power) < 1.0e-12))

    results_csv = output_dir / "fit_results_by_point.csv"
    metrics_csv = output_dir / "fit_metrics_by_power_params.csv"
    best_csv = output_dir / "best_params_by_power.csv"
    _write_csv(results_csv, rows)
    _write_csv(metrics_csv, metric_rows)
    _write_csv(best_csv, best_by_power)

    summary = {
        "case": "center0p30_uniform_fit_all_powers_cs_wire",
        "target_csv": str(args.target_csv),
        "thermal_power_kwt": float(args.thermal_power_kwt),
        "voltage_points_v": voltages,
        "output_dir": str(output_dir),
        "pressure_points_torr": pressures,
        "wire_resistance_points_ohm": wires,
        "duration_s": float(args.duration_s),
        "dt_s": float(args.dt_s),
        "progress_interval_s": float(args.progress_interval_s),
        "workers": int(args.workers),
        "total_wall_time_s": time.perf_counter() - wall_start,
        "best_by_power": best_by_power,
        "files": {
            "results_csv": str(results_csv),
            "metrics_csv": str(metrics_csv),
            "best_csv": str(best_csv),
        },
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
