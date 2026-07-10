from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

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

VALIDATION_POWERS_KW = (3.26, 3.52, 4.01, 4.49, 5.0, 5.6)
VALIDATION_POWER_LABELS = {
    3.26: "3.26",
    3.52: "3.52",
    4.01: "4.01",
    4.49: "4.49",
    5.0: "5.0",
    5.6: "5.6",
}

DEFAULT_5P6_VOLTAGE_LIST = "1.49793,1.48345,1.47103,1.45862,1.44414,1.43172,1.41517,1.39655,1.37586,1.35517,1.33448,1.31379,1.29517,1.27448,1.25793,1.23931,1.21862,1.19793,1.17724,1.15655,1.14207,1.12345,1.10483,1.08621,1.06759,1.0469,1.02828,1.00759,0.98897,0.96828,0.94966,0.93103,0.91034,0.88966,0.87103,0.85034,0.82966"
DEFAULT_5P6_CURRENT_LIST = "200.23256,208.60465,216.97674,225.34884,232.7907,241.16279,250.46512,257.90698,265.34884,272.7907,281.16279,288.60465,296.04651,303.48837,309.06977,315.5814,323.02326,330.46512,337.90698,345.34884,350.93023,357.44186,364.88372,371.39535,378.83721,387.2093,394.65116,402.09302,409.53488,417.90698,424.4186,431.86047,440.23256,447.67442,455.11628,464.4186,472.7907"
DEFAULT_5P0_VOLTAGE_LIST = "0.72143,0.73903,0.75976,0.77527,0.79286,0.81695,0.8335,0.85759,0.88189,0.91142,0.94117,0.9751,1.00589,1.05533,1.09827,1.13891,1.18164,1.21138,1.25747,1.31256,1.34105,1.37624,1.40703"
DEFAULT_5P0_CURRENT_LIST = "477.5286,467.1685,456.7141,448.3318,441.9274,435.052,427.1406,417.7223,409.34,396.5311,385.135,373.3621,361.0241,342.7525,326.4589,312.6139,297.8272,287.467,269.6664,251.3949,238.0209,218.8075,200.536"
DEFAULT_4P49_VOLTAGE_LIST = "0.56202,0.59282,0.62466,0.65,0.6764,0.71159,0.73568,0.76416,0.79726,0.8268,0.85654,0.88733,0.91917,0.94766,0.9795,1.00589,1.03564,1.06748,1.10267,1.13116,1.163,1.19169,1.21474,1.24322"
DEFAULT_4P49_CURRENT_LIST = "476.0217,462.6477,450.3097,438.5368,427.1406,413.7666,402.9356,391.0685,377.7887,368.3704,355.5615,344.1653,333.3342,321.4672,308.6582,299.805,289.9158,278.0487,268.6304,256.8575,247.9101,233.5943,219.7494,206.9404"
DEFAULT_4P01_VOLTAGE_LIST = "0.48955,0.51929,0.54338,0.57417,0.60832,0.64016,0.67744,0.70488,0.74007,0.76856,0.79831,0.83119,0.86534,0.89823,0.93132,0.96316,1.00149,1.03669,1.07523,1.10037,1.13116,1.16195,1.17955"
DEFAULT_4P01_CURRENT_LIST = "478.4705,464.6255,448.3318,437.5008,424.6919,411.7888,398.038,385.7001,373.3621,361.9659,347.1792,335.3121,323.9159,309.6001,297.2621,286.431,274.5639,263.262,250.3589,239.5278,225.6829,212.4031,204.0208"
DEFAULT_3P52_VOLTAGE_LIST = "0.40827,0.42692,0.45331,0.48515,0.50945,0.53018,0.55993,0.59177,0.62361,0.65,0.67954,0.71599,0.74343,0.77966,0.8092,0.8379,0.86869,0.90158,0.93886,0.96861,0.9994,1.02789,1.06517,1.09157"
DEFAULT_3P52_CURRENT_LIST = "473.5729,463.6837,449.8388,434.5811,423.1849,410.376,397.473,383.7222,370.8191,361.0241,347.6501,334.8412,325.8938,313.5558,302.7247,289.9158,281.5335,267.2177,255.8215,246.4032,233.1234,223.2341,212.4031,202.0429"
DEFAULT_3P26_VOLTAGE_LIST = "0.33789,0.36659,0.38858,0.41812,0.44242,0.47865,0.49939,0.53018,0.56202,0.59512,0.62026,0.6521,0.68289,0.71599,0.74217,0.77422,0.80375,0.84439,0.86974,0.90493,0.93342,0.96966,1.00045,1.02244"
DEFAULT_3P26_CURRENT_LIST = "478.4705,465.6615,450.8748,438.5368,424.221,409.34,397.473,382.6862,368.3704,355.5615,342.7525,330.4146,319.0184,305.6444,293.8715,282.0044,272.5861,258.8354,249.417,238.0209,226.7189,216.3588,204.9626,198.5581"

DEFAULT_TARGET_SERIES_BY_POWER = {
    3.26: {"voltage_list": DEFAULT_3P26_VOLTAGE_LIST, "target_current_list": DEFAULT_3P26_CURRENT_LIST},
    3.52: {"voltage_list": DEFAULT_3P52_VOLTAGE_LIST, "target_current_list": DEFAULT_3P52_CURRENT_LIST},
    4.01: {"voltage_list": DEFAULT_4P01_VOLTAGE_LIST, "target_current_list": DEFAULT_4P01_CURRENT_LIST},
    4.49: {"voltage_list": DEFAULT_4P49_VOLTAGE_LIST, "target_current_list": DEFAULT_4P49_CURRENT_LIST},
    5.0: {"voltage_list": DEFAULT_5P0_VOLTAGE_LIST, "target_current_list": DEFAULT_5P0_CURRENT_LIST},
    5.6: {"voltage_list": DEFAULT_5P6_VOLTAGE_LIST, "target_current_list": DEFAULT_5P6_CURRENT_LIST},
}

FIRST_RUN_DURATION_S = 800.0
FIRST_RUN_DT_S = 0.05
CONTINUATION_DURATION_S = 100.0
CONTINUATION_DT_S = 0.1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validation voltage scan for one thermal power, one wire resistance, and one cesium pressure.")
    parser.add_argument("--thermal-power-kwt", type=float, default=3.26)
    parser.add_argument("--wire-resistance-ohm", type=float, default=1.5e-4)
    parser.add_argument("--cesium-pressure-torr", type=float, default=1.45)
    parser.add_argument("--voltage-list", type=str, default=None)
    parser.add_argument("--target-current-list", type=str, default=None)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--force-first-run", action="store_true", help="Ignore validation restarts and run the 800 s first-run schedule.")
    return parser


def _parse_float_list(value: Any, *, allow_empty: bool = False) -> list[float]:
    if value is None:
        values: list[float] = []
    elif isinstance(value, (int, float, np.floating)):
        values = [float(value)]
    elif isinstance(value, (list, tuple)):
        values = [float(item) for item in value]
    else:
        values = [float(part.strip()) for part in str(value).split(",") if part.strip()]
    if not values and not allow_empty:
        raise ValueError("list argument must contain at least one value.")
    return values


def _power_match(power_kwt: float) -> float:
    for candidate in VALIDATION_POWERS_KW:
        if abs(float(power_kwt) - candidate) < 1.0e-9:
            return candidate
    raise ValueError(f"unsupported validation power {power_kwt}; expected one of {VALIDATION_POWERS_KW}")


def _restart_token(value: float) -> str:
    text = f"{float(value):.8g}"
    return text.replace("-", "m").replace("+", "").replace(".", "p")


def validation_directory(power_kwt: float) -> Path:
    power = _power_match(power_kwt)
    index = VALIDATION_POWERS_KW.index(power) + 1
    return CASE_DIR / f"Validation_{index}-{VALIDATION_POWER_LABELS[power]}"


def validation_restart_path(power_kwt: float, voltage_v: float) -> Path:
    power = _power_match(power_kwt)
    power_token = VALIDATION_POWER_LABELS[power].replace(".", "p")
    voltage_token = _restart_token(voltage_v)
    return validation_directory(power) / "restarts" / f"restart_{power_token}_{voltage_token}.npz"


def validation_results_path(power_kwt: float, wire_resistance_ohm: float, cesium_pressure_torr: float) -> Path:
    power = _power_match(power_kwt)
    power_token = VALIDATION_POWER_LABELS[power].replace(".", "p")
    return validation_directory(power) / "results" / (
        f"results_{power_token}_Rw{_restart_token(wire_resistance_ohm)}_Cs{_restart_token(cesium_pressure_torr)}.csv"
    )


def ensure_validation_directories() -> None:
    for power in VALIDATION_POWERS_KW:
        base = validation_directory(power)
        (base / "restarts").mkdir(parents=True, exist_ok=True)
        (base / "results").mkdir(parents=True, exist_ok=True)


def select_run_schedule(restart_path: Path, *, force_first_run: bool = False) -> dict[str, Any]:
    is_first_run = force_first_run or not restart_path.exists()
    if is_first_run:
        return {"is_first_run": True, "duration_s": FIRST_RUN_DURATION_S, "dt_s": FIRST_RUN_DT_S}
    return {"is_first_run": False, "duration_s": CONTINUATION_DURATION_S, "dt_s": CONTINUATION_DT_S}


def _thermal_restart_path(power_kwt: float) -> Path:
    power = _power_match(power_kwt)
    label = VALIDATION_POWER_LABELS[power].replace(".", "p")
    return CASE_DIR / "runs" / f"thermal_center0p30_uniform_{label}kw_1000s_dt0p05" / "preheat_restart.npz"


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


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


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
        "J_mean_A_cm2": None if j.size == 0 else float(np.mean(j)),
        "J_max_A_cm2": None if j.size == 0 else float(np.max(j)),
    }


def _run_voltage_case(task: dict[str, Any]) -> dict[str, Any]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    power_kwt = float(task["power_kwt"])
    voltage_v = float(task["voltage_v"])
    target_current_a = float(task["target_current_a"])
    wire_resistance_ohm = float(task["wire_resistance_ohm"])
    cesium_pressure_torr = float(task["cesium_pressure_torr"])
    restart_path = Path(task["restart_path"])
    thermal_restart_path = Path(task["thermal_restart_path"])
    schedule = select_run_schedule(restart_path, force_first_run=bool(task.get("force_first_run", False)))

    config = LuchauSingleTFEConfig(
        thermal_power_w=power_kwt * 1000.0,
        target_voltage_v=voltage_v,
        heater_length_m=0.30,
        cesium_pressure_torr=cesium_pressure_torr,
        i_guess_a=max(50.0, target_current_a),
        wire_resistance_ohm=wire_resistance_ohm,
    )
    build = build_luchau_single_tfe(config)
    load_path = thermal_restart_path if schedule["is_first_run"] else restart_path
    build["system"].load_global_state(str(load_path))
    build["tfe"].clear_tec_sources()
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)

    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name=f"Validation_P{power_kwt:g}_U{voltage_v:g}_Cs{cesium_pressure_torr:g}_Rw{wire_resistance_ohm:g}",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    duration_s = float(schedule["duration_s"])
    dt_s = float(schedule["dt_s"])
    steps = max(1, int(round(duration_s / dt_s)))
    previous = _temperature_vector(build)
    last_rate = math.nan
    wall_start = time.perf_counter()
    for _ in range(steps):
        build["system"].step(dt_s)
        current = _temperature_vector(build)
        last_rate = float(np.max(np.abs(current - previous)) / dt_s)
        previous = current
    wall_time_s = time.perf_counter() - wall_start

    restart_path.parent.mkdir(parents=True, exist_ok=True)
    build["system"].save_global_state(str(restart_path))

    electric = _electric_summary(thermo_model)
    iout = electric["Iout_a"]
    error = None if iout is None else float(iout - target_current_a)
    return {
        "power_kwt": power_kwt,
        "voltage_v": voltage_v,
        "target_current_a": target_current_a,
        "wire_resistance_ohm": wire_resistance_ohm,
        "cesium_pressure_torr": cesium_pressure_torr,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(cesium_pressure_torr),
        "is_first_run": bool(schedule["is_first_run"]),
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "loaded_restart_path": str(load_path),
        "saved_restart_path": str(restart_path),
        "wall_time_s": wall_time_s,
        "max_temperature_rate_k_s": last_rate,
        **_temperature_stats(build),
        **electric,
        "current_error_a": error,
        "current_abs_error_a": None if error is None else abs(error),
    }


def _build_tasks(args: argparse.Namespace) -> list[dict[str, Any]]:
    power = _power_match(float(args.thermal_power_kwt))
    default_series = DEFAULT_TARGET_SERIES_BY_POWER[power]
    voltage_list = args.voltage_list if args.voltage_list is not None else default_series["voltage_list"]
    target_current_list = (
        args.target_current_list if args.target_current_list is not None else default_series["target_current_list"]
    )
    voltages = _parse_float_list(voltage_list)
    target_currents = _parse_float_list(target_current_list, allow_empty=True)
    if target_currents and len(target_currents) != len(voltages):
        raise ValueError("target-current-list length must match voltage-list length.")
    if not target_currents:
        target_currents = [0.0] * len(voltages)

    thermal_restart = _thermal_restart_path(power)
    if not thermal_restart.exists():
        raise FileNotFoundError(f"Missing thermal restart for {power} kW: {thermal_restart}")

    tasks = []
    for voltage, target_current in zip(voltages, target_currents):
        tasks.append({
            "power_kwt": power,
            "voltage_v": float(voltage),
            "target_current_a": float(target_current),
            "wire_resistance_ohm": float(args.wire_resistance_ohm),
            "cesium_pressure_torr": float(args.cesium_pressure_torr),
            "thermal_restart_path": str(thermal_restart),
            "restart_path": str(validation_restart_path(power, float(voltage))),
            "force_first_run": bool(args.force_first_run),
        })
    return tasks


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
    ensure_validation_directories()
    tasks = _build_tasks(args)

    rows: list[dict[str, Any]] = []
    wall_start = time.perf_counter()
    if int(args.workers) == 1:
        for task in tasks:
            rows.append(_run_voltage_case(task))
    else:
        with ProcessPoolExecutor(max_workers=int(args.workers)) as executor:
            futures = [executor.submit(_run_voltage_case, task) for task in tasks]
            for future in as_completed(futures):
                rows.append(future.result())

    rows.sort(key=lambda row: float(row["voltage_v"]))
    for row in rows:
        print(
            f"{float(row['voltage_v']):.5f} "
            f"{float(row['target_current_a']):.5f} "
            f"{float(row['Iout_a']):.5f} "
            f"{float(row['max_temperature_rate_k_s']):.6e}",
            flush=True,
        )

    results_csv = validation_results_path(
        float(args.thermal_power_kwt),
        float(args.wire_resistance_ohm),
        float(args.cesium_pressure_torr),
    )
    _write_csv(results_csv, rows)

    summary = {
        "thermal_power_kwt": float(args.thermal_power_kwt),
        "wire_resistance_ohm": float(args.wire_resistance_ohm),
        "cesium_pressure_torr": float(args.cesium_pressure_torr),
        "workers": int(args.workers),
        "total_wall_time_s": time.perf_counter() - wall_start,
        "results_csv": str(results_csv),
        "n_points": len(rows),
        "n_first_run": sum(1 for row in rows if bool(row["is_first_run"])),
        "n_continuation": sum(1 for row in rows if not bool(row["is_first_run"])),
    }
    summary_path = results_csv.with_suffix(".summary.json")
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
