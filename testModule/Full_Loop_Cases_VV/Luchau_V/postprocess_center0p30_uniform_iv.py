from __future__ import annotations

import csv
import json
from dataclasses import replace
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


RUN_DIR = CASE_DIR / "runs" / "lookup_center0p30_uniform_3260W_cs3p55_1p00V_plus200s_dt0p05"
RESTART = RUN_DIR / "final_restart.npz"
VOLTAGE_POINTS_V = (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
I_GUESSES_A = (50.0, 100.0, 180.0, 300.0, 500.0, 800.0, 1000.0)


def _temperature_stats(build: dict[str, Any]) -> dict[str, float]:
    tfe = build["tfe"]
    pellet = np.asarray(tfe.solids["pellet"].T, dtype=float)
    emitter = np.asarray(tfe.solids["emitter"].T, dtype=float)
    collector = np.asarray(tfe.solids["collector"].T, dtype=float)
    coolant = np.array([float(vol.T) for vol in build["channel"].volumes], dtype=float)
    return {
        "pellet_temperature_mean_k": float(np.mean(pellet)),
        "pellet_temperature_max_k": float(np.max(pellet)),
        "emitter_temperature_mean_k": float(np.mean(emitter)),
        "emitter_temperature_max_k": float(np.max(emitter)),
        "collector_temperature_mean_k": float(np.mean(collector)),
        "collector_temperature_max_k": float(np.max(collector)),
        "coolant_outlet_temperature_k": float(coolant[-1]),
    }


def _surface_temperatures(build: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    tfe = build["tfe"]
    return (
        np.asarray(tfe.solids["emitter"].boundaries["right"].T_surface, dtype=float).reshape(1, tfe.mesh.n_axial),
        np.asarray(tfe.solids["collector"].boundaries["left"].T_surface, dtype=float).reshape(1, tfe.mesh.n_axial),
    )


def _electric_summary(model: Any) -> dict[str, Any]:
    global_results = model.get_global_results() or {}
    res = model.get_tec_results(0) or {}
    j = np.asarray(res.get("J", []), dtype=float)
    ue = np.asarray(res.get("UE", []), dtype=float)
    uc = np.asarray(res.get("UC", []), dtype=float)
    current = _finite_float(global_results.get("Iout"))
    voltage = _finite_float(global_results.get("Uout"))
    ueuc = ue - uc if ue.size and uc.size else np.array([], dtype=float)
    return {
        "Iout_a": current,
        "Uout_v": voltage,
        "electric_power_w": None if current is None or voltage is None else float(current * voltage),
        "converged": bool(global_results.get("converged", False)),
        "iteration_count": global_results.get("iteration_count"),
        "J_mean_A_cm2": None if j.size == 0 else float(np.mean(j)),
        "J_max_A_cm2": None if j.size == 0 else float(np.max(j)),
        "UE_minus_UC_mean_v": None if ueuc.size == 0 else float(np.mean(ueuc)),
    }


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _calculate_with_retries(build: dict[str, Any], base_config: LuchauSingleTFEConfig, voltage: float) -> dict[str, Any]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    emitter_t, collector_t = _surface_temperatures(build)
    attempts: list[dict[str, Any]] = []
    for guess in I_GUESSES_A:
        config = replace(base_config, target_voltage_v=float(voltage), i_guess_a=float(guess))
        model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
        configure_luchau_thermocalc(model, build, config)
        model.set_temperatures(emitter_t, collector_t)
        model.calculate(verbose=False)
        result = _electric_summary(model)
        result["i_guess_a"] = float(guess)
        attempts.append(result)
        if result["converged"] and result["Iout_a"] is not None:
            return result
    result = attempts[-1]
    result["failed_attempts"] = attempts
    return result


def main() -> int:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base_config = LuchauSingleTFEConfig(
        thermal_power_w=3260.0,
        target_voltage_v=1.0,
        heater_length_m=0.30,
        cesium_pressure_torr=3.55,
        i_guess_a=185.0,
        wire_resistance_ohm=0.0,
    )
    build = build_luchau_single_tfe(base_config)
    build["system"].load_global_state(str(RESTART))
    build["tfe"].update_neutronic_power(float(base_config.thermal_power_w), alpha=1.0)
    lookup = _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)
    stats = _temperature_stats(build)
    rows = []
    for voltage in VOLTAGE_POINTS_V:
        rows.append(
            {
                "target_voltage_v": float(voltage),
                "cesium_pressure_torr": float(base_config.cesium_pressure_torr),
                "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(base_config.cesium_pressure_torr)),
                **stats,
                **_calculate_with_retries(build, base_config, float(voltage)),
            }
        )

    csv_path = RUN_DIR / "iv_scan_0p4_to_1p5_same_thermal_state.csv"
    serializable_rows = [{k: v for k, v in row.items() if k != "failed_attempts"} for row in rows]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(serializable_rows[0].keys()))
        writer.writeheader()
        writer.writerows(serializable_rows)

    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=160)
    ax.plot(
        [row["target_voltage_v"] for row in serializable_rows],
        [row["Iout_a"] if row["converged"] else np.nan for row in serializable_rows],
        marker="o",
        lw=1.6,
    )
    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Current (A)")
    ax.set_title("Center 0.3 m uniform, 3.26 kW, Pcs=3.55 torr")
    ax.grid(True, alpha=0.3)
    plot_path = RUN_DIR / "iv_scan_current_voltage.png"
    fig.tight_layout()
    fig.savefig(plot_path)

    summary = {
        "case": "center0p30_uniform_same_state_iv_scan",
        "source_restart": str(RESTART),
        "lookup": lookup,
        "rows": rows,
        "csv": str(csv_path),
        "plot": str(plot_path),
    }
    with (RUN_DIR / "iv_scan_0p4_to_1p5_same_thermal_state_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
