from __future__ import annotations

import argparse
import csv
import json
import time
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
    build_luchau_single_tfe,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import json_default


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare thermal-only center-0.3m uniform heating restarts for Luchau single TFE.")
    parser.add_argument("--power-list-kw", type=str, required=True)
    parser.add_argument("--duration-s", type=float, default=1000.0)
    parser.add_argument("--dt-s", type=float, default=0.05)
    parser.add_argument("--force", action="store_true")
    return parser


def _parse_power_list(value: str) -> list[float]:
    powers = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not powers:
        raise ValueError("power-list-kw must contain at least one value.")
    for power in powers:
        if power <= 0.0:
            raise ValueError("powers must be positive.")
    return powers


def _temperature_vector(build: dict[str, Any]) -> np.ndarray:
    return np.concatenate(
        [np.asarray(solid.T, dtype=float).ravel() for solid in build["tfe"].solids.values() if hasattr(solid, "T")]
    )


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


def _label_for_power(power_kw: float) -> str:
    return f"thermal_center0p30_uniform_{str(power_kw).replace('.', 'p')}kw_1000s_dt0p05"


def _run_power(power_kw: float, duration_s: float, dt_s: float, force: bool) -> dict[str, Any]:
    output_dir = CASE_DIR / "runs" / _label_for_power(power_kw)
    restart_path = output_dir / "preheat_restart.npz"
    summary_path = output_dir / "summary.json"
    if restart_path.exists() and summary_path.exists() and not force:
        with summary_path.open("r", encoding="utf-8") as f:
            summary = json.load(f)
        summary["skipped_existing"] = True
        return summary

    output_dir.mkdir(parents=True, exist_ok=True)
    config = LuchauSingleTFEConfig(
        thermal_power_w=float(power_kw) * 1000.0,
        target_voltage_v=1.0,
        heater_length_m=0.30,
        cesium_pressure_torr=1.0,
        wire_resistance_ohm=0.0,
    )
    build = build_luchau_single_tfe(config)
    steps = int(round(float(duration_s) / float(dt_s)))
    previous = _temperature_vector(build)
    history = []
    wall_start = time.perf_counter()
    for step in range(steps):
        build["system"].step(float(dt_s))
        current = _temperature_vector(build)
        rate = float(np.max(np.abs(current - previous)) / float(dt_s))
        previous = current
        if step == 0 or (step + 1) % 100 == 0 or step + 1 == steps:
            history.append({
                "step": int(step + 1),
                "elapsed_s": float((step + 1) * float(dt_s)),
                "absolute_time_s": float(build["system"].global_time),
                "max_temperature_rate_k_s": rate,
                **_temperature_stats(build),
            })
    wall_time_s = time.perf_counter() - wall_start
    build["system"].save_global_state(str(restart_path))

    history_csv = output_dir / "history.csv"
    with history_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    summary = {
        "case": "center0p30_uniform_thermal_only_preheat",
        "config": asdict(config),
        "power_kw": float(power_kw),
        "duration_s": float(duration_s),
        "dt_s": float(dt_s),
        "steps": steps,
        "wall_time_s": wall_time_s,
        "final": history[-1],
        "files": {"restart": str(restart_path), "history_csv": str(history_csv)},
    }
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return summary


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summaries = [_run_power(power, float(args.duration_s), float(args.dt_s), bool(args.force)) for power in _parse_power_list(args.power_list_kw)]
    print(json.dumps({"summaries": summaries}, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
