"""V15 V71 fixed-voltage lookup continuation from the 1000 s thermal restart."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_cold_start import (  # noqa: E402
    ColdStartConfig,
    _apply_core_power,
    _write_json,
    build_cold_start_case,
    collect_metrics,
)


@dataclass(frozen=True)
class ContinuationConfig:
    output_dir: Path
    duration_s: float = 1000.0
    dt_s: float = 0.05
    record_interval_s: float = 10.0
    checkpoint_interval_s: float = 100.0
    core_power_w: float = 106000.0
    tec_voltage_v: float = 27.2
    tec_current_guess_a: float = 385.6
    thermo_update_interval_s: float = 0.8
    lookup_db: Path = Path("ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr")
    lookup_regions: Sequence[str] = ("core", "startup", "high_power", "accident")


def configure_tec(build: Dict[str, Any], config: ContinuationConfig) -> None:
    core = build["core"]
    core.tec_lookup_enabled = True
    core.tec_lookup_db = str(config.lookup_db)
    core.tec_lookup_regions = tuple(config.lookup_regions)
    core.enable_tec_coupled = True
    core._build_thermo_calc()
    core.setup_tec_circuit(
        "fixed_u",
        float(config.tec_voltage_v),
        I_guess=float(config.tec_current_guess_a),
        topology="series",
    )
    core.thermo_update_interval = float(config.thermo_update_interval_s)


def _metrics(build: Dict[str, Any], dt_s: float, solve_ms: float, solve_count: int) -> Dict[str, Any]:
    row = collect_metrics(build, dt_s)
    results = build["core"].get_tec_circuit_global_results().get("main") or {}
    current = float(results.get("Iout", np.nan))
    voltage = float(results.get("Uout", np.nan))
    row.update({
        "tec_current_A": current,
        "tec_voltage_V": voltage,
        "tec_power_W": current * voltage,
        "tec_converged": bool(results.get("converged", False)),
        "tec_iteration_count": int(results.get("iteration_count", 0)),
        "tec_zero_emission_skipped": bool(results.get("zero_emission_skipped", False)),
        "tec_calculate_ms_last": float(solve_ms),
        "tec_solve_count": int(solve_count),
    })
    return row


def run(config: ContinuationConfig, restart_path: Path) -> Dict[str, Any]:
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = ColdStartConfig(
        output_dir=out_dir,
        initial_temperature_k=723.0,
        space_temperature_k=4.0,
        core_power_w=float(config.core_power_w),
        tec_voltage_v=float(config.tec_voltage_v),
        tec_current_guess_a=float(config.tec_current_guess_a),
        lookup_db=Path(config.lookup_db),
        lookup_regions=tuple(config.lookup_regions),
    )
    build = build_cold_start_case(base)
    system = build["system"]
    system.load_global_state(str(restart_path))
    configure_tec(build, config)
    core = build["core"]
    core.post_step(0.0, float(system.global_time))

    solve_count = 0
    total_solve_ms = 0.0
    max_solve_ms = 0.0
    original_calculate = core.thermo_calc.calculate

    def timed_calculate(verbose: bool = False) -> float:
        nonlocal solve_count, total_solve_ms, max_solve_ms
        started = time.perf_counter()
        elapsed_ms = float(original_calculate(verbose=verbose))
        measured_ms = (time.perf_counter() - started) * 1000.0
        solve_count += 1
        total_solve_ms += measured_ms
        max_solve_ms = max(max_solve_ms, measured_ms)
        return elapsed_ms

    core.thermo_calc.calculate = timed_calculate
    history_path = out_dir / "coupled_history.csv"
    initial = _metrics(build, 0.0, 0.0, solve_count)
    fields = list(initial)
    with history_path.open("w", newline="", encoding="utf-8") as handle:
        csv.DictWriter(handle, fieldnames=fields).writeheader()

    start_time = float(system.global_time)
    final_time = start_time + float(config.duration_s)
    next_record = start_time
    next_checkpoint = start_time + float(config.checkpoint_interval_s)
    latest = initial
    while system.global_time < final_time - 1.0e-12:
        dt = min(float(config.dt_s), final_time - float(system.global_time))
        _apply_core_power(build, config.core_power_w)
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=100)
        _apply_core_power(build, config.core_power_w)

        if system.global_time + 1.0e-12 >= next_record or system.global_time >= final_time - 1.0e-12:
            latest = _metrics(build, dt, max_solve_ms, solve_count)
            with history_path.open("a", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=fields).writerow(latest)
            print(json.dumps(latest, sort_keys=True), flush=True)
            next_record = system.global_time + float(config.record_interval_s)

        if config.checkpoint_interval_s > 0.0 and system.global_time + 1.0e-12 >= next_checkpoint:
            system.save_global_state(str(out_dir / f"coupled_checkpoint_t{int(round(system.global_time)):04d}s.npz"))
            next_checkpoint += float(config.checkpoint_interval_s)

    restart_out = out_dir / "coupled_1000s_restart.npz"
    system.save_global_state(str(restart_out))
    summary = {
        "case": "V15_V71_cold_start_tec_lookup_fixed_u",
        "absolute_start_time_s": start_time,
        "absolute_end_time_s": float(system.global_time),
        "dt_s": float(config.dt_s),
        "tec_update_interval_s": float(config.thermo_update_interval_s),
        "tec_solve_count": int(solve_count),
        "tec_total_measured_ms": float(total_solve_ms),
        "tec_average_measured_ms": float(total_solve_ms / solve_count) if solve_count else 0.0,
        "tec_max_measured_ms": float(max_solve_ms),
        "restart_path": str(restart_out),
        "latest_metrics": latest,
    }
    _write_json(out_dir / "coupled_summary.json", summary)
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=1000.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--checkpoint-interval", type=float, default=100.0)
    parser.add_argument("--tec-update-interval", type=float, default=0.8)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    config = ContinuationConfig(
        output_dir=args.output_dir,
        duration_s=float(args.duration),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        thermo_update_interval_s=float(args.tec_update_interval),
        lookup_db=REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr",
    )
    result = run(config, args.restart)
    print(json.dumps(result["latest_metrics"], indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
