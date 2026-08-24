"""Continue the half-radiator case from its post-scram checkpoint."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "testModule").is_dir()
)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (  # noqa: E402
    DebugRunConfig,
    build_debug_case,
    collect_metrics,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1.run_v14_210kw_fixed_power_loca_1 import (  # noqa: E402
    _neutronics_metrics,
    append_postprocessing_histories,
    build_snapshot_payload,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.run_v14_helium_depressurization import (  # noqa: E402
    collect_temperature_peaks,
    find_nonfinite_model_state,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (  # noqa: E402
    ReactivityControlRunConfig,
    load_baseline_debug_config,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.heatpipe_failure_accident import (  # noqa: E402
    SUMMARY_FIELDS,
    _temperature_metrics,
    apply_heatpipe_failure,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_heatpipe_half_radiator_failure.run_v14_heatpipe_half_radiator_failure import (  # noqa: E402
    FAILURE_MAP,
    FAILURE_MODE,
)


CASE_DIR = Path(__file__).resolve().parent
SOURCE_DIR = CASE_DIR / "runs" / "cycle_parallel_20260810"
RESTART = SOURCE_DIR / "checkpoint_tplus_05300.000s.npz"
OUTPUT_DIR = CASE_DIR / "runs" / "continuation_from_tplus5300_20260813"
ORIGINAL_ACCIDENT_START_S = 19864.69999990125
TARGET_ELAPSED_S = 5668.144369
SCRAM_DOLLARS = -2.0
DT_S = 0.05
HISTORY_INTERVAL_S = 10.0
CHECKPOINT_INTERVAL_S = 100.0


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def _append_summary(path: Path, row: dict) -> None:
    new_file = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        writer.writerow(row)


def _row(build: dict, accident: dict, *, start_time: float, dt_s: float) -> dict:
    system, core = build["system"], build["core"]
    rho = SCRAM_DOLLARS * float(core.point_reactor.beta_total)
    row = {
        **collect_metrics(build, stage_index=1, dt_s=dt_s),
        **_neutronics_metrics(core, rho, SCRAM_DOLLARS),
        **_temperature_metrics(collect_temperature_peaks(core)),
        "time_s": float(system.global_time),
        "accident_elapsed_s": float(system.global_time) - ORIGINAL_ACCIDENT_START_S,
        "case_mode": "fixed_power_then_scram",
        "failure_mode": FAILURE_MODE,
        "external_heat_enabled": True,
        "external_heat_period_s": 5668.144369,
        "external_heat_phase_s": (
            (float(system.global_time) - 13864.19999998857) % 5668.144369
        ),
        "tec_electrical_calculation_enabled": True,
        "fixed_power_control_active": False,
        "scram_active": True,
        "scram_time_absolute_s": 20059.799999898427,
        "scram_elapsed_s": 195.09999999717547,
        "scram_trigger_component": "collector",
        "scram_trigger_representative": "Ring1",
        "scram_trigger_axial_position_m": 0.35906,
        "scram_trigger_actual_K": 1023.0000866203926,
        "scram_trigger_limit_K": 1023.0,
        "external_reactivity": rho,
        "external_reactivity_dollars": SCRAM_DOLLARS,
        "core_power_W": float(core.last_total_core_power),
        "stop_reason": "completed",
    }
    return row


def main() -> int:
    if not RESTART.is_file():
        raise FileNotFoundError(RESTART)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    history = OUTPUT_DIR / "history.csv"
    if history.exists():
        raise FileExistsError(history)

    source_config = json.loads((SOURCE_DIR / "run_config.json").read_text(encoding="utf-8"))
    runtime = ReactivityControlRunConfig(
        restart_in=RESTART,
        output_dir=OUTPUT_DIR,
        duration_s=TARGET_ELAPSED_S,
        dt_s=DT_S,
        record_interval_s=HISTORY_INTERVAL_S,
        checkpoint_interval_s=CHECKPOINT_INTERVAL_S,
        min_fluid_temperature_stop_k=None,
        external_heat_enabled=True,
        external_heat_period_s=float(source_config["external_heat_period_s"]),
        external_heat_time_origin_s=float(source_config["external_heat_time_origin_s"]),
    )
    debug, _ = load_baseline_debug_config(runtime)
    build = build_debug_case(debug, apply_fixed_power=False)
    system, core = build["system"], build["core"]
    if not core.has_point_reactor:
        raise RuntimeError("post-scram restart did not restore point kinetics")
    signature = apply_heatpipe_failure(build, FAILURE_MAP)
    reference_fluid = {
        key: np.asarray(getattr(system.fluid_solver, key + "_vec"), dtype=float).copy()
        for key in ("T", "P", "h", "W")
    }
    accident = {
        "reference_fluid": reference_fluid,
        "feedback_reference_total": float(core.feedback_reference_result.total),
        "tec_open_circuit_active": False,
        "tec_open_circuit_time_s": float("nan"),
    }
    start_time = float(system.global_time)
    end_time = ORIGINAL_ACCIDENT_START_S + TARGET_ELAPSED_S
    rows = []
    next_record = start_time
    next_checkpoint = start_time + CHECKPOINT_INTERVAL_S

    run_config = dict(source_config)
    run_config.update({
        "case": "V14_210kW_heatpipe_half_radiator_failure_continuation",
        "continuation": True,
        "continuation_source_restart": str(RESTART),
        "original_accident_start_time_s": ORIGINAL_ACCIDENT_START_S,
        "target_accident_elapsed_s": TARGET_ELAPSED_S,
        "scram_reactivity_dollars": SCRAM_DOLLARS,
        "history_interval_s": HISTORY_INTERVAL_S,
        "checkpoint_interval_s": CHECKPOINT_INTERVAL_S,
        "point_kinetics_enabled": True,
        "fixed_power_control_active": False,
        "scram_active": True,
        "failure_signature": signature,
    })
    _write_json(OUTPUT_DIR / "run_config.json", run_config)
    _write_json(OUTPUT_DIR / "continuation_event.json", {
        "source_checkpoint": str(RESTART),
        "source_accident_elapsed_s": float(start_time - ORIGINAL_ACCIDENT_START_S),
        "target_accident_elapsed_s": TARGET_ELAPSED_S,
        "scram_active": True,
        "external_reactivity_dollars": SCRAM_DOLLARS,
    })

    def record(row: dict) -> None:
        rho = float(row["external_reactivity"])
        payload = build_snapshot_payload(
            build, accident, start_time=ORIGINAL_ACCIDENT_START_S,
            coolant_present=True, hydraulic_solve_enabled=True,
            external_reactivity=rho, external_reactivity_dollars=SCRAM_DOLLARS,
        )
        append_postprocessing_histories(OUTPUT_DIR, payload)
        _append_summary(history, row)
        rows.append(row)
        _write_json(OUTPUT_DIR / "latest_state.json", {
            "case": run_config["case"],
            "latest_checkpoint_path": str(OUTPUT_DIR / "latest_restart.npz"),
            "latest_metrics": row,
        })

    row = _row(build, accident, start_time=ORIGINAL_ACCIDENT_START_S, dt_s=0.0)
    record(row)
    while float(system.global_time) < end_time - 1.0e-9:
        now = float(system.global_time)
        dt = min(DT_S, end_time - now, max(1.0e-12, next_record - now) if next_record > now else DT_S)
        rho = SCRAM_DOLLARS * float(core.point_reactor.beta_total)
        system.step(
            dt,
            inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=True,
            fluid_max_iter=int(debug.fluid_max_iter),
            reactivity_control=rho,
        )
        if find_nonfinite_model_state(build) is not None:
            raise RuntimeError("nonfinite state during continuation")
        now = float(system.global_time)
        if now >= next_record - 1.0e-9 or now >= end_time - 1.0e-9:
            record(_row(build, accident, start_time=ORIGINAL_ACCIDENT_START_S, dt_s=dt))
            next_record = now + HISTORY_INTERVAL_S
        if now >= next_checkpoint - 1.0e-9:
            checkpoint = OUTPUT_DIR / f"checkpoint_tplus_{now - ORIGINAL_ACCIDENT_START_S:09.3f}s.npz"
            system.save_global_state(str(checkpoint))
            system.save_global_state(str(OUTPUT_DIR / "latest_restart.npz"))
            next_checkpoint += CHECKPOINT_INTERVAL_S

    final = OUTPUT_DIR / "final_restart.npz"
    system.save_global_state(str(final))
    result = {
        "case": run_config["case"],
        "source_checkpoint": str(RESTART),
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "accident_elapsed_s": float(system.global_time - ORIGINAL_ACCIDENT_START_S),
        "scram_time_absolute_s": 20059.799999898427,
        "scram_trigger_component": "collector",
        "scram_trigger_actual_K": 1023.0000866203926,
        "external_reactivity_dollars": SCRAM_DOLLARS,
        "stop_reason": "completed",
        "final_restart_path": str(final),
        "history_path": str(history),
        "latest_metrics": rows[-1],
    }
    _write_json(OUTPUT_DIR / "run_summary.json", result)
    _write_json(OUTPUT_DIR / "latest_state.json", {
        "case": run_config["case"],
        "latest_checkpoint_path": str(final),
        "latest_metrics": rows[-1],
    })
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
