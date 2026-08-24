"""Hold the final 20% state for one external-heat period."""

from __future__ import annotations

import json
import time
from pathlib import Path

from .run_v14_low_power_fixed_i import LowPowerRunConfig, run_low_power_case


CASE_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = (
    CASE_DIR / "runs"
    / "continuation_40pct_period10s_then_40to20_record1s_Q121086_20260815"
)
PHASE2_DIR = OUTPUT_ROOT / "phase2_40to20_quintic1500_record1s"
PHASE3_DIR = OUTPUT_ROOT / "phase3_20pct_long_cooling_period10s"
EXTERNAL_HEAT_PERIOD_S = 5668.144369
Q_20_W = 121086.033950196


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _wait_for_phase2() -> dict:
    summary_path = PHASE2_DIR / "run_summary.json"
    restart_path = PHASE2_DIR / "final_restart.npz"
    while not summary_path.is_file() or not restart_path.is_file():
        if summary_path.is_file():
            summary = _read_json(summary_path)
            if summary.get("stop_reason") != "completed":
                raise RuntimeError(f"phase 2 did not complete: {summary}")
        time.sleep(30.0)
    summary = _read_json(summary_path)
    if summary.get("stop_reason") != "completed":
        raise RuntimeError(f"phase 2 did not complete: {summary}")
    return summary


def main() -> int:
    phase2 = _wait_for_phase2()
    if PHASE3_DIR.exists():
        raise FileExistsError(f"output already exists: {PHASE3_DIR}")

    current_a = float(phase2["initial_current_A"])
    flow_kg_s = float(phase2["final_flow_setpoint_kg_s"])
    result = run_low_power_case(LowPowerRunConfig(
        restart_in=PHASE2_DIR / "final_restart.npz",
        output_dir=PHASE3_DIR,
        duration_s=EXTERNAL_HEAT_PERIOD_S,
        dt_s=0.05,
        tec_update_interval_s=0.05,
        record_interval_s=10.0,
        ramp_shape="quintic",
        final_power_w=Q_20_W,
        final_flow_kg_s=flow_kg_s,
        checkpoint_interval_s=100.0,
        fixed_current_a=current_a,
        initial_power_w=Q_20_W,
        initial_flow_kg_s=flow_kg_s,
        external_heat_enabled=True,
    ))
    if result.get("stop_reason") != "completed":
        raise RuntimeError(f"20% long-cooling run did not complete: {result}")

    (OUTPUT_ROOT / "long_cooling_summary.json").write_text(
        json.dumps({
            "phase2_run": str(PHASE2_DIR),
            "phase3_long_cooling": result,
            "duration_s": EXTERNAL_HEAT_PERIOD_S,
            "power_setpoint_W": Q_20_W,
            "flow_setpoint_kg_s": flow_kg_s,
            "history_interval_s": 10.0,
            "checkpoint_interval_s": 100.0,
            "external_heat_enabled": True,
        }, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
