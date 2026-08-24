"""Advance the calibrated 20% endpoint for one more external-heat period."""

from __future__ import annotations

import json
from pathlib import Path

from .run_v14_low_power_fixed_i import LowPowerRunConfig, run_low_power_case


CASE_DIR = Path(__file__).resolve().parent
SOURCE_DIR = (
    CASE_DIR / "runs"
    / "continuation_40pct_period10s_then_40to20_record1s_Q121086_20260815"
    / "phase3_20pct_long_cooling_period10s"
)
OUTPUT_DIR = (
    CASE_DIR / "runs"
    / "continuation_20pct_additional_period_record50s_Q121086_20260816"
)
EXTERNAL_HEAT_PERIOD_S = 5668.144369
Q_20_W = 121086.033950196


def main() -> int:
    source_summary_path = SOURCE_DIR / "run_summary.json"
    source_restart = SOURCE_DIR / "final_restart.npz"
    if not source_summary_path.is_file() or not source_restart.is_file():
        raise FileNotFoundError(f"completed phase 3 restart is required: {SOURCE_DIR}")
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"output already exists: {OUTPUT_DIR}")

    source = json.loads(source_summary_path.read_text(encoding="utf-8"))
    if source.get("stop_reason") != "completed":
        raise RuntimeError(f"source phase 3 did not complete: {source}")
    current_a = float(source["initial_current_A"])
    flow_kg_s = float(source["final_flow_setpoint_kg_s"])
    result = run_low_power_case(LowPowerRunConfig(
        restart_in=source_restart,
        output_dir=OUTPUT_DIR,
        duration_s=EXTERNAL_HEAT_PERIOD_S,
        dt_s=0.05,
        tec_update_interval_s=0.05,
        record_interval_s=50.0,
        ramp_shape="quintic",
        final_power_w=Q_20_W,
        final_flow_kg_s=flow_kg_s,
        checkpoint_interval_s=0.0,
        fixed_current_a=current_a,
        initial_power_w=Q_20_W,
        initial_flow_kg_s=flow_kg_s,
        external_heat_enabled=True,
    ))
    if result.get("stop_reason") != "completed":
        raise RuntimeError(f"additional period did not complete: {result}")

    (OUTPUT_DIR / "additional_period_summary.json").write_text(
        json.dumps({
            "source_run": str(SOURCE_DIR),
            "duration_s": EXTERNAL_HEAT_PERIOD_S,
            "power_setpoint_W": Q_20_W,
            "flow_setpoint_kg_s": flow_kg_s,
            "history_interval_s": 50.0,
            "checkpoint_interval_s": 0.0,
            "external_heat_enabled": True,
            "run": result,
        }, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
