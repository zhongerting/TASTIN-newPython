"""Run the complete full-power -> 40% -> 20% fixed-I trajectory."""

from __future__ import annotations

import json
from pathlib import Path

from .run_v14_low_power_fixed_i import LowPowerRunConfig, run_low_power_case


CASE_DIR = Path(__file__).resolve().parent
SOURCE_40_RUN = (
    CASE_DIR / "runs" /
    "final_quintic1000_Q066_W100_tec005_full_orbit_record1s_20260813"
)
OUTPUT_ROOT = CASE_DIR / "runs" / "continuation_40pct_period10s_then_40to20_record1s_Q121086_20260815"
PHASE_1_DIR = OUTPUT_ROOT / "phase1_40pct_period_10s"
PHASE_2_DIR = OUTPUT_ROOT / "phase2_40to20_quintic1500_record1s"
EXTERNAL_HEAT_PERIOD_S = 5668.144369
Q_40_W = 138600.0
Q_20_W = 121086.033950196
RAMP_40_TO_20_S = 1500.0


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _require_completed(result: dict, label: str) -> None:
    if result.get("stop_reason") != "completed":
        raise RuntimeError(f"{label} did not complete: {result}")


def main() -> int:
    source_summary_path = SOURCE_40_RUN / "run_summary.json"
    source_restart = SOURCE_40_RUN / "final_restart.npz"
    if not source_summary_path.is_file() or not source_restart.is_file():
        raise FileNotFoundError(
            f"40% source result is incomplete: {SOURCE_40_RUN}"
        )
    if OUTPUT_ROOT.exists():
        raise FileExistsError(f"output already exists: {OUTPUT_ROOT}")

    source_summary = _read_json(source_summary_path)
    _require_completed(source_summary, "40% source run")
    current_a = float(source_summary["initial_current_A"])
    flow_kg_s = float(source_summary["final_flow_setpoint_kg_s"])

    phase1 = run_low_power_case(LowPowerRunConfig(
        restart_in=source_restart,
        output_dir=PHASE_1_DIR,
        duration_s=EXTERNAL_HEAT_PERIOD_S,
        dt_s=0.05,
        tec_update_interval_s=0.05,
        record_interval_s=10.0,
        ramp_shape="quintic",
        final_power_w=Q_40_W,
        final_flow_kg_s=flow_kg_s,
        checkpoint_interval_s=60.0,
        fixed_current_a=current_a,
        initial_power_w=Q_40_W,
        initial_flow_kg_s=flow_kg_s,
        external_heat_enabled=True,
    ))
    _require_completed(phase1, "40% full-period hold")

    phase1_restart = Path(phase1["restart_path"])
    phase2 = run_low_power_case(LowPowerRunConfig(
        restart_in=phase1_restart,
        output_dir=PHASE_2_DIR,
        duration_s=RAMP_40_TO_20_S,
        dt_s=0.05,
        tec_update_interval_s=0.05,
        record_interval_s=1.0,
        ramp_duration_s=RAMP_40_TO_20_S,
        ramp_shape="quintic",
        final_power_w=Q_20_W,
        final_flow_kg_s=flow_kg_s,
        checkpoint_interval_s=60.0,
        fixed_current_a=current_a,
        initial_power_w=Q_40_W,
        initial_flow_kg_s=flow_kg_s,
        external_heat_enabled=True,
    ))
    _require_completed(phase2, "40% to 20% descent")

    (OUTPUT_ROOT / "continuation_summary.json").write_text(
        json.dumps({
            "source_40_run": str(SOURCE_40_RUN),
            "current_A": current_a,
            "flow_kg_s": flow_kg_s,
            "phase1_40pct_period": phase1,
            "phase2_40to20_descent": phase2,
        }, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
