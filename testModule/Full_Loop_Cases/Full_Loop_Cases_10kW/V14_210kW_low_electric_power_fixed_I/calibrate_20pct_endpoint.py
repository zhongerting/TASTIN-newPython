"""Calibrate the 20% electric-power endpoint at one external-heat period."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

from .run_v14_low_power_fixed_i import LowPowerRunConfig, run_low_power_case


CASE_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = CASE_DIR / "runs"
SOURCE_20PCT_DIR = (
    CASE_DIR.parent / "V14_20pct_electric_power_fixed_I" / "runs"
    / "noext_iter2_Q120000_1500s"
)
FULL_POWER_RUN = (
    CASE_DIR / "runs"
    / "final_quintic1000_Q066_W100_tec005_full_orbit_record1s_20260813"
)
EXTERNAL_HEAT_PERIOD_S = 5668.144369
CANDIDATE_POWER_W = 121000.0


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-power", type=float, default=CANDIDATE_POWER_W)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    candidate_power_w = float(args.final_power)
    output_dir = args.output_dir or (
        OUTPUT_ROOT
        / f"endpoint_calibration_from20pct_Q{int(round(candidate_power_w))}_one_orbit_20260815"
    )
    source_summary_path = SOURCE_20PCT_DIR / "run_summary.json"
    source_restart = SOURCE_20PCT_DIR / "final_restart.npz"
    if not source_summary_path.is_file() or not source_restart.is_file():
        raise FileNotFoundError("formal 20% endpoint restart is required")
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")

    source = _read_json(source_summary_path)
    if source.get("stop_reason") != "completed":
        raise RuntimeError(f"formal 20% endpoint did not complete: {source}")
    full_power = _read_json(FULL_POWER_RUN / "run_summary.json")
    full_power_electric_w = float(full_power["initial_electric_power_W"])
    target_electric_w = 0.2 * full_power_electric_w
    current_a = float(source["initial_current_A"])
    flow_kg_s = float(source["final_flow_setpoint_kg_s"])

    result = run_low_power_case(LowPowerRunConfig(
        restart_in=source_restart,
        output_dir=output_dir,
        duration_s=EXTERNAL_HEAT_PERIOD_S,
        dt_s=0.05,
        tec_update_interval_s=0.05,
        record_interval_s=EXTERNAL_HEAT_PERIOD_S + 1.0,
        ramp_shape="quintic",
        final_power_w=candidate_power_w,
        final_flow_kg_s=flow_kg_s,
        checkpoint_interval_s=0.0,
        fixed_current_a=current_a,
        initial_power_w=candidate_power_w,
        initial_flow_kg_s=flow_kg_s,
        external_heat_enabled=True,
        enforce_outlet_limit=False,
    ))
    latest = result.get("latest_metrics", {})
    final_electric_w = float(latest["tec_main_electric_power_W"])
    (output_dir / "endpoint_calibration_summary.json").write_text(
        json.dumps({
            "source_20pct_endpoint": str(SOURCE_20PCT_DIR),
            "candidate_thermal_power_W": candidate_power_w,
            "full_power_electric_W": full_power_electric_w,
            "target_electric_power_W": target_electric_w,
            "final_electric_power_W": final_electric_w,
            "final_electric_ratio_to_full": final_electric_w / full_power_electric_w,
            "power_error_W": final_electric_w - target_electric_w,
            "run": result,
        }, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    if result.get("stop_reason") != "completed":
        raise RuntimeError(f"endpoint calibration did not complete: {result}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
