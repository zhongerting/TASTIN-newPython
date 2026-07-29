"""Run fixed 210 kW V14 with N18 external heat for two orbital periods."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
from typing import Optional, Sequence

import numpy as np

REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    run_debug_case,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    ReactivityControlRunConfig,
    load_baseline_debug_config,
)


CASE_DIR = Path(__file__).resolve().parent
CASE_NAME = "V14_10kW_210kW_fixed_power_external_heat_2orbits"
ORBIT_PERIOD_S = 5668.144369
TOTAL_DURATION_S = 2.0 * ORBIT_PERIOD_S
DEFAULT_RESTART = (
    CASE_DIR.parent
    / "V14_210kW_reactivity_control"
    / "反应性控制"
    / "checkpoint_t013864s.npz"
)
DEFAULT_OUTPUT_DIR = CASE_DIR / "runs" / "default"


def restart_time_s(path: Path) -> float:
    with np.load(path, allow_pickle=False) as data:
        return float(data["System/global_time"][0])


def build_run_config(args: argparse.Namespace):
    phase_origin = restart_time_s(args.restart_in)
    debug, _ = load_baseline_debug_config(ReactivityControlRunConfig(
        restart_in=args.restart_in,
        output_dir=args.output_dir,
        external_heat_enabled=True,
        external_heat_period_s=ORBIT_PERIOD_S,
        external_heat_time_origin_s=phase_origin,
    ))
    return replace(
        debug,
        output_dir=args.output_dir,
        stage_durations_s=(TOTAL_DURATION_S,),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        min_fluid_temperature_stop_k=(
            None if float(args.min_fluid_temperature_stop) <= 0.0
            else float(args.min_fluid_temperature_stop)
        ),
        power_w=210000.0,
        external_heat_enabled=True,
        external_heat_period_s=ORBIT_PERIOD_S,
        external_heat_time_origin_s=phase_origin,
        case_prefix=CASE_NAME,
    )


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--checkpoint-interval", type=float, default=600.0)
    parser.add_argument("--min-fluid-temperature-stop", type=float, default=500.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = run_debug_case(build_run_config(args))
    print(json.dumps(result["latest_metrics"], indent=2, sort_keys=True, ensure_ascii=False))
    print("Stop reason: {}".format(result["stop_reason"]))
    print("Saved outputs to: {}".format(result["output_dir"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
