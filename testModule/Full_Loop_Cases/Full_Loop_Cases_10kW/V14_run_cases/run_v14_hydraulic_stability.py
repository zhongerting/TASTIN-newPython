"""Hydraulic-only V14_10kW stability run."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np


REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    build_v14_case_a_system,
)


DEFAULT_OUTPUT_PATH = Path(__file__).with_name("v14_hydraulic_stability_result.json")


def _nonfinite_names(objects, values) -> list[str]:
    mask = ~np.isfinite(np.asarray(values, dtype=float))
    return [getattr(obj, "name", f"item_{idx}") for idx, obj in enumerate(objects) if bool(mask[idx])]


def _junction_flow(net: Any, name: str) -> float:
    for junction in net.junctions_obj:
        if getattr(junction, "name", "") == name:
            return float(getattr(junction, "W", np.nan))
    return float("nan")


def run_stability_case(
    *,
    output_path: Optional[Path] = None,
    n_steps: int = 200,
    last_window: int = 40,
    hydraulic_init_dt_s: float = 0.01,
    hydraulic_step_dt_s: float = 1.0e-4,
    hydraulic_tol_kg_s: float = 1.0e-4,
    hydraulic_max_iter: int = 1000,
    stable_flow_tol_kg_s: float = 1.0e-6,
) -> Dict[str, Any]:
    output = DEFAULT_OUTPUT_PATH if output_path is None else Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(main_tec_enabled=False),
        flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
        pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
        radiator_config=V14HeatPipeRadiatorConfig(
            ring_emissivity=0.0,
            hp_emissivity=0.0,
            fin_emissivity=0.0,
        ),
    )
    net = build["system"].fluid_solver

    hydraulic_init_converged = bool(
        net.initialize_hydraulics(
            dt=float(hydraulic_init_dt_s),
            tol=float(hydraulic_tol_kg_s),
            max_iter=int(hydraulic_max_iter),
        )
    )

    pump_a_history = []
    pump_b_history = []
    for _ in range(int(n_steps)):
        net.step_hydraulic(float(hydraulic_step_dt_s))
        pump_a_history.append(_junction_flow(net, "J_PumpA"))
        pump_b_history.append(_junction_flow(net, "J_PumpB"))

    window = max(2, min(int(last_window), len(pump_a_history)))
    pump_a_tail = np.asarray(pump_a_history[-window:], dtype=float)
    pump_b_tail = np.asarray(pump_b_history[-window:], dtype=float)
    last_window_max_change = float(np.max(np.abs(np.diff(pump_a_tail)))) if len(pump_a_tail) > 1 else float("nan")
    stable_flow_reached = bool(
        hydraulic_init_converged
        and np.all(np.isfinite(pump_a_tail))
        and np.all(np.isfinite(pump_b_tail))
        and last_window_max_change < float(stable_flow_tol_kg_s)
    )

    result: Dict[str, Any] = {
        "case": "V14_10kW_hydraulic_stability",
        "description": "Hydraulic-only V14_10kW run. No heat power, no TEC electrical calculation, no solid thermal step.",
        "case_version": str(build["case_version"]),
        "hydraulic_init_converged": hydraulic_init_converged,
        "stable_flow_reached": stable_flow_reached,
        "stable_flow_tol_kg_s": float(stable_flow_tol_kg_s),
        "n_steps": int(n_steps),
        "last_window": int(window),
        "hydraulic_step_dt_s": float(hydraulic_step_dt_s),
        "n_volumes": int(net.n_vol),
        "n_junctions": int(net.n_junc),
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "pump_single_head_pa": float(build["pump_single_head_pa"]),
        "design_total_flow_kg_s": float(build["total_flow_design_kg_s"]),
        "upper_ring_design_flow_kg_s": float(build["upper_ring_total_flow_design_kg_s"]),
        "lower_ring_design_flow_kg_s": float(build["lower_ring_total_flow_design_kg_s"]),
        "final_pump_a_flow_kg_s": float(pump_a_history[-1]),
        "final_pump_b_flow_kg_s": float(pump_b_history[-1]),
        "last_window_mean_pump_a_flow_kg_s": float(np.mean(pump_a_tail)),
        "last_window_max_abs_pump_flow_change_kg_s": last_window_max_change,
        "min_pressure_pa": float(np.min(net.P_vec)),
        "max_pressure_pa": float(np.max(net.P_vec)),
        "min_flow_kg_s": float(np.min(net.W_vec)),
        "max_flow_kg_s": float(np.max(net.W_vec)),
        "nonfinite_flow_junctions": _nonfinite_names(net.junctions_obj, net.W_vec),
        "nonfinite_pressure_volumes": _nonfinite_names(net.volumes_obj, net.P_vec),
    }

    output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")
    return result


if __name__ == "__main__":
    saved = run_stability_case()
    print(json.dumps(saved, indent=2, sort_keys=True, ensure_ascii=False))
