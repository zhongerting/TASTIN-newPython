"""Hydraulic-only V14 flow-path smoke run.

This case intentionally does not initialize thermal couplers or advance solids.
It only checks that the V14 closed hydraulic loop is connected, has a single
passive pressure reference, and can complete one hydraulic solve step.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    build_v14_case_a_system,
)


DEFAULT_OUTPUT_PATH = Path(__file__).with_name("v14_flow_path_smoke_result.json")


def _selected_junction_flows(build: Dict[str, Any]) -> Dict[str, float]:
    names = {
        "J_RadiatorInletHeader_to_HotOutletBranch_1",
        "J_HotOutletBranch_1_to_InletMix_I1",
        "J_I1_to_A1_I1_to_O1",
        "J_A1_I1_to_O1_to_O1",
        "J_OutletMix_O1_Manifold_1",
        "J_Manifold_1_to_RadiatorOutletHeader",
        "J_PumpA",
        "J_PumpB",
        "J_PumpOutletNode_to_CoreInletSegment",
        "J_CoreInletSegment_to_CoreInletConnector",
    }
    result = {}
    for junction in build["system"].fluid_solver.junctions_obj:
        name = getattr(junction, "name", "")
        if name in names:
            result[name] = float(getattr(junction, "W", np.nan))
    return result


def _selected_volume_pressures(build: Dict[str, Any]) -> Dict[str, float]:
    names = {
        "CoreInletConnector",
        "CoreOutletConnector",
        "RadiatorInletHeader_Vol_01",
        "HotOutletBranch_1_Vol_01",
        "InletMix_I1",
        "OutletMix_O1",
        "Manifold_1_Vol_01",
        "RadiatorOutletHeader_Vol_01",
        "PumpMidNode",
        "PumpOutletNode",
        "CoreInletSegment_Vol_01",
    }
    result = {}
    for volume in build["system"].fluid_solver.volumes_obj:
        name = getattr(volume, "name", "")
        if name in names:
            result[name] = float(getattr(volume, "P", np.nan))
    return result


def _nonfinite_names(objects, values) -> list[str]:
    mask = ~np.isfinite(np.asarray(values, dtype=float))
    return [getattr(obj, "name", f"item_{idx}") for idx, obj in enumerate(objects) if bool(mask[idx])]


def run_smoke_case(
    *,
    output_path: Optional[Path] = None,
    hydraulic_init_dt_s: float = 0.01,
    hydraulic_step_dt_s: float = 1.0e-4,
    hydraulic_tol_kg_s: float = 1.0e-4,
    hydraulic_max_iter: int = 1000,
) -> Dict[str, Any]:
    output = DEFAULT_OUTPUT_PATH if output_path is None else Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(main_tec_enabled=False),
        flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
        pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
        radiator_config=V14HeatPipeRadiatorConfig(
            hot_branch_n_nodes=1,
            manifold_node_counts=(1, 1, 1),
            hp_n_con=2,
            n_fin_height=3,
            ring_emissivity=0.0,
            hp_emissivity=0.0,
            fin_emissivity=0.0,
        ),
    )
    net = build["system"].fluid_solver

    pressure_reference_volumes = [
        getattr(volume, "name", "")
        for volume in net.volumes_obj
        if bool(getattr(volume, "is_pressure_reference", False))
    ]
    fixed_pressure_boundary_volumes = [
        getattr(volume, "name", "")
        for volume in net.volumes_obj
        if bool(getattr(volume, "is_pressure_boundary", False))
    ]

    hydraulic_init_converged = bool(
        net.initialize_hydraulics(
            dt=float(hydraulic_init_dt_s),
            tol=float(hydraulic_tol_kg_s),
            max_iter=int(hydraulic_max_iter),
        )
    )
    net.step_hydraulic(float(hydraulic_step_dt_s))

    result: Dict[str, Any] = {
        "case": "V14_flow_path_hydraulic_only_smoke",
        "case_version": str(build["case_version"]),
        "description": "Hydraulic-only closed-loop V14 smoke. Thermal radiation emissivities are set to zero and no SystemManager thermal coupling step is run.",
        "hydraulic_init_converged": hydraulic_init_converged,
        "hydraulic_step_completed": True,
        "hydraulic_init_dt_s": float(hydraulic_init_dt_s),
        "hydraulic_step_dt_s": float(hydraulic_step_dt_s),
        "hydraulic_tol_kg_s": float(hydraulic_tol_kg_s),
        "hydraulic_max_iter": int(hydraulic_max_iter),
        "n_volumes": int(net.n_vol),
        "n_junctions": int(net.n_junc),
        "pressure_reference_volumes": pressure_reference_volumes,
        "fixed_pressure_boundary_volumes": fixed_pressure_boundary_volumes,
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "pump_single_head_pa": float(build["pump_single_head_pa"]),
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "single_ring_branch_flow_design_kg_s": float(build["single_ring_branch_flow_design_kg_s"]),
        "min_pressure_pa_after_step": float(np.min(net.P_vec)),
        "max_pressure_pa_after_step": float(np.max(net.P_vec)),
        "max_abs_flow_kg_s_after_step": float(np.max(np.abs(net.W_vec))),
        "min_flow_kg_s_after_step": float(np.min(net.W_vec)),
        "max_flow_kg_s_after_step": float(np.max(net.W_vec)),
        "nonfinite_flow_junctions": _nonfinite_names(net.junctions_obj, net.W_vec),
        "nonfinite_pressure_volumes": _nonfinite_names(net.volumes_obj, net.P_vec),
        "selected_junction_flows_kg_s": _selected_junction_flows(build),
        "selected_volume_pressures_pa": _selected_volume_pressures(build),
    }

    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    saved = run_smoke_case()
    print(json.dumps(saved, indent=2, sort_keys=True, ensure_ascii=False))
