"""V15 pipe-fin radiator TOPAZ-II full-loop case builder."""

from typing import Any, Dict, Optional

from .common_builder import build_full_loop_common_base
from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig
from .v15_pipefin_radiator import V15PipeFinRadiatorConfig, attach_v15_pipefin_radiator


def build_v15_case_a_system(
    core_config: Optional[FullLoopCoreConfig] = None,
    flow_config: Optional[FullLoopFlowConfig] = None,
    pump_config: Optional[FullLoopPumpConfig] = None,
    radiator_config: Optional[V15PipeFinRadiatorConfig] = None,
) -> Dict[str, Any]:
    """Build V15 CaseA with common core/pump layer and local pipe-fin radiator."""
    core_cfg = FullLoopCoreConfig() if core_config is None else core_config
    flow_cfg = FullLoopFlowConfig() if flow_config is None else flow_config
    pump_cfg = FullLoopPumpConfig() if pump_config is None else pump_config
    rad_cfg = V15PipeFinRadiatorConfig() if radiator_config is None else radiator_config

    build = build_full_loop_common_base(
        core_config=core_cfg,
        flow_config=flow_cfg,
        pump_config=pump_cfg,
        radiator_connector=lambda common_build: attach_v15_pipefin_radiator(common_build, rad_cfg),
        connect_pump_outlet_to_core=False,
    )
    build["case_version"] = "v15_pipefin_radiator_full_loop"
    build["radiator_config"] = rad_cfg
    return build


def v15_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, float]:
    """Return compact topology/flow diagnostics for V15 checks."""
    return {
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "radiator_tube_count": int(len(build.get("radiator_tube_channels", []))),
        "radiator_unit_count": int(len(build.get("radiator_units", []))),
        "upper_header_count": int(len(build.get("radiator_upper_headers", []))),
        "lower_header_count": int(len(build.get("radiator_lower_headers", []))),
        "cold_return_branch_count": int(len(build.get("cold_return_branches", []))),
        "single_radiator_tube_flow_design_kg_s": float(build["single_radiator_tube_flow_design_kg_s"]),
        "cold_return_branch_flow_design_kg_s": float(build["cold_return_branch_flow_design_kg_s"]),
    }
