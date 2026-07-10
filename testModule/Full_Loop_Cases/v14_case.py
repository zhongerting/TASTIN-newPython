"""V14 heat-pipe radiator TOPAZ-II full-loop case builder."""

from typing import Any, Dict, Optional

from .common_builder import build_full_loop_common_base
from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig
from .v14_heatpipe_radiator import V14HeatPipeRadiatorConfig, attach_v14_heatpipe_radiator


def build_v14_case_a_system(
    core_config: Optional[FullLoopCoreConfig] = None,
    flow_config: Optional[FullLoopFlowConfig] = None,
    pump_config: Optional[FullLoopPumpConfig] = None,
    radiator_config: Optional[V14HeatPipeRadiatorConfig] = None,
) -> Dict[str, Any]:
    """Build V14 CaseA with the common core/pump layer and local heat-pipe radiator."""
    core_cfg = FullLoopCoreConfig() if core_config is None else core_config
    flow_cfg = FullLoopFlowConfig() if flow_config is None else flow_config
    pump_cfg = FullLoopPumpConfig() if pump_config is None else pump_config
    rad_cfg = V14HeatPipeRadiatorConfig() if radiator_config is None else radiator_config

    build = build_full_loop_common_base(
        core_config=core_cfg,
        flow_config=flow_cfg,
        pump_config=pump_cfg,
        radiator_connector=lambda common_build: attach_v14_heatpipe_radiator(common_build, rad_cfg),
    )
    build["case_version"] = "v14_heatpipe_radiator_full_loop"
    build["radiator_config"] = rad_cfg
    return build


def v14_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, float]:
    """Return compact topology/flow diagnostics for V14 smoke checks."""
    hot_to_ring = build.get("hot_outlet_to_ring_junctions", [])
    manifold_to_header = build.get("manifold_to_outlet_header_junctions", [])
    return {
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "single_ring_in_total_kg_s": float(sum(getattr(j, "W", 0.0) for j in hot_to_ring)),
        "single_ring_out_total_kg_s": float(sum(getattr(j, "W", 0.0) for j in manifold_to_header)),
        "ring_sector_count": float(len(build.get("ring_sectors", []))),
        "ring_hp_count": float(len(build.get("ring_hps", []))),
    }
