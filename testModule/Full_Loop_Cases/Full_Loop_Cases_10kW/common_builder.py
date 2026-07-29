from typing import Any, Callable, Dict, Optional

from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager

from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig
from .common_core_builder import build_full_loop_core, make_coolant
from .common_flow_builder import build_common_flow_objects


def _sync_design_flows(build: Dict[str, Any]) -> None:
    total_flow = float(build["total_flow_design_kg_s"])
    tfe_flow = float(build["single_tfe_flow_design_kg_s"])
    net = build["system"].fluid_solver
    for idx, junc in enumerate(net.junctions_obj):
        name = getattr(junc, "name", "")
        if hasattr(junc, "design_flow_kg_s"):
            value = float(junc.design_flow_kg_s)
        else:
            value = tfe_flow if name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_") or name.startswith("Chan_") else total_flow
        junc.W = value
        target_value = (
            float(build.get("pump_target_flow_kg_s", value))
            if bool(getattr(junc, "is_pump_junction", False)) and hasattr(junc, "target_W")
            else value
        )
        if hasattr(junc, "set_flow_rate"):
            junc.set_flow_rate(target_value)
        elif hasattr(junc, "target_W"):
            junc.target_W = target_value
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()
        net.W_vec[idx] = value
    if hasattr(net, "W_old"):
        net.W_old[:] = net.W_vec
    if hasattr(net, "W_iterate"):
        net.W_iterate[:] = net.W_vec
    net._refresh_cached_boundary_targets()


def build_full_loop_common_base(
    core_config: FullLoopCoreConfig,
    flow_config: FullLoopFlowConfig,
    pump_config: FullLoopPumpConfig,
    radiator_connector: Optional[Callable[[Dict[str, Any]], None]] = None,
    close_with_placeholder_bridge: bool = False,
    connect_pump_outlet_to_core: bool = True,
) -> Dict[str, Any]:
    if radiator_connector is not None and close_with_placeholder_bridge:
        raise ValueError("Use either radiator_connector or close_with_placeholder_bridge, not both.")
    if radiator_connector is None and not close_with_placeholder_bridge:
        raise ValueError("A radiator_connector or close_with_placeholder_bridge=True is required to close the loop.")

    coolant, coolant_name = make_coolant(core_config.coolant_material)
    flow = build_common_flow_objects(
        core_config,
        flow_config,
        pump_config,
        material=coolant,
        close_with_placeholder_bridge=close_with_placeholder_bridge,
        connect_pump_outlet_to_core=connect_pump_outlet_to_core,
    )
    core_build = build_full_loop_core(
        core_config,
        core_inlet_connector=flow["core_inlet_connector"],
        core_outlet_connector=flow["core_outlet_connector"],
        total_flow_kg_s=float(flow_config.total_flow_kg_s),
    )

    volumes = list(flow["volumes"]) + list(core_build["core_fluid_volumes"])
    junctions = list(flow["junctions"]) + list(core_build["core_fluid_junctions"])
    build: Dict[str, Any] = {
        "case_version": "full_loop_common_base",
        "system": None,
        "core": core_build["core"],
        "tfes": core_build["tfes"],
        "fluid_channels": core_build["fluid_channels"],
        "ring_multipliers": core_build["ring_multipliers"],
        "tec_ring_multipliers": core_build["tec_ring_multipliers"],
        "coolant_material": coolant_name,
        "total_flow_design_kg_s": float(flow_config.total_flow_kg_s),
        "single_tfe_flow_design_kg_s": core_build["single_tfe_flow_design_kg_s"],
        **{key: value for key, value in flow.items() if key not in {"volumes", "junctions"}},
    }
    if radiator_connector is not None:
        extra = radiator_connector(build)
        if extra is not None:
            raise ValueError("radiator_connector must mutate build in place and return None.")
        volumes.extend(build.get("radiator_adapter_volumes", []))
        junctions.extend(build.get("radiator_adapter_junctions", []))

    hydraulic_net = HydraulicNetwork(volumes=volumes, junctions=junctions, gravity_vector=0.0)
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core_build["core"])
    for component in build.get("radiator_adapter_components", []):
        system.add_component(component)
    build["system"] = system
    _sync_design_flows(build)
    return build
