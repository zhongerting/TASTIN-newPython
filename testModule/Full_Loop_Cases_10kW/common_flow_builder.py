from typing import Any, Dict, List

from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
    PumpJunction,
)

from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig


class FlowControlledPumpJunction(PumpJunction):
    def __init__(self, *args: Any, target_flow_kg_s: float = 0.0, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self.W = float(target_flow_kg_s)
        self.target_W = float(target_flow_kg_s)

    def set_flow_rate(self, flow_kg_s: float) -> None:
        self.target_W = float(flow_kg_s)


def extend_channel_objects(volumes: List[Any], junctions: List[Any], channel: IncompressibleFluidChannel) -> None:
    volumes.extend(channel.volumes)
    junctions.extend(channel.internal_junctions)


def make_channel(
    *,
    name: str,
    n_nodes: int,
    length: float,
    area: float,
    dh: float,
    material: Any,
    initial_p: float,
    initial_t: float,
) -> IncompressibleFluidChannel:
    return IncompressibleFluidChannel(
        name=name,
        n_nodes=max(1, int(n_nodes)),
        total_length=float(length),
        flow_area=float(area),
        hydraulic_diam=float(dh),
        initial_P=float(initial_p),
        initial_T=float(initial_t),
        material=material,
    )


def _set_initial_flow(junctions: List[Any], value: float) -> None:
    for junc in junctions:
        junc.W = float(value)
        if hasattr(junc, "set_flow_rate"):
            junc.set_flow_rate(float(value))
        elif hasattr(junc, "target_W"):
            junc.target_W = float(value)
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()


def build_common_flow_objects(
    core_config: FullLoopCoreConfig,
    flow_config: FullLoopFlowConfig,
    pump_config: FullLoopPumpConfig,
    *,
    material: Any,
    close_with_placeholder_bridge: bool,
    connect_pump_outlet_to_core: bool = True,
) -> Dict[str, Any]:
    if float(flow_config.total_flow_kg_s) <= 0.0:
        raise ValueError("total_flow_kg_s must be positive.")
    if float(pump_config.pump_total_head_pa) <= 0.0:
        raise ValueError("pump_total_head_pa must be positive.")
    if int(pump_config.pump_count) != 2:
        raise ValueError("Full loop common builder currently requires two series pumps.")

    initial_p = float(core_config.reference_pressure_pa)
    initial_t = float(core_config.inlet_temperature_k)

    core_inlet = IncompressibleFluidVolume(
        name="CoreInletConnector",
        volume=float(flow_config.core_connector_volume_m3),
        length=float(flow_config.core_connector_length_m),
        flow_area=float(flow_config.core_inlet_segment_area_m2),
        hydraulic_diam=float(flow_config.core_inlet_segment_dh_m),
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    core_outlet = IncompressibleFluidVolume(
        name="CoreOutletConnector",
        volume=float(flow_config.core_connector_volume_m3),
        length=float(flow_config.core_connector_length_m),
        flow_area=float(flow_config.radiator_header_area_m2),
        hydraulic_diam=float(flow_config.radiator_header_dh_m),
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    core_inlet.is_pressure_boundary = False
    core_inlet.is_pressure_reference = True
    core_inlet.target_P = initial_p

    core_inlet_segment = None
    if bool(connect_pump_outlet_to_core):
        core_inlet_segment = make_channel(
            name="CoreInletSegment",
            n_nodes=int(flow_config.core_inlet_segment_n_nodes),
            length=float(flow_config.core_inlet_segment_length_m),
            area=float(flow_config.core_inlet_segment_area_m2),
            dh=float(flow_config.core_inlet_segment_dh_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
    radiator_inlet_header = make_channel(
        name="RadiatorInletHeader",
        n_nodes=int(flow_config.radiator_header_n_nodes),
        length=float(flow_config.radiator_inlet_header_length_m),
        area=float(flow_config.radiator_header_area_m2),
        dh=float(flow_config.radiator_header_dh_m),
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    radiator_outlet_header = make_channel(
        name="RadiatorOutletHeader",
        n_nodes=int(flow_config.radiator_header_n_nodes),
        length=float(flow_config.radiator_outlet_header_length_m),
        area=float(flow_config.radiator_header_area_m2),
        dh=float(flow_config.radiator_header_dh_m),
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )

    pump_area = float(flow_config.radiator_header_area_m2 if pump_config.pump_area_m2 is None else pump_config.pump_area_m2)
    pump_dh = float(flow_config.radiator_header_dh_m if pump_config.pump_dh_m is None else pump_config.pump_dh_m)
    pump_mid = IncompressibleFluidVolume(
        name="PumpMidNode",
        volume=float(pump_config.pump_node_volume_m3),
        length=float(pump_config.pump_node_length_m),
        flow_area=pump_area,
        hydraulic_diam=pump_dh,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    pump_outlet = IncompressibleFluidVolume(
        name="PumpOutletNode",
        volume=float(pump_config.pump_node_volume_m3),
        length=float(pump_config.pump_node_length_m),
        flow_area=pump_area,
        hydraulic_diam=pump_dh,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )

    volumes: List[Any] = [core_inlet, core_outlet, pump_mid, pump_outlet]
    junctions: List[Any] = []
    for channel in ([core_inlet_segment] if core_inlet_segment is not None else []) + [radiator_inlet_header, radiator_outlet_header]:
        extend_channel_objects(volumes, junctions, channel)

    j_core_to_radiator = FlowJunction(
        name="J_CoreOutletConnector_to_RadiatorInletHeader",
        from_vol=core_outlet,
        to_vol=radiator_inlet_header.volumes[0],
        flow_area=float(flow_config.radiator_header_area_m2),
        k_loss=float(flow_config.connector_k_loss),
        hydraulic_diam=float(flow_config.radiator_header_dh_m),
    )
    pump_single_head = 0.5 * float(pump_config.pump_total_head_pa)
    pump_cls = FlowControlledPumpJunction if bool(pump_config.pump_flow_control) else PumpJunction
    pump_target_flow = (
        float(flow_config.total_flow_kg_s)
        if pump_config.target_flow_kg_s is None
        else float(pump_config.target_flow_kg_s)
    )
    pump_kwargs = {"target_flow_kg_s": pump_target_flow} if bool(pump_config.pump_flow_control) else {}
    pump_a = pump_cls(
        name="J_PumpA",
        from_vol=radiator_outlet_header.volumes[-1],
        to_vol=pump_mid,
        flow_area=pump_area,
        k_loss=float(flow_config.connector_k_loss),
        delta_p=pump_single_head,
        **pump_kwargs,
    )
    pump_b = pump_cls(
        name="J_PumpB",
        from_vol=pump_mid,
        to_vol=pump_outlet,
        flow_area=pump_area,
        k_loss=float(flow_config.connector_k_loss),
        delta_p=pump_single_head,
        **pump_kwargs,
    )
    j_pump_to_core_inlet_segment = None
    j_core_inlet_segment_to_core = None
    if core_inlet_segment is not None:
        j_pump_to_core_inlet_segment = FlowJunction(
            name="J_PumpOutletNode_to_CoreInletSegment",
            from_vol=pump_outlet,
            to_vol=core_inlet_segment.volumes[0],
            flow_area=float(flow_config.core_inlet_segment_area_m2),
            k_loss=float(flow_config.connector_k_loss),
            hydraulic_diam=float(flow_config.core_inlet_segment_dh_m),
        )
        j_core_inlet_segment_to_core = FlowJunction(
            name="J_CoreInletSegment_to_CoreInletConnector",
            from_vol=core_inlet_segment.volumes[-1],
            to_vol=core_inlet,
            flow_area=float(flow_config.core_inlet_segment_area_m2),
            k_loss=float(flow_config.connector_k_loss),
            hydraulic_diam=float(flow_config.core_inlet_segment_dh_m),
        )

    junctions.append(j_core_to_radiator)
    if close_with_placeholder_bridge:
        j_bridge = FlowJunction(
            name="J_RadiatorPlaceholderBridge",
            from_vol=radiator_inlet_header.volumes[-1],
            to_vol=radiator_outlet_header.volumes[0],
            flow_area=float(flow_config.radiator_header_area_m2),
            k_loss=float(flow_config.radiator_bridge_k_loss),
            custom_length=float(flow_config.placeholder_bridge_length_m),
            hydraulic_diam=float(flow_config.radiator_header_dh_m),
        )
        junctions.append(j_bridge)
    else:
        j_bridge = None
    junctions.extend([pump_a, pump_b])
    if j_pump_to_core_inlet_segment is not None and j_core_inlet_segment_to_core is not None:
        junctions.extend([j_pump_to_core_inlet_segment, j_core_inlet_segment_to_core])

    _set_initial_flow(junctions, float(flow_config.total_flow_kg_s))
    if bool(pump_config.pump_flow_control):
        pump_a.set_flow_rate(pump_target_flow)
        pump_b.set_flow_rate(pump_target_flow)

    return {
        "volumes": volumes,
        "junctions": junctions,
        "core_inlet_connector": core_inlet,
        "core_outlet_connector": core_outlet,
        "core_inlet_segment": core_inlet_segment,
        "radiator_inlet_header": radiator_inlet_header,
        "radiator_outlet_header": radiator_outlet_header,
        "pump_mid_node": pump_mid,
        "pump_outlet_node": pump_outlet,
        "pump_a": pump_a,
        "pump_b": pump_b,
        "j_core_outlet_to_radiator_inlet": j_core_to_radiator,
        "j_radiator_outlet_to_pump_a": pump_a,
        "j_pump_outlet_to_core_inlet_segment": j_pump_to_core_inlet_segment,
        "j_core_inlet_segment_to_core_inlet": j_core_inlet_segment_to_core,
        "j_pump_b_to_core_inlet": j_core_inlet_segment_to_core,
        "placeholder_bridge_junction": j_bridge,
        "pump_total_head_pa": float(pump_config.pump_total_head_pa),
        "pump_single_head_pa": pump_single_head,
        "pump_target_flow_kg_s": pump_target_flow if bool(pump_config.pump_flow_control) else float(flow_config.total_flow_kg_s),
    }
