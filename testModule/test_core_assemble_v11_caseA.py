import os
import sys
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_6segment_v9_interface as ring_cfg
from Components.BaseComponent import BaseComponent
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
    MacroFlowJunction,
    PumpJunction,
)
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager
from test_core_assemble_v8_caseA import (
    CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    V8_PHYSICAL_RING_COUNT,
    _case_a_electric_diagnostics,
    build_v8_case_a_system,
)
from test_core_assemble_v9_caseA import (
    AREA_CORE_BRANCH,
    AREA_HEADER,
    AREA_SMALL_BRANCH,
    DH_CORE_BRANCH,
    DH_HEADER,
    DH_SMALL_BRANCH,
    L_COLD_RETURN_1,
    L_COLD_RETURN_23,
    L_HOT_OUTLET,
    L_RADIATOR_INNER_HEADER_53,
    L_RADIATOR_OUTER_HEADER_52,
    _channel_delta_p,
    _channel_t_in_out,
    _extend_channel_objects,
    _make_channel,
    _rename_channel,
    _selected_old_core_junctions,
)


V11_CASE_VERSION = "v11_closed_core_collector_ring_pumped_loop"
V11_DEFAULT_INLET_TEMPERATURE_K = 753.330663091
V11_DEFAULT_REFERENCE_PRESSURE_PA = 166471.52
V11_DEFAULT_CONNECTOR_VOLUME_M3 = 1.0e-5
V11_DEFAULT_CONNECTOR_LENGTH_M = 0.02
V11_DEFAULT_PUMP_TOTAL_HEAD_PA = 6466.56
SIGMA_SB = 5.670374419e-8


class FluidChannelRadiationSink(BaseComponent):
    """Equivalent radiation sink applied directly to a fluid channel."""

    def __init__(
        self,
        name: str,
        channel: IncompressibleFluidChannel,
        emissivity: float,
        t_space: float,
        perimeter: float,
        area_scale: float = 1.0,
    ):
        super().__init__(name)
        self.channel = channel
        self.emissivity = float(emissivity)
        self.t_space = float(t_space)
        self.perimeter = float(perimeter)
        self.area_scale = float(area_scale)
        self.node_area = self.perimeter * float(channel.node_length) * self.area_scale
        self.last_radiation_distribution_w = np.zeros(channel.n_nodes, dtype=float)
        self.last_radiation_w = 0.0

    def pre_step(self, dt: float, current_time: float):
        temperatures = np.asarray(self.channel.temperature_vector, dtype=float)
        q_reject = (
            self.emissivity
            * SIGMA_SB
            * self.node_area
            * (np.maximum(temperatures, 1.0e-6) ** 4 - self.t_space ** 4)
        )
        q_reject = np.nan_to_num(q_reject, nan=0.0, posinf=0.0, neginf=0.0)
        self.last_radiation_distribution_w[:] = q_reject
        self.last_radiation_w = float(np.sum(q_reject))
        self.channel.add_heat_source_distribution(-q_reject)


def _is_core_channel_junction(name: str) -> bool:
    return (
        name.startswith("J_PlenumIn_")
        or name.startswith("J_PlenumOut_")
        or name.startswith("Chan_")
    )


def _is_hot_outlet_channel_junction(name: str) -> bool:
    return name.startswith("HotOutletBranch_") and "_Junc_" in name


def _is_ring_sector_channel_junction(name: str) -> bool:
    return name.startswith("A") and "_Channel_Junc_" in name


def _is_ring_sector_link_junction(name: str) -> bool:
    return (
        name.startswith("J_I")
        and "_to_A" in name
    ) or (
        name.startswith("J_A")
        and "_to_O" in name
    ) or (
        name.startswith("J_A")
        and "_to_I" in name
    )


def _is_manifold_channel_junction(name: str) -> bool:
    return name.startswith("Manifold_") and "_Junc_" in name


def _is_outlet_mix_to_manifold_junction(name: str) -> bool:
    return name.startswith("J_OutletMix_") and "_Manifold_" in name


def _is_cold_return_channel_junction(name: str) -> bool:
    return name.startswith("ColdReturnBranch_") and "_Junc_" in name


def _is_radiator_header_channel_junction(name: str) -> bool:
    return (
        name.startswith("RadiatorInnerHeader_53")
        or name.startswith("RadiatorOuterHeader_52")
    ) and "_Junc_" in name


def _is_preserved_ring_restart_junction(name: str) -> bool:
    return (
        _is_ring_sector_channel_junction(name)
        or _is_ring_sector_link_junction(name)
        or _is_manifold_channel_junction(name)
        or _is_outlet_mix_to_manifold_junction(name)
    )


def _build_ring_only(
    material: Any,
    initial_p: float,
    initial_t: float,
    ring_emissivity: Optional[float] = None,
    hp_emissivity: Optional[float] = None,
    fin_emissivity: Optional[float] = None,
) -> Dict[str, Any]:
    old_ring_emissivity = ring_cfg.RING_EMISSIVITY
    old_hp_emissivity = ring_cfg.cfg.HP_EMISSIVITY
    old_fin_emissivity = ring_cfg.cfg.FIN_EMISSIVITY
    if ring_emissivity is not None:
        ring_cfg.RING_EMISSIVITY = float(ring_emissivity)
    if hp_emissivity is not None:
        ring_cfg.cfg.HP_EMISSIVITY = float(hp_emissivity)
    if fin_emissivity is not None:
        ring_cfg.cfg.FIN_EMISSIVITY = float(fin_emissivity)

    inlet_mix_nodes = {
        key: ring_cfg.build_mix_node(f"InletMix_{key}", "inlet")
        for key in ring_cfg.INLET_MIX_KEYS
    }
    outlet_mix_nodes = {
        key: ring_cfg.build_mix_node(f"OutletMix_{key}", "outlet")
        for key in ring_cfg.OUTLET_MIX_KEYS
    }
    mix_nodes = {**inlet_mix_nodes, **outlet_mix_nodes}
    for vol in [*inlet_mix_nodes.values(), *outlet_mix_nodes.values()]:
        vol.material = material
        vol.P = float(initial_p)
        vol.T = float(initial_t)
        vol.h = material.enthalpy(float(initial_t))

    sectors = []
    solids = []
    ring_hps = []
    segment_links = []
    segment_entry_links = []
    segment_exit_links = []

    try:
        for sector_name, start_key, end_key, multipliers in ring_cfg.SEGMENT_SPECS:
            channel = IncompressibleFluidChannel(
                name=f"{sector_name}_Channel",
                n_nodes=ring_cfg.N_SECTOR,
                total_length=ring_cfg.L_SECTOR,
                flow_area=ring_cfg.AREA_RING,
                hydraulic_diam=ring_cfg.DH_RING,
                initial_P=float(initial_p),
                initial_T=float(initial_t),
                material=material,
            )
            solid = ring_cfg.build_sector_solid(f"{sector_name}_Solid")
            ring_hp = ring_cfg.cfg.build_ring_hp(
                name=f"{sector_name}_RingHP",
                fluid_channel=channel,
                solid_header=solid,
                hp_multipliers=multipliers,
            )
            ring_cfg.configure_ring_hp_heat_pipe_solver(ring_hp)
            sectors.append(channel)
            solids.append(solid)
            ring_hps.append(ring_hp)

            entry_link = FlowJunction(
                name=f"J_{start_key}_to_{sector_name}",
                from_vol=mix_nodes[start_key],
                to_vol=channel.volumes[0],
                flow_area=ring_cfg.AREA_RING,
                k_loss=ring_cfg.K_INLET_MIX_TO_RING_SEGMENT,
            )
            exit_link = FlowJunction(
                name=f"J_{sector_name}_to_{end_key}",
                from_vol=channel.volumes[-1],
                to_vol=mix_nodes[end_key],
                flow_area=ring_cfg.AREA_RING,
                k_loss=ring_hp.outlet_k_loss + ring_cfg.K_RING_SEGMENT_TO_OUTLET_MIX,
                dynamic_loss_params=ring_hp.outlet_dynamic_loss_params,
            )
            segment_entry_links.append(entry_link)
            segment_exit_links.append(exit_link)
            segment_links.extend([entry_link, exit_link])
    finally:
        ring_cfg.RING_EMISSIVITY = old_ring_emissivity
        ring_cfg.cfg.HP_EMISSIVITY = old_hp_emissivity
        ring_cfg.cfg.FIN_EMISSIVITY = old_fin_emissivity

    return {
        "inlet_mix_nodes": inlet_mix_nodes,
        "outlet_mix_nodes": outlet_mix_nodes,
        "sectors": sectors,
        "solids": solids,
        "ring_hps": ring_hps,
        "segment_links": segment_links,
        "segment_entry_links": segment_entry_links,
        "segment_exit_links": segment_exit_links,
    }


def build_v11_case_a_system(
    inlet_temperature_k: float = V11_DEFAULT_INLET_TEMPERATURE_K,
    total_inlet_flow_kg_s: float = CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    reference_pressure_pa: float = V11_DEFAULT_REFERENCE_PRESSURE_PA,
    pump_total_head_pa: float = V11_DEFAULT_PUMP_TOTAL_HEAD_PA,
    pipe_n_nodes: int = 8,
    external_pipe_n_nodes: int = 5,
    connector_volume_m3: float = V11_DEFAULT_CONNECTOR_VOLUME_M3,
    connector_length_m: float = V11_DEFAULT_CONNECTOR_LENGTH_M,
    solid_heat_capacity_scale: float = 1.0,
    solid_heat_capacity_scale_scope: str = "global_outer",
    global_outer_heat_capacity_scale: Optional[float] = None,
    ring_multipliers: Optional[Sequence[int]] = None,
    tec_ring_multipliers: Optional[Sequence[int]] = None,
    coolant_material: str = "SodiumPotassium78",
    ring_emissivity: Optional[float] = None,
    hp_emissivity: Optional[float] = None,
    fin_emissivity: Optional[float] = None,
    outer_header_emissivity: float = 0.0,
    outer_header_t_space_k: Optional[float] = None,
    outer_header_area_scale: float = 1.0,
) -> Dict[str, Any]:
    if total_inlet_flow_kg_s <= 0.0:
        raise ValueError("total_inlet_flow_kg_s must be positive.")

    build = build_v8_case_a_system(
        inlet_temperature_k=float(inlet_temperature_k),
        pipe_n_nodes=int(pipe_n_nodes),
        inlet_plenum_volume_m3=float(connector_volume_m3),
        outlet_plenum_volume_m3=float(connector_volume_m3),
        plenum_length_m=float(connector_length_m),
        solid_heat_capacity_scale=solid_heat_capacity_scale,
        solid_heat_capacity_scale_scope=solid_heat_capacity_scale_scope,
        global_outer_heat_capacity_scale=global_outer_heat_capacity_scale,
        coolant_material=coolant_material,
        ring_multipliers=ring_multipliers,
        tec_ring_multipliers=tec_ring_multipliers,
    )

    material = build["inlet_boundary"].material
    core = build["core"]
    tfes = build["tfes"]
    fluid_channels = build["fluid_channels"]
    core_inlet = build["inlet_plenum"]
    core_outlet = build["outlet_plenum"]
    core_inlet.name = "CoreInletConnector"
    core_outlet.name = "CoreOutletConnector"

    initial_p = float(reference_pressure_pa)
    initial_t = float(inlet_temperature_k)
    total_flow = float(total_inlet_flow_kg_s)
    pump_total_head = float(pump_total_head_pa)
    pump_single_head = 0.5 * pump_total_head
    macro_branch_flow = total_flow / 3.0
    single_ring_branch_flow = macro_branch_flow / 2.0
    total_multiplier = float(sum(build["ring_multipliers"].values()))
    tfe_single_flow = total_flow / total_multiplier

    core_inlet.is_pressure_boundary = False
    core_inlet.is_pressure_reference = True
    core_inlet.target_P = initial_p

    radiator_merge = IncompressibleFluidVolume(
        name="V10_RadiatorManifoldMerge",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    header_to_cold_split = IncompressibleFluidVolume(
        name="V11_PumpOutletDistributor_51",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    cold_return_merge = IncompressibleFluidVolume(
        name="V10_ColdReturnOutletMerge",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_CORE_BRANCH * 3.0,
        hydraulic_diam=DH_CORE_BRANCH,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    pump_mid_node = IncompressibleFluidVolume(
        name="V11_PumpMidNode",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
        material=material,
        initial_P=initial_p + pump_single_head,
        initial_T=initial_t,
    )

    _rename_channel(build["inlet_pipe_1"], "ColdReturnBranch_1")
    _rename_channel(build["inlet_pipe_23"], "ColdReturnBranch_2_3_Rep")

    radiator_inner_header_53 = _make_channel(
        name="RadiatorInnerHeader_53",
        n_nodes=max(1, int(round(external_pipe_n_nodes * L_RADIATOR_INNER_HEADER_53 / 0.40911))),
        length=L_RADIATOR_INNER_HEADER_53,
        area=AREA_HEADER,
        dh=DH_HEADER,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    radiator_outer_header_52 = _make_channel(
        name="RadiatorOuterHeader_52",
        n_nodes=1,
        length=L_RADIATOR_OUTER_HEADER_52,
        area=AREA_HEADER,
        dh=DH_HEADER,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    hot_outlet_branches = [
        _make_channel(
            name=f"HotOutletBranch_{idx}",
            n_nodes=int(pipe_n_nodes),
            length=L_HOT_OUTLET,
            area=AREA_CORE_BRANCH,
            dh=DH_CORE_BRANCH,
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        for idx in range(1, 4)
    ]
    manifolds = [
        _make_channel(
            name=f"Manifold_{idx}",
            n_nodes=ring_cfg.MANIFOLD_NODE_COUNTS[idx - 1],
            length=ring_cfg.MANIFOLD_LENGTHS[idx - 1],
            area=AREA_SMALL_BRANCH,
            dh=DH_SMALL_BRANCH,
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        for idx in range(1, 4)
    ]
    ring = _build_ring_only(
        material,
        initial_p,
        initial_t,
        ring_emissivity=ring_emissivity,
        hp_emissivity=hp_emissivity,
        fin_emissivity=fin_emissivity,
    )
    effective_ring_emissivity = float(ring_cfg.RING_EMISSIVITY if ring_emissivity is None else ring_emissivity)
    effective_hp_emissivity = float(ring_cfg.cfg.HP_EMISSIVITY if hp_emissivity is None else hp_emissivity)
    effective_fin_emissivity = float(ring_cfg.cfg.FIN_EMISSIVITY if fin_emissivity is None else fin_emissivity)
    effective_outer_header_t_space = float(ring_cfg.T_SPACE if outer_header_t_space_k is None else outer_header_t_space_k)
    outer_header_radiation_sink = None
    if float(outer_header_emissivity) > 0.0:
        outer_header_radiation_sink = FluidChannelRadiationSink(
            name="RadiatorOuterHeader_52_RadiationSink",
            channel=radiator_outer_header_52,
            emissivity=float(outer_header_emissivity),
            t_space=effective_outer_header_t_space,
            perimeter=np.pi * DH_HEADER,
            area_scale=float(outer_header_area_scale),
        )

    all_vols: List[Any] = [
        core_inlet,
        core_outlet,
        radiator_merge,
        pump_mid_node,
        header_to_cold_split,
        cold_return_merge,
    ]
    all_juncs: List[Any] = []
    for channel in fluid_channels.values():
        _extend_channel_objects(all_vols, all_juncs, channel)
    for channel in hot_outlet_branches:
        _extend_channel_objects(all_vols, all_juncs, channel)
    all_vols.extend(ring["inlet_mix_nodes"].values())
    all_vols.extend(ring["outlet_mix_nodes"].values())
    for channel in ring["sectors"]:
        _extend_channel_objects(all_vols, all_juncs, channel)
    for channel in manifolds:
        _extend_channel_objects(all_vols, all_juncs, channel)
    for channel in (radiator_inner_header_53, radiator_outer_header_52, build["inlet_pipe_1"], build["inlet_pipe_23"]):
        _extend_channel_objects(all_vols, all_juncs, channel)

    hot_outlet_entry_junctions = []
    hot_outlet_to_ring_junctions = []
    for idx, channel in enumerate(hot_outlet_branches, start=1):
        hot_outlet_entry_junctions.append(
            FlowJunction(
                name=f"J_CoreOutletConnector_to_HotOutletBranch_{idx}",
                from_vol=core_outlet,
                to_vol=channel.volumes[0],
                flow_area=AREA_CORE_BRANCH,
                k_loss=0.0,
            )
        )
        inlet_key = ring_cfg.INLET_MIX_KEYS[idx - 1]
        hot_outlet_to_ring_junctions.append(
            MacroFlowJunction(
                name=f"J_HotOutletBranch_{idx}_to_InletMix_{inlet_key}",
                from_vol=channel.volumes[-1],
                to_vol=ring["inlet_mix_nodes"][inlet_key],
                macro_vol=channel.volumes[-1],
                multiplier=2,
                flow_area=AREA_CORE_BRANCH,
                k_loss=ring_cfg.K_HOT_LEG_TO_INLET_MIX,
            )
        )

    outlet_mix_to_manifold = []
    manifold_to_merge = []
    for idx, key in enumerate(ring_cfg.OUTLET_MIX_KEYS, start=1):
        outlet_mix_to_manifold.append(
            FlowJunction(
                name=f"J_OutletMix_{key}_Manifold_{idx}",
                from_vol=ring["outlet_mix_nodes"][key],
                to_vol=manifolds[idx - 1].volumes[0],
                flow_area=AREA_SMALL_BRANCH,
                k_loss=ring_cfg.K_OUTLET_MIX_TO_MANIFOLD,
            )
        )
        manifold_to_merge.append(
            MacroFlowJunction(
                name=f"J_Manifold_{idx}_to_RadiatorMerge",
                from_vol=manifolds[idx - 1].volumes[-1],
                to_vol=radiator_merge,
                macro_vol=radiator_merge,
                multiplier=2,
                flow_area=AREA_SMALL_BRANCH,
                k_loss=ring_cfg.K_MANIFOLD_TO_OUTLET_HEADER,
            )
        )

    j_merge_to_h53 = FlowJunction(
        name="J_RadiatorMerge_to_RadiatorInnerHeader_53",
        from_vol=radiator_merge,
        to_vol=radiator_inner_header_53.volumes[0],
        flow_area=AREA_HEADER,
        k_loss=0.0,
    )
    j_h53_to_h52 = FlowJunction(
        name="J_RadiatorInnerHeader_53_to_RadiatorOuterHeader_52",
        from_vol=radiator_inner_header_53.volumes[-1],
        to_vol=radiator_outer_header_52.volumes[0],
        flow_area=AREA_HEADER,
        k_loss=0.0,
    )
    pump_a = PumpJunction(
        name="J_RadiatorOuterHeader_52_to_PumpA",
        from_vol=radiator_outer_header_52.volumes[-1],
        to_vol=pump_mid_node,
        flow_area=AREA_HEADER,
        k_loss=0.0,
        delta_p=pump_single_head,
    )
    pump_b = PumpJunction(
        name="J_PumpA_to_PumpOutletDistributor_51",
        from_vol=pump_mid_node,
        to_vol=header_to_cold_split,
        flow_area=AREA_HEADER,
        k_loss=0.0,
        delta_p=pump_single_head,
    )
    j_cold_1_in = FlowJunction(
        name="J_PumpOutletDistributor_51_to_ColdReturnBranch_1",
        from_vol=header_to_cold_split,
        to_vol=build["inlet_pipe_1"].volumes[0],
        flow_area=AREA_CORE_BRANCH,
        k_loss=0.0,
    )
    j_cold_23_in = MacroFlowJunction(
        name="J_PumpOutletDistributor_51_to_ColdReturnBranch_2_3_Rep",
        from_vol=header_to_cold_split,
        to_vol=build["inlet_pipe_23"].volumes[0],
        macro_vol=header_to_cold_split,
        multiplier=2,
        flow_area=AREA_CORE_BRANCH,
        k_loss=0.0,
    )
    j_cold_1_out = FlowJunction(
        name="J_ColdReturnBranch_1_to_ColdReturnOutletMerge",
        from_vol=build["inlet_pipe_1"].volumes[-1],
        to_vol=cold_return_merge,
        flow_area=AREA_CORE_BRANCH,
        k_loss=0.0,
    )
    j_cold_23_out = MacroFlowJunction(
        name="J_ColdReturnBranch_2_3_Rep_to_ColdReturnOutletMerge",
        from_vol=build["inlet_pipe_23"].volumes[-1],
        to_vol=cold_return_merge,
        macro_vol=cold_return_merge,
        multiplier=2,
        flow_area=AREA_CORE_BRANCH,
        k_loss=0.0,
    )
    j_cold_merge_to_core_inlet = FlowJunction(
        name="J_ColdReturnOutletMerge_to_CoreInletConnector",
        from_vol=cold_return_merge,
        to_vol=core_inlet,
        flow_area=AREA_CORE_BRANCH * 3.0,
        k_loss=0.0,
    )

    all_juncs.extend(_selected_old_core_junctions(build))
    all_juncs.extend(hot_outlet_entry_junctions)
    all_juncs.extend(hot_outlet_to_ring_junctions)
    all_juncs.extend(ring["segment_links"])
    all_juncs.extend(outlet_mix_to_manifold)
    all_juncs.extend(manifold_to_merge)
    all_juncs.extend([
        j_merge_to_h53,
        j_h53_to_h52,
        pump_a,
        pump_b,
        j_cold_1_in,
        j_cold_23_in,
        j_cold_1_out,
        j_cold_23_out,
        j_cold_merge_to_core_inlet,
    ])

    for junc in all_juncs:
        name = getattr(junc, "name", "")
        if _is_core_channel_junction(name):
            value = tfe_single_flow
        elif name.startswith("J_HotOutletBranch_") and "_to_InletMix_" in name:
            value = single_ring_branch_flow
        elif name.startswith("J_CoreOutletConnector_to_HotOutletBranch_"):
            value = macro_branch_flow
        elif (
            _is_ring_sector_channel_junction(name)
            or _is_ring_sector_link_junction(name)
            or _is_manifold_channel_junction(name)
            or _is_outlet_mix_to_manifold_junction(name)
            or name.startswith("J_Manifold_")
        ):
            value = single_ring_branch_flow
        elif _is_hot_outlet_channel_junction(name):
            value = macro_branch_flow
        elif _is_radiator_header_channel_junction(name):
            value = total_flow
        elif _is_cold_return_channel_junction(name):
            value = macro_branch_flow
        elif name in {
            "J_PumpOutletDistributor_51_to_ColdReturnBranch_1",
            "J_PumpOutletDistributor_51_to_ColdReturnBranch_2_3_Rep",
            "J_ColdReturnBranch_1_to_ColdReturnOutletMerge",
            "J_ColdReturnBranch_2_3_Rep_to_ColdReturnOutletMerge",
        }:
            value = macro_branch_flow
        elif (
            name.startswith("J_Radiator")
            or name.startswith("J_Pump")
            or name == "J_ColdReturnOutletMerge_to_CoreInletConnector"
        ):
            value = total_flow
        else:
            value = total_flow
        junc.W = value
        if hasattr(junc, "target_W"):
            junc.set_flow_rate(value)

    hydraulic_net = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)
    for ring_hp in ring["ring_hps"]:
        system.add_component(ring_hp)
    if outer_header_radiation_sink is not None:
        system.add_component(outer_header_radiation_sink)

    build.update(
        {
            "system": system,
            "case_version": V11_CASE_VERSION,
            "core": core,
            "tfes": tfes,
            "fluid_channels": fluid_channels,
            "core_inlet_connector": core_inlet,
            "core_outlet_connector": core_outlet,
            "hot_outlet_branches": hot_outlet_branches,
            "hot_outlet_entry_junctions": hot_outlet_entry_junctions,
            "hot_outlet_to_ring_junctions": hot_outlet_to_ring_junctions,
            "ring": ring,
            "inlet_mix_nodes": ring["inlet_mix_nodes"],
            "outlet_mix_nodes": ring["outlet_mix_nodes"],
            "manifolds": manifolds,
            "outlet_mix_to_manifold": outlet_mix_to_manifold,
            "manifold_to_merge": manifold_to_merge,
            "radiator_branch_merge": radiator_merge,
            "radiator_inner_header_53": radiator_inner_header_53,
            "radiator_outer_header_52": radiator_outer_header_52,
            "outer_header_radiation_sink": outer_header_radiation_sink,
            "pump_mid_node": pump_mid_node,
            "pump_outlet_distributor": header_to_cold_split,
            "cold_return_branch_1": build["inlet_pipe_1"],
            "cold_return_branch_2_3_rep": build["inlet_pipe_23"],
            "cold_return_merge": cold_return_merge,
            "pump_a": pump_a,
            "pump_b": pump_b,
            "j_cold_merge_to_core_inlet": j_cold_merge_to_core_inlet,
            "j_cold_1_in": j_cold_1_in,
            "j_cold_23_in": j_cold_23_in,
            "j_cold_1_out": j_cold_1_out,
            "j_cold_23_out": j_cold_23_out,
            "pump_total_head_pa": pump_total_head,
            "pump_single_head_pa": pump_single_head,
            "reference_pressure_pa": float(reference_pressure_pa),
            "total_flow_design_kg_s": total_flow,
            "macro_branch_flow_design_kg_s": macro_branch_flow,
            "single_ring_branch_flow_design_kg_s": single_ring_branch_flow,
            "single_tfe_flow_design_kg_s": tfe_single_flow,
            "coolant_material": coolant_material,
            "target_core_inlet_t_k": float(inlet_temperature_k),
            "physical_ring_count": V8_PHYSICAL_RING_COUNT,
            "passive_tfe_names": ["Ring3_Open"],
            "external_pipe_n_nodes": int(external_pipe_n_nodes),
            "connector_volume_m3": float(connector_volume_m3),
            "connector_length_m": float(connector_length_m),
            "ring_emissivity": effective_ring_emissivity,
            "hp_emissivity": effective_hp_emissivity,
            "fin_emissivity": effective_fin_emissivity,
            "outer_header_emissivity": float(outer_header_emissivity),
            "outer_header_t_space_k": effective_outer_header_t_space,
            "outer_header_area_scale": float(outer_header_area_scale),
        }
    )
    return build


def set_v11_pump_total_head(build: Dict[str, Any], total_head_pa: float) -> None:
    total_head = float(total_head_pa)
    single_head = 0.5 * total_head
    build["pump_a"].set_delta_p(single_head)
    build["pump_b"].set_delta_p(single_head)
    build["pump_total_head_pa"] = total_head
    build["pump_single_head_pa"] = single_head


def reset_v11_design_flows(build: Dict[str, Any], *, preserve_ring_restart_flows: bool = False) -> None:
    total_flow = float(build["total_flow_design_kg_s"])
    macro_branch_flow = float(build["macro_branch_flow_design_kg_s"])
    single_ring_branch_flow = float(build["single_ring_branch_flow_design_kg_s"])
    tfe_flow = float(build["single_tfe_flow_design_kg_s"])
    net = build["system"].fluid_solver
    for idx, junc in enumerate(net.junctions_obj):
        name = getattr(junc, "name", "")
        if preserve_ring_restart_flows and _is_preserved_ring_restart_junction(name):
            value = float(junc.W)
        elif name in {
            "J_RadiatorOuterHeader_52_to_PumpA",
            "J_PumpA_to_PumpOutletDistributor_51",
            "J_ColdReturnOutletMerge_to_CoreInletConnector",
        }:
            value = total_flow
        elif _is_core_channel_junction(name):
            value = tfe_flow
        elif name.startswith("J_HotOutletBranch_") and "_to_InletMix_" in name:
            value = single_ring_branch_flow
        elif name.startswith("J_CoreOutletConnector_to_HotOutletBranch_") or _is_hot_outlet_channel_junction(name):
            value = macro_branch_flow
        elif (
            _is_ring_sector_channel_junction(name)
            or _is_ring_sector_link_junction(name)
            or _is_manifold_channel_junction(name)
            or _is_outlet_mix_to_manifold_junction(name)
            or name.startswith("J_Manifold_")
        ):
            value = single_ring_branch_flow
        elif _is_cold_return_channel_junction(name) or "ColdReturnBranch_" in name:
            value = macro_branch_flow
        elif _is_radiator_header_channel_junction(name):
            value = total_flow
        elif name in {
            "J_RadiatorMerge_to_RadiatorInnerHeader_53",
            "J_RadiatorInnerHeader_53_to_RadiatorOuterHeader_52",
            "J_RadiatorOuterHeader_52_to_PumpA",
            "J_PumpA_to_PumpOutletDistributor_51",
            "J_ColdReturnOutletMerge_to_CoreInletConnector",
        }:
            value = total_flow
        else:
            value = total_flow
        junc.W = value
        if hasattr(junc, "target_W") and not (preserve_ring_restart_flows and _is_preserved_ring_restart_junction(name)):
            junc.set_flow_rate(value)
        net.W_vec[idx] = value
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()
    if hasattr(net, "W_old"):
        net.W_old[:] = net.W_vec
    if hasattr(net, "W_iterate"):
        net.W_iterate[:] = net.W_vec
    net._refresh_cached_boundary_targets()


def v11_flow_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    fixed_pressure = [vol.name for vol in net.volumes_obj if bool(getattr(vol, "is_pressure_boundary", False))]
    pressure_reference = [
        vol.name for vol in net.volumes_obj
        if bool(getattr(vol, "is_pressure_reference", False))
    ]
    hot_single = [float(j.W) for j in build["hot_outlet_to_ring_junctions"]]
    ring_out = [float(j.W) for j in build["outlet_mix_to_manifold"]]
    pump_flow = 0.5 * (float(build["pump_a"].W) + float(build["pump_b"].W))
    return {
        "closed_loop_flow_kg_s": float(build["j_cold_merge_to_core_inlet"].W),
        "pump_a_flow_kg_s": float(build["pump_a"].W),
        "pump_b_flow_kg_s": float(build["pump_b"].W),
        "pump_mean_flow_kg_s": pump_flow,
        "pump_total_head_pa": float(build["pump_a"].delta_p + build["pump_b"].delta_p),
        "pump_a_head_pa": float(build["pump_a"].delta_p),
        "pump_b_head_pa": float(build["pump_b"].delta_p),
        "hot_outlet_macro_flows_kg_s": [float(j.W) for j in build["hot_outlet_entry_junctions"]],
        "hot_to_single_ring_flows_kg_s": hot_single,
        "single_ring_in_total_kg_s": float(sum(hot_single)),
        "single_ring_out_total_kg_s": float(sum(ring_out)),
        "outlet_mix_to_manifold_flows_kg_s": ring_out,
        "macro_manifold_to_merge_flow_kg_s": float(sum(j.W * 2.0 for j in build["manifold_to_merge"])),
        "fixed_pressure_boundaries": fixed_pressure,
        "fixed_pressure_boundary_count": float(len(fixed_pressure)),
        "pressure_reference": pressure_reference,
        "pressure_reference_count": float(len(pressure_reference)),
        "has_hot_outlet_merge": bool(any(getattr(vol, "name", "") == "V9_HotOutletMerge" for vol in net.volumes_obj)),
        "has_pump_junction": bool(any(bool(getattr(j, "is_pump_junction", False)) for j in net.junctions_obj)),
        "pump_junction_count": float(sum(bool(getattr(j, "is_pump_junction", False)) for j in net.junctions_obj)),
    }


def v11_temperature_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    out = {
        "core_inlet_connector_t_k": float(build["core_inlet_connector"].T),
        "core_outlet_connector_t_k": float(build["core_outlet_connector"].T),
        "core_connector_delta_t_k": float(build["core_outlet_connector"].T - build["core_inlet_connector"].T),
        "pump_mid_node_t_k": float(build["pump_mid_node"].T),
        "pump_outlet_distributor_t_k": float(build["pump_outlet_distributor"].T),
        "cold_return_merge_t_k": float(build["cold_return_merge"].T),
    }
    for idx, channel in enumerate(build["hot_outlet_branches"], start=1):
        temps = _channel_t_in_out(channel)
        out[f"hot_outlet_branch_{idx}_t_in_k"] = temps["t_in_k"]
        out[f"hot_outlet_branch_{idx}_t_out_k"] = temps["t_out_k"]
    for idx, channel in enumerate(build["manifolds"], start=1):
        temps = _channel_t_in_out(channel)
        out[f"manifold_{idx}_t_in_k"] = temps["t_in_k"]
        out[f"manifold_{idx}_t_out_k"] = temps["t_out_k"]
    for key, vol in build["inlet_mix_nodes"].items():
        out[f"inlet_mix_{key}_t_k"] = float(vol.T)
    for key, vol in build["outlet_mix_nodes"].items():
        out[f"outlet_mix_{key}_t_k"] = float(vol.T)
    for key in ("radiator_inner_header_53", "radiator_outer_header_52", "cold_return_branch_1", "cold_return_branch_2_3_rep"):
        temps = _channel_t_in_out(build[key])
        out[f"{key}_t_in_k"] = temps["t_in_k"]
        out[f"{key}_t_out_k"] = temps["t_out_k"]
    return out


def _ring_wall_radiation_w(build: Dict[str, Any]) -> float:
    total = 0.0
    for solid in build.get("ring", {}).get("solids", []):
        for boundary in getattr(solid, "boundaries", {}).values():
            if boundary is None:
                continue
            if hasattr(boundary, "compute_net_flux_for_solver"):
                boundary.compute_net_flux_for_solver()
            for condition in getattr(boundary, "conditions", []):
                if hasattr(condition, "G_rad") and hasattr(condition, "q_flux"):
                    total += -float(np.sum(condition.q_flux))
    return total


def v11_radiation_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    outer_sink = build.get("outer_header_radiation_sink")
    outer_header_radiation_w = (
        float(getattr(outer_sink, "last_radiation_w", 0.0))
        if outer_sink is not None
        else 0.0
    )
    hp_fin_radiation_w = 0.0
    hp_fin_radiation_macro_w = 0.0
    for ring_hp in build.get("ring", {}).get("ring_hps", []):
        rejection = float(ring_hp.get_total_heat_rejection_scaled())
        hp_fin_radiation_w += rejection
        hp_fin_radiation_macro_w += 2.0 * rejection
    core_inlet_t = float(build["core_inlet_connector"].T)
    return {
        "target_core_inlet_t_k": float(build.get("target_core_inlet_t_k", np.nan)),
        "core_inlet_minus_target_k": core_inlet_t - float(build.get("target_core_inlet_t_k", core_inlet_t)),
        "radiator_cooling_delta_t_k": (
            float(build["core_outlet_connector"].T)
            - float(build["core_inlet_connector"].T)
        ),
        "radiator_outer_to_core_inlet_delta_t_k": (
            float(build["radiator_outer_header_52"].volumes[-1].T)
            - float(build["core_inlet_connector"].T)
        ),
        "ring_emissivity": float(build.get("ring_emissivity", np.nan)),
        "hp_emissivity": float(build.get("hp_emissivity", np.nan)),
        "fin_emissivity": float(build.get("fin_emissivity", np.nan)),
        "outer_header_emissivity": float(build.get("outer_header_emissivity", 0.0)),
        "outer_header_t_space_k": float(build.get("outer_header_t_space_k", np.nan)),
        "outer_header_radiation_w": outer_header_radiation_w,
        "ring_wall_radiation_w": _ring_wall_radiation_w(build),
        "hp_fin_radiation_w": hp_fin_radiation_w,
        "hp_fin_radiation_macro_w": hp_fin_radiation_macro_w,
    }


def v11_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    electric = _case_a_electric_diagnostics(build["core"])
    flow = v11_flow_diagnostics(build)
    temp = v11_temperature_diagnostics(build)
    radiation = v11_radiation_diagnostics(build)
    core_heat_power_w = sum(
        float(tfe.neutronic_data.total_power) * float(build["ring_multipliers"][name])
        for name, tfe in build["tfes"].items()
    )
    return {
        "case_version": V11_CASE_VERSION,
        "absolute_time_s": float(build["system"].global_time),
        "coolant_material": build["coolant_material"],
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "core_heat_power_w": core_heat_power_w,
        **electric,
        **flow,
        **temp,
        **radiation,
    }


__all__ = [
    "V11_CASE_VERSION",
    "V11_DEFAULT_PUMP_TOTAL_HEAD_PA",
    "build_v11_case_a_system",
    "reset_v11_design_flows",
    "set_v11_pump_total_head",
    "v11_basic_diagnostics",
    "v11_flow_diagnostics",
    "v11_temperature_diagnostics",
]
