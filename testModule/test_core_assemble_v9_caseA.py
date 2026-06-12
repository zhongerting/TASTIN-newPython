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

from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
    MacroFlowJunction,
)
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager
from test_core_assemble_v8_caseA import (
    CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    V8_PHYSICAL_RING_COUNT,
    V8_REPRESENTATIVE_NAMES,
    V8_RING_MAPPING,
    V8_RING_MULTIPLIERS,
    V8_TEC_RING_MULTIPLIERS,
    _case_a_electric_diagnostics,
    build_v8_case_a_system,
)


V9_CASE_VERSION = "v9_open_external_piping"
V9_DEFAULT_INLET_TEMPERATURE_K = 743.0
V9_DEFAULT_OUTLET_PRESSURE_PA = 160000.0
V9_DEFAULT_CONNECTOR_VOLUME_M3 = 1.0e-5
V9_DEFAULT_CONNECTOR_LENGTH_M = 0.02

AREA_SMALL_BRANCH = 2.5434e-4
DH_SMALL_BRANCH = 0.018
AREA_HEADER = 0.001734
DH_HEADER = 0.047
AREA_CORE_BRANCH = 5.9798e-4
DH_CORE_BRANCH = 0.0276

L_RADIATOR_OUTLET_38 = 0.40911
L_RADIATOR_OUTLET_44_50 = 1.41912
L_RADIATOR_INNER_HEADER_53 = 1.50969
L_RADIATOR_OUTER_HEADER_52 = 0.0915
L_COLD_RETURN_1 = 1.89021
L_COLD_RETURN_23 = 2.50705
L_HOT_OUTLET = 2.19632


def _rename_channel(channel: IncompressibleFluidChannel, new_name: str) -> None:
    channel.name = new_name
    for idx, vol in enumerate(channel.volumes, start=1):
        vol.name = f"{new_name}_Vol_{idx:02d}"
    for idx, junc in enumerate(channel.internal_junctions, start=1):
        junc.name = f"{new_name}_Junc_{idx}_{idx + 1}"


def _extend_channel_objects(vols: List[Any], juncs: List[Any], channel: IncompressibleFluidChannel) -> None:
    vols.extend(channel.volumes)
    juncs.extend(channel.internal_junctions)


def _make_channel(
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


def _selected_old_core_junctions(build: Dict[str, Any]) -> List[Any]:
    keep = []
    for junc in build["system"].fluid_solver.junctions_obj:
        name = getattr(junc, "name", "")
        if name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_"):
            keep.append(junc)
    return keep


def build_v9_case_a_system(
    inlet_temperature_k: float = V9_DEFAULT_INLET_TEMPERATURE_K,
    total_inlet_flow_kg_s: float = CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    outlet_pressure_pa: float = V9_DEFAULT_OUTLET_PRESSURE_PA,
    pipe_n_nodes: int = 8,
    external_pipe_n_nodes: int = 5,
    connector_volume_m3: float = V9_DEFAULT_CONNECTOR_VOLUME_M3,
    connector_length_m: float = V9_DEFAULT_CONNECTOR_LENGTH_M,
    solid_heat_capacity_scale: float = 1.0,
    solid_heat_capacity_scale_scope: str = "all",
    global_outer_heat_capacity_scale: Optional[float] = None,
    ring_multipliers: Optional[Sequence[int]] = None,
    tec_ring_multipliers: Optional[Sequence[int]] = None,
    coolant_material: str = "SodiumPotassium78",
) -> Dict[str, Any]:
    """Build V9 CaseA: V8 core plus open external piping without collector-ring heat pipes."""
    if total_inlet_flow_kg_s <= 0.0:
        raise ValueError("total_inlet_flow_kg_s must be positive.")
    if connector_volume_m3 <= 0.0 or connector_length_m <= 0.0:
        raise ValueError("Connector volume and length must be positive.")

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

    _rename_channel(build["inlet_pipe_1"], "ColdReturnBranch_1")
    _rename_channel(build["inlet_pipe_23"], "ColdReturnBranch_2_3_Rep")
    build["j_inlet_pipe_1_out"].name = "J_ColdReturnBranch_1_to_CoreInletConnector"
    build["j_inlet_pipe_23_out"].name = "J_ColdReturnBranch_2_3_Rep_to_CoreInletConnector"

    initial_p = float(outlet_pressure_pa)
    initial_t = float(inlet_temperature_k)
    single_cold_flow = float(total_inlet_flow_kg_s) / 3.0
    total_multiplier = float(sum(build["ring_multipliers"].values()))
    tfe_single_flow = float(total_inlet_flow_kg_s) / total_multiplier

    inlet_boundary = IncompressibleBoundaryVolume(
        name="V9_InletBoundary_FixedFlow",
        material=material,
        P=initial_p + 5000.0,
        T=initial_t,
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
    )
    inlet_connector = IncompressibleFluidVolume(
        name="V9_InletConnector",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    radiator_merge = IncompressibleFluidVolume(
        name="V9_RadiatorOutletBranchMerge",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    header_to_cold_split = IncompressibleFluidVolume(
        name="V9_HeaderToColdReturnSplit",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_HEADER,
        hydraulic_diam=DH_HEADER,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    hot_outlet_merge = IncompressibleFluidVolume(
        name="V9_HotOutletMerge",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_CORE_BRANCH * 3.0,
        hydraulic_diam=DH_CORE_BRANCH,
        material=material,
        initial_P=initial_p,
        initial_T=initial_t,
    )
    outlet_boundary = IncompressibleBoundaryVolume(
        name="V9_OutletBoundary_FixedPressure",
        material=material,
        P=float(outlet_pressure_pa),
        T=initial_t,
        flow_area=AREA_CORE_BRANCH * 3.0,
        hydraulic_diam=DH_CORE_BRANCH,
    )
    outlet_boundary.is_pressure_boundary = True

    radiator_outlet_38 = _make_channel(
        name="RadiatorOutletBranch_38",
        n_nodes=external_pipe_n_nodes,
        length=L_RADIATOR_OUTLET_38,
        area=AREA_SMALL_BRANCH,
        dh=DH_SMALL_BRANCH,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    radiator_outlet_44_50 = _make_channel(
        name="RadiatorOutletBranch_44_50_Rep",
        n_nodes=max(1, int(round(external_pipe_n_nodes * L_RADIATOR_OUTLET_44_50 / L_RADIATOR_OUTLET_38))),
        length=L_RADIATOR_OUTLET_44_50,
        area=AREA_SMALL_BRANCH,
        dh=DH_SMALL_BRANCH,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    radiator_inner_header_53 = _make_channel(
        name="RadiatorInnerHeader_53",
        n_nodes=max(1, int(round(external_pipe_n_nodes * L_RADIATOR_INNER_HEADER_53 / L_RADIATOR_OUTLET_38))),
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

    all_vols: List[Any] = [
        inlet_boundary,
        inlet_connector,
        radiator_merge,
        header_to_cold_split,
        core_inlet,
        core_outlet,
        hot_outlet_merge,
        outlet_boundary,
    ]
    all_juncs: List[Any] = []
    for channel in (
        radiator_outlet_38,
        radiator_outlet_44_50,
        radiator_inner_header_53,
        radiator_outer_header_52,
        build["inlet_pipe_1"],
        build["inlet_pipe_23"],
    ):
        _extend_channel_objects(all_vols, all_juncs, channel)
    for channel in fluid_channels.values():
        _extend_channel_objects(all_vols, all_juncs, channel)
    for channel in hot_outlet_branches:
        _extend_channel_objects(all_vols, all_juncs, channel)

    j_inlet = InletJunction(
        name="J_V9_InletBoundary_to_InletConnector",
        from_vol=inlet_boundary,
        to_vol=inlet_connector,
        W_initial=float(total_inlet_flow_kg_s),
    )
    j_rad_38_in = FlowJunction(
        name="J_InletConnector_to_RadiatorOutletBranch_38",
        from_vol=inlet_connector,
        to_vol=radiator_outlet_38.volumes[0],
        flow_area=AREA_SMALL_BRANCH,
        k_loss=0.0,
    )
    j_rad_44_50_in = MacroFlowJunction(
        name="J_InletConnector_to_RadiatorOutletBranch_44_50_Rep",
        from_vol=inlet_connector,
        to_vol=radiator_outlet_44_50.volumes[0],
        macro_vol=inlet_connector,
        multiplier=2,
        flow_area=AREA_SMALL_BRANCH,
        k_loss=0.0,
    )
    j_rad_38_out = FlowJunction(
        name="J_RadiatorOutletBranch_38_to_Merge",
        from_vol=radiator_outlet_38.volumes[-1],
        to_vol=radiator_merge,
        flow_area=AREA_SMALL_BRANCH,
        k_loss=0.0,
    )
    j_rad_44_50_out = MacroFlowJunction(
        name="J_RadiatorOutletBranch_44_50_Rep_to_Merge",
        from_vol=radiator_outlet_44_50.volumes[-1],
        to_vol=radiator_merge,
        macro_vol=radiator_merge,
        multiplier=2,
        flow_area=AREA_SMALL_BRANCH,
        k_loss=0.0,
    )
    j_merge_to_h53 = FlowJunction(
        name="J_Merge_to_RadiatorInnerHeader_53",
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
    j_h52_to_split = FlowJunction(
        name="J_RadiatorOuterHeader_52_to_ColdReturnSplit",
        from_vol=radiator_outer_header_52.volumes[-1],
        to_vol=header_to_cold_split,
        flow_area=AREA_HEADER,
        k_loss=0.0,
    )
    j_cold_1_in = FlowJunction(
        name="J_ColdReturnSplit_to_ColdReturnBranch_1",
        from_vol=header_to_cold_split,
        to_vol=build["inlet_pipe_1"].volumes[0],
        flow_area=AREA_CORE_BRANCH,
        k_loss=0.0,
    )
    j_cold_23_in = MacroFlowJunction(
        name="J_ColdReturnSplit_to_ColdReturnBranch_2_3_Rep",
        from_vol=header_to_cold_split,
        to_vol=build["inlet_pipe_23"].volumes[0],
        macro_vol=header_to_cold_split,
        multiplier=2,
        flow_area=AREA_CORE_BRANCH,
        k_loss=0.0,
    )

    hot_outlet_entry_junctions = []
    hot_outlet_exit_junctions = []
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
        hot_outlet_exit_junctions.append(
            FlowJunction(
                name=f"J_HotOutletBranch_{idx}_to_HotOutletMerge",
                from_vol=channel.volumes[-1],
                to_vol=hot_outlet_merge,
                flow_area=AREA_CORE_BRANCH,
                k_loss=0.0,
            )
        )
    j_hot_merge_to_outlet = FlowJunction(
        name="J_HotOutletMerge_to_OutletBoundary",
        from_vol=hot_outlet_merge,
        to_vol=outlet_boundary,
        flow_area=AREA_CORE_BRANCH * 3.0,
        k_loss=0.0,
    )

    all_juncs.extend(
        [
            j_inlet,
            j_rad_38_in,
            j_rad_44_50_in,
            j_rad_38_out,
            j_rad_44_50_out,
            j_merge_to_h53,
            j_h53_to_h52,
            j_h52_to_split,
            j_cold_1_in,
            j_cold_23_in,
            build["j_inlet_pipe_1_out"],
            build["j_inlet_pipe_23_out"],
        ]
    )
    all_juncs.extend(_selected_old_core_junctions(build))
    all_juncs.extend(hot_outlet_entry_junctions)
    all_juncs.extend(hot_outlet_exit_junctions)
    all_juncs.append(j_hot_merge_to_outlet)

    for junc in all_juncs:
        name = getattr(junc, "name", "")
        if name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_") or name.startswith("Chan_"):
            junc.W = tfe_single_flow
        elif junc in hot_outlet_entry_junctions or junc in hot_outlet_exit_junctions:
            junc.W = single_cold_flow
        elif junc in {
            j_cold_1_in,
            j_cold_23_in,
            build["j_inlet_pipe_1_out"],
            build["j_inlet_pipe_23_out"],
            j_rad_38_in,
            j_rad_44_50_in,
            j_rad_38_out,
            j_rad_44_50_out,
        }:
            junc.W = single_cold_flow
        else:
            junc.W = float(total_inlet_flow_kg_s)
        if hasattr(junc, "target_W"):
            junc.set_flow_rate(junc.W)

    hydraulic_net = HydraulicNetwork(
        volumes=all_vols,
        junctions=all_juncs,
        gravity_vector=0.0,
    )
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)

    build.update(
        {
            "system": system,
            "case_version": V9_CASE_VERSION,
            "core": core,
            "tfes": tfes,
            "fluid_channels": fluid_channels,
            "inlet_boundary": inlet_boundary,
            "outlet_boundary": outlet_boundary,
            "inlet_connector": inlet_connector,
            "radiator_branch_merge": radiator_merge,
            "header_to_cold_split": header_to_cold_split,
            "inlet_plenum": core_inlet,
            "outlet_plenum": core_outlet,
            "core_inlet_connector": core_inlet,
            "core_outlet_connector": core_outlet,
            "hot_outlet_merge": hot_outlet_merge,
            "radiator_outlet_branch_38": radiator_outlet_38,
            "radiator_outlet_branch_44_50_rep": radiator_outlet_44_50,
            "radiator_inner_header_53": radiator_inner_header_53,
            "radiator_outer_header_52": radiator_outer_header_52,
            "cold_return_branch_1": build["inlet_pipe_1"],
            "cold_return_branch_2_3_rep": build["inlet_pipe_23"],
            "hot_outlet_branches": hot_outlet_branches,
            "j_inlet": j_inlet,
            "j_rad_38_in": j_rad_38_in,
            "j_rad_44_50_in": j_rad_44_50_in,
            "j_rad_38_out": j_rad_38_out,
            "j_rad_44_50_out": j_rad_44_50_out,
            "j_merge_to_h53": j_merge_to_h53,
            "j_h53_to_h52": j_h53_to_h52,
            "j_h52_to_split": j_h52_to_split,
            "j_cold_1_in": j_cold_1_in,
            "j_cold_23_in": j_cold_23_in,
            "hot_outlet_entry_junctions": hot_outlet_entry_junctions,
            "hot_outlet_exit_junctions": hot_outlet_exit_junctions,
            "j_hot_merge_to_outlet": j_hot_merge_to_outlet,
            "total_flow_design_kg_s": float(total_inlet_flow_kg_s),
            "single_pipe_flow_design_kg_s": single_cold_flow,
            "single_cold_branch_flow_design_kg_s": single_cold_flow,
            "single_tfe_flow_design_kg_s": tfe_single_flow,
            "coolant_material": coolant_material,
            "physical_ring_count": V8_PHYSICAL_RING_COUNT,
            "passive_tfe_names": ["Ring3_Open"],
            "external_pipe_n_nodes": int(external_pipe_n_nodes),
            "connector_volume_m3": float(connector_volume_m3),
            "connector_length_m": float(connector_length_m),
        }
    )
    return build


def reset_v9_design_flows(build: Dict[str, Any]) -> None:
    total_flow = float(build["total_flow_design_kg_s"])
    branch_flow = float(build["single_pipe_flow_design_kg_s"])
    tfe_flow = float(build["single_tfe_flow_design_kg_s"])
    net = build["system"].fluid_solver
    total_names = {
        "J_V9_InletBoundary_to_InletConnector",
        "J_Merge_to_RadiatorInnerHeader_53",
        "J_RadiatorInnerHeader_53_to_RadiatorOuterHeader_52",
        "J_RadiatorOuterHeader_52_to_ColdReturnSplit",
        "J_HotOutletMerge_to_OutletBoundary",
    }
    branch_names = {
        "J_InletConnector_to_RadiatorOutletBranch_38",
        "J_InletConnector_to_RadiatorOutletBranch_44_50_Rep",
        "J_RadiatorOutletBranch_38_to_Merge",
        "J_RadiatorOutletBranch_44_50_Rep_to_Merge",
        "J_ColdReturnSplit_to_ColdReturnBranch_1",
        "J_ColdReturnSplit_to_ColdReturnBranch_2_3_Rep",
        "J_ColdReturnBranch_1_to_CoreInletConnector",
        "J_ColdReturnBranch_2_3_Rep_to_CoreInletConnector",
    }
    for junc in net.junctions_obj:
        name = getattr(junc, "name", "")
        if name in total_names:
            value = total_flow
        elif name in branch_names or name.startswith("J_CoreOutletConnector_to_HotOutletBranch_") or name.startswith("J_HotOutletBranch_"):
            value = branch_flow
        elif name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_") or name.startswith("Chan_"):
            value = tfe_flow
        elif any(
            token in name
            for token in (
                "RadiatorOutletBranch_38_Junc",
                "RadiatorOutletBranch_44_50_Rep_Junc",
                "ColdReturnBranch_1_Junc",
                "ColdReturnBranch_2_3_Rep_Junc",
                "HotOutletBranch_",
            )
        ):
            value = branch_flow
        else:
            value = total_flow
        junc.W = value
        if hasattr(junc, "target_W"):
            junc.set_flow_rate(value)
        elif hasattr(junc, "target_W"):
            junc.target_W = value
    for idx, junc in enumerate(net.junctions_obj):
        net.W_vec[idx] = junc.W
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()
    if hasattr(net, "W_old"):
        net.W_old[:] = net.W_vec
    if hasattr(net, "W_iterate"):
        net.W_iterate[:] = net.W_vec
    net._refresh_cached_boundary_targets()


def _channel_delta_p(channel: IncompressibleFluidChannel) -> float:
    return float(channel.volumes[0].P - channel.volumes[-1].P)


def _channel_t_in_out(channel: IncompressibleFluidChannel) -> Dict[str, float]:
    return {
        "t_in_k": float(channel.volumes[0].T),
        "t_out_k": float(channel.volumes[-1].T),
    }


def v9_flow_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    junc_by_name = {getattr(j, "name", ""): j for j in net.junctions_obj}
    hot_branch_flows = [
        float(j.W)
        for j in build["hot_outlet_entry_junctions"]
    ]
    tfe_single_flows = {
        name: float(channel.internal_junctions[0].W)
        for name, channel in build["fluid_channels"].items()
    }
    tfe_macro_flows = {
        name: float(tfe_single_flows[name] * build["ring_multipliers"][name])
        for name in build["fluid_channels"]
    }
    return {
        "inlet_total_flow_kg_s": float(build["j_inlet"].W),
        "radiator_branch_38_flow_kg_s": float(build["j_rad_38_in"].W),
        "radiator_branch_44_50_single_flow_kg_s": float(build["j_rad_44_50_in"].W),
        "radiator_branch_44_50_macro_flow_kg_s": float(build["j_rad_44_50_in"].W * 2.0),
        "cold_return_branch_1_flow_kg_s": float(build["j_cold_1_in"].W),
        "cold_return_branch_2_3_single_flow_kg_s": float(build["j_cold_23_in"].W),
        "cold_return_branch_2_3_macro_flow_kg_s": float(build["j_cold_23_in"].W * 2.0),
        "hot_outlet_branch_flows_kg_s": hot_branch_flows,
        "hot_outlet_total_flow_kg_s": float(sum(hot_branch_flows)),
        "tfe_single_flow_kg_s": tfe_single_flows,
        "tfe_macro_flow_kg_s": tfe_macro_flows,
        "tfe_total_macro_flow_kg_s": float(sum(tfe_macro_flows.values())),
        "core_inlet_pressure_pa": float(build["core_inlet_connector"].P),
        "core_outlet_pressure_pa": float(build["core_outlet_connector"].P),
        "core_delta_p_pa": float(build["core_inlet_connector"].P - build["core_outlet_connector"].P),
        "radiator_branch_38_delta_p_pa": _channel_delta_p(build["radiator_outlet_branch_38"]),
        "radiator_branch_44_50_delta_p_pa": _channel_delta_p(build["radiator_outlet_branch_44_50_rep"]),
        "radiator_inner_header_53_delta_p_pa": _channel_delta_p(build["radiator_inner_header_53"]),
        "radiator_outer_header_52_delta_p_pa": _channel_delta_p(build["radiator_outer_header_52"]),
        "cold_return_branch_1_delta_p_pa": _channel_delta_p(build["cold_return_branch_1"]),
        "cold_return_branch_2_3_delta_p_pa": _channel_delta_p(build["cold_return_branch_2_3_rep"]),
        "hot_outlet_branch_delta_p_pa": [
            _channel_delta_p(channel)
            for channel in build["hot_outlet_branches"]
        ],
        "junction_count": float(len(net.junctions_obj)),
        "volume_count": float(len(net.volumes_obj)),
        "fixed_pressure_boundary_count": float(
            sum(bool(getattr(vol, "is_pressure_boundary", False)) for vol in net.volumes_obj)
        ),
        "has_pump_junction": bool(any(bool(getattr(j, "is_pump_junction", False)) for j in net.junctions_obj)),
        "inlet_junction_target_w_kg_s": float(getattr(junc_by_name["J_V9_InletBoundary_to_InletConnector"], "target_W")),
    }


def v9_temperature_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    out = {
        "inlet_boundary_t_k": float(build["inlet_boundary"].T),
        "outlet_boundary_t_k": float(build["outlet_boundary"].T),
        "core_inlet_connector_t_k": float(build["core_inlet_connector"].T),
        "core_outlet_connector_t_k": float(build["core_outlet_connector"].T),
        "core_connector_delta_t_k": float(build["core_outlet_connector"].T - build["core_inlet_connector"].T),
    }
    for key in (
        "radiator_outlet_branch_38",
        "radiator_outlet_branch_44_50_rep",
        "radiator_inner_header_53",
        "radiator_outer_header_52",
        "cold_return_branch_1",
        "cold_return_branch_2_3_rep",
    ):
        temps = _channel_t_in_out(build[key])
        out[f"{key}_t_in_k"] = temps["t_in_k"]
        out[f"{key}_t_out_k"] = temps["t_out_k"]
    for idx, channel in enumerate(build["hot_outlet_branches"], start=1):
        temps = _channel_t_in_out(channel)
        out[f"hot_outlet_branch_{idx}_t_in_k"] = temps["t_in_k"]
        out[f"hot_outlet_branch_{idx}_t_out_k"] = temps["t_out_k"]
    return out


def v9_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    electric = _case_a_electric_diagnostics(build["core"])
    flow = v9_flow_diagnostics(build)
    temp = v9_temperature_diagnostics(build)
    inlet_h = float(getattr(build["inlet_boundary"], "h", np.nan))
    outlet_h = float(getattr(build["hot_outlet_merge"], "h", np.nan))
    coolant_enthalpy_rise_w = float(build["total_flow_design_kg_s"]) * (outlet_h - inlet_h)
    core_heat_power_w = sum(
        float(tfe.neutronic_data.total_power) * float(build["ring_multipliers"][name])
        for name, tfe in build["tfes"].items()
    )
    return {
        "case_version": V9_CASE_VERSION,
        "absolute_time_s": float(build["system"].global_time),
        "coolant_material": build["coolant_material"],
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "core_heat_power_w": core_heat_power_w,
        "coolant_enthalpy_rise_w": coolant_enthalpy_rise_w,
        **electric,
        **flow,
        **temp,
    }


__all__ = [
    "AREA_CORE_BRANCH",
    "CASE_A_DESIGN_TOTAL_FLOW_KG_S",
    "V8_REPRESENTATIVE_NAMES",
    "V8_RING_MAPPING",
    "V8_RING_MULTIPLIERS",
    "V8_TEC_RING_MULTIPLIERS",
    "V9_CASE_VERSION",
    "build_v9_case_a_system",
    "reset_v9_design_flows",
    "v9_basic_diagnostics",
    "v9_flow_diagnostics",
    "v9_temperature_diagnostics",
]
