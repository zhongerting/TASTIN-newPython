import os
import sys
from typing import Any, Dict, Optional, Sequence

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidVolume, PumpJunction
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager
from test_core_assemble_v12_caseA import (
    AREA_FLOW_NETWORK_MAIN,
    DH_FLOW_NETWORK_SMALL,
    V12_DEFAULT_CONNECTOR_LENGTH_M,
    V12_DEFAULT_CONNECTOR_VOLUME_M3,
    build_v12_case_a_system,
    radiator_flow_diagnostics,
    radiator_radiation_breakdown,
)
from test_core_assemble_v8_caseA import (
    CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    V8_REPRESENTATIVE_NAMES,
    _case_a_electric_diagnostics,
)


V13_CASE_VERSION = "v13_closed_core_pipefin_radiator_pumped_loop"
V13_DEFAULT_INLET_TEMPERATURE_K = 727.0
V13_DEFAULT_REFERENCE_PRESSURE_PA = 207927.58
V13_DEFAULT_PUMP_TOTAL_HEAD_PA = 7900.0


class FlowControlledPumpJunction(PumpJunction):
    """Pump-topology junction with a prescribed mass-flow target."""

    def __init__(self, *args: Any, W_initial: float = 0.0, **kwargs: Any):
        super().__init__(*args, **kwargs)
        self.W = float(W_initial)
        self.target_W = float(W_initial)

    def set_flow_rate(self, W: float) -> None:
        self.target_W = float(W)


def _set_junction_flow(junc: Any, value: float) -> None:
    junc.W = float(value)
    if hasattr(junc, "set_flow_rate"):
        junc.set_flow_rate(float(value))
    elif hasattr(junc, "target_W"):
        junc.target_W = float(value)
    if hasattr(junc, "update_velocity"):
        junc.update_velocity()


def _channel_delta_p(channel: Any) -> float:
    return float(channel.volumes[0].P - channel.volumes[-1].P)


def _channel_t_in_out(channel: Any) -> Dict[str, float]:
    return {
        "t_in_k": float(channel.volumes[0].T),
        "t_out_k": float(channel.volumes[-1].T),
    }


def build_v13_case_a_system(
    inlet_temperature_k: float = V13_DEFAULT_INLET_TEMPERATURE_K,
    total_inlet_flow_kg_s: float = CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    reference_pressure_pa: float = V13_DEFAULT_REFERENCE_PRESSURE_PA,
    pump_total_head_pa: float = V13_DEFAULT_PUMP_TOTAL_HEAD_PA,
    pipe_n_nodes: int = 8,
    connector_volume_m3: float = V12_DEFAULT_CONNECTOR_VOLUME_M3,
    connector_length_m: float = V12_DEFAULT_CONNECTOR_LENGTH_M,
    solid_heat_capacity_scale: float = 1.0,
    solid_heat_capacity_scale_scope: str = "global_outer",
    coolant_material: str = "SodiumPotassium78",
    ring_multipliers: Optional[Sequence[int]] = None,
    tec_ring_multipliers: Optional[Sequence[int]] = None,
    enable_tec_coupled: bool = True,
    n_tubes: int = 78,
    n_axial: int = 8,
    n_radial_wall: int = 1,
    n_fin_width: int = 12,
    tube_length_m: float = 1.85,
    tube_inner_diameter_m: float = 0.007,
    tube_outer_diameter_m: float = 0.008,
    upper_header_centerline_diameter_m: float = 0.824,
    lower_header_centerline_diameter_m: float = 1.346,
    header_inner_diameter_m: float = 0.020,
    fin_thickness_m: float = 0.0004,
    fin_width_upper_m: float = 0.03319,
    fin_width_lower_m: float = 0.05421,
    tube_emissivity: float = 0.80,
    fin_emissivity: float = 0.80,
    tube_area_scale: float = 1.0,
    fin_area_scale: float = 0.35,
    t_space_k: float = 3.0,
    fin_conductivity_w_m_k: float = 348.9,
    fin_view_factor: float = 1.0,
    fin_contact_resistance_m2k_w: float = 0.0,
    radiator_header_k_loss: float = 1.0,
    radiator_tube_inlet_k_loss: float = 100.0,
    radiator_tube_outlet_k_loss: float = 100.0,
    connector_k_loss: float = 0.0,
    pump_flow_control: bool = True,
    fluid_solid_coupling_scheme: str = "local_implicit",
    solid_ode_method: str = "RK45",
) -> Dict[str, Any]:
    if total_inlet_flow_kg_s <= 0.0:
        raise ValueError("total_inlet_flow_kg_s must be positive.")
    if pump_total_head_pa <= 0.0:
        raise ValueError("pump_total_head_pa must be positive.")

    base = build_v12_case_a_system(
        inlet_temperature_k=float(inlet_temperature_k),
        total_inlet_flow_kg_s=float(total_inlet_flow_kg_s),
        outlet_pressure_pa=float(reference_pressure_pa) - 9000.0,
        pipe_n_nodes=int(pipe_n_nodes),
        connector_volume_m3=float(connector_volume_m3),
        connector_length_m=float(connector_length_m),
        solid_heat_capacity_scale=float(solid_heat_capacity_scale),
        solid_heat_capacity_scale_scope=solid_heat_capacity_scale_scope,
        coolant_material=coolant_material,
        ring_multipliers=ring_multipliers,
        tec_ring_multipliers=tec_ring_multipliers,
        enable_tec_coupled=bool(enable_tec_coupled),
        n_tubes=int(n_tubes),
        n_axial=int(n_axial),
        n_radial_wall=int(n_radial_wall),
        n_fin_width=int(n_fin_width),
        tube_length_m=float(tube_length_m),
        tube_inner_diameter_m=float(tube_inner_diameter_m),
        tube_outer_diameter_m=float(tube_outer_diameter_m),
        upper_header_centerline_diameter_m=float(upper_header_centerline_diameter_m),
        lower_header_centerline_diameter_m=float(lower_header_centerline_diameter_m),
        header_inner_diameter_m=float(header_inner_diameter_m),
        fin_thickness_m=float(fin_thickness_m),
        fin_width_upper_m=float(fin_width_upper_m),
        fin_width_lower_m=float(fin_width_lower_m),
        tube_emissivity=float(tube_emissivity),
        fin_emissivity=float(fin_emissivity),
        tube_area_scale=float(tube_area_scale),
        fin_area_scale=float(fin_area_scale),
        t_space_k=float(t_space_k),
        fin_conductivity_w_m_k=float(fin_conductivity_w_m_k),
        fin_view_factor=float(fin_view_factor),
        fin_contact_resistance_m2k_w=float(fin_contact_resistance_m2k_w),
        radiator_header_k_loss=float(radiator_header_k_loss),
        radiator_tube_inlet_k_loss=float(radiator_tube_inlet_k_loss),
        radiator_tube_outlet_k_loss=float(radiator_tube_outlet_k_loss),
        connector_k_loss=float(connector_k_loss),
        fluid_solid_coupling_scheme=fluid_solid_coupling_scheme,
        solid_ode_method=solid_ode_method,
    )

    core_inlet = base["core_inlet_connector"]
    core_inlet.is_pressure_boundary = False
    core_inlet.is_pressure_reference = True
    core_inlet.P = float(reference_pressure_pa)
    core_inlet.target_P = float(reference_pressure_pa)

    old_net = base["system"].fluid_solver
    excluded_volumes = {base["inlet_boundary"], base["outlet_boundary"]}
    excluded_junctions = {base["j_inlet"], base["j_to_outlet"]}
    volumes = [vol for vol in old_net.volumes_obj if vol not in excluded_volumes]
    junctions = [
        junc for junc in old_net.junctions_obj
        if junc not in excluded_junctions
        and junc.from_vol not in excluded_volumes
        and junc.to_vol not in excluded_volumes
    ]

    material = base["core_inlet_connector"].material
    pump_single_head = 0.5 * float(pump_total_head_pa)
    pump_mid_node = IncompressibleFluidVolume(
        name="V13_PumpMidNode",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
        material=material,
        initial_P=float(reference_pressure_pa) - pump_single_head,
        initial_T=float(inlet_temperature_k),
    )
    pump_outlet_node = IncompressibleFluidVolume(
        name="V13_PumpOutletNode",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
        material=material,
        initial_P=float(reference_pressure_pa),
        initial_T=float(inlet_temperature_k),
    )
    pipe09_last = base["flow_network_cold_pipes"][-1].volumes[-1]
    pipe11_first = base["pipe11_core_inlet_header"].volumes[0]
    pump_cls = FlowControlledPumpJunction if bool(pump_flow_control) else PumpJunction
    pump_kwargs = {"W_initial": float(total_inlet_flow_kg_s)} if bool(pump_flow_control) else {}
    pump_a = pump_cls(
        name="J_Pipe09_to_V13_PumpA",
        from_vol=pipe09_last,
        to_vol=pump_mid_node,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
        custom_length=float(connector_length_m),
        delta_p=pump_single_head,
        **pump_kwargs,
    )
    pump_b = pump_cls(
        name="J_V13_PumpA_to_PumpB",
        from_vol=pump_mid_node,
        to_vol=pump_outlet_node,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
        custom_length=float(connector_length_m),
        delta_p=pump_single_head,
        **pump_kwargs,
    )
    j_pump_to_pipe11 = FlowJunction(
        name="J_V13_PumpOutlet_to_Pipe11",
        from_vol=pump_outlet_node,
        to_vol=pipe11_first,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
        custom_length=float(connector_length_m),
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
    )

    volumes.extend([pump_mid_node, pump_outlet_node])
    junctions.extend([pump_a, pump_b, j_pump_to_pipe11])

    hydraulic_net = HydraulicNetwork(volumes=volumes, junctions=junctions, gravity_vector=0.0)
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(base["core"])
    for unit in base["radiator_units"]:
        system.add_component(unit)

    base.update(
        {
            "system": system,
            "case_version": V13_CASE_VERSION,
            "inlet_boundary": None,
            "outlet_boundary": None,
            "pump_mid_node": pump_mid_node,
            "pump_outlet_node": pump_outlet_node,
            "pump_a": pump_a,
            "pump_b": pump_b,
            "j_pump_to_pipe11": j_pump_to_pipe11,
            "pump_total_head_pa": float(pump_total_head_pa),
            "pump_single_head_pa": pump_single_head,
            "pump_flow_control": bool(pump_flow_control),
            "reference_pressure_pa": float(reference_pressure_pa),
            "open_loop_source_case_version": "v12_open_core_pipefin_radiator",
        }
    )
    reset_v13_design_flows(base)
    return base


def set_v13_pump_total_head(build: Dict[str, Any], total_head_pa: float) -> None:
    total_head = float(total_head_pa)
    single_head = 0.5 * total_head
    build["pump_a"].set_delta_p(single_head)
    build["pump_b"].set_delta_p(single_head)
    build["pump_total_head_pa"] = total_head
    build["pump_single_head_pa"] = single_head


def reset_v13_design_flows(build: Dict[str, Any]) -> None:
    total_flow = float(build["total_flow_design_kg_s"])
    branch_flow = float(build["single_branch_flow_design_kg_s"])
    tfe_flow = float(build["single_tfe_flow_design_kg_s"])
    tube_flow = float(build["single_radiator_tube_flow_design_kg_s"])
    half_flow = 0.5 * total_flow
    branch_names = {
        "J_CoreInletDistribution_to_CoreInletBranch_1",
        "J_CoreInletDistribution_to_CoreInletBranch_2_3_Rep",
        "J_V12_CoreInletBranch_1_to_CoreInletConnector",
        "J_V12_CoreInletBranch_2_3_Rep_to_CoreInletConnector",
    }
    half_names = {
        "J_RadiatorInletSplit_to_UpperHeader_A",
        "J_RadiatorInletSplit_to_UpperHeader_B",
        "J_LowerHeader_A_to_RadiatorOutletMix",
        "J_LowerHeader_B_to_RadiatorOutletMix",
    }
    net = build["system"].fluid_solver
    for idx, junc in enumerate(net.junctions_obj):
        name = getattr(junc, "name", "")
        if name in branch_names or name.startswith("V12_CoreInletBranch_"):
            value = branch_flow
        elif name in half_names:
            value = half_flow
        elif name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_") or name.startswith("Chan_"):
            value = tfe_flow
        elif (
            name.startswith("J_RadiatorUpper_to_Tube_")
            or name.startswith("J_RadiatorTube_")
            or name.startswith("V12_RadiatorTubeFluid_")
        ):
            value = tube_flow
        elif name.startswith("J_RadiatorUpperRing_") or name.startswith("J_RadiatorLowerRing_"):
            value = 0.0
        else:
            value = total_flow
        _set_junction_flow(junc, value)
        net.W_vec[idx] = float(junc.W)
    if hasattr(net, "W_old"):
        net.W_old[:] = net.W_vec
    if hasattr(net, "W_iterate"):
        net.W_iterate[:] = net.W_vec
    net._refresh_cached_boundary_targets()


def v13_flow_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    fixed_pressure = [vol.name for vol in net.volumes_obj if bool(getattr(vol, "is_pressure_boundary", False))]
    pressure_reference = [
        vol.name for vol in net.volumes_obj
        if bool(getattr(vol, "is_pressure_reference", False))
    ]
    pump_flow = 0.5 * (float(build["pump_a"].W) + float(build["pump_b"].W))
    return {
        "closed_loop_flow_kg_s": float(build["j_pump_to_pipe11"].W),
        "pump_a_flow_kg_s": float(build["pump_a"].W),
        "pump_b_flow_kg_s": float(build["pump_b"].W),
        "pump_mean_flow_kg_s": pump_flow,
        "pump_total_head_pa": float(build["pump_a"].delta_p + build["pump_b"].delta_p),
        "pump_a_head_pa": float(build["pump_a"].delta_p),
        "pump_b_head_pa": float(build["pump_b"].delta_p),
        "pump_flow_control": bool(build.get("pump_flow_control", False)),
        "fixed_pressure_boundaries": fixed_pressure,
        "fixed_pressure_boundary_count": float(len(fixed_pressure)),
        "pressure_reference": pressure_reference,
        "pressure_reference_count": float(len(pressure_reference)),
        "has_pump_junction": bool(any(bool(getattr(j, "is_pump_junction", False)) for j in net.junctions_obj)),
        "pump_junction_count": float(sum(bool(getattr(j, "is_pump_junction", False)) for j in net.junctions_obj)),
    }


def v13_temperature_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    cold_pipe_temps = {}
    cold_pipe_dps = {}
    for channel in build["flow_network_cold_pipes"]:
        key = channel.name
        temps = _channel_t_in_out(channel)
        cold_pipe_temps[f"{key}_t_in_k"] = temps["t_in_k"]
        cold_pipe_temps[f"{key}_t_out_k"] = temps["t_out_k"]
        cold_pipe_dps[f"{key}_delta_p_pa"] = _channel_delta_p(channel)
    return {
        "core_inlet_connector_t_k": float(build["core_inlet_connector"].T),
        "core_outlet_connector_t_k": float(build["core_outlet_connector"].T),
        "core_connector_delta_t_k": float(build["core_outlet_connector"].T - build["core_inlet_connector"].T),
        "radiator_inlet_split_t_k": float(build["radiator_inlet_split"].T),
        "radiator_outlet_mix_t_k": float(build["radiator_outlet_mix"].T),
        "radiator_delta_t_k": float(build["radiator_inlet_split"].T - build["radiator_outlet_mix"].T),
        "pump_mid_node_t_k": float(build["pump_mid_node"].T),
        "pump_outlet_node_t_k": float(build["pump_outlet_node"].T),
        **cold_pipe_temps,
        **cold_pipe_dps,
    }


def v13_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    inlet_h = float(build["core_inlet_connector"].h)
    outlet_h = float(build["core_outlet_connector"].h)
    coolant_enthalpy_rise_w = float(v13_flow_diagnostics(build)["closed_loop_flow_kg_s"]) * (outlet_h - inlet_h)
    core_heat_power_w = sum(
        float(tfe.neutronic_data.total_power) * float(build["ring_multipliers"][name])
        for name, tfe in build["tfes"].items()
    )
    return {
        "case_version": V13_CASE_VERSION,
        "absolute_time_s": float(build["system"].global_time),
        "coolant_material": build["coolant_material"],
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "core_heat_power_w": core_heat_power_w,
        "coolant_enthalpy_rise_w": coolant_enthalpy_rise_w,
        **_case_a_electric_diagnostics(build["core"]),
        "core_inlet_pressure_pa": float(build["core_inlet_connector"].P),
        "core_outlet_pressure_pa": float(build["core_outlet_connector"].P),
        "core_delta_p_pa": float(build["core_inlet_connector"].P - build["core_outlet_connector"].P),
        "pipe11_delta_p_pa": _channel_delta_p(build["pipe11_core_inlet_header"]),
        "pipe05_delta_p_pa": _channel_delta_p(build["pipe05_core_outlet_to_radiator"]),
        "junction_count": float(len(net.junctions_obj)),
        "volume_count": float(len(net.volumes_obj)),
        "tec_coupled_enabled": bool(getattr(build["core"], "enable_tec_coupled", False)),
        **v13_temperature_diagnostics(build),
        **v13_flow_diagnostics(build),
        **radiator_flow_diagnostics(build),
        **radiator_radiation_breakdown(build),
    }


__all__ = [
    "V13_CASE_VERSION",
    "V13_DEFAULT_INLET_TEMPERATURE_K",
    "V13_DEFAULT_PUMP_TOTAL_HEAD_PA",
    "V13_DEFAULT_REFERENCE_PRESSURE_PA",
    "build_v13_case_a_system",
    "reset_v13_design_flows",
    "set_v13_pump_total_head",
    "v13_basic_diagnostics",
    "v13_flow_diagnostics",
]
