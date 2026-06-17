import math
import os
import sys
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Components.RadiatorPipeWithFin import RadiatorPipeWithFin
from Materials.Solids.WallMaterial import SS316
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
    V8_REPRESENTATIVE_NAMES,
    V8_RING_MULTIPLIERS,
    V8_TEC_RING_MULTIPLIERS,
    _case_a_electric_diagnostics,
    build_v8_case_a_system,
)


V12_CASE_VERSION = "v12_open_core_pipefin_radiator"
V12_DEFAULT_INLET_TEMPERATURE_K = 753.330663091
V12_DEFAULT_OUTLET_PRESSURE_PA = 160000.0
V12_DEFAULT_CONNECTOR_VOLUME_M3 = 1.0e-5
V12_DEFAULT_CONNECTOR_LENGTH_M = 0.02

AREA_FLOW_NETWORK_MAIN = 3.8e-4
DH_FLOW_NETWORK_SMALL = 0.014
DH_FLOW_NETWORK_HEADER = 0.0276

FLOW_NETWORK_PIPE_SPECS = {
    "Pipe11_CoreInletHeader": {
        "sys": 11,
        "length_m": 0.13,
        "area_m2": AREA_FLOW_NETWORK_MAIN,
        "dh_m": DH_FLOW_NETWORK_SMALL,
        "initial_t_k": 745.0,
    },
    "Pipe05_CoreOutletToRadiator": {
        "sys": 5,
        "length_m": 0.13,
        "area_m2": AREA_FLOW_NETWORK_MAIN,
        "dh_m": DH_FLOW_NETWORK_SMALL,
        "initial_t_k": 845.0,
    },
    "Pipe06_RadiatorOutlet": {
        "sys": 6,
        "length_m": 0.043408,
        "area_m2": AREA_FLOW_NETWORK_MAIN,
        "dh_m": DH_FLOW_NETWORK_SMALL,
        "initial_t_k": 845.0,
    },
    "Pipe07_HeatExchangerHotSide": {
        "sys": 7,
        "length_m": 0.005426,
        "area_m2": AREA_FLOW_NETWORK_MAIN,
        "dh_m": 0.0679,
        "initial_t_k": 823.0,
    },
    "Pipe08_ReturnInnerPipe": {
        "sys": 8,
        "length_m": 0.13,
        "area_m2": AREA_FLOW_NETWORK_MAIN,
        "dh_m": 0.047,
        "initial_t_k": 745.0,
    },
    "Pipe09_ValveSegment": {
        "sys": 9,
        "length_m": 0.13,
        "area_m2": AREA_FLOW_NETWORK_MAIN,
        "dh_m": 0.047,
        "initial_t_k": 745.0,
    },
}


def nak_internal_nu(reynolds: Any, prandtl: Any, _p_d_ratio: float = 1.1) -> np.ndarray:
    pe = np.maximum(np.asarray(reynolds, dtype=float) * np.asarray(prandtl, dtype=float), 1.0)
    return 7.0 + 0.025 * pe ** 0.8


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


def _make_flow_network_pipe(
    name: str,
    *,
    material: Any,
    initial_p: float,
    n_nodes: int,
) -> IncompressibleFluidChannel:
    spec = FLOW_NETWORK_PIPE_SPECS[name]
    return _make_channel(
        name=name,
        n_nodes=n_nodes,
        length=float(spec["length_m"]),
        area=float(spec["area_m2"]),
        dh=float(spec["dh_m"]),
        material=material,
        initial_p=initial_p,
        initial_t=float(spec["initial_t_k"]),
    )


def _selected_core_junctions(build: Dict[str, Any]) -> List[Any]:
    keep = []
    for junc in build["system"].fluid_solver.junctions_obj:
        name = getattr(junc, "name", "")
        if name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_"):
            keep.append(junc)
    return keep


def _iter_channel_junctions(channels: Iterable[IncompressibleFluidChannel]):
    for channel in channels:
        yield from channel.internal_junctions


def _build_pipefin_radiator(
    *,
    material: Any,
    initial_p: float,
    initial_t: float,
    n_tubes: int,
    n_axial: int,
    n_radial_wall: int,
    n_fin_width: int,
    tube_length_m: float,
    tube_inner_diameter_m: float,
    tube_outer_diameter_m: float,
    upper_header_centerline_diameter_m: float,
    lower_header_centerline_diameter_m: float,
    header_inner_diameter_m: float,
    fin_thickness_m: float,
    fin_width_upper_m: float,
    fin_width_lower_m: float,
    tube_emissivity: float,
    fin_emissivity: float,
    tube_area_scale: float,
    fin_area_scale: float,
    t_space_k: float,
    fin_conductivity_w_m_k: float,
    fin_view_factor: float,
    fin_contact_resistance_m2k_w: float,
    fluid_solid_coupling_scheme: str,
    solid_ode_method: str,
) -> Dict[str, Any]:
    wall_mat = SS316(name="V12_TOPAZ2_Radiator_SS316")
    tube_area = math.pi * float(tube_inner_diameter_m) ** 2 / 4.0
    header_area = math.pi * float(header_inner_diameter_m) ** 2 / 4.0
    upper_seg_len = math.pi * float(upper_header_centerline_diameter_m) / int(n_tubes)
    lower_seg_len = math.pi * float(lower_header_centerline_diameter_m) / int(n_tubes)

    upper_nodes = []
    lower_nodes = []
    tube_channels = []
    radiator_units = []
    for i in range(int(n_tubes)):
        upper_nodes.append(
            _make_channel(
                name=f"V12_RadiatorUpperHeader_{i + 1:02d}",
                n_nodes=1,
                length=upper_seg_len,
                area=header_area,
                dh=float(header_inner_diameter_m),
                material=material,
                initial_p=initial_p,
                initial_t=initial_t,
            )
        )
        lower_nodes.append(
            _make_channel(
                name=f"V12_RadiatorLowerHeader_{i + 1:02d}",
                n_nodes=1,
                length=lower_seg_len,
                area=header_area,
                dh=float(header_inner_diameter_m),
                material=material,
                initial_p=initial_p,
                initial_t=initial_t,
            )
        )
        channel = _make_channel(
            name=f"V12_RadiatorTubeFluid_{i + 1:02d}",
            n_nodes=int(n_axial),
            length=float(tube_length_m),
            area=tube_area,
            dh=float(tube_inner_diameter_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        radiator = RadiatorPipeWithFin(
            name=f"V12_RadiatorTube_{i + 1:02d}",
            fluid_channel=channel,
            wall_material=wall_mat,
            tube_inner_diameter=float(tube_inner_diameter_m),
            tube_outer_diameter=float(tube_outer_diameter_m),
            tube_length=float(tube_length_m),
            n_axial=int(n_axial),
            n_radial_wall=max(1, int(n_radial_wall)),
            fin_thickness=float(fin_thickness_m),
            fin_width_upper=float(fin_width_upper_m),
            fin_width_lower=float(fin_width_lower_m),
            n_fin_width=int(n_fin_width),
            correlation_func=nak_internal_nu,
            tube_emissivity=float(tube_emissivity),
            fin_emissivity=float(fin_emissivity),
            tube_area_scale=float(tube_area_scale),
            fin_area_scale=float(fin_area_scale),
            T_space=float(t_space_k),
            initial_temp=float(initial_t),
            fin_conductivity=float(fin_conductivity_w_m_k),
            fin_view_factor=float(fin_view_factor),
            contact_resistance_m2k_w=float(fin_contact_resistance_m2k_w),
            coupling_time_scheme=fluid_solid_coupling_scheme,
            solid_ode_method=solid_ode_method,
        )
        tube_channels.append(channel)
        radiator_units.append(radiator)

    return {
        "upper_nodes": upper_nodes,
        "lower_nodes": lower_nodes,
        "tube_channels": tube_channels,
        "radiator_units": radiator_units,
        "tube_flow_area": tube_area,
        "header_area": header_area,
        "upper_seg_len": upper_seg_len,
        "lower_seg_len": lower_seg_len,
        "header_dh": float(header_inner_diameter_m),
    }


def build_v12_case_a_system(
    inlet_temperature_k: float = V12_DEFAULT_INLET_TEMPERATURE_K,
    total_inlet_flow_kg_s: float = CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    outlet_pressure_pa: float = V12_DEFAULT_OUTLET_PRESSURE_PA,
    pipe_n_nodes: int = 8,
    connector_volume_m3: float = V12_DEFAULT_CONNECTOR_VOLUME_M3,
    connector_length_m: float = V12_DEFAULT_CONNECTOR_LENGTH_M,
    solid_heat_capacity_scale: float = 1.0,
    solid_heat_capacity_scale_scope: str = "global_outer",
    coolant_material: str = "SodiumPotassium78",
    ring_multipliers: Optional[Sequence[int]] = None,
    tec_ring_multipliers: Optional[Sequence[int]] = None,
    enable_tec_coupled: bool = False,
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
    fluid_solid_coupling_scheme: str = "current",
    solid_ode_method: str = "RK45",
) -> Dict[str, Any]:
    if total_inlet_flow_kg_s <= 0.0:
        raise ValueError("total_inlet_flow_kg_s must be positive.")
    if connector_volume_m3 <= 0.0 or connector_length_m <= 0.0:
        raise ValueError("connector_volume_m3 and connector_length_m must be positive.")

    multipliers = list(V8_RING_MULTIPLIERS if ring_multipliers is None else ring_multipliers)
    if enable_tec_coupled:
        tec_multipliers = list(V8_TEC_RING_MULTIPLIERS if tec_ring_multipliers is None else tec_ring_multipliers)
    else:
        tec_multipliers = [0 for _ in multipliers]
    base = build_v8_case_a_system(
        inlet_temperature_k=float(inlet_temperature_k),
        pipe_n_nodes=int(pipe_n_nodes),
        inlet_plenum_volume_m3=float(connector_volume_m3),
        outlet_plenum_volume_m3=float(connector_volume_m3),
        plenum_length_m=float(connector_length_m),
        solid_heat_capacity_scale=float(solid_heat_capacity_scale),
        solid_heat_capacity_scale_scope=solid_heat_capacity_scale_scope,
        coolant_material=coolant_material,
        ring_multipliers=multipliers,
        tec_ring_multipliers=tec_multipliers,
    )

    material = base["inlet_boundary"].material
    core = base["core"]
    core.enable_tec_coupled = bool(enable_tec_coupled)
    if not enable_tec_coupled:
        core.thermo_calc = None
    core.point_reactor = None

    core_inlet = base["inlet_plenum"]
    core_outlet = base["outlet_plenum"]
    core_inlet.name = "V12_CoreInletConnector"
    core_outlet.name = "V12_CoreOutletConnector"
    _rename_channel(base["inlet_pipe_1"], "V12_CoreInletBranch_1")
    _rename_channel(base["inlet_pipe_23"], "V12_CoreInletBranch_2_3_Rep")
    base["j_inlet_pipe_1_out"].name = "J_V12_CoreInletBranch_1_to_CoreInletConnector"
    base["j_inlet_pipe_23_out"].name = "J_V12_CoreInletBranch_2_3_Rep_to_CoreInletConnector"

    initial_p_hot = float(outlet_pressure_pa) + 5000.0
    initial_p_cold = float(outlet_pressure_pa) + 9000.0
    initial_t = float(inlet_temperature_k)
    total_flow = float(total_inlet_flow_kg_s)
    branch_flow = total_flow / 3.0
    tube_flow = total_flow / float(n_tubes)
    tfe_single_flow = total_flow / float(sum(multipliers))

    inlet_boundary = IncompressibleBoundaryVolume(
        name="V12_InletBoundary_FixedFlow",
        material=material,
        P=initial_p_cold,
        T=initial_t,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
    )
    pipe11 = _make_flow_network_pipe(
        "Pipe11_CoreInletHeader",
        material=material,
        initial_p=initial_p_cold,
        n_nodes=max(1, int(pipe_n_nodes)),
    )
    core_inlet_distribution = IncompressibleFluidVolume(
        name="V12_CoreInletDistribution",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
        material=material,
        initial_P=initial_p_cold,
        initial_T=initial_t,
    )
    pipe5 = _make_flow_network_pipe(
        "Pipe05_CoreOutletToRadiator",
        material=material,
        initial_p=initial_p_hot,
        n_nodes=max(1, int(pipe_n_nodes)),
    )
    radiator_inlet_split = IncompressibleFluidVolume(
        name="V12_RadiatorInletSplit",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
        material=material,
        initial_P=initial_p_hot,
        initial_T=823.0,
    )
    radiator_outlet_mix = IncompressibleFluidVolume(
        name="V12_RadiatorOutletMix",
        volume=float(connector_volume_m3),
        length=float(connector_length_m),
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
        material=material,
        initial_P=float(outlet_pressure_pa) + 2000.0,
        initial_T=745.0,
    )
    cold_pipes = [
        _make_flow_network_pipe(
            name,
            material=material,
            initial_p=float(outlet_pressure_pa) + 1000.0,
            n_nodes=max(1, int(pipe_n_nodes)),
        )
        for name in (
            "Pipe06_RadiatorOutlet",
            "Pipe07_HeatExchangerHotSide",
            "Pipe08_ReturnInnerPipe",
            "Pipe09_ValveSegment",
        )
    ]
    outlet_boundary = IncompressibleBoundaryVolume(
        name="V12_OutletBoundary_FixedPressure",
        material=material,
        P=float(outlet_pressure_pa),
        T=745.0,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
    )
    outlet_boundary.is_pressure_boundary = True

    radiator = _build_pipefin_radiator(
        material=material,
        initial_p=initial_p_hot,
        initial_t=823.0,
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
        fluid_solid_coupling_scheme=fluid_solid_coupling_scheme,
        solid_ode_method=solid_ode_method,
    )

    volumes: List[Any] = [
        inlet_boundary,
        core_inlet_distribution,
        core_inlet,
        core_outlet,
        radiator_inlet_split,
        radiator_outlet_mix,
        outlet_boundary,
    ]
    junctions: List[Any] = []
    for channel in [pipe11, base["inlet_pipe_1"], base["inlet_pipe_23"], pipe5, *cold_pipes]:
        _extend_channel_objects(volumes, junctions, channel)
    for channel in base["fluid_channels"].values():
        _extend_channel_objects(volumes, junctions, channel)
    for channel in [*radiator["upper_nodes"], *radiator["lower_nodes"], *radiator["tube_channels"]]:
        _extend_channel_objects(volumes, junctions, channel)

    j_inlet = InletJunction(
        name="J_V12_InletBoundary_to_Pipe11",
        from_vol=inlet_boundary,
        to_vol=pipe11.volumes[0],
        W_initial=total_flow,
    )
    j_pipe11_to_dist = FlowJunction(
        name="J_Pipe11_to_CoreInletDistribution",
        from_vol=pipe11.volumes[-1],
        to_vol=core_inlet_distribution,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
        custom_length=float(connector_length_m),
        hydraulic_diam=DH_FLOW_NETWORK_SMALL,
    )
    j_cold_1_in = FlowJunction(
        name="J_CoreInletDistribution_to_CoreInletBranch_1",
        from_vol=core_inlet_distribution,
        to_vol=base["inlet_pipe_1"].volumes[0],
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
    )
    j_cold_23_in = MacroFlowJunction(
        name="J_CoreInletDistribution_to_CoreInletBranch_2_3_Rep",
        from_vol=core_inlet_distribution,
        to_vol=base["inlet_pipe_23"].volumes[0],
        macro_vol=core_inlet_distribution,
        multiplier=2,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
    )
    j_core_to_pipe5 = FlowJunction(
        name="J_CoreOutletConnector_to_Pipe05",
        from_vol=core_outlet,
        to_vol=pipe5.volumes[0],
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
    )
    j_pipe5_to_rad = FlowJunction(
        name="J_Pipe05_to_RadiatorInletSplit",
        from_vol=pipe5.volumes[-1],
        to_vol=radiator_inlet_split,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
    )
    inlet_idx_a = 0
    inlet_idx_b = int(n_tubes) // 2
    outlet_idx_a = 0
    outlet_idx_b = int(n_tubes) // 2
    j_rad_in_a = FlowJunction(
        name="J_RadiatorInletSplit_to_UpperHeader_A",
        from_vol=radiator_inlet_split,
        to_vol=radiator["upper_nodes"][inlet_idx_a].volumes[0],
        flow_area=radiator["header_area"],
        k_loss=float(connector_k_loss),
        hydraulic_diam=radiator["header_dh"],
    )
    j_rad_in_b = FlowJunction(
        name="J_RadiatorInletSplit_to_UpperHeader_B",
        from_vol=radiator_inlet_split,
        to_vol=radiator["upper_nodes"][inlet_idx_b].volumes[0],
        flow_area=radiator["header_area"],
        k_loss=float(connector_k_loss),
        hydraulic_diam=radiator["header_dh"],
    )
    j_rad_out_a = FlowJunction(
        name="J_LowerHeader_A_to_RadiatorOutletMix",
        from_vol=radiator["lower_nodes"][outlet_idx_a].volumes[0],
        to_vol=radiator_outlet_mix,
        flow_area=radiator["header_area"],
        k_loss=float(connector_k_loss),
        hydraulic_diam=radiator["header_dh"],
    )
    j_rad_out_b = FlowJunction(
        name="J_LowerHeader_B_to_RadiatorOutletMix",
        from_vol=radiator["lower_nodes"][outlet_idx_b].volumes[0],
        to_vol=radiator_outlet_mix,
        flow_area=radiator["header_area"],
        k_loss=float(connector_k_loss),
        hydraulic_diam=radiator["header_dh"],
    )

    radiator_junctions = [j_rad_in_a, j_rad_in_b, j_rad_out_a, j_rad_out_b]
    for i in range(int(n_tubes)):
        j = (i + 1) % int(n_tubes)
        radiator_junctions.extend(
            [
                FlowJunction(
                    name=f"J_RadiatorUpperRing_{i + 1:02d}_to_{j + 1:02d}",
                    from_vol=radiator["upper_nodes"][i].volumes[0],
                    to_vol=radiator["upper_nodes"][j].volumes[0],
                    flow_area=radiator["header_area"],
                    k_loss=float(radiator_header_k_loss),
                    custom_length=radiator["upper_seg_len"],
                    hydraulic_diam=radiator["header_dh"],
                ),
                FlowJunction(
                    name=f"J_RadiatorLowerRing_{i + 1:02d}_to_{j + 1:02d}",
                    from_vol=radiator["lower_nodes"][i].volumes[0],
                    to_vol=radiator["lower_nodes"][j].volumes[0],
                    flow_area=radiator["header_area"],
                    k_loss=float(radiator_header_k_loss),
                    custom_length=radiator["lower_seg_len"],
                    hydraulic_diam=radiator["header_dh"],
                ),
                FlowJunction(
                    name=f"J_RadiatorUpper_to_Tube_{i + 1:02d}",
                    from_vol=radiator["upper_nodes"][i].volumes[0],
                    to_vol=radiator["tube_channels"][i].volumes[0],
                    flow_area=radiator["tube_flow_area"],
                    k_loss=float(radiator_tube_inlet_k_loss),
                    custom_length=0.5 * float(tube_length_m) / int(n_axial),
                    hydraulic_diam=float(tube_inner_diameter_m),
                ),
                FlowJunction(
                    name=f"J_RadiatorTube_{i + 1:02d}_to_Lower",
                    from_vol=radiator["tube_channels"][i].volumes[-1],
                    to_vol=radiator["lower_nodes"][i].volumes[0],
                    flow_area=radiator["tube_flow_area"],
                    k_loss=float(radiator_tube_outlet_k_loss),
                    custom_length=0.5 * float(tube_length_m) / int(n_axial),
                    hydraulic_diam=float(tube_inner_diameter_m),
                ),
            ]
        )

    cold_link_junctions = []
    previous = radiator_outlet_mix
    for channel in cold_pipes:
        cold_link_junctions.append(
            FlowJunction(
                name=f"J_{getattr(previous, 'name', 'Prev')}_to_{channel.name}",
                from_vol=previous,
                to_vol=channel.volumes[0],
                flow_area=channel.area,
                k_loss=float(connector_k_loss),
                hydraulic_diam=channel.d_h,
            )
        )
        previous = channel.volumes[-1]
    j_to_outlet = FlowJunction(
        name="J_Pipe09_to_OutletBoundary",
        from_vol=previous,
        to_vol=outlet_boundary,
        flow_area=AREA_FLOW_NETWORK_MAIN,
        k_loss=float(connector_k_loss),
        hydraulic_diam=DH_FLOW_NETWORK_HEADER,
    )

    junctions.extend(
        [
            j_inlet,
            j_pipe11_to_dist,
            j_cold_1_in,
            j_cold_23_in,
            base["j_inlet_pipe_1_out"],
            base["j_inlet_pipe_23_out"],
        ]
    )
    junctions.extend(_selected_core_junctions(base))
    junctions.extend([j_core_to_pipe5, j_pipe5_to_rad, *radiator_junctions, *cold_link_junctions, j_to_outlet])

    total_flow_junctions = [
        j_inlet,
        j_pipe11_to_dist,
        j_core_to_pipe5,
        j_pipe5_to_rad,
        j_to_outlet,
        *cold_link_junctions,
        *pipe11.internal_junctions,
        *pipe5.internal_junctions,
        *list(_iter_channel_junctions(cold_pipes)),
    ]
    branch_flow_junctions = [
        j_cold_1_in,
        j_cold_23_in,
        base["j_inlet_pipe_1_out"],
        base["j_inlet_pipe_23_out"],
        *base["inlet_pipe_1"].internal_junctions,
        *base["inlet_pipe_23"].internal_junctions,
    ]
    for junc in junctions:
        name = getattr(junc, "name", "")
        if junc in total_flow_junctions:
            value = total_flow
        elif junc in branch_flow_junctions:
            value = branch_flow
        elif junc in {j_rad_in_a, j_rad_in_b, j_rad_out_a, j_rad_out_b}:
            value = 0.5 * total_flow
        elif name.startswith("J_PlenumIn_") or name.startswith("J_PlenumOut_") or name.startswith("Chan_"):
            value = tfe_single_flow
        elif (
            name.startswith("J_RadiatorUpper_to_Tube_")
            or name.startswith("J_RadiatorTube_")
            or name.startswith("V12_RadiatorTubeFluid_")
        ):
            value = tube_flow
        else:
            value = 0.0
        junc.W = float(value)
        if hasattr(junc, "set_flow_rate"):
            junc.set_flow_rate(float(value))
        elif hasattr(junc, "target_W"):
            junc.target_W = float(value)
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()

    hydraulic_net = HydraulicNetwork(volumes=volumes, junctions=junctions, gravity_vector=0.0)
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)
    for unit in radiator["radiator_units"]:
        system.add_component(unit)

    base.update(
        {
            "system": system,
            "case_version": V12_CASE_VERSION,
            "core": core,
            "tfes": base["tfes"],
            "fluid_channels": base["fluid_channels"],
            "inlet_boundary": inlet_boundary,
            "outlet_boundary": outlet_boundary,
            "core_inlet_connector": core_inlet,
            "core_outlet_connector": core_outlet,
            "core_inlet_distribution": core_inlet_distribution,
            "pipe11_core_inlet_header": pipe11,
            "pipe05_core_outlet_to_radiator": pipe5,
            "flow_network_cold_pipes": cold_pipes,
            "radiator_inlet_split": radiator_inlet_split,
            "radiator_outlet_mix": radiator_outlet_mix,
            "radiator_upper_nodes": radiator["upper_nodes"],
            "radiator_lower_nodes": radiator["lower_nodes"],
            "radiator_tube_channels": radiator["tube_channels"],
            "radiator_units": radiator["radiator_units"],
            "j_inlet": j_inlet,
            "j_pipe11_to_dist": j_pipe11_to_dist,
            "j_cold_1_in": j_cold_1_in,
            "j_cold_23_in": j_cold_23_in,
            "j_core_to_pipe5": j_core_to_pipe5,
            "j_pipe5_to_rad": j_pipe5_to_rad,
            "radiator_inlet_junctions": [j_rad_in_a, j_rad_in_b],
            "radiator_outlet_junctions": [j_rad_out_a, j_rad_out_b],
            "j_to_outlet": j_to_outlet,
            "total_flow_design_kg_s": total_flow,
            "single_branch_flow_design_kg_s": branch_flow,
            "single_tfe_flow_design_kg_s": tfe_single_flow,
            "single_radiator_tube_flow_design_kg_s": tube_flow,
            "coolant_material": coolant_material,
            "ring_multipliers": {name: int(mult) for name, mult in zip(V8_REPRESENTATIVE_NAMES, multipliers)},
            "tec_ring_multipliers": {name: int(mult) for name, mult in zip(V8_REPRESENTATIVE_NAMES, tec_multipliers)},
            "passive_tfe_names": ["Ring3_Open"] if enable_tec_coupled else list(V8_REPRESENTATIVE_NAMES),
            "pipe_n_nodes": int(pipe_n_nodes),
            "radiator_geometry": {
                "n_tubes": int(n_tubes),
                "n_axial": int(n_axial),
                "n_fin_width": int(n_fin_width),
                "tube_length_m": float(tube_length_m),
                "tube_inner_diameter_m": float(tube_inner_diameter_m),
                "tube_outer_diameter_m": float(tube_outer_diameter_m),
                "header_inner_diameter_m": float(header_inner_diameter_m),
                "tube_emissivity": float(tube_emissivity),
                "fin_emissivity": float(fin_emissivity),
                "fin_area_scale": float(fin_area_scale),
                "radiator_tube_inlet_k_loss": float(radiator_tube_inlet_k_loss),
                "radiator_tube_outlet_k_loss": float(radiator_tube_outlet_k_loss),
            },
            "flow_network_pipe_specs": FLOW_NETWORK_PIPE_SPECS,
        }
    )
    return base


def reset_v12_design_flows(build: Dict[str, Any]) -> None:
    total_flow = float(build["total_flow_design_kg_s"])
    branch_flow = float(build["single_branch_flow_design_kg_s"])
    tfe_flow = float(build["single_tfe_flow_design_kg_s"])
    tube_flow = float(build["single_radiator_tube_flow_design_kg_s"])
    half_flow = 0.5 * total_flow
    net = build["system"].fluid_solver
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
    for junc in net.junctions_obj:
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
        junc.W = float(value)
        if hasattr(junc, "set_flow_rate"):
            junc.set_flow_rate(float(value))
        elif hasattr(junc, "target_W"):
            junc.target_W = float(value)
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()
    for idx, junc in enumerate(net.junctions_obj):
        net.W_vec[idx] = float(junc.W)
    if hasattr(net, "W_old"):
        net.W_old[:] = net.W_vec
    if hasattr(net, "W_iterate"):
        net.W_iterate[:] = net.W_vec
    net._refresh_cached_boundary_targets()


def radiator_tube_mass_flows(build: Dict[str, Any]) -> np.ndarray:
    flows = []
    for channel in build["radiator_tube_channels"]:
        vals = [float(junc.W) for junc in channel.internal_junctions]
        flows.append(float(np.mean(vals)) if vals else 0.0)
    return np.asarray(flows, dtype=float)


def radiator_flow_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    flows = radiator_tube_mass_flows(build)
    mean = float(np.mean(flows)) if flows.size else 0.0
    return {
        "radiator_tube_count": int(flows.size),
        "radiator_tube_total_flow_kg_s": float(np.sum(flows)) if flows.size else 0.0,
        "radiator_tube_mean_flow_kg_s": mean,
        "radiator_tube_min_flow_kg_s": float(np.min(flows)) if flows.size else 0.0,
        "radiator_tube_max_flow_kg_s": float(np.max(flows)) if flows.size else 0.0,
        "radiator_tube_flow_spread_over_mean": float((np.max(flows) - np.min(flows)) / mean)
        if flows.size and mean
        else 0.0,
        "radiator_tube_min_index": int(np.argmin(flows) + 1) if flows.size else None,
        "radiator_tube_max_index": int(np.argmax(flows) + 1) if flows.size else None,
    }


def radiator_radiation_breakdown(build: Dict[str, Any]) -> Dict[str, float]:
    q_tube = 0.0
    q_fin = 0.0
    for unit in build["radiator_units"]:
        breakdown = unit.get_heat_exchange_breakdown()
        q_tube += float(np.sum(breakdown["bare_radiation"]))
        q_fin += float(np.sum(breakdown["fin_radiation"]))
    return {
        "q_radiator_tube_w": q_tube,
        "q_radiator_fin_w": q_fin,
        "q_radiator_total_w": q_tube + q_fin,
    }


def _channel_delta_p(channel: IncompressibleFluidChannel) -> float:
    return float(channel.volumes[0].P - channel.volumes[-1].P)


def _channel_t_in_out(channel: IncompressibleFluidChannel) -> Dict[str, float]:
    return {
        "t_in_k": float(channel.volumes[0].T),
        "t_out_k": float(channel.volumes[-1].T),
    }


def v12_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    inlet_h = float(build["core_inlet_connector"].h)
    outlet_h = float(build["core_outlet_connector"].h)
    coolant_enthalpy_rise_w = float(build["total_flow_design_kg_s"]) * (outlet_h - inlet_h)
    core_heat_power_w = sum(
        float(tfe.neutronic_data.total_power) * float(build["ring_multipliers"][name])
        for name, tfe in build["tfes"].items()
    )
    cold_pipe_temps = {}
    cold_pipe_dps = {}
    for channel in build["flow_network_cold_pipes"]:
        key = channel.name
        temps = _channel_t_in_out(channel)
        cold_pipe_temps[f"{key}_t_in_k"] = temps["t_in_k"]
        cold_pipe_temps[f"{key}_t_out_k"] = temps["t_out_k"]
        cold_pipe_dps[f"{key}_delta_p_pa"] = _channel_delta_p(channel)
    return {
        "case_version": V12_CASE_VERSION,
        "absolute_time_s": float(build["system"].global_time),
        "coolant_material": build["coolant_material"],
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "core_heat_power_w": core_heat_power_w,
        "coolant_enthalpy_rise_w": coolant_enthalpy_rise_w,
        **_case_a_electric_diagnostics(build["core"]),
        "inlet_boundary_t_k": float(build["inlet_boundary"].T),
        "core_inlet_connector_t_k": float(build["core_inlet_connector"].T),
        "core_outlet_connector_t_k": float(build["core_outlet_connector"].T),
        "core_connector_delta_t_k": float(build["core_outlet_connector"].T - build["core_inlet_connector"].T),
        "radiator_inlet_split_t_k": float(build["radiator_inlet_split"].T),
        "radiator_outlet_mix_t_k": float(build["radiator_outlet_mix"].T),
        "radiator_delta_t_k": float(build["radiator_inlet_split"].T - build["radiator_outlet_mix"].T),
        "outlet_boundary_t_k": float(build["outlet_boundary"].T),
        "inlet_boundary_pressure_pa": float(build["inlet_boundary"].P),
        "core_inlet_pressure_pa": float(build["core_inlet_connector"].P),
        "core_outlet_pressure_pa": float(build["core_outlet_connector"].P),
        "outlet_boundary_pressure_pa": float(build["outlet_boundary"].P),
        "core_delta_p_pa": float(build["core_inlet_connector"].P - build["core_outlet_connector"].P),
        "pipe11_delta_p_pa": _channel_delta_p(build["pipe11_core_inlet_header"]),
        "pipe05_delta_p_pa": _channel_delta_p(build["pipe05_core_outlet_to_radiator"]),
        "junction_count": float(len(net.junctions_obj)),
        "volume_count": float(len(net.volumes_obj)),
        "fixed_pressure_boundary_count": float(
            sum(bool(getattr(vol, "is_pressure_boundary", False)) for vol in net.volumes_obj)
        ),
        "tec_coupled_enabled": bool(getattr(build["core"], "enable_tec_coupled", False)),
        **cold_pipe_temps,
        **cold_pipe_dps,
        **radiator_flow_diagnostics(build),
        **radiator_radiation_breakdown(build),
    }


__all__ = [
    "FLOW_NETWORK_PIPE_SPECS",
    "V12_CASE_VERSION",
    "V12_DEFAULT_INLET_TEMPERATURE_K",
    "V12_DEFAULT_OUTLET_PRESSURE_PA",
    "build_v12_case_a_system",
    "radiator_flow_diagnostics",
    "radiator_tube_mass_flows",
    "reset_v12_design_flows",
    "v12_basic_diagnostics",
]
