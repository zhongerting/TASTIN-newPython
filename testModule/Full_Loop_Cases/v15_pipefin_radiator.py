"""Independent V15 pipe-fin radiator adapter for Full_Loop_Cases."""

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np

from Components.RadiatorPipeWithFin import RadiatorPipeWithFin
from Materials.Solids.WallMaterial import SS316
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidVolume

from .common_flow_builder import extend_channel_objects, make_channel


def nak_internal_nu(reynolds: Any, prandtl: Any, _p_d_ratio: float = 1.1) -> np.ndarray:
    pe = np.maximum(np.asarray(reynolds, dtype=float) * np.asarray(prandtl, dtype=float), 1.0)
    return 7.0 + 0.025 * pe**0.8


@dataclass(frozen=True)
class V15PipeFinRadiatorConfig:
    """V15 pipe-fin radiator geometry and loss parameters."""

    n_tubes: int = 78
    n_axial: int = 8
    n_radial_wall: int = 1
    n_fin_width: int = 12
    tube_length_m: float = 1.85
    tube_inner_diameter_m: float = 0.007
    tube_outer_diameter_m: float = 0.008
    upper_header_centerline_diameter_m: float = 0.824
    lower_header_centerline_diameter_m: float = 1.346
    header_inner_diameter_m: float = 0.020
    fin_thickness_m: float = 0.0004
    fin_width_upper_m: float = 0.03319
    fin_width_lower_m: float = 0.05421
    tube_emissivity: float = 0.80
    fin_emissivity: float = 0.80
    tube_area_scale: float = 1.0
    fin_area_scale: float = 0.35
    t_space_k: float = 3.0
    fin_conductivity_w_m_k: float = 348.9
    fin_view_factor: float = 1.0
    fin_contact_resistance_m2k_w: float = 0.0
    radiator_header_k_loss: float = 1.0
    radiator_tube_inlet_k_loss: float = 100.0
    radiator_tube_outlet_k_loss: float = 100.0
    connector_k_loss: float = 0.0
    cold_return_branch_length_m: float = 1.89021
    cold_return_branch_area_m2: float = float(np.pi * 0.0138**2)
    cold_return_branch_dh_m: float = 0.0276
    cold_return_branch_n_nodes: int = 1
    fluid_solid_coupling_scheme: str = "current"
    solid_ode_method: str = "RK45"

    @property
    def tube_flow_area_m2(self) -> float:
        return float(np.pi * self.tube_inner_diameter_m**2 / 4.0)

    @property
    def header_flow_area_m2(self) -> float:
        return float(np.pi * self.header_inner_diameter_m**2 / 4.0)

    @property
    def upper_header_segment_length_m(self) -> float:
        return float(np.pi * self.upper_header_centerline_diameter_m / int(self.n_tubes))

    @property
    def lower_header_segment_length_m(self) -> float:
        return float(np.pi * self.lower_header_centerline_diameter_m / int(self.n_tubes))


def _validate_config(config: V15PipeFinRadiatorConfig) -> None:
    if int(config.n_tubes) != 78:
        raise ValueError("V15 requires n_tubes=78 for the first full-loop pipe-fin model.")
    if int(config.n_axial) <= 0:
        raise ValueError("n_axial must be positive.")
    if int(config.cold_return_branch_n_nodes) <= 0:
        raise ValueError("cold_return_branch_n_nodes must be positive.")


def _set_design_flow(junction: Any, flow_kg_s: float) -> None:
    value = float(flow_kg_s)
    junction.design_flow_kg_s = value
    junction.W = value
    if hasattr(junction, "set_flow_rate"):
        junction.set_flow_rate(value)
    elif hasattr(junction, "target_W"):
        junction.target_W = value
    if hasattr(junction, "update_velocity"):
        junction.update_velocity()


def _set_channel_design_flow(channel: Any, flow_kg_s: float) -> None:
    for junction in channel.internal_junctions:
        _set_design_flow(junction, flow_kg_s)


def _make_volume(
    *,
    name: str,
    volume: float,
    length: float,
    area: float,
    dh: float,
    material: Any,
    initial_p: float,
    initial_t: float,
) -> IncompressibleFluidVolume:
    return IncompressibleFluidVolume(
        name=name,
        volume=float(volume),
        length=float(length),
        flow_area=float(area),
        hydraulic_diam=float(dh),
        material=material,
        initial_P=float(initial_p),
        initial_T=float(initial_t),
    )


def _rename_channel(channel: Any, new_name: str) -> None:
    channel.name = new_name
    for idx, vol in enumerate(channel.volumes, start=1):
        vol.name = f"{new_name}_Vol_{idx:02d}"
    for idx, junc in enumerate(channel.internal_junctions, start=1):
        junc.name = f"{new_name}_Junc_{idx}_{idx + 1}"


def _make_radiator_unit(
    *,
    index: int,
    channel: Any,
    wall_material: Any,
    initial_t: float,
    config: V15PipeFinRadiatorConfig,
) -> RadiatorPipeWithFin:
    return RadiatorPipeWithFin(
        name=f"RadiatorTube_{index:02d}",
        fluid_channel=channel,
        wall_material=wall_material,
        tube_inner_diameter=float(config.tube_inner_diameter_m),
        tube_outer_diameter=float(config.tube_outer_diameter_m),
        tube_length=float(config.tube_length_m),
        n_axial=int(config.n_axial),
        n_radial_wall=max(1, int(config.n_radial_wall)),
        fin_thickness=float(config.fin_thickness_m),
        fin_width_upper=float(config.fin_width_upper_m),
        fin_width_lower=float(config.fin_width_lower_m),
        n_fin_width=int(config.n_fin_width),
        correlation_func=nak_internal_nu,
        tube_emissivity=float(config.tube_emissivity),
        fin_emissivity=float(config.fin_emissivity),
        tube_area_scale=float(config.tube_area_scale),
        fin_area_scale=float(config.fin_area_scale),
        T_space=float(config.t_space_k),
        initial_temp=float(initial_t),
        fin_conductivity=float(config.fin_conductivity_w_m_k),
        fin_view_factor=float(config.fin_view_factor),
        contact_resistance_m2k_w=float(config.fin_contact_resistance_m2k_w),
        coupling_time_scheme=str(config.fluid_solid_coupling_scheme),
        solid_ode_method=str(config.solid_ode_method),
    )


def attach_v15_pipefin_radiator(build: Dict[str, Any], config: V15PipeFinRadiatorConfig) -> None:
    """Attach V15 pipe-fin radiator and cold-return branches to a common build."""
    _validate_config(config)

    material = build["core_inlet_connector"].material
    initial_p = float(build["core_inlet_connector"].P)
    initial_t = float(build["core_inlet_connector"].T)
    total_flow = float(build["total_flow_design_kg_s"])
    half_flow = 0.5 * total_flow
    tube_flow = total_flow / float(config.n_tubes)
    cold_branch_flow = total_flow / 3.0

    _rename_channel(build["radiator_outlet_header"], "RadiatorOuterHeader")

    volumes: List[Any] = []
    junctions: List[Any] = []
    components: List[Any] = []

    flow_area_main = float(build["radiator_inlet_header"].area)
    dh_main = float(build["radiator_inlet_header"].d_h)
    connector_length = 0.02
    connector_volume = 1.0e-5

    inlet_distributor = _make_volume(
        name="RadiatorInletDistributor",
        volume=connector_volume,
        length=connector_length,
        area=flow_area_main,
        dh=dh_main,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    inner_header = _make_volume(
        name="RadiatorInnerHeader",
        volume=connector_volume,
        length=connector_length,
        area=flow_area_main,
        dh=dh_main,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    pump_outlet_distributor = _make_volume(
        name="PumpOutletDistributor",
        volume=connector_volume,
        length=connector_length,
        area=flow_area_main,
        dh=dh_main,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
    )
    volumes.extend([inlet_distributor, inner_header, pump_outlet_distributor])

    j_rad_header_to_distributor = FlowJunction(
        name="J_RadiatorInletHeader_to_RadiatorInletDistributor",
        from_vol=build["radiator_inlet_header"].volumes[-1],
        to_vol=inlet_distributor,
        flow_area=flow_area_main,
        k_loss=float(config.connector_k_loss),
        hydraulic_diam=dh_main,
    )
    j_inner_to_outer = FlowJunction(
        name="J_RadiatorInnerHeader_to_RadiatorOuterHeader",
        from_vol=inner_header,
        to_vol=build["radiator_outlet_header"].volumes[0],
        flow_area=flow_area_main,
        k_loss=float(config.connector_k_loss),
        hydraulic_diam=dh_main,
    )
    j_pump_to_distributor = FlowJunction(
        name="J_PumpOutletNode_to_PumpOutletDistributor",
        from_vol=build["pump_outlet_node"],
        to_vol=pump_outlet_distributor,
        flow_area=flow_area_main,
        k_loss=float(config.connector_k_loss),
        hydraulic_diam=dh_main,
    )
    for junc in (j_rad_header_to_distributor, j_inner_to_outer, j_pump_to_distributor):
        _set_design_flow(junc, total_flow)
    junctions.extend([j_rad_header_to_distributor, j_inner_to_outer, j_pump_to_distributor])

    wall_mat = SS316(name="V15_Radiator_SS316")
    upper_headers = []
    lower_headers = []
    tube_channels = []
    radiator_units = []
    for idx in range(1, int(config.n_tubes) + 1):
        upper = make_channel(
            name=f"RadiatorUpperHeader_{idx:02d}",
            n_nodes=1,
            length=float(config.upper_header_segment_length_m),
            area=float(config.header_flow_area_m2),
            dh=float(config.header_inner_diameter_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        lower = make_channel(
            name=f"RadiatorLowerHeader_{idx:02d}",
            n_nodes=1,
            length=float(config.lower_header_segment_length_m),
            area=float(config.header_flow_area_m2),
            dh=float(config.header_inner_diameter_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        tube = make_channel(
            name=f"RadiatorTubeFluid_{idx:02d}",
            n_nodes=int(config.n_axial),
            length=float(config.tube_length_m),
            area=float(config.tube_flow_area_m2),
            dh=float(config.tube_inner_diameter_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        extend_channel_objects(volumes, junctions, upper)
        extend_channel_objects(volumes, junctions, lower)
        extend_channel_objects(volumes, junctions, tube)
        _set_channel_design_flow(upper, half_flow)
        _set_channel_design_flow(lower, half_flow)
        _set_channel_design_flow(tube, tube_flow)
        unit = _make_radiator_unit(
            index=idx,
            channel=tube,
            wall_material=wall_mat,
            initial_t=initial_t,
            config=config,
        )
        components.append(unit)
        upper_headers.append(upper)
        lower_headers.append(lower)
        tube_channels.append(tube)
        radiator_units.append(unit)

    inlet_idx_a = 0
    inlet_idx_b = int(config.n_tubes) // 2
    for suffix, idx in (("A", inlet_idx_a), ("B", inlet_idx_b)):
        junc = FlowJunction(
            name=f"J_RadiatorInletDistributor_to_UpperHeader_{suffix}",
            from_vol=inlet_distributor,
            to_vol=upper_headers[idx].volumes[0],
            flow_area=float(config.header_flow_area_m2),
            k_loss=float(config.connector_k_loss),
            hydraulic_diam=float(config.header_inner_diameter_m),
        )
        _set_design_flow(junc, half_flow)
        junctions.append(junc)
    for suffix, idx in (("A", inlet_idx_a), ("B", inlet_idx_b)):
        junc = FlowJunction(
            name=f"J_LowerHeader_{suffix}_to_RadiatorInnerHeader",
            from_vol=lower_headers[idx].volumes[0],
            to_vol=inner_header,
            flow_area=float(config.header_flow_area_m2),
            k_loss=float(config.connector_k_loss),
            hydraulic_diam=float(config.header_inner_diameter_m),
        )
        _set_design_flow(junc, half_flow)
        junctions.append(junc)

    for i in range(int(config.n_tubes)):
        j = (i + 1) % int(config.n_tubes)
        upper_ring = FlowJunction(
            name=f"J_RadiatorUpperRing_{i + 1:02d}_to_{j + 1:02d}",
            from_vol=upper_headers[i].volumes[0],
            to_vol=upper_headers[j].volumes[0],
            flow_area=float(config.header_flow_area_m2),
            k_loss=float(config.radiator_header_k_loss),
            custom_length=float(config.upper_header_segment_length_m),
            hydraulic_diam=float(config.header_inner_diameter_m),
        )
        lower_ring = FlowJunction(
            name=f"J_RadiatorLowerRing_{i + 1:02d}_to_{j + 1:02d}",
            from_vol=lower_headers[i].volumes[0],
            to_vol=lower_headers[j].volumes[0],
            flow_area=float(config.header_flow_area_m2),
            k_loss=float(config.radiator_header_k_loss),
            custom_length=float(config.lower_header_segment_length_m),
            hydraulic_diam=float(config.header_inner_diameter_m),
        )
        to_tube = FlowJunction(
            name=f"J_RadiatorUpper_to_Tube_{i + 1:02d}",
            from_vol=upper_headers[i].volumes[0],
            to_vol=tube_channels[i].volumes[0],
            flow_area=float(config.tube_flow_area_m2),
            k_loss=float(config.radiator_tube_inlet_k_loss),
            custom_length=0.5 * float(config.tube_length_m) / int(config.n_axial),
            hydraulic_diam=float(config.tube_inner_diameter_m),
        )
        from_tube = FlowJunction(
            name=f"J_RadiatorTube_{i + 1:02d}_to_Lower",
            from_vol=tube_channels[i].volumes[-1],
            to_vol=lower_headers[i].volumes[0],
            flow_area=float(config.tube_flow_area_m2),
            k_loss=float(config.radiator_tube_outlet_k_loss),
            custom_length=0.5 * float(config.tube_length_m) / int(config.n_axial),
            hydraulic_diam=float(config.tube_inner_diameter_m),
        )
        _set_design_flow(upper_ring, half_flow)
        _set_design_flow(lower_ring, half_flow)
        _set_design_flow(to_tube, tube_flow)
        _set_design_flow(from_tube, tube_flow)
        junctions.extend([upper_ring, lower_ring, to_tube, from_tube])

    cold_return_branches = []
    for idx in range(1, 4):
        branch = make_channel(
            name=f"ColdReturnBranch_{idx}",
            n_nodes=int(config.cold_return_branch_n_nodes),
            length=float(config.cold_return_branch_length_m),
            area=float(config.cold_return_branch_area_m2),
            dh=float(config.cold_return_branch_dh_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        extend_channel_objects(volumes, junctions, branch)
        _set_channel_design_flow(branch, cold_branch_flow)
        j_in = FlowJunction(
            name=f"J_PumpOutletDistributor_to_ColdReturnBranch_{idx}",
            from_vol=pump_outlet_distributor,
            to_vol=branch.volumes[0],
            flow_area=float(config.cold_return_branch_area_m2),
            k_loss=float(config.connector_k_loss),
            hydraulic_diam=float(config.cold_return_branch_dh_m),
        )
        j_out = FlowJunction(
            name=f"J_ColdReturnBranch_{idx}_to_CoreInletConnector",
            from_vol=branch.volumes[-1],
            to_vol=build["core_inlet_connector"],
            flow_area=float(config.cold_return_branch_area_m2),
            k_loss=float(config.connector_k_loss),
            hydraulic_diam=float(config.cold_return_branch_dh_m),
        )
        _set_design_flow(j_in, cold_branch_flow)
        _set_design_flow(j_out, cold_branch_flow)
        junctions.extend([j_in, j_out])
        cold_return_branches.append(branch)

    build["radiator_adapter_volumes"] = volumes
    build["radiator_adapter_junctions"] = junctions
    build["radiator_adapter_components"] = components
    build["radiator_inlet_distributor"] = inlet_distributor
    build["radiator_inner_header"] = inner_header
    build["radiator_outer_header"] = build["radiator_outlet_header"]
    build["pump_outlet_distributor"] = pump_outlet_distributor
    build["radiator_upper_headers"] = upper_headers
    build["radiator_lower_headers"] = lower_headers
    build["radiator_tube_channels"] = tube_channels
    build["radiator_units"] = radiator_units
    build["cold_return_branches"] = cold_return_branches
    build["single_radiator_tube_flow_design_kg_s"] = tube_flow
    build["cold_return_branch_flow_design_kg_s"] = cold_branch_flow
