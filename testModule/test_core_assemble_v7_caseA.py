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

from profiler import TEASAProfiler
from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.BoundaryVolume import (
    IncompressibleBoundaryVolume,
    InletJunction,
)
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
    MacroFlowJunction,
    NonUniformIncompressibleFluidChannel,
)
from Materials.Solids.UO2 import UO2
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.ZrH import ZirconiumHydride
from Materials.Solids.BerylliumOxide import BerylliumOxide
from Materials.Solids.GasGaps import Xenon, Cesium, CarbonDioxide, Helium
from Materials.Fluids.Sodium import Sodium
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Components.TFEUnit import TFEUnit, TFEGeometry, TFEMeshParams, GapConfig
from Components.ReactorCore import (
    ReactorCore,
    GlobalAnnulusStructureConfig,
    GlobalGapStructureConfig,
)
from test_core_assemble_v5 import (
    build_axial_power_profile,
    build_global_moderator_meshes,
    build_ring_power_factors,
)
from test_core_assemble_v6 import get_time_dependent_dt_cap_v6


CASE_A_TOTAL_MULTIPLIER = 37.0
CASE_A_DESIGN_TOTAL_FLOW_KG_S = 1.30
CASE_A_DESIGN_CHANNEL_FLOW_KG_S = CASE_A_DESIGN_TOTAL_FLOW_KG_S / CASE_A_TOTAL_MULTIPLIER
CASE_A_RING_MULTIPLIERS = (1, 6, 12, 18)
CASE_A_TEC_RING_MULTIPLIERS = (1, 6, 12, 15)


class HeatCapacityScaledMaterial:
    def __init__(self, base_material, heat_capacity_scale: float):
        self.base_material = base_material
        self.heat_capacity_scale = float(heat_capacity_scale)
        self.name = f"{getattr(base_material, 'name', base_material.__class__.__name__)}_CpScale{self.heat_capacity_scale:g}"
        self.formula = getattr(base_material, "formula", "")

    def __getattr__(self, name):
        return getattr(self.base_material, name)

    def conductivity(self, T):
        return self.base_material.conductivity(T)

    def density(self, T):
        return self.base_material.density(T)

    def heat_capacity(self, T):
        return self.heat_capacity_scale * self.base_material.heat_capacity(T)

    def emissivity(self, T):
        return self.base_material.emissivity(T)


def _scale_solid_heat_capacity(material, heat_capacity_scale: float):
    if np.isclose(float(heat_capacity_scale), 1.0):
        return material
    return HeatCapacityScaledMaterial(material, heat_capacity_scale)


def _make_case_a_coolant(coolant_material: str):
    material_key = str(coolant_material).strip().lower()
    if material_key in {"sodium", "na"}:
        return Sodium(), "Sodium"
    if material_key in {
        "sodiumpotassium78",
        "sodium_potassium_78",
        "sodium-potassium78",
        "nak",
        "nak78",
        "na-k78",
    }:
        return SodiumPotassium78(), "SodiumPotassium78"
    raise ValueError(
        "coolant_material must be one of: Sodium, SodiumPotassium78, NaK, NaK78."
    )


def _extend_channel_objects(all_vols, all_juncs, channel):
    all_vols.extend(channel.volumes)
    all_juncs.extend(channel.internal_junctions)


def _pipe_delta_p(channel: IncompressibleFluidChannel) -> float:
    return float(channel.volumes[0].P - channel.volumes[-1].P)


def _remove_file_if_present(path: Optional[str]):
    if not path:
        return
    if os.path.exists(path):
        os.remove(path)


def build_v7_case_a_system(
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: Optional[float] = None,
    inlet_plenum_volume_m3: float = 5.0e-3,
    outlet_plenum_volume_m3: float = 5.0e-3,
    plenum_length_m: float = 0.05,
    pipe_n_nodes: int = 8,
    solid_heat_capacity_scale: float = 1.0,
    solid_heat_capacity_scale_scope: str = "all",
    global_outer_heat_capacity_scale: Optional[float] = None,
    ring_multipliers: Optional[Sequence[int]] = None,
    tec_ring_multipliers: Optional[Sequence[int]] = None,
    representative_names: Optional[Sequence[str]] = None,
    representative_ring_mapping: Optional[Dict[str, int]] = None,
    representative_power_factors: Optional[Dict[str, float]] = None,
    physical_ring_count: Optional[int] = None,
    core_name: str = "TASTIN_Core_V7_CaseA",
    coolant_material: str = "Sodium",
) -> Dict[str, Any]:
    if inlet_plenum_volume_m3 <= 0.0 or outlet_plenum_volume_m3 <= 0.0:
        raise ValueError("Plenum volumes must be positive.")
    if plenum_length_m <= 0.0:
        raise ValueError("plenum_length_m must be positive.")
    if pipe_n_nodes < 1:
        raise ValueError("pipe_n_nodes must be >= 1.")

    l_lower = 0.065
    l_active = 0.377
    l_upper = 0.065
    n_lower = 6
    n_active = 25
    n_upper = 6
    n_total = n_lower + n_active + n_upper
    total_height = l_lower + l_active + l_upper

    node_lengths = np.array(
        [l_lower / n_lower] * n_lower
        + [l_active / n_active] * n_active
        + [l_upper / n_upper] * n_upper,
        dtype=float,
    )
    common_y_faces = np.insert(np.cumsum(node_lengths), 0, 0.0)

    geom_params = TFEGeometry(
        r_pellet_inner=4.0e-3,
        r_pellet_outer=8.5e-3,
        r_fission_gas_outer=8.65e-3,
        r_emitter_outer=9.8e-3,
        r_collector_inner=10.3e-3,
        r_collector_outer=11.85e-3,
        r_inner_clad_inner=11.90e-3,
        r_inner_clad_outer=12.25e-3,
        r_coolant_inner=12.25e-3,
        r_coolant_outer=12.95e-3,
        r_outer_clad_outer=13.30e-3,
        r_moderator_inner=13.52e-3,
        r_moderator_outer=16.27e-3,
        height=total_height,
    )

    mesh_params = TFEMeshParams(
        n_axial=n_total,
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=1,
        n_r_outer_clad=1,
        n_r_moderator=3,
    )

    axial_power_profile = build_axial_power_profile(
        n_lower=n_lower,
        n_active=n_active,
        n_upper=n_upper,
    )

    base_solid_materials = {
        "UO2": UO2(),
        "MoNb": MoNb(),
        "Molybdenum": Molybdenum(),
        "StainlessSteel": AusteniticStainlessSteel(),
        "ZrH": ZirconiumHydride(),
        "BerylliumOxide": BerylliumOxide(),
    }
    scale_scope = str(solid_heat_capacity_scale_scope).lower()
    if scale_scope not in {"all", "global_outer", "tfe_only"}:
        raise ValueError(
            "solid_heat_capacity_scale_scope must be 'all', 'global_outer', or 'tfe_only'."
        )
    tfe_heat_capacity_scale = float(solid_heat_capacity_scale)
    global_heat_capacity_scale = (
        float(global_outer_heat_capacity_scale)
        if global_outer_heat_capacity_scale is not None
        else tfe_heat_capacity_scale
    )

    materials_dict = {
        key: _scale_solid_heat_capacity(value, tfe_heat_capacity_scale)
        if scale_scope in {"all", "tfe_only"}
        else value
        for key, value in base_solid_materials.items()
    }
    coolant, coolant_material_name = _make_case_a_coolant(coolant_material)
    materials_dict.update({
        "Sodium": coolant,
    })

    def global_material(key: str):
        if scale_scope == "global_outer":
            return _scale_solid_heat_capacity(
                base_solid_materials[key],
                tfe_heat_capacity_scale,
            )
        if scale_scope == "tfe_only":
            if global_outer_heat_capacity_scale is not None:
                return _scale_solid_heat_capacity(
                    base_solid_materials[key],
                    global_heat_capacity_scale,
                )
            return base_solid_materials[key]
        return materials_dict[key]

    sodium = materials_dict["Sodium"]

    cfg_fg = GapConfig(
        mode="simplified",
        h_eq=5678.0,
        material=Xenon(),
        emissivity_inner=0.15,
        emissivity_outer=0.15,
    )
    cfg_tec = GapConfig(
        mode="simplified",
        h_eq=29.0,
        material=Cesium(),
        emissivity_inner=0.15,
        emissivity_outer=0.60,
    )
    cfg_he = GapConfig(
        mode="simplified",
        h_eq=5678.0,
        material=Helium(),
        emissivity_inner=0.60,
        emissivity_outer=0.80,
    )
    cfg_co2 = GapConfig(
        mode="simplified",
        h_eq=53.6,
        material=CarbonDioxide(),
        emissivity_inner=0.80,
        emissivity_outer=0.80,
    )

    ring_names = (
        list(representative_names)
        if representative_names is not None
        else ["Center", "Ring1", "Ring2", "Ring3"]
    )
    if not ring_names or len(set(ring_names)) != len(ring_names):
        raise ValueError("representative_names must contain unique values.")
    multipliers = (
        list(ring_multipliers)
        if ring_multipliers is not None
        else list(CASE_A_RING_MULTIPLIERS)
    )
    if len(multipliers) != len(ring_names):
        raise ValueError("ring_multipliers must match representative_names length.")
    if any(int(mult) <= 0 for mult in multipliers):
        raise ValueError("ring_multipliers values must be positive.")
    multipliers = [int(mult) for mult in multipliers]
    tec_multipliers = (
        list(tec_ring_multipliers)
        if tec_ring_multipliers is not None
        else list(CASE_A_TEC_RING_MULTIPLIERS)
    )
    if len(tec_multipliers) != len(ring_names):
        raise ValueError("tec_ring_multipliers must match representative_names length.")
    if any(int(mult) < 0 for mult in tec_multipliers):
        raise ValueError("tec_ring_multipliers values must be non-negative.")
    tec_multipliers = [int(mult) for mult in tec_multipliers]
    for tfe_mult, tec_mult in zip(multipliers, tec_multipliers):
        if tec_mult > tfe_mult:
            raise ValueError("tec_ring_multipliers cannot exceed ring_multipliers.")
    total_multiplier = float(sum(multipliers))
    if channel_inlet_flow_kg_s is None:
        channel_inlet_flow_kg_s = CASE_A_DESIGN_TOTAL_FLOW_KG_S / total_multiplier

    p_inlet_sys = 165370.0
    p_outlet_sys = 161270.0
    w_tfe_single_design = float(channel_inlet_flow_kg_s)

    tfe_flow_area = np.pi * (geom_params.r_coolant_outer**2 - geom_params.r_coolant_inner**2)
    tfe_d_h = 2.0 * (geom_params.r_coolant_outer - geom_params.r_coolant_inner)

    total_flow_design = w_tfe_single_design * total_multiplier
    single_pipe_flow_design = total_flow_design / 3.0

    pipe_area_in = 597.9816e-6
    pipe_area_out = 5.9798e-4
    pipe_d_h = 27.6e-3
    l_inlet_pipe_1 = 1.89021
    l_inlet_pipe_23 = 2.50705
    l_outlet_pipe = 2.196

    all_fluid_vols = []
    all_fluid_juncs = []
    fluid_channels = {}

    inlet_boundary = IncompressibleBoundaryVolume(
        name="Global_Inlet_Boundary",
        material=sodium,
        P=p_inlet_sys,
        T=inlet_temperature_k,
        flow_area=pipe_area_in * 3.0,
        hydraulic_diam=pipe_d_h,
    )
    outlet_boundary = IncompressibleBoundaryVolume(
        name="Global_Outlet_Boundary",
        material=sodium,
        P=p_outlet_sys,
        T=inlet_temperature_k,
        flow_area=pipe_area_out * 3.0,
        hydraulic_diam=pipe_d_h,
    )
    inlet_boundary.is_pressure_boundary = True
    outlet_boundary.is_pressure_boundary = True

    inlet_plenum = IncompressibleFluidVolume(
        name="Global_Inlet_Plenum",
        volume=float(inlet_plenum_volume_m3),
        length=float(plenum_length_m),
        flow_area=max(float(inlet_plenum_volume_m3) / float(plenum_length_m), tfe_flow_area * total_multiplier),
        hydraulic_diam=tfe_d_h,
        material=sodium,
        initial_P=p_inlet_sys,
        initial_T=inlet_temperature_k,
    )
    outlet_plenum = IncompressibleFluidVolume(
        name="Global_Outlet_Plenum",
        volume=float(outlet_plenum_volume_m3),
        length=float(plenum_length_m),
        flow_area=max(float(outlet_plenum_volume_m3) / float(plenum_length_m), tfe_flow_area * total_multiplier),
        hydraulic_diam=tfe_d_h,
        material=sodium,
        initial_P=p_outlet_sys,
        initial_T=inlet_temperature_k,
    )

    inlet_pipe_1 = IncompressibleFluidChannel(
        name="InletPipe1",
        n_nodes=pipe_n_nodes,
        total_length=l_inlet_pipe_1,
        flow_area=pipe_area_in,
        hydraulic_diam=pipe_d_h,
        initial_P=p_inlet_sys,
        initial_T=inlet_temperature_k,
        material=sodium,
    )
    inlet_pipe_23 = IncompressibleFluidChannel(
        name="InletPipe23Rep",
        n_nodes=pipe_n_nodes,
        total_length=l_inlet_pipe_23,
        flow_area=pipe_area_in,
        hydraulic_diam=pipe_d_h,
        initial_P=p_inlet_sys,
        initial_T=inlet_temperature_k,
        material=sodium,
    )
    outlet_pipe = IncompressibleFluidChannel(
        name="OutletPipeRep",
        n_nodes=pipe_n_nodes,
        total_length=l_outlet_pipe,
        flow_area=pipe_area_out,
        hydraulic_diam=pipe_d_h,
        initial_P=p_outlet_sys,
        initial_T=inlet_temperature_k,
        material=sodium,
    )

    all_fluid_vols.extend([inlet_boundary, inlet_plenum, outlet_plenum, outlet_boundary])
    _extend_channel_objects(all_fluid_vols, all_fluid_juncs, inlet_pipe_1)
    _extend_channel_objects(all_fluid_vols, all_fluid_juncs, inlet_pipe_23)
    _extend_channel_objects(all_fluid_vols, all_fluid_juncs, outlet_pipe)

    j_inlet_pipe_1 = InletJunction(
        name="J_InletBoundary_to_InletPipe1",
        from_vol=inlet_boundary,
        to_vol=inlet_pipe_1.volumes[0],
        W_initial=single_pipe_flow_design,
    )
    j_inlet_pipe_23 = InletJunction(
        name="J_InletBoundary_to_InletPipe23Rep",
        from_vol=inlet_boundary,
        to_vol=inlet_pipe_23.volumes[0],
        W_initial=single_pipe_flow_design,
    )
    j_inlet_pipe_1_out = FlowJunction(
        name="J_InletPipe1_to_InletPlenum",
        from_vol=inlet_pipe_1.volumes[-1],
        to_vol=inlet_plenum,
        flow_area=pipe_area_in,
    )
    j_inlet_pipe_23_out = MacroFlowJunction(
        name="J_InletPipe23Rep_to_InletPlenum",
        from_vol=inlet_pipe_23.volumes[-1],
        to_vol=inlet_plenum,
        macro_vol=inlet_plenum,
        multiplier=2,
        flow_area=pipe_area_in,
    )
    j_inlet_pipe_1_out.W = single_pipe_flow_design
    j_inlet_pipe_23_out.W = single_pipe_flow_design
    all_fluid_juncs.extend([
        j_inlet_pipe_1,
        j_inlet_pipe_23,
        j_inlet_pipe_1_out,
        j_inlet_pipe_23_out,
    ])

    for name, mult in zip(ring_names, multipliers):
        chan = NonUniformIncompressibleFluidChannel(
            name=f"Chan_{name}",
            node_lengths=node_lengths,
            flow_area=tfe_flow_area,
            hydraulic_diam=tfe_d_h,
            initial_P=p_inlet_sys,
            initial_T=inlet_temperature_k,
            material=sodium,
        )
        fluid_channels[name] = chan
        _extend_channel_objects(all_fluid_vols, all_fluid_juncs, chan)

        j_macro_in = MacroFlowJunction(
            name=f"J_PlenumIn_{name}",
            from_vol=inlet_plenum,
            to_vol=chan.volumes[0],
            macro_vol=inlet_plenum,
            multiplier=mult,
            flow_area=tfe_flow_area,
        )
        j_macro_out = MacroFlowJunction(
            name=f"J_PlenumOut_{name}",
            from_vol=chan.volumes[-1],
            to_vol=outlet_plenum,
            macro_vol=outlet_plenum,
            multiplier=mult,
            flow_area=tfe_flow_area,
        )
        j_macro_in.W = w_tfe_single_design
        j_macro_out.W = w_tfe_single_design
        all_fluid_juncs.extend([j_macro_in, j_macro_out])

    j_outlet_pipe_in = MacroFlowJunction(
        name="J_OutletPlenum_to_OutletPipeRep",
        from_vol=outlet_plenum,
        to_vol=outlet_pipe.volumes[0],
        macro_vol=outlet_plenum,
        multiplier=3,
        flow_area=pipe_area_out,
    )
    j_outlet_pipe_out = FlowJunction(
        name="J_OutletPipeRep_to_OutletBoundary",
        from_vol=outlet_pipe.volumes[-1],
        to_vol=outlet_boundary,
        flow_area=pipe_area_out,
    )
    j_outlet_pipe_in.W = single_pipe_flow_design
    j_outlet_pipe_out.W = single_pipe_flow_design
    all_fluid_juncs.extend([j_outlet_pipe_in, j_outlet_pipe_out])

    tfes = {}
    for name in ring_names:
        tfes[name] = TFEUnit(
            name=name,
            geometry=geom_params,
            mesh_params=mesh_params,
            materials=materials_dict,
            coolant_channel=fluid_channels[name],
            fission_gas_config=cfg_fg,
            tec_gap_config=cfg_tec,
            he_gap_config=cfg_he,
            co2_gap_config=cfg_co2,
            power_fraction=1.0,
            axial_power_profile=axial_power_profile,
            axial_length_allocation=[l_lower, l_active, l_upper],
            axial_node_allocation=[n_lower, n_active, n_upper],
            axial_contact_resistance=0.0,
        )

    tfe_multipliers = {name: mult for name, mult in zip(ring_names, multipliers)}
    tfe_tec_multipliers = {name: mult for name, mult in zip(ring_names, tec_multipliers)}
    tfe_power_factors = (
        dict(representative_power_factors)
        if representative_power_factors is not None
        else build_ring_power_factors(ring_names, multipliers)
    )
    ring_mapping = (
        dict(representative_ring_mapping)
        if representative_ring_mapping is not None
        else {name: i for i, name in enumerate(ring_names)}
    )
    if set(ring_mapping) != set(ring_names):
        raise ValueError("representative_ring_mapping keys must match representative_names.")
    n_physical_rings = (
        int(physical_ring_count)
        if physical_ring_count is not None
        else len(ring_names)
    )
    if n_physical_rings <= 0:
        raise ValueError("physical_ring_count must be positive.")

    mod_meshes = build_global_moderator_meshes(
        inner_radius=geom_params.r_moderator_outer,
        outer_radius=60.0e-3,
        n_rings=n_physical_rings,
        y_faces=common_y_faces,
        height=geom_params.height,
        n_axial=n_total,
    )

    moderator_barrel_gap_cfg = GlobalGapStructureConfig(
        mode="simplified",
        width=5.0e-3,
        h_eq=0.0,
        emissivity_inner=0.8,
        emissivity_outer=0.8,
    )
    barrel_cfg = GlobalAnnulusStructureConfig(
        material=global_material("StainlessSteel"),
        inner_radius=65.0e-3,
        thickness=3.0e-3,
        n_radial=3,
        initial_temp=inlet_temperature_k,
        outer_surface_emissivity=0.05,
    )
    barrel_reflector_gap_cfg = GlobalGapStructureConfig(
        mode="simplified",
        width=2.0e-3,
        h_eq=0.0,
        emissivity_inner=0.8,
        emissivity_outer=0.8,
    )
    reflector_cfg = GlobalAnnulusStructureConfig(
        material=global_material("BerylliumOxide"),
        outer_radius=102.0e-3,
        n_radial=8,
        initial_temp=inlet_temperature_k,
        outer_surface_emissivity=0.2,
    )

    core = ReactorCore(
        name=core_name,
        tfe_dict=tfes,
        tfe_multipliers=tfe_multipliers,
        tec_multipliers=tfe_tec_multipliers,
        tfe_power_factors=tfe_power_factors,
        mod_meshes=mod_meshes,
        mod_material=global_material("ZrH"),
        ring_mapping=ring_mapping,
        barrel_config=barrel_cfg,
        reflector_config=reflector_cfg,
        moderator_barrel_gap_config=moderator_barrel_gap_cfg,
        barrel_reflector_gap_config=barrel_reflector_gap_cfg,
        T_space=200.0,
        alpha_tec=0.5,
        enable_tec_coupled=True,
    )

    core.setup_tec_circuit(
        mode_str="fixed_u",
        target_value=27.2,
        I_guess=284.0
    )

    hydraulic_net = HydraulicNetwork(
        volumes=all_fluid_vols,
        junctions=all_fluid_juncs,
        gravity_vector=0.0,
    )
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)

    return {
        "system": system,
        "core": core,
        "tfes": tfes,
        "ring_names": ring_names,
        "ring_multipliers": tfe_multipliers,
        "tec_ring_multipliers": tfe_tec_multipliers,
        "fluid_channels": fluid_channels,
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "inlet_plenum": inlet_plenum,
        "outlet_plenum": outlet_plenum,
        "inlet_pipe_1": inlet_pipe_1,
        "inlet_pipe_23": inlet_pipe_23,
        "outlet_pipe": outlet_pipe,
        "j_inlet_pipe_1": j_inlet_pipe_1,
        "j_inlet_pipe_23": j_inlet_pipe_23,
        "j_inlet_pipe_1_out": j_inlet_pipe_1_out,
        "j_inlet_pipe_23_out": j_inlet_pipe_23_out,
        "j_outlet_pipe_in": j_outlet_pipe_in,
        "j_outlet_pipe_out": j_outlet_pipe_out,
        "total_flow_design_kg_s": total_flow_design,
        "single_pipe_flow_design_kg_s": single_pipe_flow_design,
        "coolant_material": coolant_material_name,
        "solid_heat_capacity_scale": tfe_heat_capacity_scale,
        "solid_heat_capacity_scale_scope": scale_scope,
        "global_outer_heat_capacity_scale": global_heat_capacity_scale,
    }


def _case_a_flow_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    inlet_plenum = build["inlet_plenum"]
    outlet_plenum = build["outlet_plenum"]
    inlet_pipe_1 = build["inlet_pipe_1"]
    inlet_pipe_23 = build["inlet_pipe_23"]
    outlet_pipe = build["outlet_pipe"]
    j_inlet_pipe_1 = build["j_inlet_pipe_1"]
    j_inlet_pipe_23 = build["j_inlet_pipe_23"]
    j_inlet_pipe_1_out = build["j_inlet_pipe_1_out"]
    j_inlet_pipe_23_out = build["j_inlet_pipe_23_out"]
    j_outlet_pipe_in = build["j_outlet_pipe_in"]
    j_outlet_pipe_out = build["j_outlet_pipe_out"]
    fluid_channels = build["fluid_channels"]
    ring_multipliers = build["ring_multipliers"]

    inlet_pipe_23_macro_flow = float(j_inlet_pipe_23_out.W * 2.0)
    outlet_pipe_macro_flow = float(j_outlet_pipe_in.W * 3.0)
    tfe_single_flows = {
        name: float(channel.internal_junctions[0].W)
        for name, channel in fluid_channels.items()
    }
    tfe_macro_flows = {
        name: float(tfe_single_flows[name] * ring_multipliers[name])
        for name in fluid_channels
    }

    return {
        "inlet_pipe_1_boundary_flow_kg_s": float(j_inlet_pipe_1.W),
        "inlet_pipe_23_boundary_single_flow_kg_s": float(j_inlet_pipe_23.W),
        "inlet_pipe_1_to_plenum_flow_kg_s": float(j_inlet_pipe_1_out.W),
        "inlet_pipe_23_to_plenum_single_flow_kg_s": float(j_inlet_pipe_23_out.W),
        "inlet_pipe_23_to_plenum_macro_flow_kg_s": inlet_pipe_23_macro_flow,
        "inlet_total_macro_flow_kg_s": float(j_inlet_pipe_1_out.W + inlet_pipe_23_macro_flow),
        "outlet_pipe_single_flow_kg_s": float(j_outlet_pipe_in.W),
        "outlet_pipe_macro_flow_kg_s": outlet_pipe_macro_flow,
        "outlet_pipe_to_boundary_single_flow_kg_s": float(j_outlet_pipe_out.W),
        "inlet_plenum_pressure_pa": float(inlet_plenum.P),
        "outlet_plenum_pressure_pa": float(outlet_plenum.P),
        "core_plenum_delta_p_pa": float(inlet_plenum.P - outlet_plenum.P),
        "inlet_pipe_1_delta_p_pa": _pipe_delta_p(inlet_pipe_1),
        "inlet_pipe_23_delta_p_pa": _pipe_delta_p(inlet_pipe_23),
        "outlet_pipe_delta_p_pa": _pipe_delta_p(outlet_pipe),
        "tfe_single_flow_kg_s": tfe_single_flows,
        "tfe_macro_flow_kg_s": tfe_macro_flows,
        "tfe_total_macro_flow_kg_s": float(sum(tfe_macro_flows.values())),
    }


def _case_a_electric_diagnostics(core: ReactorCore) -> Dict[str, Any]:
    def _empty() -> Dict[str, Any]:
        return {
            "tec_total_voltage_v": None,
            "tec_total_current_a": None,
            "tec_total_electric_power_w": None,
            "tec_load_resistance_ohm": None,
            "tec_main_voltage_v": None,
            "tec_main_current_a": None,
            "tec_main_electric_power_w": None,
            "tec_reserved_parallel_voltage_v": None,
            "tec_reserved_parallel_current_a": None,
            "tec_reserved_parallel_electric_power_w": None,
        }

    def _power(results: Optional[Dict[str, Any]]) -> Optional[float]:
        if not results:
            return None
        return float(results.get("Uout", 0.0)) * float(results.get("Iout", 0.0))

    if hasattr(core, "get_tec_circuit_global_results"):
        group_results = core.get_tec_circuit_global_results()
        main_results = group_results.get("main")
        if not main_results:
            return _empty()
        reserved_results = group_results.get("reserved_parallel")
        main_voltage = float(main_results.get("Uout", 0.0))
        main_current = float(main_results.get("Iout", 0.0))
        main_power = _power(main_results)
        reserved_power = _power(reserved_results)
        total_power = float(main_power or 0.0) + float(reserved_power or 0.0)
        return {
            "tec_total_voltage_v": main_voltage,
            "tec_total_current_a": main_current,
            "tec_total_electric_power_w": total_power,
            "tec_load_resistance_ohm": float(main_results.get("Rload", 0.0)),
            "tec_main_voltage_v": main_voltage,
            "tec_main_current_a": main_current,
            "tec_main_electric_power_w": main_power,
            "tec_reserved_parallel_voltage_v": (
                None if not reserved_results else float(reserved_results.get("Uout", 0.0))
            ),
            "tec_reserved_parallel_current_a": (
                None if not reserved_results else float(reserved_results.get("Iout", 0.0))
            ),
            "tec_reserved_parallel_electric_power_w": reserved_power,
        }

    thermo_calc = getattr(core, "thermo_calc", None)
    if thermo_calc is None:
        return _empty()

    global_results = thermo_calc.get_global_results()
    if not global_results:
        return _empty()

    voltage = float(global_results.get("Uout", 0.0))
    current = float(global_results.get("Iout", 0.0))
    power = voltage * current
    return {
        "tec_total_voltage_v": voltage,
        "tec_total_current_a": current,
        "tec_total_electric_power_w": power,
        "tec_load_resistance_ohm": float(global_results.get("Rload", 0.0)),
        "tec_main_voltage_v": voltage,
        "tec_main_current_a": current,
        "tec_main_electric_power_w": power,
        "tec_reserved_parallel_voltage_v": None,
        "tec_reserved_parallel_current_a": None,
        "tec_reserved_parallel_electric_power_w": None,
    }


def _case_a_reset_design_flows_after_restart(build: Dict[str, Any]) -> None:
    """Apply current CaseA design flow targets after loading an old restart."""
    system = build["system"]
    fluid_net = system.fluid_solver
    single_pipe_flow = float(build["single_pipe_flow_design_kg_s"])
    total_multiplier = float(sum(build["ring_multipliers"].values()))
    tfe_single_flow = float(build["total_flow_design_kg_s"] / total_multiplier)

    pipe_junctions = [
        build["j_inlet_pipe_1"],
        build["j_inlet_pipe_23"],
        build["j_inlet_pipe_1_out"],
        build["j_inlet_pipe_23_out"],
        build["j_outlet_pipe_in"],
        build["j_outlet_pipe_out"],
    ]
    for channel_key in ("inlet_pipe_1", "inlet_pipe_23", "outlet_pipe"):
        pipe_junctions.extend(build[channel_key].internal_junctions)

    tfe_junctions = []
    for channel in build["fluid_channels"].values():
        tfe_junctions.extend(channel.internal_junctions)

    for junc in fluid_net.junctions_obj:
        if junc in pipe_junctions:
            junc.W = single_pipe_flow
            if hasattr(junc, "set_flow_rate"):
                junc.set_flow_rate(single_pipe_flow)
            elif hasattr(junc, "target_W"):
                junc.target_W = single_pipe_flow
        elif junc in tfe_junctions or junc.name.startswith("J_PlenumIn_") or junc.name.startswith("J_PlenumOut_"):
            junc.W = tfe_single_flow
            if hasattr(junc, "target_W"):
                junc.target_W = tfe_single_flow

    for idx, junc in enumerate(fluid_net.junctions_obj):
        fluid_net.W_vec[idx] = junc.W
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()

    if hasattr(fluid_net, "W_old"):
        fluid_net.W_old[:] = fluid_net.W_vec
    if hasattr(fluid_net, "W_iterate"):
        fluid_net.W_iterate[:] = fluid_net.W_vec
    fluid_net._refresh_cached_boundary_targets()


def run_test_v7_case_a_flow_only(
    run_duration_s: float = 0.01,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: Optional[float] = None,
    max_dt: float = 0.001,
    safety_factor: float = 20.0,
    pipe_n_nodes: int = 8,
    hydraulic_init_dt: float = 0.1,
    hydraulic_init_tol: float = 1.0e-7,
    hydraulic_init_max_iter: int = 500,
    ring_multipliers: Optional[Sequence[int]] = None,
) -> Dict[str, Any]:
    print("=== TASTIN V7 CaseA: inlet/outlet pipe flow-only test ===")

    build = build_v7_case_a_system(
        inlet_temperature_k=inlet_temperature_k,
        channel_inlet_flow_kg_s=channel_inlet_flow_kg_s,
        pipe_n_nodes=pipe_n_nodes,
        ring_multipliers=ring_multipliers,
    )
    system = build["system"]
    fluid_net = system.fluid_solver
    core = build["core"]

    core.point_reactor = None
    core.enable_tec_coupled = False
    core.update_neutronic_power(p_total=0.0, p_fiss=0.0, p_decay=0.0, alpha=1.0)

    hydraulic_initialized = fluid_net.initialize_hydraulics(
        dt=hydraulic_init_dt,
        tol=hydraulic_init_tol,
        max_iter=hydraulic_init_max_iter,
    )

    current_time = 0.0
    stop_time = current_time + float(run_duration_s)
    history_time = []
    history_diagnostics = []

    while current_time < stop_time:
        relative_time = current_time
        dt_cap = get_time_dependent_dt_cap_v6(relative_time, max_dt)
        dt = fluid_net.get_max_stable_dt(max_limit=dt_cap) / max(float(safety_factor), 1.0)
        dt = max(1.0e-4, min(dt, dt_cap))
        dt = min(dt, stop_time - current_time)
        if hasattr(fluid_net, "set_time"):
            fluid_net.set_time(current_time)
        fluid_net.step_Picard(
            dt,
            max_iter=100,
        )
        current_time += dt
        system.global_time = current_time
        history_time.append(current_time)
        history_diagnostics.append(_case_a_flow_diagnostics(build))

    final_diagnostics = _case_a_flow_diagnostics(build)
    final_summary = {
        "case_name": "V7_CaseA_FlowOnly",
        "hydraulic_initialized": bool(hydraulic_initialized),
        "final_time_s": current_time,
        "run_duration_s": float(run_duration_s),
        "core_point_reactor_is_none": core.point_reactor is None,
        "core_power_w": float(core.last_total_core_power),
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "single_pipe_flow_design_kg_s": float(build["single_pipe_flow_design_kg_s"]),
        **final_diagnostics,
    }

    print("CaseA flow-only run completed.")
    print(final_summary)

    return {
        "build": build,
        "system": system,
        "fluid_net": fluid_net,
        "core": core,
        "history_time": history_time,
        "history_diagnostics": history_diagnostics,
        "final_summary": final_summary,
    }


def run_test_v7_case_a_heated(
    run_duration_s: float = 3600.0,
    target_time_s: Optional[float] = None,
    total_power_w: float = 115000.0,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: Optional[float] = None,
    max_dt: float = 1.0,
    safety_factor: float = 20.0,
    pipe_n_nodes: int = 8,
    restart_file: Optional[str] = None,
    final_restart_file: str = "test_core_assemble_v7_caseA_heated_restart_t3600.npz",
    save_interval: float = 20.0,
    keep_only_latest_restart: bool = True,
    reset_design_flow_after_restart: bool = True,
    solid_heat_capacity_scale: float = 1.0,
    solid_heat_capacity_scale_scope: str = "all",
    global_outer_heat_capacity_scale: Optional[float] = None,
    use_time_dependent_dt_cap: bool = True,
    thermo_update_interval_s: Optional[float] = None,
    ring_multipliers: Optional[Sequence[int]] = None,
) -> Dict[str, Any]:
    print("=== TASTIN V7 CaseA: heated core long transient ===")

    build = build_v7_case_a_system(
        inlet_temperature_k=inlet_temperature_k,
        channel_inlet_flow_kg_s=channel_inlet_flow_kg_s,
        pipe_n_nodes=pipe_n_nodes,
        solid_heat_capacity_scale=solid_heat_capacity_scale,
        solid_heat_capacity_scale_scope=solid_heat_capacity_scale_scope,
        global_outer_heat_capacity_scale=global_outer_heat_capacity_scale,
        ring_multipliers=ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    tfes = build["tfes"]
    ring_names = build["ring_names"]
    if thermo_update_interval_s is not None:
        core.thermo_update_interval = float(thermo_update_interval_s)

    core.point_reactor = None
    core.enable_tec_coupled = True
    core.update_neutronic_power(
        p_total=float(total_power_w),
        p_fiss=float(total_power_w),
        p_decay=0.0,
        alpha=1.0,
    )
    system.initialize_system()
    if restart_file:
        if not os.path.exists(restart_file):
            raise FileNotFoundError(f"Restart file not found: {restart_file}")
        system.load_global_state(restart_file)
        core.point_reactor = None
        core.enable_tec_coupled = True
        core.update_neutronic_power(
            p_total=float(total_power_w),
            p_fiss=float(total_power_w),
            p_decay=0.0,
            alpha=1.0,
        )
        if reset_design_flow_after_restart:
            _case_a_reset_design_flows_after_restart(build)

    current_time = float(system.global_time)
    stop_time = float(target_time_s) if target_time_s is not None else current_time + float(run_duration_s)
    if stop_time < current_time:
        raise ValueError(
            f"target stop time {stop_time} is earlier than current time {current_time}."
        )
    next_save_time = current_time + float(save_interval) if save_interval > 0.0 else None
    latest_restart_file = restart_file

    history_time = []
    history_max_fuel = []
    history_inlet_plenum_temp = []
    history_outlet_plenum_temp = []
    history_core_delta_p = []
    history_total_macro_flow = []
    history_tfe_total_macro_flow = []
    history_center_outlet_temp = []
    history_total_voltage = []
    history_total_current = []
    history_total_electric_power = []

    while current_time < stop_time:
        relative_time = current_time
        core.update_neutronic_power(
            p_total=float(total_power_w),
            p_fiss=float(total_power_w),
            p_decay=0.0,
            alpha=1.0,
        )

        dt_cap = (
            get_time_dependent_dt_cap_v6(relative_time, max_dt)
            if use_time_dependent_dt_cap
            else float(max_dt)
        )
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=dt_cap,
            safety_factor=safety_factor,
        )
        dt = min(dt, stop_time - current_time)
        system.step(dt, inner_iter=1)
        current_time = float(system.global_time)

        flow_diag = _case_a_flow_diagnostics(build)
        max_fuel = max(float(np.max(tfe.solids["pellet"].T)) for tfe in tfes.values())
        history_time.append(current_time)
        history_max_fuel.append(max_fuel)
        history_inlet_plenum_temp.append(float(build["inlet_plenum"].T))
        history_outlet_plenum_temp.append(float(build["outlet_plenum"].T))
        history_core_delta_p.append(float(flow_diag["core_plenum_delta_p_pa"]))
        history_total_macro_flow.append(float(flow_diag["inlet_total_macro_flow_kg_s"]))
        history_tfe_total_macro_flow.append(float(flow_diag["tfe_total_macro_flow_kg_s"]))
        history_center_outlet_temp.append(float(tfes["Center"].coolant.volumes[-1].T))
        electric_diag = _case_a_electric_diagnostics(core)
        history_total_voltage.append(electric_diag["tec_total_voltage_v"])
        history_total_current.append(electric_diag["tec_total_current_a"])
        history_total_electric_power.append(electric_diag["tec_total_electric_power_w"])

        if len(history_time) % 100 == 0:
            if electric_diag["tec_total_electric_power_w"] is None:
                electric_text = "U=    N/A V | I=    N/A A | Pelec=    N/A kW"
            else:
                electric_text = (
                    f"U={electric_diag['tec_total_voltage_v']:8.3f} V | "
                    f"I={electric_diag['tec_total_current_a']:8.3f} A | "
                    f"Pelec={electric_diag['tec_total_electric_power_w'] / 1000.0:8.3f} kW"
                )
            print(
                f"t={current_time:9.3f} s | "
                f"P={float(total_power_w) / 1000.0:8.3f} kW | "
                f"TinPl={history_inlet_plenum_temp[-1]:8.2f} K | "
                f"ToutPl={history_outlet_plenum_temp[-1]:8.2f} K | "
                f"TfuelMax={max_fuel:8.2f} K | "
                f"dP={history_core_delta_p[-1]:9.2f} Pa | "
                f"{electric_text}"
            )

        if next_save_time is not None and current_time >= next_save_time:
            checkpoint_path = f"test_core_assemble_v7_caseA_heated_restart_t{int(next_save_time)}.npz"
            print(f"[Checkpoint] Saving restart file: {checkpoint_path}")
            system.save_global_state(checkpoint_path)
            if keep_only_latest_restart and latest_restart_file and latest_restart_file != checkpoint_path:
                _remove_file_if_present(latest_restart_file)
            latest_restart_file = checkpoint_path
            next_save_time += float(save_interval)

    final_diagnostics = _case_a_flow_diagnostics(build)
    final_electric_diagnostics = _case_a_electric_diagnostics(core)
    final_summary = {
        "case_name": "V7_CaseA_Heated",
        "restart_file": restart_file,
        "target_time_s": stop_time,
        "final_time_s": current_time,
        "run_duration_s": float(run_duration_s),
        "core_point_reactor_is_none": core.point_reactor is None,
        "core_power_w": float(core.last_total_core_power),
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "single_pipe_flow_design_kg_s": float(build["single_pipe_flow_design_kg_s"]),
        "solid_heat_capacity_scale": float(build["solid_heat_capacity_scale"]),
        "solid_heat_capacity_scale_scope": build["solid_heat_capacity_scale_scope"],
        "global_outer_heat_capacity_scale": float(build["global_outer_heat_capacity_scale"]),
        "use_time_dependent_dt_cap": bool(use_time_dependent_dt_cap),
        "thermo_update_interval_s": float(core.thermo_update_interval),
        "final_max_fuel_k": history_max_fuel[-1] if history_max_fuel else None,
        "final_inlet_plenum_k": float(build["inlet_plenum"].T),
        "final_outlet_plenum_k": float(build["outlet_plenum"].T),
        "final_center_outlet_k": float(tfes["Center"].coolant.volumes[-1].T),
        "final_restart_file": final_restart_file,
        **final_electric_diagnostics,
        **final_diagnostics,
    }

    if final_restart_file:
        print(f"[Final] Saving restart file: {final_restart_file}")
        system.save_global_state(final_restart_file)
        if keep_only_latest_restart and latest_restart_file and latest_restart_file != final_restart_file:
            _remove_file_if_present(latest_restart_file)
        latest_restart_file = final_restart_file

    print("CaseA heated run completed.")
    print(final_summary)

    return {
        "build": build,
        "system": system,
        "core": core,
        "tfes": tfes,
        "ring_names": ring_names,
        "history_time": history_time,
        "history_max_fuel": history_max_fuel,
        "history_inlet_plenum_temp": history_inlet_plenum_temp,
        "history_outlet_plenum_temp": history_outlet_plenum_temp,
        "history_core_delta_p": history_core_delta_p,
        "history_total_macro_flow": history_total_macro_flow,
        "history_tfe_total_macro_flow": history_tfe_total_macro_flow,
        "history_center_outlet_temp": history_center_outlet_temp,
        "history_total_voltage_v": history_total_voltage,
        "history_total_current_a": history_total_current,
        "history_total_electric_power_w": history_total_electric_power,
        "final_summary": final_summary,
    }


if __name__ == "__main__":
    run_test_v7_case_a_flow_only()
    TEASAProfiler.report()
