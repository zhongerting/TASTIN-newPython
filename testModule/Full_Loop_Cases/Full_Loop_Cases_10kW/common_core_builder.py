from typing import Any, Dict, Sequence, Tuple

import numpy as np

from Components.ReactorCore import (
    GlobalAnnulusStructureConfig,
    GlobalGapStructureConfig,
    ReactorCore,
)
from Components.TFEUnit import GapConfig, TFEGeometry, TFEMeshParams, TFEUnit
from Materials.Fluids.Sodium import Sodium
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Materials.Solids.BerylliumOxide import BerylliumOxide
from Materials.Solids.GasGaps import CarbonDioxide, Cesium, Helium, Xenon
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.UO2 import UO2
from Materials.Solids.ZrH import ZirconiumHydride
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.Components import (
    MacroFlowJunction,
    NonUniformIncompressibleFluidChannel,
)

from .common_config import FullLoopCoreConfig, validate_sequence_lengths


CORE_RING_SHARES_RAW = np.array([0.019568969, 0.120310302, 0.180465534, 0.319655034, 0.360000072], dtype=float)
AXIAL_SHAPE_COEFFS = np.array(
    [5.0392372538e-02, -3.2174418071e-02, 6.847842042e-03, -4.513204066e-03, 2.890683804e-03],
    dtype=float,
)


class HeatCapacityScaledMaterial:
    def __init__(self, base_material: Any, heat_capacity_scale: float):
        self.base_material = base_material
        self.heat_capacity_scale = float(heat_capacity_scale)
        self.name = f"{getattr(base_material, 'name', base_material.__class__.__name__)}_CpScale{self.heat_capacity_scale:g}"
        self.formula = getattr(base_material, "formula", "")

    def __getattr__(self, name: str) -> Any:
        return getattr(self.base_material, name)

    def conductivity(self, temperature):
        return self.base_material.conductivity(temperature)

    def density(self, temperature):
        return self.base_material.density(temperature)

    def heat_capacity(self, temperature):
        return self.heat_capacity_scale * self.base_material.heat_capacity(temperature)

    def emissivity(self, temperature):
        return self.base_material.emissivity(temperature)


def _scale_solid_heat_capacity(material: Any, heat_capacity_scale: float) -> Any:
    if np.isclose(float(heat_capacity_scale), 1.0):
        return material
    return HeatCapacityScaledMaterial(material, heat_capacity_scale)


def make_coolant(coolant_material: str):
    key = str(coolant_material).strip().lower()
    if key in {"sodium", "na"}:
        return Sodium(), "Sodium"
    if key in {"sodiumpotassium78", "sodium_potassium_78", "sodium-potassium78", "nak", "nak78", "na-k78"}:
        return SodiumPotassium78(), "SodiumPotassium78"
    raise ValueError("coolant_material must be Sodium or SodiumPotassium78/NaK78.")


def _coreinput_axial_shape(z: np.ndarray) -> np.ndarray:
    z2 = np.asarray(z, dtype=float) ** 2
    shape = (
        AXIAL_SHAPE_COEFFS[0]
        + AXIAL_SHAPE_COEFFS[1] * z2
        + AXIAL_SHAPE_COEFFS[2] * z2**2
        + AXIAL_SHAPE_COEFFS[3] * z2**3
        + AXIAL_SHAPE_COEFFS[4] * z2**4
    )
    return np.maximum(shape, 0.0)


def build_axial_power_profile(n_lower: int, n_active: int, n_upper: int) -> np.ndarray:
    z_centers = np.linspace(-1.0, 1.0, int(n_active))
    active_profile = _coreinput_axial_shape(z_centers)
    active_profile /= float(np.sum(active_profile))
    return np.concatenate((np.zeros(int(n_lower)), active_profile, np.zeros(int(n_upper))))


def build_ring_power_factors(
    names: Sequence[str],
    multipliers: Sequence[int],
    ring_mapping: Dict[str, int],
) -> Dict[str, float]:
    shares = CORE_RING_SHARES_RAW / float(np.sum(CORE_RING_SHARES_RAW))
    ring_totals = {}
    for name, multiplier in zip(names, multipliers):
        ring_idx = int(ring_mapping[name])
        ring_totals[ring_idx] = ring_totals.get(ring_idx, 0) + int(multiplier)
    return {
        name: float(shares[int(ring_mapping[name])] / ring_totals[int(ring_mapping[name])])
        for name in names
    }


def build_global_moderator_meshes(
    inner_radius: float,
    outer_radius: float,
    n_rings: int,
    y_faces: np.ndarray,
    height: float,
    n_axial: int,
    radial_edges=None,
):
    meshes = []
    radial_edges = (
        np.asarray(radial_edges, dtype=float)
        if radial_edges is not None
        else np.linspace(float(inner_radius), float(outer_radius), int(n_rings) + 1)
    )
    if len(radial_edges) != int(n_rings) + 1:
        raise ValueError("radial_edges must contain n_rings + 1 values.")
    for r_in, r_out in zip(radial_edges[:-1], radial_edges[1:]):
        meshes.append(
            Mesh2D(
                x_dim=float(r_out - r_in),
                n_x=3,
                y_dim=float(height),
                n_y=int(n_axial),
                y_faces=y_faces,
                geometry_type="cylindrical",
                inner_radius=float(r_in),
            )
        )
    return meshes


def build_full_loop_core(
    core_config: FullLoopCoreConfig,
    *,
    core_inlet_connector: Any,
    core_outlet_connector: Any,
    total_flow_kg_s: float,
) -> Dict[str, Any]:
    names, multipliers, tec_multipliers = validate_sequence_lengths(core_config)
    if float(total_flow_kg_s) <= 0.0:
        raise ValueError("total_flow_kg_s must be positive.")

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

    geom = TFEGeometry(
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
    mesh = TFEMeshParams(
        n_axial=n_total,
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=1,
        n_r_outer_clad=1,
        n_r_moderator=3,
    )
    axial_power_profile = build_axial_power_profile(n_lower, n_active, n_upper)

    base_solid_materials = {
        "UO2": UO2(),
        "MoNb": MoNb(),
        "Molybdenum": Molybdenum(),
        "StainlessSteel": AusteniticStainlessSteel(),
        "ZrH": ZirconiumHydride(),
        "BerylliumOxide": BerylliumOxide(),
    }
    scope = str(core_config.solid_heat_capacity_scale_scope).lower()
    if scope not in {"all", "global_outer", "tfe_only"}:
        raise ValueError("solid_heat_capacity_scale_scope must be all, global_outer, or tfe_only.")
    tfe_cp_scale = float(core_config.solid_heat_capacity_scale)
    global_cp_scale = (
        float(core_config.global_outer_heat_capacity_scale)
        if core_config.global_outer_heat_capacity_scale is not None
        else tfe_cp_scale
    )
    tfe_materials = {
        key: _scale_solid_heat_capacity(value, tfe_cp_scale) if scope in {"all", "tfe_only"} else value
        for key, value in base_solid_materials.items()
    }
    coolant, coolant_name = make_coolant(core_config.coolant_material)
    tfe_materials["Sodium"] = coolant

    def global_material(key: str):
        if scope == "global_outer":
            return _scale_solid_heat_capacity(base_solid_materials[key], global_cp_scale)
        if scope == "tfe_only":
            if core_config.global_outer_heat_capacity_scale is not None:
                return _scale_solid_heat_capacity(base_solid_materials[key], global_cp_scale)
            return base_solid_materials[key]
        return tfe_materials[key]

    tec_gap_gas = str(core_config.tec_gap_gas).strip().lower()
    if tec_gap_gas == "cesium":
        tec_gap_material = Cesium()
    elif tec_gap_gas == "helium":
        tec_gap_material = Helium()
    else:
        raise ValueError("tec_gap_gas must be Cesium or Helium.")

    cfg_fg = GapConfig("simplified", 5678.0, Xenon(), 0.15, 0.15)
    cfg_tec = GapConfig(
        "simplified", float(core_config.tec_gap_h_eq_w_m2_k), tec_gap_material, 0.15, 0.60,
    )
    cfg_he = GapConfig("simplified", 5678.0, Helium(), 0.60, 0.80)
    cfg_co2 = GapConfig("simplified", 53.6, CarbonDioxide(), 0.80, 0.80)

    flow_area = np.pi * (geom.r_coolant_outer**2 - geom.r_coolant_inner**2)
    hydraulic_d = 2.0 * (geom.r_coolant_outer - geom.r_coolant_inner)
    single_tfe_flow = float(total_flow_kg_s) / float(sum(multipliers))
    fluid_channels = {}
    fluid_volumes = []
    fluid_junctions = []
    tfes = {}

    for name, multiplier in zip(names, multipliers):
        channel = NonUniformIncompressibleFluidChannel(
            name=f"Chan_{name}",
            node_lengths=node_lengths,
            flow_area=flow_area,
            hydraulic_diam=hydraulic_d,
            initial_P=float(core_config.reference_pressure_pa),
            initial_T=float(core_config.inlet_temperature_k),
            material=coolant,
        )
        fluid_channels[name] = channel
        fluid_volumes.extend(channel.volumes)
        fluid_junctions.extend(channel.internal_junctions)

        j_in = MacroFlowJunction(
            name=f"J_PlenumIn_{name}",
            from_vol=core_inlet_connector,
            to_vol=channel.volumes[0],
            macro_vol=core_inlet_connector,
            multiplier=int(multiplier),
            flow_area=flow_area,
        )
        j_out = MacroFlowJunction(
            name=f"J_PlenumOut_{name}",
            from_vol=channel.volumes[-1],
            to_vol=core_outlet_connector,
            macro_vol=core_outlet_connector,
            multiplier=int(multiplier),
            flow_area=flow_area,
        )
        j_in.W = single_tfe_flow
        j_out.W = single_tfe_flow
        fluid_junctions.extend([j_in, j_out])

        tfes[name] = TFEUnit(
            name=name,
            geometry=geom,
            mesh_params=mesh,
            materials=tfe_materials,
            coolant_channel=channel,
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

    ring_multipliers = {name: int(mult) for name, mult in zip(names, multipliers)}
    tec_map = {name: int(mult) for name, mult in zip(names, tec_multipliers)}
    power_factors = build_ring_power_factors(names, multipliers, core_config.representative_ring_mapping)
    mod_meshes = build_global_moderator_meshes(
        inner_radius=0.0,
        outer_radius=164.0e-3,
        n_rings=int(core_config.physical_ring_count),
        y_faces=common_y_faces,
        height=geom.height,
        n_axial=n_total,
        radial_edges=[0.0, 21.0e-3, 53.25e-3, 86.75e-3, 120.5e-3, 164.0e-3],
    )
    core = ReactorCore(
        name=str(core_config.core_name),
        tfe_dict=tfes,
        tfe_multipliers=ring_multipliers,
        tec_multipliers=tec_map,
        tfe_power_factors=power_factors,
        mod_meshes=mod_meshes,
        mod_material=global_material("ZrH"),
        ring_mapping=dict(core_config.representative_ring_mapping),
        barrel_config=GlobalAnnulusStructureConfig(
            material=global_material("StainlessSteel"),
            inner_radius=164.0e-3,
            outer_radius=166.0e-3,
            n_radial=3,
            initial_temp=float(core_config.inlet_temperature_k),
            outer_surface_emissivity=0.05,
        ),
        reflector_config=GlobalAnnulusStructureConfig(
            material=global_material("BerylliumOxide"),
            outer_radius=261.0e-3,
            n_radial=8,
            initial_temp=float(core_config.inlet_temperature_k),
            outer_surface_emissivity=0.2,
        ),
        moderator_barrel_gap_config=GlobalGapStructureConfig(
            mode="simplified",
            width=0.0,
            h_eq=5678.0,
            emissivity_inner=0.8,
            emissivity_outer=0.8,
        ),
        barrel_reflector_gap_config=GlobalGapStructureConfig(
            mode="simplified",
            width=0.0,
            h_eq=5678.0,
            emissivity_inner=0.8,
            emissivity_outer=0.8,
        ),
        T_space=200.0,
        alpha_tec=0.5,
        enable_tec_coupled=bool(core_config.main_tec_enabled),
        tec_lookup_enabled=core_config.tec_lookup_enabled,
        tec_lookup_db=core_config.tec_lookup_db,
        tec_lookup_regions=core_config.tec_lookup_regions,
    )
    if bool(core_config.main_tec_enabled):
        core.setup_tec_circuit(
            mode_str=str(core_config.main_tec_mode),
            target_value=float(core_config.main_tec_target_value),
            I_guess=float(core_config.main_tec_current_guess_a),
            topology=str(core_config.main_tec_topology),
        )
        reserved = core_config.reserved_parallel_tec
        if reserved.enabled:
            core.setup_reserved_parallel_tec_circuit(
                mode_str=str(reserved.mode),
                target_value=float(reserved.target_value),
                I_guess=float(reserved.current_guess_a),
                multipliers=reserved.multipliers,
            )

    return {
        "core": core,
        "tfes": tfes,
        "fluid_channels": fluid_channels,
        "core_fluid_volumes": fluid_volumes,
        "core_fluid_junctions": fluid_junctions,
        "ring_multipliers": ring_multipliers,
        "tec_ring_multipliers": tec_map,
        "single_tfe_flow_design_kg_s": single_tfe_flow,
        "coolant": coolant,
        "coolant_material": coolant_name,
    }
