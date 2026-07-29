"""Independent V14 heat-pipe radiator adapter for Full_Loop_Cases."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

from Components.ExternalHeatSources import OrbitalTableHeatSource, load_csv_flux_table_library
from Components.RadiatorThermalShield import RadiatorThermalShield
from Components.RingHP import RingHP
from Materials.Solids.KHP import PotassiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidVolume,
)

from .common_flow_builder import extend_channel_objects, make_channel


SEGMENT_PATH: Tuple[Tuple[str, str, str], ...] = (
    ("A1_I1_to_O1", "I1", "O1"),
    ("A2_O1_to_I2", "O1", "I2"),
    ("A3_I2_to_O2", "I2", "O2"),
    ("A4_O2_to_I3", "O2", "I3"),
    ("A5_I3_to_O3", "I3", "O3"),
    ("A6_O3_to_I1", "O3", "I1"),
)
UPPER_HP_MULTIPLIERS: Tuple[Tuple[int, int, int], ...] = ((8, 9, 9),) * 5 + ((8, 8, 8),)
LOWER_HP_MULTIPLIERS: Tuple[Tuple[int, int, int], ...] = ((10, 10, 11),) * 6
EXTERNAL_HEAT_CSV = (
    next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
    / 'Components'
    / 'ExternalHeatSources'
    / 'is58p5_w0_8p12_N18_sum.csv'
)

@dataclass(frozen=True)
class V14HeatPipeRadiatorConfig:
    """V14 heat-pipe radiator geometry and loss parameters."""

    hot_branch_length_m: float = 2.19632
    hot_branch_area_m2: float = float(np.pi * 0.0138**2)
    hot_branch_dh_m: float = 0.0276
    hot_branch_n_nodes: int = 8
    hot_branch_k_loss: float = 0.0
    hot_branch_to_ring_k_loss: float = 1.0

    ring_sector_length_m: float = 0.793
    ring_sector_n_nodes: int = 3
    ring_rect_length_m: float = 0.110
    ring_rect_width_m: float = 0.040
    ring_wall_thickness_m: float = 0.002
    ring_entry_k_loss: float = 0.5
    ring_exit_k_loss: float = 0.7
    ring_emissivity: float = 0.2

    manifold_lengths_m: Tuple[float, float, float] = (0.40911, 1.41912, 1.41912)
    manifold_node_counts: Tuple[int, int, int] = (5, 17, 17)
    manifold_area_m2: float = float(np.pi * 0.009**2)
    manifold_dh_m: float = 0.018
    outlet_mix_to_manifold_k_loss: float = 1.0
    manifold_to_outlet_header_k_loss: float = 1.1

    t_space_k: float = 200.0
    hp_initial_temp_k: float = 800.0
    hp_r_vapor_m: float = 0.0075
    hp_wick_thickness_m: float = 0.0005
    hp_wall_thickness_m: float = 0.0010
    hp_l_eva_m: float = 0.0605
    hp_l_aba_m: float = 0.0415
    hp_l_con_m: float = 0.47
    hp_n_eva: int = 1
    hp_n_aba: int = 1
    hp_n_con: int = 12
    hp_n_wick: int = 1
    hp_n_wall: int = 2
    hp_porosity: float = 0.6
    fin_thickness_m: float = 0.0004
    fin_height_m: float = 0.020
    n_fin_height: int = 15
    hp_emissivity: float = 0.75
    fin_emissivity: float = 0.75
    hp_up_view_factor: float = 0.0
    upper_hp_down_view_factor: float = 0.3
    lower_hp_down_view_factor: float = 0.3
    hp_crossflow_c: float = 0.65
    hp_crossflow_k_cal: float = 1.0
    hp_crossflow_wake_factor: float = 1.0
    external_heat_enabled: bool = False
    external_heat_period_s: float = 5668.14
    external_heat_time_origin_s: float = 0.0
    external_heat_scale_factor: float = 1.0
    external_heat_absorption_efficiency: float = 0.992

    thermal_shield_enabled: bool = False
    thermal_shield_initially_active: bool | None = None
    thermal_shield_active_until_s: float | None = None
    thermal_shield_model: str = "fortran_shield2"
    thermal_shield_view_factor: float = 0.8
    thermal_shield_inner_emissivity: float = 0.8
    thermal_shield_outer_emissivity: float = 0.1
    thermal_shield_conductivity_w_m_k: float = 0.0008
    thermal_shield_thickness_m: float = 0.01

    inlet_mix_volume_factor: float = 0.10

    @property
    def ring_area_m2(self) -> float:
        return float(self.ring_rect_length_m * self.ring_rect_width_m)

    @property
    def ring_perimeter_m(self) -> float:
        return float(2.0 * (self.ring_rect_length_m + self.ring_rect_width_m))

    @property
    def ring_dh_m(self) -> float:
        return float(4.0 * self.ring_area_m2 / self.ring_perimeter_m)

    @property
    def ring_inner_radius_m(self) -> float:
        return float(self.ring_perimeter_m / (2.0 * np.pi))

    @property
    def hp_r_in_m(self) -> float:
        return float(self.hp_r_vapor_m + self.hp_wick_thickness_m)

    @property
    def hp_r_out_m(self) -> float:
        return float(self.hp_r_in_m + self.hp_wall_thickness_m)

    @property
    def fin_wrap_ratio(self) -> float:
        return float((2.0 * self.fin_thickness_m) / (2.0 * np.pi * self.hp_r_out_m))


def _validate_config(config: V14HeatPipeRadiatorConfig) -> None:
    if int(config.ring_sector_n_nodes) != 3:
        raise ValueError("V14 ring_sector_n_nodes must be 3 to match the declared segment multipliers.")
    if len(config.manifold_lengths_m) != 3 or len(config.manifold_node_counts) != 3:
        raise ValueError("V14 requires exactly three manifold lengths and node counts.")
    for label, down_view in (
        ("upper", config.upper_hp_down_view_factor),
        ("lower", config.lower_hp_down_view_factor),
    ):
        view_sum = float(config.hp_up_view_factor) + float(down_view)
        if view_sum < -1.0e-12 or view_sum > 1.0 + 1.0e-12:
            raise ValueError(f"V14 {label} heat-pipe view-factor sum must be in [0, 1], got {view_sum}.")


def _external_heat_config_for_sector(library: Any, sector_index: int,
                                     config: V14HeatPipeRadiatorConfig) -> Dict[str, Any]:
    return {
        'use_embedded_table': True,
        'table_library': library,
        'table_ids_by_node': [3 * int(sector_index) + i for i in range(3)],
        'table_scale_factor': float(config.external_heat_scale_factor),
        'external_heat_absorption_efficiency': float(config.external_heat_absorption_efficiency),
        'table_periodic': True,
        'time_origin_s': float(config.external_heat_time_origin_s),
        'wall_illumination_factor': 0.5,
        'fin_illuminated_area_scale': 1.0,
        'fin_loading_mode': 'distributed_fin_absorption',
    }

def _lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    pe = np.maximum(np.asarray(Re, dtype=float) * np.asarray(Pr, dtype=float), 1.0)
    return 7.0 + 0.025 * (pe**0.8)


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


def _make_mix_node(
    *,
    name: str,
    material: Any,
    initial_p: float,
    initial_t: float,
    area: float,
    dh: float,
    length: float,
) -> IncompressibleFluidVolume:
    return IncompressibleFluidVolume(
        name=name,
        volume=float(area) * float(length),
        length=float(length),
        flow_area=float(area),
        hydraulic_diam=float(dh),
        initial_P=float(initial_p),
        initial_T=float(initial_t),
        material=material,
    )


def _make_ring_solid(name: str, config: V14HeatPipeRadiatorConfig) -> HeatConduction2D:
    inner_radius = config.ring_inner_radius_m
    mesh = Mesh2D(
        x_dim=float(config.ring_wall_thickness_m),
        n_x=1,
        y_dim=float(config.ring_sector_length_m),
        n_y=int(config.ring_sector_n_nodes),
        geometry_type="cylindrical",
        inner_radius=inner_radius,
    )
    solid = HeatConduction2D(
        mesh=mesh,
        material=SS316(name=f"{name}_SS316"),
        name=name,
        initial_temp=float(config.hp_initial_temp_k),
    )
    bare_area_array = solid.boundaries["right"].area / 2.0
    solid.boundaries["right"].add_dynamic_radiation_condition(
        emissivity=float(config.ring_emissivity),
        bare_area_array=bare_area_array,
        T_env=float(config.t_space_k),
    )
    return solid


def _make_ring_hp(
    *,
    name: str,
    prefix: str,
    fluid_channel: Any,
    solid_header: HeatConduction2D,
    hp_multipliers: Sequence[int],
    external_heat_config: Dict[str, Any] | None,
    config: V14HeatPipeRadiatorConfig,
) -> RingHP:
    wall = SS316(name=f"{name}_Wall")
    hp_fluid = PotassiumHP(name=f"{name}_HP_Fluid_K")
    wick = WickMaterial(
        name=f"{name}_Wick_K",
        solid_mat=wall,
        fluid_mat=hp_fluid,
        porosity=float(config.hp_porosity),
        r_vapor=float(config.hp_r_vapor_m),
        r_in_wall=float(config.hp_r_in_m),
    )
    return RingHP(
        name=name,
        fluid_channel=fluid_channel,
        solid_header=solid_header,
        hp_multipliers=list(hp_multipliers),
        header_flow_area=float(config.ring_area_m2),
        header_dh=float(config.ring_dh_m),
        header_heated_perimeter=float(config.ring_perimeter_m),
        hp_r_out=float(config.hp_r_out_m),
        hp_r_in=float(config.hp_r_in_m),
        hp_r_vapor=float(config.hp_r_vapor_m),
        hp_L_eva=float(config.hp_l_eva_m),
        hp_L_con=float(config.hp_l_con_m),
        hp_L_aba=float(config.hp_l_aba_m),
        hp_n_eva=int(config.hp_n_eva),
        hp_n_con=int(config.hp_n_con),
        hp_n_aba=int(config.hp_n_aba),
        hp_n_wick=int(config.hp_n_wick),
        hp_n_wall=int(config.hp_n_wall),
        porosity_hp=float(config.hp_porosity),
        HP_initial_temp=float(config.hp_initial_temp_k),
        hp_wall_mat=wall,
        hp_fluid_mat=hp_fluid,
        hp_wick_mat=wick,
        fin_thickness=float(config.fin_thickness_m),
        fin_height=float(config.fin_height_m),
        n_fin_height=int(config.n_fin_height),
        fin_wrap_ratio=float(config.fin_wrap_ratio),
        emissivity=float(config.hp_emissivity),
        fin_emissivity=float(config.fin_emissivity),
        up_view_factor=float(config.hp_up_view_factor),
        down_view_factor=float(
            config.upper_hp_down_view_factor
            if prefix == "Upper"
            else config.lower_hp_down_view_factor
        ),
        T_space=float(config.t_space_k),
        header_correlation_func=_lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=external_heat_config,
        hp_crossflow_c=float(config.hp_crossflow_c),
        hp_crossflow_k_cal=float(config.hp_crossflow_k_cal),
        hp_crossflow_wake_factor=float(config.hp_crossflow_wake_factor),
    )


def _make_ring_set(
    *,
    prefix: str,
    mix_nodes: Dict[str, Any],
    material: Any,
    initial_p: float,
    initial_t: float,
    branch_flow: float,
    hp_multipliers: Sequence[Tuple[int, int, int]],
    external_heat_library: Any,
    config: V14HeatPipeRadiatorConfig,
):
    ring_sectors = []
    ring_solids = []
    ring_hps = []
    components: List[Any] = []
    volumes: List[Any] = []
    junctions: List[Any] = []
    entry_links = []
    exit_links = []

    for sector_index, ((sector_name, start_key, end_key), multipliers) in enumerate(
            zip(SEGMENT_PATH, hp_multipliers)):
        full_name = f"{prefix}_{sector_name}"
        channel = make_channel(
            name=f"{full_name}_Channel",
            n_nodes=int(config.ring_sector_n_nodes),
            length=float(config.ring_sector_length_m),
            area=float(config.ring_area_m2),
            dh=float(config.ring_dh_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        extend_channel_objects(volumes, junctions, channel)
        _set_channel_design_flow(channel, branch_flow)
        solid = _make_ring_solid(f"{full_name}_Solid", config)
        ring_hp = _make_ring_hp(
            name=f"{full_name}_RingHP",
            prefix=prefix,
            fluid_channel=channel,            solid_header=solid,
            hp_multipliers=multipliers,
            external_heat_config=(
                _external_heat_config_for_sector(external_heat_library, sector_index, config)
                if external_heat_library is not None else None
            ),
            config=config,
        )
        components.append(ring_hp)
        ring_sectors.append(channel)
        ring_solids.append(solid)
        ring_hps.append(ring_hp)

        entry_link = FlowJunction(
            name=f"J_{prefix}_{start_key}_to_{sector_name}",
            from_vol=mix_nodes[start_key],
            to_vol=channel.volumes[0],
            flow_area=float(config.ring_area_m2),
            k_loss=float(config.ring_entry_k_loss),
            hydraulic_diam=float(config.ring_dh_m),
        )
        exit_link = FlowJunction(
            name=f"J_{prefix}_{sector_name}_to_{end_key}",
            from_vol=channel.volumes[-1],
            to_vol=mix_nodes[end_key],
            flow_area=float(config.ring_area_m2),
            k_loss=float(config.ring_exit_k_loss) + float(ring_hp.outlet_k_loss),
            hydraulic_diam=float(config.ring_dh_m),
            dynamic_loss_params=ring_hp.outlet_dynamic_loss_params,
        )
        _set_design_flow(entry_link, branch_flow)
        _set_design_flow(exit_link, branch_flow)
        junctions.extend([entry_link, exit_link])
        entry_links.append(entry_link)
        exit_links.append(exit_link)

    return {
        "volumes": volumes,
        "junctions": junctions,
        "components": components,
        "ring_sectors": ring_sectors,
        "ring_solids": ring_solids,
        "ring_hps": ring_hps,
        "entry_links": entry_links,
        "exit_links": exit_links,
    }


def attach_v14_heatpipe_radiator(build: Dict[str, Any], config: V14HeatPipeRadiatorConfig) -> None:
    """Attach the V14_10kW explicit upper/lower heat-pipe radiator network."""
    _validate_config(config)

    external_heat_library = None
    if config.external_heat_enabled:
        if float(config.external_heat_period_s) <= 0.0:
            raise ValueError('V14 external_heat_period_s must be positive.')
        external_heat_library = load_csv_flux_table_library(
            str(EXTERNAL_HEAT_CSV), float(config.external_heat_period_s)
        )
        if external_heat_library.available_ids() != tuple(range(18)):
            raise ValueError(f'V14 external heat CSV must contain 18 flux columns: {EXTERNAL_HEAT_CSV}')

    material = build["core_inlet_connector"].material
    initial_p = float(build["core_inlet_connector"].P)
    initial_t = float(build["core_inlet_connector"].T)
    total_flow = float(build["total_flow_design_kg_s"])
    macro_branch_flow = total_flow / 3.0
    upper_hp_count = sum(sum(values) for values in UPPER_HP_MULTIPLIERS)
    lower_hp_count = sum(sum(values) for values in LOWER_HP_MULTIPLIERS)
    total_hp_count = upper_hp_count + lower_hp_count
    upper_flow = total_flow * float(upper_hp_count) / float(total_hp_count)
    lower_flow = total_flow * float(lower_hp_count) / float(total_hp_count)
    upper_branch_flow = upper_flow / 3.0
    lower_branch_flow = lower_flow / 3.0

    volumes: List[Any] = []
    junctions: List[Any] = []
    components: List[Any] = []

    hot_branches = []
    for idx in range(3):
        branch = make_channel(
            name=f"HotOutletBranch_{idx + 1}",
            n_nodes=int(config.hot_branch_n_nodes),
            length=float(config.hot_branch_length_m),
            area=float(config.hot_branch_area_m2),
            dh=float(config.hot_branch_dh_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        extend_channel_objects(volumes, junctions, branch)
        _set_channel_design_flow(branch, macro_branch_flow)
        hot_branches.append(branch)

        split = FlowJunction(
            name=f"J_RadiatorInletHeader_to_HotOutletBranch_{idx + 1}",
            from_vol=build["radiator_inlet_header"].volumes[-1],
            to_vol=branch.volumes[0],
            flow_area=float(config.hot_branch_area_m2),
            k_loss=float(config.hot_branch_k_loss),
            hydraulic_diam=float(config.hot_branch_dh_m),
        )
        _set_design_flow(split, macro_branch_flow)
        junctions.append(split)

    ring_node_length = float(config.ring_sector_length_m) / float(config.ring_sector_n_nodes)
    inlet_mix_length = max(1.0e-6, float(config.inlet_mix_volume_factor) * ring_node_length)
    outlet_mix_length = max(1.0e-6, 2.0 * float(config.manifold_dh_m))

    def make_inlet_mix_nodes(prefix: str):
        inlet = {
            key: _make_mix_node(
                name=f"{prefix}_InletMix_{key}",
                material=material,
                initial_p=initial_p,
                initial_t=initial_t,
                area=float(config.ring_area_m2),
                dh=float(config.ring_dh_m),
                length=inlet_mix_length,
            )
            for key in ("I1", "I2", "I3")
        }
        volumes.extend(inlet.values())
        return inlet

    outlet_mix_nodes = {
        key: _make_mix_node(
            name=f"OutletMix_{key}",
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
            area=float(config.manifold_area_m2),
            dh=float(config.manifold_dh_m),
            length=outlet_mix_length,
        )
        for key in ("O1", "O2", "O3")
    }
    volumes.extend(outlet_mix_nodes.values())
    upper_mix_nodes = {**make_inlet_mix_nodes("Upper"), **outlet_mix_nodes}
    lower_mix_nodes = {**make_inlet_mix_nodes("Lower"), **outlet_mix_nodes}
    hot_to_ring = []
    for idx, key in enumerate(("I1", "I2", "I3")):
        for prefix, mix_nodes, flow in (
            ("Upper", upper_mix_nodes, upper_branch_flow),
            ("Lower", lower_mix_nodes, lower_branch_flow),
        ):
            junc = FlowJunction(
                name=f"J_HotOutletBranch_{idx + 1}_to_{prefix}_InletMix_{key}",
                from_vol=hot_branches[idx].volumes[-1],
                to_vol=mix_nodes[key],
                flow_area=float(config.hot_branch_area_m2),
                k_loss=float(config.hot_branch_to_ring_k_loss),
                custom_length=float(config.hot_branch_length_m) / max(1, int(config.hot_branch_n_nodes)),
                hydraulic_diam=float(config.hot_branch_dh_m),
            )
            _set_design_flow(junc, flow)
            junctions.append(junc)
            hot_to_ring.append(junc)

    upper = _make_ring_set(
        prefix="Upper",
        mix_nodes=upper_mix_nodes,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
        branch_flow=upper_branch_flow,
        hp_multipliers=UPPER_HP_MULTIPLIERS,
        external_heat_library=external_heat_library,
        config=config,
    )
    lower = _make_ring_set(
        prefix="Lower",
        mix_nodes=lower_mix_nodes,
        material=material,
        initial_p=initial_p,
        initial_t=initial_t,
        branch_flow=lower_branch_flow,
        hp_multipliers=LOWER_HP_MULTIPLIERS,
        external_heat_library=external_heat_library,
        config=config,
    )
    for ring in (upper, lower):
        volumes.extend(ring["volumes"])
        junctions.extend(ring["junctions"])
        components.extend(ring["components"])

    manifolds = []
    outlet_mix_to_manifold = []
    manifold_to_header = []
    for idx, key in enumerate(("O1", "O2", "O3")):
        manifold = make_channel(
            name=f"Manifold_{idx + 1}",
            n_nodes=int(config.manifold_node_counts[idx]),
            length=float(config.manifold_lengths_m[idx]),
            area=float(config.manifold_area_m2),
            dh=float(config.manifold_dh_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        extend_channel_objects(volumes, junctions, manifold)
        _set_channel_design_flow(manifold, macro_branch_flow)
        manifolds.append(manifold)

        to_manifold = FlowJunction(
            name=f"J_OutletMix_{key}_to_Manifold_{idx + 1}",
            from_vol=outlet_mix_nodes[key],
            to_vol=manifold.volumes[0],
            flow_area=float(config.manifold_area_m2),
            k_loss=float(config.outlet_mix_to_manifold_k_loss),
            hydraulic_diam=float(config.manifold_dh_m),
        )
        _set_design_flow(to_manifold, macro_branch_flow)
        junctions.append(to_manifold)
        outlet_mix_to_manifold.append(to_manifold)
        to_header = FlowJunction(
            name=f"J_Manifold_{idx + 1}_to_RadiatorOutletHeader",
            from_vol=manifold.volumes[-1],
            to_vol=build["radiator_outlet_header"].volumes[0],
            flow_area=float(config.manifold_area_m2),
            k_loss=float(config.manifold_to_outlet_header_k_loss),
            custom_length=float(config.manifold_lengths_m[idx]) / max(1, int(config.manifold_node_counts[idx])),
            hydraulic_diam=float(config.manifold_dh_m),
        )
        _set_design_flow(to_header, macro_branch_flow)
        junctions.append(to_header)
        manifold_to_header.append(to_header)

    ring_sectors = upper["ring_sectors"] + lower["ring_sectors"]
    ring_solids = upper["ring_solids"] + lower["ring_solids"]
    ring_hps = upper["ring_hps"] + lower["ring_hps"]
    radiator_units = [hp for ring_hp in ring_hps for hp in ring_hp.hp_units]
    thermal_shield = None
    if config.thermal_shield_enabled:
        direct_external_heat_sources = [
            unit.fin_external_heat_source.sources[0]
            for unit in radiator_units
            if getattr(unit, "fin_external_heat_source", None) is not None
        ]
        shield_external_heat_source = (
            OrbitalTableHeatSource(
                shape=(6,),
                table_ids=(0, 3, 6, 9, 12, 15),
                table_library=external_heat_library,
                scale_factor=float(config.external_heat_scale_factor),
                periodic=True,
                time_origin_s=float(config.external_heat_time_origin_s),
            )
            if external_heat_library is not None else None
        )
        thermal_shield = RadiatorThermalShield(
            name="V14_RadiatorThermalShield",
            radiator_units=radiator_units,
            active_until_s=config.thermal_shield_active_until_s,
            background_temperature_k=float(config.t_space_k),
            shield_view_factor=float(config.thermal_shield_view_factor),
            inner_emissivity=float(config.thermal_shield_inner_emissivity),
            outer_emissivity=float(config.thermal_shield_outer_emissivity),
            conductivity_w_m_k=float(config.thermal_shield_conductivity_w_m_k),
            thickness_m=float(config.thermal_shield_thickness_m),
            model=config.thermal_shield_model,
            external_heat_source=shield_external_heat_source,
            direct_external_heat_sources=direct_external_heat_sources,
            external_heat_absorption_factor=float(config.external_heat_absorption_efficiency),
            active_override=config.thermal_shield_initially_active,
        )
        components.insert(0, thermal_shield)

    build["radiator_adapter_volumes"] = volumes
    build["radiator_adapter_junctions"] = junctions
    build["radiator_adapter_components"] = components
    build["hot_outlet_branches"] = hot_branches
    build["upper_mix_nodes"] = upper_mix_nodes
    build["lower_mix_nodes"] = lower_mix_nodes
    build["outlet_mix_nodes"] = outlet_mix_nodes
    build["upper_ring_sectors"] = upper["ring_sectors"]
    build["lower_ring_sectors"] = lower["ring_sectors"]
    build["ring_sectors"] = ring_sectors
    build["ring_solids"] = ring_solids
    build["ring_hps"] = ring_hps
    build["radiator_units"] = radiator_units
    build["radiator_thermal_shield"] = thermal_shield
    build["upper_heatpipe_count"] = upper_hp_count
    build["lower_heatpipe_count"] = lower_hp_count
    build["upper_ring_segment_entry_junctions"] = upper["entry_links"]
    build["lower_ring_segment_entry_junctions"] = lower["entry_links"]
    build["ring_segment_entry_junctions"] = upper["entry_links"] + lower["entry_links"]
    build["ring_segment_exit_junctions"] = upper["exit_links"] + lower["exit_links"]
    build["hot_outlet_to_ring_junctions"] = hot_to_ring
    build["manifolds"] = manifolds
    build["outlet_mix_to_manifold_junctions"] = outlet_mix_to_manifold
    build["manifold_to_outlet_header_junctions"] = manifold_to_header
    build["upper_ring_branch_flow_design_kg_s"] = upper_branch_flow
    build["lower_ring_branch_flow_design_kg_s"] = lower_branch_flow
    build["upper_ring_total_flow_design_kg_s"] = upper_flow
    build["lower_ring_total_flow_design_kg_s"] = lower_flow
    build["macro_hot_branch_flow_design_kg_s"] = macro_branch_flow
    build["external_heat_enabled"] = bool(config.external_heat_enabled)
    build["external_heat_csv"] = str(EXTERNAL_HEAT_CSV) if config.external_heat_enabled else None
    build["external_heat_period_s"] = float(config.external_heat_period_s)
    build["external_heat_time_origin_s"] = float(config.external_heat_time_origin_s)
