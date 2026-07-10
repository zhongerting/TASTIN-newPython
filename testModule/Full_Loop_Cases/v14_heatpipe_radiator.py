"""Independent V14 heat-pipe radiator adapter for Full_Loop_Cases."""

from dataclasses import dataclass
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

from Components.RingHP import RingHP
from Materials.Solids.KHP import PotassiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidVolume,
    MacroFlowJunction,
)

from .common_flow_builder import extend_channel_objects, make_channel


SEGMENT_SPECS: Tuple[Tuple[str, str, str, Tuple[int, int, int]], ...] = (
    ("A1_I1_to_O1", "I1", "O1", (5, 6, 6)),
    ("A2_O1_to_I2", "O1", "I2", (5, 5, 6)),
    ("A3_I2_to_O2", "I2", "O2", (5, 6, 6)),
    ("A4_O2_to_I3", "O2", "I3", (5, 5, 6)),
    ("A5_I3_to_O3", "I3", "O3", (5, 6, 6)),
    ("A6_O3_to_I1", "O3", "I1", (5, 6, 6)),
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
    hp_l_eva_m: float = 0.100
    hp_l_aba_m: float = 0.0
    hp_l_con_m: float = 0.500
    hp_n_eva: int = 1
    hp_n_aba: int = 0
    hp_n_con: int = 12
    hp_n_wick: int = 1
    hp_n_wall: int = 2
    hp_porosity: float = 0.6
    fin_thickness_m: float = 0.0004
    fin_height_m: float = 0.020
    n_fin_height: int = 15
    hp_emissivity: float = 0.75
    fin_emissivity: float = 0.75
    hp_crossflow_c: float = 0.65
    hp_crossflow_k_cal: float = 1.0
    hp_crossflow_wake_factor: float = 1.0

    inlet_mix_volume_factor: float = 0.10
    symmetric_ring_multiplier: int = 2

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
    if int(config.symmetric_ring_multiplier) != 2:
        raise ValueError("V14 uses one explicit collector ring with symmetric_ring_multiplier=2.")


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
    fluid_channel: Any,
    solid_header: HeatConduction2D,
    hp_multipliers: Sequence[int],
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
        up_view_factor=0.0,
        down_view_factor=0.3,
        T_space=float(config.t_space_k),
        header_correlation_func=_lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=None,
        hp_crossflow_c=float(config.hp_crossflow_c),
        hp_crossflow_k_cal=float(config.hp_crossflow_k_cal),
        hp_crossflow_wake_factor=float(config.hp_crossflow_wake_factor),
    )


def attach_v14_heatpipe_radiator(build: Dict[str, Any], config: V14HeatPipeRadiatorConfig) -> None:
    """Attach the V14 heat-pipe radiator network to a common full-loop build."""
    _validate_config(config)

    material = build["core_inlet_connector"].material
    initial_p = float(build["core_inlet_connector"].P)
    initial_t = float(build["core_inlet_connector"].T)
    total_flow = float(build["total_flow_design_kg_s"])
    macro_branch_flow = total_flow / 3.0
    single_ring_branch_flow = macro_branch_flow / float(config.symmetric_ring_multiplier)

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
    inlet_mix_nodes = {
        key: _make_mix_node(
            name=f"InletMix_{key}",
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
            area=float(config.ring_area_m2),
            dh=float(config.ring_dh_m),
            length=inlet_mix_length,
        )
        for key in ("I1", "I2", "I3")
    }
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
    volumes.extend(inlet_mix_nodes.values())
    volumes.extend(outlet_mix_nodes.values())
    mix_nodes = {**inlet_mix_nodes, **outlet_mix_nodes}

    hot_to_ring = []
    for idx, key in enumerate(("I1", "I2", "I3")):
        junc = MacroFlowJunction(
            name=f"J_HotOutletBranch_{idx + 1}_to_InletMix_{key}",
            from_vol=hot_branches[idx].volumes[-1],
            to_vol=inlet_mix_nodes[key],
            macro_vol=hot_branches[idx].volumes[-1],
            multiplier=int(config.symmetric_ring_multiplier),
            flow_area=float(config.hot_branch_area_m2),
            k_loss=float(config.hot_branch_to_ring_k_loss),
            custom_length=float(config.hot_branch_length_m) / max(1, int(config.hot_branch_n_nodes)),
        )
        _set_design_flow(junc, single_ring_branch_flow)
        junctions.append(junc)
        hot_to_ring.append(junc)

    ring_sectors = []
    ring_solids = []
    ring_hps = []
    segment_entry_links = []
    segment_exit_links = []
    for sector_name, start_key, end_key, multipliers in SEGMENT_SPECS:
        channel = make_channel(
            name=f"{sector_name}_Channel",
            n_nodes=int(config.ring_sector_n_nodes),
            length=float(config.ring_sector_length_m),
            area=float(config.ring_area_m2),
            dh=float(config.ring_dh_m),
            material=material,
            initial_p=initial_p,
            initial_t=initial_t,
        )
        extend_channel_objects(volumes, junctions, channel)
        _set_channel_design_flow(channel, single_ring_branch_flow)
        solid = _make_ring_solid(f"{sector_name}_Solid", config)
        ring_hp = _make_ring_hp(
            name=f"{sector_name}_RingHP",
            fluid_channel=channel,
            solid_header=solid,
            hp_multipliers=multipliers,
            config=config,
        )
        components.append(ring_hp)
        ring_sectors.append(channel)
        ring_solids.append(solid)
        ring_hps.append(ring_hp)

        entry_link = FlowJunction(
            name=f"J_{start_key}_to_{sector_name}",
            from_vol=mix_nodes[start_key],
            to_vol=channel.volumes[0],
            flow_area=float(config.ring_area_m2),
            k_loss=float(config.ring_entry_k_loss),
            hydraulic_diam=float(config.ring_dh_m),
        )
        exit_link = FlowJunction(
            name=f"J_{sector_name}_to_{end_key}",
            from_vol=channel.volumes[-1],
            to_vol=mix_nodes[end_key],
            flow_area=float(config.ring_area_m2),
            k_loss=float(config.ring_exit_k_loss) + float(ring_hp.outlet_k_loss),
            hydraulic_diam=float(config.ring_dh_m),
            dynamic_loss_params=ring_hp.outlet_dynamic_loss_params,
        )
        _set_design_flow(entry_link, single_ring_branch_flow)
        _set_design_flow(exit_link, single_ring_branch_flow)
        junctions.extend([entry_link, exit_link])
        segment_entry_links.append(entry_link)
        segment_exit_links.append(exit_link)

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
        _set_channel_design_flow(manifold, single_ring_branch_flow)
        manifolds.append(manifold)

        to_manifold = FlowJunction(
            name=f"J_OutletMix_{key}_Manifold_{idx + 1}",
            from_vol=outlet_mix_nodes[key],
            to_vol=manifold.volumes[0],
            flow_area=float(config.manifold_area_m2),
            k_loss=float(config.outlet_mix_to_manifold_k_loss),
            hydraulic_diam=float(config.manifold_dh_m),
        )
        to_header = MacroFlowJunction(
            name=f"J_Manifold_{idx + 1}_to_RadiatorOutletHeader",
            from_vol=manifold.volumes[-1],
            to_vol=build["radiator_outlet_header"].volumes[0],
            macro_vol=build["radiator_outlet_header"].volumes[0],
            multiplier=int(config.symmetric_ring_multiplier),
            flow_area=float(config.manifold_area_m2),
            k_loss=float(config.manifold_to_outlet_header_k_loss),
            custom_length=float(config.manifold_lengths_m[idx]) / max(1, int(config.manifold_node_counts[idx])),
        )
        _set_design_flow(to_manifold, single_ring_branch_flow)
        _set_design_flow(to_header, single_ring_branch_flow)
        junctions.extend([to_manifold, to_header])
        outlet_mix_to_manifold.append(to_manifold)
        manifold_to_header.append(to_header)

    build["radiator_adapter_volumes"] = volumes
    build["radiator_adapter_junctions"] = junctions
    build["radiator_adapter_components"] = components
    build["hot_outlet_branches"] = hot_branches
    build["inlet_mix_nodes"] = inlet_mix_nodes
    build["outlet_mix_nodes"] = outlet_mix_nodes
    build["ring_sectors"] = ring_sectors
    build["ring_solids"] = ring_solids
    build["ring_hps"] = ring_hps
    build["ring_segment_entry_junctions"] = segment_entry_links
    build["ring_segment_exit_junctions"] = segment_exit_links
    build["hot_outlet_to_ring_junctions"] = hot_to_ring
    build["manifolds"] = manifolds
    build["outlet_mix_to_manifold_junctions"] = outlet_mix_to_manifold
    build["manifold_to_outlet_header_junctions"] = manifold_to_header
    build["single_ring_branch_flow_design_kg_s"] = single_ring_branch_flow
    build["macro_hot_branch_flow_design_kg_s"] = macro_branch_flow
