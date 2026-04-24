import json
import math
import os
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from profiler import TEASAProfiler
from Solvers.SystemManager import SystemManager
from Solvers.Couplers import FluidSolidCouple, SolidSolidCouple2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
from Solvers.Hydrodynamics.Components import FlowJunction, NonUniformIncompressibleFluidChannel
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.BerylliumOxide import BerylliumOxide
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel


AXIAL_SHAPE_COEFFS = np.array(
    [
        6.27905178e-02,
        -7.13913811e-02,
        1.35276842e-02,
        -1.02326367e-03,
        3.90936491e-05,
    ],
    dtype=float,
)


def coreinput_axial_shape(z: np.ndarray) -> np.ndarray:
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
    z_centers = np.linspace(-1.0, 1.0, n_active)
    active_profile = coreinput_axial_shape(z_centers)
    active_profile /= np.sum(active_profile)
    return np.concatenate((np.zeros(n_lower), active_profile, np.zeros(n_upper)))


def lyon_martinelli_correlation(Re: np.ndarray, Pr: np.ndarray, p_d_ratio: float) -> np.ndarray:
    pe = np.maximum(np.asarray(Re, dtype=float) * np.asarray(Pr, dtype=float), 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


def parse_env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    return default if raw is None else int(raw)


def parse_env_float(name: str, default: float) -> float:
    raw = os.environ.get(name)
    return default if raw is None else float(raw)


def parse_env_str(name: str, default: str) -> str:
    raw = os.environ.get(name)
    return default if raw is None else raw


@dataclass(frozen=True)
class LayerSpec:
    name: str
    material_factory: Any
    thickness_m: float
    n_radial: int


@dataclass
class BenchmarkConfig:
    mode: str
    n_steps: int
    dt: float
    print_every: int
    l_lower: float
    l_active: float
    l_upper: float
    n_lower: int
    n_active: int
    n_upper: int
    inlet_temperature_k: float
    env_temperature_k: float
    p_inlet_pa: float
    p_outlet_pa: float
    inlet_flow_guess_kg_s: float
    total_power_w: float
    source_primary_amp: float
    source_secondary_amp: float
    source_primary_period_s: float
    source_secondary_period_s: float
    emissivity: float
    outlet_k_losses: Tuple[float, ...]
    channel_power_factors: Tuple[float, ...]
    channel_names: Tuple[str, ...]
    layer_specs: Tuple[LayerSpec, ...]
    coolant_inner_radius_m: float
    coolant_outer_radius_m: float
    summary_json_path: Optional[str] = None


@dataclass
class ChannelAssembly:
    name: str
    channel: NonUniformIncompressibleFluidChannel
    inlet_junction: FlowJunction
    outlet_junction: FlowJunction
    solids: Dict[str, HeatConduction2D]
    fluid_coupler: FluidSolidCouple
    solid_couplers: List[SolidSolidCouple2D]
    radiation_bc: Any
    channel_power_base_w: float
    source_phase: float
    source_weight_flat: np.ndarray
    heated_layer_name: str

    def evaluate_source_power_distribution(self, t: float) -> np.ndarray:
        scale = (
            1.0
            + self.source_primary_amp * math.sin(2.0 * math.pi * t / self.source_primary_period_s + self.source_phase)
            + self.source_secondary_amp * math.sin(
                2.0 * math.pi * t / self.source_secondary_period_s + 0.5 * self.source_phase
            )
        )
        scale = max(scale, 0.15)
        return self.channel_power_base_w * scale * self.source_weight_flat

    source_primary_amp: float = 0.0
    source_secondary_amp: float = 0.0
    source_primary_period_s: float = 1.0
    source_secondary_period_s: float = 1.0


def build_mode_config(mode: str) -> BenchmarkConfig:
    mode_key = mode.strip().lower()
    default_print_every = 2
    if mode_key == "smoke":
        n_lower, n_active, n_upper = 3, 12, 3
        layer_radial = (4, 4, 5, 6)
        n_steps = 3
        dt = 0.02
        default_print_every = 1
    elif mode_key == "stability":
        n_lower, n_active, n_upper = 6, 25, 6
        layer_radial = (5, 6, 8, 10)
        n_steps = 200
        dt = 0.03
        default_print_every = 20
    elif mode_key == "stress":
        n_lower, n_active, n_upper = 8, 35, 8
        layer_radial = (7, 8, 10, 12)
        n_steps = 12
        dt = 0.02
        default_print_every = 2
    else:
        mode_key = "baseline"
        n_lower, n_active, n_upper = 6, 25, 6
        layer_radial = (5, 6, 8, 10)
        n_steps = 8
        dt = 0.03
        default_print_every = 2

    n_steps = parse_env_int("BENCH_N_STEPS", n_steps)
    dt = parse_env_float("BENCH_DT", dt)
    print_every = parse_env_int("BENCH_PRINT_EVERY", default_print_every)

    layer_specs = (
        LayerSpec("inner_liner", AusteniticStainlessSteel, 0.35e-3, layer_radial[0]),
        LayerSpec("transition_shell", Molybdenum, 0.70e-3, layer_radial[1]),
        LayerSpec("heated_shell", MoNb, 0.90e-3, layer_radial[2]),
        LayerSpec("outer_shield", BerylliumOxide, 2.07e-3, layer_radial[3]),
    )

    return BenchmarkConfig(
        mode=mode_key,
        n_steps=n_steps,
        dt=dt,
        print_every=print_every,
        l_lower=0.065,
        l_active=0.377,
        l_upper=0.065,
        n_lower=n_lower,
        n_active=n_active,
        n_upper=n_upper,
        inlet_temperature_k=parse_env_float("BENCH_INLET_T_K", 743.0),
        env_temperature_k=parse_env_float("BENCH_ENV_T_K", 250.0),
        p_inlet_pa=parse_env_float("BENCH_PIN_PA", 165370.0),
        p_outlet_pa=parse_env_float("BENCH_POUT_PA", 161270.0),
        inlet_flow_guess_kg_s=parse_env_float("BENCH_W_IN_GUESS", 0.0351),
        total_power_w=parse_env_float("BENCH_TOTAL_POWER_W", 20000.0),
        source_primary_amp=parse_env_float("BENCH_SOURCE_PRIMARY_AMP", 0.12),
        source_secondary_amp=parse_env_float("BENCH_SOURCE_SECONDARY_AMP", 0.04),
        source_primary_period_s=parse_env_float("BENCH_SOURCE_PRIMARY_PERIOD", 1.8),
        source_secondary_period_s=parse_env_float("BENCH_SOURCE_SECONDARY_PERIOD", 0.45),
        emissivity=parse_env_float("BENCH_EMISSIVITY", 0.82),
        outlet_k_losses=(0.84, 0.92, 1.00, 1.08),
        channel_power_factors=(0.92, 1.00, 1.08, 1.16),
        channel_names=("ChanA", "ChanB", "ChanC", "ChanD"),
        layer_specs=layer_specs,
        coolant_inner_radius_m=12.25e-3,
        coolant_outer_radius_m=12.95e-3,
        summary_json_path=os.environ.get("BENCH_SUMMARY_JSON"),
    )


def build_clustered_faces(inner_radius: float, thickness: float, n_cells: int) -> np.ndarray:
    eta = np.linspace(0.0, 1.0, n_cells + 1)
    cluster = 0.5 * (1.0 - np.cos(np.pi * eta))
    return inner_radius + thickness * cluster


def build_node_lengths(config: BenchmarkConfig) -> np.ndarray:
    return np.array(
        [config.l_lower / config.n_lower] * config.n_lower
        + [config.l_active / config.n_active] * config.n_active
        + [config.l_upper / config.n_upper] * config.n_upper,
        dtype=float,
    )


def create_layer_solid(
    solid_name: str,
    inner_radius: float,
    layer_spec: LayerSpec,
    y_faces: np.ndarray,
    initial_temp: float,
) -> HeatConduction2D:
    x_faces = build_clustered_faces(inner_radius, layer_spec.thickness_m, layer_spec.n_radial)
    mesh = Mesh2D(
        x_dim=layer_spec.thickness_m,
        n_x=layer_spec.n_radial,
        y_dim=float(y_faces[-1] - y_faces[0]),
        n_y=len(y_faces) - 1,
        geometry_type="cylindrical",
        inner_radius=inner_radius,
        x_faces=x_faces,
        y_faces=y_faces,
    )
    return HeatConduction2D(
        mesh=mesh,
        material=layer_spec.material_factory(),
        name=solid_name,
        initial_temp=initial_temp,
    )


def reset_layer_boundary_conditions(solid: HeatConduction2D, t_ref: float) -> None:
    for side in ("left", "right", "top", "bottom"):
        solid.boundaries[side].clear_conditions()
    for side in ("top", "bottom"):
        solid.boundaries[side].add_resistance_condition(T_ext=t_ref, R_ext=1.0e15)


def seed_initial_temperature_field(
    solid: HeatConduction2D,
    base_temperature: float,
    channel_index: int,
    layer_index: int,
) -> None:
    mesh = solid.mesh
    t_2d = solid.T.reshape(mesh.shape_nodes)
    r_norm = (mesh.x_centers - mesh.x_faces[0]) / max(mesh.x_faces[-1] - mesh.x_faces[0], 1e-12)
    z_norm = (mesh.y_centers - mesh.y_faces[0]) / max(mesh.y_faces[-1] - mesh.y_faces[0], 1e-12)
    radial_term = 18.0 * (1.0 - r_norm[:, None])
    axial_term = 12.0 * np.cos(np.pi * (z_norm[None, :] - 0.5))
    channel_bias = 6.0 * channel_index
    layer_bias = 4.0 * layer_index
    t_2d[:] = base_temperature + radial_term + axial_term + channel_bias + layer_bias
    solid.initialize_state()


def create_heated_layer_source(
    mesh: Mesh2D,
    channel_power_base_w: float,
    axial_profile: np.ndarray,
    source_phase: float,
    config: BenchmarkConfig,
) -> Tuple[Any, np.ndarray]:
    volumes_2d = mesh.geom_data.volumes.reshape(mesh.shape_nodes)
    weights_2d = volumes_2d * axial_profile[np.newaxis, :]
    weights_2d /= np.sum(weights_2d)
    weight_flat = weights_2d.reshape(-1)

    def source_callback(t: float, t_current: np.ndarray) -> np.ndarray:
        scale = (
            1.0
            + config.source_primary_amp * math.sin(2.0 * math.pi * t / config.source_primary_period_s + source_phase)
            + config.source_secondary_amp
            * math.sin(2.0 * math.pi * t / config.source_secondary_period_s + 0.5 * source_phase)
        )
        scale = max(scale, 0.15)
        return channel_power_base_w * scale * weight_flat

    return source_callback, weight_flat


def build_parallel_benchmark_system(config: BenchmarkConfig) -> Tuple[SystemManager, List[ChannelAssembly], Dict[str, Any]]:
    node_lengths = build_node_lengths(config)
    y_faces = np.insert(np.cumsum(node_lengths), 0, 0.0)
    axial_profile = build_axial_power_profile(config.n_lower, config.n_active, config.n_upper)

    sodium = Sodium()
    flow_area = np.pi * (config.coolant_outer_radius_m**2 - config.coolant_inner_radius_m**2)
    hydraulic_diam = 2.0 * (config.coolant_outer_radius_m - config.coolant_inner_radius_m)
    heated_perimeter = 2.0 * np.pi * config.coolant_inner_radius_m

    inlet_plenum = IncompressibleBoundaryVolume(
        name="Bench_Inlet",
        material=sodium,
        P=config.p_inlet_pa,
        T=config.inlet_temperature_k,
    )
    inlet_plenum.is_pressure_boundary = True

    outlet_plenum = IncompressibleBoundaryVolume(
        name="Bench_Outlet",
        material=sodium,
        P=config.p_outlet_pa,
        T=config.inlet_temperature_k,
    )
    outlet_plenum.is_pressure_boundary = True

    all_volumes = [inlet_plenum, outlet_plenum]
    all_junctions: List[FlowJunction] = []
    assemblies: List[ChannelAssembly] = []

    power_factors = np.asarray(config.channel_power_factors, dtype=float)
    power_factors /= np.sum(power_factors)

    for channel_index, channel_name in enumerate(config.channel_names):
        channel = NonUniformIncompressibleFluidChannel(
            name=channel_name,
            node_lengths=node_lengths,
            flow_area=flow_area,
            hydraulic_diam=hydraulic_diam,
            initial_P=config.p_inlet_pa,
            initial_T=config.inlet_temperature_k,
            material=sodium,
        )
        inlet_junction = FlowJunction(
            name=f"J_In_{channel_name}",
            from_vol=inlet_plenum,
            to_vol=channel.volumes[0],
            flow_area=flow_area,
        )
        outlet_junction = FlowJunction(
            name=f"J_Out_{channel_name}",
            from_vol=channel.volumes[-1],
            to_vol=outlet_plenum,
            flow_area=flow_area,
            k_loss=config.outlet_k_losses[channel_index],
        )

        all_volumes.extend(channel.volumes)
        all_junctions.extend(channel.internal_junctions)
        all_junctions.extend([inlet_junction, outlet_junction])

        solids: Dict[str, HeatConduction2D] = {}
        solid_couplers: List[SolidSolidCouple2D] = []

        current_inner_radius = config.coolant_inner_radius_m
        for layer_index, layer_spec in enumerate(config.layer_specs):
            solid = create_layer_solid(
                solid_name=f"{channel_name}_{layer_spec.name}",
                inner_radius=current_inner_radius,
                layer_spec=layer_spec,
                y_faces=y_faces,
                initial_temp=config.inlet_temperature_k,
            )
            reset_layer_boundary_conditions(solid, config.inlet_temperature_k)
            solids[layer_spec.name] = solid
            current_inner_radius += layer_spec.thickness_m

        radiation_bc = solids["outer_shield"].boundaries["right"].add_dynamic_radiation_condition(
            emissivity=config.emissivity,
            bare_area_array=solids["outer_shield"].boundaries["right"].area.copy(),
            T_env=config.env_temperature_k,
        )

        for idx in range(len(config.layer_specs) - 1):
            solid_left = solids[config.layer_specs[idx].name]
            solid_right = solids[config.layer_specs[idx + 1].name]
            solid_couplers.append(
                SolidSolidCouple2D(
                    obj1=solid_left,
                    obj2=solid_right,
                    direction="right",
                    contact_resistance=0.0,
                )
            )

        channel_power_base_w = config.total_power_w * power_factors[channel_index]
        source_phase = 0.45 * channel_index
        heated_layer = solids["heated_shell"]
        source_callback, source_weight_flat = create_heated_layer_source(
            mesh=heated_layer.mesh,
            channel_power_base_w=channel_power_base_w,
            axial_profile=axial_profile,
            source_phase=source_phase,
            config=config,
        )
        heated_layer.set_source_term(source_callback)

        for layer_index, layer_name in enumerate(solids.keys()):
            seed_initial_temperature_field(
                solid=solids[layer_name],
                base_temperature=config.inlet_temperature_k,
                channel_index=channel_index,
                layer_index=layer_index,
            )

        fluid_coupler = FluidSolidCouple(
            name=f"FluidSolid_{channel_name}",
            fluid=channel,
            solid_boundary_region=solids["inner_liner"].boundaries["left"],
            heated_perimeter=heated_perimeter,
            correlation_func=lyon_martinelli_correlation,
            solid_node_capacitance=solids["inner_liner"].get_boundary_node_capacitance("left"),
        )

        assemblies.append(
            ChannelAssembly(
                name=channel_name,
                channel=channel,
                inlet_junction=inlet_junction,
                outlet_junction=outlet_junction,
                solids=solids,
                fluid_coupler=fluid_coupler,
                solid_couplers=solid_couplers,
                radiation_bc=radiation_bc,
                channel_power_base_w=channel_power_base_w,
                source_phase=source_phase,
                source_weight_flat=source_weight_flat,
                heated_layer_name="heated_shell",
                source_primary_amp=config.source_primary_amp,
                source_secondary_amp=config.source_secondary_amp,
                source_primary_period_s=config.source_primary_period_s,
                source_secondary_period_s=config.source_secondary_period_s,
            )
        )

    hydraulic_net = HydraulicNetwork(
        volumes=all_volumes,
        junctions=all_junctions,
        gravity_vector=0.0,
    )
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)

    for assembly in assemblies:
        for solid in assembly.solids.values():
            system.add_solid_component(solid)
        system.add_coupler(assembly.fluid_coupler)
        for coupler in assembly.solid_couplers:
            system.add_coupler(coupler)

    reference = {
        "axial_profile": axial_profile,
        "node_lengths": node_lengths,
        "y_faces": y_faces,
        "flow_area": flow_area,
        "hydraulic_diam": hydraulic_diam,
        "heated_perimeter": heated_perimeter,
    }
    return system, assemblies, reference


def compute_total_solid_energy(assemblies: List[ChannelAssembly]) -> float:
    total = 0.0
    for assembly in assemblies:
        for solid in assembly.solids.values():
            total += float(np.sum(solid.thermal_capacitance * solid.T))
    return total


def compute_total_fluid_energy(assemblies: List[ChannelAssembly]) -> float:
    total = 0.0
    for assembly in assemblies:
        channel = assembly.channel
        temperatures = channel.temperature_vector
        pressures = channel.pressure_vector
        rho = channel.density_vector
        cp = channel.material.heat_capacity(temperatures, pressures)
        total += float(np.sum(rho * cp * channel.node_volume * temperatures))
    return total


def compute_total_radiation_power(assemblies: List[ChannelAssembly]) -> float:
    total = 0.0
    for assembly in assemblies:
        q_inflow = np.sum(assembly.solids["outer_shield"].boundaries["right"].current_flux)
        total += max(0.0, float(-q_inflow))
    return total


def compute_total_fluid_enthalpy_lift(assemblies: List[ChannelAssembly]) -> float:
    total = 0.0
    for assembly in assemblies:
        channel = assembly.channel
        mass_flow = float(assembly.inlet_junction.W)
        h_in = float(channel.volumes[0].h)
        h_out = float(channel.volumes[-1].h)
        total += mass_flow * (h_out - h_in)
    return total


def compute_total_source_power(assemblies: List[ChannelAssembly], t: float) -> float:
    total = 0.0
    for assembly in assemblies:
        total += float(np.sum(assembly.evaluate_source_power_distribution(t)))
    return total


def collect_grid_audit(assemblies: List[ChannelAssembly], node_lengths: np.ndarray) -> Dict[str, float]:
    audit: Dict[str, float] = {
        "axial_dz_min_m": float(np.min(node_lengths)),
        "axial_dz_max_m": float(np.max(node_lengths)),
        "axial_dz_ratio": float(np.max(node_lengths) / np.min(node_lengths)),
    }
    reference = assemblies[0]
    for layer_name, solid in reference.solids.items():
        radial_widths = np.diff(solid.mesh.x_faces)
        audit[f"{layer_name}.dr_min_m"] = float(np.min(radial_widths))
        audit[f"{layer_name}.dr_max_m"] = float(np.max(radial_widths))
        audit[f"{layer_name}.dr_ratio"] = float(np.max(radial_widths) / np.min(radial_widths))
        audit[f"{layer_name}.left_half_min_m"] = float(np.min(solid.mesh.dx_matrix[0, :]))
        audit[f"{layer_name}.left_half_max_m"] = float(np.max(solid.mesh.dx_matrix[0, :]))
        audit[f"{layer_name}.right_half_min_m"] = float(np.min(solid.mesh.dx_matrix[-1, :]))
        audit[f"{layer_name}.right_half_max_m"] = float(np.max(solid.mesh.dx_matrix[-1, :]))
        audit[f"{layer_name}.bottom_half_min_m"] = float(np.min(solid.mesh.dy_matrix[:, 0]))
        audit[f"{layer_name}.bottom_half_max_m"] = float(np.max(solid.mesh.dy_matrix[:, 0]))
        audit[f"{layer_name}.top_half_min_m"] = float(np.min(solid.mesh.dy_matrix[:, -1]))
        audit[f"{layer_name}.top_half_max_m"] = float(np.max(solid.mesh.dy_matrix[:, -1]))
    return audit


def collect_boundary_audit(assemblies: List[ChannelAssembly]) -> Dict[str, float]:
    audit: Dict[str, float] = {}
    reference = assemblies[0]
    for layer_name, solid in reference.solids.items():
        for side in ("left", "right", "top", "bottom"):
            boundary = solid.boundaries[side]
            audit[f"{layer_name}.{side}.R_int_min"] = float(np.min(boundary.R_internal))
            audit[f"{layer_name}.{side}.R_int_max"] = float(np.max(boundary.R_internal))
            audit[f"{layer_name}.{side}.T_surface_mean"] = float(np.mean(boundary.T_surface))
    return audit


def collect_channel_summary(assemblies: List[ChannelAssembly]) -> Dict[str, Dict[str, float]]:
    summary: Dict[str, Dict[str, float]] = {}
    for assembly in assemblies:
        heated_shell = assembly.solids["heated_shell"]
        outer_shield = assembly.solids["outer_shield"]
        summary[assembly.name] = {
            "mass_flow_kg_s": float(assembly.inlet_junction.W),
            "fluid_outlet_k": float(assembly.channel.volumes[-1].T),
            "inner_wall_mean_k": float(np.mean(assembly.solids["inner_liner"].boundaries["left"].T_surface)),
            "heated_layer_max_k": float(np.max(heated_shell.T)),
            "outer_wall_mean_k": float(np.mean(outer_shield.boundaries["right"].T_surface)),
            "radiation_power_w": max(0.0, float(-np.sum(outer_shield.boundaries["right"].current_flux))),
        }
    return summary


def collect_profiler_subset() -> Dict[str, Dict[str, float]]:
    keys = (
        "SystemManager.step",
        "FluidSolidCouple.execute",
        "SolidSolidCouple2D.sync",
        "BaseHeatConduction.get_derivatives",
        "HeatConduction2D._compute_internal_resistance",
        "HeatConduction2D._update_boundaries_state",
        "HeatConduction2D._compute_fluxes",
    )
    subset: Dict[str, Dict[str, float]] = {}
    for key in keys:
        if key in TEASAProfiler.stats:
            subset[key] = {
                "count": int(TEASAProfiler.stats[key]["count"]),
                "time_s": float(TEASAProfiler.stats[key]["time"]),
            }
    return subset


def run_benchmark(config: BenchmarkConfig) -> Dict[str, Any]:
    TEASAProfiler.stats = {}
    print(f"=== HeatConduction benchmark | mode={config.mode} ===")
    print("Parameter references: test_parallel_channels topology + test_core_assemble_v6 geometry scales")

    system, assemblies, reference = build_parallel_benchmark_system(config)
    grid_audit = collect_grid_audit(assemblies, reference["node_lengths"])

    init_start = time.perf_counter()
    system.initialize_system()
    init_elapsed = time.perf_counter() - init_start

    solid_energy_0 = compute_total_solid_energy(assemblies)
    fluid_energy_0 = compute_total_fluid_energy(assemblies)

    source_integral = 0.0
    radiation_integral = 0.0
    fluid_enthalpy_integral = 0.0
    prev_source = compute_total_source_power(assemblies, system.global_time)
    prev_radiation = compute_total_radiation_power(assemblies)
    prev_fluid_lift = compute_total_fluid_enthalpy_lift(assemblies)

    step_wall_times: List[float] = []
    history_time: List[float] = []
    history_heated_max: List[float] = []
    history_radiation: List[float] = []
    history_outlet_mean: List[float] = []

    for step_idx in range(1, config.n_steps + 1):
        step_start = time.perf_counter()
        system.step(config.dt, inner_iter=1)
        step_elapsed = time.perf_counter() - step_start
        step_wall_times.append(step_elapsed)

        current_time = float(system.global_time)
        current_source = compute_total_source_power(assemblies, current_time)
        current_radiation = compute_total_radiation_power(assemblies)
        current_fluid_lift = compute_total_fluid_enthalpy_lift(assemblies)

        source_integral += 0.5 * (prev_source + current_source) * config.dt
        radiation_integral += 0.5 * (prev_radiation + current_radiation) * config.dt
        fluid_enthalpy_integral += 0.5 * (prev_fluid_lift + current_fluid_lift) * config.dt
        prev_source = current_source
        prev_radiation = current_radiation
        prev_fluid_lift = current_fluid_lift

        channel_summary = collect_channel_summary(assemblies)
        heated_max = max(item["heated_layer_max_k"] for item in channel_summary.values())
        outlet_mean = float(np.mean([item["fluid_outlet_k"] for item in channel_summary.values()]))

        history_time.append(current_time)
        history_heated_max.append(heated_max)
        history_radiation.append(current_radiation)
        history_outlet_mean.append(outlet_mean)

        if step_idx % max(config.print_every, 1) == 0:
            print(
                f"step={step_idx:03d}/{config.n_steps:03d}"
                f" | t={current_time:7.3f} s"
                f" | step_wall={step_elapsed:8.4f} s"
                f" | T_heated_max={heated_max:8.2f} K"
                f" | T_out_mean={outlet_mean:8.2f} K"
                f" | Q_rad={current_radiation:10.3f} W"
            )

    solid_energy_1 = compute_total_solid_energy(assemblies)
    fluid_energy_1 = compute_total_fluid_energy(assemblies)
    energy_residual = (
        source_integral
        - radiation_integral
        - fluid_enthalpy_integral
        - (solid_energy_1 - solid_energy_0)
        - (fluid_energy_1 - fluid_energy_0)
    )

    boundary_audit = collect_boundary_audit(assemblies)
    channel_summary = collect_channel_summary(assemblies)
    profiler_subset = collect_profiler_subset()

    result = {
        "mode": config.mode,
        "n_steps": config.n_steps,
        "dt_s": config.dt,
        "simulated_time_s": config.n_steps * config.dt,
        "initialization_time_s": init_elapsed,
        "step_wall_time_total_s": float(np.sum(step_wall_times)),
        "step_wall_time_mean_s": float(np.mean(step_wall_times)),
        "step_wall_time_max_s": float(np.max(step_wall_times)),
        "heat_source_integral_j": float(source_integral),
        "radiation_integral_j": float(radiation_integral),
        "fluid_enthalpy_integral_j": float(fluid_enthalpy_integral),
        "solid_energy_change_proxy_j": float(solid_energy_1 - solid_energy_0),
        "fluid_energy_change_proxy_j": float(fluid_energy_1 - fluid_energy_0),
        "energy_residual_proxy_j": float(energy_residual),
        "grid_audit": grid_audit,
        "boundary_audit": boundary_audit,
        "channel_summary": channel_summary,
        "profiler": profiler_subset,
        "history_tail": {
            "time_s": history_time[-3:],
            "heated_layer_max_k": history_heated_max[-3:],
            "radiation_power_w": history_radiation[-3:],
            "fluid_outlet_mean_k": history_outlet_mean[-3:],
        },
        "parameter_reference": {
            "parallel_topology": "testModule/test_parallel_channels.py",
            "geometry_scale": "testModule/test_core_assemble_v6.py",
            "axial_profile_shape": "testModule/test_core_assemble_v5.py",
        },
    }

    print("\n--- Final Summary ---")
    print(json.dumps(result, indent=2, ensure_ascii=False))

    if config.summary_json_path:
        with open(config.summary_json_path, "w", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, ensure_ascii=False)
        print(f"Summary written to: {config.summary_json_path}")

    return result


if __name__ == "__main__":
    mode = parse_env_str("BENCH_MODE", "baseline")
    benchmark_config = build_mode_config(mode)
    run_benchmark(benchmark_config)
    print("\n--- TEASA Profiler ---")
    TEASAProfiler.report()
