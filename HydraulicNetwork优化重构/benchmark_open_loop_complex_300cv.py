import os
import sys
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from profiler import TEASAProfiler
from Solvers.Hydrodynamics.BoundaryVolume import (
    IncompressibleBoundaryVolume,
    InletJunction,
)
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
)
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager


def env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    return float(raw)


def env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    return int(raw)


def env_bool(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def circular_diameter_from_area(area: float) -> float:
    return float(np.sqrt(4.0 * area / np.pi))


def normalized_profile(n_nodes: int, tilt: float, center_bulge: float) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, n_nodes)
    profile = 1.0 + tilt * x + center_bulge * np.cos(np.pi * x)
    profile = np.clip(profile, 0.2, None)
    return profile / np.sum(profile)


def seed_channel_flow(channel: IncompressibleFluidChannel, mass_flow: float) -> None:
    for junc in channel.internal_junctions:
        junc.W = mass_flow
        junc.update_velocity()


def seeded_flow_junction(
    name: str,
    from_vol,
    to_vol,
    flow_area: float,
    k_loss: float,
    w_initial: float,
    custom_length: float = None,
) -> FlowJunction:
    junc = FlowJunction(
        name=name,
        from_vol=from_vol,
        to_vol=to_vol,
        flow_area=flow_area,
        k_loss=k_loss,
        custom_length=custom_length,
    )
    junc.W = w_initial
    junc.update_velocity()
    return junc


@dataclass
class BenchmarkConfig:
    inlet_temp: float = env_float("BENCH_T_INLET", 968.0)
    init_temp: float = env_float("BENCH_T_INIT", 863.0)
    outlet_pressure: float = env_float("BENCH_P_OUTLET", 1.61e5)
    inlet_pressure_guess: float = env_float("BENCH_P_INLET", 1.66e5)
    total_flow: float = env_float("BENCH_W_TOTAL", 2.2)
    gravity: float = env_float("BENCH_GRAVITY", 0.0)
    dt: float = env_float("BENCH_DT", 0.02)
    n_steps: int = env_int("BENCH_N_STEPS", 2000)
    inner_iter: int = env_int("BENCH_INNER_ITER", 1)
    print_every: int = env_int("BENCH_PRINT_EVERY", 200)
    init_dt: float = env_float("BENCH_DT_INIT", 0.05)
    init_tol: float = env_float("BENCH_INIT_TOL", 1.0e-5)
    init_max_iter: int = env_int("BENCH_INIT_MAX_ITER", 500)
    heat_scale: float = env_float("BENCH_HEAT_SCALE", 1.0)
    cooling_scale: float = env_float("BENCH_COOLING_SCALE", 1.0)
    heat_primary_amp: float = env_float("BENCH_HEAT_PRIMARY_AMP", 0.22)
    heat_secondary_amp: float = env_float("BENCH_HEAT_SECONDARY_AMP", 0.07)

    def __post_init__(self) -> None:
        if env_bool("BENCH_DOUBLE_ITER", False):
            self.inner_iter = 2


class SyntheticFluidSourceCoupler:
    def __init__(
        self,
        name: str,
        volumetric_loads: Sequence[Tuple[str, Sequence, np.ndarray, Dict[str, float]]],
        implicit_cooling_loads: Sequence[Tuple[str, Sequence, np.ndarray, np.ndarray, float]],
        time_getter: Callable[[], float] = None,
    ) -> None:
        self.name = name
        self.volumetric_loads = list(volumetric_loads)
        self.implicit_cooling_loads = list(implicit_cooling_loads)
        self.time_getter = time_getter

    def set_time_getter(self, time_getter: Callable[[], float]) -> None:
        self.time_getter = time_getter

    def current_time(self) -> float:
        if self.time_getter is None:
            return 0.0
        return float(self.time_getter())

    @staticmethod
    def _periodic_multiplier(time_value: float, spec: Dict[str, float]) -> float:
        bias = spec.get("bias", 1.0)
        primary_amp = spec.get("primary_amp", 0.0)
        primary_period = max(spec.get("primary_period", 1.0), 1.0e-12)
        primary_phase = spec.get("primary_phase", 0.0)
        secondary_amp = spec.get("secondary_amp", 0.0)
        secondary_period = max(spec.get("secondary_period", primary_period), 1.0e-12)
        secondary_phase = spec.get("secondary_phase", 0.0)
        min_scale = spec.get("min_scale", 0.4)
        max_scale = spec.get("max_scale", 1.6)

        scale = bias
        scale += primary_amp * np.sin(2.0 * np.pi * time_value / primary_period + primary_phase)
        scale += secondary_amp * np.sin(2.0 * np.pi * time_value / secondary_period + secondary_phase)
        return float(np.clip(scale, min_scale, max_scale))

    @TEASAProfiler.profile
    def execute(self) -> None:
        time_value = self.current_time()

        for _, volumes, q_vol_array, periodic_spec in self.volumetric_loads:
            heat_scale = self._periodic_multiplier(time_value, periodic_spec)
            for vol, q_vol in zip(volumes, q_vol_array):
                vol.Q_vol += float(q_vol) * heat_scale

        for _, volumes, explicit_array, implicit_array, _ in self.implicit_cooling_loads:
            for vol, q_explicit, lam in zip(volumes, explicit_array, implicit_array):
                vol.Q_wall += float(q_explicit)
                vol.implicit_coeff += float(lam)

    def base_heating_power(self) -> float:
        return float(sum(np.sum(q_vol_array) for _, _, q_vol_array, _ in self.volumetric_loads))

    def current_heating_power(self) -> float:
        time_value = self.current_time()
        total = 0.0
        for _, _, q_vol_array, periodic_spec in self.volumetric_loads:
            total += float(np.sum(q_vol_array)) * self._periodic_multiplier(time_value, periodic_spec)
        return total

    def total_cooling_conductance(self) -> float:
        return float(sum(np.sum(implicit_array) for _, _, _, implicit_array, _ in self.implicit_cooling_loads))

    def effective_cooling_power(self) -> float:
        total = 0.0
        for _, volumes, _, implicit_array, sink_temp in self.implicit_cooling_loads:
            for vol, lam in zip(volumes, implicit_array):
                total += float(lam) * (float(vol.T) - sink_temp)
        return total


def make_plenum(
    name: str,
    area: float,
    length: float,
    initial_pressure: float,
    initial_temp: float,
    material: SodiumPotassium78,
) -> IncompressibleFluidVolume:
    return IncompressibleFluidVolume(
        name=name,
        volume=area * length,
        length=length,
        flow_area=area,
        hydraulic_diam=circular_diameter_from_area(area),
        initial_P=initial_pressure,
        initial_T=initial_temp,
        material=material,
    )


def build_branch_module(
    branch_index: int,
    config: BenchmarkConfig,
    material: SodiumPotassium78,
    branch_flow_guess: float,
    junction_length: float,
    areas: Dict[str, float],
    lengths: Dict[str, float],
    losses: Dict[str, float],
):
    inlet_pipe = IncompressibleFluidChannel(
        name=f"Branch{branch_index}_InletPipe",
        n_nodes=8,
        total_length=lengths["branch_inlet"],
        flow_area=areas["branch_pipe"],
        hydraulic_diam=circular_diameter_from_area(areas["branch_pipe"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    upper_collector = IncompressibleFluidChannel(
        name=f"Branch{branch_index}_UpperCollector",
        n_nodes=16,
        total_length=lengths["upper_collector"],
        flow_area=areas["collector"],
        hydraulic_diam=circular_diameter_from_area(areas["collector"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    lower_collector = IncompressibleFluidChannel(
        name=f"Branch{branch_index}_LowerCollector",
        n_nodes=16,
        total_length=lengths["lower_collector"],
        flow_area=areas["collector"],
        hydraulic_diam=circular_diameter_from_area(areas["collector"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    outlet_pipe = IncompressibleFluidChannel(
        name=f"Branch{branch_index}_OutletPipe",
        n_nodes=8,
        total_length=lengths["branch_outlet"],
        flow_area=areas["branch_pipe"],
        hydraulic_diam=circular_diameter_from_area(areas["branch_pipe"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )

    for channel in (inlet_pipe, upper_collector, lower_collector, outlet_pipe):
        seed_channel_flow(channel, branch_flow_guess)

    j_inlet_to_upper = seeded_flow_junction(
        name=f"J_Branch{branch_index}_InletPipe_to_UpperCollector",
        from_vol=inlet_pipe.volumes[-1],
        to_vol=upper_collector.volumes[0],
        flow_area=min(areas["branch_pipe"], areas["collector"]),
        k_loss=losses["pipe_to_upper"],
        w_initial=branch_flow_guess,
        custom_length=junction_length,
    )
    j_upper_to_lower = seeded_flow_junction(
        name=f"J_Branch{branch_index}_UpperCollector_to_LowerCollector",
        from_vol=upper_collector.volumes[-1],
        to_vol=lower_collector.volumes[0],
        flow_area=areas["collector"],
        k_loss=losses["upper_to_lower"],
        w_initial=branch_flow_guess,
        custom_length=junction_length,
    )
    j_lower_to_outlet = seeded_flow_junction(
        name=f"J_Branch{branch_index}_LowerCollector_to_OutletPipe",
        from_vol=lower_collector.volumes[-1],
        to_vol=outlet_pipe.volumes[0],
        flow_area=min(areas["branch_pipe"], areas["collector"]),
        k_loss=losses["lower_to_pipe"],
        w_initial=branch_flow_guess,
        custom_length=junction_length,
    )

    return {
        "inlet_pipe": inlet_pipe,
        "upper_collector": upper_collector,
        "lower_collector": lower_collector,
        "outlet_pipe": outlet_pipe,
        "segment_junctions": [
            j_inlet_to_upper,
            j_upper_to_lower,
            j_lower_to_outlet,
        ],
    }


def build_synthetic_source_coupler(branches, radiator_1, radiator_2, config: BenchmarkConfig) -> SyntheticFluidSourceCoupler:
    branch_total_powers = config.heat_scale * np.array([30000.0, 33500.0, 28500.0, 32000.0], dtype=float)
    upper_fraction = np.array([0.57, 0.58, 0.56, 0.58], dtype=float)
    lower_fraction = 1.0 - upper_fraction
    primary_periods = np.array([6.0, 7.4, 9.2, 11.3], dtype=float)
    branch_phases = np.array([0.0, 0.7, 1.6, 2.5], dtype=float)

    volumetric_loads = []
    for idx, branch in enumerate(branches):
        upper_total = branch_total_powers[idx] * upper_fraction[idx]
        lower_total = branch_total_powers[idx] * lower_fraction[idx]

        upper_profile = normalized_profile(16, tilt=0.10 - 0.05 * idx, center_bulge=0.18)
        lower_profile = normalized_profile(16, tilt=-0.08 + 0.04 * idx, center_bulge=0.12)

        volumetric_loads.append(
            (
                f"Branch{idx + 1}_UpperCollector_Qvol",
                branch["upper_collector"].volumes,
                upper_total * upper_profile,
                {
                    "bias": 1.0,
                    "primary_amp": config.heat_primary_amp,
                    "primary_period": primary_periods[idx],
                    "primary_phase": branch_phases[idx],
                    "secondary_amp": config.heat_secondary_amp,
                    "secondary_period": primary_periods[idx] / 2.7,
                    "secondary_phase": branch_phases[idx] + 0.35,
                    "min_scale": 0.55,
                    "max_scale": 1.45,
                },
            )
        )
        volumetric_loads.append(
            (
                f"Branch{idx + 1}_LowerCollector_Qvol",
                branch["lower_collector"].volumes,
                lower_total * lower_profile,
                {
                    "bias": 1.0,
                    "primary_amp": 0.9 * config.heat_primary_amp,
                    "primary_period": 1.15 * primary_periods[idx],
                    "primary_phase": branch_phases[idx] + 0.55,
                    "secondary_amp": 0.85 * config.heat_secondary_amp,
                    "secondary_period": primary_periods[idx] / 3.2,
                    "secondary_phase": branch_phases[idx] + 1.10,
                    "min_scale": 0.55,
                    "max_scale": 1.45,
                },
            )
        )

    sink_temperature = 835.0
    radiator_1_ha = config.cooling_scale * 520.0
    radiator_2_ha = config.cooling_scale * 380.0
    radiator_1_profile = normalized_profile(20, tilt=-0.20, center_bulge=0.08)
    radiator_2_profile = normalized_profile(20, tilt=0.10, center_bulge=0.05)

    radiator_1_implicit = radiator_1_ha * radiator_1_profile
    radiator_2_implicit = radiator_2_ha * radiator_2_profile

    implicit_cooling_loads = [
        (
            "RadiatorMain1_ImplicitCooling",
            radiator_1.volumes,
            radiator_1_implicit * sink_temperature,
            radiator_1_implicit,
            sink_temperature,
        ),
        (
            "RadiatorMain2_ImplicitCooling",
            radiator_2.volumes,
            radiator_2_implicit * sink_temperature,
            radiator_2_implicit,
            sink_temperature,
        ),
    ]

    return SyntheticFluidSourceCoupler(
        name="SyntheticFluidSourceCoupler",
        volumetric_loads=volumetric_loads,
        implicit_cooling_loads=implicit_cooling_loads,
    )


def build_benchmark_model(config: BenchmarkConfig) -> Dict[str, object]:
    material = SodiumPotassium78()

    areas = {
        "buffer": 1.25e-3,
        "hot_leg": 7.50e-4,
        "header": 1.60e-3,
        "branch_pipe": 2.25e-4,
        "collector": 3.10e-4,
        "radiator": 9.50e-4,
        "return_leg": 7.00e-4,
    }
    lengths = {
        "buffer": 1.8,
        "hot_leg": 1.2,
        "plenum": 0.18,
        "branch_inlet": 0.8,
        "upper_collector": 1.6,
        "lower_collector": 1.6,
        "branch_outlet": 0.8,
        "radiator_1": 2.0,
        "radiator_2": 2.0,
        "return_leg": 1.8,
    }

    inlet_boundary = IncompressibleBoundaryVolume(
        name="Benchmark_InletBoundary",
        material=material,
        P=config.inlet_pressure_guess,
        T=config.inlet_temp,
    )
    outlet_boundary = IncompressibleBoundaryVolume(
        name="Benchmark_OutletBoundary",
        material=material,
        P=config.outlet_pressure,
        T=config.init_temp,
    )
    outlet_boundary.is_pressure_boundary = True

    inlet_buffer = IncompressibleFluidChannel(
        name="InletBuffer",
        n_nodes=18,
        total_length=lengths["buffer"],
        flow_area=areas["buffer"],
        hydraulic_diam=circular_diameter_from_area(areas["buffer"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    hot_leg = IncompressibleFluidChannel(
        name="HotLegMain",
        n_nodes=12,
        total_length=lengths["hot_leg"],
        flow_area=areas["hot_leg"],
        hydraulic_diam=circular_diameter_from_area(areas["hot_leg"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    distributor = make_plenum(
        name="DistributorPlenum",
        area=areas["header"],
        length=lengths["plenum"],
        initial_pressure=config.outlet_pressure,
        initial_temp=config.init_temp,
        material=material,
    )
    merge_header = make_plenum(
        name="MergeHeader",
        area=areas["header"],
        length=lengths["plenum"],
        initial_pressure=config.outlet_pressure,
        initial_temp=config.init_temp,
        material=material,
    )
    radiator_1 = IncompressibleFluidChannel(
        name="RadiatorMain1",
        n_nodes=20,
        total_length=lengths["radiator_1"],
        flow_area=areas["radiator"],
        hydraulic_diam=circular_diameter_from_area(areas["radiator"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    radiator_2 = IncompressibleFluidChannel(
        name="RadiatorMain2",
        n_nodes=20,
        total_length=lengths["radiator_2"],
        flow_area=areas["radiator"],
        hydraulic_diam=circular_diameter_from_area(areas["radiator"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    return_leg = IncompressibleFluidChannel(
        name="ReturnLeg",
        n_nodes=18,
        total_length=lengths["return_leg"],
        flow_area=areas["return_leg"],
        hydraulic_diam=circular_diameter_from_area(areas["return_leg"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )
    outlet_buffer = IncompressibleFluidChannel(
        name="OutletBuffer",
        n_nodes=18,
        total_length=lengths["buffer"],
        flow_area=areas["buffer"],
        hydraulic_diam=circular_diameter_from_area(areas["buffer"]),
        initial_P=config.outlet_pressure,
        initial_T=config.init_temp,
        material=material,
    )

    mainline_flow_guess = config.total_flow
    branch_flow_guess = config.total_flow / 4.0
    for channel in (inlet_buffer, hot_leg, radiator_1, radiator_2, return_leg, outlet_buffer):
        seed_channel_flow(channel, mainline_flow_guess)

    branch_loss_data = [
        {
            "distributor_to_branch": 2.10,
            "pipe_to_upper": 0.65,
            "upper_to_lower": 0.90,
            "lower_to_pipe": 0.70,
            "branch_to_merge": 1.80,
        },
        {
            "distributor_to_branch": 2.45,
            "pipe_to_upper": 0.80,
            "upper_to_lower": 1.00,
            "lower_to_pipe": 0.85,
            "branch_to_merge": 2.00,
        },
        {
            "distributor_to_branch": 2.85,
            "pipe_to_upper": 0.95,
            "upper_to_lower": 1.10,
            "lower_to_pipe": 0.95,
            "branch_to_merge": 2.25,
        },
        {
            "distributor_to_branch": 3.15,
            "pipe_to_upper": 1.05,
            "upper_to_lower": 1.20,
            "lower_to_pipe": 1.05,
            "branch_to_merge": 2.45,
        },
    ]

    branch_modules = []
    for idx, loss_data in enumerate(branch_loss_data, start=1):
        branch_modules.append(
            build_branch_module(
                branch_index=idx,
                config=config,
                material=material,
                branch_flow_guess=branch_flow_guess,
                junction_length=0.10,
                areas=areas,
                lengths=lengths,
                losses=loss_data,
            )
        )

    all_volumes = [
        inlet_boundary,
        *inlet_buffer.volumes,
        *hot_leg.volumes,
        distributor,
    ]

    for branch in branch_modules:
        all_volumes.extend(branch["inlet_pipe"].volumes)
        all_volumes.extend(branch["upper_collector"].volumes)
        all_volumes.extend(branch["lower_collector"].volumes)
        all_volumes.extend(branch["outlet_pipe"].volumes)

    all_volumes.extend(
        [
            merge_header,
            *radiator_1.volumes,
            *radiator_2.volumes,
            *return_leg.volumes,
            *outlet_buffer.volumes,
            outlet_boundary,
        ]
    )

    inlet_junction = InletJunction(
        name="J_InletBoundary_to_InletBuffer",
        from_vol=inlet_boundary,
        to_vol=inlet_buffer.volumes[0],
        W_initial=config.total_flow,
    )
    inlet_junction.update_velocity()

    j_inlet_buffer_to_hot_leg = seeded_flow_junction(
        name="J_InletBuffer_to_HotLegMain",
        from_vol=inlet_buffer.volumes[-1],
        to_vol=hot_leg.volumes[0],
        flow_area=min(areas["buffer"], areas["hot_leg"]),
        k_loss=0.45,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )
    j_hot_leg_to_distributor = seeded_flow_junction(
        name="J_HotLegMain_to_DistributorPlenum",
        from_vol=hot_leg.volumes[-1],
        to_vol=distributor,
        flow_area=min(areas["hot_leg"], areas["header"]),
        k_loss=0.70,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )

    distributor_to_branch = []
    branch_to_merge = []
    for idx, branch in enumerate(branch_modules, start=1):
        loss_data = branch_loss_data[idx - 1]
        j_dist_to_branch = seeded_flow_junction(
            name=f"J_DistributorPlenum_to_Branch{idx}_InletPipe",
            from_vol=distributor,
            to_vol=branch["inlet_pipe"].volumes[0],
            flow_area=min(areas["header"], areas["branch_pipe"]),
            k_loss=loss_data["distributor_to_branch"],
            w_initial=branch_flow_guess,
            custom_length=0.10,
        )
        j_branch_to_merge = seeded_flow_junction(
            name=f"J_Branch{idx}_OutletPipe_to_MergeHeader",
            from_vol=branch["outlet_pipe"].volumes[-1],
            to_vol=merge_header,
            flow_area=min(areas["header"], areas["branch_pipe"]),
            k_loss=loss_data["branch_to_merge"],
            w_initial=branch_flow_guess,
            custom_length=0.10,
        )
        distributor_to_branch.append(j_dist_to_branch)
        branch_to_merge.append(j_branch_to_merge)

    j_merge_to_radiator_1 = seeded_flow_junction(
        name="J_MergeHeader_to_RadiatorMain1",
        from_vol=merge_header,
        to_vol=radiator_1.volumes[0],
        flow_area=min(areas["header"], areas["radiator"]),
        k_loss=0.85,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )
    j_radiator_1_to_radiator_2 = seeded_flow_junction(
        name="J_RadiatorMain1_to_RadiatorMain2",
        from_vol=radiator_1.volumes[-1],
        to_vol=radiator_2.volumes[0],
        flow_area=areas["radiator"],
        k_loss=0.55,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )
    j_radiator_2_to_return = seeded_flow_junction(
        name="J_RadiatorMain2_to_ReturnLeg",
        from_vol=radiator_2.volumes[-1],
        to_vol=return_leg.volumes[0],
        flow_area=min(areas["radiator"], areas["return_leg"]),
        k_loss=0.60,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )
    j_return_to_outlet_buffer = seeded_flow_junction(
        name="J_ReturnLeg_to_OutletBuffer",
        from_vol=return_leg.volumes[-1],
        to_vol=outlet_buffer.volumes[0],
        flow_area=min(areas["return_leg"], areas["buffer"]),
        k_loss=0.40,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )
    outlet_junction = seeded_flow_junction(
        name="J_OutletBuffer_to_OutletBoundary",
        from_vol=outlet_buffer.volumes[-1],
        to_vol=outlet_boundary,
        flow_area=areas["buffer"],
        k_loss=0.30,
        w_initial=mainline_flow_guess,
        custom_length=0.10,
    )

    all_junctions = [
        inlet_junction,
        *inlet_buffer.internal_junctions,
        j_inlet_buffer_to_hot_leg,
        *hot_leg.internal_junctions,
        j_hot_leg_to_distributor,
    ]

    for branch, j_dist_to_branch, j_branch_to_merge in zip(branch_modules, distributor_to_branch, branch_to_merge):
        all_junctions.append(j_dist_to_branch)
        all_junctions.extend(branch["inlet_pipe"].internal_junctions)
        all_junctions.extend(branch["segment_junctions"][:1])
        all_junctions.extend(branch["upper_collector"].internal_junctions)
        all_junctions.extend(branch["segment_junctions"][1:2])
        all_junctions.extend(branch["lower_collector"].internal_junctions)
        all_junctions.extend(branch["segment_junctions"][2:])
        all_junctions.extend(branch["outlet_pipe"].internal_junctions)
        all_junctions.append(j_branch_to_merge)

    all_junctions.extend(
        [
            j_merge_to_radiator_1,
            *radiator_1.internal_junctions,
            j_radiator_1_to_radiator_2,
            *radiator_2.internal_junctions,
            j_radiator_2_to_return,
            *return_leg.internal_junctions,
            j_return_to_outlet_buffer,
            *outlet_buffer.internal_junctions,
            outlet_junction,
        ]
    )

    network = HydraulicNetwork(all_volumes, all_junctions, gravity_vector=config.gravity)
    system = SystemManager(fluid_network=network)
    source_coupler = build_synthetic_source_coupler(branch_modules, radiator_1, radiator_2, config)
    source_coupler.set_time_getter(lambda: system.global_time)
    system.add_coupler(source_coupler)

    return {
        "config": config,
        "material": material,
        "network": network,
        "system": system,
        "source_coupler": source_coupler,
        "all_volumes": all_volumes,
        "all_junctions": all_junctions,
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "inlet_buffer": inlet_buffer,
        "hot_leg": hot_leg,
        "distributor": distributor,
        "branches": branch_modules,
        "distributor_to_branch": distributor_to_branch,
        "branch_to_merge": branch_to_merge,
        "merge_header": merge_header,
        "radiator_1": radiator_1,
        "radiator_2": radiator_2,
        "return_leg": return_leg,
        "outlet_buffer": outlet_buffer,
        "inlet_junction": inlet_junction,
        "outlet_junction": outlet_junction,
    }


def print_model_summary(model: Dict[str, object]) -> None:
    config: BenchmarkConfig = model["config"]
    source_coupler: SyntheticFluidSourceCoupler = model["source_coupler"]
    print("=" * 78)
    print("HydraulicNetwork Open-Loop Benchmark")
    print("=" * 78)
    print(f"Fluid nodes      : {len(model['all_volumes'])}")
    print(f"Flow junctions   : {len(model['all_junctions'])}")
    print(f"Time step        : {config.dt:.5f} s")
    print(f"Benchmark steps  : {config.n_steps}")
    print(f"Inner iterations : {config.inner_iter}")
    print(f"Inlet temperature: {config.inlet_temp:.2f} K")
    print(f"Initial temp     : {config.init_temp:.2f} K")
    print(f"Outlet pressure  : {config.outlet_pressure:.2f} Pa")
    print(f"Total flow       : {config.total_flow:.4f} kg/s")
    print(f"Gravity          : {config.gravity:.3f} m/s^2")
    print(f"Base Q_vol       : {source_coupler.base_heating_power():.2f} W")
    print(f"Heat forcing     : periodic, multi-frequency")
    print(f"Primary amplitude: +/- {100.0 * config.heat_primary_amp:.1f} %")
    print(f"Secondary amplitude: +/- {100.0 * config.heat_secondary_amp:.1f} %")
    print(f"Total hA         : {source_coupler.total_cooling_conductance():.2f} W/K")
    print("=" * 78)


def print_progress(step_index: int, model: Dict[str, object], step_wall_time: float) -> None:
    inlet_junction: InletJunction = model["inlet_junction"]
    distributor_to_branch: Sequence[FlowJunction] = model["distributor_to_branch"]
    outlet_buffer: IncompressibleFluidChannel = model["outlet_buffer"]
    radiator_2: IncompressibleFluidChannel = model["radiator_2"]
    source_coupler: SyntheticFluidSourceCoupler = model["source_coupler"]

    branch_flows = np.array([float(j.W) for j in distributor_to_branch], dtype=float)
    total_flow = max(abs(float(inlet_junction.W)), 1.0e-12)
    split_pct = 100.0 * branch_flows / total_flow
    print(
        f"step={step_index:4d} | "
        f"wall={step_wall_time:.4f} s | "
        f"W_in={float(inlet_junction.W):.4f} kg/s | "
        f"Q_vol={source_coupler.current_heating_power():8.1f} W | "
        f"T_rad2_out={float(radiator_2.volumes[-1].T):.2f} K | "
        f"T_out={float(outlet_buffer.volumes[-1].T):.2f} K | "
        f"split=[{split_pct[0]:5.1f}, {split_pct[1]:5.1f}, {split_pct[2]:5.1f}, {split_pct[3]:5.1f}] %"
    )


def print_final_summary(model: Dict[str, object], total_wall_time: float) -> None:
    inlet_junction: InletJunction = model["inlet_junction"]
    outlet_junction: FlowJunction = model["outlet_junction"]
    distributor_to_branch: Sequence[FlowJunction] = model["distributor_to_branch"]
    branches = model["branches"]
    source_coupler: SyntheticFluidSourceCoupler = model["source_coupler"]
    inlet_buffer: IncompressibleFluidChannel = model["inlet_buffer"]
    hot_leg: IncompressibleFluidChannel = model["hot_leg"]
    radiator_1: IncompressibleFluidChannel = model["radiator_1"]
    radiator_2: IncompressibleFluidChannel = model["radiator_2"]
    return_leg: IncompressibleFluidChannel = model["return_leg"]
    outlet_buffer: IncompressibleFluidChannel = model["outlet_buffer"]
    config: BenchmarkConfig = model["config"]

    branch_flows = np.array([float(j.W) for j in distributor_to_branch], dtype=float)
    total_flow = max(abs(float(inlet_junction.W)), 1.0e-12)
    branch_split = 100.0 * branch_flows / total_flow

    print("\n" + "=" * 78)
    print("Benchmark Summary")
    print("=" * 78)
    print(f"Total wall time      : {total_wall_time:.6f} s")
    print(f"Average step time    : {total_wall_time / max(config.n_steps, 1):.6f} s")
    print(f"Inlet mass flow      : {float(inlet_junction.W):.6f} kg/s")
    print(f"Outlet mass flow     : {float(outlet_junction.W):.6f} kg/s")
    print(
        "Branch split         : "
        f"[{branch_split[0]:.2f}, {branch_split[1]:.2f}, {branch_split[2]:.2f}, {branch_split[3]:.2f}] %"
    )
    print(f"Inlet buffer outlet T: {float(inlet_buffer.volumes[-1].T):.3f} K")
    print(f"Hot leg outlet T     : {float(hot_leg.volumes[-1].T):.3f} K")
    print(f"Radiator 1 outlet T  : {float(radiator_1.volumes[-1].T):.3f} K")
    print(f"Radiator 2 outlet T  : {float(radiator_2.volumes[-1].T):.3f} K")
    print(f"Return leg outlet T  : {float(return_leg.volumes[-1].T):.3f} K")
    print(f"Outlet buffer outlet T: {float(outlet_buffer.volumes[-1].T):.3f} K")
    print(f"Base Q_vol           : {source_coupler.base_heating_power():.3f} W")
    print(f"Final-step Q_vol     : {source_coupler.current_heating_power():.3f} W")
    print(f"Total hA             : {source_coupler.total_cooling_conductance():.3f} W/K")
    print(f"Effective cooling    : {source_coupler.effective_cooling_power():.3f} W")

    for idx, branch in enumerate(branches, start=1):
        print(
            f"Branch {idx:1d} outlet T     : "
            f"{float(branch['outlet_pipe'].volumes[-1].T):.3f} K"
        )

    print("=" * 78)


def run_benchmark(config: BenchmarkConfig = None) -> Dict[str, object]:
    if config is None:
        config = BenchmarkConfig()

    model = build_benchmark_model(config)
    system: SystemManager = model["system"]

    print_model_summary(model)

    init_start = time.perf_counter()
    system.initialize_system(
        dt_init=config.init_dt,
        tol=config.init_tol,
        max_iter=config.init_max_iter,
    )
    init_wall_time = time.perf_counter() - init_start
    print(f"Initialization time : {init_wall_time:.6f} s")

    TEASAProfiler.stats.clear()

    total_start = time.perf_counter()
    for step_index in range(1, config.n_steps + 1):
        step_start = time.perf_counter()
        system.step(dt=config.dt, inner_iter=config.inner_iter)
        step_wall_time = time.perf_counter() - step_start

        if step_index == 1 or step_index % max(config.print_every, 1) == 0 or step_index == config.n_steps:
            print_progress(step_index, model, step_wall_time)

    total_wall_time = time.perf_counter() - total_start
    print_final_summary(model, total_wall_time)
    TEASAProfiler.report()
    return model


if __name__ == "__main__":
    run_benchmark()
