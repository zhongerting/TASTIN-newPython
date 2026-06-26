import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Materials.Solids.WallMaterial import SS316
from Components.RadiatorPipeWithFin import RadiatorPipeWithFin
from Components.RadiatorThermalShield import RadiatorThermalShield
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidChannel, IncompressibleFluidVolume
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager


SIGMA = 5.670374419e-8


def nak_internal_nu(Re, Pr, _p_d_ratio=1.1):
    pe = np.maximum(np.asarray(Re, dtype=float) * np.asarray(Pr, dtype=float), 1.0)
    return 7.0 + 0.025 * pe ** 0.8


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _parse_float_list(text: str) -> List[float]:
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def _iter_channel_junctions(channels: Iterable[IncompressibleFluidChannel]):
    for channel in channels:
        yield from channel.internal_junctions


def build_model(args: argparse.Namespace) -> Dict[str, Any]:
    if getattr(args, "hydraulic_calibrated", False):
        args.tube_inlet_k_loss = 100.0
        args.tube_outlet_k_loss = 100.0

    n_tubes = int(args.n_tubes)
    n_axial = int(args.n_axial)
    inlet_temperature = float(args.inlet_temperature_k)
    initial_temperature = float(args.initial_temperature_k)
    outlet_pressure_value = float(args.outlet_pressure_pa)

    nak = SodiumPotassium78()
    wall_mat = SS316(name="TOPAZ2_Radiator_SS316")

    tube_length = float(args.tube_length_m)
    tube_inner_d = float(args.tube_inner_diameter_m)
    tube_outer_d = float(args.tube_outer_diameter_m)
    tube_flow_area = math.pi * (tube_inner_d ** 2) / 4.0

    upper_d = float(args.upper_header_centerline_diameter_m)
    lower_d = float(args.lower_header_centerline_diameter_m)
    header_inner_d = float(args.header_inner_diameter_m)
    header_area = math.pi * (header_inner_d ** 2) / 4.0
    upper_seg_len = math.pi * upper_d / n_tubes
    lower_seg_len = math.pi * lower_d / n_tubes
    header_dh = header_inner_d

    inlet_area = math.pi * (float(args.hot_leg_inner_diameter_m) ** 2) / 4.0
    outlet_area = inlet_area

    inlet_a = IncompressibleBoundaryVolume(
        name="TOPAZ2_Inlet_A",
        material=nak,
        P=outlet_pressure_value + float(args.initial_pressure_rise_pa),
        T=inlet_temperature,
        flow_area=inlet_area,
        hydraulic_diam=float(args.hot_leg_inner_diameter_m),
    )
    inlet_b = IncompressibleBoundaryVolume(
        name="TOPAZ2_Inlet_B",
        material=nak,
        P=outlet_pressure_value + float(args.initial_pressure_rise_pa),
        T=inlet_temperature,
        flow_area=inlet_area,
        hydraulic_diam=float(args.hot_leg_inner_diameter_m),
    )
    outlet_pressure = IncompressibleBoundaryVolume(
        name="TOPAZ2_OutletPressure",
        material=nak,
        P=outlet_pressure_value,
        T=initial_temperature,
        flow_area=outlet_area,
        hydraulic_diam=float(args.cold_leg_inner_diameter_m),
    )
    outlet_pressure.is_pressure_boundary = True
    outlet_mix = IncompressibleFluidVolume(
        name="TOPAZ2_OutletMix",
        volume=float(args.outlet_mix_volume_m3),
        length=float(args.outlet_mix_length_m),
        flow_area=2.0 * outlet_area,
        hydraulic_diam=float(args.cold_leg_inner_diameter_m),
        initial_P=outlet_pressure_value,
        initial_T=initial_temperature,
        material=nak,
    )

    upper_nodes = []
    lower_nodes = []
    for i in range(n_tubes):
        upper = IncompressibleFluidChannel(
            name=f"TOPAZ2_UpperHeader_{i + 1:02d}",
            n_nodes=1,
            total_length=upper_seg_len,
            flow_area=header_area,
            hydraulic_diam=header_dh,
            initial_P=outlet_pressure_value + float(args.initial_pressure_rise_pa),
            initial_T=inlet_temperature,
            material=nak,
        )
        lower = IncompressibleFluidChannel(
            name=f"TOPAZ2_LowerHeader_{i + 1:02d}",
            n_nodes=1,
            total_length=lower_seg_len,
            flow_area=header_area,
            hydraulic_diam=header_dh,
            initial_P=outlet_pressure_value,
            initial_T=initial_temperature,
            material=nak,
        )
        upper_nodes.append(upper)
        lower_nodes.append(lower)

    tube_channels: List[IncompressibleFluidChannel] = []
    radiator_units: List[RadiatorPipeWithFin] = []

    for i in range(n_tubes):
        channel = IncompressibleFluidChannel(
            name=f"TOPAZ2_RadTubeFluid_{i + 1:02d}",
            n_nodes=n_axial,
            total_length=tube_length,
            flow_area=tube_flow_area,
            hydraulic_diam=tube_inner_d,
            initial_P=outlet_pressure_value + 0.5 * float(args.initial_pressure_rise_pa),
            initial_T=initial_temperature,
            material=nak,
        )
        radiator = RadiatorPipeWithFin(
            name=f"TOPAZ2_RadTube_{i + 1:02d}",
            fluid_channel=channel,
            wall_material=wall_mat,
            tube_inner_diameter=tube_inner_d,
            tube_outer_diameter=tube_outer_d,
            tube_length=tube_length,
            n_axial=n_axial,
            n_radial_wall=max(1, int(args.n_radial_wall)),
            fin_thickness=float(args.fin_thickness_m),
            fin_width_upper=float(args.fin_width_upper_m),
            fin_width_lower=float(args.fin_width_lower_m),
            n_fin_width=int(args.n_fin_width),
            correlation_func=nak_internal_nu,
            tube_emissivity=float(args.tube_emissivity),
            fin_emissivity=float(args.fin_emissivity),
            tube_area_scale=float(args.tube_area_scale),
            fin_area_scale=float(args.fin_area_scale),
            T_space=float(args.t_space_k),
            initial_temp=initial_temperature,
            fin_conductivity=float(args.fin_conductivity_w_m_k),
            fin_view_factor=float(args.fin_view_factor),
            contact_resistance_m2k_w=float(args.fin_contact_resistance_m2k_w),
            coupling_time_scheme=args.fluid_solid_coupling_scheme,
            solid_ode_method=args.solid_ode_method,
        )
        tube_channels.append(channel)
        radiator_units.append(radiator)

    volumes = [inlet_a, inlet_b, outlet_pressure, outlet_mix]
    junctions = []

    for channels in (upper_nodes, lower_nodes, tube_channels):
        for channel in channels:
            volumes.extend(channel.volumes)
            junctions.extend(channel.internal_junctions)

    total_flow = float(args.total_mass_flow_kg_s)
    inlet_idx_a = 0
    inlet_idx_b = n_tubes // 2
    outlet_idx_a = 0
    outlet_idx_b = n_tubes // 2

    junctions.append(
        InletJunction(
            name="TOPAZ2_Inlet_A_to_Upper",
            from_vol=inlet_a,
            to_vol=upper_nodes[inlet_idx_a].volumes[0],
            W_initial=0.5 * total_flow,
        )
    )
    junctions.append(
        InletJunction(
            name="TOPAZ2_Inlet_B_to_Upper",
            from_vol=inlet_b,
            to_vol=upper_nodes[inlet_idx_b].volumes[0],
            W_initial=0.5 * total_flow,
        )
    )
    junctions.append(
        FlowJunction(
            name="TOPAZ2_Lower_to_Outlet_A",
            from_vol=lower_nodes[outlet_idx_a].volumes[0],
            to_vol=outlet_mix,
            flow_area=outlet_area,
            k_loss=float(args.outlet_k_loss),
            custom_length=float(args.cold_leg_length_m),
            hydraulic_diam=float(args.cold_leg_inner_diameter_m),
        )
    )
    junctions.append(
        FlowJunction(
            name="TOPAZ2_Lower_to_Outlet_B",
            from_vol=lower_nodes[outlet_idx_b].volumes[0],
            to_vol=outlet_mix,
            flow_area=outlet_area,
            k_loss=float(args.outlet_k_loss),
            custom_length=float(args.cold_leg_length_m),
            hydraulic_diam=float(args.cold_leg_inner_diameter_m),
        )
    )
    junctions.append(
        FlowJunction(
            name="TOPAZ2_OutletMix_to_PressureBoundary",
            from_vol=outlet_mix,
            to_vol=outlet_pressure,
            flow_area=2.0 * outlet_area,
            k_loss=float(args.outlet_k_loss),
            custom_length=float(args.outlet_mix_length_m),
            hydraulic_diam=float(args.cold_leg_inner_diameter_m),
        )
    )

    for i in range(n_tubes):
        j = (i + 1) % n_tubes
        junctions.append(
            FlowJunction(
                name=f"TOPAZ2_UpperRing_{i + 1:02d}_to_{j + 1:02d}",
                from_vol=upper_nodes[i].volumes[0],
                to_vol=upper_nodes[j].volumes[0],
                flow_area=header_area,
                k_loss=float(args.header_k_loss),
                custom_length=upper_seg_len,
                hydraulic_diam=header_dh,
            )
        )
        junctions.append(
            FlowJunction(
                name=f"TOPAZ2_LowerRing_{i + 1:02d}_to_{j + 1:02d}",
                from_vol=lower_nodes[i].volumes[0],
                to_vol=lower_nodes[j].volumes[0],
                flow_area=header_area,
                k_loss=float(args.header_k_loss),
                custom_length=lower_seg_len,
                hydraulic_diam=header_dh,
            )
        )
        junctions.append(
            FlowJunction(
                name=f"TOPAZ2_Upper_to_Tube_{i + 1:02d}",
                from_vol=upper_nodes[i].volumes[0],
                to_vol=tube_channels[i].volumes[0],
                flow_area=tube_flow_area,
                k_loss=float(args.tube_inlet_k_loss),
                custom_length=0.5 * tube_length / n_axial,
                hydraulic_diam=tube_inner_d,
            )
        )
        junctions.append(
            FlowJunction(
                name=f"TOPAZ2_Tube_{i + 1:02d}_to_Lower",
                from_vol=tube_channels[i].volumes[-1],
                to_vol=lower_nodes[i].volumes[0],
                flow_area=tube_flow_area,
                k_loss=float(args.tube_outlet_k_loss),
                custom_length=0.5 * tube_length / n_axial,
                hydraulic_diam=tube_inner_d,
            )
        )

    network = HydraulicNetwork(volumes=volumes, junctions=junctions, gravity_vector=0.0)
    system = SystemManager(fluid_network=network, start_time=0.0)
    for radiator in radiator_units:
        system.add_component(radiator)

    return {
        "args": args,
        "system": system,
        "network": network,
        "nak": nak,
        "inlets": [inlet_a, inlet_b],
        "outlets": [outlet_pressure],
        "outlet_mix": outlet_mix,
        "upper_nodes": upper_nodes,
        "lower_nodes": lower_nodes,
        "tube_channels": tube_channels,
        "radiator_units": radiator_units,
        "tube_solids": [unit.wall for unit in radiator_units],
        "couplers": [unit.coupler for unit in radiator_units],
        "external_area_per_node": radiator_units[0].tube_bare_area + radiator_units[0].fin_radiating_area,
        "tube_bare_area_per_node": radiator_units[0].tube_bare_area,
        "fin_area_per_node": radiator_units[0].fin_radiating_area,
        "tube_flow_area": tube_flow_area,
        "header_area": header_area,
        "upper_seg_len": upper_seg_len,
        "lower_seg_len": lower_seg_len,
        "tube_length": tube_length,
    }


def attach_radiator_thermal_shield(build: Dict[str, Any], args: argparse.Namespace) -> RadiatorThermalShield:
    active_until = None
    if getattr(args, "shield_active_until_s", None) is not None:
        active_until = float(build["system"].global_time) + float(args.shield_active_until_s)
    shield = RadiatorThermalShield(
        name="TOPAZ2_RadiatorThermalShield",
        radiator_units=build["radiator_units"],
        active_until_s=active_until,
        background_temperature_k=float(args.shield_background_temperature_k),
        shield_view_factor=float(args.shield_view_factor),
        inner_emissivity=float(args.shield_inner_emissivity),
        outer_emissivity=float(args.shield_outer_emissivity),
        conductivity_w_m_k=float(args.shield_conductivity_w_m_k),
        thickness_m=float(args.shield_thickness_m),
        solar_heat_flux_w_m2=float(args.shield_solar_heat_flux_w_m2),
        relaxation=float(args.shield_relaxation),
        model=str(args.shield_model),
    )
    system = build["system"]
    if shield not in system.components:
        system.components.insert(0, shield)
    build["radiator_thermal_shield"] = shield
    return shield


def radiation_shield_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    shield = build.get("radiator_thermal_shield")
    if shield is None:
        return {
            "radiation_shield_enabled": False,
            "radiation_shield_active": False,
            "radiation_shield_model": None,
            "radiation_shield_effective_background_mean_k": None,
            "radiation_shield_inner_temperature_mean_k": None,
            "radiation_shield_outer_temperature_mean_k": None,
            "radiation_shield_q_from_radiator_w": 0.0,
            "radiation_shield_q_solar_w": 0.0,
            "radiation_shield_q_to_space_w": 0.0,
            "radiation_shield_solver_failures": 0,
        }
    diagnostics = shield.get_diagnostics()
    diagnostics["radiation_shield_enabled"] = True
    return diagnostics


def initialize_flow_guess(build: Dict[str, Any]) -> None:
    args = build["args"]
    total_flow = float(args.total_mass_flow_kg_s)
    n_tubes = int(args.n_tubes)
    tube_guess = total_flow / n_tubes
    for channel in build["tube_channels"]:
        for junc in channel.internal_junctions:
            junc.W = tube_guess
            junc.update_velocity()

    # Header guesses follow the ideal symmetric distribution between inlet and symmetry plane.
    quarter = n_tubes // 4
    for junc in build["network"].junctions_obj:
        name = getattr(junc, "name", "")
        if name.startswith("TOPAZ2_UpperRing_") or name.startswith("TOPAZ2_LowerRing_"):
            junc.W = 0.0
            junc.update_velocity()
        elif name.startswith("TOPAZ2_Upper_to_Tube_") or name.startswith("TOPAZ2_Tube_"):
            junc.W = tube_guess
            junc.update_velocity()
        elif name.startswith("TOPAZ2_Lower_to_Outlet_"):
            junc.W = 0.5 * total_flow
            junc.update_velocity()

    for inlet in [j for j in build["network"].junctions_obj if isinstance(j, InletJunction)]:
        inlet.W = inlet.target_W
        inlet.update_velocity()


def tube_mass_flows(build: Dict[str, Any]) -> np.ndarray:
    flows = []
    for channel in build["tube_channels"]:
        vals = [j.W for j in channel.internal_junctions]
        if vals:
            flows.append(float(np.mean(vals)))
        else:
            flows.append(0.0)
    return np.array(flows, dtype=float)


def tube_flow_distribution(build: Dict[str, Any]) -> Dict[str, Any]:
    flows = tube_mass_flows(build)
    total = float(np.sum(flows))
    mean = float(np.mean(flows)) if flows.size else 0.0
    min_idx = int(np.argmin(flows)) if flows.size else -1
    max_idx = int(np.argmax(flows)) if flows.size else -1
    n_half = flows.size // 2
    if n_half > 0 and flows.size % 2 == 0:
        symmetry_error = float(np.max(np.abs(flows[:n_half] - flows[n_half:])))
    else:
        symmetry_error = float("nan")

    rows = []
    for i, flow in enumerate(flows, start=1):
        rows.append(
            {
                "tube": i,
                "mass_flow_kg_s": float(flow),
                "percent_of_tube_sum": float(flow / total * 100.0) if total else 0.0,
                "relative_to_mean": float(flow / mean) if mean else 0.0,
            }
        )

    return {
        "total_tube_flow_kg_s": total,
        "mean_tube_flow_kg_s": mean,
        "min_tube_flow_kg_s": float(flows[min_idx]) if flows.size else 0.0,
        "min_tube_index": min_idx + 1 if flows.size else None,
        "max_tube_flow_kg_s": float(flows[max_idx]) if flows.size else 0.0,
        "max_tube_index": max_idx + 1 if flows.size else None,
        "flow_spread_kg_s": float(np.max(flows) - np.min(flows)) if flows.size else 0.0,
        "flow_spread_over_mean": float((np.max(flows) - np.min(flows)) / mean) if mean else 0.0,
        "max_over_min": float(np.max(flows) / np.min(flows)) if flows.size and np.min(flows) else 0.0,
        "symmetry_error_kg_s": symmetry_error,
        "tube_flows": rows,
    }


def _named_junctions(build: Dict[str, Any], prefix: str) -> List[Any]:
    return [
        junc
        for junc in build["network"].junctions_obj
        if getattr(junc, "name", "").startswith(prefix)
    ]


def header_flow_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    upper = _named_junctions(build, "TOPAZ2_UpperRing_")
    lower = _named_junctions(build, "TOPAZ2_LowerRing_")

    def summarize(junctions):
        flows = np.array([float(j.W) for j in junctions], dtype=float)
        abs_flows = np.abs(flows)
        return {
            "n_segments": int(flows.size),
            "min_signed_flow_kg_s": float(np.min(flows)) if flows.size else 0.0,
            "max_signed_flow_kg_s": float(np.max(flows)) if flows.size else 0.0,
            "max_abs_flow_kg_s": float(np.max(abs_flows)) if flows.size else 0.0,
            "direction_reversal_count": int(np.sum(flows < 0.0)) if flows.size else 0,
            "segment_flows": [
                {
                    "segment": i,
                    "name": getattr(junc, "name", ""),
                    "mass_flow_kg_s": float(junc.W),
                }
                for i, junc in enumerate(junctions, start=1)
            ],
        }

    return {
        "upper_header": summarize(upper),
        "lower_header": summarize(lower),
    }


def pressure_budget_diagnostics(
    build: Dict[str, Any],
    tube_indices: Tuple[int, ...] = (1, 20, 40, 59),
) -> List[Dict[str, float]]:
    rows = []
    upper_nodes = build["upper_nodes"]
    lower_nodes = build["lower_nodes"]
    tube_channels = build["tube_channels"]
    n_tubes = len(tube_channels)
    for tube_index in tube_indices:
        if tube_index < 1 or tube_index > n_tubes:
            continue
        idx = tube_index - 1
        upper_p = float(upper_nodes[idx].volumes[0].P)
        lower_p = float(lower_nodes[idx].volumes[0].P)
        flow = float(np.mean([j.W for j in tube_channels[idx].internal_junctions]))
        rows.append(
            {
                "tube": int(tube_index),
                "upper_header_pressure_pa": upper_p,
                "lower_header_pressure_pa": lower_p,
                "drive_delta_p_pa": upper_p - lower_p,
                "mass_flow_kg_s": flow,
            }
        )
    return rows


def radiation_breakdown(build: Dict[str, Any]) -> Tuple[float, float, float]:
    tube_total = 0.0
    fin_total = 0.0
    for unit in build["radiator_units"]:
        breakdown = unit.get_heat_exchange_breakdown()
        tube_total += float(np.sum(breakdown["bare_radiation"]))
        fin_total += float(np.sum(breakdown["fin_radiation"]))
    return tube_total, fin_total, tube_total + fin_total


def mixed_outlet_state(build: Dict[str, Any]) -> Dict[str, float]:
    outlet_juncs = [
        j for j in build["network"].junctions_obj
        if getattr(j, "name", "").startswith("TOPAZ2_Lower_to_Outlet_")
    ]
    m_total = 0.0
    h_total = 0.0
    for junc in outlet_juncs:
        w = max(float(junc.W), 0.0)
        donor = junc.from_vol if junc.W >= 0.0 else junc.to_vol
        m_total += w
        h_total += w * float(donor.h)
    if m_total <= 0.0:
        return {"mass_flow_kg_s": 0.0, "enthalpy_j_kg": 0.0, "temperature_k": float("nan")}
    h_mix = h_total / m_total
    p_out = float(build["args"].outlet_pressure_pa)
    t_mix = float(build["nak"].temperature_from_enthalpy(h_mix, p_out))
    return {"mass_flow_kg_s": m_total, "enthalpy_j_kg": h_mix, "temperature_k": t_mix}


def collect_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    args = build["args"]
    flows = tube_mass_flows(build)
    mix = mixed_outlet_state(build)
    q_tube, q_fin, q_total = radiation_breakdown(build)
    fin_temps = []
    fin_root_temps = []
    fin_tip_temps = []
    fin_iteration_counts = []
    fin_max_deltas = []
    fin_warm_start_flags = []
    for unit in build["radiator_units"]:
        fin = unit.get_fin_temperature_distribution()
        if fin.size:
            fin_temps.append(fin.reshape(-1))
            fin_root_temps.append(fin[:, 0])
            fin_tip_temps.append(fin[:, -1])
        fin_iteration_counts.append(float(getattr(unit, "last_fin_iteration_count", 0)))
        fin_max_deltas.append(float(getattr(unit, "last_fin_max_delta", 0.0)))
        fin_warm_start_flags.append(bool(getattr(unit, "last_fin_used_warm_start", False)))
    if fin_temps:
        fin_all = np.concatenate(fin_temps)
        fin_root_all = np.concatenate(fin_root_temps)
        fin_tip_all = np.concatenate(fin_tip_temps)
    else:
        fin_all = np.array([], dtype=float)
        fin_root_all = np.array([], dtype=float)
        fin_tip_all = np.array([], dtype=float)
    inlet_h = float(build["nak"].enthalpy(float(args.inlet_temperature_k), float(args.outlet_pressure_pa)))
    q_fluid_drop = mix["mass_flow_kg_s"] * (inlet_h - mix["enthalpy_j_kg"])
    diagnostics = {
        "time_s": float(build["system"].global_time),
        "tube_emissivity": float(args.tube_emissivity),
        "fin_emissivity": float(args.fin_emissivity),
        "inlet_temperature_k": float(args.inlet_temperature_k),
        "mixed_outlet_temperature_k": mix["temperature_k"],
        "target_outlet_temperature_k": float(args.target_outlet_temperature_k),
        "outlet_delta_from_target_k": mix["temperature_k"] - float(args.target_outlet_temperature_k),
        "outlet_mass_flow_kg_s": mix["mass_flow_kg_s"],
        "tube_flow_min_kg_s": float(np.min(flows)) if flows.size else 0.0,
        "tube_flow_max_kg_s": float(np.max(flows)) if flows.size else 0.0,
        "tube_flow_mean_kg_s": float(np.mean(flows)) if flows.size else 0.0,
        "tube_flow_spread_kg_s": float(np.max(flows) - np.min(flows)) if flows.size else 0.0,
        "tube_flow_rel_spread": float((np.max(flows) - np.min(flows)) / np.mean(flows)) if flows.size and np.mean(flows) else 0.0,
        "tube_flow_min_index": int(np.argmin(flows) + 1) if flows.size else -1,
        "tube_flow_max_index": int(np.argmax(flows) + 1) if flows.size else -1,
        "q_tube_radiation_w": q_tube,
        "q_fin_radiation_w": q_fin,
        "q_total_radiation_w": q_total,
        "q_fluid_temperature_drop_w": q_fluid_drop,
        "energy_residual_w": q_fluid_drop - q_total,
        "total_effective_area_m2": float(int(args.n_tubes) * np.sum(build["external_area_per_node"])),
        "fin_temperature_mean_k": float(np.mean(fin_all)) if fin_all.size else float("nan"),
        "fin_temperature_min_k": float(np.min(fin_all)) if fin_all.size else float("nan"),
        "fin_root_temperature_mean_k": float(np.mean(fin_root_all)) if fin_root_all.size else float("nan"),
        "fin_tip_temperature_mean_k": float(np.mean(fin_tip_all)) if fin_tip_all.size else float("nan"),
        "fin_root_to_tip_delta_mean_k": float(np.mean(fin_root_all - fin_tip_all)) if fin_tip_all.size else float("nan"),
        "fin_iteration_mean": float(np.mean(fin_iteration_counts)) if fin_iteration_counts else 0.0,
        "fin_iteration_max": float(np.max(fin_iteration_counts)) if fin_iteration_counts else 0.0,
        "fin_max_delta_mean_k": float(np.mean(fin_max_deltas)) if fin_max_deltas else 0.0,
        "fin_max_delta_max_k": float(np.max(fin_max_deltas)) if fin_max_deltas else 0.0,
        "fin_warm_start_fraction": float(np.mean(fin_warm_start_flags)) if fin_warm_start_flags else 0.0,
        "n_fin_width": int(args.n_fin_width),
        "upper_header_arc_length_sum_m": float(int(args.n_tubes) * build["upper_seg_len"]),
        "lower_header_arc_length_sum_m": float(int(args.n_tubes) * build["lower_seg_len"]),
    }
    diagnostics.update(radiation_shield_diagnostics(build))
    return diagnostics


def write_history_row(path: Path, row: Dict[str, Any], write_header: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()), extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def run_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_model(args)
    initialize_flow_guess(build)

    system = build["system"]
    if bool(getattr(args, "enable_radiation_shield", False)):
        shield = attach_radiator_thermal_shield(build, args)
        shield.pre_step(0.0, float(system.global_time))
    system.initialize_system(dt_init=float(args.init_dt), tol=float(args.hydraulic_tol), max_iter=int(args.hydraulic_max_iter))
    if bool(getattr(args, "enable_radiation_shield", False)):
        build["radiator_thermal_shield"].pre_step(0.0, float(system.global_time))

    output_dir = Path(args.output_dir)
    history_path = output_dir / f"{args.case_prefix}_history.csv"
    latest_state_path = output_dir / f"{args.case_prefix}_latest_state.json"
    if history_path.exists() and not args.append_history:
        history_path.unlink()

    next_record = 0.0
    t_end = float(args.duration)
    max_dt = float(args.max_dt)
    write_header = not history_path.exists()
    last_record = collect_diagnostics(build)
    while system.global_time < t_end - 1.0e-12:
        dt = min(max_dt, t_end - system.global_time)
        system.step(dt=dt, inner_iter=int(args.inner_iter), convergence_tol=float(args.convergence_tol))
        if system.global_time >= next_record - 1.0e-12:
            last_record = collect_diagnostics(build)
            write_history_row(history_path, last_record, write_header)
            write_header = False
            next_record += float(args.record_interval)

    last_record = collect_diagnostics(build)
    if not history_path.exists() or (last_record.get("time_s") != system.global_time):
        write_history_row(history_path, last_record, write_header)
    latest = {
        "case_prefix": args.case_prefix,
        "history_path": history_path,
        "latest_record": last_record,
        "geometry": {
            "n_tubes": int(args.n_tubes),
            "n_axial": int(args.n_axial),
            "tube_length_m": float(args.tube_length_m),
            "tube_inner_diameter_m": float(args.tube_inner_diameter_m),
            "tube_outer_diameter_m": float(args.tube_outer_diameter_m),
            "fin_thickness_m": float(args.fin_thickness_m),
            "n_fin_width": int(args.n_fin_width),
            "upper_header_centerline_diameter_m": float(args.upper_header_centerline_diameter_m),
            "lower_header_centerline_diameter_m": float(args.lower_header_centerline_diameter_m),
            "header_inner_diameter_m": float(args.header_inner_diameter_m),
            "effective_area_m2": last_record["total_effective_area_m2"],
            "hydraulic_calibrated": bool(getattr(args, "hydraulic_calibrated", False)),
            "tube_inlet_k_loss": float(args.tube_inlet_k_loss),
            "tube_outlet_k_loss": float(args.tube_outlet_k_loss),
        },
        "radiation_shield_enabled": bool(getattr(args, "enable_radiation_shield", False)),
        "radiation_shield_model": getattr(args, "shield_model", None),
        "radiation_shield_active_until_s": getattr(args, "shield_active_until_s", None),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with latest_state_path.open("w", encoding="utf-8") as f:
        json.dump(latest, f, indent=2, ensure_ascii=False, default=_json_default)
    return latest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Standalone TOPAZ-II NaK tube-fin radiator case.")
    parser.add_argument("--output-dir", default="CoolantLoop/topaz2_tube_fin_radiator")
    parser.add_argument("--case-prefix", default="topaz2_tube_fin_radiator")
    parser.add_argument("--duration", type=float, default=50.0)
    parser.add_argument("--record-interval", type=float, default=5.0)
    parser.add_argument("--max-dt", type=float, default=0.2)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--convergence-tol", type=float, default=1.0e-3)
    parser.add_argument("--append-history", action="store_true")
    parser.add_argument("--init-dt", type=float, default=0.05)
    parser.add_argument("--hydraulic-tol", type=float, default=1.0e-6)
    parser.add_argument("--hydraulic-max-iter", type=int, default=800)

    parser.add_argument("--n-tubes", type=int, default=78)
    parser.add_argument("--n-axial", type=int, default=8)
    parser.add_argument("--n-radial-wall", type=int, default=1)
    parser.add_argument("--tube-length-m", type=float, default=1.85)
    parser.add_argument("--tube-outer-diameter-m", type=float, default=0.008)
    parser.add_argument("--tube-inner-diameter-m", type=float, default=0.007)
    parser.add_argument("--upper-header-centerline-diameter-m", type=float, default=0.824)
    parser.add_argument("--lower-header-centerline-diameter-m", type=float, default=1.346)
    parser.add_argument("--header-inner-diameter-m", type=float, default=0.020)
    parser.add_argument("--hot-leg-inner-diameter-m", type=float, default=0.030)
    parser.add_argument("--cold-leg-inner-diameter-m", type=float, default=0.030)
    parser.add_argument("--cold-leg-length-m", type=float, default=3.5)
    parser.add_argument("--fin-thickness-m", type=float, default=0.0004)
    parser.add_argument("--fin-width-upper-m", type=float, default=0.03319)
    parser.add_argument("--fin-width-lower-m", type=float, default=0.05421)
    parser.add_argument("--n-fin-width", type=int, default=12)
    parser.add_argument("--fin-conductivity-w-m-k", type=float, default=348.9)
    parser.add_argument("--fin-view-factor", type=float, default=1.0)
    parser.add_argument("--fin-contact-resistance-m2k-w", type=float, default=0.0)
    parser.add_argument("--tube-area-scale", type=float, default=1.0)
    parser.add_argument("--fin-area-scale", type=float, default=0.35)

    parser.add_argument("--total-mass-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--inlet-temperature-k", type=float, default=823.0)
    parser.add_argument("--initial-temperature-k", type=float, default=727.0)
    parser.add_argument("--target-outlet-temperature-k", type=float, default=727.0)
    parser.add_argument("--outlet-pressure-pa", type=float, default=160000.0)
    parser.add_argument("--initial-pressure-rise-pa", type=float, default=6000.0)
    parser.add_argument("--outlet-mix-volume-m3", type=float, default=1.0e-4)
    parser.add_argument("--outlet-mix-length-m", type=float, default=0.05)

    parser.add_argument("--tube-emissivity", type=float, default=0.80)
    parser.add_argument("--fin-emissivity", type=float, default=0.80)
    parser.add_argument("--scan-emissivity-grid", action="store_true")
    parser.add_argument("--tube-emissivity-values", default="0.75,0.80,0.85,0.88,0.90,0.92")
    parser.add_argument("--fin-emissivity-values", default="0.75,0.80,0.85,0.88,0.90,0.92")
    parser.add_argument("--t-space-k", type=float, default=3.0)
    parser.add_argument("--enable-radiation-shield", "--enable-radiator-shield", dest="enable_radiation_shield", action="store_true")
    parser.add_argument("--shield-active-until-s", type=float, default=None)
    parser.add_argument("--shield-inner-emissivity", type=float, default=0.8)
    parser.add_argument("--shield-outer-emissivity", type=float, default=0.1)
    parser.add_argument("--shield-conductivity-w-m-k", type=float, default=0.0008)
    parser.add_argument("--shield-thickness-m", type=float, default=0.01)
    parser.add_argument("--shield-view-factor", type=float, default=0.8)
    parser.add_argument("--shield-solar-heat-flux-w-m2", type=float, default=0.0)
    parser.add_argument("--shield-background-temperature-k", type=float, default=3.0)
    parser.add_argument("--shield-relaxation", type=float, default=1.0)
    parser.add_argument("--shield-model", choices=("segment_balance", "fortran_shield2"), default="fortran_shield2")
    parser.add_argument("--header-k-loss", type=float, default=1.0)
    parser.add_argument("--tube-inlet-k-loss", type=float, default=2.0)
    parser.add_argument("--tube-outlet-k-loss", type=float, default=2.0)
    parser.add_argument("--outlet-k-loss", type=float, default=2.0)
    parser.add_argument("--hydraulic-calibrated", action="store_true")
    parser.add_argument("--solid-ode-method", default="BDF", choices=("RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"))
    parser.add_argument("--fluid-solid-coupling-scheme", default="current", choices=("current", "local_implicit"))
    return parser


def make_default_args(**overrides) -> argparse.Namespace:
    args = build_arg_parser().parse_args([])
    for key, value in overrides.items():
        if not hasattr(args, key):
            raise AttributeError(f"Unknown TOPAZ-II radiator argument: {key}")
        setattr(args, key, value)
    return args


def parse_args() -> argparse.Namespace:
    parser = build_arg_parser()
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.scan_emissivity_grid:
        results = []
        base_prefix = args.case_prefix
        tube_values = _parse_float_list(args.tube_emissivity_values)
        fin_values = _parse_float_list(args.fin_emissivity_values)
        for tube_eps in tube_values:
            for fin_eps in fin_values:
                case_args = argparse.Namespace(**vars(args))
                case_args.scan_emissivity_grid = False
                case_args.tube_emissivity = tube_eps
                case_args.fin_emissivity = fin_eps
                case_args.case_prefix = f"{base_prefix}_tube{tube_eps:.3f}_fin{fin_eps:.3f}".replace(".", "p")
                latest = run_case(case_args)
                record = latest["latest_record"]
                results.append(record)
                print(json.dumps(record, ensure_ascii=False, default=_json_default))

        best = min(
            results,
            key=lambda row: abs(float(row["outlet_delta_from_target_k"])),
        )
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        summary_path = output_dir / f"{base_prefix}_emissivity_scan_summary.csv"
        with summary_path.open("w", newline="", encoding="utf-8-sig") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        print(json.dumps({"best": best, "summary_path": str(summary_path)}, indent=2, ensure_ascii=False, default=_json_default))
    else:
        latest = run_case(args)
        print(json.dumps(latest["latest_record"], indent=2, ensure_ascii=False, default=_json_default))


if __name__ == "__main__":
    main()
