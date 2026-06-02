import argparse
import csv
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import numpy as np

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent
for path in (str(CURRENT_DIR), str(ROOT_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from Components.TFEUnit import GapConfig, TFEGeometry, TFEMeshParams, TFEUnit
from Components.tec_electric import electric_field_from_node_potential, joule_power_from_electric_field
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.BerylliumOxide import BerylliumOxide
from Materials.Solids.GasGaps import CarbonDioxide, Cesium, Helium, Xenon
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.UO2 import UO2
from Materials.Solids.ZrH import ZirconiumHydride
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import FlowJunction, NonUniformIncompressibleFluidChannel
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Couplers import FluidSolidCouple, TECCouple2D
from Solvers.SystemManager import SystemManager


SOURCE_SNAPSHOT = "test_core_assemble_v7_caseA_faststeady_restart_t18800.npz"
SOURCE_SNAPSHOT_TIME_S = 18800.0
SOURCE_FLUID_SHAPE = (176, 179)
N_AXIAL = 37
FIXED_NUCLEAR_HEAT_W = 3489.8792830760617
FIXED_INLET_VELOCITY_M_S = 0.7533237780835761
FIXED_INLET_MASS_FLOW_KG_S = 0.03513560646568767
FIXED_OUTLET_PRESSURE_PA = 161961.33075474895
TEC_TARGET_VOLTAGE_V = 0.8
TEC_WIRE_RESISTANCE_OHM = np.zeros(4, dtype=float)
ENERGY_HISTORY_FIELDS = (
    "time_s",
    "dt_s",
    "nuclear_heat_w",
    "coolant_enthalpy_pickup_w",
    "electrical_output_w",
    "solid_storage_rate_w",
    "fluid_storage_rate_w",
    "outer_clad_heat_loss_w",
    "residual_w",
    "relative_residual",
    "inlet_mass_flow_kg_s",
    "outlet_mass_flow_kg_s",
)
AXIAL_SHAPE_COEFFS = np.array(
    [6.27905178e-02, -7.13913811e-02, 1.35276842e-02, -1.02326367e-03, 3.90936491e-05],
    dtype=float,
)

SOURCE_INLET_PLENUM_VOLUME_INDEX = 1
SOURCE_OUTLET_PLENUM_VOLUME_INDEX = 2
SOURCE_CENTER_VOLUME_SLICE = slice(28, 65)
SOURCE_CENTER_INTERNAL_JUNCTION_SLICE = slice(25, 61)
SOURCE_CENTER_MACRO_IN_JUNCTION_INDEX = 169
SOURCE_CENTER_MACRO_OUT_JUNCTION_INDEX = 170
SOURCE_TFE_PREFIX = "Macro_TASTIN_Core_V7_CaseA/TFEs/Center"

SOLID_SOURCE_PREFIXES = {
    "pellet": "Solid_Center_Pellet",
    "emitter": "Solid_Center_Emitter",
    "collector": "Solid_Center_Collector",
    "inner_clad": "Solid_Center_InnerClad",
    "outer_clad": "Solid_Center_OuterClad",
}


@dataclass
class SnapshotMapping:
    inlet_plenum_p_pa: float
    inlet_plenum_t_k: float
    inlet_plenum_h_j_kg: float
    outlet_plenum_p_pa: float
    outlet_plenum_t_k: float
    outlet_plenum_h_j_kg: float
    center_p_pa: np.ndarray
    center_t_k: np.ndarray
    center_h_j_kg: np.ndarray
    center_internal_w_kg_s: np.ndarray
    center_macro_in_w_kg_s: float
    center_macro_out_w_kg_s: float
    solid_temperatures_k: Dict[str, np.ndarray]


def _build_axial_power_profile(n_lower: int, n_active: int, n_upper: int) -> np.ndarray:
    z_squared = np.linspace(-1.0, 1.0, n_active) ** 2
    active_profile = (
        AXIAL_SHAPE_COEFFS[0]
        + AXIAL_SHAPE_COEFFS[1] * z_squared
        + AXIAL_SHAPE_COEFFS[2] * z_squared**2
        + AXIAL_SHAPE_COEFFS[3] * z_squared**3
        + AXIAL_SHAPE_COEFFS[4] * z_squared**4
    )
    active_profile = np.maximum(active_profile, 0.0)
    active_profile /= np.sum(active_profile)
    return np.concatenate((np.zeros(n_lower), active_profile, np.zeros(n_upper)))


def _resolve_path(path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute() or path.exists():
        return path
    return ROOT_DIR / path


def load_center_snapshot(snapshot_path: Path) -> SnapshotMapping:
    if not snapshot_path.exists():
        raise FileNotFoundError(f"Source snapshot not found: {snapshot_path}")

    with np.load(snapshot_path, allow_pickle=False) as data:
        fluid_shape = tuple(int(value) for value in data["Fluid/shape"])
        if fluid_shape != SOURCE_FLUID_SHAPE:
            raise ValueError(
                f"Expected source Fluid/shape {SOURCE_FLUID_SHAPE}, got {fluid_shape}."
            )
        source_time = float(data["System/global_time"][0])
        if not np.isclose(source_time, SOURCE_SNAPSHOT_TIME_S):
            raise ValueError(
                f"Expected source snapshot time {SOURCE_SNAPSHOT_TIME_S}, got {source_time}."
            )
        source_n_axial = int(data[f"{SOURCE_TFE_PREFIX}/n_axial"][0])
        if source_n_axial != N_AXIAL:
            raise ValueError(f"Expected {N_AXIAL} axial nodes, got {source_n_axial}.")

        source_power = float(data[f"{SOURCE_TFE_PREFIX}/neutronic/total_power"][0])
        if not np.isclose(source_power, FIXED_NUCLEAR_HEAT_W, rtol=0.0, atol=1.0e-9):
            raise ValueError(
                f"Expected Center TFE power {FIXED_NUCLEAR_HEAT_W}, got {source_power}."
            )

        p_vec = np.asarray(data["Fluid/P_vec"], dtype=float)
        t_vec = np.asarray(data["Fluid/T_vec"], dtype=float)
        h_vec = np.asarray(data["Fluid/h_vec"], dtype=float)
        w_vec = np.asarray(data["Fluid/W_vec"], dtype=float)
        solids = {
            name: np.asarray(data[f"{prefix}/T"], dtype=float).copy()
            for name, prefix in SOLID_SOURCE_PREFIXES.items()
        }

    center_p = p_vec[SOURCE_CENTER_VOLUME_SLICE].copy()
    center_t = t_vec[SOURCE_CENTER_VOLUME_SLICE].copy()
    center_h = h_vec[SOURCE_CENTER_VOLUME_SLICE].copy()
    center_w = w_vec[SOURCE_CENTER_INTERNAL_JUNCTION_SLICE].copy()
    for label, array in (
        ("center pressure", center_p),
        ("center temperature", center_t),
        ("center enthalpy", center_h),
    ):
        if array.shape != (N_AXIAL,):
            raise ValueError(f"{label} shape must be {(N_AXIAL,)}, got {array.shape}.")
    if center_w.shape != (N_AXIAL - 1,):
        raise ValueError(
            f"center internal flow shape must be {(N_AXIAL - 1,)}, got {center_w.shape}."
        )

    return SnapshotMapping(
        inlet_plenum_p_pa=float(p_vec[SOURCE_INLET_PLENUM_VOLUME_INDEX]),
        inlet_plenum_t_k=float(t_vec[SOURCE_INLET_PLENUM_VOLUME_INDEX]),
        inlet_plenum_h_j_kg=float(h_vec[SOURCE_INLET_PLENUM_VOLUME_INDEX]),
        outlet_plenum_p_pa=float(p_vec[SOURCE_OUTLET_PLENUM_VOLUME_INDEX]),
        outlet_plenum_t_k=float(t_vec[SOURCE_OUTLET_PLENUM_VOLUME_INDEX]),
        outlet_plenum_h_j_kg=float(h_vec[SOURCE_OUTLET_PLENUM_VOLUME_INDEX]),
        center_p_pa=center_p,
        center_t_k=center_t,
        center_h_j_kg=center_h,
        center_internal_w_kg_s=center_w,
        center_macro_in_w_kg_s=float(w_vec[SOURCE_CENTER_MACRO_IN_JUNCTION_INDEX]),
        center_macro_out_w_kg_s=float(w_vec[SOURCE_CENTER_MACRO_OUT_JUNCTION_INDEX]),
        solid_temperatures_k=solids,
    )


def _build_single_tfe(mapping: SnapshotMapping) -> Dict[str, object]:
    l_lower, l_active, l_upper = 0.065, 0.377, 0.065
    n_lower, n_active, n_upper = 6, 25, 6
    node_lengths = np.array(
        [l_lower / n_lower] * n_lower
        + [l_active / n_active] * n_active
        + [l_upper / n_upper] * n_upper,
        dtype=float,
    )
    geometry = TFEGeometry(
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
        height=l_lower + l_active + l_upper,
    )
    mesh = TFEMeshParams(
        n_axial=N_AXIAL,
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=1,
        n_r_outer_clad=1,
        n_r_moderator=3,
    )
    sodium = Sodium()
    flow_area = np.pi * (geometry.r_coolant_outer**2 - geometry.r_coolant_inner**2)
    hydraulic_diam = 2.0 * (geometry.r_coolant_outer - geometry.r_coolant_inner)

    inlet = IncompressibleBoundaryVolume(
        name="SingleTFE_Inlet",
        material=sodium,
        P=mapping.inlet_plenum_p_pa,
        T=mapping.inlet_plenum_t_k,
        flow_area=flow_area,
        hydraulic_diam=hydraulic_diam,
    )
    outlet = IncompressibleBoundaryVolume(
        name="SingleTFE_Outlet",
        material=sodium,
        P=FIXED_OUTLET_PRESSURE_PA,
        T=mapping.outlet_plenum_t_k,
        flow_area=flow_area,
        hydraulic_diam=hydraulic_diam,
    )
    outlet.is_pressure_boundary = True
    channel = NonUniformIncompressibleFluidChannel(
        name="SingleTFE_Channel",
        node_lengths=node_lengths,
        flow_area=flow_area,
        hydraulic_diam=hydraulic_diam,
        initial_P=mapping.inlet_plenum_p_pa,
        initial_T=mapping.inlet_plenum_t_k,
        material=sodium,
    )
    inlet_junction = InletJunction(
        name="SingleTFE_InletFlow",
        from_vol=inlet,
        to_vol=channel.volumes[0],
        W_initial=FIXED_INLET_MASS_FLOW_KG_S,
    )
    outlet_junction = FlowJunction(
        name="SingleTFE_OutletPressure",
        from_vol=channel.volumes[-1],
        to_vol=outlet,
        flow_area=flow_area,
    )
    network = HydraulicNetwork(
        volumes=[inlet, *channel.volumes, outlet],
        junctions=[inlet_junction, *channel.internal_junctions, outlet_junction],
        gravity_vector=0.0,
    )

    tfe = TFEUnit(
        name="SingleTFE_Center",
        geometry=geometry,
        mesh_params=mesh,
        materials={
            "UO2": UO2(),
            "MoNb": MoNb(),
            "Molybdenum": Molybdenum(),
            "StainlessSteel": AusteniticStainlessSteel(),
            "ZrH": ZirconiumHydride(),
            "BerylliumOxide": BerylliumOxide(),
        },
        coolant_channel=channel,
        fission_gas_config=GapConfig(
            mode="simplified",
            h_eq=5678.0,
            material=Xenon(),
            emissivity_inner=0.15,
            emissivity_outer=0.15,
        ),
        tec_gap_config=GapConfig(
            mode="simplified",
            h_eq=29.0,
            material=Cesium(),
            emissivity_inner=0.15,
            emissivity_outer=0.60,
        ),
        he_gap_config=GapConfig(
            mode="simplified",
            h_eq=5678.0,
            material=Helium(),
            emissivity_inner=0.60,
            emissivity_outer=0.80,
        ),
        co2_gap_config=GapConfig(
            mode="simplified",
            h_eq=53.6,
            material=CarbonDioxide(),
            emissivity_inner=0.80,
            emissivity_outer=0.80,
        ),
        axial_power_profile=_build_axial_power_profile(n_lower, n_active, n_upper),
        axial_length_allocation=[l_lower, l_active, l_upper],
        axial_node_allocation=[n_lower, n_active, n_upper],
        strict_adiabatic_single_tfe=True,
    )
    system = SystemManager(fluid_network=network, start_time=0.0)
    system.add_component(tfe)
    system.initialize_system()

    return {
        "system": system,
        "network": network,
        "tfe": tfe,
        "channel": channel,
        "inlet": inlet,
        "outlet": outlet,
        "inlet_junction": inlet_junction,
        "outlet_junction": outlet_junction,
        "flow_area_m2": flow_area,
    }


def _reset_and_map_snapshot(case: Dict[str, object], mapping: SnapshotMapping) -> None:
    system = case["system"]
    network = case["network"]
    tfe = case["tfe"]
    inlet = case["inlet"]
    outlet = case["outlet"]
    inlet_junction = case["inlet_junction"]

    network.P_vec[:] = np.concatenate(
        ([mapping.inlet_plenum_p_pa], mapping.center_p_pa, [FIXED_OUTLET_PRESSURE_PA])
    )
    network.T_vec[:] = np.concatenate(
        ([mapping.inlet_plenum_t_k], mapping.center_t_k, [mapping.outlet_plenum_t_k])
    )
    network.h_vec[:] = np.concatenate(
        ([mapping.inlet_plenum_h_j_kg], mapping.center_h_j_kg, [mapping.outlet_plenum_h_j_kg])
    )
    network.W_vec[:] = np.concatenate(
        (
            [FIXED_INLET_MASS_FLOW_KG_S],
            mapping.center_internal_w_kg_s,
            [mapping.center_macro_out_w_kg_s],
        )
    )
    inlet_junction.W = FIXED_INLET_MASS_FLOW_KG_S
    inlet_junction.target_W = FIXED_INLET_MASS_FLOW_KG_S
    inlet.P = mapping.inlet_plenum_p_pa
    inlet.T = mapping.inlet_plenum_t_k
    inlet.h = mapping.inlet_plenum_h_j_kg
    outlet.P = FIXED_OUTLET_PRESSURE_PA
    outlet.target_P = FIXED_OUTLET_PRESSURE_PA
    outlet.T = mapping.outlet_plenum_t_k
    outlet.h = mapping.outlet_plenum_h_j_kg

    network._refresh_cached_pressure_targets()
    network._refresh_cached_boundary_targets()
    network._update_fluid_properties()
    network._sync_vectors_to_objects()
    network.W_old[:] = network.W_vec
    network.W_iterate[:] = network.W_vec

    for name, solid in tfe.solids.items():
        if name not in mapping.solid_temperatures_k:
            raise KeyError(f"Unexpected single-TFE solid '{name}' in snapshot mapper.")
        source_t = mapping.solid_temperatures_k[name]
        if source_t.shape != solid.T.shape:
            raise ValueError(
                f"Solid '{name}' temperature shape mismatch: {source_t.shape} != {solid.T.shape}."
            )
        solid.T[:] = source_t
        solid.dTdt[:] = 0.0
        solid.Q_source[:] = 0.0
        solid.use_external_source_buffer = False
        solid.source_callback = None
        solid.current_time = 0.0
        for boundary in solid.boundaries.values():
            boundary.clear_boundary_conditions()
        solid.initialize_state()

    for attr in (
        "emitter_voltage",
        "emitter_resistivity",
        "collector_voltage",
        "collector_resistivity",
        "current_density",
        "emitter_joule_heat",
        "collector_joule_heat",
    ):
        getattr(tfe.electric_data, attr).fill(0.0)
    for attr in (
        "emitter_work_function",
        "collector_work_function",
        "barrier_voltage_drop",
        "emitter_temperature",
        "electron_cooling_flux",
        "electron_heating_flux",
    ):
        getattr(tfe.plasma_data, attr).fill(0.0)
    tfe.couplers["tec_couple"].set_tec_sources(
        Q_emitter=np.zeros(N_AXIAL),
        Q_collector=np.zeros(N_AXIAL),
    )
    tfe.update_neutronic_power(
        p_total=FIXED_NUCLEAR_HEAT_W,
        p_fiss=FIXED_NUCLEAR_HEAT_W,
        p_decay=0.0,
        alpha=1.0,
    )

    system.global_time = 0.0
    system._sync_solid_times_to_global()
    system._prepare_fluid_sources_for_coupling()
    system._run_couplers(current_time=0.0)
    system._refresh_solid_boundary_cache(update_flux=True, current_time=0.0)
    network.save_state()


def _assert_strict_adiabatic(tfe: TFEUnit) -> None:
    boundary = tfe.solids["outer_clad"].boundaries["right"]
    if boundary.conditions:
        raise AssertionError("Strict adiabatic outer-clad boundary still has boundary conditions.")
    flux = np.asarray(boundary.compute_net_flux_for_solver(), dtype=float)
    if not np.array_equal(flux, np.zeros_like(flux)):
        raise AssertionError(f"Strict adiabatic outer-clad flux is not exactly zero: {flux}.")


def _capture_storage_state(case: Dict[str, object]) -> Dict[str, object]:
    network = case["network"]
    tfe = case["tfe"]
    return {
        "solid": {
            name: (solid.T.copy(), solid.thermal_capacitance.copy())
            for name, solid in tfe.solids.items()
        },
        "fluid_h": network.h_vec.copy(),
        "fluid_mass": network.rho_vec * network.V_vec,
    }


def _storage_rates(case: Dict[str, object], before: Dict[str, object], dt: float) -> Dict[str, float]:
    network = case["network"]
    tfe = case["tfe"]
    channel_slice = slice(1, 1 + N_AXIAL)
    solid_storage_w = 0.0
    for name, solid in tfe.solids.items():
        old_t, old_cap = before["solid"][name]
        solid_storage_w += float(
            np.sum(0.5 * (old_cap + solid.thermal_capacitance) * (solid.T - old_t)) / dt
        )

    new_mass = network.rho_vec[channel_slice] * network.V_vec[channel_slice]
    fluid_storage_w = float(
        np.sum(
            0.5
            * (before["fluid_mass"][channel_slice] + new_mass)
            * (network.h_vec[channel_slice] - before["fluid_h"][channel_slice])
        )
        / dt
    )
    return {
        "solid_storage_rate_w": solid_storage_w,
        "fluid_storage_rate_w": fluid_storage_w,
    }


def _outer_clad_heat_loss_w(tfe: TFEUnit) -> float:
    flux_into_solid = tfe.solids["outer_clad"].boundaries["right"].compute_net_flux_for_solver()
    return -float(np.sum(flux_into_solid))


def _synchronize_latest_audit(case: Dict[str, object]) -> None:
    system = case["system"]
    current_time = float(system.global_time)
    system._prepare_fluid_sources_for_coupling()
    system._run_couplers(interface_relaxation=1.0, current_time=current_time)
    system._refresh_solid_boundary_cache(update_flux=True, current_time=current_time)


def _sum_boundary_fluxes(solid) -> Dict[str, float]:
    return {
        location: float(np.sum(boundary.current_flux))
        for location, boundary in solid.boundaries.items()
    }


def _collect_solid_energy_balance(
    case: Dict[str, object],
    before: Optional[Dict[str, object]],
    dt: Optional[float],
):
    rows = []
    max_abs_fd_residual = 0.0
    max_abs_solver_residual = 0.0
    for solid_name, solid in case["tfe"].solids.items():
        solid.get_derivatives(float(case["system"].global_time), solid.T.copy())
        boundary_fluxes = _sum_boundary_fluxes(solid)
        boundary_inflow = float(sum(boundary_fluxes.values()))
        source = float(np.sum(solid.Q_source))
        solver_storage = float(np.sum(solid.thermal_capacitance * solid.dTdt))
        solver_residual = source + boundary_inflow - solver_storage
        fd_storage = None
        fd_residual = None
        if before is not None and dt is not None:
            old_t, old_cap = before["solid"][solid_name]
            fd_storage = float(
                np.sum(0.5 * (old_cap + solid.thermal_capacitance) * (solid.T - old_t)) / dt
            )
            fd_residual = source + boundary_inflow - fd_storage
            max_abs_fd_residual = max(max_abs_fd_residual, abs(fd_residual))
        max_abs_solver_residual = max(max_abs_solver_residual, abs(solver_residual))
        rows.append({
            "solid": solid_name,
            "source_w": source,
            "boundary_inflow_w": boundary_inflow,
            "boundary_left_inflow_w": boundary_fluxes.get("left", 0.0),
            "boundary_right_inflow_w": boundary_fluxes.get("right", 0.0),
            "boundary_bottom_inflow_w": boundary_fluxes.get("bottom", 0.0),
            "boundary_top_inflow_w": boundary_fluxes.get("top", 0.0),
            "fd_storage_w": fd_storage,
            "solver_cap_dtdt_storage_w": solver_storage,
            "fd_residual_w": fd_residual,
            "solver_residual_w": solver_residual,
        })
    return rows, {
        "max_abs_fd_residual_w": max_abs_fd_residual,
        "max_abs_solver_residual_w": max_abs_solver_residual,
        "sum_fd_residual_w": float(sum(row["fd_residual_w"] or 0.0 for row in rows)),
        "sum_solver_residual_w": float(sum(row["solver_residual_w"] for row in rows)),
    }


def _collect_interface_balance(case: Dict[str, object]):
    tfe = case["tfe"]
    rows = []
    ordinary_sum_abs = 0.0
    ordinary_max_abs = 0.0
    tec_boundary_inflow_sum = 0.0
    fluid_solid_sum_abs = 0.0
    fluid_solid_max_abs = 0.0
    merged_fluid_source = np.array([
        volume.Q_wall + volume.Q_vol - volume.implicit_coeff * volume.T
        for volume in case["channel"].volumes
    ], dtype=float)
    mapped_fluid_source = np.zeros_like(merged_fluid_source)

    for coupler_name, coupler in tfe.couplers.items():
        if isinstance(coupler, FluidSolidCouple):
            lam = np.asarray(coupler._last_lambda, dtype=float)
            wall_t = np.asarray(coupler.solid_bound.T_surface, dtype=float)
            fluid_t = np.asarray(coupler.fluid.temperature_vector, dtype=float)
            fluid_source = lam * (wall_t - fluid_t)
            solid_to_fluid = -np.asarray(coupler.solid_bound.current_flux, dtype=float)
            mapped_fluid_source += fluid_source
            for axial_index, (solid_w, fluid_w) in enumerate(zip(solid_to_fluid, fluid_source)):
                residual = float(fluid_w - solid_w)
                fluid_solid_sum_abs += abs(residual)
                fluid_solid_max_abs = max(fluid_solid_max_abs, abs(residual))
                rows.append({
                    "interface_kind": "fluid_solid",
                    "coupler": coupler_name,
                    "axial_index": axial_index,
                    "side1": getattr(coupler.solid_bound, "name", "solid_boundary"),
                    "side2": getattr(coupler.fluid, "name", "fluid"),
                    "side1_inflow_w": -float(solid_w),
                    "side2_inflow_w": float(fluid_w),
                    "inflow_sum_residual_w": residual,
                })
            continue

        if not hasattr(coupler, "bound1") or not hasattr(coupler, "bound2"):
            continue
        q1 = np.asarray(coupler.bound1.current_flux, dtype=float)
        q2 = np.asarray(coupler.bound2.current_flux, dtype=float)
        for axial_index, (q1_w, q2_w) in enumerate(zip(q1, q2)):
            residual = float(q1_w + q2_w)
            if isinstance(coupler, TECCouple2D):
                tec_boundary_inflow_sum += residual
            else:
                ordinary_sum_abs += abs(residual)
                ordinary_max_abs = max(ordinary_max_abs, abs(residual))
            rows.append({
                "interface_kind": "tec_solid_solid" if isinstance(coupler, TECCouple2D) else "solid_solid",
                "coupler": coupler_name,
                "axial_index": axial_index,
                "side1": getattr(coupler.obj1, "name", ""),
                "side2": getattr(coupler.obj2, "name", ""),
                "side1_inflow_w": float(q1_w),
                "side2_inflow_w": float(q2_w),
                "inflow_sum_residual_w": residual,
            })

    merged_delta = merged_fluid_source - mapped_fluid_source
    for axial_index, (mapped_w, merged_w, residual) in enumerate(
        zip(mapped_fluid_source, merged_fluid_source, merged_delta)
    ):
        rows.append({
            "interface_kind": "fluid_source_merge",
            "coupler": "all_fluid_solid",
            "axial_index": axial_index,
            "side1": "fluid_solid_couplers",
            "side2": "channel_volume_source",
            "side1_inflow_w": float(mapped_w),
            "side2_inflow_w": -float(merged_w),
            "inflow_sum_residual_w": float(residual),
        })
    return rows, {
        "ordinary_solid_solid_sum_abs_residual_w": ordinary_sum_abs,
        "ordinary_solid_solid_max_abs_node_residual_w": ordinary_max_abs,
        "tec_boundary_inflow_sum_w": tec_boundary_inflow_sum,
        "fluid_solid_sum_abs_residual_w": fluid_solid_sum_abs,
        "fluid_solid_max_abs_node_residual_w": fluid_solid_max_abs,
        "fluid_source_merge_sum_abs_residual_w": float(np.sum(np.abs(merged_delta))),
        "fluid_source_merge_max_abs_node_residual_w": float(np.max(np.abs(merged_delta))),
    }


def _collect_fluid_volume_balance(
    case: Dict[str, object],
    before: Optional[Dict[str, object]],
    dt: Optional[float],
    applied_sources,
):
    network = case["network"]
    inflow = np.zeros(network.n_vol)
    outflow = np.zeros(network.n_vol)
    for j_idx, (_, idx_from, idx_to) in enumerate(network.junction_descriptors):
        flow = float(network.W_vec[j_idx])
        if flow >= 0.0:
            enthalpy_flow = flow * float(network.h_vec[idx_from])
            outflow[idx_from] += float(network.M_from_vec[j_idx]) * enthalpy_flow
            inflow[idx_to] += float(network.M_to_vec[j_idx]) * enthalpy_flow
        else:
            enthalpy_flow = -flow * float(network.h_vec[idx_to])
            inflow[idx_from] += float(network.M_from_vec[j_idx]) * enthalpy_flow
            outflow[idx_to] += float(network.M_to_vec[j_idx]) * enthalpy_flow

    matrix_residual = np.asarray(
        network.energy_matrix.dot(network.h_vec) - network.energy_rhs_buffer,
        dtype=float,
    )
    new_mass = network.rho_vec * network.V_vec
    rows = []
    finite_matrix_abs_sum = 0.0
    finite_matrix_max_abs = 0.0
    for index, volume in enumerate(network.volumes_obj):
        fixed_boundary = bool(getattr(volume, "is_pressure_boundary", False))
        boundary_volume = isinstance(volume, IncompressibleBoundaryVolume)
        latest_effective_source = float(volume.Q_wall + volume.Q_vol - volume.implicit_coeff * volume.T)
        q_wall, q_vol, implicit_coeff = applied_sources[index]
        effective_source = float(q_wall + q_vol - implicit_coeff * volume.T)
        fd_storage = None
        fd_residual = None
        if before is not None and dt is not None:
            fd_storage = float(
                0.5
                * (before["fluid_mass"][index] + new_mass[index])
                * (network.h_vec[index] - before["fluid_h"][index])
                / dt
            )
            fd_residual = float(inflow[index] - outflow[index] + effective_source - fd_storage)
        if not boundary_volume:
            finite_matrix_abs_sum += abs(float(matrix_residual[index]))
            finite_matrix_max_abs = max(finite_matrix_max_abs, abs(float(matrix_residual[index])))
        rows.append({
            "volume_index": index,
            "volume": getattr(volume, "name", f"volume_{index}"),
            "is_boundary_volume": boundary_volume,
            "is_fixed_pressure_boundary": fixed_boundary,
            "enthalpy_inflow_w": float(inflow[index]),
            "enthalpy_outflow_w": float(outflow[index]),
            "applied_wall_source_w": float(q_wall),
            "applied_volumetric_source_w": float(q_vol),
            "applied_implicit_coeff_w_per_k": float(implicit_coeff),
            "applied_effective_source_w": effective_source,
            "latest_synced_effective_source_w": latest_effective_source,
            "fd_storage_w": fd_storage,
            "fd_residual_w": fd_residual,
            "solver_matrix_residual_w": float(matrix_residual[index]),
        })
    return rows, {
        "finite_volume_sum_abs_matrix_residual_w": finite_matrix_abs_sum,
        "finite_volume_max_abs_matrix_residual_w": finite_matrix_max_abs,
        "finite_volume_sum_fd_residual_w": float(sum(
            row["fd_residual_w"] or 0.0
            for row in rows
            if not row["is_boundary_volume"]
        )),
    }


def _tail_residual_windows(history):
    if not history:
        return {}
    final_time = float(history[-1]["time_s"])
    windows = {}
    for width_s in (1.0, 2.0, 5.0):
        selected = [row for row in history if row["time_s"] > final_time - width_s]
        windows[f"last_{int(width_s)}s"] = {
            "row_count": len(selected),
            "mean_residual_w": float(np.mean([row["residual_w"] for row in selected])),
            "mean_relative_residual": float(np.mean([row["relative_residual"] for row in selected])),
            "mean_abs_relative_residual": float(
                np.mean([abs(row["relative_residual"]) for row in selected])
            ),
        }
    return windows


def _build_tec_model(tfe: TFEUnit):
    bindings_path = ROOT_DIR / "ThermoCalc" / "bindings.cpp"
    bindings_source = bindings_path.read_text(encoding="utf-8")
    failures = []
    if "get_scalar(data.sideAreaE, i)" in bindings_source or "get_scalar(data.sideAreaC, i)" in bindings_source:
        failures.append("bindings.cpp still collapses sideAreaE/sideAreaC to one scalar per TFE")
    for field in ("phiE", "phiC", "Vd", "joulePowerE", "joulePowerC"):
        if f'.def_readwrite("{field}"' not in bindings_source:
            failures.append(f"SingleTEC binding does not expose {field}")
    if failures:
        raise RuntimeError("TEC capability check failed: " + "; ".join(failures))

    try:
        from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

        model = ThermoCalcModel(n_elements=1, n_nodes=N_AXIAL)
        node_lengths = np.diff(tfe.common_y_faces)
        emitter_area = np.asarray(tfe.solids["emitter"].boundaries["right"].area, dtype=float)
        collector_area = np.asarray(tfe.solids["collector"].boundaries["left"].area, dtype=float)
        model._input_data.dlE = node_lengths[np.newaxis, :]
        model._input_data.dlC = node_lengths[np.newaxis, :]
        model._input_data.sideAreaE = emitter_area[np.newaxis, :]
        model._input_data.sideAreaC = collector_area[np.newaxis, :]
        model._input_data.resistanceWire = TEC_WIRE_RESISTANCE_OHM[np.newaxis, :].copy()
        model._input_data.wireU = np.array([[TEC_TARGET_VOLTAGE_V, TEC_TARGET_VOLTAGE_V, 0.0, 0.0]])
        model.setup_circuit_mode("fixed_u", TEC_TARGET_VOLTAGE_V)
        model.set_temperatures(
            tfe.solids["emitter"].boundaries["right"].T_surface[np.newaxis, :],
            tfe.solids["collector"].boundaries["left"].T_surface[np.newaxis, :],
        )
        model.calculate(verbose=False)
        results = model.get_tec_results(0)
        for field in ("phiE", "phiC", "Vd", "joulePowerE", "joulePowerC"):
            value = np.asarray(results[field], dtype=float)
            if value.shape != (N_AXIAL,):
                raise RuntimeError(f"TEC field {field} has shape {value.shape}, expected {(N_AXIAL,)}.")
            if not np.all(np.isfinite(value)):
                raise RuntimeError(f"TEC field {field} contains non-finite values.")
        if np.any(np.asarray(results["joulePowerE"]) < 0.0) or np.any(np.asarray(results["joulePowerC"]) < 0.0):
            raise RuntimeError("TEC Joule power fields must be non-negative.")
        return model
    except Exception as exc:
        raise RuntimeError(f"TEC capability check failed at runtime: {type(exc).__name__}: {exc}") from exc


def _apply_tec_sources(model, tfe: TFEUnit) -> Dict[str, object]:
    model.set_temperatures(
        tfe.solids["emitter"].boundaries["right"].T_surface[np.newaxis, :],
        tfe.solids["collector"].boundaries["left"].T_surface[np.newaxis, :],
    )
    model.calculate(verbose=False)
    results = model.get_tec_results(0)
    emitter_potential = np.asarray(results["UE"], dtype=float)
    collector_potential = np.asarray(results["UC"], dtype=float)
    current_density_a_m2 = np.asarray(results["J"], dtype=float) * 1.0e4
    emitter_temperature = np.asarray(results["TE"], dtype=float)
    phi_e = np.asarray(results["phiE"], dtype=float)
    q_e_flux = -current_density_a_m2 * (phi_e + 2.0 * 8.617e-5 * emitter_temperature)
    q_c_flux = current_density_a_m2 * (
        phi_e + 2.0 * 8.617e-5 * emitter_temperature - (emitter_potential - collector_potential)
    )
    tfe.update_electric_field_diagnostics(
        E_emit=electric_field_from_node_potential(emitter_potential, y_faces=tfe.common_y_faces),
        rho_emit=np.asarray(results["rhoE"], dtype=float),
        E_coll=electric_field_from_node_potential(collector_potential, y_faces=tfe.common_y_faces),
        rho_coll=np.asarray(results["rhoC"], dtype=float),
    )
    tfe.update_joule_power_sources(
        Q_emitter_axial=np.asarray(results["joulePowerE"], dtype=float),
        Q_collector_axial=np.asarray(results["joulePowerC"], dtype=float),
        alpha=1.0,
    )
    tfe.update_plasma_flux(q_e_flux=q_e_flux, q_c_flux=q_c_flux, alpha=1.0)
    return results


def _audit_step(
    case: Dict[str, object],
    before: Dict[str, object],
    dt: float,
    electrical_output_w: float,
) -> Dict[str, float]:
    network = case["network"]
    tfe = case["tfe"]
    storage = _storage_rates(case, before, dt)
    inlet_w = float(network.W_vec[0])
    outlet_w = float(network.W_vec[-1])
    inlet_h = float(network.h_vec[0])
    outlet_h = float(network.h_vec[N_AXIAL])
    coolant_pickup_w = outlet_w * outlet_h - inlet_w * inlet_h
    outer_loss_w = _outer_clad_heat_loss_w(tfe)
    nuclear_heat_w = float(np.sum(tfe.solids["pellet"].Q_source))
    residual_w = (
        nuclear_heat_w
        - coolant_pickup_w
        - electrical_output_w
        - storage["solid_storage_rate_w"]
        - storage["fluid_storage_rate_w"]
        - outer_loss_w
    )
    return {
        "time_s": float(case["system"].global_time),
        "dt_s": float(dt),
        "nuclear_heat_w": nuclear_heat_w,
        "coolant_enthalpy_pickup_w": coolant_pickup_w,
        "electrical_output_w": float(electrical_output_w),
        **storage,
        "outer_clad_heat_loss_w": outer_loss_w,
        "residual_w": residual_w,
        "relative_residual": residual_w / nuclear_heat_w if nuclear_heat_w else 0.0,
        "inlet_mass_flow_kg_s": inlet_w,
        "outlet_mass_flow_kg_s": outlet_w,
    }


def _write_csv(path: Path, rows, fieldnames=None) -> None:
    rows = list(rows)
    if not rows and fieldnames is None:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fieldnames or rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_profiles(output_dir: Path, case: Dict[str, object]) -> None:
    network = case["network"]
    channel = case["channel"]
    tfe = case["tfe"]
    _write_csv(
        output_dir / "fluid_profile_latest.csv",
        (
            {
                "axial_index": index,
                "z_center_m": volume.z_coordinate,
                "pressure_pa": network.P_vec[index + 1],
                "temperature_k": network.T_vec[index + 1],
                "enthalpy_j_kg": network.h_vec[index + 1],
                "mass_flow_out_kg_s": network.W_vec[min(index + 1, N_AXIAL)],
            }
            for index, volume in enumerate(channel.volumes)
        ),
    )
    rows = []
    for solid_name, solid in tfe.solids.items():
        nx, ny = solid.shape_nodes
        temperature = solid.T.reshape(nx, ny)
        for radial_index in range(nx):
            for axial_index in range(ny):
                rows.append(
                    {
                        "solid": solid_name,
                        "radial_index": radial_index,
                        "axial_index": axial_index,
                        "temperature_k": temperature[radial_index, axial_index],
                    }
                )
    _write_csv(output_dir / "solid_temperature_profile_latest.csv", rows)


def _write_tec_nodes(
    output_dir: Path,
    tfe: TFEUnit,
    results: Dict[str, object],
    global_results: Dict[str, float],
):
    area = tfe.plasma_area_diagnostics
    emitter = tfe.solids["emitter"]
    collector = tfe.solids["collector"]
    emitter_joule = np.asarray(tfe.electric_data.emitter_joule_heat, dtype=float).reshape(
        emitter.mesh.shape_nodes
    ).sum(axis=0)
    collector_joule = np.asarray(tfe.electric_data.collector_joule_heat, dtype=float).reshape(
        collector.mesh.shape_nodes
    ).sum(axis=0)
    cpp_emitter_joule = np.asarray(results["joulePowerE"], dtype=float)
    cpp_collector_joule = np.asarray(results["joulePowerC"], dtype=float)
    legacy_emitter_flat, _ = joule_power_from_electric_field(
        tfe.electric_data.emitter_voltage,
        tfe.electric_data.emitter_resistivity,
        emitter.vols_flat,
        emitter.mesh.shape_nodes,
    )
    legacy_collector_flat, _ = joule_power_from_electric_field(
        tfe.electric_data.collector_voltage,
        tfe.electric_data.collector_resistivity,
        collector.vols_flat,
        collector.mesh.shape_nodes,
    )
    legacy_emitter_joule = legacy_emitter_flat.reshape(emitter.mesh.shape_nodes).sum(axis=0)
    legacy_collector_joule = legacy_collector_flat.reshape(collector.mesh.shape_nodes).sum(axis=0)
    emitter_electron = np.asarray(area["emitter_power_emitter_area_w"], dtype=float)
    collector_electron = np.asarray(area["collector_power_emitter_area_w"], dtype=float)
    electron_boundary_diff = -emitter_electron - collector_electron
    joule = emitter_joule + collector_joule
    cpp_joule = cpp_emitter_joule + cpp_collector_joule
    legacy_joule = legacy_emitter_joule + legacy_collector_joule
    mapped_minus_cpp = joule - cpp_joule
    legacy_minus_cpp = legacy_joule - cpp_joule
    node_current = np.asarray(results["J"], dtype=float) * 1.0e4 * np.asarray(area["sideAreaE_m2"])
    weights = np.abs(node_current)
    if float(np.sum(weights)) <= 0.0:
        weights = np.ones_like(weights)
    terminal_power = float(global_results["Uout"]) * float(global_results["Iout"])
    terminal_alloc = terminal_power * weights / float(np.sum(weights))
    requested_closure = electron_boundary_diff - (terminal_alloc - joule)
    plus_joule_closure = electron_boundary_diff - (terminal_alloc + joule)
    cpp_plus_joule_closure = electron_boundary_diff - (terminal_alloc + cpp_joule)
    _write_csv(
        output_dir / "tec_node_balance_latest.csv",
        (
            {
                "axial_index": index,
                "current_density_a_m2": np.asarray(results["J"])[index] * 1.0e4,
                "node_current_a": node_current[index],
                "emitter_potential_v": np.asarray(results["UE"])[index],
                "collector_potential_v": np.asarray(results["UC"])[index],
                "phi_e_v": np.asarray(results["phiE"])[index],
                "phi_c_v": np.asarray(results["phiC"])[index],
                "vd_v": np.asarray(results["Vd"])[index],
                "emitter_electron_power_w": emitter_electron[index],
                "collector_electron_power_w": collector_electron[index],
                "electron_boundary_power_diff_w": electron_boundary_diff[index],
                "emitter_joule_power_w": emitter_joule[index],
                "collector_joule_power_w": collector_joule[index],
                "joule_power_w": joule[index],
                "cpp_emitter_joule_power_w": cpp_emitter_joule[index],
                "cpp_collector_joule_power_w": cpp_collector_joule[index],
                "cpp_joule_power_w": cpp_joule[index],
                "mapped_minus_cpp_joule_power_w": mapped_minus_cpp[index],
                "legacy_gradient_joule_power_w": legacy_joule[index],
                "legacy_gradient_minus_cpp_joule_power_w": legacy_minus_cpp[index],
                "terminal_power_alloc_w": terminal_alloc[index],
                "electron_diff_minus_terminal_less_joule_w": requested_closure[index],
                "electron_diff_minus_terminal_plus_joule_w": plus_joule_closure[index],
                "electron_diff_minus_terminal_plus_cpp_joule_w": cpp_plus_joule_closure[index],
            }
            for index in range(N_AXIAL)
        ),
    )
    return {
        "single_tec_u_field_v": float(results["U"]),
        "single_tec_current_a": float(results["I"]),
        "terminal_voltage_v": float(global_results["Uout"]),
        "terminal_current_a": float(global_results["Iout"]),
        "terminal_power_w": terminal_power,
        "emitter_electron_power_w": float(np.sum(emitter_electron)),
        "collector_electron_power_w": float(np.sum(collector_electron)),
        "electron_boundary_power_diff_w": float(np.sum(electron_boundary_diff)),
        "emitter_joule_power_w": float(np.sum(emitter_joule)),
        "collector_joule_power_w": float(np.sum(collector_joule)),
        "joule_power_w": float(np.sum(joule)),
        "cpp_emitter_joule_power_w": float(np.sum(cpp_emitter_joule)),
        "cpp_collector_joule_power_w": float(np.sum(cpp_collector_joule)),
        "cpp_joule_power_w": float(np.sum(cpp_joule)),
        "mapped_minus_cpp_joule_power_w": float(np.sum(mapped_minus_cpp)),
        "mapped_minus_cpp_joule_power_max_abs_node_w": float(np.max(np.abs(mapped_minus_cpp))),
        "legacy_gradient_joule_power_w": float(np.sum(legacy_joule)),
        "legacy_gradient_minus_cpp_joule_power_w": float(np.sum(legacy_minus_cpp)),
        "electron_diff_minus_terminal_less_joule_w": float(np.sum(requested_closure)),
        "electron_diff_minus_terminal_plus_joule_w": float(np.sum(plus_joule_closure)),
        "electron_diff_minus_terminal_plus_cpp_joule_w": float(np.sum(cpp_plus_joule_closure)),
    }


def _append_restart_metadata(restart_path: Path) -> None:
    with np.load(restart_path, allow_pickle=False) as data:
        state = {key: np.asarray(data[key]) for key in data.files}
    state["SingleTFE/source_snapshot_time_s"] = np.array([SOURCE_SNAPSHOT_TIME_S])
    state["SingleTFE/relative_time_s"] = np.array([state["System/global_time"][0]])
    np.savez_compressed(restart_path, **state)


def run_case(args) -> Dict[str, object]:
    snapshot_path = _resolve_path(args.source_snapshot)
    output_dir = _resolve_path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    mapping = load_center_snapshot(snapshot_path)
    case = _build_single_tfe(mapping)
    _reset_and_map_snapshot(case, mapping)
    tfe = case["tfe"]
    _assert_strict_adiabatic(tfe)

    density = float(Sodium().density(mapping.inlet_plenum_t_k, mapping.inlet_plenum_p_pa))
    inlet_velocity = FIXED_INLET_MASS_FLOW_KG_S / (density * case["flow_area_m2"])
    if not np.isclose(inlet_velocity, FIXED_INLET_VELOCITY_M_S, rtol=0.0, atol=1.0e-12):
        raise AssertionError(
            f"Fixed inlet velocity mismatch: expected {FIXED_INLET_VELOCITY_M_S}, got {inlet_velocity}."
        )

    tec_model = None
    tec_results = None
    if args.mode == "tec":
        tec_model = _build_tec_model(tfe)
        tec_results = _apply_tec_sources(tec_model, tfe)

    history = []
    last_before = None
    last_dt = None
    last_applied_fluid_sources = None
    stop_time = float(args.duration_s)
    while case["system"].global_time < stop_time:
        dt = case["system"].compute_adaptive_dt(
            min_dt=1.0e-6,
            max_dt=float(args.max_dt_s),
            safety_factor=0.8,
        )
        dt = min(dt, stop_time - case["system"].global_time)
        before = _capture_storage_state(case)
        last_before = before
        last_dt = dt
        if tec_model is not None:
            tec_results = _apply_tec_sources(tec_model, tfe)
        case["system"].step(dt, inner_iter=int(args.inner_iter))
        last_applied_fluid_sources = case["system"]._capture_fluid_sources()
        _assert_strict_adiabatic(tfe)
        electrical_output_w = 0.0
        if tec_model is not None:
            global_results = tec_model.get_global_results()
            electrical_output_w = float(global_results["Uout"] * global_results["Iout"])
        history.append(_audit_step(case, before, dt, electrical_output_w))

    if tec_model is not None:
        tec_results = _apply_tec_sources(tec_model, tfe)
    _synchronize_latest_audit(case)
    if last_applied_fluid_sources is None:
        last_applied_fluid_sources = case["system"]._capture_fluid_sources()
    solid_rows, solid_summary = _collect_solid_energy_balance(case, last_before, last_dt)
    interface_rows, interface_summary = _collect_interface_balance(case)
    fluid_rows, fluid_summary = _collect_fluid_volume_balance(
        case,
        last_before,
        last_dt,
        last_applied_fluid_sources,
    )

    restart_path = output_dir / "latest_restart.npz"
    case["system"].save_global_state(str(restart_path))
    _append_restart_metadata(restart_path)
    _write_csv(output_dir / "energy_balance_history.csv", history, fieldnames=ENERGY_HISTORY_FIELDS)
    _write_csv(output_dir / "solid_energy_balance_latest.csv", solid_rows)
    _write_csv(output_dir / "interface_balance_latest.csv", interface_rows)
    _write_csv(output_dir / "fluid_volume_balance_latest.csv", fluid_rows)
    _write_profiles(output_dir, case)
    tec_energy_balance = None
    if tec_results is not None:
        global_results = tec_model.get_global_results()
        tec_energy_balance = _write_tec_nodes(output_dir, tfe, tec_results, global_results)
        tec_energy_balance.update({
            "global_terminal_voltage_v": float(global_results["Uout"]),
            "global_terminal_current_a": float(global_results["Iout"]),
            "global_terminal_power_w": float(global_results["Uout"] * global_results["Iout"]),
        })

    final_audit = history[-1] if history else None
    summary = {
        "case_name": "single_tfe_energy_conservation_v7",
        "mode": args.mode,
        "source_snapshot": os.fspath(snapshot_path),
        "source_snapshot_time_s": SOURCE_SNAPSHOT_TIME_S,
        "relative_time_s": float(case["system"].global_time),
        "n_axial": N_AXIAL,
        "coolant": "Sodium",
        "nuclear_heat_w": FIXED_NUCLEAR_HEAT_W,
        "inlet_velocity_m_s": FIXED_INLET_VELOCITY_M_S,
        "inlet_mass_flow_kg_s": FIXED_INLET_MASS_FLOW_KG_S,
        "outlet_fixed_pressure_pa": FIXED_OUTLET_PRESSURE_PA,
        "tec_target_voltage_v": TEC_TARGET_VOLTAGE_V,
        "tec_wire_resistance_ohm": TEC_WIRE_RESISTANCE_OHM.tolist(),
        "strict_adiabatic_outer_clad": True,
        "outer_clad_boundary_condition_count": len(
            tfe.solids["outer_clad"].boundaries["right"].conditions
        ),
        "outer_clad_heat_loss_w": _outer_clad_heat_loss_w(tfe),
        "history_rows": len(history),
        "inner_iter": int(args.inner_iter),
        "final_audit": final_audit,
        "tail_residual_windows": _tail_residual_windows(history),
        "solid_energy_balance": solid_summary,
        "interface_balance": interface_summary,
        "fluid_volume_balance": fluid_summary,
        "tec_energy_balance": tec_energy_balance,
        "outputs": {
            "restart": os.fspath(restart_path),
            "energy_history": os.fspath(output_dir / "energy_balance_history.csv"),
            "solid_energy_balance": os.fspath(output_dir / "solid_energy_balance_latest.csv"),
            "interface_balance": os.fspath(output_dir / "interface_balance_latest.csv"),
            "fluid_volume_balance": os.fspath(output_dir / "fluid_volume_balance_latest.csv"),
            "fluid_profile": os.fspath(output_dir / "fluid_profile_latest.csv"),
            "solid_profile": os.fspath(output_dir / "solid_temperature_profile_latest.csv"),
            "tec_node_balance": (
                os.fspath(output_dir / "tec_node_balance_latest.csv")
                if tec_results is not None
                else None
            ),
        },
    }
    with (output_dir / "latest_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, ensure_ascii=True)
    return summary


def parse_args():
    parser = argparse.ArgumentParser(description="Single Center-TFE v7 energy-conservation diagnostic.")
    parser.add_argument("--mode", choices=("thermal-baseline", "tec"), default="thermal-baseline")
    parser.add_argument("--duration-s", type=float, default=10.0)
    parser.add_argument("--max-dt-s", type=float, default=0.1)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--source-snapshot", default=SOURCE_SNAPSHOT)
    parser.add_argument(
        "--output-dir",
        default="testModule/single_tfe_energy_conservation_v7",
    )
    args = parser.parse_args()
    if args.duration_s < 0.0:
        parser.error("--duration-s must be non-negative")
    if args.max_dt_s <= 0.0:
        parser.error("--max-dt-s must be positive")
    if args.inner_iter < 1:
        parser.error("--inner-iter must be >= 1")
    return args


if __name__ == "__main__":
    print(json.dumps(run_case(parse_args()), indent=2, ensure_ascii=True))
