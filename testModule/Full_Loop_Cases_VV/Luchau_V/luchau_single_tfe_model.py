from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from Components.BaseComponent import BaseComponent
from Components.TFEUnit import GapConfig, TFEGeometry, TFEMeshParams, TFEUnit
from Components.tec_electric import electric_field_from_node_potential
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
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
from Solvers.SystemManager import SystemManager


# 定义Luchau单TFFE模型的参数
LUCHAU_AXIAL_LENGTH_ALLOCATION_M = (0.065, 0.377, 0.065)
LUCHAU_AXIAL_NODE_ALLOCATION = (6, 25, 6)
LUCHAU_TOTAL_FLOW_KG_S = 1.3
LUCHAU_PHYSICAL_TFE_COUNT = 37

# 定义Luchau单TFFE模型的配置参数
@dataclass(frozen=True)
class LuchauSingleTFEConfig:
    thermal_power_w: float | None
    target_voltage_v: float | None
    inlet_temperature_k: float = 727.0
    reference_pressure_pa: float = 205000.0
    outlet_pressure_pa: float = 205000.0
    single_tfe_flow_kg_s: float = LUCHAU_TOTAL_FLOW_KG_S / LUCHAU_PHYSICAL_TFE_COUNT
    heater_length_m: float = 0.30
    cesium_pressure_torr: float = 0.4
    tec_gap_h_eq_w_m2_k: float = 29.0
    i_guess_a: float = 150.0
    wire_resistance_ohm: float = 0.0

    def __post_init__(self) -> None:
        required_positive = {
            "thermal_power_w": self.thermal_power_w,
            "target_voltage_v": self.target_voltage_v,
            "single_tfe_flow_kg_s": self.single_tfe_flow_kg_s,
            "heater_length_m": self.heater_length_m,
            "cesium_pressure_torr": self.cesium_pressure_torr,
        }
        for name, value in required_positive.items():
            if value is None:
                raise ValueError(f"{name} must be provided.")
            if float(value) <= 0.0:
                raise ValueError(f"{name} must be positive.")
        if float(self.heater_length_m) > sum(LUCHAU_AXIAL_LENGTH_ALLOCATION_M):
            raise ValueError("heater_length_m cannot exceed the full TFE length.")


def build_node_lengths() -> np.ndarray:
    lower, active, upper = LUCHAU_AXIAL_LENGTH_ALLOCATION_M
    n_lower, n_active, n_upper = LUCHAU_AXIAL_NODE_ALLOCATION
    return np.array(
        [lower / n_lower] * n_lower
        + [active / n_active] * n_active
        + [upper / n_upper] * n_upper,
        dtype=float,
    )


def build_axial_faces() -> np.ndarray:
    return np.insert(np.cumsum(build_node_lengths()), 0, 0.0)


def build_center_heater_profile(config: LuchauSingleTFEConfig) -> np.ndarray:
    faces = build_axial_faces()
    full_length = float(faces[-1])
    heater_length = float(config.heater_length_m)
    heater_start = 0.5 * (full_length - heater_length)
    heater_end = heater_start + heater_length
    overlap = np.maximum(
        0.0,
        np.minimum(faces[1:], heater_end) - np.maximum(faces[:-1], heater_start),
    )
    total = float(np.sum(overlap))
    if total <= 0.0:
        raise ValueError("center heater profile has zero heated length.")
    return overlap / total


def build_luchau_geometry() -> TFEGeometry:
    return TFEGeometry(
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
        height=sum(LUCHAU_AXIAL_LENGTH_ALLOCATION_M),
    )


def build_luchau_mesh_params() -> TFEMeshParams:
    return TFEMeshParams(
        n_axial=sum(LUCHAU_AXIAL_NODE_ALLOCATION),
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=1,
        n_r_outer_clad=1,
        n_r_moderator=3,
    )


def build_luchau_materials(coolant: SodiumPotassium78) -> dict[str, Any]:
    return {
        "UO2": UO2(),
        "MoNb": MoNb(),
        "Molybdenum": Molybdenum(),
        "StainlessSteel": AusteniticStainlessSteel(),
        "ZrH": ZirconiumHydride(),
        "BerylliumOxide": BerylliumOxide(),
        "Sodium": coolant,
    }


def build_luchau_single_tfe(config: LuchauSingleTFEConfig) -> dict[str, Any]:
    geometry = build_luchau_geometry()
    mesh = build_luchau_mesh_params()
    node_lengths = build_node_lengths()
    axial_faces = build_axial_faces()
    heater_profile = build_center_heater_profile(config)
    coolant = SodiumPotassium78()

    flow_area = math.pi * (geometry.r_coolant_outer**2 - geometry.r_coolant_inner**2)
    hydraulic_diam = 2.0 * (geometry.r_coolant_outer - geometry.r_coolant_inner)

    inlet = IncompressibleBoundaryVolume(
        name="Luchau_Inlet",
        material=coolant,
        P=float(config.reference_pressure_pa),
        T=float(config.inlet_temperature_k),
        flow_area=flow_area,
        hydraulic_diam=hydraulic_diam,
    )
    outlet = IncompressibleBoundaryVolume(
        name="Luchau_Outlet",
        material=coolant,
        P=float(config.outlet_pressure_pa),
        T=float(config.inlet_temperature_k),
        flow_area=flow_area,
        hydraulic_diam=hydraulic_diam,
    )
    outlet.is_pressure_boundary = True

    channel = NonUniformIncompressibleFluidChannel(
        name="Luchau_Channel",
        node_lengths=node_lengths,
        flow_area=flow_area,
        hydraulic_diam=hydraulic_diam,
        initial_P=float(config.reference_pressure_pa),
        initial_T=float(config.inlet_temperature_k),
        material=coolant,
    )

    inlet_junction = InletJunction(
        name="Luchau_InletFlow",
        from_vol=inlet,
        to_vol=channel.volumes[0],
        W_initial=float(config.single_tfe_flow_kg_s),
    )
    outlet_junction = FlowJunction(
        name="Luchau_OutletPressure",
        from_vol=channel.volumes[-1],
        to_vol=outlet,
        flow_area=flow_area,
    )
    outlet_junction.W = float(config.single_tfe_flow_kg_s)

    network = HydraulicNetwork(
        volumes=[inlet, *channel.volumes, outlet],
        junctions=[inlet_junction, *channel.internal_junctions, outlet_junction],
        gravity_vector=0.0,
    )

    tfe = TFEUnit(
        name="Luchau_SingleTFE",
        geometry=geometry,
        mesh_params=mesh,
        materials=build_luchau_materials(coolant),
        coolant_channel=channel,
        fission_gas_config=GapConfig("simplified", 5678.0, Xenon(), 0.15, 0.15),
        tec_gap_config=GapConfig("simplified", float(config.tec_gap_h_eq_w_m2_k), Cesium(), 0.15, 0.60),
        he_gap_config=GapConfig("simplified", 5678.0, Helium(), 0.60, 0.80),
        co2_gap_config=GapConfig("simplified", 53.6, CarbonDioxide(), 0.80, 0.80),
        power_fraction=1.0,
        axial_power_profile=heater_profile,
        axial_length_allocation=list(LUCHAU_AXIAL_LENGTH_ALLOCATION_M),
        axial_node_allocation=list(LUCHAU_AXIAL_NODE_ALLOCATION),
        axial_contact_resistance=0.0,
        strict_adiabatic_single_tfe=True,
    )
    _reset_solid_temperatures(tfe, float(config.inlet_temperature_k))

    system = SystemManager(fluid_network=network, start_time=0.0)
    system.add_component(tfe)
    system.initialize_system()
    tfe.update_neutronic_power(float(config.thermal_power_w), alpha=1.0)

    return {
        "system": system,
        "network": network,
        "tfe": tfe,
        "channel": channel,
        "inlet": inlet,
        "outlet": outlet,
        "inlet_junction": inlet_junction,
        "outlet_junction": outlet_junction,
        "node_lengths_m": node_lengths,
        "axial_faces_m": axial_faces,
        "heater_profile": heater_profile,
        "flow_area_m2": flow_area,
        "hydraulic_diam_m": hydraulic_diam,
    }


def _reset_solid_temperatures(tfe: TFEUnit, temperature_k: float) -> None:
    for solid in tfe.get_solids():
        solid.T[:] = float(temperature_k)
        if hasattr(solid, "initialize_state"):
            solid.initialize_state()


def cesium_pressure_from_tcs(tcs_k: float) -> float:
    if float(tcs_k) <= 0.0:
        raise ValueError("tcs_k must be positive.")
    tcs = float(tcs_k)
    return 2.45e8 / math.sqrt(tcs) * math.exp(-8910.0 / tcs)


def tcs_from_cesium_pressure(pcs_torr: float) -> float:
    if float(pcs_torr) <= 0.0:
        raise ValueError("pcs_torr must be positive.")
    target = float(pcs_torr)
    lo = 300.0
    hi = 1200.0
    while cesium_pressure_from_tcs(lo) > target:
        lo *= 0.8
    while cesium_pressure_from_tcs(hi) < target:
        hi *= 1.2
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if cesium_pressure_from_tcs(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def configure_luchau_thermocalc(thermo_model: Any, build: dict[str, Any], config: LuchauSingleTFEConfig) -> None:
    tfe: TFEUnit = build["tfe"]
    geom = tfe.geom
    node_lengths = np.asarray(build["node_lengths_m"], dtype=float)
    n_nodes = int(tfe.mesh.n_axial)
    shape = (1, n_nodes)

    emitter_temperature = _axial_mean_temperature(tfe.solids["emitter"], tfe.mesh.n_r_emitter, n_nodes).reshape(shape)
    collector_temperature = _axial_mean_temperature(tfe.solids["collector"], tfe.mesh.n_r_collector, n_nodes).reshape(shape)

    input_data = thermo_model._input_data
    input_data.dlE = node_lengths.reshape(shape)
    input_data.dlC = node_lengths.reshape(shape)
    input_data.sideAreaE = (2.0 * math.pi * geom.r_emitter_outer * node_lengths).reshape(shape)
    input_data.sideAreaC = (2.0 * math.pi * geom.r_collector_inner * node_lengths).reshape(shape)
    input_data.crossAreaE = np.array([math.pi * (geom.r_emitter_outer**2 - geom.r_fission_gas_outer**2)], dtype=float)
    input_data.crossAreaC = np.array([math.pi * (geom.r_collector_outer**2 - geom.r_collector_inner**2)], dtype=float)
    input_data.d_gap = np.array([(geom.r_collector_inner - geom.r_emitter_outer) * 1000.0], dtype=float)
    input_data.resistanceWire = np.full((1, 4), float(config.wire_resistance_ohm), dtype=float)
    input_data.wireU = np.array([[float(config.target_voltage_v), float(config.target_voltage_v), 0.0, 0.0]], dtype=float)

    tcs = np.full(shape, tcs_from_cesium_pressure(float(config.cesium_pressure_torr)), dtype=float)
    thermo_model.setup_circuit_mode("fixed_u", float(config.target_voltage_v), I_guess=float(config.i_guess_a))
    thermo_model.set_temperatures(emitter_temperature, collector_temperature)
    thermo_model.set_tcs(tcs)


class LuchauSingleTFEThermoCalcCoupler(BaseComponent):
    def __init__(self, name: str, tfe: TFEUnit, thermo_model: Any, alpha_tec: float = 1.0):
        super().__init__(name)
        self.tfe = tfe
        self.thermo_model = thermo_model
        self.alpha_tec = float(alpha_tec)
        self.last_global_results: dict[str, Any] | None = None

    def pre_step(self, dt: float, current_time: float):
        self.sync_thermo_electric()

    def sync_thermo_electric(self) -> None:
        tfe = self.tfe
        n_nodes = int(tfe.mesh.n_axial)
        emitter_t = np.asarray(tfe.solids["emitter"].boundaries["right"].T_surface, dtype=float).reshape(1, n_nodes)
        collector_t = np.asarray(tfe.solids["collector"].boundaries["left"].T_surface, dtype=float).reshape(1, n_nodes)
        self.thermo_model.set_temperatures(emitter_t, collector_t)
        self.thermo_model.calculate(verbose=False)
        self.last_global_results = self.thermo_model.get_global_results()

        res = self.thermo_model.get_tec_results(0)
        if res is None:
            tfe.clear_tec_sources()
            return

        zero = bool((self.last_global_results or {}).get("zero_emission_skipped", False))
        if zero and not all(key in res for key in ("joulePowerE", "joulePowerC")):
            tfe.clear_tec_sources()
            return

        ue_abs = np.asarray(res.get("UE", np.zeros(n_nodes)), dtype=float)
        uc_abs = np.asarray(res.get("UC", np.zeros(n_nodes)), dtype=float)
        rho_e = np.asarray(res.get("rhoE", np.ones(n_nodes) * 1.0e-6), dtype=float)
        rho_c = np.asarray(res.get("rhoC", np.ones(n_nodes) * 1.0e-6), dtype=float)
        joule_power_e = np.asarray(res.get("joulePowerE", np.zeros(n_nodes)), dtype=float)
        joule_power_c = np.asarray(res.get("joulePowerC", np.zeros(n_nodes)), dtype=float)

        if hasattr(tfe, "common_y_faces"):
            y_faces = np.asarray(tfe.common_y_faces, dtype=float)
            e_emit = electric_field_from_node_potential(ue_abs, y_faces=y_faces)
            e_coll = electric_field_from_node_potential(uc_abs, y_faces=y_faces)
        else:
            input_data = self.thermo_model._input_data
            e_emit = electric_field_from_node_potential(ue_abs, node_lengths=np.asarray(input_data.dlE[0], dtype=float))
            e_coll = electric_field_from_node_potential(uc_abs, node_lengths=np.asarray(input_data.dlC[0], dtype=float))

        j_density = np.asarray(res.get("J", np.zeros(n_nodes)), dtype=float) * 1.0e4
        phi_e = np.asarray(res.get("phiE", np.zeros(n_nodes)), dtype=float)
        te = np.asarray(res.get("TE", emitter_t[0]), dtype=float)
        q_e_flux = -1.0 * j_density * (phi_e + 2.0 * 8.617e-5 * te)
        q_c_flux = 1.0 * j_density * (phi_e + 2.0 * 8.617e-5 * te - (ue_abs - uc_abs))

        tfe.update_electric_field_diagnostics(
            E_emit=e_emit,
            rho_emit=rho_e,
            E_coll=e_coll,
            rho_coll=rho_c,
        )
        tfe.update_joule_power_sources(
            Q_emitter_axial=joule_power_e,
            Q_collector_axial=joule_power_c,
            alpha=self.alpha_tec,
        )
        tfe.update_plasma_flux(
            q_e_flux=q_e_flux,
            q_c_flux=q_c_flux,
            alpha=self.alpha_tec,
        )
def _axial_mean_temperature(solid: Any, n_radial: int, n_axial: int) -> np.ndarray:
    values = np.asarray(solid.T, dtype=float).reshape(int(n_radial), int(n_axial))
    return np.mean(values, axis=0)
