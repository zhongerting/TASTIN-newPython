"""Run zero-power V14/V15 cooldowns with shield or direct radiator heating."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Components.ExternalHeatSources import (  # noqa: E402
    OrbitalTableHeatSource,
    W0_8P12_ORBIT_PERIOD_S,
)
from Components.ExternalHeatSources.embedded_flux_tables import load_csv_flux_table_library  # noqa: E402
from Components.RadiatorThermalShield import RadiatorThermalShield  # noqa: E402
from Solvers.Couplers import FluidSolidCouple  # noqa: E402
from testModule.Full_Loop_Cases import (  # noqa: E402
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    V15PipeFinRadiatorConfig,
    build_v14_case_a_system,
    build_v15_case_a_system,
)

INITIAL_TEMPERATURE_K = 350.0
SPACE_TEMPERATURE_K = 4.0
TARGET_FLOW_KG_S = 0.26
STOP_TEMPERATURE_K = 260.0
DIRECT_EXTERNAL_HEAT_SCALE_FACTOR = 1.0
WIRE_RESISTANCE_BASE_OHM = np.asarray([0.001552, 0.001024, 0.000336, 0.000608])
OPERATING_PARAMETERS = {
    "v14": {"pump_head_pa": 8516.44886986068, "emissivity": 0.84, "wire_scale": 1.5},
    "v15": {"pump_head_pa": 40083.288387952, "emissivity": 0.815, "wire_scale": 1.5},
}
SHIELD_PARAMETERS = {
    "background_temperature_k": SPACE_TEMPERATURE_K,
    "shield_view_factor": 0.8,
    "inner_emissivity": 0.8,
    "outer_emissivity": 0.1,
    "conductivity_w_m_k": 0.0008,
    "thickness_m": 0.01,
}


def aggregate_flux_to_six(values: Iterable[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size not in (18, 78):
        raise ValueError("External heat must contain 18 or 78 circumferential values.")
    return values.reshape(6, values.size // 6).mean(axis=1)


def expand_six_flux(values: Iterable[float], unit_count: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.shape != (6,) or int(unit_count) % 6:
        raise ValueError("Six shield sectors must map evenly to radiator units.")
    return np.repeat(values, int(unit_count) // 6)


def _set_uniform_initial_state(build: Dict[str, Any], temperature_k: float) -> None:
    system = build["system"]
    net = system.fluid_solver
    for volume in net.volumes_obj:
        volume.T = float(temperature_k)
        volume.h = float(volume.material.enthalpy(volume.T, volume.P))
        volume.update_properties(volume.material)
    net._initialize_state_from_objects()
    net._update_fluid_properties()
    net._sync_vectors_to_objects()
    for solid in system.solid_components.values():
        solid.T[...] = float(temperature_k)
        if hasattr(solid, "dTdt"):
            solid.dTdt[...] = 0.0
        solid.current_time = float(system.global_time)
        solid._update_properties()
        solid._update_boundaries_state(current_time=float(system.global_time))
        solid.set_ode_method("implicit_euler")
    for unit in _radiator_units(build):
        if hasattr(unit, "last_fin_temperature"):
            unit.last_fin_temperature[...] = float(temperature_k)
        if hasattr(unit, "last_fin_effective_temperature_distribution"):
            unit.last_fin_effective_temperature_distribution[...] = float(temperature_k)


def _radiator_units(build: Dict[str, Any]) -> list[Any]:
    if "radiator_units" in build:
        return list(build["radiator_units"])
    return [hp for ring in build["ring_hps"] for hp in ring.hp_units]


def _set_v14_area_multipliers(build: Dict[str, Any]) -> None:
    symmetric = float(build["radiator_config"].symmetric_ring_multiplier)
    for ring in build["ring_hps"]:
        for _, hp, multiplier in ring._iter_present_hp_units_with_multiplier():
            hp.radiation_area_multiplier = symmetric * float(multiplier)


def _build_case(
    case: str,
    external_heat_target: str = "shield",
    initial_temperature_k: float = INITIAL_TEMPERATURE_K,
    pump_head_pa: Optional[float] = None,
    emissivity: Optional[float] = None,
) -> Dict[str, Any]:
    if external_heat_target not in ("shield", "radiator"):
        raise ValueError("external_heat_target must be 'shield' or 'radiator'.")
    direct_external_heat = external_heat_target == "radiator"
    parameters = OPERATING_PARAMETERS[case]
    pump_head_pa = parameters["pump_head_pa"] if pump_head_pa is None else float(pump_head_pa)
    emissivity = parameters["emissivity"] if emissivity is None else float(emissivity)
    core = FullLoopCoreConfig(
        inlet_temperature_k=float(initial_temperature_k),
        main_tec_enabled=False,
    )
    flow = FullLoopFlowConfig(total_flow_kg_s=TARGET_FLOW_KG_S)
    pump = FullLoopPumpConfig(
        pump_total_head_pa=float(pump_head_pa),
        pump_flow_control=True,
        target_flow_kg_s=TARGET_FLOW_KG_S,
    )
    if case == "v14":
        build = build_v14_case_a_system(
            core_config=core,
            flow_config=flow,
            pump_config=pump,
            radiator_config=V14HeatPipeRadiatorConfig(
                t_space_k=SPACE_TEMPERATURE_K,
                hp_initial_temp_k=float(initial_temperature_k),
                hp_emissivity=float(emissivity),
                fin_emissivity=float(emissivity),
                external_heat_enabled=direct_external_heat,
                external_heat_scale_factor=DIRECT_EXTERNAL_HEAT_SCALE_FACTOR,
            ),
        )
        _set_v14_area_multipliers(build)
        csv_path = REPO_ROOT / "Components" / "ExternalHeatSources" / "is58p5_w0_8p12_N18_sum.csv"
    else:
        build = build_v15_case_a_system(
            core_config=core,
            flow_config=flow,
            pump_config=pump,
            radiator_config=V15PipeFinRadiatorConfig(
                t_space_k=SPACE_TEMPERATURE_K,
                tube_emissivity=float(emissivity),
                fin_emissivity=float(emissivity),
                external_heat_enabled=direct_external_heat,
                external_heat_scale_factor=DIRECT_EXTERNAL_HEAT_SCALE_FACTOR,
                solid_ode_method="implicit_euler",
            ),
        )
        csv_path = REPO_ROOT / "Components" / "ExternalHeatSources" / "is58p5_w0_8p12_N78_sum.csv"

    build["core"].point_reactor = None
    build["core"].enable_tec_coupled = False
    build["core"].update_neutronic_power(p_total=0.0, p_fiss=0.0, p_decay=0.0, alpha=1.0)
    _set_uniform_initial_state(build, float(initial_temperature_k))

    units = _radiator_units(build)
    library = load_csv_flux_table_library(str(csv_path), W0_8P12_ORBIT_PERIOD_S)
    source_scale_factor = (
        SHIELD_PARAMETERS["outer_emissivity"]
        if external_heat_target == "shield"
        else DIRECT_EXTERNAL_HEAT_SCALE_FACTOR
    )
    source = OrbitalTableHeatSource(
        shape=(len(library.available_ids()),),
        table_ids=library.available_ids(),
        table_library=library,
        scale_factor=source_scale_factor,
        periodic=True,
    )
    first_table = library.get_table(library.available_ids()[0])
    build["diagnostic_external_heat_source"] = source
    build["external_heat_csv"] = str(csv_path)
    build["external_heat_period_s"] = float(first_table.time[-1])
    build["external_heat_target"] = external_heat_target
    build["prestart_pump_head_pa"] = float(pump_head_pa)
    build["prestart_emissivity"] = float(emissivity)
    if external_heat_target == "shield":
        shield = RadiatorThermalShield(
            name=f"{case.upper()}_RadiatorThermalShield",
            radiator_units=units,
            model=RadiatorThermalShield.MODEL_SEGMENT_BALANCE,
            solar_heat_flux_w_m2=np.zeros(len(units)),
            **SHIELD_PARAMETERS,
        )
        build["system"].components.insert(0, shield)
        build["radiator_thermal_shield"] = shield
        build["shield_external_heat_source"] = source
        build["shield_external_heat_csv"] = str(csv_path)
        build["shield_sector_count"] = 6
    return build


def _set_local_implicit_coupling(system) -> int:
    count = 0
    for coupler in system.couplers:
        if isinstance(coupler, FluidSolidCouple) and coupler.solid_node_capacitance is not None:
            coupler.set_coupling_time_scheme("local_implicit")
            count += 1
    return count


def _update_external_heat(build: Dict[str, Any]) -> np.ndarray:
    source = build["diagnostic_external_heat_source"]
    six_flux = aggregate_flux_to_six(source.get_heat_flux(build["system"].global_time))
    if build["external_heat_target"] == "shield":
        shield = build["radiator_thermal_shield"]
        shield.solar_heat_flux_w_m2 = expand_six_flux(six_flux, len(shield.radiator_units))
    return six_flux


def _update_shield_heat(build: Dict[str, Any]) -> np.ndarray:
    if build["external_heat_target"] != "shield":
        raise ValueError("Shield heat update requires external_heat_target='shield'.")
    return _update_external_heat(build)


def _radiator_rejection_w(build: Dict[str, Any]) -> float:
    if "ring_hps" in build:
        return _v14_heatpipe_rejection_w(build) + _v14_ring_wall_rejection_w(build)
    return float(sum(
        np.sum(unit.get_heat_exchange_breakdown()["gross_rejection"])
        for unit in build["radiator_units"]
    ))


def _v14_heatpipe_rejection_w(build: Dict[str, Any]) -> float:
    symmetric = float(build["radiator_config"].symmetric_ring_multiplier)
    return symmetric * float(sum(
        ring.get_total_heat_rejection_scaled()
        for ring in build["ring_hps"]
    ))


def _v14_ring_wall_rejection_w(build: Dict[str, Any]) -> float:
    symmetric = float(build["radiator_config"].symmetric_ring_multiplier)
    single_ring = sum(
        max(0.0, -float(np.sum(solid.boundaries["right"].current_flux)))
        for solid in build["ring_solids"]
    )
    return symmetric * float(single_ring)


def _direct_external_heat_w(build: Dict[str, Any]) -> float:
    if build["external_heat_target"] != "radiator":
        return 0.0
    if "ring_hps" in build:
        symmetric = float(build["radiator_config"].symmetric_ring_multiplier)
        return symmetric * float(sum(
            ring.get_total_external_heat_absorption_scaled(build["system"].global_time)
            for ring in build["ring_hps"]
        ))
    return float(sum(
        np.sum(
            unit.get_external_heat_absorption_distribution(
                build["system"].global_time
            )[2]
        )
        for unit in build["radiator_units"]
    ))


def collect_metrics(build: Dict[str, Any], six_flux: np.ndarray) -> Dict[str, Any]:
    system = build["system"]
    net = system.fluid_solver
    solids = [np.asarray(solid.T, dtype=float) for solid in system.solid_components.values()]
    metrics = {
        "time_s": float(system.global_time),
        "min_fluid_T_K": float(np.min(net.T_vec)),
        "max_fluid_T_K": float(np.max(net.T_vec)),
        "min_solid_T_K": float(min(np.min(value) for value in solids)),
        "max_solid_T_K": float(max(np.max(value) for value in solids)),
        "core_inlet_T_K": float(build["core_inlet_connector"].T),
        "core_outlet_T_K": float(build["core_outlet_connector"].T),
        "pump_flow_kg_s": float(build["pump_a"].W),
        "pressure_drop_Pa": float(np.max(net.P_vec) - np.min(net.P_vec)),
        "radiator_rejection_W": _radiator_rejection_w(build),
        "radiator_external_heat_W": _direct_external_heat_w(build),
    }
    if "ring_hps" in build:
        metrics.update({
            "radiator_heatpipe_rejection_W": _v14_heatpipe_rejection_w(build),
            "radiator_ring_wall_rejection_W": _v14_ring_wall_rejection_w(build),
        })
    if build["external_heat_target"] == "shield":
        shield = build["radiator_thermal_shield"]
        metrics.update({
            "shield_external_heat_W": float(shield.last_q_solar_w),
            "shield_to_space_W": float(shield.last_q_to_space_w),
            "shield_from_radiator_W": float(shield.last_q_from_radiator_w),
            "shield_effective_background_K": float(shield.last_effective_background_mean_k),
            "shield_inner_temperature_K": float(shield.last_inner_temperature_mean_k),
            "shield_outer_temperature_K": float(shield.last_outer_temperature_mean_k),
            **{f"shield_sector_{index + 1}_W_m2": float(value) for index, value in enumerate(six_flux)},
        })
    else:
        metrics.update({
            f"external_heat_sector_{index + 1}_W_m2": float(value)
            for index, value in enumerate(six_flux)
        })
    return metrics


class StateChunkWriter:
    def __init__(
        self,
        output_dir: Path,
        system: Any,
        samples_per_chunk: int,
        radiator_units: Iterable[Any],
    ):
        self.output_dir = output_dir / "state_history"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.system = system
        self.radiator_units = list(radiator_units)
        self.samples_per_chunk = max(1, int(samples_per_chunk))
        self.samples: list[Dict[str, Any]] = []
        self.chunk_index = 0

    def append(self, metrics: Dict[str, Any]) -> None:
        net = self.system.fluid_solver
        fin_background = []
        for unit in self.radiator_units:
            background = getattr(unit, "radiation_background_temperature", None)
            if background is None:
                background = np.full(
                    np.asarray(unit.last_fin_radiation_distribution).shape,
                    float(unit.T_space),
                )
            fin_background.append(np.asarray(background, dtype=float))
        self.samples.append({
            "metrics": dict(metrics),
            "fluid_P": np.array(net.P_vec, copy=True),
            "fluid_T": np.array(net.T_vec, copy=True),
            "fluid_h": np.array(net.h_vec, copy=True),
            "fluid_rho": np.array(net.rho_vec, copy=True),
            "junction_W": np.array(net.W_vec, copy=True),
            "fluid_Q_wall": np.asarray([volume.Q_wall for volume in net.volumes_obj], dtype=float),
            "fluid_Q_vol": np.asarray([volume.Q_vol for volume in net.volumes_obj], dtype=float),
            "solids": [np.array(solid.T, copy=True) for solid in self.system.solid_components.values()],
            "fin_temperature": np.stack([
                np.asarray(unit.last_fin_temperature, dtype=float)
                for unit in self.radiator_units
            ]),
            "fin_radiation": np.stack([
                np.asarray(unit.last_fin_radiation_distribution, dtype=float)
                for unit in self.radiator_units
            ]),
            "fin_absorption": np.stack([
                np.asarray(unit.last_fin_absorption_distribution, dtype=float)
                for unit in self.radiator_units
            ]),
            "fin_net_from_root": np.stack([
                np.asarray(unit.last_fin_net_from_root_distribution, dtype=float)
                for unit in self.radiator_units
            ]),
            "radiation_background": np.stack(fin_background),
        })
        if len(self.samples) >= self.samples_per_chunk:
            self.flush()

    def flush(self) -> None:
        if not self.samples:
            return
        net = self.system.fluid_solver
        solids = list(self.system.solid_components.items())
        payload: Dict[str, Any] = {
            "volume_names": np.asarray([volume.name for volume in net.volumes_obj]),
            "junction_names": np.asarray([junction.name for junction in net.junctions_obj]),
            "solid_names": np.asarray([name for name, _ in solids]),
            "radiator_unit_names": np.asarray([
                unit.name for unit in self.radiator_units
            ]),
            "metrics_json": np.asarray([json.dumps(sample["metrics"], sort_keys=True) for sample in self.samples]),
        }
        for key in ("fluid_P", "fluid_T", "fluid_h", "fluid_rho", "junction_W", "fluid_Q_wall", "fluid_Q_vol"):
            payload[key] = np.stack([sample[key] for sample in self.samples])
        for index, _ in enumerate(solids):
            payload[f"solid_{index:03d}_T"] = np.stack([sample["solids"][index] for sample in self.samples])
        for key in (
            "fin_temperature",
            "fin_radiation",
            "fin_absorption",
            "fin_net_from_root",
            "radiation_background",
        ):
            payload[key] = np.stack([sample[key] for sample in self.samples])
        start = int(round(self.samples[0]["metrics"]["time_s"]))
        end = int(round(self.samples[-1]["metrics"]["time_s"]))
        path = self.output_dir / f"state_{self.chunk_index:04d}_t{start:05d}_to_t{end:05d}.npz"
        np.savez_compressed(path, **payload)
        self.chunk_index += 1
        self.samples.clear()


def run(
    case: str,
    output_dir: Path,
    duration_s: float,
    max_dt_s: float,
    record_interval_s: float = 1.0,
    restart_interval_s: float = 60.0,
    restart_in: Optional[Path] = None,
    external_heat_target: str = "shield",
    stop_temperature_k: float = STOP_TEMPERATURE_K,
    initial_temperature_k: float = INITIAL_TEMPERATURE_K,
    pump_head_pa: Optional[float] = None,
    emissivity: Optional[float] = None,
    wire_scale: Optional[float] = None,
    record_full_state: bool = True,
) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    build = _build_case(
        case,
        external_heat_target=external_heat_target,
        initial_temperature_k=initial_temperature_k,
        pump_head_pa=pump_head_pa,
        emissivity=emissivity,
    )
    system = build["system"]
    if restart_in is None:
        system.initialize_system(dt_init=0.01, tol=1.0e-5, max_iter=1000)
    else:
        system.load_global_state(str(restart_in))
        for solid in system.solid_components.values():
            solid.set_ode_method("implicit_euler")
    _set_local_implicit_coupling(system)
    start_time = float(system.global_time)
    target_time = start_time + float(duration_s)
    history_path = output_dir / "history.csv"
    latest_restart = output_dir / "latest_restart.npz"
    writer = (
        StateChunkWriter(
            output_dir,
            system,
            round(60.0 / record_interval_s),
            _radiator_units(build),
        )
        if record_full_state
        else None
    )
    next_record = start_time
    next_restart = start_time + float(restart_interval_s)
    fields = None
    stop_reason = "duration_limit"
    latest: Dict[str, Any] = {}

    while system.global_time < target_time - 1.0e-10:
        six_flux = _update_external_heat(build)
        build["core"].update_neutronic_power(p_total=0.0, p_fiss=0.0, p_decay=0.0, alpha=1.0)
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=max_dt_s,
            safety_factor=0.8,
            respect_fluid_cfl=False,
        )
        dt = min(dt, max_dt_s, next_record - system.global_time if next_record > system.global_time else max_dt_s, target_time - system.global_time)
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=300)

        if system.global_time >= next_record - 1.0e-10:
            six_flux = _update_external_heat(build)
            latest = collect_metrics(build, six_flux)
            if not all(math.isfinite(float(value)) for value in latest.values()):
                raise FloatingPointError(f"Non-finite state at t={system.global_time}: {latest}")
            write_header = fields is None and not history_path.exists()
            fields = list(latest)
            with history_path.open("a", newline="", encoding="utf-8") as handle:
                csv_writer = csv.DictWriter(handle, fieldnames=fields)
                if write_header:
                    csv_writer.writeheader()
                csv_writer.writerow(latest)
            if writer is not None:
                writer.append(latest)
            (output_dir / "latest_state.json").write_text(
                json.dumps(latest, indent=2, sort_keys=True), encoding="utf-8"
            )
            print(json.dumps(latest, sort_keys=True), flush=True)
            next_record += float(record_interval_s)

        if system.global_time >= next_restart - 1.0e-10:
            checkpoint = output_dir / f"restart_t{int(round(system.global_time)):05d}s.npz"
            system.save_global_state(str(checkpoint))
            system.save_global_state(str(latest_restart))
            next_restart += float(restart_interval_s)

        if float(np.min(system.fluid_solver.T_vec)) <= float(stop_temperature_k):
            stop_reason = "minimum_coolant_temperature"
            break

    if not latest or latest.get("time_s") != float(system.global_time):
        six_flux = _update_external_heat(build)
        latest = collect_metrics(build, six_flux)
        if writer is not None:
            writer.append(latest)
    if writer is not None:
        writer.flush()
    system.save_global_state(str(latest_restart))
    summary = {
        "case": case,
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "stop_reason": stop_reason,
        "stop_temperature_k": float(stop_temperature_k),
        "target_flow_kg_s": TARGET_FLOW_KG_S,
        "pump_total_head_pa": float(build["prestart_pump_head_pa"]),
        "radiator_emissivity": float(build["prestart_emissivity"]),
        "wire_resistance_scale": float(
            OPERATING_PARAMETERS[case]["wire_scale"] if wire_scale is None else wire_scale
        ),
        "initial_temperature_k": float(initial_temperature_k),
        "space_temperature_k": SPACE_TEMPERATURE_K,
        "core_power_w": 0.0,
        "tec_enabled": False,
        "full_state_history_enabled": bool(record_full_state),
        "external_heat_target": external_heat_target,
        "external_heat_csv": build["external_heat_csv"],
        "external_heat_period_s": build["external_heat_period_s"],
        "external_heat_source_scale_factor": float(
            build["diagnostic_external_heat_source"].scale_factor
        ),
        "latest_metrics": latest,
        "latest_restart": str(latest_restart),
    }
    summary["wire_resistance_ohm"] = (
        WIRE_RESISTANCE_BASE_OHM * summary["wire_resistance_scale"]
    ).tolist()
    if external_heat_target == "shield":
        summary["shield_parameters"] = SHIELD_PARAMETERS
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("v14", "v15"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=10000.0)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--record-interval", type=float, default=1.0)
    parser.add_argument("--restart-interval", type=float, default=60.0)
    parser.add_argument("--restart-in", type=Path)
    parser.add_argument("--external-heat-target", choices=("shield", "radiator"), default="shield")
    parser.add_argument("--stop-temperature", type=float, default=STOP_TEMPERATURE_K)
    parser.add_argument("--initial-temperature", type=float, default=INITIAL_TEMPERATURE_K)
    parser.add_argument("--pump-head", type=float)
    parser.add_argument("--emissivity", type=float)
    parser.add_argument("--wire-scale", type=float)
    parser.add_argument(
        "--record-full-state",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    args = parser.parse_args()
    run(
        args.case,
        args.output_dir,
        args.duration,
        args.max_dt,
        args.record_interval,
        args.restart_interval,
        args.restart_in,
        args.external_heat_target,
        args.stop_temperature,
        args.initial_temperature,
        args.pump_head,
        args.emissivity,
        args.wire_scale,
        args.record_full_state,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
