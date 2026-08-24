"""Run the three sequential V14 cold-start stages."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.common_config import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.run_v14_shield_radiator_startup import (
    _set_uniform_temperature,
    append_v14_system_history,
    capture_v14_history_reference,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.v14_case import (
    build_v14_case_a_system,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.v14_heatpipe_radiator import (
    V14HeatPipeRadiatorConfig,
)
from testModule.run_v8_caseA_common import WIRE_RESISTANCE_OHM
from testModule.v13_startup_control import apply_tec_gap_h_eq


CASE_DIR = Path(__file__).resolve().parent
INITIAL_TEMPERATURE_K = 300.0
INITIAL_POWER_W = 21.0
STAGE1_TARGET_POWER_W = 60_000.0
FINAL_POWER_W = 210_000.0
STAGE1_RAMP_W_S = 600.0
STAGE2_RAMP_W_S = 80.0
STAGE1_END_S = (STAGE1_TARGET_POWER_W - INITIAL_POWER_W) / STAGE1_RAMP_W_S
FINAL_POWER_TIME_S = STAGE1_END_S + (
    FINAL_POWER_W - STAGE1_TARGET_POWER_W
) / STAGE2_RAMP_W_S
CESIUM_TEC_TIME_S = 2800.0
END_TIME_S = 10_000.0
FLOW_KG_S = 2.46
HELIUM_GAP_H_W_M2_K = 5678.0
CESIUM_GAP_H_W_M2_K = 29.0
LOAD_RESISTANCE_OHM = 0.003
CURRENT_LIMIT_A = 216.0
WIRE_RESISTANCE_SCALE = 0.335
SPACE_TEMPERATURE_K = 4.0
RADIATOR_EMISSIVITY = 0.7475
LOWER_HP_DOWN_VIEW_FACTOR = 0.4
EXTERNAL_HEAT_PERIOD_S = 5668.144369
TEC_LOOKUP_DB = "ThermoCalc/emission_runtime_db_v2/pcs_0p02_5torr"

STAGES = {
    1: (
        0.0,
        STAGE1_END_S,
        "stage_1_ramp_21w_to_60kw",
    ),
    2: (
        STAGE1_END_S,
        CESIUM_TEC_TIME_S,
        "stage_2_ramp_60kw_to_210kw_hold_2800s",
    ),
    3: (
        CESIUM_TEC_TIME_S,
        END_TIME_S,
        "stage_3_cesium_tec_to_10000s",
    ),
}


def prescribed_power_w(time_s: float) -> float:
    time_s = float(time_s)
    if time_s <= STAGE1_END_S:
        return min(
            STAGE1_TARGET_POWER_W,
            INITIAL_POWER_W + STAGE1_RAMP_W_S * time_s,
        )
    if time_s <= FINAL_POWER_TIME_S:
        return min(
            FINAL_POWER_W,
            STAGE1_TARGET_POWER_W
            + STAGE2_RAMP_W_S * (time_s - STAGE1_END_S),
        )
    return FINAL_POWER_W


def power_ramp_rate_w_s(time_s: float) -> float:
    if float(time_s) < STAGE1_END_S - 1.0e-9:
        return STAGE1_RAMP_W_S
    if float(time_s) < FINAL_POWER_TIME_S - 1.0e-9:
        return STAGE2_RAMP_W_S
    return 0.0


def _apply_power(core: Any, power_w: float) -> None:
    core.update_neutronic_power(
        p_total=float(power_w),
        p_fiss=float(power_w),
        p_decay=0.0,
        alpha=1.0,
    )


def build_case():
    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(
            inlet_temperature_k=INITIAL_TEMPERATURE_K,
            tec_gap_h_eq_w_m2_k=HELIUM_GAP_H_W_M2_K,
            tec_gap_gas="Helium",
            main_tec_enabled=True,
            main_tec_mode="fixed_r",
            main_tec_target_value=LOAD_RESISTANCE_OHM,
            main_tec_current_guess_a=CURRENT_LIMIT_A,
            main_tec_topology="series",
            tec_lookup_enabled=True,
            tec_lookup_db=TEC_LOOKUP_DB,
            tec_lookup_regions=(
                "core", "startup", "high_power", "accident"
            ),
        ),
        flow_config=FullLoopFlowConfig(total_flow_kg_s=FLOW_KG_S),
        pump_config=FullLoopPumpConfig(
            pump_total_head_pa=1.0,
            pump_flow_control=True,
            target_flow_kg_s=FLOW_KG_S,
        ),
        radiator_config=V14HeatPipeRadiatorConfig(
            t_space_k=SPACE_TEMPERATURE_K,
            hp_initial_temp_k=INITIAL_TEMPERATURE_K,
            hp_emissivity=RADIATOR_EMISSIVITY,
            fin_emissivity=RADIATOR_EMISSIVITY,
            lower_hp_down_view_factor=LOWER_HP_DOWN_VIEW_FACTOR,
            external_heat_enabled=True,
            external_heat_period_s=EXTERNAL_HEAT_PERIOD_S,
            external_heat_time_origin_s=0.0,
            external_heat_absorption_efficiency=0.992,
            thermal_shield_enabled=False,
        ),
    )
    core = build["core"]
    gap = next(iter(build["tfes"].values())).couplers["tec_couple"]
    gap_h_w_m2_k = float(gap.k_gas) / float(gap.gap)
    core.point_reactor = None
    core.enable_tec_coupled = False
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    _apply_power(core, INITIAL_POWER_W)
    return build


def _apply_wire_resistance_without_calculate(core: Any) -> None:
    thermo_calc = core.thermo_calc
    wire = (
        np.asarray(WIRE_RESISTANCE_OHM, dtype=float)
        * WIRE_RESISTANCE_SCALE
    )
    thermo_calc._input_data.resistanceWire = np.tile(
        wire, (thermo_calc.N_elem, 1)
    )
    thermo_calc.build()


def _configure_tec(build, mode: str) -> None:
    core = build["core"]
    gap = next(iter(build["tfes"].values())).couplers["tec_couple"]
    gap_h_w_m2_k = float(gap.k_gas) / float(gap.gap)
    target = CURRENT_LIMIT_A if mode == "fixed_i" else LOAD_RESISTANCE_OHM
    core.setup_tec_circuit(
        mode,
        target,
        I_guess=CURRENT_LIMIT_A,
        topology="series",
    )
    _apply_wire_resistance_without_calculate(core)
    core.enable_tec_coupled = True
    core.set_thermo_update_time(
        float(build["system"].global_time)
        - float(core.thermo_update_interval)
    )


def activate_cesium_tec(build, mode: str = "fixed_r") -> None:
    updated = apply_tec_gap_h_eq(
        build["core"], CESIUM_GAP_H_W_M2_K
    )
    if updated != len(build["tfes"]):
        raise RuntimeError(
            f"Updated {updated} TEC gaps; expected {len(build['tfes'])}."
        )
    _configure_tec(build, mode)


def _tec_metrics(core: Any) -> dict:
    if not core.enable_tec_coupled:
        return {
            "tec_current_a": 0.0,
            "tec_voltage_v": 0.0,
            "tec_electric_power_w": 0.0,
            "tec_converged": True,
            "tec_generating": False,
        }
    results = core.get_tec_circuit_global_results().get("main") or {}
    current = float(results.get("Iout", 0.0))
    voltage = float(results.get("Uout", 0.0))
    converged = bool(results.get("converged", False))
    if not converged:
        current = 0.0
        voltage = 0.0
    power = current * voltage
    return {
        "tec_current_a": current,
        "tec_voltage_v": voltage,
        "tec_electric_power_w": power,
        "tec_converged": converged,
        "tec_generating": bool(converged and power > 1.0e-9),
    }


def _collect(build, stage: int, stage_start_s: float, switch_time_s) -> dict:
    system = build["system"]
    core = build["core"]
    gap = next(iter(build["tfes"].values())).couplers["tec_couple"]
    gap_h_w_m2_k = float(gap.k_gas) / float(gap.gap)
    solids = [
        np.asarray(solid.T, dtype=float)
        for solid in system.solid_components.values()
    ]
    tec = _tec_metrics(core)
    row = {
        "time_s": float(system.global_time),
        "stage_elapsed_s": float(system.global_time - stage_start_s),
        "stage": int(stage),
        "power_control_mode": "prescribed_total_power",
        "prescribed_power_w": prescribed_power_w(system.global_time),
        "power_ramp_rate_w_s": power_ramp_rate_w_s(system.global_time),
        "target_flow_kg_s": FLOW_KG_S,
        "pump_a_flow_kg_s": float(build["pump_a"].W),
        "pump_b_flow_kg_s": float(build["pump_b"].W),
        "fluid_min_k": float(np.min(system.fluid_solver.T_vec)),
        "fluid_max_k": float(np.max(system.fluid_solver.T_vec)),
        "solid_min_k": float(min(np.min(value) for value in solids)),
        "solid_max_k": float(max(np.max(value) for value in solids)),
        "tec_gap_gas": (
            "Cesium"
            if math.isclose(gap_h_w_m2_k, CESIUM_GAP_H_W_M2_K)
            else "Helium"
        ),
        "tec_gap_h_w_m2_k": gap_h_w_m2_k,
        "tec_enabled": bool(core.enable_tec_coupled),
        "tec_control_mode": (
            str(core.tec_circuit_mode)
            if core.enable_tec_coupled
            else "off"
        ),
        "tec_current_limit_switch_time_s": (
            "" if switch_time_s is None else float(switch_time_s)
        ),
        "tec_series_element_count": int(core.total_virtual_elements),
        **tec,
    }
    for key, value in row.items():
        if isinstance(value, (str, bool)):
            continue
        if not math.isfinite(float(value)):
            raise FloatingPointError(f"Non-finite {key}: {value}")
    return row


def _prepare_stage(stage: int, restart_in: Path | None, resume: bool = False):
    stage_start_s, stage_end_s, _ = STAGES[stage]
    build = build_case()
    system = build["system"]
    core = build["core"]
    gap = next(iter(build["tfes"].values())).couplers["tec_couple"]
    gap_h_w_m2_k = float(gap.k_gas) / float(gap.gap)
    if stage == 1:
        if restart_in is not None:
            raise ValueError("Stage 1 must start from the 300 K cold state.")
        _set_uniform_temperature(build, INITIAL_TEMPERATURE_K)
        system.initialize_system(
            dt_init=0.01,
            tol=2.0e-5,
            max_iter=1000,
        )
    else:
        if restart_in is None:
            raise ValueError(f"Stage {stage} requires the previous restart.")
        system.load_global_state(str(restart_in))
        for solid in system.solid_components.values():
            solid.set_ode_method("implicit_euler")
    current_time_s = float(system.global_time)
    valid_start = math.isclose(
        current_time_s, stage_start_s, rel_tol=0.0, abs_tol=2.0e-6
    )
    valid_resume = (
        resume and stage_start_s < current_time_s < stage_end_s
    )
    if not (valid_start or valid_resume):
        raise ValueError(
            f"Stage {stage} expected {stage_start_s} or an in-stage resume, "
            f"got {system.global_time}."
        )
    build["pump_a"].set_flow_rate(FLOW_KG_S)
    apply_tec_gap_h_eq(core, HELIUM_GAP_H_W_M2_K)
    core.enable_tec_coupled = False
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    switch_time_s = None
    if stage == 3:
        loaded_mode = str(getattr(core, "tec_circuit_mode", "fixed_r"))
        mode = "fixed_i" if loaded_mode == "fixed_i" else "fixed_r"
        activate_cesium_tec(build, mode)
        if mode == "fixed_i":
            switch_time_s = float(system.global_time)
    _apply_power(core, prescribed_power_w(system.global_time))
    return build, switch_time_s


def _load_history_reference(output_dir: Path, build) -> dict:
    net = build["system"].fluid_solver
    volume_rows = {}
    junction_rows = {}
    with (output_dir / "history_coolant.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        for row in csv.DictReader(handle):
            if row["entity_type"] == "volume":
                volume_rows.setdefault(row["name"], row)
            else:
                junction_rows.setdefault(row["name"], row)
            if (
                len(volume_rows) == len(net.volumes_obj)
                and len(junction_rows) == len(net.junctions_obj)
            ):
                break
    if (
        len(volume_rows) != len(net.volumes_obj)
        or len(junction_rows) != len(net.junctions_obj)
    ):
        raise ValueError("Existing coolant history lacks reference entities.")

    with (output_dir / "history_reactivity.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        first_reactivity = next(csv.DictReader(handle))
    return {
        "feedback_reference_total": float(
            first_reactivity["feedback_total_absolute"]
        ),
        "reference_fluid": {
            "T": np.asarray([
                float(volume_rows[volume.name]["reference_temperature_K"])
                for volume in net.volumes_obj
            ]),
            "P": np.asarray([
                float(volume_rows[volume.name]["reference_pressure_Pa"])
                for volume in net.volumes_obj
            ]),
            "h": np.asarray([
                float(volume_rows[volume.name]["reference_enthalpy_J_kg"])
                for volume in net.volumes_obj
            ]),
            "W": np.asarray([
                float(junction_rows[junction.name]["reference_mass_flow_kg_s"])
                for junction in net.junctions_obj
            ]),
        },
        "tec_open_circuit_active": False,
        "tec_open_circuit_time_s": float("nan"),
    }


def run_stage(
    stage: int,
    *,
    restart_in: Path | None,
    output_dir: Path,
    max_dt_s: float = 0.2,
    record_interval_s: float = 1.0,
    checkpoint_interval_s: float = 50.0,
    resume: bool = False,
    history_reference_restart: Path | None = None,
) -> Path:
    if output_dir.exists() and not resume:
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if resume and not output_dir.exists():
        raise FileNotFoundError(f"Resume output directory is missing: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=resume)
    stage_start_s, stage_end_s, _ = STAGES[stage]
    build, switch_time_s = _prepare_stage(stage, restart_in, resume=resume)
    system = build["system"]
    core = build["core"]
    gap = next(iter(build["tfes"].values())).couplers["tec_couple"]
    gap_h_w_m2_k = float(gap.k_gas) / float(gap.gap)
    if resume:
        history_reference = _load_history_reference(output_dir, build)
    else:
        history_reference = capture_v14_history_reference(build)
    elapsed_s = float(system.global_time) - stage_start_s
    next_record_s = stage_start_s + (
        math.floor(elapsed_s / record_interval_s + 1.0e-8) + 1
    ) * record_interval_s
    next_checkpoint_s = stage_start_s + (
        math.floor(elapsed_s / checkpoint_interval_s + 1.0e-8) + 1
    ) * checkpoint_interval_s
    checkpoints = [str(path) for path in sorted(output_dir.glob("checkpoint_t*.npz"))]
    switch_restart_path = None

    latest = _collect(build, stage, stage_start_s, switch_time_s)
    if not resume:
        append_v14_system_history(
            output_dir, latest, build, history_reference, stage_start_s
        )
    print(json.dumps(latest, sort_keys=True), flush=True)

    while system.global_time < stage_end_s - 1.0e-10:
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=float(max_dt_s),
            safety_factor=0.8,
            respect_fluid_cfl=False,
        )
        dt = min(
            float(dt),
            stage_end_s - float(system.global_time),
            next_record_s - float(system.global_time),
            next_checkpoint_s - float(system.global_time),
        )
        next_time = float(system.global_time) + dt
        _apply_power(core, prescribed_power_w(next_time))
        system.step(
            dt,
            inner_iter=1,
            fail_on_fluid_nonconvergence=False,
            fluid_max_iter=1,
        )

        tec = _tec_metrics(core)
        if (
            stage == 3
            and core.tec_circuit_mode == "fixed_r"
            and tec["tec_current_a"] >= CURRENT_LIMIT_A
        ):
            switch_time_s = float(system.global_time)
            switch_restart = output_dir / (
                f"checkpoint_tec_switch_t{switch_time_s:09.3f}s.npz"
            )
            system.save_global_state(str(switch_restart))
            switch_restart_path = str(switch_restart)
            checkpoints.append(switch_restart_path)
            _configure_tec(build, "fixed_i")
            print(
                f"TEC switched to fixed_i=216 A at "
                f"t={switch_time_s:.6f} s after saving {switch_restart}",
                flush=True,
            )

        if system.global_time >= next_record_s - 1.0e-10:
            latest = _collect(
                build, stage, stage_start_s, switch_time_s
            )
            append_v14_system_history(
                output_dir,
                latest,
                build,
                history_reference,
                stage_start_s,
            )
            print(json.dumps(latest, sort_keys=True), flush=True)
            while next_record_s <= system.global_time + 1.0e-10:
                next_record_s += record_interval_s

        if (
            system.global_time >= next_checkpoint_s - 1.0e-10
            and system.global_time < stage_end_s - 1.0e-10
        ):
            checkpoint = output_dir / (
                f"checkpoint_t{system.global_time:09.3f}s.npz"
            )
            checkpoint = output_dir / (
                f"checkpoint_t{system.global_time:09.3f}s.npz"
            )
            system.save_global_state(str(checkpoint))
            checkpoints.append(str(checkpoint))
            while next_checkpoint_s <= system.global_time + 1.0e-10:
                next_checkpoint_s += checkpoint_interval_s

    if float(latest["time_s"]) < float(system.global_time) - 1.0e-10:
        latest = _collect(build, stage, stage_start_s, switch_time_s)
        append_v14_system_history(
            output_dir, latest, build, history_reference, stage_start_s
        )
        print(json.dumps(latest, sort_keys=True), flush=True)

    final_restart = output_dir / "final_restart.npz"
    system.save_global_state(str(final_restart))
    latest = _collect(build, stage, stage_start_s, switch_time_s)
    summary = {
        "stage": stage,
        "stage_start_s": stage_start_s,
        "stage_end_s": float(system.global_time),
        "restart_in": None if restart_in is None else str(restart_in),
        "restart_out": str(final_restart),
        "checkpoint_paths": checkpoints,
        "record_interval_s": record_interval_s,
        "checkpoint_interval_s": checkpoint_interval_s,
        "tec_current_limit_switch_time_s": switch_time_s,
        "tec_current_limit_switch_restart": switch_restart_path,
        "latest": latest,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return final_restart


def self_test() -> None:
    assert prescribed_power_w(0.0) == INITIAL_POWER_W
    assert math.isclose(
        prescribed_power_w(STAGE1_END_S), STAGE1_TARGET_POWER_W
    )
    assert math.isclose(
        prescribed_power_w(FINAL_POWER_TIME_S), FINAL_POWER_W
    )
    assert prescribed_power_w(CESIUM_TEC_TIME_S) == FINAL_POWER_W
    assert STAGES[1][1] == STAGES[2][0]
    assert STAGES[2][1] == STAGES[3][0]
    assert STAGES[3][1] == 10_000.0
    assert WIRE_RESISTANCE_SCALE == 0.335
    assert SPACE_TEMPERATURE_K == 4.0
    assert RADIATOR_EMISSIVITY == 0.7475
    assert LOWER_HP_DOWN_VIEW_FACTOR == 0.4
    assert sum((1, 6, 9, 18, 24)) == 58
    print("V14 startup stage checks passed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("1", "2", "3", "all"), default="all"
    )
    parser.add_argument("--base-dir", type=Path, default=CASE_DIR)
    parser.add_argument("--restart-in", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--history-reference-restart", type=Path)
    parser.add_argument("--max-dt", type=float, default=0.2)
    parser.add_argument("--record-interval", type=float, default=1.0)
    parser.add_argument(
        "--checkpoint-interval", type=float, default=50.0
    )
    parser.add_argument("--stage-start", type=float)
    parser.add_argument("--stage-end", type=float)
    parser.add_argument("--output-name")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0

    if any(
        value is not None
        for value in (args.stage_start, args.stage_end, args.output_name)
    ):
        if args.stage == "all":
            parser.error("Stage overrides require one explicit --stage.")
        stage = int(args.stage)
        default_start, default_end, default_name = STAGES[stage]
        stage_start = default_start if args.stage_start is None else args.stage_start
        stage_end = default_end if args.stage_end is None else args.stage_end
        if stage_end <= stage_start:
            parser.error("--stage-end must be greater than --stage-start.")
        STAGES[stage] = (
            float(stage_start),
            float(stage_end),
            default_name if args.output_name is None else args.output_name,
        )

    stages = (
        (1, 2, 3) if args.stage == "all" else (int(args.stage),)
    )
    restart = args.restart_in
    for stage in stages:
        _, _, name = STAGES[stage]
        restart = run_stage(
            stage,
            restart_in=restart,
            output_dir=args.base_dir / name,
            max_dt_s=args.max_dt,
            resume=args.resume,
            history_reference_restart=args.history_reference_restart,
            record_interval_s=args.record_interval,
            checkpoint_interval_s=args.checkpoint_interval,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
