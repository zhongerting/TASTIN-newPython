import argparse
import os
import sys
import traceback
import warnings
from dataclasses import dataclass

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_6segment_geometry100hp_potassium_mixed as base
from Components.basicComponents.HeatPipe2D import HeatPipe2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D


DEFAULT_RESTART = os.path.join(
    base.current_dir,
    "collector_ring_6segment_geometry100hp_potassium_mixed_500s_resume_restart.npz",
)


@dataclass(frozen=True)
class WarningTestCase:
    key: str
    anisotropic: bool
    conductivity_cap: float | None
    frozen_properties: bool
    integrator: str = "bdf"
    theta: float = 0.6
    implicit_boundary: bool = False

    @property
    def label(self) -> str:
        cap_label = "uncapped" if self.conductivity_cap is None else f"cap{self.conductivity_cap:.1e}".replace(".", "p")
        mode = "aniso" if self.anisotropic else "iso"
        frozen = "frozen" if self.frozen_properties else "liveprops"
        boundary = "ibc" if self.implicit_boundary else "explicitbc"
        integrator_label = self.integrator
        if self.integrator == "theta_implicit":
            integrator_label = f"theta{self.theta:.2f}".replace(".", "p")
        return f"{self.key}_{mode}_{cap_label}_{frozen}_{integrator_label}_{boundary}"


CASES = {
    "A": WarningTestCase("A", anisotropic=False, conductivity_cap=1.2e6, frozen_properties=False),
    "B": WarningTestCase("B", anisotropic=True, conductivity_cap=None, frozen_properties=False),
    "C": WarningTestCase("C", anisotropic=True, conductivity_cap=1.2e6, frozen_properties=False),
    "D": WarningTestCase("D", anisotropic=True, conductivity_cap=1.2e6, frozen_properties=True),
    "E": WarningTestCase("E", anisotropic=True, conductivity_cap=1.2e6, frozen_properties=False, integrator="implicit_euler"),
    "F": WarningTestCase("F", anisotropic=True, conductivity_cap=None, frozen_properties=False, integrator="implicit_euler"),
    "G": WarningTestCase("G", anisotropic=True, conductivity_cap=1.2e6, frozen_properties=False, integrator="theta_implicit", theta=0.6),
    "H": WarningTestCase("H", anisotropic=True, conductivity_cap=None, frozen_properties=False, integrator="theta_implicit", theta=0.6),
    "I": WarningTestCase("I", anisotropic=True, conductivity_cap=None, frozen_properties=False, integrator="theta_implicit", theta=0.6, implicit_boundary=True),
    "J": WarningTestCase("J", anisotropic=True, conductivity_cap=None, frozen_properties=False, integrator="implicit_euler", implicit_boundary=True),
}


def iter_heat_pipe_solids(model):
    for ring_hp in model["ring_hps"]:
        for hp_unit in ring_hp.hp_units:
            yield hp_unit.hp


def configure_heat_pipes(model, test_case: WarningTestCase):
    count = 0
    for hp in iter_heat_pipe_solids(model):
        hp.set_wick_conductivity_mode(test_case.anisotropic)
        hp.wick_mat.set_conductivity_cap(test_case.conductivity_cap)
        hp.enable_frozen_property_correction = test_case.frozen_properties
        hp.set_time_integrator(test_case.integrator)
        hp.set_implicit_boundary_linearization(test_case.implicit_boundary)
        if test_case.integrator == "theta_implicit":
            hp.set_theta_implicit_value(test_case.theta)
        count += 1
    return count


def format_array_range(obj, attr_name):
    if not hasattr(obj, attr_name):
        return ""
    arr = np.asarray(getattr(obj, attr_name), dtype=float)
    if arr.size == 0:
        return f"{attr_name}=empty"
    return (
        f"{attr_name}[min,max]=({np.nanmin(arr):.6e},{np.nanmax(arr):.6e}) "
        f"finite={bool(np.isfinite(arr).all())}"
    )


def build_solid_summary(solid, dt, test_case: WarningTestCase):
    temperature = np.asarray(solid.T, dtype=float)
    lines = [
        f"case={test_case.label}",
        f"dt={float(dt):.12g}",
        f"name={getattr(solid, 'name', type(solid).__name__)}",
        f"type={type(solid).__name__}",
        f"current_time={float(getattr(solid, 'current_time', np.nan)):.12g}",
        (
            "temperature[min,max,dT]="
            f"({float(np.nanmin(temperature)):.9f},"
            f"{float(np.nanmax(temperature)):.9f},"
            f"{float(np.nanmax(temperature) - np.nanmin(temperature)):.9f})"
        ),
    ]

    for attr_name in ("thermal_capacitance", "G_x_inner", "G_y_inner", "Q_source", "dTdt"):
        text = format_array_range(solid, attr_name)
        if text:
            lines.append(text)

    if isinstance(solid, HeatPipe2D):
        temp_2d = temperature.reshape(solid.shape_nodes)
        wick_temperature = temp_2d[: solid.n_wick, :]
        k_axial = np.asarray(solid.wick_mat.conductivity_axial(wick_temperature), dtype=float)
        k_structural = np.asarray(solid.wick_mat.conductivity_radial(wick_temperature), dtype=float)
        if solid.use_anisotropic_wick_conductivity:
            k_radial_model = k_structural
        else:
            k_radial_model = k_axial
        lines.extend(
            [
                f"wick_mode={'anisotropic' if solid.use_anisotropic_wick_conductivity else 'isotropic'}",
                f"frozen_property_correction={bool(solid.enable_frozen_property_correction)}",
                f"implicit_boundary_linearization={bool(solid.implicit_boundary_linearization)}",
                f"wick_temperature[min,max]=({float(np.nanmin(wick_temperature)):.9f},{float(np.nanmax(wick_temperature)):.9f})",
                f"k_axial[min,max]=({float(np.nanmin(k_axial)):.6e},{float(np.nanmax(k_axial)):.6e})",
                f"k_radial_model[min,max]=({float(np.nanmin(k_radial_model)):.6e},{float(np.nanmax(k_radial_model)):.6e})",
                f"k_structural[min,max]=({float(np.nanmin(k_structural)):.6e},{float(np.nanmax(k_structural)):.6e})",
            ]
        )

    return "\n".join(lines)


def write_first_warning(path, summary, exc):
    with open(path, "w", encoding="utf-8") as f:
        f.write(summary)
        f.write("\n\ntraceback:\n")
        traceback.print_exception(type(exc), exc, exc.__traceback__, limit=10, file=f)


def run_stage(
    *,
    test_case: WarningTestCase,
    restart_from,
    t_end,
    case_name,
    csv_path,
    restart_save_path,
    warning_as_error,
    first_warning_path,
):
    original_build_model = base.build_model
    original_hc_step = HeatConduction2D.step
    original_hp_solve = HeatPipe2D._solve_ivp_step
    first_warning_written = False

    def wrapped_build_model():
        model = original_build_model()
        count = configure_heat_pipes(model, test_case)
        print(
                f"Configured case {test_case.label}: heat_pipes={count}, "
                f"anisotropic={test_case.anisotropic}, cap={test_case.conductivity_cap}, "
                f"frozen_properties={test_case.frozen_properties}, integrator={test_case.integrator}, "
                f"theta={test_case.theta}, implicit_boundary={test_case.implicit_boundary}"
            )
        return model

    def record_warning_once(solid, dt, exc):
        nonlocal first_warning_written
        if not first_warning_written:
            summary = build_solid_summary(solid, dt, test_case)
            print("FIRST_WARNING")
            print(summary)
            write_first_warning(first_warning_path, summary, exc)
            first_warning_written = True

    def patched_hc_step(self, dt, method="BDF", **kwargs):
        try:
            return original_hc_step(self, dt, method=method, **kwargs)
        except RuntimeWarning as exc:
            record_warning_once(self, dt, exc)
            raise

    def patched_hp_solve(self, dt, method="BDF", **kwargs):
        try:
            return original_hp_solve(self, dt, method=method, **kwargs)
        except RuntimeWarning as exc:
            record_warning_once(self, dt, exc)
            raise

    base.build_model = wrapped_build_model
    HeatConduction2D.step = patched_hc_step
    HeatPipe2D._solve_ivp_step = patched_hp_solve
    try:
        with warnings.catch_warnings(record=True) as warning_records:
            if warning_as_error:
                warnings.simplefilter("error", RuntimeWarning)
            else:
                warnings.simplefilter("always", RuntimeWarning)
            model, history = base.run_case(
                case_name=case_name,
                t_end=t_end,
                min_dt=1.0e-3,
                max_dt=0.5,
                safety_factor=1.0,
                inner_iter=2,
                print_every_time=5.0,
                csv_path=csv_path,
                restart_from=restart_from,
                restart_save_path=restart_save_path,
                restart_save_every=999.0,
            )
        bdf_count = sum(
            "invalid value encountered in subtract" in str(record.message)
            for record in warning_records
        )
        return {
            "ok": True,
            "history": history,
            "warning_count": len(warning_records),
            "bdf_warning_count": bdf_count,
            "first_warning_written": first_warning_written,
        }
    except RuntimeWarning as exc:
        return {
            "ok": False,
            "history": [],
            "warning_count": None,
            "bdf_warning_count": None,
            "first_warning_written": first_warning_written,
            "exception": repr(exc),
        }
    finally:
        base.build_model = original_build_model
        HeatConduction2D.step = original_hc_step
        HeatPipe2D._solve_ivp_step = original_hp_solve


def run_test_case(test_case: WarningTestCase, restart_from, transition_end, validation_end):
    transition_name = f"warning_test_{test_case.label}_{int(round(transition_end))}s"
    validation_name = f"warning_test_{test_case.label}_{int(round(validation_end))}s"
    transition_restart = os.path.join(base.current_dir, f"{transition_name}_restart.npz")
    validation_restart = os.path.join(base.current_dir, f"{validation_name}_restart.npz")
    first_warning_path = os.path.join(base.current_dir, f"{validation_name}_first_warning.txt")

    transition = run_stage(
        test_case=test_case,
        restart_from=restart_from,
        t_end=transition_end,
        case_name=transition_name,
        csv_path=os.path.join(base.current_dir, f"{transition_name}_history.csv"),
        restart_save_path=transition_restart,
        warning_as_error=False,
        first_warning_path=os.path.join(base.current_dir, f"{transition_name}_first_warning.txt"),
    )
    if not transition["ok"]:
        print(f"TRANSITION_FAILED {test_case.label}: {transition.get('exception')}")
        return transition

    last_transition = transition["history"][-1] if transition["history"] else {}
    print(
        f"TRANSITION_DONE {test_case.label}: warnings={transition['warning_count']} "
        f"bdf={transition['bdf_warning_count']} "
        f"Tout={float(last_transition.get('T_out_avg', np.nan)):.6f}"
    )

    validation = run_stage(
        test_case=test_case,
        restart_from=transition_restart,
        t_end=validation_end,
        case_name=validation_name,
        csv_path=os.path.join(base.current_dir, f"{validation_name}_history.csv"),
        restart_save_path=validation_restart,
        warning_as_error=True,
        first_warning_path=first_warning_path,
    )
    if validation["ok"]:
        last = validation["history"][-1]
        delta_t = float(last["T_inlet_buffer_out"]) - float(last["T_out_avg"])
        print(
            f"VALIDATION_OK {test_case.label}: final_time={float(last['time']):.6f} "
            f"Tout={float(last['T_out_avg']):.6f} dT={delta_t:.6f}"
        )
    else:
        print(f"VALIDATION_WARNING {test_case.label}: {validation.get('exception')}")
        print(f"FIRST_WARNING_PATH {first_warning_path}")
    return validation


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run warning sensitivity matrix for the 6-segment potassium collector-ring case."
    )
    parser.add_argument("--case", choices=sorted(CASES), required=True, help="Case key to run.")
    parser.add_argument("--restart-from", default=DEFAULT_RESTART, help="Input restart path for transition stage.")
    parser.add_argument("--transition-end", type=float, default=510.0, help="Transition stage absolute end time [s].")
    parser.add_argument("--validation-end", type=float, default=530.0, help="Validation stage absolute end time [s].")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_test_case(
        CASES[args.case],
        restart_from=args.restart_from,
        transition_end=args.transition_end,
        validation_end=args.validation_end,
    )
