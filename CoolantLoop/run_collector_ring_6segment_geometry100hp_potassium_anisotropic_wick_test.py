import argparse
import os
import sys
import warnings

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_6segment_geometry100hp_potassium_mixed as base


DEFAULT_RESTART = os.path.join(
    base.current_dir,
    "collector_ring_6segment_geometry100hp_potassium_mixed_500s_resume_restart.npz",
)


def iter_heat_pipe_solids(model):
    for ring_hp in model["ring_hps"]:
        for hp_unit in ring_hp.hp_units:
            yield hp_unit.hp


def apply_wick_conductivity_mode(model, anisotropic):
    count = 0
    for hp in iter_heat_pipe_solids(model):
        hp.set_wick_conductivity_mode(anisotropic)
        count += 1
    return count


def collect_wick_diagnostics(model, anisotropic, top_n=8):
    rows = []
    for hp in iter_heat_pipe_solids(model):
        temperature = np.asarray(hp.T, dtype=float).reshape(hp.shape_nodes)
        wick_temperature = temperature[:hp.n_wick, :]
        k_axial = np.asarray(hp.wick_mat.conductivity_axial(wick_temperature), dtype=float)
        k_structural = np.asarray(hp.wick_mat.conductivity_radial(wick_temperature), dtype=float)
        k_radial_model = k_structural if anisotropic else k_axial
        rows.append(
            {
                "name": hp.name,
                "T_wick_min": float(np.min(wick_temperature)),
                "T_wick_max": float(np.max(wick_temperature)),
                "k_axial_max": float(np.max(k_axial)),
                "k_axial_mean": float(np.mean(k_axial)),
                "k_radial_model_max": float(np.max(k_radial_model)),
                "k_radial_model_mean": float(np.mean(k_radial_model)),
                "k_structural_max": float(np.max(k_structural)),
                "k_structural_mean": float(np.mean(k_structural)),
                "axial_delta_T": float(np.max(temperature) - np.min(temperature)),
            }
        )
    rows.sort(key=lambda row: row["k_axial_max"], reverse=True)
    return rows[:top_n]


def print_wick_diagnostics(model, anisotropic, top_n=8):
    rows = collect_wick_diagnostics(model, anisotropic=anisotropic, top_n=top_n)
    print("Top wick-conductivity diagnostics:")
    for row in rows:
        print(
            "  {name}: T_wick=[{T_wick_min:.3f}, {T_wick_max:.3f}] K | "
            "k_axial_max={k_axial_max:.6e} W/m/K | "
            "k_radial_model_max={k_radial_model_max:.6e} W/m/K | "
            "k_structural_max={k_structural_max:.6e} W/m/K | "
            "dT_hp={axial_delta_T:.6f} K".format(**row)
        )


def build_mode_case(base_build_model, anisotropic):
    model = base_build_model()
    count = apply_wick_conductivity_mode(model, anisotropic)
    mode = "anisotropic" if anisotropic else "isotropic"
    print(f"Wick conductivity mode: {mode}; configured heat-pipe solids: {count}")
    return model


def run_case(
    *,
    anisotropic=True,
    t_end=501.0,
    restart_from=DEFAULT_RESTART,
    case_name=None,
    csv_path=None,
    restart_save_path=None,
    restart_save_every=25.0,
    min_dt=1.0e-3,
    max_dt=0.5,
    safety_factor=1.0,
    inner_iter=2,
    print_every_time=1.0,
    diagnostics_top_n=8,
):
    mode = "anisowick" if anisotropic else "isowick_control"
    if case_name is None:
        case_name = f"collector_ring_6segment_geometry100hp_potassium_{mode}_{int(round(t_end))}s"
    if csv_path is None:
        csv_path = os.path.join(base.current_dir, f"{case_name}_history.csv")
    if restart_save_path is None:
        restart_save_path = os.path.join(base.current_dir, f"{case_name}_restart.npz")

    original_build_model = base.build_model

    def wrapped_build_model():
        return build_mode_case(original_build_model, anisotropic)

    base.build_model = wrapped_build_model
    try:
        with warnings.catch_warnings(record=True) as warning_records:
            warnings.simplefilter("always", RuntimeWarning)
            model, history = base.run_case(
                case_name=case_name,
                t_end=t_end,
                min_dt=min_dt,
                max_dt=max_dt,
                safety_factor=safety_factor,
                inner_iter=inner_iter,
                print_every_time=print_every_time,
                csv_path=csv_path,
                restart_from=restart_from,
                restart_save_path=restart_save_path,
                restart_save_every=restart_save_every,
            )
    finally:
        base.build_model = original_build_model

    bdf_warning_count = sum(
        "invalid value encountered in subtract" in str(record.message)
        for record in warning_records
    )
    print("=" * 70)
    print("Anisotropic wick test summary")
    print("=" * 70)
    print(f"Mode              : {'anisotropic' if anisotropic else 'isotropic control'}")
    print(f"RuntimeWarnings   : {len(warning_records)}")
    print(f"BDF invalid warns : {bdf_warning_count}")
    if history:
        last = history[-1]
        delta_t = float(last["T_inlet_buffer_out"]) - float(last["T_out_avg"])
        print(f"Final time        : {float(last['time']):.6f} s")
        print(f"T_out_avg         : {float(last['T_out_avg']):.6f} K")
        print(f"Inlet-outlet dT   : {delta_t:.6f} K")
        print(f"W_in/W_out        : {float(last['W_in_total']):.9f} / {float(last['W_out_total']):.9f} kg/s")
        print(f"W_ring_in/out     : {float(last['W_ring_in_total']):.9f} / {float(last['W_ring_out_total']):.9f} kg/s")
    print_wick_diagnostics(model, anisotropic=anisotropic, top_n=diagnostics_top_n)
    print("=" * 70)
    return model, history, warning_records


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the 6-segment 100-KHP collector-ring case with optional anisotropic wick conductivity."
    )
    parser.add_argument("--isotropic-control", action="store_true", help="Keep isotropic wick conductivity for A/B control.")
    parser.add_argument("--t-end", type=float, default=501.0, help="Absolute end time [s].")
    parser.add_argument("--restart-from", default=DEFAULT_RESTART, help="Input restart path.")
    parser.add_argument("--case-name", default=None, help="Output case name prefix.")
    parser.add_argument("--csv-path", default=None, help="Output history CSV path.")
    parser.add_argument("--restart-save-path", default=None, help="Output final restart path.")
    parser.add_argument("--restart-save-every", type=float, default=25.0, help="Checkpoint interval [s].")
    parser.add_argument("--min-dt", type=float, default=1.0e-3, help="Minimum adaptive timestep [s].")
    parser.add_argument("--max-dt", type=float, default=0.5, help="Maximum adaptive timestep [s].")
    parser.add_argument("--safety-factor", type=float, default=1.0, help="Adaptive timestep safety factor.")
    parser.add_argument("--inner-iter", type=int, default=2, help="SystemManager inner coupling iterations.")
    parser.add_argument("--print-every-time", type=float, default=1.0, help="Progress print interval [s].")
    parser.add_argument("--diagnostics-top-n", type=int, default=8, help="Number of heat pipes in wick diagnostics.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_case(
        anisotropic=not args.isotropic_control,
        t_end=args.t_end,
        restart_from=args.restart_from,
        case_name=args.case_name,
        csv_path=args.csv_path,
        restart_save_path=args.restart_save_path,
        restart_save_every=args.restart_save_every,
        min_dt=args.min_dt,
        max_dt=args.max_dt,
        safety_factor=args.safety_factor,
        inner_iter=args.inner_iter,
        print_every_time=args.print_every_time,
        diagnostics_top_n=args.diagnostics_top_n,
    )
