import argparse
import csv
import os
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_history(npz_path: str):
    data = np.load(npz_path)
    return {
        "time": np.asarray(data["time"], dtype=float),
        "y_centers": np.asarray(data["y_centers"], dtype=float),
        "wick_k_history": np.asarray(data["wick_k_history"], dtype=float),
        "temperature_history": np.asarray(data["temperature_history"], dtype=float),
        "input_power": float(np.asarray(data["input_power"], dtype=float)[0]),
        "condenser_h": float(np.asarray(data["condenser_h"], dtype=float)[0]),
        "condenser_temp": float(np.asarray(data["condenser_temp"], dtype=float)[0]),
        "initial_temp": float(np.asarray(data["initial_temp"], dtype=float)[0]),
    }


def sample_late_stage(history: dict, interval: float, start_time: float):
    time = history["time"]
    y_centers = history["y_centers"]
    wick_k_history = history["wick_k_history"][:, 0, :]
    wick_temperature_history = history["temperature_history"][:, 0, :]
    t_final = float(time[-1])

    target_times = np.arange(start_time, t_final + 0.5 * interval, interval, dtype=float)

    sampled_times = []
    sampled_k_curves = []
    sampled_t_curves = []

    for target_t in target_times:
        idx = int(np.argmin(np.abs(time - target_t)))
        sampled_times.append(float(time[idx]))
        sampled_k_curves.append(np.array(wick_k_history[idx], copy=True))
        sampled_t_curves.append(np.array(wick_temperature_history[idx], copy=True))

    return (
        y_centers,
        np.asarray(sampled_times, dtype=float),
        np.asarray(sampled_k_curves, dtype=float),
        np.asarray(sampled_t_curves, dtype=float),
    )


def build_reference_wick_material():
    startup_script_path = Path(CURRENT_DIR) / "test_single_hp_power_convection_startup.py"
    spec = importlib.util.spec_from_file_location("single_hp_power_convection_startup", startup_script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load startup script: {startup_script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    hp = module.build_startup_hp(
        initial_temp=293.15,
        input_power=1800.0,
        condenser_h=147.67,
        condenser_temp=293.15,
    )
    hp.wick_mat._ensure_lookup_table()
    return hp.wick_mat


def write_sample_csv(
    csv_path: str,
    y_centers: np.ndarray,
    sampled_times: np.ndarray,
    sampled_k_curves: np.ndarray,
    sampled_t_curves: np.ndarray,
):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    fieldnames = ["axial_position_m"]
    for t in sampled_times:
        tag = f"{int(round(t))}s"
        fieldnames.extend([f"T_wick_{tag}", f"k_wick_{tag}"])

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, y in enumerate(y_centers):
            row = {"axial_position_m": float(y)}
            for j, t in enumerate(sampled_times):
                tag = f"{int(round(t))}s"
                row[f"T_wick_{tag}"] = float(sampled_t_curves[j, i])
                row[f"k_wick_{tag}"] = float(sampled_k_curves[j, i])
            writer.writerow(row)


def plot_late_stage(
    y_centers: np.ndarray,
    sampled_times: np.ndarray,
    sampled_k_curves: np.ndarray,
    output_path: str,
    title: str,
):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=180)
    colors = plt.cm.plasma(np.linspace(0.0, 1.0, len(sampled_times)))
    x_mm = y_centers * 1000.0

    for color, t, curve in zip(colors, sampled_times, sampled_k_curves):
        ax.plot(
            x_mm,
            curve,
            color=color,
            linewidth=1.8,
            marker="o",
            markersize=3.2,
            label=f"{t:.0f} s",
        )

    ax.set_xlabel("Axial Position [mm]")
    ax.set_ylabel("Wick Effective Conductivity [W/m/K]")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(ncol=2, fontsize=8, frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot late-stage wick conductivity distribution.")
    parser.add_argument(
        "--npz-path",
        default=os.path.join(CURRENT_DIR, "single_hp_power_convection_startup_3600s.npz"),
    )
    parser.add_argument(
        "--output-path",
        default=os.path.join(CURRENT_DIR, "single_hp_power_convection_wick_k_late_stage.png"),
    )
    parser.add_argument(
        "--csv-path",
        default=os.path.join(CURRENT_DIR, "single_hp_power_convection_wick_k_late_stage.csv"),
    )
    parser.add_argument("--interval", type=float, default=300.0)
    parser.add_argument("--start-time", type=float, default=2400.0)
    return parser.parse_args()


def main():
    args = parse_args()
    history = load_history(args.npz_path)
    y_centers, sampled_times, sampled_k_curves, sampled_t_curves = sample_late_stage(
        history,
        interval=args.interval,
        start_time=args.start_time,
    )
    write_sample_csv(args.csv_path, y_centers, sampled_times, sampled_k_curves, sampled_t_curves)

    title = (
        f"Late-Stage Wick Conductivity Every {args.interval:.0f} s "
        f"({args.start_time:.0f}-{sampled_times[-1]:.0f} s)"
    )
    plot_late_stage(
        y_centers=y_centers,
        sampled_times=sampled_times,
        sampled_k_curves=sampled_k_curves,
        output_path=args.output_path,
        title=title,
    )

    wick_mat = build_reference_wick_material()
    print(f"PNG written to: {args.output_path}")
    print(f"Sample CSV written to: {args.csv_path}")
    print(f"Sampled times: {sampled_times.tolist()}")
    print(
        "High-nonlinearity temperature band [K]: "
        f"{wick_mat._lookup_high_nonlinear_min:.6f} to {wick_mat._lookup_high_nonlinear_max:.6f}"
    )
    print("Conductivity clipping upper bound [W/m/K]: 1000000.000000")


if __name__ == "__main__":
    main()
