import argparse
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_history(npz_path: str):
    data = np.load(npz_path)
    return {
        "time": np.asarray(data["time"], dtype=float),
        "y_centers": np.asarray(data["y_centers"], dtype=float),
        "outer_wall_temperature_history": np.asarray(data["outer_wall_temperature_history"], dtype=float),
        "initial_temp": float(np.asarray(data["initial_temp"], dtype=float)[0]),
    }


def sample_wall_temperature_every(history: dict, interval: float, include_t0: bool = True):
    time = history["time"]
    y_centers = history["y_centers"]
    wall_history = history["outer_wall_temperature_history"]
    t_final = float(time[-1])

    target_times = np.arange(interval, t_final + 0.5 * interval, interval, dtype=float)

    sampled_times = []
    sampled_curves = []

    if include_t0:
        sampled_times.append(0.0)
        sampled_curves.append(np.full_like(y_centers, history["initial_temp"], dtype=float))

    for target_t in target_times:
        idx = int(np.argmin(np.abs(time - target_t)))
        sampled_times.append(float(time[idx]))
        sampled_curves.append(np.array(wall_history[idx], copy=True))

    return y_centers, np.asarray(sampled_times, dtype=float), np.asarray(sampled_curves, dtype=float)


def write_sample_csv(csv_path: str, y_centers: np.ndarray, sampled_times: np.ndarray, sampled_curves: np.ndarray):
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    fieldnames = ["axial_position_m"] + [f"T_wall_{int(round(t))}s" for t in sampled_times]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, y in enumerate(y_centers):
            row = {"axial_position_m": float(y)}
            for j, t in enumerate(sampled_times):
                row[f"T_wall_{int(round(t))}s"] = float(sampled_curves[j, i])
            writer.writerow(row)


def plot_wall_temperature(y_centers: np.ndarray,
                          sampled_times: np.ndarray,
                          sampled_curves: np.ndarray,
                          output_path: str,
                          title: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=180)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(sampled_times)))
    x_mm = y_centers * 1000.0

    for color, t, curve in zip(colors, sampled_times, sampled_curves):
        ax.plot(
            x_mm,
            curve,
            color=color,
            linewidth=1.8,
            marker="o",
            markersize=3.2,
            label=f"{t:.0f} s"
        )

    ax.set_xlabel("Axial Position [mm]")
    ax.set_ylabel("Outer Wall Temperature [K]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8, frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Plot outer wall temperature every N seconds.")
    parser.add_argument(
        "--npz-path",
        default=os.path.join(CURRENT_DIR, "single_hp_cold_start_1000s.npz"),
    )
    parser.add_argument(
        "--output-path",
        default=os.path.join(CURRENT_DIR, "single_hp_cold_start_wall_temperature_every100s.png"),
    )
    parser.add_argument(
        "--csv-path",
        default=os.path.join(CURRENT_DIR, "single_hp_cold_start_wall_temperature_every100s.csv"),
    )
    parser.add_argument("--interval", type=float, default=100.0)
    parser.add_argument("--no-t0", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    history = load_history(args.npz_path)
    y_centers, sampled_times, sampled_curves = sample_wall_temperature_every(
        history,
        interval=args.interval,
        include_t0=not args.no_t0,
    )

    title = f"Single HP Cold-Start Outer Wall Temperature Every {args.interval:.0f} s"
    plot_wall_temperature(
        y_centers=y_centers,
        sampled_times=sampled_times,
        sampled_curves=sampled_curves,
        output_path=args.output_path,
        title=title,
    )
    write_sample_csv(args.csv_path, y_centers, sampled_times, sampled_curves)

    print(f"PNG written to: {args.output_path}")
    print(f"Sample CSV written to: {args.csv_path}")
    print(f"Sampled times: {sampled_times.tolist()}")


if __name__ == "__main__":
    main()
