import argparse
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CURRENT_DIR = Path(__file__).resolve().parent
L_EVA = 0.265
L_ABA = 0.235
L_CON = 0.500


def load_case(npz_path: Path, label: str):
    data = np.load(npz_path)

    y_centers = np.asarray(data["y_centers"], dtype=float)
    outer_wall_temperature = np.asarray(data["outer_wall_temperature_history"][-1], dtype=float)
    wick_k = np.asarray(data["wick_k_history"][-1, 0, :], dtype=float)

    n_eva = int(np.asarray(data["n_eva"], dtype=int)[0])
    n_aba = int(np.asarray(data["n_aba"], dtype=int)[0])
    n_con = int(np.asarray(data["n_con"], dtype=int)[0])

    q_eva = np.asarray(data["q_eva_input_distribution_history"][-1], dtype=float)
    q_aba = np.asarray(data["q_aba_loss_distribution_history"][-1], dtype=float)
    q_con = np.asarray(data["q_con_rejection_distribution_history"][-1], dtype=float)

    dy_eva = L_EVA / n_eva
    dy_aba = L_ABA / n_aba
    dy_con = L_CON / n_con
    q_linear = np.concatenate([
        q_eva / dy_eva,
        -q_aba / dy_aba,
        -q_con / dy_con,
    ])

    return {
        "label": label,
        "npz_path": str(npz_path),
        "y_centers": y_centers,
        "outer_wall_temperature": outer_wall_temperature,
        "wick_k": wick_k,
        "q_linear": q_linear,
        "n_eva": n_eva,
        "n_aba": n_aba,
        "n_con": n_con,
    }


def plot_curves(cases, key: str, ylabel: str, title: str, output_path: Path, *, yscale: str = "linear"):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=180)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(cases)))

    for color, case in zip(colors, cases):
        ax.plot(
            case["y_centers"] * 1000.0,
            case[key],
            color=color,
            linewidth=1.8,
            marker="o",
            markersize=2.8,
            label=(
                f"{case['label']} "
                f"({case['n_eva']}/{case['n_aba']}/{case['n_con']})"
            ),
        )

    for x_mm in (L_EVA * 1000.0, (L_EVA + L_ABA) * 1000.0):
        ax.axvline(x_mm, color="gray", linewidth=1.0, linestyle="--", alpha=0.5)

    ax.set_xlabel("Axial Position [mm]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_yscale(yscale)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, ncol=2, frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def write_csv(csv_path: Path, cases, key: str, value_name: str):
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "n_eva", "n_aba", "n_con", "axial_position_m", value_name])
        for case in cases:
            for y, value in zip(case["y_centers"], case[key]):
                writer.writerow([
                    case["label"],
                    case["n_eva"],
                    case["n_aba"],
                    case["n_con"],
                    float(y),
                    float(value),
                ])


def parse_case_specs(case_items):
    if not case_items:
        raise ValueError("At least one --case LABEL=PATH item is required.")

    case_specs = []
    for item in case_items:
        if "=" not in item:
            raise ValueError(f"Invalid --case item '{item}'. Expected LABEL=PATH.")
        label, raw_path = item.split("=", 1)
        label = label.strip()
        raw_path = raw_path.strip()
        if not label or not raw_path:
            raise ValueError(f"Invalid --case item '{item}'. Expected LABEL=PATH.")
        case_specs.append((label, Path(raw_path)))
    return case_specs


def parse_args():
    parser = argparse.ArgumentParser(description="Plot steady-state grid-comparison results for single HP startup.")
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help="Case spec in the form LABEL=PATH_TO_NPZ. Repeat for multiple cases.",
    )
    parser.add_argument(
        "--title-prefix",
        default="Steady-State",
        help="Prefix added to all figure titles.",
    )
    parser.add_argument(
        "--temperature-png",
        default=str(CURRENT_DIR / "single_hp_power_convection_steady_temperature_grid_comparison.png"),
    )
    parser.add_argument(
        "--heat-png",
        default=str(CURRENT_DIR / "single_hp_power_convection_steady_axial_heat_grid_comparison.png"),
    )
    parser.add_argument(
        "--wick-k-png",
        default=str(CURRENT_DIR / "single_hp_power_convection_steady_wick_k_grid_comparison.png"),
    )
    parser.add_argument(
        "--temperature-csv",
        default=str(CURRENT_DIR / "single_hp_power_convection_steady_temperature_grid_comparison.csv"),
    )
    parser.add_argument(
        "--heat-csv",
        default=str(CURRENT_DIR / "single_hp_power_convection_steady_axial_heat_grid_comparison.csv"),
    )
    parser.add_argument(
        "--wick-k-csv",
        default=str(CURRENT_DIR / "single_hp_power_convection_steady_wick_k_grid_comparison.csv"),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    case_specs = parse_case_specs(args.case)

    cases = [load_case(path, label) for label, path in case_specs]

    plot_curves(
        cases=cases,
        key="outer_wall_temperature",
        ylabel="Outer Wall Temperature [K]",
        title=f"{args.title_prefix} Outer Wall Temperature Distribution",
        output_path=Path(args.temperature_png),
    )
    plot_curves(
        cases=cases,
        key="q_linear",
        ylabel="External Heat Exchange per Axial Length [W/m]",
        title=f"{args.title_prefix} Axial Heat Exchange Distribution",
        output_path=Path(args.heat_png),
    )
    plot_curves(
        cases=cases,
        key="wick_k",
        ylabel="Wick Effective Conductivity [W/m/K]",
        title=f"{args.title_prefix} Wick Conductivity Distribution",
        output_path=Path(args.wick_k_png),
        yscale="log",
    )

    write_csv(Path(args.temperature_csv), cases, "outer_wall_temperature", "outer_wall_temperature_K")
    write_csv(Path(args.heat_csv), cases, "q_linear", "external_heat_exchange_W_per_m")
    write_csv(Path(args.wick_k_csv), cases, "wick_k", "wick_effective_conductivity_W_per_m_K")

    print(f"Temperature PNG written to: {args.temperature_png}")
    print(f"Heat PNG written to: {args.heat_png}")
    print(f"Wick-k PNG written to: {args.wick_k_png}")
    print(f"Temperature CSV written to: {args.temperature_csv}")
    print(f"Heat CSV written to: {args.heat_csv}")
    print(f"Wick-k CSV written to: {args.wick_k_csv}")


if __name__ == "__main__":
    main()
