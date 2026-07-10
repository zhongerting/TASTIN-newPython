"""Single HeatPipe2D cold-start run with an 800 K resistance heater.

This script is intentionally independent from unittest. It creates a refined
single-heat-pipe fixture, runs a transient, and writes outer-wall axial-node
temperature histories plus plots.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Components.basicComponents.HeatPipe2D import HeatPipe2D
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Solvers.HeatConduction.Mesh import Mesh2D


def build_single_heat_pipe(
    *,
    heater_temperature_k: float,
    heater_ha_w_per_k: float,
    condenser_emissivity: float,
    condenser_background_k: float,
    initial_temperature_k: float,
    n_eva: int = 15,
    n_con: int = 135,
) -> HeatPipe2D:
    """Build a refined sodium heat pipe for the 800 K heater startup case."""

    l_total = 0.60
    l_eva = 0.06
    l_aba = 0.0
    l_con = 0.54
    r_vapor = 8.5e-3
    r_in_wall = 9.0e-3
    r_out_wall = 11.0e-3
    porosity = 0.675
    n_wick = 2
    n_wall = 2
    n_aba = 0

    if n_eva <= 0 or n_con <= 0:
        raise ValueError("n_eva and n_con must be positive.")

    x_faces = np.array([r_vapor, 8.75e-3, r_in_wall, 0.010, r_out_wall], dtype=float)
    y_faces = np.concatenate(
        (
            np.linspace(0.0, l_eva, n_eva + 1),
            l_eva + np.linspace(0.0, l_con, n_con + 1)[1:],
        )
    )

    hp = HeatPipe2D(
        mesh=Mesh2D(
            x_dim=r_out_wall - r_vapor,
            n_x=n_wick + n_wall,
            y_dim=l_total,
            n_y=n_eva + n_aba + n_con,
            geometry_type="cylindrical",
            inner_radius=r_vapor,
            x_faces=x_faces,
            y_faces=y_faces,
        ),
        solid1=SS316(name="Run800K_SS316"),
        solid2=SodiumHP(name="Run800K_SodiumHP"),
        solid3=SS316(),
        n_wick=n_wick,
        porosity=porosity,
        n_eva=n_eva,
        n_aba=n_aba,
        n_con=n_con,
        name="Single_HP_800K_Heater",
        initial_temp=initial_temperature_k,
    )

    hp.boundaries["outer_eva"].clear_conditions()
    hp.boundaries["outer_con"].clear_conditions()
    hp.boundaries["outer_aba"].clear_conditions()

    if heater_ha_w_per_k <= 0.0:
        raise ValueError("heater_ha_w_per_k must be positive.")
    hp.boundaries["outer_eva"].add_resistance_condition(
        T_ext=heater_temperature_k,
        R_ext=1.0 / heater_ha_w_per_k,
    )

    outer_con = hp.boundaries["outer_con"]
    outer_con.add_dynamic_radiation_condition(
        emissivity=condenser_emissivity,
        bare_area_array=np.array(outer_con.area, dtype=float),
        T_env=condenser_background_k,
    )

    hp.set_time_integrator("theta_implicit")
    hp.set_theta_implicit_value(0.7)
    hp.set_face_conductance_mode("resistance_split_full")
    hp.set_wick_conductivity_mode(True)
    hp.enable_frozen_property_correction = True
    hp.max_outer_property_corrections = 3
    hp.outer_property_tol = 1.0e-4
    hp.set_implicit_boundary_linearization(True)

    return hp


def representative_nodes(hp: HeatPipe2D) -> dict[str, tuple[int, int]]:
    """Return a few auxiliary nodes for compact summaries."""

    n_con_start = hp.n_eva + hp.n_aba
    return {
        "wick_eva_mid": (0, hp.n_eva // 2),
        "wall_eva_outer_mid": (hp.shape_nodes[0] - 1, hp.n_eva // 2),
        "wall_con_start": (hp.shape_nodes[0] - 1, n_con_start),
        "wick_con_start": (0, n_con_start),
        "wall_con_mid": (hp.shape_nodes[0] - 1, n_con_start + hp.n_con // 2),
        "wick_con_mid": (0, n_con_start + hp.n_con // 2),
        "wall_con_end": (hp.shape_nodes[0] - 1, hp.shape_nodes[1] - 1),
        "wick_con_end": (0, hp.shape_nodes[1] - 1),
    }


def outer_wall_temperature_columns(hp: HeatPipe2D) -> list[str]:
    return [
        f"wall_outer_z{j:03d}_m{float(z):.5f}"
        for j, z in enumerate(np.asarray(hp.mesh.y_centers, dtype=float))
    ]


def sample_temperatures(hp: HeatPipe2D, nodes: dict[str, tuple[int, int]]) -> dict[str, float]:
    t_2d = hp.T.reshape(hp.shape_nodes)
    row = {name: float(t_2d[i, j]) for name, (i, j) in nodes.items()}
    row.update(
        {
            name: float(t_2d[-1, j])
            for j, name in enumerate(outer_wall_temperature_columns(hp))
        }
    )

    n_con_start = hp.n_eva + hp.n_aba
    row["mean_outer_eva"] = float(np.mean(t_2d[-1, : hp.n_eva]))
    row["mean_outer_con"] = float(np.mean(t_2d[-1, n_con_start:]))
    row["max_temperature"] = float(np.max(t_2d))
    row["min_temperature"] = float(np.min(t_2d))
    return row


def write_csv(csv_path: Path, rows: list[dict[str, float]]) -> None:
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_outer_wall_plot(plot_path: Path, rows: list[dict[str, float]], hp: HeatPipe2D) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    time_s = np.array([row["time_s"] for row in rows], dtype=float)
    z_centers = np.asarray(hp.mesh.y_centers, dtype=float)
    wall_columns = outer_wall_temperature_columns(hp)

    fig, ax = plt.subplots(figsize=(12, 7), dpi=150)
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=float(z_centers[0]), vmax=float(z_centers[-1]))
    for z, name in zip(z_centers, wall_columns):
        ax.plot(
            time_s,
            [row[name] for row in rows],
            color=cmap(norm(z)),
            linewidth=0.8,
            alpha=0.75,
        )
    ax.plot(time_s, [row["mean_outer_eva"] for row in rows], "k--", linewidth=2.0, label="mean_outer_eva")
    ax.plot(time_s, [row["mean_outer_con"] for row in rows], "k:", linewidth=2.0, label="mean_outer_con")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Outer Wall Axial Node Temperatures")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Axial position z [m]")
    fig.tight_layout()
    fig.savefig(plot_path)
    plt.close(fig)


def write_outer_wall_heatmap(plot_path: Path, rows: list[dict[str, float]], hp: HeatPipe2D) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    time_s = np.array([row["time_s"] for row in rows], dtype=float)
    z_centers = np.asarray(hp.mesh.y_centers, dtype=float)
    wall_columns = outer_wall_temperature_columns(hp)
    temperatures = np.array([[row[name] for name in wall_columns] for row in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
    mesh = ax.pcolormesh(time_s, z_centers, temperatures.T, shading="auto", cmap="inferno")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Axial position z [m]")
    ax.set_title("Outer Wall Temperature Field")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Temperature [K]")
    fig.tight_layout()
    fig.savefig(plot_path)
    plt.close(fig)


def run_case(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    hp = build_single_heat_pipe(
        heater_temperature_k=args.heater_temperature_k,
        heater_ha_w_per_k=args.heater_ha_w_per_k,
        condenser_emissivity=args.condenser_emissivity,
        condenser_background_k=args.condenser_background_k,
        initial_temperature_k=args.initial_temperature_k,
        n_eva=args.n_eva,
        n_con=args.n_con,
    )
    nodes = representative_nodes(hp)

    rows: list[dict[str, float]] = []
    n_steps = int(round(args.duration_s / args.dt_s))
    record_every = max(1, int(round(args.record_interval_s / args.dt_s)))

    row0 = {"time_s": 0.0}
    row0.update(sample_temperatures(hp, nodes))
    rows.append(row0)

    for step_index in range(1, n_steps + 1):
        ok = hp.step(args.dt_s)
        if not ok:
            raise RuntimeError(f"HeatPipe2D step failed at t={step_index * args.dt_s:.6g} s")

        t_2d = hp.T.reshape(hp.shape_nodes)
        if not np.all(np.isfinite(t_2d)):
            raise RuntimeError(f"Non-finite temperature at t={step_index * args.dt_s:.6g} s")
        if float(np.min(t_2d)) < 0.0:
            raise RuntimeError(f"Negative temperature at t={step_index * args.dt_s:.6g} s")

        if step_index % record_every == 0 or step_index == n_steps:
            row = {"time_s": step_index * args.dt_s}
            row.update(sample_temperatures(hp, nodes))
            rows.append(row)

    csv_path = output_dir / "single_hp_cold_start_800k_history.csv"
    outer_wall_plot_path = output_dir / "single_hp_cold_start_800k_outer_wall_axial_temperatures.png"
    heatmap_path = output_dir / "single_hp_cold_start_800k_outer_wall_heatmap.png"
    summary_path = output_dir / "single_hp_cold_start_800k_summary.json"
    write_csv(csv_path, rows)
    write_outer_wall_plot(outer_wall_plot_path, rows, hp)
    write_outer_wall_heatmap(heatmap_path, rows, hp)

    initial_heat_w = (args.heater_temperature_k - args.initial_temperature_k) * args.heater_ha_w_per_k
    summary = {
        "duration_s": args.duration_s,
        "dt_s": args.dt_s,
        "record_interval_s": args.record_interval_s,
        "heater_temperature_k": args.heater_temperature_k,
        "heater_ha_w_per_k": args.heater_ha_w_per_k,
        "heater_resistance_k_per_w": 1.0 / args.heater_ha_w_per_k,
        "initial_heater_power_w": initial_heat_w,
        "condenser_emissivity": args.condenser_emissivity,
        "condenser_background_k": args.condenser_background_k,
        "initial_temperature_k": args.initial_temperature_k,
        "n_eva": args.n_eva,
        "n_aba": hp.n_aba,
        "n_con": args.n_con,
        "n_axial": hp.shape_nodes[1],
        "history_csv": str(csv_path),
        "plot_png": str(outer_wall_plot_path),
        "outer_wall_plot_png": str(outer_wall_plot_path),
        "outer_wall_heatmap_png": str(heatmap_path),
        "final": rows[-1],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary["summary_json"] = str(summary_path)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration-s", type=float, default=1800.0)
    parser.add_argument("--dt-s", type=float, default=0.1)
    parser.add_argument("--record-interval-s", type=float, default=5.0)
    parser.add_argument("--heater-temperature-k", type=float, default=800.0)
    parser.add_argument("--heater-ha-w-per-k", type=float, default=0.3)
    parser.add_argument("--condenser-emissivity", type=float, default=0.03)
    parser.add_argument("--condenser-background-k", type=float, default=4.0)
    parser.add_argument("--initial-temperature-k", type=float, default=300.0)
    parser.add_argument("--n-eva", type=int, default=15)
    parser.add_argument("--n-con", type=int, default=135)
    parser.add_argument(
        "--output-dir",
        default="testModule/single_hp_cold_start_800k_1800s_refined_wall",
    )
    return parser.parse_args()


def main() -> None:
    summary = run_case(parse_args())
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
