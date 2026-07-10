import argparse
import json
import time
import os
import sys
from typing import Dict, Any, Callable

import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from Materials.Base import SolidMaterial
from Solvers.HeatConduction.HeatConduction import HeatConduction1D, HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh1D, Mesh2D


class ConstantSolidMaterial(SolidMaterial):
    def __init__(self, name: str, k: float, rho: float, cp: float):
        super().__init__(name=name, formula="constant")
        self._k = float(k)
        self._rho = float(rho)
        self._cp = float(cp)

    def _broadcast(self, value: float, T) -> np.ndarray:
        arr = np.asarray(T, dtype=float)
        if arr.shape == ():
            return float(value)
        return np.full_like(arr, value, dtype=float)

    def conductivity(self, T):
        return self._broadcast(self._k, T)

    def density(self, T):
        return self._broadcast(self._rho, T)

    def heat_capacity(self, T):
        return self._broadcast(self._cp, T)


def build_1d_case(
    n_nodes: int = 4000,
    length_m: float = 0.50,
    source_w: float = 1500.0,
) -> HeatConduction1D:
    mesh = Mesh1D(total_dim=length_m, n_volumes=n_nodes, geometry_type="cartesian", height=1.0)
    solid = HeatConduction1D(
        mesh=mesh,
        material=ConstantSolidMaterial("bench_const", k=18.0, rho=7800.0, cp=500.0),
        initial_temp=600.0,
    )

    for bc in solid.boundaries.values():
        bc.clear_conditions()

    solid.boundaries["inner"].add_resistance_condition(T_ext=500.0, R_ext=0.12)
    solid.boundaries["outer"].add_convection_condition(T_fluid=320.0, h_coeff=60.0)
    solid.link_source_buffer(np.full(solid.N, source_w, dtype=float))
    solid.set_ode_method("implicit_euler")

    return solid


def build_2d_case(
    n_x: int = 90,
    n_y: int = 90,
    x_dim_m: float = 0.05,
    y_dim_m: float = 0.20,
    source_w: float = 900.0,
) -> HeatConduction2D:
    mesh = Mesh2D(
        x_dim=x_dim_m,
        n_x=n_x,
        y_dim=y_dim_m,
        n_y=n_y,
        geometry_type="cartesian",
    )
    solid = HeatConduction2D(
        mesh=mesh,
        material=ConstantSolidMaterial("bench_const", k=18.0, rho=7800.0, cp=500.0),
        initial_temp=600.0,
    )

    for bc in solid.boundaries.values():
        bc.clear_conditions()

    solid.boundaries["left"].add_resistance_condition(T_ext=500.0, R_ext=0.08)
    solid.boundaries["right"].add_convection_condition(T_fluid=320.0, h_coeff=55.0)
    solid.boundaries["top"].add_resistance_condition(T_ext=400.0, R_ext=0.20)
    solid.boundaries["bottom"].add_resistance_condition(T_ext=620.0, R_ext=0.20)
    solid.link_source_buffer(np.full(solid.N, source_w, dtype=float))
    solid.set_ode_method("implicit_euler")

    return solid


def run_case(name: str, builder: Callable[[], Any], dt: float, n_steps: int) -> Dict[str, float]:
    solid = builder()
    t0 = time.perf_counter()
    for _ in range(n_steps):
        ok = solid.step(dt)
        if not ok:
            raise RuntimeError(f"{name}: implicit_euler step failed at configured n_steps={n_steps}")
    wall_s = time.perf_counter() - t0
    temperature = np.asarray(solid.T, dtype=float)
    return {
        "name": name,
        "ode_method": solid.ode_method,
        "n_nodes": int(solid.N),
        "dt_s": float(dt),
        "steps": int(n_steps),
        "wall_time_s": float(wall_s),
        "wall_time_per_step_s": float(wall_s / max(int(n_steps), 1)),
        "temp_min_k": float(temperature.min()),
        "temp_max_k": float(temperature.max()),
        "temp_mean_k": float(temperature.mean()),
        "temp_sum_k": float(temperature.sum()),
        "q_source_sum_w": float(np.sum(solid.Q_source)),
    }


def run_baseline_case(name: str, builder: Callable[[], Any], dt: float, n_steps: int, method: str) -> Dict[str, float]:
    solid = builder()
    if method not in HeatConduction1D.VALID_ODE_METHODS:
        raise ValueError(f"Unsupported ODE method '{method}'")
    solid.set_ode_method(method)
    t0 = time.perf_counter()
    for _ in range(n_steps):
        ok = solid.step(dt)
        if not ok:
            raise RuntimeError(f"{name}: baseline method {method} failed at configured n_steps={n_steps}")
    wall_s = time.perf_counter() - t0
    temperature = np.asarray(solid.T, dtype=float)
    return {
        "name": name,
        "ode_method": solid.ode_method,
        "n_nodes": int(solid.N),
        "dt_s": float(dt),
        "steps": int(n_steps),
        "wall_time_s": float(wall_s),
        "wall_time_per_step_s": float(wall_s / max(int(n_steps), 1)),
        "temp_min_k": float(temperature.min()),
        "temp_max_k": float(temperature.max()),
        "temp_mean_k": float(temperature.mean()),
        "temp_sum_k": float(temperature.sum()),
        "q_source_sum_w": float(np.sum(solid.Q_source)),
    }


def print_case_report(result: Dict[str, float]) -> None:
    print(
        f"{result['name']:>18s} | {result['ode_method']:<12s}"
        f" | nodes={result['n_nodes']:>7d}"
        f" | dt={result['dt_s']:.4f}s"
        f" | steps={result['steps']:>4d}"
        f" | wall={result['wall_time_s']:.4f}s"
        f" | wall/step={result['wall_time_per_step_s']:.6f}s"
        f" | Tmin={result['temp_min_k']:.4f}"
        f" | Tmax={result['temp_max_k']:.4f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Micro-benchmark implicit_euler heat conduction stepping for fixed 1D/2D configurations."
    )
    parser.add_argument("--steps-1d", type=int, default=40)
    parser.add_argument("--steps-2d", type=int, default=20)
    parser.add_argument("--dt-1d", type=float, default=0.02)
    parser.add_argument("--dt-2d", type=float, default=0.02)
    parser.add_argument("--n1d", type=int, default=4000)
    parser.add_argument("--n2d-x", type=int, default=90)
    parser.add_argument("--n2d-y", type=int, default=90)
    parser.add_argument("--baseline", action="store_true")
    parser.add_argument("--baseline-method", default="BDF", choices=("BDF", "RK45"))
    parser.add_argument("--summary-json", default="")
    args = parser.parse_args()

    if args.n1d < 8:
        raise ValueError("n1d must be >= 8")
    if args.n2d_x < 4 or args.n2d_y < 4:
        raise ValueError("n2d_x and n2d_y must be >= 4")

    results = []
    now_stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"HeatConduction implicit_euler micro-benchmark @ {now_stamp}")
    print(
        f"Config: 1D n={args.n1d}, dt={args.dt_1d}, steps={args.steps_1d}; "
        f"2D nx={args.n2d_x}, ny={args.n2d_y}, dt={args.dt_2d}, steps={args.steps_2d}"
    )

    case_1d = run_case(
        "1D implicit_euler",
        lambda: build_1d_case(n_nodes=args.n1d),
        dt=args.dt_1d,
        n_steps=args.steps_1d,
    )
    print_case_report(case_1d)
    results.append(case_1d)

    case_2d = run_case(
        "2D implicit_euler",
        lambda: build_2d_case(n_x=args.n2d_x, n_y=args.n2d_y),
        dt=args.dt_2d,
        n_steps=args.steps_2d,
    )
    print_case_report(case_2d)
    results.append(case_2d)

    if args.baseline:
        base_1d = run_baseline_case(
            "1D baseline",
            lambda: build_1d_case(n_nodes=args.n1d),
            dt=args.dt_1d,
            n_steps=args.steps_1d,
            method=args.baseline_method,
        )
        print_case_report(base_1d)
        results.append(base_1d)

        base_2d = run_baseline_case(
            "2D baseline",
            lambda: build_2d_case(n_x=args.n2d_x, n_y=args.n2d_y),
            dt=args.dt_2d,
            n_steps=args.steps_2d,
            method=args.baseline_method,
        )
        print_case_report(base_2d)
        results.append(base_2d)

    if args.summary_json:
        with open(args.summary_json, "w", encoding="utf-8") as stream:
            json.dump(results, stream, ensure_ascii=False, indent=2)
        print(f"Summary written to {args.summary_json}")

    print("DONE")


if __name__ == "__main__":
    main()
