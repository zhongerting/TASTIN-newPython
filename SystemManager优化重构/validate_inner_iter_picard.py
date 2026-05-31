import csv
import os
import sys

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from Solvers.SystemManager import SystemManager


class ValidationVolume:
    def __init__(self, temperature):
        self.T = float(temperature)
        self.Q_wall = 0.0
        self.Q_vol = 0.0
        self.implicit_coeff = 0.0

    def add_coupling_source(self, explicit_part, implicit_factor):
        self.Q_wall += float(explicit_part)
        self.implicit_coeff += float(implicit_factor)


class ValidationFluid:
    def __init__(self, volume, heat_capacity):
        self.volumes_obj = [volume]
        self.n_vol = 1
        self.n_junc = 0
        self.heat_capacity = float(heat_capacity)
        self.solve_log = []

    @property
    def T_vec(self):
        return np.array([self.volumes_obj[0].T], dtype=float)

    def initialize_hydraulics(self, dt=0.1, tol=1e-5, max_iter=500):
        return True

    def step_Picard(self, dt, max_iter=20, tol=1e-4):
        volume = self.volumes_obj[0]
        q_applied = volume.Q_wall + volume.Q_vol - volume.implicit_coeff * volume.T
        iteration = len(self.solve_log) + 1
        source_wall_T = (
            volume.Q_wall / volume.implicit_coeff
            if abs(volume.implicit_coeff) > 1.0e-30
            else np.nan
        )
        self.solve_log.append(
            {
                "iteration": iteration,
                "T_before": volume.T,
                "Q_wall": volume.Q_wall,
                "implicit_coeff": volume.implicit_coeff,
                "source_wall_T": source_wall_T,
                "q_applied": q_applied,
            }
        )
        volume.T += dt * q_applied / self.heat_capacity
        return True

    def get_max_stable_dt(self, max_limit=0.5):
        return max_limit

    def save_state(self):
        self._backup_T = self.volumes_obj[0].T

    def load_state(self):
        self.volumes_obj[0].T = self._backup_T


class ValidationBoundary:
    def __init__(self, solid):
        self.solid = solid
        self.T_surface = np.array([solid.T[0]], dtype=float)
        self.current_flux = np.zeros(1, dtype=float)

    def compute_net_flux_for_solver(self):
        self.current_flux[:] = self.solid.boundary_G * (
            self.solid.boundary_T_ext - self.solid.T[0]
        )
        return self.current_flux


class ValidationSolid:
    def __init__(self, temperature, heat_capacity):
        self.name = "validation_solid"
        self.N = 1
        self.T = np.array([float(temperature)], dtype=float)
        self.current_time = 0.0
        self.heat_capacity = float(heat_capacity)
        self.boundary_T_ext = float(temperature)
        self.boundary_G = 0.0
        self.boundaries = {"left": ValidationBoundary(self)}
        self.solve_log = []

    def save_state(self):
        self._backup_T = self.T.copy()
        self._backup_time = self.current_time

    def load_state(self):
        self.T[:] = self._backup_T
        self.current_time = self._backup_time
        self._update_boundaries_state(current_time=self.current_time)

    def step(self, dt):
        q_applied = self.boundary_G * (self.boundary_T_ext - self.T[0])
        iteration = len(self.solve_log) + 1
        self.solve_log.append(
            {
                "iteration": iteration,
                "T_before": float(self.T[0]),
                "T_ext": self.boundary_T_ext,
                "G": self.boundary_G,
                "q_applied": q_applied,
            }
        )
        self.T[0] += dt * q_applied / self.heat_capacity
        self.current_time += dt
        self._update_boundaries_state(current_time=self.current_time)
        return True

    def _update_properties(self):
        pass

    def _compute_internal_resistance(self):
        pass

    def _update_boundaries_state(self, current_time=None):
        self.boundaries["left"].T_surface[:] = self.T[0]

    def _compute_fluxes(self, current_time):
        return self.boundaries["left"].compute_net_flux_for_solver()


class ValidationCoupler:
    def __init__(self, fluid, solid, conductance):
        self.name = "validation_interface"
        self.fluid = fluid
        self.solid = solid
        self.G = float(conductance)
        self.execute_log = []

    def execute(self):
        volume = self.fluid.volumes_obj[0]
        T_f = float(volume.T)
        T_s = float(self.solid.T[0])
        q_to_fluid_at_coupler_state = self.G * (T_s - T_f)
        q_to_solid_at_coupler_state = self.G * (T_f - T_s)

        legacy_avg_q = np.nan
        if self.execute_log:
            legacy_avg_q = 0.5 * (
                self.execute_log[0]["q_to_fluid_at_coupler_state"]
                + q_to_fluid_at_coupler_state
            )

        volume.add_coupling_source(self.G * T_s, self.G)
        self.solid.boundary_T_ext = T_f
        self.solid.boundary_G = self.G

        self.execute_log.append(
            {
                "iteration": len(self.execute_log) + 1,
                "T_fluid_coupler": T_f,
                "T_solid_coupler": T_s,
                "q_to_fluid_at_coupler_state": q_to_fluid_at_coupler_state,
                "q_to_solid_at_coupler_state": q_to_solid_at_coupler_state,
                "legacy_avg_q": legacy_avg_q,
            }
        )


def run_validation():
    volume = ValidationVolume(temperature=600.0)
    fluid = ValidationFluid(volume=volume, heat_capacity=1000.0)
    solid = ValidationSolid(temperature=650.0, heat_capacity=1500.0)
    coupler = ValidationCoupler(fluid=fluid, solid=solid, conductance=100.0)

    manager = SystemManager(fluid_network=fluid)
    manager.solid_components[solid.name] = solid
    manager.add_coupler(coupler)
    manager.step(dt=0.5, inner_iter=4, convergence_tol=0.0)

    rows = []
    for c_row, f_row, s_row in zip(coupler.execute_log, fluid.solve_log, solid.solve_log):
        rows.append(
            {
                **c_row,
                "fluid_T_before_solve": f_row["T_before"],
                "fluid_Q_wall": f_row["Q_wall"],
                "fluid_implicit_coeff": f_row["implicit_coeff"],
                "fluid_source_wall_T": f_row["source_wall_T"],
                "fluid_q_applied": f_row["q_applied"],
                "solid_T_before_solve": s_row["T_before"],
                "solid_T_ext": s_row["T_ext"],
                "solid_G": s_row["G"],
                "solid_q_applied": s_row["q_applied"],
            }
        )

    csv_path = os.path.join(CURRENT_DIR, "inner_iter_picard_validation.csv")
    png_path = os.path.join(CURRENT_DIR, "inner_iter_picard_validation.png")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    iterations = np.array([row["iteration"] for row in rows], dtype=int)
    q_state = np.array([row["q_to_fluid_at_coupler_state"] for row in rows], dtype=float)
    q_state_solid = np.array([row["q_to_solid_at_coupler_state"] for row in rows], dtype=float)
    q_fluid_applied = np.array([row["fluid_q_applied"] for row in rows], dtype=float)
    q_solid_applied = np.array([row["solid_q_applied"] for row in rows], dtype=float)
    legacy_avg = np.array([row["legacy_avg_q"] for row in rows], dtype=float)
    implicit_coeff = np.array([row["fluid_implicit_coeff"] for row in rows], dtype=float)

    plt.figure(figsize=(11, 8))

    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(iterations, [row["T_fluid_coupler"] for row in rows], "o-", label="fluid T used by solid BC")
    ax1.plot(iterations, [row["T_solid_coupler"] for row in rows], "s-", label="solid T used by fluid source")
    ax1.set_xlabel("Picard iteration")
    ax1.set_ylabel("Temperature [K]")
    ax1.set_title("Coupler state written to both sides")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(iterations, q_state, "o-", label="coupler q to fluid")
    ax2.plot(iterations, q_state_solid, "s-", label="coupler q to solid")
    ax2.plot(iterations, q_state + q_state_solid, "k--", label="coupler balance residual")
    ax2.set_xlabel("Picard iteration")
    ax2.set_ylabel("Heat rate [W]")
    ax2.set_title("Same-state interface pair")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(iterations, q_fluid_applied, "o-", label="new fluid source after rollback")
    ax3.plot(iterations, q_solid_applied, "s-", label="new solid BC source after rollback")
    ax3.plot(iterations, legacy_avg, "x--", label="old one-sided fluid avg (reference)")
    ax3.set_xlabel("Picard iteration")
    ax3.set_ylabel("Heat rate [W]")
    ax3.set_title("No _fluid_total_Q_backup / avg_Q rewrite")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(iterations, implicit_coeff, "o-", label="fluid implicit coeff")
    ax4.plot(iterations, [row["solid_G"] for row in rows], "s-", label="solid boundary conductance")
    ax4.set_xlabel("Picard iteration")
    ax4.set_ylabel("Conductance [W/K]")
    ax4.set_title("Both sides keep current coupler conductance")
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    plt.tight_layout()
    plt.savefig(png_path, dpi=160)
    plt.close()

    max_state_balance = float(np.max(np.abs(q_state + q_state_solid)))
    avg_reference_delta = float(np.nanmax(np.abs(q_fluid_applied[1:] - legacy_avg[1:])))
    min_implicit = float(np.min(implicit_coeff))

    print(f"CSV: {csv_path}")
    print(f"PNG: {png_path}")
    print(f"iterations: {len(rows)}")
    print(f"max coupler-state balance residual: {max_state_balance:.6e} W")
    print(f"min fluid implicit coeff during solves: {min_implicit:.6e} W/K")
    print(f"max delta vs old one-sided avg reference: {avg_reference_delta:.6e} W")
    print(f"final fluid T: {volume.T:.6f} K")
    print(f"final solid T: {solid.T[0]:.6f} K")


if __name__ == "__main__":
    run_validation()
