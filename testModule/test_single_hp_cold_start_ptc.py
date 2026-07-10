import unittest

import numpy as np

from Components.basicComponents.HeatPipe2D import HeatPipe2D
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.Mesh import Mesh2D


class TestSodiumWickPTC(unittest.TestCase):
    @staticmethod
    def _build_hp(
        n_wick: int = 2,
        n_wall: int = 2,
        structural_only: bool = False,
        enable_ptc: bool = False,
    ) -> HeatPipe2D:
        L_total = 0.60
        L_eva = 0.06
        L_aba = 0.0
        L_con = 0.54
        r_vapor = 8.5e-3
        r_in_wall = 9.0e-3
        r_out_wall = 11.0e-3
        porosity = 0.675
        n_eva = 3
        n_aba = 0
        n_con = 18

        n_y = n_eva + n_aba + n_con
        if n_y != 21:
            raise AssertionError("Unexpected axial section count for cold-start fixture.")

        x_faces = np.array(
            [r_vapor, 8.75e-3, r_in_wall, 0.010, r_out_wall], dtype=float
        )
        y_faces = np.concatenate(
            (
                np.linspace(0.0, L_eva, n_eva + 1),
                L_eva + np.linspace(0.0, L_con, n_con + 1)[1:],
            )
        )

        mesh = Mesh2D(
            x_dim=r_out_wall - r_vapor,
            n_x=n_wick + n_wall,
            y_dim=L_total,
            n_y=n_y,
            geometry_type="cylindrical",
            inner_radius=r_vapor,
            x_faces=x_faces,
            y_faces=y_faces,
        )

        mat_wall = SS316(name="PTC_SS316")
        mat_fluid = SodiumHP(name="PTC_SodiumHP")

        hp = HeatPipe2D(
            mesh=mesh,
            solid1=mat_wall,
            solid2=mat_fluid,
            solid3=SS316(),
            n_wick=n_wick,
            porosity=porosity,
            n_eva=n_eva,
            n_aba=n_aba,
            n_con=n_con,
            name="PTC_Single_HP",
            initial_temp=300.0,
        )

        hp.boundaries["outer_eva"].clear_conditions()
        hp.boundaries["outer_con"].clear_conditions()
        hp.boundaries["outer_aba"].clear_conditions()

        area_eva = np.asarray(hp.boundaries["outer_eva"].area, dtype=float)
        total_power = 150.0
        q_eva = total_power * area_eva / np.sum(area_eva)
        hp.boundaries["outer_eva"].add_flux_condition(q_flux=q_eva)

        outer_con = hp.boundaries["outer_con"]
        outer_con.add_dynamic_radiation_condition(
            emissivity=0.03,
            bare_area_array=np.array(outer_con.area, dtype=float),
            T_env=300.0,
        )

        hp.set_time_integrator("theta_implicit")
        hp.set_theta_implicit_value(0.7)
        hp.set_face_conductance_mode("resistance_split_full")
        hp.set_wick_conductivity_mode(enable_ptc)
        hp.enable_frozen_property_correction = True
        hp.max_outer_property_corrections = 3
        hp.outer_property_tol = 1.0e-4

        if structural_only:
            # Keep explicit structural-only axial conductivity while preserving
            # wall/wick geometric/material interfaces from HeatPipe2D construction.
            hp.wick_mat.conductivity_axial = hp.wick_mat.conductivity_structural
            hp._property_cache_initialized = False
            hp._wick_temperature_cache[:] = np.nan
            hp._wall_temperature_cache[:] = np.nan
            hp._update_properties()

        return hp

    def _run_transient(
        self,
        hp: HeatPipe2D,
        *,
        t_end: float = 5.0,
        dt: float = 0.05,
        stop_outer_con_delta: float | None = None,
    ):
        n_steps = int(round(t_end / dt))
        n_eva = hp.n_eva
        n_con_start = hp.n_eva + hp.n_aba

        initial_2d = hp.T.reshape(hp.shape_nodes)
        initial = {
            "T_outer_eva": float(np.mean(initial_2d[-1, :n_eva])),
            "T_outer_con": float(np.mean(initial_2d[-1, n_con_start:])),
            "T_wick_con": float(np.mean(initial_2d[:hp.n_wick, n_con_start:])),
            "delta_eva_to_con": float(
                np.mean(initial_2d[-1, :n_eva]) - np.mean(initial_2d[-1, n_con_start:])
            ),
        }

        elapsed_time = 0.0
        threshold_reached = stop_outer_con_delta is None
        for step_index in range(n_steps):
            ok = hp.step(dt)
            self.assertTrue(isinstance(ok, bool), "hp.step() did not return a boolean.")
            self.assertTrue(ok, "Transient step failed.")

            elapsed_time = (step_index + 1) * dt
            T_2d = hp.T.reshape(hp.shape_nodes)
            self.assertTrue(np.all(np.isfinite(T_2d)), "Non-finite temperature detected.")
            self.assertGreaterEqual(np.nanmin(T_2d), 0.0)

            if stop_outer_con_delta is not None:
                outer_con_delta = float(np.mean(T_2d[-1, n_con_start:])) - initial["T_outer_con"]
                if outer_con_delta >= stop_outer_con_delta:
                    threshold_reached = True
                    break

        final_2d = hp.T.reshape(hp.shape_nodes)
        final = {
            "T_outer_eva": float(np.mean(final_2d[-1, :n_eva])),
            "T_outer_con": float(np.mean(final_2d[-1, n_con_start:])),
            "T_wick_con": float(np.mean(final_2d[:hp.n_wick, n_con_start:])),
            "delta_eva_to_con": float(
                np.mean(final_2d[-1, :n_eva]) - np.mean(final_2d[-1, n_con_start:])
            ),
            "delta_outer_con": 0.0,
            "delta_wick_con": 0.0,
            "elapsed_time": elapsed_time,
        }
        final["delta_outer_con"] = final["T_outer_con"] - initial["T_outer_con"]
        final["delta_wick_con"] = final["T_wick_con"] - initial["T_wick_con"]

        if not threshold_reached:
            self.fail(
                f"Condenser outer wall warmed by {final['delta_outer_con']:.3f} K, "
                f"below required {stop_outer_con_delta:.3f} K after {t_end:.3f} s."
            )

        return {"initial": initial, "final": final}
    def test_wick_conductivity_decomposition(self):
        mat_wick = WickMaterial(
            name="PTC_Wick_Composite",
            solid_mat=SS316(),
            fluid_mat=SodiumHP(),
            porosity=0.675,
            r_vapor=8.5e-3,
            r_in_wall=9.0e-3,
        )

        temperatures = np.array([300.0, 360.0, 371.0, 450.0, 600.0, 800.0, 1000.0])

        k_axial = np.asarray(mat_wick.conductivity_axial(temperatures), dtype=float)
        k_radial = np.asarray(mat_wick.conductivity_radial(temperatures), dtype=float)
        k_structural = np.asarray(mat_wick.conductivity_structural(temperatures), dtype=float)
        k_pseudo = np.asarray(mat_wick.conductivity_pseudothermal(temperatures), dtype=float)

        self.assertTrue(np.all(np.isfinite(k_axial)))
        self.assertTrue(np.all(np.isfinite(k_radial)))
        self.assertTrue(np.all(k_axial >= 0.0))
        self.assertTrue(np.all(k_radial >= 0.0))

        self.assertTrue(np.allclose(k_axial, k_structural + k_pseudo, rtol=1e-10, atol=1e-8))
        self.assertGreater(float(k_axial[-1]), 10.0 * max(float(k_axial[0]), 1.0))

    def test_single_hp_cold_start_ptc_vs_structural_only(self):
        min_condenser_outer_delta = 20.0
        ptc_result = self._run_transient(
            self._build_hp(enable_ptc=True, structural_only=False),
            t_end=240.0,
            dt=0.05,
            stop_outer_con_delta=min_condenser_outer_delta,
        )
        ptc_elapsed_time = ptc_result["final"]["elapsed_time"]
        struct_result = self._run_transient(
            self._build_hp(enable_ptc=False, structural_only=True),
            t_end=ptc_elapsed_time,
            dt=0.05,
        )

        ptc_final = ptc_result["final"]
        struct_final = struct_result["final"]

        self.assertGreater(ptc_final["T_outer_eva"], ptc_final["T_outer_con"])
        self.assertGreater(struct_final["T_outer_eva"], struct_final["T_outer_con"])

        self.assertGreaterEqual(ptc_final["delta_outer_con"], min_condenser_outer_delta)
        self.assertLessEqual(ptc_final["elapsed_time"], 240.0)
        self.assertGreater(ptc_final["delta_outer_con"], struct_final["delta_outer_con"])
        self.assertGreater(ptc_final["delta_wick_con"], struct_final["delta_wick_con"])
        self.assertLess(ptc_final["delta_eva_to_con"], struct_final["delta_eva_to_con"])

        self.assertTrue(np.isfinite(ptc_final["T_outer_eva"]))
        self.assertTrue(np.isfinite(ptc_final["T_outer_con"]))
        self.assertTrue(np.isfinite(struct_final["T_outer_eva"]))
        self.assertTrue(np.isfinite(struct_final["T_outer_con"]))

        self.assertLess(
            ptc_final["delta_eva_to_con"],
            struct_final["delta_eva_to_con"] - 1.0,
        )
if __name__ == "__main__":
    unittest.main()
