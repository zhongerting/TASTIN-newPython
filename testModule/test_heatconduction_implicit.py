import warnings
import unittest
from unittest.mock import patch

import numpy as np

from Materials.Base import SolidMaterial
from Solvers.HeatConduction.Boundary import DynamicRadiationResistanceBC
from Solvers.HeatConduction.HeatConduction import HeatConduction1D
from Solvers.HeatConduction.Mesh import Mesh1D


class ConstantSolid(SolidMaterial):
    """Stable material with constant properties for compact regression tests."""

    def __init__(self,
                 k: float = 12.0,
                 rho: float = 7800.0,
                 cp: float = 600.0):
        super().__init__(name="Constant_Solid", formula="const")
        self._k = float(k)
        self._rho = float(rho)
        self._cp = float(cp)

    def _broadcast_scalar(self, value: float, T) -> np.ndarray:
        arr = np.asarray(T, dtype=float)
        if arr.shape == ():
            return float(value)
        return np.full_like(arr, value, dtype=float)

    def conductivity(self, T):
        return self._broadcast_scalar(self._k, T)

    def density(self, T):
        return self._broadcast_scalar(self._rho, T)

    def heat_capacity(self, T):
        return self._broadcast_scalar(self._cp, T)


def _build_test_solid(*,
                      method: str = "BDF",
                      with_dynamic_radiation: bool = False,
                      initial_temp: float = 600.0,
                      n_nodes: int = 12,
                      total_length: float = 0.25) -> HeatConduction1D:
    material = ConstantSolid()
    mesh = Mesh1D(total_dim=total_length, n_volumes=n_nodes)
    solid = HeatConduction1D(mesh=mesh, material=material, initial_temp=initial_temp)

    # Remove implicit default BC then apply controlled conditions.
    for boundary in solid.boundaries.values():
        boundary.clear_conditions()

    if with_dynamic_radiation:
        solid.boundaries["inner"].add_dynamic_radiation_condition(
            emissivity=0.85,
            bare_area_array=solid.boundaries["inner"].area.copy(),
            T_env=600.0,
        )
    else:
        solid.boundaries["inner"].add_resistance_condition(T_ext=600.0, R_ext=0.12)

    solid.boundaries["outer"].add_convection_condition(T_fluid=300.0, h_coeff=35.0)
    solid.link_source_buffer(np.full(solid.N, 1200.0))

    solid.set_ode_method(method)
    return solid


class HeatConductionImplicitEulerTests(unittest.TestCase):
    def _run_constant_steps(self, solid: HeatConduction1D, dt: float, n_steps: int):
        for _ in range(n_steps):
            if not solid.step(dt):
                return False
        return True

    def test_generic_solid_defaults_to_implicit_euler(self):
        solid = HeatConduction1D(
            mesh=Mesh1D(total_dim=0.1, n_volumes=2),
            material=ConstantSolid(),
        )
        self.assertEqual(solid.ode_method, "implicit_euler")

    def test_set_ode_method_accepts_implicit_and_keeps_existing(self):
        solid = _build_test_solid()
        for method in ("RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"):
            solid.set_ode_method(method)
            self.assertEqual(solid.ode_method, method)
            self.assertTrue(solid.step(1.0e-3))

        solid.set_ode_method("implicit_euler")
        self.assertEqual(solid.ode_method, "implicit_euler")

    def test_implicit_euler_matches_solve_ivp_baseline(self):
        reference = _build_test_solid(method="BDF")
        implicit = _build_test_solid(method="implicit_euler")

        ok_ref = self._run_constant_steps(reference, dt=0.25, n_steps=8)
        ok_implicit = self._run_constant_steps(implicit, dt=0.25, n_steps=8)

        self.assertTrue(ok_ref, "Reference solve_ivp step failed.")
        self.assertTrue(ok_implicit, "Implicit Euler step failed.")
        max_abs_delta = float(np.max(np.abs(reference.T - implicit.T)))
        self.assertLess(
            max_abs_delta,
            0.5,
            msg=f"Implicit Euler deviates from solve_ivp baseline by {max_abs_delta:.3f} K",
        )

    def test_dynamic_radiation_boundary_is_implicitly_linearized(self):
        reference = _build_test_solid(method="BDF", with_dynamic_radiation=True)
        implicit = _build_test_solid(method="implicit_euler", with_dynamic_radiation=True)

        ok_ref = self._run_constant_steps(reference, dt=0.20, n_steps=6)
        ok_implicit = self._run_constant_steps(implicit, dt=0.20, n_steps=6)

        self.assertTrue(ok_ref, "Reference solve_ivp step failed.")
        self.assertTrue(ok_implicit, "Implicit Euler step failed.")

        max_abs_delta = float(np.max(np.abs(reference.T - implicit.T)))
        self.assertLess(
            max_abs_delta,
            0.5,
            msg=f"Dynamic-radiation implicit result deviates by {max_abs_delta:.3f} K",
        )

        dyn_bc = next(
            (
                condition
                for condition in implicit.boundaries["inner"].conditions
                if isinstance(condition, DynamicRadiationResistanceBC)
            ),
            None,
        )
        self.assertIsNotNone(dyn_bc, "Expected a dynamic radiation BC on inner boundary.")
        self.assertTrue(np.all(np.isfinite(dyn_bc.h_rad)))
        self.assertTrue(np.all(np.array(dyn_bc.h_rad) > 0.0))
        self.assertTrue(np.all(np.array(dyn_bc.R_ext) > 0.0))

    def test_implicit_failure_falls_back_to_solve_ivp_with_warning(self):
        solid = _build_test_solid(method="implicit_euler")
        self.assertTrue(
            hasattr(solid, "_implicit_euler_step"),
            "Expected internal implicit Euler step hook for fallback testing.",
        )

        final_state = solid.T + 0.4
        fake_solution = type(
            "FakeSolveResult",
            (),
            {
                "success": True,
                "message": "ok",
                "y": np.column_stack((solid.T.copy(), final_state)),
            },
        )()

        with patch.object(solid, "_implicit_euler_step", return_value=False), \
                patch("Solvers.HeatConduction.HeatConduction.solve_ivp", return_value=fake_solution) as solve_ivp_mock, \
                patch("warnings.warn") as warn_mock, \
                patch("builtins.print") as print_mock, \
                warnings.catch_warnings(record=True) as captured_warnings:
            warnings.simplefilter("always")
            ok = solid.step(0.5)

            self.assertTrue(ok, "Implicit fallback step should report success.")
            self.assertEqual(solve_ivp_mock.call_count, 1)

            fallback_warned = any(
                "fallback" in str(item.message).lower() for item in captured_warnings
            ) or any(
                "fallback" in str(call.args[0]).lower()
                for call in print_mock.call_args_list
                if call.args
            )

            # Support either explicit warnings.warn usage or printed warning message.
            warn_called = warn_mock.called or fallback_warned
            self.assertTrue(warn_called, "Expected warning or warning-like fallback message.")
            self.assertTrue(np.allclose(solid.T, final_state))


if __name__ == "__main__":
    unittest.main()
