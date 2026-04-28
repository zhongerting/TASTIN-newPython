import os
import sys
import unittest

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from Components.BaseComponent import BaseComponent
from Solvers.Couplers import FluidSolidCouple
from Solvers.SystemManager import SystemManager


class FakeVolume:
    def __init__(self, temperature=300.0):
        self.T = float(temperature)
        self.Q_wall = 0.0
        self.Q_vol = 0.0
        self.implicit_coeff = 0.0

    def add_coupling_source(self, explicit_part, implicit_factor):
        self.Q_wall += float(explicit_part)
        self.implicit_coeff += float(implicit_factor)


class FakeFluid:
    def __init__(self, volumes=None, converged=True, stable_dt=1.0):
        self.volumes_obj = volumes or [FakeVolume()]
        self.n_vol = len(self.volumes_obj)
        self.n_junc = 0
        self.converged = converged
        self.stable_dt = stable_dt
        self.step_calls = 0

    @property
    def T_vec(self):
        return np.array([vol.T for vol in self.volumes_obj], dtype=float)

    def initialize_hydraulics(self, dt=0.1, tol=1e-5, max_iter=500):
        return True

    def step_Picard(self, dt, max_iter=20, tol=1e-4):
        self.step_calls += 1
        for vol in self.volumes_obj:
            q_total = vol.Q_wall + vol.Q_vol - vol.implicit_coeff * vol.T
            vol.T += dt * q_total
        return self.converged

    def get_max_stable_dt(self, max_limit=0.5):
        return min(self.stable_dt, max_limit)

    def save_state(self):
        self._backup = [vol.T for vol in self.volumes_obj]

    def load_state(self):
        for vol, temp in zip(self.volumes_obj, self._backup):
            vol.T = temp

    def load_state_dict(self, data, prefix):
        pass


class FakeBoundary:
    def __init__(self, solid, log=None):
        self.solid = solid
        self.log = log
        self.T_surface = np.array([solid.T[0]], dtype=float)
        self.current_flux = np.zeros(1, dtype=float)

    def compute_net_flux_for_solver(self):
        self.current_flux[:] = self.T_surface * 2.0
        if self.log is not None:
            self.log.append("boundary_flux")
        return self.current_flux


class FakeSolid:
    def __init__(self, name="solid", temperature=400.0, fail=False, log=None):
        self.name = name
        self.N = 1
        self.T = np.array([temperature], dtype=float)
        self.current_time = 0.0
        self.fail = fail
        self.log = log if log is not None else []
        self.boundaries = {"left": FakeBoundary(self, self.log)}
        self.nuclear_power = None

    def save_state(self):
        self._backup_T = self.T.copy()
        self._backup_time = self.current_time

    def load_state(self):
        self.T[:] = self._backup_T
        self.current_time = self._backup_time

    def step(self, dt):
        if self.fail:
            return False
        self.T += dt
        self.current_time += dt
        return True

    def _update_properties(self):
        self.log.append("properties")

    def _compute_internal_resistance(self):
        self.log.append("resistance")

    def _update_boundaries_state(self, current_time=None):
        self.log.append("refresh")
        self.boundaries["left"].T_surface[:] = self.T[0]
        if current_time is not None:
            self.boundary_time = current_time

    def _compute_fluxes(self, current_time):
        self.log.append("solid_flux")
        return self.boundaries["left"].compute_net_flux_for_solver()

    def set_nuclear_power(self, p_fiss, p_decay, p_total):
        self.nuclear_power = (p_fiss, p_decay, p_total)


class HookComponent(BaseComponent):
    def __init__(self, name="hook", volume=None, pre_source=0.0):
        super().__init__(name)
        self.volume = volume
        self.pre_source = pre_source
        self.pre_times = []
        self.post_times = []

    def pre_step(self, dt, current_time):
        self.pre_times.append(current_time)
        if self.volume is not None:
            self.volume.Q_wall += self.pre_source

    def post_step(self, dt, current_time):
        self.post_times.append(current_time)


class StepStateComponent(BaseComponent):
    def __init__(self, name="step-state"):
        super().__init__(name)
        self.value = 0.0
        self.save_calls = 0
        self.load_calls = 0

    def save_step_state(self):
        self.save_calls += 1
        return {"value": self.value}

    def load_step_state(self, state):
        self.load_calls += 1
        self.value = state["value"]

    def pre_step(self, dt, current_time):
        self.value += 10.0


class HandledIterationComponent(BaseComponent):
    def __init__(self):
        super().__init__("handled-iteration")
        self.iteration_index = None

    def save_step_state(self):
        return {"iteration_index": self.iteration_index}

    def load_step_state(self, state):
        self.iteration_index = state["iteration_index"]

    def advance_neutronics(self, dt, reactivity_control=0.0, iteration_index=0):
        self.iteration_index = iteration_index
        return True

    def commit_neutronics(self):
        return True


class LogCoupler:
    def __init__(self, name, log, volume=None, boundary=None):
        self.name = name
        self.log = log
        self.volume = volume
        self.boundary = boundary
        self.execute_count = 0
        self.seen_surface_temperatures = []

    def sync(self):
        self.log.append(f"{self.name}:sync")

    def execute(self):
        self.execute_count += 1
        self.log.append(f"{self.name}:execute")
        if self.boundary is not None:
            self.seen_surface_temperatures.append(float(self.boundary.T_surface[0]))
        if self.volume is not None:
            self.volume.Q_wall += 100.0 * self.execute_count
            self.volume.implicit_coeff += 2.0 + self.execute_count


class ExecuteOnlyCoupler:
    def __init__(self, name, log, volume=None, boundary=None):
        self.name = name
        self.log = log
        self.volume = volume
        self.boundary = boundary
        self.execute_count = 0
        self.seen_surface_temperatures = []

    def execute(self):
        self.execute_count += 1
        self.log.append(f"{self.name}:execute")
        if self.boundary is not None:
            self.seen_surface_temperatures.append(float(self.boundary.T_surface[0]))
        if self.volume is not None:
            self.volume.Q_wall += 100.0 * self.execute_count
            self.volume.implicit_coeff += 2.0 + self.execute_count


class DiagnosticCoupler:
    def __init__(self, residuals):
        self.name = "diagnostic"
        self.residuals = list(residuals)
        self.execute_count = 0

    def execute(self):
        self.execute_count += 1

    def get_coupling_diagnostics(self):
        index = min(max(self.execute_count - 1, 0), len(self.residuals) - 1)
        return {
            "name": self.name,
            "interface_residual": self.residuals[index],
        }


class BoundaryTimeCoupler:
    def __init__(self, boundary, solid):
        self.boundary = boundary
        self.solid = solid
        self.boundary_times = []
        self.solid_times = []

    def execute(self):
        self.boundary_times.append(self.solid.boundary_time)
        self.solid_times.append(self.solid.current_time)


class SyncOnlyCoupler:
    def __init__(self, name, log):
        self.name = name
        self.log = log

    def sync(self):
        self.log.append(f"{self.name}:sync")


class ResetTrackingCoupler:
    def __init__(self):
        self.reset_count = 0
        self.execute_count = 0

    def reset_interface_relaxation(self):
        self.reset_count += 1

    def execute(self):
        self.execute_count += 1


class FakeMaterial:
    def viscosity(self, T, P):
        return np.ones_like(T, dtype=float)

    def conductivity(self, T, P):
        return np.ones_like(T, dtype=float)

    def heat_capacity(self, T, P):
        return np.ones_like(T, dtype=float)

    def prandtl_number(self, T, P):
        return np.ones_like(T, dtype=float)


class SequenceCorrelation:
    def __init__(self, values):
        self.values = list(values)
        self.calls = 0

    def __call__(self, Re, Pr, P_D_ratio):
        value = self.values[min(self.calls, len(self.values) - 1)]
        self.calls += 1
        return np.full_like(Re, value, dtype=float)


class FakeResistanceBC:
    def __init__(self, T_ext, R_ext, shape):
        self.T_ext = np.asarray(T_ext, dtype=float).reshape(shape).copy()
        self.R_ext = np.asarray(R_ext, dtype=float).reshape(shape).copy()


class FakeFluidSolidBoundary:
    def __init__(self, wall_temperatures):
        self.shape = (1,)
        self.T_surface = np.array([wall_temperatures[0]], dtype=float)
        self.wall_temperatures = list(wall_temperatures)
        self.flux_calls = 0
        self.solid_bc = None

    def add_resistance_condition(self, T_ext, R_ext):
        self.solid_bc = FakeResistanceBC(T_ext, R_ext, self.shape)
        return self.solid_bc

    def compute_net_flux_for_solver(self):
        index = min(self.flux_calls, len(self.wall_temperatures) - 1)
        self.T_surface[:] = self.wall_temperatures[index]
        self.flux_calls += 1
        return np.zeros(self.shape, dtype=float)


class FakeFluidSolidChannel(FakeFluid):
    def __init__(self, wall_temperatures=(400.0, 800.0)):
        super().__init__([FakeVolume(temperature=300.0)])
        self.volumes = self.volumes_obj
        self.n_nodes = 1
        self.node_length = 1.0
        self.d_h = 1.0
        self.material = FakeMaterial()
        self.pressure = np.array([1.0], dtype=float)
        self.density = np.array([1.0], dtype=float)
        self.velocity = np.array([1.0], dtype=float)
        self.wall_temperatures = wall_temperatures

    @property
    def temperature_vector(self):
        return np.array([vol.T for vol in self.volumes], dtype=float)

    @property
    def pressure_vector(self):
        return self.pressure.copy()

    @property
    def density_vector(self):
        return self.density.copy()

    @property
    def velocity_vector(self):
        return self.velocity.copy()

    def add_coupling_source_distribution(self, explicit_arr, implicit_arr):
        for vol, explicit, implicit in zip(self.volumes, explicit_arr, implicit_arr):
            vol.add_coupling_source(explicit, implicit)

    def step_Picard(self, dt, max_iter=20, tol=1e-4):
        self.step_calls += 1
        for vol in self.volumes:
            vol.T = 500.0
        return self.converged


class FakePointReactor:
    def __init__(self):
        self.step_calls = 0
        self.commit_calls = 0
        self.fission_power = 10.0
        self.decay_power = 1.0
        self.total_power = 11.0

    def step(self, dt, reactivity_control, reactivity_feedback):
        self.step_calls += 1
        self.fission_power += 1.0
        self.total_power = self.fission_power + self.decay_power
        return True

    def commit(self):
        self.commit_calls += 1
        return True

    def save_step_state(self):
        return {
            "step_calls": self.step_calls,
            "commit_calls": self.commit_calls,
            "fission_power": self.fission_power,
            "decay_power": self.decay_power,
            "total_power": self.total_power,
        }

    def load_step_state(self, state):
        self.step_calls = state["step_calls"]
        self.commit_calls = state["commit_calls"]
        self.fission_power = state["fission_power"]
        self.decay_power = state["decay_power"]
        self.total_power = state["total_power"]


class NeutronicsComponent(BaseComponent):
    def __init__(self, handled=False, committed=False):
        super().__init__("neutronics")
        self.handled = handled
        self.committed = committed
        self.advance_calls = 0
        self.commit_calls = 0

    def advance_neutronics(self, dt, reactivity_control=0.0, iteration_index=0):
        self.advance_calls += 1
        return self.handled

    def commit_neutronics(self):
        self.commit_calls += 1
        return self.committed


def make_manager(fluid=None):
    return SystemManager(fluid or FakeFluid())


class SystemManagerLifecycleTests(unittest.TestCase):
    def test_persistent_and_pre_step_fluid_sources_survive_clear_without_accumulating(self):
        vol = FakeVolume(temperature=300.0)
        fluid = FakeFluid([vol])
        manager = make_manager(fluid)
        component = HookComponent(volume=vol, pre_source=5.0)
        manager.components.append(component)
        manager.add_persistent_fluid_source(lambda _manager: vol.add_coupling_source(10.0, 0.0))

        manager.step(1.0)
        manager.step(1.0)

        self.assertAlmostEqual(vol.T, 330.0)
        self.assertAlmostEqual(vol.Q_wall, 15.0)
        self.assertEqual(component.pre_times, [0.0, 1.0])
        self.assertEqual(component.post_times, [1.0, 2.0])

    def test_solid_failure_rolls_back_and_skips_post_step_and_neutronics_commit(self):
        vol = FakeVolume(temperature=300.0)
        fluid = FakeFluid([vol])
        manager = make_manager(fluid)
        manager.solid_components["bad"] = FakeSolid(fail=True)
        hook = HookComponent(volume=vol, pre_source=7.0)
        stateful = StepStateComponent()
        point = FakePointReactor()
        manager.components.append(hook)
        manager.components.append(stateful)
        manager.add_point_reactor(point)

        with self.assertRaisesRegex(RuntimeError, "Solid 'bad' integration failed"):
            manager.step(1.0)

        self.assertAlmostEqual(manager.global_time, 0.0)
        self.assertAlmostEqual(vol.T, 300.0)
        self.assertAlmostEqual(vol.Q_wall, 0.0)
        self.assertAlmostEqual(stateful.value, 0.0)
        self.assertEqual(stateful.save_calls, 1)
        self.assertEqual(stateful.load_calls, 1)
        self.assertEqual(hook.pre_times, [0.0])
        self.assertEqual(hook.post_times, [])
        self.assertEqual(point.step_calls, 0)
        self.assertEqual(point.commit_calls, 0)
        self.assertAlmostEqual(point.fission_power, 10.0)
        self.assertEqual(manager.last_step_diagnostics["status"], "failed")

    def test_fluid_nonconvergence_default_warns_but_strict_mode_fails(self):
        manager = make_manager(FakeFluid(converged=False))
        manager.step(0.5)
        self.assertAlmostEqual(manager.global_time, 0.5)

        strict_manager = make_manager(FakeFluid(converged=False))
        with self.assertRaisesRegex(RuntimeError, "Fluid solver NOT converged"):
            strict_manager.step(0.5, fail_on_fluid_nonconvergence=True)
        self.assertAlmostEqual(strict_manager.global_time, 0.0)

    def test_post_step_receives_advanced_time_and_solid_time_is_synchronized(self):
        manager = make_manager(FakeFluid())
        manager.global_time = 2.0
        solid = FakeSolid()
        hook = HookComponent()
        manager.solid_components["solid"] = solid
        manager.components.append(hook)

        manager.step(0.25)

        self.assertEqual(hook.pre_times, [2.0])
        self.assertEqual(hook.post_times, [2.25])
        self.assertAlmostEqual(solid.current_time, 2.25)

    def test_couplers_run_sync_before_execute_and_refresh_surface_cache(self):
        log = []
        vol = FakeVolume()
        manager = make_manager(FakeFluid([vol]))
        solid = FakeSolid(temperature=455.0, log=log)
        manager.solid_components["solid"] = solid
        execute_coupler = ExecuteOnlyCoupler("execute", log, boundary=solid.boundaries["left"])
        sync_coupler = SyncOnlyCoupler("sync", log)
        manager.add_coupler(execute_coupler)
        manager.add_coupler(sync_coupler)

        manager.step(0.0)

        self.assertLess(log.index("sync:sync"), log.index("execute:execute"))
        self.assertEqual(execute_coupler.seen_surface_temperatures, [455.0])

    def test_inner_picard_keeps_coupler_implicit_source_not_one_sided_average(self):
        vol = FakeVolume()
        manager = make_manager(FakeFluid([vol]))
        coupler = ExecuteOnlyCoupler("fluid-solid", [], volume=vol)
        manager.add_coupler(coupler)

        manager.step(0.0, inner_iter=2, convergence_tol=0.0)

        self.assertEqual(coupler.execute_count, 2)
        self.assertAlmostEqual(vol.Q_wall, 200.0)
        self.assertAlmostEqual(vol.implicit_coeff, 4.0)
        self.assertFalse(hasattr(manager, "_fluid_total_Q_backup"))

    def test_picard_correction_refreshes_solid_cache_at_predicted_time(self):
        manager = make_manager(FakeFluid())
        manager.global_time = 2.0
        solid = FakeSolid()
        solid.current_time = 2.0
        manager.solid_components["solid"] = solid
        coupler = BoundaryTimeCoupler(solid.boundaries["left"], solid)
        manager.add_coupler(coupler)

        manager.step(0.25, inner_iter=2, convergence_tol=0.0)

        self.assertEqual(coupler.boundary_times, [2.0, 2.25])
        self.assertEqual(coupler.solid_times, [2.0, 2.25])

    def test_picard_thermal_rollback_keeps_latest_component_neutronics_trial(self):
        manager = make_manager(FakeFluid())
        component = HandledIterationComponent()
        manager.components.append(component)

        manager.step(0.0, inner_iter=2, convergence_tol=0.0)

        self.assertEqual(component.iteration_index, 1)

    def test_interface_relaxation_state_resets_once_per_global_step(self):
        manager = make_manager(FakeFluid())
        coupler = ResetTrackingCoupler()
        manager.add_coupler(coupler)

        manager.step(0.0, inner_iter=3, convergence_tol=0.0, interface_relaxation=0.5)

        self.assertEqual(coupler.reset_count, 1)
        self.assertEqual(coupler.execute_count, 3)

    def test_interface_relaxation_damps_fluid_source_and_solid_robin_together(self):
        fluid = FakeFluidSolidChannel()
        boundary = FakeFluidSolidBoundary(wall_temperatures=(400.0, 800.0))
        manager = make_manager(fluid)
        coupler = FluidSolidCouple(
            name="fluid-solid",
            fluid=fluid,
            solid_boundary_region=boundary,
            heated_perimeter=1.0,
            correlation_func=SequenceCorrelation((10.0, 30.0)),
        )
        manager.add_coupler(coupler)

        manager.step(0.0, inner_iter=2, convergence_tol=0.0, interface_relaxation=0.5)

        vol = fluid.volumes[0]
        self.assertAlmostEqual(vol.Q_wall, 12000.0)
        self.assertAlmostEqual(vol.implicit_coeff, 20.0)
        self.assertAlmostEqual(coupler.solid_bc.T_ext[0], 400.0)
        self.assertAlmostEqual(coupler.solid_bc.R_ext[0], 0.05)
        diagnostics = coupler.get_coupling_diagnostics()
        self.assertTrue(diagnostics["relaxed"])
        self.assertGreater(diagnostics["interface_residual"], 0.0)
        self.assertGreater(diagnostics["source_residual"], 0.0)
        self.assertEqual(manager.last_step_diagnostics["couplers"][0]["name"], "fluid-solid")

    def test_interface_relaxation_one_keeps_fluid_solid_coupler_default_behavior(self):
        fluid = FakeFluidSolidChannel()
        boundary = FakeFluidSolidBoundary(wall_temperatures=(400.0, 800.0))
        manager = make_manager(fluid)
        coupler = FluidSolidCouple(
            name="fluid-solid",
            fluid=fluid,
            solid_boundary_region=boundary,
            heated_perimeter=1.0,
            correlation_func=SequenceCorrelation((10.0, 30.0)),
        )
        manager.add_coupler(coupler)

        manager.step(0.0, inner_iter=2, convergence_tol=0.0)

        vol = fluid.volumes[0]
        self.assertAlmostEqual(vol.Q_wall, 24000.0)
        self.assertAlmostEqual(vol.implicit_coeff, 30.0)
        self.assertAlmostEqual(coupler.solid_bc.T_ext[0], 500.0)
        self.assertAlmostEqual(coupler.solid_bc.R_ext[0], 1.0 / 30.0)

    def test_last_step_diagnostics_records_picard_convergence(self):
        manager = make_manager(FakeFluid())

        manager.step(0.0, inner_iter=2, convergence_tol=1.0)

        diagnostics = manager.last_step_diagnostics
        self.assertEqual(diagnostics["status"], "completed")
        self.assertEqual(diagnostics["iterations"], 2)
        self.assertTrue(diagnostics["converged"])
        self.assertTrue(diagnostics["temperature_converged"])
        self.assertEqual(diagnostics["fluid_converged_by_iteration"], [True, True])
        self.assertAlmostEqual(diagnostics["t_end"], 0.0)

    def test_interface_convergence_uses_public_coupler_diagnostics(self):
        manager = make_manager(FakeFluid())
        coupler = DiagnosticCoupler(residuals=[10.0, 0.1])
        manager.add_coupler(coupler)

        manager.step(
            0.0,
            inner_iter=3,
            convergence_tol=1.0,
            interface_convergence_tol=0.5,
        )

        diagnostics = manager.last_step_diagnostics
        self.assertEqual(diagnostics["iterations"], 2)
        self.assertTrue(diagnostics["converged"])
        self.assertTrue(diagnostics["interface_converged"])
        self.assertAlmostEqual(diagnostics["interface_residual"], 0.1)

    def test_interface_convergence_requires_coupler_residual_diagnostics(self):
        manager = make_manager(FakeFluid())
        manager.add_coupler(ExecuteOnlyCoupler("execute", []))

        with self.assertRaisesRegex(RuntimeError, "interface_convergence_tol requires"):
            manager.step(
                0.0,
                inner_iter=2,
                convergence_tol=1.0,
                interface_convergence_tol=1.0,
            )

        self.assertEqual(manager.last_step_diagnostics["status"], "failed")

    def test_interface_convergence_requires_multiple_picard_iterations(self):
        manager = make_manager(FakeFluid())

        with self.assertRaisesRegex(ValueError, "inner_iter >= 2"):
            manager.step(0.0, inner_iter=1, interface_convergence_tol=1.0)

    def test_neutronics_fallback_uses_return_value_not_method_presence(self):
        manager = make_manager(FakeFluid())
        component = NeutronicsComponent(handled=False, committed=False)
        point = FakePointReactor()
        manager.components.append(component)
        manager.add_point_reactor(point)

        manager.step(0.0)

        self.assertEqual(component.advance_calls, 1)
        self.assertEqual(component.commit_calls, 1)
        self.assertEqual(point.step_calls, 1)
        self.assertEqual(point.commit_calls, 1)

    def test_neutronics_component_truthy_result_suppresses_fallback(self):
        manager = make_manager(FakeFluid())
        component = NeutronicsComponent(handled=True, committed=True)
        point = FakePointReactor()
        manager.components.append(component)
        manager.add_point_reactor(point)

        manager.step(0.0)

        self.assertEqual(point.step_calls, 0)
        self.assertEqual(point.commit_calls, 0)

    def test_compute_adaptive_dt_does_not_raise_physical_limit_to_min_dt(self):
        manager = make_manager(FakeFluid(stable_dt=1.0e-6))

        dt = manager.compute_adaptive_dt(min_dt=1.0e-4, max_dt=1.0, safety_factor=1.0)

        self.assertAlmostEqual(dt, 1.0e-6)

    def test_compute_adaptive_dt_shrinks_after_nonconverged_step_diagnostics(self):
        manager = make_manager(FakeFluid(stable_dt=1.0))
        manager.last_step_diagnostics = {
            "fluid_converged_by_iteration": [False],
            "converged": False,
            "iterations": 1,
            "inner_iter_limit": 1,
        }

        dt = manager.compute_adaptive_dt(min_dt=1.0e-4, max_dt=1.0, safety_factor=1.0)

        self.assertAlmostEqual(dt, 0.5)

    def test_invalid_coupler_and_duplicate_component_registration_fail_fast(self):
        manager = make_manager(FakeFluid())
        component = HookComponent()
        manager.add_component(component)

        with self.assertRaises(TypeError):
            manager.add_coupler(object())
        with self.assertRaises(TypeError):
            manager.add_persistent_fluid_source(object())
        with self.assertRaisesRegex(ValueError, "already registered"):
            manager.add_component(component)

    def test_load_global_state_synchronizes_solid_time_and_boundary_cache(self):
        manager = make_manager(FakeFluid())
        solid = FakeSolid(temperature=612.0)
        manager.solid_components["solid"] = solid

        restart = os.path.join(os.path.dirname(__file__), "_tmp_system_manager_restart.npz")
        try:
            np.savez_compressed(restart, **{"System/global_time": np.array([5.0])})
            manager.load_global_state(restart)
        finally:
            if os.path.exists(restart):
                os.remove(restart)

        self.assertAlmostEqual(manager.global_time, 5.0)
        self.assertAlmostEqual(solid.current_time, 5.0)
        self.assertAlmostEqual(solid.boundaries["left"].T_surface[0], 612.0)


if __name__ == "__main__":
    unittest.main()
