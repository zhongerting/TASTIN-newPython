import logging
from typing import Callable, Dict, List, Union

import numpy as np

from Components.BaseComponent import BaseComponent
from Solvers.Couplers import FluidSolidCouple, SolidSolidCouple1D, SolidSolidCouple2D
from Solvers.HeatConduction.HeatConduction import BaseHeatConduction
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from profiler import TEASAProfiler


logger = logging.getLogger(__name__)


class SystemManager:
    """
    协调液压网络、固体求解器、耦合器、宏观组件以及可选的点堆中子动力学，
    用于全局瞬态计算步。
    """

    def __init__(self, fluid_network: HydraulicNetwork, start_time: float = 0.0):
        self.global_time = float(start_time)
        self.fluid_solver = fluid_network
        if hasattr(self.fluid_solver, 'set_time'):
            self.fluid_solver.set_time(self.global_time)

        self.solid_components: Dict[str, BaseHeatConduction] = {}
        self.couplers: List[Union[FluidSolidCouple, SolidSolidCouple1D, SolidSolidCouple2D]] = []
        self.point_reactor = None
        self.components: List[BaseComponent] = []

        self._persistent_fluid_sources: List[Callable[['SystemManager'], None]] = []
        self.last_step_diagnostics = None
        self._component_step_states = []
        self._point_reactor_step_state = None
        logger.info(f"SystemManager 在 t={self.global_time}s 时初始化")
        logger.info(
            f"  - 流体求解器已连接: {self.fluid_solver.n_vol} 个控制体, {self.fluid_solver.n_junc} 个连接点"
        )

    # =========================================================================
    # Registration API
    # =========================================================================

    def add_point_reactor(self, reactor):
        self.point_reactor = reactor
        logger.info("  [Register] Neutronics: PointReactor attached to SystemManager.")

    def add_persistent_fluid_source(self, source_func: Callable[['SystemManager'], None]):
        """
        注册一个在瞬态源被清除后重新应用的流体源。

        source_func 接收 SystemManager；调用者在 callable 中显式写入
        需要的源项逻辑。
        """
        if not callable(source_func):
            raise TypeError("persistent fluid source must be callable")

        self._persistent_fluid_sources.append(source_func)
        return source_func

    def add_solid_component(self, component: BaseHeatConduction):
        """
        将固体传热组件注册到系统管理器中。

        参数:
            component: 继承自 BaseHeatConduction 的固体传热求解器实例

        异常:
            TypeError: 如果组件不是 BaseHeatConduction 的子类
            ValueError: 如果同名组件已存在
        """
        # 类型检查
        if not isinstance(component, BaseHeatConduction):
            raise TypeError(f"Component must inherit from BaseHeatConduction, got {type(component)}")

        # 获取或生成组件名称
        name = getattr(component, 'name', None)
        if not name:
            name = f"Solid_{len(self.solid_components)}"
            logger.warning(f"Component {component} has no 'name' attribute. Assigned: {name}")

        # 检查名称唯一性
        if name in self.solid_components:
            raise ValueError(f"Solid component with name '{name}' already exists in SystemManager!")

        # 注册组件
        self.solid_components[name] = component
        logger.info(f"  [Register] Solid: '{name}' (Type: {type(component).__name__}, Nodes: {component.N})")

    def add_coupler(self, coupler: Union[FluidSolidCouple, SolidSolidCouple1D, SolidSolidCouple2D]):
        if not hasattr(coupler, 'execute') and not hasattr(coupler, 'sync'):
            raise TypeError(f"Registered object {coupler} is not a valid coupler (missing execute/sync).")

        self.couplers.append(coupler)
        c_name = getattr(coupler, 'name', f"Coupler_{len(self.couplers)}")
        logger.info(f"  [Register] Coupler: '{c_name}' (Type: {type(coupler).__name__})")

    def add_component(self, component: BaseComponent):
        if component in self.components:
            raise ValueError(f"Component '{getattr(component, 'name', component)}' is already registered.")

        self.components.append(component)

        for solid in component.get_solids():
            self.add_solid_component(solid)

        for coupler in component.get_couplers():
            self.add_coupler(coupler)

        logger.info(
            f"  [Register] Component: '{component.name}' added. "
            f"Extracted {len(component.get_solids())} solids and {len(component.get_couplers())} couplers."
        )

    # =========================================================================
    # Initialization
    # =========================================================================

    def initialize_system(self, dt_init: float = 0.1, tol: float = 1e-5, max_iter: int = 500):
        """
        初始化系统，包括流体求解器的水力学初始化和耦合器的同步。

        参数:
            dt_init: 初始化时间步长，默认为0.1秒
            tol: 收敛容差，默认为1e-5
            max_iter: 最大迭代次数，默认为500

        异常:
            RuntimeError: 当水力学初始化失败时抛出
        """
        logger.info(f"System Initialization Started (dt={dt_init}, tol={tol})...")

        logger.info("  >> [Step 1] Initializing Hydraulics...")
        success = self.fluid_solver.initialize_hydraulics(dt=dt_init, tol=tol, max_iter=max_iter)
        if not success:
            logger.error("Hydraulic initialization FAILED! Check topology or boundary conditions.")
            raise RuntimeError("System initialization failed at hydraulic stage.")

        logger.info("  >> [Step 2] Synchronizing Couplers...")
        self._prepare_fluid_sources_for_coupling()
        self._run_couplers()

        logger.info("System Initialization Completed Successfully.\n")

    # =========================================================================
    # Fluid source lifecycle
    # =========================================================================

    def _iter_fluid_volumes(self):
        if hasattr(self.fluid_solver, 'volumes_obj'):
            return list(self.fluid_solver.volumes_obj)
        if hasattr(self.fluid_solver, 'volumes'):
            return list(self.fluid_solver.volumes)
        return []

    def _clear_fluid_sources(self):
        for vol in self._iter_fluid_volumes():
            if hasattr(vol, 'Q_wall'):
                vol.Q_wall = 0.0
            if hasattr(vol, 'Q_vol'):
                vol.Q_vol = 0.0
            if hasattr(vol, 'implicit_coeff'):
                vol.implicit_coeff = 0.0

    def _apply_persistent_fluid_sources(self):
        for source in self._persistent_fluid_sources:
            source(self)

    def _prepare_fluid_sources_for_coupling(self):
        self._clear_fluid_sources()
        self._apply_persistent_fluid_sources()

    def _capture_fluid_sources(self):
        snapshot = []
        for vol in self._iter_fluid_volumes():
            snapshot.append((
                getattr(vol, 'Q_wall', 0.0),
                getattr(vol, 'Q_vol', 0.0),
                getattr(vol, 'implicit_coeff', 0.0),
            ))
        return snapshot

    def _restore_fluid_sources(self, snapshot):
        for vol, (q_wall, q_vol, implicit_coeff) in zip(self._iter_fluid_volumes(), snapshot):
            if hasattr(vol, 'Q_wall'):
                vol.Q_wall = q_wall
            if hasattr(vol, 'Q_vol'):
                vol.Q_vol = q_vol
            if hasattr(vol, 'implicit_coeff'):
                vol.implicit_coeff = implicit_coeff

    # =========================================================================
    # Coupler scheduling and boundary cache lifecycle
    # =========================================================================

    def _refresh_solid_boundary_cache(self, update_flux: bool = False, current_time: float = None):
        cache_time = self.global_time if current_time is None else current_time
        for solid in self.solid_components.values():
            if hasattr(solid, '_update_properties'):
                solid._update_properties()
            if hasattr(solid, '_compute_internal_resistance'):
                solid._compute_internal_resistance()
            if hasattr(solid, '_update_boundaries_state'):
                solid._update_boundaries_state(current_time=cache_time)

            if update_flux:
                if hasattr(solid, '_compute_fluxes'):
                    solid._compute_fluxes(cache_time)
                elif hasattr(solid, 'boundaries'):
                    for boundary in solid.boundaries.values():
                        if hasattr(boundary, 'compute_net_flux_for_solver'):
                            boundary.compute_net_flux_for_solver()

    def _run_couplers(self, interface_relaxation: float = 1.0, current_time: float = None):
        self._refresh_solid_boundary_cache(update_flux=False, current_time=current_time)

        for coupler in self.couplers:
            if hasattr(coupler, 'sync'):
                coupler.sync()

        self._refresh_solid_boundary_cache(update_flux=False, current_time=current_time)

        for coupler in self.couplers:
            if isinstance(coupler, FluidSolidCouple):
                coupler.execute(interface_relaxation=interface_relaxation)
            elif hasattr(coupler, 'execute'):
                coupler.execute()

    # =========================================================================
    # Transient stepping
    # =========================================================================

    @TEASAProfiler.profile
    def step(
        self,
        dt: float,
        inner_iter: int = 1,
        convergence_tol: float = 1e-3,
        reactivity_control: float = 0.0,
        fail_on_fluid_nonconvergence: bool = False,
        interface_relaxation: float = 1.0,
        interface_convergence_tol: float = None,
    ):
        if inner_iter < 1:
            raise ValueError("inner_iter 必须 >= 1")
        if not (0.0 < interface_relaxation <= 1.0):
            raise ValueError("interface_relaxation 必须在 (0, 1] 范围内。")
        if interface_convergence_tol is not None and interface_convergence_tol < 0.0:
            raise ValueError("interface_convergence_tol must be non-negative or None.")
        if interface_convergence_tol is not None and inner_iter < 2:
            raise ValueError("interface_convergence_tol requires inner_iter >= 2.")

        t_start = self.global_time
        entry_fluid_sources = self._capture_fluid_sources()
        self._save_system_state(include_components=True, include_point_reactor=True)
        diagnostics = self._new_step_diagnostics(
            t_start=t_start,
            dt=dt,
            inner_iter=inner_iter,
            convergence_tol=convergence_tol,
            interface_relaxation=interface_relaxation,
            interface_convergence_tol=interface_convergence_tol,
        )
        self.last_step_diagnostics = diagnostics

        T_f_prev = None
        T_s_prev = {}
        component_neutronics_handled = False

        try:
            for coupler in self.couplers:
                reset_relaxation = getattr(coupler, 'reset_interface_relaxation', None)
                if callable(reset_relaxation):
                    reset_relaxation()

            self._prepare_fluid_sources_for_coupling()
            for comp in self.components:
                if hasattr(comp, 'pre_step'):
                    comp.pre_step(dt, t_start)

            base_fluid_sources = self._capture_fluid_sources()

            for k in range(inner_iter):
                self._restore_fluid_sources(base_fluid_sources)
                coupling_time = t_start if k == 0 else t_start + dt
                self._run_couplers(
                    interface_relaxation=interface_relaxation,
                    current_time=coupling_time,
                )
                coupler_diagnostics = self._collect_coupler_diagnostics()
                diagnostics["coupler_diagnostics_by_iteration"].append(coupler_diagnostics)
                if coupler_diagnostics:
                    diagnostics["couplers"] = coupler_diagnostics

                handled, fallback_power = self._advance_neutronics_for_iteration(
                    dt=dt,
                    reactivity_control=reactivity_control,
                    iteration_index=k,
                )
                component_neutronics_handled = component_neutronics_handled or handled

                if k > 0:
                    self._rollback_system_state(include_components=False, include_point_reactor=False)
                    self._apply_pending_nuclear_power(fallback_power)
                else:
                    self._apply_pending_nuclear_power(fallback_power)

                if hasattr(self.fluid_solver, 'set_time'):
                    self.fluid_solver.set_time(coupling_time)

                fluid_converged = self.fluid_solver.step_Picard(
                    dt,
                    max_iter=20 if inner_iter > 1 else 100,
                )
                diagnostics["fluid_converged_by_iteration"].append(bool(fluid_converged))

                if not fluid_converged:
                    message = f"Fluid solver NOT converged at t={t_start:.4f}s"
                    diagnostics["warnings"].append(message)
                    if fail_on_fluid_nonconvergence:
                        raise RuntimeError(message)
                    logger.warning(message)

                for name, solid in self.solid_components.items():
                    success = solid.step(dt)
                    if not success:
                        raise RuntimeError(f"Solid '{name}' integration failed at t={t_start:.4f}s")

                diagnostics["iterations"] = k + 1

                if inner_iter > 1:
                    T_f_curr = self._fluid_temperature_vector()
                    T_s_curr = self._solid_temperature_snapshot()
                    if k > 0:
                        err_f = self._max_abs_delta(T_f_curr, T_f_prev)
                        err_s = 0.0
                        for name, T_arr in T_s_curr.items():
                            if name in T_s_prev:
                                err_s = max(err_s, self._max_abs_delta(T_arr, T_s_prev[name]))

                        diagnostics["err_fluid_temperature"] = err_f
                        diagnostics["err_solid_temperature"] = err_s
                        temperature_converged = err_f < convergence_tol and err_s < convergence_tol

                        if interface_convergence_tol is None:
                            interface_converged = True
                        else:
                            interface_residual = self._interface_residual_from_diagnostics(coupler_diagnostics)
                            diagnostics["interface_residual"] = interface_residual
                            interface_converged = interface_residual < interface_convergence_tol

                        diagnostics["temperature_converged"] = temperature_converged
                        diagnostics["interface_converged"] = interface_converged

                        if temperature_converged and interface_converged:
                            diagnostics["converged"] = all(diagnostics["fluid_converged_by_iteration"])
                            break

                    T_f_prev = T_f_curr
                    T_s_prev = T_s_curr
                else:
                    diagnostics["temperature_converged"] = True
                    diagnostics["interface_converged"] = interface_convergence_tol is None
                    diagnostics["converged"] = bool(fluid_converged)

            self._commit_neutronics(component_neutronics_handled)

            self.global_time = t_start + dt
            self._sync_solid_times_to_global()
            self._refresh_solid_boundary_cache(update_flux=True)

            for comp in self.components:
                if hasattr(comp, 'post_step'):
                    comp.post_step(dt, self.global_time)

            diagnostics["status"] = "completed"
            diagnostics["t_end"] = self.global_time
            diagnostics["component_neutronics_handled"] = component_neutronics_handled

        except Exception as exc:
            self._rollback_system_state(include_components=True, include_point_reactor=True)
            self._restore_fluid_sources(entry_fluid_sources)
            self.global_time = t_start
            self._sync_solid_times_to_global()
            self._refresh_solid_boundary_cache(update_flux=True)
            diagnostics["status"] = "failed"
            diagnostics["t_end"] = t_start
            diagnostics["exception"] = f"{type(exc).__name__}: {exc}"
            raise

    def _new_step_diagnostics(
        self,
        *,
        t_start: float,
        dt: float,
        inner_iter: int,
        convergence_tol: float,
        interface_relaxation: float,
        interface_convergence_tol: float,
    ):
        return {
            "status": "running",
            "t_start": t_start,
            "t_end": None,
            "dt": dt,
            "inner_iter_limit": inner_iter,
            "iterations": 0,
            "convergence_tol": convergence_tol,
            "interface_relaxation": interface_relaxation,
            "interface_convergence_tol": interface_convergence_tol,
            "converged": False,
            "temperature_converged": False,
            "interface_converged": interface_convergence_tol is None,
            "err_fluid_temperature": None,
            "err_solid_temperature": None,
            "interface_residual": None,
            "fluid_converged_by_iteration": [],
            "warnings": [],
            "couplers": [],
            "coupler_diagnostics_by_iteration": [],
            "component_neutronics_handled": False,
            "exception": None,
        }

    def _collect_coupler_diagnostics(self):
        diagnostics = []
        for coupler in self.couplers:
            getter = getattr(coupler, 'get_coupling_diagnostics', None)
            if not callable(getter):
                continue

            data = getter()
            if data is None:
                continue

            entry = dict(data)
            entry.setdefault("name", getattr(coupler, "name", type(coupler).__name__))
            entry.setdefault("type", type(coupler).__name__)
            diagnostics.append(entry)

        return diagnostics

    @staticmethod
    def _interface_residual_from_diagnostics(coupler_diagnostics) -> float:
        residuals = []
        for diagnostics in coupler_diagnostics:
            residual = diagnostics.get("interface_residual")
            if residual is None:
                continue
            residuals.append(float(residual))

        if not residuals:
            raise RuntimeError(
                "interface_convergence_tol requires at least one coupler "
                "diagnostic entry with 'interface_residual'."
            )

        return max(residuals)

    def _advance_neutronics_for_iteration(
        self,
        *,
        dt: float,
        reactivity_control: float,
        iteration_index: int,
    ):
        component_handled = False
        for comp in self.components:
            if hasattr(comp, 'advance_neutronics'):
                handled = comp.advance_neutronics(
                    dt=dt,
                    reactivity_control=reactivity_control,
                    iteration_index=iteration_index,
                )
                component_handled = component_handled or bool(handled)

        if self.point_reactor is None or component_handled:
            return component_handled, None

        total_feedback = 0.0
        for solid in self.solid_components.values():
            if hasattr(solid, 'get_reactivity_feedback'):
                total_feedback += solid.get_reactivity_feedback()

        self.point_reactor.step(dt, reactivity_control, total_feedback)
        return False, (
            self.point_reactor.fission_power,
            self.point_reactor.decay_power,
            self.point_reactor.total_power,
        )

    def _apply_pending_nuclear_power(self, fallback_power):
        """
        将待处理的中子功率应用到所有固体组件。

        参数:
            fallback_power: 包含裂变功率、衰变功率和总功率的元组 (p_fiss, p_decay, p_total)，
                           如果为 None 则不执行任何操作
        """
        if fallback_power is None:
            return

        p_fiss, p_decay, p_total = fallback_power
        for solid in self.solid_components.values():
            if hasattr(solid, 'set_nuclear_power'):
                solid.set_nuclear_power(p_fiss, p_decay, p_total)

    def _commit_neutronics(self, component_neutronics_handled: bool):
        component_committed = False
        for comp in self.components:
            if hasattr(comp, 'commit_neutronics'):
                component_committed = component_committed or bool(comp.commit_neutronics())

        if self.point_reactor is not None and not component_neutronics_handled and not component_committed:
            self.point_reactor.commit()

    def _fluid_temperature_vector(self) -> np.ndarray:
        if hasattr(self.fluid_solver, 'T_vec'):
            return np.asarray(self.fluid_solver.T_vec, dtype=float).copy()
        return np.array([getattr(vol, 'T', 0.0) for vol in self._iter_fluid_volumes()], dtype=float)

    def _solid_temperature_snapshot(self):
        snapshot = {}
        for name, solid in self.solid_components.items():
            if hasattr(solid, 'T'):
                snapshot[name] = np.asarray(solid.T, dtype=float).copy()
        return snapshot

    @staticmethod
    def _max_abs_delta(current, previous) -> float:
        if current is None or previous is None:
            return float('inf')
        current_arr = np.asarray(current, dtype=float)
        previous_arr = np.asarray(previous, dtype=float)
        if current_arr.size == 0 and previous_arr.size == 0:
            return 0.0
        return float(np.max(np.abs(current_arr - previous_arr)))

    # =========================================================================
    # Adaptive time-step control
    # =========================================================================

    def compute_adaptive_dt(
        self,
        min_dt: float = 1e-4,
        max_dt: float = 0.5,
        safety_factor: float = 0.8,
        use_last_step_diagnostics: bool = True,
    ) -> float:
        raw_dt_fluid = self.fluid_solver.get_max_stable_dt(max_limit=max_dt)
        dt_fluid = raw_dt_fluid * safety_factor

        dt_coupler = max_dt
        for coupler in self.couplers:
            if hasattr(coupler, 'get_max_stable_dt'):
                dt_c = coupler.get_max_stable_dt(safety_factor=safety_factor, max_limit=max_dt)
                dt_coupler = min(dt_coupler, dt_c)

        dt_target = min(dt_fluid, dt_coupler)
        if dt_target < min_dt:
            logger.warning(
                "Adaptive dt target %.6e is below min_dt %.6e; returning the physical target.",
                dt_target,
                min_dt,
            )
            dt_clamped = dt_target
        else:
            dt_clamped = min(dt_target, max_dt)

        if use_last_step_diagnostics and self.last_step_diagnostics is not None:
            diag = self.last_step_diagnostics
            fluid_ok = all(diag.get("fluid_converged_by_iteration", [True]))
            picard_ok = bool(diag.get("converged", True))
            iterations = int(diag.get("iterations", 0) or 0)
            iteration_limit = int(diag.get("inner_iter_limit", 1) or 1)

            if not fluid_ok or not picard_ok:
                dt_clamped *= 0.5
            elif iteration_limit > 1 and iterations >= iteration_limit:
                dt_clamped *= 0.8

            if dt_clamped < min_dt:
                logger.warning(
                    "Adaptive dt target %.6e is below min_dt %.6e after convergence control.",
                    dt_clamped,
                    min_dt,
                )

        if hasattr(self, '_last_dt'):
            dt_final = min(dt_clamped, 1.2 * self._last_dt)
        else:
            dt_final = dt_clamped

        self._last_dt = dt_final
        return dt_final

    # =========================================================================
    # State management for inner iteration and restart
    # =========================================================================

    def _save_system_state(self, include_components: bool = False, include_point_reactor: bool = False):
        if hasattr(self.fluid_solver, 'save_state'):
            self.fluid_solver.save_state()

        for solid in self.solid_components.values():
            if hasattr(solid, 'save_state'):
                solid.save_state()

        if include_components:
            self._component_step_states = []
            for comp in self.components:
                save_step_state = getattr(comp, 'save_step_state', None)
                if callable(save_step_state):
                    self._component_step_states.append((comp, save_step_state()))

        if include_point_reactor and self.point_reactor is not None:
            save_step_state = getattr(self.point_reactor, 'save_step_state', None)
            self._point_reactor_step_state = save_step_state() if callable(save_step_state) else None

    def _rollback_system_state(self, include_components: bool = False, include_point_reactor: bool = False):
        if hasattr(self.fluid_solver, 'load_state'):
            self.fluid_solver.load_state()

        for solid in self.solid_components.values():
            if hasattr(solid, 'load_state'):
                solid.load_state()

        if include_components:
            for comp, state in reversed(self._component_step_states):
                load_step_state = getattr(comp, 'load_step_state', None)
                if callable(load_step_state):
                    load_step_state(state)

        if include_point_reactor and self.point_reactor is not None:
            load_step_state = getattr(self.point_reactor, 'load_step_state', None)
            if callable(load_step_state):
                load_step_state(self._point_reactor_step_state)

    def _sync_solid_times_to_global(self):
        for solid in self.solid_components.values():
            if hasattr(solid, 'current_time'):
                solid.current_time = self.global_time

    def save_global_state(self, filepath: str):
        logger.info(f"Saving global state to {filepath} ...")
        global_state = {'System/global_time': np.array([self.global_time])}

        if hasattr(self, '_last_dt'):
            global_state['System/last_dt'] = np.array([self._last_dt])

        if hasattr(self.fluid_solver, 'get_state_dict'):
            global_state.update(self.fluid_solver.get_state_dict(prefix="Fluid"))

        for name, solid in self.solid_components.items():
            if hasattr(solid, 'get_state_dict'):
                global_state.update(solid.get_state_dict(prefix=f"Solid_{name}"))

        for comp in self.components:
            if hasattr(comp, 'get_state_dict'):
                global_state.update(comp.get_state_dict(prefix=f"Macro_{comp.name}"))

        if self.point_reactor is not None and hasattr(self.point_reactor, 'get_state_dict'):
            global_state.update(self.point_reactor.get_state_dict(prefix="PointReactor"))

        np.savez_compressed(filepath, **global_state)
        logger.info(f"Global state successfully saved to {filepath}.")

    def load_global_state(self, filepath: str):
        import os

        logger.info(f"Loading global state from {filepath} ...")
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Restart file not found: {filepath}")

        with np.load(filepath, allow_pickle=False) as data:
            data_dict = dict(data)

            if 'System/global_time' in data_dict:
                self.global_time = float(data_dict['System/global_time'][0])
            if 'System/last_dt' in data_dict:
                self._last_dt = float(data_dict['System/last_dt'][0])

            if hasattr(self.fluid_solver, 'load_state_dict'):
                self.fluid_solver.load_state_dict(data_dict, prefix="Fluid")

            for name, solid in self.solid_components.items():
                if hasattr(solid, 'load_state_dict'):
                    solid.load_state_dict(data_dict, prefix=f"Solid_{name}")

            for comp in self.components:
                if hasattr(comp, 'load_state_dict'):
                    comp.load_state_dict(data_dict, prefix=f"Macro_{comp.name}")

            if self.point_reactor is not None and hasattr(self.point_reactor, 'load_state_dict'):
                self.point_reactor.load_state_dict(data_dict, prefix="PointReactor")

        self._sync_solid_times_to_global()
        self._refresh_solid_boundary_cache(update_flux=True)
        self._prepare_fluid_sources_for_coupling()
        self._run_couplers()

        logger.info(f"Global state successfully loaded. Resuming from t={self.global_time}s.")
