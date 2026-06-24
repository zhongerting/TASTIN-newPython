import numpy as np
import warnings
from typing import Dict, Optional, Union, Any, Callable
from dataclasses import dataclass
from abc import ABC, abstractmethod
from scipy.integrate import solve_ivp
from scipy.sparse import spdiags, coo_matrix
from scipy.sparse.linalg import spsolve

from Materials.Base import SolidMaterial
from Solvers.HeatConduction.Mesh import Mesh2D, Mesh1D
from Solvers.HeatConduction.Boundary import BoundaryRegion
from profiler import TEASAProfiler


# =================================================================
#  1. 导热求解器基类 (BaseHeatConduction)
#  提取 1D 和 2D 共用的状态管理、物性更新和热源接口
# =================================================================
class BaseHeatConduction(ABC):
    VALID_ODE_METHODS = {"RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA", "implicit_euler", "theta_implicit"}
    IMPLICIT_ALGEBRAIC_METHODS = {"implicit_euler", "theta_implicit"}

    def __init__(self,
                 mesh: Union[Mesh1D, Mesh2D],
                 material: SolidMaterial,
                 name: str = "Unnamed_Solid",  # [修改] 新增 name 参数
                 initial_temp: float = 298.15):
        self.mesh = mesh
        self.material = material
        self.name = name
        self.ode_method = "BDF"
        self.theta_implicit_value = 1.0
        self.max_implicit_property_iterations = 2
        self.implicit_property_tol = 1.0e-3
        self.implicit_fallback_method = "BDF"
        self.implicit_fallback_to_solve_ivp = True
        self.minimum_physical_temperature = 0.0
        self.maximum_physical_temperature = np.inf
        self.N = mesh.N  # 节点总数 (1D: n_vols, 2D: nx*ny)

        # --- 状态变量 (Flattened 1D Arrays) ---
        # 无论几何维度如何，为了 ODE 求解器效率，核心状态总是 1D 数组
        self.T = np.full(self.N, initial_temp, dtype=float)
        self.dTdt = np.zeros(self.N, dtype=float)

        # --- 物性缓存 ---
        self.k_node = np.zeros(self.N, dtype=float)
        self.rho_node = np.zeros(self.N, dtype=float)
        self.cp_node = np.zeros(self.N, dtype=float)
        self.thermal_capacitance = np.zeros(self.N, dtype=float)

        # --- 热源 [W] ---
        self.Q_source = np.zeros(self.N, dtype=float)
        self.source_callback: Optional[Callable[[float, np.ndarray], np.ndarray]] = None
        self.use_external_source_buffer = False

        # --- 边界容器 ---
        # 键通常为 'inner'/'outer' (1D) 或 'left'/'right'/'top'/'bottom' (2D)
        self.boundaries: Dict[str, BoundaryRegion] = {}

        # 内部仿真时间记录
        self.current_time = 0.0
        self.last_step_success = True
        self.last_step_failure_message = ""
        self.last_step_failure_time: Optional[float] = None
        self.last_attempted_step_end_time: Optional[float] = None
        self.last_trial_temperature_min = float(initial_temp)
        self.last_trial_temperature_max = float(initial_temp)
        self.last_trial_temperature_time: Optional[float] = None
        self.last_step_used_fallback = False
        self.last_step_fallback_message = ""
        self.last_implicit_outer_iterations = 0

        # 初始化
        self._update_properties()

    def set_ode_method(self, method: str):
        method = str(method)
        if method not in self.VALID_ODE_METHODS:
            valid = ", ".join(sorted(self.VALID_ODE_METHODS))
            raise ValueError(f"Unsupported ODE method '{method}'. Valid methods: {valid}")
        self.ode_method = method

    def _reset_step_diagnostics(self, attempted_end_time: Optional[float] = None):
        current_min = float(np.min(self.T))
        current_max = float(np.max(self.T))
        self.last_step_success = False
        self.last_step_failure_message = ""
        self.last_step_failure_time = None
        self.last_attempted_step_end_time = attempted_end_time
        self.last_trial_temperature_min = current_min
        self.last_trial_temperature_max = current_max
        self.last_trial_temperature_time = self.current_time
        self.last_step_used_fallback = False
        self.last_step_fallback_message = ""
        self.last_implicit_outer_iterations = 0

    def _record_trial_state(self, t: float, T_current: np.ndarray):
        trial_min = float(np.min(T_current))
        trial_max = float(np.max(T_current))

        if trial_min < self.last_trial_temperature_min:
            self.last_trial_temperature_min = trial_min
            self.last_trial_temperature_time = t

        if trial_max > self.last_trial_temperature_max:
            self.last_trial_temperature_max = trial_max

    def _mark_step_success(self):
        self.last_step_success = True
        self.last_step_failure_message = ""
        self.last_step_failure_time = None

    def _mark_step_failure(self, message: str, failure_time: Optional[float] = None):
        self.last_step_success = False
        self.last_step_failure_message = message
        self.last_step_failure_time = self.current_time if failure_time is None else failure_time

    def save_state(self):
        """保存当前温度场，用于全局 Picard 回滚和隐式 fallback。"""
        self.T_backup = self.T.copy()
        self.current_time_backup = self.current_time

    def load_state(self):
        """回滚温度场，并刷新与温度相关的缓存。"""
        self.T[:] = self.T_backup
        self.current_time = self.current_time_backup
        self._update_properties()
        self._compute_internal_resistance()
        self._update_boundaries_state(current_time=self.current_time)

    # 新增的初始化函数
    def initialize_state(self):
        """
        [新增] 强制初始化内部状态

        目的:
        在仿真开始前 (t=0)，确保物性、热阻和边界状态已经基于 initial_temp 完成了计算。
        防止耦合器在第一步读取到 0 度或未初始化的边界数据。
        """
        # 1. 基于初始 T 更新物性 (k, rho, cp)
        self._update_properties()

        # 2. 计算几何热阻 (依赖于刚刚更新的 k)
        self._compute_internal_resistance()

        # 3. 将 T 和 R 推送给 BoundaryRegion
        self._update_boundaries_state(current_time=self.current_time)

        # 4. 让边界区域自我更新一次 T_surface
        self._compute_fluxes(0.0)
        # 如果 BoundaryRegion 有 compute_surface_temperature 之类的方法，应触发它
        # 目前我们主要保证 T_adj_node 和 R_internal 正确即可

    # --- 通用接口: 物性更新 ---
    @TEASAProfiler.profile
    def _update_properties(self):
        """通用物性更新逻辑: 调用 Material 向量化接口"""
        self.k_node[:] = self.material.conductivity(self.T)
        self.rho_node[:] = self.material.density(self.T)
        self.cp_node[:] = self.material.heat_capacity(self.T)

        # 更新热容 (rho * cp * V)
        # 注意: mesh.volumes 在 1D 和 2D 中都应提供 flattened 数组
        if hasattr(self.mesh, 'volumes'):
            # Mesh1D
            vols = self.mesh.volumes
        else:
            # Mesh2D (假设 geom_data.volumes 是 flat 的)
            vols = self.mesh.geom_data.volumes

        self.thermal_capacitance[:] = self.rho_node * self.cp_node * vols

    # --- 通用接口: 热源管理 ---
    def set_source_term(self, source_func: Optional[Callable[[float, np.ndarray], np.ndarray]] = None):
        """设置内热源函数"""
        self.source_callback = source_func
        self.use_external_source_buffer = False

    def link_source_buffer(self, external_buffer: np.ndarray):
        """内存绑定外部热源数组 (高性能耦合用)"""
        if external_buffer.shape != (self.N,):
            raise ValueError(f"External buffer shape {external_buffer.shape} mismatch with mesh N={self.N}")
        self.Q_source = external_buffer
        self.use_external_source_buffer = True
        self.source_callback = None

    def _update_sources(self, t: float):
        """更新热源项"""
        if self.source_callback:
            self.Q_source[:] = self.source_callback(t, self.T)
        elif self.use_external_source_buffer:
            pass  # 外部直接修改内存，无需操作
        else:
            self.Q_source[:] = 0.0

    def get_state_dict(self, prefix: str) -> dict:
        """
        [持久化接口] 将固体核心状态序列化为 NumPy 兼容的扁平字典
        """
        state = {
            # 1. 网格指纹 (用于加载时的绝对安全校验)
            f"{prefix}/shape": np.array(self.T.shape),

            # 2. 核心状态与源项
            f"{prefix}/T": self.T,
            f"{prefix}/dTdt": self.dTdt,
            f"{prefix}/Q_source": self.Q_source,

            # 3. 逻辑开关 (转为 1D 数组存储)
            f"{prefix}/use_ext_buffer": np.array([self.use_external_source_buffer], dtype=bool)
        }

        # 4. 边界状态展开 (将边界对象中的数值数组打散存储)
        for loc, bnd in self.boundaries.items():
            # 这里提取边界中最关键的缓存：表面温度与热流密度
            if hasattr(bnd, 'T_surface'):
                state[f"{prefix}/boundaries/{loc}/T_surface"] = getattr(bnd, 'T_surface')
            if hasattr(bnd, 'q_net'):
                state[f"{prefix}/boundaries/{loc}/q_net"] = getattr(bnd, 'q_net')
            if hasattr(bnd, 'h_conv'):
                state[f"{prefix}/boundaries/{loc}/h_conv"] = getattr(bnd, 'h_conv')

        return state

    def load_state_dict(self, data: dict, prefix: str):
        """
        [持久化接口] 从扁平字典中恢复固体状态，并进行防御性校验
        """
        # --- 1. 致命级网格校验 (Shape Validation) ---
        shape_key = f"{prefix}/shape"
        if shape_key not in data:
            raise KeyError(f"Missing shape fingerprint for {self.name} in restart file.")

        saved_shape = tuple(data[shape_key])
        if saved_shape != self.T.shape:
            raise ValueError(f"CRITICAL: Mesh mismatch for Solid '{self.name}'. "
                             f"File has {saved_shape}, current memory has {self.T.shape}. "
                             f"Cannot restart!")

        # --- 2. 恢复核心状态与源项 ---
        self.T[:] = data[f"{prefix}/T"]
        self.dTdt[:] = data[f"{prefix}/dTdt"]

        # 处理外部源缓冲机制
        use_ext = bool(data[f"{prefix}/use_ext_buffer"][0])
        self.use_external_source_buffer = use_ext
        if not use_ext:
            self.Q_source[:] = data[f"{prefix}/Q_source"]
        # 注：如果 use_ext 为 True，说明 Q_source 被外部系统(如 Pellet)引用，
        # 等待外部宏观组件自己去恢复它，固体底层不应覆盖共享内存。

        # --- 3. 恢复边界状态 ---
        for loc, bnd in self.boundaries.items():
            t_surf_key = f"{prefix}/boundaries/{loc}/T_surface"
            q_net_key = f"{prefix}/boundaries/{loc}/q_net"
            h_conv_key = f"{prefix}/boundaries/{loc}/h_conv"

            if t_surf_key in data and hasattr(bnd, 'T_surface'):
                bnd.T_surface[:] = data[t_surf_key]
            if q_net_key in data and hasattr(bnd, 'q_net'):
                bnd.q_net[:] = data[q_net_key]
            if h_conv_key in data and hasattr(bnd, 'h_conv'):
                bnd.h_conv[:] = data[h_conv_key]

        # --- 4. 强制状态同步 ---
        # 重载了温度后，必须强制刷新底层的物性和内部几何热阻
        self._update_properties()
        self._compute_internal_resistance()
        self._update_boundaries_state(current_time=self.current_time)

    # --- 抽象接口: 核心计算 ---
    @abstractmethod
    def _compute_internal_resistance(self):
        """计算几何热阻 (子类实现)"""
        pass

    @abstractmethod
    def _compute_fluxes(self, t: float) -> np.ndarray:
        """计算净热流 (Q_net) (子类实现)"""
        pass

    @abstractmethod
    def _update_boundaries_state(self, current_time: Optional[float] = None):
        """推送状态到边界对象 (子类实现)"""
        pass

    # --- ODE Solver 标准接口 ---
    @TEASAProfiler.profile
    def get_derivatives(self, t: float, T_current: np.ndarray) -> np.ndarray:
        """计算 dT/dt"""
        # 1. 接收状态
        self.T[:] = T_current
        self._record_trial_state(t, T_current)

        # 2. 更新物性
        self._update_properties()

        # 3. 计算几何热阻 (依赖 k)
        self._compute_internal_resistance()

        # 4. 更新边界状态
        self._update_boundaries_state(current_time=t)

        # 5. 计算净热流 (Flux_in - Flux_out)
        Q_net_conduction = self._compute_fluxes(t)

        # 6. 更新热源
        self._update_sources(t)

        # 7. 能量平衡: dT/dt = (Q_cond + Q_source) / C
        total_Q = Q_net_conduction + self.Q_source

        with np.errstate(divide='ignore'):
            self.dTdt = total_Q / self.thermal_capacitance

        return self.dTdt

    def _solve_ivp_step(self, dt: float, method: str = None, **kwargs) -> bool:
        if method is None or method in self.IMPLICIT_ALGEBRAIC_METHODS:
            method = self.implicit_fallback_method
        method = str(method)
        if method not in self.VALID_ODE_METHODS or method in self.IMPLICIT_ALGEBRAIC_METHODS:
            valid = ", ".join(sorted(self.VALID_ODE_METHODS - self.IMPLICIT_ALGEBRAIC_METHODS))
            raise ValueError(f"Unsupported solve_ivp method '{method}'. Valid solve_ivp methods: {valid}")

        t_span = (self.current_time, self.current_time + dt)
        self._reset_step_diagnostics(t_span[1])

        # 动态检测是否存在预计算的雅可比稀疏矩阵方法 (兼容 1D 和 2D)
        if hasattr(self, 'get_jac_sparsity') and method in ['BDF', 'Radau']:
            # 仅在隐式且支持的情况下注入稀疏模式
            kwargs['jac_sparsity'] = self.get_jac_sparsity()

        # solve_ivp 要求 fun(t, y) -> dy/dt
        # 我们的 get_derivatives 签名正好符合 (self 绑定后)
        sol = solve_ivp(
            fun=self.get_derivatives,
            t_span=t_span,
            y0=self.T.copy(),
            method=method,
            **kwargs
        )

        if sol.success:
            # 取最后一个时间点的解作为当前状态
            self.T[:] = sol.y[:, -1]
            self.current_time += dt

            # 保持步进结束后的边界缓存与最终时刻一致，
            # 这样时间相关边界在步后查询时不会停留在旧时刻。
            self._update_properties()
            self._compute_internal_resistance()
            self._update_boundaries_state(current_time=self.current_time)
            self._compute_fluxes(self.current_time)
            self._mark_step_success()

            # 可选：可以在这里调用一个后处理钩子
            # self.post_step_update()
            return True
        else:
            self._mark_step_failure(sol.message, failure_time=t_span[1])
            print(f"HeatConduction step failed at t={self.current_time}: {sol.message}")
            return False

    def _is_temperature_state_physical(self, temperature_field: np.ndarray) -> bool:
        return (
            np.all(np.isfinite(temperature_field))
            and float(np.min(temperature_field)) >= float(self.minimum_physical_temperature)
            and float(np.max(temperature_field)) <= float(self.maximum_physical_temperature)
        )

    def _fallback_solve_ivp_after_implicit_failure(self, dt: float, message: str, **kwargs) -> bool:
        self.load_state()
        self.last_step_used_fallback = True
        self.last_step_fallback_message = message
        warnings.warn(
            f"Solid '{self.name}' implicit heat solve failed; falling back to {self.implicit_fallback_method}: {message}",
            RuntimeWarning,
            stacklevel=2,
        )
        success = self._solve_ivp_step(dt, method=self.implicit_fallback_method, **kwargs)
        self.last_step_used_fallback = True
        self.last_step_fallback_message = message
        return success

    def _implicit_euler_step(self, dt: float, **kwargs) -> bool:
        return self._theta_implicit_step(dt, theta=1.0, method_name="Implicit Euler", **kwargs)

    def _theta_implicit_step(self, dt: float, theta: float = None, method_name: str = "Theta implicit", **kwargs) -> bool:
        if theta is None:
            theta = self.theta_implicit_value
        theta = float(theta)
        if theta <= 0.0 or theta > 1.0:
            raise ValueError("theta must be in (0, 1]")

        attempted_end_time = self.current_time + dt
        self._reset_step_diagnostics(attempted_end_time)
        self.save_state()

        time_start = self.current_time
        T_start = self.T.copy()
        guess_T = T_start.copy()
        accepted_T = None
        max_outer_iters = max(int(self.max_implicit_property_iterations), 1)
        failure_message = ""

        try:
            for outer_iter in range(max_outer_iters):
                self.last_implicit_outer_iterations = outer_iter + 1
                self.T[:] = guess_T
                self.current_time = time_start
                self._update_properties()
                self._compute_internal_resistance()
                self._update_boundaries_state(current_time=time_start)
                self._update_sources(time_start)

                rhs = self._assemble_theta_implicit_rhs(dt, theta, T_start)
                matrix = self._assemble_theta_implicit_matrix(dt, theta)
                candidate_T = np.asarray(spsolve(matrix, rhs), dtype=float)
                self._record_trial_state(attempted_end_time, candidate_T)

                if not self._is_temperature_state_physical(candidate_T):
                    failure_message = f"{method_name} returned non-finite or unphysical temperature state"
                    accepted_T = None
                    break

                accepted_T = candidate_T
                delta = float(np.max(np.abs(candidate_T - guess_T)))
                scale = max(1.0, float(np.max(np.abs(candidate_T))))
                if delta <= self.implicit_property_tol * scale:
                    break
                guess_T = candidate_T.copy()
        except Exception as exc:
            failure_message = f"{type(exc).__name__}: {exc}"
            accepted_T = None

        if accepted_T is None:
            if self.implicit_fallback_to_solve_ivp:
                return self._fallback_solve_ivp_after_implicit_failure(dt, failure_message or f"{method_name} failed", **kwargs)
            self.load_state()
            self._mark_step_failure(failure_message or f"{method_name} failed", failure_time=attempted_end_time)
            return False

        self.T[:] = accepted_T
        self.current_time = attempted_end_time
        self._update_properties()
        self._compute_internal_resistance()
        self._update_boundaries_state(current_time=self.current_time)
        self._compute_fluxes(self.current_time)
        self._mark_step_success()
        return True

    def _boundary_theta_terms_from_mapping(self, mapping):
        conductance = np.zeros(self.N, dtype=float)
        source = np.zeros(self.N, dtype=float)
        old_flux = np.zeros(self.N, dtype=float)
        explicit_flux = np.zeros(self.N, dtype=float)

        for boundary, flat_indices in mapping:
            indices = np.asarray(flat_indices, dtype=int).reshape(-1)
            if indices.size == 0:
                continue
            boundary_flux = np.asarray(boundary.compute_net_flux_for_solver(), dtype=float).reshape(-1)
            if boundary_flux.size != indices.size:
                raise ValueError(
                    f"Boundary flux size {boundary_flux.size} does not match mapped node count {indices.size}"
                )

            resistance_mask = getattr(boundary, "_resistance_mask", None)
            if resistance_mask is None:
                G_sum = getattr(boundary, "G_sum", None)
                if G_sum is None:
                    np.add.at(explicit_flux, indices, boundary_flux)
                    continue
                resistance_mask = np.asarray(G_sum, dtype=float).reshape(-1) > 1.0e-20
            else:
                resistance_mask = np.asarray(resistance_mask, dtype=bool).reshape(-1)

            if resistance_mask.size != indices.size:
                raise ValueError("Boundary resistance mask size does not match mapped node count")

            if np.any(resistance_mask):
                R_eff = np.asarray(boundary.R_eff, dtype=float).reshape(-1)
                R_internal = np.asarray(boundary.R_internal, dtype=float).reshape(-1)
                T_eff = np.asarray(boundary.T_eff, dtype=float).reshape(-1)
                R_total = R_eff + R_internal
                with np.errstate(divide='ignore', invalid='ignore'):
                    G_total = np.divide(
                        1.0,
                        R_total,
                        out=np.zeros_like(R_total, dtype=float),
                        where=resistance_mask,
                    )
                G_total = np.nan_to_num(G_total, nan=0.0, posinf=1.0e20, neginf=0.0)
                active_indices = indices[resistance_mask]
                active_G = G_total[resistance_mask]
                np.add.at(conductance, active_indices, active_G)
                np.add.at(source, active_indices, active_G * T_eff[resistance_mask])
                np.add.at(old_flux, active_indices, boundary_flux[resistance_mask])

            if np.any(~resistance_mask):
                np.add.at(explicit_flux, indices[~resistance_mask], boundary_flux[~resistance_mask])

        return conductance, source, old_flux, explicit_flux

    def _build_csc_template_from_entries(self, rows: np.ndarray, cols: np.ndarray):
        template = coo_matrix(
            (np.zeros(rows.size, dtype=float), (rows, cols)),
            shape=(self.N, self.N),
        ).tocsc()
        position_by_entry = {}
        for col in range(self.N):
            start = template.indptr[col]
            end = template.indptr[col + 1]
            for pos in range(start, end):
                position_by_entry[(int(template.indices[pos]), col)] = pos
        entry_positions = np.array(
            [position_by_entry[(int(row), int(col))] for row, col in zip(rows, cols)],
            dtype=int,
        )
        return template, entry_positions

    def _assemble_theta_implicit_rhs(self, dt: float, theta: float, T_start: np.ndarray) -> np.ndarray:
        inv_dt = 1.0 / max(float(dt), 1.0e-30)
        boundary_conductance, boundary_source, boundary_old_flux, boundary_explicit_flux = self._compute_boundary_theta_terms()
        rhs = self.thermal_capacitance * inv_dt * T_start
        if theta < 1.0:
            rhs = rhs + (1.0 - theta) * self._compute_internal_flux_for_temperature(T_start)
            rhs = rhs + (1.0 - theta) * boundary_old_flux
        rhs = rhs + theta * boundary_source + boundary_explicit_flux + self.Q_source
        self._last_boundary_conductance_for_matrix = boundary_conductance
        return rhs

    @abstractmethod
    def _compute_boundary_theta_terms(self):
        pass

    @abstractmethod
    def _compute_internal_flux_for_temperature(self, T_flat: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def _assemble_theta_implicit_matrix(self, dt: float, theta: float):
        pass

    # 时间步进接口，用于推进物体温度场
    @TEASAProfiler.profile
    def step(self, dt: float, method: str = None, **kwargs) -> bool:
        """
        执行一个时间步长的积分，推进固体温度场。
        """
        if method is None:
            method = self.ode_method
        else:
            method = str(method)
            if method not in self.VALID_ODE_METHODS:
                valid = ", ".join(sorted(self.VALID_ODE_METHODS))
                raise ValueError(f"Unsupported ODE method '{method}'. Valid methods: {valid}")

        if method == "implicit_euler":
            self.save_state()
            success = self._implicit_euler_step(dt, **kwargs)
            if not success and self.implicit_fallback_to_solve_ivp and not self.last_step_used_fallback:
                return self._fallback_solve_ivp_after_implicit_failure(dt, "implicit_euler returned False", **kwargs)
            return success
        if method == "theta_implicit":
            self.save_state()
            success = self._theta_implicit_step(dt, theta=self.theta_implicit_value, **kwargs)
            if not success and self.implicit_fallback_to_solve_ivp and not self.last_step_used_fallback:
                return self._fallback_solve_ivp_after_implicit_failure(dt, "theta_implicit returned False", **kwargs)
            return success
        return self._solve_ivp_step(dt, method=method, **kwargs)


# =================================================================
#  2. 一维导热求解器 (HeatConduction1D)
#  修正说明：BoundaryRegion 初始化移除 name 参数
# =================================================================

class HeatConduction1D(BaseHeatConduction):
    """
    一维瞬态导热求解器
    """

    def __init__(self,
                 mesh: Mesh1D,
                 material: SolidMaterial,
                 initial_temp: float = 298.15,
                 name: str = "Unnamed_Solid"):  # [修改] 新增 name 参数):
        super().__init__(mesh=mesh, material=material, name=name, initial_temp=initial_temp)

        # --- 1D 特有结构 ---
        # 几何导热热阻 R_geom (Size N+1)
        # 包含: 0=左边界, 1..N-1=内部界面, N=右边界
        self.R_geom = np.zeros(self.N + 1, dtype=float)

        # 界面热流 (Size N+1)
        self.interface_flux = np.zeros(self.N + 1, dtype=float)

        # --- 边界初始化 (修正处) ---
        # 严格参照 Boundary.py，__init__ 只有 shape 和 area_array

        # Left / Inner Boundary (Index 0)
        self.boundaries['inner'] = BoundaryRegion(
            shape=(1,),
            area_array=np.array([self.mesh.face_areas[0]])
        )

        # Right / Outer Boundary (Index N)
        self.boundaries['outer'] = BoundaryRegion(
            shape=(1,),
            area_array=np.array([self.mesh.face_areas[-1]])
        )

        # [新增] 初始化边界状态！
        # 这一步至关重要，确保边界持有正确的 T_init 和 R_int
        self.initialize_state()

    def attach_boundary(self, location: str, boundary_region: BoundaryRegion):
        """
        [兼容接口] 允许外部替换边界对象
        """
        if location not in ['inner', 'outer']:
            raise ValueError("Location must be 'inner' or 'outer'")

        if boundary_region.shape != (1,):
            raise ValueError("BoundaryRegion for 1D solver must have shape (1,)")

        self.boundaries[location] = boundary_region

    def _compute_internal_resistance(self):
        """
        [Step 2] 计算内部导热热阻 R_geom
        **完全保留您原代码中关于圆柱/笛卡尔坐标的处理逻辑**
        """
        # 1. 基础数据准备
        geo_type = self.mesh.geometry_type
        H = self.mesh.H  # 高度/深度

        # 导热系数 (节点间调和平均)
        k_L = self.k_node[:-1]
        k_R = self.k_node[1:]

        with np.errstate(divide='ignore', invalid='ignore'):
            k_inter = 2 * k_L * k_R / (k_L + k_R)
            k_inter = np.nan_to_num(k_inter, nan=0.0)

        # ==========================================================
        # 分支 A: 圆柱坐标 (使用精确对数热阻)
        # ==========================================================
        if geo_type == 'cylindrical':
            r_nodes = self.mesh.node_centers
            r_faces = self.mesh.face_locations

            r_i = r_nodes[:-1]  # 左节点
            r_next = r_nodes[1:]  # 右节点

            # --- 1. 内部节点间热阻 ---
            r_i_safe = np.maximum(r_i, 1e-12)
            self.R_geom[1:-1] = np.log(r_next / r_i_safe) / (2.0 * np.pi * k_inter * H)

            # --- 2. 边界热阻 (Face -> Node) ---
            # Inner Boundary (Face 0 -> Node 0)
            r_f0 = max(r_faces[0], 1e-12)
            r_n0 = max(r_nodes[0], 1e-12)

            if abs(r_n0 - r_f0) < 1e-12:
                self.R_geom[0] = 0.0
            else:
                self.R_geom[0] = np.abs(np.log(r_n0 / r_f0)) / (2.0 * np.pi * self.k_node[0] * H)

            # Outer Boundary (Node N-1 -> Face N)
            r_fN = max(r_faces[-1], 1e-12)
            r_nN = max(r_nodes[-1], 1e-12)
            self.R_geom[-1] = np.abs(np.log(r_fN / r_nN)) / (2.0 * np.pi * self.k_node[-1] * H)

        # ==========================================================
        # 分支 B: 笛卡尔坐标 (使用线性热阻)
        # ==========================================================
        else:
            # --- 1. 内部节点间热阻 ---
            dx = self.mesh.dr_node_to_node
            area = self.mesh.face_areas[1:-1]
            self.R_geom[1:-1] = dx / (k_inter * area)

            # --- 2. 边界热阻 ---
            # Inner (Left)
            dr_0 = self.mesh.dr_node_to_face[0, 0]
            A_0 = self.mesh.face_areas[0]
            self.R_geom[0] = dr_0 / (self.k_node[0] * A_0)

            # Outer (Right)
            dr_N = self.mesh.dr_node_to_face[-1, 1]
            A_N = self.mesh.face_areas[-1]
            self.R_geom[-1] = dr_N / (self.k_node[-1] * A_N)

    def _update_boundaries_state(self, current_time: Optional[float] = None):
        """
        [Step 3] 将内部状态推送给边界对象
        适配：传递 (1,) 维度的数组
        """
        # Inner Boundary
        if self.boundaries['inner']:
            self.boundaries['inner'].update_internal_state(
                T_node=self.T[0:1],  # 切片保持维度 (1,)
                R_int=self.R_geom[0:1],
                current_time=current_time
            )

        # Outer Boundary
        if self.boundaries['outer']:
            self.boundaries['outer'].update_internal_state(
                T_node=self.T[-1:],
                R_int=self.R_geom[-1:],
                current_time=current_time
            )

    def _compute_fluxes(self, t: float) -> np.ndarray:
        """
        [Step 4] 计算所有界面的热流
        """
        # --- A. 内部通量 (i -> i+1) ---
        delta_T = self.T[:-1] - self.T[1:]

        with np.errstate(divide='ignore', invalid='ignore'):
            flux = delta_T / self.R_geom[1:-1]
            # 处理可能的无效值
            self.interface_flux[1:-1] = np.nan_to_num(flux)

        # --- B. 边界通量 ---
        # 符号约定: interface_flux 正方向为 x+ (Left -> Right)

        # Inner (Left): Boundary 返回 "流入节点" (Ext -> Node 0, 即 x+)
        if self.boundaries['inner']:
            q_in = self.boundaries['inner'].compute_net_flux_for_solver()
            self.interface_flux[0] = q_in[0]  # 取标量
        else:
            self.interface_flux[0] = 0.0

        # Outer (Right): Boundary 返回 "流入节点" (Ext -> Node N-1, 即 x-)
        # interface_flux[-1] 定义为流出 (Node N-1 -> Ext, 即 x+)
        # 所以 Flux_out = - Flux_in
        if self.boundaries['outer']:
            q_in_outer = self.boundaries['outer'].compute_net_flux_for_solver()
            self.interface_flux[-1] = -q_in_outer[0]
        else:
            self.interface_flux[-1] = 0.0

        # 返回净导热项 Q_net [W]
        # Q_i = Flux_in - Flux_out
        return self.interface_flux[:-1] - self.interface_flux[1:]

    def _compute_boundary_theta_terms(self):
        mapping = []
        if 'inner' in self.boundaries:
            mapping.append((self.boundaries['inner'], np.array([0], dtype=int)))
        if 'outer' in self.boundaries:
            mapping.append((self.boundaries['outer'], np.array([self.N - 1], dtype=int)))
        return self._boundary_theta_terms_from_mapping(mapping)

    def _compute_internal_flux_for_temperature(self, T_flat: np.ndarray) -> np.ndarray:
        T_arr = np.asarray(T_flat, dtype=float).reshape(-1)
        flux = np.zeros(self.N + 1, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            flux[1:-1] = np.nan_to_num((T_arr[:-1] - T_arr[1:]) / self.R_geom[1:-1])
        return flux[:-1] - flux[1:]

    def _assemble_theta_implicit_matrix(self, dt: float, theta: float):
        inv_dt = 1.0 / max(float(dt), 1.0e-30)
        boundary_conductance = getattr(
            self,
            "_last_boundary_conductance_for_matrix",
            np.zeros(self.N, dtype=float),
        )
        if not hasattr(self, "_implicit_matrix_cache") or self._implicit_matrix_cache.get("N") != self.N:
            if self.N > 1:
                left = np.arange(self.N - 1, dtype=int)
                right = left + 1
                diag_idx = np.arange(self.N, dtype=int)
                rows = np.concatenate((left, right, diag_idx))
                cols = np.concatenate((right, left, diag_idx))
                data = np.empty(rows.size, dtype=float)
                matrix, entry_positions = self._build_csc_template_from_entries(rows, cols)
                cache = {
                    "N": self.N,
                    "rows": rows,
                    "cols": cols,
                    "data": data,
                    "matrix": matrix,
                    "entry_positions": entry_positions,
                    "left": left,
                    "right": right,
                    "lr_slice": slice(0, self.N - 1),
                    "rl_slice": slice(self.N - 1, 2 * (self.N - 1)),
                    "diag_slice": slice(2 * (self.N - 1), 2 * (self.N - 1) + self.N),
                    "diag_work": np.empty(self.N, dtype=float),
                    "conductance_work": np.empty(self.N - 1, dtype=float),
                }
            else:
                rows = np.array([0], dtype=int)
                cols = np.array([0], dtype=int)
                matrix, entry_positions = self._build_csc_template_from_entries(rows, cols)
                cache = {
                    "N": self.N,
                    "rows": rows,
                    "cols": cols,
                    "data": np.empty(1, dtype=float),
                    "matrix": matrix,
                    "entry_positions": entry_positions,
                    "left": np.array([], dtype=int),
                    "right": np.array([], dtype=int),
                    "lr_slice": slice(0, 0),
                    "rl_slice": slice(0, 0),
                    "diag_slice": slice(0, 1),
                    "diag_work": np.empty(1, dtype=float),
                    "conductance_work": np.empty(0, dtype=float),
                }
            self._implicit_matrix_cache = cache

        cache = self._implicit_matrix_cache
        diag = cache["diag_work"]
        diag[:] = self.thermal_capacitance * inv_dt
        diag += float(theta) * np.asarray(boundary_conductance, dtype=float)
        data = cache["data"]

        if self.N > 1:
            cache["conductance_work"].fill(0.0)
            with np.errstate(divide='ignore', invalid='ignore'):
                conductances = np.divide(
                    1.0,
                    self.R_geom[1:-1],
                    out=cache["conductance_work"],
                    where=np.abs(self.R_geom[1:-1]) > 1.0e-300,
                )
            conductances = cache["conductance_work"]
            np.nan_to_num(conductances, copy=False, nan=0.0, posinf=1.0e20, neginf=0.0)
            conductances *= float(theta)
            np.add.at(diag, cache["left"], conductances)
            np.add.at(diag, cache["right"], conductances)
            data[cache["lr_slice"]] = -conductances
            data[cache["rl_slice"]] = -conductances

        data[cache["diag_slice"]] = diag
        cache["matrix"].data[cache["entry_positions"]] = data
        return cache["matrix"]

    # 在 HeatConduction1D 类中添加
    def get_boundary_node_capacitance(self, location: str) -> np.ndarray:
        """
        [接口] 获取指定边界处的节点热容数组 [J/K]
        """
        if location == 'inner':
            # 第一个节点
            return self.thermal_capacitance[0:1].copy()
        elif location == 'outer':
            # 最后一个节点
            return self.thermal_capacitance[-1:].copy()
        else:
            raise ValueError(f"Invalid boundary location for 1D: {location}")

    def post_step_update(self):
        """时间步后处理，保留空接口"""
        pass


# =================================================================
#  3. 二维导热求解器 (HeatConduction2D)
#  严格适配无 name 参数的 BoundaryRegion
# =================================================================

class HeatConduction2D(BaseHeatConduction):
    """
    二维瞬态导热求解器 (支持 Cartesian & Cylindrical via Mesh2D)
    """

    def __init__(self,
                 mesh: Mesh2D,
                 material: SolidMaterial,
                 name: str = "Unnamed_Solid",  # [修改] 新增 name 参数
                 initial_temp: float = 298.15):
        super().__init__(mesh=mesh, material=material, name=name, initial_temp=initial_temp)

        if name is not None:
            self.name = name

        # --- 2D 特有缓存 ---
        # 缓存网格形状视图，避免重复调用 reshape
        self.shape_nodes = self.mesh.shape_nodes  # (nx, ny)

        # [新增] 预分配 _compute_fluxes 的累加数组，避免每秒发生数万次内存分配与垃圾回收
        self.Q_net_2d_buffer = np.zeros(self.shape_nodes, dtype=float)

        # 内部界面热阻缓存 (在 _compute_internal_resistance 中更新)
        # R_x_inner: X方向内部界面热阻, Size (nx-1, ny)
        # R_y_inner: Y方向内部界面热阻, Size (nx, ny-1)
        # self.R_x_inner = None
        # self.R_y_inner = None

        # [修改] 内部界面热导缓存 (使用 Conductance = 1 / Resistance)
        # G_x_inner: X方向内部界面热导, Size (nx-1, ny)
        # G_y_inner: Y方向内部界面热导, Size (nx, ny-1)
        nx, ny = self.shape_nodes
        self.G_x_inner = np.empty((max(nx - 1, 0), ny), dtype=float)
        self.G_y_inner = np.empty((nx, max(ny - 1, 0)), dtype=float)
        self._conductance_work_x = np.empty_like(self.G_x_inner)
        self._conductance_work_y = np.empty_like(self.G_y_inner)
        self._flux_x_buffer = np.empty_like(self.G_x_inner)
        self._flux_y_buffer = np.empty_like(self.G_y_inner)
        self._jac_sparsity_cache = None

        # --- 边界初始化 (修正: 无 name 参数) ---
        # 必须确保 BoundaryRegion 接收正确的 shape 和 area_array
        # 使用 .copy() 防止意外修改网格原始数据

        # 1. Left / Inner Boundary (x=0) -> Shape: (ny,)
        self.boundaries['left'] = BoundaryRegion(
            shape=(self.mesh.n_y,),
            area_array=self.mesh.area_x_matrix[0, :].copy()
        )

        # 2. Right / Outer Boundary (x=nx) -> Shape: (ny,)
        self.boundaries['right'] = BoundaryRegion(
            shape=(self.mesh.n_y,),
            area_array=self.mesh.area_x_matrix[-1, :].copy()
        )

        # 3. Bottom Boundary (y=0) -> Shape: (nx,)
        self.boundaries['bottom'] = BoundaryRegion(
            shape=(self.mesh.n_x,),
            area_array=self.mesh.area_y_matrix[:, 0].copy()
        )

        # 4. Top Boundary (y=ny) -> Shape: (nx,)
        self.boundaries['top'] = BoundaryRegion(
            shape=(self.mesh.n_x,),
            area_array=self.mesh.area_y_matrix[:, -1].copy()
        )

        # [新增] 初始化边界状态！
        # 这一步至关重要，确保边界持有正确的 T_init 和 R_int
        self.initialize_state()

    def get_jac_sparsity(self):
        """
        [TEASA 优化接口] 获取 2D 导热雅可比矩阵的稀疏模式 (Jacobian Sparsity)。
        根据自变量 T 扁平化重塑为 (nx, ny) 的内存布局，
        生成五点差分格式的稀疏对角矩阵。

        用于隐式求解器 (如 BDF, Radau) 时，将极其昂贵的 O(N^3) 稠密矩阵操作
        以及 O(N^2) 的有限差分扰动次数，降维至 O(N) 级别。

        动态生成差分格式的稀疏对角矩阵，完美适配降维退化 (1D/0D)。
        """
        if self._jac_sparsity_cache is not None:
            return self._jac_sparsity_cache

        N = self.N
        nx, ny = self.shape_nodes

        offsets = []
        data_rows = []

        # 1. 永远存在的主对角线 (自身节点)
        offsets.append(0)
        data_rows.append(np.ones(N))

        # 2. 只有当 Y 方向 (Axial) 网格数 > 1 时，才存在 Y 方向的真实物理邻居 (偏移量为 ±1)
        if ny > 1:
            offsets.extend([-1, 1])
            data_rows.extend([np.ones(N), np.ones(N)])

        # 3. 只有当 X 方向 (Radial) 网格数 > 1 时，才存在 X 方向的真实物理邻居 (偏移量为 ±ny)
        if nx > 1:
            offsets.extend([-ny, ny])
            data_rows.extend([np.ones(N), np.ones(N)])

        # 将 list 转换为 numpy 数组结构，满足 spdiags 的输入要求
        offsets = np.array(offsets)
        data = np.vstack(data_rows)

        # 使用 scipy.sparse 生成稀疏矩阵 (csc 格式在 ODE 求解器中效率较高)
        self._jac_sparsity_cache = spdiags(data, offsets, N, N, format='csc')

        return self._jac_sparsity_cache

        # 旧版，无法处理 1D/0D 的导热网格模型
        # N = self.N
        # nx, ny = self.shape_nodes
        #
        # # 2D 网格中，一个节点 (i, j) 最多依赖 5 个节点：
        # # 自身 (0), Y方向相邻 (±1), X方向相邻 (±ny)
        # offsets = np.array([0, -1, 1, -ny, ny])
        #
        # # 构造对角线数据，用 1 标记可能非零的依赖关系
        # data = np.ones((5, N))
        #
        # # 使用 scipy.sparse 生成稀疏矩阵 (csc 格式在 ODE 求解器中效率较高)
        # jac_sparsity = spdiags(data, offsets, N, N, format='csc')
        #
        # return jac_sparsity

    @TEASAProfiler.profile
    def _compute_internal_resistance(self):
        """
        [Step 2] 计算内部几何热阻
        R = distance / (k_face * area)
        [更新版本] 直接计算热导，解决 k = 0 时的除零错误
        """
        # 创建 2D 视图
        k_2d = self.k_node.reshape(self.shape_nodes)

        # === X 方向 (Radial/Horizontal) ===
        # 内部界面: indices 1 to nx-1
        k_left = k_2d[:-1, :]
        k_right = k_2d[1:, :]

        # 界面导热系数 (调和平均)
        # [优化点] 分母加上 1e-30，绝对防止 0/0 报错，无需使用 np.errstate 和 np.nan_to_num
        np.add(k_left, k_right, out=self._conductance_work_x)
        self._conductance_work_x += 1e-30
        np.multiply(k_left, k_right, out=self.G_x_inner)
        self.G_x_inner *= 2.0
        np.divide(self.G_x_inner, self._conductance_work_x, out=self.G_x_inner)

        dx_inner = self.mesh.dx_matrix[1:-1, :]
        # [优化点] 防止 dx 为 0 导致除零
        area_x_inner = self.mesh.area_x_matrix[1:-1, :]

        # 计算热导 G_x (乘法优先，最后做一次除法)
        np.maximum(dx_inner, 1e-30, out=self._conductance_work_x)
        np.multiply(self.G_x_inner, area_x_inner, out=self.G_x_inner)
        np.divide(self.G_x_inner, self._conductance_work_x, out=self.G_x_inner)

        # === Y 方向 (Axial/Vertical) ===
        # 内部界面: indices 1 to ny-1
        k_bottom = k_2d[:, :-1]
        k_top = k_2d[:, 1:]

        # 界面导热系数 (调和平均)
        np.add(k_bottom, k_top, out=self._conductance_work_y)
        self._conductance_work_y += 1e-30
        np.multiply(k_bottom, k_top, out=self.G_y_inner)
        self.G_y_inner *= 2.0
        np.divide(self.G_y_inner, self._conductance_work_y, out=self.G_y_inner)

        dy_inner = self.mesh.dy_matrix[:, 1:-1]
        # [优化点] 防止 dy 为 0 导致除零
        area_y_inner = self.mesh.area_y_matrix[:, 1:-1]

        # 计算热导 G_y
        np.maximum(dy_inner, 1e-30, out=self._conductance_work_y)
        np.multiply(self.G_y_inner, area_y_inner, out=self.G_y_inner)
        np.divide(self.G_y_inner, self._conductance_work_y, out=self.G_y_inner)

    def _update_boundaries_state(self, current_time: Optional[float] = None):
        """
        [Step 3] 推送状态到边界
        注意：需要临时计算边界处的导热热阻 R_int = d / (k * A)
        """
        # 视图引用
        T_2d = self.T.reshape(self.shape_nodes)
        k_2d = self.k_node.reshape(self.shape_nodes)

        # 获取几何数据
        dx_mat = self.mesh.dx_matrix  # (nx+1, ny)
        dy_mat = self.mesh.dy_matrix  # (nx, ny+1)
        ax_mat = self.mesh.area_x_matrix  # (nx+1, ny)
        ay_mat = self.mesh.area_y_matrix  # (nx, ny+1)

        # [辅助函数] 快速计算边界热阻，加 1e-30 避免除以零，丢弃 with np.errstate
        def calc_R_int_safe(d, k, A):
            denom = k * A
            denom = np.maximum(denom, 1e-30)
            return d / denom

        # --- 1. Left (x=0) ---
        self.boundaries['left'].update_internal_state(
            T_2d[0, :],
            calc_R_int_safe(dx_mat[0, :], k_2d[0, :], ax_mat[0, :]),
            current_time=current_time
        )

        # --- 2. Right (x=nx) ---
        self.boundaries['right'].update_internal_state(
            T_2d[-1, :],
            calc_R_int_safe(dx_mat[-1, :], k_2d[-1, :], ax_mat[-1, :]),
            current_time=current_time
        )

        # --- 3. Bottom (y=0) ---
        self.boundaries['bottom'].update_internal_state(
            T_2d[:, 0],
            calc_R_int_safe(dy_mat[:, 0], k_2d[:, 0], ay_mat[:, 0]),
            current_time=current_time
        )

        # --- 4. Top (y=ny) ---
        self.boundaries['top'].update_internal_state(
            T_2d[:, -1],
            calc_R_int_safe(dy_mat[:, -1], k_2d[:, -1], ay_mat[:, -1]),
            current_time=current_time
        )

    @TEASAProfiler.profile
    def _compute_fluxes(self, t: float) -> np.ndarray:
        """
        [Step 4] 计算净热流 Q_net
        [优化版本] 无内存分配，无 np.errstate，使用热导(G)将除法变为乘法。
        """
        # [优化点 1] 复用预分配的内存，代替高频调用的 np.zeros()
        Q_net_2d = self.Q_net_2d_buffer
        Q_net_2d.fill(0.0)

        T_2d = self.T.reshape(self.shape_nodes)

        # === A. 内部通量 (Internal Fluxes) ===
        # [优化点 2] 使用热导 (G_x_inner) 进行乘法运算，不再有除以 0 的风险

        # X 方向
        np.subtract(T_2d[:-1, :], T_2d[1:, :], out=self._flux_x_buffer)
        np.multiply(self._flux_x_buffer, self.G_x_inner, out=self._flux_x_buffer)
        Q_net_2d[:-1, :] -= self._flux_x_buffer
        Q_net_2d[1:, :] += self._flux_x_buffer

        # Y 方向
        np.subtract(T_2d[:, :-1], T_2d[:, 1:], out=self._flux_y_buffer)
        np.multiply(self._flux_y_buffer, self.G_y_inner, out=self._flux_y_buffer)
        Q_net_2d[:, :-1] -= self._flux_y_buffer
        Q_net_2d[:, 1:] += self._flux_y_buffer

        # === B. 边界通量 (Boundary Fluxes) ===
        # BoundaryRegion 返回的是 "流入节点" (Inflow to Node)
        Q_net_2d[0, :] += self.boundaries['left'].compute_net_flux_for_solver()
        Q_net_2d[-1, :] += self.boundaries['right'].compute_net_flux_for_solver()
        Q_net_2d[:, 0] += self.boundaries['bottom'].compute_net_flux_for_solver()
        Q_net_2d[:, -1] += self.boundaries['top'].compute_net_flux_for_solver()

        # 返回 Flatten 后的净热流数组
        # 注意: flatten() 默认会返回一个深拷贝 (copy)，所以它不会破坏 buffer 的内存
        return Q_net_2d.flatten()

    def _compute_boundary_theta_terms(self):
        nx, ny = self.shape_nodes
        mapping = []
        if 'left' in self.boundaries:
            mapping.append((self.boundaries['left'], np.arange(ny, dtype=int)))
        if 'right' in self.boundaries:
            mapping.append((self.boundaries['right'], (nx - 1) * ny + np.arange(ny, dtype=int)))
        if 'bottom' in self.boundaries:
            mapping.append((self.boundaries['bottom'], np.arange(nx, dtype=int) * ny))
        if 'top' in self.boundaries:
            mapping.append((self.boundaries['top'], np.arange(nx, dtype=int) * ny + (ny - 1)))
        return self._boundary_theta_terms_from_mapping(mapping)

    def _compute_internal_flux_for_temperature(self, T_flat: np.ndarray) -> np.ndarray:
        Q_net_2d = self.Q_net_2d_buffer
        Q_net_2d.fill(0.0)
        T_2d = np.asarray(T_flat, dtype=float).reshape(self.shape_nodes)

        if self.G_x_inner.size > 0:
            np.subtract(T_2d[:-1, :], T_2d[1:, :], out=self._flux_x_buffer)
            np.multiply(self._flux_x_buffer, self.G_x_inner, out=self._flux_x_buffer)
            Q_net_2d[:-1, :] -= self._flux_x_buffer
            Q_net_2d[1:, :] += self._flux_x_buffer

        if self.G_y_inner.size > 0:
            np.subtract(T_2d[:, :-1], T_2d[:, 1:], out=self._flux_y_buffer)
            np.multiply(self._flux_y_buffer, self.G_y_inner, out=self._flux_y_buffer)
            Q_net_2d[:, :-1] -= self._flux_y_buffer
            Q_net_2d[:, 1:] += self._flux_y_buffer

        return Q_net_2d.reshape(-1).copy()

    def _assemble_theta_implicit_matrix(self, dt: float, theta: float):
        nx, ny = self.shape_nodes
        inv_dt = 1.0 / max(float(dt), 1.0e-30)
        boundary_conductance = getattr(
            self,
            "_last_boundary_conductance_for_matrix",
            np.zeros(self.N, dtype=float),
        )
        cache_sig = (self.N, nx, ny)
        if not hasattr(self, "_implicit_matrix_cache") or self._implicit_matrix_cache.get("sig") != cache_sig:
            row_parts = []
            col_parts = []
            data_sizes = []
            cache = {
                "sig": cache_sig,
                "x_left": np.array([], dtype=int),
                "x_right": np.array([], dtype=int),
                "y_bottom": np.array([], dtype=int),
                "y_top": np.array([], dtype=int),
            }

            if self.G_x_inner.size > 0:
                i_idx, j_idx = np.meshgrid(np.arange(nx - 1), np.arange(ny), indexing='ij')
                x_left = (i_idx * ny + j_idx).reshape(-1)
                x_right = x_left + ny
                cache["x_left"] = x_left
                cache["x_right"] = x_right
                row_parts.extend((x_left, x_right))
                col_parts.extend((x_right, x_left))
                data_sizes.extend((x_left.size, x_right.size))

            if self.G_y_inner.size > 0:
                i_idx, j_idx = np.meshgrid(np.arange(nx), np.arange(ny - 1), indexing='ij')
                y_bottom = (i_idx * ny + j_idx).reshape(-1)
                y_top = y_bottom + 1
                cache["y_bottom"] = y_bottom
                cache["y_top"] = y_top
                row_parts.extend((y_bottom, y_top))
                col_parts.extend((y_top, y_bottom))
                data_sizes.extend((y_bottom.size, y_top.size))

            diag_idx = np.arange(self.N, dtype=int)
            row_parts.append(diag_idx)
            col_parts.append(diag_idx)
            data_sizes.append(self.N)

            rows = np.concatenate(row_parts)
            cols = np.concatenate(col_parts)
            data = np.empty(rows.size, dtype=float)
            matrix, entry_positions = self._build_csc_template_from_entries(rows, cols)
            cursor = 0
            slice_names = []
            if self.G_x_inner.size > 0:
                slice_names.extend(("x_lr_slice", "x_rl_slice"))
            if self.G_y_inner.size > 0:
                slice_names.extend(("y_bt_slice", "y_tb_slice"))
            slice_names.append("diag_slice")
            for name, size in zip(slice_names, data_sizes):
                cache[name] = slice(cursor, cursor + size)
                cursor += size
            cache.update({
                "rows": rows,
                "cols": cols,
                "data": data,
                "matrix": matrix,
                "entry_positions": entry_positions,
                "diag_work": np.empty(self.N, dtype=float),
                "x_conductance_work": np.empty(self.G_x_inner.size, dtype=float),
                "y_conductance_work": np.empty(self.G_y_inner.size, dtype=float),
            })
            self._implicit_matrix_cache = cache

        cache = self._implicit_matrix_cache
        diag = cache["diag_work"]
        diag[:] = self.thermal_capacitance * inv_dt
        diag += float(theta) * np.asarray(boundary_conductance, dtype=float)
        data = cache["data"]

        if self.G_x_inner.size > 0:
            conductance = cache["x_conductance_work"]
            conductance[:] = self.G_x_inner.reshape(-1)
            conductance *= float(theta)
            np.add.at(diag, cache["x_left"], conductance)
            np.add.at(diag, cache["x_right"], conductance)
            data[cache["x_lr_slice"]] = -conductance
            data[cache["x_rl_slice"]] = -conductance

        if self.G_y_inner.size > 0:
            conductance = cache["y_conductance_work"]
            conductance[:] = self.G_y_inner.reshape(-1)
            conductance *= float(theta)
            np.add.at(diag, cache["y_bottom"], conductance)
            np.add.at(diag, cache["y_top"], conductance)
            data[cache["y_bt_slice"]] = -conductance
            data[cache["y_tb_slice"]] = -conductance

        data[cache["diag_slice"]] = diag
        cache["matrix"].data[cache["entry_positions"]] = data
        return cache["matrix"]

    def get_boundary_node_capacitance(self, location: str) -> np.ndarray:
        """
        [接口] 获取指定边界处的节点热容数组 [J/K]
        用于耦合器 (FluidSolidCouple) 计算稳定性时间步长。
        """
        # 1. 检查 location 是否合法
        if location not in ['left', 'right', 'top', 'bottom']:
            raise ValueError(f"Invalid boundary location for 2D: {location}")

        # 2. 将扁平化的热容数组重塑为 2D 视图 (nx, ny)
        # 注意: self.thermal_capacitance 是 BaseHeatConduction 维护的
        cap_2d = self.thermal_capacitance.reshape(self.shape_nodes)

        # 3. 根据位置切片提取
        # 返回 copy() 以防止外部修改影响内部状态
        if location == 'left':
            # x=0 这一列 (对应 mesh.area_x_matrix[0, :])
            return cap_2d[0, :].copy()

        elif location == 'right':
            # x=nx-1 这一列
            return cap_2d[-1, :].copy()

        elif location == 'bottom':
            # y=0 这一行
            return cap_2d[:, 0].copy()

        elif location == 'top':
            # y=ny-1 这一行
            return cap_2d[:, -1].copy()

    def save_state(self):
        """保存当前温度场"""
        # 使用 copy() 深拷贝，防止引用被修改
        self.T_backup = self.T.copy()
        self.current_time_backup = self.current_time

    def load_state(self):
        """回滚温度场"""
        self.T[:] = self.T_backup
        self.current_time = self.current_time_backup

        # 回滚后，立即基于旧温度刷新物性和边界状态
        self._update_properties()
        self._compute_internal_resistance()
        self._update_boundaries_state(current_time=self.current_time)
