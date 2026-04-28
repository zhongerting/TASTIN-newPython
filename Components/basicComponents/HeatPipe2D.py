import numpy as np
from scipy.integrate import solve_ivp

from Materials.Base import SolidMaterial
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.Boundary import BoundaryRegion
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D


class HeatPipe2D(HeatConduction2D):
    """
    二维柱坐标热管纯导热求解器。
    - 径向区分 wick / wall
    - 轴向外壁切分为蒸发段、绝热段、冷凝段
    """

    def __init__(self,
                 mesh: Mesh2D,
                 solid1: SolidMaterial,
                 solid2: SolidMaterial,
                 solid3: SolidMaterial,
                 n_wick: int,
                 porosity: float,
                 n_eva: int,
                 n_aba: int,
                 n_con: int,
                 name: str = "Unnamed_Solid",
                 emissivity: float = 0.8,
                 up_view_factor: float = 1.0,
                 down_view_factor: float = 1.0,
                 initial_temp: float = 298.15):
        self.n_wick = n_wick
        self.wall_mat = solid1

        if self.n_wick >= mesh.n_x or self.n_wick <= 0:
            raise ValueError(
                f"Wick radial nodes (n_wick={n_wick}) must be strictly between 0 and total radial nodes (mesh.n_x={mesh.n_x})"
            )

        r_vapor = mesh.inner_radius
        r_in_wall = mesh.x_faces[self.n_wick]
        self.wick_mat = WickMaterial(
            name="HP_Wick_Composite",
            solid_mat=solid3,
            fluid_mat=solid2,
            porosity=porosity,
            r_vapor=r_vapor,
            r_in_wall=r_in_wall
        )

        self.n_eva = n_eva
        self.n_aba = n_aba
        self.n_con = n_con
        self.emissivity = emissivity
        self.up_view_factor = up_view_factor
        self.down_view_factor = down_view_factor

        if self.n_eva + self.n_aba + self.n_con != mesh.n_y:
            raise ValueError(
                f"Axial sections sum ({n_eva}+{n_aba}+{n_con}) must equal total mesh n_y ({mesh.n_y})"
            )

        self.shape_nodes = mesh.shape_nodes
        self.name = name
        self._idx_eva = self.n_eva
        self._idx_aba = self.n_eva + self.n_aba
        self._idx_con = mesh.n_y
        self._slice_eva = slice(0, self._idx_eva)
        self._slice_aba = slice(self._idx_eva, self._idx_aba)
        self._slice_con = slice(self._idx_aba, self._idx_con)

        self._property_update_tol = 5.0e-2
        self._property_cache_initialized = False
        self._k_2d_view = None
        self._k_axial_2d_view = None
        self._k_radial_2d_view = None
        self._rho_2d_view = None
        self._cp_2d_view = None
        self._cap_2d_view = None
        self._volumes_2d_view = None
        self._wick_temperature_cache = np.full((self.n_wick, mesh.n_y), np.nan, dtype=float)
        self._wall_temperature_cache = np.full((mesh.n_x - self.n_wick, mesh.n_y), np.nan, dtype=float)

        self._R_int_left_buffer = np.zeros(mesh.n_y, dtype=float)
        self._R_int_out_buffer = np.zeros(mesh.n_y, dtype=float)
        self._R_int_bottom_buffer = np.zeros(mesh.n_x, dtype=float)
        self._R_int_top_buffer = np.zeros(mesh.n_x, dtype=float)
        self._boundary_work_y = np.zeros(mesh.n_y, dtype=float)
        self._boundary_work_x = np.zeros(mesh.n_x, dtype=float)
        self._face_resistance_work_x = np.empty((max(mesh.n_x - 1, 0), mesh.n_y), dtype=float)
        self._face_resistance_work_y = np.empty((mesh.n_x, max(mesh.n_y - 1, 0)), dtype=float)

        self.enable_frozen_property_correction = False
        self.max_outer_property_corrections = 2
        self.outer_property_tol = 1.0e-3
        self._properties_frozen = False
        self.minimum_physical_temperature = 0.0
        self.maximum_physical_temperature = np.inf
        self.use_anisotropic_wick_conductivity = False
        self.face_conductance_mode = "legacy_harmonic"

        super().__init__(mesh, self.wall_mat, name=name, initial_temp=initial_temp)

        self._setup_virtual_boundaries()
        self._update_boundaries_state(current_time=0.0)

    def _setup_virtual_boundaries(self):
        outer_surface_areas = self.mesh.area_x_matrix[-1, :].copy()
        self.boundaries['outer_eva'] = BoundaryRegion(
            shape=(self.n_eva,),
            area_array=outer_surface_areas[self._slice_eva]
        )
        self.boundaries['outer_aba'] = BoundaryRegion(
            shape=(self.n_aba,),
            area_array=outer_surface_areas[self._slice_aba]
        )
        self.boundaries['outer_con'] = BoundaryRegion(
            shape=(self.n_con,),
            area_array=outer_surface_areas[self._slice_con]
        )

        if 'right' in self.boundaries:
            del self.boundaries['right']

    def _ensure_property_views(self):
        if self._k_2d_view is not None:
            return

        self._k_2d_view = self.k_node.reshape(self.shape_nodes)
        self._k_axial_2d_view = self._k_2d_view
        self._k_radial_2d_view = np.empty(self.shape_nodes, dtype=float)
        self._rho_2d_view = self.rho_node.reshape(self.shape_nodes)
        self._cp_2d_view = self.cp_node.reshape(self.shape_nodes)
        self._cap_2d_view = self.thermal_capacitance.reshape(self.shape_nodes)
        self._volumes_2d_view = self.mesh.geom_data.volumes.reshape(self.shape_nodes)

    @staticmethod
    def _fill_boundary_resistance(target: np.ndarray,
                                  distance: np.ndarray,
                                  conductivity: np.ndarray,
                                  area: np.ndarray,
                                  work: np.ndarray):
        np.multiply(conductivity, area, out=work)
        np.maximum(work, 1.0e-30, out=work)
        np.divide(distance, work, out=target)

    def _refresh_current_state(self, current_time: float):
        self._properties_frozen = False
        self._update_properties()
        self._compute_internal_resistance()
        self._update_boundaries_state(current_time=current_time)
        self._compute_fluxes(current_time)

    def set_wick_conductivity_mode(self, anisotropic: bool):
        self.use_anisotropic_wick_conductivity = bool(anisotropic)
        self._property_cache_initialized = False

    def set_face_conductance_mode(self, mode: str):
        valid_modes = {
            "legacy_harmonic",
            "resistance_split_axial",
            "resistance_split_xy",
            "resistance_split_full",
        }
        if mode not in valid_modes:
            raise ValueError(f"Unsupported face conductance mode: {mode}")
        self.face_conductance_mode = mode

    def _reset_boundary_condition_diagnostics(self):
        for boundary in self.boundaries.values():
            for condition in getattr(boundary, "conditions", []):
                if hasattr(condition, "reset_diagnostics"):
                    condition.reset_diagnostics()

    def _collect_boundary_condition_diagnostics(self):
        diagnostics = []
        for location, boundary in self.boundaries.items():
            for condition in getattr(boundary, "conditions", []):
                if getattr(condition, "invalid_temperature_detected", False):
                    diagnostics.append(
                        f"{location}_trial_min={float(condition.min_trial_temperature):.6f}K"
                    )
        return diagnostics

    def _is_temperature_state_physical(self, temperature_field: np.ndarray) -> bool:
        return (
            np.all(np.isfinite(temperature_field))
            and float(np.min(temperature_field)) >= self.minimum_physical_temperature
            and float(np.max(temperature_field)) <= self.maximum_physical_temperature
        )

    def _compose_failure_message(self, base_message: str) -> str:
        """组装失败信息，包含温度状态及边界条件的诊断信息"""
        details = []
        # 检查最小试算温度是否非有限或超出物理下限
        if not np.isfinite(self.last_trial_temperature_min):
            details.append("trial_temperature_nonfinite")
        elif self.last_trial_temperature_min < self.minimum_physical_temperature:
            details.append(f"trial_temperature_min={self.last_trial_temperature_min:.6f}K")

        # 检查最大试算温度是否非有限或超出物理上限
        if not np.isfinite(self.last_trial_temperature_max):
            details.append("trial_temperature_max_nonfinite")
        elif self.last_trial_temperature_max > self.maximum_physical_temperature:
            details.append(f"trial_temperature_max={self.last_trial_temperature_max:.6f}K")

        # 收集边界条件诊断信息
        details.extend(self._collect_boundary_condition_diagnostics())
        if not details:
            return base_message
        return f"{base_message}; " + "; ".join(details)

    def _solve_ivp_step(self, dt: float, method: str = 'BDF', **kwargs):
        t_span = (self.current_time, self.current_time + dt)
        solve_kwargs = dict(kwargs)

        if hasattr(self, 'get_jac_sparsity') and method in ['BDF', 'Radau']:
            solve_kwargs['jac_sparsity'] = self.get_jac_sparsity()

        sol = solve_ivp(
            fun=self.get_derivatives,
            t_span=t_span,
            y0=self.T.copy(),
            method=method,
            **solve_kwargs
        )

        if not sol.success:
            failure_message = self._compose_failure_message(sol.message)
            self._mark_step_failure(failure_message, failure_time=t_span[1])
            print(f"HeatConduction step failed at t={self.current_time}: {failure_message}")
            return False, None

        candidate_T = sol.y[:, -1]
        if not self._is_temperature_state_physical(candidate_T):
            failure_message = self._compose_failure_message(
                "Unphysical temperature state returned by solver"
            )
            self._mark_step_failure(failure_message, failure_time=t_span[1])
            print(f"HeatConduction step failed at t={self.current_time}: {failure_message}")
            return False, None

        return True, candidate_T

    def _update_properties(self):
        if self._properties_frozen and self._property_cache_initialized:
            return

        self._ensure_property_views()

        T_2d = self.T.reshape(self.shape_nodes)
        k_axial_2d = self._k_axial_2d_view
        k_radial_2d = self._k_radial_2d_view
        rho_2d = self._rho_2d_view
        cp_2d = self._cp_2d_view

        T_wick = T_2d[:self.n_wick, :]
        k_wick_axial = k_axial_2d[:self.n_wick, :]
        k_wick_radial = k_radial_2d[:self.n_wick, :]
        rho_wick = rho_2d[:self.n_wick, :]
        cp_wick = cp_2d[:self.n_wick, :]

        T_wall = T_2d[self.n_wick:, :]
        k_wall_axial = k_axial_2d[self.n_wick:, :]
        k_wall_radial = k_radial_2d[self.n_wick:, :]
        rho_wall = rho_2d[self.n_wick:, :]
        cp_wall = cp_2d[self.n_wick:, :]

        properties_changed = not self._property_cache_initialized

        if not self._property_cache_initialized:
            if self.use_anisotropic_wick_conductivity:
                k_wick_axial[:] = self.wick_mat.conductivity_axial(T_wick)
                k_wick_radial[:] = self.wick_mat.conductivity_radial(T_wick)
            else:
                k_wick_axial[:] = self.wick_mat.conductivity(T_wick)
                k_wick_radial[:] = k_wick_axial
            rho_wick[:] = self.wick_mat.density(T_wick)
            cp_wick[:] = self.wick_mat.heat_capacity(T_wick)
            self._wick_temperature_cache[:] = T_wick
        else:
            wick_update_mask = np.abs(T_wick - self._wick_temperature_cache) > self._property_update_tol
            wick_update_mask |= self.wick_mat.is_high_nonlinearity_temperature(T_wick)
            wick_update_mask |= self.wick_mat.is_high_nonlinearity_temperature(self._wick_temperature_cache)

            if np.any(wick_update_mask):
                T_wick_update = T_wick[wick_update_mask]
                if self.use_anisotropic_wick_conductivity:
                    k_wick_axial[wick_update_mask] = self.wick_mat.conductivity_axial(T_wick_update)
                    k_wick_radial[wick_update_mask] = self.wick_mat.conductivity_radial(T_wick_update)
                else:
                    k_wick_axial[wick_update_mask] = self.wick_mat.conductivity(T_wick_update)
                    k_wick_radial[wick_update_mask] = k_wick_axial[wick_update_mask]
                rho_wick[wick_update_mask] = self.wick_mat.density(T_wick_update)
                cp_wick[wick_update_mask] = self.wick_mat.heat_capacity(T_wick_update)
                self._wick_temperature_cache[wick_update_mask] = T_wick_update
                properties_changed = True

        if T_wall.size > 0:
            if not self._property_cache_initialized:
                k_wall_axial[:] = self.wall_mat.conductivity(T_wall)
                k_wall_radial[:] = k_wall_axial
                rho_wall[:] = self.wall_mat.density(T_wall)
                cp_wall[:] = self.wall_mat.heat_capacity(T_wall)
                self._wall_temperature_cache[:] = T_wall
            else:
                wall_update_mask = np.abs(T_wall - self._wall_temperature_cache) > self._property_update_tol
                if np.any(wall_update_mask):
                    T_wall_update = T_wall[wall_update_mask]
                    k_wall_axial[wall_update_mask] = self.wall_mat.conductivity(T_wall_update)
                    k_wall_radial[wall_update_mask] = k_wall_axial[wall_update_mask]
                    rho_wall[wall_update_mask] = self.wall_mat.density(T_wall_update)
                    cp_wall[wall_update_mask] = self.wall_mat.heat_capacity(T_wall_update)
                    self._wall_temperature_cache[wall_update_mask] = T_wall_update
                    properties_changed = True

        if properties_changed:
            np.multiply(rho_2d, cp_2d, out=self._cap_2d_view)
            self._cap_2d_view *= self._volumes_2d_view

        self._property_cache_initialized = True

    def _compute_radial_conductance_legacy(self):
        k_radial = self._k_radial_2d_view
        k_left = k_radial[:-1, :]
        k_right = k_radial[1:, :]
        np.add(k_left, k_right, out=self._conductance_work_x)
        self._conductance_work_x += 1.0e-30
        np.multiply(k_left, k_right, out=self.G_x_inner)
        self.G_x_inner *= 2.0
        np.divide(self.G_x_inner, self._conductance_work_x, out=self.G_x_inner)

        dx_inner = self.mesh.dx_matrix[1:-1, :]
        area_x_inner = self.mesh.area_x_matrix[1:-1, :]
        np.maximum(dx_inner, 1.0e-30, out=self._conductance_work_x)
        np.multiply(self.G_x_inner, area_x_inner, out=self.G_x_inner)
        np.divide(self.G_x_inner, self._conductance_work_x, out=self.G_x_inner)

    def _compute_axial_conductance_legacy(self):
        k_axial = self._k_axial_2d_view
        k_bottom = k_axial[:, :-1]
        k_top = k_axial[:, 1:]
        np.add(k_bottom, k_top, out=self._conductance_work_y)
        self._conductance_work_y += 1.0e-30
        np.multiply(k_bottom, k_top, out=self.G_y_inner)
        self.G_y_inner *= 2.0
        np.divide(self.G_y_inner, self._conductance_work_y, out=self.G_y_inner)

        dy_inner = self.mesh.dy_matrix[:, 1:-1]
        area_y_inner = self.mesh.area_y_matrix[:, 1:-1]
        np.maximum(dy_inner, 1.0e-30, out=self._conductance_work_y)
        np.multiply(self.G_y_inner, area_y_inner, out=self.G_y_inner)
        np.divide(self.G_y_inner, self._conductance_work_y, out=self.G_y_inner)

    def _compute_axial_conductance_resistance_split(self):
        k_axial = self._k_axial_2d_view
        k_bottom = k_axial[:, :-1]
        k_top = k_axial[:, 1:]
        area_y_inner = self.mesh.area_y_matrix[:, 1:-1]

        dy_bottom_half = (self.mesh.y_faces[1:-1] - self.mesh.y_centers[:-1])[np.newaxis, :]
        dy_top_half = (self.mesh.y_centers[1:] - self.mesh.y_faces[1:-1])[np.newaxis, :]

        np.multiply(k_bottom, area_y_inner, out=self._conductance_work_y)
        np.maximum(self._conductance_work_y, 1.0e-30, out=self._conductance_work_y)
        np.divide(dy_bottom_half, self._conductance_work_y, out=self._face_resistance_work_y)

        np.multiply(k_top, area_y_inner, out=self._conductance_work_y)
        np.maximum(self._conductance_work_y, 1.0e-30, out=self._conductance_work_y)
        self._face_resistance_work_y += dy_top_half / self._conductance_work_y

        np.maximum(self._face_resistance_work_y, 1.0e-30, out=self._face_resistance_work_y)
        np.divide(1.0, self._face_resistance_work_y, out=self.G_y_inner)

    def _compute_radial_conductance_resistance_split(self):
        k_radial = self._k_radial_2d_view
        k_left = k_radial[:-1, :]
        k_right = k_radial[1:, :]

        if self.mesh.geometry_type == "cylindrical":
            radial_length = self.mesh.dy[np.newaxis, :]
            r_left = self.mesh.x_centers[:-1][:, np.newaxis]
            r_face = self.mesh.x_faces[1:-1][:, np.newaxis]
            r_right = self.mesh.x_centers[1:][:, np.newaxis]

            log_left = np.log(np.maximum(r_face, 1.0e-30) / np.maximum(r_left, 1.0e-30))
            log_right = np.log(np.maximum(r_right, 1.0e-30) / np.maximum(r_face, 1.0e-30))

            np.multiply(k_left, radial_length, out=self._conductance_work_x)
            self._conductance_work_x *= 2.0 * np.pi
            np.maximum(self._conductance_work_x, 1.0e-30, out=self._conductance_work_x)
            np.divide(log_left, self._conductance_work_x, out=self._face_resistance_work_x)

            np.multiply(k_right, radial_length, out=self._conductance_work_x)
            self._conductance_work_x *= 2.0 * np.pi
            np.maximum(self._conductance_work_x, 1.0e-30, out=self._conductance_work_x)
            self._face_resistance_work_x += log_right / self._conductance_work_x
        else:
            area_x_inner = self.mesh.area_x_matrix[1:-1, :]
            dx_left_half = (self.mesh.x_faces[1:-1] - self.mesh.x_centers[:-1])[:, np.newaxis]
            dx_right_half = (self.mesh.x_centers[1:] - self.mesh.x_faces[1:-1])[:, np.newaxis]

            np.multiply(k_left, area_x_inner, out=self._conductance_work_x)
            np.maximum(self._conductance_work_x, 1.0e-30, out=self._conductance_work_x)
            np.divide(dx_left_half, self._conductance_work_x, out=self._face_resistance_work_x)

            np.multiply(k_right, area_x_inner, out=self._conductance_work_x)
            np.maximum(self._conductance_work_x, 1.0e-30, out=self._conductance_work_x)
            self._face_resistance_work_x += dx_right_half / self._conductance_work_x

        np.maximum(self._face_resistance_work_x, 1.0e-30, out=self._face_resistance_work_x)
        np.divide(1.0, self._face_resistance_work_x, out=self.G_x_inner)

    def _compute_internal_resistance(self):
        self._ensure_property_views()

        mode = self.face_conductance_mode
        if mode == "legacy_harmonic":
            self._compute_radial_conductance_legacy()
            self._compute_axial_conductance_legacy()
            return

        if mode == "resistance_split_axial":
            self._compute_radial_conductance_legacy()
            self._compute_axial_conductance_resistance_split()
            return

        if mode in {"resistance_split_xy", "resistance_split_full"}:
            self._compute_radial_conductance_resistance_split()
            self._compute_axial_conductance_resistance_split()
            return

        raise ValueError(f"Unsupported face conductance mode: {mode}")

    def _fill_boundary_halfcell_resistance_radial(self,
                                                  target: np.ndarray,
                                                  conductivity: np.ndarray,
                                                  radius_center: float,
                                                  radius_face: float,
                                                  axial_lengths: np.ndarray):
        np.multiply(conductivity, axial_lengths, out=self._boundary_work_y)
        self._boundary_work_y *= 2.0 * np.pi
        np.maximum(self._boundary_work_y, 1.0e-30, out=self._boundary_work_y)
        log_ratio = np.log(np.maximum(radius_face, radius_center) / np.maximum(min(radius_face, radius_center), 1.0e-30))
        target[:] = log_ratio / self._boundary_work_y

    def _fill_boundary_halfcell_resistance_axial(self,
                                                 target: np.ndarray,
                                                 distances: float,
                                                 conductivity: np.ndarray,
                                                 area: np.ndarray):
        np.multiply(conductivity, area, out=self._boundary_work_x)
        np.maximum(self._boundary_work_x, 1.0e-30, out=self._boundary_work_x)
        target[:] = distances / self._boundary_work_x

    def _update_boundaries_state(self, current_time: float = None):
        self._ensure_property_views()

        T_2d = self.T.reshape(self.shape_nodes)
        k_axial_2d = self._k_axial_2d_view
        k_radial_2d = self._k_radial_2d_view
        dx_mat = self.mesh.dx_matrix
        dy_mat = self.mesh.dy_matrix
        ax_mat = self.mesh.area_x_matrix
        ay_mat = self.mesh.area_y_matrix

        if not self._properties_frozen:
            if self.face_conductance_mode == "resistance_split_full":
                if self.mesh.geometry_type == "cylindrical":
                    axial_lengths = self.mesh.dy
                    self._fill_boundary_halfcell_resistance_radial(
                        self._R_int_left_buffer,
                        k_radial_2d[0, :],
                        float(self.mesh.x_centers[0]),
                        float(self.mesh.x_faces[0]),
                        axial_lengths
                    )
                    self._fill_boundary_halfcell_resistance_radial(
                        self._R_int_out_buffer,
                        k_radial_2d[-1, :],
                        float(self.mesh.x_centers[-1]),
                        float(self.mesh.x_faces[-1]),
                        axial_lengths
                    )
                else:
                    self._fill_boundary_resistance(
                        self._R_int_left_buffer,
                        dx_mat[0, :],
                        k_radial_2d[0, :],
                        ax_mat[0, :],
                        self._boundary_work_y
                    )
                    self._fill_boundary_resistance(
                        self._R_int_out_buffer,
                        dx_mat[-1, :],
                        k_radial_2d[-1, :],
                        ax_mat[-1, :],
                        self._boundary_work_y
                    )

                self._fill_boundary_halfcell_resistance_axial(
                    self._R_int_bottom_buffer,
                    float(self.mesh.y_centers[0] - self.mesh.y_faces[0]),
                    k_axial_2d[:, 0],
                    ay_mat[:, 0]
                )
                self._fill_boundary_halfcell_resistance_axial(
                    self._R_int_top_buffer,
                    float(self.mesh.y_faces[-1] - self.mesh.y_centers[-1]),
                    k_axial_2d[:, -1],
                    ay_mat[:, -1]
                )
            else:
                self._fill_boundary_resistance(
                    self._R_int_left_buffer,
                    dx_mat[0, :],
                    k_radial_2d[0, :],
                    ax_mat[0, :],
                    self._boundary_work_y
                )
                self._fill_boundary_resistance(
                    self._R_int_out_buffer,
                    dx_mat[-1, :],
                    k_radial_2d[-1, :],
                    ax_mat[-1, :],
                    self._boundary_work_y
                )
                self._fill_boundary_resistance(
                    self._R_int_bottom_buffer,
                    dy_mat[:, 0],
                    k_axial_2d[:, 0],
                    ay_mat[:, 0],
                    self._boundary_work_x
                )
                self._fill_boundary_resistance(
                    self._R_int_top_buffer,
                    dy_mat[:, -1],
                    k_axial_2d[:, -1],
                    ay_mat[:, -1],
                    self._boundary_work_x
                )

        if 'left' in self.boundaries:
            self.boundaries['left'].update_internal_state(
                T_2d[0, :],
                self._R_int_left_buffer,
                current_time=current_time
            )

        if 'bottom' in self.boundaries:
            self.boundaries['bottom'].update_internal_state(
                T_2d[:, 0],
                self._R_int_bottom_buffer,
                current_time=current_time
            )

        if 'top' in self.boundaries:
            self.boundaries['top'].update_internal_state(
                T_2d[:, -1],
                self._R_int_top_buffer,
                current_time=current_time
            )

        T_out = T_2d[-1, :]

        if 'outer_eva' in self.boundaries:
            self.boundaries['outer_eva'].update_internal_state(
                T_out[self._slice_eva],
                self._R_int_out_buffer[self._slice_eva],
                current_time=current_time
            )

        if 'outer_aba' in self.boundaries:
            self.boundaries['outer_aba'].update_internal_state(
                T_out[self._slice_aba],
                self._R_int_out_buffer[self._slice_aba],
                current_time=current_time
            )

        if 'outer_con' in self.boundaries:
            self.boundaries['outer_con'].update_internal_state(
                T_out[self._slice_con],
                self._R_int_out_buffer[self._slice_con],
                current_time=current_time
            )

    def _compute_fluxes(self, t: float) -> np.ndarray:
        Q_net_2d = self.Q_net_2d_buffer
        Q_net_2d.fill(0.0)
        T_2d = self.T.reshape(self.shape_nodes)

        np.subtract(T_2d[:-1, :], T_2d[1:, :], out=self._flux_x_buffer)
        np.multiply(self._flux_x_buffer, self.G_x_inner, out=self._flux_x_buffer)
        Q_net_2d[:-1, :] -= self._flux_x_buffer
        Q_net_2d[1:, :] += self._flux_x_buffer

        np.subtract(T_2d[:, :-1], T_2d[:, 1:], out=self._flux_y_buffer)
        np.multiply(self._flux_y_buffer, self.G_y_inner, out=self._flux_y_buffer)
        Q_net_2d[:, :-1] -= self._flux_y_buffer
        Q_net_2d[:, 1:] += self._flux_y_buffer

        if 'left' in self.boundaries:
            Q_net_2d[0, :] += self.boundaries['left'].compute_net_flux_for_solver()

        if 'bottom' in self.boundaries:
            Q_net_2d[:, 0] += self.boundaries['bottom'].compute_net_flux_for_solver()

        if 'top' in self.boundaries:
            Q_net_2d[:, -1] += self.boundaries['top'].compute_net_flux_for_solver()

        if 'outer_eva' in self.boundaries:
            Q_net_2d[-1, self._slice_eva] += self.boundaries['outer_eva'].compute_net_flux_for_solver()

        if 'outer_aba' in self.boundaries:
            Q_net_2d[-1, self._slice_aba] += self.boundaries['outer_aba'].compute_net_flux_for_solver()

        if 'outer_con' in self.boundaries:
            Q_net_2d[-1, self._slice_con] += self.boundaries['outer_con'].compute_net_flux_for_solver()

        return Q_net_2d.reshape(-1).copy()

    def get_derivatives(self, t: float, T_current: np.ndarray) -> np.ndarray:
        if not self._properties_frozen:
            return super().get_derivatives(t, T_current)

        self.T[:] = T_current
        self._record_trial_state(t, T_current)
        self._update_boundaries_state(current_time=t)
        Q_net_conduction = self._compute_fluxes(t)
        self._update_sources(t)

        total_Q = Q_net_conduction + self.Q_source
        self.dTdt.fill(0.0)
        np.divide(
            total_Q,
            self.thermal_capacitance,
            out=self.dTdt,
            where=self.thermal_capacitance > 1.0e-30
        )
        return self.dTdt

    def step(self, dt: float, method: str = 'BDF', **kwargs) -> bool:
        attempted_end_time = self.current_time + dt
        self._reset_step_diagnostics(attempted_end_time)
        self._reset_boundary_condition_diagnostics()
        self.save_state()

        if not self.enable_frozen_property_correction:
            success = super().step(dt, method=method, **kwargs)
            if not success:
                failure_message = self.last_step_failure_message or self._compose_failure_message(
                    "Solver step failed"
                )
                self.load_state()
                self._mark_step_failure(failure_message, failure_time=attempted_end_time)
                return False

            if not self._is_temperature_state_physical(self.T):
                failure_message = self._compose_failure_message(
                    "Accepted state violates physical temperature bound"
                )
                self.load_state()
                self._mark_step_failure(failure_message, failure_time=attempted_end_time)
                return False

            self._mark_step_success()
            return True

        time_start = self.current_time
        T_start = self.T.copy()
        guess_T = T_start.copy()
        accepted_T = None
        max_outer_iters = max(int(self.max_outer_property_corrections), 1)

        for _ in range(max_outer_iters):
            self.T[:] = guess_T
            self.current_time = time_start
            self._properties_frozen = False
            self._update_properties()
            self._compute_internal_resistance()
            self._update_boundaries_state(current_time=time_start)
            self._compute_fluxes(time_start)

            self.T[:] = T_start
            self.current_time = time_start
            self._properties_frozen = True

            success, candidate_T = self._solve_ivp_step(dt, method=method, **kwargs)
            if not success:
                failure_message = self.last_step_failure_message or self._compose_failure_message(
                    "Frozen-property correction step failed"
                )
                self.load_state()
                self._mark_step_failure(failure_message, failure_time=attempted_end_time)
                return False

            accepted_T = candidate_T
            delta = float(np.max(np.abs(candidate_T - guess_T)))
            scale = max(1.0, float(np.max(np.abs(candidate_T))))
            if delta <= self.outer_property_tol * scale:
                break

            guess_T = candidate_T.copy()

        if accepted_T is None or not self._is_temperature_state_physical(accepted_T):
            failure_message = self._compose_failure_message(
                "Frozen-property correction produced unphysical state"
            )
            self.load_state()
            self._mark_step_failure(failure_message, failure_time=attempted_end_time)
            return False

        self.T[:] = accepted_T
        self.current_time = time_start + dt
        self._refresh_current_state(current_time=self.current_time)
        self._mark_step_success()
        return True

    def get_boundary_node_capacitance(self, location: str) -> np.ndarray:
        cap_2d = self.thermal_capacitance.reshape(self.shape_nodes)

        if location == 'outer_eva':
            return cap_2d[-1, self._slice_eva].copy()
        if location == 'outer_aba':
            return cap_2d[-1, self._slice_aba].copy()
        if location == 'outer_con':
            return cap_2d[-1, self._slice_con].copy()

        return super().get_boundary_node_capacitance(location)
