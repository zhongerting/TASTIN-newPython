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

        self.enable_frozen_property_correction = False
        self.max_outer_property_corrections = 2
        self.outer_property_tol = 1.0e-3
        self._properties_frozen = False

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
            print(f"HeatConduction step failed at t={self.current_time}: {sol.message}")
            return False, None

        return True, sol.y[:, -1]

    def _update_properties(self):
        if self._properties_frozen and self._property_cache_initialized:
            return

        self._ensure_property_views()

        T_2d = self.T.reshape(self.shape_nodes)
        k_2d = self._k_2d_view
        rho_2d = self._rho_2d_view
        cp_2d = self._cp_2d_view

        T_wick = T_2d[:self.n_wick, :]
        k_wick = k_2d[:self.n_wick, :]
        rho_wick = rho_2d[:self.n_wick, :]
        cp_wick = cp_2d[:self.n_wick, :]

        T_wall = T_2d[self.n_wick:, :]
        k_wall = k_2d[self.n_wick:, :]
        rho_wall = rho_2d[self.n_wick:, :]
        cp_wall = cp_2d[self.n_wick:, :]

        properties_changed = not self._property_cache_initialized

        if not self._property_cache_initialized:
            k_wick[:] = self.wick_mat.conductivity(T_wick)
            rho_wick[:] = self.wick_mat.density(T_wick)
            cp_wick[:] = self.wick_mat.heat_capacity(T_wick)
            self._wick_temperature_cache[:] = T_wick
        else:
            wick_update_mask = np.abs(T_wick - self._wick_temperature_cache) > self._property_update_tol
            wick_update_mask |= self.wick_mat.is_high_nonlinearity_temperature(T_wick)
            wick_update_mask |= self.wick_mat.is_high_nonlinearity_temperature(self._wick_temperature_cache)

            if np.any(wick_update_mask):
                T_wick_update = T_wick[wick_update_mask]
                k_wick[wick_update_mask] = self.wick_mat.conductivity(T_wick_update)
                rho_wick[wick_update_mask] = self.wick_mat.density(T_wick_update)
                cp_wick[wick_update_mask] = self.wick_mat.heat_capacity(T_wick_update)
                self._wick_temperature_cache[wick_update_mask] = T_wick_update
                properties_changed = True

        if T_wall.size > 0:
            if not self._property_cache_initialized:
                k_wall[:] = self.wall_mat.conductivity(T_wall)
                rho_wall[:] = self.wall_mat.density(T_wall)
                cp_wall[:] = self.wall_mat.heat_capacity(T_wall)
                self._wall_temperature_cache[:] = T_wall
            else:
                wall_update_mask = np.abs(T_wall - self._wall_temperature_cache) > self._property_update_tol
                if np.any(wall_update_mask):
                    T_wall_update = T_wall[wall_update_mask]
                    k_wall[wall_update_mask] = self.wall_mat.conductivity(T_wall_update)
                    rho_wall[wall_update_mask] = self.wall_mat.density(T_wall_update)
                    cp_wall[wall_update_mask] = self.wall_mat.heat_capacity(T_wall_update)
                    self._wall_temperature_cache[wall_update_mask] = T_wall_update
                    properties_changed = True

        if properties_changed:
            np.multiply(rho_2d, cp_2d, out=self._cap_2d_view)
            self._cap_2d_view *= self._volumes_2d_view

        self._property_cache_initialized = True

    def _update_boundaries_state(self, current_time: float = None):
        self._ensure_property_views()

        T_2d = self.T.reshape(self.shape_nodes)
        k_2d = self._k_2d_view
        dx_mat = self.mesh.dx_matrix
        dy_mat = self.mesh.dy_matrix
        ax_mat = self.mesh.area_x_matrix
        ay_mat = self.mesh.area_y_matrix

        if not self._properties_frozen:
            self._fill_boundary_resistance(
                self._R_int_left_buffer,
                dx_mat[0, :],
                k_2d[0, :],
                ax_mat[0, :],
                self._boundary_work_y
            )
            self._fill_boundary_resistance(
                self._R_int_out_buffer,
                dx_mat[-1, :],
                k_2d[-1, :],
                ax_mat[-1, :],
                self._boundary_work_y
            )
            self._fill_boundary_resistance(
                self._R_int_bottom_buffer,
                dy_mat[:, 0],
                k_2d[:, 0],
                ay_mat[:, 0],
                self._boundary_work_x
            )
            self._fill_boundary_resistance(
                self._R_int_top_buffer,
                dy_mat[:, -1],
                k_2d[:, -1],
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
        if not self.enable_frozen_property_correction:
            return super().step(dt, method=method, **kwargs)

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
                self.T[:] = T_start
                self.current_time = time_start
                self._refresh_current_state(current_time=time_start)
                return False

            accepted_T = candidate_T
            delta = float(np.max(np.abs(candidate_T - guess_T)))
            scale = max(1.0, float(np.max(np.abs(candidate_T))))
            if delta <= self.outer_property_tol * scale:
                break

            guess_T = candidate_T.copy()

        self.T[:] = accepted_T
        self.current_time = time_start + dt
        self._refresh_current_state(current_time=self.current_time)
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
