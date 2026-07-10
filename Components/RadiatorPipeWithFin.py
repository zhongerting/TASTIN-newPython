import numpy as np

from Components.BaseComponent import BaseComponent
from Solvers.Couplers import FluidSolidCouple
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D


class RadiatorPipeWithFin(BaseComponent):
    """
    NaK radiator pipe with a reduced-order fin radiation branch.

    Differences from HPwithFin:
    - the pipe wall is a normal cylindrical HeatConduction2D solid, not HeatPipe2D;
    - the inner surface couples to an external fluid channel through FluidSolidCouple;
    - the bare outer tube radiates directly to space;
    - the copper bands are solved as quasi-steady 1D fin strips per axial slice.
    """

    sigma = 5.670374419e-8

    def __init__(
        self,
        name: str,
        fluid_channel,
        wall_material,
        tube_inner_diameter: float,
        tube_outer_diameter: float,
        tube_length: float,
        n_axial: int,
        n_radial_wall: int,
        fin_thickness: float,
        fin_width_upper: float,
        fin_width_lower: float,
        n_fin_width: int,
        correlation_func,
        tube_emissivity: float = 0.8,
        fin_emissivity: float = 0.8,
        tube_area_scale: float = 1.0,
        fin_area_scale: float = 1.0,
        T_space: float = 3.0,
        initial_temp: float = 727.0,
        fin_conductivity: float = 348.9,
        fin_view_factor: float = 1.0,
        contact_resistance_m2k_w: float = 0.0,
        coupling_time_scheme: str = "current",
        solid_ode_method: str = "BDF",
    ):
        super().__init__(name)

        self.fluid_channel = fluid_channel
        self.tube_inner_diameter = float(tube_inner_diameter)
        self.tube_outer_diameter = float(tube_outer_diameter)
        self.tube_length = float(tube_length)
        self.n_axial = int(n_axial)
        self.n_radial_wall = int(n_radial_wall)
        self.fin_thickness = float(fin_thickness)
        self.fin_width_upper = float(fin_width_upper)
        self.fin_width_lower = float(fin_width_lower)
        self.n_fin_width = int(n_fin_width)
        self.tube_emissivity = float(tube_emissivity)
        self.fin_emissivity = float(fin_emissivity)
        self.tube_area_scale = float(tube_area_scale)
        self.fin_area_scale = float(fin_area_scale)
        self.T_space = float(T_space)
        self._default_radiation_background_k = float(T_space)
        self.radiation_background_temperature = np.full(self.n_axial, float(T_space), dtype=float)
        self.fin_conductivity = float(fin_conductivity)
        self.fin_view_factor = float(fin_view_factor)
        self.contact_resistance_m2k_w = float(contact_resistance_m2k_w)

        if self.n_axial <= 0:
            raise ValueError("n_axial must be positive")
        if self.n_fin_width <= 0:
            raise ValueError("n_fin_width must be positive")
        if self.fin_thickness <= 0.0:
            raise ValueError("fin_thickness must be positive")
        if self.tube_outer_diameter <= self.tube_inner_diameter:
            raise ValueError("tube_outer_diameter must be larger than tube_inner_diameter")

        tube_wall_thickness = 0.5 * (self.tube_outer_diameter - self.tube_inner_diameter)
        self.mesh = Mesh2D(
            x_dim=tube_wall_thickness,
            n_x=max(1, self.n_radial_wall),
            y_dim=self.tube_length,
            n_y=self.n_axial,
            geometry_type="cylindrical",
            inner_radius=0.5 * self.tube_inner_diameter,
        )
        self.wall = HeatConduction2D(
            mesh=self.mesh,
            material=wall_material,
            name=f"{self.name}_Wall",
            initial_temp=initial_temp,
        )
        self.wall.set_ode_method(solid_ode_method)

        self.node_length = self.tube_length / self.n_axial
        self.tube_bare_area = (
            np.array(self.wall.boundaries["right"].area, dtype=float)
            * self.tube_area_scale
        )

        self.fin_pitch_width = np.linspace(
            self.fin_width_upper,
            self.fin_width_lower,
            self.n_axial,
        )
        self.fin_net_width = np.maximum(self.fin_pitch_width - self.tube_outer_diameter, 0.0)
        self.fin_height = 0.5 * self.fin_net_width
        # Two identical half-fins are represented as one doubled-width 1D strip.
        self.fin_strip_width = np.full(self.n_axial, 2.0 * self.node_length)
        self.fin_radiating_area = (
            2.0
            * self.fin_strip_width
            * self.fin_height
            * self.fin_area_scale
            * self.fin_view_factor
        )
        self.fin_root_area = self.fin_strip_width * self.fin_thickness
        with np.errstate(divide="ignore", invalid="ignore"):
            self.contact_resistance = self.contact_resistance_m2k_w / self.fin_root_area
        self.contact_resistance = np.nan_to_num(self.contact_resistance, nan=0.0, posinf=0.0, neginf=0.0)

        self.bc_tube_radiation = self.wall.boundaries["right"].add_dynamic_radiation_condition(
            emissivity=self.tube_emissivity,
            bare_area_array=self.tube_bare_area,
            T_env=self.T_space,
        )
        self.bc_fin = self.wall.boundaries["right"].add_resistance_condition(
            T_ext=self.T_space,
            R_ext=np.full(self.n_axial, 1.0e15, dtype=float),
        )

        self.coupler = FluidSolidCouple(
            name=f"{self.name}_FluidSolid",
            fluid=self.fluid_channel,
            solid_boundary_region=self.wall.boundaries["left"],
            heated_perimeter=np.pi * self.tube_inner_diameter,
            correlation_func=correlation_func,
            solid_node_capacitance=self.wall.get_boundary_node_capacitance("left"),
            coupling_time_scheme=coupling_time_scheme,
        )

        self.last_fin_temperature = np.zeros((self.n_axial, self.n_fin_width), dtype=float)
        self.last_fin_radiation_distribution = np.zeros(self.n_axial, dtype=float)
        self.last_fin_net_from_root_distribution = np.zeros(self.n_axial, dtype=float)
        self.last_fin_conductance_distribution = np.full(self.n_axial, 1.0e-15, dtype=float)
        self.last_fin_effective_temperature_distribution = np.full(self.n_axial, self.T_space, dtype=float)
        self.last_fin_equivalent_resistance_distribution = np.full(self.n_axial, 1.0e15, dtype=float)
        self.last_tube_radiation_distribution = np.zeros(self.n_axial, dtype=float)
        self.last_fin_iteration_count = 0
        self.last_fin_max_delta = 0.0
        self.last_fin_used_warm_start = False
        self._has_valid_fin_temperature = False
        self._fin_scratch_shape = None
        self._fin_scratch = {}

        self.wall.initialize_state()

    def set_radiation_background_temperature(self, value):
        """Update the equivalent radiation background used by tube and fin radiation."""
        if np.isscalar(value):
            background = np.full(self.n_axial, float(value), dtype=float)
        else:
            background = np.asarray(value, dtype=float)
            if background.shape != (self.n_axial,):
                background = np.broadcast_to(background, (self.n_axial,))
                background = np.array(background, dtype=float, copy=True)
        background = np.nan_to_num(
            background,
            nan=self._default_radiation_background_k,
            posinf=1.0e6,
            neginf=self._default_radiation_background_k,
        )
        background = np.maximum(background, 1.0e-3)
        self.radiation_background_temperature[:] = background
        self.T_space = float(np.mean(background))
        self.bc_tube_radiation.update_params(T_env=background)

    def restore_default_radiation_background(self):
        self.set_radiation_background_temperature(self._default_radiation_background_k)

    def get_radiation_surface_temperature(self):
        return np.asarray(self.wall.boundaries["right"].T_surface, dtype=float).copy()

    @staticmethod
    def _solve_tridiagonal_inplace(a, b, c, d, c_prime, d_prime, solution):
        nh = b.shape[1]
        denom_floor = 1.0e-12

        b0 = b[:, 0]
        b0_safe = np.copysign(
            np.maximum(np.abs(b0), denom_floor),
            np.where(b0 == 0.0, 1.0, b0),
        )
        c_prime[:, 0] = c[:, 0] / b0_safe
        d_prime[:, 0] = d[:, 0] / b0_safe

        for j in range(1, nh):
            denom = b[:, j] - a[:, j] * c_prime[:, j - 1]
            denom_safe = np.copysign(
                np.maximum(np.abs(denom), denom_floor),
                np.where(denom == 0.0, 1.0, denom),
            )
            c_prime[:, j] = c[:, j] / denom_safe
            d_prime[:, j] = (d[:, j] - a[:, j] * d_prime[:, j - 1]) / denom_safe

        solution[:, -1] = d_prime[:, -1]
        for j in range(nh - 2, -1, -1):
            solution[:, j] = d_prime[:, j] - c_prime[:, j] * solution[:, j + 1]

    def _get_fin_scratch(self, n_active: int, n_width: int):
        shape = (int(n_active), int(n_width))
        if self._fin_scratch_shape != shape:
            self._fin_scratch_shape = shape
            self._fin_scratch = {
                "a": np.zeros(shape),
                "b": np.zeros(shape),
                "c": np.zeros(shape),
                "d": np.zeros(shape),
                "c_prime": np.zeros(shape),
                "d_prime": np.zeros(shape),
                "T_new": np.zeros(shape),
                "rad_term": np.zeros(shape),
                "s_a": np.zeros(shape),
                "s_b": np.zeros(shape),
                "s_c": np.zeros(shape),
                "s_d": np.zeros(shape),
                "s_c_prime": np.zeros(shape),
                "s_d_prime": np.zeros(shape),
                "sensitivity": np.zeros(shape),
                "rad_deriv": np.zeros(shape),
            }
        return self._fin_scratch

    def _active_fin_mask(self, T_root: np.ndarray) -> np.ndarray:
        return (
            np.isfinite(T_root)
            & (self.fin_height > 0.0)
            & (self.fin_strip_width > 0.0)
            & (self.fin_radiating_area > 0.0)
        )

    def _solve_fin_quasi_steady(self, T_root: np.ndarray):
        T_root = np.asarray(T_root, dtype=float)
        na = len(T_root)
        nw = self.n_fin_width
        root_initial = np.repeat(T_root[:, np.newaxis], nw, axis=1)
        warm_start_valid = (
            self._has_valid_fin_temperature
            and getattr(self, "last_fin_temperature", None) is not None
            and self.last_fin_temperature.shape == (na, nw)
            and np.all(np.isfinite(self.last_fin_temperature))
        )
        if warm_start_valid:
            T = np.array(self.last_fin_temperature, dtype=float, copy=True)
            invalid_rows = ~np.isfinite(T_root)
            if np.any(invalid_rows):
                T[invalid_rows, :] = root_initial[invalid_rows, :]
        else:
            T = root_initial
        Q_fin = np.zeros_like(T_root)
        active = self._active_fin_mask(T_root)
        if not np.any(active):
            self.last_fin_iteration_count = 0
            self.last_fin_max_delta = 0.0
            self.last_fin_used_warm_start = bool(warm_start_valid)
            return T, Q_fin

        T_active = T[active].copy()
        T_root_active = T_root[active]
        height_active = self.fin_height[active]
        width_active = self.fin_strip_width[active]
        area_scale_active = self.fin_area_scale * self.fin_view_factor
        dx = height_active / nw
        Ac = width_active * self.fin_thickness
        P_rad = 2.0 * width_active * area_scale_active
        G = self.fin_conductivity * Ac / np.maximum(dx, 1.0e-30)
        G_base = 2.0 * G

        scratch = self._get_fin_scratch(len(T_root_active), nw)
        a = scratch["a"]
        b = scratch["b"]
        c = scratch["c"]
        d = scratch["d"]
        c_prime = scratch["c_prime"]
        d_prime = scratch["d_prime"]
        T_new = scratch["T_new"]
        rad_term = scratch["rad_term"]

        max_iter = 50
        tol = 1.0e-4
        iteration_count = 0
        max_delta = 0.0
        for iteration in range(max_iter):
            background_active = self.radiation_background_temperature[active]
            rad_term[:, :] = self.fin_emissivity * self.sigma * (
                T_active ** 2 + background_active[:, np.newaxis] ** 2
            ) * (T_active + background_active[:, np.newaxis])
            rad_term *= P_rad[:, np.newaxis] * dx[:, np.newaxis]

            a.fill(0.0)
            b.fill(0.0)
            c.fill(0.0)
            d.fill(0.0)

            b[:, 0] = G_base + G + rad_term[:, 0]
            c[:, 0] = -G
            d[:, 0] = G_base * T_root_active + rad_term[:, 0] * background_active

            if nw > 1:
                a[:, 1:-1] = -G[:, np.newaxis]
                b[:, 1:-1] = 2.0 * G[:, np.newaxis] + rad_term[:, 1:-1]
                c[:, 1:-1] = -G[:, np.newaxis]
                d[:, 1:-1] = rad_term[:, 1:-1] * background_active[:, np.newaxis]

                a[:, -1] = -G
                b[:, -1] = G + rad_term[:, -1]
                d[:, -1] = rad_term[:, -1] * background_active

            self._solve_tridiagonal_inplace(a, b, c, d, c_prime, d_prime, T_new)
            err = float(np.max(np.abs(T_new - T_active)))
            iteration_count = iteration + 1
            max_delta = err
            T_active[:, :] = T_new
            if err < tol:
                break

        T[active] = T_active
        Q_fin[active] = np.sum(
            self.fin_emissivity
            * self.sigma
            * P_rad[:, np.newaxis]
            * dx[:, np.newaxis]
            * (T_active ** 4 - background_active[:, np.newaxis] ** 4),
            axis=1,
        )
        self.last_fin_iteration_count = int(iteration_count)
        self.last_fin_max_delta = float(max_delta)
        self.last_fin_used_warm_start = bool(warm_start_valid)
        self._has_valid_fin_temperature = True
        return T, Q_fin

    def _compute_fin_tangent_conductance(self, T_root: np.ndarray, T_fin: np.ndarray, active_mask: np.ndarray):
        T_root = np.asarray(T_root, dtype=float)
        T_fin = np.asarray(T_fin, dtype=float)
        active = np.asarray(active_mask, dtype=bool).copy()
        active &= np.isfinite(T_root)
        active &= np.all(np.isfinite(T_fin), axis=1)
        conductance = np.full_like(T_root, np.nan, dtype=float)
        if not np.any(active):
            return conductance

        nw = self.n_fin_width
        T_active = T_fin[active]
        height_active = self.fin_height[active]
        width_active = self.fin_strip_width[active]
        area_scale_active = self.fin_area_scale * self.fin_view_factor
        dx = height_active / nw
        Ac = width_active * self.fin_thickness
        P_rad = 2.0 * width_active * area_scale_active
        G = self.fin_conductivity * Ac / np.maximum(dx, 1.0e-30)
        G_base = 2.0 * G

        scratch = self._get_fin_scratch(int(np.sum(active)), nw)
        a = scratch["s_a"]
        b = scratch["s_b"]
        c = scratch["s_c"]
        d = scratch["s_d"]
        c_prime = scratch["s_c_prime"]
        d_prime = scratch["s_d_prime"]
        sensitivity = scratch["sensitivity"]
        rad_deriv = scratch["rad_deriv"]

        a.fill(0.0)
        b.fill(0.0)
        c.fill(0.0)
        d.fill(0.0)
        rad_deriv[:, :] = (
            4.0
            * self.fin_emissivity
            * self.sigma
            * P_rad[:, np.newaxis]
            * dx[:, np.newaxis]
            * (T_active ** 3)
        )

        b[:, 0] = G_base + G + rad_deriv[:, 0]
        c[:, 0] = -G
        d[:, 0] = G_base

        if nw > 1:
            a[:, 1:-1] = -G[:, np.newaxis]
            b[:, 1:-1] = 2.0 * G[:, np.newaxis] + rad_deriv[:, 1:-1]
            c[:, 1:-1] = -G[:, np.newaxis]
            a[:, -1] = -G
            b[:, -1] = G + rad_deriv[:, -1]

        self._solve_tridiagonal_inplace(a, b, c, d, c_prime, d_prime, sensitivity)
        conductance[active] = G_base * (1.0 - sensitivity[:, 0])
        return conductance

    def pre_step(self, dt: float, current_time: float):
        T_root, _ = self.wall.boundaries["right"].get_coupling_surface_snapshot()
        T_root = np.asarray(T_root, dtype=float)
        T_fin, Q_fin = self._solve_fin_quasi_steady(T_root)
        active = self._active_fin_mask(T_root)
        lambda_raw = self._compute_fin_tangent_conductance(T_root, T_fin, active)

        T_space_safe = np.maximum(self.radiation_background_temperature, 1.0e-3)
        T_seg_safe = np.maximum(T_fin, 1.0e-3)
        nw = self.n_fin_width
        dx_arr = np.zeros(self.n_axial, dtype=float)
        mask = self.fin_height > 0.0
        dx_arr[mask] = self.fin_height[mask] / nw
        P_rad = 2.0 * self.fin_strip_width * self.fin_area_scale * self.fin_view_factor
        A_seg = P_rad[:, np.newaxis] * dx_arr[:, np.newaxis]
        lambda_fallback = np.sum(
            self.fin_emissivity
            * self.sigma
            * A_seg
            * (T_seg_safe + T_space_safe[:, np.newaxis])
            * (T_seg_safe ** 2 + T_space_safe[:, np.newaxis] ** 2),
            axis=1,
        )
        lambda_fallback = np.clip(
            np.nan_to_num(lambda_fallback, nan=1.0e-12, posinf=1.0e12, neginf=1.0e-12),
            1.0e-12,
            1.0e12,
        )
        lambda_valid = np.isfinite(lambda_raw) & (lambda_raw > 1.0e-12)
        lambda_fin = np.where(lambda_valid, lambda_raw, lambda_fallback)
        lambda_fin = np.clip(lambda_fin, 1.0e-12, 1.0e12)

        R_fin = self.contact_resistance + 1.0 / lambda_fin
        R_fin = np.nan_to_num(R_fin, nan=1.0e15, posinf=1.0e15, neginf=1.0e15)
        T_fin_eff = T_root - Q_fin * R_fin
        T_fin_eff = np.nan_to_num(T_fin_eff, nan=float(np.mean(T_space_safe)), posinf=1.0e12, neginf=-1.0e12)
        self.bc_fin.update_params(T_ext=T_fin_eff, R_ext=R_fin)

        self.last_fin_temperature = np.array(T_fin, copy=True)
        self.last_fin_radiation_distribution = np.array(Q_fin, copy=True)
        self.last_fin_net_from_root_distribution = np.array(Q_fin, copy=True)
        self.last_fin_conductance_distribution = np.array(lambda_fin, copy=True)
        self.last_fin_effective_temperature_distribution = np.array(T_fin_eff, copy=True)
        self.last_fin_equivalent_resistance_distribution = np.array(R_fin, copy=True)
        self.last_tube_radiation_distribution = (
            self.tube_emissivity
            * self.sigma
            * self.tube_bare_area
            * (
                np.asarray(self.wall.boundaries["right"].T_surface, dtype=float) ** 4
                - self.radiation_background_temperature ** 4
            )
        )

    def get_solids(self) -> list:
        return [self.wall]

    def get_couplers(self) -> list:
        return [self.coupler]

    def get_heat_exchange_breakdown(self) -> dict:
        T_wall = np.asarray(self.wall.boundaries["right"].T_surface, dtype=float)
        tube_radiation = (
            self.tube_emissivity
            * self.sigma
            * self.tube_bare_area
            * (T_wall ** 4 - self.radiation_background_temperature ** 4)
        )
        return {
            "bare_radiation": tube_radiation,
            "fin_radiation": np.array(self.last_fin_radiation_distribution, copy=True),
            "fin_net_from_root": np.array(self.last_fin_net_from_root_distribution, copy=True),
            "fin_conductance": np.array(self.last_fin_conductance_distribution, copy=True),
            "fin_effective_temperature": np.array(self.last_fin_effective_temperature_distribution, copy=True),
            "fin_equivalent_resistance": np.array(self.last_fin_equivalent_resistance_distribution, copy=True),
            "fin_iteration_count": self.last_fin_iteration_count,
            "fin_max_delta": self.last_fin_max_delta,
            "fin_used_warm_start": self.last_fin_used_warm_start,
            "gross_rejection": tube_radiation + self.last_fin_radiation_distribution,
            "net_rejection": tube_radiation + self.last_fin_net_from_root_distribution,
        }

    def get_temperature_distribution(self):
        return self.wall.T.reshape(self.wall.shape_nodes)

    def get_fin_temperature_distribution(self):
        return np.array(self.last_fin_temperature, copy=True)
