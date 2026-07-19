from __future__ import annotations

from typing import Iterable, Optional

import numpy as np

from Components.BaseComponent import BaseComponent


class RadiatorThermalShield(BaseComponent):
    """
    Quasi-steady equivalent thermal shield for pipe-fin radiator startup cases.

    The v1 model does not add a solid ODE. It updates each radiator unit's
    equivalent radiation background temperature so existing tube/fin radiation
    paths can reuse their normal boundary machinery.
    """

    sigma = 5.670374419e-8
    MODEL_SEGMENT_BALANCE = "segment_balance"
    MODEL_FORTRAN_SHIELD2 = "fortran_shield2"

    XSS_BASE = np.array([0.041496, 0.086026, 0.076388, 0.077237, 0.076389, 0.086026], dtype=float)
    XS_U = np.full(6, 0.056333, dtype=float)
    XS_D = np.full(6, 0.036094, dtype=float)
    XS_R1_BASE = np.array([0.168781211, 0.0176293, 0.008384135, 0.011616, 0.008383836, 0.0176293], dtype=float)
    XS_R2_BASE = np.array([0.19733, 0.01343104, 0.0019838, 0.0034083, 0.001984, 0.013432], dtype=float)
    XU_R = np.array([3.27e-3] * 6 + [2.51e-3] * 6, dtype=float)
    XD_R = np.array([2.37e-2, 2.37e-2, 2.36e-2, 2.37e-2, 2.37e-2, 2.37e-2] + [1.00e-1] * 6, dtype=float)
    XUD = 0.039436
    KSU = 2.739424494
    KSD = 1.126253671
    KUD = 0.411984
    KR1S = 4.9703
    KR2S = 4.1224
    XR1_U = 0.0044079
    XR1_D = 0.10374
    XR2_U = 0.0037477
    XR2_D = 0.30589

    def __init__(
            self,
            name: str,
            radiator_units: Iterable,
            active_until_s: Optional[float] = None,
            background_temperature_k: float = 3.0,
            shield_view_factor: float = 0.8,
            inner_emissivity: float = 0.8,
            outer_emissivity: float = 0.8,
            conductivity_w_m_k: float = 1.0,
            thickness_m: float = 0.002,
            solar_heat_flux_w_m2: float = 0.0,
            relaxation: float = 1.0,
            model: str = MODEL_SEGMENT_BALANCE,
            qsss_w_m2: Optional[Iterable[float]] = None,
            strict_fortran: bool = False,
            shield2_initial_temperature_k: float = 200.0,
            shield2_tol_k: float = 1.0e-3,
            shield2_max_iterations: int = 200):
        super().__init__(name)
        self.radiator_units = list(radiator_units)
        self.active_until_s = None if active_until_s is None else float(active_until_s)
        self.background_temperature_k = float(background_temperature_k)
        self.shield_view_factor = float(np.clip(shield_view_factor, 0.0, 1.0))
        self.inner_emissivity = max(float(inner_emissivity), 1.0e-12)
        self.outer_emissivity = max(float(outer_emissivity), 1.0e-12)
        self.conductivity_w_m_k = max(float(conductivity_w_m_k), 1.0e-12)
        self.thickness_m = max(float(thickness_m), 1.0e-12)
        solar_flux = np.asarray(solar_heat_flux_w_m2, dtype=float)
        if solar_flux.ndim > 1 or (solar_flux.ndim == 1 and solar_flux.size not in (1, len(self.radiator_units))):
            raise ValueError("solar_heat_flux_w_m2 must be scalar or match radiator_units.")
        self.solar_heat_flux_w_m2 = solar_flux
        self.relaxation = float(np.clip(relaxation, 0.0, 1.0))
        self.model = str(model)
        if self.model not in {self.MODEL_SEGMENT_BALANCE, self.MODEL_FORTRAN_SHIELD2}:
            raise ValueError(f"Unsupported radiation shield model: {self.model}")
        if qsss_w_m2 is None:
            self.qsss_w_m2 = np.zeros(8, dtype=float)
        else:
            self.qsss_w_m2 = np.asarray(qsss_w_m2, dtype=float)
            if self.qsss_w_m2.shape != (8,):
                raise ValueError("qsss_w_m2 must have shape (8,).")
        self.strict_fortran = bool(strict_fortran)
        self.shield2_initial_temperature_k = max(float(shield2_initial_temperature_k), 1.0e-3)
        self.shield2_tol_k = max(float(shield2_tol_k), 1.0e-12)
        self.shield2_max_iterations = max(1, int(shield2_max_iterations))

        self.last_active = False
        self.last_inner_temperature_mean_k = self.background_temperature_k
        self.last_outer_temperature_mean_k = self.background_temperature_k
        self.last_effective_background_mean_k = self.background_temperature_k
        self.last_q_from_radiator_w = 0.0
        self.last_q_solar_w = 0.0
        self.last_q_to_space_w = 0.0
        self.last_energy_residual_w = 0.0
        self.last_energy_residual_rel = 0.0
        self.last_solver_failures = 0
        self.last_background_by_unit = []
        self.last_shield2_tc4_rad_12 = np.zeros(12, dtype=float)
        self.last_shield2_qrr_12 = np.zeros(12, dtype=float)
        self.last_shield2_applied_qrr_12 = np.zeros(12, dtype=float)
        self.last_shield2_qrr_weight_12 = np.ones(12, dtype=float)
        self.last_shield2_qrr_bg4_12 = np.zeros(12, dtype=float)
        self.last_shield2_qrr_bgT_12 = np.zeros(12, dtype=float)
        self.last_shield2_inner_t4_8 = np.zeros(8, dtype=float)
        self.last_shield2_inner_temperature_8 = np.full(8, self.background_temperature_k, dtype=float)
        self.last_shield2_outer_temperature_8 = np.full(8, self.background_temperature_k, dtype=float)
        self.last_shield2_iteration_count = 0
        self.last_shield2_converged = False
        self.last_shield2_residual_k = 0.0
        self._shield2_xs_s = self._rotated_rows(self.XSS_BASE)
        self._shield2_xsr1 = self._rotated_rows(self.XS_R1_BASE)
        self._shield2_xsr2 = self._rotated_rows(self.XS_R2_BASE)
        self._shield2_xsr = np.hstack([self._shield2_xsr1, self._shield2_xsr2])
        self._shield2_matrix = self._build_shield2_matrix()

    @staticmethod
    def _rotated_rows(base: Iterable[float]) -> np.ndarray:
        base_arr = np.asarray(base, dtype=float)
        base12 = np.concatenate([base_arr, base_arr])
        rows = []
        for i in range(6):
            start = 6 - i
            rows.append(base12[start:start + 6])
        return np.asarray(rows, dtype=float)

    def _build_shield2_matrix(self) -> np.ndarray:
        matrix = np.zeros((8, 8), dtype=float)
        matrix[0:6, 0:6] = -self._shield2_xs_s
        for i in range(6):
            matrix[i, i] = 1.0 - self._shield2_xs_s[i, i]
            matrix[i, 6] = -self.XS_U[i]
            matrix[i, 7] = -self.XS_D[i]
            matrix[6, i] = -self.XS_U[i] * self.KSU
            matrix[7, i] = -self.XS_D[i] * self.KSD
        matrix[6, 6] = 1.0
        matrix[6, 7] = -self.XUD
        matrix[7, 6] = -self.XUD * self.KUD
        matrix[7, 7] = 1.0
        return matrix

    def _is_active(self, current_time: float) -> bool:
        if self.active_until_s is None:
            return True
        return float(current_time) <= self.active_until_s

    def _unit_radiation_area(self, unit, shape) -> np.ndarray:
        tube_area = np.asarray(getattr(unit, "tube_bare_area", np.ones(shape)), dtype=float)
        fin_area = np.asarray(getattr(unit, "fin_radiating_area", np.zeros(shape)), dtype=float)
        area = tube_area + fin_area
        area *= float(getattr(unit, "radiation_area_multiplier", 1.0))
        if area.shape != tuple(shape):
            area = np.broadcast_to(area, shape)
        area = np.nan_to_num(area, nan=0.0, posinf=0.0, neginf=0.0)
        return np.maximum(area, 0.0)

    def _surface_temperature(self, unit) -> np.ndarray:
        if hasattr(unit, "get_radiation_surface_temperature"):
            return np.asarray(unit.get_radiation_surface_temperature(), dtype=float)
        wall = getattr(unit, "wall", None)
        if wall is not None and "right" in getattr(wall, "boundaries", {}):
            return np.asarray(wall.boundaries["right"].T_surface, dtype=float)
        raise AttributeError(f"Radiator unit {getattr(unit, 'name', unit)!r} lacks radiation surface temperature.")

    def _set_unit_background(self, unit, background: np.ndarray) -> None:
        if hasattr(unit, "set_radiation_background_temperature"):
            unit.set_radiation_background_temperature(background)
            return
        raise AttributeError(f"Radiator unit {getattr(unit, 'name', unit)!r} lacks radiation background setter.")

    def _shield2_tube_t4(self, unit) -> float:
        fin_temperature = getattr(unit, "last_fin_temperature", None)
        has_valid_fin = bool(getattr(unit, "_has_valid_fin_temperature", True))
        if fin_temperature is not None and has_valid_fin:
            fin = np.asarray(fin_temperature, dtype=float)
            finite = np.isfinite(fin)
            if fin.size and np.any(finite) and float(np.nanmax(fin)) > 1.0:
                return float(np.mean(np.maximum(fin[finite], 1.0e-3) ** 4))
        surface = self._surface_temperature(unit)
        return float(np.mean(np.maximum(surface, 1.0e-3) ** 4))

    def _shield2_nonheating_limit_t4(self, unit, surface: np.ndarray) -> float:
        limits = [float(np.min(np.maximum(surface, 1.0e-3) ** 4))]
        fin_temperature = getattr(unit, "last_fin_temperature", None)
        has_valid_fin = bool(getattr(unit, "_has_valid_fin_temperature", True))
        if fin_temperature is not None and has_valid_fin:
            fin = np.asarray(fin_temperature, dtype=float)
            finite = np.isfinite(fin)
            if fin.size and np.any(finite) and float(np.nanmax(fin)) > 1.0:
                limits.append(float(np.min(np.maximum(fin[finite], 1.0e-3) ** 4)))
        return max(min(limits), 1.0e-12)

    @staticmethod
    def _shield2_sector_indices(n_tubes: int) -> np.ndarray:
        if n_tubes <= 0:
            return np.zeros(0, dtype=int)
        sectors = np.floor(np.arange(n_tubes, dtype=float) * 12.0 / float(n_tubes)).astype(int)
        return np.clip(sectors, 0, 11)

    def _shield2_tc4_rad(self) -> tuple[np.ndarray, np.ndarray]:
        n_tubes = len(self.radiator_units)
        sectors = self._shield2_sector_indices(n_tubes)
        tube_t4 = np.asarray([self._shield2_tube_t4(unit) for unit in self.radiator_units], dtype=float)
        tc4_rad = np.zeros(12, dtype=float)
        for sector in range(12):
            values = tube_t4[sectors == sector]
            tc4_rad[sector] = float(np.mean(values)) if values.size else 0.0
        return tc4_rad, sectors

    def _shield2_outer_temperature(self, inner_temperature: float, qsss: float) -> float:
        t4_coeff = 1.0 / (self.thickness_m / self.conductivity_w_m_k * self.outer_emissivity * self.sigma)
        inner_temperature = max(float(inner_temperature), 1.0e-3)
        qsss = float(qsss)

        def residual(outer_temperature: float) -> float:
            return outer_temperature ** 4 + t4_coeff * outer_temperature - qsss / self.sigma - inner_temperature * t4_coeff

        lower = 1.0e-9
        upper = max(inner_temperature, self.background_temperature_k, self.shield2_initial_temperature_k, 300.0)
        while residual(upper) < 0.0:
            upper *= 2.0
            if upper > 1.0e6:
                break
        for _ in range(100):
            mid = 0.5 * (lower + upper)
            if residual(mid) > 0.0:
                upper = mid
            else:
                lower = mid
        return float(0.5 * (lower + upper))

    def _solve_fortran_shield2(self, tc4_rad: np.ndarray, qsss: np.ndarray) -> dict:
        tc4 = np.asarray(tc4_rad, dtype=float)
        if tc4.shape != (12,):
            raise ValueError("tc4_rad must have shape (12,).")
        qsss = np.asarray(qsss, dtype=float)
        if qsss.shape != (8,):
            raise ValueError("qsss must have shape (8,).")

        inner_temperature = np.full(8, self.shield2_initial_temperature_k, dtype=float)
        y_inner = inner_temperature ** 4
        outer_temperature = np.array(inner_temperature, copy=True)
        converged = False
        residual_k = 0.0
        iteration_count = 0

        base_b = np.zeros(8, dtype=float)
        for i in range(6):
            base_b[i] = float(np.sum(self._shield2_xsr[i, :] * tc4))
        base_b[6] = float(np.sum(self.XU_R * tc4))
        base_b[7] = float(np.sum(self.XD_R * tc4))

        for iteration in range(1, self.shield2_max_iterations + 1):
            qqs = np.zeros(8, dtype=float)
            for i in range(7):
                outer_temperature[i] = self._shield2_outer_temperature(inner_temperature[i], qsss[i])
                qqs[i] = outer_temperature[i] ** 4 - qsss[i] / self.sigma
            qqs[7] = 0.0
            b = base_b - qqs
            y_inner = np.linalg.solve(self._shield2_matrix, b)
            y_inner = np.nan_to_num(y_inner, nan=1.0e-12, posinf=1.0e24, neginf=1.0e-12)
            y_inner = np.maximum(y_inner, 1.0e-12)
            next_inner_temperature = y_inner ** 0.25
            residual_k = float(np.sum(np.abs(next_inner_temperature - inner_temperature)))
            inner_temperature = next_inner_temperature
            iteration_count = iteration
            if residual_k < self.shield2_tol_k:
                converged = True
                break

        for i in range(7):
            outer_temperature[i] = self._shield2_outer_temperature(inner_temperature[i], qsss[i])
        outer_temperature[7] = inner_temperature[7]

        qrr = np.zeros(12, dtype=float)
        qrr_weight = np.zeros(12, dtype=float)
        for j in range(12):
            for mm in range(6):
                if j < 6:
                    qrr[j] += self._shield2_xsr[mm, j] * y_inner[mm] * self.KR1S
                    qrr_weight[j] += self._shield2_xsr[mm, j] * self.KR1S
                else:
                    qrr[j] += self._shield2_xsr[mm, j] * y_inner[mm] * self.KR2S
                    qrr_weight[j] += self._shield2_xsr[mm, j] * self.KR2S
            if j < 6:
                qrr[j] += self.XR1_U * y_inner[6] + self.XR1_D * y_inner[7]
                qrr_weight[j] += self.XR1_U + self.XR1_D
            else:
                upper_factor = self.XR1_U if self.strict_fortran else self.XR2_U
                qrr[j] += upper_factor * y_inner[6] + self.XR2_D * y_inner[7]
                qrr_weight[j] += upper_factor + self.XR2_D
        qrr = np.maximum(np.nan_to_num(qrr, nan=1.0e-12, posinf=1.0e24, neginf=1.0e-12), 1.0e-12)
        qrr_weight = np.maximum(
            np.nan_to_num(qrr_weight, nan=1.0e-12, posinf=1.0e24, neginf=1.0e-12),
            1.0e-12,
        )
        qrr_bg4 = np.maximum(qrr / qrr_weight, 1.0e-12)
        qrr_bgT = qrr_bg4 ** 0.25

        return {
            "tc4_rad": tc4,
            "qrr": qrr,
            "qrr_weight": qrr_weight,
            "qrr_bg4": qrr_bg4,
            "qrr_bgT": qrr_bgT,
            "y_inner": y_inner,
            "inner_temperature": inner_temperature,
            "outer_temperature": outer_temperature,
            "iteration_count": iteration_count,
            "converged": converged,
            "residual_k": residual_k,
        }

    def _solve_segment_balance(
            self,
            radiator_temperature: float,
            area: float,
            solar_heat_flux_w_m2: float) -> tuple[float, float, float, float, float, float, bool]:
        area = float(area)
        background = max(self.background_temperature_k, 1.0e-3)
        radiator_temperature = max(float(radiator_temperature), 1.0e-3)
        if area <= 0.0:
            return background, background, background, 0.0, 0.0, 0.0, True

        view = self.shield_view_factor
        q_solar = float(solar_heat_flux_w_m2) * area
        if view <= 0.0:
            outer_t4 = background ** 4 + q_solar / (self.outer_emissivity * self.sigma * area)
            outer_temperature = np.power(max(outer_t4, 1.0e-12), 0.25)
            return background, outer_temperature, background, 0.0, q_solar, q_solar, True

        conductance = self.conductivity_w_m_k * area / self.thickness_m

        def values(inner_temperature: float) -> tuple[float, float, float, float]:
            inner_temperature = max(float(inner_temperature), 1.0e-3)
            q_in = (
                self.inner_emissivity
                * self.sigma
                * area
                * view
                * (radiator_temperature ** 4 - inner_temperature ** 4)
            )
            outer_temperature = inner_temperature - q_in / conductance
            outer_t4 = max(outer_temperature, 1.0e-3) ** 4
            q_out = self.outer_emissivity * self.sigma * area * (outer_t4 - background ** 4)
            residual = q_in + q_solar - q_out
            return q_in, outer_temperature, q_out, residual

        lower = 1.0e-3
        solar_outer_t = np.power(
            max(background ** 4 + max(q_solar, 0.0) / (self.outer_emissivity * self.sigma * area), 1.0e-12),
            0.25,
        )
        upper = max(radiator_temperature, background, solar_outer_t, 300.0) + 10.0
        _, _, _, f_lower = values(lower)
        _, _, _, f_upper = values(upper)
        expansions = 0
        while f_upper > 0.0 and expansions < 20:
            upper *= 1.5
            _, _, _, f_upper = values(upper)
            expansions += 1

        solved = np.isfinite(f_lower) and np.isfinite(f_upper) and f_lower >= 0.0 and f_upper <= 0.0
        if solved:
            inner_temperature = upper
            for _ in range(80):
                mid = 0.5 * (lower + upper)
                q_in, outer_temperature, q_out, residual = values(mid)
                inner_temperature = mid
                if abs(residual) <= 1.0e-10 * max(abs(q_in) + abs(q_solar) + abs(q_out), 1.0):
                    break
                if residual > 0.0:
                    lower = mid
                else:
                    upper = mid
            q_in, outer_temperature, q_out, residual = values(inner_temperature)
        else:
            inner_temperature = np.clip(radiator_temperature, background, max(radiator_temperature, background))
            q_in, outer_temperature, q_out, residual = values(inner_temperature)

        effective_t4 = (1.0 - view) * background ** 4 + view * max(inner_temperature, 1.0e-3) ** 4
        effective_background = np.power(max(effective_t4, 1.0e-12), 0.25)
        return inner_temperature, outer_temperature, effective_background, q_in, q_solar, q_out, solved

    def _compute_effective_background(
            self,
            surface_temperature: np.ndarray,
            area: np.ndarray,
            solar_heat_flux_w_m2: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
        surface = np.maximum(np.asarray(surface_temperature, dtype=float), 1.0e-3)
        area = np.maximum(np.asarray(area, dtype=float), 0.0)
        inner = np.empty_like(surface)
        outer = np.empty_like(surface)
        effective = np.empty_like(surface)
        q_in = np.empty_like(surface)
        q_solar = np.empty_like(surface)
        q_out = np.empty_like(surface)
        failures = 0
        for index in np.ndindex(surface.shape):
            (
                inner[index],
                outer[index],
                effective[index],
                q_in[index],
                q_solar[index],
                q_out[index],
                solved,
            ) = self._solve_segment_balance(surface[index], area[index], solar_heat_flux_w_m2)
            if not solved:
                failures += 1

        if self.last_background_by_unit:
            previous = self.last_background_by_unit.pop(0)
            if np.shape(previous) != np.shape(effective):
                previous = None
            if previous is not None and self.relaxation < 1.0:
                effective = previous + self.relaxation * (effective - previous)
        return inner, outer, effective, q_in, q_solar, q_out, failures

    def pre_step(self, dt: float, current_time: float):
        active = self._is_active(current_time)
        self.last_active = bool(active)
        if self.model == self.MODEL_FORTRAN_SHIELD2:
            self._pre_step_fortran_shield2(active)
            return
        next_backgrounds = []
        inner_temperatures = []
        outer_temperatures = []
        effective_backgrounds = []
        q_from = 0.0
        q_solar = 0.0
        q_to_space = 0.0
        solver_failures = 0

        previous_backgrounds = list(self.last_background_by_unit)
        self.last_background_by_unit = previous_backgrounds

        solar_flux = np.broadcast_to(self.solar_heat_flux_w_m2, (len(self.radiator_units),))
        for unit_index, unit in enumerate(self.radiator_units):
            surface = self._surface_temperature(unit)
            area = self._unit_radiation_area(unit, surface.shape)
            if active:
                (
                    inner_temperature,
                    outer_temperature,
                    effective_background,
                    q_in_distribution,
                    q_solar_distribution,
                    q_out_distribution,
                    failures,
                ) = self._compute_effective_background(surface, area, float(solar_flux[unit_index]))
                solver_failures += failures
            else:
                effective_background = np.full_like(surface, max(self.background_temperature_k, 1.0e-3))
                inner_temperature = effective_background
                outer_temperature = effective_background
                q_in_distribution = np.zeros_like(surface)
                q_solar_distribution = np.zeros_like(surface)
                q_out_distribution = np.zeros_like(surface)

            self._set_unit_background(unit, effective_background)
            next_backgrounds.append(np.array(effective_background, copy=True))
            inner_temperatures.append(inner_temperature)
            outer_temperatures.append(outer_temperature)
            effective_backgrounds.append(effective_background)
            q_from += float(np.sum(q_in_distribution))
            q_solar += float(np.sum(q_solar_distribution))
            q_to_space += float(np.sum(q_out_distribution))

        self.last_background_by_unit = next_backgrounds
        if inner_temperatures:
            inner_all = np.concatenate([np.ravel(arr) for arr in inner_temperatures])
            outer_all = np.concatenate([np.ravel(arr) for arr in outer_temperatures])
            effective_all = np.concatenate([np.ravel(arr) for arr in effective_backgrounds])
            self.last_inner_temperature_mean_k = float(np.mean(inner_all))
            self.last_outer_temperature_mean_k = float(np.mean(outer_all))
            self.last_effective_background_mean_k = float(np.mean(effective_all))
        else:
            self.last_inner_temperature_mean_k = self.background_temperature_k
            self.last_outer_temperature_mean_k = self.background_temperature_k
            self.last_effective_background_mean_k = self.background_temperature_k
        self.last_q_from_radiator_w = q_from
        self.last_q_solar_w = q_solar
        self.last_q_to_space_w = q_to_space
        self.last_energy_residual_w = q_from + q_solar - q_to_space
        denominator = max(abs(q_from) + abs(q_solar) + abs(q_to_space), 1.0)
        self.last_energy_residual_rel = abs(self.last_energy_residual_w) / denominator
        self.last_solver_failures = int(solver_failures)

    def _pre_step_fortran_shield2(self, active: bool) -> None:
        if not active:
            default_background = max(self.background_temperature_k, 1.0e-3)
            for unit in self.radiator_units:
                self._set_unit_background(unit, np.full_like(self._surface_temperature(unit), default_background))
            self.last_inner_temperature_mean_k = default_background
            self.last_outer_temperature_mean_k = default_background
            self.last_effective_background_mean_k = default_background
            self.last_q_from_radiator_w = 0.0
            self.last_q_solar_w = 0.0
            self.last_q_to_space_w = 0.0
            self.last_energy_residual_w = 0.0
            self.last_energy_residual_rel = 0.0
            self.last_solver_failures = 0
            self.last_shield2_iteration_count = 0
            self.last_shield2_converged = False
            self.last_shield2_residual_k = 0.0
            return

        tc4_rad, sectors = self._shield2_tc4_rad()
        result = self._solve_fortran_shield2(tc4_rad, self.qsss_w_m2)
        qrr = result["qrr"]
        qrr_bg4 = result["qrr_bg4"]
        qrr_bgT = result["qrr_bgT"]

        q_from = 0.0
        for unit_index, unit in enumerate(self.radiator_units):
            sector = int(sectors[unit_index])
            surface = self._surface_temperature(unit)
            background = np.full_like(surface, qrr_bgT[sector], dtype=float)
            self._set_unit_background(unit, background)
            area = self._unit_radiation_area(unit, surface.shape)
            q_from += float(np.sum(
                self.inner_emissivity
                * self.sigma
                * area
                * (np.maximum(surface, 1.0e-3) ** 4 - qrr_bg4[sector])
            ))

        self.last_shield2_tc4_rad_12 = np.array(result["tc4_rad"], copy=True)
        self.last_shield2_qrr_12 = np.array(qrr, copy=True)
        self.last_shield2_qrr_weight_12 = np.array(result["qrr_weight"], copy=True)
        self.last_shield2_qrr_bg4_12 = np.array(qrr_bg4, copy=True)
        self.last_shield2_qrr_bgT_12 = np.array(qrr_bgT, copy=True)
        self.last_shield2_applied_qrr_12 = np.array(qrr_bg4, copy=True)
        self.last_shield2_inner_t4_8 = np.array(result["y_inner"], copy=True)
        self.last_shield2_inner_temperature_8 = np.array(result["inner_temperature"], copy=True)
        self.last_shield2_outer_temperature_8 = np.array(result["outer_temperature"], copy=True)
        self.last_shield2_iteration_count = int(result["iteration_count"])
        self.last_shield2_converged = bool(result["converged"])
        self.last_shield2_residual_k = float(result["residual_k"])
        self.last_inner_temperature_mean_k = float(np.mean(self.last_shield2_inner_temperature_8))
        self.last_outer_temperature_mean_k = float(np.mean(self.last_shield2_outer_temperature_8))
        self.last_effective_background_mean_k = float(np.mean(qrr_bgT))
        self.last_q_from_radiator_w = q_from
        self.last_q_solar_w = 0.0
        self.last_q_to_space_w = q_from
        self.last_energy_residual_w = 0.0
        self.last_energy_residual_rel = 0.0
        self.last_solver_failures = 0 if self.last_shield2_converged else 1

    def get_diagnostics(self) -> dict:
        return {
            "radiation_shield_model": self.model,
            "radiation_shield_active": bool(self.last_active),
            "radiation_shield_effective_background_mean_k": float(self.last_effective_background_mean_k),
            "radiation_shield_inner_temperature_mean_k": float(self.last_inner_temperature_mean_k),
            "radiation_shield_outer_temperature_mean_k": float(self.last_outer_temperature_mean_k),
            "radiation_shield_q_from_radiator_w": float(self.last_q_from_radiator_w),
            "radiation_shield_q_solar_w": float(self.last_q_solar_w),
            "radiation_shield_q_to_space_w": float(self.last_q_to_space_w),
            "radiation_shield_energy_residual_w": float(self.last_energy_residual_w),
            "radiation_shield_energy_residual_rel": float(self.last_energy_residual_rel),
            "radiation_shield_solver_failures": int(self.last_solver_failures),
            "radiation_shield_view_factor": float(self.shield_view_factor),
            "radiation_shield_conductivity_w_m_k": float(self.conductivity_w_m_k),
            "radiation_shield_thickness_m": float(self.thickness_m),
            "radiation_shield_tc4_rad_12": self.last_shield2_tc4_rad_12.tolist(),
            "radiation_shield_qrr_12": self.last_shield2_qrr_12.tolist(),
            "radiation_shield_applied_qrr_12": self.last_shield2_applied_qrr_12.tolist(),
            "radiation_shield_qrr_weight_12": self.last_shield2_qrr_weight_12.tolist(),
            "radiation_shield_qrr_bg4_12": self.last_shield2_qrr_bg4_12.tolist(),
            "radiation_shield_qrr_bgT_12": self.last_shield2_qrr_bgT_12.tolist(),
            "radiation_shield_inner_t4_8": self.last_shield2_inner_t4_8.tolist(),
            "radiation_shield_inner_temperature_8_k": self.last_shield2_inner_temperature_8.tolist(),
            "radiation_shield_outer_temperature_8_k": self.last_shield2_outer_temperature_8.tolist(),
            "radiation_shield_iteration_count": int(self.last_shield2_iteration_count),
            "radiation_shield_converged": bool(self.last_shield2_converged),
            "radiation_shield_residual_k": float(self.last_shield2_residual_k),
            "radiation_shield_strict_fortran": bool(self.strict_fortran),
        }
