import numpy as np
from typing import Union

from Materials.Base import SolidMaterial

Numeric = Union[float, np.ndarray]


class WickMaterial(SolidMaterial):
    """
    热管吸液芯复合材料。
    - 密度取工质密度
    - 比热按体积热容加权
    - 导热系数 = 结构等效导热 + 赝热导
    """

    def __init__(self,
                 name: str,
                 solid_mat: SolidMaterial,
                 fluid_mat: SolidMaterial,
                 porosity: float,
                 r_vapor: float,
                 r_in_wall: float):
        super().__init__(name=name)

        self.solid = solid_mat
        self.fluid = fluid_mat
        self.porosity = porosity
        self.r_vapor = r_vapor
        self.r_in_wall = r_in_wall

        self.R_gas = 8.314
        self.k_boltzmann = 1.380649e-23
        self.conductivity_cap = 5.0e7

        self.use_lookup_table = True
        self.lookup_table_temp_min = 250.0
        self.lookup_table_temp_max = min(float(getattr(self.fluid, "T_crit", 2503.7)), 2500.0)
        self.lookup_table_num_points = 8192

        self._lookup_table_ready = False
        self._lookup_temperature_grid = None
        self._lookup_log_p_sat = None
        self._lookup_mu_v = None
        self._lookup_h_fg = None
        self._lookup_log_k_total = None
        self._lookup_high_nonlinear_min = None
        self._lookup_high_nonlinear_max = None

    def density(self, T: Numeric) -> Numeric:
        return self.fluid.density(T)

    def heat_capacity(self, T: Numeric) -> Numeric:
        phi = self.porosity

        rho_f = self.fluid.density(T)
        cp_f = self.fluid.heat_capacity(T)
        rho_s = self.solid.density(T)
        cp_s = self.solid.heat_capacity(T)

        capa_wf = rho_f * cp_f
        capa_s = rho_s * cp_s
        capa_eff = phi * capa_wf + (1.0 - phi) * capa_s

        return capa_eff / rho_f

    def conductivity(self, T: Numeric) -> Numeric:
        if not self.use_lookup_table:
            return self._conductivity_direct(T)

        self._ensure_lookup_table()
        return self._lookup_k_total(T)

    def is_high_nonlinearity_temperature(self, T: Numeric) -> Numeric:
        self._ensure_lookup_table()

        T_arr = np.asarray(T, dtype=float)
        if self._lookup_high_nonlinear_min is None or self._lookup_high_nonlinear_max is None:
            result = np.zeros_like(T_arr, dtype=bool)
        else:
            result = (
                (T_arr >= self._lookup_high_nonlinear_min)
                & (T_arr <= self._lookup_high_nonlinear_max)
            )

        if np.isscalar(T):
            return bool(result)
        return result

    def _conductivity_direct(self, T: Numeric) -> Numeric:
        phi = self.porosity

        k_f = self.fluid.conductivity(T)
        k_s = self.solid.conductivity(T)

        with np.errstate(divide='ignore', invalid='ignore'):
            numerator = k_s + k_f - (1.0 - phi) * (k_f - k_s)
            denominator = k_f + k_s + (1.0 - phi) * (k_f - k_s)
            k_self = k_f * (numerator / denominator)
            k_self = np.nan_to_num(k_self, nan=0.0)

        P_sat = self.fluid.saturation_pressure(T)
        mu_v = self.fluid.vapor_viscosity(T)
        h_fg = self.fluid.latent_heat(T)
        pse1 = self._pse1_conductivity(T, P_sat, mu_v, h_fg)

        T_safe = np.maximum(T, 1.0e-3)
        P_safe = np.maximum(P_sat, 1.0e-12)
        v_ave = np.sqrt((8.0 * self.R_gas * T_safe) / (np.pi * self.fluid.molar_mass))

        with np.errstate(divide='ignore', invalid='ignore'):
            parameter1 = (mu_v * v_ave) / (self.r_vapor * P_safe)
            parameter2 = 1.0 + (8.0 / 3.0) * parameter1
            k_pse_total = pse1 / parameter2
            k_pse_total = np.nan_to_num(k_pse_total, nan=0.0)

        k_total = k_self + k_pse_total
        return np.clip(k_total, a_min=0.0, a_max=self.conductivity_cap)

    def _ensure_lookup_table(self):
        if self._lookup_table_ready:
            return

        t_min = float(self.lookup_table_temp_min)
        t_max = float(self.lookup_table_temp_max)
        if t_max <= t_min:
            raise ValueError(f"Invalid WickMaterial lookup range: [{t_min}, {t_max}]")

        T_grid = np.linspace(t_min, t_max, int(self.lookup_table_num_points), dtype=float)

        P_sat = np.asarray(self.fluid.saturation_pressure(T_grid), dtype=float)
        mu_v = np.asarray(self.fluid.vapor_viscosity(T_grid), dtype=float)
        h_fg = np.asarray(self.fluid.latent_heat(T_grid), dtype=float)
        k_total = np.asarray(self._conductivity_direct(T_grid), dtype=float)

        self._lookup_temperature_grid = T_grid
        self._lookup_log_p_sat = np.log(np.maximum(P_sat, 1.0e-30))
        self._lookup_mu_v = np.maximum(mu_v, 1.0e-30)
        self._lookup_h_fg = np.maximum(h_fg, 1.0e-30)
        self._lookup_log_k_total = np.log(np.maximum(k_total, 1.0e-30))

        dlogk_dT = np.gradient(self._lookup_log_k_total, T_grid)
        grad_max = float(np.max(np.abs(dlogk_dT))) if dlogk_dT.size > 0 else 0.0
        if grad_max > 0.0:
            nonlinear_mask = np.abs(dlogk_dT) >= 0.05 * grad_max
            if np.any(nonlinear_mask):
                self._lookup_high_nonlinear_min = float(T_grid[np.argmax(nonlinear_mask)])
                self._lookup_high_nonlinear_max = float(
                    T_grid[len(T_grid) - 1 - np.argmax(nonlinear_mask[::-1])]
                )

        self._lookup_table_ready = True

    def _interp_table(self, T: Numeric, table: np.ndarray, *, log_space: bool = False) -> Numeric:
        T_arr = np.asarray(T, dtype=float)
        T_clipped = np.clip(T_arr, self.lookup_table_temp_min, self.lookup_table_temp_max)
        interp_flat = np.interp(T_clipped.reshape(-1), self._lookup_temperature_grid, table)
        if log_space:
            interp_flat = np.exp(interp_flat)

        result = interp_flat.reshape(T_arr.shape)
        if np.isscalar(T):
            return float(result)
        return result

    def _lookup_k_total(self, T: Numeric) -> Numeric:
        return self._interp_table(T, self._lookup_log_k_total, log_space=True)

    def _pse1_conductivity(self, T: Numeric, P_sat: Numeric, mu_v: Numeric, h_fg: Numeric) -> Numeric:
        T_safe = np.maximum(T, 1.0e-3)
        mu_safe = np.maximum(mu_v, 1.0e-12)

        Rv = self.r_vapor
        Rw = self.r_in_wall
        Mg = self.fluid.molar_mass
        R_gas = self.R_gas

        parameter1 = (Rv ** 4) / (Rw ** 2 - Rv ** 2)

        with np.errstate(divide='ignore', invalid='ignore'):
            parameter2 = (h_fg * Mg * P_sat) ** 2
            parameter3 = 4.0 * mu_safe * (R_gas ** 2) * (T_safe ** 3)
            pse1_cond_temp = parameter1 * parameter2 / parameter3
            pse1_cond_temp = np.nan_to_num(pse1_cond_temp, nan=0.0, posinf=self.conductivity_cap)

        return pse1_cond_temp
