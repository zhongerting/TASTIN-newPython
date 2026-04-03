import numpy as np
from typing import Union
from Materials.Base import SolidMaterial

# 定义类型别名以支持向量化提示
Numeric = Union[float, np.ndarray]


class PotassiumHP(SolidMaterial):
    """
    热管液态金属工质：钾 (Potassium, K)
    引入等效热容法 (Apparent Heat Capacity Method) 处理固液相变 (糊状区 Mushy Zone)
    """

    def __init__(self, name: str = "Potassium_Wick_Fluid"):
        super().__init__(name=name)
        # 基础物理常数
        self.molar_mass = 0.039098  # 摩尔质量 [kg/mol]

        # 相变(融化)参数
        self.T_melt = 336.8  # 融化温度 [K]
        self.H_sf = 5.932e4  # 熔化潜热 [J/kg]
        self.dt_thaw = 1.0  # 糊状区半宽 [K] (总宽度为 2K)

    def _get_mushy_fraction(self, T: np.ndarray) -> np.ndarray:
        """计算糊状区插值因子 f, 范围 [0, 1]"""
        f = (T - self.T_melt + self.dt_thaw) / (2.0 * self.dt_thaw)
        return np.clip(f, 0.0, 1.0)

    def density(self, T: Numeric) -> Numeric:
        """密度 [kg/m^3]"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.15

        rho_solid = 864.0 - 0.24162 * TC
        rho_liquid = 949.1192 - 0.24463 * T

        f = self._get_mushy_fraction(T)
        rho_mushy = rho_solid + f * (rho_liquid - rho_solid)

        return np.where(T < self.T_melt - self.dt_thaw, rho_solid,
                        np.where(T > self.T_melt + self.dt_thaw, rho_liquid, rho_mushy))

    def conductivity(self, T: Numeric) -> Numeric:
        """导热系数 [W/(m*K)]"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.15

        cond_solid = 126.0 - 6.028e-2 * T

        T_safe = np.maximum(T, 1e-3)
        cond_liquid = 43.8 - 2.22e-2 * TC + 3.95e3 / T_safe

        f = self._get_mushy_fraction(T)
        cond_mushy = cond_solid + f * (cond_liquid - cond_solid)

        return np.where(T < self.T_melt - self.dt_thaw, cond_solid,
                        np.where(T > self.T_melt + self.dt_thaw, cond_liquid, cond_mushy))

    def heat_capacity(self, T: Numeric) -> Numeric:
        """定压比热容 [J/(kg*K)]"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.15

        capa_solid = 537.9 + 0.8002 * T
        capa_liquid = 1000.0 * (0.8389 - 3.6741e-4 * TC + 4.592e-7 * TC ** 2)

        f = self._get_mushy_fraction(T)
        base_capa_mushy = capa_solid + f * (capa_liquid - capa_solid)
        latent_capa = self.H_sf / (2.0 * self.dt_thaw)
        capa_mushy = base_capa_mushy + latent_capa

        return np.where(T < self.T_melt - self.dt_thaw, capa_solid,
                        np.where(T > self.T_melt + self.dt_thaw, capa_liquid, capa_mushy))

    # =====================================================================
    # 赝热导率专用扩展接口
    # =====================================================================

    def saturation_pressure(self, T: Numeric) -> Numeric:
        """饱和蒸汽压 [Pa]"""
        T_safe = np.maximum(np.asarray(T, dtype=float), 1e-3)
        return 4.0168e11 * (10.0 ** (-4625.3 / T_safe)) / (T_safe ** 0.7)

    def vapor_viscosity(self, T: Numeric) -> Numeric:
        """蒸汽动力粘度 [Pa*s]"""
        T = np.asarray(T, dtype=float)
        return 4.86372e-6 + 1.15683e-8 * T

    def latent_heat(self, T: Numeric) -> Numeric:
        """汽化潜热 [J/kg]"""
        T = np.asarray(T, dtype=float)
        return np.full_like(T, 1983.7e3)
