import numpy as np
from typing import Union
from Materials.Base import SolidMaterial

# 定义类型别名以支持向量化提示
Numeric = Union[float, np.ndarray]


class SodiumHP(SolidMaterial):
    """
    热管液态金属工质：钠 (Sodium, Na)
    引入等效热容法 (Apparent Heat Capacity Method) 处理固液相变 (糊状区 Mushy Zone)
    """

    def __init__(self, name: str = "Sodium_Wick_Fluid"):
        super().__init__(name=name)
        # 基础物理常数
        self.molar_mass = 0.02299  # 摩尔质量 [kg/mol] (精确补充)

        # 相变(融化)参数
        self.T_melt = 371.0  # 融化温度 [K]
        self.H_sf = 114.7e3  # 熔化潜热 [J/kg]
        self.dt_thaw = 1.0  # 糊状区半宽 [K] (总宽度为 2K)

        self.T_crit = 2503.7  # 临界温度 [K]

    def _get_mushy_fraction(self, T: np.ndarray) -> np.ndarray:
        """计算糊状区插值因子 f, 范围 [0, 1]"""
        f = (T - self.T_melt + self.dt_thaw) / (2.0 * self.dt_thaw)
        return np.clip(f, 0.0, 1.0)

    def density(self, T: Numeric) -> Numeric:
        """密度 [kg/m^3] (固液两相及糊状区插值)"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.15

        # 1. 固态密度计算
        rho_solid = np.where(
            T < 273.15,
            0.9725,
            0.9725 - 20.11e-5 * TC - 1.5e-7 * TC ** 2
        ) * 1e3

        # 2. 液态密度计算
        T_safe = np.minimum(T, self.T_crit)
        rho_liquid = 219.0 + 275.32 * (1.0 - T_safe / self.T_crit) + 511.58 * np.sqrt(1.0 - T_safe / self.T_crit)

        # 3. 糊状区平滑过渡
        f = self._get_mushy_fraction(T)
        rho = rho_solid + f * (rho_liquid - rho_solid)

        # 4. 根据温度区间返回
        return np.where(T < self.T_melt - self.dt_thaw, rho_solid,
                        np.where(T > self.T_melt + self.dt_thaw, rho_liquid, rho))

    def conductivity(self, T: Numeric) -> Numeric:
        """导热系数 [W/(m*K)] (固液两相及糊状区插值)"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.15

        # 1. 固态导热
        cond_solid = np.where(
            T <= 273.15,
            135.0,
            (1.356 - 1.67e-3 * TC) * 1e2
        )

        # 2. 液态导热
        cond_liquid = 124.67 - 0.11381 * T + 5.5226e-5 * T ** 2 - 1.1842e-8 * T ** 3

        # 3. 糊状区插值 (cond_s + f * (cond_l - cond_s))
        f = self._get_mushy_fraction(T)
        cond_mushy = cond_solid + f * (cond_liquid - cond_solid)

        return np.where(T < self.T_melt - self.dt_thaw, cond_solid,
                        np.where(T > self.T_melt + self.dt_thaw, cond_liquid, cond_mushy))

    def heat_capacity(self, T: Numeric) -> Numeric:
        """定压比热容 [J/(kg*K)] (固液相变潜热等效热容法)"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.15

        # 1. 固态显热容
        capa_solid = np.where(
            T <= 273.15,
            0.2865 * 4.1868,
            (0.2865 + 1.551e-4 * TC + 0.2516e-5 * TC ** 2) * 4.1868
        ) * 1e3

        # 2. 液态显热容
        capa_liquid = 1630.138954 - 0.833500608 * T + 0.000462815872 * T ** 2

        # 3. 糊状区等效热容 = 潜热释放等效部分 + 显热部分平滑插值
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
        T_safe = np.maximum(T, 1e-3)
        return 2.29e11 * (10.0 ** (-5567.0 / T_safe)) / (T_safe ** 0.5)

    def vapor_viscosity(self, T: Numeric) -> Numeric:
        """蒸汽动力粘度 [Pa*s]"""
        T = np.asarray(T, dtype=float)
        return 6.083e-9 * T + 1.2606e-5

    def latent_heat(self, T: Numeric) -> Numeric:
        """汽化潜热 [J/kg]"""
        T_safe = np.minimum(np.asarray(T, dtype=float), self.T_crit)
        Vhfg = 393.37 * (1.0 - T_safe / self.T_crit) + 4398.6 * ((1.0 - T_safe / self.T_crit) ** 0.29302)
        return 1000.0 * Vhfg

