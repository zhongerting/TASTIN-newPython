import numpy as np
from typing import Union
from Materials.Base import SolidMaterial

# 定义类型别名以支持向量化提示
Numeric = Union[float, np.ndarray]


class SS321(SolidMaterial):
    """321 不锈钢 (SS321)"""

    def __init__(self, name: str = "SS321"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        """密度 [kg/m^3]"""
        T = np.asarray(T, dtype=float)
        return np.full_like(T, 8090.0)

    def conductivity(self, T: Numeric) -> Numeric:
        """导热系数 [W/(m*K)]"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.0
        return 14.5 + 0.015 * TC

    def heat_capacity(self, T: Numeric) -> Numeric:
        """定压比热容 [J/(kg*K)]"""
        T = np.asarray(T, dtype=float)
        return np.full_like(T, 500.0)


class SS316(SolidMaterial):
    """316L 不锈钢 (SS316L)"""

    def __init__(self, name: str = "SS316"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        """密度 [kg/m^3]"""
        T = np.asarray(T, dtype=float)
        return np.full_like(T, 7900.0)

    def conductivity(self, T: Numeric) -> Numeric:
        """导热系数 [W/(m*K)]"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.0
        return 15.61342 + 0.01324 * TC

    def heat_capacity(self, T: Numeric) -> Numeric:
        """定压比热容 [J/(kg*K)]"""
        T = np.asarray(T, dtype=float)
        TC = T - 273.0
        return 509.77978 + 0.14008 * TC


class SS316H(SolidMaterial):
    """316H 不锈钢 (SS316H)"""

    def __init__(self, name: str = "SS316H"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        """密度 [kg/m^3]"""
        T = np.asarray(T, dtype=float)
        return np.full_like(T, 7900.0)

    def conductivity(self, T: Numeric) -> Numeric:
        """导热系数 [W/(m*K)]"""
        T = np.asarray(T, dtype=float)
        return 9.2 + 0.0175 * T - 0.000002 * T ** 2

    def heat_capacity(self, T: Numeric) -> Numeric:
        """定压比热容 [J/(kg*K)]"""
        T = np.asarray(T, dtype=float)
        # 防除零保护
        T_safe = np.maximum(T, 1e-3)
        return 472.0 + 0.136 * T - 2820000.0 / (T_safe ** 2)


class Haynes(SolidMaterial):
    """Haynes 合金"""

    def __init__(self, name: str = "Haynes"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        """密度 [kg/m^3]"""
        T = np.asarray(T, dtype=float)
        return np.full_like(T, 8180.0)

    def conductivity(self, T: Numeric) -> Numeric:
        """导热系数 [W/(m*K)]"""
        T = np.asarray(T, dtype=float)

        cond_low = 0.017 * (T - 273.15) + 9.033
        cond_high = 6e-8 * T ** 3 - 0.0002 * T ** 2 + 0.2423 * T - 73.001

        return np.where(T < 773.15, cond_low, cond_high)

    def heat_capacity(self, T: Numeric) -> Numeric:
        """定压比热容 [J/(kg*K)]"""
        T = np.asarray(T, dtype=float)

        capa_low = 0.145 * (T - 273.15) + 453.7
        capa_high = 0.112 * T + 443.64

        return np.where(T < 773.15, capa_low, capa_high)
