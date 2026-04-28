import numpy as np
from typing import Union

from Materials.Base import SolidMaterial

Numeric = Union[float, np.ndarray]


def _clip_temperature_input(T: Numeric) -> np.ndarray:
    return np.maximum(np.asarray(T, dtype=float), 1.0e-3)


class SS321(SolidMaterial):
    """321 stainless steel."""

    def __init__(self, name: str = "SS321"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        return np.full_like(T, 8090.0)

    def conductivity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        temp_c = T - 273.0
        return 14.5 + 0.015 * temp_c

    def heat_capacity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        return np.full_like(T, 500.0)


class SS316(SolidMaterial):
    """316L stainless steel."""

    def __init__(self, name: str = "SS316"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        return np.full_like(T, 7900.0)

    def conductivity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        temp_c = T - 273.0
        conductivity = 15.61342 + 0.01324 * temp_c
        return np.maximum(conductivity, 1.0e-6)

    def heat_capacity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        temp_c = T - 273.0
        return 509.77978 + 0.14008 * temp_c


class SS316H(SolidMaterial):
    """316H stainless steel."""

    def __init__(self, name: str = "SS316H"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        return np.full_like(T, 7900.0)

    def conductivity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        conductivity = 9.2 + 0.0175 * T - 0.000002 * T ** 2
        return np.maximum(conductivity, 1.0e-6)

    def heat_capacity(self, T: Numeric) -> Numeric:
        T_safe = _clip_temperature_input(T)
        return 472.0 + 0.136 * T_safe - 2820000.0 / (T_safe ** 2)


class Haynes(SolidMaterial):
    """Haynes alloy."""

    def __init__(self, name: str = "Haynes"):
        super().__init__(name=name)

    def density(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        return np.full_like(T, 8180.0)

    def conductivity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        cond_low = 0.017 * (T - 273.15) + 9.033
        cond_high = 6e-8 * T ** 3 - 0.0002 * T ** 2 + 0.2423 * T - 73.001
        return np.where(T < 773.15, cond_low, cond_high)

    def heat_capacity(self, T: Numeric) -> Numeric:
        T = _clip_temperature_input(T)
        capa_low = 0.145 * (T - 273.15) + 453.7
        capa_high = 0.112 * T + 443.64
        return np.where(T < 773.15, capa_low, capa_high)
