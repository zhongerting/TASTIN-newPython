import numpy as np
from Materials.Base import SolidMaterial, Numeric


class ZirconiumHydride(SolidMaterial):
    """
    慢化剂材料 (Zirconium Hydride, ZrH) 物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 导热系数为常数 (20.0 W/mK)
    - 密度为常数 (5615.0 kg/m^3)
    - 比热容为温度的线性函数
    """

    def __init__(self):
        super().__init__(name="ZrH", formula="ZrH")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCZrH
        公式: k = 20.0
        """
        # [向量化关键] 
        # 虽然是常数，但为了支持 HeatConduction 中的切片操作 (k[:-1])，
        # 必须返回与 T 维度一致的数组。
        # 使用 0.0 * T + 20.0 是一种简便的广播方法。
        return 20.0 + 0.0 * T

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENZrH
        公式: rho = 5615.0
        """
        # 密度通常用于乘法运算 (rho * Cp)，标量可以自动广播，因此返回常数是安全的。
        return 5615.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPZrH

        注意: 原代码中除以 1000.0，说明原单位为 kJ/(kg·K)。
        此处移除除法，直接返回 SI 单位 J/(kg·K)。
        公式: Cp = 187.0 + 0.745 * T
        """
        # 线性公式，天然支持向量化
        # (187.0 + 0.745*T) J/kgK
        cp = 187.0 + 0.745 * T
        return cp
