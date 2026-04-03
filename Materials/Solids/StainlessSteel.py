import numpy as np
from Materials.Base import SolidMaterial, Numeric


class AusteniticStainlessSteel(SolidMaterial):
    """
    奥氏体不锈钢 (Austenitic Stainless Steel) 物性模型

    来源: PROPSODIUM.FOR
    适用范围: 0 - 1700 K (根据原代码注释)
    """

    def __init__(self):
        super().__init__(name="Austenitic_SS", formula="Fe-Cr-Ni")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCAUSTENTIC
        公式: k = 9.2 + 0.0175*T - 2.0e-6*T^2
        """
        T = np.asarray(T)
        # 多项式，天然支持向量化
        k = 9.2 + 0.0175 * T - 2.0e-6 * T ** 2
        return k

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENAUSTENTIC
        公式: rho = 8084 - 0.4209*T - 3.894e-5*T^2
        """
        T = np.asarray(T)
        rho = 8084.0 - 0.4209 * T - 3.894e-5 * T ** 2
        return rho

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPAUSTENTIC

        注意：原 Fortran 代码返回 kJ/(kg·K)，此处已转换为 J/(kg·K)。
        公式: Cp = 472 + 0.136*T - 2.82e6/T^2
        """
        T = np.asarray(T)

        # [向量化修复关键点]
        # 1. 构造安全温度数组，避免 T=0 时除以零报错 (仅用于计算公式中间值)
        #    对于 T<=0 的位置，暂时用 1.0 代替
        T_safe = np.where(T <= 0, 1.0, T)

        # 2. 计算公式值
        cp_calc = 472.0 + 0.136 * T_safe - 2.82e6 / (T_safe ** 2)

        # 3. 使用 np.where 根据真实温度 T 组合最终结果
        #    如果 T <= 0，返回 472.0；否则返回计算值
        return np.where(T <= 0, 472.0, cp_calc)

    def thermal_expansion_coeff(self, T: Numeric) -> Numeric:
        """
        线膨胀系数 alpha [1/K]
        对应 Fortran: EXPAUSTENTIC
        公式: alpha = 1.789e-5 + 2.398e-9*T + 3.269e-13*T^2
        """
        T = np.asarray(T)
        alpha = 1.789e-5 + 2.398e-9 * T + 3.269e-13 * T ** 2
        return alpha
