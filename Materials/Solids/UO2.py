import numpy as np
from Materials.Base import SolidMaterial, Numeric


class UO2(SolidMaterial):
    """
    二氧化铀 (Uranium Dioxide) 燃料物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度设定为常数 (10400 kg/m^3)
    - 导热和比热容包含针对特定温度段的修正逻辑(原代码中对低温段进行了截断处理)
    """

    def __init__(self):
        super().__init__(name="UO2", formula="UO2")
        # 经验公式中的理论密度因子，原代码中硬编码为 0.95
        self.td_fraction = 0.95

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCUO2

        逻辑复现:
        - T > 1500K: 使用 Harding & Martin 形式的经验公式
        - T <= 1500K: 锁定为 1500K 时的导热系数 (原代码逻辑)
        """
        T = np.asarray(T)

        # [向量化] 使用 np.maximum 实现截断逻辑: max(T, 1500.0)
        T_calc = np.maximum(T, 1500.0)

        # 摄氏温度 t
        t_c = T_calc - 273.15

        # 常用项预计算
        term1_denom = t_c + 2360.0 - 2120.0 * (self.td_fraction ** 3.0)
        term1 = 4575.0 / term1_denom

        term2_coeff = 5.67e-10 * (self.td_fraction ** 4.0)
        term2_base = t_c - 860.0
        term2 = term2_coeff * (term2_base ** 3.0)

        k = term1 + term2
        return k

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENUO2

        逻辑复现: 返回常数 10400.0
        """
        return 10400.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPUO2

        逻辑复现:
        - T < 700K: 锁定为 700K 时的比热容
        - T >= 700K: 使用多项式公式
        """
        T = np.asarray(T)

        # [向量化] 使用 np.maximum 实现截断逻辑: max(T, 700.0)
        T_calc = np.maximum(T, 700.0)

        t_c = T_calc - 273.15

        # 公式: Cp = 297 + 0.025*t - 6.12e6 / t^2
        term1 = 297.0
        term2 = 2.5e-2 * t_c
        term3 = 6.12e6 / (t_c ** 2.0)

        cp = term1 + term2 - term3
        return cp
