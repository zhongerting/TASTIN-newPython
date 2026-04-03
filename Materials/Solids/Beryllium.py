import numpy as np
from Materials.Base import SolidMaterial, Numeric


class Beryllium(SolidMaterial):
    """
    金属铍 (Be) 反射层材料物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度为常数 (1830.0 kg/m^3)
    - 导热系数和比热容均为分段函数 (低温/高温段取边界定值，中间段为多项式)
    """

    def __init__(self):
        super().__init__(name="Be", formula="Be")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCBe
        """
        # 确保输入为数组
        T = np.asarray(T)

        # 逻辑复现:
        # T <= 198.1: 取 T=198.1
        # T >= 1556.0: 取 T=1556.0
        # 中间段: 保持原值
        # 使用 np.clip 一步完成截断
        t = np.clip(T, 198.1, 1556.0)

        # 公式: 304.362 - 0.472*T + 3.38e-4*T^2 - 9.0729e-8*T^3
        k = 304.362 - 0.472 * t + 3.38035e-4 * (t ** 2.0) - 9.0729e-8 * (t ** 3.0)
        return k

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENBe
        公式: rho = 1830.0
        """
        return 1830.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPBe
        """
        T = np.asarray(T)

        # 1. 计算中间段的多项式值
        cp_poly = 256.23 + 6.58 * T - 5.262e-3 * (T ** 2.0) + 1.54667e-6 * (T ** 3.0)

        # 2. 定义分段条件和对应值
        conditions = [T <= 223.1, T >= 1556.0]
        choices = [142.3, 357.7]

        # 3. 使用 np.select 组合结果 (若都不满足则使用 default=cp_poly)
        cp = np.select(conditions, choices, default=cp_poly)

        return cp
