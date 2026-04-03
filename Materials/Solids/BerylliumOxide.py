import numpy as np
from Materials.Base import SolidMaterial, Numeric


class BerylliumOxide(SolidMaterial):
    """
    氧化铍 (BeO) 反射层材料物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度为常数 (2800.0 kg/m^3)
    - 导热系数与温度平方成反比
    - 比热容为线性关系
    """

    def __init__(self):
        super().__init__(name="BeO", formula="BeO")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCBeO
        公式: k = 4.0e7 / T^2
        """
        T = np.asarray(T)

        # 向量化处理防除零逻辑
        # 1. 创建安全计算的温度数组 (避免底数为0报错)
        #    如果 T <= 0, 暂时用 1.0 代替进行计算
        T_safe = np.where(T <= 0, 1.0, T)

        # 2. 应用公式
        k_calc = 4.0e7 / (T_safe ** 2.0)

        # 3. 还原特殊情况 (如果 T <= 0, 返回 4.0e7)
        return np.where(T <= 0, 4.0e7, k_calc)

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENBeO
        公式: rho = 2800.0
        """
        return 2800.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPBeO

        注意: 原代码 (1343.0 + 1.8469*T)/1000.0 单位为 kJ/kgK。
        此处转换为 J/kgK。
        """
        T = np.asarray(T)
        # 线性公式，天然支持向量化
        return 1343.0 + 1.8469 * T
