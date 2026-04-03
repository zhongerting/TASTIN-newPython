from Materials.Base import SolidMaterial, Numeric


class Molybdenum(SolidMaterial):
    """
    接收极材料 (Molybdenum, Mo) 物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度为常数 (10200 kg/m^3)
    - 导热系数和比热容为温度的线性函数
    """

    def __init__(self):
        super().__init__(name="Mo", formula="Mo")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCMo
        公式: k = (110.0 - 0.015 * T) * 0.9
        """
        # 原代码: TCMo=(110.0-0.015*TT)*0.9
        # 线性公式，天然支持向量化
        k = (110.0 - 0.015 * T) * 0.9
        return k

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENMo
        公式: rho = 10200.0
        """
        return 10200.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPMo

        注意: 原代码中数值量级为 e-3 (如 101.8e-3)，单位推断为 kJ/(kg·K)。
        此处乘以 1000 转换为 SI 单位 J/(kg·K)。
        公式: Cp = 101.8 + 0.025 * T
        """
        # 线性公式，天然支持向量化
        # 转换: (0.1018 + 0.000025*T) kJ/kgK * 1000 -> 101.8 + 0.025*T J/kgK
        cp = (101.8e-3 + 0.025e-3 * T) * 1000.0
        return cp

