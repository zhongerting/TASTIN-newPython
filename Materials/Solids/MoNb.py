from Materials.Base import SolidMaterial, Numeric


class MoNb(SolidMaterial):
    """
    发射极材料 (Mo-Nb 合金) 物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度为常数 (11506 kg/m^3)
    - 导热系数和比热容为温度的线性函数
    """

    def __init__(self):
        super().__init__(name="MoNb", formula="Mo-Nb")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCMoNb
        公式: k = 120.0 - 0.014 * T
        """
        # 线性公式，天然支持向量化
        k = 120.0 - 0.014 * T
        return k

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENMoNb
        公式: rho = 11506.0
        """
        return 11506.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPMoNb

        注意: 原代码中数值量级为 1e-3，根据上下文推断单位为 kJ/(kg·K)。
        此处乘以 1000 转换为 SI 单位 J/(kg·K)。
        公式: Cp = 65.05 + 0.0133 * T
        """
        # 线性公式，天然支持向量化
        # (0.06505 + 0.0000133*T) kJ/kgK * 1000 -> J/kgK
        cp = (65.05e-3 + 0.0133e-3 * T) * 1000.0
        return cp
