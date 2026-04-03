from Materials.Base import SolidMaterial


class WickStructure(SolidMaterial):
    """
    吸液芯骨架 (Wick Structure) 物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 这是一个等效材料模型，数值设定具有一定的奇异性（可能是占位符）。
    - 密度: 1.0 (原代码值)
    - 导热: 0.15
    - 比热: 100.0 (原代码值，未除以1000，但在 HPTsCal 中调用时乘了 1000)
    """

    def __init__(self):
        # 假设吸液芯材质，这里原代码未指明，通常为丝网或烧结粉末
        super().__init__(name="Wick", formula="Composite")

    def conductivity(self, T: float) -> float:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCWICK
        公式: k = 0.15
        """
        return 0.15

    def density(self, T: float) -> float:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENWICK
        公式: rho = 1.0
        注意: 这个值极小，可能是作为孔隙率修正的基准或相对值，
        但在重构阶段我们必须忠实复现。
        """
        return 1.0

    def heat_capacity(self, T: float) -> float:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPWICK

        注意: 原代码返回 100.0。
        在 HPTsCal.for 中，该值被乘以 1.0e3 使用。
        为了配合标准接口，我们这里直接返回 100.0 * 1000.0 = 100000.0
        """
        # 原代码: CPWICK=100.0
        return 100.0
