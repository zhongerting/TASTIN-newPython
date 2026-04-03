import numpy as np
from Materials.Base import SolidMaterial, Numeric


class CarbonDioxide(SolidMaterial):
    """
    二氧化碳 (CO2) 气隙材料模型

    来源: PROPSODIUM.FOR
    应用: 通常用于外部绝缘气隙
    物理特征:
    - 密度: 1.1 kg/m^3 (常数)
    - 导热: 随温度线性变化
    - 比热: 100.0 * 1000 (原代码值)
    """

    def __init__(self):
        super().__init__(name="CO2", formula="CO2")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCCO2
        公式: k = 0.0064 + 14.4e-5 * T
        """
        # 线性公式，自动支持向量化
        return 0.0064 + 1.44e-4 * T

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENCO2
        """
        return 1.1

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPCO2
        """
        return 0.88 * 1000.0


class Helium(SolidMaterial):
    """
    氦气 (Helium, He) 气隙材料模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度: 4.0/22.4 ≈ 0.178 kg/m^3 (理想气体估算)
    - 导热: 随温度线性变化
    - 比热: 10000.0 * 1000 (极大的值，可能作为数值稳定处理)
    """

    def __init__(self):
        super().__init__(name="Helium", formula="He")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCHE
        公式: k = 0.00031 * T + 0.059
        """
        return 3.1e-4 * T + 0.059

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENHE
        """
        return 4.0 / 22.4

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: HECP
        """
        return 5.19 * 1000.0


class Xenon(SolidMaterial):
    """
    氙气 (Xenon, Xe) / 裂变气体模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 导热系数固定为 0.3 W/mK (注: 真实 Xe 约为 0.005，这可能是包含辐射的等效值)
    - 密度: 1.0
    - 比热: 100.0 * 1000
    """

    def __init__(self):
        super().__init__(name="FissionGas_Xe", formula="Xe")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCXe
        """
        return 0.3

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENXe
        """
        return 1.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPXe
        """
        return 0.158*1000


class Cesium(SolidMaterial):
    """
    铯蒸汽 (Cesium Vapor, Cs) 间隙模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 导热系数固定为 0.15 W/mK
    - 密度: 1.0
    - 比热: 100.0 * 1000
    """

    def __init__(self):
        super().__init__(name="CesiumVapor", formula="Cs")

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCCs
        """
        return 0.15

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENCs
        """
        return 1.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPCs
        """
        return 100.0
