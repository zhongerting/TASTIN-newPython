import numpy as np
from Materials.Base import SolidMaterial, Numeric


class BoronCarbide(SolidMaterial):
    """
    碳化硼 (Boron Carbide, B4C) 控制棒材料物性模型

    来源: PROPSODIUM.FOR
    物理特征:
    - 密度设定为常数 (2300 kg/m^3)
    - 导热系数和比热容采用查表线性插值计算
    """

    def __init__(self):
        super().__init__(name="B4C", formula="B4C")

        # ==========================================
        # 导热系数查找表 (来源: TCB4C 函数)
        # ==========================================
        self._k_temp = np.array([
            400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
            750.0, 800.0, 850.0, 900.0, 950.0, 1000.0
        ])
        self._k_vals = np.array([
            16.10, 14.67, 13.69, 13.62, 12.81, 10.91, 12.17,
            12.32, 11.48, 9.196, 9.903, 10.21, 8.434
        ])

        # ==========================================
        # 比热容查找表 (来源: CPB4C 函数)
        # ==========================================
        self._cp_temp = np.array([
            300.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0,
            1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0,
            1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0
        ])
        # 原值为 kJ/kgK，此处直接存为 J/kgK (*1000)
        self._cp_vals = np.array([
            0.959, 1.206, 1.341, 1.468, 1.589, 1.704, 1.812, 1.913,
            2.008, 2.096, 2.178, 2.253, 2.321, 2.383, 2.438, 2.487,
            2.529, 2.564, 2.593, 2.615, 2.631, 2.640
        ]) * 1000.0

    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: TCB4C
        """
        # np.interp 自动支持标量和数组输入，无需额外处理
        return np.interp(T, self._k_temp, self._k_vals)

    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: DENB4C
        """
        # 返回常数。虽然是标量，但在 HeatConduction 的 numpy 运算中会被自动广播，符合 Numeric 定义。
        return 2300.0

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        比热容 Cp [J/(kg·K)]
        对应 Fortran: CPB4C
        """
        # np.interp 自动支持标量和数组输入
        return np.interp(T, self._cp_temp, self._cp_vals)
