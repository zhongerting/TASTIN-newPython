import numpy as np
from Materials.Base import TwoPhaseFluidMaterial, Numeric
from scipy import optimize

class NaK78(TwoPhaseFluidMaterial):
    """
    NaK78 合金流体物性模型

    来源: 高温热管专著V4
    适用范围: 265.0 K - 2503.7 K (临界温度)
    修改状态: 已向量化 (Vectorized for NumPy)
    """
    def __init__(self):
        super().__init__(name="Sodium", formula="Na")
        self.T_crit = 2503.7  # 临界温度 [K]
        self.P_crit = 25.64e6  # 临界压力

    def density(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """密度 rho [kg/m^3]"""
        v = self.specific_volume(T, P)
        return 1.0 / np.maximum(v, 1e-10)



