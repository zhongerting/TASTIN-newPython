from abc import ABC, abstractmethod
from typing import Union
import numpy as np

# 定义通用的数值类型：可以是单个浮点数，也可以是 numpy 数组
Numeric = Union[float, np.ndarray]


# 物性的最底层抽象基类
class Material(ABC):
    def __init__(self, name: str, formula: str = ""):
        self.name = name
        self.formula = formula

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.name}>"


# 固体材料基类
class SolidMaterial(Material):
    def __init__(self, name: str, formula: str = ""):
        super().__init__(name, formula)
        self.molar_mass = None

    @abstractmethod
    def conductivity(self, T: Numeric) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: MaterCond
        :param T: 温度 [K] (标量或数组)
        """
        pass

    @abstractmethod
    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        定压比热容 Cp [J/(kg·K)]
        对应 Fortran: MaterCp
        :param T: 温度 [K] (标量或数组)
        """
        pass

    @abstractmethod
    def density(self, T: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: MaterDens
        :param T: 温度 [K] (标量或数组)
        """
        pass

    def emissivity(self, T: Numeric) -> Numeric:
        """
        表面发射率 epsilon [-] (可选，默认返回 0.8)
        对应 Fortran: MaterEmiss
        """
        # 如果 T 是数组，返回同维度的 0.8；如果是标量，返回 0.8
        if np.ndim(T) > 0:
            return np.full_like(T, 0.8, dtype=float)
        return 0.8

    def diffusivity(self, T: Numeric) -> Numeric:
        """
        热扩散系数 alpha = k / (rho * cp) [m^2/s]
        (派生属性，通常不需要重写)
        """
        return self.conductivity(T) / (self.density(T) * self.heat_capacity(T))

    def saturation_pressure(self, T):
        pass

    def vapor_viscosity(self, T):
        pass

    def latent_heat(self, T):
        pass


# 流体材料基类
class FluidMaterial(Material):
    @abstractmethod
    def density(self, T: Numeric, P: Numeric) -> Numeric:
        """
        密度 rho [kg/m^3]
        对应 Fortran: SDENPHSO (钠), PDENTPO (钾)
        """
        pass

    @abstractmethod
    def viscosity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        动力粘度 mu [Pa·s]
        对应 Fortran: VISCLSO
        """
        pass

    @abstractmethod
    def heat_capacity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        定压比热容 Cp [J/(kg·K)]
        对应 Fortran: CPSATLTPO
        """
        pass

    @abstractmethod
    def conductivity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        导热系数 k [W/(m·K)]
        对应 Fortran: THCONLSO
        """
        pass

    @abstractmethod
    def enthalpy(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        比焓 h [J/kg]
        对应 Fortran: SFHPSO
        """
        pass

    @abstractmethod
    def temperature_from_enthalpy(self, h: Numeric, P: Numeric) -> Numeric:
        """
        根据焓和压力反求温度 T = f(h, P)
        用于能量方程求解 T
        对应 Fortran: TPHSO
        """
        pass

    def prandtl_number(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        普朗特数 Pr = mu * Cp / k [-]
        对应 Fortran: PRFSO
        """
        mu = self.viscosity(T, P)
        cp = self.heat_capacity(T, P)
        k = self.conductivity(T, P)
        return mu * cp / k

        # ==========================================================================
        # 2. 两相流接口 (Two-Phase Properties)
        #    注：若子类是不可压缩单相流体，可实现为抛出异常或返回固定值
        # ==========================================================================

    @abstractmethod
    def saturation_pressure(self, T: Numeric) -> Numeric:
        """饱和蒸汽压 P_sat [Pa]"""
        pass

    @abstractmethod
    def saturation_temperature(self, P: Numeric) -> Numeric:
        """饱和温度 T_sat [K]"""
        pass

    @abstractmethod
    def latent_heat(self, T: Numeric) -> Numeric:
        """汽化潜热 h_fg [J/kg]"""
        pass

    # ==========================================================================
    # 3. 兼容性层 (Compatibility Layer)
    #    用于桥接 Correlations.py 中的调用与 Sodium.py 等具体实现
    # ==========================================================================

    def liquid_density(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """[兼容接口] 获取液相密度"""
        # 优先尝试调用具体的饱和液比体积接口 (如 Sodium.py 中定义)
        if hasattr(self, 'specific_volume_liquid_sat'):
            vol = getattr(self, 'specific_volume_liquid_sat')(T)
            # 防止除零 (vol 可能是数组)
            return 1.0 / np.maximum(vol, 1e-10)
        # 回退到通用 density
        return self.density(T, P)

    def vapor_density(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """[兼容接口] 获取气相密度"""
        if hasattr(self, 'specific_volume_vapor_sat'):
            vol = getattr(self, 'specific_volume_vapor_sat')(T)
            return 1.0 / np.maximum(vol, 1e-10)
        raise NotImplementedError(f"{self.name} 未定义气相密度接口 (specific_volume_vapor_sat)")

    def liquid_viscosity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """[兼容接口] 获取液相粘度"""
        if hasattr(self, 'viscosity_liquid'):
            return getattr(self, 'viscosity_liquid')(T)
        # 回退到通用 viscosity
        return self.viscosity(T, P)

    def vapor_viscosity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """[兼容接口] 获取气相粘度"""
        if hasattr(self, 'viscosity_vapor'):
            return getattr(self, 'viscosity_vapor')(T)
        raise NotImplementedError(f"{self.name} 未定义气相粘度接口 (viscosity_vapor)")

    def liquid_enthalpy_sat(self, T: Numeric) -> Numeric:
        """[兼容接口] 饱和液焓"""
        if hasattr(self, 'enthalpy_saturated_liquid'):
            return getattr(self, 'enthalpy_saturated_liquid')(T)
        # 尝试通过通用 enthalpy 估算 (假设 P 为该温度下的饱和压)
        P_sat = self.saturation_pressure(T)
        return self.enthalpy(T, P_sat)

    def vapor_enthalpy_sat(self, T: Numeric) -> Numeric:
        """[兼容接口] 饱和气焓"""
        if hasattr(self, 'enthalpy_saturated_vapor'):
            return getattr(self, 'enthalpy_saturated_vapor')(T)
        # 尝试通过液焓 + 潜热
        return self.liquid_enthalpy_sat(T) + self.latent_heat(T)

    # ==========================================================================
    # 4. 派生属性与导数
    # ==========================================================================

    @abstractmethod
    def liquid_density_derivative_P(self, P: Numeric) -> Numeric:
        """密度对压力的导数 drho_dp"""
        pass

    @abstractmethod
    def liquid_density_derivative_T(self, T: Numeric) -> Numeric:
        """密度对温度的导数 drho_dT"""
        pass


# 为了保持向后兼容（如果之前的代码引用了 TwoPhaseFluidMaterial），可以定义别名
TwoPhaseFluidMaterial = FluidMaterial
