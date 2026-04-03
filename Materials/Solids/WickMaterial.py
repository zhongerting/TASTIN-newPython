from abc import ABC

import numpy as np
from typing import Union
from Materials.Base import SolidMaterial

# 定义类型别名以支持向量化提示
Numeric = Union[float, np.ndarray]


class WickMaterial(SolidMaterial):
    """
    热管吸液芯复合材料 (Heat Pipe Wick Composite Material)

    物理特性:
    1. 由固体骨架 (如不锈钢) 和工作工质 (如液态金属钠) 组成。
    2. 等效热容基于孔隙率 (Porosity) 按照体积热容加权计算。
    3. 导热系数由“结构等效导热系数”和由气液相变引起的“赝导热系数”叠加而成。
    """

    def __init__(self,
                 name: str,
                 solid_mat: SolidMaterial,
                 fluid_mat: SolidMaterial,  # 工质现已作为 SolidMaterial 的衍生类
                 porosity: float,
                 r_vapor: float,
                 r_in_wall: float):
        """
        初始化吸液芯材料

        :param name: 材料名称描述
        :param solid_mat: 固体骨架结构材料 (例如 SS316 实例)
        :param fluid_mat: 工作工质材料对象 (例如 SodiumHP 实例)
        :param porosity: 孔隙率 PHI [-]
        :param r_vapor: 蒸汽腔外半径 Rv [m] (用于赝热导率计算)
        :param r_in_wall: 管壁内半径 Rw [m]
        """
        super().__init__(name=name)

        # 1. 保存物理对象与几何参数
        self.solid = solid_mat
        self.fluid = fluid_mat
        self.porosity = porosity
        self.r_vapor = r_vapor
        self.r_in_wall = r_in_wall

        # 2. 物理常数缓存
        self.R_gas = 8.314  # 理想气体常数 R [J/(mol·K)]
        self.k_boltzmann = 1.380649e-23  # 玻尔兹曼常数

    # =====================================================================
    # 1. 基础物性：孔隙率加权计算 (重写 SolidMaterial 基类方法)
    # =====================================================================

    def density(self, T: Numeric) -> Numeric:
        """
        [向量化接口] 等效密度 [kg/m^3]
        按照物理与数值技巧约定：为了使得 (density * heat_capacity) 完美等于等效体积热容，
        这里将密度简单固定为工质的密度。
        """
        return self.fluid.density(T)

    def heat_capacity(self, T: Numeric) -> Numeric:
        """
        [向量化接口] 等效定压比热容 [J/(kg·K)]
        公式: capa[i] = PHI * capawf[i] + (1 - PHI) * scapa; (体积热容加权)
        """
        phi = self.porosity

        # 1. 获取流体与固体的密度和质量比热容
        rho_f = self.fluid.density(T)
        cp_f = self.fluid.heat_capacity(T)
        rho_s = self.solid.density(T)
        cp_s = self.solid.heat_capacity(T)

        # 2. 计算各自的体积热容 [J/(m^3*K)]
        capa_wf = rho_f * cp_f
        s_capa = rho_s * cp_s

        # 3. 计算等效体积热容 capa_eff
        capa_eff = phi * capa_wf + (1.0 - phi) * s_capa

        # 4. 除以 self.density(T) 的返回值 (即 rho_f)，得到供求解器使用的等效质量热容
        return capa_eff / rho_f

    def conductivity(self, T: Numeric) -> Numeric:
        """
        [向量化接口] 总等效导热系数 [W/(m·K)]
        最终的导热系数 = 等效的材料导热系数 + 附加赝导热系数
        """
        phi = self.porosity

        # 1. 获取本征导热系数
        k_f = self.fluid.conductivity(T)
        k_s = self.solid.conductivity(T)

        # 2. 计算结构等效导热系数 (k_self)
        with np.errstate(divide='ignore', invalid='ignore'):
            numerator = k_s + k_f - (1.0 - phi) * (k_f - k_s)
            denominator = k_f + k_s + (1.0 - phi) * (k_f - k_s)
            k_self = k_f * (numerator / denominator)
            k_self = np.nan_to_num(k_self, nan=0.0)

        # 3. 获取相变流体特有物性 (通过扩展接口调用)
        P_sat = self.fluid.saturation_pressure(T)
        mu_v = self.fluid.vapor_viscosity(T)
        h_fg = self.fluid.latent_heat(T)

        # 4. 计算附加赝导热系数
        pse1 = self._pse1_conductivity(T, P_sat, mu_v, h_fg)

        # 采用我们之前探讨的物理修正项：总的相变导热能力应当受到分子稀薄效应(Kn)的阻力折减
        T_safe = np.maximum(T, 1e-3)
        P_safe = np.maximum(P_sat, 1e-12)
        vAve = np.sqrt((8.0 * self.R_gas * T_safe) / (np.pi * self.fluid.molar_mass))

        with np.errstate(divide='ignore', invalid='ignore'):
            parameter1 = (mu_v * vAve) / (self.r_vapor * P_safe)
            parameter2 = 1.0 + (8.0 / 3.0) * parameter1

            # 【物理修正】：相变总导热 = 连续流导热(pse1) / 阻力修正(parameter2)
            k_pse_total = pse1 / parameter2
            k_pse_total = np.nan_to_num(k_pse_total, nan=0.0)

        # 5. 总导热系数 = 结构等效导热 + 相变赝导热
        k_total = k_self + k_pse_total

        # 限制上限防刚性溢出 (对应原 C++ 的 > 1e6 截断)
        return np.clip(k_total, a_min=0.0, a_max=1e6)

    # =====================================================================
    # 2. 内部微观物理方法：赝热导率计算
    # =====================================================================

    def _pse1_conductivity(self, T: Numeric, P_sat: Numeric, mu_v: Numeric, h_fg: Numeric) -> Numeric:
        """
        [内部向量化方法] 连续流体制下的一阶赝热导率 (Pseudo-thermal conductivity order 1)
        """
        T_safe = np.maximum(T, 1e-3)
        mu_safe = np.maximum(mu_v, 1e-12)

        Rv = self.r_vapor
        Rw = self.r_in_wall
        Mg = self.fluid.molar_mass
        R_gas = self.R_gas

        # 参数 1: 几何项
        parameter1 = (Rv ** 4) / (Rw ** 2 - Rv ** 2)

        with np.errstate(divide='ignore', invalid='ignore'):
            # 参数 2: 热力学驱动项
            parameter2 = (h_fg * Mg * P_sat) ** 2

            # 参数 3: 阻力项
            parameter3 = 4.0 * mu_safe * (R_gas ** 2) * (T_safe ** 3)

            # 计算并清理异常值
            pse1CondTemp = parameter1 * parameter2 / parameter3
            pse1CondTemp = np.nan_to_num(pse1CondTemp, nan=0.0, posinf=1e6)

        return pse1CondTemp
