import numpy as np
from typing import Union
from Materials.Base import SolidMaterial

# 定义类型别名以支持向量化提示
Numeric = Union[float, np.ndarray]


class PotassiumHP(SolidMaterial):
    """
    热管液态金属工质：钾 (Potassium, K)
    引入等效热容法 (Apparent Heat Capacity Method) 处理固液相变 (糊状区 Mushy Zone)
    """

    # NOTE(KHP): molar_mass, T_melt, T_crit, H_sf, and dt_thaw have
    # been updated from user-provided values on 2026-06-04/2026-06-05. For derived
    # property formulas below, each formula line is annotated with its
    # source/range when user-verified, or marked as an unverified
    # pre-existing repository formula.
    def __init__(self, name: str = "Potassium_Wick_Fluid"):
        super().__init__(name=name)
        # 基础物理常数
        # User-provided value refined on 2026-06-05 from NIST molecular weight.
        self.molar_mass = 3.90983e-2  # 摩尔质量 [kg/mol]
        self.universal_gas_constant = 8.314462618  # [J/(mol*K)]
        self.vapor_specific_gas_constant = self.universal_gas_constant / self.molar_mass
        self.vapor_gamma = 5.0 / 3.0

        # 相变(融化)参数
        # User-provided values, 2026-06-04. Reference/range not yet supplied.
        self.T_melt = 336.35  # 融化温度 [K]
        self.H_sf_molar_cal = 553.8  # 熔化潜热 [cal/mol]
        self.H_sf = self.H_sf_molar_cal * 4.184 / self.molar_mass  # 熔化潜热 [J/kg]
        self.dt_thaw = 1.0  # 糊状区半宽 [K] (总宽度为 2K)

        self.T_crit = 2287.0  # 临界温度 [K]

    def _get_mushy_fraction(self, T: np.ndarray) -> np.ndarray:
        """计算糊状区插值因子 f, 范围 [0, 1]"""
        f = (T - self.T_melt + self.dt_thaw) / (2.0 * self.dt_thaw)
        return np.clip(f, 0.0, 1.0)

    @staticmethod
    def _smoothstep(x: np.ndarray) -> np.ndarray:
        x = np.clip(x, 0.0, 1.0)
        return x * x * (3.0 - 2.0 * x)

    def density(self, T: Numeric) -> Numeric:
        # Solid potassium density formula verified from user-provided data.
        # Source validity note: recommended for 270-320 K; using it near
        # the mushy interval is an extrapolation until a melt-range formula is provided.
        # Liquid potassium density formula verified from user-provided data.
        # Source validity note: recommended for 623-1030 K, approximately extendable
        # to 1123 K if the system pressure maintains liquid or saturated-liquid behavior.
        """密度 [kg/m^3]"""
        T = np.asarray(T, dtype=float)
        # Solid K density [kg/m3].
        # Formula: rho_s = 857.6 / (1 + 2.39e-4*(T - 300)).
        # Range: 270-320 K; normal pressure, near-room-temperature solid K.
        # Reference: Schouten, D. R. and Swenson, C. A.,
        # "Linear-thermal-expansion measurements on potassium metal from
        # 2 to 320 K," Physical Review B, Vol. 10, No. 6, pp. 2175-2185, 1974.
        rho_solid = 857.6 / (1.0 + 2.39e-4 * (T - 300.0))
        # Saturated/liquid K density [kg/m3].
        # Formula: rho_l_sat = 890.29 - 2.113e-1*T.
        # Range: 623-1123 K; above the normal boiling point this should be
        # interpreted as saturated liquid at the corresponding saturation pressure.
        # Reference: Engineering ToolBox, "Potassium -- Thermophysical
        # Properties vs. Temperature," saturated liquid density data from
        # 350 C to 850 C; Hoffman and Cox, ORNL-TM-2126/NASA-CR-96741, 1968.
        rho_liquid = 890.29 - 2.113e-1 * T

        # Mushy-zone interpolation, apparent heat-capacity model support.
        # Range: [T_melt-dt_thaw, T_melt+dt_thaw]. This interpolation is a
        # repository model assumption, not a separately referenced property formula.
        f = self._get_mushy_fraction(T)
        rho_mushy = rho_solid + f * (rho_liquid - rho_solid)

        return np.where(T < self.T_melt - self.dt_thaw, rho_solid,
                        np.where(T > self.T_melt + self.dt_thaw, rho_liquid, rho_mushy))

    def conductivity(self, T: Numeric) -> Numeric:
        # Solid potassium conductivity formula verified from user-provided data.
        # Source validity note: recommended for 270-336.8 K.
        # Solid and liquid potassium conductivity formulas verified from user-provided data.
        """导热系数 [W/(m*K)]"""
        T = np.asarray(T, dtype=float)

        # Solid K thermal conductivity [W/(m*K)].
        # Formula: lambda_s = 102.5 - 1.04e-1*(T - 298.2).
        # Range: 270-336.8 K; normal pressure solid K from room temperature
        # to near melting. Do not extrapolate to the low-temperature nonlinear zone.
        # Reference: Ho, C. Y., Powell, R. W., and Liley, P. E.,
        # "Thermal Conductivity of the Elements," Journal of Physical and
        # Chemical Reference Data, Vol. 1, No. 2, pp. 279-421, 1972.
        cond_solid = 102.5 - 1.04e-1 * (T - 298.2)

        # Liquid K thermal conductivity [W/(m*K)].
        # Formula: lambda_l = 66.09 - 3.579e-2*T.
        # Range: 336.8-1000 K; normal pressure liquid K from near melting to
        # before the normal boiling point. Use a higher-order correlation or
        # table interpolation for high-pressure liquid K or near-critical conditions.
        # Reference: Ho, C. Y., Powell, R. W., and Liley, P. E.,
        # "Thermal Conductivity of the Elements," Journal of Physical and
        # Chemical Reference Data, Vol. 1, No. 2, pp. 279-421, 1972.
        cond_liquid = 66.09 - 3.579e-2 * T

        # Mushy-zone interpolation, repository model assumption.
        # Range: [T_melt-dt_thaw, T_melt+dt_thaw].
        f = self._get_mushy_fraction(T)
        cond_mushy = cond_solid + f * (cond_liquid - cond_solid)

        return np.where(T < self.T_melt - self.dt_thaw, cond_solid,
                        np.where(T > self.T_melt + self.dt_thaw, cond_liquid, cond_mushy))

    def heat_capacity(self, T: Numeric) -> Numeric:
        # Solid and liquid potassium heat-capacity formulas verified from user-provided data.
        """定压比热容 [J/(kg*K)]"""
        T = np.asarray(T, dtype=float)
        # Solid K constant-pressure heat capacity [J/(kg*K)].
        # Shomate temperature variable: tau = T/1000.
        # Range: 298-336.35 K; normal pressure solid K near melting.
        # Reference: Chase, M. W., Jr., NIST-JANAF Thermochemical Tables,
        # Fourth Edition, Journal of Physical and Chemical Reference Data,
        # Monograph 9, pp. 1-1951, 1998.
        tau = T / 1000.0
        capa_solid = (
            -63.47410
            - 3226.340 * tau
            + 14644.60 * tau ** 2
            - 16229.50 * tau ** 3
            + 16.29410 * tau ** -2
        ) / 3.90983e-2
        # Liquid K constant-pressure heat capacity [J/(kg*K)].
        # Shomate temperature variable: tau = T/1000.
        # Range: 336.35-1039.54 K; normal pressure liquid K from near melting
        # to near normal boiling. Do not use for solid K, near-critical states,
        # or high-pressure non-saturated liquid without further validation.
        # Reference: Chase, M. W., Jr., NIST-JANAF Thermochemical Tables,
        # Fourth Edition, Journal of Physical and Chemical Reference Data,
        # Monograph 9, pp. 1-1951, 1998.
        capa_liquid = (
            40.27113
            - 30.54542 * tau
            + 26.49505 * tau ** 2
            - 5.727854 * tau ** 3
            - 0.063477 * tau ** -2
        ) / 3.90983e-2

        # Mushy-zone interpolation plus latent heat over a 2*dt_thaw interval.
        # Range: [T_melt-dt_thaw, T_melt+dt_thaw]. This is the repository's
        # apparent heat-capacity treatment for melting.
        f = self._get_mushy_fraction(T)
        base_capa_mushy = capa_solid + f * (capa_liquid - capa_solid)
        latent_capa = self.H_sf / (2.0 * self.dt_thaw)
        capa_mushy = base_capa_mushy + latent_capa

        return np.where(T < self.T_melt - self.dt_thaw, capa_solid,
                        np.where(T > self.T_melt + self.dt_thaw, capa_liquid, capa_mushy))

    # =====================================================================
    # 赝热导率专用扩展接口
    # =====================================================================

    def saturation_pressure(self, T: Numeric) -> Numeric:
        # Potassium saturation-pressure formula verified from user-provided data.
        """饱和蒸汽压 [Pa]"""
        T_safe = np.maximum(np.asarray(T, dtype=float), 1e-3)
        # Saturation pressure [Pa].
        # Formula: ln(p_sat) = 25.109 - 10488/T - 0.448*ln(T), in Pa.
        # Range: 350-1000 K; saturated potassium vapor pressure from near
        # liquid-potassium melting to near normal boiling. This equation is for
        # the vapor-liquid coexistence line, not superheated vapor or compressed
        # subcooled liquid. Above 900 K, the pressure is total saturated pressure;
        # monomer/dimer partial pressures require additional vapor chemistry.
        # References: Bystrov, Kagan, et al., Liquid-Metal Heat Pipes:
        # Thermophysical Properties and Processes, Nauka, 1988; Fink and
        # Leibowitz, Thermophysical Properties of Alkali Metals,
        # ANL-CEN-RSD-82-4, Argonne National Laboratory, 1982; Faghri,
        # Heat Pipe Science and Technology, Second Edition, 2016.
        log_p_sat = 25.109 - 10488.0 / T_safe - 0.448 * np.log(T_safe)
        return np.exp(log_p_sat)

    def saturated_liquid_density(self, T: Numeric) -> Numeric:
        # Potassium saturated-liquid density formula verified from user-provided data.
        """Saturated liquid density [kg/m3]."""
        T = np.asarray(T, dtype=float)
        # Saturated K liquid density [kg/m3].
        # Formula: rho_l_sat = 890.29 - 2.113e-1*T.
        # Range: 623-1123 K; saturated liquid potassium. Above the normal
        # boiling point, interpret this as the saturated liquid at the
        # corresponding saturation pressure. Do not use for solid K,
        # near-critical states, or high-accuracy equation-of-state work.
        # Reference: Engineering ToolBox, "Potassium -- Thermophysical
        # Properties vs. Temperature," saturated liquid density data from
        # 350 C to 850 C; Hoffman and Cox, A Preliminary Collation of the
        # Thermodynamic and Transport Properties of Potassium,
        # ORNL-TM-2126/NASA-CR-96741, Oak Ridge National Laboratory, 1968.
        return 890.29 - 2.113e-1 * T

    def specific_volume_liquid_sat(self, T: Numeric) -> Numeric:
        """Saturated liquid specific volume [m3/kg]."""
        rho = self.saturated_liquid_density(T)
        return 1.0 / np.maximum(rho, 1.0e-30)

    def vapor_density(self, T: Numeric) -> Numeric:
        # Potassium saturated-vapor density formula verified from user-provided data.
        """Saturated vapor density [kg/m3]."""
        T_safe = np.maximum(np.asarray(T, dtype=float), 1.0e-3)
        # Saturated K vapor density [kg/m3].
        # Formula: rho_v_sat = 2.398e3*exp(-8.698e3/T).
        # Range: 623-1123 K; saturated potassium vapor. Above the normal
        # boiling point, interpret this as the saturated vapor at the
        # corresponding saturation pressure. Do not use for solid K,
        # near-critical states, high-accuracy equation-of-state work, or cases
        # with significant dimer/polymer effects.
        # Reference: Engineering ToolBox, "Potassium -- Thermophysical
        # Properties vs. Temperature," saturated vapor density data from
        # 350 C to 850 C; Hoffman and Cox, A Preliminary Collation of the
        # Thermodynamic and Transport Properties of Potassium,
        # ORNL-TM-2126/NASA-CR-96741, Oak Ridge National Laboratory, 1968.
        return 2.398e3 * np.exp(-8.698e3 / T_safe)

    def specific_volume_vapor_sat(self, T: Numeric) -> Numeric:
        """Saturated vapor specific volume [m3/kg]."""
        rho = self.vapor_density(T)
        return 1.0 / np.maximum(rho, 1.0e-30)

    def vapor_viscosity(self, T: Numeric) -> Numeric:
        # Potassium vapor-viscosity formula verified from user-provided data.
        """蒸汽动力粘度 [Pa*s]"""
        T = np.asarray(T, dtype=float)
        # Vapor dynamic viscosity [Pa*s].
        # Formula: mu_v = 5.450e-6 + 2.830e-8*T - 5.600e-12*T^2.
        # Range: 350-1000 K; saturated potassium vapor viscosity for low-pressure
        # gas-phase flow resistance and momentum transport from post-melting
        # cold start to medium/high-temperature operation. Do not use for liquid
        # potassium. Above 900 K, K2 dimerization can lower the actual vapor
        # viscosity relative to a monatomic ideal-gas estimate; this polynomial
        # is a globally smoothed engineering fit. Strong nonequilibrium plasma,
        # ionized, or ultrahigh-pressure dense states require collision-integral
        # corrections such as Chapman-Enskog treatment.
        # References: Fink and Leibowitz, Thermophysical Properties of Alkali
        # Metals, ANL-CEN-RSD-82-4, Argonne National Laboratory, 1982; Bystrov,
        # Kagan, et al., Liquid-Metal Heat Pipes: Thermophysical Properties and
        # Processes, Nauka, 1988; Vargaftik, Handbook of Physical Properties of
        # Liquids and Gases, Hemisphere Publishing, 1975.
        return 5.450e-6 + 2.830e-8 * T - 5.600e-12 * T ** 2

    def vapor_heat_capacity(self, T: Numeric) -> Numeric:
        # Potassium vapor heat-capacity formula verified from user-provided data.
        """Vapor constant-pressure heat capacity [J/(kg*K)]."""
        T = np.asarray(T, dtype=float)
        tau = T / 1000.0
        # Vapor K constant-pressure heat capacity [J/(kg*K)].
        # Shomate temperature variable: tau = T/1000.
        # Formula: cp_v = (20.66122 + 0.391869*tau - 0.417344*tau^2
        #                 + 0.145582*tau^3 + 0.003764*tau^-2)/3.90983e-2.
        # The numerator is in J/(mol*K); 3.90983e-2 kg/mol is the molar mass.
        # Range: 1039.54-1800 K; potassium vapor gas-phase heat capacity.
        # Low-pressure or medium/low-density vapor may also be approximated as
        # monatomic ideal gas cp ~= 5*R_u/(2*M_K) ~= 531.6 J/(kg*K).
        # Do not use for near-critical states, high-pressure dense vapor, or
        # cases with significant dimer/polymer effects.
        # Reference: Chase, M. W., Jr., NIST-JANAF Thermochemical Tables,
        # Fourth Edition, Journal of Physical and Chemical Reference Data,
        # Monograph 9, 1998; NIST Chemistry WebBook, SRD 69, Potassium,
        # Gas Phase Heat Capacity Shomate Equation.
        return (
            20.66122
            + 0.391869 * tau
            - 0.417344 * tau ** 2
            + 0.145582 * tau ** 3
            + 0.003764 * tau ** -2
        ) / 3.90983e-2

    def vapor_gas_constant(self) -> float:
        # Potassium vapor ideal-gas constant verified from user-provided data.
        """Specific gas constant for potassium vapor [J/(kg*K)]."""
        # Formula: R_K = R_u/M_K = 8.314462618/3.90983e-2 = 212.65 J/(kg*K).
        # Range/use: engineering checks for low-pressure or medium/low-density
        # potassium vapor. Do not use as a real-gas equation of state in
        # near-critical, high-pressure dense-vapor, liquid, or two-phase regions.
        # Reference: NIST Chemistry WebBook, SRD 69, Potassium, molecular
        # weight 39.0983; NIST Atomic Weights and Isotopic Compositions for
        # Potassium, standard atomic weight 39.0983(1).
        return self.vapor_specific_gas_constant

    def vapor_heat_capacity_ratio(self) -> float:
        # Monatomic ideal-gas heat-capacity ratio for potassium vapor.
        """Ideal-gas heat-capacity ratio gamma [-]."""
        # Formula: gamma_K = 5/3 = 1.6667.
        # Range/use: low-pressure or medium/low-density monatomic potassium
        # vapor checks. Do not use for liquid K, two-phase flow, near-critical
        # states, high-pressure dense vapor, or cases with significant
        # dimer/polymer effects.
        # Reference: Chase, M. W., Jr., NIST-JANAF Thermochemical Tables,
        # Fourth Edition, Journal of Physical and Chemical Reference Data,
        # Monograph 9, 1998.
        return self.vapor_gamma

    def vapor_ideal_cp(self) -> float:
        """Monatomic ideal-gas cp check value [J/(kg*K)]."""
        # Formula: cp_v,ideal = 5/2*R_K ~= 531.6 J/(kg*K).
        return 2.5 * self.vapor_gas_constant()

    def vapor_ideal_cv(self) -> float:
        """Monatomic ideal-gas cv check value [J/(kg*K)]."""
        # Formula: cv_v,ideal = 3/2*R_K ~= 319.0 J/(kg*K).
        return 1.5 * self.vapor_gas_constant()

    def vapor_sound_speed(self, T: Numeric) -> Numeric:
        # Potassium vapor ideal-gas sound-speed formula verified from user-provided data.
        """Ideal-gas vapor sound speed [m/s]."""
        T = np.asarray(T, dtype=float)
        # Formula: a_v = sqrt(gamma_K*R_K*T) ~= 18.83*sqrt(T).
        # Range: 700-1800 K; engineering check for low-pressure or
        # medium/low-density potassium vapor. This is not valid for liquid K,
        # two-phase regions, near-critical states, or high-pressure dense vapor.
        # Saturated two-phase flow requires a two-phase sound-speed model.
        # References: NIST Chemistry WebBook, SRD 69, Potassium, molecular
        # weight 39.0983; NIST Atomic Weights and Isotopic Compositions for
        # Potassium, standard atomic weight 39.0983(1); Chase, M. W., Jr.,
        # NIST-JANAF Thermochemical Tables, Fourth Edition, 1998.
        return np.sqrt(self.vapor_heat_capacity_ratio() * self.vapor_gas_constant() * T)

    def liquid_viscosity(self, T: Numeric) -> Numeric:
        # Potassium liquid-viscosity formula verified from user-provided data.
        """Liquid dynamic viscosity [Pa*s]."""
        T = np.asarray(T, dtype=float)
        T_safe = np.maximum(T, 1.0e-3)
        # Liquid K dynamic viscosity [Pa*s].
        # Formula: mu_l = 4.293e-5 * exp(1017.72/T).
        # Range: 623-1123 K; liquid potassium in the medium/high-temperature
        # range. Above the normal boiling point, the system pressure must
        # maintain liquid or saturated-liquid behavior. Do not use for solid K,
        # near-critical states, or high-accuracy viscosity standards.
        # Reference: Engineering ToolBox, "Potassium -- Thermophysical
        # Properties vs. Temperature," liquid viscosity data from 350 C to 850 C.
        return 4.293e-5 * np.exp(1017.72 / T_safe)

    def viscosity(self, T: Numeric) -> Numeric:
        """Alias for liquid dynamic viscosity [Pa*s]."""
        return self.liquid_viscosity(T)

    def surface_tension(self, T: Numeric) -> Numeric:
        # Potassium liquid surface-tension formula verified from user-provided data.
        """Liquid surface tension [N/m]."""
        T = np.asarray(T, dtype=float)
        # Liquid K surface tension [N/m].
        # Formula: sigma_l = 1.3794e-1 - 6.927e-5*T.
        # Range: 623-1123 K; liquid potassium in the medium/high-temperature
        # range. Above the normal boiling point, the system pressure must
        # maintain liquid or saturated-liquid behavior. Do not use for solid K,
        # near-critical states, or high-accuracy interface-stability analysis.
        # Reference: Engineering ToolBox, "Potassium -- Thermophysical
        # Properties vs. Temperature," liquid surface-tension data from 350 C to 850 C.
        return 1.3794e-1 - 6.927e-5 * T

    def latent_heat(self, T: Numeric) -> Numeric:
        # Potassium latent-heat formula verified from user-provided data.
        """汽化潜热 [J/kg]"""
        T = np.asarray(T, dtype=float)
        # Latent heat of vaporization [J/kg].
        # Formula: h_fg = 2.487169e6 - 396.5976*T - 0.102412*T^2.
        # Range: 350-900 K; polynomial fit to saturated-line potassium latent
        # heat data. It is intended for post-melting cold start through
        # medium/high-temperature heat-pipe operation, and is not valid for
        # superheated vapor, subcooled liquid, solid sublimation, near-critical
        # states, or dense supercritical potassium.
        # Watson check form:
        # h_fg_Watson = 1.966838e6*((2223.0 - T)/1191.0)^0.38,
        # with T_b = 1032.0 K, T_c = 2223.0 K, and
        # h_fg,b = 76.9e3/3.90983e-2 = 1.966838e6 J/kg.
        # References: NIST Chemistry WebBook, SRD 69, Potassium thermochemical
        # data; Faghri, Heat Pipe Science and Technology, Second Edition, 2016;
        # Bystrov, Kagan, et al., Liquid-Metal Heat Pipes, Nauka, 1988; Fink
        # and Leibowitz, Thermophysical Properties of Alkali Metals,
        # ANL-CEN-RSD-82-4, Argonne National Laboratory, 1982.
        return 2.487169e6 - 396.5976 * T - 0.102412 * T ** 2
