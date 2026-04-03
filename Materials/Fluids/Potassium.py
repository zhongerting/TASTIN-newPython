import numpy as np
from Materials.Base import TwoPhaseFluidMaterial, Numeric


class Potassium(TwoPhaseFluidMaterial):
    """
    钾 (Potassium) 流体物性模型

    来源: PropPotassium.for
    数据源: 钱增源《低熔点金属的热物性》
    熔点: 63.2 C (336.35 K)
    """

    def __init__(self):
        super().__init__(name="Potassium", formula="K")
        self.T_melt = 336.35
        # 能量单位换算系数 (cal -> J)，源自 Fortran 代码中的隐式常数
        self._J_per_cal = 4.1858518

    # ==========================================================================
    # 1. PSATTPO: 饱和蒸气压 (Pressure)
    # ==========================================================================
    def saturation_pressure(self, T: Numeric) -> Numeric:
        """
        饱和蒸汽压 P_sat [Pa]
        对应 Fortran: PSATTPO
        公式: P(MPa) = 4.0168E11 * 10^(-4625.3/T) * T^(-0.7)
        """
        T = np.asarray(T)
        # 避免 T=0 导致除零或幂运算错误
        T_safe = np.where(T <= 0, 1.0, T)

        # P_MPa = 4.0168e11 * 10**(-4625.3/T) * T**(-0.7)
        val_mpa = 4.0168e11 * np.power(10.0, -4625.3 / T_safe) / (T_safe ** 0.7)

        # T <= 0 返回 0.0
        return np.where(T <= 0, 0.0, val_mpa * 1.0e6)

    # ==========================================================================
    # 2. TSATPPO: 饱和温度 (Temperature)
    # ==========================================================================
    def saturation_temperature(self, P: Numeric) -> Numeric:
        """
        饱和温度 T_sat [K]
        对应 Fortran: TSATPPO
        实现方式: 向量化的牛顿迭代法反求 PSATTPO
        """
        P_pa = np.asarray(P)

        # 1. 处理低压边界 (P <= 1e-5 -> T_melt)
        mask_low = P_pa <= 1e-5

        # 2. 初始猜测
        t_guess = np.full_like(P_pa, 800.0, dtype=float)
        t_guess[mask_low] = self.T_melt

        # 目标 ln(P_MPa)，防止 log(0)
        P_safe_mpa = np.maximum(P_pa, 1e-10) / 1.0e6
        target_log_p = np.log(P_safe_mpa)

        ln10 = np.log(10.0)
        B = 4625.3

        # 3. 向量化迭代 (固定 20 步，对于牛顿法通常足够收敛)
        for _ in range(20):
            # 仅对未达到低压边界的点进行迭代
            # 这里为了保持数组形状一致，对所有点计算，最后用 mask 还原即可，或者全计算无妨

            # 计算当前 P_calc (MPa) 对应的压力
            # P_MPa = A * 10^(-B/T) * T^(-0.7)
            # 为了防止 t_guess 跑飞到 0导致报错
            t_iter = np.maximum(t_guess, 10.0)

            # 直接复用公式计算 p_calc_mpa (避免递归调用 saturation_pressure 的额外开销)
            # ln(P) = ln(A) - (B*ln10)/T - 0.7*ln(T)
            # 我们直接在对数域迭代更稳定: f(T) = ln(P_calc) - ln(P_target) = 0

            # log(4.0168e11) approx 26.719
            ln_p_calc = np.log(4.0168e11) - (B * ln10) / t_iter - 0.7 * np.log(t_iter)

            f_val = ln_p_calc - target_log_p

            # 导数 f'(T) = d(ln P)/dT = (B * ln10) / T^2 - 0.7 / T
            df_dt = (B * ln10) / (t_iter ** 2) - 0.7 / t_iter

            # 避免导数为0
            df_dt = np.where(df_dt == 0, 1.0, df_dt)

            delta = f_val / df_dt

            # 阻尼防止发散 (限制步长不超过当前温度的 20%)
            max_step = 0.2 * t_iter
            delta = np.clip(delta, -max_step, max_step)

            t_guess = t_guess - delta

            # 检查收敛 (可选，若追求极致速度可提前跳出，但向量化下全跑完往往更快)
            if np.all(np.abs(delta) < 0.001):
                break

        # 还原低压边界值
        return np.where(mask_low, self.T_melt, t_guess)

    # ==========================================================================
    # 3. PHSFPO: 熔化潜热
    # ==========================================================================
    def latent_heat_fusion(self, T: Numeric) -> Numeric:
        val_kj = 14.17 * self._J_per_cal
        return val_kj * 1000.0

    # ==========================================================================
    # 4. PHFGPO: 汽化潜热
    # ==========================================================================
    def latent_heat(self, T: Numeric) -> Numeric:
        val_kj = 473.9 * self._J_per_cal
        return val_kj * 1000.0

    # ==========================================================================
    # 5. HSATLTPO: 饱和液体焓
    # ==========================================================================
    def enthalpy_saturated_liquid(self, T: Numeric) -> Numeric:
        T = np.asarray(T)
        t_c = T - 273.15

        # 分段函数向量化
        # Range 1: 63.2 <= t_c <= 800.0
        # Range 2: else

        cond1 = (t_c >= 63.2) & (t_c <= 800.0)

        h1 = 56.179 + 0.84074 * t_c - 1.5844e-4 * (t_c ** 2) + 1.0499e-7 * (t_c ** 3)
        h2 = 0.7104 * t_c + 1.0385e-3 * (t_c ** 2)

        h_kj = np.where(cond1, h1, h2)
        return h_kj * 1000.0

    # ==========================================================================
    # 6. HSATGTPO: 饱和气体焓
    # ==========================================================================
    def enthalpy_saturated_vapor(self, T: Numeric) -> Numeric:
        return self.enthalpy_saturated_liquid(T) + self.latent_heat(T)

    # ==========================================================================
    # 7. SSATSTPO: 饱和固态熵
    # ==========================================================================
    def entropy_saturated_solid(self, T: Numeric) -> Numeric:
        T = np.asarray(T)
        t_c = T - 273.15

        # 范围: 0 <= t_c <= 63.2
        cond = (t_c >= 0.0) & (t_c <= 63.2)

        # 计算正常值
        # 防止 log10(0) 或负数
        T_safe = np.maximum(T, 1e-10)
        s_calc = 0.32938 * np.log10(T_safe) + 2.0770e-8 * T_safe - 1.36986

        # 计算边界值 (63.2C -> 336.35K)
        T_bound = 336.35
        s_bound = 0.32938 * np.log10(T_bound) + 2.0770e-8 * T_bound - 1.36986

        s_kj = np.where(cond, s_calc, s_bound)
        return s_kj * 1000.0

    # ==========================================================================
    # 8. SSATLTPO: 饱和液态熵
    # ==========================================================================
    def entropy_saturated_liquid(self, T: Numeric) -> Numeric:
        T = np.asarray(T)
        T_safe = np.maximum(T, 1e-10)
        s_kj = 2.18919 * np.log10(T_safe) + 4.88620e-4 * T_safe + 1.5718e-7 * T_safe - 5.04665
        return s_kj * 1000.0

    # ==========================================================================
    # 9. VSATGTPO: 饱和气体比体积
    # ==========================================================================
    def specific_volume_vapor_sat(self, T: Numeric) -> Numeric:
        T = np.asarray(T)

        # 1. P (Pa)
        P_pa = self.saturation_pressure(T)

        # 2. dP/dT
        B = 4625.3
        ln10 = np.log(10.0)
        T_safe = np.maximum(T, 1e-10)
        dlnP_dT = (B * ln10) / (T_safe ** 2) - 0.7 / T_safe
        dp_dt = P_pa * dlnP_dT

        # 3. rho_f -> v_f
        rho_f = self.density(T)
        # 避免除以0
        rho_safe = np.where(rho_f > 0, rho_f, 1.0)
        v_f = np.where(rho_f > 0, 1.0 / rho_safe, 0.001)

        # 4. h_fg
        h_fg = self.latent_heat(T)

        # 5. Clapeyron
        # v_g = v_f + h_fg / (T * dp_dt)
        # 保护分母
        denom = T * dp_dt
        # 只有当 dp_dt > 0 且 P > 1e-5 时计算有效
        valid_mask = (dp_dt > 1e-20) & (P_pa > 1e-5)

        term2 = np.divide(h_fg, denom, out=np.zeros_like(h_fg), where=valid_mask)
        v_g = v_f + term2

        return np.where(valid_mask, v_g, 0.0)

    # ==========================================================================
    # 10. PDENTPO: 密度
    # ==========================================================================
    def density(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        T = np.asarray(T)
        t_c = T - 273.15

        # 分段逻辑: Solid < 63.2 <= Liquid < 756.5 <= Gas < 2100

        # Solid
        rho_sol = 864.0 - 0.24162 * t_c

        # Liquid
        rho_liq = 841.5 - 0.2127 * t_c - 2.7e-5 * (t_c ** 2) + 4.77e-9 * (t_c ** 3)

        # Gas
        T_safe = np.maximum(T, 1e-10)
        rho_gas = 1.884233822e9 * np.power(10.0, -4625.3 / T_safe) / (T_safe ** 1.7)

        # 组合
        cond_sol = t_c < 63.2
        cond_liq = t_c < 756.5
        cond_gas = t_c < 2100.0

        # np.select(conditions, choices, default)
        return np.select(
            [cond_sol, cond_liq, cond_gas],
            [rho_sol, rho_liq, rho_gas],
            default=0.0
        )

    # ==========================================================================
    # 11. PSUPTDPO: 过热蒸气压
    # ==========================================================================
    def pressure_superheated(self, T: Numeric, rho: Numeric) -> Numeric:
        T = np.asarray(T)
        rho = np.asarray(rho)

        # VV = 39098 / rho
        rho_safe = np.where(rho > 0, rho, 1.0)
        vv = 39098.0 / rho_safe

        T_safe = np.maximum(T, 1.0)
        con1 = -0.3812e-5 * T * np.exp(7007.0 / T_safe)
        con2 = 6.193e-2 * np.exp(8168.0 / T_safe)
        con3 = -0.4615 * np.exp(10059.0 / T_safe)

        term = 1.0 + con1 / vv + con2 / (vv ** 2) + con3 / (vv ** 3)
        p_mpa = (term / vv) * T * 8.314

        res = p_mpa * 1.0e6
        return np.where(rho > 0, res, 0.0)

    # ==========================================================================
    # 12. CEXPPO: 体膨胀系数
    # ==========================================================================
    def thermal_expansion_cubic(self, T: Numeric) -> Numeric:
        T = np.asarray(T)
        t_c = T - 273.15

        # 插值数据
        tw_data = np.array([
            336.65, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0,
            1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0
        ]) + 273.15
        ex_data = np.array([
            30.86, 29.61, 29.06, 29.60, 30.74, 32.17, 33.74, 35.39,
            37.16, 39.28, 42.20, 46.74, 54.26, 66.83, 87.49, 120.5
        ])

        # 1. 基础插值
        val_interp = np.interp(T, tw_data, ex_data)

        # 2. 覆盖逻辑 (63.2 <= t_c <= 1250.0)
        val_poly = 1.0 + 2.581e-4 * t_c + 1.308e-7 * (t_c ** 2) + 1.98e-12 * (t_c ** 3)

        cond_poly = (t_c >= 63.2) & (t_c <= 1250.0)
        return np.where(cond_poly, val_poly, val_interp)

    # ==========================================================================
    # 13 & 14. 压缩系数 (插值)
    # ==========================================================================
    def adiabatic_compressibility_liquid(self, T: Numeric) -> Numeric:
        tt_data = np.array([100, 200, 300, 400, 500, 600, 700, 800]) + 273.15
        com_data = np.array([35.033, 38.186, 41.725, 45.720, 50.236, 55.367, 61.229, 67.940])
        val_inv_mpa = np.interp(T, tt_data, com_data) * 1.0e-5
        return val_inv_mpa / 1.0e6

    def isothermal_compressibility_liquid(self, T: Numeric) -> Numeric:
        tt_data = np.array([100, 200, 300, 400, 500, 600, 700, 800]) + 273.15
        com_data = np.array([39.171, 44.044, 49.617, 55.962, 63.171, 71.317, 80.570, 80.865])
        val_inv_mpa = np.interp(T, tt_data, com_data) * 1.0e-5
        return val_inv_mpa / 1.0e6

    # ==========================================================================
    # 15. 占位符
    # ==========================================================================
    def specific_volume_gas_working(self, P: Numeric, T: Numeric) -> Numeric:
        return 0.0

    def temperature_from_pressure_enthalpy(self, P: Numeric, h: Numeric) -> Numeric:
        return 0.0

    def entropy(self, P: Numeric, T: Numeric) -> Numeric:
        return 0.0

    # ==========================================================================
    # 16. VISCVPO: 气体动力粘度
    # ==========================================================================
    def viscosity_vapor(self, T: Numeric) -> Numeric:
        # 公式一直有效，简单返回
        return 4.86372e-6 + 1.15683e-8 * T

    # ==========================================================================
    # 17. VISCLPO: 液体动力粘度
    # ==========================================================================
    def viscosity_liquid(self, T: Numeric) -> Numeric:
        T = np.asarray(T)
        rho = self.density(T)
        # 保护除法
        T_safe = np.maximum(T, 1.0)

        # Base formula: A * rho^(1/3) * exp(B * rho / T)
        # Range 1: 336.35 < T < 653.15
        # Coeffs: A=1.131E-5, B=0.68
        v1 = 1.131e-5 * (rho ** (1.0 / 3.0)) * np.exp(0.68 * rho / T_safe)

        # Range 2: else (High temp)
        # Coeffs: A=7.99E-6, B=0.978
        v2 = 7.99e-6 * (rho ** (1.0 / 3.0)) * np.exp(0.978 * rho / T_safe)

        cond1 = (T > 336.35) & (T < 653.15)
        return np.where(cond1, v1, v2)

    def viscosity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        return self.viscosity_liquid(T)

    # ==========================================================================
    # 18. THCONSPO: 热导率
    # ==========================================================================
    def conductivity_solid(self, T: Numeric) -> Numeric:
        # Formula valid for 2.3 < T < 336.35
        return (0.301 - 1.44e-4 * T) * 418.68

    def conductivity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        T = np.asarray(T)
        T_safe = np.maximum(T, 1.0)
        # Liquid formula
        return 43.8 - 2.22e-2 * T + 3.95e3 / T_safe

    def conductivity_vapor(self, T: Numeric) -> Numeric:
        tt_data = np.array([700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0])
        th_data = np.array([118.0, 148.0, 175.0, 200.0, 225.0, 248.0, 269.0])
        k_val = np.interp(T, tt_data, th_data)
        return k_val * 1.0e-3

    def conductivity_gas(self, T: Numeric, P: Numeric) -> Numeric:
        return 0.0

    # ==========================================================================
    # 21 & 22. 普朗特数
    # ==========================================================================
    def prandtl_liquid(self, T: Numeric) -> Numeric:
        tw_k = np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0]) + 273.15
        pr_vals = np.array([729.0, 542.0, 455.0, 405.0, 367.0, 349.0, 346.0, 355.0, 374.0]) * 1.0e-5
        return np.interp(T, tw_k, pr_vals)

    def prandtl_gas(self, T: Numeric) -> Numeric:
        return 0.8

    def prandtl_number(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        return self.prandtl_liquid(T)

    # ==========================================================================
    # 23-25. 比热容
    # ==========================================================================
    def heat_capacity_solid(self, T: Numeric) -> Numeric:
        # 原逻辑是 100 < T < 336 使用公式，否则复用。
        # 向量化下直接返回公式即可，除非有严格定义
        return 537.9 + 0.8002 * T

    def heat_capacity_liquid(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        t_c = T - 273.15
        # 原逻辑范围 336 < T < 1523.15，否则复用。直接公式。
        val_kj = 0.8389 - 3.6741e-4 * t_c + 4.592e-7 * (t_c ** 2)
        return val_kj * 1000.0

    def heat_capacity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        return self.heat_capacity_liquid(T, P)

    def heat_capacity_saturated_liquid(self, T: Numeric) -> Numeric:
        T = np.asarray(T)

        # Range 1: 336.35 < T < 1523.15 -> Use Liquid Formula
        cp_liq = self.heat_capacity_liquid(T)

        # Range 2: Else -> Linear formula
        val_kj = (0.1285 + 1.9116e-4 * T) * 4.18
        cp_lin = val_kj * 1000.0

        cond = (T > 336.35) & (T < 1523.15)
        return np.where(cond, cp_liq, cp_lin)

    # ==========================================================================
    # 26. SOUNDVLPO: 液态声速
    # ==========================================================================
    def sound_velocity(self, T: Numeric) -> Numeric:
        T = np.asarray(T)
        t_c = T - 273.15

        # 1. 插值
        wt_data = np.array([
            100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650,
            700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250
        ], dtype=float)
        c_data = np.array([
            1868, 1841, 1815, 1788, 1761, 1734, 1707, 1680, 1653, 1626,
            1599, 1571, 1545, 1518, 1491, 1484, 1437, 1410, 1383, 1356,
            1329, 1302, 1276, 1249
        ], dtype=float)
        c_interp = np.interp(t_c, wt_data, c_data)

        # 2. 覆盖 (100 <= t_c <= 800)
        c_poly = 1922.3 - 0.539 * t_c

        cond = (t_c >= 100.0) & (t_c <= 800.0)
        return np.where(cond, c_poly, c_interp)

    # ==========================================================================
    # 27. STLPO: 表面张力
    # ==========================================================================
    def surface_tension(self, T: Numeric) -> Numeric:
        t_c = T - 273.15
        sigma_mn = 115.7 - 0.064 * t_c
        return sigma_mn / 1000.0

    # ==========================================================================
    # 28. 补充通用接口 (基于 Sodium.py 逻辑移植)
    # ==========================================================================

    def enthalpy(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        通用焓接口 h [J/kg]
        逻辑:
          如果 T <= T_sat(P): 返回饱和液体焓 (忽略过冷度对焓的微小影响)
          如果 T >  T_sat(P): 返回饱和气体焓 (暂时忽略过热度，或可视作含气率x=1)
        """
        T = np.asarray(T)
        P = np.asarray(P)

        # 获取当前压力下的饱和温度
        T_sat = self.saturation_temperature(P)

        # 分别计算液/气焓
        h_liq = self.enthalpy_saturated_liquid(T)
        h_gas = self.enthalpy_saturated_vapor(T)

        # 根据温度判断相态
        return np.where(T <= T_sat, h_liq, h_gas)

    def temperature_from_enthalpy(self, h: Numeric, P: Numeric) -> Numeric:
        """
        根据焓值反求温度 T [K]
        实现: 向量化牛顿迭代法 (Newton-Raphson)
        f(T) = h_liq(T) - h_target = 0
        f'(T) = Cp_liq(T)
        """
        h_target = np.asarray(h)

        # 初始猜测 (假设在液相区，取中间温度 800K)
        T_guess = np.full_like(h_target, 800.0, dtype=float)

        # 迭代求解
        for _ in range(10):  # 通常 5-10 步收敛
            h_calc = self.enthalpy_saturated_liquid(T_guess)
            cp = self.heat_capacity_liquid(T_guess)

            # f(T) / f'(T)
            diff = h_calc - h_target
            # 保护除零
            delta = diff / np.maximum(cp, 1.0)

            # 更新
            T_guess = T_guess - delta

            # 简单的收敛判定 (可选)
            if np.all(np.abs(delta) < 0.01):
                break

        return T_guess

    def _finite_diff_P(self, func, P, eps=100.0):
        """
        向量化的有限差分辅助函数 (针对压力 P)
        eps 取 100 Pa 以保证数值稳定性
        """
        v1 = func(P)
        v2 = func(P - eps)
        return (v1 - v2) / eps

    def liquid_density_derivative_P(self, P: Numeric) -> Numeric:
        """
        饱和液密度对压力的导数 drho/dP | sat
        逻辑: d(rho_sat)/dP
        """
        # 定义辅助 lambda: 给定 P -> T_sat -> rho_sat
        func = lambda p: self.density(self.saturation_temperature(p))
        return self._finite_diff_P(func, P)

    def liquid_density_derivative_T(self, T: Numeric) -> Numeric:
        """
        液密度对温度的导数 drho/dT [kg/(m^3·K)]
        实现: 解析求导 (Analytical Derivative)

        原公式 (Density):
        rho = 841.5 - 0.2127*t_c - 2.7e-5*t_c^2 + 4.77e-9*t_c^3

        导数 d(rho)/d(t_c):
        drho = -0.2127 - 5.4e-5*t_c + 1.431e-8*t_c^2
        """
        T = np.asarray(T)
        t_c = T - 273.15

        # 解析导数公式
        drho_dt = -0.2127 - 5.4e-5 * t_c + 1.431e-8 * (t_c ** 2)

        # 仅对液相范围有效，其他范围设为 0 或对应导数
        # (这里主要关注液相，直接返回公式计算值)
        return drho_dt
