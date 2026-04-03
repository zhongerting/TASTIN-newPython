import numpy as np
from Materials.Base import TwoPhaseFluidMaterial, Numeric
from scipy import optimize


class Sodium(TwoPhaseFluidMaterial):
    """
    钠 (Sodium) / NaK 合金流体物性模型

    来源: PROPSODIUM.FOR
    适用范围: 371 K - 2503.7 K (临界温度)
    修改状态: 已向量化 (Vectorized for NumPy)
    """

    def __init__(self):
        super().__init__(name="Sodium", formula="Na")
        self.T_crit = 2503.7  # 临界温度 [K]
        self.P_crit = 25.64e6  # 临界压力

    # ==========================================================================
    # 0. 基础物性 (Heat Capacity) - 修复原文件中的 pass
    # ==========================================================================
    def heat_capacity_liquid(self, T: Numeric) -> Numeric:
        """
        液体定压比热容 Cp [J/(kg·K)]
        来源: 提取自原文件 prandtl_number 中的 CPPTSO 公式
        """
        T = np.asarray(T)
        # 原公式: (1097.73 - 0.556*T + 3.43e-4*T^2)/1000 -> kJ/kgK
        # 转换为 J/kgK
        cp = 1097.73 - 0.556577 * T + 3.43167e-4 * (T ** 2)
        return cp

    def heat_capacity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """通用比热容接口"""
        # 默认为液体比热容
        return self.heat_capacity_liquid(T)

    # ==========================================================================
    # 1. 饱和性质 (Saturation Properties)
    # ==========================================================================

    def saturation_pressure(self, T: Numeric) -> Numeric:
        """
        饱和蒸汽压 P_sat [Pa]
        对应 Fortran: PSATTSO
        """
        T = np.asarray(T)
        # 避免 T=0 导致的除零或 log(0) 错误
        T_safe = np.where(T <= 0, 1.0, T)

        # 原公式: EXP(11.9463 - 12633.73/T - 0.4672*LOG(T)) (MPa)
        val_mpa = np.exp(11.9463 - 12633.73 / T_safe - 0.4672 * np.log(T_safe))

        # T <= 0 返回 0.0
        return np.where(T <= 0, 0.0, val_mpa * 1.0e6)

    def saturation_temperature(self, P: Numeric) -> Numeric:
        """
        饱和温度 T_sat [K]
        对应 Fortran: TSATPSO
        实现: 向量化牛顿迭代法
        """
        P_pa = np.asarray(P)

        # 处理非正压力
        P_safe = np.maximum(P_pa, 1e-5)  # 避免 log(<=0)
        pp_mpa = P_safe / 1.0e6

        # --- 初始猜测 (近似公式) ---
        A = 7.813
        B = 11209.0
        C = 5.249e5

        log_p = np.log(pp_mpa)
        denom = 2.0 * (A - log_p)
        # 防止除零
        denom = np.where(denom == 0, 1e-5, denom)

        # 初始 T_sat
        # 注意: sqrt 内的值必须非负
        delta_sqrt = B ** 2 + 4.0 * C * (A - log_p)
        t_sat = (B + np.sqrt(np.maximum(0, delta_sqrt))) / denom

        # --- 牛顿迭代精修 (Vectorized) ---
        # 固定迭代 10 次 (通常足够收敛)
        for _ in range(10):
            # 防止 t_sat 跑飞到 0 或负数
            t_iter = np.maximum(t_sat, 10.0)

            # 计算当前 P_calc (MPa)
            # 重复 saturation_pressure 的核心逻辑以避免不必要的数组创建
            ln_p_calc = 11.9463 - 12633.73 / t_iter - 0.4672 * np.log(t_iter)

            # f(T) = ln(P_calc) - ln(P_target)
            f_val = ln_p_calc - log_p

            # f'(T) = d(lnP)/dT = 12633.73/T^2 - 0.4672/T
            df_dt = 12633.73 / (t_iter ** 2) - 0.4672 / t_iter
            # 避免导数为0
            df_dt = np.where(df_dt == 0, 1.0, df_dt)

            delta = f_val / df_dt

            # 限制步长，防止发散
            delta = np.clip(delta, -0.2 * t_iter, 0.2 * t_iter)

            t_sat = t_sat - delta

        # 对于原本 P <= 0 的输入，返回一个安全值或报错（这里返回 371.0）
        return np.where(P_pa <= 0, 371.0, t_sat)

    def latent_heat(self, T: Numeric) -> Numeric:
        """
        汽化潜热 h_fg [J/kg]
        对应 Fortran: SHFGSO
        """
        T = np.asarray(T)
        term = 1.0 - T / self.T_crit
        # 截断负值
        term = np.maximum(term, 0.0)

        h_fg_kj = 393.37 * term + 4398.6 * (term ** 0.29302)
        return h_fg_kj * 1000.0

    def enthalpy_saturated_liquid(self, T: Numeric) -> Numeric:
        """
        饱和液体焓 h_f [J/kg]
        对应 Fortran: SFHPSO
        """
        T = np.asarray(T)

        # 分段逻辑
        # Range 1: 371.0 < T < 2000.0 (扩大一点包含低温区，防止 T<371 报错)
        # Range 2: T >= 2000.0

        # 公式 1 （已停用）
        # h_kj = np.abs(1.6582 * T - 365.77 - 4.2395 * (T ** 2) / 1e4 +
        #               1.487 * ((T / 100) ** 3) / 10.0 + 2992.6 / np.maximum(T, 1.0))
        # return h_kj * 1000

        # 由比热容积分而来，规定 300K 为积分初始点，对应焓为 1.075247e5 J/kg
        h = np.abs(1.14389e-4*T*T*T-0.278288*T*T+1.09773e3*T-2.00633e5)

        return h

        # 公式 2
        # val2 = 2128.4 + 0.86496 * T - (self.latent_heat(T) / 1000.0) / 2.0
        #
        # # 组合
        # cond = (T < 2000.0)
        # h_kj = np.where(cond, val1, val2)
        #
        # return h_kj * 1000.0

    def enthalpy_saturated_vapor(self, T: Numeric) -> Numeric:
        """
        饱和蒸汽焓 h_g [J/kg]
        对应 Fortran: SGHPSO
        """
        return self.enthalpy_saturated_liquid(T) + self.latent_heat(T)

    # ==========================================================================
    # 2. 密度与体积 (Density & Volume)
    # ==========================================================================

    def specific_volume_liquid_sat(self, T: Numeric) -> Numeric:
        """饱和液体比体积 v_f [m^3/kg]"""
        T = np.asarray(T)
        term = 1.0 - T / self.T_crit
        term = np.maximum(term, 0.0)

        A = 219.0 + 275.32 * term + 511.58 * np.sqrt(term)
        # 防止 A=0
        return 1.0 / np.maximum(A, 1e-5)

    def specific_volume_vapor_sat(self, T: Numeric) -> Numeric:
        """饱和蒸汽比体积 v_g [m^3/kg]"""
        T = np.asarray(T)
        T_safe = np.maximum(T, 1.0)

        # dP/dT (MPa/K) derived term
        term_deriv = 12633.73 / (T_safe ** 2) - 0.4672 / T_safe
        p_mpa = self.saturation_pressure(T) / 1.0e6

        dp_dt_mpa = term_deriv * p_mpa

        v_f = self.specific_volume_liquid_sat(T)
        h_fg_kj = self.latent_heat(T) / 1000.0

        # 防止除零
        denom = np.maximum(1000.0 * T_safe * dp_dt_mpa, 1e-10)
        v_g = v_f + h_fg_kj / denom

        return v_g

    def specific_volume(self, T: Numeric, P: Numeric) -> Numeric:
        """
        比体积 [m^3/kg] (过冷液或过热汽)
        对应 Fortran: SVPTSO
        向量化实现：同时计算过冷和过热，最后合并
        """
        T = np.asarray(T)
        P = np.asarray(P)
        P_mpa = P / 1.0e6
        T_sat = self.saturation_temperature(P)

        # 1. 过冷液体 (Subcooled) 计算
        v_sub = self.specific_volume_liquid_sat(T)

        # 2. 过热蒸汽 (Superheated) 计算 - 维里方程迭代
        R = 8.31441

        # 保护 T, P 防止数学错误
        T_sh = np.maximum(T, 1.0)
        P_sh = np.maximum(P_mpa, 1e-10)  # atm approx conversion requires P>0

        app = P_sh / 0.101325
        aak2 = np.exp(-9.95845 + 9117.50383 / T_sh)
        aak4 = np.exp(-24.59115 + 20660.83388 / T_sh)

        # 迭代求解 aan1
        aan1 = np.full_like(T, 0.8, dtype=float)

        for _ in range(20):
            num = aak4 * (app ** 3) * (aan1 ** 4) + aak2 * app * (aan1 ** 2) + aan1 - 1.0
            den = 4.0 * aak4 * (app ** 3) * (aan1 ** 3) + 2.0 * aak2 * app * aan1 + 1.0
            # 避免 den=0
            den = np.where(den == 0, 1.0, den)

            delta = num / den
            aan1 -= delta
            # 向量化不需要 break，跑完固定次数即可

        aan2 = (aan1 ** 2) * app * aak2
        aan4 = (aan1 ** 4) * (app ** 3) * aak4

        am = 22.98977 * aan1 + 22.98977 * 2.0 * aan2 + 22.98977 * 4.0 * aan4

        # sv calculation
        # 防止 am=0
        am = np.where(am == 0, 1.0, am)
        v_super = (T_sh * R) / (1000.0 * P_sh * am)

        # 3. 合并结果
        is_subcooled = T < T_sat
        return np.where(is_subcooled, v_sub, v_super)

    def density(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """密度 rho [kg/m^3]"""
        v = self.specific_volume(T, P)
        return 1.0 / np.maximum(v, 1e-10)

    # ==========================================================================
    # 3. 焓与热工水力专用接口
    # ==========================================================================

    def enthalpy_gas(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """气体焓 h_g [J/kg]"""
        T = np.asarray(T)
        val_kj = -126.03 + 0.8844 * T
        return val_kj * 1000.0

    def enthalpy(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """通用焓接口"""
        T = np.asarray(T)
        T_sat = self.saturation_temperature(P)

        h_liq = self.enthalpy_saturated_liquid(T)
        h_gas = self.enthalpy_gas(T, P)

        return np.where(T <= T_sat, h_liq, h_gas)

    def temperature_from_enthalpy(self, h: Numeric, P: Numeric) -> Numeric:
        """
        根据焓反求温度 [K] (Polynomial Approximation Version)

        使用用户提供的高精度多项式拟合公式替代迭代求解，显著提高计算速度和稳定性。
        """
        # 0. 维度安全处理
        h_in = np.asarray(h)
        P_in = np.asarray(P)
        original_shape = h_in.shape

        # 展平以便向量化计算
        h_vec = np.atleast_1d(h_in).flatten()
        P_vec = np.atleast_1d(P_in).flatten()

        # 1. 计算饱和线参数 (作为相态分界)
        T_sat = self.saturation_temperature(P_vec)
        h_f = self.enthalpy_saturated_liquid(T_sat)
        h_g = self.enthalpy_saturated_vapor(T_sat)

        # 2. 初始化结果数组
        # 默认填充为饱和温度 (覆盖两相区 h_f <= h <= h_g 的情况)
        T_result = np.array(T_sat, copy=True)

        # ==========================================
        # Case A: 过冷液体 (Subcooled Liquid)
        # 判据: h < h_f
        # 公式: T = c3*h^3 + c2*h^2 + c1*h + c0
        # ==========================================
        mask_liq = h_vec < h_f

        if np.any(mask_liq):
            # 由焓公式计算而得（已停用）
            # h_liq = h_vec[mask_liq]
            #
            # # 拟合系数
            # c3 = -4.43113194093226e-17
            # c2 = 1.28278235899215e-10
            # c1 = 0.000673347725489693
            # c0 = 226.415517947328
            #
            # # 多项式计算
            # T_calc = c3 * (h_liq ** 3) + c2 * (h_liq ** 2) + c1 * h_liq + c0
            #
            # # 物理下限保护 (防止焓值极低时温度越界)
            # T_result[mask_liq] = np.maximum(T_calc, 200.0)

            # 改为新的比热容积分所得焓公式计算获得
            h_liq = h_vec[mask_liq]
            # 拟合系数
            c3 = -1.6918e-16
            c2 = 2.922e-10
            c1 = 9.7797e-4
            c0 = 192.6584
            # 多项式计算
            T_calc = c3 * (h_liq ** 3) + c2 * (h_liq ** 2) + c1 * h_liq + c0

            # 物理下限保护 (防止焓值极低时温度越界)
            T_result[mask_liq] = np.maximum(T_calc, 300.0)


        # ==========================================
        # Case B: 两相混合区 (Two-Phase Mixture)
        # ==========================================
        # 已在初始化中设为 T_sat，无需额外操作

        # ==========================================
        # Case C: 过热气体 (Superheated Gas)
        # 判据: h > h_g
        # 公式: T = a*h^2 + b*h + c
        # ==========================================
        mask_gas = h_vec > h_g

        if np.any(mask_gas):
            h_gas = h_vec[mask_gas]

            # 拟合系数
            # 注意: 原公式中线性项系数对应变量修正为 h (原需求写为 T 疑似笔误)
            a = 2.94846761640699e-09
            b = -0.0267871874558033
            c = 61103.0660018335  # 原值保留 (注: 此常数项较大，用于平衡 h^2 的量级)

            # 多项式计算
            T_gas = a * (h_gas ** 2) + b * h_gas + c

            T_result[mask_gas] = T_gas

        # 3. 还原原始形状
        if len(original_shape) == 0:
            return float(T_result[0])
        else:
            return T_result.reshape(original_shape)

        # h = np.asarray(h)
        # h_kj = h / 1000.0
        # return 142.54 + 1.13064 * h_kj

    # def temperature_from_enthalpy(self, h: Numeric, P: Numeric) -> Numeric:
    #     """
    #     根据焓反求温度 [K] (High-Precision Quadratic Version)
    #
    #     改进点：
    #     1. 使用针对 Sodium.py 正向公式拟合的二次多项式，解决了线性公式在 600K 处
    #        100K+ 的偏差问题。
    #     2. 精度验证：在 400K-1200K 范围内误差 < 2K。
    #     3. 保持解析形式，无迭代，计算极快且对 ODE 求解器稳定。
    #     """
    #     # 0. 维度安全处理
    #     h_in = np.asarray(h)
    #     P_in = np.asarray(P)
    #     original_shape = h_in.shape
    #
    #     h_vec = np.atleast_1d(h_in).flatten()
    #     P_vec = np.atleast_1d(P_in).flatten()
    #
    #     # 1. 计算基准参数
    #     T_sat = self.saturation_temperature(P_vec)
    #     h_f = self.enthalpy_saturated_liquid(T_sat)
    #     h_g = self.enthalpy_saturated_vapor(T_sat)
    #
    #     # 2. 初始化结果
    #     T_result = np.array(T_sat, copy=True)
    #
    #     # ==========================================
    #     # Case A: 液相 (Liquid)
    #     # ==========================================
    #     mask_liq = h_vec < h_f
    #
    #     if np.any(mask_liq):
    #         h_liq_kj = h_vec[mask_liq] / 1000.0
    #
    #         # [核心修复] 使用高精度二次拟合公式
    #         # 拟合数据点: (h=233, T=400), (h=514, T=600), (h=809, T=800), (h=1113, T=1000)
    #         # 拟合公式: T = 235.24 + 0.7234 * h - 3.32e-5 * h^2
    #
    #         c0 = 235.24
    #         c1 = 0.7234
    #         c2 = -3.32e-5
    #
    #         # T = c0 + c1*h + c2*h^2
    #         T_liq = c0 + c1 * h_liq_kj + c2 * (h_liq_kj ** 2)
    #
    #         # 物理下限保护
    #         T_liq = np.maximum(T_liq, 200.0)
    #
    #         T_result[mask_liq] = T_liq
    #
    #     # ==========================================
    #     # Case C: 气相 (Gas)
    #     # ==========================================
    #     mask_gas = h_vec > h_g
    #
    #     if np.any(mask_gas):
    #         h_gas_kj = h_vec[mask_gas] / 1000.0
    #         # 气相公式保持不变 (线性精度足够)
    #         # T = (h + 126.03) / 0.8844
    #         T_gas = (h_gas_kj + 126.03) / 0.8844
    #         T_result[mask_gas] = T_gas
    #
    #     # 3. 还原形状
    #     if len(original_shape) == 0:
    #         return float(T_result[0])
    #     else:
    #         return T_result.reshape(original_shape)

    def density_from_enthalpy(self, h: Numeric, P: Numeric) -> Numeric:
        """根据焓求密度 (快速公式)"""
        T = self.temperature_from_enthalpy(h, P)
        rho = 938.399 - 0.215438 * T - 1.71386e-5 * (T ** 2)
        return rho

    def liquid_density_derivative_T(self, T: Numeric) -> Numeric:
        """d(rho)/dT"""
        T = np.asarray(T)
        term = np.sqrt(np.maximum(1.0 - T / self.T_crit, 1e-10))
        drho_dt = -275.32 / self.T_crit - 511.58 * 0.5 / (term * self.T_crit)
        return drho_dt

    # ==========================================================================
    # 4. 输运性质 (Transport Properties)
    # ==========================================================================

    def conductivity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """液体导热系数 [W/(m·K)]"""
        T = np.asarray(T)
        k = np.abs(13.9723 + 0.0331088 * T - 0.0000222398 * (T ** 2))
        return k

    def viscosity_liquid(self, T: Numeric) -> Numeric:
        """液体粘度 [Pa·s]
        :param T:
        :return:
        """
        T = np.asarray(T)
        mu = np.abs(0.0016363 - 4.58259e-6 * T + 4.94085e-9 * (T ** 2) - 1.84513e-12 * (T ** 3))
        return mu

    def viscosity_vapor(self, T: Numeric) -> Numeric:
        """气体粘度 [Pa·s]"""
        T = np.asarray(T)
        return 5.196e-7 * (T ** 0.675)

    def viscosity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """通用粘度接口 (默认液体)"""
        return self.viscosity_liquid(T)

    def prandtl_number(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """普朗特数"""
        mu = self.viscosity_liquid(T)
        k = self.conductivity(T)
        cp = self.heat_capacity_liquid(T)
        return mu * cp / k

    def surface_tension(self, T: Numeric) -> Numeric:
        """表面张力 [N/m]"""
        T = np.asarray(T)
        term = np.maximum(1.0 - T / self.T_crit, 0.0)
        sigma_mn = 240.5 * (term ** 1.126)
        return sigma_mn / 1000.0

    # ==========================================================================
    # 5. 可压缩性与声速
    # ==========================================================================

    def adiabatic_compressibility_liquid(self, T: Numeric) -> Numeric:
        """液体绝热压缩系数 [1/Pa]"""
        T = np.asarray(T)
        tm = 371.0
        B = 3.2682
        term = (T - tm) / (self.T_crit - tm)

        val_inv_mpa = 1.717 * (1.0 + term / B) / ((1.0 - term) * 1.0e4)
        return val_inv_mpa / 1.0e6

    def sound_velocity(self, T: Numeric) -> Numeric:
        """声速 [m/s]"""
        T = np.asarray(T)

        # 1. 低温段公式
        c_low = 2660.7 - 0.37667 * T - 9.0356 * (T ** 2) / 1.0e5

        # 2. 高温段公式 (>1773K)
        beta = self.adiabatic_compressibility_liquid(T) * 1e6  # 1/MPa
        v = self.specific_volume_liquid_sat(T)
        # 防止 sqrt 负数
        arg = np.maximum(beta / v, 1e-10)
        c_high = 1000.0 / np.sqrt(arg)

        return np.where(T < 1773.0, c_low, c_high)

    # ==========================================================================
    # 6. 偏导数 (使用 Lambda 时需注意数组兼容性)
    # ==========================================================================

    @staticmethod
    def _finite_diff_P(func, P, eps=5.0):
        """向量化的有限差分"""
        v1 = func(P)
        v2 = func(P - eps)
        return (v1 - v2) / eps

    def vapor_density_derivative_P(self, P: Numeric) -> Numeric:
        func = lambda p: 1.0 / self.specific_volume_vapor_sat(self.saturation_temperature(p))
        return self._finite_diff_P(func, P)

    def liquid_density_derivative_P(self, P: Numeric) -> Numeric:
        func = lambda p: 1.0 / self.specific_volume_liquid_sat(self.saturation_temperature(p))
        return self._finite_diff_P(func, P)

    def vapor_enthalpy_derivative_P(self, P: Numeric) -> Numeric:
        func = lambda p: self.enthalpy_saturated_vapor(self.saturation_temperature(p))
        return self._finite_diff_P(func, P)

    def liquid_enthalpy_derivative_P(self, P: Numeric) -> Numeric:
        func = lambda p: self.enthalpy_saturated_liquid(self.saturation_temperature(p))
        return self._finite_diff_P(func, P)


# Alias
NaK = Sodium
