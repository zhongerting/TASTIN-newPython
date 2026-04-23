import numpy as np
from Materials.Base import TwoPhaseFluidMaterial, Numeric
from scipy import optimize

class SodiumPotassium78(TwoPhaseFluidMaterial):
    """
    钠钾合金 (Sodium-Potassium78) / NaK 合金流体物性模型

    本系统对应工质为：NaK78，并且只考虑液态条件下的物性。

    主要来源为：《Liquid Metals Handbook》
               《Sodium-NaK Engineering Handbook》
    适用范围: 371 K - 2503.7 K (临界温度)
    修改状态: 已向量化 (Vectorized for NumPy)
    """

    #TODO：添加气态条件下的物性模型，并将其与液态条件下的模型设置通用接口
    #TODO：添加固态物性模型

    def __init__(self):
        super().__init__(name="Sodium", formula="Na")
        """
        来源于 Kopp 定律扩展：在常温下，固体化合物或固溶体合金的摩尔等压热容 ($C_p$) 近似等于构成该物质的各单质或简单化合物的摩尔热容之和。
        Kubaschewski, O., Alcock, C. B., & Spencer, P. J. (1993). Materials Thermochemistry (6th Edition).
        """
        self.T_crit = 2503.7  # 临界温度 [K]
        self.P_crit = 18.00e6  # 临界压力

        self.N_Na = 0.324   # 金属钠摩尔分数
        self.N_K = 0.676    # 金属钾摩尔分数

        self.m_Na = 0.22    # 金属钠质量分数
        self.m_K = 0.78    # 金属钾质量分数

    # ==========================================================================
    # 1. 密度 (Density) 和 比体积 (Specific Volume)
    # ==========================================================================
    def density(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：密度 rho [kg/m³]
        输入：温度 T [K]，压力 P [Pa]
        来源: 采用质量加权计算，$$rho_{NaK} = w_{Na} \rho_{Na} + w_K \rho_{K}$$
        """
        return 1.0 / self.specific_volume(T, P)
    
    def specific_volume(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：比体积 V [m³/kg]
        输入：温度 T [K]，压力 P [Pa]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """

        # 液态 NaK78 物性
        T = np.asarray(T)

        # 摄氏度和开尔文换算
        T_ssd = T - 273.15  # 转换为摄氏度

        # 液态 K 比体积
        rho_K = 1000. * (0.8415 - 2.172e-4 * T_ssd - 2.7e-8 * T_ssd ** 2 + 4.77e-12 * T_ssd ** 3)
        specific_volume_K = 1.0 / rho_K

        #液态 Na 比体积
        rho_Na = 1000. * (0.9501 - 2.2976e-4 * T_ssd - 1.46e-8 * T_ssd ** 2 + 5.638e-12 * T_ssd ** 3)
        specific_volume_Na = 1.0 / rho_Na

        # 计算 NaK 比体积
        specific_volume_NaK = self.N_Na * specific_volume_Na + self.N_K * specific_volume_K

        return specific_volume_NaK * 1.003 # 考虑 NaK78 比体积的微小变化
    
    def specific_volume_liquid_sat(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        return self.specific_volume(T, P)
    
    # ==========================================================================
    # 2. 动力粘度系数 (Dynamic Viscosity)
    # ==========================================================================
    def viscosity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：动力粘度 mu [Pa·s]
        输入：温度 T [K]，压力 P [Pa]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """
        T = np.asarray(T)

        # 1. 获取密度 (kg/m^3) 并转换为公式所需的 g/cm^3
        rho_kg = self.density(T, P)
        rho_g = rho_kg / 1000.0

        # 2. 预先计算整个温区的高、低公式独立值
        # Eq. 1.18 (T < 400°C 时的基准)
        mu_low = 0.116 * (rho_g**(1/3)) * np.exp((688.0 * rho_g) / T)

        # Eq. 1.19 (T > 400°C 时的基准)
        mu_high = 0.082 * (rho_g**(1/3)) * np.exp((979.0 * rho_g) / T)

        # 3. 定义平滑过渡权重函数
        T_threshold = 673.15  # 切换点 400°C (K)
        k = 0.1               # 过渡带宽系数，控制平滑区域的宽度

        # 计算权重系数 w，利用 tanh 使其在 T_threshold 附近平滑地从 0 变到 1
        w = 0.5 * (1.0 + np.tanh(k * (T - T_threshold)))

        # 4. 加权融合计算最终粘度 (cP)
        eta_cP = (1.0 - w) * mu_low + w * mu_high

        # 5. 单位换算：厘泊 (cP) 转换为国际标准单位 Pa·s
        mu = eta_cP * 1e-3

        # 标量兼容处理
        if mu.ndim == 0:
            return float(mu)
            
        return mu
    
    # ==========================================================================
    # 3. 声速 (Sound Speed) 和 绝热可压缩系数 (Isentropic Compressibility)
    # ==========================================================================
    def sound_velocity(self, T: Numeric) -> Numeric:
        """
        输出：声速 c [m/s]
        输入：温度 T [K]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """
        T = np.asarray(T)
        T_ssd = T - 273.15  # 转换为摄氏度

        # 液态Na声速
        c_Na = 2578.0 - 0.52 * T_ssd
        # 液态K声速
        c_K = 1922.0 - 0.54 * T_ssd

        return self.m_Na * c_Na + self.m_K * c_K
    
    def adiabatic_compressibility_liquid(self, T: Numeric) -> Numeric:
        """
        输出：绝热可压缩系数 (Isentropic Compressibility)
        输入：温度 T [K]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """
        T = np.asarray(T)
        T_ssd = T - 273.15  # 转换为摄氏度
        
        # 计算NaK密度
        rho_kg = self.density(T)
        # 计算NaK声速
        c = self.sound_velocity(T)
        return 1 / (rho_kg * c**2)
    
    # ==========================================================================
    # 5. 表面张力 (Surface Tension)
    # ==========================================================================
    def surface_tension(self, T: Numeric) -> Numeric:
        """
        输出：表面张力 σ [N/m]
        输入：温度 T [K]
        来源: 仅做占位拟合，没有合适的具体数据
        """
        T = np.asarray(T)
        T_ssd = T - 273.15  # 转换为摄氏度

        return 0.119 - 8e-5 * T
    
    # ==========================================================================
    # 6. 导热系数 (Conductivity)
    # ==========================================================================
    def conductivity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：导热系数 k [W/(m·K)]
        输入：温度 T [K]，压力 P [Pa]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """
        T = np.asarray(T)
        
        T_ssd = T - 273.15  # 转换为摄氏度

        k = 21.4 + 2.07e-2 * T_ssd - 2.2e-5 * T_ssd**2

        return k
    
    # ==========================================================================
    # 7. 比热容 (Specific Heat Capacity)
    # ==========================================================================
    def heat_capacity(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：比热容 cp [J/(kg·K)]
        输入：温度 T [K]，压力 P [Pa]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """
        T = np.asarray(T)
        T_ssd = T - 273.15  # 转换为摄氏度
        
        cp = 970.688 - 0.3690288 * T_ssd + 3.43088e-4 * T_ssd**2
        
        return cp

    # ==========================================================================
    # 8. 普朗特数 (Prandtl Number)
    # ==========================================================================
    def prandtl_number(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：普朗特数 Pr
        输入：温度 T [K]，压力 P [Pa]
        来源: 《SODIUM-NaK ENGINEERING HANDBOOK. VOLUME I. SODIUM CHEMISTRY AND PHYSICAL PROPERTIES》
        """
        mu = self.viscosity(T, P)
        k = self.conductivity(T, P)
        cp = self.heat_capacity(T, P)
        return mu * cp / k
    
    # ==========================================================================
    # 8. 焓值 (Enthalpy) 和 反求温度
    # ==========================================================================
    def enthalpy(self, T: Numeric, P: Numeric = 0.0) -> Numeric:
        """
        输出：焓值 h [J/kg/K]
        输入：温度 T [K]，压力 P [Pa]
        来源: 比热容积分而得
        """
        T = np.asarray(T)
        T_ssd = T - 273.15  # 转换为摄氏度

        h = 145.607 + 970.688 * T_ssd - 0.1845144 * T_ssd**2 + 1.1436267e-4 * T_ssd**3
        
        return h
    
    def temperature_from_enthalpy(self, h: Numeric, P: Numeric) -> Numeric:
        """
        输出：焓值 h [J/kg/K]
        输入：温度 T [K]，压力 P [Pa]
        来源: 拟合函数的反函数
        """

        # 将输入转化为 numpy 数组以支持向量化运算
        h = np.asarray(h)
        
        # 引入缩放因子，防止高阶拟合时的计算溢出或精度损失
        H = h / 1.0e5
        
        # 6阶拟合高精度系数 (已乘入原始输出中的 1.0e+02 乘子)
        c6 =  0.0000207246601
        c5 =  0.0002980337881
        c4 = -0.0149643035931
        c3 = -0.0246563914243
        c2 =  1.9697000874028
        c1 =  103.0537331204521
        c0 =  272.9895037621607

        # 采用秦九韶算法（嵌套乘法）高效计算多项式
        # 等价于: T = c6*H^6 + c5*H^5 + c4*H^4 + c3*H^3 + c2*H^2 + c1*H + c0
        T = (((((c6 * H + c5) * H + c4) * H + c3) * H + c2) * H + c1) * H + c0
        
        return T
    
    # ==========================================================================
    # 9. 饱和温度 (Saturation Temperature) 和 饱和蒸汽压 (Saturation Pressure)
    # ==========================================================================
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

    # ==========================================================================
    # 10. 汽化潜热 (Latent Heat of Vaporization)
    # ==========================================================================

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
    
    # ==========================================================================
    # 11. 偏导数 (使用 Lambda 时需注意数组兼容性)
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
    
    def liquid_density_derivative_T(self, T: Numeric) -> Numeric:
        """
        液体密度对温度偏导 rho_dT [kg/m³/K]
        由密度公式计算拟合获得
        """
        T = np.asarray(T)
        # 将 MATLAB 控制台输出的系数填入此处
        # 提示：由于密度随温度升高而降低，计算出的偏导数应全为负值
        a = 1.39777730e-08
        b = -5.69419325e-05
        c = -2.04847715e-01
        
        # 提速写法：d_rho_dT = a*T^2 + b*T + c
        d_rho_dT = (a * T + b) * T + c
        return d_rho_dT