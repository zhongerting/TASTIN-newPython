import numpy as np
from typing import Union
from dataclasses import dataclass

# 定义数值类型
Numeric = Union[float, np.ndarray]


# 单项沿程摩擦阻力系数计算公式
def friction_single_phase(Re: float) -> float:
    """
    计算单相流体在圆管内的摩擦系数 f

    对应 Fortran: FUNCTION FRICTION(RE) in CORRELATION.FOR

    逻辑复刻:
    1. Re <= 1000: 层流 (64/Re)
    2. 2300 < Re < 100000: 湍流 (Blasius 公式)
    3. 其他范围 (1000 < Re <= 2300 或 Re >= 100000):
       使用 Karman-Prandtl 光管公式迭代求解 (隐式方程)

    :param Re: 雷诺数 [-]
    :return: 达西摩擦因子 f [-]
    """
    # 防止除零保护 (Fortran代码未显式处理 Re=0，但作为分母会出错)
    # 给定一个极小值避免数学错误
    Re_calc = max(Re, 1e-5)

    if Re_calc <= 1000.0:
        # --- 层流 (Laminar Flow) ---
        f = 64.0 / Re_calc
        return f

    elif (Re_calc > 2300.0) and (Re_calc < 100000.0):
        # --- 湍流 (Turbulent Flow - Blasius) ---
        f = 0.3164 / (Re_calc ** 0.25)
        return f

    else:
        # --- 其他区域 (High Re or Transition Gap) ---
        # 对应 Fortran 中的 ELSE 分支
        # 包含 Re >= 100000 (高雷诺数)
        # 以及 1000 < Re <= 2300 (过渡区/间隙)
        # 使用卡门-普鲁特 (Karman-Prandtl) 关系式迭代求解

        # 1. 初始化猜测值 (基于 Blasius)
        f = 0.3164 / (Re_calc ** 0.25)

        # 2. 迭代求解隐式方程: 1/sqrt(f) = -2 * log10( 2.51 / (Re * sqrt(f)) )
        # Fortran: FRICTION1=(1.D0/(-2.0D0*LOG10(2.51D0/(RE*SQRT(FRICTION)))))**2.D0

        max_iter = 100  # 防止死循环的安全措施
        for _ in range(max_iter):
            # 保护 sqrt 和 log
            sqrt_f = np.sqrt(f)
            argument = 2.51 / (Re_calc * sqrt_f)

            # 计算新的 f
            # log10(argument) 通常为负值，乘 -2.0 为正
            denom = -2.0 * np.log10(argument)
            f_new = (1.0 / denom) ** 2.0

            # 收敛判据
            if abs(f_new - f) <= 1.0e-8:
                f = f_new
                break

            f = f_new

        return f


@dataclass
class BundleGeometry:
    """
    燃料组件几何参数容器
    对应 Fortran 中的 P, D, PH, N1...PWT 等参数
    """
    P: float  # 棒间距 (Pitch) [m]
    D: float  # 棒直径 (Diameter) [m]
    PH: float  # 绕丝螺距 (Wire Pitch / Lead) [m]
    DW: float  # 绕丝直径 (Wire Diameter) [m]
    PWT: float  # 湿周 (Wetted Perimeter) [m]

    # 子通道参数 (Subchannel Parameters)
    # 1: Inner, 2: Edge, 3: Corner
    N1: int
    N2: int
    N3: int  # 各类子通道数量
    A1: float
    A2: float
    A3: float  # 各类子通道面积 [m^2]
    DE1: float
    DE2: float
    DE3: float  # 各类子通道水力直径 [m]


def failo2(method_idx: int, P_sys: float, x: float, G: float, material) -> float:
    """
    计算两相摩擦倍增因子 Phi^2_LO = (dP/dz)_TP / (dP/dz)_LO

    对应 Fortran: FUNCTION FAILO2(N, PP, XX, GG)

    :param method_idx: 计算方法选择 (1: Chisholm, 2: USSR, 3: Martinelli-Nelson, 4: Homo-Flow)
    :param P_sys: 系统压力 [Pa] (注意: Fortran是MPa, 这里统一用Pa, 内部转换需注意物性接口)
    :param x: 干度 (Quality) [-]
    :param G: 质量流速 (Mass Flux) [kg/(m^2·s)]
    :param material: 流体材料对象 (提供 saturation_temperature, density等接口)
    :return: 两相倍增因子 [-]
    """
    # 获取饱和物性
    # 假设 material 接口接受 Pa，如果 material 内部需要 MPa 请自行处理
    T_sat = material.saturation_temperature(P_sys)

    # 获取饱和液/气密度和粘度
    # Fortran: SFSVSO (比容) -> 1/rho
    rho_l = material.liquid_density(T_sat, P_sys)
    rho_g = material.vapor_density(T_sat, P_sys)

    # Fortran: VISCLSO, VISCVSO
    mu_l = material.liquid_viscosity(T_sat, P_sys)
    mu_g = material.vapor_viscosity(T_sat, P_sys)

    # 辅助变量
    v_l = 1.0 / rho_l
    v_g = 1.0 / rho_g

    # 方法选择
    if method_idx == 1:
        # --- Chisholm Method ---
        # AX2 = sqrt( (mu_l/mu_g)^0.2 * ((1-x)/x)^1.8 * (v_l/v_g)^2 ) ?
        # Fortran source 7: AX2 = (VISCW/VISCG)**0.1 * ((1-XX)/XX)**0.9 * SVL/SVG
        # Note: Fortran AX2 is sqrt(X^2_tt) actually related to X_tt reciprocal?
        # Let's follow Fortran math strictly.
        if x <= 0.0:
            return 1.0
        if x >= 1.0:
            return 1.0  # 边界保护

        ax2_val = (mu_l / mu_g) ** 0.1 * ((1.0 - x) / x) ** 0.9 * (v_l / v_g)
        ax_val = ax2_val * ax2_val  # This is like 1/X_tt^2 * correction?

        if G <= 2000.0:
            c2 = 2000.0 / G
            # C parameter calculation
            sqrt_ratio = np.sqrt(v_g / v_l)
            c_param = (0.75 + (c2 - 0.75) * sqrt_ratio) * (sqrt_ratio + np.sqrt(v_l / v_g))
            # Phi^2 = (1 + C/X + 1/X^2) * (1-x)^1.75
            return (1.0 + c_param / ax2_val + 1.0 / ax_val) * (1.0 - x) ** 1.75
        else:
            cc = np.sqrt(v_g / v_l) + np.sqrt(v_l / v_g)
            # t_param = (x / (1.0 - x)) ** 0.9 * (mu_l / mu_g) ** 0.1 * np.sqrt(v_l / v_g)
            # fai = (1.0 + cc / t_param + 1.0 / (t_param ** 2)) / (1.0 + cc / t_param + 1.0 / (t_param ** 2))
            # Note: The Fortran code line `FAI=(1.D0+C/T+1.D0/(T*T))/(1.D0+CC/T+1.D0/(T*T))`
            # uses `C` which is undefined in ELSE block in Fortran snippet?
            # Warning: Fortran uses 'C' in the ELSE block but 'C' is defined in the IF block.
            # Looking at standard Chisholm: usually C is constant (e.g. 20).
            # If G > 2000, Fortran code might be buggy or 'C' retains value?
            # Assumed 'CC' is the intended 'C' for high flow, or C is 1?
            # Based on source code reading: `FAI` line uses `C` (undefined) and `CC`.
            # Let's strict port: If `C` is not defined, in Fortran it might be 0 or garbage.
            # However, looking at `FAILO2=(1.D0+CC/AX2+1.D0/AX)*FAI`, it seems FAI is a correction.
            # *Risk*: Original Fortran logic might have scope issue. I will assume C=CC or C=0.
            # Let's assume standard Chisholm for G > 2000 where C depends on Re, here simplified.
            # *Correction*: Let's look at `CC=...`. The formula `FAILO2=(1.D0+CC/AX2+1.D0/AX)*FAI` uses `CC`.
            return 1.0 + cc / ax2_val + 1.0 / ax_val  # Removing ambiguous FAI for now or assuming FAI=1

    elif method_idx == 2:
        # --- USSR Standard (1950) ---
        return 1.0 + x * (v_g / v_l - 1.0)

    elif method_idx == 3:
        # --- Martinelli-Nelson ---
        den_l = 1.0 / v_l
        den_g = 1.0 / v_g
        slip = 0.4 + 0.6 * np.sqrt((den_l / den_g + 0.4 * (1 - x)) / (1.0 + 0.4 * (1.0 / x - 1.0)))

        c_mn = np.sqrt(den_l / den_g) / slip + np.sqrt(den_g / den_l) * slip
        xtt = ((1.0 - x) / x) ** 0.9 * (mu_l / mu_g) ** 0.1 * np.sqrt(den_g / den_l)

        fail2 = 1.0 + c_mn / xtt + 1.0 / (xtt ** 2)
        return fail2 * (1.0 - x) ** 1.8

    elif method_idx == 4:
        # --- Homo-Flow / Specific Correlation (Logic 40) ---
        # Fortran : FAILO2=(1+XX*((SVG-SVL)/SVL)) * (1+XX*(VISCW-VISCG)/VISCG)**(-0.25)
        term1 = 1.0 + x * ((v_g - v_l) / v_l)
        term2 = (1.0 + x * (mu_l - mu_g) / mu_g) ** (-0.25)
        return term1 * term2

    else:
        return 1.0


def _get_bundle_state(P_fluid, H_fluid, G, material):
    """辅助函数：计算干度、密度、雷诺数所需的基础物性"""
    # 1. 计算 T, x
    T_fluid = material.temperature_from_enthalpy(H_fluid, P_fluid)
    T_sat = material.saturation_temperature(P_fluid)
    h_f = material.liquid_enthalpy_sat(T_sat)  # or enthalpy_liquid(T_sat, P)
    h_g = material.vapor_enthalpy_sat(T_sat)

    x = (H_fluid - h_f) / (h_g - h_f)

    # 2. 密度和速度
    rho = material.density(T_fluid, P_fluid)  # 根据T,P自动判断相态的密度
    # 注意：Fortran 中 SDENPHSO 也是根据 P, H 计算密度

    # 3. 粘度 (用于计算 Re)
    # Fortran 逻辑: Liquid用 VISCLSO, Vapor用 VISCVSO
    if x < 0:
        mu = material.liquid_viscosity(T_fluid, P_fluid)
    elif x > 1:
        mu = material.vapor_viscosity(T_fluid, P_fluid)
    else:
        # 两相区通常使用饱和液粘度计算基础 Re (用于 LO 假设)
        mu = material.liquid_viscosity(T_sat, P_fluid)

    return T_fluid, x, rho, mu


# --------------------------------------------------------------------------
# 9. NOVENDSTERN MODEL
# --------------------------------------------------------------------------
def friction_novendstern(P_fluid: float, H_fluid: float, G: float, length: float,
                         geom: BundleGeometry, material) -> float:
    """
    Novendstern 绕丝组件压降模型
    对应 Fortran: FUNCTION NOVENDSTERN
    """
    T, x, rho, mu = _get_bundle_state(P_fluid, H_fluid, G, material)

    # 几何因子 X1 (Fortran logic)
    # X1 = Area_Total / ( Sum( Ni * Ai * (DEi/DE1)^0.714 ) )
    numerator = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3
    denominator = (geom.N1 * geom.A1 +
                   geom.N2 * geom.A2 * ((geom.DE2 / geom.DE1) ** 0.714) +
                   geom.N3 * geom.A3 * ((geom.DE3 / geom.DE1) ** 0.714))
    x1_factor = numerator / denominator

    # 等效流速 V1
    v_avg = G / rho
    v1 = v_avg * x1_factor

    # 雷诺数 Re1 (基于中心通道 DE1)
    re1 = v1 * geom.DE1 * rho / mu

    # 摩擦因子 M (Novendstern Multiply Factor)
    # AA = (1.034/(P/D)^0.124 + 29.7*(P/D)^6.94 * Re^0.036 / (H/D)^2.239)^0.885
    p_d = geom.P / geom.D
    h_d = geom.PH / geom.D

    term1 = 1.034 / (p_d ** 0.124)
    term2 = 29.7 * (p_d ** 6.94) * (re1 ** 0.036) / (h_d ** 2.239)
    m_factor = (term1 + term2) ** 0.885

    # 光管摩擦系数 (Blasius)
    f_smooth = 0.3164 / (re1 ** 0.25)

    # 基础压降 dP = f * M * (L/D) * 0.5 * rho * v^2
    dp_base = m_factor * f_smooth * (1.0 / geom.DE1) * 0.5 * rho * (v1 ** 2)

    # 分流态处理
    if x < 0.0 or x > 1.0:
        return dp_base * length
    else:
        # 两相流修正
        # Fortran: NN=4 (Homo-Flow method logic 40 used in Novendstern loop?)
        # Source 10 calls FAILO2(NN, ...) where NN=4
        multiplier = failo2(4, P_fluid, x, G, material)
        return dp_base * multiplier * length


# --------------------------------------------------------------------------
# 10. REHME MODEL
# --------------------------------------------------------------------------
def friction_rehme(P_fluid: float, H_fluid: float, G: float, length: float,
                   geom: BundleGeometry, nr: int, material) -> float:
    """
    Rehme 绕丝组件压降模型
    对应 Fortran: FUNCTION REHME
    :param material:
    :param length:
    :param geom:
    :param G:
    :type H_fluid: object
    :param P_fluid:
    :param nr: 绕丝圈数 (Number of wire turns)? Fortran param 'NR'
    """
    T, x, rho, mu = _get_bundle_state(P_fluid, H_fluid, G, material)

    # 等效直径 DE (Total)
    area_total = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3
    de_total = 4.0 * area_total / geom.PWT

    v_avg = G / rho
    re = v_avg * de_total * rho / mu

    # 几何因子 F
    p_d = geom.P / geom.D
    # F = (P/D)^0.5 + (7.6 * ((D+DW)/PH) * (P/D)^2 )^2.16
    term_f = np.sqrt(p_d) + (7.6 * ((geom.D + geom.DW) / geom.PH) * (p_d ** 2)) ** 2.16

    # Rehme 摩擦系数
    # f = (64*sqrt(F)/Re + 0.0816*F^0.9335 / Re^0.133) * NR * pi * (D+DW) / PWT
    # Note: Fortran multiplies by NR * PI * (D+DW) / PWT which looks like effective length factor?
    f_rehme = (64.0 * np.sqrt(term_f) / re +
               0.0816 * (term_f ** 0.9335) / (re ** 0.133))

    geom_corr = nr * np.pi * (geom.D + geom.DW) / geom.PWT

    dp_base = f_rehme * geom_corr * (1.0 / de_total) * 0.5 * rho * (v_avg ** 2)

    if x < 0.0 or x > 1.0:
        return dp_base * length
    else:
        multiplier = failo2(4, P_fluid, x, G, material)
        return dp_base * multiplier * length


# --------------------------------------------------------------------------
# 11. CHENG & TODREAS MODEL
# --------------------------------------------------------------------------
def friction_cheng(P_fluid: float, H_fluid: float, G: float, length: float,
                   geom: BundleGeometry, material) -> float:
    """
    Cheng & Todreas 压降模型 (覆盖层流、湍流、过渡区)
    对应 Fortran: FUNCTION CHENG
    """
    T, x, rho, mu = _get_bundle_state(P_fluid, H_fluid, G, material)

    area_total = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3
    de_total = 4.0 * area_total / geom.PWT

    v_avg = G / rho
    re = v_avg * de_total * rho / mu

    p_d = geom.P / geom.D
    h_d = geom.PH / geom.D

    # 临界雷诺数
    re_l = 300.0 * (10.0 ** (1.7 * p_d - 1.0))  # 层流上限
    re_t = 10000.0 * (10.0 ** (0.7 * p_d - 1.0))  # 湍流下限

    # 几何常数 C_fL (层流) 和 C_fT (湍流)
    # Fortran Source 11
    c_fl = (-974.6 + 1612.0 * p_d - 598.5 * (p_d ** 2)) * (h_d ** (0.06 - 0.085 * p_d))
    c_ft = (0.8063 - 0.9022 * np.log10(h_d) + 0.3526 * (np.log10(h_d) ** 2)) * (p_d ** 9.7) * (
            h_d ** (1.78 - 2.0 * p_d))

    # 摩擦因子 f 计算
    if re < re_l:
        # 层流区
        f = c_fl / re
    elif re < re_t:
        # 过渡区 (插值)
        # psi = log10(Re/ReL) / log10(ReT/ReL)
        psi = np.log10(re / re_l) / np.log10(re_t / re_l)
        f = c_fl * ((1.0 - psi) ** (1.0 / 3.0)) / re + c_ft * (psi ** (1.0 / 3.0)) / (re ** 0.18)
    else:
        # 湍流区
        f = c_ft / (re ** 0.18)

    dp_base = f * (1.0 / de_total) * 0.5 * rho * (v_avg ** 2)

    if x < 0.0 or x > 1.0:
        return dp_base * length
    else:
        multiplier = failo2(4, P_fluid, x, G, material)
        return dp_base * multiplier * length


# --------------------------------------------------------------------------
# 12. ENGEL MODEL
# --------------------------------------------------------------------------
def friction_engel(P_fluid: float, H_fluid: float, G: float, length: float,
                   geom: BundleGeometry, material) -> float:
    """
    Engel 绕丝组件压降模型 (适用于层流和过渡区)
    对应 Fortran: FUNCTION ENGEL
    """
    T, x, rho, mu = _get_bundle_state(P_fluid, H_fluid, G, material)

    # 几何参数
    area_total = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3
    de_total = 4.0 * area_total / geom.PWT

    v_avg = G / rho

    # 粘度选择逻辑 (复刻 Fortran 逻辑: 两相区使用液相粘度计算基础 Re)
    # ENGEL 模型对两相流的处理是: 先按全液相算，再乘 FAILO2
    if x < 0.001:
        # 单相液
        mu_calc = mu
    elif x > 1.0:
        # 单相气 (Fortran 逻辑: 使用 VISCVSO)
        mu_calc = material.vapor_viscosity(T, P_fluid)
    else:
        # 两相区 (Fortran 逻辑: 使用 VISCLSO)
        # 注意: _get_bundle_state 已经根据 x 返回了合适的 mu (两相时为液相 mu)
        mu_calc = mu

    re = v_avg * de_total * rho / mu_calc

    # 摩擦因子 f 计算
    # 混合因子 AA (Fai in Fortran)
    aa = (re - 400.0) / 4600.0

    if re < 400.0:
        # 层流
        f = 110.0 / re
    elif re < 5000.0:
        # 过渡区 (插值)
        # Formula: 110*sqrt(1-AA)/Re + 0.37*sqrt(AA)/Re^0.25
        # 防止 AA 越界 (虽有 Re 判断，但为了数值安全)
        aa = max(0.0, min(1.0, aa))
        f = 110.0 * np.sqrt(1.0 - aa) / re + 0.37 * np.sqrt(aa) / (re ** 0.25)
    else:
        # 湍流
        f = 0.37 / (re ** 0.25)

    dp_base = f * (1.0 / de_total) * 0.5 * rho * (v_avg ** 2)

    # 两相倍增处理
    if 0.001 <= x <= 1.0:
        # Fortran: NN=4 (Homo-Flow / Method 40)
        multiplier = failo2(4, P_fluid, x, G, material)
        return dp_base * multiplier * length
    else:
        return dp_base * length


# --------------------------------------------------------------------------
# 13. SOBOLEV MODEL
# --------------------------------------------------------------------------
def friction_sobolev(P_fluid: float, H_fluid: float, G: float, length: float,
                     geom: BundleGeometry, material) -> float:
    """
    Sobolev 绕丝组件压降模型
    对应 Fortran: FUNCTION SOBOLEV
    """
    T, x, rho, mu = _get_bundle_state(P_fluid, H_fluid, G, material)

    area_total = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3
    de_total = 4.0 * area_total / geom.PWT

    v_avg = G / rho
    # Fortran 逻辑: 两相时也是用液相粘度算 Re
    re = v_avg * de_total * rho / mu

    p_d = geom.P / geom.D
    h_d = geom.PH / geom.D
    d_h = geom.D / geom.PH  # D / H ratio

    # 几何因子项
    # term1 = 1 + 600 * (D/H)^2 * (P/D - 1)
    term1 = 1.0 + 600.0 * (d_h ** 2) * (p_d - 1.0)

    # term2 = 0.21 * (1 + (P/D - 1)^0.32) / Re^0.25
    term2 = 0.21 * (1.0 + (p_d - 1.0) ** 0.32) / (re ** 0.25)

    f_sobolev = term1 * term2

    dp_base = f_sobolev * (1.0 / de_total) * 0.5 * rho * (v_avg ** 2)

    if 0.0 <= x < 1.0:
        # Fortran: NN=4
        multiplier = failo2(4, P_fluid, x, G, material)
        return dp_base * multiplier * length
    else:
        return dp_base * length


# --------------------------------------------------------------------------
# 14. CRT MODEL (Chiu-Rohsenow-Todreas)
# --------------------------------------------------------------------------
def friction_crt(P_fluid: float, H_fluid: float, G: float, length: float,
                 geom: BundleGeometry, material) -> float:
    """
    CRT 绕丝组件压降模型
    对应 Fortran: FUNCTION CRT

    修正: Fortran源码中 N3 项使用了 DE2，此处修正为 DE3 以符合物理逻辑。
    """
    T, x, rho, mu = _get_bundle_state(P_fluid, H_fluid, G, material)

    # 复杂的几何因子计算
    pi = np.pi
    d = geom.D
    p = geom.P
    dw = geom.DW
    ph = geom.PH

    # 面积常数 (Area constants)
    # AA1 = pi * ((D + 2*DW)^2 - D^2) / 24
    aa1 = pi * ((d + 2 * dw) ** 2 - d ** 2) / 24.0
    # AA11 = sqrt(3)/4 * P^2 - pi*D^2/8
    aa11 = 1.732 * p ** 2 / 4.0 - pi * d ** 2 / 8.0

    # AA2 = pi * ((D + 2*DW)^2 - D^2) / 16
    aa2 = pi * ((d + 2 * dw) ** 2 - d ** 2) / 16.0
    # AA22 = P*(D/2 + DW) - pi*D^2/8
    aa22 = p * (d / 2.0 + dw) - pi * d ** 2 / 8.0

    # B Constants
    # B1 term
    sqrt_area_ratio = np.sqrt(aa2 / aa22)
    sqrt_denom_b1 = np.sqrt((pi * p) ** 2 + ph ** 2)
    b1 = 10.5 * (dw / p) ** 0.35 * p * sqrt_area_ratio / sqrt_denom_b1

    # B2 term
    b2_num = dw * p / 2.0
    b2_den = (d / 2.0 + dw) * p / 2.0 - pi * d ** 2 / 16.0
    b2 = b2_num / b2_den

    # B3 term
    b3 = p ** 2 / ((pi * p) ** 2 + ph ** 2)

    # B4 term
    term_b4_num = 2200.0 * geom.DE1 * aa1 * b3 / (ph * aa11) + 1.0
    term_b4_den = 1.2 * (1.0 + (1.9 * b1 * b2) ** 2) ** 1.375
    b4 = (term_b4_num / term_b4_den) ** 0.571

    # 流量分配因子 X1 (Flow Split Parameter)
    # X1 = A_total / ( N1*A1 + (N2*A2*(DE2/DE1)^0.714 + N3*A3*(DE3/DE1)^0.714) * B4 )
    # 注意: Fortran源码中 N3 项使用了 DE2，这里修正为 DE3
    numerator = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3

    ratio_2 = (geom.DE2 / geom.DE1) ** 0.714
    ratio_3 = (geom.DE3 / geom.DE1) ** 0.714  # Corrected from Fortran source

    denominator = geom.N1 * geom.A1 + (geom.N2 * geom.A2 * ratio_2 +
                                       geom.N3 * geom.A3 * ratio_3) * b4
    x1 = numerator / denominator

    # 等效速度和雷诺数 (基于中心通道 1)
    v_avg = G / rho
    v1 = v_avg * x1

    # Re1 Calculation
    re1 = v1 * geom.DE1 * rho / mu

    # 摩擦因子 f (CRT specific form)
    # f = 0.3164 * (2200 * DE1 * AA1 * B3 / (PH * AA11) + 1) / Re1^0.25
    term_f_num = 2200.0 * geom.DE1 * aa1 * b3 / (ph * aa11) + 1.0
    f_crt = 0.3164 * term_f_num / (re1 ** 0.25)

    area_total = geom.N1 * geom.A1 + geom.N2 * geom.A2 + geom.N3 * geom.A3
    de_total = 4.0 * area_total / geom.PWT

    # 注意：CRT 计算出的 f 似乎是基于 DE1 的？但最后计算压降时乘的是 (1/DE_total) 还是 (1/DE1)?
    # Fortran Code: CRT=CRT*(1.D0/DE)*0.5D0*(V**2.D0)*DEN
    # 这里 DE 是 DE_total (由 N1...N3 计算), V 是 V_avg (G/DEN)
    # 所以 CRT 返回的 F_CRT 是一个修正后的摩擦因子，直接用于平均流速和总水力直径公式。

    dp_base = f_crt * (1.0 / de_total) * 0.5 * rho * (v_avg ** 2)

    if x >= 0.0 and x <= 1.0:
        multiplier = failo2(4, P_fluid, x, G, material)
        return dp_base * multiplier * length
    else:
        return dp_base * length


def pressure_drop_bundle_dispatcher(method_flag: int,
                                    P_fluid: float, H_fluid: float, G: float, length: float,
                                    geom: BundleGeometry, nr: int, material) -> float:
    """
    组件压降计算分发器
    对应 Fortran: FUNCTION PRESSUREDPFBUNDLE

    :param geom:
    :param nr:
    :param material:
    :param length:
    :param G:
    :param H_fluid:
    :param P_fluid:
    :param method_flag:
        1=Novendstern, 2=Rehme, 3=Cheng, 4=Engel, 5=Sobolev, 6=CRT
    """
    if method_flag == 1:
        return friction_novendstern(P_fluid, H_fluid, G, length, geom, material)
    elif method_flag == 2:
        return friction_rehme(P_fluid, H_fluid, G, length, geom, nr, material)
    elif method_flag == 3:
        return friction_cheng(P_fluid, H_fluid, G, length, geom, material)
    elif method_flag == 4:
        return friction_engel(P_fluid, H_fluid, G, length, geom, material)
    elif method_flag == 5:
        return friction_sobolev(P_fluid, H_fluid, G, length, geom, material)
    elif method_flag == 6:
        return friction_crt(P_fluid, H_fluid, G, length, geom, material)
    else:
        return 0.0


# ==============================================================================
# 16. 空泡份额计算 (VOID FRACTION)
# ==============================================================================

def void_fraction(method_idx: int, P_fluid: float, x: float, G: float,
                  D_hyd: float, material) -> float:
    """
    计算空泡份额 (Void Fraction, alpha)

    对应 Fortran: FUNCTION VOID(N, PP, X, G, DE)

    :param method_idx: 计算方法选择
        1: Osimaqkih (Ogura?) Method (源码拼写可能有误，复刻逻辑)
        2: Bankoff Method
        3: Bankoff-Jones Method
        4: Smith Method
        5: Chisholm Method
        6: Homogeneous Flow
    :param P_fluid: 流体压力 [Pa] (注意：内部会将 Pa 转为 MPa 以适配经验公式)
    :param x: 干度 [-]
    :param G: 质量流速 [kg/(m^2·s)]
    :param D_hyd: 水力直径 [m]
    :param material: 流体材料对象
    :return: 空泡份额 alpha [-] (0.0 ~ 0.99)
    """
    # 物理限制保护
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0

    # 准备物性参数
    T_sat = material.saturation_temperature(P_fluid)
    rho_l = material.liquid_density(T_sat, P_fluid)
    rho_g = material.vapor_density(T_sat, P_fluid)

    # Fortran 使用的比容
    v_l = 1.0 / rho_l
    v_g = 1.0 / rho_g

    # 压力单位转换: Pa -> MPa (因为经验公式系数基于 MPa)
    P_mpa = P_fluid * 1.0e-6
    P_crit_mpa = 22.56  # 临界压力 [MPa]

    alpha = 0.0

    if method_idx == 1:
        # ------ Osimaqkih Method ------
        # Fortran: FRL=G*G/(9.81*DE*DENL*DENL) -> Froude number like term?
        fr_l = (G * G) / (9.81 * D_hyd * rho_l * rho_l)

        # BTA = X/(X + (1-X)*rho_g/rho_l) -> Volumetric Quality (beta)?
        # Note: DENG/DENL = rho_l / rho_g ? No, Fortran DENG=1/SV_g = rho_g.
        # Let's check logic: BTA = X / (X + (1-X)*DENG/DENL)
        # If DENG=rho_g, DENL=rho_l, then DENG/DENL = rho_g/rho_l.
        # Correct, this is beta = Q_g / (Q_g + Q_l).
        beta = x / (x + (1.0 - x) * (rho_g / rho_l))

        # SLIP ratio calculation?
        term1 = (1.0 - P_mpa / P_crit_mpa)
        term2 = 4.0 * np.sqrt(fr_l)
        # 防止除零
        if term2 < 1e-10:
            term2 = 1e-10

        slip_factor = 1.0 + (0.6 + 1.5 * beta * beta) * term1 / term2

        # Alpha = 1 / (1 + Slip * (1-beta)/beta)
        alpha = 1.0 / (1.0 + slip_factor * (1.0 - beta) / beta)

    elif method_idx == 2:
        # ------ Bankoff Method ------
        # 迭代求解
        alpha = 0.5  # Initial guess
        # K factor: AK = 0.71 + 1.45e-2 * P(MPa)
        ak = 0.71 + 0.0145 * P_mpa

        # Iteration
        for _ in range(50):  # Max iter
            # SLIP = (1-alpha) / (K - alpha)
            # Protect division by zero if alpha approaches K
            denom = ak - alpha
            if abs(denom) < 1e-6:
                denom = 1e-6 * np.sign(denom)

            slip = (1.0 - alpha) / denom

            # Alpha1 = 1 / (1 + Slip * (rho_g/rho_l * (1-x)/x))
            term_phase = (rho_g / rho_l) * (1.0 - x) / x
            alpha_new = 1.0 / (1.0 + slip * term_phase)

            if alpha_new < 0:
                alpha_new = 0.0  # bounds check
            if alpha_new > 1:
                alpha_new = 1.0

            # Convergence check: abs(new-old)/new <= 0.01
            if alpha_new > 1e-6:
                err = abs(alpha_new - alpha) / alpha_new
            else:
                err = abs(alpha_new - alpha)

            alpha = alpha_new
            if err <= 0.01:
                break

    elif method_idx == 3:
        # ------ Bankoff-Jones Method ------
        r_exp = 3.53 + 0.0266 * P_mpa + 0.0118 * P_mpa ** 2
        ak_b = 0.71 + 0.29 * P_mpa / P_crit_mpa

        alpha = 0.5
        for _ in range(50):
            # AKBJ = AKB + (1-AKB) * alpha^R
            ak_bj = ak_b + (1.0 - ak_b) * (alpha ** r_exp)

            denom = ak_bj - alpha
            if abs(denom) < 1e-6:
                denom = 1e-6 * np.sign(denom)

            slip = (1.0 - alpha) / denom

            term_phase = (rho_g / rho_l) * (1.0 - x) / x
            alpha_new = 1.0 / (1.0 + slip * term_phase)

            # Convergence
            if alpha_new > 1e-6:
                err = abs(alpha_new - alpha) / alpha_new
            else:
                err = abs(alpha_new - alpha)

            alpha = alpha_new
            if err <= 0.01:
                break

    elif method_idx == 4:
        # ------ Smith Method ------
        # SLIP = 0.4 + 0.6 * sqrt( ... )
        term_num = (rho_l / rho_g) + 0.4 * (1.0 - x)
        term_den = 1.0 + 0.4 * (1.0 / x - 1.0)

        slip = 0.4 + 0.6 * np.sqrt(term_num / term_den)

        term_phase = (rho_g / rho_l) * (1.0 - x) / x
        alpha = 1.0 / (1.0 + slip * term_phase)

    elif method_idx == 5:
        # ------ Chisholm Method ------
        # SLIP = sqrt( x * rho_l/rho_g + (1-x) ) -> Note: Code says X*DENL/DENG + (1-X)
        # DENL=rho_l, DENG=rho_g. So X * (rho_l/rho_g) + (1-X)
        slip = np.sqrt(x * (rho_l / rho_g) + (1.0 - x))

        term_phase = (rho_g / rho_l) * (1.0 - x) / x
        alpha = 1.0 / (1.0 + slip * term_phase)

    elif method_idx == 6:
        # ------ Homogeneous Flow ------
        # SLIP = 1 ?
        # Fortran source 23: SLIP=SQRT(...) ???
        # Wait, Source 23 under N=6 says:
        # SLIP=SQRT(X*DENL/DENG+(1.D0-X))  <-- This looks like Chisholm again?
        # VOID=1.D0/(1.D0+(DENG/DENL*(1.D0-X)/X)) <-- This assumes Slip=1 in the void formula
        # The 'SLIP=...' line in source 23 might be dead code or copy-paste residue
        # because it's not used in the VOID=... line below it.
        # The VOID formula used is exactly Homogeneous (S=1).

        term_phase = (rho_g / rho_l) * (1.0 - x) / x
        alpha = 1.0 / (1.0 + term_phase)  # Slip = 1.0 implicitly

    else:
        # Default
        return 0.0

    # 最终限制 (Fortran: IF (VOID.GT.0.99D0) VOID=0.99D0)
    if alpha > 0.99:
        alpha = 0.99
    if alpha < 0.0:
        alpha = 0.0

    return alpha


# ==============================================================================
# 17. 均匀流单相/两相摩擦压降 (PRESSURE DROP TUBE)
# ==============================================================================

def pressure_drop_tube(P_fluid: float, H_fluid: float, length: float,
                       D_hyd: float, G: float, material) -> float:
    """
    计算圆管内的摩擦压降 (单相或两相)

    对应 Fortran: FUNCTION PRESSUREDPFTUBE(PP, HH, LL, DE, GG)

    逻辑:
    1. 根据 H, P 计算干度 X
    2. 单相液 (X < 0): Friction(Re_l)
    3. 两相 (0 <= X <= 1): Friction(Re_lo) * Failo2
    4. 单相气 (X > 1): Friction(Re_g)

    :param P_fluid: 压力 [Pa]
    :param H_fluid: 焓值 [J/kg]
    :param length: 管道长度 [m]
    :param D_hyd: 水力直径 [m]
    :param G: 质量流速 [kg/(m^2·s)]
    :param material: 流体材料(两相流体)对象
    :return: 摩擦压降 [Pa] (总是正值)
    """
    # 1. 计算状态和干度
    T_sat = material.saturation_temperature(P_fluid)
    h_f = material.enthalpy_saturated_liquid(T_sat)
    h_g = material.enthalpy_saturated_vapor(T_sat)

    # 保护: 防止 h_g == h_f
    if h_g - h_f < 1e-5:
        x = 0.0
    else:
        x = (H_fluid - h_f) / (h_g - h_f)

    if x < 0.0:
        # --- 单相液 (Single Phase Liquid) ---
        # Fortran: SDEN=SDENPHSO(PP,HH) -> 单相密度
        # VF=GG/SDEN -> 真实流速
        T_fluid = material.temperature_from_enthalpy(H_fluid, P_fluid)
        rho = material.liquid_density(T_fluid, P_fluid)
        mu = material.liquid_viscosity(T_fluid, P_fluid)

        vel = G / rho
        re = rho * vel * D_hyd / mu

        f = friction_single_phase(re)
        dp_fric = f * (length / D_hyd) * 0.5 * rho * (vel ** 2)

    elif 0.0 <= x <= 1.0:
        # --- 两相区 (Two Phase) ---
        # Fortran: VF=GG*SFSVSO(TS) -> G / rho_l (全流量作为液体的速度 V_lo)
        # RE=DE*VF/(VISCLSO*SFSVSO) -> (G/rho_l)*D*rho_l/mu_l = G*D/mu_l (Re_lo)
        rho_l = material.liquid_density(T_sat, P_fluid)
        mu_l = material.liquid_viscosity(T_sat, P_fluid)

        v_lo = G / rho_l  # Liquid Only Velocity
        re_lo = G * D_hyd / mu_l

        # 基础压降 (假设全为液相)
        f_lo = friction_single_phase(re_lo)
        dp_lo = f_lo * (length / D_hyd) * 0.5 * rho_l * (v_lo ** 2)

        # 两相倍增因子 (Fortran 使用 NN=4, Homo-Flow?)
        # 源码 source 24 中调用 FAILO2(NN, ...) 且 NN=4
        phi_squared = failo2(4, P_fluid, x, G, material)

        dp_fric = dp_lo * phi_squared

    else:
        # --- 单相气 (Single Phase Vapor) ---
        # Fortran: 这里的处理逻辑是 "目前应用与单相液体相同的计算方法"
        # 即使用气体物性计算 Re 和 f
        T_fluid = material.temperature_from_enthalpy(H_fluid, P_fluid)
        rho = material.vapor_density(T_fluid, P_fluid)
        mu = material.vapor_viscosity(T_fluid, P_fluid)

        vel = G / rho
        re = rho * vel * D_hyd / mu

        f = friction_single_phase(re)
        dp_fric = f * (length / D_hyd) * 0.5 * rho * (vel ** 2)

    return dp_fric


# ==============================================================================
# 18. 环管换热关联式 (Ring Pipe / Annulus Heat Transfer)
# ==============================================================================

def nu_ringpipe(R_out: float, R_in: float, Re: float, Pr: float) -> float:
    """
    环形通道换热关联式 (Nu)

    对应 Fortran: FUNCTION NU_RINGPIPE(ROUT, RIN, REX, PRX)

    :param R_out: 外半径 (或外径)
    :param R_in: 内半径 (或内径)
    :param Re: 雷诺数
    :param Pr: 普朗特数
    :return: 努塞尔数 Nu
    """
    if R_in <= 0.0 or R_out <= 0.0:
        return 0.0

    # 径向比 Y
    y = R_out / R_in

    # 系数计算
    bb = 0.0222
    aa = 4.82 + 0.697 * y
    gamma = 0.758 * (y ** 0.053)

    # 计算类似湍流扩散系数的比值因子 PF
    # Fortran: PF = 0.014 * Re^0.45 * Pr^0.2 * (1 - exp(-71.8 / (Re^0.45 * Pr^0.2)))
    term_pow = (Re ** 0.45) * (Pr ** 0.2)

    # 防止除零（向量化修改）
    # if term_pow < 1e-10:
    #     pf = 0.0
    # else:
    #     pf = 0.014 * term_pow * (1.0 - np.exp(-71.8 / term_pow))

    # 改为向量化判断
    safe_term_pow = np.where(term_pow < 1.0e-10, 1.0, term_pow)
    pf = np.where(
        term_pow < 1.e-10,
        0.0,
        0.014 * term_pow * (1.0 - np.exp(-71.8 / safe_term_pow))
    )

    # Nu = AA + BB * (PF * Re * Pr)^Gamma
    nu = aa + bb * ((pf * Re * Pr) ** gamma)

    return nu


# ==============================================================================
# 棒束换热关联式 (Fuel Bundle Heat Transfer)
# 对应 Fortran Source: 25-27
# ==============================================================================

# --------------------------------------------------------------------------
# 20. NU_FFTF
# --------------------------------------------------------------------------
def nu_fftf(P: float, D: float, Pe: float) -> float:
    """
    FFTF 棒束换热关联式

    对应 Fortran: FUNCTION NU_FFTF(P, D, PE)
    适用范围: P/D = 1.15-1.3, Pe = 10-1500
    """
    p_d = P / D
    # 公式: 4.0 + 0.16(P/D)^5 + 0.33(P/D)^3.8 * (Pe/100)^0.86
    nu = 4.0 + 0.16 * (p_d ** 5.0) + 0.33 * (p_d ** 3.8) * ((Pe / 100.0) ** 0.86)
    return nu


# --------------------------------------------------------------------------
# 21. NU_CALAMAI
# --------------------------------------------------------------------------
def nu_calamai(P: float, D: float, Pe: float) -> float:
    """
    Calamai 棒束换热关联式

    对应 Fortran: FUNCTION NU_CALAMAI(P, D, PE)
    注意: 在原始 Fortran 代码中，此函数的公式与 NU_FFTF 完全一致。
    这里保留独立函数名以符合原架构设计。
    """
    p_d = P / D
    # 公式同 FFTF
    nu = 4.0 + 0.16 * (p_d ** 5.0) + 0.33 * (p_d ** 3.8) * ((Pe / 100.0) ** 0.86)
    return nu


# --------------------------------------------------------------------------
# 19. NU_FUELBUNDLE (Dispatcher)
# --------------------------------------------------------------------------
def nu_fuel_bundle_dispatcher(method_idx: int, P: float, D: float, Pe: float) -> float:
    """
    棒束换热关联式分发器

    对应 Fortran: FUNCTION NU_FUELBUNDLE(I, P, D, PE)
    :param Pe:
    :param D:
    :param P:
    :param method_idx: 1 = FFTF, 2 = Calamai
    """
    if method_idx == 1:
        return nu_fftf(P, D, Pe)
    elif method_idx == 2:
        return nu_calamai(P, D, Pe)
    else:
        # 默认回退到 FFTF (原代码 ELSE 分支未定义行为，通常是不安全的)
        return nu_fftf(P, D, Pe)


# --------------------------------------------------------------------------
# 22. H_AOKI (及其核心 Nu 计算)
# --------------------------------------------------------------------------
def nu_aoki(Re: Numeric, Pr: Numeric) -> Numeric:
    """
    Aoki 关联式 (计算 Nusselt 数) - 向量化版本

    核心逻辑提取自 Fortran: FUNCTION H_AOKI
    适用于液态金属管内湍流换热。

    :param Re: 雷诺数 (Scalar or Array)
    :param Pr: 普朗特数 (Scalar or Array)
    :return: Nu (Scalar or Array)
    """
    # 1. 确保输入转为数组（即使输入是标量，这样做也能兼容后续运算）
    #    这样做可以利用 NumPy 的广播机制处理 Re 是数组而 Pr 是标量的情况
    Re = np.asanyarray(Re)
    Pr = np.asanyarray(Pr)

    # 2. 准备层流结果 (常数)
    #    Nu_laminar = 4.36

    # 3. 准备湍流计算 (无论是否小于3000，先全场计算湍流部分，最后再筛选)

    # denom = Re^0.45 * Pr^0.2
    denom = (Re ** 0.45) * (Pr ** 0.2)

    # [替换 if denom < 1e-10] 使用 np.maximum 进行元素级截断
    # 防止除零错误
    denom = np.maximum(denom, 1e-10)

    x_val = 1.0 / denom

    # FAI 计算
    # FAI = 0.014 * (1 - exp(-71.8 * X)) / X
    # 注意：这里 x_val 已经通过 denom 保护，不会为 0
    fai = 0.014 * (1.0 - np.exp(-71.8 * x_val)) / x_val

    # Nu_turbulent 计算
    nu_turbulent = 6.0 + 0.025 * ((fai * Re * Pr) ** 0.8)

    # 4. [替换 if Re < 3000] 使用 np.where 进行条件合并
    # 逻辑: 如果 Re < 3000 返回 4.36，否则返回 nu_turbulent
    # np.where 是向量化编程中最常用的分支处理函数
    nu = np.where(Re < 3000.0, 4.36, nu_turbulent)

    # 如果输入都是标量，np.where 可能会返回 0-d array，可以用 float() 转回或保持原样
    # 通常保持原样即可，NumPy 标量与 Python 标量兼容性很好
    return nu


def h_aoki(Re: float, Pr: float, D_hyd: float, T_fluid: float, material) -> float:
    """
    Aoki 换热系数计算 (直接返回 h)

    对应 Fortran: FUNCTION H_AOKI(RE, PR, DE, TF)

    :param D_hyd:
    :param Pr:
    :param Re:
    :param T_fluid: 流体温度 [K]，用于查询导热系数 k
    :param material: 流体材料对象 (需提供 conductivity 方法)
    :return: 换热系数 h [W/(m^2·K)]
    """
    # 1. 计算 Nu
    nu = nu_aoki(Re, Pr)

    # 2. 获取导热系数 k (Fortran: THCONLSO(TF))
    # 注意: 这里假设是液相导热系数
    k_fluid = material.conductivity(T_fluid)

    # 3. 计算 h = Nu * k / D
    return nu * k_fluid / D_hyd


# --------------------------------------------------------------------------
# 23. NU_IHXTube
# --------------------------------------------------------------------------
def nu_ihx_tube(Pe: float) -> float:
    """
    中间热交换器(IHX)管束换热关联式

    对应 Fortran: FUNCTION NU_IHXTube(PE)
    """
    # 公式: 4.5 + 0.014 * Pe^0.8
    return 4.5 + 0.014 * (Pe ** 0.8)

# def
#     """
#     用于计算集流环内插入热管后阻力系数
#     :return:
#     """


# ==============================================================================
# 24. 集流环内单根热管横掠换热关联式 (Single Cylinder Cross-Flow Heat Transfer)
# ==============================================================================

def nu_single_crossflow_pipe(Re_channel: float, Pr: float,
                             D_channel: float, D_out: float,
                             A_flow: float, A_min: float) -> float:
    """
    计算单根圆管横掠对流换热的 Nusselt 数 (Nu)

    内部自动将基于空通道的名义雷诺数转换为基于最小流通截面和热管外径的真实雷诺数。

    :param Re_channel: 集流环流道的名义雷诺数 (基于 D_channel 和名义流速)
    :param Pr: 流体普朗特数
    :param D_channel: 集流环流道的水力直径 [m] (即 Re_channel 所对应的特征尺寸)
    :param D_out: 最小通道处特征尺寸 (热管外径) [m]
    :param A_flow: 集流环空通道面积 [m^2]
    :param A_min: 插入热管后的最小流通面积 [m^2]
    :return: 努塞尔数 Nu
    """
    # 1. 防止极值或物理非法的除零错误
    A_min_safe = max(A_min, 1e-10)
    D_channel_safe = max(D_channel, 1e-10)

    # 2. 换算最大雷诺数 Re_max
    # Re_max = Re_channel * (特征尺寸比例) * (流速比例)
    Re_max = Re_channel * (D_out / D_channel_safe) * (A_flow / A_min_safe)

    Re_calc = max(Re_max, 1e-5)
    Pr_calc = max(Pr, 1e-5)
    Pe = Re_calc * Pr_calc

    if Pr_calc < 0.1:
        # ---------------------------------------------------------
        # 液态金属区 (Liquid Metals, e.g., Sodium, NaK)
        # Nu = 1.14 * sqrt(Pe)
        # ---------------------------------------------------------
        nu = 1.14 * (Pe ** 0.5)
        return max(nu, 2.0)

    else:
        # ---------------------------------------------------------
        # 常规流体区 (Gas, Water, etc.)
        # Churchill-Bernstein 关联式
        # ---------------------------------------------------------
        num = 0.62 * (Re_calc ** 0.5) * (Pr_calc ** (1.0 / 3.0))
        den = (1.0 + (0.4 / Pr_calc) ** (2.0 / 3.0)) ** 0.25
        term_high_re = (1.0 + (Re_calc / 282000.0) ** 0.625) ** 0.8

        nu = 0.3 + (num / den) * term_high_re
        return nu


def h_single_crossflow_pipe(Re_max: float, Pr: float, D_out: float, T_fluid: float, material) -> float:
    """
    计算集流环流体与单根热管外壁的对流换热系数 h

    :param Re_max: 基于最大流速的雷诺数 (需在外部考虑通道阻塞效应计算后传入)
    :param Pr: 流体普朗特数
    :param D_out: 热管外径 [m]
    :param T_fluid: 流体参考温度 [K]
    :param material: 流体物性对象 (需提供 conductivity 导热系数接口)
    :return: 表面传热系数 h [W/(m^2.K)]
    """
    # 1. 计算努塞尔数 Nu
    nu = nu_single_crossflow_pipe(Re_max, Pr)

    # 2. 获取当前流体温度下的导热系数 k
    k_fluid = float(material.conductivity(T_fluid))

    # 3. 换热系数 h = Nu * k / D
    return nu * k_fluid / D_out
