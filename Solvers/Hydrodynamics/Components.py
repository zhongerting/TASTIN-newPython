import math
import numpy as np
from typing import List, Optional, Tuple, TYPE_CHECKING, Union
from dataclasses import dataclass, field
from Correlations.Correlations import friction_single_phase

# 防止循环导入，仅用于类型提示
if TYPE_CHECKING:
    from Materials.Base import FluidMaterial


# ==============================================================================
# 1. 基础物理单元：流体控制体 (FluidVolume)
# ==============================================================================

class FluidVolume:
    """
    流体控制体 (Control Volume) - 标量存储容器

    物理责任:
    1. 存储状态: P (压力), h (焓), T (温度), rho (密度)
    2. 质量守恒 -> 贡献 dP/dt
    3. 能量守恒 -> 贡献 dh/dt
    """

    def __init__(self,
                 name: str,
                 volume: float,
                 length: float,
                 flow_area: float,
                 hydraulic_diam: float,
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None):
        """
        初始化控制体
        :param name: 唯一标识名
        :param volume: 节点体积 [m^3]
        :param length: 节点长度 [m] (用于计算梯度和动量方程中的惯性长度)
        :param flow_area: 流通面积 [m^2]
        :param hydraulic_diam: 水力直径 [m] (用于Re数计算)
        :param initial_P: 初始压力 [Pa]
        :param initial_T: 初始温度 [K]
        :param material: 关联的流体物性对象 (必须在计算前绑定)
        """
        self.name = name

        # --- 几何参数 ---
        self.vol = volume
        self.len = length
        self.area = flow_area

        self.d_h = hydraulic_diam

        # --- 状态变量 (State Variables) ---
        self.P = initial_P
        # 初始焓值需要物性库计算，这里暂时给初值，必须在 initialize() 中更新
        self.h = 0.0
        self.T = initial_T
        self.rho = 0.0
        self.mu = 0.0  # 动力粘度

        # --- 导数缓存 (Derivatives) ---
        self.dydt_P = 0.0
        self.dydt_h = 0.0

        # --- 拓扑连接 (Topology) ---
        # 存储连接到本节点的 Junction 对象
        self.inlet_junctions: List['FlowJunction'] = []
        self.outlet_junctions: List['FlowJunction'] = []

        # --- 外部源项接口 ---
        # Q_wall: 来自壁面/燃料棒的热流 [W] (正值代表加热流体)
        self.Q_wall = 0.0
        # Q_volumetric: 体积热源 [W] (如伽马释热)
        self.Q_vol = 0.0

        # [新增] 隐式源项系数 (Implicit Source Coefficient) [W/K]
        # 对应方程中的 -lambda * T 项。此系数应为非负值 (lambda = h*A)。
        self.implicit_coeff = 0.0

        # 新增材料
        self.material = material

        self.initialize(self.material)

        # 添加z坐标空置值（默认为0）
        # 后续FluidChannel调用使用
        self.z_coordinate = 0.

    def initialize(self, material: 'FluidMaterial'):
        """根据初始 P, T 初始化 h, rho 等导出量"""
        self.h = float(material.enthalpy(self.T, self.P))
        self.update_properties(material)

    def update_properties(self, material: 'FluidMaterial'):
        """
        [Step A] 根据当前 P, h 更新 T, rho, mu 等物性
        必须在每个时间步开始时调用
        """
        # 1. h, P -> T
        self.T = float(material.temperature_from_enthalpy(self.h, self.P))

        # 2. T, P -> rho, mu
        self.rho = float(material.density(self.T, self.P))
        self.mu = float(material.viscosity(self.T, self.P))

    def get_volume_derivatives(self, material: 'FluidMaterial') -> Tuple[float, float]:
        """
        [Step B] 求解质量与能量守恒方程，返回 (dP/dt, dh/dt)

        方程组 (2x2):
        [ d_rho/dP    d_rho/dh ] [ dP/dt ]   [ (Sum(W_in) - Sum(W_out)) / V ]
        [    0           rho   ] [ dh/dt ] = [ (Sum(Wh_net) + Q_total) / V  ]

        注意：此处采用简化能量方程 (忽略 dP/dt 做功项)，因此矩阵是上三角的，可解耦。
        """
        # --- 1. 计算源项 (RHS) ---

        # 质量通量 [kg/s]
        sum_W_in = sum(j.W for j in self.inlet_junctions)
        sum_W_out = sum(j.W for j in self.outlet_junctions)
        mass_balance_rhs = (sum_W_in - sum_W_out) / self.vol  # d(rho)/dt+

        # 新增：设定虚拟声速（对于不可压缩流体，虚拟声速）
        # 对于热离子钠钾冷却空间堆，冷却剂流速最大约为 1.0m/s（稳态运行工况）
        # 设置虚拟声速一般要求远大于实际流速
        c_art = 200.0  # [m/s] 虚拟可调声速

        # 人工可压缩
        drho_dP_artificial = 1.0 / (c_art * c_art)

        # 能量通量 [W] -> [W/m^3]
        # 施主元胞策略 (Donor Cell):
        # 流入节点的焓 = 上游节点的焓 (如果W>0) 或 本节点焓 (如果倒流)

        flux_h_net = 0.0

        # 处理入口连接
        for j in self.inlet_junctions:
            # j.to_vol is self
            h_donor = j.from_vol.h if j.W >= 0 else self.h
            flux_h_net += j.W * h_donor

        # 处理出口连接
        for j in self.outlet_junctions:
            # j.from_vol is self
            h_donor = self.h if j.W >= 0 else j.to_vol.h
            flux_h_net -= j.W * h_donor

        source_total = self.Q_wall + self.Q_vol - self.implicit_coeff * self.T
        energy_balance_rhs = (flux_h_net + source_total) / self.vol
        # energy_balance_rhs = (flux_h_net + self.Q_wall + self.Q_vol) / self.vol

        # --- 2. 准备物性偏导数 (系数矩阵 A) ---

        # d(rho)/dP (等温压缩性)
        # 如果物性库没有直接提供，可用 finite_diff 计算，或者对于液态金属取常数
        # 这里假设 Sodium.py 实现了 liquid_density_derivative_P
        try:
            drho_dP = float(material.liquid_density_derivative_P(self.P))
        except AttributeError:
            drho_dP = 1e-7  # Fallback for slightly compressible

        # d(rho)/dh = (d(rho)/dT) * (dT/dh)
        # dT/dh 约为 1/Cp
        dT_dh = 1.0 / float(material.heat_capacity(self.T, self.P))
        drho_dT = float(material.liquid_density_derivative_T(self.T))
        drho_dh = drho_dT * dT_dh

        # --- 3. 求解线性方程组 ---
        # Eq 2 (Energy): rho * dh/dt = RHS_energy
        # => dh/dt = RHS_energy / rho
        self.dydt_h = energy_balance_rhs / max(self.rho, 1e-5)

        # Eq 1 (Mass): drho/dP * dP/dt + drho/dh * dh/dt = RHS_mass
        # => dP/dt = (RHS_mass - drho/dh * dh/dt) / drho/dP

        term_correction = drho_dh * self.dydt_h

        # 防止除零 (对于不可压缩流体，drho_dP 极小，需要人工压缩性或极小值截断)
        denom_P = max(drho_dP, 5e-4)

        # 忽略之前的 drho_dP 以及 denom_P 的计算
        # 直接使用设定的人工压缩性 drho_dP_artificial
        # denom_P = drho_dP_artificial
        denom_P = min(drho_dP_artificial, denom_P)

        self.dydt_P = (mass_balance_rhs - term_correction) / denom_P

        return self.dydt_P, self.dydt_h

    def add_coupling_source(self, explicit_part: float, implicit_factor: float):
        """
        添加耦合源项: Q = explicit_part - implicit_factor * T_fluid
        :param explicit_part: 显式部分 (如 h * A * T_wall) [W]
        :param implicit_factor: 隐式系数 (如 h * A) [W/K]
        """
        self.Q_wall += explicit_part
        self.implicit_coeff += implicit_factor


# ==============================================================================
# 1.1. 基础物理单元：不可压缩流体控制体 (IncompressibleFluidVolume)
#      用于计算不可压缩冷却剂的流动换热，将压力-流量进行解耦，
#      需要扩展一个 流网类 用于重新计算压力分布
# ==============================================================================
class IncompressibleFluidVolume(FluidVolume):
    """
    不可压缩流体控制体 (Incompressible Control Volume)

    设计目的：
    用于 '参考点法/代数压力法' 求解策略。

    关键特性：
    1. 状态向量缩减：get_volume_derivatives 仅返回焓导数 (dydt_h)，不返回压力导数。
       这意味着压力 P 不再是 ODE 求解器的状态变量。
    2. 压力代数更新：必须由外部管理器调用 update_pressure_algebraic() 来更新 P。
    3. 能量方程：采用非守恒输运形式，消除质量不平衡带来的误差。
    """

    def __init__(self,
                 name: str,
                 volume: float,
                 length: float,
                 flow_area: float,
                 hydraulic_diam: float,
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None):
        """
        初始化不可压缩控制体

        完全复用父类 FluidVolume 的初始化逻辑。
        """
        super().__init__(name, volume, length, flow_area, hydraulic_diam,
                         initial_P, initial_T, material)

        # 可以在这里添加额外的标志位，例如：
        self.is_compressible = False

        # 添加一个压力定点标志，
        self.is_pressurize = False

    def get_volume_derivatives(self, material: 'FluidMaterial') -> float:
        """
        计算能量方程导数 dh/dt

        注意：
        - 返回值仅为一个 float (dydt_h)。
        - 调用方 (SystemManager) 在构建 dydt 向量时需注意这一点，不要期待 tuple。
        """

        # ---------------------------------------------------------
        # 1. 计算净焓通量 (Net Enthalpy Flux) - 输运形式
        #    rho * V * dh/dt = Sum(W_in * (h_in - h)) + Sum(W_out * (h_out - h)) + Q
        # ---------------------------------------------------------
        net_enthalpy_flux = 0.0

        # --- 处理入口连接 ---
        for junc in self.inlet_junctions:
            w = junc.W
            if w > 0:
                # 正向流入: 携带上游焓 (h_in) 修改当前焓
                # 贡献项: W * (h_donor - h_self)
                h_in = junc.from_vol.h
                net_enthalpy_flux += w * (h_in - self.h)
            else:
                # 反向流出 (倒流):
                # 相当于流出项，且 h_flux = h_self
                # 输运项 W * (h_self - h_self) = 0，无贡献
                pass

        # --- 处理出口连接 ---
        for junc in self.outlet_junctions:
            w = junc.W
            if w > 0:
                # 正向流出:
                # 输运项 -W * (h_self - h_self) = 0，无贡献
                pass
            else:
                # 反向流入 (倒流): w < 0
                # 相当于入口，携带下游焓 (h_downstream)
                # 贡献项: -W * (h_donor - h_self)  (注意 -W 是正值)
                h_in = junc.to_vol.h
                net_enthalpy_flux -= w * (h_in - self.h)

        # ---------------------------------------------------------
        # 2. 处理源项 (Heat Source & Coupling)
        # ---------------------------------------------------------
        # Q_total = Q_explicit + Q_implicit * T
        # 兼容基类可能存在的源项属性
        q_explicit = getattr(self, 'source_explicit', 0.0)
        q_implicit = getattr(self, 'source_implicit', 0.0)

        # 注意：对于半隐式源项 Q = A + B*T，这里的 T 是当前时刻温度
        Q_total = q_explicit + q_implicit * self.T

        # ---------------------------------------------------------
        # 3. 计算导数 dh/dt
        # ---------------------------------------------------------
        # 保护除零异常 (虽然初始化后 rho/vol 不应为0)
        if self.rho > 1e-9 and self.vol > 1e-9:
            self.dydt_h = (net_enthalpy_flux + Q_total) / (self.rho * self.vol)
        else:
            self.dydt_h = 0.0

        # 只返回焓导数，这就把 P 从 ODE 系统中“物理隔离”了
        return self.dydt_h

    def update_pressure_algebraic(self, new_P: float):
        """
        代数压力更新接口

        该方法应在 dydt 计算之前，由 SystemManager 的 '压力分布更新算法' 调用。
        """
        self.P = new_P

        # [关键]: 立即更新物性
        # 虽然 P 不在 ODE 中演化，但它通过影响 rho 和 T 间接影响动量和能量方程
        if self.material:
            # 1. 更新温度 (T = f(P, h))
            self.T = self.material.temperature_from_enthalpy(self.h, self.P)

            # 2. 更新密度 (rho = f(P, T))
            # 这对重力压降和惯性项至关重要
            self.rho = self.material.density(self.T, self.P)

            # 3. 更新粘度 (影响摩擦系数)
            self.mu = self.material.viscosity(self.T, self.P)


# ==============================================================================
# 1.2. 用于封闭不可压缩网络的被动压力参考控制体
# ==============================================================================

class PressurizerVolume(IncompressibleFluidVolume):
    """
    封闭液压网络的被动绝对压力参考点。

    该控制体在压力求解后锚定绝对压力水平，
    但它不是固定压力边界，也不作为质量源/汇。
    """

    is_pressure_reference = True
    is_pressure_boundary = False

    def __init__(self,
                 name: str,
                 volume: float,
                 length: float,
                 flow_area: float,
                 hydraulic_diam: float,
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None):
        super().__init__(name, volume, length, flow_area, hydraulic_diam,
                         initial_P, initial_T, material)
        self.is_pressure_reference = True
        self.is_pressure_boundary = False
        self.target_P = float(initial_P)
        self._fixed_target_P = float(initial_P)
        self.current_time = 0.0
        self._pressure_table_times = None
        self._pressure_table_values = None

    def set_pressure(self, pressure: float):
        """设置固定的绝对压力目标值 [Pa]。"""
        self._fixed_target_P = float(pressure)
        if self._pressure_table_times is None:
            self.target_P = self._fixed_target_P

    def set_time(self, time: float):
        """更新本地时间以及当前的目标压力。"""
        self.current_time = float(time)
        self.target_P = self.compute_target_pressure(self.current_time)

    def set_pressure_table(self, times, pressures):
        """
        设置时间-绝对压力表格。

        时间单位为秒，压力单位为Pa。超出表格范围的值将通过np.interp钳位到最近的端点。
        """
        times_arr = np.asarray(times, dtype=float)
        pressures_arr = np.asarray(pressures, dtype=float)

        if times_arr.ndim != 1 or pressures_arr.ndim != 1:
            raise ValueError("稳压器压力表的时间和压力必须为一维数组。")
        if times_arr.size == 0:
            raise ValueError("稳压器压力表不能为空。")
        if times_arr.size != pressures_arr.size:
            raise ValueError("稳压器压力表的时间和压力必须具有相同的长度。")
        if not np.all(np.isfinite(times_arr)) or not np.all(np.isfinite(pressures_arr)):
            raise ValueError("稳压器压力表的值必须为有限值。")
        if np.any(np.diff(times_arr) <= 0.0):
            raise ValueError("稳压器压力表的时间必须严格递增。")

        self._pressure_table_times = times_arr.copy()
        self._pressure_table_values = pressures_arr.copy()
        self.target_P = self.compute_target_pressure(self.current_time)

    def clear_pressure_table(self):
        """Return to the fixed absolute pressure target."""
        self._pressure_table_times = None
        self._pressure_table_values = None
        self.target_P = self._fixed_target_P

    def compute_target_pressure(self, time: Optional[float] = None) -> float:
        """Return the current absolute pressure target [Pa]."""
        if self._pressure_table_times is None:
            return self._fixed_target_P

        t_eval = self.current_time if time is None else float(time)
        return float(np.interp(t_eval, self._pressure_table_times, self._pressure_table_values))


# ==============================================================================
# 1.2. Fixed-pressure pump volume for incompressible pressure networks
# ==============================================================================

class FixedPressurePumpVolume(IncompressibleFluidVolume):
    """
    Incompressible control volume with a fixed pressure rise.

    Positive delta_p raises pressure in the channel index direction
    (upstream -> pump -> downstream). The pressure rise is only consumed by
    algebraic incompressible pressure-distribution updates; the inherited
    enthalpy equation is unchanged.
    """

    def __init__(self,
                 name: str,
                 volume: float,
                 length: float,
                 flow_area: float,
                 hydraulic_diam: float,
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None,
                 delta_p: float = 0.0):
        super().__init__(name, volume, length, flow_area, hydraulic_diam,
                         initial_P, initial_T, material)
        self.delta_p = float(delta_p)

    def set_delta_p(self, delta_p: float):
        """Update the fixed pressure rise [Pa]."""
        self.delta_p = float(delta_p)

    def get_pressure_rise(self) -> float:
        """Return the fixed pressure rise [Pa]."""
        return self.delta_p


# ==============================================================================
# 2. 基础物理单元：流动连接 (FlowJunction)
# ==============================================================================

class FlowJunction:
    """
    流动连接 (Junction) - 矢量传输通道

    物理责任:
    1. 存储状态: W (质量流量), vel (流速)
    2. 动量守恒 -> 贡献 dW/dt
    """

    def __init__(self,
                 name: str,
                 from_vol: 'FluidVolume',
                 to_vol: 'FluidVolume',
                 flow_area: Optional[float] = None,
                 k_loss: float = 0.0,
                 custom_length: Optional[float] = None,
                 hydraulic_diam: Optional[float] = None,
                 dynamic_loss_params: Optional[dict] = None):
        """
        初始化连接
        :param name: 唯一标识名
        :param from_vol: 上游控制体
        :param to_vol: 下游控制体
        :param flow_area: 连接截面积 [m^2] (默认取两端较小值)
        :param k_loss: 形状阻力系数 (弯头/阀门)
        :param custom_length: 自定义惯性长度 [m] (默认取两端长度的一半之和)
        :param hydraulic_diam: 连接件自身水力直径 [m]；默认回退到两端控制体直径
        :param dynamic_loss_params: 当 k_loss < 0 时，所需的动态阻力参数配置字典
        """
        self.name = name
        self.from_vol = from_vol
        self.to_vol = to_vol

        # --- 几何参数 ---
        if flow_area is None:
            self.area = min(from_vol.area, to_vol.area)
        else:
            self.area = flow_area

        # 惯性长度 L_inertia
        if custom_length is None:
            self.length = 0.5 * (from_vol.len + to_vol.len)
        else:
            self.length = custom_length

        self.hydraulic_diameter = None if hydraulic_diam is None else float(hydraulic_diam)
        self.d_h = 0.0 if self.hydraulic_diameter is None else self.hydraulic_diameter
        self.k_loss = k_loss

        # [新增] 保存动态计算所需的参数配置
        self.dynamic_loss_params = dynamic_loss_params if dynamic_loss_params is not None else {}

        # --- 状态变量 ---
        self.W = 0.0  # [kg/s]
        self.vel = 0.0  # [m/s]
        self.dydt_W = 0.0

        # --- 自动注册拓扑 ---
        self.from_vol.outlet_junctions.append(self)
        self.to_vol.inlet_junctions.append(self)

    # 新增的函数，用于动态计算局部阻力系数
    def _compute_dynamic_k_loss(self) -> float:
        """
        动态计算局部阻力系数 (当给定 k_loss < 0 时被触发)

        当前支持模型:
        1. 'crossflow_tube_bank': 集流环内横掠管束模型
        """
        loss_model = self.dynamic_loss_params.get('model', 'none')

        if loss_model == 'single_crossflow_pipe':
            # --- 1. 获取几何参数 ---
            D_out = self.dynamic_loss_params.get('D_out', 0.015)  # 热管外径 [m]
            L_pipe = self.dynamic_loss_params.get('L_pipe', self.length)  # 热管迎风长度 [m]

            # 流道名义面积 (空管道截面积)
            A_flow = self.area
            # 热管的迎风投影面积
            A_proj = D_out * L_pipe

            # 物理保护：防止热管完全堵死管道导致除零错误 (保留 1% 极小间隙)
            if A_proj >= A_flow * 0.99:
                return 1e5  # 返回极大的阻力系数，符合物理上的“堵死

            # 计算流体实际通过的最小截面积
            A_min = A_flow - A_proj

            # --- 2. 获取上游物性与名义流速 ---
            rho = self.from_vol.rho
            mu = self.from_vol.mu
            v_nominal = abs(self.vel)

            if v_nominal < 1e-10 or rho < 1e-10 or mu < 1e-10:
                return 0.0

            # --- 3. 计算真实最大流速 V_max (基于连续性方程) ---
            # 流体在最窄处的流速
            v_max = v_nominal * (A_flow / A_min)

            # --- 4. 计算雷诺数 Re (基于管外径和最大流速) ---
            Re_D = (rho * v_max * D_out) / mu

            # --- 5. 计算单根圆管绕流阻力系数 C_D (经典工程拟合公式) ---
            if Re_D < 1.0:
                C_D = 24.0 / Re_D
            elif Re_D < 1000.0:
                C_D = 24.0 / Re_D + 3.0 / (Re_D ** 0.5) + 0.34
            elif Re_D < 2e5:
                C_D = 1.2
            else:
                C_D = 0.3

            # --- 6. 换算为 FlowJunction 适用的等效 K_eq ---
            # 真实的物理压降: dP = C_D * (A_proj / A_flow) * (0.5 * rho * V_max^2)
            # 求解器期望的公式: dP = K_eq * (0.5 * rho * V_nominal^2)
            # 因为 V_max = V_nominal * (A_flow / A_min)
            # 所以: K_eq = C_D * (A_proj / A_flow) * (A_flow / A_min)^2

            area_ratio = A_flow / A_min
            K_eq = C_D * (A_proj / A_flow) * (area_ratio ** 2)

            return K_eq

        if loss_model == 'inline_tube_bank_euler':
            eps = 1.0e-10
            D_out = float(self.dynamic_loss_params.get('D_out', 0.015))
            L_pipe = float(self.dynamic_loss_params.get('L_pipe', self.length))
            N_rows = float(self.dynamic_loss_params.get('N_rows', 1.0))

            A_flow = max(float(self.area), eps)
            A_proj = D_out * L_pipe
            if A_proj >= A_flow * 0.99:
                return 1.0e5

            A_min = A_flow - A_proj
            if A_min <= eps:
                return 1.0e5

            pitch_ratio = self.dynamic_loss_params.get('pitch_ratio', None)
            if pitch_ratio is None:
                if L_pipe <= eps:
                    return 0.0
                pitch_ratio = (A_flow / L_pipe) / D_out
            pitch_ratio = float(pitch_ratio)
            if pitch_ratio <= 1.0:
                return 1.0e5

            rho = max(float(self.from_vol.rho), eps)
            mu = max(float(self.from_vol.mu), eps)
            v_nominal = abs(float(self.vel))
            if v_nominal < eps or N_rows <= 0.0:
                return 0.0

            area_ratio = A_flow / A_min
            v_max = v_nominal * area_ratio
            Re_D = max((rho * v_max * D_out) / mu, eps)
            Eu = 0.67 * (pitch_ratio - 1.0) ** (-0.5) * Re_D ** (-0.15)
            return Eu * N_rows * (area_ratio ** 2)

        return 0.0

        # TODO: 接入实际的局部阻力计算逻辑（如我们之前讨论的 k_factor_local）
        # 目前先预留为空，默认返回 0.0

    def update_velocity(self):
        """更新流速 (用于计算阻力)"""
        # 使用施主密度 (Donor Cell Density)
        rho = self.from_vol.rho if self.W >= 0 else self.to_vol.rho
        # 避免除以零
        self.vel = self.W / (rho * self.area) if rho > 1e-9 else 0.0

    def get_mass_flow_for(self, target_vol: 'FluidVolume') -> float:
        """
        [接口] 获取针对特定控制体的质量流量。
        对于普通的 1:1 连接，无论谁问，流量都是真实的单管流量 W。
        """
        return self.W

    def calculate_friction_pressure_drop(self, velocity: float, density: float, viscosity: float) -> float:
        """
                计算摩擦压降 (Delta P_friction) 的标量值 (Magnitude)

                物理含义: dP = f * (L/D) * 0.5 * rho * v^2
                调用: Correlations.friction_single_phase(Re)
                返回: dP_fric [Pa] (总是 >= 0)
                """
        v_abs = abs(velocity)
        if v_abs < 1e-10:
            return 0.0

        # 获取特征尺寸 D_eff (水力直径)。连接件显式给定时优先使用自身水力直径，
        # 否则保持旧逻辑，按施主/相邻控制体直径兜底。
        D_eff = self.hydraulic_diameter if self.hydraulic_diameter is not None else self.from_vol.d_h

        # 如果上游是边界/大容器(D=0)，则使用下游管道的直径作为特征尺寸
        if D_eff < 1e-9:
            D_eff = self.to_vol.d_h

        # 如果还是 0 (两头都是无穷大?), 则无法计算摩擦，返回 0
        if D_eff < 1e-9:
            return 0.0

        # 1. 计算雷诺数 Re
        # 防止粘度为0或极小
        mu_safe = max(viscosity, 1e-10)
        Re = (density * v_abs * D_eff) / mu_safe

        # 2. 调用外部关联式库计算摩擦因子 f
        # 注意：这里调用的是单相摩擦系数公式
        f = friction_single_phase(Re)

        # 3. 计算压降 (Magnitude)
        # dP = f * (L/D) * 0.5 * rho * v^2
        term_geom = self.length / D_eff
        term_dynamic = 0.5 * density * (v_abs ** 2)

        dP_fric = f * term_geom * term_dynamic

        return dP_fric

    def calculate_form_loss_pressure_drop(self, velocity: float, density: float) -> float:
        """
        计算形阻压降 (局部损失) 的标量值 (Magnitude)
        dP = K * 0.5 * rho * v^2
        """
        v_abs = abs(velocity)

        # 1. 基础局部阻力 (来自固定的几何突变、弯头等)
        has_dynamic_loss = (
            self.dynamic_loss_params
            and self.dynamic_loss_params.get('model', 'none') != 'none'
        )
        total_k_loss = self.k_loss
        if has_dynamic_loss and total_k_loss < 0.0:
            total_k_loss = 0.0

        # 2. 如果配置了动态阻力模型，则计算并叠加动态阻力
        if has_dynamic_loss:
            total_k_loss += self._compute_dynamic_k_loss()

        # 3. 计算最终压降
        dP_form = total_k_loss * 0.5 * density * (v_abs ** 2)
        return dP_form

    def compute_pump_head(self, time: Optional[float] = None) -> float:
        """默认没有泵，返回 0.0"""
        # TODO: 后续添加泵模型
        return 0.0

    def get_momentum_derivative(self, material: 'FluidMaterial') -> float:
        """
        [Step C] 求解动量守恒方程，返回 dW/dt
        """
        # 1. 压差驱动
        delta_P = self.from_vol.P - self.to_vol.P

        # 2. 重力项 (预留)
        term_grav = 0.0

        # 准备物性参数 (Donor Properties)
        rho = self.from_vol.rho if self.W >= 0 else self.to_vol.rho

        # 注意：需要从 material 获取粘度 mu
        # 这里假设 FluidVolume.update_properties 已经计算了 mu 并存储，或者实时计算
        # 如果 FluidVolume 没有 mu 属性，需要在此处调用 material.viscosity
        # 建议 FluidVolume 增加 self.mu 缓存，或如下调用:
        T_donor = self.from_vol.T if self.W >= 0 else self.to_vol.T
        P_donor = self.from_vol.P if self.W >= 0 else self.to_vol.P
        mu = material.viscosity(T_donor, P_donor)  # 使用 Material 接口

        # ==========================================================
        # 3. 调用阻力接口 (Friction & Form Loss)
        # ==========================================================

        dP_fric_mag = self.calculate_friction_pressure_drop(self.vel, rho, mu)
        dP_form_mag = self.calculate_form_loss_pressure_drop(self.vel, rho)

        dP_resist_total = dP_fric_mag + dP_form_mag

        # 确定阻力方向：总是与流速方向相反
        # -1.0 * sign(v) * dP
        term_resistance = -1.0 * np.sign(self.vel) * dP_resist_total

        # ==========================================================

        # 4. 泵扬程
        term_pump = self.compute_pump_head()

        # 5. 总和求解 dW/dt
        inertia = self.length / self.area

        self.dydt_W = (delta_P + term_grav + term_pump + term_resistance) / inertia

        return self.dydt_W

# ==============================================================================
# 2.1 基础物理单元：魔法接口 (MacroFlowJunction)
# ==============================================================================

class PumpJunction(FlowJunction):
    """
    Main pump junction that adds a pressure rise to the hydraulic momentum source.

    Positive delta_p drives flow in the from_vol -> to_vol direction. The pressure
    rise affects only the hydraulic momentum equation; energy transport remains
    the inherited FlowJunction behavior.
    """

    is_pump_junction = True

    def __init__(self,
                 name: str,
                 from_vol: 'FluidVolume',
                 to_vol: 'FluidVolume',
                 flow_area: Optional[float] = None,
                 k_loss: float = 0.0,
                 custom_length: Optional[float] = None,
                 dynamic_loss_params: Optional[dict] = None,
                 delta_p: float = 0.0):
        super().__init__(
            name=name,
            from_vol=from_vol,
            to_vol=to_vol,
            flow_area=flow_area,
            k_loss=k_loss,
            custom_length=custom_length,
            dynamic_loss_params=dynamic_loss_params,
        )
        self.delta_p = float(delta_p)
        self.current_time = 0.0
        self._pressure_table_times = None
        self._pressure_table_values = None

    def set_delta_p(self, delta_p: float):
        """设置固定的泵压升 [Pa]"""
        self.delta_p = float(delta_p)

    def set_time(self, time: float):
        """更新用于压力表插值的本地时间"""
        self.current_time = float(time)

    def set_pressure_table(self, times, pressures):
        """
        设置手动时间-压升表格。

        时间单位为秒，压力单位为Pa。超出表格范围的值将通过np.interp钳位到最近的端点。
        """
        times_arr = np.asarray(times, dtype=float)
        pressures_arr = np.asarray(pressures, dtype=float)

        if times_arr.ndim != 1 or pressures_arr.ndim != 1:
            raise ValueError("泵压力表时间和压力必须为一维数组。")
        if times_arr.size == 0:
            raise ValueError("泵压力表不能为空。")
        if times_arr.size != pressures_arr.size:
            raise ValueError("泵压力表时间和压力必须具有相同的长度。")
        if not np.all(np.isfinite(times_arr)) or not np.all(np.isfinite(pressures_arr)):
            raise ValueError("泵压力表值必须为有限值。")
        if np.any(np.diff(times_arr) <= 0.0):
            raise ValueError("泵压力表时间必须严格递增。")

        self._pressure_table_times = times_arr.copy()
        self._pressure_table_values = pressures_arr.copy()

    def clear_pressure_table(self):
        """清除压力表，恢复为固定的泵压升。"""
        self._pressure_table_times = None
        self._pressure_table_values = None

    def compute_pump_head(self, time: Optional[float] = None) -> float:
        """返回当前的泵压升 [Pa]。"""
        if self._pressure_table_times is None:
            return self.delta_p

        t_eval = self.current_time if time is None else float(time)
        return float(np.interp(t_eval, self._pressure_table_times, self._pressure_table_values))


class MacroFlowJunction(FlowJunction):
    """
    宏观乘子接管 (Multiplier Junction / Area-Scaled Junction)

    物理责任:
    用于连接代表单根管子的“微观控制体 (Micro Volume)”与代表全堆/全局的“宏观控制体 (Macro Volume)”。
    在计算动量方程 (压降、流速、阻力) 时，严格保持单管 (1x) 的物理纯洁性；
    在分发质量和能量源项时，对宏观端应用乘子 (Mx)，实现广延量的非对称缩放。
    """

    def __init__(self,
                 name: str,
                 from_vol: 'FluidVolume',
                 to_vol: 'FluidVolume',
                 macro_vol: 'FluidVolume',
                 multiplier: int,
                 flow_area: Optional[float] = None,
                 k_loss: float = 0.0,
                 custom_length: Optional[float] = None,
                 dynamic_loss_params: Optional[dict] = None):
        """
        初始化宏观乘子接管

        :param macro_vol: 必须是 from_vol 或 to_vol 中的一个，指定哪一端是“被放大的宏观端” (例如 Inlet/Outlet 联箱)
        :param multiplier: 放大倍数 (代表该微观通道实际上是由多少根相同的管子并联而成的)
        """

        # 1. 正常初始化基类 (动量计算参数与普通接管完全一致)
        super().__init__(name=name,
                         from_vol=from_vol,
                         to_vol=to_vol,
                         flow_area=flow_area,
                         k_loss=k_loss,
                         custom_length=custom_length,
                         dynamic_loss_params=dynamic_loss_params)

        self.macro_vol = macro_vol
        self.multiplier = float(multiplier)

        # 2. 参数校验与乘子分配
        if macro_vol not in (from_vol, to_vol):
            raise ValueError(f"[{self.name}] macro_vol 必须是 from_vol 或 to_vol 中的一个！")

        # 预先判断并缓存供体和受体的乘子，避免在极高频的计算循环中做 if 判断
        self.multiplier_from = self.multiplier if from_vol == macro_vol else 1.0
        self.multiplier_to = self.multiplier if to_vol == macro_vol else 1.0

    def get_mass_flow_for(self, target_vol: 'FluidVolume') -> float:
        """
        [重写接口] 质量与能量守恒的魔法发生地。

        当矩阵求解器/能量方程求解器向这个接管索要流量时，
        它会根据索要者的身份，给出经过缩放的“欺骗性”流量。
        """
        if target_vol == self.from_vol:
            return self.W * self.multiplier_from
        elif target_vol == self.to_vol:
            return self.W * self.multiplier_to
        else:
            # 防御性编程：如果不相关的 Volume 来问，返回 0
            return 0.0


# ==============================================================================
# 3. 复合组件：流体通道 (FluidChannel)
# ==============================================================================

class FluidChannel:
    """
    流体通道 (Fluid Channel) - 复合组件

    物理责任:
    1. 自动化创建并连接一系列 FluidVolume (网格划分)。
    2. 提供向量化的状态访问接口 (用于耦合器和求解器)。
    3. 管理通道级别的几何参数 (如总长、入口标高)。
    """

    def __init__(self,
                 name: str,
                 n_nodes: int,
                 total_length: float,
                 flow_area: float,
                 hydraulic_diam: float,
                 active_height_start: float = 0.0,  # 用于轴向功率分布映射
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None):
        """
        初始化流体通道

        :param name: 通道名称
        :param n_nodes: 轴向网格节点数 (Control Volumes)
        :param total_length: 通道总长度 [m]
        :param flow_area: 流通面积 [m^2]
        :param hydraulic_diam: 水力直径 [m]
        :param active_height_start: 通道起点在堆芯活性区的高度坐标 [m] (便于后续物理耦合)
        :param initial_P: 初始压力 [Pa] (通常指出口或参考点压力，初始化时会推导全场)
        :param initial_T: 初始温度 [K]
        :param material: 流体材料对象
        """
        self.name = name
        self.n_nodes = n_nodes
        self.total_length = total_length
        self.area = flow_area
        self.d_h = hydraulic_diam
        self.material = material

        # 网格几何
        self.node_length = total_length / n_nodes
        self.node_volume = flow_area * self.node_length
        self.active_height_start = active_height_start

        # --- 子对象容器 ---
        self.volumes: List[FluidVolume] = []
        self.internal_junctions: List[FlowJunction] = []

        # --- 1. 创建控制体 (Volumes) ---
        for i in range(n_nodes):
            vol_name = f"{name}_Vol_{i + 1:02d}"
            vol = FluidVolume(
                name=vol_name,
                volume=self.node_volume,
                length=self.node_length,
                flow_area=flow_area,
                hydraulic_diam=hydraulic_diam,
                initial_P=initial_P,  # 简化的平坦分布，后续可用 initialize_pressure_drop 修正
                initial_T=initial_T,
                material=material
            )
            # 标记网格索引和物理坐标 (便于调试和数据映射)
            vol.mesh_index = i
            vol.z_coordinate = active_height_start + (i + 0.5) * self.node_length

            self.volumes.append(vol)

        # --- 2. 创建内部连接 (Internal Junctions) ---
        # N 个节点之间有 N-1 个内部连接
        for i in range(n_nodes - 1):
            junc_name = f"{name}_Junc_{i + 1}_{i + 2}"
            junc = FlowJunction(
                name=junc_name,
                from_vol=self.volumes[i],
                to_vol=self.volumes[i + 1],
                flow_area=flow_area,
                k_loss=0.0,  # 内部连接通常无局部形阻
                custom_length=self.node_length  # 内部连接的惯性长度 = dx/2 + dx/2 = dx
            )
            self.internal_junctions.append(junc)

        # --- 3. 外部接口 (Boundary Hooks) ---
        # 这些属性将在建立系统连接时被赋值 (连接到 Plenum 或 Pump)
        self.inlet_junction: Optional[FlowJunction] = None
        self.outlet_junction: Optional[FlowJunction] = None

        # 计算焓值

    def initialize_state(self):
        """对内部所有节点执行初始化 (h, rho, mu)"""
        if self.material is None:
            raise ValueError(f"Channel {self.name} has no material assigned!")

        for vol in self.volumes:
            vol.initialize(self.material)

        for junc in self.internal_junctions:
            junc.update_velocity()

    # ==========================================================================
    # 向量化数据访问接口 (Vectorized Accessors)
    # 这些属性对耦合器 (Coupler) 非常重要，避免了在 Python 中写显式循环
    # ==========================================================================

    @property
    def temperature_vector(self) -> np.ndarray:
        """返回所有节点的温度数组 [T1, T2, ..., Tn]"""
        return np.array([v.T for v in self.volumes])

    @property
    def pressure_vector(self) -> np.ndarray:
        """返回所有节点的压力数组"""
        return np.array([v.P for v in self.volumes])

    @property
    def density_vector(self) -> np.ndarray:
        """返回所有节点的密度数组"""
        return np.array([v.rho for v in self.volumes])

    @property
    def velocity_vector(self) -> np.ndarray:
        """
        返回所有节点的特征流速数组 (用于换热关联式计算 Re)
        注意: Volume 本身没有速度，这里取其相连 Junction 的平均值或主要流动方向值
        """
        # 简化策略: 取节点入口和出口的平均流速
        # 需注意边界节点的 inlet/outlet 可能为空
        vels = []
        for vol in self.volumes:
            v_in = sum(j.vel for j in vol.inlet_junctions) / max(len(vol.inlet_junctions), 1)
            v_out = sum(j.vel for j in vol.outlet_junctions) / max(len(vol.outlet_junctions), 1)
            vels.append((v_in + v_out) / 2.0)
        return np.array(vels)

    # ==========================================================================
    # 批量源项操作接口 (Batch Source Term Operations)
    # ==========================================================================

    def add_coupling_source_distribution(self, explicit_arr: np.ndarray, implicit_arr: np.ndarray):
        """
        向量化添加耦合源项
        """
        for i, (exp, imp) in enumerate(zip(explicit_arr, implicit_arr)):
            self.volumes[i].add_coupling_source(exp, imp)

    def add_heat_source_distribution(self, heat_dist: np.ndarray):
        """
        将空间分布的热源 [W] 施加到各个节点
        :param heat_dist: 长度为 n_nodes 的数组，单位 [W]
        """
        if len(heat_dist) != self.n_nodes:
            raise ValueError(f"Heat source length {len(heat_dist)} != nodes {self.n_nodes}")

        for i, Q in enumerate(heat_dist):
            # 累加到节点的 Q_wall 或 Q_vol (根据物理意义，这里假设是来自燃料棒的壁面热流)
            self.volumes[i].Q_wall += Q

    def clear_sources(self):
        """在每个时间步开始前清除源项"""
        for vol in self.volumes:
            vol.Q_wall = 0.0
            vol.Q_vol = 0.0
            vol.implicit_coeff = 0.0  # [新增] 清除隐式系数


# ==============================================================================
# 3.1. 复合组件：不可压缩流体通道 (IncompressibleFluidChannel)
# ==============================================================================
class IncompressibleFluidChannel(FluidChannel):
    """
    不可压缩流体通道 (Incompressible Fluid Channel)

    功能:
    1. 自动构建 IncompressibleFluidVolume 网格。
    2. 提供代数压力分布更新方法 (参考点法)，不依赖 z_coordinate 属性。
    """

    def __init__(self,
                 name: str,
                 n_nodes: int,
                 total_length: float,
                 flow_area: float,
                 hydraulic_diam: float,
                 active_height_start: float = 0.0,
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None):
        """
        初始化不可压缩通道
        :param name:
        :param n_nodes:
        :param total_length:
        :param flow_area:
        :param hydraulic_diam:
        :param active_height_start:
        :param initial_P:
        :param initial_T:
        :param material:
        """
        # 复刻 FluidChannel 的初始化逻辑，不调用 super().__init__ 以避免创建错误的 Volume 类型
        self.name = name
        self.n_nodes = n_nodes
        self.total_length = total_length
        self.area = flow_area
        self.d_h = hydraulic_diam
        self.material = material

        # 网格几何
        self.node_length = total_length / n_nodes
        self.node_volume = flow_area * self.node_length
        self.active_height_start = active_height_start

        # --- 子对象容器 ---
        self.volumes: List[IncompressibleFluidVolume] = []
        self.internal_junctions: List[FlowJunction] = []

        # --- 1. 创建控制体 (使用 IncompressibleFluidVolume) ---
        for i in range(n_nodes):
            vol_name = f"{name}_Vol_{i + 1:02d}"

            # 使用不可压缩控制体
            vol = IncompressibleFluidVolume(
                name=vol_name,
                volume=self.node_volume,
                length=self.node_length,
                flow_area=flow_area,
                hydraulic_diam=hydraulic_diam,
                initial_P=initial_P,
                initial_T=initial_T,
                material=material
            )

            # 兼容性设置: 虽然不用于计算，但保持与 FluidChannel 一致的属性赋值
            vol.mesh_index = i
            # z_coordinate 在这里赋值，但后续计算我们不依赖它，以防万一
            vol.z_coordinate = active_height_start + (i + 0.5) * self.node_length

            self.volumes.append(vol)

        # --- 2. 创建内部连接 (Internal Junctions) ---
        for i in range(n_nodes - 1):
            junc_name = f"{name}_Junc_{i + 1}_{i + 2}"
            junc = FlowJunction(
                name=junc_name,
                from_vol=self.volumes[i],
                to_vol=self.volumes[i + 1],
                flow_area=flow_area,
                k_loss=0.0,
                custom_length=self.node_length  # 连接长度 = dx
            )
            self.internal_junctions.append(junc)

        # --- 3. 外部接口 ---
        self.inlet_junction: Optional[FlowJunction] = None
        self.outlet_junction: Optional[FlowJunction] = None

    # ==========================================================================
    # 压力分布更新算法 (代数递推)
    # 修正: 使用 junc.length 计算重力压降，不调用 z_coordinate
    # ==========================================================================

    def update_pressure_distribution_downstream(self, P_inlet: float):
        """
        [正向递推]: 已知入口压力 P_inlet，推算全场压力。
        假设流道为垂直向上布置 (Flow from i to i+1 is UP)。
        """
        current_P = P_inlet

        # 1. 更新第一个节点压力 (近似认为入口界面压力传递到节点中心)
        # 注意: 严格来说 P_node = P_inlet - 0.5 * dP，这里简化处理
        self.volumes[0].update_pressure_algebraic(current_P)

        # 2. 遍历内部连接
        for i, junc in enumerate(self.internal_junctions):
            # junc 连接 volumes[i] -> volumes[i+1]
            junc.update_velocity()

            # 使用上游物性
            rho = self.volumes[i].rho
            mu = self.volumes[i].mu

            # 计算摩擦和形阻压降 (标量，恒为正)
            dp_fric = junc.calculate_friction_pressure_drop(junc.vel, rho, mu)
            dp_form = junc.calculate_form_loss_pressure_drop(junc.vel, rho)

            # 空间条件不考虑重力加速度
            g = 0

            # 计算重力压降 (假设垂直向上流，高度差为连接长度)
            dp_grav = rho * g * junc.length
            # 如果是垂直向上，h > 0，重力导致压力降低
            # dp_grav = rho * 9.81 * junc.length

            # 确定流动方向对阻力的影响
            # 如果流速为正 (i->i+1)，阻力使压力降低 (-dp)
            # 如果流速为负 (i<-i+1)，阻力使压力升高 (+dp，即 -(-dp))
            # 动量方程稳态形式: P_in - P_out = Resistance + Gravity
            # P_out = P_in - (Resistance_term + Gravity_term)

            sign_flow = np.sign(junc.W) if abs(junc.W) > 1e-10 else 1.0

            # P_{i+1} = P_i - (Sign * (Fric + Form) + Grav)
            # 注意: Gravity 始终向下，如果 junc 定义为向上(i->i+1)，则 Gravity 始终导致压降
            pressure_rise = 0.0
            next_vol = self.volumes[i + 1]
            if hasattr(next_vol, "get_pressure_rise"):
                pressure_rise = float(next_vol.get_pressure_rise())

            next_P = current_P - (sign_flow * (dp_fric + dp_form) + dp_grav) + pressure_rise

            self.volumes[i + 1].update_pressure_algebraic(next_P)
            current_P = next_P

    def update_pressure_distribution_upstream(self, P_outlet: float):
        """
        [反向递推]: 已知出口压力 P_outlet，推算全场压力。
        """
        current_P = P_outlet

        # 更新最后一个节点
        self.volumes[-1].update_pressure_algebraic(current_P)

        # 空间条件不考虑重力加速度
        g = 0

        # 逆序遍历
        # 索引 i 对应连接 volumes[i] -> volumes[i+1]
        for i in range(len(self.internal_junctions) - 1, -1, -1):
            junc = self.internal_junctions[i]
            vol_curr = self.volumes[i + 1]  # 当前已知压力的节点 (下游)
            vol_prev = self.volumes[i]  # 待求压力的节点 (上游)

            junc.update_velocity()

            # 使用当前节点(下游)物性做近似
            rho = vol_curr.rho
            mu = vol_curr.mu

            dp_fric = junc.calculate_friction_pressure_drop(junc.vel, rho, mu)
            dp_form = junc.calculate_form_loss_pressure_drop(junc.vel, rho)

            # 重力压降
            dp_grav = rho * g * junc.length

            sign_flow = np.sign(junc.W) if abs(junc.W) > 1e-10 else 1.0

            # 反向公式: P_i = P_{i+1} + (Sign * (Fric + Form) + Grav)
            pressure_rise = 0.0
            if hasattr(vol_curr, "get_pressure_rise"):
                pressure_rise = float(vol_curr.get_pressure_rise())

            prev_P = current_P + (sign_flow * (dp_fric + dp_form) + dp_grav) - pressure_rise

            vol_prev.update_pressure_algebraic(prev_P)
            current_P = prev_P

class NonUniformIncompressibleFluidChannel(IncompressibleFluidChannel):
    """
    非均匀不可压缩流体通道 (Non-Uniform Incompressible Fluid Channel)

    架构特点:
    1. 通过接收节点长度数组 (node_lengths) 实现真正的非均匀轴向网格。
    2. 对外暴露数组化属性 (self.node_length)，完美向下兼容 FluidSolidCouple 的广播机制。
    3. 对内采用严格的标量切片实例化，保护底层 0D 控制体 (FluidVolume) 的 ODE 求解器。
    4. 自动继承并无缝使用父类的代数压力分布求解算法。
    """

    def __init__(self,
                 name: str,
                 node_lengths: Union[List[float], np.ndarray],
                 flow_area: float,
                 hydraulic_diam: float,
                 active_height_start: float = 0.0,
                 initial_P: float = 1.013e5,
                 initial_T: float = 600.0,
                 material: Optional['FluidMaterial'] = None):
        """
        初始化非均匀不可压缩通道

        :param name: 通道名称
        :param node_lengths: 包含所有网格节点轴向长度的数组 [m]
        :param flow_area: 统一的流通面积 [m^2]
        :param hydraulic_diam: 统一的水力直径 [m]
        """

        # 不调用 super().__init__，因为父类期望的是标量总长，这里我们完全接管初始化逻辑
        self.name = name

        # --- 1. 处理非均匀几何数组 (鸭子类型接口兼容) ---
        # 强制转换为 numpy 数组，保障后续 Couplers.py 中的向量化运算绝对安全
        self.node_length = np.asarray(node_lengths, dtype=float)
        self.n_nodes = len(self.node_length)
        self.total_length = float(np.sum(self.node_length))

        self.area = flow_area
        self.d_h = hydraulic_diam
        self.material = material
        self.active_height_start = active_height_start

        # 暴露给 Couplers 的体积数组属性
        self.node_volume = self.area * self.node_length

        # [物理修正] 计算节点中心绝对标高 (z_coordinate)
        # 利用累加 (cumsum) 算出每个界面起点，再加本身长度的一半
        faces_start = np.insert(np.cumsum(self.node_length), 0, 0.0)[:-1]
        z_centers = self.active_height_start + faces_start + self.node_length / 2.0

        # --- 子对象容器 ---
        self.volumes: List[IncompressibleFluidVolume] = []
        self.internal_junctions: List[FlowJunction] = []

        # --- 2. 创建控制体 (降维保护) ---
        for i in range(self.n_nodes):
            vol_name = f"{name}_Vol_{i + 1:02d}"

            # 严格使用 [i] 切片传入标量，保护 FluidVolume 的内部数据结构不被数组污染
            vol = IncompressibleFluidVolume(
                name=vol_name,
                volume=self.node_volume[i],
                length=self.node_length[i],
                flow_area=self.area,
                hydraulic_diam=self.d_h,
                initial_P=initial_P,
                initial_T=initial_T,
                material=self.material
            )

            vol.mesh_index = i
            vol.z_coordinate = float(z_centers[i])

            self.volumes.append(vol)

        # --- 3. 创建内部连接 (物理修正) ---
        for i in range(self.n_nodes - 1):
            junc_name = f"{name}_Junc_{i + 1}_{i + 2}"

            # [关键物理修正] 相邻两节点的交界面特征长度 (惯性长度/重力压降长度)
            # 严格等于上游节点长度的一半 + 下游节点长度的一半
            inertia_length = (self.node_length[i] + self.node_length[i + 1]) / 2.0

            junc = FlowJunction(
                name=junc_name,
                from_vol=self.volumes[i],
                to_vol=self.volumes[i + 1],
                flow_area=self.area,
                k_loss=0.0,
                custom_length=float(inertia_length)  # 将修正后的长度传入
            )
            self.internal_junctions.append(junc)

        # --- 4. 外部接口 ---
        self.inlet_junction: Optional[FlowJunction] = None
        self.outlet_junction: Optional[FlowJunction] = None

    # ==========================================================================
    # 压力分布算法 (update_pressure_distribution) 自动继承自 IncompressibleFluidChannel！
    # 因为我们在初始化 FlowJunction 时，已经赋予了修正过的 custom_length。
    # 父类的算法在执行 dp_grav = rho * g * junc.length 时，天然就会拿到正确的物理长度。
    # ==========================================================================
