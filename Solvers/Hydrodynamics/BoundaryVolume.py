import numpy as np
from typing import List, Optional, Tuple, TYPE_CHECKING
from dataclasses import dataclass, field
from Correlations.Correlations import friction_single_phase
from Solvers.Hydrodynamics.Components import FluidVolume, FlowJunction, IncompressibleFluidVolume

# 防止循环导入，仅用于类型提示
if TYPE_CHECKING:
    from Materials.Base import FluidMaterial

class BoundaryVolume(FluidVolume):
    """
    边界控制体：用于设定系统的压力边界和温度/焓边界。
    状态 (P, T, h) 由用户指定，不通过微分方程求解。
    """
    """
    边界控制体：用于设定系统的压力边界和温度/焓边界。
    [修改后]: 支持瞬态变化的压力边界，采用松弛法驱动。
    """

    def __init__(self, name: str, material: 'FluidMaterial',
                 P: float, T: float, flow_area: float = 1.0,
                 hydraulic_diam: float = 0.0):
        super().__init__(name, volume=1.0, length=0.0,
                         flow_area=flow_area,
                         hydraulic_diam=hydraulic_diam,
                         initial_P=P, initial_T=T, material=material)

        self.material = material
        self.update_properties(material)

        self.T = T
        self.P = P

        # [新增] 目标压力，初始化为当前压力
        self.target_P = P

        # [新增] 目标温度 (可选，若需变温边界也可类似处理)
        self.target_T = T

    def set_state(self, P: Optional[float] = None, T: Optional[float] = None):
        """
        [修改]: 更新目标状态，而不是直接强制修改当前状态。
        """
        if P is not None:
            self.target_P = P  # 更新目标
        if T is not None:
            self.target_T = T

        # 注意：这里不再调用 update_properties，因为 P 和 T 还没变，
        # 它们会随着 dydt 的积分在后续时间步逐渐变化。

    def get_volume_derivatives(self, material: 'FluidMaterial') -> Tuple[float, float]:
        """
        [修改]: 返回驱动导数，迫使 P 追随 target_P。
        物理方程: dP/dt = K * (P_target - P_current)
        """
        # 松弛增益 K [1/s]
        # K=100.0 表示约 0.01s 的响应时间，足够快且对 BDF 求解器稳定
        gain_k = 50.0

        # 1. 压力导数
        self.dydt_P = (self.target_P - self.P) * gain_k

        # 2. 焓导数 (简化处理)
        # 如果温度恒定，焓随压力变化很小 (dh ~ vdP)。
        # 这里为了简单，我们假设焓也跟随目标焓松弛，或者直接设为 0 (如果忽略 T 变化)
        # 更严谨的做法是：计算目标焓 h_target = h(T_target, P_target)
        # 然后 d_h/dt = K * (h_target - h_current)
        if self.material:
            h_target = self.material.enthalpy(self.target_T, self.target_P)
            self.dydt_h = (h_target - self.h) * gain_k
        else:
            self.dydt_h = 0.0

        return self.dydt_P, self.dydt_h

    def update_properties(self, material: 'FluidMaterial'):
        """保持原样"""
        if self.material:
            self.h = self.material.enthalpy(self.T, self.P)
            self.rho = self.material.density(self.T, self.P)
            self.mu = self.material.viscosity(self.T, self.P)


class IncompressibleBoundaryVolume(IncompressibleFluidVolume):
    """
    不可压缩流体边界控制体 (Incompressible Boundary Volume)

    设计目的：
    1. 作为不可压缩流体回路的 '压力参考点' (Pressure Anchor)。
    2. 提供固定的或随时间变化的压力/温度边界条件。
    3. 支持两种热工模式：
       - 无限大热池 (Infinite Sink/Source): 温度恒定，不随流入流出改变。
       - 有限体积混合 (Mixing Plenum): 温度根据能量方程动态计算 (参与 ODE)。
    """

    def __init__(self,
                 name: str,
                 material: 'FluidMaterial',
                 P: float,
                 T: float,
                 flow_area: float = 1e12,
                 hydraulic_diam: float = 0.0,
                 mixing_enabled: bool = False,
                 vol_mixing: float = 1.0):
        """
        初始化边界

        :param P: 边界压力 [Pa] (系统的参考压力)
        :param T: 边界初始温度 [K]
        :param mixing_enabled: 是否启用物理混合计算?
                               False(默认) = 无限大热池, dh/dt=0
                               True = 有限体积集管, dh/dt由能量平衡决定
        :param vol_mixing: 混合体积 [m^3] (仅当 mixing_enabled=True 时有效)
        """
        # 确定节点体积：如果是无限大热池，给一个巨大的虚拟体积防止除零；如果是混合模式，使用真实体积
        volume_val = vol_mixing if mixing_enabled else 1.0e9

        # 调用 IncompressibleFluidVolume 的构造函数
        super().__init__(name,
                         volume=volume_val,
                         length=0.0,
                         flow_area=flow_area,
                         hydraulic_diam=hydraulic_diam,
                         initial_P=P,
                         initial_T=T,
                         material=material)

        self.material = material
        self.mixing_enabled = mixing_enabled

        # 外部强制变温速率 [J/kg/s]，用于实现预设的温度曲线 (Ramp/Sine等)
        self.forcing_dydt_h = 0.0

        # 标记为压力边界 (SystemManager 可据此识别从哪里开始推算压力)
        self.is_pressure_boundary = False

        self.target_P = P

        # 立即更新物性
        self.update_properties(material)

        # 默认进出口边界流动长度为0，并且水力学直径是一个很大的值
        self.length = 0.
        self.d_h = 1e12

    def get_volume_derivatives(self, material: 'FluidMaterial') -> float:
        """
        计算状态导数 (支持 混合模式 和 强制变温模式)

        Return:
            dydt_h (float): 焓随时间的变化率 [J/kg/s]
        """

        # --- 1. 物理混合计算 (Finite Volume Mixing) ---
        # 如果启用了混合模式 (例如作为闭式回路的 Plenum，或者有回流)，
        # 则调用父类 IncompressibleFluidVolume 的逻辑计算基于流量的能量平衡。
        if self.mixing_enabled:
            # 调用父类方法计算物理导数 (Flux based)
            # 此时会根据流入流出的焓差计算 dh/dt
            phys_dydt = super().get_volume_derivatives(material)
        else:
            # 无限大热池模式，物理导数为 0 (忽略流入流出的影响)
            phys_dydt = 0.0

        # --- 2. 外部强制驱动 (External Forcing) ---
        # 如果需要实现预设的变温曲线 (如 ramp, sine)，可以从外部设置 self.forcing_dydt_h
        # Total dh/dt = Physics + Forcing

        self.dydt_h = phys_dydt + self.forcing_dydt_h

        return self.dydt_h

    def set_boundary_state(self, P: float = None, T: float = None):
        """
        [瞬态控制接口] 动态设置边界状态

        在瞬态计算中 (如 dydt_wrapper)，调用此方法修改边界条件。

        :param P: 新的压力值 [Pa] (用于模拟压力激增)
        :param T: 新的温度值 [K] (直接重置温度，慎用，通常建议用 forcing_dydt_h 平滑过渡)
        """
        state_changed = False

        if P is not None:
            self.P = P
            self.target_P = P
            state_changed = True

        if T is not None:
            self.T = T
            # 如果修改了 T，需要同步更新 h
            if self.material:
                self.h = self.material.enthalpy(self.T, self.P)
            state_changed = True

        # 如果状态改变，必须立即更新物性 (rho, mu)
        # 这样在动量方程求解时，连接该边界的 Junction 才能取到正确的密度
        if state_changed and self.material:
            self.update_properties(self.material)
            # 如果是压力边界，这里的 self.P 改变会通过 SystemManager 的压力更新算法传播到整个回路

    def update_pressure_algebraic(self, new_P: float):
        """
        覆盖基类方法。

        对于边界节点，压力通常由用户设定 (set_boundary_state) 或作为固定参考点。
        但在某些特殊的环路求解逻辑中（如从出口倒推回来），可能需要临时更新它。

        默认行为：只更新 P 和物性，不覆盖用户设定的 target_P (如果有的话)。
        """
        super().update_pressure_algebraic(new_P)


class InletJunction(FlowJunction):
    """
    流量边界连接：用于强制规定流经该连接的质量流量。
    不求解动量方程。
    """

    def __init__(self, name: str, from_vol: FluidVolume, to_vol: FluidVolume,
                 W_initial: float = 0.0):
        # 调用基类初始化
        super().__init__(name, from_vol, to_vol)
        # 设定初始流量
        self.W = W_initial
        self.target_W = W_initial  # 目标流量

    def set_flow_rate(self, W: float):
        """设置目标流量 (可随时间变化)"""
        self.target_W = W

    def get_momentum_derivative(self, material: 'FluidMaterial') -> float:
        """
        覆盖基类方法。
        采用 [方案 B (松弛法)]：
        不再强制赋值 self.W，而是返回一个驱动导数，指示求解器去改变 W。

        物理方程变为: dW/dt = K * (W_target - W_current)
        """

        # [关键修改 1] 移除这行强制赋值！
        # self.W = self.target_W  <-- 删除这行！绝对不能在导数函数里修改状态变量！

        # [关键修改 2] 定义松弛增益 K (Relaxation Gain) [1/s]
        # K 代表了控制的"硬度"。
        # K = 10.0  意味着约 0.1s 的响应滞后 (比较柔和)
        # K = 100.0 意味着约 0.01s 的响应滞后 (非常灵敏，但也更刚性)
        # 推荐范围: 50.0 ~ 200.0
        gain_k = 0.1

        # [关键修改 3] 计算偏差驱动的导数
        # 如果当前流量 < 目标，导数为正 (加速)
        # 如果当前流量 > 目标，导数为负 (减速)
        self.dydt_W = (self.target_W - self.W) * gain_k

        return self.dydt_W
