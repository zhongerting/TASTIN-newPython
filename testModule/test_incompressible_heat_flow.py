import numpy as np
import matplotlib.pyplot as plt
import logging
import sys
import time
from typing import List, Any, Optional

# --- 引入项目组件 ---
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Solvers.Hydrodynamics.Components import (FluidChannel, FlowJunction, FluidVolume, IncompressibleFluidChannel,
                                              IncompressibleFluidVolume)
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume, InletJunction, IncompressibleBoundaryVolume
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.Couplers import FluidSolidCouple
from MathSolvers.solver_module import NuclearODESolver
from Correlations.Correlations import nu_aoki

from Solvers.HeatConduction.Boundary import BoundaryRegion

from MathSolvers.optimization_utils import FluidJacobianBlockLayout

class SystemManager:
    """
    系统管理器 (System Manager) - 核心求解调度器

    功能:
    1. 状态向量管理 (构建 y, 分发 y)
    2. 压力场求解 (代数递推 / 矩阵求解)
    3. 导数计算 (dydt)
    4. 耦合调度
    """

    def __init__(self,
                 channels: List[FluidChannel] = None,
                 junctions: List[FlowJunction] = None,
                 boundaries: List[IncompressibleBoundaryVolume] = None,
                 couples: List[FluidSolidCouple] = None,
                 components: List[Any] = None):  # components 指固体组件

        self.channels = channels if channels else []
        self.junctions = junctions if junctions else []
        self.boundaries = boundaries if boundaries else []
        self.couples = couples if couples else []
        self.components = components if components else []

        # --- 1. 展平所有控制体 (用于索引) ---
        self.volumes: List[FluidVolume] = []
        # 先加入边界 (便于处理 ghost cells)
        self.volumes.extend(self.boundaries)
        # 再加入通道内的节点
        for ch in self.channels:
            self.volumes.extend(ch.volumes)

        self.n_vol = len(self.volumes)
        self.n_junc = len(self.junctions)

        # --- 2. 识别系统类型 (不可压缩 vs 可压缩) ---
        # 检查是否包含 IncompressibleFluidVolume
        self.is_incompressible = False
        for vol in self.volumes:
            if isinstance(vol, IncompressibleFluidVolume):
                self.is_incompressible = True
                break

        if self.is_incompressible:
            print("🔧 SystemManager: Detected Incompressible Architecture. Pressure DOF removed from ODE.")

        # --- 3. 计算维度 & 建立索引 ---
        self.dim = 0

        # A. 流体部分
        if self.is_incompressible:
            # 不可压缩: y = [h...h, W...W]
            # P 不在状态向量中
            self.idx_h_start = 0
            self.dim += self.n_vol

            self.idx_W_start = self.dim
            self.dim += self.n_junc
        else:
            # 可压缩: y = [P...P, h...h, W...W]
            self.idx_P_start = 0
            self.dim += self.n_vol

            self.idx_h_start = self.dim
            self.dim += self.n_vol

            self.idx_W_start = self.dim
            self.dim += self.n_junc

        # B. 固体部分
        self.idx_solid_start = self.dim
        for comp in self.components:
            # 假设组件有 mesh.n_volumes 属性
            if hasattr(comp, 'mesh'):
                self.dim += comp.mesh.n_volumes
            elif hasattr(comp, 'T'):
                self.dim += comp.T.size

        # --- 4. 压力求解器预处理 (仅不可压缩) ---
        if self.is_incompressible:
            self._setup_pressure_solver()

    def _setup_pressure_solver(self):
        """
        [不可压缩专用] 预分析压力传播路径
        寻找压力边界(Anchor)，构建计算顺序。
        目前实现：基于边界的简单启发式搜索。
        """
        self.pressure_anchors = []
        for bc in self.boundaries:
            if isinstance(bc, IncompressibleBoundaryVolume) and getattr(bc, 'is_pressure_boundary', False):
                self.pressure_anchors.append(bc)

        if not self.pressure_anchors:
            print("⚠️ Warning: Incompressible system but no Pressure Boundary found! Pressure might drift.")

    def get_initial_state(self) -> np.ndarray:
        """构建初始状态向量 y"""
        y = np.zeros(self.dim)

        # 1. Fluid Volumes
        if self.is_incompressible:
            # 只取 h
            hs = [v.h for v in self.volumes]
            y[self.idx_h_start: self.idx_h_start + self.n_vol] = hs
        else:
            # 取 P, h
            Ps = [v.P for v in self.volumes]
            hs = [v.h for v in self.volumes]
            y[self.idx_P_start: self.idx_P_start + self.n_vol] = Ps
            y[self.idx_h_start: self.idx_h_start + self.n_vol] = hs

        # 2. Junctions W
        Ws = [j.W for j in self.junctions]
        y[self.idx_W_start: self.idx_W_start + self.n_junc] = Ws

        # 3. Solid T
        current_idx = self.idx_solid_start
        for comp in self.components:
            n = comp.T.size
            y[current_idx: current_idx + n] = comp.T.flatten()
            current_idx += n

        return y

    def update_system_state(self, y: np.ndarray):
        """将 y 向量写回物理对象 (不触发压力计算，仅赋值状态)"""
        # 1. Fluid
        if self.is_incompressible:
            hs = y[self.idx_h_start: self.idx_h_start + self.n_vol]
            for i, vol in enumerate(self.volumes):
                vol.h = hs[i]
                # P 不变 (保留上一步的值或初始值)
        else:
            Ps = y[self.idx_P_start: self.idx_P_start + self.n_vol]
            hs = y[self.idx_h_start: self.idx_h_start + self.n_vol]
            for i, vol in enumerate(self.volumes):
                vol.P = Ps[i]
                vol.h = hs[i]

        # 2. Junctions
        Ws = y[self.idx_W_start: self.idx_W_start + self.n_junc]
        for i, junc in enumerate(self.junctions):
            junc.W = Ws[i]
            # 更新流速 vel (这对于后续压力计算至关重要)
            junc.update_velocity()

        # 3. Solids
        current_idx = self.idx_solid_start
        for comp in self.components:
            n = comp.T.size
            T_vals = y[current_idx: current_idx + n]
            comp.T[:] = T_vals
            current_idx += n

            # 触发固体内部更新
            if hasattr(comp, '_update_properties'): comp._update_properties()
            if hasattr(comp, '_compute_internal_resistance'): comp._compute_internal_resistance()
            if hasattr(comp, '_update_boundaries_state'): comp._update_boundaries_state()

    def solve_pressure_distribution(self):
        """
        [核心算法] 代数压力分布求解 (Reference Point Method)
        仅在不可压缩模式下调用。
        """
        if not self.is_incompressible:
            return

        # 策略：从 Pressure Anchors 出发，向连接的 Channel 扩散
        # 注意：这需要 SystemManager 知道 Channel 的存在。

        # 1. 确保所有连接的流速已更新 (update_system_state 已做，但为了保险)
        # 并且更新边界的物性 (rho, mu)，因为压降计算依赖它们
        for bc in self.boundaries:
            bc.update_properties(bc.material)

        # 2. 遍历所有压力锚点
        for anchor in self.pressure_anchors:
            anchor_P = anchor.P

            # 寻找连接到该 Anchor 的通道
            # 情况 A: Anchor -> Junction -> Channel (Inlet)
            for junc in anchor.outlet_junctions:
                # junc.from_vol 是 anchor
                # junc.to_vol 应该是 Channel 的第一个节点

                # 找到所属的 Channel (反查)
                target_channel = self._find_channel_by_volume(junc.to_vol)
                if target_channel:
                    # 计算跨过 Inlet Junction 的压降
                    # P_node_0 = P_anchor - dP_junc
                    rho = anchor.rho
                    mu = anchor.mu
                    dp_fric = junc.calculate_friction_pressure_drop(junc.vel, rho, mu)
                    dp_form = junc.calculate_form_loss_pressure_drop(junc.vel, rho)
                    # 重力 (假设 junc.length 为高差)
                    # 需判断方向: 这里简单假设 junc 是短连接，忽略重力，或 junc.length 很小

                    sign_flow = np.sign(junc.W) if abs(junc.W) > 1e-10 else 1.0
                    p_drop = sign_flow * (dp_fric + dp_form)

                    # 得到通道入口节点的压力
                    p_channel_inlet = anchor_P - p_drop

                    # 调用通道的递推方法
                    if hasattr(target_channel, 'update_pressure_distribution_downstream'):
                        target_channel.update_pressure_distribution_downstream(p_channel_inlet)

            # 情况 B: Channel (Outlet) -> Junction -> Anchor
            for junc in anchor.inlet_junctions:
                # junc.to_vol 是 anchor
                # junc.from_vol 是 Channel 的最后一个节点

                target_channel = self._find_channel_by_volume(junc.from_vol)
                if target_channel:
                    # P_node_last = P_anchor + dP_junc
                    rho = anchor.rho  # 近似使用 anchor 物性，或 junc.from_vol
                    mu = anchor.mu
                    dp_fric = junc.calculate_friction_pressure_drop(junc.vel, rho, mu)
                    dp_form = junc.calculate_form_loss_pressure_drop(junc.vel, rho)

                    sign_flow = np.sign(junc.W) if abs(junc.W) > 1e-10 else 1.0
                    p_drop = sign_flow * (dp_fric + dp_form)

                    p_channel_outlet = anchor_P + p_drop

                    if hasattr(target_channel, 'update_pressure_distribution_upstream'):
                        target_channel.update_pressure_distribution_upstream(p_channel_outlet)

    def _find_channel_by_volume(self, vol: FluidVolume) -> Optional[FluidChannel]:
        """辅助：查找 Volume 属于哪个 Channel"""
        for ch in self.channels:
            if vol in ch.volumes:
                return ch
        return None

    def dydt_function(self, t, y):
        """
        ODE 右端项函数
        """
        # [Step 0] 清除源项
        self.clear_all_sources()

        # [Step 1] 更新状态 (赋值 h, W, T_solid)
        self.update_system_state(y)

        # [Step 2] 压力分布求解 (关键差异点)
        if self.is_incompressible:
            self.solve_pressure_distribution()

        # [Step 3] 更新流体物性 (T, rho, mu)
        # 必须在 P 更新后进行
        for vol in self.volumes:
            if isinstance(vol, BoundaryVolume):
                # 边界可能在 set_boundary_state 中已更新，但再调一次无妨
                pass
            else:
                vol.update_properties(vol.material)

        # [Step 4] 执行耦合 (计算源项)
        for couple in self.couples:
            couple.execute()

        # [Step 5] 计算导数
        dydt = np.zeros(self.dim)

        # 5.1 Fluid Volumes
        if self.is_incompressible:
            # 只取 dh/dt
            for i, vol in enumerate(self.volumes):
                # get_volume_derivatives 此时只返回 dh (float)
                dh = vol.get_volume_derivatives(vol.material)
                dydt[self.idx_h_start + i] = dh
        else:
            # 取 dP/dt, dh/dt
            for i, vol in enumerate(self.volumes):
                dP, dh = vol.get_volume_derivatives(vol.material)
                dydt[self.idx_P_start + i] = dP
                dydt[self.idx_h_start + i] = dh

        # 5.2 Junctions
        for i, junc in enumerate(self.junctions):
            mat = junc.from_vol.material if junc.from_vol.material else junc.to_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[self.idx_W_start + i] = dW

        # 5.3 Solids
        current_idx = self.idx_solid_start
        for comp in self.components:
            n = comp.T.size
            dT = comp.get_derivatives(t, comp.T)
            dydt[current_idx: current_idx + n] = dT.flatten()
            current_idx += n

        # 引入不可压缩流体修正
        if self.is_incompressible:
            for ch in self.channels:
                # 1. 收集该通道所有相关的连接 (Inlet -> Internal -> Outlet)
                channel_junctions = []
                # 找入口连接 (连接到 ch.volumes[0] 的)
                for j in self.junctions:
                    if j.to_vol == ch.volumes[0]:
                        channel_junctions.append(j)

                # 加内部连接
                channel_junctions.extend(ch.internal_junctions)

                # 找出口连接 (从 ch.volumes[-1] 连出去的)
                for j in self.junctions:
                    if j.from_vol == ch.volumes[-1]:
                        channel_junctions.append(j)

                if not channel_junctions:
                    continue

                # 2. 确定“领队” (Leader)
                # 通常选第一个连接，特别是当它是 InletJunction 时
                leader_junc = channel_junctions[0]

                # 获取领队在全局 y 向量中的索引
                # self.junctions 列表的顺序对应 y 中 dW 的顺序
                leader_idx_global = self.junctions.index(leader_junc)

                # 读取领队算出的加速度
                leader_acc = dydt[self.idx_W_start + leader_idx_global]

                # 3. 广播加速度 (Followers)
                for junc in channel_junctions[1:]:
                    follower_idx = self.junctions.index(junc)
                    # [核心操作] 强制覆盖！
                    dydt[self.idx_W_start + follower_idx] = leader_acc

        # print(t)

        if t > 5:
            print("t = ", t, " Start debug")

        return dydt

    def clear_all_sources(self):
        for vol in self.volumes:
            vol.Q_vol = 0.0
            vol.Q_wall = 0.0
            if hasattr(vol, 'implicit_coeff'): vol.implicit_coeff = 0.0

        # 如果 Channel 提供了清除接口，也调用一下（通常 Channel 是操作 volumes 的代理）
        for ch in self.channels:
            ch.clear_sources()


class SimpleFlowTest:
    """
    开式管道纯流动计算算例 (修复坐标获取方式)
    """

    def __init__(self):
        # 1. 几何与工况参数
        self.L = 0.375  # [m]
        self.D_inner = 0.03  # [m]
        self.flow_area = np.pi * (self.D_inner / 2.0) ** 2
        self.hydraulic_diam = self.D_inner

        self.n_nodes = 20  # 网格数

        # 边界条件
        self.W_in = 0.33  # [kg/s] 入口流量
        self.P_out = 116000.0  # [Pa] 出口压力 (压力锚点)
        self.T_init = 743.0  # [K]

        # 材料
        self.mat_na = Sodium()

        # 系统对象
        self.manager = None

    def build_system(self):
        print("\n🔧 Building System...")

        # --- A. 建立边界 ---
        self.inlet_bc = IncompressibleBoundaryVolume(
            name="InletBC",
            material=self.mat_na,
            P=120000.0,
            T=self.T_init
        )

        self.outlet_bc = IncompressibleBoundaryVolume(
            name="OutletBC",
            material=self.mat_na,
            P=self.P_out,
            T=self.T_init
        )
        self.outlet_bc.is_pressure_boundary = True

        # --- B. 建立通道 ---
        self.channel = IncompressibleFluidChannel(
            name="TestChannel",
            n_nodes=self.n_nodes,
            total_length=self.L,
            flow_area=self.flow_area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=self.P_out,
            initial_T=self.T_init,
            material=self.mat_na
        )

        for junc in self.channel.internal_junctions:
            junc.W = 0.1

        # --- C. 建立连接 ---
        self.junc_in = InletJunction(
            name="JuncIn",
            from_vol=self.inlet_bc,
            to_vol=self.channel.volumes[0],
            W_initial=0.1
        )
        self.junc_in.set_flow_rate(self.W_in)

        self.junc_out = FlowJunction(
            name="JuncOut",
            from_vol=self.channel.volumes[-1],
            to_vol=self.outlet_bc,
            flow_area=self.flow_area
        )
        self.junc_out.W = 0.1

        # --- D. 组装 SystemManager ---
        self.manager = SystemManager(
            channels=[self.channel],
            junctions=[self.junc_in] + self.channel.internal_junctions + [self.junc_out],
            boundaries=[self.inlet_bc, self.outlet_bc]
        )

        print("✅ System Built Successfully.")

    def run_simulation(self, t_end=20.0):
        print(f"\n🚀 Starting Simulation (0.0s -> {t_end}s)...")

        y0 = self.manager.get_initial_state()
        fun_dydt = self.manager.dydt_function

        # 使用自定义求解器 (RK45 适合非刚性流动)
        solver = NuclearODESolver(method='BDF', rtol=1e-5, atol=1e-8)

        start_time = time.time()
        res = solver.solve(fun=fun_dydt, t_span=(0.0, t_end), y0=y0) # , max_step=0.05)

        print(f"✅ Simulation Done in {time.time() - start_time:.3f}s")
        return res

    def plot_results(self, res):
        """后处理绘图"""
        if not res['success']:
            print("Skipping plots due to solver failure.")
            return

        t = res['t']
        y = res['y']

        # 1. 流量趋势
        idx_W_start = self.manager.idx_W_start
        W_in_trace = y[idx_W_start, :]

        # 2. 还原最终时刻状态 (获取 P 分布)
        y_final = y[:, -1]
        self.manager.update_system_state(y_final)
        self.manager.solve_pressure_distribution()

        # --- [修复] 手动计算坐标，不依赖 .z_coordinate ---
        dx = self.L / self.n_nodes
        # 节点中心坐标: 0.5*dx, 1.5*dx, ...
        z_coords = [(i + 0.5) * dx for i in range(self.n_nodes)]

        p_dist = [v.P for v in self.channel.volumes]

        # --- 绘图 ---
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # 图1: 流量启动
        ax1.plot(t, W_in_trace, 'b-', linewidth=2, label='Inlet Flow')
        ax1.axhline(self.W_in, color='r', linestyle='--', label='Target Flow')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Mass Flow (kg/s)')
        ax1.set_title('Flow Rate Startup')
        ax1.legend()
        ax1.grid(True)

        # 图2: 压力分布
        # 添加出口点以便绘图完整 (假设出口在 L 处)
        z_plot = z_coords + [self.L]
        p_plot = p_dist + [self.outlet_bc.P]

        ax2.plot(z_plot, np.array(p_plot) / 1000.0, 'o-', color='darkorange', label='Pressure')
        ax2.set_xlabel('Axial Position (m)')
        ax2.set_ylabel('Pressure (kPa)')
        ax2.set_title(f'Pressure Distribution at t={t[-1]:.2f}s')

        total_dp = p_dist[0] - self.outlet_bc.P
        ax2.text(0.1, 0.1, f"Total dP ≈ {total_dp:.1f} Pa",
                 transform=ax2.transAxes, bbox=dict(facecolor='white', alpha=0.9))

        ax2.grid(True)
        ax2.legend()

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    test = SimpleFlowTest()
    test.build_system()
    result = test.run_simulation(t_end=20.0)
    test.plot_results(result)
