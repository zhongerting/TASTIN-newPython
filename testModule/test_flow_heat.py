import numpy as np
import matplotlib.pyplot as plt
import logging
import sys
import time
from typing import List

# --- 引入项目组件 ---
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, FluidVolume
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume, InletJunction
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.Couplers import FluidSolidCouple
from MathSolvers.solver_module import NuclearODESolver
from Correlations.Correlations import nu_aoki

from Solvers.HeatConduction.Boundary import BoundaryRegion

from MathSolvers.optimization_utils import FluidJacobianBlockLayout
USE_JAC = True


# 创建系统管理器
class SystemManager:
    def __init__(self, volumes: List[FluidVolume] = None,
                 junctions: List[FlowJunction] = None,
                 couples: List[FluidSolidCouple] = None,
                 components: List[HeatConduction2D] = None):
        self.couples = couples
        self.volumes = volumes if volumes else []
        self.junctions = junctions if junctions else []
        self.components = components if components else []

        self.n_vol = len(self.volumes)
        self.n_junc = len(self.junctions)
        self.n_comp = len(self.components)

        # --- 计算维度 ---
        # 流体部分: P, h (volumes) + W (junctions)
        self.dim_fluid = 2 * self.n_vol + self.n_junc

        # 固体部分: T (components nodes)
        self.dim_solid = 0
        for comp in self.components:
            self.dim_solid += comp.mesh.n_volumes

        # 总维度
        self.dim = self.dim_fluid + self.dim_solid

        # 自动寻找入口连接
        self.inlet_junc: InletJunction = None
        for j in self.junctions:
            if isinstance(j, InletJunction):
                self.inlet_junc = j
                break
        # 注意：如果是闭式回路可能没有入口连接，此处视具体需求抛出异常或忽略
        if self.inlet_junc is None and self.n_junc > 0:
            print("Warning: No InletJunction found in system.")
            pass

    def get_initial_state(self) -> np.ndarray:
        """收集所有组件的初始状态构建向量 y"""
        y = np.zeros(self.dim)

        # 1. Volumes: P
        y[0: self.n_vol] = [v.P for v in self.volumes]

        # 2. Volumes: h
        y[self.n_vol: 2 * self.n_vol] = [v.h for v in self.volumes]

        # 3. Junctions: W
        y[2 * self.n_vol: 2 * self.n_vol + self.n_junc] = [j.W for j in self.junctions]

        # 4. Components: T
        if self.components:
            offset = 2 * self.n_vol + self.n_junc
            for comp in self.components:
                n_nodes = comp.mesh.n_volumes
                # 扁平化赋值
                y[offset: offset + n_nodes] = comp.T
                offset += n_nodes

        return y

    def clear_all_sources(self):
        """清除流体组件的源项，防止累加"""
        if self.volumes:
            for vol in self.volumes:
                vol.Q_vol = 0.0
                vol.Q_wall = 0.0
                vol.implicit_coeff = 0.0

    def update_system_state(self, y: np.ndarray):
        """
        将状态向量 y 反向写入到物理对象中
        (用于重启动、输出结果或调试)
        """
        # --- 1. 解包流体部分 ---
        Ps = y[0: self.n_vol]
        hs = y[self.n_vol: 2 * self.n_vol]
        Ws = y[2 * self.n_vol: 2 * self.n_vol + self.n_junc]

        # 固体部分起始索引
        offset_solid = 2 * self.n_vol + self.n_junc

        # --- 2. 更新流体控制体 ---
        for i, vol in enumerate(self.volumes):
            vol.P = Ps[i]
            vol.h = hs[i]
            if isinstance(vol, BoundaryVolume):
                vol.update_properties(vol.material)
            else:
                if vol.material:
                    vol.T = vol.material.temperature_from_enthalpy(vol.h, vol.P)
                    vol.rho = vol.material.density(vol.T, vol.P)
                    vol.mu = vol.material.viscosity(vol.T, vol.P)

        # --- 3. 更新连接管 ---
        for i, junc in enumerate(self.junctions):
            junc.W = Ws[i]
            junc.update_velocity()

        # --- 4. 更新固体组件 ---
        if self.components:
            current_idx = 0
            for comp in self.components:
                n_nodes = comp.mesh.n_volumes
                # 获取当前组件的温度片段
                T_slice = y[offset_solid + current_idx: offset_solid + current_idx + n_nodes]

                # 更新组件状态
                comp.T[:] = T_slice
                comp._update_properties()

                current_idx += n_nodes

        # print(f"✅ System State Updated from vector (Dim={len(y)})")

    def dydt_function_flow_init(self, t, y):
        """
        流动初始化专用的 ODE 右端项函数
        只计算流动导数，忽略固体导热（固体温度保持恒定或随动）
        """
        # [Step 1] 解包状态并更新组件属性
        # 注意：这里我们复用 update_system_state 的逻辑手动展开，以提高效率
        Ps = y[0: self.n_vol]
        hs = y[self.n_vol: 2 * self.n_vol]
        Ws = y[2 * self.n_vol: 2 * self.n_vol + self.n_junc]

        # 更新流体
        for i, vol in enumerate(self.volumes):
            vol.P = Ps[i]
            vol.h = hs[i]
            if isinstance(vol, BoundaryVolume):
                vol.update_properties(vol.material)
            else:
                vol.T = vol.material.temperature_from_enthalpy(vol.h, vol.P)
                vol.rho = vol.material.density(vol.T, vol.P)
                vol.mu = vol.material.viscosity(vol.T, vol.P)

        # 更新连接
        for i, junc in enumerate(self.junctions):
            junc.W = Ws[i]
            junc.update_velocity()

        # [Step 2] 计算物理导数
        dydt = np.zeros(self.dim)

        # 流体 P, h 导数
        for i, vol in enumerate(self.volumes):
            dP, dh = vol.get_volume_derivatives(vol.material)
            dydt[i] = dP
            dydt[self.n_vol + i] = dh

        # 流体 W 导数
        for i, junc in enumerate(self.junctions):
            mat = junc.from_vol.material
            if mat is None:
                mat = junc.to_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[2 * self.n_vol + i] = dW

        # 固体导数保持为 0.0 (默认初始化)

        return dydt

    def dydt_function(self, t, y):
        """
        主 ODE 右端项函数 (流热耦合)
        """
        # [Step 0] 必须最先清除源项！
        self.clear_all_sources()

        # =================================================================
        # [Step 1] 状态解包与更新 (CRITICAL: 必须在 Couple Execute 之前)
        # =================================================================
        Ps = y[0: self.n_vol]
        hs = y[self.n_vol: 2 * self.n_vol]
        Ws = y[2 * self.n_vol: 2 * self.n_vol + self.n_junc]

        offset_solid = 2 * self.n_vol + self.n_junc

        # 1.1 更新流体属性 (P, h -> T, rho)
        for i, vol in enumerate(self.volumes):
            vol.P = Ps[i]
            vol.h = hs[i]
            if isinstance(vol, BoundaryVolume):
                vol.update_properties(vol.material)
            else:
                vol.T = vol.material.temperature_from_enthalpy(vol.h, vol.P)
                vol.rho = vol.material.density(vol.T, vol.P)
                vol.mu = vol.material.viscosity(vol.T, vol.P)

        # 1.2 更新连接流速
        for i, junc in enumerate(self.junctions):
            junc.W = Ws[i]
            junc.update_velocity()

        # 1.3 更新固体温度 (CRITICAL)
        # 必须先更新固体温度，并推送到边界(BoundaryRegion)，FluidSolidCouple 才能读到正确的壁温
        if self.components:
            current_idx = 0
            for comp in self.components:
                n_nodes = comp.mesh.n_volumes
                # 提取当前组件的温度切片
                T_slice = y[offset_solid + current_idx: offset_solid + current_idx + n_nodes]

                # 赋值并强制更新内部状态
                comp.T[:] = T_slice

                # [核心修正] 显式调用内部更新，确保边界条件(BoundaryRegion)获得最新温度
                # 如果不调用这些，couples.execute 读取的壁面温度将是上一步的旧值 -> 震荡
                comp._update_properties()
                comp._compute_internal_resistance()
                comp._update_boundaries_state()

                current_idx += n_nodes

        # =================================================================
        # [Step 2] 执行耦合计算 (Source Terms)
        # =================================================================
        # 此时 FluidVolume 和 Solid Boundary 都已经拥有了时刻 t 的状态
        if self.couples:
            for couple in self.couples:
                couple.execute()

        # =================================================================
        # [Step 3] 计算导数 (Derivatives)
        # =================================================================
        dydt = np.zeros(self.dim)

        # 3.1 流体 P, h
        for i, vol in enumerate(self.volumes):
            dP, dh = vol.get_volume_derivatives(vol.material)
            dydt[i] = dP
            dydt[self.n_vol + i] = dh

        # 3.2 流体 W
        for i, junc in enumerate(self.junctions):
            mat = junc.from_vol.material
            if mat is None:
                mat = junc.to_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[2 * self.n_vol + i] = dW

        # 3.3 固体 T
        if self.components:
            current_idx = 0
            for comp in self.components:
                n_nodes = comp.mesh.n_volumes

                # 计算 dT/dt
                # 注意：这里 get_derivatives 内部会再次调用 update，但这是必要的以计算 Flux
                dTdt = comp.get_derivatives(t, comp.T)

                # 正确赋值到 dydt 对应的切片位置
                target_slice = slice(offset_solid + current_idx, offset_solid + current_idx + n_nodes)
                dydt[target_slice] = dTdt

                current_idx += n_nodes

        # if t > 20:
        print("t = ", t)

        return dydt


# 建模
class HeatFlowSystem:
    def __init__(self, n_axial=20, n_radial=3):
        # =================================================================
        # 1. 几何与工况参数
        # =================================================================
        self.L = 0.375  # 长度 [m]
        self.D_inner = 0.03  # 内径 [m]
        self.R_inner = self.D_inner / 2.0
        self.thk = 0.01  # 壁厚 [m]

        self.W_in = 0.33  # 进口流量 [kg/s]
        self.P_out = 116000.0  # 出口压力 [Pa]
        self.T_init = 743.0  # 初始温度 [K]
        self.q_flux = 100.0  # 外壁热流 [W/m2] (加热器功率)

        self.n_axial = n_axial
        self.n_radial = n_radial

        # 辅助几何计算
        self.area_flow = np.pi * self.R_inner ** 2
        self.heat_area = np.pi * self.D_inner * self.L
        self.wet_perimeter = 2 * np.pi * self.R_inner
        self.heat_perimeter = self.wet_perimeter

        self.P_in_guess = self.P_out + 70.0  # 初始压力猜测

        # =================================================================
        # 2. 材料初始化
        # =================================================================
        self.mat_na = Sodium()
        self.mat_ss = AusteniticStainlessSteel()

        # =================================================================
        # 3. 流动建模 (Hydrodynamics)
        # =================================================================
        # 3.1. 进出口边界
        self.inlet = BoundaryVolume("inlet", self.mat_na, self.P_in_guess, self.T_init, self.area_flow, self.D_inner)
        self.outlet = BoundaryVolume("outlet", self.mat_na, self.P_out, self.T_init, self.area_flow, self.D_inner)

        # 3.2. 主流道
        self.channel = FluidChannel("channel", self.n_axial, self.L, self.area_flow, self.D_inner,
                                    initial_P=self.P_out, initial_T=self.T_init, material=self.mat_na)

        # 3.3. 连接 (Junctions)
        self.junc_in = InletJunction("junc_in", from_vol=self.inlet, to_vol=self.channel.volumes[0],
                                     W_initial=self.W_in)
        self.junc_in.k_loss = 0.5

        self.junc_out = FlowJunction("junc_out", from_vol=self.channel.volumes[-1], to_vol=self.outlet,
                                     flow_area=self.area_flow, k_loss=0.5)

        # 3.4. 流动初始化 (Propagate properties)
        # 确保所有组件在 t=0 时有合理的物性
        self.inlet.update_properties(self.inlet.material)
        self.outlet.update_properties(self.outlet.material)
        self.channel.initialize_state()
        for junc in self.channel.internal_junctions:
            junc.W = self.W_in
            junc.update_velocity()

        # 3.5. 收集流体对象
        self.all_vols = [self.inlet] + self.channel.volumes + [self.outlet]
        self.all_juncs = [self.junc_in] + self.channel.internal_junctions + [self.junc_out]

        # =================================================================
        # 4. 固体导热建模 (Solid Heat Conduction) - [核心修改]
        # =================================================================
        # 4.1. 建立固体网格
        # 几何类型: 圆柱 (Cylindrical)
        # X方向: 径向 (厚度 self.thk) -> n_radial
        # Y方向: 轴向 (长度 self.L)   -> n_axial
        # inner_radius: 必须设置为管道内径，这样 X=0 处就是流固交界面
        self.solid_mesh = Mesh2D(
            x_dim=self.thk,  # 径向厚度 (R_out - R_in)
            n_x=self.n_radial,
            y_dim=self.L,  # 轴向长度
            n_y=self.n_axial,
            geometry_type='cylindrical',
            inner_radius=self.R_inner  # 关键参数
        )

        # 4.2. 建立固体求解器
        self.wall_conduction = HeatConduction2D(
            mesh=self.solid_mesh,
            material=self.mat_ss,
            initial_temp=self.T_init
        )

        # 4.3. 设置固体边界条件

        # [Left Boundary] (内表面 x=0):
        # 获取引用，用于传给 FluidSolidCouple。
        # 此时它还是绝热的，Couple 会在初始化时自动添加 ResistanceBC。
        self.heat_bound_flow = self.wall_conduction.boundaries['left']

        # [Right Boundary] (外表面 x=thk):
        # 施加固定热流边界 (模拟外部电加热器)
        # 构造一个形状为 (n_axial,) 的数组，表示沿轴向均匀加热
        q_flux_array = np.full((self.n_axial,), self.q_flux)
        T_wall_out = np.full((self.n_axial,), 820.0)
        R_wall_out = np.full((self.n_axial,), 1e-5)

        # 根据 Boundary.py，方法名为 add_flux_condition
        # self.wall_conduction.boundaries['right'].add_flux_condition(q_flux_array)
        self.wall_conduction.boundaries['right'].add_resistance_condition(T_wall_out, R_wall_out)

        # 4.4. 注册组件列表 (用于 SystemManager)
        self.components = [self.wall_conduction]

        # =================================================================
        # 5. 流固耦合 (Coupling)
        # =================================================================
        def correlation_robust(Re, Pr, ratio):
            Nu = nu_aoki(Re, Pr)
            # 限制 Nu 范围以保证数值稳定性 (特别是在流速极低或刚启动时)
            return np.clip(Nu, 4.36, 50.0)

        self.couple1 = FluidSolidCouple(
            name="couple1",
            fluid=self.channel,
            solid_boundary_region=self.heat_bound_flow,  # 传入固体的内边界对象
            heated_perimeter=self.heat_perimeter,
            correlation_func=correlation_robust
        )

        self.couples = [self.couple1]

        # 占位变量
        self.res_flow_init = None
        self.res_init = None
        self.manager = None
        self.y0 = None

    def set_state(self):
        """初始化系统管理器并构建初始状态向量 y0"""
        self.manager = SystemManager(
            volumes=self.all_vols,
            junctions=self.all_juncs,
            couples=self.couples,
            components=self.components  # [关键] 必须传入固体组件
        )
        self.y0 = self.manager.get_initial_state()

    def flow_init(self):
        """
        流动初始化：先粗算后精算
        此时固体温度保持初始值不变，流体场适应几何和边界
        """
        print("🌊 Starting Flow Initialization...")
        # 1. 粗算 (大容差，快速收敛)
        solver_flow_init = NuclearODESolver(method='BDF', rtol=1e-1, atol=1e-1)
        self.res_flow_init = solver_flow_init.solve(
            self.manager.dydt_function_flow_init,
            t_span=(0.0, 10.0),
            y0=self.y0
        )

        if self.res_flow_init['success']:
            self.y0 = self.res_flow_init['y'][:, -1]
            self.manager.update_system_state(self.y0)
            print("✅ Coarse flow init successful.")
        else:
            print("❌ Coarse flow init failed.")
            return

        # 2. 精算 (小容差，消除残差震荡)
        # 推荐开启，特别是对于耦合计算，初值越平滑越不容易炸
        print("🌊 Refining Flow State...")
        solver_flow_init_fined = NuclearODESolver(method='BDF', rtol=1e-4, atol=1e-6)
        res_fined = solver_flow_init_fined.solve(
            self.manager.dydt_function_flow_init,
            t_span=(0.0, 10.0),
            y0=self.y0
        )
        if res_fined['success']:
            self.y0 = res_fined['y'][:, -1]
            self.manager.update_system_state(self.y0)
            print("✅ Fine flow init successful.")

    def save_restart_file(self, filename: str):
        """保存当前状态向量 y0 到文件"""
        if self.y0 is None:
            print("⚠️ No state to save!")
            return
        try:
            if not filename.endswith('.npy'):
                filename += '.npy'
            np.save(filename, self.y0)
            print(f"💾 Checkpoint saved to: {filename}")
        except Exception as e:
            print(f"❌ Save failed: {e}")

    def load_restart_file(self, filename: str):
        """从文件加载状态向量，并更新物理对象"""
        try:
            if not filename.endswith('.npy'):
                filename += '.npy'
            y_loaded = np.load(filename)

            if len(y_loaded) != self.manager.dim:
                print(f"❌ Dimension mismatch! File: {len(y_loaded)}, System: {self.manager.dim}")
                return

            self.y0 = y_loaded
            self.manager.update_system_state(self.y0)
            print(f"📂 Checkpoint loaded from: {filename}")

        except FileNotFoundError:
            print(f"❌ File not found: {filename}")
        except Exception as e:
            print(f"❌ Load failed: {e}")

    def calculate(self, step, starttime, endtime, rtol, atol, max_step=1e20, method='BDF'):
        """主计算循环"""
        # 针对 t=0 时刻的剧烈瞬态，强制第一步使用极小步长
        current_first_step = 1e-6 if starttime == 0.0 else None

        solver_init = NuclearODESolver(method=method, rtol=rtol, atol=atol)

        print(f"🚀 Calculation Start: {starttime}s -> {endtime}s")

        # 运行求解
        self.res_init = solver_init.solve(
            self.manager.dydt_function,
            t_span=(starttime, endtime),
            y0=self.y0,
            max_step=max_step,
            first_step=current_first_step
        )

        # 更新状态
        if self.res_init['success'] and self.res_init['y'].shape[1] > 0:
            y_final = self.res_init['y'][:, -1]
            # 更新 y0，以便连续计算
            self.y0 = y_final
            # 将数值结果反向注入回物理对象 (P, T, W 等)
            self.manager.update_system_state(y_final)
            print(f"💾 Internal state updated to time t={self.res_init['t'][-1]:.4f}s")
        else:
            print(f"❌ Calculation Failed or No Steps: {self.res_init.get('message', 'Unknown')}")

        return self.res_init

# 纯流动+雅可比计算
class HydroOnlyHeatCase:
    def __init__(self, n_axial=20):
        # =================================================================
        # 1. 几何与工况参数 (严格参考 test_flow_heat.py)
        # =================================================================
        self.L = 0.375  # [m]
        self.D_inner = 0.03  # [m]
        self.R_inner = self.D_inner / 2.0

        self.W_in = 0.33  # [kg/s] (入口流量边界)
        self.P_out = 116000.0  # [Pa]   (出口压力边界)
        self.T_init = 743.0  # [K]

        self.n_axial = n_axial

        # 辅助几何
        self.area_flow = np.pi * self.R_inner ** 2

        # 初始压力猜测 (假设只有微小压降)
        self.P_in_guess = self.P_out + 100.0

        # 材料
        self.mat_na = Sodium()

        self.manager = None
        self.y0 = None

    def build_system(self):
        print(f"\n🔧 Building Hydro System (Based on HeatFlowSystem Geometry)")
        print(f"   L={self.L}m, D={self.D_inner}m, Nodes={self.n_axial}")
        print(f"   BC: Win={self.W_in} kg/s, Pout={self.P_out} Pa")

        # =================================================================
        # 2. 流动建模 (Hydrodynamics) - 复刻 HeatFlowSystem
        # =================================================================
        # 2.1. 进出口边界
        # 入口边界: 这里的 P_in_guess 只是初始值，实际压力由 InletJunction 强行推入流量决定
        self.inlet = BoundaryVolume("inlet", self.mat_na, self.P_in_guess, self.T_init, self.area_flow, self.D_inner)
        self.outlet = BoundaryVolume("outlet", self.mat_na, self.P_out, self.T_init, self.area_flow, self.D_inner)

        # 2.2. 主流道
        self.channel = FluidChannel("channel", self.n_axial, self.L, self.area_flow, self.D_inner,
                                    initial_P=self.P_out, initial_T=self.T_init, material=self.mat_na)

        # 2.3. 连接 (Junctions)
        # [关键] 使用 InletJunction，这意味着流量 W 是强制边界条件，而不是由压差计算出来的
        self.junc_in = InletJunction("junc_in", from_vol=self.inlet, to_vol=self.channel.volumes[0],
                                     W_initial=self.W_in)
        self.junc_in.k_loss = 0.5

        self.junc_out = FlowJunction("junc_out", from_vol=self.channel.volumes[-1], to_vol=self.outlet,
                                     flow_area=self.area_flow, k_loss=0.5)

        # 2.4. 初始化状态
        self.inlet.update_properties(self.inlet.material)
        self.outlet.update_properties(self.outlet.material)
        self.channel.initialize_state()

        # 初始时刻，给内部所有节点赋予流量初值，加速收敛
        for junc in self.channel.internal_junctions:
            junc.W = self.W_in
            junc.update_velocity()
        self.junc_out.W = self.W_in

        # 2.5. 收集组件
        all_vols = [self.inlet] + self.channel.volumes + [self.outlet]
        all_juncs = [self.junc_in] + self.channel.internal_junctions + [self.junc_out]

        # =================================================================
        # 3. SystemManager 包装
        # =================================================================
        # 仅传入流体组件，couples 和 components 留空
        self.manager = SystemManager(
            volumes=all_vols,
            junctions=all_juncs,
            couples=[],
            components=[]
        )

        self.y0 = self.manager.get_initial_state()
        print(f"   System Dim: {self.manager.dim}")

    def run_simulation(self):
        # 1. 准备 Jacobian
        jac_sparsity = None
        if USE_JAC:
            jac_builder = FluidJacobianBlockLayout(self.manager.volumes, self.manager.junctions)
            jac_sparsity = jac_builder.get_sparsity_matrix()

        # 2. 设置求解器
        # BDF 或 Radau 均可。因为这里是强制流量入口，刚性可能不如压差驱动那么强，但仍建议用隐式
        solver = NuclearODESolver(method='BDF', rtol=1e-1, atol=1e-2)

        res = solver.solve(
            fun=self.manager.dydt_function_flow_init,
            t_span=(0.0, 2.0),
            y0=self.y0,
            jac_sparsity=jac_sparsity
        )

        print(f"\n🚀 Starting Hydro Simulation: 0.0s -> 2.0s")
        start_time = time.time()

        solver.rtol = 1e-3
        solver.atol = 1e-6

        # 3. 调用求解
        # 使用 flow_init 模式的导数函数 (只计算流体)
        res = solver.solve(
            fun=self.manager.dydt_function_flow_init,
            t_span=(0.0, 2.0),
            y0=self.y0,
            jac_sparsity=jac_sparsity
        )

        end_time = time.time()

        if res['success']:
            print(f"✅ Simulation Complete in {end_time - start_time:.3f}s")

            y_final = res['y'][:, -1]
            self.manager.update_system_state(y_final)

            # 结果分析
            # 对于 InletJunction，流量应该是恒定的 0.33
            # 我们主要关注压力分布
            P_inlet_node = self.channel.volumes[0].P
            P_outlet_node = self.channel.volumes[-1].P
            dP_total = P_inlet_node - self.P_out

            print(f"\n📊 Result Analysis (Steady State):")
            print(f"   Imposed Flow : {self.junc_in.W:.4f} kg/s")
            print(f"   Outlet Flow  : {self.junc_out.W:.4f} kg/s (Should match inlet)")
            print(f"   Inlet Node P : {P_inlet_node:.1f} Pa")
            print(f"   Outlet BC P  : {self.P_out:.1f} Pa")
            print(f"   Total dP     : {dP_total:.1f} Pa")

            self.plot_pressure_profile(res)
        else:
            print(f"❌ Simulation Failed: {res['message']}")

    def plot_pressure_profile(self, res):
        # 提取最终时刻的压力分布
        # SystemManager 布局: [P...P, h...h, W...W]
        # P 的索引是 0 到 n_vol-1

        # 注意: manager.volumes 包含了 boundary volumes (inlet, outlet)
        # 顺序: [InletVol, ChVol1, ..., ChVolN, OutletVol]

        y_final = res['y'][:, -1]
        n_vol = self.manager.n_vol
        P_all = y_final[0:n_vol]

        # 提取坐标 (仅针对 Channel 内部节点)
        # InletVol 和 OutletVol 是虚拟边界，没有明确的空间坐标用于绘图，或者设为 0 和 L

        z_coords = []
        P_channel = []

        # 遍历 manager.volumes，只取 FluidChannel 的部分
        for i, vol in enumerate(self.manager.volumes):
            if vol in self.channel.volumes:
                # 假设 Channel 已经计算了 z_coordinate (HeatFlowSystem 中未显式计算，这里手动算一下)
                # self.L = 0.375, n=20, dx = 0.375/20
                idx_in_channel = self.channel.volumes.index(vol)
                dx = self.channel.node_length
                z = (idx_in_channel + 0.5) * dx
                z_coords.append(z)
                P_channel.append(P_all[i])

        plt.figure(figsize=(8, 5))
        plt.plot(z_coords, P_channel, 'o-', label='Pressure Profile')

        # 标记出口边界
        plt.plot(self.L, self.P_out, 'rx', label='Outlet BC')

        plt.title(f'Pressure Distribution (Win={self.W_in}kg/s)')
        plt.xlabel('Axial Position (m)')
        plt.ylabel('Pressure (Pa)')
        plt.grid(True)
        plt.legend()
        plt.show()

# 压力突变计算
class PressureSurgeSimulation:
    def __init__(self):
        # ==========================================
        # 1. 几何与初始工况 (参考 HeatFlowSystem)
        # ==========================================
        self.L = 0.375  # [m]
        self.D_inner = 0.03  # [m]
        self.area_flow = np.pi * (self.D_inner / 2) ** 2

        self.n_nodes = 20

        # 初始边界条件
        self.W_target = 0.33  # [kg/s] 目标流量
        self.P_out_initial = 116000.0  # [Pa] 初始出口压力
        self.T_init = 743.0  # [K]

        # 瞬态扰动设定
        self.surge_time = 1.0  # [s] 触发时间
        self.P_out_surge = 160000.0  # [Pa] 激增后的出口压力 (+44kPa)

        self.mat_na = Sodium()
        self.manager = None
        self.y0 = None

    def build_system(self):
        print(f"\n🔧 Building System for Pressure Surge Test")

        # 1. 建立边界 (使用新版 BoundaryVolume)
        # 注意：这里我们传入初始 P，内部会自动设置 target_P = initial_P
        self.inlet_bc = BoundaryVolume("InletBC", self.mat_na,
                                       P=self.P_out_initial + 1000,  # 初始猜测值
                                       T=self.T_init,
                                       flow_area=1.0)  # 大面积以模拟集管

        self.outlet_bc = BoundaryVolume("OutletBC", self.mat_na,
                                        P=self.P_out_initial,
                                        T=self.T_init,
                                        flow_area=1.0)

        # 2. 建立通道
        self.channel = FluidChannel("Channel", self.n_nodes, self.L, self.area_flow, self.D_inner,
                                    initial_P=self.P_out_initial,
                                    initial_T=self.T_init,
                                    material=self.mat_na)

        # 3. 建立连接
        # [关键] 使用新版 InletJunction，设定 W_initial
        self.junc_in = InletJunction("InletJunc", self.inlet_bc, self.channel.volumes[0],
                                     W_initial=self.W_target)
        self.junc_in.k_loss = 0.5

        # 普通出口连接
        self.junc_out = FlowJunction("OutletJunc", self.channel.volumes[-1], self.outlet_bc,
                                     flow_area=self.area_flow, k_loss=0.5)

        # 4. 初始化状态
        self.inlet_bc.update_properties(self.mat_na)
        self.outlet_bc.update_properties(self.mat_na)
        self.channel.initialize_state()

        # 给内部节点赋初始流量，帮助收敛
        for junc in self.channel.internal_junctions:
            junc.W = self.W_target
            junc.update_velocity()
        self.junc_out.W = self.W_target

        # 5. 包装
        # 只需要流体部分
        all_vols = [self.inlet_bc] + self.channel.volumes + [self.outlet_bc]
        all_juncs = [self.junc_in] + self.channel.internal_junctions + [self.junc_out]

        self.manager = SystemManager(volumes=all_vols, junctions=all_juncs)
        self.y0 = self.manager.get_initial_state()

    def dydt_wrapper(self, t, y):
        """
        包装导数函数，用于在运行过程中动态修改边界条件
        """
        # [核心逻辑]：在 t > surge_time 时修改出口边界的目标压力
        # 注意：不要每一步都 set_state，只要在跨越时间点时设置一次即可
        # 但为了简单，这里每次调用 set_state 开销很小（只是赋值 target_P），是安全的

        if t >= self.surge_time:
            # 施加压力激增：设定新的目标压力
            self.outlet_bc.set_state(P=self.P_out_surge)
        else:
            # 保持初始压力
            self.outlet_bc.set_state(P=self.P_out_initial)

        # 调用原有的导数计算
        # 使用 flow_init 版本只计算流体，效率高
        # return self.manager.dydt_function_flow_init(t, y)
        return self.manager.dydt_function(t,y)
    def run(self):
        # 1. 雅可比
        jac_sparsity = None
        if USE_JAC:
            # 新版 BoundaryVolume 的导数只依赖自身 (dydt = K*(Target - P))
            # 现有的 FluidJacobianBlockLayout 包含对角项，所以完全兼容
            jac_builder = FluidJacobianBlockLayout(self.manager.volumes, self.manager.junctions)
            jac_sparsity = jac_builder.get_sparsity_matrix()

        # 2. 求解器
        # 刚性求解器 (Radau 或 BDF)
        solver = NuclearODESolver(method='Radau', rtol=1e-1, atol=1e-1)

        # 3. 求解
        # res = solver.solve(
        #     fun=self.dydt_wrapper,  # 使用包装后的函数
        #     t_span=(0.0, 4.0),  # 总时长
        #     y0=self.y0,
        #     jac_sparsity=jac_sparsity,
        #     max_step=0.05  # 限制最大步长以确保捕捉到扰动瞬间
        # )

        print(f"🚀 Starting Pressure Surge Simulation")
        print(f"   0.0s -> {self.surge_time}s : Steady State (P_out = {self.P_out_initial / 1000:.1f} kPa)")
        print(f"   {self.surge_time}s -> 3.0s : SURGE (Target P_out -> {self.P_out_surge / 1000:.1f} kPa)")

        start_time = time.time()

        solver.rtol = 1e-3
        solver.rtol = 1e-3

        # 3. 求解
        res = solver.solve(
            fun=self.dydt_wrapper,  # 使用包装后的函数
            t_span=(0.0, 3.0),  # 总时长
            y0=self.y0,
            jac_sparsity=jac_sparsity,
            max_step=0.05  # 限制最大步长以确保捕捉到扰动瞬间
        )

        print(f"✅ Done in {time.time() - start_time:.3f}s")

        if res['success']:
            self.plot_results(res)
        else:
            print(f"❌ Failed: {res['message']}")

    def plot_results(self, res):
        t = res['t']
        y = res['y']

        # 提取数据 (SystemManager 布局: P...P, h...h, W...W)
        n_vol = self.manager.n_vol

        # 1. 出口边界压力 (OutletBC, index = -1)
        # 注意: manager.volumes 最后一个是 OutletBC
        P_outlet_hist = y[n_vol - 1, :]

        # 2. 入口通道压力 (Channel Node 0, index = 1)
        # index 0 是 InletBC，index 1 是 Channel[0]
        P_ch_inlet_hist = y[1, :]

        # 3. 入口流量 (InletJunc, index = 2*n_vol)
        idx_W_start = 2 * n_vol
        W_inlet_hist = y[idx_W_start, :]

        # 绘图
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # 上图：压力响应
        ax1.plot(t, P_outlet_hist / 1000, 'r--', label='Outlet BC Pressure (State)')
        ax1.plot(t, P_ch_inlet_hist / 1000, 'b-', label='Channel Inlet Pressure')

        # 画出目标设定值的阶跃线
        t_set = [0, self.surge_time, self.surge_time, 3.0]
        p_set = [self.P_out_initial, self.P_out_initial, self.P_out_surge, self.P_out_surge]
        ax1.plot(t_set, np.array(p_set) / 1000, 'g:', alpha=0.5, label='Outlet Target Setting')

        ax1.set_ylabel('Pressure (kPa)')
        ax1.set_title('Pressure Surge Response')
        ax1.legend()
        ax1.grid(True)

        # 下图：流量响应
        ax2.plot(t, W_inlet_hist, 'k-', label='Inlet Mass Flow')
        ax2.set_ylabel('Mass Flow (kg/s)')
        ax2.set_xlabel('Time (s)')
        ax2.set_title('Inlet Flow Regulation (Target=0.33 kg/s)')
        ax2.legend()
        ax2.grid(True)

        plt.tight_layout()
        plt.show()

# 压力-流量突变
class NaturalReductionTest:
    def __init__(self):
        # 1. 几何参数 (参考 test_flow_heat)
        self.L = 0.375
        self.D_inner = 0.03
        self.area_flow = np.pi * (self.D_inner / 2) ** 2
        self.n_nodes = 20

        # 2. 初始工况
        self.W_target = 0.33  # [kg/s] 初始目标流量
        self.P_out_initial = 116000.0  # [Pa] 初始出口压力
        self.T_init = 743.0

        # 3. 瞬态参数
        self.surge_time = 1.0  # [s] 扰动开始时间
        # [注意]：由于该通道流阻很小(~200Pa)，如果出口压力增加太多，流量会瞬间倒流。
        # 这里演示一个温和的压力上升 (+500 Pa)，以便观察流量下降过程。
        self.surge_delta_P = 500.0

        self.mat_na = Sodium()
        self.manager = None
        self.y0 = None

    def build_system(self, inlet_pressure_guess):
        """
        构建系统，入口压力由参数指定
        """
        # 1. 入口压力边界 (Inlet Plenum)
        # 这是一个无穷大容器，压力恒定为 inlet_pressure_guess
        self.inlet_bc = BoundaryVolume("InletPlenum", self.mat_na,
                                       P=inlet_pressure_guess,
                                       T=self.T_init, flow_area=1.0)

        # 2. 出口压力边界 (Outlet Plenum)
        self.outlet_bc = BoundaryVolume("OutletPlenum", self.mat_na,
                                        P=self.P_out_initial,
                                        T=self.T_init, flow_area=1.0)

        # 3. 通道
        self.channel = FluidChannel("Channel", self.n_nodes, self.L, self.area_flow, self.D_inner,
                                    initial_P=self.P_out_initial,
                                    initial_T=self.T_init,
                                    material=self.mat_na)

        # 4. 连接 (全部使用普通 FlowJunction)
        # [核心修改] 这里不再使用 InletJunction，而是普通连接
        # 流量将完全由 (P_inlet_bc - P_channel_inlet) 决定
        self.junc_in = FlowJunction("JuncIn", self.inlet_bc, self.channel.volumes[0],
                                    flow_area=self.area_flow, k_loss=0.5)

        self.junc_out = FlowJunction("JuncOut", self.channel.volumes[-1], self.outlet_bc,
                                     flow_area=self.area_flow, k_loss=0.5)

        # 5. 初始化
        self.inlet_bc.update_properties(self.mat_na)
        self.outlet_bc.update_properties(self.mat_na)
        self.channel.initialize_state()

        # 赋予初始流量猜想值 (仅用于加速收敛，不强制)
        for junc in self.channel.internal_junctions:
            junc.W = self.W_target
        self.junc_in.W = self.W_target
        self.junc_out.W = self.W_target

        # 6. 包装
        all_vols = [self.inlet_bc] + self.channel.volumes + [self.outlet_bc]
        all_juncs = [self.junc_in] + self.channel.internal_junctions + [self.junc_out]

        self.manager = SystemManager(volumes=all_vols, junctions=all_juncs)
        self.y0 = self.manager.get_initial_state()

    def find_steady_state_pressure(self):
        """
        预计算：寻找维持 0.33 kg/s 所需的入口压力
        """
        print("🔍 Searching for required Inlet Pressure...")

        # 估算流阻: dP ~ K * 0.5 * rho * v^2
        # v = 0.55 m/s, rho=850, K_tot ~ 1.5 -> dP ~ 200 Pa
        estimated_dP = 200.0
        current_P_in = self.P_out_initial + estimated_dP

        # 简单迭代几次找到准确的 P_in
        for i in range(3):
            self.build_system(current_P_in)
            # 跑一小段稳态
            solver = NuclearODESolver(method='BDF', rtol=1e-3, atol=1e-3)
            res = solver.solve(self.manager.dydt_function_flow_init, (0, 1.0), self.y0)

            # 获取计算出的稳态流量
            y_final = res['y'][:, -1]
            self.manager.update_system_state(y_final)
            w_result = self.junc_in.W

            print(f"   Iter {i}: P_in={current_P_in:.1f} Pa -> Flow={w_result:.4f} kg/s")

            # 根据误差调整 P_in (简单的比例调整)
            # Flow ~ sqrt(dP) -> dP ~ Flow^2
            # New_dP = Old_dP * (Target_W / Current_W)^2
            current_dP = current_P_in - self.P_out_initial
            new_dP = current_dP * (self.W_target / w_result) ** 2
            current_P_in = self.P_out_initial + new_dP

        print(f"✅ Found Inlet Pressure: {current_P_in:.2f} Pa (dP = {new_dP:.2f} Pa)")
        return current_P_in

    def dydt_transient_wrapper(self, t, y):
        # 在 t > surge_time 时，提高出口边界的目标压力
        if t >= self.surge_time:
            # 模拟下游背压升高
            target_P = self.P_out_initial + self.surge_delta_P
            self.outlet_bc.set_state(P=target_P)
        else:
            self.outlet_bc.set_state(P=self.P_out_initial)

        # return self.manager.dydt_function_flow_init(t, y)
        return self.manager.dydt_function(t, y)

    def run(self):
        # 1. 确定入口压力
        required_P_in = self.find_steady_state_pressure()

        # 2. 重建最终系统
        self.build_system(required_P_in)

        # 3. 设置求解器
        jac_sparsity = None
        if USE_JAC:
            jac = FluidJacobianBlockLayout(self.manager.volumes, self.manager.junctions)
            jac_sparsity = jac.get_sparsity_matrix()

        solver = NuclearODESolver(method='Radau', rtol=1e-5, atol=1e-5)

        print(f"\n🚀 Starting Natural Reduction Test")
        print(f"   Time 0.0 -> {self.surge_time}s: Steady State")
        print(f"   Time {self.surge_time} -> 3.0s: Outlet Pressure +{self.surge_delta_P:.1f} Pa")
        print(f"   (Inlet Pressure held constant at {required_P_in:.1f} Pa)")

        # 4. 运行
        res = solver.solve(
            fun=self.dydt_transient_wrapper,
            t_span=(0.0, 3.0),
            y0=self.y0,
            jac_sparsity=jac_sparsity
        )

        if res['success']:
            self.plot_results(res)
        else:
            print("❌ Failed.")

    def plot_results(self, res):
        t = res['t']
        y = res['y']
        n_vol = self.manager.n_vol

        # 提取流量 (JuncIn)
        idx_W = 2 * n_vol  # 第一个 Junction 就是 JuncIn
        W_trace = y[idx_W, :]

        # 提取压力 (OutletBC)
        P_outlet_trace = y[n_vol - 1, :]

        fig, ax1 = plt.subplots(figsize=(10, 6))

        # 绘制压力
        color = 'tab:red'
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Outlet Pressure (Pa)', color=color)
        ax1.plot(t, P_outlet_trace, color=color, linestyle='--', label='Outlet P (Surge)')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True)

        # 绘制流量
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('Mass Flow (kg/s)', color=color)
        ax2.plot(t, W_trace, color=color, linewidth=2, label='Mass Flow Response')
        ax2.tick_params(axis='y', labelcolor=color)

        # 标注入门流量线
        ax2.axhline(self.W_target, color='green', linestyle=':', alpha=0.5, label='Initial Flow')

        plt.title('Natural Flow Reduction due to Back-Pressure Surge')
        fig.tight_layout()
        plt.show()

#

class HeatTransferTest:
    def __init__(self):
        # 1. 几何参数
        self.L = 0.375  # [m]
        self.D_inner = 0.03  # [m]
        self.area_flow = np.pi * (self.D_inner / 2) ** 2
        self.perimeter = np.pi * self.D_inner
        self.n_nodes = 20

        # 2. 工况参数
        self.W_in = 0.33  # [kg/s] 入口流量
        self.P_out = 116000.0  # [Pa] 出口压力
        self.T_inlet = 743.0  # [K] 入口流体温度

        # 3. 热边界参数 (恒温壁面)
        self.T_wall_const = 850.0  # [K]

        self.mat_na = Sodium()
        self.manager = None
        self.y0 = None

    def build_system(self):
        print(f"\n🔧 Building Heat Transfer System (Constant Wall T={self.T_wall_const}K)")

        # --- A. 流体系统建模 ---

        # 1. 边界容积
        # 入口：压力给一个猜测值，实际由流量决定；温度固定为 T_inlet
        self.inlet_bc = BoundaryVolume("InletPlenum", self.mat_na,
                                       P=self.P_out + 500, T=self.T_inlet, flow_area=1.0)

        self.outlet_bc = BoundaryVolume("OutletPlenum", self.mat_na,
                                        P=self.P_out, T=self.T_inlet, flow_area=1.0)

        # 2. 通道 (初始温度设为入口温度)
        self.channel = FluidChannel("Channel", self.n_nodes, self.L, self.area_flow, self.D_inner,
                                    initial_P=self.P_out,
                                    initial_T=self.T_inlet,  # 初始全场冷态
                                    material=self.mat_na)

        # 3. 连接 (使用 FlowJunction, 我们通过寻找压力来控制流量，或简化直接给初值)
        # 为了简化测试，这里我们假设已经通过之前的步骤找到了稳态流量对应的入口压力
        # 或者我们使用 InletJunction (如果BoundaryVolume支持) 来强迫流量
        # 这里为了物理自洽，我们使用普通连接，并在 dydt 中保持入口恒定压力

        # 预估压降 ~200Pa
        P_in_est = self.P_out + 200.0
        self.inlet_bc.P = P_in_est
        self.inlet_bc.target_P = P_in_est

        self.junc_in = FlowJunction("JuncIn", self.inlet_bc, self.channel.volumes[0],
                                    flow_area=self.area_flow, k_loss=0.5)
        self.junc_out = FlowJunction("JuncOut", self.channel.volumes[-1], self.outlet_bc,
                                     flow_area=self.area_flow, k_loss=0.5)

        # 初始化流体状态
        self.inlet_bc.update_properties(self.mat_na)
        self.outlet_bc.update_properties(self.mat_na)
        self.channel.initialize_state()

        # 赋流量初值
        for junc in self.channel.internal_junctions: junc.W = self.W_in
        self.junc_in.W = self.W_in
        self.junc_out.W = self.W_in

        # --- B. 热边界建模 (虚拟固体) ---

        # 1. 创建边界区域 (代表管壁内表面)
        # 计算每个节点的换热面积
        node_len = self.L / self.n_nodes
        area_per_node = self.perimeter * node_len
        area_array = np.full(self.n_nodes, area_per_node)

        self.wall_boundary = BoundaryRegion(shape=(self.n_nodes,), area_array=area_array)

        # [关键]：设定恒定壁温
        # 直接修改 T_surface。因为没有固体求解器去覆盖它，它将保持恒定。
        self.wall_boundary.T_surface[:] = self.T_wall_const

        # --- C. 耦合器建模 ---

        def correlation_robust(Re, Pr, ratio):
            Nu = nu_aoki(Re, Pr)
            # 限制 Nu 范围以保证数值稳定性 (特别是在流速极低或刚启动时)
            return np.clip(Nu, 4.36, 50.0)

        self.coupler = FluidSolidCouple(
            name="ChannelWallCouple",
            fluid=self.channel,
            solid_boundary_region=self.wall_boundary,
            heated_perimeter=self.perimeter,
            correlation_func=correlation_robust  # 注入您的关联式
        )

        # --- D. 系统包装 ---

        all_vols = [self.inlet_bc] + self.channel.volumes + [self.outlet_bc]
        all_juncs = [self.junc_in] + self.channel.internal_junctions + [self.junc_out]

        # 注意：couples 列表加入我们的耦合器
        self.manager = SystemManager(volumes=all_vols, junctions=all_juncs, couples=[self.coupler])

        self.y0 = self.manager.get_initial_state()

    # ==========================================================================
    # 3. 核心导数函数 (含能量源项处理)
    # ==========================================================================

    def dydt_coupled_flow(self, t, y):
        print(t)

        # 1. 将求解器向量 y 更新到对象状态 (P, h, W)
        self.manager.update_system_state(y)

        # 2. [关键] 清除上一时间步累积的源项
        # 必须调用，否则 Q_wall 会一直累加导致温度爆炸
        self.channel.clear_sources()

        # 3. [关键] 执行耦合计算
        # 这会自动：
        #   a. 计算 Re, Pr, Nu, h
        #   b. 读取 wall_boundary.T_surface (恒定850K)
        #   c. 计算 Q = hA(Tw - Tf)
        #   d. 调用 fluid.add_coupling_source 注入源项
        self.coupler.execute()

        # 4. 计算导数
        # 使用 standard dydt_function (它会调用 get_volume_derivatives)
        # get_volume_derivatives 内部使用了 energy_balance_rhs (包含 Q_wall)

        # 由于 SystemManager 默认 dydt 可能包含了一些不需要的固体逻辑(如果有组件的话)
        # 这里我们手动写一个纯流体+耦合的导数流程，或者复用 manager 的逻辑
        # 为了安全，这里手动构建，确保逻辑透明

        dydt = np.zeros(len(y))
        n_vol = self.manager.n_vol
        n_junc = self.manager.n_junc

        # A. Fluid Volumes (P, h)
        for i, vol in enumerate(self.manager.volumes):
            if isinstance(vol, BoundaryVolume):
                # 调用边界容积的导数 (支持松弛驱动或为0)
                dP, dh = vol.get_volume_derivatives(vol.material)
            else:
                # 内部容积 (包含耦合源项)
                dP, dh = vol.get_volume_derivatives(vol.material)

            dydt[i] = dP
            dydt[n_vol + i] = 0.0

        # B. Fluid Junctions (W)
        for i, junc in enumerate(self.manager.junctions):
            mat = junc.from_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[2 * n_vol + i] = dW

        return dydt

    def dydt_coupled(self, t, y):
        print(t)

        # 1. 将求解器向量 y 更新到对象状态 (P, h, W)
        self.manager.update_system_state(y)

        # 2. [关键] 清除上一时间步累积的源项
        # 必须调用，否则 Q_wall 会一直累加导致温度爆炸
        self.channel.clear_sources()

        # 3. [关键] 执行耦合计算
        # 这会自动：
        #   a. 计算 Re, Pr, Nu, h
        #   b. 读取 wall_boundary.T_surface (恒定850K)
        #   c. 计算 Q = hA(Tw - Tf)
        #   d. 调用 fluid.add_coupling_source 注入源项
        self.coupler.execute()

        # 4. 计算导数
        # 使用 standard dydt_function (它会调用 get_volume_derivatives)
        # get_volume_derivatives 内部使用了 energy_balance_rhs (包含 Q_wall)

        # 由于 SystemManager 默认 dydt 可能包含了一些不需要的固体逻辑(如果有组件的话)
        # 这里我们手动写一个纯流体+耦合的导数流程，或者复用 manager 的逻辑
        # 为了安全，这里手动构建，确保逻辑透明

        dydt = np.zeros(len(y))
        n_vol = self.manager.n_vol
        n_junc = self.manager.n_junc

        # A. Fluid Volumes (P, h)
        for i, vol in enumerate(self.manager.volumes):
            if isinstance(vol, BoundaryVolume):
                # 调用边界容积的导数 (支持松弛驱动或为0)
                dP, dh = vol.get_volume_derivatives(vol.material)
            else:
                # 内部容积 (包含耦合源项)
                dP, dh = vol.get_volume_derivatives(vol.material)

            dydt[i] = 0.0
            dydt[n_vol + i] = dh

        # B. Fluid Junctions (W)
        for i, junc in enumerate(self.manager.junctions):
            mat = junc.from_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[2 * n_vol + i] = dW

        return dydt

    def run(self):
        # 雅可比矩阵 (流体部分)
        jac_sparsity = None
        if USE_JAC:
            jac = FluidJacobianBlockLayout(self.manager.volumes, self.manager.junctions)
            jac_sparsity = jac.get_sparsity_matrix()

        # 求解器设置
        # 涉及能量方程，温度变化会导致密度变化，进而影响压力，建议使用隐式
        # 建立初始流动
        solver = NuclearODESolver(method='Radau', rtol=1e-3, atol=1e-3)

        res_init = solver.solve(
            fun=self.dydt_coupled_flow,
            t_span=(0.0, 5.0),
            y0=self.y0,
            jac_sparsity=jac_sparsity
        )

        y_final = res_init['y'][:, -1]

        self.y0 = y_final

        solver.rtol = 1e-5
        solver.atol = 1e-6

        print(f"🚀 Starting Coupled Heat Transfer Simulation")
        print(f"   Inlet T: {self.T_inlet} K")
        print(f"   Wall  T: {self.T_wall_const} K (Constant)")

        start_time = time.time()

        # 运行 2.0 秒，足够建立热平衡
        # 初始阶段液体是 743K，壁面是 850K，会有剧烈换热
        res = solver.solve(
            fun=self.dydt_coupled,
            t_span=(0.0, 50.0),
            y0=self.y0,
            jac_sparsity=jac_sparsity
        )

        print(f"✅ Done in {time.time() - start_time:.3f}s")

        if res['success']:
            self.analyze_and_plot(res)
        else:
            print(f"❌ Failed: {res['message']}")

    def analyze_and_plot(self, res):
        # 提取最终状态
        y_final = res['y'][:, -1]
        self.manager.update_system_state(y_final)

        # 提取沿程数据
        z_coords = []
        T_fluid = []
        T_wall = []
        h_coeffs = []  # 我们想看看换热系数分布

        # 重新执行一次 couple 以获取当前的 h (因为 h 不在状态向量里)
        self.channel.clear_sources()
        self.coupler.execute()

        # 获取 Nu/h 数据 (需要一点黑客手段，或者我们在 couple 里 print)
        # 这里我们手动算一下用于绘图
        # execute() 已经把 R_ext 写进了 solid_bc
        # R_conv = 1 / (hA) -> h = 1 / (R_conv * A)

        R_ext_vals = self.coupler.solid_bc.R_ext
        areas = self.coupler.node_areas
        h_calc = 1.0 / (R_ext_vals * areas)

        for i, vol in enumerate(self.channel.volumes):
            # 几何中心坐标
            z = (i + 0.5) * self.channel.node_length
            z_coords.append(z)
            T_fluid.append(vol.T)
            T_wall.append(self.T_wall_const)  # 恒定
            h_coeffs.append(h_calc[i])

        # 绘图
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # 上图：温度分布
        ax1.plot(z_coords, T_wall, 'r-', linewidth=2, label='Wall Temperature (850K)')
        ax1.plot(z_coords, T_fluid, 'b-o', label='Fluid Temperature')
        ax1.set_ylabel('Temperature (K)')
        ax1.set_title('Axial Temperature Profile (Steady State)')
        ax1.legend()
        ax1.grid(True)

        # 下图：换热系数
        ax2.plot(z_coords, h_coeffs, 'g-s', label='Heat Transfer Coeff (h)')
        ax2.set_ylabel('h (W/m^2-K)')
        ax2.set_xlabel('Axial Position (m)')
        ax2.set_title('Heat Transfer Coefficient Distribution')
        ax2.legend()
        ax2.grid(True)

        plt.tight_layout()
        plt.show()

        # 能量守恒检查
        # Q_total_in = m * (h_out - h_in)
        # Q_wall_total = sum(Q_wall)

        h_in = self.inlet_bc.h
        h_out = self.outlet_bc.h  # 或者是 channel[-1].h
        # 注意: outlet_bc 的 h 是通过 update_properties 算的，取决于 T_outlet
        # 实际上我们要看流出流道的焓
        h_exit = self.channel.volumes[-1].h
        m_flow = self.junc_in.W  # 稳态下入口流量

        Q_fluid_gain = m_flow * (h_exit - h_in)

        print(f"\n📊 Energy Balance Check:")
        print(f"   Fluid Inlet T : {self.channel.volumes[0].T:.2f} K")
        print(f"   Fluid Exit  T : {self.channel.volumes[-1].T:.2f} K")
        print(f"   Fluid Heat Gain: {Q_fluid_gain / 1000:.2f} kW")
        # 由于我们没有直接累加 Q_wall，这里略过对比，但直观上 T 升高说明加热成功。


# ==============================================================================
# 主执行流程
# ==============================================================================
if __name__ == "__main__":
    test = HeatTransferTest()
    test.build_system()
    test.run()

    # test = NaturalReductionTest()
    # test.run()

    # sim = PressureSurgeSimulation()
    # sim.build_system()
    # sim.run()

    # test = HydroOnlyHeatCase(n_axial=20)
    # test.build_system()
    # test.run_simulation()

    # # 1. 正常实例化 (还是原来的配方)
    # test1 = HeatFlowSystem(n_axial=20, n_radial=3)
    #
    # # 【关键一步】清空组件列表！
    # # 这告诉求解器："不要计算固体的温度变化，把它当成死的背景"
    # test1.components = []
    #
    # # 2. 初始化 (此时状态向量 y 只包含流体，不含固体)
    # test1.set_state()
    # test1.flow_init()
    #
    # # 3. 施加阶跃突增
    # print("\n⚡ 触发温度突增：内壁面直接变为 820K")
    # # 直接修改固体组件内部的温度数组
    # test1.wall_conduction.T[:] = 820.0
    #
    # # ⚠️ 重要：手动刷新一下边界状态，确保流固耦合器能读到这个新温度
    # test1.wall_conduction._update_boundaries_state()
    #
    # test1.flow_init()
    #
    # # 4. 运行瞬态计算
    # # 此时固体温度会死死固定在 820K，只测试流体压力是否爆炸
    # # 使用严格的 rtol 来捕捉刚性响应
    # test1.calculate(step=10, starttime=0.0, endtime=0.1, rtol=1e-6, atol=1e-1)
    #
    # # 5. 打印结果验证
    # print("\n📊 最终压力分布 (Pa):")
    # # 打印前5个节点的压力
    # print(test1.y0[0:5])


    # 1. 初始化系统
    # test1 = HeatFlowSystem(n_axial=40, n_radial=3)
    # test1.set_state()
    #
    # # 2. 流动初始化
    # test1.flow_init()
    # # 手动赋值结果
    # if test1.res_flow_init['success']:
    #     test1.y0 = test1.res_flow_init['y'][:, -1]
    #     test1.manager.update_system_state(test1.y0)
    #
    # # 3. 耦合瞬态（0~100s，初步计算）并保存状态
    # # test1.calculate(step=20, starttime=0.0, endtime=100.0, rtol=0.1, atol=0.1)
    # # test1.save_restart_file("restart_t100s")
    #
    # # 4. 继续计算（10~20s，精细计算）
    # # 考虑分段的 atol
    #
    # # 4.1. 定义各物理量的数量级估算
    # N_vol = test1.manager.n_vol  # 流体节点数
    # N_junc = test1.manager.n_junc  # 连接数
    # N_solid = test1.manager.dim_solid  # 固体节点数
    #
    # # 4.2. 设定针对性的绝对容限 (atol)
    # # 压力 P: 量级 ~1e5 Pa -> 允许误差 ~1 Pa
    # atol_P = [1.0] * N_vol
    #
    # # 焓 h: 量级 ~1e6 J/kg -> 允许误差 ~10 J/kg
    # atol_h = [10.0] * N_vol
    #
    # # 流量 W: 量级 ~0.3 kg/s -> 允许误差 ~1e-4 kg/s
    # atol_W = [1e-4] * N_junc
    #
    # # 温度 T (固体): 量级 ~800 K -> 允许误差 1e-2 K
    # atol_T = [1e-2] * N_solid
    #
    # # 拼接成完整的 atol 向量
    # atol_custom = atol_P + atol_h + atol_W + atol_T
    #
    # test1.calculate(step=20, starttime=0.0, endtime=1.0, rtol=1e-1, atol=atol_custom, max_step=1)
    # test1.calculate(step=20, starttime=1.0, endtime=10.0, rtol=1e-1, atol=atol_custom, max_step=1)
    # test1.calculate(step=20, starttime=10.0, endtime=100.0, rtol=1e-1, atol=atol_custom, max_step=1)
    # test1.calculate(step=20, starttime=1000.0, endtime=2500.0, rtol=1e-1, atol=atol_custom, max_step=1)
    #

    print("end")
