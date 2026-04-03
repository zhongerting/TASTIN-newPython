import numpy as np
import matplotlib.pyplot as plt
from typing import List
import logging
import sys

# 物性组件
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel

# 流体组件
from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, FluidVolume
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume, InletJunction

# 二维导热组件
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D

# 求解器组件
from MathSolvers.solver_module import NuclearODESolver

# 耦合组件
from Solvers.Couplers import SolidSolidCouple2D
from Solvers.Couplers import FluidSolidCouple
from Correlations.Correlations import nu_aoki

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("TestRun")


class CoupledSystem2D:
    """
    流体-2D固体耦合系统容器
    负责物理对象的实例化和几何对齐
    """

    def __init__(self):
        self.total_dim = None
        self.sl_s_T = None
        self.idx_bv_h = None
        self.idx_bv_P = None
        self.idx_j_out = None
        self.idx_j_in = None
        self.sl_f_h = None
        self.sl_j_int = None
        self.n_solid_nodes = None
        self.sl_f_P = None
        self.n_junc_int = None
        self.n_fluid_vol = None
        self.coupler = None
        self.bc_outer = None
        logger.info("Initializing CoupledSystem2D...")

        # [新增] 热力求解开关
        # True: 正常计算温度变化; False: 强制温度不变(用于纯水力初始化)
        self.enable_thermal_solve = True

        # ==========================================
        # 1. 全局参数定义
        # ==========================================
        # 几何参数
        self.L_channel = 0.6  # 轴向长度 [m]
        self.D_inner = 0.03  # 流体直径 / 固体内径 [m]
        self.R_inner = self.D_inner / 2.0
        self.wall_thickness = 0.005  # 壁厚 [m]

        # 网格划分
        self.N_axial = 20  # 轴向控制体数量 (流体节点数 = 固体轴向层数)
        self.N_radial = 3  # 固体径向控制体层数

        # 初始条件
        self.T_init = 743.0  # 初始温度 [K]
        self.P_init = 116000.0  # 初始压力 [Pa] (116 kPa)
        self.W_init = 0.33  # 初始流量 [kg/s]

        # ==========================================
        # 2. 实例化材料
        # ==========================================
        self.mat_na = Sodium()
        self.mat_ss = AusteniticStainlessSteel()
        logger.info(f"Materials loaded: {self.mat_na.name}, {self.mat_ss.name}")

        # ==========================================
        # 3. 构建流体系统 (1D)
        # ==========================================
        # 计算流通面积
        area_flow = np.pi * (self.R_inner ** 2)

        # 3.1 边界控制体 (Boundary Volumes)
        # 进口: 用于设定 T_in 和 P_in (通过松弛控制)
        self.bv_in = BoundaryVolume(
            name="Inlet_Plenum",
            material=self.mat_na,
            P=self.P_init + 70.0,
            T=self.T_init
        )

        # 出口: 用于设定 P_out
        self.bv_out = BoundaryVolume(
            name="Outlet_Plenum",
            material=self.mat_na,
            P=self.P_init,  # 初始设为与进口一致，防止激波
            T=self.T_init
        )

        # 3.2 流体通道 (Fluid Channel)
        self.channel = FluidChannel(
            name="Test_Section_Na",
            n_nodes=self.N_axial,
            total_length=self.L_channel,
            flow_area=area_flow,
            hydraulic_diam=self.D_inner,
            initial_P=self.P_init,
            initial_T=self.T_init,
            material=self.mat_na
        )

        # 初始化内部物性 (h, rho, mu)
        self.channel.initialize_state()

        # 3.3 连接组件 (Junctions)
        # 进口连接: 强制流量边界 (InletJunction)
        self.junc_in = InletJunction(
            name="Junc_Inlet",
            from_vol=self.bv_in,
            to_vol=self.channel.volumes[0],
            W_initial=self.W_init
        )
        # 将连接注册到通道入口槽 (Hook)
        self.channel.inlet_junction = self.junc_in

        # 出口连接: 普通连接，由压差驱动流出
        self.junc_out = FlowJunction(
            name="Junc_Outlet",
            from_vol=self.channel.volumes[-1],
            to_vol=self.bv_out
        )
        # 初始流量猜测
        self.junc_out.W = self.W_init
        # 将连接注册到通道出口槽
        self.channel.outlet_junction = self.junc_out

        logger.info(f"Fluid system constructed. Channel nodes: {self.channel.n_nodes}")

        # ==========================================
        # 4. 构建固体系统 (2D Cylindrical)
        # ==========================================
        # 4.1 生成 2D 网格
        # X 维度 -> 径向 (Radial), Y 维度 -> 轴向 (Axial)
        # x_dim: 壁厚, y_dim: 长度
        self.mesh_solid = Mesh2D(
            x_dim=self.wall_thickness,
            n_x=self.N_radial,
            y_dim=self.L_channel,
            n_y=self.N_axial,
            geometry_type='cylindrical',
            inner_radius=self.R_inner
        )

        # 4.2 实例化导热求解器
        self.solid = HeatConduction2D(
            mesh=self.mesh_solid,
            material=self.mat_ss,
            initial_temp=self.T_init
        )

        logger.info(f"Solid system constructed. Mesh: 2D Cylindrical")
        logger.info(f"  - Radial (X): {self.mesh_solid.n_x} layers, Thickness={self.wall_thickness}m")
        logger.info(f"  - Axial  (Y): {self.mesh_solid.n_y} layers, Length={self.L_channel}m")
        logger.info(f"  - Total Solid Nodes: {self.solid.N}")

        # 验证拓扑对齐 (Sanity Check)
        if self.channel.n_nodes != self.mesh_solid.n_y:
            raise ValueError(
                f"Topology Mismatch: Fluid Nodes ({self.channel.n_nodes}) != Solid Axial Nodes ({self.mesh_solid.n_y})")

        logger.info("System Object Construction Completed.")

    def setup_boundaries_and_couplers(self):
        """
        [Step 2] 设置边界条件与耦合关系
        """
        logger.info("Setting up boundaries and couplers...")

        # ==========================================
        # 1. 固体外壁面边界 (Outer Wall)
        # ==========================================
        # 对应 Mesh2D 的 'right' 边界 (X/R 方向的终点)
        # 初始状态：绝热 (q_flux = 0.0)，稍后在计算脚本中开启加热
        self.bc_outer = self.solid.boundaries['right'].add_flux_condition(q_flux=0.0)
        logger.info("Solid Outer Boundary (Right): FluxBC initialized (q=0).")

        # ==========================================
        # 2. 流固耦合器 (Inner Wall <-> Fluid)
        # ==========================================
        # 定义换热关联式包装器
        # 签名需匹配: func(Re, Pr, P_D_ratio) -> Nu
        def correlation_wrapper(Re, Pr, ratio):
            # 液态金属 Aoki 关联式，不依赖 P/D (设为None或默认)
            return nu_aoki(Re, Pr)

        # 实例化耦合器
        # 注意: solid_boundary_region 使用 'left' (内表面)
        self.coupler = FluidSolidCouple(
            name="Coup_Na_SS_2D",
            fluid=self.channel,
            solid_boundary_region=self.solid.boundaries['left'],
            heated_perimeter=np.pi * self.D_inner,
            correlation_func=correlation_wrapper
        )
        logger.info(f"Coupler '{self.coupler.name}' initialized.")

        # ==========================================
        # 3. 构建状态向量索引映射 (Indexing)
        # ==========================================
        # 我们需要将所有微分变量打包成一个扁平向量 y

        # --- 计数 ---
        self.n_fluid_vol = self.channel.n_nodes
        self.n_junc_int = len(self.channel.internal_junctions)
        self.n_solid_nodes = self.solid.N  # Flattened 2D nodes (3 * 20 = 60)

        # --- 切片生成 (Slices) ---
        idx = 0

        # 1. 流体控制体 (P, h)
        self.sl_f_P = slice(idx, idx + self.n_fluid_vol)
        idx += self.n_fluid_vol
        self.sl_f_h = slice(idx, idx + self.n_fluid_vol)
        idx += self.n_fluid_vol

        # 2. 连接流量 (W)
        # 内部连接
        self.sl_j_int = slice(idx, idx + self.n_junc_int)
        idx += self.n_junc_int
        # 进口连接 (InletJunction)
        self.idx_j_in = idx
        idx += 1
        # 出口连接 (FlowJunction)
        self.idx_j_out = idx
        idx += 1

        # 3. 边界容积状态 (P, h) - 用于松弛边界
        # 顺序: [BV_In_P, BV_Out_P, BV_In_h, BV_Out_h]
        self.idx_bv_P = slice(idx, idx + 2)
        idx += 2
        self.idx_bv_h = slice(idx, idx + 2)
        idx += 2

        # 4. 固体温度 (T) - 2D Flattened
        self.sl_s_T = slice(idx, idx + self.n_solid_nodes)
        idx += self.n_solid_nodes

        self.total_dim = idx
        logger.info(f"State Vector Built. Total Dimensions: {self.total_dim}")
        logger.info(f"  - Fluid P/h: {self.n_fluid_vol}*2")
        logger.info(f"  - Junctions: {self.n_junc_int + 2}")
        logger.info(f"  - Solid T:   {self.n_solid_nodes} (3x20)")

    def pack_state(self) -> np.ndarray:
        """[Helper] 将对象状态打包为向量 y"""
        y = np.zeros(self.total_dim)
        # Fluid Volumes
        y[self.sl_f_P] = self.channel.pressure_vector
        y[self.sl_f_h] = [v.h for v in self.channel.volumes]

        # Junctions
        y[self.sl_j_int] = [j.W for j in self.channel.internal_junctions]
        y[self.idx_j_in] = self.junc_in.W
        y[self.idx_j_out] = self.junc_out.W

        # Boundary Volumes
        y[self.idx_bv_P] = [self.bv_in.P, self.bv_out.P]
        y[self.idx_bv_h] = [self.bv_in.h, self.bv_out.h]

        # Solid T (直接拷贝扁平数组)
        y[self.sl_s_T] = self.solid.T

        return y

    def unpack_state(self, y: np.ndarray):
        """[Helper] 将向量 y 分发回对象"""
        # 1. Fluid Volumes
        P_vec = y[self.sl_f_P]
        h_vec = y[self.sl_f_h]
        for i, v in enumerate(self.channel.volumes):
            v.P = P_vec[i]
            v.h = h_vec[i]
            # [关键] 必须立即更新导出物性 (T, rho, mu)，因为后续计算依赖它们
            v.update_properties(self.mat_na)

        # 2. Junctions
        W_int = y[self.sl_j_int]
        for i, j in enumerate(self.channel.internal_junctions):
            j.W = W_int[i]
            j.update_velocity()  # 更新流速

        self.junc_in.W = y[self.idx_j_in]
        self.junc_in.update_velocity()

        self.junc_out.W = y[self.idx_j_out]
        self.junc_out.update_velocity()

        # 3. Boundary Volumes
        bps = y[self.idx_bv_P]
        bhs = y[self.idx_bv_h]

        self.bv_in.P = bps[0]
        self.bv_in.h = bhs[0]
        self.bv_in.update_properties(self.mat_na)

        self.bv_out.P = bps[1]
        self.bv_out.h = bhs[1]
        self.bv_out.update_properties(self.mat_na)

        # 4. Solid T
        self.solid.T[:] = y[self.sl_s_T]

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        [Core] ODE 右端项函数 dy/dt = f(t, y)
        """
        # --- A. 状态同步 ---
        self.unpack_state(y)

        # --- B. 物理耦合计算 ---
        self.channel.clear_sources()

        # --- C. 计算导数 ---
        dydt = np.zeros(self.total_dim)

        # 1. 总是计算：流体动量 (dW/dt) 和 质量 (dP/dt)
        #    因为 P 和 W 是水力初始化的核心

        # dP/dt (Mass Balance)
        dP_list, dh_list = [], []
        for v in self.channel.volumes:
            dp, dh = v.get_volume_derivatives(self.mat_na)
            dP_list.append(dp)
            dh_list.append(dh)
        dydt[self.sl_f_P] = dP_list

        # dW/dt (Momentum)
        dW_list = []
        for j in self.channel.internal_junctions:
            dW_list.append(j.get_momentum_derivative(self.mat_na))
        dydt[self.sl_j_int] = dW_list
        dydt[self.idx_j_in] = self.junc_in.get_momentum_derivative(self.mat_na)
        dydt[self.idx_j_out] = self.junc_out.get_momentum_derivative(self.mat_na)

        # 边界容积 P
        dp_in, dh_in = self.bv_in.get_volume_derivatives(self.mat_na)
        dp_out, dh_out = self.bv_out.get_volume_derivatives(self.mat_na)
        dydt[self.idx_bv_P] = [dp_in, dp_out]

        # 2. 条件计算：热力方程 (dh/dt, dT/dt)
        if self.enable_thermal_solve:
            # --- 正常求解模式 ---
            dydt[self.sl_f_h] = dh_list
            dydt[self.idx_bv_h] = [dh_in, dh_out]

            # 固体导热
            dT_solid = self.solid.get_derivatives(t, self.solid.T)
            dydt[self.sl_s_T] = dT_solid
        else:
            # --- 纯水力初始化模式 (锁定温度) ---
            # 强制能量导数为 0
            dydt[self.sl_f_h] = 0.0
            dydt[self.idx_bv_h] = 0.0
            dydt[self.sl_s_T] = 0.0

            # 注意：get_volume_derivatives 内部已经计算了 dP，
            # 即使 dh/dt 被强行置零，dP/dt 依然有效（因为矩阵解耦或弱耦合）。

        return dydt


# ==============================================================================
# 主执行流程
# ==============================================================================
if __name__ == "__main__":
    # --- 1. 系统实例化 ---
    logger.info(">>> [Init] Instantiating System...")
    sys_coupled = CoupledSystem2D()
    sys_coupled.setup_boundaries_and_couplers()

    # 实例化求解器
    solver = NuclearODESolver(method='BDF', rtol=1e-3, atol=1e-6)

    # ==========================================================================
    # --- 2. Step 1: 纯水力初始化 (Hydraulic Init) ---
    # ==========================================================================
    logger.info(">>> [Step 1] Initializing Hydraulics (Thermal Locked)...")

    # [关键] 锁定热力计算，只算 P 和 W
    sys_coupled.enable_thermal_solve = False

    # 设定边界：零热流，但有流量
    sys_coupled.junc_in.set_flow_rate(0.33)
    sys_coupled.bc_outer.update_params(q_flux=0.0)

    y0 = sys_coupled.pack_state()

    # 计算 5 秒，足以让流量从 0.33 (初始猜测) 稳定下来，或消除压力的非物理震荡
    # 由于不涉及温度计算，这步通常非常快且稳定
    solver.rtol = 1e-3
    res1 = solver.solve(sys_coupled.rhs, (0.0, 5.0), y0, t_eval=[5.0])

    if not res1['success']:
        logger.error("Hydraulic init failed")
        sys.exit(1)

    y_hydro_ready = res1['y'][:, -1]
    logger.info("   -> Hydraulics established.")

    # ==========================================================================
    # --- 3. Step 2: 耦合加热计算 (Coupled Heating) ---
    # ==========================================================================
    logger.info(">>> [Step 2] Starting Coupled Heating (5s - 60s)...")

    # [关键] 解锁热力计算，开启加热
    sys_coupled.enable_thermal_solve = True
    sys_coupled.bc_outer.update_params(q_flux=10000.0)

    # 恢复高精度容差
    solver.rtol = 1e-6
    solver.atol = 1e-8

    # 从第 5 秒开始接续计算
    t_span = (5.0, 60.0)
    t_eval = np.linspace(5.0, 60.0, 56)

    res_heat = solver.solve(sys_coupled.rhs, t_span, y_hydro_ready, t_eval=t_eval)

    if not res_heat['success']:
        logger.error("Coupled heating calculation failed.")
        sys.exit(1)

    logger.info("   -> Coupled simulation completed.")

    # ==========================================================================
    # --- 4. Step 3: 数据提取与守恒验证 (Verification) ---
    # ==========================================================================
    logger.info(">>> [Step 3] Verification & Analysis")

    # 提取最终时刻状态
    y_final = res_heat['y'][:, -1]
    sys_coupled.unpack_state(y_final)

    # --- A. 质量守恒 ---
    W_in = sys_coupled.junc_in.W
    W_out = sys_coupled.junc_out.W
    mass_err = abs(W_in - W_out)

    # --- B. 能量守恒 ---
    # 1. 输入热量 (Wall Input)
    # 获取外壁面面积总和
    # Solid 'right' boundary 对应外表面
    area_outer = np.sum(sys_coupled.solid.boundaries['right'].area)
    Q_input = 10000.0 * area_outer  # [W]

    # 2. 输出热量 (Fluid Removal)
    # Q_fluid = W * (h_out - h_in)
    # 取最后一节流体控制体的焓 vs 进口边界焓
    h_in = sys_coupled.bv_in.h
    h_out = sys_coupled.channel.volumes[-1].h
    Q_fluid = W_out * (h_out - h_in)

    # 3. 误差
    energy_err_abs = Q_input - Q_fluid
    energy_err_rel = abs(energy_err_abs) / Q_input * 100.0

    print("-" * 60)
    print(f"VERIFICATION REPORT (t=60s)")
    print("-" * 60)
    print(f"Mass Balance:")
    print(f"  Win  : {W_in:.6f} kg/s")
    print(f"  Wout : {W_out:.6f} kg/s")
    print(f"  Error: {mass_err:.2e} kg/s")
    print("-" * 60)
    print(f"Energy Balance:")
    print(f"  Heat Input (Solid Wall) : {Q_input:.2f} W")
    print(f"  Heat Removal (Fluid)    : {Q_fluid:.2f} W")
    print(f"  Abs Error               : {energy_err_abs:.2f} W")
    print(f"  Rel Error               : {energy_err_rel:.2f} %")
    print("-" * 60)

    if energy_err_rel < 1.0:
        logger.info("VERIFICATION PASSED: Energy balance < 1%")
    else:
        logger.warning(
            f"VERIFICATION WARNING: Energy balance error is {energy_err_rel:.2f}% (System might not be fully steady yet)")

    # ==========================================================================
    # --- 5. 后处理绘图 (Plotting) ---
    # ==========================================================================

    # 提取数据
    z_coords = [v.z_coordinate for v in sys_coupled.channel.volumes]
    T_fluid = sys_coupled.channel.temperature_vector
    P_fluid = sys_coupled.channel.pressure_vector / 1000.0  # kPa

    # 提取 2D 固体温度
    # 扁平化逻辑: k = i * ny + j (i=Radial Layer, j=Axial Node)
    # Ny = 20
    Ny = sys_coupled.N_axial

    # 内壁面 (Layer 0, i=0) -> Indices 0 to 19
    T_solid_inner = sys_coupled.solid.T[0:Ny]

    # 外壁面 (Layer N_rad-1, i=2) -> Indices 40 to 59
    # 更通用的写法:
    start_idx = (sys_coupled.N_radial - 1) * Ny
    end_idx = start_idx + Ny
    T_solid_outer = sys_coupled.solid.T[start_idx:end_idx]

    # 绘图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Subplot 1: Temperatures
    ax1.plot(z_coords, T_solid_outer, 'r-s', label='Solid Outer Wall (Heated)', markersize=5)
    ax1.plot(z_coords, T_solid_inner, 'm-^', label='Solid Inner Wall (Coupled)', markersize=5)
    ax1.plot(z_coords, T_fluid, 'b-o', label='Sodium Fluid', linewidth=2)

    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Axial Temperature Distribution (t=60s)')
    ax1.legend()
    ax1.grid(True)

    # Subplot 2: Pressure & Radial Delta T
    ax2_r = ax2.twinx()
    ax2.plot(z_coords, P_fluid, 'g--', label='Fluid Pressure')
    ax2.set_ylabel('Pressure (kPa)', color='g')
    ax2.tick_params(axis='y', labelcolor='g')

    # 计算径向温差
    delta_T_wall = T_solid_outer - T_solid_inner
    ax2_r.bar(z_coords, delta_T_wall, width=0.01, alpha=0.3, color='orange', label='Wall Radial Delta T')
    ax2_r.set_ylabel('Radial Delta T (K)', color='orange')

    ax2.set_xlabel('Axial Position (m)')
    ax2.set_title('Pressure and Wall Radial Temperature Gradient')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('result_2d_coupled.png')
    logger.info("Plot saved to 'result_2d_coupled.png'")

    plt.show()
