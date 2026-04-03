import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

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

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("TestRun")

class CoupledSystem:
    def __init__(self, n_axial=20, n_radial=3):
        logger.info(f"Initializing System: Axial Nodes={n_axial}, Radial Layers={n_radial}")

        # 1. 几何与工况参数
        self.L = 0.6  # 长度 [m]
        self.D_inner = 0.03  # 内径 [m]
        self.R_inner = self.D_inner / 2.0
        self.thk = 0.005  # 壁厚 [m]

        self.W_in = 0.33  # 进口流量 [kg/s]
        self.P_out = 116000.0  # 出口压力 [Pa]
        self.T_init = 743.0  # 初始温度 [K]
        self.q_flux = 100.0  # 外壁热流 [W/m2]

        self.P_in_guess = self.P_out + 70.0

        self.n_axial = n_axial
        self.n_radial = n_radial

        # 2. 材料
        self.mat_na = Sodium()
        self.mat_ss = AusteniticStainlessSteel()

        # 3. 流体侧构建
        area_flow = np.pi * self.R_inner ** 2

        # 边界容积
        self.bv_in = BoundaryVolume("Inlet", self.mat_na, P=self.P_in_guess, T=self.T_init)  # P稍高以提供初始压差
        self.bv_out = BoundaryVolume("Outlet", self.mat_na, P=self.P_out, T=self.T_init)

        # 通道
        self.channel = FluidChannel(
            name="Pipe_Na", n_nodes=n_axial, total_length=self.L,
            flow_area=area_flow, hydraulic_diam=self.D_inner,
            initial_P=self.P_out, initial_T=self.T_init, material=self.mat_na
        )
        self.channel.initialize_state()

        p_initial = np.linspace(self.bv_in.P, self.bv_out.P, self.channel.n_nodes)

        for i, vol in enumerate(self.channel.volumes):
            vol.P = p_initial[i]

        for junc in self.channel.internal_junctions:
            junc.W = self.W_in

        # 连接
        self.junc_in = InletJunction("J_In", self.bv_in, self.channel.volumes[0], W_initial=self.W_in)
        self.channel.inlet_junction = self.junc_in

        self.junc_out = FlowJunction("J_Out", self.channel.volumes[-1], self.bv_out)
        self.junc_out.W = self.W_in
        self.channel.outlet_junction = self.junc_out

        # 4. 固体侧构建 (2D 圆柱坐标)
        # Mesh2D: x=Radial, y=Axial
        self.mesh = Mesh2D(
            x_dim=self.thk, n_x=n_radial,
            y_dim=self.L, n_y=n_axial,
            geometry_type='cylindrical', inner_radius=self.R_inner
        )

        self.solid = HeatConduction2D(self.mesh, self.mat_ss, initial_temp=self.T_init)

        # 5. 边界与耦合
        self.setup_boundaries_and_couplers()

        # 6. 状态向量索引映射
        self.build_indexing()

    def setup_boundaries_and_couplers(self):
        # A. 固体外壁面 (Right Boundary) -> 固定热流
        # 在 Mesh2D Cylindrical 中，'right' 对应外半径方向
        self.bc_outer = self.solid.boundaries['right'].add_flux_condition(q_flux=self.q_flux)

        # B. 流固耦合 (Inner Wall <-> Fluid)
        # 关联式包装器：增加限幅防止刚性震荡
        def correlation_robust(Re, Pr, ratio):
            Nu = nu_aoki(Re, Pr)
            # 液态金属 Nu 通常在 5-20 之间，过高会导致数值不稳定
            # 限制 Nu 最大为 500 (强湍流)，最小 4.36 (层流)
            return np.clip(Nu, 4.36, 500.0)

        self.coupler = FluidSolidCouple(
            name="Couple_Na_SS",
            fluid=self.channel,
            solid_boundary_region=self.solid.boundaries['left'],  # 内半径
            heated_perimeter=np.pi * self.D_inner,
            correlation_func=correlation_robust
        )

    def build_indexing(self):
        # 计数
        self.n_f_vol = self.channel.n_nodes
        self.n_junc_int = len(self.channel.internal_junctions)
        self.n_solid = self.solid.N

        idx = 0
        # Fluid State: P, h
        self.sl_f_P = slice(idx, idx + self.n_f_vol)
        idx += self.n_f_vol
        self.sl_f_h = slice(idx, idx + self.n_f_vol)
        idx += self.n_f_vol

        # Junctions: W (Internal + In + Out)
        self.sl_j_int = slice(idx, idx + self.n_junc_int)
        idx += self.n_junc_int
        self.idx_j_in = idx
        idx += 1
        self.idx_j_out = idx
        idx += 1

        # Boundary Volumes: P, h (Relaxation variables)
        self.idx_bv_P = slice(idx, idx + 2)
        idx += 2
        self.idx_bv_h = slice(idx, idx + 2)
        idx += 2

        # Solid State: T
        self.sl_s_T = slice(idx, idx + self.n_solid)
        idx += self.n_solid

        self.total_dim = idx

    def pack_state(self):
        y = np.zeros(self.total_dim)
        y[self.sl_f_P] = self.channel.pressure_vector
        y[self.sl_f_h] = [v.h for v in self.channel.volumes]
        y[self.sl_j_int] = [j.W for j in self.channel.internal_junctions]
        y[self.idx_j_in] = self.junc_in.W
        y[self.idx_j_out] = self.junc_out.W
        y[self.idx_bv_P] = [self.bv_in.P, self.bv_out.P]
        y[self.idx_bv_h] = [self.bv_in.h, self.bv_out.h]
        y[self.sl_s_T] = self.solid.T
        return y

    def unpack_state(self, y):
        # Fluid P, h
        P_vec = y[self.sl_f_P]
        h_vec = y[self.sl_f_h]
        for i, v in enumerate(self.channel.volumes):
            v.P = P_vec[i]
            v.h = h_vec[i]
            v.update_properties(self.mat_na)  # Update T, rho

        # Junction W
        W_int = y[self.sl_j_int]
        for i, j in enumerate(self.channel.internal_junctions):
            j.W = W_int[i]
            j.update_velocity()

        self.junc_in.W = y[self.idx_j_in]
        self.junc_in.update_velocity()
        self.junc_out.W = y[self.idx_j_out]
        self.junc_out.update_velocity()

        # Boundary
        bps = y[self.idx_bv_P]
        bhs = y[self.idx_bv_h]
        self.bv_in.P = bps[0]
        self.bv_in.h = bhs[0]
        self.bv_in.update_properties(self.mat_na)
        self.bv_out.P = bps[1]
        self.bv_out.h = bhs[1]
        self.bv_out.update_properties(self.mat_na)

        # Solid T
        self.solid.T[:] = y[self.sl_s_T]
        self.solid._update_properties()

    def rhs(self, t, y):
        self.unpack_state(y)

        # --- 物理耦合 ---
        self.channel.clear_sources()
        self.coupler.execute()

        dydt = np.zeros(self.total_dim)

        # 1. Fluid Mass & Energy
        dP_list, dh_list = [], []
        for v in self.channel.volumes:
            dp, dh = v.get_volume_derivatives(self.mat_na)
            dP_list.append(dp)
            dh_list.append(dh)
        dydt[self.sl_f_P] = dP_list
        dydt[self.sl_f_h] = dh_list

        # 2. Fluid Momentum
        dW_list = []
        for j in self.channel.internal_junctions:
            dW_list.append(j.get_momentum_derivative(self.mat_na))
        dydt[self.sl_j_int] = dW_list
        dydt[self.idx_j_in] = self.junc_in.get_momentum_derivative(self.mat_na)
        dydt[self.idx_j_out] = self.junc_out.get_momentum_derivative(self.mat_na)

        # 3. Boundary Relaxation
        dp_in, dh_in = self.bv_in.get_volume_derivatives(self.mat_na)
        dp_out, dh_out = self.bv_out.get_volume_derivatives(self.mat_na)
        dydt[self.idx_bv_P] = [dp_in, dp_out]
        dydt[self.idx_bv_h] = [dh_in, dh_out]

        # 4. Solid Conduction
        dydt[self.sl_s_T] = self.solid.get_derivatives(t, self.solid.T)

        return dydt


# ==============================================================================
# 执行脚本
# ==============================================================================
if __name__ == "__main__":
    # 实例化系统
    sys_sim = CoupledSystem(n_axial=20, n_radial=3)

    # ---------------------------------------------------------
    # Phase 1: 粗计算 (0 - 20s)
    # 目的: 快速建立流场和初步温度梯度，允许较大误差
    # ---------------------------------------------------------
    logger.info(">>> Phase 1: High Tolerance Calculation (0-20s)...")
    solver_coarse = NuclearODESolver(method='RK45', rtol=1e-1, atol=1e-1)

    y0 = sys_sim.pack_state()
    t_span1 = (0.0, 20.0)

    # max_step=0.1 依然建议加上，防止 RK45 在初期步长过大导致耦合滞后
    res1 = solver_coarse.solve(sys_sim.rhs, t_span1, y0, max_step=0.1)

    if not res1['success']:
        logger.error("Phase 1 Failed!")
        sys.exit(1)

    y_interim = res1['y'][:, -1]
    logger.info(">>> Phase 1 Completed.")

    # ---------------------------------------------------------
    # Phase 2: 精细续算 (20 - 40s)
    # 目的: 在已有解的基础上，收敛到高精度稳态
    # ---------------------------------------------------------
    logger.info(">>> Phase 2: Low Tolerance Calculation (20-40s)...")
    solver_fine = NuclearODESolver(method='RK45', rtol=1e-4, atol=1e-6)

    t_span2 = (20.0, 40.0)
    t_eval2 = np.linspace(20.0, 40.0, 101)  # 0.2s 间隔输出

    res2 = solver_fine.solve(sys_sim.rhs, t_span2, y_interim, t_eval=t_eval2, max_step=0.05)

    if not res2['success']:
        logger.error("Phase 2 Failed!")
        sys.exit(1)

    logger.info(">>> Phase 2 Completed.")

    # ==============================================================================
    # 结果分析与绘图
    # ==============================================================================
    y_final = res2['y'][:, -1]
    sys_sim.unpack_state(y_final)

    # 1. 守恒验证
    W_in = sys_sim.junc_in.W
    W_out = sys_sim.junc_out.W

    # 获取外表面积 (Right Boundary)
    area_outer = np.sum(sys_sim.solid.boundaries['right'].area)
    Q_input = sys_sim.q_flux * area_outer

    h_in = sys_sim.bv_in.h
    h_out = sys_sim.channel.volumes[-1].h
    Q_fluid = W_out * (h_out - h_in)

    err_Q = abs(Q_input - Q_fluid) / Q_input * 100

    print("-" * 50)
    print(f"VERIFICATION (t=40s)")
    print("-" * 50)
    print(f"Flow Rate: In={W_in:.4f}, Out={W_out:.4f} kg/s")
    print(f"Heat Input (Wall): {Q_input:.2f} W")
    print(f"Heat Fluid (Out-In): {Q_fluid:.2f} W")
    print(f"Energy Balance Error: {err_Q:.2f} %")

    # 2. 绘图
    z_coords = [v.z_coordinate for v in sys_sim.channel.volumes]
    T_fluid = sys_sim.channel.temperature_vector

    # 提取 2D 固体温度
    # Mesh2D flattened: Layer by Layer.
    # Layer 0 (Indices 0-19) -> Inner
    # Layer 1 (Indices 20-39) -> Middle
    # Layer 2 (Indices 40-59) -> Outer
    nx = sys_sim.n_radial
    ny = sys_sim.n_axial

    T_solid = sys_sim.solid.T
    T_inner = T_solid[0:ny]  # 第一层 (内壁)
    T_outer = T_solid[(nx - 1) * ny: nx * ny]  # 最后一层 (外壁)

    plt.figure(figsize=(10, 6))
    plt.plot(z_coords, T_outer, 'r-s', label='Solid Outer (Heated)', markersize=4)
    plt.plot(z_coords, T_inner, 'm-^', label='Solid Inner (Interface)', markersize=4)
    plt.plot(z_coords, T_fluid, 'b-o', label='Fluid Sodium', linewidth=2)

    plt.title(f'Sodium Pipe Flow (2D Solid Conduction)\nTime=40s, Q"={sys_sim.q_flux}W/m2')
    plt.xlabel('Axial Position (m)')
    plt.ylabel('Temperature (K)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig('test_result_20_40s.png')
    plt.show()
