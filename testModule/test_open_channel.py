import numpy as np
import matplotlib.pyplot as plt
from typing import List

# 引入你的模块 (假设路径正确)
from Materials.Fluids.Sodium import Sodium
from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, FluidVolume
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume
from MathSolvers.solver_module import NuclearODESolver

# ==============================================================================
# 1. 定义压力控制曲线 (Sigmoid Ramp)
# ==============================================================================
# 初始设置
P_OUT_FIXED = 115000.0  # [Pa] 出口恒定压力
DELTA_P_START = 500.0  # [Pa] 初始压差
P_IN_START = P_OUT_FIXED + DELTA_P_START

# 瞬态设置
T_RAMP_START = 15.0  # [s] 开始变化的时间
RAMP_DURATION = 20.0  # [s] 变化持续时间
DELTA_P_END = DELTA_P_START * 2.0  # [Pa] 最终压差 (翻倍)


def pressure_schedule(t: float) -> float:
    """
    计算当前时刻的进口压力 P_in(t)。
    在 T_RAMP_START 后，使用 Sigmoid 函数将压差平滑翻倍。
    """
    # 目标：从 P_in_start 增加到 P_in_end
    p_start = P_IN_START
    p_end = P_OUT_FIXED + DELTA_P_END

    # Sigmoid 参数
    # 中心点设在变化区间的中间
    t_center = T_RAMP_START + (RAMP_DURATION / 2.0)
    # 宽度参数: 20s 的区间，width=3.5 能保证两端平滑且覆盖约 95% 的变化
    width = 3.5

    # 计算
    arg = (t - t_center) / width
    arg = np.clip(arg, -50.0, 50.0)  # 防止溢出

    # Sigmoid: Start + (End - Start) / (1 + exp(-x))
    # 注意这里 arg 是 (t-center)，所以用正号在分母里是 1 + exp(-arg) 写法，或者如下：
    p_current = p_start + (p_end - p_start) / (1.0 + np.exp(-arg))

    return p_current


# ==============================================================================
# 2. 系统管理器 (支持变压力边界)
# ==============================================================================
class SystemManager:
    def __init__(self, volumes: List[FluidVolume], junctions: List[FlowJunction]):
        self.volumes = volumes
        self.junctions = junctions
        self.n_vol = len(volumes)
        self.n_junc = len(junctions)
        self.dim = 2 * self.n_vol + self.n_junc

        # [关键] 自动寻找入口边界体积 (Inlet Boundary Volume)
        # 我们假设它是列表中第一个 BoundaryVolume，或者通过名字/位置识别
        # 在此算例中，我们明确知道它是 volumes[0] (因为构建顺序)
        self.bc_inlet: BoundaryVolume = None

        if isinstance(volumes[0], BoundaryVolume):
            self.bc_inlet = volumes[0]
        else:
            # 搜索 fallback
            for v in volumes:
                if isinstance(v, BoundaryVolume) and "Inlet" in v.name:
                    self.bc_inlet = v
                    break

        if self.bc_inlet is None:
            raise ValueError("SystemManager 无法找到 Inlet BoundaryVolume!")

    def get_initial_state(self) -> np.ndarray:
        y = np.zeros(self.dim)
        y[0:self.n_vol] = [v.P for v in self.volumes]
        y[self.n_vol:2 * self.n_vol] = [v.h for v in self.volumes]
        y[2 * self.n_vol:] = [j.W for j in self.junctions]
        return y

    def dydt_function(self, t, y):
        """ODE 右端项"""

        # --- [Step 0] 动态更新边界压力 ---
        # 计算当前时刻的目标进口压力
        target_P_in = pressure_schedule(t)

        # 强制更新边界对象的状态 (这是 BoundaryVolume.py 中 set_state 的作用)
        # 注意: 这不涉及 dydt，而是直接修改代数方程的源项
        self.bc_inlet.set_state(P=target_P_in)

        # --- [Step 1] 解包状态 ---
        Ps = y[0:self.n_vol]
        hs = y[self.n_vol:2 * self.n_vol]
        Ws = y[2 * self.n_vol:]

        # 更新 Volumes
        for i, vol in enumerate(self.volumes):
            vol.P = Ps[i]
            vol.h = hs[i]
            if isinstance(vol, BoundaryVolume):
                vol.update_properties(vol.material)
            else:
                vol.T = vol.material.temperature_from_enthalpy(vol.h, vol.P)
                vol.rho = vol.material.density(vol.T, vol.P)
                vol.mu = vol.material.viscosity(vol.T, vol.P)

        # 更新 Junctions
        for i, junc in enumerate(self.junctions):
            junc.W = Ws[i]
            junc.update_velocity()

        # --- [Step 2] 计算导数 ---
        dydt = np.zeros(self.dim)

        for i, vol in enumerate(self.volumes):
            # 获取材料 (BoundaryVolume 自带，FluidVolume 需确保已分配)
            mat = vol.material
            dP, dh = vol.get_volume_derivatives(mat)
            dydt[i] = dP
            dydt[self.n_vol + i] = dh

        for i, junc in enumerate(self.junctions):
            mat = junc.from_vol.material
            if mat is None: mat = junc.to_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[2 * self.n_vol + i] = dW

        return dydt


# ==============================================================================
# 3. 主运行程序 (三阶段串联)
# ==============================================================================
def run_pressure_ramp_test():
    print("🚀 启动模拟: 变压差瞬态测试 (dP -> 2*dP)")

    # --- A. 构建物理模型 ---
    sodium = Sodium()
    L_pipe = 0.5
    D_pipe = 0.03
    Area = np.pi * (D_pipe / 2) ** 2
    N_nodes = 10
    T_in = 743.0

    # 通道
    channel = FluidChannel("Pipe", N_nodes, L_pipe, Area, D_pipe,
                           initial_P=P_OUT_FIXED, initial_T=T_in, material=sodium)

    # 边界 (Inlet P 初始设为 P_IN_START)
    bc_inlet = BoundaryVolume("Inlet_BC", sodium, P=P_IN_START, T=T_in, flow_area=Area, hydraulic_diam=D_pipe)
    bc_outlet = BoundaryVolume("Outlet_BC", sodium, P=P_OUT_FIXED, T=T_in, flow_area=Area, hydraulic_diam=D_pipe)

    # 线性初始化管道内部压力
    P_dist = np.linspace(P_IN_START, P_OUT_FIXED, N_nodes)
    for i, v in enumerate(channel.volumes):
        v.P = P_dist[i]
        v.update_properties(sodium)

    # 连接 (使用普通 FlowJunction，因为流量由压差驱动)
    junc_in = FlowJunction("J_In", bc_inlet, channel.volumes[0], flow_area=Area)
    junc_out = FlowJunction("J_Out", channel.volumes[-1], bc_outlet, flow_area=Area)

    # 估算初值流量 (用于加速收敛)
    rho_est = sodium.density(T_in, P_OUT_FIXED)
    # dP = f * L/D * 0.5 * rho * v^2 -> v ~ sqrt(dP)
    # 粗略估计 v ~ 1.0 m/s 对应 dP ~ 500Pa (非常粗略)
    W_guess = 2.0  # 给一个非零初值即可
    for j in [junc_in] + channel.internal_junctions + [junc_out]:
        j.W = W_guess

    # 管理器
    all_vols = [bc_inlet] + channel.volumes + [bc_outlet]
    all_juncs = [junc_in] + channel.internal_junctions + [junc_out]
    manager = SystemManager(all_vols, all_juncs)

    # 容差配置
    atol_fine = [1.0] * manager.n_vol + [10.0] * manager.n_vol + [1e-4] * manager.n_junc

    # ==========================================================================
    # Phase 1: 粗略初始化 (0 - 5s)
    # ==========================================================================
    print("\n[Phase 1] 粗略初始化 (t=0->5s, rtol=1e-1)...")
    solver_p1 = NuclearODESolver(method='BDF', rtol=1e-1, atol=atol_fine)  # atol保持细致以免W跑飞
    y0 = manager.get_initial_state()

    res1 = solver_p1.solve(manager.dydt_function, (0.0, 5.0), y0)
    if not res1['success']:
        return print("❌ Phase 1 失败")

    # ==========================================================================
    # Phase 2: 高精度稳态 (5 - 15s)
    # ==========================================================================
    print("[Phase 2] 高精度稳态 (t=5->15s, rtol=1e-4)...")
    solver_p2 = NuclearODESolver(method='BDF', rtol=1e-4, atol=atol_fine)
    y1_final = res1['y'][:, -1]

    res2 = solver_p2.solve(manager.dydt_function, (5.0, 15.0), y1_final)
    if not res2['success']:
        return print("❌ Phase 2 失败")

    # 记录稳态流量用于后续对比
    W_steady_start = res2['y'][2 * manager.n_vol, -1]
    print(f"   >> 初始稳态流量: {W_steady_start:.4f} kg/s")

    # ==========================================================================
    # Phase 3: 压力瞬态 (15 - 45s)
    # ==========================================================================
    print("[Phase 3] 压力瞬态 (t=15->45s, rtol=1e-4)...")
    # 此时 dydt_function 会自动根据 t > 15 调用 pressure_schedule 改变压力
    solver_p3 = NuclearODESolver(method='BDF', rtol=1e-4, atol=atol_fine)
    y2_final = res2['y'][:, -1]

    # 稍微延长一点时间以便观察新稳态
    t_eval_p3 = np.linspace(15.0, 60.0, 300)
    res3 = solver_p3.solve(manager.dydt_function, (15.0, 60.0), y2_final, t_eval=t_eval_p3)
    if not res3['success']: return print("❌ Phase 3 失败")

    W_steady_end = res3['y'][2 * manager.n_vol, -1]
    print(f"   >> 最终稳态流量: {W_steady_end:.4f} kg/s")
    ratio = W_steady_end / W_steady_start
    print(f"   >> 流量倍率: {ratio:.3f} (理论值 sqrt(2) ≈ 1.414)")

    # ==========================================================================
    # 数据拼接与绘图
    # ==========================================================================
    t_all = np.concatenate((res1['t'], res2['t'], res3['t']))
    y_all = np.concatenate((res1['y'], res2['y'], res3['y']), axis=1)

    plot_pressure_transient(t_all, y_all, manager)


def plot_pressure_transient(t, y, manager):
    # 提取数据
    n_vol = manager.n_vol
    idx_p_in = 0  # Inlet BC Index
    idx_w_in = 2 * n_vol

    P_boundary_in = y[idx_p_in, :]
    W_in = y[idx_w_in, :]

    # 理论压差曲线 (用于核对)
    Delta_P_curve = P_boundary_in - P_OUT_FIXED

    # 绘图
    plt.figure(figsize=(10, 10))

    # 1. 压差输入
    plt.subplot(2, 1, 1)
    plt.plot(t, Delta_P_curve, 'k-', label='Boundary Delta P', linewidth=2)
    plt.axvline(x=15.0, color='gray', linestyle='--')
    plt.axvline(x=35.0, color='gray', linestyle='--')
    plt.text(7.5, DELTA_P_START, "Steady 1", ha='center', color='blue')
    plt.text(25.0, DELTA_P_START * 1.5, "Ramp (20s)", ha='center', color='blue')
    plt.text(40.0, DELTA_P_END, "Steady 2", ha='center', color='blue')
    plt.ylabel('Delta P [Pa]')
    plt.title('Forcing Function: Pressure Drop Ramp')
    plt.grid(True)
    plt.legend()

    # 2. 流量响应
    plt.subplot(2, 1, 2)
    plt.plot(t, W_in, 'b-', label='Mass Flow Rate', linewidth=2)

    # 绘制理论倍率线
    W_start = np.mean(W_in[(t > 5) & (t < 15)])
    W_theoretical_end = W_start * 1.486 # np.sqrt(2.0)
    plt.axhline(y=W_theoretical_end, color='r', linestyle=':', label='Theoretical Limit (x1.414)')

    plt.axvline(x=15.0, color='gray', linestyle='--')
    plt.ylabel('Flow [kg/s]')
    plt.title('System Response: Mass Flow Rate')
    plt.xlabel('Time [s]')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_pressure_ramp_test()