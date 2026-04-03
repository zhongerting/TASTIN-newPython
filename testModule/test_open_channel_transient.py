import numpy as np
import matplotlib.pyplot as plt
from typing import List

# 确保这些模块在你的路径中
from Materials.Fluids.Sodium import Sodium
from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, FluidVolume
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume, InletJunction
from MathSolvers.solver_module import NuclearODESolver

# ==============================================================================
# 1. 定义流量控制曲线 (平滑版 - 消除压力震荡)
# ==============================================================================
T_TRIP = 20.0  # 定义动作开始的大致时间区间


def flow_schedule(t: float) -> float:
    """
    使用 Sigmoid 函数生成的平滑流量曲线。
    特点：在 T_TRIP 附近平滑过渡，消除导数突变引起的压力尖峰。
    """
    # --- 1. 流量参数 ---
    w_start = 0.45  # [kg/s] 初始稳态流量
    w_end = 0.225  # [kg/s] 最终惰转流量

    # return 0.33

    # --- 2. 形状参数 (关键) ---
    # t_center: 下降的中点时刻 (流量降到一半的时间)
    # 你原本是 20s 开始降，tau=10s (很慢)。
    # 为了模拟类似的缓慢过程，我们将中点设在 20 + 8 = 28s 左右
    t_center = T_TRIP + 20.0

    # width: 曲线的"宽度"或平缓程度 (类似时间常数)
    # width 越大，下降越平缓，拐弯越圆滑。
    # 为了匹配你原来 tau=10 的慢速衰减，这里设为 3.0 ~ 4.0
    width = 3.5

    # --- 3. Sigmoid 计算公式 ---
    # 公式: f(t) = End + (Start - End) / (1 + exp( (t - t_center) / width ))

    # [安全措施] 防止 exp 计算溢出 (当 t 很大或很小时)
    # 将指数的输入限制在 [-50, 50] 之间，这对结果精度没有影响
    arg = (t - t_center) / width
    arg = np.clip(arg, -50.0, 50.0)

    w_target = w_end + (w_start - w_end) / (1.0 + np.exp(arg))

    return w_target


# ==============================================================================
# 2. 系统管理器
# ==============================================================================
class SystemManager:
    def __init__(self, volumes: List[FluidVolume], junctions: List[FlowJunction]):
        self.volumes = volumes
        self.junctions = junctions
        self.n_vol = len(volumes)
        self.n_junc = len(junctions)
        self.dim = 2 * self.n_vol + self.n_junc

        # 自动寻找入口连接
        self.inlet_junc: InletJunction = None
        for j in junctions:
            if isinstance(j, InletJunction):
                self.inlet_junc = j
                break
        if self.inlet_junc is None:
            raise ValueError("系统必须包含一个 InletJunction 用于控制流量！")

    def get_initial_state(self) -> np.ndarray:
        y = np.zeros(self.dim)
        y[0:self.n_vol] = [v.P for v in self.volumes]
        y[self.n_vol:2 * self.n_vol] = [v.h for v in self.volumes]
        y[2 * self.n_vol:] = [j.W for j in self.junctions]
        return y

    def dydt_function(self, t, y):
        """ODE 右端项函数"""
        # [Step 0] 更新控制目标
        target_w = flow_schedule(t)
        self.inlet_junc.set_flow_rate(target_w)

        # [Step 1] 解包求解器状态
        Ps = y[0:self.n_vol]
        hs = y[self.n_vol:2 * self.n_vol]
        Ws = y[2 * self.n_vol:]

        # 更新控制体
        for i, vol in enumerate(self.volumes):
            vol.P = Ps[i]
            vol.h = hs[i]
            if isinstance(vol, BoundaryVolume):
                vol.update_properties(vol.material)
            else:
                vol.T = vol.material.temperature_from_enthalpy(vol.h, vol.P)
                vol.rho = vol.material.density(vol.T, vol.P)
                vol.mu = vol.material.viscosity(vol.T, vol.P)

        # 更新连接管
        for i, junc in enumerate(self.junctions):
            junc.W = Ws[i]
            junc.update_velocity()

        # [Step 2] 计算物理导数
        dydt = np.zeros(self.dim)

        for i, vol in enumerate(self.volumes):
            dP, dh = vol.get_volume_derivatives(vol.material)
            dydt[i] = dP
            dydt[self.n_vol + i] = dh

        for i, junc in enumerate(self.junctions):
            mat = junc.from_vol.material
            if mat is None: mat = junc.to_vol.material
            dW = junc.get_momentum_derivative(mat)
            dydt[2 * self.n_vol + i] = dW

        return dydt


# ==============================================================================
# 3. 主运行程序 (两阶段计算)
# ==============================================================================
def run_transient_test():
    print("🚀 启动模拟: 钠回路泵惰转瞬态 (Two-Phase Calculation)")

    # --- A. 参数与组件 ---
    sodium = Sodium()
    L_pipe = 0.6
    D_pipe = 0.03
    Area = np.pi * (D_pipe / 2) ** 2
    N_nodes = 20
    W_nominal = 0.33
    P_out = 116000.0
    T_in = 743.0

    channel = FluidChannel("Pipe", N_nodes, L_pipe, Area, D_pipe,
                           initial_P=P_out, initial_T=T_in, material=sodium)

    # 线性压力初始化
    P_in_guess = P_out + 70.0
    P_dist = np.linspace(P_in_guess, P_out, N_nodes)
    for i, v in enumerate(channel.volumes):
        v.P = P_dist[i]
        v.update_properties(sodium)

    bc_inlet = BoundaryVolume("Inlet", sodium, P=P_in_guess, T=T_in, flow_area=Area, hydraulic_diam=D_pipe)
    bc_outlet = BoundaryVolume("Outlet", sodium, P=P_out, T=T_in, flow_area=Area, hydraulic_diam=D_pipe)

    junc_in = InletJunction("J_In", bc_inlet, channel.volumes[0], W_initial=W_nominal)
    junc_out = FlowJunction("J_Out", channel.volumes[-1], bc_outlet, flow_area=Area)

    for j in channel.internal_junctions: j.W = W_nominal
    junc_out.W = W_nominal

    # --- B. 系统管理器 ---
    all_vols = [bc_inlet] + channel.volumes + [bc_outlet]
    all_juncs = [junc_in] + channel.internal_junctions + [junc_out]
    manager = SystemManager(all_vols, all_juncs)

    # 构造 y0
    y0 = manager.get_initial_state()
    # 强制 y0 入口流量对齐
    y0[2 * manager.n_vol] = flow_schedule(0.0)

    # 容差设置 (列表式)
    atol_list = [1.0] * manager.n_vol + [10.0] * manager.n_vol + [1e-4] * manager.n_junc

    # ==========================================================================
    # Phase 1: 稳态初始化 (Initialization Run)
    # ==========================================================================
    print(f"\n[Phase 1] 稳态初始化 (t=0.0 -> {T_TRIP}s)...")
    print(">> 使用大容差 (rtol=1e-1) 消除初值不平衡")

    # 使用大容差求解器
    solver_init = NuclearODESolver(method='BDF', rtol=1e-1, atol=atol_list)

    # 运行
    res_init = solver_init.solve(manager.dydt_function, (0.0, T_TRIP), y0)

    if not res_init['success']:
        print(f"❌ 初始化失败: {res_init['message']}")
        return

    # ==========================================================================
    # Phase 2: 瞬态计算 (Transient Run)
    # ==========================================================================
    print(f"\n[Phase 2] 瞬态惰转 (t={T_TRIP}s -> 15.0s)...")
    print(">> 使用标准容差 (rtol=1e-4) 捕捉流量下降")

    # 提取 Phase 1 的最终状态作为 Phase 2 的起点
    y_steady = res_init['y'][:, -1]

    # 稍微收紧一点容差 (可选，保持 1e-1 也是安全的，这里演示收紧到 1e-2)
    solver_trans = NuclearODESolver(method='BDF', rtol=1e-3, atol=atol_list)

    t_end = 100.0
    t_eval_trans = np.linspace(T_TRIP, t_end, 200)

    res_trans = solver_trans.solve(manager.dydt_function, (T_TRIP, t_end), y_steady, t_eval=t_eval_trans)

    if res_trans['success']:
        print("✅ 瞬态计算成功!")

        # 拼接数据用于绘图
        # 注意: res_init 的时间点可能很少，res_trans 我们指定了 t_eval
        t_all = np.concatenate((res_init['t'], res_trans['t']))
        y_all = np.concatenate((res_init['y'], res_trans['y']), axis=1)

        t_figure = t_all[10:]
        y_figure = y_all[:,10:]

        plot_results(t_figure, y_figure, manager)
    else:
        print(f"❌ 瞬态计算失败: {res_trans['message']}")


def plot_results(t, y, manager):
    # ==========================================================================
    # [关键修改] 数据切片：只保留 t > 10.0s 的数据
    # ==========================================================================
    start_time_display = 1

    # 找到时间数组中大于 start_time_display 的所有索引
    mask = t >= start_time_display

    # 使用布尔掩码 (Mask) 过滤时间和状态数据
    t_plot = t[mask]
    y_plot = y[:, mask]  # 注意 y 是二维数组，列对应时间

    # 如果数据为空 (例如模拟时长不足10s)，则回退到显示全部
    if len(t_plot) == 0:
        print(f"⚠️ 警告: 模拟总时长不足 {start_time_display}s，显示全部数据。")
        t_plot = t
        y_plot = y

    t = t_plot
    y = y_plot

    # --- 1. 提取流量 ---
    # n_vol 是包含边界在内的总节点数
    n_vol = manager.n_vol

    # InletJunction 是第一个连接 (连接 BC_Inlet 和 Pipe_Vol_1)
    idx_w_in = 2 * n_vol
    # OutletJunction 是最后一个连接
    idx_w_out = 2 * n_vol + manager.n_junc - 1

    W_in = y[idx_w_in, :]
    W_out = y[idx_w_out, :]
    W_set = [flow_schedule(ti) for ti in t]

    # --- 2. [关键修正] 提取压差 ---
    # 我们要看的是管道(Channel)两端的压差，而不是边界(Boundary)的压差
    # y 的前 n_vol 项是压力
    # 结构: [BC_Inlet, Pipe_1, Pipe_2, ..., Pipe_N, BC_Outlet]

    # 索引 0: BC_Inlet (恒定压力)
    # 索引 1: Pipe_Vol_0 (管道入口，压力随流量变化) <-- 选这个！
    idx_p_pipe_in = 1

    # 索引 n_vol-1: BC_Outlet (恒定压力)
    # 索引 n_vol-2: Pipe_Vol_N-1 (管道出口) <-- 选这个！
    idx_p_pipe_out = n_vol - 2

    P_pipe_in = y[idx_p_pipe_in, :]
    P_pipe_out = y[idx_p_pipe_out, :]

    # 计算管道本体的压差
    dP_channel = P_pipe_in - P_pipe_out

    # (可选) 边界压差，用于对比验证它是直线的
    P_bc_in = y[0, :]
    P_bc_out = y[n_vol - 1, :]
    dP_boundary = P_bc_in - P_bc_out

    # --- 3. 绘图 ---
    plt.figure(figsize=(10, 10))

    # 子图1: 流量
    plt.subplot(3, 1, 1)
    plt.plot(t, W_set, 'k--', label='Set Point', alpha=0.6)
    plt.plot(t, W_in, 'b-', label='Inlet Flow', linewidth=2)
    plt.plot(t, W_out, 'r:', label='Outlet Flow', linewidth=2)
    plt.axvline(x=T_TRIP, color='gray', linestyle='--')
    plt.ylabel('Flow [kg/s]')
    plt.title('Flow Coast-down')
    plt.legend()
    plt.grid(True)

    # 子图2: 管道压差 (应随流量下降)
    plt.subplot(3, 1, 2)
    plt.plot(t, dP_channel, 'g-', label='Channel dP (Dynamic)', linewidth=2)
    plt.axvline(x=T_TRIP, color='gray', linestyle='--')
    plt.ylabel('Channel dP [Pa]')
    plt.title('Channel Pressure Drop Response (Correct)')
    plt.legend()
    plt.grid(True)

    # 子图3: 边界/节点压力绝对值 (用于诊断)
    plt.subplot(3, 1, 3)
    plt.plot(t, P_pipe_in, 'b-', label='Pipe Inlet Node P')
    plt.plot(t, P_bc_in, 'k--', label='Boundary Inlet P (Fixed)')
    plt.ylabel('Absolute Pressure [Pa]')
    plt.title('Inlet Pressure: Node vs Boundary')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_transient_test()