import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import logging

# 引入核心组件
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Solvers.Hydrodynamics.Components import FluidChannel
from Solvers.HeatConduction.Mesh import Mesh1D
from Solvers.HeatConduction.HeatConduction import HeatConduction1D
from Solvers.Couplers import FluidSolidCouple
from MathSolvers.solver_module import NuclearODESolver

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


class SimpleLumpedSystem:
    def __init__(self):
        self.h_old = 0.0

        # --- 1. 准备物性 ---
        self.mat_na = Sodium()
        self.mat_ss = AusteniticStainlessSteel()

        # --- 2. 建立流体 (静止的钠池) ---
        # 1个节点, 长度1m, 面积0.01m2 -> 体积 0.01 m3
        self.fluid = FluidChannel(
            name="SodiumTank", n_nodes=1, total_length=1.0,
            flow_area=0.01, hydraulic_diam=0.1,
            initial_P=116000, initial_T=600.0, material=self.mat_na
        )
        self.fluid.initialize_state()

        # 赋予一个虚拟流速以便计算 Re 数 (假设有搅拌)
        # Hack: 直接给 junction 列表塞一个假对象，或者在关联式里写死 Nu
        # 这里我们选择在耦合器的 correlation_func 里写死 Nu=10，模拟自然对流

        # --- 3. 建立固体 (不锈钢块) ---
        # 1个节点 (集总参数), 厚度 1cm
        mesh_solid = Mesh1D(total_dim=0.01, n_volumes=1, geometry_type='cartesian')
        self.solid = HeatConduction1D(mesh_solid, self.mat_ss, initial_temp=650.0)

        # --- 4. 建立耦合 ---
        # 将流体与固体的 'inner' (左) 边界连接
        # 换热面积 = 1.0 m (长度) * 0.4 m (假设湿周) = 0.4 m2
        self.heated_perimeter = 0.4
        self.heat_area = self.heated_perimeter * self.fluid.node_length  # 0.4 m2

        # 强制设置固体边界的面积，确保物理量守恒计算正确
        self.solid.boundaries['inner'].area[:] = self.heat_area

        # 定义关联式 (模拟常数换热系数)
        # 返回 Nu。 h = Nu * k / D. 设 Nu=50
        def const_nu_correlation(Re, Pr, ratio):
            return 50

        self.coupler = FluidSolidCouple(
            name="TankCoupler",
            fluid=self.fluid,
            solid_boundary_region=self.solid.boundaries['inner'],
            heated_perimeter=self.heated_perimeter,
            correlation_func=const_nu_correlation
        )

        # 打印初始状态
        # 固体质量 = rho * V = 7900 * (0.01 * 0.4) approx 30kg
        # 流体质量 = rho * V = 850 * 0.01 approx 8.5kg
        print(f"Init Fluid T: {self.fluid.volumes[0].T:.2f} K")
        print(f"Init Solid T: {self.solid.T[0]:.2f} K")

    def rhs(self, t, y):
        """
        状态向量 y: [h_fluid, T_solid]
        """
        h_fluid_curr = y[0]
        T_solid_curr = y[1]

        # if self.h_old < h_fluid_curr:
        #     print(t, "!!!0")
        # else:
        #     print(t, "!!!1")

        self.h_old = h_fluid_curr

        # 1. 更新对象状态
        # 流体: 更新 h, 然后反推 T, rho, mu
        vol = self.fluid.volumes[0]
        vol.h = h_fluid_curr
        vol.update_properties(self.mat_na)  # 关键：更新 vol.T

        # 固体: 更新 T
        self.solid.T[0] = T_solid_curr
        self.solid._update_properties()

        # 2. 执行耦合 (计算 Q)
        # 清除旧源项
        self.fluid.clear_sources()
        # 执行耦合逻辑 (计算 h, T_surf, Q_wall)
        self.coupler.execute()

        # print(t)

        # if vol.T > self.solid.T[0]:
        #     print(t)
        # else:
        #     print("2")
        # if abs(vol.T - self.solid.T[0]) < 5.0:
        #     print("3")

        # 3. 计算导数

        # --- 流体能量方程 ---
        # rho * V * dh/dt = Q_total
        # Q_total = Q_wall (Explicit) - Implicit * T_fluid
        # 注意: execute() 已经把 explicit 和 implicit 填入 vol 了
        Q_fluid_total = vol.Q_wall + vol.Q_vol - vol.implicit_coeff * vol.T

        mass_fluid = vol.rho * vol.vol
        dydt_P, dh_dt = vol.get_volume_derivatives(vol.material)

        # if dh_dt < 0:
            # print("4")

        # --- 固体导热方程 ---
        # m * Cp * dT/dt = -Q_out
        # HeatConduction1D.get_derivatives 会自动处理边界通量
        # 我们这里只有一个节点，可以直接利用 boundary flux
        # get_derivatives 内部调用:
        #   compute_internal_resistance()
        #   update_boundaries_state() -> 将 T_solid 推送给 Boundary
        #   compute_fluxes() -> 从 Boundary 获取 Q (由耦合器写入的)

        dT_solid_dt = self.solid.get_derivatives(t, np.array([T_solid_curr]))

        return [dh_dt, dT_solid_dt[0]]


def run_simulation():
    system = SimpleLumpedSystem()

    # 初始状态 y0 = [h_fluid, T_solid]
    y0 = [system.fluid.volumes[0].h, system.solid.T[0]]

    # 时间积分 0 -> 200s
    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 1001)

    print("Starting integration...")
    # solver = NuclearODESolver(method='BDF', rtol=1e-3, atol=1e-6)
    sol = solve_ivp(system.rhs, t_span, y0, t_eval=t_eval, method='RK45', max_step=0.1)
    # sol = solver.solve(system.rhs, t_span, y0, t_eval)

    if sol.success:
        print("Integration successful!")

        # --- 后处理与绘图 ---
        h_res = sol.y[0]
        T_solid_res = sol.y[1]
        t_res = sol.t

        # 将流体焓 h 转换为温度 T 以便绘图
        T_fluid_res = []
        for h in h_res:
            T_val = system.mat_na.temperature_from_enthalpy(h, 1.0e5)
            T_fluid_res.append(T_val)

        plt.figure(figsize=(8, 6))
        plt.plot(t_res, T_solid_res, 'r-', label='Solid (Steel) Temp', linewidth=2)
        plt.plot(t_res, T_fluid_res, 'b-', label='Fluid (Sodium) Temp', linewidth=2)

        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (K)')
        plt.title('Transient Coupled Heating: Steel Block in Sodium Tank')
        plt.legend()
        plt.grid(True)
        plt.savefig('simple_transient_result.png')
        print("Result saved to 'simple_transient_result.png'")
        plt.show()
    else:
        print("Integration failed:", sol.message)


if __name__ == "__main__":
    run_simulation()
