import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import logging
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

# 引入核心组件
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume, InletJunction
from Solvers.HeatConduction.Mesh import Mesh1D
from Solvers.HeatConduction.HeatConduction import HeatConduction1D
from Solvers.HeatConduction.Boundary import BoundaryRegion
from Solvers.Couplers import FluidSolidCouple
from MathSolvers.solver_module import NuclearODESolver

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


class SimpleLumpedSystem:
    def __init__(self, sync_boundary_before_couple: bool = False):
        self.h_old = 0.0
        self.sync_boundary_before_couple = sync_boundary_before_couple

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
        print(f"Boundary sync before coupler: {self.sync_boundary_before_couple}")

    def _sync_solid_boundary(self, t: float):
        """Push the latest solid state to BoundaryRegion before coupler reads it."""
        self.solid._update_properties()
        self.solid._compute_internal_resistance()
        self.solid._update_boundaries_state(current_time=t)

    def _set_state(self, h_fluid: float, T_solid: float, t: float):
        vol = self.fluid.volumes[0]
        vol.h = h_fluid
        vol.update_properties(self.mat_na)
        self.solid.T[0] = T_solid
        self.solid._update_properties()
        if self.sync_boundary_before_couple:
            self._sync_solid_boundary(t)
        return vol

    def evaluate_interface_audit(self, t: float, h_fluid: float, T_solid: float):
        """Evaluate instantaneous fluid/solid interface energy balance."""
        vol = self._set_state(h_fluid, T_solid, t)
        self.fluid.clear_sources()
        self.coupler.execute()

        Q_fluid = vol.Q_wall + vol.Q_vol - vol.implicit_coeff * vol.T
        dT_solid_dt = self.solid.get_derivatives(t, np.array([T_solid]))[0]
        Q_solid_storage = self.solid.thermal_capacitance[0] * dT_solid_dt
        Q_solid_boundary = self.solid.boundaries['inner'].current_flux[0]

        return {
            "Q_fluid": float(Q_fluid),
            "Q_solid_boundary": float(Q_solid_boundary),
            "Q_solid_storage": float(Q_solid_storage),
            "interface_residual": float(Q_fluid + Q_solid_boundary),
            "storage_residual": float(Q_fluid + Q_solid_storage),
        }

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
        vol = self._set_state(h_fluid_curr, T_solid_curr, t)

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


def _as_float(value):
    return float(np.asarray(value).reshape(-1)[0])


class FlowThroughEnthalpyAudit:
    """Flow-through channel audit: Q_wall = W*(h_out-h_in) + storage."""

    def __init__(
            self,
            n_nodes: int = 12,
            length: float = 1.0,
            diameter: float = 0.03,
            mass_flow: float = 0.20,
            inlet_temp: float = 600.0,
            wall_temp: float = 650.0,
            pressure: float = 116000.0):
        self.n_nodes = n_nodes
        self.length = length
        self.diameter = diameter
        self.mass_flow = mass_flow
        self.inlet_temp = inlet_temp
        self.wall_temp = wall_temp
        self.pressure = pressure

        self.mat_na = Sodium()
        self.flow_area = np.pi * (0.5 * diameter) ** 2
        self.perimeter = np.pi * diameter

        self.inlet = BoundaryVolume(
            name="Audit_Inlet",
            material=self.mat_na,
            P=pressure,
            T=inlet_temp,
            flow_area=self.flow_area,
            hydraulic_diam=diameter,
        )
        self.outlet = BoundaryVolume(
            name="Audit_Outlet",
            material=self.mat_na,
            P=pressure,
            T=inlet_temp,
            flow_area=self.flow_area,
            hydraulic_diam=diameter,
        )

        self.channel = FluidChannel(
            name="Audit_Channel",
            n_nodes=n_nodes,
            total_length=length,
            flow_area=self.flow_area,
            hydraulic_diam=diameter,
            initial_P=pressure,
            initial_T=inlet_temp,
            material=self.mat_na,
        )
        self.channel.initialize_state()

        self.junc_in = InletJunction(
            name="Audit_Junc_In",
            from_vol=self.inlet,
            to_vol=self.channel.volumes[0],
            W_initial=mass_flow,
        )
        self.junc_out = FlowJunction(
            name="Audit_Junc_Out",
            from_vol=self.channel.volumes[-1],
            to_vol=self.outlet,
            flow_area=self.flow_area,
        )
        self.junc_out.W = mass_flow
        self.channel.inlet_junction = self.junc_in
        self.channel.outlet_junction = self.junc_out

        for junc in self.channel.internal_junctions:
            junc.W = mass_flow
            junc.update_velocity()
        self.junc_in.update_velocity()
        self.junc_out.update_velocity()

        area_nodes = np.full(n_nodes, self.perimeter * self.channel.node_length)
        self.wall_boundary = BoundaryRegion(shape=(n_nodes,), area_array=area_nodes)

        def constant_nu(Re, Pr, ratio):
            return np.full_like(np.asarray(Re, dtype=float), 20.0)

        self.coupler = FluidSolidCouple(
            name="Audit_Flow_Wall_Coupler",
            fluid=self.channel,
            solid_boundary_region=self.wall_boundary,
            heated_perimeter=self.perimeter,
            correlation_func=constant_nu,
        )

        self.h_inlet = float(self.inlet.h)

    def _sync_fixed_wall(self, t: float):
        wall_nodes = np.full(self.n_nodes, self.wall_temp)
        zero_resistance = np.zeros(self.n_nodes)
        self.wall_boundary.update_internal_state(
            T_node=wall_nodes,
            R_int=zero_resistance,
            current_time=t,
        )

    def _set_state(self, h_vec: np.ndarray, t: float):
        for vol, h in zip(self.channel.volumes, h_vec):
            vol.P = self.pressure
            vol.h = float(h)
            vol.update_properties(self.mat_na)

        self.inlet.P = self.pressure
        self.inlet.T = self.inlet_temp
        self.inlet.h = self.h_inlet
        self.inlet.update_properties(self.mat_na)

        for junc in self.channel.internal_junctions:
            junc.W = self.mass_flow
            junc.update_velocity()
        self.junc_in.W = self.mass_flow
        self.junc_in.target_W = self.mass_flow
        self.junc_in.update_velocity()
        self.junc_out.W = self.mass_flow
        self.junc_out.update_velocity()

        self._sync_fixed_wall(t)

    def _apply_coupling(self, t: float):
        self.channel.clear_sources()
        self._sync_fixed_wall(t)
        self.coupler.execute()

    def rhs(self, t, h_vec):
        self._set_state(h_vec, t)
        self._apply_coupling(t)
        return np.array(
            [vol.get_volume_derivatives(self.mat_na)[1] for vol in self.channel.volumes],
            dtype=float,
        )

    def audit_state(self, t: float, h_vec: np.ndarray):
        self._set_state(h_vec, t)
        self._apply_coupling(t)

        dhdt = np.array(
            [vol.get_volume_derivatives(self.mat_na)[1] for vol in self.channel.volumes],
            dtype=float,
        )
        q_wall_nodes = np.array(
            [vol.Q_wall + vol.Q_vol - vol.implicit_coeff * vol.T for vol in self.channel.volumes],
            dtype=float,
        )
        storage_nodes = np.array(
            [vol.rho * vol.vol * d_h for vol, d_h in zip(self.channel.volumes, dhdt)],
            dtype=float,
        )

        q_wall = float(np.sum(q_wall_nodes))
        q_storage = float(np.sum(storage_nodes))
        h_out = float(self.channel.volumes[-1].h)
        q_enthalpy = self.mass_flow * (h_out - self.h_inlet)
        q_balance_residual = q_wall - q_enthalpy - q_storage

        return {
            "Q_wall": q_wall,
            "Q_enthalpy": float(q_enthalpy),
            "Q_storage": q_storage,
            "Q_balance_residual": float(q_balance_residual),
            "h_in": self.h_inlet,
            "h_out": h_out,
            "T_out": float(self.channel.volumes[-1].T),
            "T_mean": float(np.mean([vol.T for vol in self.channel.volumes])),
        }


def run_case(sync_boundary_before_couple: bool):
    label = "scheduler-synced" if sync_boundary_before_couple else "manual-unsynced"
    print("\n" + "=" * 72)
    print(f"Running fluid-solid audit case: {label}")
    print("=" * 72)

    system = SimpleLumpedSystem(sync_boundary_before_couple=sync_boundary_before_couple)
    y0 = [system.fluid.volumes[0].h, system.solid.T[0]]
    t_span = (0.0, 50.0)
    t_eval = np.linspace(t_span[0], t_span[1], 1001)

    sol = solve_ivp(
        system.rhs,
        t_span,
        y0,
        t_eval=t_eval,
        method='RK45',
        max_step=0.1,
        rtol=1.0e-7,
        atol=1.0e-8,
    )
    if not sol.success:
        print(f"Integration failed for {label}: {sol.message}")
        return None

    h_res = sol.y[0]
    T_solid_res = sol.y[1]
    t_res = sol.t

    T_fluid_res = []
    audits = []
    for t, h, T_solid in zip(t_res, h_res, T_solid_res):
        T_fluid_res.append(_as_float(system.mat_na.temperature_from_enthalpy(h, 1.0e5)))
        audits.append(system.evaluate_interface_audit(t, h, T_solid))

    audit_arrays = {
        key: np.array([item[key] for item in audits], dtype=float)
        for key in audits[0]
    }
    max_q = max(np.max(np.abs(audit_arrays["Q_fluid"])), 1.0)
    max_idx = int(np.argmax(np.abs(audit_arrays["interface_residual"])))
    max_interface_residual = np.max(np.abs(audit_arrays["interface_residual"]))
    max_storage_residual = np.max(np.abs(audit_arrays["storage_residual"]))

    print(f"Final T_fluid       = {T_fluid_res[-1]:.6f} K")
    print(f"Final T_solid       = {T_solid_res[-1]:.6f} K")
    print(f"Final Q_fluid       = {audit_arrays['Q_fluid'][-1]: .6e} W")
    print(f"Final Q_solid_bound = {audit_arrays['Q_solid_boundary'][-1]: .6e} W")
    print(f"Final Q_solid_store = {audit_arrays['Q_solid_storage'][-1]: .6e} W")
    print(f"Final interface res = {audit_arrays['interface_residual'][-1]: .6e} W")
    print(f"Max interface res   = {max_interface_residual: .6e} W")
    print(f"Max residual time   = {t_res[max_idx]: .6e} s")
    print(f"Max relative res    = {max_interface_residual / max_q: .6e}")
    print(f"Max storage res     = {max_storage_residual: .6e} W")

    return {
        "label": label,
        "t": t_res,
        "T_fluid": np.array(T_fluid_res, dtype=float),
        "T_solid": T_solid_res,
        "audit": audit_arrays,
    }


def run_flow_through_enthalpy_audit():
    print("\n" + "=" * 72)
    print("Running flow-through enthalpy audit")
    print("=" * 72)

    system = FlowThroughEnthalpyAudit()
    y0 = np.array([vol.h for vol in system.channel.volumes], dtype=float)
    t_span = (0.0, 80.0)
    t_eval = np.linspace(t_span[0], t_span[1], 801)

    sol = solve_ivp(
        system.rhs,
        t_span,
        y0,
        t_eval=t_eval,
        method="BDF",
        max_step=0.2,
        rtol=1.0e-8,
        atol=1.0e-5,
    )
    if not sol.success:
        print(f"Flow-through audit failed: {sol.message}")
        return None

    audits = [system.audit_state(t, sol.y[:, i]) for i, t in enumerate(sol.t)]
    audit_arrays = {
        key: np.array([item[key] for item in audits], dtype=float)
        for key in audits[0]
    }

    final = audits[-1]
    max_q = max(np.max(np.abs(audit_arrays["Q_wall"])), 1.0)
    max_balance_res = np.max(np.abs(audit_arrays["Q_balance_residual"]))
    steady_gap = final["Q_wall"] - final["Q_enthalpy"]

    print(f"Final T_out             = {final['T_out']:.6f} K")
    print(f"Final Q_wall_total      = {final['Q_wall']: .6e} W")
    print(f"Final W*(h_out-h_in)    = {final['Q_enthalpy']: .6e} W")
    print(f"Final storage term      = {final['Q_storage']: .6e} W")
    print(f"Final steady gap        = {steady_gap: .6e} W")
    print(f"Final full residual     = {final['Q_balance_residual']: .6e} W")
    print(f"Max full residual       = {max_balance_res: .6e} W")
    print(f"Max relative residual   = {max_balance_res / max_q: .6e}")

    plt.figure(figsize=(8, 6))
    plt.plot(sol.t, audit_arrays["Q_wall"], label="Q_wall_total")
    plt.plot(sol.t, audit_arrays["Q_enthalpy"], "--", label="W*(h_out-h_in)")
    plt.plot(sol.t, audit_arrays["Q_storage"], ":", label="storage")
    plt.plot(sol.t, audit_arrays["Q_enthalpy"] + audit_arrays["Q_storage"],
             "-.", label="enthalpy + storage")
    plt.xlabel("Time (s)")
    plt.ylabel("Power (W)")
    plt.title("Flow-Through Enthalpy Balance Audit")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("simple_flow_enthalpy_audit.png", dpi=160)
    print("Flow audit saved to 'simple_flow_enthalpy_audit.png'")

    return {
        "t": sol.t,
        "audit": audit_arrays,
        "solution": sol,
    }


def run_simulation_v2():
    original = run_case(sync_boundary_before_couple=False)
    synced = run_case(sync_boundary_before_couple=True)
    flow_audit = run_flow_through_enthalpy_audit()
    results = [item for item in (original, synced) if item is not None]
    if not results:
        return

    plt.figure(figsize=(8, 6))
    for item in results:
        plt.plot(item["t"], item["T_solid"], label=f"{item['label']} solid")
        plt.plot(item["t"], item["T_fluid"], linestyle="--", label=f"{item['label']} fluid")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    plt.title("Fluid-Solid Coupling Temperature Audit")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("simple_transient_result.png", dpi=160)
    print("Result saved to 'simple_transient_result.png'")

    plt.figure(figsize=(8, 6))
    for item in results:
        plt.semilogy(
            item["t"],
            np.abs(item["audit"]["interface_residual"]) + 1.0e-30,
            label=f"{item['label']}: |Q_fluid + Q_solid_boundary|",
        )
        plt.semilogy(
            item["t"],
            np.abs(item["audit"]["storage_residual"]) + 1.0e-30,
            linestyle="--",
            label=f"{item['label']}: |Q_fluid + C*dT/dt|",
        )
    plt.xlabel("Time (s)")
    plt.ylabel("Energy residual (W)")
    plt.title("Fluid-Solid Interface Energy Residual")
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.savefig("simple_fluid_solid_energy_audit.png", dpi=160)
    print("Audit saved to 'simple_fluid_solid_energy_audit.png'")


if __name__ == "__main__":
    run_simulation_v2()
