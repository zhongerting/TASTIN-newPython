import numpy as np
from Materials.Base import SolidMaterial, Numeric

from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Couplers import GapCouple2D, ActiveGapCouple2D, TECCouple2D


class VirtualMaterial(SolidMaterial):
    """
    [验证专用] 虚拟线性材料

    特性:
    1. 物性 (k, rho, cp) 为常数，不随温度变化。
    2. 用于消除非线性干扰，验证求解器和耦合器的数值逻辑。
    """

    def __init__(self, k: float = 10.0, rho: float = 1.0, cp: float = 1.0, name: str = "Virtual_Mat"):
        """
        :param k: 导热系数 [W/mK] (默认 10, 方便口算)
        :param rho: 密度 [kg/m3] (默认 1)
        :param cp: 比热容 [J/kgK] (默认 1)
        """
        super().__init__(name=name, formula="N/A")
        self.k_val = float(k)
        self.rho_val = float(rho)
        self.cp_val = float(cp)

    def conductivity(self, T: Numeric) -> Numeric:
        # 无论 T 是标量还是数组，都返回常数
        if np.ndim(T) > 0:
            return np.full_like(T, self.k_val, dtype=float)
        return self.k_val

    def density(self, T: Numeric) -> Numeric:
        if np.ndim(T) > 0:
            return np.full_like(T, self.rho_val, dtype=float)
        return self.rho_val

    def heat_capacity(self, T: Numeric) -> Numeric:
        if np.ndim(T) > 0:
            return np.full_like(T, self.cp_val, dtype=float)
        return self.cp_val

    # 显式定义热扩散系数 (alpha = k / (rho*cp))，方便验证脚本调用
    # @property
    # def diffusivity(self) -> float:
    #     return self.k_val / (self.rho_val * self.cp_val)

def setup_test_bench():
    print("=" * 60)
    print("STEP 1 & 2: 建立虚拟物性与单网格测试台")
    print("=" * 60)

    # 1. 实例化虚拟材料 (k=10, rho=1, cp=1)
    #    选择这些值是为了让计算结果一目了然 (如 Flux = k * dT / dx = 10 * dT / dx)
    mat = VirtualMaterial(k=10.0, rho=1.0, cp=1.0, name="Test_Linear_Mat")
    print(f"[Material] Created {mat.name}: k={mat.k_val}, rho={mat.rho_val}, cp={mat.cp_val}")

    # 2. 创建网格 (1cm x 1cm, 单网格)
    #    使用您 Mesh2D 的构造函数签名: (x_dim, n_x, y_dim, n_y, ...)
    Lx = 0.01  # [m]
    Ly = 0.01  # [m]

    # Emitter 网格
    mesh1 = Mesh2D(x_dim=Lx, n_x=1, y_dim=Ly, n_y=1, geometry_type='cartesian')

    # Collector 网格
    mesh2 = Mesh2D(x_dim=Lx, n_x=1, y_dim=Ly, n_y=1, geometry_type='cartesian')

    print(f"[Mesh] Created two 1x1 meshes (Lx={Lx}m, Ly={Ly}m).")

    # 3. 实例化导热对象 (Emitter 和 Collector)
    #    初始温度设为 300 K
    T_init = 300.0
    emitter = HeatConduction2D(mesh=mesh1, material=mat, initial_temp=T_init)
    collector = HeatConduction2D(mesh=mesh2, material=mat, initial_temp=T_init)

    # 标记名称以便调试 (Python动态属性)
    emitter.name = "Emitter_Obj"
    collector.name = "Collector_Obj"

    print(f"[Solver] Instantiated 'Emitter'   @ T={emitter.T[0]:.2f} K")
    print(f"[Solver] Instantiated 'Collector' @ T={collector.T[0]:.2f} K")

    # 4. 验证初始化状态 (Mass & Capacitance Check)
    #    体积 V = Lx * Ly (单位深度假设? Mesh2D代码中 cartesian volume = dx * dy)
    #    Check Mesh.py: volumes_2d = dx_2d * dy_2d -> 0.01 * 0.01 = 1e-4 m^3
    expected_vol = Lx * Ly
    actual_vol = emitter.mesh.geom_data.volumes[0]

    #    热容 C = rho * cp * V = 1 * 1 * 1e-4 = 1e-4 J/K
    expected_C = mat.rho_val * mat.cp_val * expected_vol
    actual_C = emitter.thermal_capacitance[0]

    print("-" * 30)
    print(f"验证几何与热容:")
    print(f"  Volume Expected: {expected_vol:.2e} m^3 | Actual: {actual_vol:.2e} m^3")
    print(f"  ThermC Expected: {expected_C:.2e} J/K   | Actual: {actual_C:.2e} J/K")

    if abs(actual_C - expected_C) < 1e-12:
        print("  >> Status: PASS ✅")
    else:
        print("  >> Status: FAIL ❌ (请检查 Mesh2D 体积计算逻辑)")

    print("=" * 60)
    return emitter, collector

def calculate_internal_resistance(obj):
    """
    [辅助] 计算单网格物体的内部导热热阻 (Node -> Surface)
    R''_int = (Lx/2) / k_solid
    """
    Lx = obj.mesh.x_dim
    k = obj.material.conductivity(300.0) # Virtual mat is const
    # 距离是长度的一半
    dist = Lx / 2.0
    return dist / k

def create_solver_pair(T1, T2):
    """
    [辅助函数] 创建一对标准的 1x1 测试求解器
    :param T1: Emitter (Obj1) 温度 [K]
    :param T2: Collector (Obj2) 温度 [K]
    """
    # 1. 虚拟材料 (k=10)
    mat = VirtualMaterial(k=10.0, rho=1.0, cp=1.0)

    # 2. 单网格 1cm x 1cm
    Lx, Ly = 0.01, 0.01
    mesh1 = Mesh2D(x_dim=Lx, n_x=1, y_dim=Ly, n_y=1, geometry_type='cartesian')
    mesh2 = Mesh2D(x_dim=Lx, n_x=1, y_dim=Ly, n_y=1, geometry_type='cartesian')

    # 3. 求解器
    obj1 = HeatConduction2D(mesh=mesh1, material=mat, initial_temp=T1)
    obj2 = HeatConduction2D(mesh=mesh2, material=mat, initial_temp=T2)

    obj1.name = "Emitter"
    obj2.name = "Collector"

    return obj1, obj2

def verify_radiation_only():
    print("\n" + "=" * 60)
    print("【测试 1】GapCouple2D - 仅辐射验证 (Radiation Only)")
    print("=" * 60)

    T_hot = 400.0
    T_cold = 300.0
    obj1, obj2 = create_solver_pair(T_hot, T_cold)

    # 实例化耦合器
    coupler = GapCouple2D(
        obj1=obj1, obj2=obj2, direction='right',
        gap_width=1e-3, gas_conductivity=0.0,
        emissivity1=1.0, emissivity2=1.0
    )

    coupler.sync()

    # 获取求解器热流
    q_in_solver = obj1.boundaries['right'].compute_net_flux_for_solver()
    q_in_solver = obj2.boundaries['left'].compute_net_flux_for_solver()
    flux_solver = -q_in_solver[0] / obj1.boundaries['right'].area[0]

    # --- 修正后的理论验证方法 ---
    # 为了验证 Stefan-Boltzmann 定律，我们应该使用求解器实际算出的“表面温度”
    # 这样可以排除内部导热压降的影响，纯粹验证辐射公式。

    # 获取 Solver 实际算出的表面温度
    T_s1 = obj1.boundaries['right'].T_surface[0]
    T_s2 = obj2.boundaries['left'].T_surface[0]

    sigma = 5.670374419e-8
    flux_theory_pure = sigma * (T_s1 ** 4 - T_s2 ** 4)

    print(f"  [诊断数据]")
    print(f"    T_node_hot  = {T_hot} K")
    print(f"    T_surf_hot  = {T_s1:.4f} K (因内部热阻降低)")
    print(f"    T_surf_cold = {T_s2:.4f} K (因内部热阻升高)")
    print("-" * 40)
    print(f"  Flux (Solver)         : {flux_solver:.4f} W/m^2")
    print(f"  Flux (Theory @ Surf)  : {flux_theory_pure:.4f} W/m^2")
    print(f"  (验证公式: Flux = sigma * (T_surf1^4 - T_surf2^4))")

    # 允许 0.1% 的误差 (由于 h_rad 的线性化近似)
    if abs(flux_solver - flux_theory_pure) / flux_solver < 1e-3:
        print("  >> Status: PASS ✅")
    else:
        print("  >> Status: FAIL ❌")

def verify_conduction_only():
    print("\n" + "=" * 60)
    print("【测试 2】GapCouple2D - 仅气体导热验证 (Gas Conduction Only)")
    print("=" * 60)

    T_hot = 400.0
    T_cold = 300.0
    k_solid = 10.0
    Lx = 0.01

    obj1, obj2 = create_solver_pair(T_hot, T_cold)

    gap_val = 0.001
    k_gas_val = 1.0

    coupler = GapCouple2D(
        obj1=obj1, obj2=obj2, direction='right',
        gap_width=gap_val, gas_conductivity=k_gas_val,
        emissivity1=0.0, emissivity2=0.0
    )

    coupler.sync()

    q_in_solver = obj1.boundaries['right'].compute_net_flux_for_solver()
    flux_solver = -q_in_solver[0] / obj1.boundaries['right'].area[0]

    # --- 修正后的理论真值计算 (考虑串联热阻) ---
    # 1. 内部热阻 R_int = (L/2) / k_solid
    R_int_val = (Lx / 2.0) / k_solid

    # 2. 间隙热阻 R_gap = gap / k_gas
    R_gap_val = gap_val / k_gas_val

    # 3. 总热阻 (串联)
    R_total = R_int_val + R_gap_val + R_int_val

    flux_theory = (T_hot - T_cold) / R_total

    print(f"  [热阻分析 Unit Area]")
    print(f"    R_int (Solid) : {R_int_val:.5f} K·m^2/W (x2)")
    print(f"    R_gap (Gas)   : {R_gap_val:.5f} K·m^2/W")
    print(f"    R_total       : {R_total:.5f} K·m^2/W")
    print("-" * 40)
    print(f"  Flux (Solver) : {flux_solver:.4f} W/m^2")
    print(f"  Flux (Theory) : {flux_theory:.4f} W/m^2")

    if abs(flux_solver - flux_theory) < 1e-5:
        print("  >> Status: PASS ✅")
    else:
        print("  >> Status: FAIL ❌")

def verify_active_gap_couple_logic():
    print("\n" + "=" * 60)
    print("【最终验证】ActiveGapCouple2D (现有类) - 主动热流源分配验证")
    print("=" * 60)

    # --- 参数设置 ---
    T_base = 300.0  # 无温差环境
    k_solid = 10.0  # 固体导热系数
    Lx = 0.01  # 固体厚度

    gap_val = 0.001  # 间隙 1mm
    k_gas_val = 0.1  # 气体导热系数 (较小，使得 R_gap 较大)

    Q_source_set = 100.0  # 设定的主动热源 [W]

    # --- 1. 创建求解器 ---
    obj1, obj2 = create_solver_pair(T_base, T_base)

    # --- 2. 实例化项目中的 ActiveGapCouple2D ---
    # 给定极小的发射率以忽略辐射，专注于导热+源项
    coupler = ActiveGapCouple2D(
        obj1=obj1, obj2=obj2, direction='right',
        gap_width=gap_val, gas_conductivity=k_gas_val,
        emissivity1=1e-10, emissivity2=1e-10
    )

    # --- 3. 施加源项 ---
    coupler.set_active_heat_source(Q_source_set)
    print(f"  Conditions:")
    print(f"    T_init   = {T_base} K")
    print(f"    Q_source = {Q_source_set} W")
    print("-" * 40)

    # --- 4. 执行同步 ---
    coupler.sync()

    # --- 5. 获取 Solver 结果 ---
    # Obj1 (Emitter): 应该流出热量
    q_in_1 = obj1.boundaries['right'].compute_net_flux_for_solver()
    flux_1 = -q_in_1[0] / obj1.boundaries['right'].area[0]  # W/m2
    power_1 = flux_1 * obj1.boundaries['right'].area[0]  # W

    # --- 6. 理论真值计算 ---
    # 物理模型: 电流源 Q 并联 R_gap => 等效电压源 V = Q * R_gap 串联
    # 分流公式: P_out = Q * R_gap / (R_int1 + R_gap + R_int2)

    # R_int (Unit Area): (L/2)/k
    R_int_unit = (Lx / 2.0) / k_solid
    # R_gap (Unit Area): gap/k_gas
    R_gap_unit = gap_val / k_gas_val

    R_total_unit = R_int_unit + R_gap_unit + R_int_unit

    # 理论预期功率
    ratio = R_gap_unit / R_total_unit
    power_expected = Q_source_set * ratio

    print(f"  [热阻分析 Unit Area]")
    print(f"    R_int (Solid) x2 : {2 * R_int_unit:.5f} K·m^2/W")
    print(f"    R_gap (Gas)      : {R_gap_unit:.5f} K·m^2/W")
    print(f"    Total Circuit R  : {R_total_unit:.5f} K·m^2/W")
    print(f"    Scaling Ratio    : {ratio:.4%}")
    print("-" * 40)
    print(f"  [结果对比]")
    print(f"    Set Source Q     : {Q_source_set:.4f} W")
    print(f"    Expected Output  : {power_expected:.4f} W (由于固体内部热阻分压)")
    print(f"    Solver Output    : {power_1:.4f} W")

    # --- 7. 判定 ---
    if abs(power_1 - power_expected) < 1e-3:
        print(f"  >> Status: PASS ✅")
    else:
        print(f"  >> Status: FAIL ❌")
        print(f"     Diff: {abs(power_1 - power_expected):.4f} W")
        print("     (如果失败，请检查 ActiveGapCouple2D.sync 中是否使用了 delta_T = Q * R_gap)")

def verify_tec_couple_logic():
    print("\n" + "=" * 60)
    print("【最终验证】TECCouple2D (现有类) - 热电转换与能量守恒验证")
    print("=" * 60)

    # --- 基础参数 ---
    T_base = 1000.0  # 高温环境
    k_solid = 10.0
    Lx = 0.01

    gap_val = 0.001
    k_gas_val = 1.0  # 使得 R_gap 较小，接近真实 TEC 情况

    # --- 理论热阻计算 (Unit Area) ---
    R_int_unit = (Lx / 2.0) / k_solid  # 0.0005
    R_gap_unit = gap_val / k_gas_val  # 0.0010
    R_total_unit = 2 * R_int_unit + R_gap_unit  # 0.0020

    # Emitter 侧外部热阻 (R_gap + R_int_coll)
    R_ext_emit = R_gap_unit + R_int_unit  # 0.0015

    print(f"  [系统热阻分析]")
    print(f"    R_int (Solid) : {R_int_unit:.5f}")
    print(f"    R_gap (Gas)   : {R_gap_unit:.5f}")
    print(f"    R_total       : {R_total_unit:.5f}")
    print("-" * 40)

    # ==========================================================================
    # Case 1: 纯电子冷却 (Pure Cooling)
    # ==========================================================================
    print("\n[Case 1] 纯电子冷却测试 (Q_emit=-100W, Q_coll=0)")
    obj1, obj2 = create_solver_pair(T_base, T_base)

    # 实例化 (关闭辐射以验证纯导热耦合逻辑)
    coupler = TECCouple2D(
        obj1=obj1, obj2=obj2, direction='right',
        gap_width=gap_val, gas_conductivity=k_gas_val,
        emissivity1=1e-10, emissivity2=1e-10
    )

    Q1_set = -100.0  # 冷却
    Q2_set = 0.0
    coupler.set_tec_sources(Q_emitter=Q1_set, Q_collector=Q2_set)
    coupler.sync()

    # 获取结果
    q_in_1 = obj1.boundaries['right'].compute_net_flux_for_solver()
    flux_1 = -q_in_1[0] / obj1.boundaries['right'].area[0]
    p1_solver = flux_1 * obj1.boundaries['right'].area[0]

    # 理论预测: Emitter 流出功率 = -Q1 * (R_ext / R_total)
    # 说明: 源项 Q1 (流出) 在界面处分流，一部分由 Emitter 内部供给(流出增加)，一部分由 Gap 供给。
    p1_theory = -Q1_set * (R_ext_emit / R_total_unit)

    print(f"    Solver Outflow: {p1_solver:.4f} W")
    print(f"    Theory Outflow: {p1_theory:.4f} W")

    if abs(p1_solver - p1_theory) < 1e-4:
        print("    >> Status: PASS ✅")
    else:
        print("    >> Status: FAIL ❌")

    # ==========================================================================
    # Case 2: 发电模式 (Power Generation & Energy Balance)
    # ==========================================================================
    print("\n[Case 2] 发电模式测试 (Q_emit=-100W, Q_coll=+60W)")
    print("         预期提取电功率: 40W")

    obj1, obj2 = create_solver_pair(T_base, T_base)

    coupler = TECCouple2D(
        obj1=obj1, obj2=obj2, direction='right',
        gap_width=gap_val, gas_conductivity=k_gas_val,
        emissivity1=1e-10, emissivity2=1e-10
    )

    Q1_set = -100.0
    Q2_set = 60.0
    coupler.set_tec_sources(Q_emitter=Q1_set, Q_collector=Q2_set)
    coupler.sync()

    # Solver 1 Outflow (Emitter)
    q_in_1 = obj1.boundaries['right'].compute_net_flux_for_solver()
    p1_solver = -q_in_1[0]  # 注意: compute_net_flux 返回的是流入，取反为流出

    # Solver 2 Inflow (Collector)
    q_in_2 = obj2.boundaries['left'].compute_net_flux_for_solver()
    p2_solver = q_in_2[0]  # 流入

    # 理论预测 (叠加原理)
    # P1(Out) = Cooling_Term - Back_Heating_Term
    #         = (-Q1 * R_ext/R_tot) - (Q2 * R_int/R_tot)
    term_cool = -Q1_set * (R_ext_emit / R_total_unit)
    term_heat = - (Q2_set * R_int_unit / R_total_unit)
    p1_theory = term_cool + term_heat

    # P2(In) = Heating_Term + Cooling_Drag_Term
    #        = (Q2 * R_ext_coll/R_tot) + (Q1 * R_int_emit/R_tot)
    #        (注: Q1为负，表示Emitter冷却会减少Collector的流入)
    #        R_ext_coll (Collector看向外部) = R_gap + R_int_emit = R_ext_emit
    term_heat_coll = Q2_set * (R_ext_emit / R_total_unit)
    term_cool_drag = Q1_set * (R_int_unit / R_total_unit)
    p2_theory = term_heat_coll + term_cool_drag

    print(f"  [Emitter 流出]")
    print(f"    Solver: {p1_solver:.4f} W  | Theory: {p1_theory:.4f} W")

    print(f"  [Collector 流入]")
    print(f"    Solver: {p2_solver:.4f} W  | Theory: {p2_theory:.4f} W")

    # 能量守恒检查
    # 热系统损失 = 流出 - 流入
    # 应等于 = 电功率输出
    heat_loss = p1_solver - p2_solver
    elec_power = -(Q1_set + Q2_set)

    print("-" * 40)
    print(f"  热系统净损失 (P_out - P_in): {heat_loss:.4f} W")
    print(f"  应等于电功率 (-(Q1+Q2))    : {elec_power:.4f} W")

    if abs(heat_loss - elec_power) < 1e-4:
        print("  >> Energy Balance: PASS ✅ (System correctly models non-conservative power extraction)")
    else:
        print("  >> Energy Balance: FAIL ❌")


if __name__ == "__main__":
    # setup_test_bench()
    # verify_radiation_only()
    # verify_conduction_only()
    # verify_active_gap_couple_logic()
    verify_tec_couple_logic()