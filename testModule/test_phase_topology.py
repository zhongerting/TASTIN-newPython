import sys
import numpy as np

from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork

# --- 1. Mock Classes (模拟物理组件) ---
class MockVolume:
    def __init__(self, name, P=1.0e5, T=600.0, z=0.0, length=1.0, flow_area=0.01):
        self.name = name
        self.P = P
        self.T = T
        self.h = T * 1000.0 # 假定 Cp=1000
        self.z = z
        self.length = length
        self.flow_area = flow_area
        # 物性 (人为设定以便手算验证)
        self.rho = 800.0
        self.mu = 1.0e-4

class MockJunction:
    def __init__(self, name, from_vol, to_vol, W=0.0):
        self.name = name
        self.from_vol = from_vol
        self.to_vol = to_vol
        self.W = W
        # 几何
        self.length = 0.5 # 连接本身的长度 (通常不用，用 inertial length)
        self.flow_area = 0.01
        self.hydraulic_diameter = 0.1
        self.k_loss = 0.0 # 局部阻力


def test_topology_phase1():
    print("🧪 Testing Phase 1: Topology Building...")

    # 1. 创建物理对象
    inlet = MockVolume("Inlet", P=1.2e5, is_pressure_boundary=True)  # 稳压器1
    vol1 = MockVolume("Vol_1")
    vol2 = MockVolume("Vol_2")
    outlet = MockVolume("Outlet", P=1.0e5, is_pressure_boundary=True)  # 稳压器2

    j1 = MockJunction("J_In", inlet, vol1)
    j2 = MockJunction("J_1_2", vol1, vol2)
    j3 = MockJunction("J_Out", vol2, outlet)

    volumes = [inlet, vol1, vol2, outlet]
    junctions = [j1, j2, j3]

    # 2. 初始化网络
    net = HydraulicNetwork(volumes, junctions)

    # 3. 验证节点映射
    print("   -> Checking Nodes...")
    assert len(net.vol_to_idx) == 4
    assert net.vol_to_idx[inlet] == 0
    assert net.vol_to_idx[outlet] == 3

    # 4. 验证稳压器识别
    print("   -> Checking Pressurizers...")
    assert len(net.fixed_pressure_indices) == 2
    assert 0 in net.fixed_pressure_indices
    assert 3 in net.fixed_pressure_indices
    assert 1 not in net.fixed_pressure_indices

    # 5. 验证连接解析
    print("   -> Checking Connectivity...")
    assert len(net.junction_descriptors) == 3
    # 检查 J_In: Inlet(0) -> Vol_1(1)
    desc_j1 = net.junction_descriptors[0]
    assert desc_j1[0] == j1
    assert desc_j1[1] == 0
    assert desc_j1[2] == 1

    # 检查 J_Out: Vol_2(2) -> Outlet(3)
    desc_j3 = net.junction_descriptors[2]
    assert desc_j3[1] == 2
    assert desc_j3[2] == 3

    # 6. 打印拓扑供人工检查
    net.debug_topology()

    print("✅ Phase 1 Verification Passed!")

def test_friction_factor_logic():
    print("\n[Test 1] Friction Factor Calculation...")

    # A. 层流 (Re=100) -> f = 64/100 = 0.64
    f_lam = HydraulicNetwork._calc_friction_factor_static(100.0)
    assert abs(f_lam - 0.64) < 1e-6, f"Laminar error: {f_lam}"
    print(f"  ✅ Laminar (Re=100): f={f_lam:.4f}")

    # B. Blasius 湍流 (Re=10000) -> f = 0.3164 / 10000^0.25 = 0.3164 / 10 = 0.03164
    f_turb = HydraulicNetwork._calc_friction_factor_static(10000.0)
    assert abs(f_turb - 0.03164) < 1e-6, f"Blasius error: {f_turb}"
    print(f"  ✅ Blasius (Re=10000): f={f_turb:.5f}")

    # C. 高雷诺数迭代 (Re=1e6)
    # Karman-Prandtl: 1/sqrt(f) = 2*log10(Re*sqrt(f)) - 0.8
    # 这是一个隐式方程，只要函数能收敛并返回合理值即可
    f_high = HydraulicNetwork._calc_friction_factor_static(1.0e6)
    print(f"  ✅ High Re (Re=1e6): f={f_high:.6f} (Iterative Solved)")
    assert 0.01 < f_high < 0.02, "High Re friction factor out of range"


def test_gravity_term():
    print("\n[Test 2] Gravity Source Term (Vertical Pipe)...")

    # 场景: 垂直管，z轴向上为正。
    # Inlet (Bottom, z=0) -> Outlet (Top, z=10)
    # 这是一个上升流动。重力应该阻碍流动，或者在静止时产生静压差。

    vol_in = MockVolume("Bot", z=0.0)
    vol_out = MockVolume("Top", z=10.0)
    junc = MockJunction("J1", vol_in, vol_out, W=0.0)  # 静止

    net = HydraulicNetwork([vol_in, vol_out], [junc], gravity_vector=9.81)

    # 手动触发系数计算
    dt = 0.1
    net._calc_momentum_coeffs(dt)

    # 检查源项 b_j
    # 理论公式: b_j = a_j * (I*W_old + rho*g*(z_in - z_out))
    # W_old = 0 -> b_j = a_j * rho * g * (-10)
    # 预期 b_j 应该是负值 (重力作为阻力/负压头)

    idx = 0  # 第一个连接
    a_j = net.A_coeffs[idx]
    b_j = net.B_coeffs[idx]

    rho = 800.0
    expected_grav_head = rho * 9.81 * (0.0 - 10.0)  # -78480 Pa
    expected_b = a_j * expected_grav_head

    print(f"  Gravity Head: {expected_grav_head:.1f} Pa")
    print(f"  Calculated b_j: {b_j:.4e}, a_j: {a_j:.4e}")

    # 验证 b_j = a_j * (-78480)
    assert abs(b_j - a_j * expected_grav_head) < 1e-3
    assert b_j < 0, "Gravity term should be negative for upward pipe"
    print("  ✅ Gravity term sign and magnitude correct.")


def test_full_momentum_coeff():
    print("\n[Test 3] Full Momentum Coefficients (Horizontal Flow)...")

    # 场景: 水平流动，Re=10000 (湍流)
    rho = 800.0
    mu = 1.0e-3  # 调整粘度使 Re=10000
    W_flow = 1.0  # kg/s
    A = 0.01
    D = 0.1
    L = 2.0  # 两个节点各长 2.0m -> 惯性长度 = 2.0

    vol_in = MockVolume("In", z=0, length=2.0, flow_area=A)
    vol_in.rho = rho
    vol_in.mu = mu

    vol_out = MockVolume("Out", z=0, length=2.0, flow_area=A)
    vol_out.rho = rho
    vol_out.mu = mu

    junc = MockJunction("Link", vol_in, vol_out, W=W_flow)
    junc.flow_area = A
    junc.hydraulic_diameter = D
    junc.k_loss = 0.5  # 增加局部阻力

    net = HydraulicNetwork([vol_in, vol_out], [junc])

    # 手算验证
    dt = 0.1

    # 1. 速度与 Re
    vel = W_flow / (rho * A)  # 1 / 8 = 0.125 m/s
    Re = (rho * vel * D) / mu  # (800 * 0.125 * 0.1) / 1e-3 = 10000

    # 2. 阻力系数
    # Blasius: f = 0.3164 / 10000^0.25 = 0.03164
    f1 = HydraulicNetwork._calc_friction_factor_static(Re)
    f = 0.03164

    print("f1 = ", f1, " f = ", f)

    # 3. 线性化阻力 K_linear
    # term_geom = f*L/D + K_loss = 0.03164 * 2.0 / 0.1 + 0.5 = 0.6328 + 0.5 = 1.1328
    # K_linear = term_geom * |W| / (2 * rho * A^2)
    #          = 1.1328 * 1.0 / (2 * 800 * 1e-4) = 1.1328 / 0.16 = 7.08
    term_geom = (f * 2.0 / 0.1) + 0.5
    K_linear_expected = (term_geom * abs(W_flow)) / (2.0 * rho * A ** 2)

    # 4. 惯性项 I_term
    # I = L / (A * dt) = 2.0 / (0.01 * 0.1) = 2000.0
    I_term_expected = 2.0 / (A * dt)

    # 5. 导纳 a_j
    # a_j = 1 / (I + K)
    a_j_expected = 1.0 / (I_term_expected + K_linear_expected)

    # 6. 源项 b_j (水平，无加速，无重力)
    # b_j = a_j * (I * W_old)
    b_j_expected = a_j_expected * (I_term_expected * W_flow)

    # 执行代码计算
    net._calc_momentum_coeffs(dt)

    a_j_calc = net.A_coeffs[0]
    b_j_calc = net.B_coeffs[0]

    print(f"  Hand calc: a_j={a_j_expected:.4e}, b_j={b_j_expected:.4f}")
    print(f"  Code calc: a_j={a_j_calc:.4e}, b_j={b_j_calc:.4f}")

    assert abs(a_j_calc - a_j_expected) / a_j_expected < 1e-3
    assert abs(b_j_calc - b_j_expected) / b_j_expected < 1e-3
    print("  ✅ Momentum coefficients match hand calculation.")


def test_upwind_property():
    print("\n[Test 4] Upwind Property Scheme...")

    # 场景: 倒流 (W < 0), 应该取下游物性
    vol_in = MockVolume("In")
    vol_in.rho = 100.0  # 上游密度低

    vol_out = MockVolume("Out")
    vol_out.rho = 1000.0  # 下游密度高

    # 倒流
    junc = MockJunction("RevLink", vol_in, vol_out, W=-1.0)

    net = HydraulicNetwork([vol_in, vol_out], [junc])

    # 强制更新一下内部 rho_vec (虽然 Mock 没实现 update_properties，但 Network 初始化会读一次)
    net.rho_vec[0] = 100.0
    net.rho_vec[1] = 1000.0

    net._calc_momentum_coeffs(0.1)

    # 检查计算中使用的 Re (间接检查 rho)
    # 如果用了下游密度 1000，惯性阻力会不同。
    # 我们直接侵入检查 _get_upwind_properties 方法
    rho_used, _ = net._get_upwind_properties(0, 1, -1.0)

    print(f"  Flow W=-1.0, Upstream rho=100, Downstream rho=1000")
    print(f"  Property used by solver: rho={rho_used}")

    assert rho_used == 1000.0, "Should pick downstream property for reverse flow"
    print("  ✅ Upwind scheme works correctly.")

def test_phase2():
    test_upwind_property()
    test_gravity_term()
    test_full_momentum_coeff()
    test_friction_factor_logic()


if __name__ == "__main__":
    # test_topology_phase1()
    test_phase2()
