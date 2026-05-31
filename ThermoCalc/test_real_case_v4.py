import numpy as np
import te_solver
import time

print(f"正在载入模块: {te_solver.__file__}")

# ==============================================================================
# 1. 原始数据准备 (与 test.cpp 保持一致)
# ==============================================================================
# TE_data = [
#     1491.29, 1500.11, 1516, 1537.34, 1562.78, 1591.18, 1621.63, 1653.35, 1685.71, 1718.16,
#     1750.24, 1781.57, 1811.82, 1840.71, 1868, 1893.49, 1917, 1938.4, 1957.56, 1974.4,
#     1988.83, 2000.8, 2010.27, 2017.21, 2021.62, 2023.48, 2022.83, 2019.69, 2014.09, 2006.08,
#     1995.73, 1983.11, 1968.29, 1951.39, 1932.5, 1911.77, 1889.32, 1865.34, 1840.01, 1813.57,
#     1786.27, 1758.44, 1730.45, 1702.77, 1675.95, 1650.68, 1627.8, 1608.39, 1593.75, 1585.55
# ]
#
# TC_data = [
#     760.774, 762.065, 763.823, 765.81, 767.951, 770.213, 772.576, 775.027, 777.551, 780.139,
#     782.779, 785.462, 788.178, 790.92, 793.677, 796.445, 799.211, 801.974, 804.72, 807.451,
#     810.151, 812.824, 815.452, 818.044, 820.578, 823.068, 825.488, 827.858, 830.144, 832.379,
#     834.517, 836.603, 838.58, 840.509, 842.317, 844.084, 845.715, 847.323, 848.779, 850.236,
#     851.526, 852.851, 853.995, 855.227, 856.262, 857.463, 858.453, 859.715, 860.732, 862.062
# ]

TE_data = [1800.0] * 50
TC_data = [800.0] * 50

TE_data_initial = [600.] * len(TE_data)
TC_data_initial = [600.] * len(TE_data)

N_elem = 6
n_node = len(TE_data)

# ==============================================================================
# 2. 填充 InputData
# ==============================================================================
input_data = te_solver.InputData()
input_data.N_elements = N_elem
input_data.n_axi = n_node

# 物理场分布 (复制 6 份)
input_data.Temitter = np.tile(TE_data, (N_elem, 1))
input_data.Tcollector = np.tile(TC_data, (N_elem, 1))
input_data.Temitter = np.tile(TE_data_initial, (N_elem, 1))
input_data.Tcollector = np.tile(TC_data_initial, (N_elem, 1))
input_data.Tcs = np.full((N_elem, n_node), 610.0)
input_data.V_init = np.full((N_elem, n_node), 0.2)

# 几何分布
dl_val = 0.377 / float(n_node)
input_data.dlE = np.full((N_elem, n_node), dl_val)
input_data.dlC = np.full((N_elem, n_node), dl_val)

# 标量参数
input_data.crossAreaE = np.full(N_elem, 6.667e-5)
input_data.crossAreaC = np.full(N_elem, 1.0786e-4)
sideE_val = 0.00092855424159680002 * 25.0 / float(n_node)
input_data.sideAreaE = np.full((N_elem, n_node), sideE_val)
sideC_val = 0.00097592945800480005 * 25.0 / float(n_node)
input_data.sideAreaC = np.full((N_elem, n_node), sideC_val)

input_data.U_init = np.full(N_elem, 1.6)
input_data.d_gap = np.full(N_elem, 0.5)
input_data.Itarget = np.full(N_elem, 200.0)

# 导线电阻与电压
input_data.resistanceWire = np.zeros((N_elem, 4))
wireU_single = np.array([0.8, 0.8, 0.0, 0.0])
input_data.wireU = np.tile(wireU_single, (N_elem, 1))

# 全局控制
input_data.mode = te_solver.CalculationMode.FixedResistance
input_data.target_val = 0.0044
input_data.I_total_init = 150.0

# ==============================================================================
# 3. 创建与计算
# ==============================================================================
print("正在调用 C++ 创建电路...")
circuit = te_solver.create_circuit(input_data)

# 设置计算条件
circuit.Utarget = 0.77 * 4.0 / 6.0
circuit.Uout = 0.77 * 4.0 / 6.0
circuit.Iout = 150.0
circuit.Rload = 0.0044
circuit.isFixedR = True
circuit.isFixedU = False

for i in range(N_elem):
    # 直接赋值给属性！Pybind11 会自动搬运数据到 C++ vector
    # 使用 [:] 原地修改或者直接属性赋值均可，属性赋值更安全
    circuit.TECs[i].Temitter = np.array(TE_data)
    circuit.TECs[i].Tcollector = np.array(TC_data)

print("正在进行计算 (Fixed Resistance Mode)...")
t0 = time.time()
circuit.calc()
t1 = time.time()
print(f"计算耗时: {(t1 - t0) * 1000:.2f} ms")

t0 = time.time()
for i in range(50):
    circuit.calc()
t1 = time.time()
print(f"50次计算耗时: {(t1 - t0) * 1000:.2f} ms")

t0 = time.time()
for i in range(100):
    circuit.calc()
t1 = time.time()
print(f"100次计算耗时: {(t1 - t0) * 1000:.2f} ms")

# ==============================================================================
# 4. 详细结果验证 (新增部分)
# ==============================================================================
print("\n" + "=" * 60)
print("详细物理场结果验证报告")
print("=" * 60)

# --- A. 验证全局电路结果 ---
print(f"[全局电路结果 CircuitTECs]")
print(f"  Iout (总电流)    : {circuit.Iout:.6f} A")
print(f"  Uout (总电压)    : {circuit.Uout:.6f} V")
print(f"  Rload (负载)     : {circuit.Rload:.6f} Ohm")
print("-" * 60)

# --- B. 验证单个元件结果 ---
# 我们选取第 0 个和第 2 个元件作为代表进行检查
check_indices = [0]


def print_vector_info(name, vec_data):
    """辅助函数：打印向量的统计信息"""
    # Pybind11 返回的 vector 可以直接作为 list 使用，也可以转 numpy
    v = np.array(vec_data)
    print(f"    {name:<12} | Min: {v.min():.4e} | Max: {v.max():.4e} | Mean: {v.mean():.4e}")
    # 打印 首端-中点-末端 的值，方便快速核对分布趋势
    mid = len(v) // 2
    print(f"                 | [0]: {v[0]:.4e}   | [{mid}]: {v[mid]:.4e}   | [-1]: {v[-1]:.4e}")


for idx in check_indices:
    tec = circuit.TECs[idx]
    print(f"\n[元件 TECs[{idx}] 详细数据]")

    # 1. 标量结果
    print(f"  1. 标量状态:")
    print(f"    Current I    : {tec.I:.6f} A")
    print(f"    Voltage U    : {tec.U:.6f} V")

    # 2. 关键物理场分布
    print(f"  2. 物理场分布 (Vectors, Length={len(tec.J)}):")

    # J: 电流密度分布
    print_vector_info("J (A/cm2)", tec.J)

    # V: 极板间电压分布
    print_vector_info("V (Volts)", tec.V)

    # UE / UC: 电极电势分布
    print_vector_info("UE (Volts)", tec.UE)
    print_vector_info("UC (Volts)", tec.UC)

    # rhoE / rhoC: 电阻率分布
    print_vector_info("rhoE (Ohm.m)", tec.rhoE)
    print_vector_info("rhoC (Ohm.m)", tec.rhoC)

    # IEsec / ICsec: 截面电流分布 (您新增的字段)
    print_vector_info("IEsecSingle", tec.IEsecSingle)
    print_vector_info("ICsecSingle", tec.ICsecSingle)

    print("-" * 30)

# --- C. 单独测试 calc_current 方法 (验证 debug 接口) ---
print("\n[Debug 接口测试]")
print("尝试手动调用 TECs[0].calc_current()...")
# 注意：直接调用可能会改变内部状态，这里仅演示调用是否成功
delta_I = circuit.TECs[0].calc_current()
print(f"calc_current() 返回值 (delta I): {delta_I:.6e}")
print("=" * 60)

t0 = time.time()
TE_data = [x + 50. for x in TE_data]
real_TE_array = np.array(TE_data)
for i in range(N_elem):
    # 直接赋值给属性！Pybind11 会自动搬运数据到 C++ vector
    # 使用 [:] 原地修改或者直接属性赋值均可，属性赋值更安全
    circuit.TECs[i].Temitter = real_TE_array
circuit.calc()
t1 = time.time()
print(f"计算耗时: {(t1 - t0) * 1000:.2f} ms")

# ==============================================================================
# 4. 详细结果验证 (新增部分)
# ==============================================================================
print("\n" + "=" * 60)
print("详细物理场结果验证报告")
print("=" * 60)

# --- A. 验证全局电路结果 ---
print(f"[全局电路结果 CircuitTECs]")
print(f"  Iout (总电流)    : {circuit.Iout:.6f} A")
print(f"  Uout (总电压)    : {circuit.Uout:.6f} V")
print(f"  Rload (负载)     : {circuit.Rload:.6f} Ohm")
print("-" * 60)

# --- B. 验证单个元件结果 ---
# 我们选取第 0 个和第 2 个元件作为代表进行检查
check_indices = [0]


def print_vector_info(name, vec_data):
    """辅助函数：打印向量的统计信息"""
    # Pybind11 返回的 vector 可以直接作为 list 使用，也可以转 numpy
    v = np.array(vec_data)
    print(f"    {name:<12} | Min: {v.min():.4e} | Max: {v.max():.4e} | Mean: {v.mean():.4e}")
    # 打印 首端-中点-末端 的值，方便快速核对分布趋势
    mid = len(v) // 2
    print(f"                 | [0]: {v[0]:.4e}   | [{mid}]: {v[mid]:.4e}   | [-1]: {v[-1]:.4e}")


for idx in check_indices:
    tec = circuit.TECs[idx]
    print(f"\n[元件 TECs[{idx}] 详细数据]")

    # 1. 标量结果
    print(f"  1. 标量状态:")
    print(f"    Current I    : {tec.I:.6f} A")
    print(f"    Voltage U    : {tec.U:.6f} V")

    # 2. 关键物理场分布
    print(f"  2. 物理场分布 (Vectors, Length={len(tec.J)}):")

    # J: 电流密度分布
    print_vector_info("J (A/cm2)", tec.J)

    # V: 极板间电压分布
    print_vector_info("V (Volts)", tec.V)

    # UE / UC: 电极电势分布
    print_vector_info("UE (Volts)", tec.UE)
    print_vector_info("UC (Volts)", tec.UC)

    # rhoE / rhoC: 电阻率分布
    print_vector_info("rhoE (Ohm.m)", tec.rhoE)
    print_vector_info("rhoC (Ohm.m)", tec.rhoC)

    # IEsec / ICsec: 截面电流分布 (您新增的字段)
    print_vector_info("IEsecSingle", tec.IEsecSingle)
    print_vector_info("ICsecSingle", tec.ICsecSingle)

    # phiE / phiC: 功函数分布 (您新增的字段)
    print_vector_info("phiE", tec.phiE)
    print_vector_info("phiC", tec.phiC)

    # Vd: 电弧降分布 (您新增的字段)
    print_vector_info("Vd", tec.Vd)

    print("-" * 30)

# --- C. 单独测试 calc_current 方法 (验证 debug 接口) ---
print("\n[Debug 接口测试]")
print("尝试手动调用 TECs[0].calc_current()...")
# 注意：直接调用可能会改变内部状态，这里仅演示调用是否成功
delta_I = circuit.TECs[0].calc_current()
print(f"calc_current() 返回值 (delta I): {delta_I:.6e}")
print("=" * 60)
