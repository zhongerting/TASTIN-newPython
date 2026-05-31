import numpy as np
import te_solver
import time

print(f"正在载入模块: {te_solver.__file__}")

# ==============================================================================
# 1. 原始数据准备 (来自 test.cpp - testCircuit)
# ==============================================================================

# 发射极温度 (50个点)
TE_data = [
    1491.29, 1500.11, 1516, 1537.34, 1562.78, 1591.18, 1621.63, 1653.35, 1685.71, 1718.16, 
    1750.24, 1781.57, 1811.82, 1840.71, 1868, 1893.49, 1917, 1938.4, 1957.56, 1974.4, 
    1988.83, 2000.8, 2010.27, 2017.21, 2021.62, 2023.48, 2022.83, 2019.69, 2014.09, 2006.08, 
    1995.73, 1983.11, 1968.29, 1951.39, 1932.5, 1911.77, 1889.32, 1865.34, 1840.01, 1813.57, 
    1786.27, 1758.44, 1730.45, 1702.77, 1675.95, 1650.68, 1627.8, 1608.39, 1593.75, 1585.55
]

# 接收极温度 (50个点)
TC_data = [
    760.774, 762.065, 763.823, 765.81, 767.951, 770.213, 772.576, 775.027, 777.551, 780.139, 
    782.779, 785.462, 788.178, 790.92, 793.677, 796.445, 799.211, 801.974, 804.72, 807.451, 
    810.151, 812.824, 815.452, 818.044, 820.578, 823.068, 825.488, 827.858, 830.144, 832.379, 
    834.517, 836.603, 838.58, 840.509, 842.317, 844.084, 845.715, 847.323, 848.779, 850.236, 
    851.526, 852.851, 853.995, 855.227, 856.262, 857.463, 858.453, 859.715, 860.732, 862.062
]

# 维度参数
N_elem = 6
n_node = len(TE_data) # 50

# ==============================================================================
# 2. 填充 InputData
# ==============================================================================
input_data = te_solver.InputData()
input_data.N_elements = N_elem
input_data.n_axi = n_node

# --- 物理场分布 (复制 6 份) ---
# np.tile 将 (50,) 的数据在行方向复制 6 次，变成 (6, 50)
input_data.Temitter   = np.tile(TE_data, (N_elem, 1))
input_data.Tcollector = np.tile(TC_data, (N_elem, 1))
input_data.Tcs        = np.full((N_elem, n_node), 600.0)
input_data.V_init     = np.full((N_elem, n_node), 0.2) # test.cpp 中是 vinput(..., 0.2) 或 1.2

# --- 几何分布 ---
# dlC = 0.377 / 50.0
dl_val = 0.377 / float(n_node)
input_data.dlE = np.full((N_elem, n_node), dl_val)
input_data.dlC = np.full((N_elem, n_node), dl_val)

# --- 标量参数 ---
# crossAreaE = 6.6669999999999997e-05
input_data.crossAreaE = np.full(N_elem, 6.667e-5)
# crossAreaC = 0.00010786
input_data.crossAreaC = np.full(N_elem, 1.0786e-4)

# sideAreaE = 0.00092855424159680002 * 25 / 50.0
sideE_val = 0.00092855424159680002 * 25.0 / float(n_node)
input_data.sideAreaE  = np.full((N_elem, n_node), sideE_val)

# sideAreaC = 0.00097592945800480005 * 25 / 50.0
sideC_val = 0.00097592945800480005 * 25.0 / float(n_node)
input_data.sideAreaC  = np.full((N_elem, n_node), sideC_val)

# U_init = 1.6, d_gap = 0.5 (from 'other' vector)
input_data.U_init     = np.full(N_elem, 1.6)
input_data.d_gap      = np.full(N_elem, 0.5)

# Itarget = 200.0
input_data.Itarget    = np.full(N_elem, 200.0)

# --- 导线电阻 (REC) ---
# test.cpp 中初始化为 {0,0,0,0}，且 T1, T2, T3 也都被设为 0
input_data.resistanceWire = np.zeros((N_elem, 4)) 

# --- 导线电压 (uwire) ---
# {0.8, 0.8, 0., 0.}
wireU_single = np.array([0.8, 0.8, 0.0, 0.0])
input_data.wireU = np.tile(wireU_single, (N_elem, 1))

# --- 全局控制模式 (这些只是初始值，稍后会覆盖) ---
input_data.mode = te_solver.CalculationMode.FixedResistance
input_data.target_val = 0.0044
input_data.I_total_init = 150.0

# ==============================================================================
# 3. 创建与计算
# ==============================================================================
print("正在调用 C++ 创建电路...")
circuit = te_solver.create_circuit(input_data)

# --- 设置与 testCircuit 一致的参数 ---
# C1->Utarget = 0.77 * 4;
# C1->Uout = 0.77 * 4;
# C1->Iout = 150.;
# C1->Rload = 0.0044;
circuit.Utarget = 0.77 * 4.0
circuit.Uout    = 0.77 * 4.0
circuit.Iout    = 150.0
circuit.Rload   = 0.0044
circuit.isFixedR = True
circuit.isFixedU = False

# [注意] test.cpp 中设置了 T1->isTail=true 和 T6->isTail=true
# 目前 binding 未暴露 isTail，如果计算结果偏差大，需要添加绑定后重新编译。
print("正在进行计算 (Fixed Resistance Mode)...")
t0 = time.time()
circuit.calc()
t1 = time.time()

print("="*40)
print(f"计算完成，耗时: {(t1-t0)*1000:.2f} ms")
print(f"C++ 计算结果:")
print(f"Uout: {circuit.Uout}")
print(f"Iout: {circuit.Iout}")
print("="*40)
