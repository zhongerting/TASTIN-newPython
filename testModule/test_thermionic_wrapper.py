import numpy as np
import time
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, '..'))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from profiler import TEASAProfiler

# 1. 获取当前测试脚本所在的绝对路径 (.../TASTIN_Project/Tests)
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. 推导出 ThermoCalc 模块所在的绝对路径
thermo_calc_dir = os.path.join(current_dir, '..', 'ThermoCalc')
thermo_calc_dir = os.path.abspath(thermo_calc_dir)

# 3. 将该路径加入系统搜索路径
if thermo_calc_dir not in sys.path:
    sys.path.insert(0, thermo_calc_dir)

# 导入我们刚刚编写的封装类
from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

# ==============================================================================
# 1. 原始数据准备 (直接复制自 test_real_case_v3.py)
# ==============================================================================
TE_data = [
    1491.29, 1500.11, 1516, 1537.34, 1562.78, 1591.18, 1621.63, 1653.35, 1685.71, 1718.16,
    1750.24, 1781.57, 1811.82, 1840.71, 1868, 1893.49, 1917, 1938.4, 1957.56, 1974.4,
    1988.83, 2000.8, 2010.27, 2017.21, 2021.62, 2023.48, 2022.83, 2019.69, 2014.09, 2006.08,
    1995.73, 1983.11, 1968.29, 1951.39, 1932.5, 1911.77, 1889.32, 1865.34, 1840.01, 1813.57,
    1786.27, 1758.44, 1730.45, 1702.77, 1675.95, 1650.68, 1627.8, 1608.39, 1593.75, 1585.55
]

TC_data = [
    760.774, 762.065, 763.823, 765.81, 767.951, 770.213, 772.576, 775.027, 777.551, 780.139,
    782.779, 785.462, 788.178, 790.92, 793.677, 796.445, 799.211, 801.974, 804.72, 807.451,
    810.151, 812.824, 815.452, 818.044, 820.578, 823.068, 825.488, 827.858, 830.144, 832.379,
    834.517, 836.603, 838.58, 840.509, 842.317, 844.084, 845.715, 847.323, 848.779, 850.236,
    851.526, 852.851, 853.995, 855.227, 856.262, 857.463, 858.453, 859.715, 860.732, 862.062
]

N_elem = 6
n_node = len(TE_data)

# ==============================================================================
# 2. 模型初始化与配置
# ==============================================================================
print("正在通过 ThermoCalcModel 初始化系统...")
model = ThermoCalcModel(n_elements=N_elem, n_nodes=n_node)

# 设置为定电阻模式: Rload = 0.0044 Ohm, 初始猜测电流 = 150.0 A
model.setup_circuit_mode(mode_str='fixed_r', target_value=0.0044, I_guess=150.0)

# 转换数据格式为 (N_elem, n_node) 的 NumPy 矩阵
TE_mat = np.tile(TE_data, (N_elem, 1))
TC_mat = np.tile(TC_data, (N_elem, 1))

# 下发温度边界
model.set_temperatures(TE_mat, TC_mat)

# ==============================================================================
# 3. 初始状态计算
# ==============================================================================
print("正在进行初始稳态计算...")
model.calculate(verbose=True)


# 验证函数提取
def print_report(mdl, title="详细物理场结果验证报告"):
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)

    # A. 全局电路结果
    global_res = mdl.get_global_results()
    print(f"[全局电路结果]")
    print(f"  Iout (总电流)    : {global_res['Iout']:.6f} A")
    print(f"  Uout (总电压)    : {global_res['Uout']:.6f} V")
    print(f"  Rload (负载)     : {global_res['Rload']:.6f} Ohm")
    print("-" * 60)

    # B. 单个元件结果验证
    check_indices = [0, 2]

    def print_vector_info(name, v):
        print(f"    {name:<12} | Min: {v.min():.4e} | Max: {v.max():.4e} | Mean: {v.mean():.4e}")
        mid = len(v) // 2
        print(f"                 | [0]: {v[0]:.4e}   | [{mid}]: {v[mid]:.4e}   | [-1]: {v[-1]:.4e}")

    for idx in check_indices:
        tec_res = mdl.get_tec_results(idx)
        print(f"\n[元件 TECs[{idx}] 详细数据]")
        print(f"  1. 标量状态:")
        print(f"    Current I    : {tec_res['I']:.6f} A")
        print(f"    Voltage U    : {tec_res['U']:.6f} V")
        print(f"  2. 物理场分布:")
        print_vector_info("J (A/cm2)", tec_res['J'])
        print_vector_info("V (Volts)", tec_res['V'])
        print_vector_info("UE (Volts)", tec_res['UE'])
        print_vector_info("UC (Volts)", tec_res['UC'])
        print_vector_info("rhoE (Ohm.m)", tec_res['rhoE'])
        print_vector_info("rhoC (Ohm.m)", tec_res['rhoC'])
        print_vector_info("IEsecSingle", tec_res['IEsecSingle'])
        print_vector_info("ICsecSingle", tec_res['ICsecSingle'])
        print("-" * 30)


print_report(model, "初始状态验证报告 (Base Temperature)")

# ==============================================================================
# 4. 模拟瞬态步进: 发射极温度整体抬升 50K
# ==============================================================================
print("\n>>> 模拟瞬态扰动：发射极温度整体增加 50K...")
TE_mat_new = TE_mat + 50.0

# 只需要两行代码即可完成更新和重新计算！
model.set_temperatures(TE_mat_new, TC_mat)

model.calculate(verbose=True)

for i in range(50):
    model.calculate(verbose=True)

print_report(model, "瞬态扰动后结果验证报告 (TE + 50K)")
print("=" * 60)
print("测试完毕！")

TEASAProfiler.report()
