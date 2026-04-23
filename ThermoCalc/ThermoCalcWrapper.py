import os
import sys
import numpy as np
import time
import logging

from profiler import TEASAProfiler

logger = logging.getLogger(__name__)

# =========================================================
# [核心修复] 动态环境变量注入 (Dynamic Path Injection)
# 让 Wrapper 自己找到并加载同一目录下的 te_solver.pyd
# =========================================================
# 获取当前 ThermoCalcWrapper.py 所在的绝对路径
current_dir = os.path.dirname(os.path.abspath(__file__))

# 将该目录置于系统搜索路径的最高优先级 (Index 0)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# =========================================================
# 宽容导入机制 (Graceful Import)
# =========================================================
try:
    import te_solver  # 现在 Python 会在 sys.path[0] 也就是 current_dir 里找到它
    HAS_TE_SOLVER = True
except ImportError as e:
    HAS_TE_SOLVER = False
    logger.warning(f"⚠️ Module 'te_solver' (C++ PyBind11) not found in {current_dir}! "
                   f"ThermoCalcModel instantiation will fail if TEC coupling is enabled. "
                   f"Detail: {e}")


class ThermoCalcModel:
    """
    热离子能量转换器 (TEC) 的纯 Python 封装外壳
    负责处理几何参数分配、状态初始化、底层 C++ 对象通信以及结果提取。
    """

    def __init__(self, n_elements: int, n_nodes: int):
        self.N_elem = n_elements
        self.n_node = n_nodes

        # 保存底层的 C++ 计算核心对象
        self._circuit = None
        self._input_data = te_solver.InputData()

        # 1. 初始化基础维度参数
        self._input_data.N_elements = self.N_elem
        self._input_data.n_axi = self.n_node

        # 2. 提供默认几何与物理参数矩阵 (用户可通过后续方法覆盖)
        self._input_data.Tcs = np.full((self.N_elem, self.n_node), 600.0)
        self._input_data.V_init = np.full((self.N_elem, self.n_node), 0.2)

        # 默认几何分配 (以单根 TFE 的经验值为基准，后续提供自定义接口)
        dl_val = 0.507 / float(self.n_node)
        self._input_data.dlE = np.full((self.N_elem, self.n_node), dl_val)
        self._input_data.dlC = np.full((self.N_elem, self.n_node), dl_val)

        self._input_data.crossAreaE = np.full(self.N_elem, 6.667e-5)
        self._input_data.crossAreaC = np.full(self.N_elem, 1.0786e-4)

        sideE_val = 0.00092855424159680002 * 25.0 / float(self.n_node)
        self._input_data.sideAreaE = np.full(self.N_elem, sideE_val)

        sideC_val = 0.00097592945800480005 * 25.0 / float(self.n_node)
        self._input_data.sideAreaC = np.full(self.N_elem, sideC_val)

        self._input_data.U_init = np.full(self.N_elem, 1.6)
        self._input_data.d_gap = np.full(self.N_elem, 0.5)
        self._input_data.Itarget = np.full(self.N_elem, 200.0)

        # 默认导线连接电阻/电压
        self._input_data.resistanceWire = np.zeros((self.N_elem, 4))
        wireU_single = np.array([0.8, 0.8, 0.0, 0.0])
        self._input_data.wireU = np.tile(wireU_single, (self.N_elem, 1))

        # 电路模式默认值 (初始以定电阻模式为主)
        # self._input_data.mode = te_solver.CalculationMode.FixedResistance
        self._input_data.mode = te_solver.CalculationMode.FixedVoltage
        # self._input_data.target_val = 0.0044
        self._input_data.target_val = 0.89 * 34
        self._input_data.I_total_init = 284.0

        # --- 温度场占位符 ---
        # 初始默认处于 600K 均匀状态
        self._T_emitter = np.full((self.N_elem, self.n_node), 600.0)
        self._T_collector = np.full((self.N_elem, self.n_node), 600.0)

    def set_temperatures(self, T_em: np.ndarray, T_co: np.ndarray):
        """
        更新发射极与接收极温度场
        :param T_em: 形状为 (N_elem, n_node) 的 numpy 数组
        :param T_co: 形状为 (N_elem, n_node) 的 numpy 数组
        """
        if T_em.shape != (self.N_elem, self.n_node) or T_co.shape != (self.N_elem, self.n_node):
            raise ValueError(f"温度场维度错误。预期 ({self.N_elem}, {self.n_node})")

        self._T_emitter = T_em.copy()
        self._T_collector = T_co.copy()

        # 如果底层电路已经创建，直接通过 Pybind11 更新属性
        if self._circuit is not None:
            for i in range(self.N_elem):
                self._circuit.TECs[i].Temitter = self._T_emitter[i, :]
                self._circuit.TECs[i].Tcollector = self._T_collector[i, :]

    def setup_circuit_mode(self, mode_str: str, target_value: float, I_guess: float = 150.0):
        """
        设置全局电路求解模式
        :param mode_str: 'fixed_R' (定电阻) 或 'fixed_U' (定电压) 或 'fixed_I' (定电流)
        :param target_value: 目标阻值(Ohm) 或 电压(V) 或 电流(A)
        :param I_guess: 迭代猜测初始电流 (A)
        """
        if mode_str.lower() == 'fixed_r':
            self._input_data.mode = te_solver.CalculationMode.FixedResistance
        elif mode_str.lower() == 'fixed_u':
            self._input_data.mode = te_solver.CalculationMode.FixedVoltage
        elif mode_str.lower() == 'fixed_i':
            self._input_data.mode = te_solver.CalculationMode.FixedCurrent
        else:
            raise ValueError(f"不支持的模式: {mode_str}")

        self._input_data.target_val = target_value
        self._input_data.I_total_init = I_guess

    def build(self):
        """
        执行 C++ 层的实例化，准备开始计算。
        仅在首次运行或几何、连接方式发生改变时调用。
        """
        # 将最新的温度场推入 InputData
        self._input_data.Temitter = self._T_emitter
        self._input_data.Tcollector = self._T_collector

        # 构建电路
        self._circuit = te_solver.create_circuit(self._input_data)

        # 根据所选模式，初始化底层 circuit 的运行标记
        if self._input_data.mode == te_solver.CalculationMode.FixedResistance:
            self._circuit.isFixedR = True
            self._circuit.isFixedU = False
            self._circuit.Rload = self._input_data.target_val
            # 初始猜测
            self._circuit.Iout = self._input_data.I_total_init
            self._circuit.Uout = self._circuit.Iout * self._circuit.Rload
            self._circuit.Utarget = self._circuit.Uout

        elif self._input_data.mode == te_solver.CalculationMode.FixedVoltage:
            self._circuit.isFixedR = False
            self._circuit.isFixedU = True
            self._circuit.Utarget = self._input_data.target_val
            self._circuit.Uout = self._input_data.target_val
            self._circuit.Iout = self._input_data.I_total_init

    @TEASAProfiler.profile
    def calculate(self, verbose: bool = False) -> float:
        """
        触发底层计算 (瞬态单步更新时调用)
        :return: 计算耗时 (ms)
        """
        if self._circuit is None:
            self.build()

        t0 = time.time()
        self._circuit.calc()
        t1 = time.time()

        dt_ms = (t1 - t0) * 1000.0
        if verbose:
            print(f"[ThermoCalc] TEC System converged in {dt_ms:.2f} ms")
        return dt_ms

    def get_global_results(self):
        """ 获取系统级电路宏观结果 """
        if self._circuit is None:
            return None
        return {
            "Iout": self._circuit.Iout,
            "Uout": self._circuit.Uout,
            "Rload": self._circuit.Rload
        }

    def get_tec_results(self, idx: int):
        """
        提取特定 TFE 元件的详细物理场结果
        :param idx: 元件索引 (0 到 N_elem-1)
        :return: 包含各项数据的字典
        """
        if self._circuit is None or idx < 0 or idx >= self.N_elem:
            return None

        tec = self._circuit.TECs[idx]
        return {
            "I": tec.I,
            "U": tec.U,
            "J": np.array(tec.J),
            "V": np.array(tec.V),
            "UE": np.array(tec.UE),
            "UC": np.array(tec.UC),
            "rhoE": np.array(tec.rhoE),
            "rhoC": np.array(tec.rhoC),
            "IEsecSingle": np.array(tec.IEsecSingle),
            "ICsecSingle": np.array(tec.ICsecSingle),
            "phiE": np.array(tec.phiE),
            "phiC": np.array(tec.phiC),
            "Vd": np.array(tec.Vd),
            "TE": np.array(tec.Temitter),
            "TC": np.array(tec.Tcollector)
        }

    def set_tcs(self, tcs_val):
        """
        [通信接口 -> C++] 动态设置/更新铯池温度 (Tcs)

        :param tcs_val: 可以是全局统一的标量温度 [K]，也可以是与网格匹配的 (N_elem, n_node) 的 NumPy 数组
        """
        # --- 1. 数据形状安全对齐 (Safety Alignment) ---
        if np.isscalar(tcs_val):
            # 如果传入的是标量，直接利用 NumPy 的广播机制填满整个数组
            self._input_data.Tcs[:] = float(tcs_val)
        else:
            # 如果传入的是数组，必须进行严格的维度校验，防止 C++ 底层越界崩溃
            tcs_arr = np.asarray(tcs_val, dtype=float)
            if tcs_arr.shape != (self.N_elem, self.n_node):
                raise ValueError(
                    f"[ThermoCalcWrapper] Tcs array shape mismatch! Expected {(self.N_elem, self.n_node)}, got {tcs_arr.shape}")
            self._input_data.Tcs[:] = tcs_arr

        # --- 2. C++ 运行时热更新 (Hot Update) ---
        if self._circuit is not None:
            # TEASA 防御性编程：根据你底层 PyBind11 的具体绑定方式进行尝试
            if hasattr(self._circuit, 'set_tcs'):
                self._circuit.set_tcs(self._input_data.Tcs)
            elif hasattr(self._circuit, 'Tcs'):
                # 如果 C++ 端直接将 Tcs 暴露为可写的 numpy array / memoryview 属性
                self._circuit.Tcs = self._input_data.Tcs
            else:
                raise AttributeError(
                    "[ThermoCalcWrapper] The underlying C++ object lacks 'set_tcs' method or 'Tcs' property. Cannot update Tcs at runtime!")

    def set_rload(self, rload_val: float):
        """
        [通信接口 -> C++] 动态设置/更新系统的外部负载电阻

        :param rload_val: 外部负载电阻 [Ohm]
        """
        # --- 1. 更新初始化结构体 (兜底保护) ---
        # 兼容你 _input_data 中可能存在的各种命名习惯
        r_val = float(rload_val)
        if hasattr(self._input_data, 'Rload'):
            self._input_data.Rload = r_val
        elif hasattr(self._input_data, 'R_load'):
            self._input_data.R_load = r_val

        # --- 2. C++ 运行时热更新 (Hot Update) ---
        if self._circuit is not None:
            # 同样进行防御性探测
            if hasattr(self._circuit, 'set_rload'):
                self._circuit.set_rload(r_val)
            elif hasattr(self._circuit, 'Rload'):
                # 根据您代码中 get_global_results() 的实现，_circuit 应该暴露了 Rload 属性
                self._circuit.Rload = r_val
            else:
                raise AttributeError(
                    "[ThermoCalcWrapper] The underlying C++ object lacks 'set_rload' method or 'Rload' property. Cannot update load resistance at runtime!")
