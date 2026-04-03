import logging
from typing import List, Tuple

import numpy as np

from Components.BaseComponent import BaseComponent
from Components.basicComponents.TECPair import TECPair
from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel
from profiler import TEASAProfiler

logger = logging.getLogger(__name__)


class TECCircuitManager(BaseComponent):
    """
    TASTIN 系统层 - 热离子电热耦合电路管理器 (The Thermo-Electric Coordinator)

    采用【纯显式算子分裂】范式：
    在每个物理时间步的最开始 (pre_step)，根据当前时刻的电极表面温度，
    瞬间求解全局电路，并将等离子体热流与焦耳热强加给底层固体网格。
    """

    def __init__(self,
                 name: str,
                 tfe_list: List[TECPair],
                 Tcs_init: float = 600.0,
                 V_init: float = 0.2):
        """
        初始化电路管理器

        :param name: 组件名称 (SystemManager 识别用)
        :param tfe_list: 包含若干 TECPair 对象的列表
        :param Tcs_init: 铯池初始给定温度 [K]
        :param V_init: 初始给定电压估值 [V]
        """
        super().__init__(name)

        if not tfe_list:
            raise ValueError("TFE list cannot be empty!")

        self.tfe_list = tfe_list
        self.N_elem = len(tfe_list)
        self.n_node = tfe_list[0].n_node

        # 检查所有元件的轴向节点数是否一致
        for i, tfe in enumerate(self.tfe_list):
            if tfe.n_node != self.n_node:
                raise ValueError(f"TFE[{i}] '{tfe.name}' n_node ({tfe.n_node}) mismatches base ({self.n_node})!")

        # 实例化底层的 C++ 封装外壳
        self.circuit = ThermoCalcModel(n_elements=self.N_elem, n_nodes=self.n_node)

        # 初始化宏观参数矩阵 (切片赋值保证内存安全)
        self.circuit._input_data.Tcs[:] = Tcs_init
        self.circuit._input_data.V_init[:] = V_init

        # 预分配温度收集矩阵，避免每次步进时触发 Python 垃圾回收
        self._T_emit_matrix = np.zeros((self.N_elem, self.n_node))
        self._T_coll_matrix = np.zeros((self.N_elem, self.n_node))

    # =========================================================
    # SystemManager 组件接管接口 (Duck Typing)
    # =========================================================

    def get_solids(self) -> list:
        """ [接口] 榨取所有电极实体，交由 SystemManager 进行导热积分 """
        solids = []
        for tec in self.tfe_list:
            solids.extend([tec.emitter, tec.collector])
        return solids

    def get_couplers(self) -> list:
        """ [接口] 榨取所有间隙耦合器，交由 SystemManager 在 Picard 迭代中同步边界 """
        return [tec.tec_gap for tec in self.tfe_list]

    def pre_step(self, dt: float, current_time: float):
        """
        [核心时间步钩子]
        SystemManager.step() 在内部循环开始前会自动调用此函数。
        在这里执行显式电热计算，提供本时间步的所有源项。
        """
        # 如果需要基于 current_time 修改外部负载 R_load 等，可以在这里构造宏观参数字典
        # macro_params = {'R_load': new_r_load, 'Tcs': new_tcs}
        macro_params = {}

        # 执行电学快照计算与源项分发
        self.sync_thermo_electric(alpha=1.0, macro_params=macro_params)

    # =========================================================
    # 物理计算与通信逻辑
    # =========================================================

    def _compute_plasma_fluxes(self, tec_results: dict) -> Tuple[np.ndarray, np.ndarray]:
        """
        [物理核心预留接口]
        根据底层返回的电学结果，计算电子冷却和加热的面热流密度。
        注意：此处返回值必须是基于 Emitter 外表面积的面热流密度 [W/m^2]！
        """
        J_density = tec_results.get('J', np.zeros(self.n_node))  # [A/cm^2]
        J_density *= 1e4    # [A/m^2]

        phiE = tec_results.get('phiE', np.zeros(self.n_node))   # [V]
        phiC = tec_results.get('phiC', np.zeros(self.n_node))   # [V]
        Vd = tec_results.get('Vd', np.zeros(self.n_node))       # [V]

        TE = tec_results.get('TE', np.zeros(self.n_node))       # [K]

        # -----------------------------------------------------------------
        # 注意符号约定: Emitter失去热量用负值，Collector获得热量用正值。
        # -----------------------------------------------------------------

        # [临时占位符]
        # q_e_flux = -1.0 * J_density * 2.0
        # q_c_flux = 1.0 * J_density * 2.0

        # 对于发射极表面热流，需要获取发射极边界温度、发射极功函数、电流密度
        q_e_flux = -1.0 * J_density * (phiE + 2.0 * 8.617e-5 * TE)
        q_c_flux = 1.0 * J_density * (phiC + 2.0 * 8.617e-5 * TE + Vd)

        return q_e_flux, q_c_flux

    @TEASAProfiler.profile
    def sync_thermo_electric(self, alpha: float = 1.0, macro_params: dict = None):
        """
        执行底层 C++ 计算并将热源散发到组件中。
        :param macro_params:
        :param alpha: 亚松弛因子 (1.0 为纯显式。如果发生理查森震荡，可调小此值)
        """

        # --- 1. 瞬态宏观边界更新 ---
        if macro_params:
            if 'Tcs' in macro_params and hasattr(self.circuit, 'set_tcs'):
                self.circuit.set_tcs(macro_params['Tcs'])
            if 'R_load' in macro_params and hasattr(self.circuit, 'set_rload'):
                self.circuit.set_rload(macro_params['R_load'])

        # --- 2. Gather (收集边界温度场) ---
        for i, tec_pair in enumerate(self.tfe_list):
            T_e, T_c = tec_pair.get_gap_surface_temperatures()
            self._T_emit_matrix[i, :] = T_e
            self._T_coll_matrix[i, :] = T_c

        # 采用内置函数进行温度更改操作
        if hasattr(self.circuit, 'set_temperatures'):
            self.circuit.set_temperatures(self._T_emit_matrix, self._T_coll_matrix)
        else:
            logger.warning("Method 'set_temperatures' not found in ThermoCalcModel. Temperatures not updated.")

        # --- 3. 执行一次底层电路计算，防止不存在circuit时直接取值
        self.circuit.calculate(verbose=False)

        # 必须使用 [:] 保证写入底层 PyBind11 绑定的 C++ 内存块！(错误调用)
        # self.circuit._input_data._T_emitter[:] = self._T_emit_matrix
        # self.circuit._input_data._T_collector[:] = self._T_coll_matrix

        # --- 4. Scatter (提取结果与物理量转换，下发源项) ---
        for i, tec_pair in enumerate(self.tfe_list):
            res = self.circuit.get_tec_results(i)
            if res is None:
                logger.error(f"Failed to fetch TEC results for element {i}")
                continue

            # 1. 提取绝对电势分布 [V] 和 电阻率 [Ohm*m]
            UE_abs = res.get('UE', np.zeros(self.n_node))
            UC_abs = res.get('UC', np.zeros(self.n_node))
            rho_e = res.get('rhoE', np.ones(self.n_node) * 1e-6)
            rho_c = res.get('rhoC', np.ones(self.n_node) * 1e-6)

            # 2. 计算每个轴向节点的电势差 dU [V]
            # np.gradient 计算相邻节点的差值，取绝对值保证发热为正
            dU_e = np.abs(np.gradient(UE_abs))
            dU_c = np.abs(np.gradient(UC_abs))

            tec_pair.set_joule_heating(dU_emit=dU_e, rho_emit=rho_e,
                                       dU_coll=dU_c, rho_coll=rho_c,
                                       alpha=alpha)

            # B. 等离子体面热流下发
            q_e, q_c = self._compute_plasma_fluxes(res)
            tec_pair.update_plasma_heat_flux(q_e_flux=q_e, q_c_flux=q_c, alpha=alpha)
