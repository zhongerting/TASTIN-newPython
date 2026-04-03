import numpy as np
import logging
from typing import Dict, List, Any, Optional

from Components.BaseComponent import BaseComponent
from Components.TFEUnit import TFEUnit
from Solvers.Couplers import SolidSolidCouple2D
from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D

logger = logging.getLogger(__name__)


class ReactorCore(BaseComponent):
    """
    TASTIN 堆芯宏观容器 (Super-Component / Core Vessel)
    集成多物理场调度机制：包含点堆动力学和热电耦合电路 (TEC) 的虚拟复制与降维映射。
    """

    def __init__(self,
                 name: str,
                 tfe_dict: Dict[str, TFEUnit],
                 tfe_multipliers: Dict[str, int],
                 # === 全局慢化剂基体参数 ===
                 mod_meshes: Optional[List[Mesh2D]] = None,  # 每一圈慢化剂的网格列表
                 mod_material: Optional[Any] = None,  # 慢化剂统一物性 (如 ZrH)
                 ring_mapping: Optional[Dict[str, int]] = None,  # TFE 到全局慢化剂圈层的映射字典
                 T_space: float = 250.0,
                 alpha_tec: float = 0.5,
                 enable_tec_coupled: bool = True):  # [新增] 电计算总开关
        """
        初始化反应堆堆芯组件

        :param name: 堆芯名称
        :param tfe_dict: 物理计算实体字典 (如 {'TFE_Center': tfe0, 'TFE_Ring1': tfe1})
        :param tfe_multipliers: 虚拟分身映射字典 (如 {'TFE_Center': 1, 'TFE_Ring1': 6})
                                指示在电路求解中该实体代表几根相同的组件串联
        :param alpha_tec: 电热耦合亚松弛因子 (避免理查森发射引起的温度震荡，推荐 0.2~0.8)
        :param enable_tec_coupled: 是否开启热电耦合计算。若为 False，则关闭所有电路反馈。
        :param mod_meshes: 慢化剂网格列表
        :param mod_material: 慢化剂材料列表。
        :param ring_mapping: 映射字典，指示哪个 TFE 对应慢化剂的哪一圈 (0, 1, 2, 3)
                             例如: {'TFE_Center': 0, 'TFE_Ring1': 1, 'TFE_Ring2': 2, 'TFE_Ring3': 3}
        """
        super().__init__(name)

        self.tfes = tfe_dict
        self.tfe_multipliers = tfe_multipliers
        self.alpha_tec = alpha_tec
        self.enable_tec_coupled = enable_tec_coupled

        # --- 多物理场模型 ---
        self.point_reactor = None

        # =========================================================
        # [核心] 构建 ThermoCalcModel 的虚拟阵列
        # =========================================================
        # 校验字典键的一致性
        for key in tfe_dict.keys():
            if key not in tfe_multipliers:
                raise ValueError(f"TFE '{key}' 必须在 tfe_multipliers 中指明复制数量！")

        self.total_virtual_elements = sum(tfe_multipliers.values())
        self.n_nodes = list(self.tfes.values())[0].mesh.n_axial

        if self.total_virtual_elements > 0:
            self.thermo_calc = ThermoCalcModel(n_elements=self.total_virtual_elements,
                                               n_nodes=self.n_nodes)
            status = "ENABLED" if self.enable_tec_coupled else "DISABLED"
            logger.info(f"ReactorCore '{name}': Built ThermoCalc circuit with "
                        f"{len(self.tfes)} physical TFE(s) representing "
                        f"{self.total_virtual_elements} virtual element(s) in series. (TEC is {status})")
        else:
            self.thermo_calc = None

        self.extra_solids = {}
        self.extra_couplers = {}

        # =========================================================
        # 2. 构建全局慢化剂分层基体 (Global Moderator Rings)
        # =========================================================
        self.ring_mapping = ring_mapping
        self.has_global_moderator = False
        self.mod_rings: List[HeatConduction2D] = []

        if mod_meshes is not None and mod_material is not None and ring_mapping is not None:
            self.has_global_moderator = True

            # 2.1 逐层实例化全局慢化剂环
            for i, mesh in enumerate(mod_meshes):
                ring = HeatConduction2D(
                    name=f"{self.name}_ModRing_{i}",
                    mesh=mesh,
                    material=mod_material,
                    initial_temp=743.0  # 全局稳态初始温度
                )
                self.mod_rings.append(ring)
                # 注册到大字典中，SystemManager 会自动接管并求解它们
                self.extra_solids[f'mod_ring_{i}'] = ring

            logger.info(f"ReactorCore '{name}': Built {len(self.mod_rings)} Global Moderator Rings.")

            # 2.2 构建层间径向固体导热耦合 (俄罗斯套娃串联)
            for i in range(len(self.mod_rings) - 1):
                coupler_name = f'mod_couple_{i}_{i + 1}'
                self.extra_couplers[coupler_name] = SolidSolidCouple2D(
                    obj1=self.mod_rings[i],
                    obj2=self.mod_rings[i + 1],
                    direction='right'
                )
            logger.info(
                f"ReactorCore '{name}': Established {len(self.mod_rings) - 1} Solid-Solid Couplers between rings.")

            # 2.3 [外部边界] 为最外圈慢化剂挂载太空辐射散热边界
            outermost_ring = self.mod_rings[-1]
            outermost_boundary = outermost_ring.boundaries['right']

            # 直接使用你 Boundary.py 中的 DynamicRadiationFluxBC
            outermost_boundary.add_dynamic_radiation_condition(
                emissivity=0.05,  # 反应堆外壳的典型黑度
                bare_area_array=outermost_boundary.area,  # 取最外层圆柱面积
                T_env=T_space
            )
            logger.info(f"ReactorCore '{name}': Attached Space Radiation BC (T_env={T_space}K) to the outermost ring.")

        # =========================================================
        # 3. [性能优化] 热电计算的更新频率控制
        # =========================================================
        self.thermo_update_interval = 0.8  # [s] 设定每 0.5 秒更新一次
        self._last_thermo_update_time = -999.0  # 初始化为负数，确保 t=0 时绝对会触发第一次计算

    def setup_tec_circuit(self, mode_str: str, target_value: float, I_guess: float = 150.0):
        """[封装接口] 统一设置底层的热电计算电路模式"""
        if self.thermo_calc is not None:
            self.thermo_calc.setup_circuit_mode(mode_str, target_value, I_guess)

    def get_solids(self) -> list:
        all_solids = list(self.extra_solids.values())
        for tfe in self.tfes.values():
            all_solids.extend(tfe.get_solids())
        return all_solids

    def get_couplers(self) -> list:
        all_couplers = list(self.extra_couplers.values())
        for tfe in self.tfes.values():
            all_couplers.extend(tfe.get_couplers())
        return all_couplers

    # =========================================================================
    # 多物理场通信核心：Pre-Step (计算并散发) & Post-Step (收集并广播)
    # =========================================================================

    def pre_step(self, dt: float, current_time: float):
        """
        [求解前钩子] Solve & Reduce: 
        求解电学状态，提取分身的代表结果，转换为源项施加到各物理 TFE 上。
        """
        # --- 1. 热电耦合电路计算 ---
        # 1.1 判断是否到达热电计算更新周期

        # 只有当开关开启且电路模型存在时才执行计算
        if self.enable_tec_coupled and self.thermo_calc is not None:
            time_since_last_update = current_time - self._last_thermo_update_time
            if time_since_last_update >= self.thermo_update_interval or self.has_global_moderator < 0:
                # 执行底层 C++ 求解
                self.thermo_calc.calculate(verbose=False)
                # 更新时间戳，记录这次执行的时刻
                self._last_thermo_update_time = current_time
            else:
                # 命中缓存：跳过计算
                # 这里必须确保你的组件能够正确地“沿用”上一次计算产生的源项（如电子冷却热流、焦耳热等）
                pass

            # Reduce：提取并下发源项
            idx = 0
            for tfe_name, mult in self.tfe_multipliers.items():
                tfe = self.tfes[tfe_name]

                # 由于串联且温度一致，同一物理批次的分身计算结果绝对一致
                # 我们只需要提取该批次的第一个分身(idx)的结果代表全体
                res = self.thermo_calc.get_tec_results(idx)
                idx += mult  # 跳过这些分身，游标指向下一批次的第一个

                if res is None:
                    continue

                # 提取焦耳热相关变量
                UE_abs = res.get('UE', np.zeros(self.n_nodes))
                UC_abs = res.get('UC', np.zeros(self.n_nodes))
                rho_e = res.get('rhoE', np.ones(self.n_nodes) * 1e-6)
                rho_c = res.get('rhoC', np.ones(self.n_nodes) * 1e-6)

                dU_e = np.abs(np.gradient(UE_abs))
                dU_c = np.abs(np.gradient(UC_abs))

                # 提取等离子体热流相关变量
                J_density = res.get('J', np.zeros(self.n_nodes)) * 1e4  # [A/cm^2] -> [A/m^2]
                phiE = res.get('phiE', np.zeros(self.n_nodes))
                phiC = res.get('phiC', np.zeros(self.n_nodes))
                Vd = res.get('Vd', np.zeros(self.n_nodes))
                TE = res.get('TE', np.zeros(self.n_nodes))

                q_e_flux = -1.0 * J_density * (phiE + 2.0 * 8.617e-5 * TE)
                q_c_flux = 1.0 * J_density * (phiC + 2.0 * 8.617e-5 * TE + Vd)

                # 下发给 TFEUnit 内部的源项接口
                tfe.update_electric_fields(dU_emit=dU_e, rho_emit=rho_e,
                                           dU_coll=dU_c, rho_coll=rho_c,
                                           alpha=self.alpha_tec)

                tfe.update_plasma_flux(q_e_flux=q_e_flux, q_c_flux=q_c_flux,
                                       alpha=self.alpha_tec)

        # --- 2. 提取虚拟慢化剂边界热流，注入全局慢化剂 ---
        if self.has_global_moderator:
            # 清空全局慢化剂的热源
            for ring in self.mod_rings:
                ring.use_external_source_buffer = True
                ring.Q_source[:] = 0.0
            for tfe_name, tfe in self.tfes.items():
                if tfe_name in self.ring_mapping:
                    ring_idx = self.ring_mapping[tfe_name]
                    mult = self.tfe_multipliers[tfe_name]
                    ring = self.mod_rings[ring_idx]

                    # A. 提取虚拟慢化剂最外层流出的热流 [W]
                    virtual_mod = tfe.solids['moderator']
                    # BoundaryRegion.current_flux 是“流入”节点的热流，所以流出为负号
                    q_flux_out = -virtual_mod.boundaries['right'].current_flux

                    # B. 转换为 M 根 TFE 排出的总功率 [W]
                    # area_out = virtual_mod.boundaries['right'].area
                    q_watts_total = q_flux_out * mult

                    # C. 转换为全局慢化剂对应圈层的体积热源 [W/m^3] 并注入
                    nx, ny = ring.shape_nodes
                    vols_2d = ring.mesh.geom_data.volumes.reshape(nx, ny)
                    vol_axial = np.sum(vols_2d, axis=0)

                    # 1. 算出该高度层均匀的体积发热率 [W/m^3]
                    q_vol_1d = q_watts_total / np.maximum(vol_axial, 1e-12)
                    q_vol_2d = np.ones((nx, 1)) * q_vol_1d[np.newaxis, :]

                    # 2. 【核心修正】乘以每个单独网格的体积，转化为绝对功率 [W]
                    q_watts_2d = q_vol_2d * vols_2d

                    # 3. 注入求解器源项
                    ring.Q_source += q_watts_2d.flatten()

        # --- 3. 向下传递钩子 ---
        for tfe in self.tfes.values():
            if hasattr(tfe, 'pre_step'):
                tfe.pre_step(dt, current_time)

    def post_step(self, dt: float, current_time: float):
        """
        [求解后钩子] Gather & Broadcast:
        从收敛的固体网格中收集表面温度，广播复制给所有分身，为下一时间步做准备。
        """
        # --- 1. 先向下传递收尾工作 ---
        for tfe in self.tfes.values():
            if hasattr(tfe, 'post_step'):
                tfe.post_step(dt, current_time)

        # --- 2. 收集温度并推流给电路模型 ---
        # [修改] 只有当开关开启且电路模型存在时才进行温度同步
        if self.enable_tec_coupled and self.thermo_calc is not None:
            T_em_matrix = np.zeros((self.total_virtual_elements, self.n_nodes))
            T_co_matrix = np.zeros((self.total_virtual_elements, self.n_nodes))

            idx = 0
            for tfe_name, mult in self.tfe_multipliers.items():
                tfe = self.tfes[tfe_name]

                # 提取真实元件的发射极外表面温度和接收极内表面温度
                T_e = tfe.solids['emitter'].boundaries['right'].T_surface
                T_c = tfe.solids['collector'].boundaries['left'].T_surface

                # 广播 (Broadcast) 给其名下的所有虚拟分身
                T_em_matrix[idx: idx + mult, :] = T_e
                T_co_matrix[idx: idx + mult, :] = T_c

                idx += mult

            # 更新到 C++ 层
            self.thermo_calc.set_temperatures(T_em_matrix, T_co_matrix)

        # --- 3. 先向下传递收尾工作 ---
        if self.has_global_moderator:
            for tfe_name, tfe in self.tfes.items():
                if tfe_name in self.ring_mapping:
                    ring_idx = self.ring_mapping[tfe_name]
                    ring = self.mod_rings[ring_idx]

                    # 提取该圈全局慢化剂的 2D 温度场，取径向平均
                    nx, ny = ring.shape_nodes
                    T_2d = ring.T.reshape(nx, ny)
                    T_avg_axial = np.mean(T_2d, axis=0)

                    # 赋给 TFE 的边界温度缓存，供下一个 pre_step 使用
                    tfe.boundary_data.moderator_temperature[:] = T_avg_axial

    # ==========================================
    # 断点续算路由 (Restart Routing)
    # ==========================================

    def get_state_dict(self, prefix: str) -> dict:
        state = {}
        # [新增] 保存系统当前的宏观开关状态
        state[f"{prefix}/enable_tec_coupled"] = np.array([self.enable_tec_coupled], dtype=bool)

        for tfe_name, tfe in self.tfes.items():
            if hasattr(tfe, 'get_state_dict'):
                state.update(tfe.get_state_dict(prefix=f"{prefix}/TFEs/{tfe_name}"))

        # 存入上次更新时间
        state[f"{prefix}/_last_thermo_update_time"] = np.array([self._last_thermo_update_time])

        return state

    def load_state_dict(self, data: dict, prefix: str):
        # [新增] 恢复系统当前的宏观开关状态
        flag_key = f"{prefix}/enable_tec_coupled"
        if flag_key in data:
            self.enable_tec_coupled = bool(data[flag_key][0])

        for tfe_name, tfe in self.tfes.items():
            if hasattr(tfe, 'load_state_dict'):
                tfe.load_state_dict(data, prefix=f"{prefix}/TFEs/{tfe_name}")

        key = f"{prefix}/_last_thermo_update_time"
        if key in data:
            self._last_thermo_update_time = float(data[key][0])
