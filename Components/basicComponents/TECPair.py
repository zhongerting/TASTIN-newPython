import numpy as np
from typing import Optional

from Solvers.HeatConduction.Mesh import Mesh2D
from Components.basicComponents.Electord import Emitter, Collector
from Solvers.Couplers import TECCouple2D
from Materials.Base import SolidMaterial


class TECPair:
    """
    TASTIN 组件层 - 热离子电极对 (Thermionic Electrode Pair)

    物理职责：
    1. 封装 Emitter 和 Collector 的 2D 圆柱导热模型。
    2. 维护极间隙的传热/辐射耦合 (TECCouple2D)。
    3. 提供与底层放电模型 (ThermoCalcModel) 交互的标准化接口。
    4. 处理面积守恒缩放与焦耳体积热源的映射。
    """

    def __init__(self,
                 name: str,
                 L_node: float,
                 n_node: int,
                 R_e_in: float,
                 delta_e: float,
                 delta_gap: float,
                 delta_c: float,
                 n_rad_e: int,
                 n_rad_c: int,
                 mat_emitter: Optional[SolidMaterial] = None,
                 mat_collector: Optional[SolidMaterial] = None,
                 T_init_e: float = 1600.0,
                 T_init_c: float = 800.0,
                 k_gap_gas: float = 0.0,  # 间隙气体导热系数 [W/mK] (如 Cs 蒸气)
                 emissivity_e: float = 0.8,  # 发射极表面发射率
                 emissivity_c: float = 0.8):  # 接收极表面发射率

        self.name = name
        self.n_node = n_node
        self.L_node = L_node
        self.total_length = L_node * n_node
        self.n_rad_e = n_rad_e
        self.n_rad_c = n_rad_c

        # --- 1. 几何尺寸推算 ---
        self.R_e_out = R_e_in + delta_e
        self.R_c_in = self.R_e_out + delta_gap
        self.R_c_out = self.R_c_in + delta_c

        # --- 2. 严格遵循 Mesh2D 构建 2D 圆柱网格 ---
        # x_dim 为径向厚度，y_dim 为轴向总长
        mesh_e = Mesh2D(x_dim=delta_e, n_x=n_rad_e,
                        y_dim=self.total_length, n_y=n_node,
                        geometry_type='cylindrical', inner_radius=R_e_in)

        mesh_c = Mesh2D(x_dim=delta_c, n_x=n_rad_c,
                        y_dim=self.total_length, n_y=n_node,
                        geometry_type='cylindrical', inner_radius=self.R_c_in)

        # --- 3. 实例化电极实体 ---
        self.emitter = Emitter(name=f"{name}_Emitter", mesh=mesh_e,
                               material=mat_emitter, initial_temp=T_init_e)
        self.collector = Collector(name=f"{name}_Collector", mesh=mesh_c,
                                   material=mat_collector, initial_temp=T_init_c)

        # --- 4. 建立间隙耦合器 (TECCouple2D) ---
        # 严格遵守 TECCouple2D 构造函数:
        # obj1='right' 连 obj2='left'，所以 direction='right'
        self.tec_gap = TECCouple2D(
            obj1=self.emitter,
            obj2=self.collector,
            direction='right',
            gap_width=delta_gap,
            gas_conductivity=k_gap_gas,
            emissivity1=emissivity_e,
            emissivity2=emissivity_c
        )

        # 亚松弛历史状态缓存 (面热流密度 W/m^2)
        self._q_e_old = np.zeros(n_node)
        self._q_c_old = np.zeros(n_node)

        # --- 5. 开放外部边界引用 ---
        self.inner_boundary = self.emitter.boundaries['left']  # 供燃料发热耦合
        self.outer_boundary = self.collector.boundaries['right']  # 供冷却剂耦合

        # --- 6. 计算中间变量 ---
        # 6.1. 上一时刻电极表面热流密度
        self._q_joule_e_old = np.zeros(self.n_rad_e * self.n_node)
        self._q_joule_c_old = np.zeros(self.n_rad_e * self.n_node)

    def get_gap_surface_temperatures(self):
        """
        [通信接口 -> C++] 提取电极靠近间隙表面的当前温度分布
        """
        T_emit_surf, _ = self.emitter.boundaries['right'].get_coupling_surface_snapshot()
        T_coll_surf, _ = self.collector.boundaries['left'].get_coupling_surface_snapshot()
        return T_emit_surf, T_coll_surf

    def set_joule_heating(self, dU_emit: np.ndarray, rho_emit: np.ndarray,
                          dU_coll: np.ndarray, rho_coll: np.ndarray, alpha: float = 1.0):
        """
        [通信接口 <- C++] 根据底层电学结果，分配焦耳热体积源项
        利用公式: Q_vol = (dU)^2 / (rho * L^2)
        """
        safe_rho_e = np.maximum(rho_emit, 1e-12)
        safe_rho_c = np.maximum(rho_coll, 1e-12)

        # 1. 计算 1D 轴向体积发热率 [W/m^3]
        q_vol_e_1d = (dU_emit ** 2) / (safe_rho_e * self.L_node ** 2)
        q_vol_c_1d = (dU_coll ** 2) / (safe_rho_c * self.L_node ** 2)

        # 2. 将 1D 数组广播为 2D 矩阵 (x:径向, y:轴向)
        q_vol_e_2d = np.ones((self.n_rad_e, 1)) * q_vol_e_1d[np.newaxis, :]
        q_vol_c_2d = np.ones((self.n_rad_c, 1)) * q_vol_c_1d[np.newaxis, :]

        # 3. 展平并乘以体积得到总功率 [W]，严格适配 Electord 的要求
        q_watts_e_flat = q_vol_e_2d.flatten() * self.emitter.vols_flat
        q_watts_c_flat = q_vol_c_2d.flatten() * self.collector.vols_flat

        # 增加亚松弛保护
        q_watts_e_new = alpha * q_watts_e_flat + (1.0 - alpha) * self._q_joule_e_old
        q_watts_c_new = alpha * q_watts_c_flat + (1.0 - alpha) * self._q_joule_c_old

        self._q_joule_e_old = q_watts_e_new
        self._q_joule_c_old = q_watts_c_new

        self.emitter.set_joule_heating(q_watts_e_new)
        self.collector.set_joule_heating(q_watts_c_new)

    def update_plasma_heat_flux(self, q_e_flux: np.ndarray, q_c_flux: np.ndarray, alpha: float = 1.0):
        """
        [通信接口 <- C++] 施加等离子体放电产生的面热流密度
        参数 q_e_flux 和 q_c_flux 是基于 Emitter 外表面积的面热流 [W/m^2]
        """

        # 0. 传入的参数是 面热流密度 [W/m^2]

        # 1. 获取网格面的面积 (Emitter 的外表面积)
        A_emit_surf = self.emitter.boundaries['right'].area

        # 2. 亚松弛更新面热流密度 [W/m^2]
        q_e_new = alpha * q_e_flux + (1.0 - alpha) * self._q_e_old
        q_c_new = alpha * q_c_flux + (1.0 - alpha) * self._q_c_old

        self._q_e_old = q_e_new
        self._q_c_old = q_c_new

        # 3. 转换为纯热流功率 [W]
        # 因为 q_c_flux 也是基于 Emitter 面积算出来的，
        # 所以它们乘以 Emitter 面积后，就是完美的守恒传输功率 (Watts)。
        # 根据 TECCouple2D 定义: Emitter侧流出输入负值，Collector侧流入输入正值。
        Q_e_watts = q_e_new * A_emit_surf
        Q_c_watts = q_c_new * A_emit_surf

        # 4. 严格调用 TECCouple2D 提供的专属设值接口
        self.tec_gap.set_tec_sources(Q_emitter=Q_e_watts, Q_collector=Q_c_watts)

    def sync_and_step(self, dt: float):
        """时间步积分接口"""
        self.tec_gap.sync()
        self.emitter.step(dt)
        self.collector.step(dt)
