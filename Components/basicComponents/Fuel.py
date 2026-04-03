import numpy as np
from typing import Optional, List

# 导入 Layer 2 的网格与基础导热求解器
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Base import SolidMaterial

from Materials.Solids.UO2 import UO2


class Fuel(HeatConduction2D):
    """
    TASTIN 组件层 (Layer 3) - 燃料芯块 (Fuel Pellet)

    本质：一个带有内热源和温度反馈的导热空心圆柱。
    架构：继承自底层 HeatConduction2D。它通过接受 cylindrical 网格来表达圆柱几何，
         并通过扩展鸭子类型接口，实现与 PointReactor (点堆) 的双向闭环。
    """

    def __init__(self,
                 name: str,
                 mesh: Mesh2D,
                 material: Optional[SolidMaterial] = None,
                 initial_temp: float = 600.0,
                 power_fraction: float = 1.0,
                 axial_power_profile: Optional[np.ndarray] = None,
                 # --- [TEASA 极简接口] 将几何分配剥离，仅接收交界面索引 ---
                 contact_resistance_interfaces: Optional[List[int]] = None,
                 axial_contact_resistance: float = 0.0):
        """
        初始化燃料组件

        :param contact_resistance_interfaces: [新增] 存在接触热阻的轴向交界面索引列表。例如 [4, 24]。
        :param axial_contact_resistance: [新增] 接触界面的总热阻值 [K/W]。
        """
        # 1. 缺省材料设置
        if material is None:
            material = UO2()

        # 2. 缓存热阻接口与数值
        self.contact_resistance_interfaces = contact_resistance_interfaces if contact_resistance_interfaces is not None else []
        self.axial_contact_resistance = axial_contact_resistance

        # 3. 提取体积参数，防止瞬态中反复调用
        self.nx, self.ny = mesh.n_x, mesh.n_y
        self.vols_flat = mesh.geom_data.volumes
        self.vols_2d = self.vols_flat.reshape(self.nx, self.ny)

        # 4. 挂载底层导热求解器 (首先调用，确立 self.mesh 等基础属性)
        super().__init__(name=name, mesh=mesh, material=material, initial_temp=initial_temp)
        self.power_fraction = power_fraction

        # [Bug修复点]：预存总体积标量，用于反应性反馈加权平均
        self.total_vol = np.sum(self.vols_flat)
        self.vols_layer = np.sum(self.vols_2d, axis=0)

        # 记录体积发热率分布 (W/m^3)
        self.Q_vol = np.zeros(self.N, dtype=float)
        self.power_allocation_weights = np.zeros(self.N, dtype=float)

        # 5. 设置轴向功率分布
        # [Bug修复点]：取消重复的冗余逻辑，直接调用单一入口，解决 len(None) 崩溃风险
        self.set_axial_power_profile(axial_power_profile)

    # --- [TEASA 极简重载] 拦截交界面热阻 ---
    def _compute_internal_resistance(self):
        """
        拦截父类的内部热阻计算，在指定的轴向交界面索引处，串联接触热阻。
        """
        super()._compute_internal_resistance()

        # 性能保护：不使用该功能时直接退出
        if self.axial_contact_resistance <= 1e-12 or not self.contact_resistance_interfaces:
            return

        # 提取面积并计算等效热阻系数
        A_y_inner = self.mesh.area_y_matrix[:, 1:-1]
        A_total = np.sum(A_y_inner[:, 0])
        r_c = self.axial_contact_resistance * A_total

        # 遍历外部传入的所有交界面索引，精准施加接触热阻
        for j_interface in self.contact_resistance_interfaces:
            j_int = int(j_interface)
            if 0 <= j_int < (self.ny - 1):
                G_orig = self.G_y_inner[:, j_int]
                R_orig = 1.0 / np.maximum(G_orig, 1e-30)

                A_ring = A_y_inner[:, j_int]
                R_contact_ring = r_c / np.maximum(A_ring, 1e-30)

                self.G_y_inner[:, j_int] = 1.0 / (R_orig + R_contact_ring)

    def set_axial_power_profile(self, profile_array: Optional[np.ndarray] = None):
        """
        设置/更新轴向功率分布，并预计算每个网格的功率分配权重矩阵。
        """
        if profile_array is None:
            # 缺省为均匀分布，每层分得 1.0 / ny
            profile_array = np.full(self.ny, 1.0 / self.ny)
        else:
            if len(profile_array) != self.ny:
                raise ValueError(f"轴向功率分布数组长度 ({len(profile_array)}) 必须与网格轴向节点数 ({self.ny}) 一致！")

        # 强制归一化 (防呆设计，防止总功率被错误缩放)
        total_profile = np.sum(profile_array)
        if total_profile > 1e-12:
            profile_array = profile_array / total_profile
        else:
            profile_array = np.zeros(self.ny)

        # 计算每个节点的权重
        weights_2d = np.zeros((self.nx, self.ny), dtype=float)
        for j in range(self.ny):
            if self.vols_layer[j] > 1e-12:
                weights_2d[:, j] = profile_array[j] * (self.vols_2d[:, j] / self.vols_layer[j])

        self.power_allocation_weights = weights_2d.flatten()

        # 严谨性校验：所有节点权重之和必须为 1.0
        if total_profile > 1e-12:
            assert np.isclose(np.sum(self.power_allocation_weights), 1.0), "功率分配权重归一化失败！"

    def set_nuclear_power(self, p_fiss: float, p_decay: float, p_total: float):
        """
        [多物理场输入接口] 接收点堆核功率，并下发至计算网格。
        """
        component_power = p_total * self.power_fraction

        # 极致优化：直接使用预计算权重进行分配，实现 O(1) 向量运算
        self.Q_source = component_power * self.power_allocation_weights

        # 加上微小值防止除0警告 (理论上体积不可能为0)
        self.Q_vol = self.Q_source / (self.vols_flat + 1e-20)

    def get_reactivity_feedback(self) -> float:
        """
        [多物理场输出接口] 计算并返回当前燃料组件引入的反应性反馈。
        """
        # [Bug修复点]：使用总体积 self.total_vol，避免标量除以数组引发维度爆炸
        T_avg = np.sum(self.T * self.vols_flat) / self.total_vol

        # 燃料温度反馈多项式
        T = T_avg
        rho_fb = (0.001360811
                  - 6.47927757e-6 * T
                  + 2.321231e-9 * (T ** 2)
                  - 3.52e-13 * (T ** 3))

        return rho_fb

