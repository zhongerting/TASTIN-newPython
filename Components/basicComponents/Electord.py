import numpy as np
from typing import Optional

# 导入底层导热求解器与网格
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Base import SolidMaterial

# 引入发射极缺省材料 MoNb (钼铌合金)
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum


class Electrode(HeatConduction2D):
    """
    TASTIN 组件层 (Layer 3) - 电极基类 (Electrode)

    本质：一个带有内热源(焦耳热)的导热空心圆柱管。
    作用：作为 Emitter 和 Collector 的通用基类，处理共有的焦耳热空间映射逻辑。
          具体的材料缺省值和反应性反馈关联式由继承它的子类实现。
    """

    def __init__(self,
                 name: str,
                 mesh: Mesh2D,
                 material: SolidMaterial,
                 initial_temp: float = 600.0):
        """
        初始化电极基类

        :param name: 组件名称
        :param mesh: 二维网格对象 (需使用 geometry_type='cylindrical')
        :param material: 电极物性对象 (必须由子类指定或由用户传入)
        :param initial_temp: 初始均匀温度 [K]
        """
        # 1. 挂载底层导热求解器
        super().__init__(mesh=mesh, material=material, initial_temp=initial_temp)
        self.name = name

        # 2. 提取体积参数，防止瞬态中反复调用
        self.vols_flat = self.mesh.geom_data.volumes
        self.total_vol = np.sum(self.vols_flat)
        self.nx, self.ny = self.mesh.shape_nodes

        # 3. 记录体积发热率分布 (W/m^3)，提供给后处理提取
        self.Q_vol = np.zeros(self.N, dtype=float)

    def set_joule_heating(self, q_joule_array: np.ndarray):
        """
        [多物理场输入接口] 接收空间分布的焦耳热源。
        由 TFE_Unit 在调用完 ThermoCalcModel 计算后统一调用。

        :param q_joule_array: 展平的 1D 数组，包含每个网格节点分配到的焦耳热功率 [W]
                              (长度必须等于 self.N)
        """
        if len(q_joule_array) != self.N:
            raise ValueError(
                f"焦耳热数组长度 ({len(q_joule_array)}) 必须与电极 '{self.name}' 的网格节点数 ({self.N}) 一致！")

        # 1. 使用 [:] 原地切片赋值，避免改变 NumPy 数组的底层内存指针
        self.Q_source[:] = q_joule_array

        # 2. [核心修复点] 激活外部热源标记！
        # 告诉 BaseHeatConduction：这个热源由外部显式维护，绝对不要在 _update_sources 中将其清零！
        self.use_external_source_buffer = True

        # 3. 同步更新体积热源项 Q_vol [W/m^3] (工程可视化用)
        self.Q_vol = self.Q_source / (self.vols_flat + 1e-20)

        # # 极致优化：直接将焦耳热数组赋给底层微分方程的热源项
        # self.Q_source = q_joule_array
        #
        # # 同步更新体积热源项 Q_vol [W/m^3] (工程可视化用)
        # self.Q_vol = self.Q_source / (self.vols_flat + 1e-20)

    def get_reactivity_feedback(self) -> float:
        """
        [多物理场输出接口] 计算并返回当前电极组件引入的反应性反馈。

        TEASA 架构限制：电极基类没有明确的物理反馈多项式，必须由子类(Emitter/Collector)实现。
        如果在子类中忘记实现，系统将在这里抛出错误以防止物理失真。
        """
        raise NotImplementedError(
            f"电极基类 '{self.name}' 不能直接提供反应性反馈，必须由 Emitter 或 Collector 子类实现！")


# =============================================================================
# 2. 发射极子类 (Emitter)
# =============================================================================
class Emitter(Electrode):
    """
    TASTIN 组件层 (Layer 3) - 发射极 (Emitter)

    继承自 Electrode 基类。
    特性：
        1. 缺省材料为 MoNb (钼铌合金)。
        2. 实现发射极专属的温度-反应性反馈多项式。
    """

    def __init__(self,
                 name: str,
                 mesh: Mesh2D,
                 material: Optional[SolidMaterial] = None,
                 initial_temp: float = 600.0):
        """
        初始化发射极组件
        """
        # 缺省材料设置：如果没有外部传入材料，则使用 MoNb
        if material is None:
            material = MoNb()

        # 调用父类 (Electrode) 的初始化
        super().__init__(name=name, mesh=mesh, material=material, initial_temp=initial_temp)

    def get_reactivity_feedback(self) -> float:
        """
        [多物理场输出接口] 计算并返回当前发射极组件引入的反应性反馈。
        """
        # 1. 采用预存的总体积进行快速体积加权平均
        T_avg = np.sum(self.T * self.vols_flat) / self.total_vol

        # 2. 发射极温度反馈多项式 (关联式提取自 TOPAZ 反应堆原 Fortran 源码)
        T = T_avg
        rho_fb = 1.e-4 * (3.46455e-6 * (T ** 2) - 0.03232167 * T + 0.74202216)

        return rho_fb

# =============================================================================
# 3. 接收极子类 (Collector)
# =============================================================================
class Collector(Electrode):
    """
    TASTIN 组件层 (Layer 3) - 接收极 (Collector)

    继承自 Electrode 基类。
    特性：
        1. 缺省材料为 Molybdenum (钼)。
        2. 实现发射极专属的温度-反应性反馈多项式。
    """

    def __init__(self,
                 name: str,
                 mesh: Mesh2D,
                 material: Optional[SolidMaterial] = None,
                 initial_temp: float = 600.0):
        """
        初始化发射极组件
        """
        # 缺省材料设置：如果没有外部传入材料，则使用 MoNb
        if material is None:
            material = Molybdenum()

        # 调用父类 (Electrode) 的初始化
        super().__init__(name=name, mesh=mesh, material=material, initial_temp=initial_temp)

    def get_reactivity_feedback(self) -> float:
        """
        [多物理场输出接口] 计算并返回当前接收极组件引入的反应性反馈。
        """

        # TODO 修改接收极反应性反馈函数
        # 1. 采用预存的总体积进行快速体积加权平均
        T_avg = np.sum(self.T * self.vols_flat) / self.total_vol

        # 2. 发射极温度反馈多项式 (关联式提取自 TOPAZ 反应堆原 C++ 源码)
        T = T_avg
        rho_fb = 1.e-4 * (3.46455e-6 * (T ** 2) - 0.03232167 * T + 0.74202216)

        return rho_fb
