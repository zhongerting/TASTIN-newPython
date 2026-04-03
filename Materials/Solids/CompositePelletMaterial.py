import numpy as np
from typing import Union
from Materials.Base import SolidMaterial


class CompositePelletMaterial(SolidMaterial):
    """
    复合燃料芯块材料类 (Composite Pellet Material)

    采用“拦截与重组”策略：
    对外伪装成一种单一材料，接收扁平化的温度场数组；
    对内根据网格的轴向切片（y方向），将温度场分发给 UO2 和 BeO 材料对象，
    计算完成后再拼装成完整的扁平数组返回给底层导热求解器。
    """

    def __init__(self,
                 fuel_mat: SolidMaterial,
                 reflector_mat: SolidMaterial,
                 shape_nodes: tuple,
                 n_lower: int,
                 n_active: int,
                 n_upper: int):
        """
        :param fuel_mat: 活性区燃料材料实例 (如 UO2)
        :param reflector_mat: 上下反射层材料实例 (如 BerylliumOxide)
        :param shape_nodes: 导热网格的形状 (nx, ny)，其中 nx 为径向节点，ny 为轴向节点
        :param n_lower: 下反射层占用的轴向网格数
        :param n_active: 活性区占用的轴向网格数
        :param n_upper: 上反射层占用的轴向网格数
        """
        # 动态生成名称
        super().__init__(name=f"Composite_{fuel_mat.name}_{reflector_mat.name}")

        self.fuel_mat = fuel_mat
        self.reflector_mat = reflector_mat
        self.shape_nodes = shape_nodes
        self.nx, self.ny = shape_nodes

        # 逻辑校验 1：确保分段网格数与总网格数严格匹配
        if n_lower + n_active + n_upper != self.ny:
            raise ValueError(f"[CompositePelletMaterial] 轴向节点数不匹配: "
                             f"{n_lower} + {n_active} + {n_upper} != 设定的总轴向节点 {self.ny}")

        # 轴向切片的分界索引
        self.idx_1 = n_lower  # 下反射层结束，活性区开始
        self.idx_2 = n_lower + n_active  # 活性区结束，上反射层开始

        # 内存优化：预分配输出缓冲区，避免在隐式积分器（每秒调用上万次）中产生垃圾回收开销
        self._out_buffer = np.zeros(self.shape_nodes, dtype=float)

    def _compute_property(self, T: Union[float, np.ndarray], prop_name: str) -> np.ndarray:
        """
        内部通用的物性分段计算引擎
        """
        # 逻辑校验 2：兼容标量输入（HeatConduction 在 initial_temp 初始化时可能会传入标量）
        if np.isscalar(T):
            T_2d = np.full(self.shape_nodes, T, dtype=float)
        else:
            T_2d = np.asarray(T).reshape(self.shape_nodes)

        out_2d = self._out_buffer
        out_2d.fill(0.0)

        # 1. 下反射层 (Lower Reflector) -> 轴向切片: [:, 0 : idx_1]
        if self.idx_1 > 0:
            func = getattr(self.reflector_mat, prop_name)
            out_2d[:, :self.idx_1] = func(T_2d[:, :self.idx_1])

        # 2. 活性区 (Active Fuel) -> 轴向切片: [:, idx_1 : idx_2]
        if self.idx_2 > self.idx_1:
            func = getattr(self.fuel_mat, prop_name)
            out_2d[:, self.idx_1:self.idx_2] = func(T_2d[:, self.idx_1:self.idx_2])

        # 3. 上反射层 (Upper Reflector) -> 轴向切片: [:, idx_2 : ny]
        if self.idx_2 < self.ny:
            func = getattr(self.reflector_mat, prop_name)
            out_2d[:, self.idx_2:] = func(T_2d[:, self.idx_2:])

        # 返回扁平化数组，匹配 BaseHeatConduction 中 k_node[:] 的赋值需求
        # 注意：flatten() 默认返回深拷贝(copy)，所以不会破坏 _out_buffer 的内存安全
        return out_2d.flatten()

    def conductivity(self, T: Union[float, np.ndarray]) -> np.ndarray:
        return self._compute_property(T, 'conductivity')

    def heat_capacity(self, T: Union[float, np.ndarray]) -> np.ndarray:
        return self._compute_property(T, 'heat_capacity')

    def density(self, T: Union[float, np.ndarray]) -> np.ndarray:
        return self._compute_property(T, 'density')
