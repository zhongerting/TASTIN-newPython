import numpy as np
from abc import ABC, abstractmethod
from typing import List, Optional, Tuple, Union

# ==========================================
# 1. 物理条件层 (Physical Condition Layer)
# ==========================================

class BaseBoundaryCondition(ABC):
    """
    边界条件抽象基类
    """

    def __init__(self, shape: Tuple[int, ...]):
        """
        :param shape: 边界数据的维度形状.
                      1D 问题: (1,)
                      2D 问题: (Ny,) [左右边界] 或 (Nx,) [上下边界]
        """
        self.shape = shape
        self.bc_type = "base"
        # 附加热阻 (R_add)，用户可用于模拟污垢或接触热阻
        # 初始化为零向量
        self.R_add = np.zeros(shape, dtype=float)

    def set_additional_resistance(self, R_val: Union[float, np.ndarray]):
        """
        设置附加热阻
        :param R_val: 标量(统一热阻) 或 数组(分布热阻)
        """
        if np.isscalar(R_val):
            self.R_add.fill(R_val)
        else:
            if R_val.shape != self.shape:
                raise ValueError(f"R_val shape {R_val.shape} mismatch with BC shape {self.shape}")
            self.R_add[:] = R_val

    @abstractmethod
    def compute_flux_from_node(self, T_node: np.ndarray, R_internal: np.ndarray) -> np.ndarray:
        """
        [向量化计算] 根据内部节点温度和内部热阻计算通过边界的热流
        :param T_node: 边界相邻控制体的中心温度 (Vector)
        :param R_internal: 控制体中心到表面的导热热阻 (Vector)
        :return: 热流 [W] (Vector, 正值代表 T_ext -> T_node)
        """
        pass

    @abstractmethod
    def update_params(self, **kwargs):
        """用于在瞬态过程中更新边界参数 (如随时间变化的 T_ext)"""
        pass

class ResistanceBC(BaseBoundaryCondition):
    """
    通用热阻型边界条件 (Robin / Dirichlet / Neumann 的统一形式)
    Flux = (T_ext - T_node) / (R_internal + R_ext + R_add)
    """

    def __init__(self,
                 T_ext: Union[float, np.ndarray],
                 R_ext: Union[float, np.ndarray],
                 shape: Tuple[int, ...]):
        super().__init__(shape)
        self.bc_type = "resistance"

        # 处理 T_ext (标量 -> 广播, 数组 -> 直接存储)
        if np.isscalar(T_ext):
            self.T_ext = np.full(shape, T_ext, dtype=float)
        else:
            self.T_ext = np.array(T_ext, dtype=float)

        # 处理 R_ext
        if np.isscalar(R_ext):
            self.R_ext = np.full(shape, R_ext, dtype=float)
        else:
            self.R_ext = np.array(R_ext, dtype=float)

    def update_params(self, T_ext: Optional[Union[float, np.ndarray]] = None,
                      R_ext: Optional[Union[float, np.ndarray]] = None):
        if T_ext is not None:
            if np.isscalar(T_ext):
                self.T_ext.fill(T_ext)
            else:
                self.T_ext[:] = T_ext
        if R_ext is not None:
            if np.isscalar(R_ext):
                self.R_ext.fill(R_ext)
            else:
                self.R_ext[:] = R_ext

    def compute_flux_from_node(self, T_node: np.ndarray, R_internal: np.ndarray) -> np.ndarray:
        # 总热阻 = 内部导热 + 外部边界 + 附加(污垢)
        # 利用 NumPy 数组运算直接处理整个边界
        R_total = R_internal + self.R_ext + self.R_add

        # 避免除以零 (处理绝热边界 R -> inf 的情况)
        # 如果 R_total 非常大，Flux 趋近于 0
        with np.errstate(divide='ignore', invalid='ignore'):
            flux = (self.T_ext - T_node) / R_total
            # 将 nan/inf 修正为 0 (通常发生在绝热 R=inf)
            flux = np.nan_to_num(flux, nan=0.0, posinf=0.0, neginf=0.0)

        return flux

class FluxBC(BaseBoundaryCondition):
    """
    固定热流边界条件 (Neumann)
    注意：在热阻网络法中，这通常作为源项处理，但在边界框架下，
    我们可以反推等效表面温度，或者直接强制返回 Flux。
    此处直接返回设定 Flux。
    """
    def __init__(self, q_flux: Union[float, np.ndarray], shape: Tuple[int, ...]):
        super().__init__(shape)
        self.bc_type = "flux"
        if np.isscalar(q_flux):
            self.q_flux = np.full(shape, q_flux, dtype=float)
        else:
            self.q_flux = np.array(q_flux, dtype=float)

    def compute_flux_from_node(self, T_node: np.ndarray, R_internal: np.ndarray) -> np.ndarray:
        return self.q_flux

    def update_params(self, q_flux: Union[float, np.ndarray]):
        if np.isscalar(q_flux):
            self.q_flux.fill(q_flux)
        else:
            self.q_flux[:] = q_flux

class DynamicRadiationResistanceBC(ResistanceBC):
    """
    动态非线性辐射热流边界条件 (Dynamic Non-linear Radiation Flux BC)
    继承自 FluxBC，利用直线法隐式求解器的实时状态注入，在每个残差评估步动态计算 T^4 辐射热流。
    """

    def __init__(self,
                 emissivity: Union[float, np.ndarray],
                 area_array: np.ndarray,
                 T_env: Union[float, np.ndarray],
                 shape: Tuple[int, ...]):
        """
        :param emissivity: 表面发射率
        :param area_array: 实际参与辐射的裸露表面积数组 [m^2]
        :param T_env: 环境辐射背景温度 [K]
        :param shape: 边界形状
        """
        # 初始化父类 FluxBC，初始热流设为 0
        if np.isscalar(T_env):
            t_env_arr = np.full(shape, T_env, dtype=float)
        else:
            t_env_arr = np.array(T_env, dtype=float)

        super().__init__(T_ext=t_env_arr, R_ext=1.0e15, shape=shape)

        # 标量/数组安全处理，兼容性更好
        if np.isscalar(emissivity):
            self.emissivity = np.full(shape, emissivity, dtype=float)
        else:
            self.emissivity = np.array(emissivity, dtype=float)

        self.area_array = np.array(area_array, dtype=float)

        if np.isscalar(T_env):
            self.T_env = np.full(shape, T_env, dtype=float)
        else:
            self.T_env = np.array(T_env, dtype=float)

        self.sigma = 5.670374419e-8
        self.T_ext[:] = self.T_env
        self.h_rad = np.zeros(shape, dtype=float)
        self.G_rad = np.zeros(shape, dtype=float)
        self.T_ref = np.full(shape, 1.0e-3, dtype=float)
        self.q_flux = np.zeros(shape, dtype=float)
        self.last_q_flux = self.q_flux
        self.invalid_temperature_detected = False
        self.min_trial_temperature = np.inf

    def reset_diagnostics(self):
        self.invalid_temperature_detected = False
        self.min_trial_temperature = np.inf

    def update_state(self, T_node: np.ndarray, T_surface: Optional[np.ndarray] = None):
        """
        [关键动态钩子]
        接收最新试探温度，动态更新 self.q_flux。
        由于继承自 FluxBC，BoundaryRegion 在计算时会自动获取更新后的 q_flux。
        """
        if T_surface is None:
            t_ref_raw = np.asarray(T_node, dtype=float)
        else:
            t_surface_arr = np.asarray(T_surface, dtype=float)
            t_node_arr = np.asarray(T_node, dtype=float)
            t_ref_raw = np.where(t_surface_arr > 1.0e-3, t_surface_arr, t_node_arr)

        if t_ref_raw.size == 0:
            return

        raw_min = float(np.min(t_ref_raw))
        if raw_min < 0.0:
            self.invalid_temperature_detected = True
            self.min_trial_temperature = min(self.min_trial_temperature, raw_min)

        T_safe = np.maximum(t_ref_raw, 1.0e-3)
        T_env_safe = np.maximum(self.T_env, 1.0e-3)

        self.T_ref[:] = T_safe
        self.h_rad[:] = (
            self.emissivity
            * self.sigma
            * (T_safe + T_env_safe)
            * (T_safe ** 2 + T_env_safe ** 2)
        )
        self.G_rad[:] = self.h_rad * self.area_array
        self.G_rad[:] = np.nan_to_num(self.G_rad, nan=0.0, posinf=0.0, neginf=0.0)
        self.G_rad[:] = np.maximum(self.G_rad, 1.0e-30)

        with np.errstate(divide='ignore', invalid='ignore'):
            self.R_ext[:] = 1.0 / self.G_rad
        self.R_ext[:] = np.nan_to_num(self.R_ext, nan=1.0e30, posinf=1.0e30, neginf=1.0e30)
        self.T_ext[:] = self.T_env

        self.q_flux[:] = self.G_rad * (self.T_env - T_safe)

    def compute_flux_from_node(self, T_node: np.ndarray, R_internal: np.ndarray) -> np.ndarray:
        """重写父类方法以支持未来可能的直接函数式调用"""
        self.update_state(T_node)
        return super().compute_flux_from_node(T_node, R_internal)

    def update_params(self,
                      T_env: Optional[Union[float, np.ndarray]] = None,
                      emissivity: Optional[Union[float, np.ndarray]] = None,
                      **kwargs):
        """支持瞬态参数更新 (如飞行器进入太空阴影区导致 T_env 突变)"""
        if T_env is not None:
            if np.isscalar(T_env):
                self.T_env.fill(T_env)
            else:
                self.T_env[:] = T_env
            self.T_ext[:] = self.T_env
        if emissivity is not None:
            if np.isscalar(emissivity):
                self.emissivity.fill(emissivity)
            else:
                self.emissivity[:] = emissivity


# Backward-compatible class name. The implementation is now resistance-based.
DynamicRadiationFluxBC = DynamicRadiationResistanceBC


# ==========================================
# 2. 边界区域管理器 (Manager Layer)
# ==========================================
class BoundaryRegion:
    """
    边界区域管理器 (原 BoundaryFace 的升级版)

    2026.1.27添加了边界条件的叠加；
    核心原理：诺顿等效 (Parallel Conductance Superposition)
    G_total = Sum(1/R_i)
    J_total = Sum(T_i/R_i) + Sum(Q_flux_j)
    """

    def __init__(self, shape: Tuple[int, ...], area_array: np.ndarray):
        """
        :param shape: 边界数据的维度形状.
        :param area_array: 边界上每个网格面的面积数组 (用于计算 h -> R)
        """
        self.shape = shape
        # 替代原本的 set_geometry，在初始化时强制绑定
        self.area = area_array

        if self.area.shape != self.shape:
            raise ValueError(f"Area array shape {self.area.shape} mismatch with BC shape {self.shape}")

        # 状态变量
        self.T_surface = np.zeros(shape)
        self.current_flux = np.zeros(shape)     # 每个边界流出的能量[W]（并非面热流密度）
        self.R_internal = np.ones(shape)
        self.T_adj_node = np.zeros(shape)

        self.conditions: List[BaseBoundaryCondition] = []

        # 2026.1.27 新增，叠加计算缓存
        # G_sum: 总外部热导 (W/K for total, or 1/R_unit if normalized)
        # 这里假设 R 是绝对热阻 [K/W] 或 单位热阻 [K m^2/W]。
        # 为了通用性，我们直接对 R 取倒数。
        self.G_sum = np.zeros(shape)

        # J_sum: 总驱动源项 (相当于诺顿等效电流源)
        # 对于电阻型 BC: J = T_ext / R_ext
        # 对于热流型 BC: J = Q_flux
        self.J_sum = np.zeros(shape)
        self.Q_sum_flux = np.zeros(shape)

        # 有效等效参数
        self.R_eff = np.full(shape, np.inf)
        self.T_eff = np.zeros(shape)
        self._g_work = np.zeros(shape)
        self._mix_work = np.zeros(shape)
        self._resistance_mask = np.zeros(shape, dtype=bool)
        # 默认绝热
        self.add_resistance_condition(T_ext=300.0, R_ext=1e15)

    # 2026.1.27 新增方法，清除原有的热导G和热流J
    # 已被废弃，不再调用！
    def clear_boundary_conditions(self):
        """
        [关键方法] 在时间步开始时调用。
        重置累加器，并重新装载静态边界条件 (self.conditions 列表中的条件)。
        这样可以清除上一时间步由 Coupler 动态施加的临时 BC，保留固定的物理 BC。
        """
        self.G_sum.fill(0.0)
        self.J_sum.fill(0.0)

        # 重新应用所有静态添加的条件 (如固定的对流换热、固定热流)
        for bc in self.conditions:
            self._accumulate_bc(bc)

    # 2026.1.27 新增方法，将已有的边界条件叠加到当前状态
    # 已被废弃，不再调用！
    def _accumulate_bc(self, bc: BaseBoundaryCondition):
        """内部辅助：将一个 BC 对象叠加到当前状态"""
        if isinstance(bc, ResistanceBC):
            # R_total_ext = R_ext + R_add (不包含 R_internal，因为那是 Solver 侧的事)
            R_ext_total = bc.R_ext + bc.R_add

            with np.errstate(divide='ignore'):
                G = 1.0 / R_ext_total

            # 处理绝热情况 (G=0)
            G = np.nan_to_num(G, posinf=0.0, neginf=0.0)

            # 累加
            self.G_sum += G
            self.J_sum += bc.T_ext * G

        elif isinstance(bc, FluxBC):
            # 固定热流直接作为源项叠加
            self.J_sum += bc.q_flux

    # 2026.1.27 新增方法，叠加动态边界条件
    # 已被废弃，不再调用！
    def update_params(self, T_ext: Union[float, np.ndarray], R_ext: Union[float, np.ndarray]):
        """
        [新增] 供外部 Coupler 调用，叠加动态边界条件。
        此方法不会覆盖旧值，而是执行并联叠加。

        :param T_ext: 外部温度
        :param R_ext: 外部热阻 (不包含固体内部热阻)
        """
        # 转换输入
        if np.isscalar(T_ext):
            T = np.full(self.shape, T_ext)
        else:
            T = np.asarray(T_ext)

        if np.isscalar(R_ext):
            R = np.full(self.shape, R_ext)
        else:
            R = np.asarray(R_ext)

        # 计算电导 G = 1/R
        with np.errstate(divide='ignore'):
            G_new = 1.0 / R

        # 处理 R=0 (Dirichlet) 或 R=inf (绝热)
        # R -> 0, G -> inf: 数值上用极大值代替，或者在 Solver 中特殊处理
        # 这里做安全截断防止溢出破坏叠加
        G_new = np.nan_to_num(G_new, posinf=1e20)

        # 叠加到总状态
        self.G_sum += G_new
        self.J_sum += T * G_new

    # --- [辅助方法] 为了兼容性和灵活性 ---

    def update_geometry(self, new_area_array: np.ndarray):
        """
        (原 set_geometry 的替代)
        如果网格变形或重新划分，可调用此方法更新面积
        """
        if new_area_array.shape != self.shape:
            raise ValueError("New geometry shape mismatch")
        self.area[:] = new_area_array

    def get_coupling_snapshot(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        (已恢复) 获取用于耦合的快照状态
        :return: (T_adj_node, R_total_internal)
                 这里的 R_total_internal 指的是从节点中心到表面的热阻
        """
        # return self.T_adj_node, self.R_internal
        return self.T_adj_node, self.R_internal

    def get_coupling_surface_snapshot(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        (已恢复) 获取用于耦合的快照状态
        :return: (T_surface, R_total_internal)
                 这里的 R_total_internal 指的是从节点中心到表面的热阻
        """
        # return self.T_surface, self.R_internal
        return self.T_surface, self.R_internal

    # --- [核心逻辑] ---

    def update_internal_state(self, T_node: np.ndarray, R_int: np.ndarray, current_time: float = None):
        """更新边界"看到"的内部状态"""
        self.T_adj_node[:] = T_node
        self.R_internal[:] = R_int

        # ====== [新增] 动态唤醒逻辑 ======
        # 遍历所有条件，如果是动态边界，则实时更新其内部热流
        for cond in self.conditions:
            if hasattr(cond, 'update_state'):
                if current_time is not None and hasattr(cond, 'current_time'):
                    cond.current_time = current_time
                if isinstance(cond, DynamicRadiationResistanceBC):
                    T_surface_ref = np.where(self.T_surface > 1.0e-3, self.T_surface, self.T_adj_node)
                    cond.update_state(T_node, T_surface=T_surface_ref)
                else:
                    cond.update_state(T_node)

    # 已经经过验证
    # 2026.1.27 修改，增加了多边界条件同时存在的处理逻辑；
    # 已经过测试
    def compute_net_flux_for_solver(self) -> np.ndarray:
        """
        计算叠加后的净热流，并更新表面温度。
        逻辑：
        1. 遍历 ResistanceBC，累加电导 G_sum 和 诺顿源 J_sum_R。
        2. 遍历 FluxBC，累加纯热流 Q_sum_flux。
        3. 若存在导热通道(ResistanceBC)，执行戴维南等效叠加 (T_eff = T_eq_R + Q * R_eq)。
        4. 若无导热通道(纯FluxBC)，直接注入热流。
        """
        if not self.conditions:
            # 无任何条件 -> 绝热
            self.current_flux.fill(0.0)
            self.T_surface[:] = self.T_adj_node[:]
            return self.current_flux

        # --- 1. 初始化累加器 ---
        self.G_sum.fill(0.0)
        self.J_sum.fill(0.0)
        self.Q_sum_flux.fill(0.0)

        # --- 2. 循环遍历并分类累加 ---
        for bc in self.conditions:
            if bc.bc_type == "resistance":
                # 计算总外部热阻 (不含 R_internal)
                np.add(bc.R_ext, bc.R_add, out=self._mix_work)

                # G = 1/R
                with np.errstate(divide='ignore', invalid='ignore'):
                    np.divide(1.0, self._mix_work, out=self._g_work)
                    # 绝热保护 (R=inf -> G=0)
                    np.nan_to_num(self._g_work, copy=False, posinf=0.0)
                    # Dirichlet保护 (R=0 -> G=inf, 截断防止溢出)
                    np.nan_to_num(self._g_work, copy=False, posinf=1e20)

                # 累加
                self.G_sum += self._g_work
                np.multiply(bc.T_ext, self._g_work, out=self._mix_work)
                self.J_sum += self._mix_work

            elif bc.bc_type == "flux":
                # 纯热流直接累加 (单位: W)
                self.Q_sum_flux += bc.q_flux

        # --- 3. 混合逻辑计算 ---
        # 使用掩码区分 "有电阻回路的区域" 和 "纯热流/绝热区域"
        # 阈值 1e-20 用于判定是否存在有效的 ResistanceBC
        np.greater(self.G_sum, 1e-20, out=self._resistance_mask)

        # 预填充纯热流结果 (应对没有 ResistanceBC 的情况)
        # 此时 Flux = Q_flux (全部注入，无处分流)
        self.current_flux[:] = self.Q_sum_flux
        self.R_eff.fill(np.inf)
        self.T_eff.fill(0.0)

        # --- 针对有 ResistanceBC 的区域进行戴维南等效计算 ---
        if np.any(self._resistance_mask):
            np.divide(1.0, self.G_sum, out=self.R_eff, where=self._resistance_mask)
            np.divide(self.J_sum, self.G_sum, out=self.T_eff, where=self._resistance_mask)
            np.multiply(self.Q_sum_flux, self.R_eff, out=self._mix_work)
            np.add(self.T_eff, self._mix_work, out=self.T_eff, where=self._resistance_mask)
            np.add(self.R_eff, self.R_internal, out=self._mix_work)
            np.subtract(self.T_eff, self.T_adj_node, out=self._g_work)
            np.divide(self._g_work, self._mix_work, out=self.current_flux, where=self._resistance_mask)

        # --- 4. 更新表面温度 ---
        # T_surf = T_node + Flux * R_internal
        # (注: current_flux 为流向节点方向? 不，compute_flux_from_node 定义通常是 (Text-Tnode)/R，即流入为正)
        # 若 current_flux 为流入节点，则节点得到热量，
        # 但这里的 R_internal 压降意味着: T_surf = T_node + Flux_in * R_int
        # (热从高温流向低温: 如果 Flux>0 (流入), 则 T_surf > T_node)
        np.multiply(self.current_flux, self.R_internal, out=self.T_surface)
        np.add(self.T_adj_node, self.T_surface, out=self.T_surface)

        for bc in self.conditions:
            if isinstance(bc, DynamicRadiationResistanceBC):
                bc.q_flux[:] = bc.G_rad * (bc.T_env - self.T_surface)

        return self.current_flux

    # 2026.1.27 修改的方法，同时处理多个边界条件（已禁用）
    # def compute_net_flux_for_solver(self) -> np.ndarray:
    #     """
    #     计算进入求解器的净热流 (基于叠加后的有效参数)
    #     """
    #     # 1. 计算有效参数 (Effective Parameters)
    #     # R_eff = 1 / G_sum
    #     with np.errstate(divide='ignore', invalid='ignore'):
    #         self.R_eff = 1.0 / self.G_sum
    #         # T_eff = J_sum / G_sum = J_sum * R_eff
    #         self.T_eff = self.J_sum * self.R_eff
    #
    #     # 清洗数据：如果 G_sum 为 0 (全绝热)，R_eff -> inf, T_eff -> NaN
    #     # 此时 Flux 应该为 0
    #     self.R_eff = np.nan_to_num(self.R_eff, posinf=1e20, nan=1e20)
    #     self.T_eff = np.nan_to_num(self.T_eff, nan=0.0)
    #
    #     # 2. 计算与内部节点的交换热流
    #     # Flux = (T_eff - T_node) / (R_eff + R_internal)
    #     # 注意：这里我们重新引入了 R_internal，因为它是串联在这一侧的
    #     R_total_path = self.R_eff + self.R_internal
    #
    #     with np.errstate(divide='ignore', invalid='ignore'):
    #         self.current_flux = (self.T_eff - self.T_adj_node) / R_total_path
    #         self.current_flux = np.nan_to_num(self.current_flux, nan=0.0)
    #
    #     # 3. 更新表面温度
    #     # T_surf = T_node + Flux * R_internal
    #     self.T_surface = self.T_adj_node + self.current_flux * self.R_internal
    #
    #     return self.current_flux

    # --- 工厂方法 (保持不变) ---
    def add_resistance_condition(self, T_ext, R_ext, R_add=0.0):
        bc = ResistanceBC(T_ext, R_ext, self.shape)
        bc.set_additional_resistance(R_add)
        self.conditions.append(bc)
        # 2026.1.27 新增，立即叠加一次，以免初始化时状态为空
        self._accumulate_bc(bc)
        return bc

    def add_convection_condition(self, T_fluid, h_coeff):
        with np.errstate(divide='ignore'):
            R_conv = 1.0 / (h_coeff * self.area)
            R_conv = np.nan_to_num(R_conv, nan=1e15, posinf=1e15)
        bc = ResistanceBC(T_fluid, R_conv, self.shape)
        self.conditions.append(bc)
        # 2026.1.27 新增，立即叠加一次，以免初始化时状态为空
        self._accumulate_bc(bc)
        return bc

    def add_flux_condition(self, q_flux):
        bc = FluxBC(q_flux, self.shape)
        self.conditions.append(bc)
        # 2026.1.27 新增，立即叠加一次，以免初始化时状态为空
        self._accumulate_bc(bc)
        return bc

    def add_dynamic_radiation_condition(self, emissivity, bare_area_array, T_env=3.0):
        """挂载动态辐射边界"""
        bc = DynamicRadiationResistanceBC(emissivity, bare_area_array, T_env, self.shape)
        self.conditions.append(bc)
        # 立即叠加一次，防止初始化阶段状态为空
        self._accumulate_bc(bc)
        return bc

    def clear_conditions(self):
        self.conditions.clear()
        self.clear_boundary_conditions()
