import numpy as np
from typing import Literal, Union, Optional, Callable, TYPE_CHECKING
from profiler import TEASAProfiler

if TYPE_CHECKING:
    from Solvers.Hydrodynamics.Components import FluidChannel
    from Solvers.HeatConduction.Boundary import BoundaryRegion, ResistanceBC


class SolidSolidCouple1D:
    """
    一维固体-固体导热耦合器
    适配基于 BoundaryRegion 的 HeatConduction1D
    """

    def __init__(self,
                 obj1,  # Type: HeatConduction1D
                 obj2,  # Type: HeatConduction1D
                 direction: Literal['right', 'left', 'outer', 'inner'],
                 contact_resistance: float = 0.0):
        """
        初始化耦合器
        :param obj1: 主物体
        :param obj2: 相对物体
        :param direction: obj2 在 obj1 的哪个方向?
                          'right'/'outer': obj1(Outer) <-> obj2(Inner)
                          'left'/'inner':  obj1(Inner) <-> obj2(Outer)
        :param contact_resistance: 接触热阻 [K/W] (总热阻，非单位面积)
        """
        self.obj1 = obj1
        self.obj2 = obj2
        self.R_contact = contact_resistance

        # --- 1. 解析方向并确定接口位置 ---
        # 1D 求解器统一使用 'inner' (左/小半径) 和 'outer' (右/大半径) 作为键
        direction = direction.lower()

        if direction in ['right', 'outer']:
            # Obj1 右侧(Outer) 连接 Obj2 左侧(Inner)
            self.loc1 = 'outer'
            self.loc2 = 'inner'
        elif direction in ['left', 'inner']:
            # Obj1 左侧(Inner) 连接 Obj2 右侧(Outer)
            self.loc1 = 'inner'
            self.loc2 = 'outer'
        else:
            raise ValueError(f"Unknown direction: {direction}")

        # --- 2. 获取现有的边界区域对象 ---
        self.bound1 = self.obj1.boundaries[self.loc1]
        self.bound2 = self.obj2.boundaries[self.loc2]

        if self.bound1 is None or self.bound2 is None:
            raise ValueError("HeatConduction1D boundary not initialized (is None).")

        # --- 3. 挂载热阻边界条件 (ResistanceBC) ---
        # 初始状态：设为环境温度 300K，R_ext=0 (稍后 sync 会覆盖)
        # 将接触热阻 R_contact 平分到两侧 (0.5 * R)，也可全加在一侧

        self.bc1: ResistanceBC = self.bound1.add_resistance_condition(
            T_ext=300.0, R_ext=0.0, R_add=self.R_contact
        )

        self.bc2: ResistanceBC = self.bound2.add_resistance_condition(
            T_ext=300.0, R_ext=0.0, R_add=self.R_contact
        )

    def sync(self):
        """
        【核心方法】同步状态
        在每个时间步计算 Flux 之前调用 (通常在 get_derivatives 内部或外部调度器中)
        """
        # 1. 获取快照 (T_surface, R_total_internal)
        # 注意：对于 1D，返回的是 shape=(1,) 的数组
        T1_surf, R1_int = self.bound1.get_coupling_snapshot()
        T2_surf, R2_int = self.bound2.get_coupling_snapshot()

        # 2. 交叉更新参数
        # Obj1 的边界条件：外部温度是 Obj2 的表面温度，外部热阻是 Obj2 的内部热阻
        self.bc1.update_params(T_ext=T2_surf, R_ext=R2_int)

        # Obj2 的边界条件：外部温度是 Obj1 的表面温度，外部热阻是 Obj1 的内部热阻
        self.bc2.update_params(T_ext=T1_surf, R_ext=R1_int)


class SolidSolidCouple2D:
    """
    二维固体-固体导热耦合器
    连接两个 HeatConduction2D 物体的边界 (Line-to-Line Coupling)
    """

    def __init__(self,
                 obj1,  # Type: HeatConduction2D
                 obj2,  # Type: HeatConduction2D
                 direction: Literal['right', 'left', 'top', 'bottom'],
                 contact_resistance: float = 0.0):
        """
        :param obj1: 主物体
        :param obj2: 邻居物体
        :param direction: obj2 在 obj1 的哪个方向?
               (例如 'right' 意味着 obj1.right 与 obj2.left 耦合)
        :param contact_resistance: 接触热阻 [K/W] (注意：如果是接触热阻系数 m^2K/W，需除以面积)
                                   此处简化为直接接受总热阻或单位面积热阻，
                                   视 BoundaryRegion.add_resistance_condition 的 R_add 定义而定。
                                   通常 ResistanceBC 里的 R 是绝对热阻 [K/W]。
        """
        self.obj1 = obj1
        self.obj2 = obj2
        self.R_contact = contact_resistance

        # --- 1. 确定接口映射关系 ---
        # 键: obj1 的方向
        # 值: (obj1 的边界键名, obj2 的边界键名)
        self.dir_map = {
            'right': ('right', 'left'),  # Obj1 右连 Obj2 左
            'left': ('left', 'right'),  # Obj1 左连 Obj2 右
            'top': ('top', 'bottom'),  # Obj1 上连 Obj2 下
            'bottom': ('bottom', 'top')  # Obj1 下连 Obj2 上
        }

        direction = direction.lower()
        if direction not in self.dir_map:
            raise ValueError(f"Unknown direction for 2D couple: {direction}")

        self.loc1, self.loc2 = self.dir_map[direction]

        # --- 2. 获取边界区域对象 ---
        self.bound1 = self.obj1.boundaries[self.loc1]
        self.bound2 = self.obj2.boundaries[self.loc2]

        # --- 3. 校验维度匹配 (关键) ---
        # 必须保证连接的两个边界拥有相同数量的网格节点
        if self.bound1.shape != self.bound2.shape:
            raise ValueError(
                f"Boundary shape mismatch in 2D coupling: "
                f"Obj1({self.loc1}) shape {self.bound1.shape} != "
                f"Obj2({self.loc2}) shape {self.bound2.shape}"
            )

        # --- 4. 挂载热阻边界条件 ---
        # 支持向量化: 如果 contact_resistance 是标量，会自动广播

        self.bc1: ResistanceBC = self.bound1.add_resistance_condition(
            T_ext=300.0, R_ext=0.0, R_add=self.R_contact
        )

        self.bc2: ResistanceBC = self.bound2.add_resistance_condition(
            T_ext=300.0, R_ext=0.0, R_add=self.R_contact
        )

    def sync(self):
        """
        同步状态 (向量化)
        无需循环，直接交换整个数组
        """
        # 1. 获取快照 (T_surf, R_int) -> 这些都是 NumPy 数组
        T1_surf, R1_int = self.bound1.get_coupling_snapshot()
        T2_surf, R2_int = self.bound2.get_coupling_snapshot()

        # 2. 交叉更新
        # 数组直接传递，ResistanceBC 内部会自动处理
        self.bc1.update_params(T_ext=T2_surf, R_ext=R2_int)
        self.bc2.update_params(T_ext=T1_surf, R_ext=R1_int)


class FluidSolidCouple:
    """
    流体-固体热耦合器 (Fluid-Solid Thermal Coupler)

    物理功能:
    1. 计算对流换热: Nu -> h
    2. 更新固体边界: 利用 ResistanceBC 设置 Robin 边界 (T_ext=T_fluid, R_ext=1/hA)
    3. 更新流体源项: 利用 add_coupling_source_distribution 注入半隐式源项 (Q = hA*Tw - hA*Tf)

    设计原则:
    - 严格调用 Components.py 和 Boundary.py 的现有 API。
    - 预留 Nu 计算接口 (correlation_func)。
    """

    def __init__(self,
                 name: str,
                 fluid: 'FluidChannel',
                 solid_boundary_region: 'BoundaryRegion',
                 heated_perimeter: float,
                 correlation_func: Callable[[float, float, float], float],
                 solid_node_capacitance: Optional[np.ndarray] = None):
        """
        初始化耦合器

        :param name: 耦合器名称
        :param fluid: 流体通道对象 (必须包含 temperature_vector, add_coupling_source_distribution 等)
        :param solid_boundary_region: 固体边界区域对象 (必须包含 T_surface, add_resistance_condition)
        :param heated_perimeter: 湿周 [m]
        :param correlation_func: Nu数计算函数接口, 签名需为: func(Re, Pr, P_D_ratio) -> Nu
        """
        self.name = name
        self.fluid = fluid
        self.solid_bound = solid_boundary_region
        self.perimeter = heated_perimeter
        self.correlation_func = correlation_func

        # [新增] 固体热容缓存 (用于稳定性计算)
        self.solid_node_capacitance = solid_node_capacitance

        # [新增] 缓存上一次计算的耦合系数 (Lambda = h * A)
        self._last_lambda = None
        self._interface_relaxation_previous = None
        self._last_coupling_diagnostics = None

        # --- 1. 校验网格一致性 ---
        # 确保流体节点数与固体边界的节点数一致
        # BoundaryRegion.shape 是一个元组，通常 (N,) 表示节点数
        if self.solid_bound.shape[0] != self.fluid.n_nodes:
            raise ValueError(f"[{self.name}] Mesh mismatch: "
                             f"Fluid({self.fluid.n_nodes}) vs Solid({self.solid_bound.shape[0]})")

        # 如果提供了固体热容，也要校验形状
        if self.solid_node_capacitance is not None:
            if self.solid_node_capacitance.shape[0] != self.fluid.n_nodes:
                raise ValueError(
                    f"[{self.name}] Solid capacitance shape mismatch: {self.solid_node_capacitance.shape}")

        # --- 2. 几何参数缓存 ---
        # 计算每个控制体的换热面积 A_node [m^2]
        # 假设流体节点长度均匀，使用 fluid.node_length
        self.node_areas = self.perimeter * self.fluid.node_length

        # --- 3. 初始化固体边界条件 ---
        # 利用 BoundaryRegion.add_resistance_condition 创建 Robin 边界
        # 初始状态：T_ext = T_fluid, R_ext = 总热阻 [K/W]

        T_fluid_init = self.fluid.temperature_vector
        # 给定一个非零初始热阻防止除零错误 (例如 hA=1000 W/K 时 R=1e-3)
        R_init = np.ones(self.fluid.n_nodes) * 1.0e-3

        # [关键调用] 引用现有的 boundary region 方法创建 BC 对象
        # 我们持有 self.solid_bc 引用，以便在 execute 中更新它
        self.solid_bc = self.solid_bound.add_resistance_condition(
            T_ext=T_fluid_init,
            R_ext=R_init
        )

    def set_solid_capacitance(self, capacitance_arr: np.ndarray):
        """[辅助] 更新/设置固体节点热容"""
        if capacitance_arr.shape[0] != self.fluid.n_nodes:
            raise ValueError(f"Capacitance shape mismatch: {capacitance_arr.shape}")
        self.solid_node_capacitance = capacitance_arr

    def reset_interface_relaxation(self):
        self._interface_relaxation_previous = None
        self._last_coupling_diagnostics = None

    @staticmethod
    def _interface_state_shapes_match(previous_state, current_state):
        if previous_state is None:
            return False
        for key, value in current_state.items():
            if key not in previous_state:
                return False
            if np.asarray(previous_state[key]).shape != np.asarray(value).shape:
                return False
        return True

    @staticmethod
    def _max_abs_delta_or_none(current, previous):
        if previous is None:
            return None
        current_arr = np.asarray(current, dtype=float)
        previous_arr = np.asarray(previous, dtype=float)
        if current_arr.shape != previous_arr.shape:
            return None
        if current_arr.size == 0:
            return 0.0
        return float(np.max(np.abs(current_arr - previous_arr)))

    @staticmethod
    def _max_existing(*values):
        finite_values = [float(value) for value in values if value is not None]
        if not finite_values:
            return None
        return max(finite_values)

    def get_coupling_diagnostics(self):
        if self._last_coupling_diagnostics is None:
            return None
        return dict(self._last_coupling_diagnostics)

    @TEASAProfiler.profile
    def execute(self, interface_relaxation: float = 1.0):
        """
        执行耦合计算步骤:
        1. 获取状态 (Fluid T, Solid T, Flow properties)
        2. 计算 Nu 和 h
        3. 更新固体 BoundaryCondition (R_ext, T_ext)
        4. 更新流体 CouplingSource (Explicit, Implicit)
        """
        # === A. 获取流体物性与状态 (Vectorized) ===
        # 使用 FluidChannel 的向量化属性
        if not (0.0 < interface_relaxation <= 1.0):
            raise ValueError("interface_relaxation must be in (0, 1].")

        T_f = np.asarray(self.fluid.temperature_vector, dtype=float)
        P_f = self.fluid.pressure_vector
        rho = self.fluid.density_vector
        vel = self.fluid.velocity_vector  # 特征流速

        # 获取物性 (调用 Material 方法，需传入数组)
        # 注意：Components.py 中 FluidVolume 使用 material.viscosity(T, P)
        mu = self.fluid.material.viscosity(T_f, P_f)
        k_f = self.fluid.material.conductivity(T_f, P_f)
        Cp_f = self.fluid.material.heat_capacity(T_f, P_f)

        # 防止粘度为 0 (除零保护)
        mu = np.maximum(mu, 1e-10)

        # === B. 计算无量纲数 ===
        # Re = (rho * v * D_h) / mu
        Re = (rho * np.abs(vel) * self.fluid.d_h) / mu

        # Pr = mu * Cp / k
        Pr = self.fluid.material.prandtl_number(T_f, P_f)

        # === C. 计算换热系数 (Correlation) ===
        # 调用预留的接口计算 Nu
        # 假设 P/D ratio 为 1.1 (或后续从 self.fluid.geom_params 获取)
        Nu = self.correlation_func(Re, Pr, 1.1)

        # h = Nu * k / D_h [W/(m^2 K)]
        h_vals = Nu * k_f / self.fluid.d_h

        # === D. 计算耦合系数 (Lambda) ===
        # Lambda = h * A [W/K]
        # 这是半隐式耦合的核心系数
        lambda_vals = h_vals * self.node_areas
        safe_lambda = np.maximum(lambda_vals, 1e-8)

        # === E. 更新固体边界 (Robin BC) ===
        # 固体方程看见的边界: Flux = (T_ext - T_surf) / R_ext
        # 物理对流公式: Flux = h*A * (T_fluid - T_surf)
        # 对比可知: T_ext = T_fluid, R_ext = 1 / (h*A) = 1 / lambda

        R_convection = 1.0 / safe_lambda

        # 直接更新 ResistanceBC 对象的属性 (数组切片更新)
        self.solid_bc.R_ext[:] = R_convection
        self.solid_bc.T_ext[:] = T_f

        # === F. 更新流体源项 (Semi-Implicit Source) ===
        # 流体方程看见的源项: Q = Explicit - Implicit * T_fluid
        # 物理对流公式: Q = h*A * (T_wall - T_fluid)
        #                 = (h*A * T_wall) - (h*A) * T_fluid
        # 对比可知:
        #   Explicit_Part = lambda * T_wall
        #   Implicit_Factor = lambda

        # 更新一下固体壁面温度
        if hasattr(self.solid_bound, "compute_net_flux_for_solver"):
            self.solid_bound.compute_net_flux_for_solver()

        # [关键数据获取] 从 BoundaryRegion 获取当前壁面温度
        T_wall = np.asarray(self.solid_bound.T_surface, dtype=float)

        current_state = {
            "lambda": safe_lambda,
            "T_f": T_f,
            "T_wall": T_wall,
        }
        previous_state = self._interface_relaxation_previous
        previous_shapes_match = self._interface_state_shapes_match(previous_state, current_state)

        if (
            interface_relaxation < 1.0
            and previous_shapes_match
        ):
            relaxed_state = {
                key: (
                    interface_relaxation * value
                    + (1.0 - interface_relaxation) * previous_state[key]
                )
                for key, value in current_state.items()
            }
        else:
            relaxed_state = current_state

        relaxed_lambda = np.maximum(relaxed_state["lambda"], 1e-8)
        relaxed_T_f = relaxed_state["T_f"]
        relaxed_T_wall = relaxed_state["T_wall"]

        self.solid_bc.R_ext[:] = 1.0 / relaxed_lambda
        self.solid_bc.T_ext[:] = relaxed_T_f

        explicit_arr = relaxed_lambda * relaxed_T_wall
        implicit_arr = relaxed_lambda

        # [关键调用] 使用 FluidChannel 提供的半隐式接口
        self.fluid.add_coupling_source_distribution(explicit_arr, implicit_arr)

        if previous_shapes_match:
            previous_explicit = previous_state["lambda"] * previous_state["T_wall"]
            previous_implicit = previous_state["lambda"]
            delta_lambda = self._max_abs_delta_or_none(relaxed_lambda, previous_state["lambda"])
            delta_T_f = self._max_abs_delta_or_none(relaxed_T_f, previous_state["T_f"])
            delta_T_wall = self._max_abs_delta_or_none(relaxed_T_wall, previous_state["T_wall"])
            delta_explicit = self._max_abs_delta_or_none(explicit_arr, previous_explicit)
            delta_implicit = self._max_abs_delta_or_none(implicit_arr, previous_implicit)
        else:
            delta_lambda = None
            delta_T_f = None
            delta_T_wall = None
            delta_explicit = None
            delta_implicit = None

        self._last_coupling_diagnostics = {
            "name": self.name,
            "type": type(self).__name__,
            "interface_relaxation": float(interface_relaxation),
            "relaxed": bool(interface_relaxation < 1.0 and previous_shapes_match),
            "has_previous_state": bool(previous_shapes_match),
            "interface_residual": self._max_existing(delta_lambda, delta_T_f, delta_T_wall),
            "source_residual": self._max_existing(delta_explicit, delta_implicit),
            "delta_lambda": delta_lambda,
            "delta_T_f": delta_T_f,
            "delta_T_wall": delta_T_wall,
            "delta_explicit": delta_explicit,
            "delta_implicit": delta_implicit,
            "max_lambda": float(np.max(relaxed_lambda)) if relaxed_lambda.size else 0.0,
            "max_explicit": float(np.max(np.abs(explicit_arr))) if explicit_arr.size else 0.0,
        }

        # 缓存当前时间步的耦合系数，用于稳定性分析和下次迭代的松弛计算
        self._last_lambda = relaxed_lambda.copy()
        self._interface_relaxation_previous = {
            "lambda": relaxed_lambda.copy(),
            "T_f": np.asarray(relaxed_T_f, dtype=float).copy(),
            "T_wall": np.asarray(relaxed_T_wall, dtype=float).copy(),
        }

    # =========================================================================
    # Phase 6: Time Step Control (Coupling Stability)
    # =========================================================================

    def get_max_stable_dt(self, safety_factor: float = 0.8, max_limit: float = 1.0) -> float:
        """
        [自适应步长] 基于显式流固耦合稳定性条件计算最大允许时间步长
        Criterion: dt < C_solid / (h * A) = C_solid / lambda

        :param safety_factor: 安全系数 (建议 0.5 - 0.8)
        :param max_limit: 最大上限 [s]
        :return: 建议的时间步长
        """
        # 如果未提供固体热容信息，无法计算限制，返回上限
        if self.solid_node_capacitance is None:
            return max_limit

        # 如果尚未执行过 execute，没有 lambda 数据，返回上限
        if self._last_lambda is None:
            return max_limit

        # 分别处理固体热容规定的时间步长和流体热容规定的时间步长
        # 1. 计算固体热容
        C_solid = self.solid_node_capacitance
        # 2. 计算流体热容
        T_f = self.fluid.temperature_vector
        P_f = self.fluid.pressure_vector
        rho = self.fluid.density_vector
        # 调用物性库计算当前温度下的流体比热容
        Cp_f = self.fluid.material.heat_capacity(T_f, P_f)
        # 提取流体节点的流通面积 (兼容不同的属性命名，默认为 1.0 防止报错)
        flow_area = getattr(self.fluid, 'area', getattr(self.fluid, 'flow_area', 1.0))
        V_fluid = flow_area * self.fluid.node_length
        C_fluid = rho * V_fluid * Cp_f

        # 防止热容极小导致除零或溢出
        denom = np.maximum(self._last_lambda, 1e-10)
        C_solid_safe = np.maximum(C_solid, 1e-10)
        C_fluid_safe = np.maximum(C_fluid, 1e-10)

        # 方法1：分别计算最小时间步，随后选取较小步长 (不适用)
        # dt_fluid = C_fluid_safe / denom
        # dt_solid = C_solid_safe / denom
        # dt_vec = min(dt_fluid, dt_solid)

        # 方法2：计算调和平均值
        C_eff = (C_solid_safe * C_fluid_safe) / (C_solid_safe + C_fluid_safe)
        dt_vec = C_eff / denom

        # 计算局部稳定性限制
        # dt_local = C_node / (h * A)
        # 防止除零
        # denom = np.maximum(self._last_lambda, 1e-10)
        # dt_vec = self.solid_node_capacitance / denom

        # 取全场最小值
        dt_stable = np.min(dt_vec) * safety_factor

        return min(dt_stable, max_limit)


class GapCouple2D(SolidSolidCouple2D):
    """
    二维间隙耦合器 (Gap Coupler)
    同时考虑：
    1. 表面对表面辐射 (Surface-to-Surface Radiation)
    2. 间隙气体导热 (Gap Gas Conduction)

    二者为并联关系。
    """

    def __init__(self,
                 obj1, obj2,
                 direction: Literal['right', 'left', 'top', 'bottom'],
                 gap_width: float,  # [新增] 间隙宽度 delta [m]
                 gas_conductivity: float,  # [新增] 气体导热系数 k_gas [W/mK]
                 emissivity1: float = 0.8,
                 emissivity2: float = 0.8):
        super().__init__(obj1, obj2, direction, contact_resistance=0.0)
        self.eps1 = emissivity1
        self.eps2 = emissivity2
        self.sigma = 5.670374419e-8

        # 物理参数
        self.gap = gap_width
        self.k_gas = gas_conductivity

        # 计算参数(总间隙等效热阻)
        self.R_gap_total = 0.

    @TEASAProfiler.profile
    def sync(self):
        # 1. 获取表面温度 (用于辐射)
        T1_surf, _ = self.bound1.get_coupling_surface_snapshot()
        T2_surf, _ = self.bound2.get_coupling_surface_snapshot()

        # 2. 获取面积
        A1 = np.maximum(self.bound1.area, 1e-12)
        A2 = np.maximum(self.bound2.area, 1e-12)

        # 3. 几何判定 (Inner/Outer)
        is_1_inner = A1 < A2
        A_in = np.where(is_1_inner, A1, A2)
        A_out = np.where(is_1_inner, A2, A1)
        eps_in = np.where(is_1_inner, self.eps1, self.eps2)
        eps_out = np.where(is_1_inner, self.eps2, self.eps1)

        # ==========================================================
        # 核心逻辑：计算并联电导 (Conductance = 1/Resistance)
        # ==========================================================

        # --- A. 辐射电导 (G_rad = 1/R_rad) ---
        # h_rad_star = sigma * (T1^2+T2^2)(T1+T2)
        T_sum = T1_surf + T2_surf
        h_rad_star = self.sigma * (T1_surf ** 2 + T2_surf ** 2) * T_sum
        h_rad_star = np.maximum(h_rad_star, 1e-20)

        with np.errstate(divide='ignore', invalid='ignore'):
            # 计算 1/eps，如果 eps=0 会得到 inf，这是我们要的（电阻无穷大）
            inv_eps_in = 1.0 / eps_in
            inv_eps_out = 1.0 / eps_out

            # 计算辐射分母
            # denom = (1/e1) + (1/e2 - 1) * (A1/A2)
            term2 = (inv_eps_out - 1.0) * (A_in / A_out)
            denom_rad = inv_eps_in + term2

            # G_rad = (A * h) / denom
            # 如果 denom 是 inf (因为 eps=0)，则 G_rad 自动为 0
            G_rad = (A_in * h_rad_star) / denom_rad

            # 清洗 NaN (0/0情况) 和 Inf，确保数值安全
            G_rad = np.nan_to_num(G_rad, nan=0.0, posinf=0.0)

        # # Denom for radiation
        # denom_rad = (1.0 / eps_in) + ((1.0 - eps_out) / eps_out) * (A_in / A_out)
        #
        # # G_rad [W/K] = (A_in * h_star) / D
        # G_rad = (A_in * h_rad_star) / denom_rad

        # --- B. 气体导热电导 (G_cond = 1/R_cond) ---
        # 假设主要通过较小的面积 A_in 进行一维导热 (同心圆柱近似)
        # G_cond [W/K] = k_gas * A_in / delta
        # 如果需要更精确的圆柱导热公式: 2*pi*k*L / ln(r2/r1)，但薄间隙下 A/delta 足够精确
        G_cond = (self.k_gas * A_in) / self.gap

        # --- C. 总等效热阻 ---
        # R_total = 1 / (G_rad + G_cond)
        G_total = G_rad + G_cond
        R_gap_total = 1.0 / np.maximum(G_total, 1e-20)

        self.R_gap_total = R_gap_total

        # ==========================================================
        # 更新边界 (串联物理热阻)
        # ==========================================================
        # 获取内部节点状态
        T1_node, R1_int = self.bound1.get_coupling_snapshot()
        T2_node, R2_int = self.bound2.get_coupling_snapshot()

        # Obj1 边界: R_ext = R_gap_total + R_int_2
        self.bc1.update_params(T_ext=T2_node, R_ext=R_gap_total + R2_int)

        # Obj2 边界: R_ext = R_gap_total + R_int_1
        self.bc2.update_params(T_ext=T1_node, R_ext=R_gap_total + R1_int)


class ActiveGapCouple2D(GapCouple2D):
    """
    带源项的二维间隙耦合器 (Active Gap Coupler)

    继承自 GapCouple2D，保留其所有功能（辐射+气体导热）。
    额外增加：随时间变化的固定热流源（如电子冷却/加热）。

    物理实现:
    利用戴维南等效 (Thevenin Equivalent) 将并联热流源转化为串联温差源。
    这避免了当温差为0时计算等效热阻产生的数值奇异点。
    """

    def __init__(self,
                 obj1, obj2,
                 direction: Literal['right', 'left', 'top', 'bottom'],
                 gap_width: float,
                 gas_conductivity: float,
                 emissivity1: float = 0.8,
                 emissivity2: float = 0.8):
        super().__init__(obj1, obj2, direction, gap_width, gas_conductivity,
                         emissivity1, emissivity2)

        # [修改] 初始化热流源数组 [W]
        # 形状必须与边界一致。我们从 bound1 获取形状。
        self.shape = self.bound1.shape
        self._current_Q_source = np.zeros(self.shape, dtype=float)

    def set_active_heat_source(self, Q_array: Union[float, np.ndarray]):
        """
        [接口] 设置当前时刻的分布式热流源

        :param Q_array: 节点热流数组 [W]。
                        正值 (+): 热量从 Obj1 流出，流入 Obj2 (Obj1 冷却, Obj2 加热)
                        负值 (-): 热量从 Obj2 流出，流入 Obj1

                        如果是标量，则自动广播为均匀热流 (每个节点一样)。
                        如果是数组，形状必须与 boundary.shape 一致。
        """
        if np.isscalar(Q_array):
            self._current_Q_source.fill(Q_array)
        else:
            Q_arr = np.asarray(Q_array, dtype=float)
            if Q_arr.shape != self.shape:
                raise ValueError(f"Heat source shape mismatch: expected {self.shape}, got {Q_arr.shape}")
            self._current_Q_source[:] = Q_arr

    def sync(self):
        """
        同步状态：
        1. 调用父类计算基础热阻 (辐射 // 气体导热)。
        2. 应用戴维南等效修正，注入固定热流源。
        """
        # Step 1: 执行基础物理计算
        # 父类会更新 bc1.R_ext 和 bc1.T_ext (此时 T_ext = T_neighbor_node)
        super().sync()

        # 快速检查：如果全零则跳过
        if np.all(np.abs(self._current_Q_source) < 1e-12):
            return

        # Step 2: 获取父类计算好的总外部回路热阻
        # 2.1. 间隙导热热阻
        R_gap = self.R_gap_total

        # 2.2. 内部导热热阻
        # R_total_1 = R_gap + R_internal_2 (这是 Obj1 看到的总阻力)
        # R_total_2 = R_gap + R_internal_1 (这是 Obj2 看到的总阻力)
        R_total_1 = self.bc1.R_ext
        R_total_2 = self.bc2.R_ext
        #
        # # 2.3. 获取父类设置的基础外部温度 (对方的节点温度)
        T1_node, R1_int = self.bound1.get_coupling_snapshot()
        T2_node, R2_int = self.bound2.get_coupling_snapshot()
        #
        # # 2.4. 计算回路总热阻
        # R_loop_total = R1_int + R2_int

        # Step 3: 应用戴维南等效修正 (Source Transformation)
        # 原理: Q_total = (T_self - T_ext_new) / R_total
        # 目标: Q_total = (T_self - T_neighbor)/R_total + Q_source
        # 推导: T_ext_new = T_neighbor - Q_source * R_total
        # 注意: set_active_heat_source 定义 Q 正方向为 1 -> 2
        # Obj1 (Source端): 为了多流出 Q，对方看起来要更冷 -> T_ext 减小
        # Obj2 (Sink端):   为了接收 Q，对方看起来要更热 -> T_ext 增大
        # 计算戴维南温差修正量
        delta_T = self._current_Q_source * R_gap

        # Step 4: 应用修正到边界
        # 获取基础外部温度 (此时 T_ext 正好等于对方的 T_node)
        T_ext_base_1 = T2_node
        T_ext_base_2 = T1_node

        # 应用修正
        self.bc1.update_params(T_ext=T_ext_base_1 - delta_T, R_ext=R_total_1)
        self.bc2.update_params(T_ext=T_ext_base_2 + delta_T, R_ext=R_total_2)


class TECCouple2D(GapCouple2D):
    """
    热电转换耦合器 (TEC Couple)
    继承自 GapCouple2D，支持非保守的非对称热流注入。

    功能：
    1. 计算间隙的辐射与气体导热 (沿用 GapCouple2D)。
    2. 允许在间隙两侧分别施加独立的热流源 (Source Terms)，以模拟：
       - 发射极 (Side 1): 电子冷却 (Electron Cooling) -> 热流流出
       - 接收极 (Side 2): 电子加热 (Electron Heating) -> 热流流入

    符号约定 (Sign Convention):
    - Q > 0: 热量流入固体 (Heating/Inflow)
    - Q < 0: 热量流出固体 (Cooling/Outflow)

    典型用法:
    - Emitter Side (Side 1): 输入负值 (e.g., -10 W/cm2)
    - Collector Side (Side 2): 输入正值 (e.g., +6 W/cm2)
    """

    def __init__(self,
                 obj1, obj2,
                 direction: Literal['right', 'left', 'top', 'bottom'],
                 gap_width: float,
                 gas_conductivity: float,
                 emissivity1: float = 0.8,
                 emissivity2: float = 0.8):
        super().__init__(obj1, obj2, direction, gap_width, gas_conductivity,
                         emissivity1, emissivity2)

        # 初始化两侧热流源缓存 (单位: W)
        # 形状与边界一致 (假设 1:1 映射)
        self.shape = self.bound1.shape
        self.Q_source_1 = np.zeros(self.shape, dtype=float)  # Emitter 侧
        self.Q_source_2 = np.zeros(self.shape, dtype=float)  # Collector 侧

    def set_tec_sources(self, Q_emitter: np.ndarray, Q_collector: np.ndarray):
        """
        设置热电转换产生的表面热流源。

        :param Q_emitter: 作用于 Obj1 (Emitter) 表面的热流 [W]。
                          电子冷却通常导致热量流失，应为负值 (-)。
        :param Q_collector: 作用于 Obj2 (Collector) 表面的热流 [W]。
                            电子加热导致热量获得，应为正值 (+)。
        """
        # 1. 转换并校验 Side 1 (Emitter)
        Q1 = np.asarray(Q_emitter, dtype=float)
        if Q1.shape != self.shape:
            # 尝试广播标量
            if Q1.ndim == 0:
                Q1 = np.full(self.shape, Q1)
            else:
                raise ValueError(f"[TECCouple2D] Emitter source shape mismatch: {Q1.shape} vs {self.shape}")

        # 2. 转换并校验 Side 2 (Collector)
        Q2 = np.asarray(Q_collector, dtype=float)
        if Q2.shape != self.shape:
            if Q2.ndim == 0:
                Q2 = np.full(self.shape, Q2)
            else:
                raise ValueError(f"[TECCouple2D] Collector source shape mismatch: {Q2.shape} vs {self.shape}")

        # 3. 存储
        self.Q_source_1[:] = Q1
        self.Q_source_2[:] = Q2

    def sync(self):
        super().sync()

        if np.all(self.Q_source_1 == 0) and np.all(self.Q_source_2 == 0):
            return

        R_ext_1 = self.bc1.R_ext

        R_ext_2 = self.bc2.R_ext

        _, R_int_1 = self.bound1.get_coupling_snapshot()
        _, R_int_2 = self.bound2.get_coupling_snapshot()

        T_node_2 = self.bc1.T_ext  # Obj1 面对的是 Obj2 的节点
        T_node_1 = self.bc2.T_ext  # Obj2 面对的是 Obj1 的节点

        T_ext_1_new = T_node_2 + (self.Q_source_1 * R_ext_1) + (self.Q_source_2 * R_int_2)
        T_ext_2_new = T_node_1 + (self.Q_source_2 * R_ext_2) + (self.Q_source_1 * R_int_1)

        self.bc1.update_params(T_ext=T_ext_1_new, R_ext=R_ext_1)
        self.bc2.update_params(T_ext=T_ext_2_new, R_ext=R_ext_2)
