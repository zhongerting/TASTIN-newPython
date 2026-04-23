import logging
from typing import Dict, List, Optional, Union
import numpy as np

# 导入核心组件的类型提示 (Type Hinting)
# 注意：使用 TYPE_CHECKING 避免循环导入，或者确保路径正确
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.HeatConduction.HeatConduction import BaseHeatConduction
from Solvers.Couplers import FluidSolidCouple, SolidSolidCouple1D, SolidSolidCouple2D
from Components.BaseComponent import BaseComponent

# 导入剖析器，用于统计各个函数消耗时间
from profiler import TEASAProfiler

# 配置日志
logger = logging.getLogger(__name__)


class SystemManager:
    """
    系统管理器 (System Manager) - The Conductor

    核心职责:
    1. 资源管理: 持有流体网络、固体组件和耦合器的引用。
    2. 调度控制: 协调流体求解器和固体求解器的交错步进 (Staggered Stepping)。
    3. 数据同步: 在物理场之间传递边界条件和源项。
    """

    def __init__(self, fluid_network: HydraulicNetwork, start_time: float = 0.0):
        """
        初始化系统管理器

        :param fluid_network: 必须注入一个已经实例化好的 HydraulicNetwork 对象。
                              它是系统的流体骨架。
        :param start_time: 仿真起始时间 [s]
        """
        # --- 1. 全局时钟 ---
        # self._last_dt = None
        self.global_time = start_time

        # --- 2. 流体求解器 (核心骨架) ---
        # 依赖注入: SystemManager 不负责创建流体网络，只负责持有和调用它
        self.fluid_solver = fluid_network

        # --- 3. 固体组件容器 (插件式架构) ---
        # 使用字典存储，方便通过名字检索 (例如: 'fuel_rod_01', 'emitter_02')
        # Key: Component Name, Value: BaseHeatConduction Instance
        self.solid_components: Dict[str, BaseHeatConduction] = {}

        # --- 4. 耦合器容器 ---
        # 存储所有的流固耦合器 (FluidSolidCouple) 和 固固耦合器 (SolidSolidCouple)
        # 它们负责在 step() 的中间阶段搬运数据
        self.couplers: List[Union[FluidSolidCouple, SolidSolidCouple1D, SolidSolidCouple2D]] = []

        # --- 5. 点堆中子动力学 ---
        # [新增] 点堆动力学求解器引用
        self.point_reactor = None

        # --- 6. 组件模型 ---
        # [新增] 宏观组件容器引用
        self.components: List[BaseComponent] = []

        logger.info(f"SystemManager initialized at t={self.global_time}s")

        logger.info(f"SystemManager initialized at t={self.global_time}s")
        logger.info(
            f"  - Fluid Solver attached: {self.fluid_solver.n_vol} volumes, {self.fluid_solver.n_junc} junctions")

    # =========================================================================
    # A. 注册接口 (Registration API)
    # =========================================================================

    def add_point_reactor(self, reactor):
        """
        [新增] 注册点堆动力学模块
        """
        self.point_reactor = reactor
        logger.info("  [Register] Neutronics: PointReactor attached to SystemManager.")

    def add_solid_component(self, component: BaseHeatConduction):
        """
        注册固体导热组件到系统管理器。

        :param component: 继承自 BaseHeatConduction 的实例
                          (例如: HeatConduction2D, FuelRod, Emitter 等)
        """
        # 1. 类型安全检查
        if not isinstance(component, BaseHeatConduction):
            raise TypeError(f"Component must inherit from BaseHeatConduction, got {type(component)}")

        # 2. 命名冲突检查
        # 组件名称将作为字典的 Key，必须唯一
        name = getattr(component, 'name', None)
        if not name:
            # 如果组件没有 name 属性，尝试生成一个默认名 (虽然我们在设计中未强制要求 name)
            # 建议在 BaseHeatConduction 或其子类中添加 self.name
            name = f"Solid_{len(self.solid_components)}"
            logger.warning(f"Component {component} has no 'name' attribute. Assigned: {name}")

        if name in self.solid_components:
            raise ValueError(f"Solid component with name '{name}' already exists in SystemManager!")

        # 3. 注册
        self.solid_components[name] = component

        # 4. 记录日志
        # 打印组件的具体类型 (比如 'FuelRod' 而不仅仅是 'HeatConduction2D')
        logger.info(f"  [Register] Solid: '{name}' (Type: {type(component).__name__}, Nodes: {component.N})")

    def add_coupler(self, coupler: Union[FluidSolidCouple, SolidSolidCouple1D, SolidSolidCouple2D]):
        """
        注册耦合器 (流-固 或 固-固)。
        """
        # 1. 鸭子类型检查 (Duck Typing Check)
        # 确保耦合器具有我们需要的方法 (execute 或 sync)
        if not hasattr(coupler, 'execute') and not hasattr(coupler, 'sync'):
            logger.warning(
                f"⚠️ Registered object {coupler} does not look like a valid coupler (no execute/sync method).")

        # 2. 注册
        self.couplers.append(coupler)

        # 3. 记录日志
        # 尝试获取 name，如果没有则用类名
        c_name = getattr(coupler, 'name', f"Coupler_{len(self.couplers)}")
        logger.info(f"  [Register] Coupler: '{c_name}' (Type: {type(coupler).__name__})")

    def add_component(self, component: BaseComponent):
        """
        [新增] 核心收件逻辑：接收宏观组件，自动拆包并注册到底层清单
        """
        # 1. 将组件本身注册到宏观名册
        self.components.append(component)

        # 2. 自动榨取底层固体，编入导热大军
        for solid in component.get_solids():
            self.add_solid_component(solid)

        # 3. 自动榨取底层耦合器，编入换热搬运工大军
        for coupler in component.get_couplers():
            self.add_coupler(coupler)

        logger.info(f"  [Register] Component: '{component.name}' added. "
                    f"Extracted {len(component.get_solids())} solids and {len(component.get_couplers())} couplers.")

    # =========================================================================
    # B. 初始化调度 (Initialization)
    # =========================================================================

    def initialize_system(self, dt_init: float = 0.1, tol: float = 1e-5, max_iter: int = 500):
        """
        系统初始化：负责将流体网络和固体组件带入一个物理一致的初始状态。

        流程:
        1. 液力找平: 在冻结温度场的情况下，迭代计算流场，直到 dW/dt = 0。
        2. 耦合预热: 执行一次耦合器，根据初始流速计算 h，并同步边界条件。

        :param dt_init: 用于液力初始化的虚拟时间步长 [s]
        :param tol: 液力收敛容差
        :param max_iter: 最大迭代次数
        """
        logger.info(f"System Initialization Started (dt={dt_init}, tol={tol})...")

        # --- 1. 流体网络：纯液力初始化 (Hydraulic Steady-State) ---
        # 这一步是为了消除初始流速猜测值不准导致的剧烈非物理震荡
        logger.info("  >> [Step 1] Initializing Hydraulics...")
        success = self.fluid_solver.initialize_hydraulics(dt=dt_init, tol=tol, max_iter=max_iter)

        if not success:
            logger.error("Hydraulic initialization FAILED! Check topology or boundary conditions.")
            # 根据策略，这里可以选择抛出异常，或者让用户决定是否继续
            raise RuntimeError("System initialization failed at hydraulic stage.")

        # --- 2. 耦合器：初始状态同步 (Coupling Pre-Sync) ---
        # 此时流速 W 已经稳定，我们需要利用这个 W 计算初始的换热系数 h 和 Re
        # 并将流体温度 T_f 推送给固体边界，将固体壁温 T_w 推送给流体源项
        logger.info("  >> [Step 2] Synchronizing Couplers...")
        self._sync_solid_boundaries_for_coupling()

        for coupler in self.couplers:
            if hasattr(coupler, 'execute'):
                # 流固耦合：计算 h, Nu, 更新 BC
                coupler.execute()
            elif hasattr(coupler, 'sync'):
                # 固固耦合：同步接触面温度
                coupler.sync()

        # --- 3. (可选) 固体热工检查 ---
        # 目前固体组件仅使用 initial_temp 初始化。
        # 如果需要稳态热工初始化 (Steady-State Thermal)，需要求解 d/dt = 0 的热传导方程。
        # 这是一个高阶功能，目前我们假设用户给定的初始温度是合理的，或者前几秒作为瞬态建立过程。
        pass

        logger.info("System Initialization Completed Successfully.\n")

    # =========================================================================
    # C. 瞬态步进 (Transient Stepping)
    # =========================================================================

    def _clear_fluid_sources(self):
        """
        [辅助] 清空所有流体节点的源项 (Q_wall, implicit_coeff)

        必须在执行耦合器 (Coupler.execute) 之前调用。
        因为耦合器是累加式写入 (+=)，如果不清空，源项会无限累积。
        """
        # 遍历流体网络中的所有节点
        for vol in self.fluid_solver.volumes_obj:
            # 清除壁面热流 (Explicit Source)
            if hasattr(vol, 'Q_wall'):
                vol.Q_wall = 0.0

            # 清除体积热源 (注意：如果有固定的核热/衰变热，需要在外部主循环中重新施加)
            if hasattr(vol, 'Q_vol'):
                vol.Q_vol = 0.0

            # 清除隐式耦合系数 (Implicit Source Coeff)
            if hasattr(vol, 'implicit_coeff'):
                vol.implicit_coeff = 0.0

    def _sync_solid_boundaries_for_coupling(self):
        """Push current solid states into BoundaryRegion caches before couplers read them."""
        for solid in self.solid_components.values():
            if hasattr(solid, '_update_properties'):
                solid._update_properties()
            if hasattr(solid, '_compute_internal_resistance'):
                solid._compute_internal_resistance()
            if hasattr(solid, '_update_boundaries_state'):
                solid._update_boundaries_state(current_time=self.global_time)

    @TEASAProfiler.profile
    def step(self, dt: float, inner_iter: int = 1, convergence_tol: float = 1e-3, reactivity_control: float = 0.0):
        """
        [修正版 6.0] 执行一个时间步长的系统演化 (Global Time Step)

        改进点:
        1. 采用 "Total Flux Averaging" (总热流平均) 策略。
        2. 利用预测温度 T* 计算校正步的衰减项，实现伪隐式 (Pseudo-Implicit) 稳定性。
        3. 解决高换热系数下的数值震荡问题。
        """

        # --- 1. 触发所有组件的 pre_step 方法
        for comp in self.components:
            if hasattr(comp, 'pre_step'):
                comp.pre_step(dt, self.global_time)

        # --- 2. 保存现场 (Snapshot at t_n) ---
        if inner_iter > 1:
            self._save_system_state()
            # [新增] 初始化总热流备份列表
            self._fluid_total_Q_backup = []

        # 用于收敛性检查
        T_f_prev = None
        T_s_prev = {}

        # --- 3. 内部迭代循环 ---
        for k in range(inner_iter):

            # [Step A] 准备源项 (Coupling Phase)
            # -------------------------------------------------------
            self._clear_fluid_sources()
            self._sync_solid_boundaries_for_coupling()

            # 执行耦合，计算当前状态下的 Q_wall (Explicit) 和 implicit_coeff
            # k=0: 基于 t_n
            # k=1: 基于 t_n+1 (Predictor)
            for coupler in self.couplers:
                if hasattr(coupler, 'execute'):
                    coupler.execute()
                elif hasattr(coupler, 'sync'):
                    coupler.sync()

            # [Step B] 源项稳定化处理 (Stabilization)
            # -------------------------------------------------------
            if inner_iter > 1:
                # 获取流体网络的所有节点
                volumes = self.fluid_solver.volumes_obj

                # 计算当前时刻的“总净热流” (Total Net Heat Flux)
                # Formula: Q_total = Q_explicit - lambda * T_current
                # 注意：
                #  - k=0 时，T_current 是 T^n
                #  - k=1 时，T_current 是 T* (预测温度) <--- 关键点！

                current_total_Q_list = []
                for i, vol in enumerate(volumes):
                    q_exp = getattr(vol, 'Q_wall', 0.0)
                    lam = getattr(vol, 'implicit_coeff', 0.0)
                    # 直接从对象读取当前温度 (HydraulicNetwork 保证了 vol.T 是最新的)
                    t_curr = vol.T

                    # 计算净热流 [W]
                    q_net = q_exp - lam * t_curr
                    current_total_Q_list.append(q_net)

                if k == 0:
                    # 备份预测步的净热流
                    self._fluid_total_Q_backup = current_total_Q_list
                else:
                    # 校正步：取平均 (Trapezoidal Rule)
                    # Q_final = 0.5 * (Q_net^n + Q_net^*)
                    for i, vol in enumerate(volumes):
                        avg_Q = 0.5 * (self._fluid_total_Q_backup[i] + current_total_Q_list[i])

                        # [关键操作]
                        # 1. 将平均后的总热流作为纯显式源项写入 Q_wall
                        vol.Q_wall = avg_Q

                        # 2. 强制将隐式系数置零
                        #    原因：我们已经在 avg_Q 中包含了两端的 -lambda*T 效应。
                        #    如果不置零，求解器在 step 时会计算 avg_Q - lambda * T^n，导致重复扣除。
                        vol.implicit_coeff = 0.0

            # [Step C] 状态回滚 (Rollback)
            # -------------------------------------------------------
            # 带着计算好的平均源项 (avg_Q)，回到 t_n 重新出发
            if k > 0:
                self._rollback_system_state()

            # [新增 Step C.2] 中子动力学积分与功率下发
            # -------------------------------------------------------
            # 回滚后，在真正的热工水力求解之前，先由宏观组件更新中子状态与功率源项。
            component_neutronics_handled = False
            for comp in self.components:
                if hasattr(comp, 'advance_neutronics'):
                    comp.advance_neutronics(
                        dt=dt,
                        reactivity_control=reactivity_control,
                        iteration_index=k
                    )
                    component_neutronics_handled = True

            # 向后兼容旧链路：若没有组件接管中子学，仍允许使用 SystemManager.point_reactor。
            if self.point_reactor is not None and not component_neutronics_handled:
                total_feedback = 0.0
                for name, solid in self.solid_components.items():
                    if hasattr(solid, 'get_reactivity_feedback'):
                        total_feedback += solid.get_reactivity_feedback()

                self.point_reactor.step(dt, reactivity_control, total_feedback)

                p_fiss = self.point_reactor.fission_power
                p_decay = self.point_reactor.decay_power
                p_total = self.point_reactor.total_power

                for name, solid in self.solid_components.items():
                    if hasattr(solid, 'set_nuclear_power'):
                        solid.set_nuclear_power(p_fiss, p_decay, p_total)

            # [Step D] 求解器步进 (Solving)
            # -------------------------------------------------------
            # 1. 流体步进
            fluid_converged = self.fluid_solver.step_Picard(dt, max_iter=20 if inner_iter > 1 else 100)

            if not fluid_converged and inner_iter == 1:
                logger.warning(f"Fluid solver NOT converged at t={self.global_time:.4f}s")

            # 2. 固体步进
            for name, solid in self.solid_components.items():
                success = solid.step(dt)
                if not success:
                    logger.error(f"Solid '{name}' integration failed")

            # [Step E] 收敛性检查
            # -------------------------------------------------------
            if inner_iter > 1:
                T_f_curr = self.fluid_solver.T_vec.copy()
                T_s_curr = {name: solid.T.copy() for name, solid in self.solid_components.items()}

                if k > 0:
                    err_f = np.max(np.abs(T_f_curr - T_f_prev))

                    err_s = 0.0
                    for name, T_arr in T_s_curr.items():
                        if name in T_s_prev:
                            diff = np.max(np.abs(T_arr - T_s_prev[name]))
                            if diff > err_s:
                                err_s = diff

                    if err_f < convergence_tol and err_s < convergence_tol:
                        break

                T_f_prev = T_f_curr
                T_s_prev = T_s_curr

        # --- 4. 固化中子先驱核与系统时间推进 ---
        component_neutronics_committed = False
        for comp in self.components:
            if hasattr(comp, 'commit_neutronics'):
                comp.commit_neutronics()
                component_neutronics_committed = True

        if self.point_reactor is not None and not component_neutronics_committed:
            # 当 Picard 迭代跳出（收敛）后，固化当前时间步的真实状态
            self.point_reactor.commit()

        # --- 5. 触发所有组件的 post_step 方法
        for comp in self.components:
            if hasattr(comp, 'post_step'):
                comp.post_step(dt, self.global_time)

        # --- 6. 全局时间推进 ---
        self.global_time += dt

    # def step(self, dt: float, inner_iter: int = 1, convergence_tol: float = 1e-3):
    #     """
    #     [修正版 6.0] 执行一个时间步长的系统演化 (Global Time Step)
    #
    #     改进点:
    #     1. 采用 "Total Flux Averaging" (总热流平均) 策略。
    #     2. 利用预测温度 T* 计算校正步的衰减项，实现伪隐式 (Pseudo-Implicit) 稳定性。
    #     3. 解决高换热系数下的数值震荡问题。
    #     """
    #
    #     # --- 1. 保存现场 (Snapshot at t_n) ---
    #     if inner_iter > 1:
    #         self._save_system_state()
    #         # [新增] 初始化总热流备份列表
    #         self._fluid_total_Q_backup = []
    #
    #     # 用于收敛性检查
    #     T_f_prev = None
    #     T_s_prev = {}
    #
    #     # --- 2. 内部迭代循环 ---
    #     for k in range(inner_iter):
    #
    #         # [Step A] 准备源项 (Coupling Phase)
    #         # -------------------------------------------------------
    #         self._clear_fluid_sources()
    #
    #         # 执行耦合，计算当前状态下的 Q_wall (Explicit) 和 implicit_coeff
    #         # k=0: 基于 t_n
    #         # k=1: 基于 t_n+1 (Predictor)
    #         for coupler in self.couplers:
    #             if hasattr(coupler, 'execute'):
    #                 coupler.execute()
    #             elif hasattr(coupler, 'sync'):
    #                 coupler.sync()
    #
    #         # [Step B] 源项稳定化处理 (Stabilization)
    #         # -------------------------------------------------------
    #         if inner_iter > 1:
    #             # 获取流体网络的所有节点
    #             volumes = self.fluid_solver.volumes_obj
    #
    #             # 计算当前时刻的“总净热流” (Total Net Heat Flux)
    #             # Formula: Q_total = Q_explicit - lambda * T_current
    #             # 注意：
    #             #  - k=0 时，T_current 是 T^n
    #             #  - k=1 时，T_current 是 T* (预测温度) <--- 关键点！
    #
    #             current_total_Q_list = []
    #             for i, vol in enumerate(volumes):
    #                 q_exp = getattr(vol, 'Q_wall', 0.0)
    #                 lam = getattr(vol, 'implicit_coeff', 0.0)
    #                 # 直接从对象读取当前温度 (HydraulicNetwork 保证了 vol.T 是最新的)
    #                 t_curr = vol.T
    #
    #                 # 计算净热流 [W]
    #                 q_net = q_exp - lam * t_curr
    #                 current_total_Q_list.append(q_net)
    #
    #             if k == 0:
    #                 # 备份预测步的净热流
    #                 self._fluid_total_Q_backup = current_total_Q_list
    #             else:
    #                 # 校正步：取平均 (Trapezoidal Rule)
    #                 # Q_final = 0.5 * (Q_net^n + Q_net^*)
    #                 for i, vol in enumerate(volumes):
    #                     avg_Q = 0.5 * (self._fluid_total_Q_backup[i] + current_total_Q_list[i])
    #
    #                     # [关键操作]
    #                     # 1. 将平均后的总热流作为纯显式源项写入 Q_wall
    #                     vol.Q_wall = avg_Q
    #
    #                     # 2. 强制将隐式系数置零
    #                     #    原因：我们已经在 avg_Q 中包含了两端的 -lambda*T 效应。
    #                     #    如果不置零，求解器在 step 时会计算 avg_Q - lambda * T^n，导致重复扣除。
    #                     vol.implicit_coeff = 0.0
    #
    #         # [Step C] 状态回滚 (Rollback)
    #         # -------------------------------------------------------
    #         # 带着计算好的平均源项 (avg_Q)，回到 t_n 重新出发
    #         if k > 0:
    #             self._rollback_system_state()
    #
    #         # [Step D] 求解器步进 (Solving)
    #         # -------------------------------------------------------
    #         # 1. 流体步进
    #         fluid_converged = self.fluid_solver.step_Picard(dt, max_iter=20 if inner_iter > 1 else 100)
    #
    #         if not fluid_converged and inner_iter == 1:
    #             logger.warning(f"Fluid solver NOT converged at t={self.global_time:.4f}s")
    #
    #         # 2. 固体步进
    #         for name, solid in self.solid_components.items():
    #             success = solid.step(dt)
    #             if not success:
    #                 logger.error(f"Solid '{name}' integration failed")
    #
    #         # [Step E] 收敛性检查
    #         # -------------------------------------------------------
    #         if inner_iter > 1:
    #             T_f_curr = self.fluid_solver.T_vec.copy()
    #             T_s_curr = {name: solid.T.copy() for name, solid in self.solid_components.items()}
    #
    #             if k > 0:
    #                 err_f = np.max(np.abs(T_f_curr - T_f_prev))
    #
    #                 err_s = 0.0
    #                 for name, T_arr in T_s_curr.items():
    #                     if name in T_s_prev:
    #                         diff = np.max(np.abs(T_arr - T_s_prev[name]))
    #                         if diff > err_s: err_s = diff
    #
    #                 if err_f < convergence_tol and err_s < convergence_tol:
    #                     break
    #
    #             T_f_prev = T_f_curr
    #             T_s_prev = T_s_curr
    #
    #     # --- 3. 全局时间推进 ---
    #     self.global_time += dt

    def compute_adaptive_dt(self,
                            min_dt: float = 1e-4,
                            max_dt: float = 0.5,
                            safety_factor: float = 0.8) -> float:
        """
        [自适应步长调度] 汇总所有物理限制，计算下一时刻的最佳时间步长

        :param min_dt: 强制下限，防止计算停滞 [s]
        :param max_dt: 强制上限，防止跨度过大导致物理失真 [s]
        :param safety_factor: 全局安全系数 (建议 0.5-0.8)
        :return: 经过修正的最佳时间步长
        """
        # 1. 收集各组件的限制建议 (内部已含各自的安全系数)
        # 流体 CFL 限制
        raw_dt_fluid = self.fluid_solver.get_max_stable_dt(max_limit=max_dt)
        dt_fluid = raw_dt_fluid * safety_factor

        # 耦合稳定性限制
        dt_coupler = max_dt
        for coupler in self.couplers:
            if hasattr(coupler, 'get_max_stable_dt'):
                dt_c = coupler.get_max_stable_dt(safety_factor=safety_factor, max_limit=max_dt)
                dt_coupler = min(dt_coupler, dt_c)

        # [测试用]查看哪个是限制时间步长的根本
        # 若流体时间步长限制起作用，则停止，用于调试
        if abs(dt_fluid) < abs(dt_coupler):
            dt_fluid = dt_fluid

        # 2. 取交集 (最小值)
        dt_target = min(dt_fluid, dt_coupler)

        # 3. 平滑与限幅策略
        # A. 限制在 [min, max] 区间
        dt_clamped = max(min_dt, min(dt_target, max_dt))

        # B. 限制单步增幅，防止步长跳变剧烈导致求解器失稳
        # 允许步长根据物理需要快速减小，但增加时不能超过上一步的 1.2 倍
        if hasattr(self, '_last_dt'):
            dt_final = min(dt_clamped, 1.2 * self._last_dt)
        else:
            dt_final = dt_clamped

        self._last_dt = dt_final
        return dt_final

    # =========================================================================
    # D. 状态管理接口 (State Management for Inner Iteration)
    # =========================================================================

    def _save_system_state(self):
        """
        [接口] 保存所有子系统的当前状态 (Snapshot at t_n)
        用于内部迭代的回滚操作。
        """
        # 1. 保存流体状态
        if hasattr(self.fluid_solver, 'save_state'):
            self.fluid_solver.save_state()

        # 2. 保存固体状态
        for solid in self.solid_components.values():
            if hasattr(solid, 'save_state'):
                solid.save_state()

    def _rollback_system_state(self):
        """
        [接口] 将所有子系统回滚到上一次保存的状态 (Reset to t_n)
        """
        # 1. 回滚流体
        if hasattr(self.fluid_solver, 'load_state'):
            self.fluid_solver.load_state()

        # 2. 回滚固体
        for solid in self.solid_components.values():
            if hasattr(solid, 'load_state'):
                solid.load_state()

    # =========================================================================
    # E. 全局断点续算 (Global Restart Capability)
    # =========================================================================

    def save_global_state(self, filepath: str):
        """
        [全局接口] 保存全系统断点续算状态为高压缩比 .npz 文件

        工作机制：
        采用扁平化收割策略，分别提取流体、固体网格、宏观组件特殊属性以及点堆的状态字典。
        """
        logger.info(f"Saving global state to {filepath} ...")
        global_state = {'System/global_time': np.array([self.global_time])}

        # 1. 记录系统级元数据 (Metadata Header)
        if hasattr(self, '_last_dt'):
            global_state['System/last_dt'] = np.array([self._last_dt])

        # 2. 收集流体网络状态
        if hasattr(self.fluid_solver, 'get_state_dict'):
            global_state.update(self.fluid_solver.get_state_dict(prefix="Fluid"))

        # 3. 收集所有固体组件状态 (网格指纹与核心温度场)
        # 这里的 solid_components 已经包含了从 TFEUnit 等宏观组件中拆包出来的所有底层实体
        for name, solid in self.solid_components.items():
            if hasattr(solid, 'get_state_dict'):
                global_state.update(solid.get_state_dict(prefix=f"Solid_{name}"))

        # 4. 收集宏观组件的特殊状态 (如 TFEUnit 的电场/等离子体缓存数组)
        for comp in self.components:
            if hasattr(comp, 'get_state_dict'):
                global_state.update(comp.get_state_dict(prefix=f"Macro_{comp.name}"))

        # 5. 收集点堆状态 (如果有)
        if self.point_reactor is not None and hasattr(self.point_reactor, 'get_state_dict'):
            global_state.update(self.point_reactor.get_state_dict(prefix="PointReactor"))

        # 一键压缩并落盘 (allow_pickle=False 保证安全性)
        import os
        np.savez_compressed(filepath, **global_state)
        logger.info(f"Global state successfully saved to {filepath}.")

    def load_global_state(self, filepath: str):
        """
        [全局接口] 读取并恢复全系统断点续算状态

        工作机制：
        读取 .npz 文件并将其转化为扁平字典，然后精确路由分发给各个下属对象进行网格校验与恢复。
        """
        import os
        logger.info(f"Loading global state from {filepath} ...")

        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Restart file not found: {filepath}")

        # 使用 allow_pickle=False 强制只读取纯数值数组，防止任意代码执行漏洞
        with np.load(filepath, allow_pickle=False) as data:
            # 将 NpzFile 对象转为标准 Python 字典供各层路由读取
            data_dict = dict(data)

            # 1. 恢复系统级元数据
            if 'System/global_time' in data_dict:
                self.global_time = float(data_dict['System/global_time'][0])
            if 'System/last_dt' in data_dict:
                self._last_dt = float(data_dict['System/last_dt'][0])

            # 2. 路由：恢复流体网络状态
            if hasattr(self.fluid_solver, 'load_state_dict'):
                self.fluid_solver.load_state_dict(data_dict, prefix="Fluid")

            # 3. 路由：恢复固体组件状态
            for name, solid in self.solid_components.items():
                if hasattr(solid, 'load_state_dict'):
                    solid.load_state_dict(data_dict, prefix=f"Solid_{name}")

            # 4. 路由：恢复宏观组件的特殊状态
            for comp in self.components:
                if hasattr(comp, 'load_state_dict'):
                    comp.load_state_dict(data_dict, prefix=f"Macro_{comp.name}")

            # 5. 路由：恢复点堆状态
            if self.point_reactor is not None and hasattr(self.point_reactor, 'load_state_dict'):
                self.point_reactor.load_state_dict(data_dict, prefix="PointReactor")

        logger.info(f"Global state successfully loaded. Resuming from t={self.global_time}s.")
