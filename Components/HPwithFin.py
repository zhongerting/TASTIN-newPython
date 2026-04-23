import numpy as np

from Components.BaseComponent import BaseComponent
from Components.basicComponents.HeatPipe2D import HeatPipe2D


class HPwithFin(BaseComponent):
    """
    带降维翅片的热管散热器组件。

    这个类的核心思想是：
    1. 热管本体仍然由 ``HeatPipe2D`` 显式求解。
    2. 翅片不再作为独立二维/三维网格离散，而是降维成“沿翅高方向的一维准稳态导热问题”。
    3. 翅片求解完成后，被等效为挂在冷凝段外壁边界上的一个变热阻支路。

    这样做的好处是：
    1. 保留翅片沿高度方向温度衰减的主要物理效应。
    2. 不破坏现有 ``BoundaryRegion`` 的边界叠加机制。
    3. 方便在同一套框架中同时支持“简化折算外热流”和“翅片直接受照外热流”两种建模。
    """

    def __init__(self,
                 name: str,
                 r_out_wall: float, r_in_wall: float, r_vapor: float,
                 L_eva: float, L_aba: float, L_con: float,
                 n_eva: int, n_aba: int, n_con: int,
                 n_wick: int, n_wall: int,
                 wall_mat, fluid_mat, wick_struct_mat, porosity: float,
                 fin_thickness: float, fin_height: float, n_fin_height: int,
                 fin_wrap_ratio: float,
                 emissivity: float = 0.8,
                 up_view_factor: float = 1.0,
                 down_view_factor: float = 1.0,
                 T_env: float = 3.0,
                 initial_temp: float = 298.15):
        super().__init__(name)

        # ===== 1. 保存翅片与环境参数 =====
        self.fin_thickness = fin_thickness
        self.fin_height = fin_height
        self.n_fin_height = int(n_fin_height)
        self.fin_wrap_ratio = fin_wrap_ratio

        # 接触热阻按单位根部面积给出，这里先保留接口，默认取理想接触。
        # 单位面积接触热阻 [K·m²/W]
        self.R_contact_sq = 0.0

        self.emissivity = emissivity
        self.T_space = T_env
        self.k_fin_mat = wall_mat
        self.sigma = 5.670374419e-8

        # 上下视角因子折算到一个平均有效发射率。
        self.up_view_factor = up_view_factor
        self.down_view_factor = down_view_factor
        self.inner_view_factor = self.up_view_factor + self.down_view_factor
        if self.inner_view_factor > 1.0 + 1.0e-12:
            raise ValueError(
                "Invalid HPwithFin inner-side view factors: "
                "up_view_factor + down_view_factor must be <= 1.0. "
                f"Got up={self.up_view_factor}, down={self.down_view_factor}."
            )
        if self.inner_view_factor < -1.0e-12:
            raise ValueError(
                "Invalid HPwithFin inner-side view factors: "
                "up_view_factor + down_view_factor must be >= 0.0. "
                f"Got up={self.up_view_factor}, down={self.down_view_factor}."
            )
        # 外侧半周完全向太空辐射，内侧半周按 up/down 角系数之和折算。
        self.F_avg = (1.0 + self.inner_view_factor) / 2.0
        self.effective_emissivity = self.emissivity * self.F_avg

        # ===== 2. 构造热管二维网格 =====
        # 轴向由蒸发段、绝热段、冷凝段拼接而成。
        faces_eva = np.linspace(0.0, L_eva, n_eva + 1)
        faces_aba = np.linspace(L_eva, L_eva + L_aba, n_aba + 1)[1:]
        faces_con = np.linspace(L_eva + L_aba, L_eva + L_aba + L_con, n_con + 1)[1:]
        y_faces_custom = np.concatenate([faces_eva, faces_aba, faces_con])

        # 径向由吸液芯区和管壁区拼接而成。
        faces_wick = np.linspace(r_vapor, r_in_wall, n_wick + 1)
        faces_wall = np.linspace(r_in_wall, r_out_wall, n_wall + 1)[1:]
        x_faces_custom = np.concatenate([faces_wick, faces_wall])

        from Solvers.HeatConduction.Mesh import Mesh2D

        n_x = n_wick + n_wall
        n_y = n_eva + n_aba + n_con

        self.hp_mesh = Mesh2D(
            x_dim=0.0,
            n_x=n_x,
            y_dim=0.0,
            n_y=n_y,
            geometry_type='cylindrical',
            inner_radius=r_vapor,
            y_faces=y_faces_custom,
            x_faces=x_faces_custom
        )

        # ===== 3. 实例化热管本体 =====
        self.hp = HeatPipe2D(
            name=f"{self.name}_HP_inner",
            mesh=self.hp_mesh,
            solid1=wall_mat,
            solid2=fluid_mat,
            solid3=wick_struct_mat,
            n_wick=n_wick,
            porosity=porosity,
            n_eva=n_eva,
            n_aba=n_aba,
            n_con=n_con,
            emissivity=emissivity,
            up_view_factor=up_view_factor,
            down_view_factor=down_view_factor,
            initial_temp=initial_temp
        )

        area_aba = self.hp.boundaries['outer_aba'].area
        area_con = self.hp.boundaries['outer_con'].area

        # ===== 4. 挂载外边界 =====
        # 绝热段和冷凝段裸壁对外辐射都直接挂在 BoundaryRegion 上。
        self.bc_rad_aba = self.hp.boundaries['outer_aba'].add_dynamic_radiation_condition(
            emissivity=self.effective_emissivity,
            bare_area_array=area_aba,
            T_env=self.T_space
        )

        self.bc_rad_con = self.hp.boundaries['outer_con'].add_dynamic_radiation_condition(
            emissivity=self.effective_emissivity,
            bare_area_array=area_con,
            T_env=self.T_space
        )

        # 翅片支路不是显式网格，而是等效成一条“外部热阻边界”。
        # pre_step() 每一步都会重新计算这条支路的等效热阻。
        # 翅片支路不是显式网格，而是等效成一条"外部热阻边界"。
        # pre_step() 每一步都会重新计算这条支路的等效热阻。
        # 这里初始化时设置一个极大的热阻(1e15 K/W)，相当于"开路"状态，
        # 确保在第一次pre_step()执行前不会有热量通过翅片支路。
        # 这个值会在“pre_step()”中根据翅片准稳态求解结果被动态修正。

        # TODO: 翅片支路的热阻计算需要根据实际物理模型进行。
        # 当前这里简单地设置为一个极大的热阻(1e15 K/W)，相当于"开路"状态，
        # 确保在第一次pre_step()执行前不会有热量通过翅片支路。
        # 这个值会在“pre_step()”中根据翅片准稳态求解结果被动态修正。
        self.bc_fin_con = self.hp.boundaries['outer_con'].add_resistance_condition(
            T_ext=self.T_space,
            R_ext=1e15
        )

        # ===== 5. 翅片几何折算 =====
        # area_fin_root: 每个轴向切片上，真正长出翅片的根部面积。
        # area_con是数组，fin_wrap_ratio是实数，结果area_fin_root是数组
        area_fin_root = area_con * self.fin_wrap_ratio

        # 降维后的一维翅片方程需要一个“翅片宽度”，
        # 这里由根部面积 / 厚度得到等效宽度。
        # fin_width_array 是降维后一维翅片模型的"特征宽度"。
        # 注意这里不是几何意义上的"翅片宽度"，而是为了保持导热截面积
        # 和辐射周长计算正确而引入的等效参数。
        #
        # 对于矩形截面翅片，导热截面积 A_c = width * thickness，
        # 辐射周长 P = 2 * width（双面辐射）。
        # 因此 width = A_c / thickness = area_fin_root / thickness。
        #
        # 这里 area_fin_root 是根部实际面积（考虑 fin_wrap_ratio 后），
        # 除以厚度得到等效宽度，确保后续一维导热方程的截面积和周长计算
        # 与真实三维几何保持热等效，而不是直接使用某个几何宽度。
        self.fin_width_array = area_fin_root / self.fin_thickness

        # 这里定义的是“一侧受照的投影面积”，用于轨道外热流加载，
        # 对当前简化模型来说，比直接使用双面辐射面积更合适。
        self.fin_illuminated_area_array = self.fin_width_array * self.fin_height

        # 保留 R_contact_abs 接口以支持未来扩展：
        # 1. 允许用户通过 set_contact_resistance() 等方法动态设置接触热阻
        # 2. 支持非理想接触情况下的热阻叠加计算
        # 3. 保持与现有边界条件接口的兼容性，便于后续模型升级
        # 当前默认理想接触（R_contact_sq = 0），但接口保留供未来使用
        # R_contact_sq 是单位面积接触热阻 [K·m²/W]
        # area_fin_root 是根部面积 [m²]
        # 因此 R_contact_abs 的单位是 [K/W]，符合热阻的量纲要求
        with np.errstate(divide='ignore', invalid='ignore'):
            self.R_contact_abs = self.R_contact_sq / area_fin_root
            self.R_contact_abs = np.nan_to_num(self.R_contact_abs, posinf=0.0)

        # ===== 6. 翅片直接受照模式所需状态 =====
        # 如果不设置 fin_external_heat_source，就退化成“翅片不直接吸热”的旧模型。
        self.fin_external_heat_source = None
        self.fin_external_area_scale = 1.0

        # ===== 7. 后处理核算所需状态 =====
        # 这一组量只用于统计“吸收了多少轨道热流”，
        # 不会反向影响主求解。
        self.external_heat_accounting_source = None
        self.wall_external_absorption_area = np.zeros_like(area_con)
        self.fin_external_absorption_area = np.zeros_like(area_con)

        # ===== 8. 缓存上一步翅片求解结果，供后处理直接使用 =====
        self.last_fin_temperature = np.zeros((len(area_con), self.n_fin_height))
        self.last_fin_radiation_distribution = np.zeros_like(area_con)
        self.last_fin_absorption_distribution = np.zeros_like(area_con)
        self.last_fin_net_from_root_distribution = np.zeros_like(area_con)
        self.last_fin_conductance_distribution = np.zeros_like(area_con)
        self.last_fin_effective_temperature_distribution = np.full_like(area_con, self.T_space)
        self.last_fin_equivalent_resistance_distribution = np.full_like(area_con, 1e15)

    def set_fin_external_heat_source(self, heat_source, illuminated_area_scale: float = 1.0):
        """
        为降维翅片模型挂载“直接受照”的外热流源。

        这个接口对应更严格的建模方法：
        轨道外热流不再全部折算到冷凝段根部边界，而是允许其中一部分直接加到翅片方程里。
        """
        if heat_source is None:
            self.fin_external_heat_source = None
            self.fin_external_area_scale = 1.0
            self.last_fin_absorption_distribution.fill(0.0)
            return

        expected_shape = self.hp.boundaries['outer_con'].shape
        if heat_source.shape != expected_shape:
            raise ValueError(
                f"Fin external heat source shape {heat_source.shape} mismatch with condenser boundary shape {expected_shape}"
            )

        self.fin_external_heat_source = heat_source
        self.fin_external_area_scale = float(illuminated_area_scale)

    def configure_external_heat_accounting(self,
                                           heat_source,
                                           wall_area_array: np.ndarray,
                                           fin_area_array: np.ndarray = None):
        """
        配置轨道外热流吸收核算信息。

        这里的目标不是施加边界，而是记录：
        1. 管壁按多大面积吸热。
        2. 翅片按多大面积吸热。
        这样后处理时就能把“总散热量”和“吸收的轨道热流”拆开看。
        """
        if heat_source is None:
            self.external_heat_accounting_source = None
            self.wall_external_absorption_area.fill(0.0)
            self.fin_external_absorption_area.fill(0.0)
            return

        # 获取冷凝段外边界的期望形状，用于验证输入热流源的维度一致性
        expected_shape = self.hp.boundaries['outer_con'].shape
        # 检查热流源形状是否与期望形状匹配，若不匹配则抛出异常
        if heat_source.shape != expected_shape:
            raise ValueError(
                f"Accounting heat source shape {heat_source.shape} mismatch with condenser boundary shape {expected_shape}"
            )

        wall_area = np.array(wall_area_array, dtype=float, copy=True)
        if wall_area.shape != expected_shape:
            wall_area = np.broadcast_to(wall_area, expected_shape)
            wall_area = np.array(wall_area, dtype=float, copy=True)

        if fin_area_array is None:
            fin_area = np.zeros(expected_shape, dtype=float)
        else:
            fin_area = np.array(fin_area_array, dtype=float, copy=True)
            if fin_area.shape != expected_shape:
                fin_area = np.broadcast_to(fin_area, expected_shape)
                fin_area = np.array(fin_area, dtype=float, copy=True)

        self.external_heat_accounting_source = heat_source
        self.wall_external_absorption_area[:] = wall_area
        self.fin_external_absorption_area[:] = fin_area

    def get_fin_illuminated_area_array(self, illuminated_area_scale: float = 1.0) -> np.ndarray:
        """
        返回翅片受照投影面积。

        illuminated_area_scale 用于后续扩展：
        例如只照到一面、部分遮挡、姿态修正等，都可以通过这个系数折算。
        """
        return np.array(self.fin_illuminated_area_array * float(illuminated_area_scale), copy=True)

    def get_external_heat_absorption_distribution(self, current_time: float):
        """
        计算当前时刻轨道外热流在管壁、翅片和总体上的吸收功率。

        返回值单位都是 W，而不是 W/m^2。
        """
        if self.external_heat_accounting_source is None:
            zero = np.zeros_like(self.wall_external_absorption_area)
            return zero.copy(), zero.copy(), zero.copy()

        q_density = np.asarray(
            self.external_heat_accounting_source.get_heat_flux(current_time),
            dtype=float
        )
        wall_absorption = q_density * self.wall_external_absorption_area
        fin_absorption = q_density * self.fin_external_absorption_area
        total_absorption = wall_absorption + fin_absorption
        return wall_absorption, fin_absorption, total_absorption

    def _solve_fin_quasi_steady(self, T_root: np.ndarray, q_fin_abs_density: np.ndarray):
        """
        对每个冷凝段轴向切片求解降维后的一维翅片问题。

        返回:
        - T_fin: 翅片温度场，形状为 (n_con, n_fin_height)
        - Q_fin_radiation: 翅片向空间辐射的热量 [W]
        - Q_fin_absorption: 翅片直接吸收的热量 [W]
        - Q_fin_net_from_root: 从热管根部抽取的净热量 [W]

        Q_fin_net_from_root 的符号约定为：当热量通过翅片根部离开热管时为正。
        """
        T_root = np.asarray(T_root, dtype=float)
        q_fin_abs_density = np.asarray(q_fin_abs_density, dtype=float)

        Nc = len(T_root)
        Nh = self.n_fin_height
        if Nh <= 0:
            raise ValueError("n_fin_height must be positive for HPwithFin")

        dx = self.fin_height / Nh
        Ac = self.fin_width_array * self.fin_thickness
        P = 2.0 * self.fin_width_array

        T = np.repeat(T_root[:, np.newaxis], Nh, axis=1)
        Q_fin_radiation_array = np.zeros_like(T_root)
        Q_fin_absorption_array = (
            q_fin_abs_density
            * self.fin_illuminated_area_array
            * self.fin_external_area_scale
        )

        active = (
            np.isfinite(T_root)
            & np.isfinite(q_fin_abs_density)
            & (self.fin_height > 0.0)
            & (Ac > 0.0)
            & (P > 0.0)
        )

        if not np.any(active):
            Q_fin_net_from_root_array = Q_fin_radiation_array - Q_fin_absorption_array
            return T, Q_fin_radiation_array, Q_fin_absorption_array, Q_fin_net_from_root_array

        T_active = T[active].copy()
        T_root_active = T_root[active]
        q_abs_density_active = q_fin_abs_density[active]
        Ac_active = Ac[active]
        P_active = P[active]
        k_fin_array = np.full_like(T_root_active, 348.9)

        q_fin_abs_segment = np.repeat(
            (
                q_abs_density_active
                * self.fin_width_array[active]
                * dx
                * self.fin_external_area_scale
            )[:, np.newaxis],
            Nh,
            axis=1
        )

        G = k_fin_array * Ac_active / dx
        G_base = 2.0 * G

        max_iter = 50
        tol = 1e-4
        Na = len(T_root_active)

        for _ in range(max_iter):
            h_rad = self.effective_emissivity * self.sigma * (
                T_active ** 2 + self.T_space ** 2
            ) * (T_active + self.T_space)
            rad_term = h_rad * P_active[:, np.newaxis] * dx

            a = np.zeros((Na, Nh))
            b = np.zeros((Na, Nh))
            c = np.zeros((Na, Nh))
            d = np.zeros((Na, Nh))

            b[:, 0] = G_base + G + rad_term[:, 0]
            c[:, 0] = -G
            d[:, 0] = G_base * T_root_active + rad_term[:, 0] * self.T_space + q_fin_abs_segment[:, 0]

            if Nh > 1:
                a[:, 1:-1] = -G[:, np.newaxis]
                b[:, 1:-1] = 2.0 * G[:, np.newaxis] + rad_term[:, 1:-1]
                c[:, 1:-1] = -G[:, np.newaxis]
                d[:, 1:-1] = rad_term[:, 1:-1] * self.T_space + q_fin_abs_segment[:, 1:-1]

                a[:, -1] = -G
                b[:, -1] = G + rad_term[:, -1]
                d[:, -1] = rad_term[:, -1] * self.T_space + q_fin_abs_segment[:, -1]

            c_prime = np.zeros((Na, Nh))
            d_prime = np.zeros((Na, Nh))

            c_prime[:, 0] = c[:, 0] / b[:, 0]
            d_prime[:, 0] = d[:, 0] / b[:, 0]

            for j in range(1, Nh):
                denom = b[:, j] - a[:, j] * c_prime[:, j - 1]
                c_prime[:, j] = c[:, j] / denom
                d_prime[:, j] = (d[:, j] - a[:, j] * d_prime[:, j - 1]) / denom

            T_new = np.zeros((Na, Nh))
            T_new[:, -1] = d_prime[:, -1]
            for j in range(Nh - 2, -1, -1):
                T_new[:, j] = d_prime[:, j] - c_prime[:, j] * T_new[:, j + 1]

            if np.max(np.abs(T_new - T_active)) < tol:
                T_active = T_new
                break
            T_active = T_new

        T[active] = T_active
        Q_fin_radiation_array[active] = np.sum(
            self.effective_emissivity
            * self.sigma
            * P_active[:, np.newaxis]
            * dx
            * (T_active ** 4 - self.T_space ** 4),
            axis=1
        )

        Q_fin_net_from_root_array = Q_fin_radiation_array - Q_fin_absorption_array
        return T, Q_fin_radiation_array, Q_fin_absorption_array, Q_fin_net_from_root_array

    def pre_step(self, dt: float, current_time: float):
        """
        在每个时间步开始前，求解翅片的准稳态一维导热问题。

        求解结果不会直接替代热管二维求解，而是被转换成一组等效热阻，
        再通过 ``self.bc_fin_con`` 挂回冷凝段边界。

        如果设置了 ``fin_external_heat_source``：
        1. 翅片会直接吸收轨道外热流。
        2. 这个吸热项以沿翅高方向分布源项的形式进入翅片方程。

        如果没有设置：
        1. 翅片只作为向深空辐射散热的支路存在。
        2. 这就对应之前的简化模型。
        """

        T_surf, _ = self.hp.boundaries['outer_con'].get_coupling_surface_snapshot()

        # 初始化翅片外热流密度为零，作为默认状态
        q_fin_abs_density = np.zeros_like(T_surf)
        # 如果设置了翅片外热流源，则获取当前时刻的热流密度
        if self.fin_external_heat_source is not None:
            q_fin_abs_density = np.asarray(
                self.fin_external_heat_source.get_heat_flux(current_time),
                dtype=float
            )

        # 求解翅片准稳态一维导热问题，获取翅片温度场及各热流分量
        (
            T,                              # 翅片温度分布，形状为 (n_con, n_fin_height)
            Q_fin_radiation_array,          # 翅片向空间辐射的总热量 [W]
            Q_fin_absorption_array,         # 翅片直接吸收的轨道外热流 [W]
            Q_fin_net_from_root_array,      # 从热管根部抽取的净热量 [W]
        ) = self._solve_fin_quasi_steady(T_surf, q_fin_abs_density)

        dT = np.maximum(0.1, 1.0e-4 * np.maximum(T_surf, 1.0))
        _, _, _, Q_fin_net_plus = self._solve_fin_quasi_steady(
            T_surf + dT,
            q_fin_abs_density
        )

        with np.errstate(divide='ignore', invalid='ignore'):
            lambda_raw = (Q_fin_net_plus - Q_fin_net_from_root_array) / dT

        lambda_min = 1.0e-12
        lambda_max = 1.0e12

        Nh = self.n_fin_height
        dx = self.fin_height / Nh
        P = 2.0 * self.fin_width_array
        T_space_safe = max(float(self.T_space), 1.0e-3)
        T_seg_safe = np.maximum(T, 1.0e-3)
        A_seg = P[:, np.newaxis] * dx
        lambda_fallback = np.sum(
            self.effective_emissivity
            * self.sigma
            * A_seg
            * (T_seg_safe + T_space_safe)
            * (T_seg_safe ** 2 + T_space_safe ** 2),
            axis=1
        )
        lambda_fallback = np.nan_to_num(lambda_fallback, nan=lambda_min, posinf=lambda_max, neginf=lambda_min)
        lambda_fallback = np.clip(lambda_fallback, lambda_min, lambda_max)

        lambda_valid = np.isfinite(lambda_raw) & (lambda_raw > lambda_min)
        lambda_fin = np.where(lambda_valid, lambda_raw, lambda_fallback)
        lambda_fin = np.clip(lambda_fin, lambda_min, lambda_max)

        R_fin_eq_array = self.R_contact_abs + 1.0 / lambda_fin
        R_fin_eq_array = np.nan_to_num(R_fin_eq_array, nan=1e15, posinf=1e15, neginf=1e15)
        R_fin_eq_array = np.maximum(R_fin_eq_array, 1.0 / lambda_max)

        T_fin_eff_array = T_surf - Q_fin_net_from_root_array * R_fin_eq_array
        T_fin_eff_array = np.nan_to_num(
            T_fin_eff_array,
            nan=float(self.T_space),
            posinf=1.0e12,
            neginf=-1.0e12
        )

        self.bc_fin_con.update_params(T_ext=T_fin_eff_array, R_ext=R_fin_eq_array)

        self.last_fin_temperature = np.array(T, copy=True)
        self.last_fin_radiation_distribution = np.array(Q_fin_radiation_array, copy=True)
        self.last_fin_absorption_distribution = np.array(Q_fin_absorption_array, copy=True)
        self.last_fin_net_from_root_distribution = np.array(Q_fin_net_from_root_array, copy=True)
        self.last_fin_conductance_distribution = np.array(lambda_fin, copy=True)
        self.last_fin_effective_temperature_distribution = np.array(T_fin_eff_array, copy=True)
        self.last_fin_equivalent_resistance_distribution = np.array(R_fin_eq_array, copy=True)

    def get_solids(self) -> list:
        """返回系统管理器需要积分的固体对象。"""
        return [self.hp]

    def get_couplers(self) -> list:
        """当前组件内部没有额外的显式耦合器。"""
        return []

    def get_heat_rejection_distribution(self):
        """
        返回绝热段与冷凝段的总散热分布。

        注意这里的冷凝段结果取的是 gross rejection，也就是：
        裸壁辐射 + 翅片辐射。

        之所以这样定义，是因为我们当前比较两种翅片外热流建模时，
        重点关注的是“总向外散了多少热”。
        """

        T_aba, _ = self.hp.boundaries['outer_aba'].get_coupling_surface_snapshot()
        area_aba = self.hp.boundaries['outer_aba'].area
        q_aba_dist = self.effective_emissivity * self.sigma * area_aba * (T_aba ** 4 - self.T_space ** 4)

        T_con, _ = self.hp.boundaries['outer_con'].get_coupling_surface_snapshot()
        area_con = self.hp.boundaries['outer_con'].area
        q_con_bare_dist = self.effective_emissivity * self.sigma * area_con * (T_con ** 4 - self.T_space ** 4)
        q_con_fin_dist = np.array(self.last_fin_radiation_distribution, copy=True)

        q_con_total_dist = q_con_bare_dist + q_con_fin_dist
        return q_aba_dist, q_con_total_dist

    def get_heat_exchange_breakdown(self):
        """
        返回更细的热量拆分结果，便于调试和验证。

        返回字典中包含：
        1. bare_radiation: 裸壁辐射
        2. fin_radiation: 翅片总辐射
        3. fin_absorption: 翅片直接吸收的轨道外热流
        4. fin_net_from_root: 真正通过根部从壁面抽走的净热量
        5. gross_rejection: 裸壁辐射 + 翅片辐射
        6. net_rejection: 裸壁辐射 + 翅片净抽热
        """
        T_con, _ = self.hp.boundaries['outer_con'].get_coupling_surface_snapshot()
        area_con = self.hp.boundaries['outer_con'].area
        q_con_bare_dist = self.effective_emissivity * self.sigma * area_con * (T_con ** 4 - self.T_space ** 4)
        q_con_fin_rad_dist = np.array(self.last_fin_radiation_distribution, copy=True)
        q_con_fin_abs_dist = np.array(self.last_fin_absorption_distribution, copy=True)
        q_con_fin_net_dist = np.array(self.last_fin_net_from_root_distribution, copy=True)
        q_con_fin_conductance_dist = np.array(self.last_fin_conductance_distribution, copy=True)
        q_con_fin_effective_temp_dist = np.array(self.last_fin_effective_temperature_distribution, copy=True)
        q_con_fin_equiv_resistance_dist = np.array(self.last_fin_equivalent_resistance_distribution, copy=True)
        return {
            'bare_radiation': q_con_bare_dist,
            'fin_radiation': q_con_fin_rad_dist,
            'fin_absorption': q_con_fin_abs_dist,
            'fin_net_from_root': q_con_fin_net_dist,
            'fin_conductance': q_con_fin_conductance_dist,
            'fin_effective_temperature': q_con_fin_effective_temp_dist,
            'fin_equivalent_resistance': q_con_fin_equiv_resistance_dist,
            'gross_rejection': q_con_bare_dist + q_con_fin_rad_dist,
            'net_rejection': q_con_bare_dist + q_con_fin_net_dist
        }

    def get_temperature_distribution(self):
        """返回热管二维温度场，形状为 (径向节点数, 轴向节点数)。"""
        return self.hp.T.reshape(self.hp.shape_nodes)
