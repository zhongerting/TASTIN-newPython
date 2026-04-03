import numpy as np
from typing import Union

from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.Boundary import BoundaryRegion
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Base import SolidMaterial


class HeatPipe2D(HeatConduction2D):
    """
    二维柱坐标热管导热求解器 (2D Cylindrical Heat Pipe Solver)

    特殊功能:
    1. 支持吸液芯的变物性与赝热导率 (将通过专用的 WickMaterial 实现)。
    2. 实现轴向边界切片: 将统一的外壁面 ('right') 拆分为蒸发段 ('outer_eva')、
       绝热段 ('outer_aba') 和冷凝段 ('outer_con')，以便独立挂载不同耦合器。
    """

    def __init__(self,
                 mesh: Mesh2D,
                 solid1: SolidMaterial,  # 管壁材料 (Wall Material)
                 solid2: SolidMaterial,  # 热管工质 (Fluid Material)
                 solid3: SolidMaterial,  # 吸液芯结构材料 (Wick Structure Material)
                 n_wick: int,  # 径向属于吸液芯的网格数
                 porosity: float,  # 吸液芯孔隙率
                 n_eva: int,
                 n_aba: int,
                 n_con: int,
                 name: str = "Unnamed_Solid",
                 emissivity: float = 0.8,
                 up_view_factor: float = 1.0,
                 down_view_factor: float = 1.0,
                 initial_temp: float = 298.15):

        # 1. 记录物理参数与几何分界
        self.n_wick = n_wick
        self.wall_mat = solid1

        # 校验径向网格
        if self.n_wick >= mesh.n_x or self.n_wick <= 0:
            raise ValueError(
                f"Wick radial nodes (n_wick={n_wick}) must be strictly between 0 and total radial nodes (mesh.n_x={mesh.n_x})")

        # 2. 自动从网格中提取吸液芯的几何边界
        # 柱坐标系下，mesh.inner_radius 即为蒸汽腔外半径 Rv
        # mesh.x_faces 记录了所有的径向界面位置，索引 n_wick 处即为吸液芯与管壁的交界 Rw
        r_vapor = mesh.inner_radius
        r_in_wall = mesh.x_faces[self.n_wick]

        # 3. 实例化吸液芯复合材料 (必须在 super().__init__ 之前完成！)
        self.wick_mat = WickMaterial(
            name="HP_Wick_Composite",
            solid_mat=solid3,
            fluid_mat=solid2,
            porosity=porosity,
            r_vapor=r_vapor,
            r_in_wall=r_in_wall
        )

        # ==================== 记录切片参数与辐射参数 ====================
        self.n_eva = n_eva
        self.n_aba = n_aba
        self.n_con = n_con

        self.emissivity = emissivity
        self.up_view_factor = up_view_factor
        self.down_view_factor = down_view_factor

        if self.n_eva + self.n_aba + self.n_con != mesh.n_y:
            raise ValueError(f"Axial sections sum ({n_eva}+{n_aba}+{n_con}) "
                             f"must equal total mesh n_y ({mesh.n_y})")

        # ==================== 提前定义 shape_nodes 网格形状 避免报错 ====================
        self.shape_nodes = mesh.shape_nodes

        # ==================== 调用父类初始化 ====================
        # 注意：这里传入 wall_mat 作为基类的占位符。
        # 父类初始化时会调用 self.initialize_state() -> self._update_properties()
        # 由于我们重写了 _update_properties，真实的物性将由我们的分层逻辑接管！
        super().__init__(mesh, self.wall_mat, name=name, initial_temp=initial_temp)

        # 构建虚拟切片边界，覆盖父类的统一 right 边界
        self._setup_virtual_boundaries()

    def _setup_virtual_boundaries(self):
        """
        [新增核心方法] 构建边界切片
        逻辑：
        屏蔽或移除父类原有的 'right' 边界对象。
        创建三个新的 BoundaryRegion 对象：'outer_eva', 'outer_aba', 'outer_con'。
        并从原始的 area_x_matrix 中切片提取对应段的表面积，赋给这三个新边界。
        """

        # 1. 获取外壁面 (X方向最外侧) 每个控制体的径向对外表面积
        # mesh.area_x_matrix 的 shape 为 (nx+1, ny)
        # 索引 -1 代表最右侧界面，包含了 ny 个控制体的外表面积
        outer_surface_areas = self.mesh.area_x_matrix[-1, :].copy()

        # 2. 计算切片索引 (Slicing Indices)
        idx_eva = self.n_eva
        idx_aba = self.n_eva + self.n_aba
        idx_con = self.n_eva + self.n_aba + self.n_con  # 应该等于 ny

        # 3. 提取对应段的面积子数组
        # 如果 Y 向使用了非均匀网格，这里的 area_eva 等数组内部的元素值自然是不同的，无需额外处理
        area_eva = outer_surface_areas[0: idx_eva]
        area_aba = outer_surface_areas[idx_eva: idx_aba]
        area_con = outer_surface_areas[idx_aba: idx_con]

        # 4. 实例化虚拟切片的 BoundaryRegion，并挂载到 boundaries 字典
        self.boundaries['outer_eva'] = BoundaryRegion(shape=(self.n_eva,), area_array=area_eva)
        self.boundaries['outer_aba'] = BoundaryRegion(shape=(self.n_aba,), area_array=area_aba)
        self.boundaries['outer_con'] = BoundaryRegion(shape=(self.n_con,), area_array=area_con)

        # 5. 【极其关键】移除原始的全局 'right' 边界
        # 销毁该对象可以防止外部 Couplers 或基础类误调用一个横跨三段物理区域的统一边界
        if 'right' in self.boundaries:
            del self.boundaries['right']

    def _update_properties(self):
        """
        [重写父类方法] 更新物性
        逻辑：
        在这里处理管壁（Solid）和吸液芯（WickMaterial）的径向异构性。
        根据径向坐标或索引，对 k_node, rho_node, cp_node 数组进行分层赋值。
        （特别是要在此时计算并调用 pse1Cond, pse2Cond 等赝热导率）。
        """
        # 1. 获取全局温度场的 2D 视图 (nx, ny)
        T_2d = self.T.reshape(self.shape_nodes)

        # 创建用于存放 2D 物性的临时数组
        k_2d = np.zeros_like(T_2d)
        rho_2d = np.zeros_like(T_2d)
        cp_2d = np.zeros_like(T_2d)

        # 2. 吸液芯区域赋值 (径向索引: 0 到 n_wick-1)
        # 注意切片是视图，T_wick 的形状是 (n_wick, ny)
        T_wick = T_2d[:self.n_wick, :]
        k_2d[:self.n_wick, :] = self.wick_mat.conductivity(T_wick)
        rho_2d[:self.n_wick, :] = self.wick_mat.density(T_wick)
        cp_2d[:self.n_wick, :] = self.wick_mat.heat_capacity(T_wick)

        # 3. 管壁区域赋值 (径向索引: n_wick 到 nx-1)
        T_wall = T_2d[self.n_wick:, :]
        k_2d[self.n_wick:, :] = self.wall_mat.conductivity(T_wall)
        rho_2d[self.n_wick:, :] = self.wall_mat.density(T_wall)
        cp_2d[self.n_wick:, :] = self.wall_mat.heat_capacity(T_wall)

        # 4. 展平并覆盖基类的 1D 物性缓存数组
        self.k_node[:] = k_2d.flatten()
        self.rho_node[:] = rho_2d.flatten()
        self.cp_node[:] = cp_2d.flatten()

        # 5. 计算热容 (rho * cp * V)
        # 注意: mesh.geom_data.volumes 本身就是一个展开的 1D 数组 (长度 N)
        vols = self.mesh.geom_data.volumes
        self.thermal_capacitance[:] = self.rho_node * self.cp_node * vols

    def _update_boundaries_state(self, current_time: float = None):
        """
        [重写父类方法] 推送内部状态到边界
        逻辑：
        调用父类方法更新 'left', 'top', 'bottom'。
        对于外壁面，需要将 T 数组和 R_int 数组沿着轴向 (Y向) 切成三段，
        分别推送给 'outer_eva', 'outer_aba', 'outer_con'。
        """
        # 1. 获取全局状态的 2D 视图 (nx, ny)
        T_2d = self.T.reshape(self.shape_nodes)
        k_2d = self.k_node.reshape(self.shape_nodes)

        # 2. 获取底层的网格几何数据
        # ax_mat: 径向界面面积, dx_mat: 节点到界面的距离
        dx_mat = self.mesh.dx_matrix
        dy_mat = self.mesh.dy_matrix
        ax_mat = self.mesh.area_x_matrix
        ay_mat = self.mesh.area_y_matrix

        # ==========================================
        # A. 常规保留边界推送 (安全判断模式)
        # ==========================================

        # --- 1. Left (左侧内部边界 x=0) ---
        if 'left' in self.boundaries:
            d_left = dx_mat[0, :]
            A_left = ax_mat[0, :]
            k_left = k_2d[0, :]
            with np.errstate(divide='ignore', invalid='ignore'):
                R_int_left = d_left / (k_left * A_left)
                R_int_left = np.nan_to_num(R_int_left, posinf=1e10)
            self.boundaries['left'].update_internal_state(
                T_2d[0, :], R_int_left, current_time=current_time
            )

        # --- 2. Bottom (底部轴向边界 y=0) ---
        if 'bottom' in self.boundaries:
            d_bot = dy_mat[:, 0]
            A_bot = ay_mat[:, 0]
            k_bot = k_2d[:, 0]
            with np.errstate(divide='ignore', invalid='ignore'):
                R_int_bot = d_bot / (k_bot * A_bot)
                R_int_bot = np.nan_to_num(R_int_bot, posinf=1e10)
            self.boundaries['bottom'].update_internal_state(
                T_2d[:, 0], R_int_bot, current_time=current_time
            )

        # --- 3. Top (顶部轴向边界 y=ny-1) ---
        if 'top' in self.boundaries:
            d_top = dy_mat[:, -1]
            A_top = ay_mat[:, -1]
            k_top = k_2d[:, -1]
            with np.errstate(divide='ignore', invalid='ignore'):
                R_int_top = d_top / (k_top * A_top)
                R_int_top = np.nan_to_num(R_int_top, posinf=1e10)
            self.boundaries['top'].update_internal_state(
                T_2d[:, -1], R_int_top, current_time=current_time
            )

        # ==========================================
        # B. 核心切片：右侧外壁面边界推送
        # ==========================================

        # 4. 获取最外侧一圈网格节点 (x = nx-1) 的全量数据
        d_out = dx_mat[-1, :]
        A_out = ax_mat[-1, :]
        k_out = k_2d[-1, :]
        T_out = T_2d[-1, :]

        with np.errstate(divide='ignore', invalid='ignore'):
            # 计算外壁面的基础传导热阻
            R_int_out = d_out / (k_out * A_out)
            R_int_out = np.nan_to_num(R_int_out, nan=1e10, posinf=1e10)

        # 5. 准备轴向截断索引
        idx_eva = self.n_eva
        idx_aba = self.n_eva + self.n_aba
        idx_con = self.mesh.n_y  # 必然等于 n_eva + n_aba + n_con

        # 6. 分段推送给切片边界
        if 'outer_eva' in self.boundaries:
            self.boundaries['outer_eva'].update_internal_state(
                T_out[0: idx_eva],
                R_int_out[0: idx_eva],
                current_time=current_time
            )

        if 'outer_aba' in self.boundaries:
            self.boundaries['outer_aba'].update_internal_state(
                T_out[idx_eva: idx_aba],
                R_int_out[idx_eva: idx_aba],
                current_time=current_time
            )

        if 'outer_con' in self.boundaries:
            self.boundaries['outer_con'].update_internal_state(
                T_out[idx_aba: idx_con],
                R_int_out[idx_aba: idx_con],
                current_time=current_time
            )

    def _compute_fluxes(self, t: float) -> np.ndarray:
        """
        [重写父类方法] 计算净热流
        逻辑：
        先调用父类方法（或复用其内部通量计算代码）。
        但对于外边界的热流流入，不再从 'right' 获取。
        而是分别从 'outer_eva', 'outer_aba', 'outer_con' 获取切片热流，
        拼接成一个完整的轴向数组后，再累加到 Q_net_2d 的外壁面位置。
        """
        # 1. 准备 2D 累加器 (最终需要 flatten 返回)
        Q_net_2d = np.zeros(self.shape_nodes)
        T_2d = self.T.reshape(self.shape_nodes)

        # ==========================================
        # A. 内部通量 (Internal Fluxes)
        # ==========================================

        # --- X 方向 (径向): Flux i -> i+1 ---
        # 对应 Q_net 的: i 处流出(-), i+1 处流入(+)
        # 使用热导 G_x_inner (不是热阻 R_x_inner)
        flux_x = (T_2d[:-1, :] - T_2d[1:, :]) * self.G_x_inner

        Q_net_2d[:-1, :] -= flux_x  # 流出内侧节点
        Q_net_2d[1:, :] += flux_x  # 流入外侧节点

        # --- Y 方向 (轴向): Flux j -> j+1 ---
        # 使用热导 G_y_inner (不是热阻 R_y_inner)
        flux_y = (T_2d[:, :-1] - T_2d[:, 1:]) * self.G_y_inner

        Q_net_2d[:, :-1] -= flux_y  # 流出下方节点
        Q_net_2d[:, 1:] += flux_y  # 流入上方节点

        # ==========================================
        # B. 常规边界通量 (Boundary Fluxes)
        # ==========================================
        # BoundaryRegion 返回的是 "流入节点" 的热流 (Inflow to Node)，直接 += 即可

        if 'left' in self.boundaries:
            Q_net_2d[0, :] += self.boundaries['left'].compute_net_flux_for_solver()

        if 'bottom' in self.boundaries:
            Q_net_2d[:, 0] += self.boundaries['bottom'].compute_net_flux_for_solver()

        if 'top' in self.boundaries:
            Q_net_2d[:, -1] += self.boundaries['top'].compute_net_flux_for_solver()

        # ==========================================
        # C. 核心拼接：右侧外壁面边界热流 (Slicing Merge)
        # ==========================================

        idx_eva = self.n_eva
        idx_aba = self.n_eva + self.n_aba
        idx_con = self.mesh.n_y

        # 分别向最外层一圈节点 (x = nx-1) 注入对应轴向区段的热流
        if 'outer_eva' in self.boundaries:
            Q_net_2d[-1, 0: idx_eva] += self.boundaries['outer_eva'].compute_net_flux_for_solver()

        if 'outer_aba' in self.boundaries:
            Q_net_2d[-1, idx_eva: idx_aba] += self.boundaries['outer_aba'].compute_net_flux_for_solver()

        if 'outer_con' in self.boundaries:
            Q_net_2d[-1, idx_aba: idx_con] += self.boundaries['outer_con'].compute_net_flux_for_solver()

        # 最终展平为 1D 数组供 ODE 求解器步进使用
        return Q_net_2d.flatten()

    def get_boundary_node_capacitance(self, location: str) -> np.ndarray:
        """
        [重写父类方法] 获取边界热容 (用于 FluidSolidCouple 自适应步长)
        逻辑：
        扩展父类方法，当 location 为 'outer_eva' 等虚拟切片时，
        返回对应轴向切片范围内的节点热容。
        """
        cap_2d = self.thermal_capacitance.reshape(self.shape_nodes)

        idx_eva = self.n_eva
        idx_aba = self.n_eva + self.n_aba
        idx_con = self.mesh.n_y

        # 拦截切片位置请求并返回外壁面对应该段的热容
        if location == 'outer_eva':
            return cap_2d[-1, 0: idx_eva].copy()
        elif location == 'outer_aba':
            return cap_2d[-1, idx_eva: idx_aba].copy()
        elif location == 'outer_con':
            return cap_2d[-1, idx_aba: idx_con].copy()

        # 其它常规边界沿用父类逻辑
        return super().get_boundary_node_capacitance(location)
