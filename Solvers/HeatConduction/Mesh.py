import numpy as np
from dataclasses import dataclass
from typing import Optional, Union, Literal, Tuple

@dataclass
class GeometricData:
    """用于封装计算出的几何数据的容器"""
    volumes: np.ndarray             # 控制体体积 [m^3] (Size: N)
    face_areas: np.ndarray          # 界面面积 [m^2] (Size: N+1)
    node_centers: np.ndarray        # 节点中心坐标 [m] (Size: N)
    face_locations: np.ndarray      # 界面坐标 [m] (Size: N+1)
    dr_node_to_node: np.ndarray     # 相邻节点间的距离 [m] (Size: N-1)
    dr_node_to_face: np.ndarray     # 节点到自身界面的距离 [m] (Size: N, 2) -> [left, right]


class Mesh1D:
    """
    一维有限体积网格类
    支持: 笛卡尔坐标 (平板), 圆柱坐标 (燃料棒/管道)
    """

    def __init__(self,
                 total_dim: float,
                 n_volumes: int,
                 geometry_type: Literal['cartesian', 'cylindrical'] = 'cartesian',
                 inner_radius: float = 0.0,
                 height: float = 1.0):
        """
        标准均匀网格构造函数

        :param total_dim: 总长度(笛卡尔) 或 外径(圆柱) [m]
        :param n_volumes: 控制体数量 N
        :param geometry_type: 'cartesian' 或 'cylindrical'
        :param inner_radius: 内径 [m] (仅圆柱有效，实心圆柱为0，管道 >0)·
        :param height: 高度/深度 [m] (用于计算绝对体积和面积，默认为1.0)
        """
        self.geometry_type = geometry_type.lower()
        self.n_volumes = n_volumes
        self.L_or_R = total_dim
        self.r_in = inner_radius
        self.H = height
        self.N = n_volumes

        # 1. 生成界面坐标 (Face Locations) - 默认均匀划分
        if self.geometry_type == 'cartesian':
            # 平板：从 0 到 L
            self.face_locations = np.linspace(0, total_dim, n_volumes + 1)
        elif self.geometry_type == 'cylindrical':
            # 圆柱：从 r_in 到 r_out
            self.face_locations = np.linspace(inner_radius, total_dim, n_volumes + 1)
        else:
            raise ValueError(f"Unknown geometry type: {geometry_type}")

        # 2. 计算其余几何参数
        self._compute_geometry()

    @classmethod
    def from_custom_faces(cls, face_locations: np.ndarray,
                          geometry_type: str = 'cartesian',
                          height: float = 1.0):
        """
        [工厂方法] 从自定义界面坐标创建非均匀网格
        用于处理非均匀网格（如燃料棒边缘加密）
        """
        instance = cls.__new__(cls)
        instance.face_locations = np.array(face_locations, dtype=float)
        instance.n_volumes = len(face_locations) - 1
        instance.geometry_type = geometry_type.lower()
        instance.H = height

        # 自动推断边界
        instance.r_in = face_locations[0]
        instance.L_or_R = face_locations[-1]

        instance._compute_geometry()
        return instance

    def _compute_geometry(self):
        """核心几何计算逻辑"""
        faces = self.face_locations
        N = self.n_volumes

        # --- A. 节点中心 (Node Centers) ---
        # 简单算术平均 (对于圆柱坐标，FVM通常也使用几何中点作为节点代表温度点，
        # 尽管对数分布更精确，但算术平均在工程上通用性更好)
        self.node_centers = 0.5 * (faces[:-1] + faces[1:])

        # --- B. 几何属性计算 ---
        self.volumes = np.zeros(N)
        self.face_areas = np.zeros(N + 1)

        if self.geometry_type == 'cartesian':
            # --- 笛卡尔坐标 ---
            # 假设横截面积为 H (即深度 * 宽度1)
            # Volume = dx * Area
            dx = faces[1:] - faces[:-1]
            cross_section_area = self.H

            self.volumes = dx * cross_section_area
            self.face_areas[:] = cross_section_area  # 面积处处相等

        elif self.geometry_type == 'cylindrical':
            # --- 圆柱坐标 ---
            # Volume = pi * (r_out^2 - r_in^2) * H
            r_inner = faces[:-1]
            r_outer = faces[1:]
            self.volumes = np.pi * (r_outer ** 2 - r_inner ** 2) * self.H

            # Area = 2 * pi * r * H
            self.face_areas = 2.0 * np.pi * faces * self.H

        # --- C. 导热距离 (Distances for Resistance) ---
        # 1. 节点到节点距离 (用于计算内部界面热阻 R[1]...R[N-1])
        # size: N-1
        self.dr_node_to_node = self.node_centers[1:] - self.node_centers[:-1]

        # 2. 节点到自身边界的距离 (用于计算边界热阻或修正)
        # size: N (每个节点有两个半宽: [左半宽, 右半宽])
        # 这对于处理边界条件 R_internal = dr / (k*A) 非常重要
        self.dr_node_to_face = np.zeros((N, 2))
        self.dr_node_to_face[:, 0] = self.node_centers - faces[:-1]  # Node - LeftFace
        self.dr_node_to_face[:, 1] = faces[1:] - self.node_centers  # RightFace - Node

    def get_geometric_resistance_factor(self, conductivity_array: np.ndarray) -> np.ndarray:
        """
        [辅助方法] 计算几何热阻因子 (Geometry/Area 部分)
        R = geometry_factor / conductivity

        注意：这只是几何部分，实际热阻还需要除以导热系数 k。
        这里预留接口帮助 HeatConduction 类简化计算。

        返回数组大小: N+1 (对应所有界面)
        """
        R_geom = np.zeros(self.n_volumes + 1)

        # 这里需要复杂的逻辑来处理 k 的调和平均位置
        # 为了保持 Mesh 类的纯粹性，我们只提供 dr/Area 的几何比率
        # 具体的 k 组合逻辑留给 HeatConduction 类

        return R_geom

    def __repr__(self):
        return (f"<Mesh1D type={self.geometry_type} N={self.n_volumes} "
                f"range=[{self.face_locations[0]:.4f}, {self.face_locations[-1]:.4f}]>")

@dataclass
class GeometricData2D:
    """
        二维网格几何数据容器
        为了求解器效率，核心数据存储为 1D 数组，但逻辑上代表 2D 场。
        """
    # --- 节点属性 (Size: N = Nx * Ny) ---
    volumes: np.ndarray  # 控制体体积 [m^3]
    node_centers_x: np.ndarray  # 节点中心 X/R 坐标 [m]
    node_centers_y: np.ndarray  # 节点中心 Y/Z 坐标 [m]

    # --- 界面属性 (用于计算 Flux) ---
    # Area X: 垂直于 X/R 轴的界面面积。
    # 逻辑 Shape: (Nx+1, Ny) -> Flattened Size: (Nx+1)*Ny
    area_x: np.ndarray

    # Area Y: 垂直于 Y/Z 轴的界面面积。
    # 逻辑 Shape: (Nx, Ny+1) -> Flattened Size: Nx*(Ny+1)
    area_y: np.ndarray

    # --- 导热距离 (Distances) ---
    # 用于计算梯度 dT/dx, dT/dy
    # 逻辑 Shape 与 Area 相同，表示跨过界面的两个节点间距离
    dx_node_to_node: np.ndarray  # 跨 X 界面的距离 (Nx+1, Ny)
    dy_node_to_node: np.ndarray  # 跨 Y 界面的距离 (Nx, Ny+1)

    # 节点到自身界面的距离 (用于边界处理)
    # Shape: (N, 4) -> [x_left, x_right, y_bottom, y_top]
    # 或者为了内存紧凑，可以拆分为 dx_node_to_face (N, 2) 和 dy...
    # 这里为了通用性，存储跨越界面的总距离通常更直接。
    # 但为了边界条件 (Boundary Condition)，我们需要节点中心到边界面的距离。
    dist_to_boundaries: np.ndarray  # Shape (N, 4): [x_in, x_out, y_in, y_out]


class Mesh2D:
    """
    二维有限体积网格类
    支持: 笛卡尔坐标 (x, y) 和 圆柱坐标 (r, z)
    存储策略: 内部数据尽量扁平化(1D)以优化求解器交互，但提供2D视图逻辑。
    """
    def __init__(self,
                 x_dim: float, n_x: int,
                 y_dim: float, n_y: int,
                 geometry_type: Literal['cartesian', 'cylindrical'] = 'cartesian',
                 inner_radius: float = 0.0,
                 y_faces: Optional[np.ndarray] = None,
                 x_faces: Optional[np.ndarray] = None):
        """
        :param x_dim: X方向长度 或 R方向外径(若inner_radius=0) 或 径向厚度
                      注意：如果是圆柱坐标，这代表 (R_outer - R_inner)
        :param n_x:   X/R 方向网格数
        :param y_dim: Y方向长度 或 Z方向高度
        :param n_y:   Y/Z 方向网格数
        :param geometry_type: 'cartesian' or 'cylindrical'
        :param inner_radius: 内径 (仅圆柱坐标有效)
        :param y_faces: [新增] 自定义轴向界面坐标数组。若提供，则生成轴向非均匀网格。
        """

        self.geometry_type = geometry_type

        # --- [新增] 处理自定义非均匀网格逻辑 (X/径向方向) ---
        if x_faces is not None:
            self.x_faces_custom = np.asarray(x_faces, dtype=float)
            self.n_x = len(self.x_faces_custom) - 1
            # 自动推断径向总厚度
            self.x_dim = self.x_faces_custom[-1] - self.x_faces_custom[0]

            if self.geometry_type == 'cylindrical':
                # 在圆柱坐标下，非均匀网格的第一个面坐标就是内径
                self.inner_radius = self.x_faces_custom[0]
            else:
                self.inner_radius = inner_radius
        else:
            self.x_faces_custom = None
            self.x_dim = x_dim
            self.n_x = n_x
            self.inner_radius = inner_radius

        # --- 处理自定义非均匀网格逻辑 (Y/轴向方向) ---
        if y_faces is not None:
            self.y_faces_custom = np.asarray(y_faces, dtype=float)
            self.n_y = len(self.y_faces_custom) - 1
            # 自动推断轴向总长度
            self.y_dim = self.y_faces_custom[-1] - self.y_faces_custom[0]
        else:
            self.y_faces_custom = None
            self.y_dim = y_dim
            self.n_y = n_y

        # 校验维度
        if self.n_y <= 0 or self.n_x <= 0:
            raise ValueError("Mesh2D requires n_x > 0 and n_y > 0. Provide valid dims and n, or faces arrays.")

        self.N = self.n_x * self.n_y
        self.n_volumes = self.N

        # 生成网格坐标
        self._generate_grid()

        # 计算几何参数 (体积, 面积, 距离)
        self.geom_data = self._compute_geometry()

        # 缓存常用的形状视图，避免重复 reshape
        self.shape_nodes = (self.n_x, self.n_y)
        self.shape_faces_x = (self.n_x + 1, self.n_y)
        self.shape_faces_y = (self.n_x, self.n_y + 1)

        # 下面是旧版网格文件生成，无法生成非均匀轴向网格
        # 为适应热管轴向网格划分，需要增加非均匀轴向网格

        # self.x_dim = x_dim
        # self.n_x = n_x
        # self.y_dim = y_dim
        # self.n_y = n_y
        # self.N = n_x * n_y
        # self.n_volumes = self.N
        # self.geometry_type = geometry_type
        # self.inner_radius = inner_radius
        #
        # # 生成网格坐标
        # self._generate_grid()
        #
        # # 计算几何参数 (体积, 面积, 距离)
        # self.geom_data = self._compute_geometry()
        #
        # # 缓存常用的形状视图，避免重复 reshape
        # self.shape_nodes = (n_x, n_y)
        # self.shape_faces_x = (n_x + 1, n_y)
        # self.shape_faces_y = (n_x, n_y + 1)

    def _generate_grid(self):
        """生成界面和节点坐标"""
        # --- Dimension 1 (X or R) ---
        # [修改] 支持非均匀径向网格
        if getattr(self, 'x_faces_custom', None) is not None:
            self.x_faces = self.x_faces_custom
        else:
            # 物理坐标范围
            r_start = 0.0 if self.geometry_type == 'cartesian' else self.inner_radius
            r_end = self.x_dim if self.geometry_type == 'cartesian' else (self.inner_radius + self.x_dim)
            self.x_faces = np.linspace(r_start, r_end, self.n_x + 1)

        self.x_centers = 0.5 * (self.x_faces[:-1] + self.x_faces[1:])
        self.dx = self.x_faces[1:] - self.x_faces[:-1]  # 向量

        # --- Dimension 2 (Y or Z) ---
        # 如果传入了自定义网格面，则直接使用；否则生成均匀网格
        if getattr(self, 'y_faces_custom', None) is not None:
            self.y_faces = self.y_faces_custom
        else:
            self.y_faces = np.linspace(0.0, self.y_dim, self.n_y + 1)

        self.y_centers = 0.5 * (self.y_faces[:-1] + self.y_faces[1:])
        self.dy = self.y_faces[1:] - self.y_faces[:-1]

        # 下面是旧版网格文件生成，无法生成非均匀轴向网格
        # 为适应热管轴向网格划分，需要增加非均匀轴向网格

        # # --- Dimension 1 (X or R) ---
        # # 物理坐标范围
        # r_start = 0.0 if self.geometry_type == 'cartesian' else self.inner_radius
        # r_end = self.x_dim if self.geometry_type == 'cartesian' else (self.inner_radius + self.x_dim)
        #
        # self.x_faces = np.linspace(r_start, r_end, self.n_x + 1)
        # self.x_centers = 0.5 * (self.x_faces[:-1] + self.x_faces[1:])
        # self.dx = self.x_faces[1:] - self.x_faces[:-1]  # 向量，虽然均匀网格是常数
        #
        # # --- Dimension 2 (Y or Z) ---
        # self.y_faces = np.linspace(0.0, self.y_dim, self.n_y + 1)
        # self.y_centers = 0.5 * (self.y_faces[:-1] + self.y_faces[1:])
        # self.dy = self.y_faces[1:] - self.y_faces[:-1]

    def _compute_geometry(self) -> GeometricData2D:
        """
        计算所有几何参数。
        注意：尽量利用广播机制 (Broadcasting) 计算，然后 Flatten 存储。
        """
        nx, ny = self.n_x, self.n_y

        # ==========================================
        # 1. 体积 (Volumes)
        # ==========================================
        # 笛卡尔: V = dx * dy (假设单位深度)
        # 圆柱:   V = pi * (r_out^2 - r_in^2) * dz

        # 准备广播: X在行(dim 0), Y在列(dim 1) -> shape (nx, ny)
        dx_2d = self.dx[:, np.newaxis]  # shape (nx, 1)
        dy_2d = self.dy[np.newaxis, :]  # shape (1, ny)

        if self.geometry_type == 'cartesian':
            # V_ij = dx_i * dy_j
            volumes_2d = dx_2d * dy_2d

            # X-Centers (2D array)
            xc_2d = self.x_centers[:, np.newaxis] * np.ones((1, ny))
            yc_2d = np.ones((nx, 1)) * self.y_centers[np.newaxis, :]

        else:  # cylindrical
            # V_ij = Area_ring_i * dz_j
            # Area_ring_i = pi * (r_out^2 - r_in^2)
            r_in = self.x_faces[:-1]
            r_out = self.x_faces[1:]
            area_ring = np.pi * (r_out ** 2 - r_in ** 2)  # shape (nx,)

            volumes_2d = area_ring[:, np.newaxis] * dy_2d  # (nx, ny)

            xc_2d = self.x_centers[:, np.newaxis] * np.ones((1, ny))  # r centers
            yc_2d = np.ones((nx, 1)) * self.y_centers[np.newaxis, :]  # z centers

        # ==========================================
        # 2. 界面面积 (Face Areas)
        # ==========================================

        # --- Area X (垂直于X/R轴的面) ---
        # Shape (nx+1, ny)
        # Cartesian: dy * 1
        # Cylindrical: 2 * pi * r_face * dz

        if self.geometry_type == 'cartesian':
            # 所有 i 的面积都是 dy_j
            # shape (nx+1, ny)
            ax_2d = np.ones((nx + 1, 1)) * dy_2d
        else:
            # Area = 2 * pi * r * dz
            # r 来自 x_faces (size nx+1)
            r_faces_col = self.x_faces[:, np.newaxis]  # (nx+1, 1)
            ax_2d = 2.0 * np.pi * r_faces_col * dy_2d  # (nx+1, ny)

        # --- Area Y (垂直于Y/Z轴的面) ---
        # Shape (nx, ny+1)
        # Cartesian: dx * 1
        # Cylindrical: pi * (r_out^2 - r_in^2) (即圆环面积 area_ring)

        if self.geometry_type == 'cartesian':
            # 所有 j 的面积都是 dx_i
            ay_2d = dx_2d * np.ones((1, ny + 1))
        else:
            # 这里的 AreaY 是 Z方向的流通面积，即径向的圆环截面
            # 它不随 Z (dim 1) 变化
            # area_ring shape is (nx,)
            r_in = self.x_faces[:-1]
            r_out = self.x_faces[1:]
            area_ring_y = np.pi * (r_out ** 2 - r_in ** 2)  # shape (nx,)

            ay_2d = area_ring_y[:, np.newaxis] * np.ones((1, ny + 1))

        # ==========================================
        # 3. 导热距离 (Distances)
        # ==========================================

        # --- Node to Node Distances ---
        # 用于计算 Flux = (T1-T2) / (dist/kA)
        # 这里的 dist 是跨越界面的两个节点中心的距离

        # X-direction (nx+1, ny)
        # 内部: center[i] - center[i-1]
        # 边界: center[0] - face[0] (左), face[end] - center[end] (右)
        # 为了通用，我们先构建完整的 (nx+1) 距离数组

        # 计算 X 方向跨界面的距离
        # 内部 (indices 1 to nx-1):
        d_nodes_x = self.x_centers[1:] - self.x_centers[:-1]
        # 边界 (index 0 and nx):
        d_bound_left = self.x_centers[0] - self.x_faces[0]
        d_bound_right = self.x_faces[-1] - self.x_centers[-1]

        # 拼接成 (nx+1,)
        dist_x_1d = np.concatenate(([d_bound_left], d_nodes_x, [d_bound_right]))
        dx_node_to_node_2d = dist_x_1d[:, np.newaxis] * np.ones((1, ny))

        # Y-direction (nx, ny+1)
        d_nodes_y = self.y_centers[1:] - self.y_centers[:-1]
        d_bound_bottom = self.y_centers[0] - self.y_faces[0]
        d_bound_top = self.y_faces[-1] - self.y_centers[-1]

        dist_y_1d = np.concatenate(([d_bound_bottom], d_nodes_y, [d_bound_top]))
        dy_node_to_node_2d = np.ones((nx, 1)) * dist_y_1d[np.newaxis, :]

        # --- Node to Boundary Distances (Store specifically for BCs) ---
        # Shape (nx, ny, 4) -> flatten to (N, 4)
        # order: x_inner, x_outer, y_inner, y_outer
        # 这是一个简单的几何距离，用于计算节点到表面的热阻

        dist_to_bounds = np.zeros((nx, ny, 4))
        # x_inner (left) distance: center_x - face_x_left
        dist_to_bounds[:, :, 0] = (self.x_centers - self.x_faces[:-1])[:, np.newaxis]
        # x_outer (right) distance: face_x_right - center_x
        dist_to_bounds[:, :, 1] = (self.x_faces[1:] - self.x_centers)[:, np.newaxis]
        # y_inner (bottom)
        dist_to_bounds[:, :, 2] = (self.y_centers - self.y_faces[:-1])[np.newaxis, :]
        # y_outer (top)
        dist_to_bounds[:, :, 3] = (self.y_faces[1:] - self.y_centers)[np.newaxis, :]

        # ==========================================
        # 4. Flatten and Return
        # ==========================================
        return GeometricData2D(
            volumes=volumes_2d.flatten(),
            # C-style (row-major): x changes last? No, numpy default is last index changes fastest.
            # Wait, numpy flatten 'C': last axis changes fastest.
            # Our arrays are (nx, ny). So 'ny' (y-axis) changes fastest.
            # k = i * ny + j. Correct.
            node_centers_x=xc_2d.flatten(),
            node_centers_y=yc_2d.flatten(),
            area_x=ax_2d.flatten(),
            area_y=ay_2d.flatten(),
            dx_node_to_node=dx_node_to_node_2d.flatten(),
            dy_node_to_node=dy_node_to_node_2d.flatten(),
            dist_to_boundaries=dist_to_bounds.reshape(-1, 4)
        )

    # ==========================================
    # 辅助方法: 索引转换
    # ==========================================

    def flatten_index(self, i: int, j: int) -> int:
        """(i, j) -> k.  i: x/r index, j: y/z index"""
        return i * self.n_y + j

    def unravel_index(self, k: int) -> Tuple[int, int]:
        """k -> (i, j)"""
        return k // self.n_y, k % self.n_y

    # ==========================================
    # 视图获取器 (View Getters)
    # 允许以 (Nx, Ny) 的形状访问数据，无内存拷贝
    # ==========================================

    @property
    def volumes_matrix(self):
        return self.geom_data.volumes.reshape(self.shape_nodes)

    @property
    def area_x_matrix(self):
        return self.geom_data.area_x.reshape(self.shape_faces_x)

    @property
    def area_y_matrix(self):
        return self.geom_data.area_y.reshape(self.shape_faces_y)

    @property
    def dx_matrix(self):
        return self.geom_data.dx_node_to_node.reshape(self.shape_faces_x)

    @property
    def dy_matrix(self):
        return self.geom_data.dy_node_to_node.reshape(self.shape_faces_y)

