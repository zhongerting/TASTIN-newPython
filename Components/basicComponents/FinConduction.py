import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from typing import Union

from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Base import SolidMaterial
from profiler import TEASAProfiler


class FinConduction(HeatConduction2D):
    """
    准稳态二维翅片导热求解器 (Quasi-Steady Fin Conduction)

    物理假设与降维策略:
    1. 准稳态 (Quasi-Steady): 忽略翅片极小的热容，每个全局时间步内瞬时达到热平衡 (dT/dt = 0)。
    2. 一维传热条 (1D Strips): 忽略翅片沿热管轴向 (Y向) 的导热，将 2D 网格退化为 N_y 个独立的 1D 方程组。
    3. 表面辐射体积化: 翅片暴露在正反两侧巨大的散热面积被等效为体积热源/热漏。边缘散热则由标准边界条件处理。
    """

    def __init__(self,
                 mesh: Mesh2D,
                 material: SolidMaterial,
                 fin_thickness: float,
                 emissivity: float = 0.8,
                 up_view_factor: float = 1.0,  # [新增] 上表面角系数
                 down_view_factor: float = 1.0,  # [新增] 下表面角系数
                 T_env: Union[float, np.ndarray] = 3.0,
                 initial_temp: float = 298.15):

        super().__init__(mesh, material, initial_temp)

        self.fin_thickness = fin_thickness
        self.emissivity = emissivity
        self.up_view_factor = up_view_factor  # [新增] 保存上角系数
        self.down_view_factor = down_view_factor  # [新增] 保存下角系数
        self.T_env = T_env
        self.sigma = 5.670374419e-8  # 斯蒂芬-玻尔兹曼常数

    def _get_boundary_linear_terms(self, boundary_key: str):
        """
        [内部辅助] 从边界区域提取戴维南等效线性系数
        将复杂的外部边界 (Couplers 挂载的热阻和热流) 转换为稀疏矩阵所需的 G_bound 和 S_bound
        Flux = G_bound * T_eff - G_bound * T_node + Q_flux_only
        """
        boundary = self.boundaries[boundary_key]
        shape = boundary.shape

        G_bound = np.zeros(shape)
        S_bound = np.zeros(shape)

        # 重建诺顿等效源
        G_sum = np.zeros(shape)
        J_sum = np.zeros(shape)
        Q_flux = np.zeros(shape)

        for bc in boundary.conditions:
            if bc.bc_type == "resistance":
                R_ext_total = bc.R_ext + bc.R_add
                with np.errstate(divide='ignore', invalid='ignore'):
                    G_i = 1.0 / R_ext_total
                    G_i = np.nan_to_num(G_i, posinf=0.0)
                G_sum += G_i
                J_sum += bc.T_ext * G_i
            elif bc.bc_type == "flux":
                Q_flux += bc.q_flux

        # 处理带有热阻连接的节点
        has_res = G_sum > 1e-20
        if np.any(has_res):
            R_eq = 1.0 / G_sum[has_res]
            # 戴维南等效温度: T_eff = T_norton + Q_flux * R_eq
            T_eff = J_sum[has_res] * R_eq + Q_flux[has_res] * R_eq
            # 串联上内部热阻
            R_tot = R_eq + boundary.R_internal[has_res]
            G_bound[has_res] = 1.0 / R_tot
            S_bound[has_res] = G_bound[has_res] * T_eff

        # 处理纯热流(Neumann)节点
        S_bound[~has_res] += Q_flux[~has_res]

        return G_bound, S_bound

    @TEASAProfiler.profile
    def step(self, dt: float, max_iter: int = 50, tol: float = 1e-4, relaxation: float = 0.8, **kwargs) -> bool:
        """
        [核心重写] 拦截基类的 ODE 求解器，转而执行非线性稳态 Picard 迭代。
        """
        nx, ny = self.shape_nodes
        N = self.N

        err = 1.0e6

        for k_iter in range(max_iter):
            T_old = self.T.copy()

            # 1. 更新当前温度下的物性和几何热阻
            self._update_properties()
            self._compute_internal_resistance()
            self._update_boundaries_state()

            # === 开始组装稀疏矩阵 A * T_new = b ===
            # 使用对角线 D (主对角) 和 off-diagonals (上下对角)
            D = np.zeros(N)
            b = np.zeros(N)

            # 2. 内部 X方向导热 (径向高度)
            # 获取调和平均后的热导 G_x (对应 R_x_inner)
            G_x = 1.0 / self.R_x_inner  # shape: (nx-1, ny)
            G_x_flat = G_x.flatten()

            # 计算相邻节点的平铺索引 (Flattened Index: k = i * ny + j)
            i_idx, j_idx = np.meshgrid(np.arange(nx - 1), np.arange(ny), indexing='ij')
            k1_flat = (i_idx * ny + j_idx).flatten()
            k2_flat = ((i_idx + 1) * ny + j_idx).flatten()

            # 填充主对角线 (流出项为正)
            np.add.at(D, k1_flat, G_x_flat)
            np.add.at(D, k2_flat, G_x_flat)

            # 填充非对角线 (流入项为负)
            rows = np.concatenate([k1_flat, k2_flat])
            cols = np.concatenate([k2_flat, k1_flat])
            vals = np.concatenate([-G_x_flat, -G_x_flat])

            # [关键降维]: 不添加 Y 方向的导热，矩阵自动解耦为 N_y 个独立的块！

            # 3. 翅片正反面辐射散热 (Surface Radiation)
            # 局部线性化 (Picard): h_rad = eps * sigma * (T^2 + T_env^2)*(T + T_env)
            h_rad = self.emissivity * self.sigma * (T_old ** 2 + self.T_env ** 2) * (T_old + self.T_env)

            # [修改] 结合上下角系数计算有效辐射面积
            # 单面暴露面积 = volume / thickness
            # 真实有效辐射面积 = (F_up + F_down) * 单面面积
            total_view_factor = self.up_view_factor + self.down_view_factor
            A_rad_eff = total_view_factor * (self.mesh.geom_data.volumes / self.fin_thickness)

            G_rad = h_rad * A_rad_eff

            D += G_rad
            b += G_rad * self.T_env

            # 4. 处理四条边界条件
            # Left (x=0) -> 翅片根部，对接热管冷凝段
            G_l, S_l = self._get_boundary_linear_terms('left')
            idx_left = np.arange(ny)  # i=0
            D[idx_left] += G_l
            b[idx_left] += S_l

            # Right (x=nx-1) -> 翅片远端边缘
            G_r, S_r = self._get_boundary_linear_terms('right')
            idx_right = (nx - 1) * ny + np.arange(ny)
            D[idx_right] += G_r
            b[idx_right] += S_r

            # Bottom (y=0) & Top (y=ny-1) -> 轴向两端边缘
            G_b, S_b = self._get_boundary_linear_terms('bottom')
            idx_bottom = np.arange(nx) * ny
            D[idx_bottom] += G_b
            b[idx_bottom] += S_b

            G_t, S_t = self._get_boundary_linear_terms('top')
            idx_top = np.arange(nx) * ny + ny - 1
            D[idx_top] += G_t
            b[idx_top] += S_t

            # 5. 添加外部体积热源 (如有)
            self._update_sources(self.current_time)
            b += self.Q_source

            # 6. 生成完整的稀疏矩阵并求解
            # 加入主对角线
            rows = np.concatenate([rows, np.arange(N)])
            cols = np.concatenate([cols, np.arange(N)])
            vals = np.concatenate([vals, D])

            A_sparse = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()

            # 使用直接求解器 (针对块三对角矩阵极其迅速)
            T_new = spla.spsolve(A_sparse, b)

            # 7. 收敛检查与松弛更新 (防止辐射引起的非线性震荡)
            err = np.max(np.abs(T_new - T_old))
            self.T = relaxation * T_new + (1 - relaxation) * T_old

            if err < tol:
                break

        if err >= tol:
            print(f"Warning: FinConduction step did not fully converge. Max Error: {err:.4e}")

        # 更新时间戳 (虽然是稳态，但维持系统时钟同步)
        self.current_time += dt
        return True
