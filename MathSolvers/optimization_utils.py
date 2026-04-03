import numpy as np
from scipy.sparse import coo_matrix
from typing import List, Dict

# 引入类型 (假设在路径中)
from Solvers.Hydrodynamics.Components import FluidVolume, FlowJunction


class FluidJacobianBlockLayout:
    """
    适配 SystemManager 块状布局的雅可比矩阵构建器

    状态向量 y 的布局:
    [ P_0, ..., P_N-1,  h_0, ..., h_N-1,  W_0, ..., W_M-1 ]
    """

    def __init__(self, volumes: List[FluidVolume], junctions: List[FlowJunction]):
        self.volumes = volumes
        self.junctions = junctions

        self.n_vol = len(volumes)
        self.n_junc = len(junctions)
        self.n_vars = 2 * self.n_vol + self.n_junc

        # 映射
        self.vol_map: Dict[FluidVolume, int] = {v: i for i, v in enumerate(volumes)}
        self.junc_map: Dict[FlowJunction, int] = {j: i for i, j in enumerate(junctions)}

    # --- 块状布局索引计算 (核心修改) ---

    def _idx_P(self, vol_idx: int) -> int:
        # P 位于向量的最前面 [0, n_vol)
        return vol_idx

    def _idx_h(self, vol_idx: int) -> int:
        # h 位于 P 之后 [n_vol, 2*n_vol)
        return self.n_vol + vol_idx

    def _idx_W(self, junc_idx: int) -> int:
        # W 位于 h 之后 [2*n_vol, end)
        return 2 * self.n_vol + junc_idx

    def get_sparsity_matrix(self) -> coo_matrix:
        rows = []
        cols = []

        def add(r, c):
            rows.append(r)
            cols.append(c)

        # 1. 遍历 Volume 方程 (P 和 h)
        for i, vol in enumerate(self.volumes):
            row_P = self._idx_P(i)
            row_h = self._idx_h(i)

            # A. 自相关 (d/dt 依赖于自身的 P, h)
            # dP/dt <- P, h
            add(row_P, self._idx_P(i));
            add(row_P, self._idx_h(i))
            # dh/dt <- P, h
            add(row_h, self._idx_P(i));
            add(row_h, self._idx_h(i))

            # B. 流量耦合 (通过 Junction)
            connected = vol.inlet_junctions + vol.outlet_junctions
            for junc in connected:
                if junc not in self.junc_map: continue

                j_idx = self.junc_map[junc]
                col_W = self._idx_W(j_idx)

                # 质量/能量方程依赖于流量 W
                add(row_P, col_W)
                add(row_h, col_W)

                # C. 施主元胞 (Donor Cell) - 能量平流
                # 依赖于邻居的 h (可能还有 P)
                neighbor = junc.from_vol if junc.to_vol == vol else junc.to_vol
                if neighbor in self.vol_map:
                    n_idx = self.vol_map[neighbor]
                    add(row_h, self._idx_h(n_idx))  # 主要是 h

        # 2. 遍历 Junction 方程 (W)
        for j, junc in enumerate(self.junctions):
            row_W = self._idx_W(j)

            # A. 自相关 (摩擦力)
            add(row_W, self._idx_W(j))

            # B. 压差驱动 (上下游 P)
            # 上游
            if junc.from_vol in self.vol_map:
                f_idx = self.vol_map[junc.from_vol]
                add(row_W, self._idx_P(f_idx))  # dW/dt <- P_in
                add(row_W, self._idx_h(f_idx))  # dW/dt <- h_in (密度)

            # 下游
            if junc.to_vol in self.vol_map:
                t_idx = self.vol_map[junc.to_vol]
                add(row_W, self._idx_P(t_idx))  # dW/dt <- P_out
                add(row_W, self._idx_h(t_idx))  # dW/dt <- h_out (密度)

        data = np.ones(len(rows), dtype=int)
        return coo_matrix((data, (rows, cols)), shape=(self.n_vars, self.n_vars))