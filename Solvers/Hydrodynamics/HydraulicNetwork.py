import math
import numpy as np
import logging
from typing import List, Any, Tuple, Dict, Set
from scipy.sparse import lil_matrix, csr_matrix, coo_matrix
from scipy.sparse.linalg import spsolve

from profiler import TEASAProfiler

from Solvers.Hydrodynamics.Components import FluidVolume, FlowJunction
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume

# 配置日志
logger = logging.getLogger(__name__)


class _FrozenBoundaryIndexSet(set):
    """带运行时保护的定压边界索引集合。"""

    def __init__(self, owner: "HydraulicNetwork"):
        super().__init__()
        self._owner = owner

    def _raise_if_locked(self):
        if self._owner._fixed_pressure_indices_locked:
            raise RuntimeError(
                "Runtime modification of fixed_pressure_indices is not supported. "
                "The cached matrix structures assume a fixed pressure-boundary set. "
                "Please rebuild the HydraulicNetwork if the fixed-pressure boundary membership must change."
            )

    def add(self, element):
        self._raise_if_locked()
        return super().add(element)

    def clear(self):
        self._raise_if_locked()
        return super().clear()

    def discard(self, element):
        self._raise_if_locked()
        return super().discard(element)

    def difference_update(self, *others):
        self._raise_if_locked()
        return super().difference_update(*others)

    def intersection_update(self, *others):
        self._raise_if_locked()
        return super().intersection_update(*others)

    def pop(self):
        self._raise_if_locked()
        return super().pop()

    def remove(self, element):
        self._raise_if_locked()
        return super().remove(element)

    def symmetric_difference_update(self, other):
        self._raise_if_locked()
        return super().symmetric_difference_update(other)

    def update(self, *others):
        self._raise_if_locked()
        return super().update(*others)

    def __iand__(self, other):
        self._raise_if_locked()
        return super().__iand__(other)

    def __ior__(self, other):
        self._raise_if_locked()
        return super().__ior__(other)

    def __isub__(self, other):
        self._raise_if_locked()
        return super().__isub__(other)

    def __ixor__(self, other):
        self._raise_if_locked()
        return super().__ixor__(other)


class HydraulicNetwork:
    """
    半隐式液力网络求解器 (Semi-Implicit Hydraulic Network Solver) - Phase 1 & 2 Skeleton

    核心职责：
    1. [Phase 1] 拓扑映射：将 FluidVolume 和 FlowJunction 映射为线性方程组的索引。
    2. [Phase 2] 系数计算：计算动量方程的导纳 (Admittance) 和源项 (Source)。
    3. [Phase 3] 矩阵求解：组装并求解 A*P = B。
    4. [Phase 4] 能量耦合：显式推进焓/温度方程。
    """

    def __init__(self, volumes: List[Any], junctions: List[Any], gravity_vector: float = 9.81):
        """
        初始化网络求解器

        :param volumes: 系统中所有的控制体列表 (包含 Channel 内部节点和 BoundaryPlenum)
        :param junctions: 系统中所有的连接列表 (包含 Inlet/Outlet 和 Internal Junctions)
        :param gravity_vector: 重力加速度 [m/s^2]。
                               注意：我们假设 z 轴垂直向上为正。
                               因此重力总是指向 z 轴负方向。
                               默认值为 9.81 (标量)，计算时会自动处理符号。
        """
        self.volumes_obj = volumes
        self.junctions_obj = junctions
        self.g = gravity_vector  # [Fix 3] 初始化重力加速度l.z_coordinate

        # 矩阵维度
        self.n_vol = len(volumes)
        self.n_junc = len(junctions)

        # --- Phase 1: 拓扑数据 ---
        # 1. 映射表 (Object -> Matrix Index)
        self.vol_to_idx: Dict[Any, int] = {}

        # 2. 邻接表 (Connectivity)
        # 存储格式: List of (junction_obj, idx_from, idx_to)
        self.junction_descriptors: List[Tuple[Any, int, int]] = []

        # 3. 边界条件索引 (Fixed Pressure Boundary Indices)
        # 存储那些被标记为“稳压器”的节点的索引
        self._fixed_pressure_indices_locked = False
        self._fixed_pressure_indices = _FrozenBoundaryIndexSet(self)
        self.fixed_pressure_mask = np.zeros(self.n_vol, dtype=bool)
        self.normal_pressure_mask = np.ones(self.n_vol, dtype=bool)
        self.fixed_pressure_idx_arr = np.empty(0, dtype=np.int32)
        self.normal_pressure_idx_arr = np.empty(0, dtype=np.int32)
        self.fixed_target_P_vec = np.empty(0)

        # --- Phase 1: 几何常数 (Geometric Constants) ---
        # 存储节点几何常数，防止计算过程中反复调用
        self.V_vec = np.zeros(self.n_vol)  # vol.vol
        self.L_node_vec = np.zeros(self.n_vol)  # vol.len
        self.A_node_vec = np.zeros(self.n_vol)  # vol.area
        self.Dh_node_vec = np.zeros(self.n_vol)  # vol.d_h
        self.z_vec = np.zeros(self.n_vol)  # vol.z_coordinate
        # --- 连接几何常数 ---
        self.A_junc_vec = np.zeros(self.n_junc)  # junc.area
        self.L_junc_vec = np.zeros(self.n_junc)  # junc.length
        self.K_loss_vec = np.zeros(self.n_junc)  # junc.k_loss
        self.idx_from_vec = np.zeros(self.n_junc, dtype=np.int32)
        self.idx_to_vec = np.zeros(self.n_junc, dtype=np.int32)
        self.Dh_from_vec = np.zeros(self.n_junc)
        self.Dh_to_vec = np.zeros(self.n_junc)
        self.L_from_half_vec = np.zeros(self.n_junc)
        self.L_to_half_vec = np.zeros(self.n_junc)
        self.A_from_node_vec = np.zeros(self.n_junc)
        self.A_to_node_vec = np.zeros(self.n_junc)
        self.is_inlet_junction_mask = np.zeros(self.n_junc, dtype=bool)
        self.target_W_vec = np.zeros(self.n_junc)
        self.inlet_junction_indices = np.empty(0, dtype=np.int32)
        # --- 乘子向量，用于网络拓扑的非对称质量/能量缩放 ---
        self.M_from_vec = np.ones(self.n_junc)
        self.M_to_vec = np.ones(self.n_junc)

        # --- Phase 2: 物性缓存数组 (Properties) ---
        # 存储节点几何常数，防止计算过程中反复调用
        self.cp_vec = np.ones(self.n_vol) * 1000.0  # 热容
        self.drho_dp_vec = np.zeros(self.n_vol)  # 压缩性
        self.drho_dt_vec = np.zeros(self.n_vol)  # 热膨胀率
        self.Q_expl_vec = np.zeros(self.n_vol)  # 显式热源 (Q_wall + Q_vol)
        self.lam_imp_vec = np.zeros(self.n_vol)  # 隐式换热系数 (implicit_coeff)

        # --- Phase 3: 状态向量 (State Vectors) ---
        # 使用 numpy 数组存储，方便矩阵运算
        self.P_vec = np.zeros(self.n_vol)  # Pressure [Pa]
        self.T_vec = np.zeros(self.n_vol)  # Temperature [K]
        self.h_vec = np.zeros(self.n_vol)  # Enthalpy [J/kg]
        self.W_vec = np.zeros(self.n_junc)  # Mass Flow [kg/s]
        # 用于迭代的流量变量
        self.W_old = np.zeros(self.n_junc)
        self.W_iterate = np.zeros(self.n_junc)
        self.W_residual = np.zeros(self.n_junc)

        # --- Phase 4: 物理系数 (Physics Coefficients) ---
        # 动量方程系数: W^{n+1} = a_j * (P_in - P_out) + b_j
        self.A_coeffs = np.zeros(self.n_junc)  # 导纳 a_j [kg*s/m^2] (量纲取决于推导，此处为简化表示)
        self.B_coeffs = np.zeros(self.n_junc)  # 源项 b_j [kg/s]

        # 辅助物理场 (用于系数计算)
        self.rho_vec = np.zeros(self.n_vol)  # Density [kg/m^3]
        self.mu_vec = np.zeros(self.n_vol)  # Viscosity [Pa.s]
        self.inlet_relaxation_gain = 50.0
        self.zero_source_vec = np.zeros(self.n_vol)
        self.pressure_rhs_buffer = np.zeros(self.n_vol)
        self.pressure_matrix = csr_matrix((self.n_vol, self.n_vol), dtype=float)
        self.pressure_diag_ptrs = np.empty(0, dtype=np.int32)
        self.pressure_junction_slot_ptrs = np.empty(0, dtype=np.int32)
        self.pressure_junction_slot_junc_idx = np.empty(0, dtype=np.int32)
        self.pressure_junction_slot_scale = np.empty(0)
        self.pressure_rhs_rows = np.empty(0, dtype=np.int32)
        self.pressure_rhs_junc_idx = np.empty(0, dtype=np.int32)
        self.pressure_rhs_scale = np.empty(0)
        self.primary_material = next(
            (getattr(vol, 'material', None) for vol in self.volumes_obj if getattr(vol, 'material', None) is not None),
            None
        )
        self.material_has_drho_dp = hasattr(self.primary_material, 'liquid_density_derivative_P')
        self.material_has_drho_dt = hasattr(self.primary_material, 'liquid_density_derivative_T')
        self.material_has_temp_from_h = hasattr(self.primary_material, 'temperature_from_enthalpy')

        self.energy_matrix = csr_matrix((self.n_vol, self.n_vol), dtype=float)
        self.energy_diag_ptrs = np.empty(0, dtype=np.int32)
        self.energy_from_diag_ptrs = np.empty(0, dtype=np.int32)
        self.energy_from_offdiag_ptrs = np.empty(0, dtype=np.int32)
        self.energy_to_diag_ptrs = np.empty(0, dtype=np.int32)
        self.energy_to_offdiag_ptrs = np.empty(0, dtype=np.int32)
        self.energy_from_junc_idx = np.empty(0, dtype=np.int32)
        self.energy_to_junc_idx = np.empty(0, dtype=np.int32)
        self.energy_from_multiplier = np.empty(0)
        self.energy_to_multiplier = np.empty(0)
        self.energy_mass_buffer = np.zeros(self.n_vol)
        self.energy_cp_buffer = np.ones(self.n_vol)
        self.energy_lam_h_buffer = np.zeros(self.n_vol)
        self.energy_q_mod_buffer = np.zeros(self.n_vol)
        self.energy_inertia_buffer = np.zeros(self.n_vol)
        self.energy_diag_buffer = np.zeros(self.n_vol)
        self.energy_rhs_buffer = np.zeros(self.n_vol)
        self.energy_aux_buffer = np.zeros(self.n_vol)
        self.energy_T_old_buffer = np.zeros(self.n_vol)
        self.energy_h_old_buffer = np.zeros(self.n_vol)
        self.enthalpy_derivative_buffer = np.zeros(self.n_vol)
        self.thermal_source_buffer = np.zeros(self.n_vol)
        self.enthalpy_mass_buffer = np.zeros(self.n_vol)
        self.enthalpy_q_total_buffer = np.zeros(self.n_vol)
        self.enthalpy_flux_buffer = np.zeros(self.n_junc)
        self.enthalpy_from_flux_buffer = np.zeros(self.n_junc)
        self.enthalpy_to_flux_buffer = np.zeros(self.n_junc)
        self.enthalpy_donor_idx_buffer = np.zeros(self.n_junc, dtype=np.int32)
        self.drho_dh_buffer = np.zeros(self.n_vol)

        # --- 执行构建 ---
        self._build_topology()
        self._build_pressure_system_cache()
        self._build_energy_system_cache()
        self._fixed_pressure_indices_locked = True
        self._initialize_state_from_objects()

        # 初始化物性 (防止后续计算系数时 rho=0 导致除零)
        # 这一步依赖于 _initialize_state_from_objects 已正确读取了 T/h
        # 我们稍后在 Phase 2 的代码中实现 _update_fluid_properties，这里先留个占位或者手动调用的逻辑
        self._update_fluid_properties()

        logger.info(f"HydraulicNetwork initialized: {self.n_vol} Nodes, {self.n_junc} Links. Gravity={self.g}")

    @property
    def fixed_pressure_indices(self) -> Set[int]:
        """定压边界索引集合。初始化完成后禁止修改成员关系。"""
        return self._fixed_pressure_indices

    def _build_topology(self):
        """
        [Phase 1 核心] 构建拓扑映射关系
        将物理对象转译为矩阵行号，并识别稳压器。
        """
        # Step 1: 建立节点索引映射
        for idx, vol in enumerate(self.volumes_obj):
            self.vol_to_idx[vol] = idx

            # --- [性能优化] 提取节点恒定几何参数 ---
            # 对于明确存在的属性，直接访问。这比底层的字典寻址要快得多。
            self.V_vec[idx] = vol.vol
            self.L_node_vec[idx] = vol.len
            self.A_node_vec[idx] = vol.area
            self.Dh_node_vec[idx] = vol.d_h

            # z_coordinate 在 FluidVolume 初始化时给了 0.0，但在特殊边界中可能遗漏，安全起见用 getattr 兜底
            self.z_vec[idx] = getattr(vol, 'z_coordinate', 0.0)

            # 检查是否为稳压器 (Pressure Boundary)
            # 使用 getattr 安全获取属性，默认为 False
            if getattr(vol, 'is_pressure_boundary', False):
                self.fixed_pressure_indices.add(idx)
                logger.debug(f"Node [{idx}] identified as Fixed Pressure Boundary.")

        if not self.fixed_pressure_indices:
            logger.warning("⚠️ No Pressure Boundary (Pressurizer) detected! "
                           "Matrix might be singular for closed loops.")

        # Step 2: 解析连接关系
        for j_idx, junc in enumerate(self.junctions_obj):
            vol_from = junc.from_vol
            vol_to = junc.to_vol

            # 完整性检查
            if vol_from not in self.vol_to_idx:
                raise ValueError(f"Junction '{getattr(junc, 'name', 'Unknown')}' connects FROM unknown volume.")
            if vol_to not in self.vol_to_idx:
                raise ValueError(f"Junction '{getattr(junc, 'name', 'Unknown')}' connects TO unknown volume.")

            idx_from = self.vol_to_idx[vol_from]
            idx_to = self.vol_to_idx[vol_to]
            self.idx_from_vec[j_idx] = idx_from
            self.idx_to_vec[j_idx] = idx_to

            # --- [性能优化] 提取连接恒定几何参数 ---
            # FlowJunction 类已明确定义了这些属性
            self.A_junc_vec[j_idx] = junc.area
            self.L_junc_vec[j_idx] = junc.length
            self.K_loss_vec[j_idx] = junc.k_loss
            self.Dh_from_vec[j_idx] = vol_from.d_h
            self.Dh_to_vec[j_idx] = vol_to.d_h
            self.L_from_half_vec[j_idx] = 0.5 * vol_from.len
            self.L_to_half_vec[j_idx] = 0.5 * vol_to.len
            self.A_from_node_vec[j_idx] = vol_from.area
            self.A_to_node_vec[j_idx] = vol_to.area

            # --- 提取宏观/微观乘子 (普通接管默认返回 1.0) ---
            self.M_from_vec[j_idx] = getattr(junc, 'multiplier_from', 1.0)
            self.M_to_vec[j_idx] = getattr(junc, 'multiplier_to', 1.0)
            self.is_inlet_junction_mask[j_idx] = hasattr(junc, 'target_W')
            if self.is_inlet_junction_mask[j_idx]:
                self.target_W_vec[j_idx] = junc.target_W

            # 存储描述符
            self.junction_descriptors.append((junc, idx_from, idx_to))

        self.inlet_junction_indices = np.flatnonzero(self.is_inlet_junction_mask).astype(np.int32, copy=False)
        logger.debug("Topology mapping built successfully.")

    def _initialize_state_from_objects(self):
        """
        从物理对象同步初始状态到求解器向量
        """
        for idx, vol in enumerate(self.volumes_obj):
            self.P_vec[idx] = vol.P
            self.h_vec[idx] = vol.h
            # 尝试读取温度，如果没有则根据焓反算 (或者留0，由后续 update_properties 处理)
            if hasattr(vol, 'T'):
                self.T_vec[idx] = vol.T
            else:
                # 如果没有 T 属性，这里可以暂存 0，后续依靠 h 计算
                pass

        for idx, junc in enumerate(self.junctions_obj):
            self.W_vec[idx] = junc.W

    def _refresh_cached_boundary_targets(self):
        """刷新入口边界的目标流量缓存。"""
        for j_idx in self.inlet_junction_indices:
            self.target_W_vec[j_idx] = self.junctions_obj[j_idx].target_W

    @staticmethod
    def _locate_csr_entries(indptr: np.ndarray, indices: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """定位给定坐标在 CSR 数据区中的位置。"""
        ptrs = np.empty(len(rows), dtype=np.int32)
        for i, (row, col) in enumerate(zip(rows, cols)):
            start = indptr[row]
            end = indptr[row + 1]
            local = np.searchsorted(indices[start:end], col)
            ptr = start + local
            if ptr >= end or indices[ptr] != col:
                raise KeyError(f"CSR entry ({row}, {col}) not found in pressure matrix cache.")
            ptrs[i] = ptr
        return ptrs

    def _refresh_cached_pressure_targets(self):
        """刷新定压边界的目标压力缓存。"""
        for local_idx, vol_idx in enumerate(self.fixed_pressure_idx_arr):
            vol = self.volumes_obj[vol_idx]
            self.fixed_target_P_vec[local_idx] = getattr(vol, 'target_P', self.P_vec[vol_idx])

    def _build_pressure_system_cache(self):
        """预构建固定拓扑下的压力矩阵结构和 RHS 映射。"""
        self.fixed_pressure_mask.fill(False)
        if self.fixed_pressure_indices:
            self.fixed_pressure_mask[list(self.fixed_pressure_indices)] = True

        np.logical_not(self.fixed_pressure_mask, out=self.normal_pressure_mask)
        self.fixed_pressure_idx_arr = np.flatnonzero(self.fixed_pressure_mask).astype(np.int32, copy=False)
        self.normal_pressure_idx_arr = np.flatnonzero(self.normal_pressure_mask).astype(np.int32, copy=False)
        self.fixed_target_P_vec = np.zeros(self.fixed_pressure_idx_arr.size)
        self._refresh_cached_pressure_targets()

        diag_rows = np.arange(self.n_vol, dtype=np.int32)
        diag_cols = diag_rows.copy()

        from_active = ~self.fixed_pressure_mask[self.idx_from_vec]
        to_active = ~self.fixed_pressure_mask[self.idx_to_vec]
        from_junc_idx = np.flatnonzero(from_active).astype(np.int32, copy=False)
        to_junc_idx = np.flatnonzero(to_active).astype(np.int32, copy=False)

        slot_rows_parts = []
        slot_cols_parts = []
        slot_junc_idx_parts = []
        slot_scale_parts = []

        if from_junc_idx.size > 0:
            slot_rows_parts.extend([
                self.idx_from_vec[from_junc_idx],
                self.idx_from_vec[from_junc_idx]
            ])
            slot_cols_parts.extend([
                self.idx_from_vec[from_junc_idx],
                self.idx_to_vec[from_junc_idx]
            ])
            slot_junc_idx_parts.extend([from_junc_idx, from_junc_idx])
            slot_scale_parts.extend([
                self.M_from_vec[from_junc_idx],
                -self.M_from_vec[from_junc_idx]
            ])

        if to_junc_idx.size > 0:
            slot_rows_parts.extend([
                self.idx_to_vec[to_junc_idx],
                self.idx_to_vec[to_junc_idx]
            ])
            slot_cols_parts.extend([
                self.idx_to_vec[to_junc_idx],
                self.idx_from_vec[to_junc_idx]
            ])
            slot_junc_idx_parts.extend([to_junc_idx, to_junc_idx])
            slot_scale_parts.extend([
                self.M_to_vec[to_junc_idx],
                -self.M_to_vec[to_junc_idx]
            ])

        if slot_rows_parts:
            slot_rows = np.concatenate(slot_rows_parts).astype(np.int32, copy=False)
            slot_cols = np.concatenate(slot_cols_parts).astype(np.int32, copy=False)
            slot_junc_idx = np.concatenate(slot_junc_idx_parts).astype(np.int32, copy=False)
            slot_scale = np.concatenate(slot_scale_parts)
        else:
            slot_rows = np.empty(0, dtype=np.int32)
            slot_cols = np.empty(0, dtype=np.int32)
            slot_junc_idx = np.empty(0, dtype=np.int32)
            slot_scale = np.empty(0)

        struct_rows = np.concatenate([diag_rows, slot_rows])
        struct_cols = np.concatenate([diag_cols, slot_cols])
        struct_data = np.ones(struct_rows.size)
        self.pressure_matrix = coo_matrix((struct_data, (struct_rows, struct_cols)), shape=(self.n_vol, self.n_vol)).tocsr()
        self.pressure_matrix.data[:] = 0.0

        self.pressure_diag_ptrs = self._locate_csr_entries(
            self.pressure_matrix.indptr,
            self.pressure_matrix.indices,
            diag_rows,
            diag_cols
        )

        if slot_rows.size > 0:
            self.pressure_junction_slot_ptrs = self._locate_csr_entries(
                self.pressure_matrix.indptr,
                self.pressure_matrix.indices,
                slot_rows,
                slot_cols
            )
            self.pressure_junction_slot_junc_idx = slot_junc_idx
            self.pressure_junction_slot_scale = slot_scale
        else:
            self.pressure_junction_slot_ptrs = np.empty(0, dtype=np.int32)
            self.pressure_junction_slot_junc_idx = np.empty(0, dtype=np.int32)
            self.pressure_junction_slot_scale = np.empty(0)

        rhs_rows_parts = []
        rhs_junc_idx_parts = []
        rhs_scale_parts = []

        if from_junc_idx.size > 0:
            rhs_rows_parts.append(self.idx_from_vec[from_junc_idx])
            rhs_junc_idx_parts.append(from_junc_idx)
            rhs_scale_parts.append(-self.M_from_vec[from_junc_idx])

        if to_junc_idx.size > 0:
            rhs_rows_parts.append(self.idx_to_vec[to_junc_idx])
            rhs_junc_idx_parts.append(to_junc_idx)
            rhs_scale_parts.append(self.M_to_vec[to_junc_idx])

        if rhs_rows_parts:
            self.pressure_rhs_rows = np.concatenate(rhs_rows_parts).astype(np.int32, copy=False)
            self.pressure_rhs_junc_idx = np.concatenate(rhs_junc_idx_parts).astype(np.int32, copy=False)
            self.pressure_rhs_scale = np.concatenate(rhs_scale_parts)
        else:
            self.pressure_rhs_rows = np.empty(0, dtype=np.int32)
            self.pressure_rhs_junc_idx = np.empty(0, dtype=np.int32)
            self.pressure_rhs_scale = np.empty(0)

    def _build_energy_system_cache(self):
        """预构建固定拓扑下的能量矩阵结构。"""
        diag_rows = np.arange(self.n_vol, dtype=np.int32)
        diag_cols = diag_rows.copy()

        self.energy_from_junc_idx = np.flatnonzero(~self.fixed_pressure_mask[self.idx_from_vec]).astype(np.int32, copy=False)
        self.energy_to_junc_idx = np.flatnonzero(~self.fixed_pressure_mask[self.idx_to_vec]).astype(np.int32, copy=False)
        self.energy_from_multiplier = self.M_from_vec[self.energy_from_junc_idx].copy()
        self.energy_to_multiplier = self.M_to_vec[self.energy_to_junc_idx].copy()

        rows_parts = [diag_rows]
        cols_parts = [diag_cols]

        if self.energy_from_junc_idx.size > 0:
            from_rows = self.idx_from_vec[self.energy_from_junc_idx]
            rows_parts.extend([from_rows, from_rows])
            cols_parts.extend([
                self.idx_from_vec[self.energy_from_junc_idx],
                self.idx_to_vec[self.energy_from_junc_idx]
            ])

        if self.energy_to_junc_idx.size > 0:
            to_rows = self.idx_to_vec[self.energy_to_junc_idx]
            rows_parts.extend([to_rows, to_rows])
            cols_parts.extend([
                self.idx_to_vec[self.energy_to_junc_idx],
                self.idx_from_vec[self.energy_to_junc_idx]
            ])

        struct_rows = np.concatenate(rows_parts).astype(np.int32, copy=False)
        struct_cols = np.concatenate(cols_parts).astype(np.int32, copy=False)
        struct_data = np.ones(struct_rows.size)
        self.energy_matrix = coo_matrix((struct_data, (struct_rows, struct_cols)), shape=(self.n_vol, self.n_vol)).tocsr()
        self.energy_matrix.data[:] = 0.0

        self.energy_diag_ptrs = self._locate_csr_entries(
            self.energy_matrix.indptr,
            self.energy_matrix.indices,
            diag_rows,
            diag_cols
        )

        if self.energy_from_junc_idx.size > 0:
            from_rows = self.idx_from_vec[self.energy_from_junc_idx]
            self.energy_from_diag_ptrs = self._locate_csr_entries(
                self.energy_matrix.indptr,
                self.energy_matrix.indices,
                from_rows,
                self.idx_from_vec[self.energy_from_junc_idx]
            )
            self.energy_from_offdiag_ptrs = self._locate_csr_entries(
                self.energy_matrix.indptr,
                self.energy_matrix.indices,
                from_rows,
                self.idx_to_vec[self.energy_from_junc_idx]
            )
        else:
            self.energy_from_diag_ptrs = np.empty(0, dtype=np.int32)
            self.energy_from_offdiag_ptrs = np.empty(0, dtype=np.int32)

        if self.energy_to_junc_idx.size > 0:
            to_rows = self.idx_to_vec[self.energy_to_junc_idx]
            self.energy_to_diag_ptrs = self._locate_csr_entries(
                self.energy_matrix.indptr,
                self.energy_matrix.indices,
                to_rows,
                self.idx_to_vec[self.energy_to_junc_idx]
            )
            self.energy_to_offdiag_ptrs = self._locate_csr_entries(
                self.energy_matrix.indptr,
                self.energy_matrix.indices,
                to_rows,
                self.idx_from_vec[self.energy_to_junc_idx]
            )
        else:
            self.energy_to_diag_ptrs = np.empty(0, dtype=np.int32)
            self.energy_to_offdiag_ptrs = np.empty(0, dtype=np.int32)

    def debug_topology(self):
        """
        [调试工具] 打印网络拓扑结构
        """
        print("\n=== Hydraulic Network Topology ===")
        print(f"Nodes ({self.n_vol}):")
        for vol, idx in self.vol_to_idx.items():
            type_str = " [FIXED]" if idx in self.fixed_pressure_indices else ""
            name = getattr(vol, 'name', f'Vol_{idx}')
            print(f"  [{idx}] {name}{type_str}")

        print(f"\nLinks ({self.n_junc}):")
        for j_desc in self.junction_descriptors:
            junc, i_from, i_to = j_desc
            name = getattr(junc, 'name', 'Link')
            print(f"  Link: [{i_from}] -> [{i_to}] via '{name}'")
        print("==================================\n")

    # =========================================================================
    # Phase 2: Coefficient Calculation (Physics Engine)
    # =========================================================================
    @TEASAProfiler.profile
    def _update_fluid_properties(self):
        """
        更新所有节点的流体物性 (rho, mu)
        依赖: vol.material, vol.T (或 vol.h)
        """

        if self.n_vol == 0:
            return

        # =========================================================================
        # 1. 向量化物性计算 (核心性能提升区)
        #    系统使用单一流体回路，统一流网采用相同流体
        # =========================================================================
        mat = self.primary_material

        if mat is not None:
            # 直接传入整个系统的温度向量 T_vec 和压力向量 P_vec
            # 单次 C/NumPy 级别调用，取代 N 次 Python 级别调用
            self.cp_vec[:] = mat.heat_capacity(self.T_vec, self.P_vec)
            raw_rho = mat.density(self.T_vec, self.P_vec)
            raw_mu = mat.viscosity(self.T_vec, self.P_vec)

            # 偏导数计算 (用于 Jacobian 和压缩性矩阵构建)
            if self.material_has_drho_dp:
                raw_drho_dp = mat.liquid_density_derivative_P(self.P_vec)
            else:
                self.drho_dp_vec.fill(1e-9)
                raw_drho_dp = self.drho_dp_vec

            if self.material_has_drho_dt:
                self.drho_dt_vec[:] = mat.liquid_density_derivative_T(self.T_vec)
            else:
                self.drho_dt_vec.fill(0.0)

            # --- 数值安全与物理边界保护 (全数组并行处理) ---
            np.maximum(raw_rho, 1e-1, out=self.rho_vec)  # 避免除零
            np.maximum(raw_mu, 1e-10, out=self.mu_vec)  # 避免极小粘度导致雷诺数溢出
            np.maximum(raw_drho_dp, 1e-11, out=self.drho_dp_vec)  # 提供微弱的人工可压缩性兜底

        else:
            # 兜底：如果未指定物性库，赋予安全的恒定默认值
            self.cp_vec.fill(1000.0)
            self.drho_dp_vec.fill(1e-9)
            self.drho_dt_vec.fill(0.0)
            # rho 和 mu 保持初始/当前值，但施加最小限制
            np.maximum(self.rho_vec, 1e-1, out=self.rho_vec)
            np.maximum(self.mu_vec, 1e-10, out=self.mu_vec)

        q_expl = self.Q_expl_vec
        lam_imp = self.lam_imp_vec
        for i, vol in enumerate(self.volumes_obj):
            vol.rho = self.rho_vec[i]
            vol.mu = self.mu_vec[i]
            q_expl[i] = vol.Q_wall + vol.Q_vol
            lam_imp[i] = vol.implicit_coeff

    @staticmethod
    def _calc_friction_factor_static(Re: float) -> float:
        """
        计算单相摩擦系数 f (静态方法，复刻 Fortran 逻辑)
        """
        Re_calc = max(Re, 1e-5)

        if Re_calc <= 1000.0:
            return 64.0 / Re_calc
        elif 2300.0 < Re_calc < 100000.0:
            return 0.3164 / (Re_calc ** 0.25)
        else:
            # Karman-Prandtl 隐式迭代
            f = 0.3164 / (Re_calc ** 0.25)
            max_iter = 50
            for _ in range(max_iter):
                sqrt_f = math.sqrt(f)
                argument = 2.51 / (Re_calc * sqrt_f)
                denom = -2.0 * math.log10(argument)
                f_new = (1.0 / denom) ** 2.0
                if abs(f_new - f) <= 1.0e-8:
                    return f_new
                f = f_new
            return f

    @staticmethod
    def _calc_friction_factor_static_vec(Re: np.ndarray) -> np.ndarray:
        """Vectorized equivalent of `_calc_friction_factor_static`."""
        Re_calc = np.maximum(np.asarray(Re, dtype=float), 1e-5)
        f = np.empty_like(Re_calc)

        laminar_mask = Re_calc <= 1000.0
        blasius_mask = (Re_calc > 2300.0) & (Re_calc < 100000.0)
        turbulent_mask = ~(laminar_mask | blasius_mask)

        f[laminar_mask] = 64.0 / Re_calc[laminar_mask]
        f[blasius_mask] = 0.3164 / np.power(Re_calc[blasius_mask], 0.25)

        if np.any(turbulent_mask):
            re_turb = Re_calc[turbulent_mask]
            f_turb = 0.3164 / np.power(re_turb, 0.25)

            for _ in range(50):
                sqrt_f = np.sqrt(f_turb)
                argument = 2.51 / (re_turb * sqrt_f)
                denom = -2.0 * np.log10(argument)
                f_new = np.power(1.0 / denom, 2.0)
                if np.max(np.abs(f_new - f_turb)) <= 1.0e-8:
                    f_turb = f_new
                    break
                f_turb = f_new

            f[turbulent_mask] = f_turb

        return f

    def _get_upwind_properties(self, idx_in: int, idx_out: int, W_flow: float) -> tuple[np.ndarray[tuple[Any, ...], Any], np.ndarray[tuple[Any, ...], Any]]:
        """
        [辅助] 根据流向获取迎风侧物性 (rho, mu)
        增加迎风侧平滑，可在极低流量下保证计算平话
        """
        # 平滑阈值
        eps = 1e-6

        # 计算平滑因子
        # 当 W_iterate >> eps 时，alpha -> 1 (完全取 in 侧)
        # 当 W_iterate << -eps 时，alpha -> 0 (完全取 out 侧)
        alpha = 0.5 * (1.0 + np.tanh(W_flow / eps))

        # 加权混合
        rho_face = alpha * self.rho_vec[idx_in] + (1.0 - alpha) * self.rho_vec[idx_out]
        mu_face = alpha * self.mu_vec[idx_in] + (1.0 - alpha) * self.mu_vec[idx_out]

        return rho_face, mu_face

        # 如果 W >= 0 (in -> out)，取上游(in)物性
        # 如果 W < 0 (out -> in)，取下游(out)物性 （已简化）
        # if W_flow >= 0:
        #     return self.rho_vec[idx_in], self.mu_vec[idx_in]
        # else:
        #     return self.rho_vec[idx_out], self.mu_vec[idx_out]

    def _get_inertial_length(self, idx_in: int, idx_out: int) -> float:
        """
        [辅助] 计算连接的惯性特征长度
        L_junc = 0.5 * L_upstream + 0.5 * L_downstream
        """
        # 获取 Volume 对象
        vol_in = self.volumes_obj[idx_in]
        vol_out = self.volumes_obj[idx_out]

        # 获取长度，如果没有 length 属性 (如 Plenum)，默认为 0.0
        L_in = getattr(vol_in, 'length', 0.0)
        L_out = getattr(vol_out, 'length', 0.0)

        # 惯性长度取平均
        L_total = 0.5 * (L_in + L_out)

        # 防止长度为 0 (如两个 Plenum 直接相连)，给极小值避免除零
        return max(L_total, 1e-6)

    @TEASAProfiler.profile
    def _calc_momentum_coeffs(self, dt: float, W_old_frozen: np.ndarray = None):
        """
        [核心] 计算动量方程系数 a_j, b_j

        方程: W^{n+1} = a_j * (P_in^{n+1} - P_out^{n+1}) + b_j

        物理模型:
        1. 惯性项: (L/A) * dW/dt
        2. 摩擦阻力: Darcy-Weisbach (显式线性化)
        3. 局部阻力: Form Loss (显式线性化)
        4. 重力项: rho * g * dz (Source)
        5. 空间加速: Bernoulli 效应 (Source)
        """
        return self._calc_momentum_coeffs_fast(dt, W_old_frozen=W_old_frozen)

        # 作为系数，需要考虑 3 种流量变量：
        # 1. 此次迭代中的流量变量，即本次迭代待求出的变量：W_vec
        # 2. 上一时间步的流量变量，用于构建源项中的惯性项系数和加速压降：W_old
        # 3. 本时间步上一次迭代的流量变量，用于构建系数矩阵中的对角线系数部分：W_iterate

        # 极小值保护
        EPS = 1e-10

        # 冻结的旧状态向量，如果没有提供历史流量，则采用当前流量
        # 先保留
        if W_old_frozen is None:
            W_old_frozen = self.W_vec

        for j_idx, (junc, idx_in, idx_out) in enumerate(self.junction_descriptors):
            # --- 0. 获取解耦的三种流量变量 ---
            W_old = W_old_frozen[j_idx]
            W_iterate = self.W_vec[j_idx]

            # ==========================================================
            # 【新增修复】：强行拦截 InletJunction (定流量边界)
            # ==========================================================
            if hasattr(junc, 'target_W'):
                # 获取松弛增益 (建议与 BoundaryVolume.py 中保持一致)
                gain_k = 50.0
                target_W = junc.target_W

                # 计算受控的新流量 (显式推进)
                W_next = W_old + dt * gain_k * (target_W - W_old)

                # 切断压力耦合 (a_j=0)，将其彻底变为强制质量源 (b_j=W_next)
                self.A_coeffs[j_idx] = 0.0
                self.B_coeffs[j_idx] = W_next

                # 跳过后续常规的摩擦、压降和惯性项计算！
                continue
            # ==========================================================

            # 阻力计算系数，采用 W_iterate 的绝对值
            W_iterate_abs = abs(W_iterate)

            # --- 1. 几何与状态参数 ---
            # 几何参数
            A_flow = max(getattr(junc, 'area', 1e-4), EPS)
            D_h_upwind = junc.from_vol.d_h
            D_h_downwind = junc.to_vol.d_h
            L_upwind = 0.5 * junc.from_vol.len
            L_downwind = 0.5 * junc.to_vol.len
            L_inertial = L_upwind + L_downwind
            # D_h = max(getattr(junc, 'hydraulic_diameter', min(junc.from_vol.hydraulic_diam),junc.to_vol.hydraulic_diam), EPS)
            # L_inertial = max(getattr(junc, 'custom_length', self._get_inertial_length(idx_in, idx_out)), EPS)
            # L_inertial = self._get_inertial_length(idx_in, idx_out)

            # 上一时刻状态
            # W_old = W_old_frozen[j_idx]
            # W_old = self.W_vec[j_idx]
            # W_abs = abs(W_old)

            # 获取迎风物性 (已经在调用的函数中修改过)
            rho_j, mu_j = self._get_upwind_properties(idx_in, idx_out, W_iterate)

            # --- 2. 计算雷诺数 ---
            # vel = W / (rho * A)
            vel = W_iterate_abs / (rho_j * A_flow)
            Re_upwind = (rho_j * vel * D_h_upwind) / mu_j
            Re_downwind = (rho_j * vel * D_h_downwind) / mu_j

            # --- 3. 计算阻力系数 (Resistance) ---
            # 3.1 摩擦系数 f
            f_upwind = self._calc_friction_factor_static(Re_upwind)
            f_downwind = self._calc_friction_factor_static(Re_downwind)

            # 3.2 局部阻力系数 K_form
            K_form = getattr(junc, 'k_loss', 0.0)

            # 3.3 总线性化阻力系数 K_linear [1/(m*s)]
            # 动量方程阻力项: dP = K_linear * W
            # K_linear = (f*L/D + K_form) * |W| / (2*rho*A^2)
            # 注意：这里的 L 取决于摩擦发生的长度。
            # 通常对于 Junction，如果它代表一段管子，应取 vol.length 的一部分。
            # 这里为简便且保守，使用惯性长度 L_inertial 作为摩擦长度。
            # 需要采用 W_iterate_abs 进行线性化
            term_geom = (f_upwind * L_upwind / D_h_upwind + f_downwind * L_downwind / D_h_downwind) + K_form
            K_linear = (term_geom * W_iterate_abs) / (2.0 * rho_j * A_flow ** 2)

            # --- 4. 计算惯性项系数 (Inertia) ---
            # I = L / (A * dt) [1/(m*s)]
            I_term = L_inertial / (A_flow * dt)

            # --- 5. 计算源项 (Gravity & Acceleration) ---

            # 5.1 重力项 (Gravity)
            # z 轴向上为正。流体从高(z_in)流向低(z_out)受重力驱动(正源项)。
            # dP_grav = rho * g * (z_in - z_out)
            z_in = getattr(self.volumes_obj[idx_in], 'z', 0.0)
            z_out = getattr(self.volumes_obj[idx_out], 'z', 0.0)
            dP_grav = rho_j * self.g * (z_in - z_out)

            # 5.2 空间加速项 (Spatial Acceleration)
            # dP_acc = 0.5 * W^2 * (1/(rho_in*A_in^2) - 1/(rho_out*A_out^2))
            # 同样使用迎风思想，或者简单的中心差分思想。
            # 获取节点流通面积 (Plenum 默认为无穷大->1e12, 1/A^2 -> 0)
            A_in_node = getattr(self.volumes_obj[idx_in], 'area', 1e12)
            A_out_node = getattr(self.volumes_obj[idx_out], 'area', 1e12)

            # 防御性保护：防止面积为0
            A_in_node = max(A_in_node, EPS)
            A_out_node = max(A_out_node, EPS)

            # 仍然保持显式，作为源项直接作用，使用 W_old 防止Picard 迭代震荡
            # ====== [TEASA 修正]: 必须考虑乘子导致的宏观节点真实流速放大 ======
            M_in = self.M_from_vec[j_idx]
            M_out = self.M_to_vec[j_idx]
            term_acc_in = (M_in * M_in) / (self.rho_vec[idx_in] * A_in_node ** 2)
            term_acc_out = (M_out * M_out) / (self.rho_vec[idx_out] * A_out_node ** 2)
            dP_acc = 0.5 * (W_old ** 2) * (term_acc_in - term_acc_out)

            # 5.3 泵扬程 (预留)
            dP_pump = 0.0

            # --- 6. 汇总计算最终系数 ---
            # 离散方程: (I + K_linear) * W^{n+1} = (P_in - P_out) + (I * W^n + Sources)
            # 形式: W^{n+1} = a_j * (P_in - P_out) + b_j

            # 导纳 a_j
            # 必须保证分母不为0 (I_term通常很大，安全)
            self.A_coeffs[j_idx] = 1.0 / (I_term + K_linear)

            # 源项 b_j
            source_total = (I_term * W_old) + dP_grav + dP_acc + dP_pump
            self.B_coeffs[j_idx] = self.A_coeffs[j_idx] * source_total

    def _calc_momentum_coeffs_fast(self, dt: float, W_old_frozen: np.ndarray = None):
        """基于缓存数组的动量系数快速路径。"""
        if self.n_junc == 0:
            return

        eps = 1.0e-10
        dt_safe = max(dt, eps)

        self._refresh_cached_boundary_targets()

        if W_old_frozen is None:
            W_old_frozen = self.W_vec

        idx_from = self.idx_from_vec
        idx_to = self.idx_to_vec
        inlet_mask = self.is_inlet_junction_mask
        normal_mask = ~inlet_mask

        self.A_coeffs.fill(0.0)
        self.B_coeffs.fill(0.0)

        if np.any(inlet_mask):
            W_old_inlet = np.asarray(W_old_frozen[inlet_mask], dtype=float)
            target_W = self.target_W_vec[inlet_mask]
            self.B_coeffs[inlet_mask] = (
                W_old_inlet
                + dt_safe * self.inlet_relaxation_gain * (target_W - W_old_inlet)
            )

        if not np.any(normal_mask):
            return

        W_iterate = self.W_vec[normal_mask]
        W_old = np.asarray(W_old_frozen[normal_mask], dtype=float)

        idx_from_n = idx_from[normal_mask]
        idx_to_n = idx_to[normal_mask]

        A_flow = np.maximum(self.A_junc_vec[normal_mask], eps)
        D_h_upwind = np.maximum(self.Dh_from_vec[normal_mask], eps)
        D_h_downwind = np.maximum(self.Dh_to_vec[normal_mask], eps)
        L_upwind = self.L_from_half_vec[normal_mask]
        L_downwind = self.L_to_half_vec[normal_mask]
        L_inertial = np.maximum(L_upwind + L_downwind, eps)

        rho_in = np.maximum(self.rho_vec[idx_from_n], eps)
        rho_out = np.maximum(self.rho_vec[idx_to_n], eps)
        mu_in = np.maximum(self.mu_vec[idx_from_n], eps)
        mu_out = np.maximum(self.mu_vec[idx_to_n], eps)

        alpha = 0.5 * (1.0 + np.tanh(W_iterate / 1.0e-6))
        rho_j = np.maximum(alpha * rho_in + (1.0 - alpha) * rho_out, eps)
        mu_j = np.maximum(alpha * mu_in + (1.0 - alpha) * mu_out, eps)

        W_iterate_abs = np.abs(W_iterate)
        vel = W_iterate_abs / (rho_j * A_flow)
        Re_upwind = (rho_j * vel * D_h_upwind) / mu_j
        Re_downwind = (rho_j * vel * D_h_downwind) / mu_j

        f_upwind = self._calc_friction_factor_static_vec(Re_upwind)
        f_downwind = self._calc_friction_factor_static_vec(Re_downwind)

        term_geom = (
            f_upwind * L_upwind / D_h_upwind
            + f_downwind * L_downwind / D_h_downwind
            + self.K_loss_vec[normal_mask]
        )
        K_linear = (term_geom * W_iterate_abs) / (2.0 * rho_j * A_flow ** 2)
        I_term = L_inertial / (A_flow * dt_safe)

        dz = self.z_vec[idx_from_n] - self.z_vec[idx_to_n]
        dP_grav = rho_j * self.g * dz

        A_in_node = np.maximum(self.A_from_node_vec[normal_mask], eps)
        A_out_node = np.maximum(self.A_to_node_vec[normal_mask], eps)
        M_in = self.M_from_vec[normal_mask]
        M_out = self.M_to_vec[normal_mask]
        term_acc_in = (M_in * M_in) / (rho_in * A_in_node ** 2)
        term_acc_out = (M_out * M_out) / (rho_out * A_out_node ** 2)
        dP_acc = 0.5 * (W_old ** 2) * (term_acc_in - term_acc_out)

        denom = np.maximum(I_term + K_linear, eps)
        A_normal = 1.0 / denom
        source_total = (I_term * W_old) + dP_grav + dP_acc

        self.A_coeffs[normal_mask] = A_normal
        self.B_coeffs[normal_mask] = A_normal * source_total

    def debug_coefficients(self):
        """[调试工具] 打印计算出的系数"""
        print("\n=== Momentum Coefficients (Phase 2) ===")
        print(f"{'Link':<10} {'W_old':<10} {'Re':<10} {'a_j (Adm)':<15} {'b_j (Src)':<15} {'dP_grav':<10}")
        for j_idx, (junc, i_in, i_out) in enumerate(self.junction_descriptors):
            W = self.W_vec[j_idx]
            # 简单估算 Re 用于显示
            A = max(getattr(junc, 'flow_area', 1e-4), 1e-4)
            D = getattr(junc, 'hydraulic_diameter', 0.03)
            mu = 0.5 * (self.mu_vec[i_in] + self.mu_vec[i_out])
            Re = abs(W) * D / (A * mu) if mu > 0 else 0

            # 获取重力压降用于调试
            rho = 0.5 * (self.rho_vec[i_in] + self.rho_vec[i_out])
            z_in = getattr(self.volumes_obj[i_in], 'z', 0.0)
            z_out = getattr(self.volumes_obj[i_out], 'z', 0.0)
            dp_g = rho * self.g * (z_in - z_out)

            name = getattr(junc, 'name', f'Link_{j_idx}')
            print(
                f"{name:<10} {W:<10.4f} {Re:<10.1e} {self.A_coeffs[j_idx]:<15.2e} {self.B_coeffs[j_idx]:<15.2e} "
                f"{dp_g:<10.1f}")
        print("=======================================\n")

    # =========================================================================
    # Phase 3: Matrix Assembly & Solver (Hydraulic Engine)
    # =========================================================================

    @TEASAProfiler.profile
    def _assemble_pressure_system(
        self,
        dt: float,
        include_thermal_expansion: bool = True,
        S_thermal_vec: np.ndarray = None
    ) -> Tuple[csr_matrix, np.ndarray]:
        """
        [Phase 3] 组装压力矩阵。

        在固定拓扑下，矩阵非零结构预构建并复用，这里只更新 data 和 RHS。
        """
        if S_thermal_vec is None:
            if include_thermal_expansion:
                S_thermal_vec = self._calc_thermal_expansion_source(dt)
            else:
                S_thermal_vec = self.zero_source_vec

        dt_safe = max(dt, 1.0e-10)
        normal_idx = self.normal_pressure_idx_arr
        fixed_idx = self.fixed_pressure_idx_arr

        data = self.pressure_matrix.data
        data.fill(0.0)

        if self.pressure_junction_slot_ptrs.size > 0:
            np.add.at(
                data,
                self.pressure_junction_slot_ptrs,
                self.pressure_junction_slot_scale * self.A_coeffs[self.pressure_junction_slot_junc_idx]
            )

        if normal_idx.size > 0:
            compressibility = (self.V_vec[normal_idx] / dt_safe) * np.maximum(self.drho_dp_vec[normal_idx], 1e-9)
            data[self.pressure_diag_ptrs[normal_idx]] += compressibility
        else:
            compressibility = self.zero_source_vec[:0]

        if fixed_idx.size > 0:
            self._refresh_cached_pressure_targets()
            data[self.pressure_diag_ptrs[fixed_idx]] = 1.0

        B = self.pressure_rhs_buffer
        B.fill(0.0)

        if normal_idx.size > 0:
            B[normal_idx] = compressibility * self.P_vec[normal_idx] + S_thermal_vec[normal_idx]

        if fixed_idx.size > 0:
            B[fixed_idx] = self.fixed_target_P_vec

        if self.pressure_rhs_rows.size > 0:
            np.add.at(
                B,
                self.pressure_rhs_rows,
                self.pressure_rhs_scale * self.B_coeffs[self.pressure_rhs_junc_idx]
            )

        return self.pressure_matrix, B

        # # 初始化矩阵
        # # 采用 LIL 格式的稀疏矩阵
        # A = lil_matrix((self.n_vol, self.n_vol), dtype=float)
        # B = np.zeros(self.n_vol)
        #
        # # --- 流量项填充 ---
        # for j_idx, (_, idx_in, idx_out) in enumerate(self.junction_descriptors):
        #     a_j = self.A_coeffs[j_idx]
        #     b_j = self.B_coeffs[j_idx]
        #
        #     # 上游 (Out)
        #     A[idx_in, idx_in] += a_j
        #     A[idx_in, idx_out] -= a_j
        #     B[idx_in] -= b_j
        #
        #     # 下游 (In)
        #     A[idx_out, idx_in] -= a_j
        #     A[idx_out, idx_out] += a_j
        #     B[idx_out] += b_j
        #
        # # --- 压缩性与源项填充 ---
        # for i, vol in enumerate(self.volumes_obj):
        #     V_cell = getattr(vol, 'vol', 1.0)
        #
        #     # 压缩性 d(rho)/dP
        #     drho_dP = 1e-7  # 人工给定极小值作为默认值
        #     if hasattr(vol, 'material') and vol.material is not None:
        #         try:
        #             drho_dP = max(vol.material.liquid_density_derivative_P(vol.P), 1e-9)
        #         except AttributeError:
        #             pass
        #
        #     # 主对角线: C = (V/dt) * drho/dP
        #     compressibility = (V_cell / dt) * drho_dP
        #     A[i, i] += compressibility
        #     B[i] += compressibility * self.P_vec[i]
        #
        #     # 添加热膨胀源项 (如果在初始化阶段，此项为0)
        #     B[i] += S_thermal_vec[i]
        #
        # # --- 边界条件处理 ---
        # for idx in self.fixed_pressure_indices:
        #     vol = self.volumes_obj[idx]
        #     p_target = getattr(vol, 'target_P', vol.P)
        #
        #     A[idx, :] = 0.0
        #     A[idx, idx] = 1.0
        #     B[idx] = p_target
        #
        # # 将 LIL 格式稀疏矩阵转换为 CSR 格式
        # A_csr = A.tocsr()
        #
        # return A_csr, B

    # 旧版完整矩阵求解构建方法
    # def _assemble_pressure_system(self, dt: float, include_thermal_expansion: bool = True) -> Tuple[
    #     np.ndarray, np.ndarray]:
    #     """
    #     [Phase 3 更新版] 组装压力矩阵
    #     :param dt: 时间步长
    #     :param include_thermal_expansion: 是否包含热膨胀源项 (S_thermal)。
    #                                       瞬态计算应为 True，稳态初始化应为 False。
    #     """
    #     # 1. 计算热膨胀源项 (仅当需要时)
    #     if include_thermal_expansion:
    #         S_thermal_vec = self._calc_thermal_expansion_source(dt)
    #     else:
    #         S_thermal_vec = np.zeros(self.n_vol)
    #
    #     # 初始化矩阵
    #     A = np.zeros((self.n_vol, self.n_vol))
    #     B = np.zeros(self.n_vol)
    #
    #     # --- 流量项填充 (同前) ---
    #     for j_idx, (_, idx_in, idx_out) in enumerate(self.junction_descriptors):
    #         a_j = self.A_coeffs[j_idx]
    #         b_j = self.B_coeffs[j_idx]
    #
    #         # 上游 (Out)
    #         A[idx_in, idx_in] += a_j
    #         A[idx_in, idx_out] -= a_j
    #         B[idx_in] -= b_j
    #
    #         # 下游 (In)
    #         A[idx_out, idx_in] -= a_j
    #         A[idx_out, idx_out] += a_j
    #         B[idx_out] += b_j
    #
    #     # --- 压缩性与源项填充 ---
    #     for i, vol in enumerate(self.volumes_obj):
    #         V_cell = getattr(vol, 'vol', 1.0)
    #
    #         # 压缩性 d(rho)/dP
    #         drho_dP = 1e-7  # Artificial default
    #         if hasattr(vol, 'material') and vol.material is not None:
    #             try:
    #                 drho_dP = max(vol.material.liquid_density_derivative_P(vol.P), 1e-9)
    #             except AttributeError:
    #                 pass
    #
    #         # 主对角线: C = (V/dt) * drho/dP
    #         compressibility = (V_cell / dt) * drho_dP
    #         A[i, i] += compressibility
    #         B[i] += compressibility * self.P_vec[i]
    #
    #         # 添加热膨胀源项 (如果在初始化阶段，此项为0)
    #         B[i] += S_thermal_vec[i]
    #
    #     # --- 边界条件处理 ---
    #     for idx in self.fixed_pressure_indices:
    #         vol = self.volumes_obj[idx]
    #         p_target = getattr(vol, 'target_P', vol.P)
    #
    #         A[idx, :] = 0.0
    #         A[idx, idx] = 1.0
    #         B[idx] = p_target
    #
    #     return A, B

    @TEASAProfiler.profile
    def _solve_linear_system(self, A_sparse, B: np.ndarray):
        """
        [Phase 3 稀疏更新版] 求解压力线性方程组 A * P = B
        Group E 评估后保留 `spsolve` 直解路径：在当前固定拓扑、约 300 节点 benchmark 下，
        `spsolve` 比 `splu`、`bicgstab`、`gmres` 和 dense solve 更快且更稳定。
        """

        try:
            self.P_vec = spsolve(A_sparse, B)
        except Exception as e:
            logger.error(f"Sparse Direct Solver failed! Matrix is likely singular. "
                         f"Check for missing pressure boundaries or disconnected nodes. Error: {e}")
            raise RuntimeError(f"Pressure system solver failed: {e}")

    # 旧版：采用直接 LU 分解求解压力方程
    # def _solve_linear_system(self, A: np.ndarray, B: np.ndarray):
    #     """
    #     调用线性求解器解出 P_new
    #     """
    #     try:
    #         # 使用 numpy 标准求解器 (对于稠密矩阵)
    #         self.P_vec = np.linalg.solve(A, B)
    #     except np.linalg.LinAlgError:
    #         logger.error("Matrix Singularity Detected! Check if pressure boundary exists or system is closed.")
    #         raise

    @TEASAProfiler.profile
    def _update_flow_rates(self):
        """
        回代计算流量 W_new
        W_j = a_j * (P_in - P_out) + b_j
        """
        if self.n_junc == 0:
            return

        pressure_drop = self.P_vec[self.idx_from_vec] - self.P_vec[self.idx_to_vec]
        self.W_vec[:] = self.A_coeffs * pressure_drop + self.B_coeffs

    def step_hydraulic(self, dt: float):
        """
        [Phase 3 主接口] 执行一步纯液力计算

        流程:
        1. 更新物性 (rho, mu)
        2. 计算动量系数 (a_j, b_j)
        3. 组装压力矩阵 (A, B)
        4. 求解压力 (P_new)
        5. 更新流量 (W_new)
        """
        # 1. 准备物理场
        self._update_fluid_properties()

        # 2. 计算系数 (Phase 2)
        self._calc_momentum_coeffs(dt)

        # 3. 组装矩阵 (Phase 3)
        A, B = self._assemble_pressure_system(dt)

        # 4. 求解
        self._solve_linear_system(A, B)

        # 5. 更新状态
        self._update_flow_rates()
        self._sync_vectors_to_objects(
            sync_pressure=True,
            sync_flow=True,
            sync_energy=False,
            sync_properties=False
        )

    # =========================================================================
    # Phase 4: Energy Coupling & Main Stepping (Thermal Engine)
    # =========================================================================

    def _calc_enthalpy_time_derivative_explicit(self) -> np.ndarray:
        """
        [辅助] 基于当前流量 W_old 计算焓的时间导数 dh/dt (预测值)
        用于估算热膨胀源项 S_thermal

        方程: M * dh/dt = Sum(W*h)_in - Sum(W*h)_out + Q_total
        """
        dh_dt_vec = self.enthalpy_derivative_buffer
        dh_dt_vec.fill(0.0)

        if self.n_junc > 0:
            donor_idx = self.enthalpy_donor_idx_buffer
            np.copyto(donor_idx, self.idx_from_vec)
            reverse_mask = self.W_vec < 0.0
            donor_idx[reverse_mask] = self.idx_to_vec[reverse_mask]

            h_flux_vec = self.enthalpy_flux_buffer
            np.multiply(self.W_vec, self.h_vec[donor_idx], out=h_flux_vec)
            np.multiply(-self.M_from_vec, h_flux_vec, out=self.enthalpy_from_flux_buffer)
            np.multiply(self.M_to_vec, h_flux_vec, out=self.enthalpy_to_flux_buffer)
            np.add.at(dh_dt_vec, self.idx_from_vec, self.enthalpy_from_flux_buffer)
            np.add.at(dh_dt_vec, self.idx_to_vec, self.enthalpy_to_flux_buffer)

        # 2. 加上源项并除以质量
        # 获取质量 M = rho * V (直接读取缓存数组)
        mass_vec = self.enthalpy_mass_buffer
        np.multiply(self.rho_vec, self.V_vec, out=mass_vec)

        # 合成净热源: Q_total = Q_explicit - imp_coeff * T_curr
        Q_total_vec = self.enthalpy_q_total_buffer
        np.multiply(self.lam_imp_vec, self.T_vec, out=Q_total_vec)
        np.subtract(self.Q_expl_vec, Q_total_vec, out=Q_total_vec)

        # 将热源累加到导数上
        dh_dt_vec += Q_total_vec

        # 防护除零，进行向量化安全相除
        # 相当于: if mass > 1e-6: dh_dt = val / mass else: dh_dt = 0.0
        np.divide(dh_dt_vec, mass_vec, out=dh_dt_vec, where=mass_vec > 1e-6)
        dh_dt_vec[mass_vec <= 1e-6] = 0.0

        return dh_dt_vec

        # for i, vol in enumerate(self.volumes_obj):
        #     # 获取质量 M = rho * V
        #     mass = self.rho_vec[i] * getattr(vol, 'vol', 1.0)
        #
        #     # A. 获取显式部分
        #     # Q_wall 在 Components.py 中被累加为: explicit_part (即 h * A * T_wall)
        #     # Q_vol 是体积热源
        #     Q_explicit = getattr(vol, 'Q_wall', 0.0) + getattr(vol, 'Q_vol', 0.0)
        #
        #     # B. 获取隐式系数
        #     # implicit_coeff 在 Components.py 中被累加为: implicit_factor (即 h * A)
        #     imp_coeff = getattr(vol, 'implicit_coeff', 0.0)
        #
        #     # C. 获取当前流体温度
        #     # 使用 SystemManager 维护的状态向量 self.T_vec[i]，保证与当前时间步一致
        #     T_curr = self.T_vec[i]
        #
        #     # D. 合成净热源
        #     # 公式: Q_net = (hA * Tw) - (hA) * Tf
        #     # 这引入了负反馈：Tf 升高会导致 Q_total 减小，从而抑制温度爆炸
        #     Q_total = Q_explicit - imp_coeff * T_curr
        #
        #     # 获取热源 Q (Wall + Volumetric) （错误方法，不调用）
        #     # 假设 FluidVolume 有 Q_wall, Q_vol 属性
        #     # Q_total = getattr(vol, 'Q_wall', 0.0) + getattr(vol, 'Q_vol', 0.0)
        #
        #     # 计算最终 dh/dt
        #     # 保护 mass 防止除零
        #     if mass > 1e-6:
        #         dh_dt_vec[i] = (dh_dt_vec[i] + Q_total) / mass
        #     else:
        #         dh_dt_vec[i] = 0.0

        # return dh_dt_vec

    def _calc_thermal_expansion_source(self, dt: float) -> np.ndarray:
        """
        计算热膨胀源项 S_thermal [kg/s] (归一化到质量方程右端项)

        S_thermal = - V * (d_rho/dt)
                  = - V * (d_rho/dh * dh/dt)
        """

        # 获取 dh/dt
        dh_dt_vec = self._calc_enthalpy_time_derivative_explicit()

        # 2. 纯数组计算 d_rho/dh = (d_rho/dT) / Cp
        # 保护 cp 防止除零
        safe_cp = self.energy_cp_buffer
        np.maximum(self.cp_vec, 1.0e-3, out=safe_cp)
        np.divide(self.drho_dt_vec, safe_cp, out=self.drho_dh_buffer)
        drho_dh_vec = self.drho_dh_buffer

        # 3. 计算 drho/dt
        drho_dt_vec = self.energy_aux_buffer
        np.multiply(drho_dh_vec, dh_dt_vec, out=drho_dt_vec)

        # 4. 计算最终源项 S = - V * drho_dt
        S_vec = self.thermal_source_buffer
        np.multiply(-self.V_vec, drho_dt_vec, out=S_vec)

        return S_vec

        # S_vec = np.zeros(self.n_vol)
        #
        # # 1. 估算 dh/dt
        # dh_dt_vec = self._calc_enthalpy_time_derivative_explicit()
        #
        # for i, vol in enumerate(self.volumes_obj):
        #     # 2. 获取 d_rho/dh
        #     # d_rho/dh = (d_rho/dT) / Cp
        #     drho_dh = 0.0
        #     if hasattr(vol, 'material') and vol.material is not None:
        #         try:
        #             # 获取当前温度下的物性导数
        #             T = self.T_vec[i]
        #             P = self.P_vec[i]
        #             drho_dT = vol.material.liquid_density_derivative_T(T)
        #             cp = vol.material.heat_capacity(T, P)
        #             if cp > 1e-3:
        #                 drho_dh = drho_dT / cp
        #         except AttributeError:
        #             pass
        #
        #     # 3. 计算源项 S = - V * drho/dt
        #     # drho/dt = drho_dh * dh_dt
        #     # 物理检查: 加热 -> dh/dt > 0. 通常 drho/dT < 0 (膨胀).
        #     # -> drho/dt < 0 (密度降低).
        #     # -> S = - V * (-val) > 0.
        #     # -> 正源项，相当于节点内生成了体积，会把流体挤出去。符合物理。
        #
        #     V_cell = getattr(vol, 'vol', 1.0)
        #     drho_dt = drho_dh * dh_dt_vec[i]
        #
        #     # 注意: 压力方程构建是基于 d(rho*V)/dt + ... = 0
        #     # 离散化: V * (rho_new - rho_old)/dt + ... = 0
        #     # V * (drho_dP * dP + drho_dh * dh) / dt ...
        #     # (V/dt * drho_dP) * dP = - (V/dt * drho_dh * dh) ...
        #     # 所以右端项增加的是: - V * drho_dh * (dh_new - dh_old)/dt
        #     # 即 - V * drho_dt
        #
        #     S_vec[i] = - V_cell * drho_dt
        #
        # return S_vec

    @TEASAProfiler.profile
    def _step_energy_implicit(self, dt: float):
        """
        [Phase 4] 全隐式能量方程 (Enthalpy-Based Fully Implicit Solver)
        采用纯焓构建雅可比矩阵，彻底打破对流 CFL 限制，严格保证能量守恒。
        """

        # ==========================================
        # 1. 向量化预计算对角线和右端项系数 (纯 NumPy 数组运算)
        # ==========================================
        # 质量保护: M = max(rho * V, 1e-9)
        mass_vec = self.energy_mass_buffer
        np.multiply(self.rho_vec, self.V_vec, out=mass_vec)
        np.maximum(mass_vec, 1e-9, out=mass_vec)

        # Cp 保护: 从预计算缓存中读取
        cp_vec = self.energy_cp_buffer
        np.maximum(self.cp_vec, 1.0, out=cp_vec)

        # 隐式抗震荡系数: lam_h = lam / cp
        lam_h_vec = self.energy_lam_h_buffer
        np.divide(self.lam_imp_vec, cp_vec, out=lam_h_vec)

        # 修正后的显式源项: Q_mod = Q_expl - lam * T_old + lam_h * h_old
        Q_expl_mod_vec = self.energy_q_mod_buffer
        np.multiply(self.lam_imp_vec, self.T_vec, out=Q_expl_mod_vec)
        np.subtract(self.Q_expl_vec, Q_expl_mod_vec, out=Q_expl_mod_vec)
        np.multiply(lam_h_vec, self.h_vec, out=self.energy_aux_buffer)
        Q_expl_mod_vec += self.energy_aux_buffer

        # 惯性项系数: inertia = M / dt
        inertia_vec = self.energy_inertia_buffer
        np.divide(mass_vec, max(dt, 1.0e-12), out=inertia_vec)

        diag_data = self.energy_diag_buffer
        np.add(inertia_vec, lam_h_vec, out=diag_data)
        if self.fixed_pressure_idx_arr.size > 0:
            diag_data[self.fixed_pressure_idx_arr] = 1.0

        B_enth = self.energy_rhs_buffer
        np.multiply(inertia_vec, self.h_vec, out=B_enth)
        B_enth += Q_expl_mod_vec
        if self.fixed_pressure_idx_arr.size > 0:
            B_enth[self.fixed_pressure_idx_arr] = self.h_vec[self.fixed_pressure_idx_arr]

        data = self.energy_matrix.data
        data.fill(0.0)
        data[self.energy_diag_ptrs] = diag_data

        if self.energy_from_junc_idx.size > 0:
            W_from = self.W_vec[self.energy_from_junc_idx]
            from_pos = W_from > 0.0
            from_neg = W_from < 0.0
            if np.any(from_pos):
                np.add.at(
                    data,
                    self.energy_from_diag_ptrs[from_pos],
                    self.energy_from_multiplier[from_pos] * W_from[from_pos]
                )
            if np.any(from_neg):
                np.add.at(
                    data,
                    self.energy_from_offdiag_ptrs[from_neg],
                    self.energy_from_multiplier[from_neg] * W_from[from_neg]
                )

        if self.energy_to_junc_idx.size > 0:
            W_to = self.W_vec[self.energy_to_junc_idx]
            to_pos = W_to > 0.0
            to_neg = W_to < 0.0
            if np.any(to_pos):
                np.add.at(
                    data,
                    self.energy_to_offdiag_ptrs[to_pos],
                    -self.energy_to_multiplier[to_pos] * W_to[to_pos]
                )
            if np.any(to_neg):
                np.add.at(
                    data,
                    self.energy_to_diag_ptrs[to_neg],
                    -self.energy_to_multiplier[to_neg] * W_to[to_neg]
                )

        try:
            h_new_vec = spsolve(self.energy_matrix, B_enth)
        except Exception as e:
            logger.error(f"Enthalpy Matrix Solve Failed (Singular Matrix). Error: {e}")
            raise RuntimeError(f"Energy system solver failed: {e}")

        # ==========================================
        # 5. 更新状态并反算真实温度 (绝对核心提速点)
        # [优化]: 利用物性库的向量化接口，单次底层调用反算全部温度
        # ==========================================
        np.copyto(self.energy_T_old_buffer, self.T_vec)
        np.copyto(self.energy_h_old_buffer, self.h_vec)

        # 1. 直接全量更新缓存向量
        self.h_vec[:] = h_new_vec

        # 2. 批量调用 EoS 反算温度，如果调用失败才使用备用的 Cp 降级计算
        mat = self.primary_material

        if mat is not None and self.material_has_temp_from_h:
            try:
                # 单次 NumPy 向量调用，无需 for 循环和 try-except 损耗
                self.T_vec[:] = mat.temperature_from_enthalpy(h_new_vec, self.P_vec)
            except Exception as e:
                logger.warning(
                    f"Vectorized EoS temperature back-calculation failed: {e}. Falling back to Cp method.")
                self.T_vec[:] = self.energy_T_old_buffer + (h_new_vec - self.energy_h_old_buffer) / cp_vec
        else:
            self.T_vec[:] = self.energy_T_old_buffer + (h_new_vec - self.energy_h_old_buffer) / cp_vec

    def _step_energy(self, dt: float):
        """
        [Phase 4 半隐式版本] 推进能量方程 (Temperature-Based Implicit Solver)

        策略:
        1. 求解温度 T (解决源项刚性问题)
        2. 强制同步焓 h (解决物性拟合误差问题)
        """
        # 1. 计算显式部分的焓通量 (Advection)
        #    Flux_net = Sum(W*h)_in - Sum(W*h)_out
        flux_net_vec = np.zeros(self.n_vol)

        for j_idx, (_, idx_in, idx_out) in enumerate(self.junction_descriptors):
            W = self.W_vec[j_idx]

            # Upwind scheme: 输运的是“焓”
            if W >= 0:
                h_flow = W * self.h_vec[idx_in]
            else:
                h_flow = W * self.h_vec[idx_out]

            flux_net_vec[idx_in] -= h_flow
            flux_net_vec[idx_out] += h_flow

        # 2. 逐节点隐式求解温度
        for i, vol in enumerate(self.volumes_obj):
            mass = self.rho_vec[i] * getattr(vol, 'vol', 1.0)

            # 获取 Cp (用于构建温度方程系数)
            cp = 1000.0
            if hasattr(vol, 'material') and vol.material is not None:
                try:
                    cp = vol.material.heat_capacity(self.T_vec[i], self.P_vec[i])
                except:
                    pass

            if mass < 1e-9:
                continue

            # A. 构建方程系数: (M*Cp/dt + lambda) * T_new = ...
            inertia_coef = (mass * cp) / dt

            Q_expl = getattr(vol, 'Q_wall', 0.0) + getattr(vol, 'Q_vol', 0.0)
            lam = getattr(vol, 'implicit_coeff', 0.0)

            # B. 求解 T_new (无条件稳定)
            rhs = inertia_coef * self.T_vec[i] + flux_net_vec[i] + Q_expl
            lhs_inv = 1.0 / (inertia_coef + lam)
            T_new = rhs * lhs_inv

            # C. 更新温度 (这是我们的“真理”)
            self.T_vec[i] = T_new
            vol.T = T_new

            # D. 同步更新焓 h (消除拟合误差)
            # 优先使用物性库的绝对焓函数，保证 h 和 T 始终在物性曲线上
            if hasattr(vol, 'material') and hasattr(vol.material, 'enthalpy'):
                try:
                    # 直接由 T 算 h，避免反算误差
                    h_new = vol.material.enthalpy(T_new, self.P_vec[i])
                    self.h_vec[i] = h_new
                    vol.h = h_new
                except:
                    # Fallback: 如果调用失败，使用 Cp 增量
                    self.h_vec[i] += cp * (T_new - self.T_vec[i])  # 注意此时 T_vec[i] 已经被覆盖了，这里逻辑需微调
                    # 修正 Fallback 逻辑: 无法获取旧 T，只能用 differential 近似
                    # 但通常只要有 material 就能调通，这里仅做极其保守的防崩溃处理
                    pass
            else:
                # Fallback: 如果没有 enthalpy 接口，使用增量更新
                # Q_net_actual = Flux + Q - lam * T
                Q_net_actual = flux_net_vec[i] + Q_expl - lam * T_new
                dh = (Q_net_actual / mass) * dt
                self.h_vec[i] += dh
                vol.h = self.h_vec[i]

    # def _step_energy(self, dt: float):
    #     """
    #     [Phase 4] 显式推进能量方程 (Corrector)
    #     使用更新后的流量 W^{n+1} 更新焓 h^{n+1}
    #     """
    #     # 计算导数 (使用更新后的 W_vec)
    #     # 注意: 此时 W_vec 已经是 W_new 了
    #     dh_dt_vec = self._calc_enthalpy_time_derivative_explicit()
    #
    #     # 显式积分
    #     self.h_vec += dh_dt_vec * dt
    #
    #     # 同步更新温度 T 和 对象属性
    #     for i, vol in enumerate(self.volumes_obj):
    #         vol.h = self.h_vec[i]
    #
    #         # 反算温度
    #         if hasattr(vol, 'material') and vol.material is not None:
    #             # T = T(h, P)
    #             try:
    #                 new_T = vol.material.temperature_from_enthalpy(vol.h, vol.P)
    #                 self.T_vec[i] = new_T
    #                 vol.T = new_T
    #             except AttributeError:
    #                 pass
    #         else:
    #             # 如果没有物性库，简单近似 dh = Cp * dT (假设 Cp=1000)
    #             vol.T += dh_dt_vec[i] * dt / 1000.0
    #             self.T_vec[i] = vol.T

    def step(self, dt: float):
        """
        [主接口] 执行一个完整的时间步 (Hydraulic + Thermal)

        Sequence:
        1. Prop Update: P_old, h_old -> rho, mu
        2. Momentum Coeffs: W_old -> a_j, b_j (Matrix coeffs)
        3. Pressure Solve: Assemble A, B (with S_thermal) -> Solve P_new
        4. Flow Update: P_new -> W_new
        5. Energy Step: W_new, h_old -> h_new
        """
        # 1. 更新物性
        self._update_fluid_properties()

        # 2. 计算动量系数
        self._calc_momentum_coeffs(dt)

        # 3. 组装并求解压力 (包含热膨胀预测)
        A, B = self._assemble_pressure_system(dt)
        self._solve_linear_system(A, B)

        # 4. 更新流量
        self._update_flow_rates()
        self._sync_vectors_to_objects(
            sync_pressure=True,
            sync_flow=True,
            sync_energy=False,
            sync_properties=False
        )

        # 5. 更新能量 (焓/温度)
        self._step_energy(dt)
        self._sync_vectors_to_objects(
            sync_pressure=False,
            sync_flow=False,
            sync_energy=True,
            sync_properties=False
        )

    @TEASAProfiler.profile
    def step_Picard(self, dt: float, max_iter: int = 20, tol: float = 1e-4) -> bool:
        """
        [主接口] 执行一个完整的时间步 (Hydraulic + Thermal)
        包含 Picard 外部迭代，以解决非线性阻力系数问题。

        :param dt: 时间步长 [s]
        :param max_iter: 最大水力迭代次数
        :param tol: 流量收敛容差 [kg/s]
        :return: 是否收敛 (True/False)
        """
        # 1. 准备步骤：更新物性 (基于 t^n 的温度)
        #    注意：通常物性在时间步内变化较慢，可以放在迭代循环外
        self._update_fluid_properties()

        converged = False
        S_thermal_vec = self._calc_thermal_expansion_source(dt)

        # --- 2. 外部迭代循环 (Hydraulic Iteration) ---
        # 目标：解决 W 与 阻力系数(W) 之间的非线性耦合

        # 进入迭代前，冻结历史流量
        np.copyto(self.W_old, self.W_vec)

        diff = 1e6

        for k in range(max_iter):
            # 2.1 备份当前流量 (用于收敛检查)
            np.copyto(self.W_iterate, self.W_vec)

            # 2.2 计算动量系数 (Phase 2)
            #     这里会使用最新的 self.W_vec 来计算 Re 和 f
            #     如果是第0次迭代，使用的是 W^n
            #     如果是第k次迭代，使用的是 W^{n+1, k-1}
            # 传入冻结的历史流量，保证惯性项守恒
            self._calc_momentum_coeffs(dt, W_old_frozen=self.W_old)

            # 2.3 组装压力矩阵 (Phase 3)
            #     时间步内冻结热膨胀源项，避免 Picard 内重复计算
            A, B = self._assemble_pressure_system(
                dt,
                include_thermal_expansion=False,
                S_thermal_vec=S_thermal_vec
            )

            # 2.4 求解压力场
            try:
                self._solve_linear_system(A, B)
            except np.linalg.LinAlgError:
                logger.error(f"Solver failed at step t={dt}")
                return False

            # 2.5 更新流量 (回代)
            self._update_flow_rates()

            # 2.6 收敛性检查
            #     检查流量变化的最大范数
            np.subtract(self.W_vec, self.W_iterate, out=self.W_residual)
            np.abs(self.W_residual, out=self.W_residual)
            diff = float(np.max(self.W_residual))

            logger.debug(f"Iter {k}: Max dW = {diff:.6e}")

            if diff < tol:
                converged = True
                logger.debug(f"Hydraulic step converged in {k+1} iterations.")
                break

            # self._step_energy_implicit(dt)

        if not converged:
            logger.warning(f"Hydraulic step NOT converged after {max_iter} iterations. Max residual: {diff:.6e}")

        # --- 3. 能量方程步进 (Corrector) ---
        # 只有当流场收敛(或尽力收敛)后，才推进温度
        # 这样保证了能量输运使用的是正确的流速
        # self._step_energy(dt)
        self._step_energy_implicit(dt)
        self._sync_vectors_to_objects(
            sync_pressure=True,
            sync_flow=True,
            sync_energy=True,
            sync_properties=False
        )

        return converged

    def step_Picard_over(self, dt: float, max_iter: int = 20, tol: float = 1e-4) -> bool:
        """
        [主接口] 执行一个完整的时间步 (Hydraulic + Thermal)
        包含 Picard 外部迭代，以解决非线性阻力系数问题。

        :param dt: 时间步长 [s]
        :param max_iter: 最大水力迭代次数
        :param tol: 流量收敛容差 [kg/s]
        :return: 是否收敛 (True/False)
        """
        # 1. 准备步骤：更新物性 (基于 t^n 的温度)
        #    注意：通常物性在时间步内变化较慢，可以放在迭代循环外
        self._update_fluid_properties()

        converged = False
        # converged = True
        S_thermal_vec = self._calc_thermal_expansion_source(dt)

        # --- 2. 外部迭代循环 (Hydraulic Iteration) ---
        # 目标：解决 W 与 阻力系数(W) 之间的非线性耦合

        np.copyto(self.W_old, self.W_vec)
        diff = 1e6

        for k in range(max_iter):
            # 2.1 备份当前流量 (用于收敛检查)
            np.copyto(self.W_iterate, self.W_vec)

            # 2.2 计算动量系数 (Phase 2)
            #     这里会使用最新的 self.W_vec 来计算 Re 和 f
            #     如果是第0次迭代，使用的是 W^n
            #     如果是第k次迭代，使用的是 W^{n+1, k-1}
            # 传入冻结的历史流量，保证惯性项守恒
            self._calc_momentum_coeffs(dt, W_old_frozen=self.W_old)

            # 2.3 组装压力矩阵 (Phase 3)
            #     时间步内冻结热膨胀源项，避免 Picard 内重复计算
            A, B = self._assemble_pressure_system(
                dt,
                include_thermal_expansion=False,
                S_thermal_vec=S_thermal_vec
            )

            # 2.4 求解压力场
            try:
                self._solve_linear_system(A, B)
            except np.linalg.LinAlgError:
                logger.error(f"Solver failed at step t={dt}")
                return False

            # 2.5 更新流量 (回代)
            self._update_flow_rates()

            # 2.6 收敛性检查
            #     检查流量变化的最大范数
            np.subtract(self.W_vec, self.W_iterate, out=self.W_residual)
            np.abs(self.W_residual, out=self.W_residual)
            diff = float(np.max(self.W_residual))

            logger.debug(f"Iter {k}: Max dW = {diff:.6e}")

            if diff < tol:
                converged = True
                logger.debug(f"Hydraulic step converged in {k+1} iterations.")
                break

            # self._step_energy_implicit(dt)

        if not converged:
            logger.warning(f"Hydraulic step NOT converged after {max_iter} iterations. Max residual: {diff:.6e}")

        # --- 3. 能量方程步进 (Corrector) ---
        # 只有当流场收敛(或尽力收敛)后，才推进温度
        # 这样保证了能量输运使用的是正确的流速
        # self._step_energy(dt)
        self._step_energy_implicit(dt)
        self._sync_vectors_to_objects(
            sync_pressure=True,
            sync_flow=True,
            sync_energy=True,
            sync_properties=False
        )

        return converged

    def initialize_hydraulics(self, dt: float = 0.1, tol: float = 1e-5, max_iter: int = 500, omega: float = 0.5) -> bool:
        """
        [初始化工具] 纯液力稳态初始化 (Hydraulic Only Initialization)

        功能：
        在冻结温度场(不更新能量)的情况下，迭代求解压力和流量，直到达到稳态。
        这用于消除初始条件的非物理震荡，建立符合当前阻力特性的初始流场。

        :param dt: 虚拟时间步长 [s]。建议取 0.1~1.0。
                   dt 越小，阻尼越大，越稳定但收敛越慢；
                   dt 越大，惯性项越小，越接近直接牛顿法，但也越容易震荡。
        :param tol: 收敛容差 [kg/s]
        :param max_iter: 最大迭代次数
        :return: 是否收敛
        """
        logger.info(f"Starting Hydraulic Initialization (Max Iter={max_iter}, Tol={tol})...")

        # 1. 初始物性更新 (基于初始温度场)
        self._update_fluid_properties()

        converged = False

        for k in range(max_iter):
            # 备份旧流量用于检查收敛
            np.copyto(self.W_iterate, self.W_vec)

            # 2. 计算动量系数
            #    注意：这里 dt 起到了松弛因子(Relaxation)的作用
            #    I = L/Adt。当收敛时 W_new = W_old，惯性项自然消失。
            self._calc_momentum_coeffs(dt, W_old_frozen=self.W_iterate)

            # 3. 组装压力矩阵 关闭热膨胀源项
            A, B = self._assemble_pressure_system(dt, include_thermal_expansion=False)

            # 4. 求解压力
            try:
                self._solve_linear_system(A, B)
            except np.linalg.LinAlgError:
                logger.error(f"Initialization solver failed at iter={k}")
                return False

            # 5. 更新流量
            self._update_flow_rates()

            # 6. 选用显式亚松驰（稳态计算）
            self.W_vec[:] = omega * self.W_vec + (1.0 - omega) * self.W_iterate

            # 6. 检查收敛 (L_inf norm)
            np.subtract(self.W_vec, self.W_iterate, out=self.W_residual)
            np.abs(self.W_residual, out=self.W_residual)
            diff = float(np.max(self.W_residual))

            if diff < tol:
                converged = True
                logger.info(f"Hydraulic initialization CONVERGED at iter {k + 1}. Residual={diff:.2e}")
                break

        if not converged:
            logger.warning(f"Initial: Hydraulic initialization NOT converged after {max_iter} iters. Residual={diff:.2e}")

        self._sync_vectors_to_objects(
            sync_pressure=True,
            sync_flow=True,
            sync_energy=False,
            sync_properties=False
        )

        return converged

    # =========================================================================
    # Phase 6: Time Step Control (Adaptive DT)
    # =========================================================================
    def get_max_stable_dt(self, max_limit: float = 0.5) -> float:
        """
        [自适应步长] 基于 CFL 条件计算最大允许时间步长
        Criterion: dt < L_node / v_max

        :param max_limit: 当流速接近0时的最大时间步上限 [s]
        :return: 建议的时间步长 (未乘安全系数)
        """
        min_dt = max_limit

        for vol in self.volumes_obj:
            # 1. 获取节点特征长度
            # 注意: BoundaryVolume (如Plenum) 通常 length 很大或未定义，
            # 且流速极低，通常不是 CFL 瓶颈。这里做个防御性检查。
            L = getattr(vol, 'len', 0.0)

            # 忽略长度为0或极小的虚拟节点
            if L < 1e-6:
                continue

            # 2. 获取节点特征流速
            # 保守策略: 取连接到该节点的所有接口中流速最大的那个
            v_max = 0.0

            # 检查入口
            if hasattr(vol, 'inlet_junctions'):
                for j in vol.inlet_junctions:
                    v_max = max(v_max, abs(j.vel))

            # 检查出口
            if hasattr(vol, 'outlet_junctions'):
                for j in vol.outlet_junctions:
                    v_max = max(v_max, abs(j.vel))

            # 3. 计算 CFL 限制
            # dt < L / v
            if v_max > 1e-6:
                dt_local = L / v_max
                if dt_local < min_dt:
                    min_dt = dt_local
            else:
                # 流速为 0 时，不受 CFL 限制，维持当前的 min_dt (即 max_limit)
                pass

        return min_dt

    def save_state(self):
        self.P_backup = self.P_vec.copy()
        self.h_backup = self.h_vec.copy()
        self.W_backup = self.W_vec.copy()
        self.T_backup = self.T_vec.copy()
        self.rho_backup = self.rho_vec.copy()
        self.mu_backup = self.mu_vec.copy()
        # ... 以及 T_vec, rho_vec 等

    def load_state(self):
        """
        [接口] 恢复到上一次保存的状态
        """
        # 1. 恢复向量数据
        self.P_vec[:] = self.P_backup
        self.T_vec[:] = self.T_backup
        self.h_vec[:] = self.h_backup
        self.W_vec[:] = self.W_backup
        self.rho_vec[:] = self.rho_backup
        self.mu_vec[:] = self.mu_backup

        # 2. 将恢复后的向量同步回对象实例
        self._sync_vectors_to_objects()

    def _sync_vectors_to_objects(
        self,
        sync_pressure: bool = True,
        sync_flow: bool = True,
        sync_energy: bool = True,
        sync_properties: bool = True
    ):
        """
        [辅助] 将求解器向量同步回物理对象 (Volume, Junction)
        确保对象属性 (vol.P, junc.W 等) 与求解器向量保持一致。
        """
        # 1. 同步节点状态 (Volume)
        for i, vol in enumerate(self.volumes_obj):
            if sync_pressure:
                vol.P = self.P_vec[i]
            if sync_energy:
                vol.T = self.T_vec[i]
                vol.h = self.h_vec[i]
            if sync_properties:
                vol.rho = self.rho_vec[i]
                vol.mu = self.mu_vec[i]

        # 2. 同步连接状态 (Junction)
        if sync_flow:
            for i, junc in enumerate(self.junctions_obj):
                junc.W = self.W_vec[i]

                # [关键] 必须同时更新流速 (vel)，因为阻力计算依赖 junc.vel
                # 如果只恢复 W 而不更新 vel，下一次阻力计算可能会用错速度
                if hasattr(junc, 'update_velocity'):
                    junc.update_velocity()

    # =========================================================================
    # Phase 7: 断点续算与状态持久化 (Restart & Persistence)
    # =========================================================================

    def get_state_dict(self, prefix: str) -> dict:
        """
        [持久化接口] 将流体网络状态序列化为扁平字典

        策略：主体采用高性能的全局向量存储，辅以边界控制体的局部参数存储。
        """
        state = {
            # 1. 拓扑指纹 (防御性校验用: [n_vol, n_junc])
            f"{prefix}/shape": np.array([self.n_vol, self.n_junc]),
            f"{prefix}/fixed_pressure_idx": np.array(sorted(self.fixed_pressure_indices), dtype=np.int32),

            # 2. 核心状态向量 (极速高压缩比)
            f"{prefix}/P_vec": self.P_vec,
            f"{prefix}/T_vec": self.T_vec,
            f"{prefix}/h_vec": self.h_vec,
            f"{prefix}/W_vec": self.W_vec,
        }

        # 3. 提取具有瞬态记忆的边界控制目标 (如阀门、稳压器的目标值)
        # 注意：只提取存在的属性，转为 1D NumPy 数组以兼容 .npz
        for i, vol in enumerate(self.volumes_obj):
            if hasattr(vol, 'target_P'):
                state[f"{prefix}/Vols/{i}/target_P"] = np.array([vol.target_P], dtype=float)

        for j, junc in enumerate(self.junctions_obj):
            if hasattr(junc, 'target_W'):
                state[f"{prefix}/Juncs/{j}/target_W"] = np.array([junc.target_W], dtype=float)

        return state

    def load_state_dict(self, data: dict, prefix: str):
        """
        [持久化接口] 从字典恢复流体网络状态
        """
        # --- 1. 致命级拓扑校验 (Topology Validation) ---
        shape_key = f"{prefix}/shape"
        if shape_key not in data:
            raise KeyError(f"Missing topology fingerprint for HydraulicNetwork in restart file.")

        saved_shape = data[shape_key]
        if saved_shape[0] != self.n_vol or saved_shape[1] != self.n_junc:
            raise ValueError(
                f"CRITICAL: Topology mismatch! Restart file has ({saved_shape[0]} Vols, {saved_shape[1]} Juncs), "
                f"but current memory has ({self.n_vol} Vols, {self.n_junc} Juncs).")

        fixed_idx_key = f"{prefix}/fixed_pressure_idx"
        if fixed_idx_key in data:
            saved_fixed_idx = np.asarray(data[fixed_idx_key], dtype=np.int32).reshape(-1)
            current_fixed_idx = np.array(sorted(self.fixed_pressure_indices), dtype=np.int32)
            if saved_fixed_idx.shape != current_fixed_idx.shape or not np.array_equal(saved_fixed_idx, current_fixed_idx):
                raise ValueError(
                    "CRITICAL: Fixed-pressure boundary set mismatch! "
                    f"Restart file has {saved_fixed_idx.tolist()}, "
                    f"but current model has {current_fixed_idx.tolist()}. "
                    "Changing fixed-pressure boundary membership across restart is not supported "
                    "with the cached matrix structures."
                )
        else:
            logger.warning(
                "Restart file is missing fixed-pressure boundary fingerprint. "
                "Skipping fixed-boundary-set validation for backward compatibility."
            )

        # --- 2. 恢复核心状态向量 ---
        self.P_vec[:] = data[f"{prefix}/P_vec"]
        self.T_vec[:] = data[f"{prefix}/T_vec"]
        self.h_vec[:] = data[f"{prefix}/h_vec"]
        self.W_vec[:] = data[f"{prefix}/W_vec"]

        # --- 3. 恢复边界控制目标 ---
        for i, vol in enumerate(self.volumes_obj):
            t_p_key = f"{prefix}/Vols/{i}/target_P"
            if t_p_key in data and hasattr(vol, 'target_P'):
                vol.target_P = float(data[t_p_key][0])

        for j, junc in enumerate(self.junctions_obj):
            t_w_key = f"{prefix}/Juncs/{j}/target_W"
            if t_w_key in data and hasattr(junc, 'target_W'):
                junc.target_W = float(data[t_w_key][0])

        self._refresh_cached_pressure_targets()
        self._refresh_cached_boundary_targets()

        # --- 4. 强制状态同步与物性刷新 ---
        # 必须把恢复到全局向量的数值同步给底层的 Physical Objects
        self._sync_vectors_to_objects()

        # 必须立即刷新一遍物性 (rho_vec, mu_vec 等)，否则接下来的第一步动量计算会基于错误的密度崩溃
        self._update_fluid_properties()
