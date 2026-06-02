import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np

from Components.BaseComponent import BaseComponent
from Components.TFEUnit import TFEUnit
from Components.tec_electric import electric_field_from_node_potential
from Solvers.Neutronics.PointReactor import PointReactor
from Solvers.Couplers import GapCouple2D, SolidSolidCouple2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

logger = logging.getLogger(__name__)


@dataclass
class GlobalAnnulusStructureConfig:
    """
    堆芯外部全局环形固体层配置。

    适用对象：
    1. 筒体层
    2. 反射层

    设计原则：
    1. 优先允许用户直接传入现成的 ``Mesh2D``
    2. 若未传入网格，则由 ``ReactorCore`` 根据半径参数自动生成空心圆柱网格
    3. 轴向网格默认与 TFE 组件完全一致，便于后续温度同步和反应性反馈统计
    """

    material: Any
    outer_radius: Optional[float] = None
    thickness: Optional[float] = None
    inner_radius: Optional[float] = None
    n_radial: int = 8
    initial_temp: float = 743.0
    outer_surface_emissivity: Optional[float] = None
    mesh: Optional[Mesh2D] = None


@dataclass
class GlobalGapStructureConfig:
    """
    堆芯外部环形间隙配置。

    支持两种模式：
    1. ``simplified``: 使用 ``GapCouple2D``，适合先快速搭好结构
    2. ``meshed``:     将间隙本身也建成独立的 ``HeatConduction2D``，便于后续扩展
    """

    mode: Literal['simplified', 'meshed'] = 'simplified'
    width: float = 0.0
    h_eq: float = 0.0
    emissivity_inner: float = 0.8
    emissivity_outer: float = 0.8
    material: Optional[Any] = None
    n_radial: int = 1
    initial_temp: float = 743.0
    mesh: Optional[Mesh2D] = None


@dataclass
class ReactivityFeedbackTemperatureSummary:
    """
    反应性反馈用的结构平均温度摘要。

    所有温度默认单位均为 K。
    若某类结构当前尚未建模，则对应字段为 ``None``。
    """

    fuel: Optional[float] = None
    emitter: Optional[float] = None
    collector: Optional[float] = None
    moderator: Optional[float] = None
    reflector: Optional[float] = None
    barrel: Optional[float] = None


@dataclass
class ReactivityFeedbackResult:
    """
    反应性反馈计算结果摘要。

    ``total`` 为各项反馈之和，便于后续直接传给点堆模块。
    """

    temperatures: ReactivityFeedbackTemperatureSummary
    fuel: float = 0.0
    electrode: float = 0.0
    moderator: float = 0.0
    reflector: float = 0.0
    total: float = 0.0


class ReactorCore(BaseComponent):
    """
    TASTIN 堆芯宏观容器。

    当前这一层负责两件事：
    1. 组织所有 TFE 以及堆芯外部全局结构的热工建模
    2. 负责 TEC 相关的温度收集与热流下发

    说明：
    1. 中子学反馈框架后续会继续收敛到 ReactorCore 层，不再依赖各个 Fuel 单独返回反馈
    2. 本次先补齐“慢化剂-间隙-筒体-间隙-反射层”的全局结构建模
    """

    def __init__(self,
                 name: str,
                 tfe_dict: Dict[str, TFEUnit],
                 tfe_multipliers: Dict[str, int],
                 tec_multipliers: Optional[Dict[str, int]] = None,
                 tfe_power_factors: Optional[Dict[str, float]] = None,
                 mod_meshes: Optional[List[Mesh2D]] = None,
                 mod_material: Optional[Any] = None,
                 ring_mapping: Optional[Dict[str, int]] = None,
                 barrel_config: Optional[GlobalAnnulusStructureConfig] = None,
                 reflector_config: Optional[GlobalAnnulusStructureConfig] = None,
                 moderator_barrel_gap_config: Optional[GlobalGapStructureConfig] = None,
                 barrel_reflector_gap_config: Optional[GlobalGapStructureConfig] = None,
                 T_space: float = 250.0,
                 alpha_tec: float = 0.5,
                 enable_tec_coupled: bool = True):
        """
        初始化堆芯组件。

        参数说明：
        1. ``tfe_power_factors`` 用于定义代表元件对应的单根燃料元件功率份额
        2. ``mod_meshes/mod_material/ring_mapping`` 用于全局慢化剂环建模
        3. ``barrel_config/reflector_config`` 用于在慢化剂外继续建立筒体与反射层
        4. ``moderator_barrel_gap_config/barrel_reflector_gap_config`` 用于建立对应层间的环形间隙
        """
        super().__init__(name)

        self.tfes = tfe_dict
        self.tfe_multipliers = tfe_multipliers
        self.tec_multipliers = dict(tec_multipliers) if tec_multipliers is not None else dict(tfe_multipliers)
        self.alpha_tec = alpha_tec
        self.enable_tec_coupled = enable_tec_coupled
        self.ring_mapping = ring_mapping
        self.tfe_power_factors: Dict[str, float] = {}
        self.tfe_group_power_shares: Dict[str, float] = {}
        self.power_factor_weighted_sum = 0.0
        self.last_total_core_power = 0.0
        self.ring_to_tfe_names: Dict[int, List[str]] = {}
        self.ring_total_multipliers: Dict[int, int] = {}
        self.last_feedback_temperatures = ReactivityFeedbackTemperatureSummary()
        self.last_feedback_result = ReactivityFeedbackResult(
            temperatures=ReactivityFeedbackTemperatureSummary()
        )
        self.feedback_reference_result = ReactivityFeedbackResult(
            temperatures=ReactivityFeedbackTemperatureSummary()
        )
        self.last_effective_reactivity_feedback = 0.0
        self.last_reactivity_control = 0.0

        # 多物理场求解器句柄。
        self.point_reactor = None

        # 统一收集由 ReactorCore 自己构建出来的额外固体和耦合器。
        self.extra_solids: Dict[str, HeatConduction2D] = {}
        self.extra_couplers: Dict[str, Any] = {}

        # 为后续反馈与后处理预留的全局结构引用。
        self.has_global_moderator = False
        self.mod_rings: List[HeatConduction2D] = []
        self.barrel: Optional[HeatConduction2D] = None
        self.reflector: Optional[HeatConduction2D] = None
        self.global_gap_solids: Dict[str, HeatConduction2D] = {}
        self.global_outer_solid: Optional[HeatConduction2D] = None
        self.barrel_outer_surface_emissivity: Optional[float] = None
        self.reflector_outer_surface_emissivity: Optional[float] = None

        self._validate_inputs()
        self._initialize_power_factors(tfe_power_factors=tfe_power_factors)
        self._initialize_ring_groups(mod_meshes=mod_meshes)
        self._build_thermo_calc()
        self._build_global_moderator_rings(mod_meshes=mod_meshes, mod_material=mod_material)
        self._build_global_outer_structures(
            barrel_config=barrel_config,
            reflector_config=reflector_config,
            moderator_barrel_gap_config=moderator_barrel_gap_config,
            barrel_reflector_gap_config=barrel_reflector_gap_config
        )
        self._attach_outer_radiation_boundary(T_space=T_space)

        # TEC 更新频率控制。
        self.thermo_update_interval = 0.8
        self._last_thermo_update_time = -999.0

    # =========================================================================
    # 初始化辅助方法
    # =========================================================================

    def _validate_inputs(self):
        """检查 TFE 字典与复制系数字典的一致性。"""
        if not self.tfes:
            raise ValueError("ReactorCore 至少需要包含一个 TFEUnit。")

        for key in self.tfes.keys():
            if key not in self.tfe_multipliers:
                raise ValueError(f"TFE '{key}' 必须在 tfe_multipliers 中给出复制数量。")
            if key not in self.tec_multipliers:
                raise ValueError(f"TFE '{key}' 必须在 tec_multipliers 中给出复制数量。")

            tfe_multiplier = int(self.tfe_multipliers[key])
            tec_multiplier = int(self.tec_multipliers[key])
            if tfe_multiplier <= 0:
                raise ValueError(f"TFE '{key}' 的 tfe_multipliers 必须为正整数。")
            if tec_multiplier < 0:
                raise ValueError(f"TFE '{key}' 的 tec_multipliers 不能为负数。")
            if tec_multiplier > tfe_multiplier:
                raise ValueError(
                    f"TFE '{key}' 的 tec_multipliers ({tec_multiplier}) "
                    f"不能大于热工/水力 tfe_multipliers ({tfe_multiplier})。"
                )

        for key in self.tec_multipliers.keys():
            if key not in self.tfes:
                raise ValueError(f"tec_multipliers 中存在未知的 TFE 键: '{key}'。")

    def _initialize_power_factors(self, tfe_power_factors: Optional[Dict[str, float]]):
        """
        初始化代表元件功率因子。

        约定：
        1. ``tfe_power_factors[name]`` 表示该代表元件所代表的“单根真实燃料元件”功率份额
        2. 全堆验算采用 ``sum(power_factor * multiplier)``
        3. 若未显式提供，则默认所有真实燃料元件等功率
        """
        total_virtual_elements = float(sum(self.tfe_multipliers.values()))
        if total_virtual_elements <= 0.0:
            raise ValueError("tfe_multipliers 的总和必须大于 0。")

        if tfe_power_factors is None:
            uniform_factor = 1.0 / total_virtual_elements
            tfe_power_factors = {name: uniform_factor for name in self.tfes.keys()}

        self.tfe_power_factors = {}
        self.tfe_group_power_shares = {}

        for key in self.tfes.keys():
            if key not in tfe_power_factors:
                raise ValueError(f"TFE '{key}' 必须在 tfe_power_factors 中给出功率因子。")

        for key, factor in tfe_power_factors.items():
            if key not in self.tfes:
                raise ValueError(f"tfe_power_factors 中存在未知的 TFE 键: '{key}'。")

            factor_value = float(factor)
            if factor_value < 0.0:
                raise ValueError(f"TFE '{key}' 的功率因子不能为负数。")

            self.tfe_power_factors[key] = factor_value
            self.tfe_group_power_shares[key] = factor_value * float(self.tfe_multipliers[key])

        self.power_factor_weighted_sum = float(sum(self.tfe_group_power_shares.values()))
        if not np.isclose(self.power_factor_weighted_sum, 1.0, atol=1e-6, rtol=1e-6):
            raise ValueError(
                "功率因子验算失败：sum(tfe_power_factor * multiplier) "
                f"= {self.power_factor_weighted_sum:.12f}，应接近 1.0。"
            )

        logger.info(
            f"ReactorCore '{self.name}': Power-factor check passed, "
            f"sum(factor * multiplier) = {self.power_factor_weighted_sum:.12f}"
        )

    def _initialize_ring_groups(self, mod_meshes: Optional[List[Mesh2D]]):
        """
        初始化“全局慢化剂圈 -> 代表元件列表”的映射。

        设计目的：
        1. 显式支持同一圈存在多个代表元件
        2. 后续反应性反馈统计可以直接复用这一层分组关系
        3. 同圈统计时统一按 multiplier 加权，而不是默认取某一个代表件代替整圈
        """
        self.ring_to_tfe_names = {}
        self.ring_total_multipliers = {}

        if self.ring_mapping is None:
            return

        n_rings = None if mod_meshes is None else len(mod_meshes)

        for tfe_name in self.tfes.keys():
            if tfe_name not in self.ring_mapping:
                raise ValueError(f"TFE '{tfe_name}' 缺少 ring_mapping 定义。")

            ring_idx = int(self.ring_mapping[tfe_name])
            if ring_idx < 0:
                raise ValueError(f"TFE '{tfe_name}' 的 ring_mapping 不能为负数。")
            if n_rings is not None and ring_idx >= n_rings:
                raise ValueError(
                    f"TFE '{tfe_name}' 的 ring_mapping={ring_idx} 超出全局慢化剂圈数 {n_rings}。"
                )

            if ring_idx not in self.ring_to_tfe_names:
                self.ring_to_tfe_names[ring_idx] = []
                self.ring_total_multipliers[ring_idx] = 0

            self.ring_to_tfe_names[ring_idx].append(tfe_name)
            self.ring_total_multipliers[ring_idx] += int(self.tfe_multipliers[tfe_name])

        if self.ring_to_tfe_names:
            logger.info(
                f"ReactorCore '{self.name}': Initialized {len(self.ring_to_tfe_names)} moderator ring group(s)."
            )

    def get_ring_member_names(self, ring_idx: int) -> List[str]:
        """返回指定慢化剂圈层中的全部代表元件名称。"""
        return list(self.ring_to_tfe_names.get(int(ring_idx), []))

    def get_ring_total_multiplier(self, ring_idx: int) -> int:
        """返回指定慢化剂圈层所代表的真实燃料元件总数。"""
        return int(self.ring_total_multipliers.get(int(ring_idx), 0))

    def compute_ring_weighted_average(self, ring_idx: int, value_getter) -> Any:
        """
        计算同圈多个代表元件的 multiplier 加权平均值。

        参数：
        1. ``ring_idx``: 目标慢化剂圈编号
        2. ``value_getter``: 回调函数，签名为 ``value_getter(tfe_name, tfe) -> scalar/array``

        返回：
        1. 若该圈无代表元件，返回 ``None``
        2. 若有多个代表元件，则按 ``multiplier`` 做严格加权平均
        """
        member_names = self.get_ring_member_names(ring_idx)
        if not member_names:
            return None

        weighted_sum = None
        total_weight = 0.0

        for tfe_name in member_names:
            tfe = self.tfes[tfe_name]
            weight = float(self.tfe_multipliers[tfe_name])
            value = np.asarray(value_getter(tfe_name, tfe), dtype=float)

            if weighted_sum is None:
                weighted_sum = np.zeros_like(value, dtype=float)

            weighted_sum += value * weight
            total_weight += weight

        if total_weight <= 0.0:
            return None

        return weighted_sum / total_weight

    def _compute_representative_structure_average_temperature(self, structure_name: str) -> Optional[float]:
        """
        计算代表元件内部某类结构的全堆体积平均温度。

        计算规则：
        1. 只统计当前存在于 ``tfe.solids`` 中的目标结构
        2. 每个代表元件按 ``multiplier`` 还原成真实元件数量
        3. 每个代表元件内部采用严格体积加权平均
        """
        numerator = 0.0
        denominator = 0.0

        for tfe_name, tfe in self.tfes.items():
            if structure_name not in tfe.solids:
                continue

            solid = tfe.solids[structure_name]
            volumes = np.asarray(solid.mesh.geom_data.volumes, dtype=float)
            multiplier = float(self.tfe_multipliers[tfe_name])

            numerator += multiplier * float(np.sum(solid.T * volumes))
            denominator += multiplier * float(np.sum(volumes))

        if denominator <= 0.0:
            return None

        return numerator / denominator

    def _compute_global_structure_average_temperature(self, solid: Optional[HeatConduction2D]) -> Optional[float]:
        """
        计算全局结构的体积平均温度。

        适用对象：
        1. 全局慢化剂环
        2. 全局筒体
        3. 全局反射层
        """
        if solid is None:
            return None

        volumes = np.asarray(solid.mesh.geom_data.volumes, dtype=float)
        denominator = float(np.sum(volumes))
        if denominator <= 0.0:
            return None

        return float(np.sum(solid.T * volumes)) / denominator

    def _compute_global_moderator_average_temperature(self) -> Optional[float]:
        """
        计算全局慢化剂的体积平均温度。

        优先级：
        1. 若已建立全局慢化剂环，则直接对全局慢化剂环求平均
        2. 否则回退到各代表元件内部 moderator，并按 multiplier 加权
        """
        if self.has_global_moderator and self.mod_rings:
            numerator = 0.0
            denominator = 0.0
            for ring in self.mod_rings:
                volumes = np.asarray(ring.mesh.geom_data.volumes, dtype=float)
                numerator += float(np.sum(ring.T * volumes))
                denominator += float(np.sum(volumes))

            if denominator > 0.0:
                return numerator / denominator

        return self._compute_representative_structure_average_temperature('moderator')

    def collect_reactivity_feedback_temperatures(self) -> ReactivityFeedbackTemperatureSummary:
        """
        采集反应性反馈所需的全堆结构平均温度。

        本方法只负责“温度统计”，不负责代入反馈公式。
        后续若要替换反馈关联式，只需改公式层，不必改统计层。
        """
        summary = ReactivityFeedbackTemperatureSummary(
            fuel=self._compute_representative_structure_average_temperature('pellet'),
            emitter=self._compute_representative_structure_average_temperature('emitter'),
            collector=self._compute_representative_structure_average_temperature('collector'),
            moderator=self._compute_global_moderator_average_temperature(),
            reflector=self._compute_global_structure_average_temperature(self.reflector),
            barrel=self._compute_global_structure_average_temperature(self.barrel),
        )
        self.last_feedback_temperatures = summary
        return summary

    @staticmethod
    def _compute_fuel_reactivity_feedback(temperature: Optional[float]) -> float:
        """燃料反馈多项式，公式来自旧版 TASTIN Fortran。"""
        if temperature is None:
            return 0.0

        T = float(temperature)
        return (
            0.001360811
            - 6.47927757e-6 * T
            + 2.321231e-9 * (T ** 2)
            - 3.52e-13 * (T ** 3)
        )

    @staticmethod
    def _compute_electrode_reactivity_feedback(emitter_temperature: Optional[float],
                                               collector_temperature: Optional[float]) -> float:
        """
        电极反馈多项式，公式来自旧版 TASTIN Fortran。

        注意：
        1. 发射极与接收极是两项之和
        2. 接收极这里采用 Fortran 版本中的独立公式，不再沿用旧 Python 中错误复用发射极公式的写法
        """
        rho_emitter = 0.0
        rho_collector = 0.0

        if emitter_temperature is not None:
            T = float(emitter_temperature)
            rho_emitter = 1.0e-4 * (3.46455e-6 * (T ** 2) - 0.03232167 * T + 0.74202216)

        if collector_temperature is not None:
            T = float(collector_temperature)
            rho_collector = 1.0e-4 * (1.7094427e-6 * (T ** 2) - 0.017584414 * T - 0.7409044027)

        return rho_emitter + rho_collector

    @staticmethod
    def _compute_moderator_reactivity_feedback(temperature: Optional[float]) -> float:
        """慢化剂反馈多项式，公式来自旧版 TASTIN Fortran。"""
        if temperature is None:
            return 0.0

        T = float(temperature)
        return -1.7251979e-8 * (T ** 2) + 5.9009e-5 * T - 0.0148233537

    @staticmethod
    def _compute_reflector_reactivity_feedback(temperature: Optional[float]) -> float:
        """反射层反馈多项式，公式来自旧版 TASTIN Fortran。"""
        if temperature is None:
            return 0.0

        T = float(temperature)
        return -3.861086e-9 * (T ** 2) + 1.27026e-5 * T - 0.0035761

    def compute_reactivity_feedback(self) -> ReactivityFeedbackResult:
        """
        计算当前时刻的全堆反应性反馈。

        当前纳入的反馈项：
        1. 燃料反馈
        2. 电极反馈（发射极 + 接收极）
        3. 慢化剂反馈
        4. 侧反射层反馈

        说明：
        1. 筒体当前只统计温度，不参与反馈公式
        2. 后续若要新增更多结构反馈，只需要在这里扩展
        """
        temperatures = self.collect_reactivity_feedback_temperatures()

        fuel_feedback = self._compute_fuel_reactivity_feedback(temperatures.fuel)
        electrode_feedback = self._compute_electrode_reactivity_feedback(
            emitter_temperature=temperatures.emitter,
            collector_temperature=temperatures.collector
        )
        moderator_feedback = self._compute_moderator_reactivity_feedback(temperatures.moderator)
        reflector_feedback = self._compute_reflector_reactivity_feedback(temperatures.reflector)

        result = ReactivityFeedbackResult(
            temperatures=temperatures,
            fuel=fuel_feedback,
            electrode=electrode_feedback,
            moderator=moderator_feedback,
            reflector=reflector_feedback,
            total=fuel_feedback + electrode_feedback + moderator_feedback + reflector_feedback
        )
        self.last_feedback_result = result
        return result

    def get_reactivity_feedback(self) -> float:
        """
        返回当前全堆总反馈，便于后续直接接入点堆模块。
        """
        return self.compute_reactivity_feedback().total

    def calibrate_reactivity_feedback_reference(self,
                                                reference_result: Optional[ReactivityFeedbackResult] = None
                                                ) -> ReactivityFeedbackResult:
        """
        标定反应性反馈基准点。

        设计目的：
        1. 点堆方程应接收“相对初始稳态”的反馈增量，而不是绝对反馈值
        2. 这样在 ``reactivity_control = 0`` 的情况下，初始化工况可以保持临界平衡
        """
        if reference_result is None:
            reference_result = self.compute_reactivity_feedback()

        self.feedback_reference_result = ReactivityFeedbackResult(
            temperatures=ReactivityFeedbackTemperatureSummary(
                fuel=reference_result.temperatures.fuel,
                emitter=reference_result.temperatures.emitter,
                collector=reference_result.temperatures.collector,
                moderator=reference_result.temperatures.moderator,
                reflector=reference_result.temperatures.reflector,
                barrel=reference_result.temperatures.barrel,
            ),
            fuel=float(reference_result.fuel),
            electrode=float(reference_result.electrode),
            moderator=float(reference_result.moderator),
            reflector=float(reference_result.reflector),
            total=float(reference_result.total)
        )
        self.last_effective_reactivity_feedback = 0.0
        return self.feedback_reference_result

    def get_effective_reactivity_feedback(self) -> float:
        """
        返回当前相对反馈增量。

        定义：
        ``effective_feedback = current_feedback - reference_feedback``
        """
        current_feedback = self.compute_reactivity_feedback()
        self.last_effective_reactivity_feedback = (
            current_feedback.total - self.feedback_reference_result.total
        )
        return self.last_effective_reactivity_feedback

    def _build_thermo_calc(self):
        """构建 TEC 电路的虚拟分身映射。"""
        self.total_virtual_elements = sum(self.tec_multipliers.values())
        self.n_nodes = self._get_reference_tfe().mesh.n_axial

        if self.total_virtual_elements <= 0 or not self.enable_tec_coupled:
            self.thermo_calc = None
            return

        try:
            self.thermo_calc = ThermoCalcModel(
                n_elements=self.total_virtual_elements,
                n_nodes=self.n_nodes
            )
            self._configure_thermo_calc_geometry()
            logger.info(
                f"ReactorCore '{self.name}': Built ThermoCalc circuit with "
                f"{len(self.tfes)} physical TFE(s) representing "
                f"{self.total_virtual_elements} virtual element(s) in series. (TEC is ENABLED)"
            )
        except Exception as exc:
            logger.warning(
                f"ReactorCore '{self.name}': failed to build ThermoCalc circuit, "
                f"disabling TEC coupling. Detail: {exc}"
            )
            self.thermo_calc = None
            self.enable_tec_coupled = False

    def _configure_thermo_calc_geometry(self):
        """
        保持 ThermoCalc 轴向长度和侧面积与 TFE 热网格一致。

        TEC 求解器接受逐节点的轴向长度和逐节点的侧面积，
        因此这两个数组都从 TFE 热网格复制。这使得反射层节点
        和有源区节点各自使用对应的物理段面积，而不是采用
        单一的有源区标量近似。
        """
        if self.thermo_calc is None:
            return

        dl_e = np.zeros((self.total_virtual_elements, self.n_nodes), dtype=float)
        dl_c = np.zeros((self.total_virtual_elements, self.n_nodes), dtype=float)
        side_area_e = np.zeros((self.total_virtual_elements, self.n_nodes), dtype=float)
        side_area_c = np.zeros((self.total_virtual_elements, self.n_nodes), dtype=float)

        idx = 0
        for tfe_name, multiplier in self.tec_multipliers.items():
            multiplier = int(multiplier)
            if multiplier <= 0:
                continue
            tfe = self.tfes[tfe_name]
            if hasattr(tfe, 'common_y_faces'):
                node_lengths = np.diff(np.asarray(tfe.common_y_faces, dtype=float))
            else:
                node_lengths = np.diff(np.asarray(tfe.solids['emitter'].mesh.y_faces, dtype=float))
            if node_lengths.shape != (self.n_nodes,):
                raise ValueError(
                    f"TFE '{tfe_name}' axial node length shape {node_lengths.shape} "
                    f"does not match ThermoCalc n_nodes {(self.n_nodes,)}."
                )

            emitter_area = np.asarray(
                tfe.solids['emitter'].boundaries['right'].area,
                dtype=float
            )
            collector_area = np.asarray(
                tfe.solids['collector'].boundaries['left'].area,
                dtype=float
            )
            if emitter_area.shape != (self.n_nodes,) or collector_area.shape != (self.n_nodes,):
                raise ValueError(
                    f"TFE '{tfe_name}' electrode side area shape "
                    f"{emitter_area.shape}/{collector_area.shape} does not match "
                    f"ThermoCalc n_nodes {(self.n_nodes,)}."
                )

            side_area_e[idx: idx + multiplier, :] = emitter_area
            side_area_c[idx: idx + multiplier, :] = collector_area
            dl_e[idx: idx + multiplier, :] = node_lengths
            dl_c[idx: idx + multiplier, :] = node_lengths
            idx += multiplier

        self.thermo_calc._input_data.dlE = dl_e
        self.thermo_calc._input_data.dlC = dl_c
        self.thermo_calc._input_data.sideAreaE = side_area_e
        self.thermo_calc._input_data.sideAreaC = side_area_c

    def _get_reference_tfe(self) -> TFEUnit:
        """返回一个代表性 TFE，用于读取统一的轴向网格。"""
        return list(self.tfes.values())[0]

    def _get_reference_y_faces(self) -> np.ndarray:
        """
        返回与 TFE 完全一致的轴向界面坐标。

        这样全局结构与 TFE 在轴向上严格对齐，后续做温度同步和反馈统计时会更自然。
        """
        ref_tfe = self._get_reference_tfe()
        if hasattr(ref_tfe, 'common_y_faces'):
            return np.asarray(ref_tfe.common_y_faces, dtype=float)

        # 向后兼容：如果未来某个 TFE 未显式保存 common_y_faces，就退回到 pellet 网格。
        return np.asarray(ref_tfe.solids['pellet'].mesh.y_faces, dtype=float)

    def _register_solid(self, key: str, solid: HeatConduction2D):
        """统一注册 ReactorCore 自己创建的附加固体。"""
        self.extra_solids[key] = solid

    def _register_coupler(self, key: str, coupler: Any):
        """统一注册 ReactorCore 自己创建的附加耦合器。"""
        self.extra_couplers[key] = coupler

    def _get_solid_inner_outer_radius(self, solid: HeatConduction2D) -> Tuple[float, float]:
        """读取圆柱坐标二维固体的内外半径。"""
        mesh = solid.mesh
        return float(mesh.x_faces[0]), float(mesh.x_faces[-1])

    def _create_cylindrical_annulus_mesh(self,
                                         inner_radius: float,
                                         outer_radius: float,
                                         n_radial: int) -> Mesh2D:
        """
        为全局环形结构生成二维空心圆柱网格。

        轴向划分严格复用 TFE 的 ``y_faces``，保证整体网格在轴向上一一对应。
        """
        y_faces = self._get_reference_y_faces()
        return Mesh2D(
            x_dim=outer_radius - inner_radius,
            n_x=n_radial,
            y_dim=y_faces[-1] - y_faces[0],
            n_y=len(y_faces) - 1,
            y_faces=y_faces,
            geometry_type='cylindrical',
            inner_radius=inner_radius
        )

    def _resolve_annulus_config(self,
                                config: GlobalAnnulusStructureConfig,
                                previous_outer_radius: float,
                                gap_width: float) -> Tuple[Mesh2D, float, float]:
        """
        解析环形固体层配置。

        优先级：
        1. 若显式传入 mesh，则直接使用 mesh
        2. 否则按 ``inner_radius/outer_radius/thickness`` 自动生成
        """
        if config.mesh is not None:
            mesh = config.mesh
            return mesh, float(mesh.x_faces[0]), float(mesh.x_faces[-1])

        inner_radius = config.inner_radius
        if inner_radius is None:
            inner_radius = previous_outer_radius + max(gap_width, 0.0)

        if config.outer_radius is not None:
            outer_radius = config.outer_radius
        elif config.thickness is not None:
            outer_radius = inner_radius + config.thickness
        else:
            raise ValueError("全局环形结构配置必须提供 outer_radius、thickness 或 mesh。")

        if outer_radius <= inner_radius:
            raise ValueError("全局环形结构的外半径必须大于内半径。")

        mesh = self._create_cylindrical_annulus_mesh(
            inner_radius=inner_radius,
            outer_radius=outer_radius,
            n_radial=config.n_radial
        )
        return mesh, float(inner_radius), float(outer_radius)

    # =========================================================================
    # 全局慢化剂、筒体、反射层建模
    # =========================================================================

    def _build_global_moderator_rings(self,
                                      mod_meshes: Optional[List[Mesh2D]],
                                      mod_material: Optional[Any]):
        """建立全局慢化剂环，并与TFE虚拟慢化剂建立热耦合。"""
        if mod_meshes is None or mod_material is None or self.ring_mapping is None:
            return

        self.has_global_moderator = True

        for i, mesh in enumerate(mod_meshes):
            ring = HeatConduction2D(
                name=f"{self.name}_ModRing_{i}",
                mesh=mesh,
                material=mod_material,
                initial_temp=743.0
            )
            self.mod_rings.append(ring)
            self._register_solid(f'mod_ring_{i}', ring)

        logger.info(f"ReactorCore '{self.name}': Built {len(self.mod_rings)} global moderator rings.")

        # 慢化剂环之间做径向固固导热耦合。
        for i in range(len(self.mod_rings) - 1):
            self._register_coupler(
                key=f'mod_couple_{i}_{i + 1}',
                coupler=SolidSolidCouple2D(
                    obj1=self.mod_rings[i],
                    obj2=self.mod_rings[i + 1],
                    direction='right'
                )
            )

    def _build_global_outer_structures(self,
                                       barrel_config: Optional[GlobalAnnulusStructureConfig],
                                       reflector_config: Optional[GlobalAnnulusStructureConfig],
                                       moderator_barrel_gap_config: Optional[GlobalGapStructureConfig],
                                       barrel_reflector_gap_config: Optional[GlobalGapStructureConfig]):
        """
        建立堆芯外部全局结构。

        当前支持的层次顺序为：
        慢化剂最外层 -> 间隙 -> 筒体 -> 间隙 -> 反射层

        说明：
        1. 如果只建到筒体，则最外表面辐射边界会挂到筒体外表面
        2. 如果继续建反射层，则辐射边界会自动挂到反射层外表面
        """
        if barrel_config is None and reflector_config is None:
            return

        if not self.has_global_moderator or not self.mod_rings:
            raise ValueError("建立全局筒体/反射层前，必须先提供全局慢化剂环。")

        previous_solid = self.mod_rings[-1]
        _, previous_outer_radius = self._get_solid_inner_outer_radius(previous_solid)

        structure_plan = [
            ('barrel', barrel_config, moderator_barrel_gap_config),
            ('reflector', reflector_config, barrel_reflector_gap_config),
        ]

        for layer_name, layer_config, gap_config in structure_plan:
            if layer_config is None:
                continue

            gap_width = 0.0 if gap_config is None else gap_config.width
            layer_mesh, layer_inner_radius, layer_outer_radius = self._resolve_annulus_config(
                config=layer_config,
                previous_outer_radius=previous_outer_radius,
                gap_width=gap_width
            )

            # 若用户显式给出内半径，则额外做一次几何一致性检查，避免几何关系悄悄出错。
            expected_min_inner = previous_outer_radius + max(gap_width, 0.0)
            if layer_inner_radius + 1e-12 < previous_outer_radius:
                raise ValueError(f"{layer_name} 的内半径不能小于前一层的外半径。")
            if gap_width > 0.0 and layer_inner_radius + 1e-12 < expected_min_inner:
                raise ValueError(f"{layer_name} 的内半径与前一层间隙宽度不一致。")

            layer_solid = HeatConduction2D(
                name=f"{self.name}_{layer_name.capitalize()}",
                mesh=layer_mesh,
                material=layer_config.material,
                initial_temp=layer_config.initial_temp
            )
            self._register_solid(layer_name, layer_solid)

            if layer_name == 'barrel':
                self.barrel = layer_solid
                self.barrel_outer_surface_emissivity = layer_config.outer_surface_emissivity
            elif layer_name == 'reflector':
                self.reflector = layer_solid
                self.reflector_outer_surface_emissivity = layer_config.outer_surface_emissivity

            actual_gap_width = layer_inner_radius - previous_outer_radius

            self._connect_outer_layers(
                left_name='mod_outer' if previous_solid is self.mod_rings[-1] else previous_solid.name,
                left_solid=previous_solid,
                right_name=layer_name,
                right_solid=layer_solid,
                gap_config=gap_config,
                gap_width=actual_gap_width
            )

            # 对用户显式声明的间隙宽度做严格双向校验，避免“配置宽度”和“实际几何宽度”悄悄不一致。
            if gap_config is None:
                if abs(actual_gap_width) > 1e-12:
                    raise ValueError(
                        f"{layer_name} 与前一层之间实际存在 {actual_gap_width:.6e} m 间隙，"
                        "但未提供对应的间隙配置。"
                    )
            else:
                if abs(actual_gap_width - gap_config.width) > 1e-12:
                    raise ValueError(
                        f"{layer_name} 的实际间隙宽度为 {actual_gap_width:.6e} m，"
                        f"与配置值 {gap_config.width:.6e} m 不一致。"
                    )

            previous_solid = layer_solid
            previous_outer_radius = layer_outer_radius

        self.global_outer_solid = previous_solid

    def _connect_outer_layers(self,
                              left_name: str,
                              left_solid: HeatConduction2D,
                              right_name: str,
                              right_solid: HeatConduction2D,
                              gap_config: Optional[GlobalGapStructureConfig],
                              gap_width: float):
        """
        连接相邻两层全局结构。

        三种情况：
        1. 无间隙：直接使用 ``SolidSolidCouple2D``
        2. ``simplified`` 间隙：使用 ``GapCouple2D``
        3. ``meshed`` 间隙：间隙本身建成独立固体，同时保留纯辐射并联通道
        """
        if gap_width <= 1e-12:
            self._register_coupler(
                key=f'{left_name}_to_{right_name}',
                coupler=SolidSolidCouple2D(
                    obj1=left_solid,
                    obj2=right_solid,
                    direction='right'
                )
            )
            return

        if gap_config is None:
            raise ValueError(
                f"{left_name} 与 {right_name} 之间存在 {gap_width:.6e} m 的径向间隙，"
                "但未提供对应的间隙配置。"
            )

        if gap_config.mode == 'simplified':
            self._register_coupler(
                key=f'{left_name}_to_{right_name}_gap',
                coupler=GapCouple2D(
                    obj1=left_solid,
                    obj2=right_solid,
                    direction='right',
                    gap_width=gap_width,
                    gas_conductivity=gap_config.h_eq * gap_width,
                    emissivity1=gap_config.emissivity_inner,
                    emissivity2=gap_config.emissivity_outer
                )
            )
            return

        if gap_config.mode != 'meshed':
            raise ValueError(f"未知的全局间隙模式: {gap_config.mode}")

        if gap_config.material is None:
            raise ValueError("网格化全局间隙模式必须提供 material。")

        if gap_config.mesh is not None:
            gap_mesh = gap_config.mesh
        else:
            _, left_outer = self._get_solid_inner_outer_radius(left_solid)
            gap_mesh = self._create_cylindrical_annulus_mesh(
                inner_radius=left_outer,
                outer_radius=left_outer + gap_width,
                n_radial=gap_config.n_radial
            )

        gap_name = f'{left_name}_to_{right_name}_gap_solid'
        gap_solid = HeatConduction2D(
            name=f"{self.name}_{gap_name}",
            mesh=gap_mesh,
            material=gap_config.material,
            initial_temp=gap_config.initial_temp
        )
        self.global_gap_solids[gap_name] = gap_solid
        self._register_solid(gap_name, gap_solid)

        # 网格化间隙承担气体导热。
        self._register_coupler(
            key=f'{left_name}_to_{gap_name}',
            coupler=SolidSolidCouple2D(
                obj1=left_solid,
                obj2=gap_solid,
                direction='right'
            )
        )
        self._register_coupler(
            key=f'{gap_name}_to_{right_name}',
            coupler=SolidSolidCouple2D(
                obj1=gap_solid,
                obj2=right_solid,
                direction='right'
            )
        )

        # 与 TFE 中 meshed gap 的处理保持一致：再额外挂一条纯辐射并联通道。
        self._register_coupler(
            key=f'{left_name}_to_{right_name}_radiation',
            coupler=GapCouple2D(
                obj1=left_solid,
                obj2=right_solid,
                direction='right',
                gap_width=gap_width,
                gas_conductivity=0.0,
                emissivity1=gap_config.emissivity_inner,
                emissivity2=gap_config.emissivity_outer
            )
        )

    def _attach_outer_radiation_boundary(self, T_space: float):
        """
        给最外层结构挂接对外辐射边界。

        优先级：
        1. 反射层
        2. 筒体
        3. 最外层慢化剂环
        """
        target_solid = self.global_outer_solid
        if target_solid is None and self.mod_rings:
            target_solid = self.mod_rings[-1]

        if target_solid is None:
            return

        if target_solid is self.reflector:
            emissivity = self.reflector_outer_surface_emissivity
            if emissivity is None:
                emissivity = 0.6
        elif target_solid is self.barrel:
            emissivity = self.barrel_outer_surface_emissivity
            if emissivity is None:
                emissivity = 0.05
        else:
            emissivity = 0.05

        outermost_boundary = target_solid.boundaries['right']
        outermost_boundary.add_dynamic_radiation_condition(
            emissivity=emissivity,
            bare_area_array=outermost_boundary.area,
            T_env=T_space
        )
        logger.info(
            f"ReactorCore '{self.name}': Attached space radiation BC "
            f"(T_env={T_space} K, emissivity={emissivity}) to '{target_solid.name}'."
        )

    # =========================================================================
    # 对外公共接口
    # =========================================================================

    def setup_tec_circuit(self, mode_str: str, target_value: float, I_guess: float = 150.0):
        """封装 TEC 电路设置接口。"""
        if self.thermo_calc is not None:
            self.thermo_calc.setup_circuit_mode(mode_str, target_value, I_guess)

    def attach_point_reactor(self, point_reactor: Optional[PointReactor] = None):
        """
        将点堆动力学模块挂接到 ReactorCore。

        若未显式传入实例，则自动创建默认 ``PointReactor``。
        """
        self.point_reactor = point_reactor if point_reactor is not None else PointReactor()
        return self.point_reactor

    def initialize_point_reactor(self, total_power_initial: float, point_reactor: Optional[PointReactor] = None):
        """
        初始化 ReactorCore 内部的点堆动力学模块，并将初始功率下发到燃料内热源。
        """
        reactor = self.attach_point_reactor(point_reactor=point_reactor)
        reactor.initialize_steady_state(total_power_initial=total_power_initial)
        self.calibrate_reactivity_feedback_reference()

        # 初始化后立即同步一次功率，确保 t=0 时燃料内热源与点堆状态一致。
        self.update_neutronic_power(
            p_total=reactor.total_power,
            p_fiss=reactor.fission_power,
            p_decay=reactor.decay_power,
            alpha=1.0
        )
        return reactor

    @property
    def has_point_reactor(self) -> bool:
        """返回当前 ReactorCore 是否已经挂接点堆模块。"""
        return self.point_reactor is not None

    def update_neutronic_power(self,
                               p_total: float,
                               p_fiss: float = 0.0,
                               p_decay: float = 0.0,
                               alpha: float = 1.0):
        """
        将全堆总功率按 ReactorCore 级别的功率因子分配到各个代表 TFE。

        说明：
        1. 代表件收到的是“单根代表件功率”
        2. 全堆总功率守恒验算应按 ``代表件功率 × multiplier`` 回推
        """
        self.last_total_core_power = float(p_total)

        for tfe_name, tfe in self.tfes.items():
            power_factor = self.tfe_power_factors[tfe_name]
            representative_power = self.last_total_core_power * power_factor
            representative_fiss_power = float(p_fiss) * power_factor
            representative_decay_power = float(p_decay) * power_factor

            tfe.update_neutronic_power(
                p_total=representative_power,
                p_fiss=representative_fiss_power,
                p_decay=representative_decay_power,
                alpha=alpha
            )

    def get_power_distribution_summary(self) -> Dict[str, Any]:
        """
        返回当前功率分配摘要，便于调试和守恒验算。
        """
        representative_powers = {
            name: self.last_total_core_power * factor
            for name, factor in self.tfe_power_factors.items()
        }
        reconstructed_total_real_power = sum(
            representative_powers[name] * self.tfe_multipliers[name]
            for name in representative_powers
        )
        return {
            'tfe_power_factors': dict(self.tfe_power_factors),
            'tfe_group_power_shares': dict(self.tfe_group_power_shares),
            'power_factor_weighted_sum': self.power_factor_weighted_sum,
            'last_total_core_power': self.last_total_core_power,
            'representative_powers': representative_powers,
            'reconstructed_total_real_power': reconstructed_total_real_power,
        }

    def advance_neutronics(self,
                           dt: float,
                           reactivity_control: float = 0.0,
                           iteration_index: int = 0) -> Optional[ReactivityFeedbackResult]:
        """
        在当前热工状态下推进点堆动力学，并把新的功率映射到各个代表元件。

        调用时机：
        1. 由 SystemManager 在每次内迭代中调用
        2. 必须发生在热工水力求解之前，以便新功率进入本次求解
        """
        if self.point_reactor is None:
            return None

        feedback_result = self.compute_reactivity_feedback()
        self.last_reactivity_control = float(reactivity_control)
        effective_feedback = feedback_result.total - self.feedback_reference_result.total
        self.last_effective_reactivity_feedback = float(effective_feedback)

        self.point_reactor.step(
            dt=dt,
            reactivity_control=self.last_reactivity_control,
            reactivity_feedback=self.last_effective_reactivity_feedback
        )
        self.update_neutronic_power(
            p_total=self.point_reactor.total_power,
            p_fiss=self.point_reactor.fission_power,
            p_decay=self.point_reactor.decay_power,
            alpha=1.0
        )
        return feedback_result

    def commit_neutronics(self):
        """
        固化本时间步的点堆状态。
        """
        if self.point_reactor is not None:
            self.point_reactor.commit()
            return True
        return False

    @staticmethod
    def _capture_solid_source_state(solid):
        state = {}
        if hasattr(solid, 'Q_source'):
            state['Q_source_ref'] = solid.Q_source
            state['Q_source_value'] = np.asarray(solid.Q_source, dtype=float).copy()
        if hasattr(solid, 'use_external_source_buffer'):
            state['use_external_source_buffer'] = bool(solid.use_external_source_buffer)
        if hasattr(solid, 'source_callback'):
            state['source_callback'] = solid.source_callback
        return state

    @staticmethod
    def _restore_solid_source_state(solid, state):
        if not state:
            return
        if 'Q_source_ref' in state:
            solid.Q_source = state['Q_source_ref']
            solid.Q_source[:] = state['Q_source_value']
        if 'use_external_source_buffer' in state:
            solid.use_external_source_buffer = state['use_external_source_buffer']
        if 'source_callback' in state:
            solid.source_callback = state['source_callback']

    def save_step_state(self):
        tfe_states = {}
        for name, tfe in self.tfes.items():
            save_step_state = getattr(tfe, 'save_step_state', None)
            if callable(save_step_state):
                tfe_states[name] = save_step_state()

        point_reactor_state = None
        if self.point_reactor is not None:
            save_step_state = getattr(self.point_reactor, 'save_step_state', None)
            if callable(save_step_state):
                point_reactor_state = save_step_state()

        solid_source_states = [
            (solid, self._capture_solid_source_state(solid))
            for solid in self.get_solids()
        ]

        return {
            '_last_thermo_update_time': self._last_thermo_update_time,
            'last_total_core_power': self.last_total_core_power,
            'last_effective_reactivity_feedback': self.last_effective_reactivity_feedback,
            'last_reactivity_control': self.last_reactivity_control,
            'tfe_states': tfe_states,
            'point_reactor_state': point_reactor_state,
            'solid_source_states': solid_source_states,
        }

    def load_step_state(self, state):
        if state is None:
            return

        self._last_thermo_update_time = state['_last_thermo_update_time']
        self.last_total_core_power = state['last_total_core_power']
        self.last_effective_reactivity_feedback = state['last_effective_reactivity_feedback']
        self.last_reactivity_control = state['last_reactivity_control']

        for name, tfe_state in state['tfe_states'].items():
            tfe = self.tfes.get(name)
            if tfe is None:
                continue
            load_step_state = getattr(tfe, 'load_step_state', None)
            if callable(load_step_state):
                load_step_state(tfe_state)

        if self.point_reactor is not None and state['point_reactor_state'] is not None:
            load_step_state = getattr(self.point_reactor, 'load_step_state', None)
            if callable(load_step_state):
                load_step_state(state['point_reactor_state'])

        for solid, source_state in state['solid_source_states']:
            self._restore_solid_source_state(solid, source_state)

    def get_solids(self) -> list:
        all_solids = list(self.extra_solids.values())
        for tfe in self.tfes.values():
            all_solids.extend(tfe.get_solids())
        return all_solids

    def get_couplers(self) -> list:
        all_couplers = list(self.extra_couplers.values())
        for tfe in self.tfes.values():
            all_couplers.extend(tfe.get_couplers())
        return all_couplers

    # =========================================================================
    # 生命周期钩子
    # =========================================================================

    def pre_step(self, dt: float, current_time: float):
        """
        求解前钩子。

        当前主要处理：
        1. TEC 电路计算并将热流/焦耳热下发到各个 TFE
        2. 从虚拟慢化剂抽取热量并注入全局慢化剂环
        """
        if self.enable_tec_coupled and self.thermo_calc is not None:
            time_since_last_update = current_time - self._last_thermo_update_time
            if time_since_last_update >= self.thermo_update_interval:
                self.thermo_calc.calculate(verbose=False)
                self._last_thermo_update_time = current_time

            idx = 0
            for tfe_name, thermal_mult in self.tfe_multipliers.items():
                tec_mult = int(self.tec_multipliers.get(tfe_name, 0))
                if tec_mult <= 0:
                    continue
                tfe = self.tfes[tfe_name]
                tec_idx = idx
                res = self.thermo_calc.get_tec_results(tec_idx)
                idx += tec_mult

                if res is None:
                    continue

                UE_abs = res.get('UE', np.zeros(self.n_nodes))
                UC_abs = res.get('UC', np.zeros(self.n_nodes))
                rho_e = res.get('rhoE', np.ones(self.n_nodes) * 1e-6)
                rho_c = res.get('rhoC', np.ones(self.n_nodes) * 1e-6)
                joule_power_e = np.asarray(res['joulePowerE'], dtype=float)
                joule_power_c = np.asarray(res['joulePowerC'], dtype=float)

                if hasattr(tfe, 'common_y_faces'):
                    E_e = electric_field_from_node_potential(
                        UE_abs,
                        y_faces=np.asarray(tfe.common_y_faces, dtype=float),
                    )
                    E_c = electric_field_from_node_potential(
                        UC_abs,
                        y_faces=np.asarray(tfe.common_y_faces, dtype=float),
                    )
                else:
                    input_data = self.thermo_calc._input_data
                    E_e = electric_field_from_node_potential(
                        UE_abs,
                        node_lengths=np.asarray(input_data.dlE[tec_idx], dtype=float),
                    )
                    E_c = electric_field_from_node_potential(
                        UC_abs,
                        node_lengths=np.asarray(input_data.dlC[tec_idx], dtype=float),
                    )

                J_density = res.get('J', np.zeros(self.n_nodes)) * 1e4
                phiE = res.get('phiE', np.zeros(self.n_nodes))
                TE = res.get('TE', np.zeros(self.n_nodes))

                q_e_flux = -1.0 * J_density * (phiE + 2.0 * 8.617e-5 * TE)
                q_c_flux = 1.0 * J_density * (
                    phiE + 2.0 * 8.617e-5 * TE - (UE_abs - UC_abs)
                )

                electric_alpha = self.alpha_tec * float(tec_mult) / float(thermal_mult)

                tfe.update_electric_field_diagnostics(
                    E_emit=E_e,
                    rho_emit=rho_e,
                    E_coll=E_c,
                    rho_coll=rho_c
                )
                tfe.update_joule_power_sources(
                    Q_emitter_axial=joule_power_e,
                    Q_collector_axial=joule_power_c,
                    alpha=electric_alpha
                )
                tfe.update_plasma_flux(
                    q_e_flux=q_e_flux,
                    q_c_flux=q_c_flux,
                    alpha=electric_alpha
                )

        for tfe in self.tfes.values():
            if hasattr(tfe, 'pre_step'):
                tfe.pre_step(dt, current_time)

        if self.has_global_moderator:
            for ring in self.mod_rings:
                ring.use_external_source_buffer = True
                ring.Q_source[:] = 0.0

            for ring_idx, ring in enumerate(self.mod_rings):
                member_names = self.get_ring_member_names(ring_idx)
                if not member_names:
                    continue

                # 同圈可能存在多个代表元件，必须把它们全部计入，并按 multiplier 还原成真实总热流。
                q_watts_total = np.zeros(ring.shape_nodes[1], dtype=float)
                for tfe_name in member_names:
                    tfe = self.tfes[tfe_name]
                    mult = self.tfe_multipliers[tfe_name]
                    virtual_mod = tfe.solids['moderator']

                    # TFEUnit.pre_step() may update the moderator outer-boundary
                    # temperature. Refresh this cache before mapping its flux to
                    # the global moderator ring, including immediately after a
                    # restart.
                    virtual_mod._update_properties()
                    virtual_mod._compute_internal_resistance()
                    virtual_mod._update_boundaries_state(current_time=current_time)
                    virtual_mod._compute_fluxes(current_time)

                    q_flux_out = -virtual_mod.boundaries['right'].current_flux
                    q_watts_total += q_flux_out * mult

                nx, ny = ring.shape_nodes
                vols_2d = ring.mesh.geom_data.volumes.reshape(nx, ny)
                vol_axial = np.sum(vols_2d, axis=0)

                q_vol_1d = q_watts_total / np.maximum(vol_axial, 1e-12)
                q_vol_2d = np.ones((nx, 1)) * q_vol_1d[np.newaxis, :]
                q_watts_2d = q_vol_2d * vols_2d
                ring.Q_source += q_watts_2d.flatten()

    def post_step(self, dt: float, current_time: float):
        """
        求解后钩子。

        当前主要处理：
        1. 将收敛后的表面温度同步到 TEC 电路
        2. 将全局慢化剂环的轴向平均温度回填到各个 TFE 外侧边界
        """
        for tfe in self.tfes.values():
            if hasattr(tfe, 'post_step'):
                tfe.post_step(dt, current_time)

        if self.enable_tec_coupled and self.thermo_calc is not None:
            T_em_matrix = np.zeros((self.total_virtual_elements, self.n_nodes))
            T_co_matrix = np.zeros((self.total_virtual_elements, self.n_nodes))

            idx = 0
            for tfe_name, mult in self.tec_multipliers.items():
                mult = int(mult)
                if mult <= 0:
                    continue
                tfe = self.tfes[tfe_name]
                T_e = tfe.solids['emitter'].boundaries['right'].T_surface
                T_c = tfe.solids['collector'].boundaries['left'].T_surface
                T_em_matrix[idx: idx + mult, :] = T_e
                T_co_matrix[idx: idx + mult, :] = T_c
                idx += mult

            self.thermo_calc.set_temperatures(T_em_matrix, T_co_matrix)

        if self.has_global_moderator:
            for ring_idx, ring in enumerate(self.mod_rings):
                member_names = self.get_ring_member_names(ring_idx)
                if not member_names:
                    continue

                nx, ny = ring.shape_nodes
                T_2d = ring.T.reshape(nx, ny)
                T_avg_axial = np.mean(T_2d, axis=0)

                # 同圈多个代表元件共享同一个慢化剂温度场。
                for tfe_name in member_names:
                    self.tfes[tfe_name].boundary_data.moderator_temperature[:] = T_avg_axial

    # =========================================================================
    # 断点续算
    # =========================================================================

    def get_state_dict(self, prefix: str) -> dict:
        """保存 ReactorCore 自身的宏观状态。"""
        # 构建状态字典，用于断点续算
        # 包含TEC耦合标志、热更新时间、功率和反应性反馈等宏观状态
        state = {
            # TEC耦合使能标志
            f"{prefix}/enable_tec_coupled": np.array([self.enable_tec_coupled], dtype=bool),
            # 上次热更新时间
            f"{prefix}/_last_thermo_update_time": np.array([self._last_thermo_update_time]),
            # 上次总堆芯功率
            f"{prefix}/last_total_core_power": np.array([self.last_total_core_power], dtype=float),
            # 上次有效反应性反馈
            f"{prefix}/last_effective_reactivity_feedback": np.array(
                [self.last_effective_reactivity_feedback], dtype=float
            ),
            # 上次反应性控制值
            f"{prefix}/last_reactivity_control": np.array([self.last_reactivity_control], dtype=float),
            # 反应性反馈参考值 - 燃料温度系数贡献
            f"{prefix}/feedback_reference/fuel": np.array([self.feedback_reference_result.fuel], dtype=float),
            # 反应性反馈参考值 - 电极温度系数贡献
            f"{prefix}/feedback_reference/electrode": np.array([self.feedback_reference_result.electrode], dtype=float),
            # 反应性反馈参考值 - 慢化剂温度系数贡献
            f"{prefix}/feedback_reference/moderator": np.array([self.feedback_reference_result.moderator], dtype=float),
            # 反应性反馈参考值 - 反射层温度系数贡献
            f"{prefix}/feedback_reference/reflector": np.array([self.feedback_reference_result.reflector], dtype=float),
            # 反应性反馈参考值 - 总反应性
            f"{prefix}/feedback_reference/total": np.array([self.feedback_reference_result.total], dtype=float),
            # 反应性反馈参考温度 - 燃料平均温度（None时用nan表示）
            f"{prefix}/feedback_reference_temperature/fuel": np.array(
                [np.nan if self.feedback_reference_result.temperatures.fuel is None
                 else self.feedback_reference_result.temperatures.fuel],
                dtype=float
            ),
            # 反应性反馈参考温度 - 发射极平均温度（None时用nan表示）
            f"{prefix}/feedback_reference_temperature/emitter": np.array(
                [np.nan if self.feedback_reference_result.temperatures.emitter is None
                 else self.feedback_reference_result.temperatures.emitter],
                dtype=float
            ),
            # 反应性反馈参考温度 - 收集极平均温度（None时用nan表示）
            f"{prefix}/feedback_reference_temperature/collector": np.array(
                [np.nan if self.feedback_reference_result.temperatures.collector is None
                 else self.feedback_reference_result.temperatures.collector],
                dtype=float
            ),
            # 反应性反馈参考温度 - 慢化剂平均温度（None时用nan表示）
            f"{prefix}/feedback_reference_temperature/moderator": np.array(
                [np.nan if self.feedback_reference_result.temperatures.moderator is None
                 else self.feedback_reference_result.temperatures.moderator],
                dtype=float
            ),
            # 反应性反馈参考温度 - 反射层平均温度（None时用nan表示）
            f"{prefix}/feedback_reference_temperature/reflector": np.array(
                [np.nan if self.feedback_reference_result.temperatures.reflector is None
                 else self.feedback_reference_result.temperatures.reflector],
                dtype=float
            ),
            # 反应性反馈参考温度 - 筒体平均温度（None时用nan表示）
            f"{prefix}/feedback_reference_temperature/barrel": np.array(
                [np.nan if self.feedback_reference_result.temperatures.barrel is None
                 else self.feedback_reference_result.temperatures.barrel],
                dtype=float
            ),
        }

        for tfe_name, tfe in self.tfes.items():
            if hasattr(tfe, 'get_state_dict'):
                state.update(tfe.get_state_dict(prefix=f"{prefix}/TFEs/{tfe_name}"))

        if self.point_reactor is not None and hasattr(self.point_reactor, 'get_state_dict'):
            state.update(self.point_reactor.get_state_dict(prefix=f"{prefix}/PointReactor"))

        return state

    def load_state_dict(self, data: dict, prefix: str):
        """恢复 ReactorCore 自身的宏观状态。"""
        flag_key = f"{prefix}/enable_tec_coupled"
        if flag_key in data:
            self.enable_tec_coupled = bool(data[flag_key][0])

        key = f"{prefix}/_last_thermo_update_time"
        if key in data:
            self._last_thermo_update_time = float(data[key][0])

        key = f"{prefix}/last_total_core_power"
        if key in data:
            self.last_total_core_power = float(data[key][0])

        key = f"{prefix}/last_effective_reactivity_feedback"
        if key in data:
            self.last_effective_reactivity_feedback = float(data[key][0])

        key = f"{prefix}/last_reactivity_control"
        if key in data:
            self.last_reactivity_control = float(data[key][0])

        reference_prefix = f"{prefix}/feedback_reference"
        reference_temperature_prefix = f"{prefix}/feedback_reference_temperature"

        def _load_optional_temperature(field_name: str) -> Optional[float]:
            temp_key = f"{reference_temperature_prefix}/{field_name}"
            if temp_key not in data:
                return None

            value = float(data[temp_key][0])
            return None if np.isnan(value) else value

        total_key = f"{reference_prefix}/total"
        if total_key in data:
            self.feedback_reference_result = ReactivityFeedbackResult(
                temperatures=ReactivityFeedbackTemperatureSummary(
                    fuel=_load_optional_temperature('fuel'),
                    emitter=_load_optional_temperature('emitter'),
                    collector=_load_optional_temperature('collector'),
                    moderator=_load_optional_temperature('moderator'),
                    reflector=_load_optional_temperature('reflector'),
                    barrel=_load_optional_temperature('barrel'),
                ),
                fuel=float(data.get(f"{reference_prefix}/fuel", np.array([0.0]))[0]),
                electrode=float(data.get(f"{reference_prefix}/electrode", np.array([0.0]))[0]),
                moderator=float(data.get(f"{reference_prefix}/moderator", np.array([0.0]))[0]),
                reflector=float(data.get(f"{reference_prefix}/reflector", np.array([0.0]))[0]),
                total=float(data[total_key][0])
            )

        for tfe_name, tfe in self.tfes.items():
            if hasattr(tfe, 'load_state_dict'):
                tfe.load_state_dict(data, prefix=f"{prefix}/TFEs/{tfe_name}")

        point_reactor_prefix = f"{prefix}/PointReactor"
        point_reactor_key = f"{point_reactor_prefix}/y_committed"
        if point_reactor_key in data:
            if self.point_reactor is None:
                self.attach_point_reactor()
            self.point_reactor.load_state_dict(data, prefix=point_reactor_prefix)
