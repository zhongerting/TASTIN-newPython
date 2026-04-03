import copy
from typing import Any, Callable, Dict, List, Optional

import numpy as np

from Components.BaseComponent import BaseComponent
from Components.HPwithFin import HPwithFin
from Components.ExternalHeatSources import (
    AlbedoHeatSource,
    CompositeHeatSource,
    EarthIRHeatSource,
    ExternalHeatFluxBC,
    OrbitalHeatSource,
    OrbitalTableHeatSource,
)
from Solvers.Couplers import FluidSolidCouple


class SingleVolumeProxy:
    """
    单节点流体代理。

    `FluidSolidCouple` 期望接收的是一个“通道式”的流体对象，
    但在 `RingHP` 中，每个控制体位置只需要和一根代表性热管蒸发段换热。

    这个代理类的作用就是：
    1. 把一个标量 `FluidVolume` 伪装成单节点向量化流体对象。
    2. 保持换热面积仍然按热管蒸发段长度 `L_eva` 计算。
    3. 保持流体等效热容仍与原始控制体体积一致。
    """

    def __init__(self, vol, channel, vol_idx: int, L_eva: float, N_hp: float):
        self.vol = vol
        self.channel = channel
        self.idx = vol_idx
        self.n_nodes = 1
        self.N_hp = N_hp

        # 用蒸发段长度定义换热特征长度 A = perimeter * L。
        self.node_length = L_eva

        # 通过 area * node_length = volume 反推单节点等效流通面积，
        # 这样既满足耦合器接口，又不会改变流体本体热容。
        self.area = vol.vol / L_eva

        self.d_h = channel.d_h
        self.material = channel.material

    @property
    def temperature_vector(self):
        return np.array([self.vol.T])

    @property
    def pressure_vector(self):
        return np.array([self.vol.P])

    @property
    def density_vector(self):
        return np.array([self.vol.rho])

    @property
    def velocity_vector(self):
        return np.array([self.channel.velocity_vector[self.idx]])

    def add_coupling_source_distribution(self, explicit_arr, implicit_arr):
        # Coupler 返回的是“单根代表热管”的源项，
        # 这里按该控制体实际代表的热管数量 N_hp 折回总量。
        total_exp = explicit_arr[0] * self.N_hp
        total_imp = implicit_arr[0] * self.N_hp
        self.vol.add_coupling_source(total_exp, total_imp)


class RingHP(BaseComponent):
    """
    集流环 + 代表性热管阵列宏观组件。

    组件内部包含两部分：
    1. 集流环流体与集流环壁面的流固换热。
    2. 每个流体控制体位置对应的一根代表性热管及其蒸发段耦合。

    从 2026-04-03 起，本组件还支持将轨道外热流配置直接传给内部构造的 `HPwithFin`，
    使 `RingHP` 与单独使用 `HPwithFin` 时的外热流能力保持一致。
    """

    def __init__(self,
                 name: str,
                 fluid_channel,
                 solid_header,
                 hp_multipliers: List[float],
                 header_flow_area: float,
                 header_dh: float,
                 header_heated_perimeter: float,
                 hp_r_out: float, hp_r_in: float, hp_r_vapor: float,
                 hp_L_eva: float, hp_L_con: float,
                 hp_n_eva: int, hp_n_con: int,
                 hp_n_wick: int, hp_n_wall: int,
                 porosity_hp: float, HP_initial_temp: float,
                 fin_thickness: float, fin_height: float, n_fin_height: int,
                 fin_wrap_ratio: float, emissivity: float, up_view_factor: float, down_view_factor: float,
                 T_space: float,
                 hp_wall_mat, hp_fluid_mat, hp_wick_mat,
                 header_correlation_func: Callable,
                 hp_crossflow_base_func: Callable,
                 C_D: float = 1.0,
                 external_heat_config: Optional[Dict[str, Any]] = None):
        super().__init__(name)

        n_nodes = fluid_channel.n_nodes
        if len(hp_multipliers) != n_nodes:
            raise ValueError(
                f"[{name}] hp_multipliers 长度 ({len(hp_multipliers)}) 必须与流体节点数 ({n_nodes}) 一致。"
            )

        self.fluid_channel = fluid_channel
        self.solid_header = solid_header
        self.external_heat_config = copy.deepcopy(external_heat_config) if external_heat_config is not None else None

        # 额外保存几何参数，便于后处理或后续扩展。
        self.header_flow_area = float(header_flow_area)
        self.header_dh = float(header_dh)
        self.hp_r_out = float(hp_r_out)
        self.hp_L_eva = float(hp_L_eva)
        self.hp_n_con = int(hp_n_con)
        self.n_header_nodes = int(n_nodes)

        # 防呆：确保集流环壁温度场、边界和热容信息都已初始化。
        self.solid_header.initialize_state()

        # ===== A. 集流环内壁流固换热 =====
        cap_header = self.solid_header.get_boundary_node_capacitance('left')
        self.coupler_header = FluidSolidCouple(
            name=f"{name}_coupler_header",
            fluid=self.fluid_channel,
            solid_boundary_region=self.solid_header.boundaries['left'],
            heated_perimeter=header_heated_perimeter,
            correlation_func=header_correlation_func,
            solid_node_capacitance=cap_header
        )

        self.hp_units: List[HPwithFin] = []
        self.coupler_hps: List[FluidSolidCouple] = []
        self._outlet_k_loss = 0.0
        self._hp_k_loss_distribution = np.zeros(n_nodes, dtype=float)
        self._hp_presence_mask = np.zeros(n_nodes, dtype=bool)
        self._hp_external_heat_enabled = np.zeros(n_nodes, dtype=bool)

        # ===== B. 为每个流体控制体构造代表热管 =====
        for i in range(n_nodes):
            N_hp = hp_multipliers[i]
            K_loss_val = 0.0

            if N_hp > 0:
                self._hp_presence_mask[i] = True

                hp_name = f"{name}_HP_node{i}"
                hp = HPwithFin(
                    name=hp_name,
                    r_out_wall=hp_r_out,
                    r_in_wall=hp_r_in,
                    r_vapor=hp_r_vapor,
                    L_eva=hp_L_eva,
                    L_aba=0.0,
                    L_con=hp_L_con,
                    n_eva=hp_n_eva,
                    n_aba=0,
                    n_con=hp_n_con,
                    n_wick=hp_n_wick,
                    n_wall=hp_n_wall,
                    wall_mat=hp_wall_mat,
                    fluid_mat=hp_fluid_mat,
                    wick_struct_mat=hp_wick_mat,
                    porosity=porosity_hp,
                    fin_thickness=fin_thickness,
                    fin_height=fin_height,
                    n_fin_height=n_fin_height,
                    fin_wrap_ratio=fin_wrap_ratio,
                    emissivity=emissivity,
                    up_view_factor=up_view_factor,
                    down_view_factor=down_view_factor,
                    T_env=T_space,
                    initial_temp=HP_initial_temp
                )
                hp.hp.initialize_state()

                # 如果用户提供了外热流配置，则在构造阶段直接挂接到这根代表热管上。
                node_external_heat_config = self._resolve_node_external_heat_config(i)
                if node_external_heat_config is not None:
                    self._attach_external_heat_to_hp(hp, node_external_heat_config)
                    self._hp_external_heat_enabled[i] = True

                self.hp_units.append(hp)

                # ===== C. 估算热管阵列对集流环流动的阻塞形阻 =====
                A_proj_single = 2.0 * hp_r_out * hp_L_eva
                A_proj_total = A_proj_single

                if A_proj_total >= header_flow_area * 0.99:
                    raise ValueError(f"[{name}] 节点 {i} 的代表热管投影面积将集流环通道几乎完全堵塞。")

                sigma = (header_flow_area - A_proj_total) / header_flow_area
                if sigma <= 0.0:
                    raise ValueError(f"[{name}] 节点 {i} 的孔隙率 sigma 非法：{sigma}。")

                # 单排简化阻力模型。
                N_eff = 1 + 0.3 * (N_hp - 1)
                # 基于形阻系数、投影面积比和有效排数计算局部阻力损失系数
                K_loss_val = C_D * (A_proj_single / header_flow_area) * N_eff

                # ===== D. 建立流体控制体与热管蒸发段耦合 =====
                vol = self.fluid_channel.volumes[i]
                proxy = SingleVolumeProxy(vol, self.fluid_channel, i, hp_L_eva, N_hp)

                # 这里仍保留当前工程中的简化换热模型。
                # 如后续需要，可替换成基于横掠管束的 Nu 相关式。
                def constant_h_corr(Re, Pr, dummy=1.1):
                    k_f = proxy.material.conductivity(np.array([proxy.vol.T]), np.array([proxy.vol.P]))[0]
                    target_h = 10000.0
                    return (target_h * proxy.d_h) / k_f

                hp_peri_single = 2.0 * np.pi * hp_r_out
                cap_hp = hp.hp.get_boundary_node_capacitance('outer_eva')

                coupler_hp = FluidSolidCouple(
                    name=f"{name}_coupler_hp_{i}",
                    fluid=proxy,
                    solid_boundary_region=hp.hp.boundaries['outer_eva'],
                    heated_perimeter=hp_peri_single,
                    correlation_func=constant_h_corr,
                    solid_node_capacitance=cap_hp
                )
                self.coupler_hps.append(coupler_hp)

            # 不论该节点是否有热管，都要把出口局部阻力分配到对应 Junction。
            if i < n_nodes - 1:
                self.fluid_channel.internal_junctions[i].k_loss = K_loss_val
            else:
                self._outlet_k_loss = K_loss_val

            self._hp_k_loss_distribution[i] = K_loss_val

    @staticmethod
    def _build_external_heat_source(shape: tuple, external_heat_config: Dict[str, Any]) -> CompositeHeatSource:
        """构造与 HPwithFin 测试脚本一致的轨道外热流对象。"""
        composite_source = CompositeHeatSource(shape)

        if external_heat_config.get('use_embedded_table', False):
            composite_source.add_source(
                OrbitalTableHeatSource(
                    shape=shape,
                    table_ids=external_heat_config.get('table_ids', 1),
                    scale_factor=external_heat_config.get('table_scale_factor', 1.0),
                    offset=external_heat_config.get('table_offset', 0.0),
                    periodic=external_heat_config.get('table_periodic', True)
                )
            )
            return composite_source

        if external_heat_config.get('add_solar', False):
            composite_source.add_source(
                OrbitalHeatSource(
                    shape=shape,
                    solar_constant=external_heat_config.get('solar_constant', 1361.0),
                    orbit_height=external_heat_config.get('orbit_height', 800.0),
                    orbit_period=external_heat_config.get('orbit_period', 7644.0),
                    orbit_inclination=external_heat_config.get('orbit_inclination', 0.0),
                    surface_normal_angles=external_heat_config.get('surface_normal_angles', (0.0, 0.0))
                )
            )

        if external_heat_config.get('add_albedo', False):
            composite_source.add_source(
                AlbedoHeatSource(
                    shape=shape,
                    albedo_factor=external_heat_config.get('albedo_factor', 0.3),
                    solar_constant=external_heat_config.get('solar_constant', 1361.0)
                )
            )

        if external_heat_config.get('add_earth_ir', False):
            composite_source.add_source(
                EarthIRHeatSource(
                    shape=shape,
                    earth_ir_flux=external_heat_config.get('earth_ir_flux', 237.0)
                )
            )

        return composite_source

    def _resolve_node_external_heat_config(self, node_index: int) -> Optional[Dict[str, Any]]:
        """
        解析某个流体节点对应的热管外热流配置。

        支持三种传参方式：
        1. `external_heat_config=None`：整个 RingHP 不加载轨道外热流。
        2. 直接给一个统一配置字典：所有代表热管共用。
        3. 在统一配置基础上，通过 `node_configs` 或 `*_by_node` 做节点级覆盖。
        """
        if self.external_heat_config is None:
            return None

        config = copy.deepcopy(self.external_heat_config)

        node_configs = config.pop('node_configs', None)
        if node_configs is not None:
            if isinstance(node_configs, dict):
                node_override = node_configs.get(node_index)
            else:
                if len(node_configs) != self.n_header_nodes:
                    raise ValueError(
                        f"[{self.name}] node_configs 长度 ({len(node_configs)}) 必须与流体节点数 ({self.n_header_nodes}) 一致。"
                    )
                node_override = node_configs[node_index]
            if node_override is not None:
                config.update(copy.deepcopy(node_override))

        # 允许把某些参数按节点给成列表/数组，在这里抽取当前节点的那一份。
        per_node_keys = [
            'table_ids',
            'surface_normal_angles',
            'solar_constant',
            'orbit_height',
            'orbit_period',
            'orbit_inclination',
            'albedo_factor',
            'earth_ir_flux',
            'wall_illumination_factor',
            'fin_illuminated_area_scale',
            'fin_loading_mode',
            'table_scale_factor',
            'table_offset',
            'table_periodic',
            'add_solar',
            'add_albedo',
            'add_earth_ir',
            'use_embedded_table',
        ]
        for key in per_node_keys:
            by_node_key = f'{key}_by_node'
            if by_node_key in config:
                value_by_node = config.pop(by_node_key)
                if isinstance(value_by_node, dict):
                    if node_index in value_by_node:
                        config[key] = copy.deepcopy(value_by_node[node_index])
                else:
                    if len(value_by_node) != self.n_header_nodes:
                        raise ValueError(
                            f"[{self.name}] {by_node_key} 长度 ({len(value_by_node)}) 必须与流体节点数 ({self.n_header_nodes}) 一致。"
                        )
                    config[key] = copy.deepcopy(value_by_node[node_index])

        if not config:
            return None
        return config

    def _attach_external_heat_to_hp(self, hp: HPwithFin, external_heat_config: Dict[str, Any]):
        """
        将轨道外热流挂接到某根内部代表热管上。

        这部分逻辑与 `test_HP_with_external_heat_source.py` 保持一致，
        以保证单独测试 `HPwithFin` 与通过 `RingHP` 构造 `HPwithFin` 时行为一致。
        """
        outer_con_boundary = hp.hp.boundaries['outer_con']
        area_con = outer_con_boundary.area
        shape = outer_con_boundary.shape

        composite_source = self._build_external_heat_source(shape, external_heat_config)

        wall_illumination_factor = external_heat_config.get('wall_illumination_factor', 0.5)
        fin_illuminated_area_scale = external_heat_config.get('fin_illuminated_area_scale', 1.0)
        fin_loading_mode = external_heat_config.get('fin_loading_mode', 'lumped_root_area')

        wall_absorption_area = area_con * wall_illumination_factor
        fin_absorption_area = hp.get_fin_illuminated_area_array(fin_illuminated_area_scale)

        if fin_loading_mode == 'lumped_root_area':
            effective_boundary_area = wall_absorption_area + fin_absorption_area
            external_bc = ExternalHeatFluxBC(
                heat_source=composite_source,
                area_array=effective_boundary_area
            )
            outer_con_boundary.conditions.append(external_bc)
        elif fin_loading_mode == 'distributed_fin_absorption':
            external_bc = ExternalHeatFluxBC(
                heat_source=composite_source,
                area_array=wall_absorption_area
            )
            outer_con_boundary.conditions.append(external_bc)
            hp.set_fin_external_heat_source(
                composite_source,
                illuminated_area_scale=fin_illuminated_area_scale
            )
        else:
            raise ValueError(
                f"[{self.name}] 不支持的 fin_loading_mode={fin_loading_mode!r}。"
            )

        hp.configure_external_heat_accounting(
            composite_source,
            wall_area_array=wall_absorption_area,
            fin_area_array=fin_absorption_area
        )

    def get_solids(self) -> list:
        solids = [self.solid_header]
        for hp in self.hp_units:
            solids.extend(hp.get_solids())
        return solids

    def get_couplers(self) -> list:
        couplers = [self.coupler_header] + self.coupler_hps
        for hp in self.hp_units:
            couplers.extend(hp.get_couplers())
        return couplers

    def pre_step(self, dt: float, current_time: float):
        for hp in self.hp_units:
            if hasattr(hp, 'pre_step'):
                hp.pre_step(dt, current_time)

    def post_step(self, dt: float, current_time: float):
        for hp in self.hp_units:
            if hasattr(hp, 'post_step'):
                hp.post_step(dt, current_time)

    @property
    def outlet_k_loss(self) -> float:
        """返回最后一个控制体出口对应的局部阻力系数。"""
        return self._outlet_k_loss

    @property
    def hp_k_loss_distribution(self) -> np.ndarray:
        """返回各流体节点对应的热管阵列局部阻力系数分布。"""
        return np.array(self._hp_k_loss_distribution, copy=True)

    def get_total_heat_rejection(self) -> float:
        """返回所有代表热管当前的总向外散热量。"""
        total_q = 0.0
        for hp in self.hp_units:
            _, q_con_dist = hp.get_heat_rejection_distribution()
            total_q += float(np.sum(q_con_dist))
        return total_q

    def get_total_external_heat_absorption(self, current_time: float) -> float:
        """返回所有代表热管当前吸收的轨道外热流总功率。"""
        total_q = 0.0
        for hp in self.hp_units:
            _, _, q_abs_dist = hp.get_external_heat_absorption_distribution(current_time)
            total_q += float(np.sum(q_abs_dist))
        return total_q

    def get_hp_status_summary(self, current_time: Optional[float] = None) -> Dict[str, Any]:
        """
        返回 RingHP 内部热管阵列的汇总信息，便于调试和后处理。
        """
        summary: Dict[str, Any] = {
            'n_header_nodes': self.n_header_nodes,
            'n_hp_units': len(self.hp_units),
            'hp_presence_mask': np.array(self._hp_presence_mask, copy=True),
            'hp_external_heat_enabled': np.array(self._hp_external_heat_enabled, copy=True),
            'hp_k_loss_distribution': np.array(self._hp_k_loss_distribution, copy=True),
            'gross_heat_rejection': self.get_total_heat_rejection(),
        }
        if current_time is not None:
            summary['external_heat_absorption'] = self.get_total_external_heat_absorption(current_time)
            summary['net_heat_rejection'] = summary['gross_heat_rejection'] - summary['external_heat_absorption']
        return summary
