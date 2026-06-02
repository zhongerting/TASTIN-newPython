import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Literal, List
import math

from Components.tec_electric import joule_power_from_electric_field


@dataclass
class TFEGeometry:
    """热离子燃料元件几何参数 (单位: 国际单位制 m)"""
    # 半径参数 (从内到外)
    r_pellet_inner: float  # 芯块内半径 (如果是实心则为0)
    r_pellet_outer: float  # 芯块外半径
    r_fission_gas_outer: float  # 裂变气隙外半径 (氙气)
    r_emitter_outer: float  # 发射极外半径
    # 跨越电极间隙
    r_collector_inner: float  # 接收极内半径 (注意这里跳过了电极间隙)
    r_collector_outer: float  # 接收极外半径
    # 跨越氦气隙 内套管
    r_inner_clad_inner: float  # 内套管内径 ( = r_collector_outer + 氦气隙厚度)
    r_inner_clad_outer: float
    # 冷却剂流道
    r_coolant_inner: float  # 冷却剂内半径
    r_coolant_outer: float  # 冷却剂外半径
    # 外套管
    r_outer_clad_outer: float  # 外套管外半径
    # 跨越 CO2 间隙
    r_moderator_inner: float
    r_moderator_outer: float  # 等效慢化剂外半径 (跳过了CO2间隙)

    # 轴向参数
    height: float  # 燃料元件有效高度

@dataclass
class TFEMeshParams:
    """热离子燃料元件径向与轴向网格划分参数"""
    n_axial: int               # 轴向节点数
    # 径向节点数分配
    n_r_pellet: int
    n_r_emitter: int
    n_r_collector: int
    n_r_inner_clad: int
    n_r_outer_clad: int
    n_r_moderator: int

@dataclass
class NeutronicData:
    """中子学与热源数据域 (现改为标量参数)"""
    total_power: float = 0.0          # 当前接收到的总功率 [W]
    _total_power_old: float = 0.0     # 上一时间步/迭代步的总功率 (用于亚松弛)


@dataclass
class ElectricFieldData:
    """电场与焦耳热计算数据域 (长度均为 n_axial)"""
    # 电场分布
    emitter_voltage: np.ndarray = field(default_factory=lambda: np.array([]))  # [V]
    emitter_resistivity: np.ndarray = field(default_factory=lambda: np.array([]))  # [Ohm*m]
    collector_voltage: np.ndarray = field(default_factory=lambda: np.array([]))  # [V]
    collector_resistivity: np.ndarray = field(default_factory=lambda: np.array([]))  # [Ohm*m]

    # 局部电流密度 [A/m^2]
    current_density: np.ndarray = field(default_factory=lambda: np.array([]))

    # [内部计算结果缓存] 映射后计算得到的总焦耳热密度 [W/m^3]
    emitter_joule_heat: np.ndarray = field(default_factory=lambda: np.array([]))
    collector_joule_heat: np.ndarray = field(default_factory=lambda: np.array([]))

@dataclass
class PlasmaCouplingData:
    """等离子体与界面热流数据域 (长度均为 n_axial)"""
    # TEC 物理状态参数 (用于计算热流或作为外部电模块反馈)
    emitter_work_function: np.ndarray = field(default_factory=lambda: np.array([]))  # 发射极功函数 [eV]
    collector_work_function: np.ndarray = field(default_factory=lambda: np.array([]))  # 接收极功函数 [eV]
    barrier_voltage_drop: np.ndarray = field(default_factory=lambda: np.array([]))  # 阻塞电压 Vd [V]
    emitter_temperature: np.ndarray = field(default_factory=lambda: np.array([]))  # 发射极温度 TE [K]

    # 外部计算传入：发射极表面的电子冷却热流 [W/m^2]
    electron_cooling_flux: np.ndarray = field(default_factory=lambda: np.array([]))

    # 外部计算传入：接收极表面的电子加热热流 [W/m^2]
    electron_heating_flux: np.ndarray = field(default_factory=lambda: np.array([]))

    # (注意：radiation_flux 已被移除，由内部 TECCouple2D 负责计算)

@dataclass
class BoundaryConditionData:
    """外部强加的边界条件数据域 (长度均为 n_axial)"""
    # 等效慢化剂的外部温度分布 [K]，用于驱动 GapCouple2D 或变定温边界
    moderator_temperature: np.ndarray = field(default_factory=lambda: np.array([]))


@dataclass
class GapConfig:
    """
    单一气隙的独立建模配置类
    """
    mode: Literal['simplified', 'meshed']

    # [新增] 间隙等效换热系数 h [W/(m^2·K)]
    # 用于 'simplified' 模式，底层将自动转换为气隙等效导热系数 k = h * d
    h_eq: float = 0.0

    # 气隙填充材料对象 (用于 'meshed' 模式)
    material: Optional[Any] = None

    # 辐射换热参数
    emissivity_inner: float = 0.8
    emissivity_outer: float = 0.8

    n_radial_nodes: int = 0

    def __post_init__(self):
        if self.mode == 'meshed' and self.n_radial_nodes <= 0:
            raise ValueError("气隙配置错误：meshed 模式必须提供 n_radial_nodes！")
        if self.mode == 'meshed' and self.material is None:
            raise ValueError("气隙配置错误：meshed 模式必须提供真实的 material 物性对象！")


from Components.BaseComponent import BaseComponent
from Components.basicComponents.Fuel import Fuel
from Components.basicComponents.Electord import Emitter, Collector
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Couplers import SolidSolidCouple2D, GapCouple2D, TECCouple2D, FluidSolidCouple
from Correlations.Correlations import nu_ringpipe
from Materials.Solids.CompositePelletMaterial import CompositePelletMaterial


class TFEUnit(BaseComponent):
    """
    热离子燃料元件 (Thermionic Fuel Element) 组件
    集成芯块、发射极、接收极、套管、慢化剂及内部各种热传递与电热耦合机制。
    """

    def __init__(self,
                 name: str,
                 geometry: TFEGeometry,
                 mesh_params: TFEMeshParams,
                 materials: Dict[str, Any],  # 包含 'UO2', 'Mo', 'Nb', 'ZrH' 等物性对象
                 coolant_channel: Any,
                 # --- [新增] 独立的气隙配置参数 ---
                 fission_gas_config: GapConfig,  # 裂变气隙配置 (Pellet - Emitter)
                 tec_gap_config: GapConfig,  # 铯气隙/极间隙配置 (Emitter - Collector)
                 he_gap_config: GapConfig,  # 氦气隙配置 (Collector - InnerClad)
                 co2_gap_config: GapConfig,  # CO2气隙配置 (OuterClad - Moderator)
                 # [新增] 功率分配与轴向分布参数
                 power_fraction: float = 1.0,
                 axial_power_profile: Optional[np.ndarray] = None,
                 # --- [TEASA 新增] 轴向网格非均匀配置与接触热阻 ---
                 axial_length_allocation: Optional[List[float]] = None,
                 axial_node_allocation: Optional[List[int]] = None,
                 axial_contact_resistance: float = 0.0,
                 strict_adiabatic_single_tfe: bool = False):  # 外部传入的冷却剂流体实例

        super().__init__(name)

        # 0. 准备输入接口与参数保存
        self.geom = geometry
        self.mesh = mesh_params
        self.materials = materials
        self.coolant = coolant_channel

        self.power_fraction = power_fraction
        self.axial_power_profile = axial_power_profile
        self.axial_contact_resistance = axial_contact_resistance
        self.strict_adiabatic_single_tfe = bool(strict_adiabatic_single_tfe)

        # 将四个间隙的配置存入一个内部字典，方便后续按名字调用和迭代
        self.gap_configs = {
            'fission_gas': fission_gas_config,
            'tec_gap': tec_gap_config,
            'he_gap': he_gap_config,
            'co2_gap': co2_gap_config
        }

        # 获取轴向节点数，用于初始化数组长度
        nz = self.mesh.n_axial
        H = self.geom.height

        # =========================================================
        # [核心架构升级 0 & 2] 生成公用 y_faces 与 交界面索引
        # =========================================================
        if axial_length_allocation is not None and axial_node_allocation is not None:
            if len(axial_length_allocation) != 3 or len(axial_node_allocation) != 3:
                raise ValueError(
                    "[TFEUnit] axial_length_allocation 和 axial_node_allocation 必须为长度为 3 的数组 [下, 活性, 上]")

            l_lower, l_act, l_upper = axial_length_allocation
            n_lower, n_act, n_upper = axial_node_allocation

            if n_lower + n_act + n_upper != nz:
                raise ValueError(
                    f"[TFEUnit] 分段节点总和 ({n_lower + n_act + n_upper}) 必须等于网格总轴向节点数 nz ({nz})")

            # 生成三段局部分布的界面坐标 (利用 concatenate 和切片 [1:] 去除重复的拼接点)
            y_faces_lower = np.linspace(0, l_lower, n_lower + 1)
            y_faces_active = np.linspace(l_lower, l_lower + l_act, n_act + 1)[1:]
            y_faces_upper = np.linspace(l_lower + l_act, H, n_upper + 1)[1:]

            self.common_y_faces = np.concatenate([y_faces_lower, y_faces_active, y_faces_upper])
            self.n_lower, self.n_active, self.n_upper = n_lower, n_act, n_upper
            self.has_reflectors = (n_lower > 0 or n_upper > 0)
        else:
            # 向下兼容：如果没有传入非均匀分配参数，回退到完全均匀网格
            self.common_y_faces = np.linspace(0.0, H, nz + 1)
            self.n_lower, self.n_active, self.n_upper = 0, nz, 0
            self.has_reflectors = False

        # 提前计算存在接触热阻的交界面索引 (将传入极简 Fuel)
        self.contact_resistance_interfaces = []
        if self.has_reflectors and self.axial_contact_resistance > 1e-12:
            if self.n_lower > 0:
                self.contact_resistance_interfaces.append(self.n_lower - 1)
            if self.n_upper > 0:
                self.contact_resistance_interfaces.append(self.n_lower + self.n_active - 1)

        # 功率分布闭环保护：如果有反射层但未提供功率分布数组，则强制仅在活性区发热
        if axial_power_profile is None and self.has_reflectors:
            profile = np.zeros(nz)
            if self.n_active > 0:
                profile[self.n_lower: self.n_lower + self.n_active] = 1.0
            self.axial_power_profile = profile
        else:
            self.axial_power_profile = axial_power_profile

        # 实例化内部数据存储 (初始化为零数组，保证内存连续与类型安全)
        self.neutronic_data = NeutronicData()
        self.electric_data = ElectricFieldData(
            emitter_voltage=np.zeros(nz),
            emitter_resistivity=np.zeros(nz),
            collector_voltage=np.zeros(nz),
            collector_resistivity=np.zeros(nz),
            current_density=np.zeros(nz),
            emitter_joule_heat=np.zeros(nz),
            collector_joule_heat=np.zeros(nz)
        )
        self.plasma_data = PlasmaCouplingData(
            emitter_work_function=np.zeros(nz),
            collector_work_function=np.zeros(nz),
            barrier_voltage_drop=np.zeros(nz),
            emitter_temperature=np.zeros(nz),
            electron_cooling_flux=np.zeros(nz),
            electron_heating_flux=np.zeros(nz)
        )
        self.boundary_data = BoundaryConditionData(
            moderator_temperature=np.zeros(nz)
        )

        # 内部组件字典，便于管理
        self.solids = {}
        self.couplers = {}

        # 顺序执行初始化构建
        self._build_solids()
        self._build_couplers()

    def _build_solids(self):
        """1. 建立固体实例 (内部实例化)"""

        H = self.geom.height
        Ny = self.mesh.n_axial

        # 先进行所有固体域建模
        # ---------------------------------------------------------
        # (1) 芯块 (Pellet) - 使用 Fuel 类
        # ---------------------------------------------------------
        delta_pellet = self.geom.r_pellet_outer - self.geom.r_pellet_inner
        # 追加 y_faces 统一网格
        mesh_pellet = Mesh2D(x_dim=delta_pellet, n_x=self.mesh.n_r_pellet,
                             y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                             geometry_type='cylindrical', inner_radius=self.geom.r_pellet_inner)

        # 动态装配复合芯块材料
        if self.has_reflectors:
            beo_mat = self.materials.get('BerylliumOxide') or self.materials.get('BeO')
            if beo_mat is None:
                raise ValueError("缺少反射层物性 'BerylliumOxide'！请在 materials 字典中提供。")

            pellet_mat = CompositePelletMaterial(
                fuel_mat=self.materials.get('UO2'),
                reflector_mat=beo_mat,
                shape_nodes=(self.mesh.n_r_pellet, Ny),
                n_lower=self.n_lower,
                n_active=self.n_active,
                n_upper=self.n_upper
            )
        else:
            pellet_mat = self.materials.get('UO2')

        self.solids['pellet'] = Fuel(
            name=f"{self.name}_Pellet",
            mesh=mesh_pellet,
            material=pellet_mat,
            initial_temp=743.15,
            power_fraction=self.power_fraction,
            axial_power_profile=self.axial_power_profile,
            # 传入之前算好的热阻接口与数值
            contact_resistance_interfaces=self.contact_resistance_interfaces,
            axial_contact_resistance=self.axial_contact_resistance
        )

        # ---------------------------------------------------------
        # (2) 发射极 (Emitter) - 使用 Emitter 类
        # ---------------------------------------------------------
        # 发射极内半径紧贴裂变气隙外围
        delta_emitter = self.geom.r_emitter_outer - self.geom.r_fission_gas_outer
        mesh_emitter = Mesh2D(x_dim=delta_emitter, n_x=self.mesh.n_r_emitter,
                              y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                              geometry_type='cylindrical', inner_radius=self.geom.r_fission_gas_outer)

        self.solids['emitter'] = Emitter(
            name=f"{self.name}_Emitter",
            mesh=mesh_emitter,
            material=self.materials.get('MoNb'),
            initial_temp=743.15
        )

        # ---------------------------------------------------------
        # (3) 接收极 (Collector) - 使用 Collector 类
        # ---------------------------------------------------------
        # 跨越电极间隙，内半径由 geom.r_collector_inner 给出
        delta_collector = self.geom.r_collector_outer - self.geom.r_collector_inner
        mesh_collector = Mesh2D(x_dim=delta_collector, n_x=self.mesh.n_r_collector,
                                y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                                geometry_type='cylindrical', inner_radius=self.geom.r_collector_inner)

        self.solids['collector'] = Collector(
            name=f"{self.name}_Collector",
            mesh=mesh_collector,
            material=self.materials.get('Molybdenum'),
            initial_temp=743.15
        )

        # ---------------------------------------------------------
        # (4) 冷却剂内套管 (Inner Clad)
        # ---------------------------------------------------------
        # 严格使用 r_inner_clad_inner，保留氦气隙的物理空间跳跃
        delta_iclad = self.geom.r_inner_clad_outer - self.geom.r_inner_clad_inner
        mesh_iclad = Mesh2D(x_dim=delta_iclad, n_x=self.mesh.n_r_inner_clad,
                            y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                            geometry_type='cylindrical',
                            inner_radius=self.geom.r_inner_clad_inner)  # [修正]

        self.solids['inner_clad'] = HeatConduction2D(
            name=f"{self.name}_InnerClad",
            mesh=mesh_iclad,
            material=self.materials.get('StainlessSteel'),
            initial_temp=743.0
        )

        # ---------------------------------------------------------
        # (6) 冷却剂外套管 (Outer Clad)
        # ---------------------------------------------------------
        delta_oclad = self.geom.r_outer_clad_outer - self.geom.r_coolant_outer
        mesh_oclad = Mesh2D(x_dim=delta_oclad, n_x=self.mesh.n_r_outer_clad,
                            y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                            geometry_type='cylindrical',
                            inner_radius=self.geom.r_coolant_outer)

        self.solids['outer_clad'] = HeatConduction2D(
            name=f"{self.name}_OuterClad",
            mesh=mesh_oclad,
            material=self.materials.get('StainlessSteel'),
            initial_temp=743.0
        )

        if not self.strict_adiabatic_single_tfe:
            # ---------------------------------------------------------
            # (7) 等效慢化剂 (Moderator)
            # ---------------------------------------------------------
            # 严格使用 r_moderator_inner，保留 CO2 间隙的物理空间跳跃
            delta_mod = self.geom.r_moderator_outer - self.geom.r_moderator_inner
            mesh_mod = Mesh2D(x_dim=delta_mod, n_x=self.mesh.n_r_moderator,
                              y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                              geometry_type='cylindrical',
                              inner_radius=self.geom.r_moderator_inner)  # [修正]

            self.solids['moderator'] = HeatConduction2D(
                name=f"{self.name}_Moderator",
                mesh=mesh_mod,
                material=self.materials.get('ZrH'),
                initial_temp=743.0
            )

        # ==========================================
        # 动态保真度气隙层 (当 mode == 'meshed' 时构建实体)
        # ==========================================

        # [A] 裂变气隙
        cfg_fg = self.gap_configs['fission_gas']
        if cfg_fg.mode == 'meshed':
            d_fg = self.geom.r_fission_gas_outer - self.geom.r_pellet_outer
            mesh_fg = Mesh2D(x_dim=d_fg, n_x=cfg_fg.n_radial_nodes, y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                             geometry_type='cylindrical', inner_radius=self.geom.r_pellet_outer)
            self.solids['fission_gas_solid'] = HeatConduction2D(
                name=f"{self.name}_FissionGasMesh", mesh=mesh_fg, material=cfg_fg.material, initial_temp=1550.0)

        # [B] TEC极间隙 (铯)
        cfg_tec = self.gap_configs['tec_gap']
        if cfg_tec.mode == 'meshed':
            d_tec = self.geom.r_collector_inner - self.geom.r_emitter_outer
            mesh_tec = Mesh2D(x_dim=d_tec, n_x=cfg_tec.n_radial_nodes, y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                              geometry_type='cylindrical', inner_radius=self.geom.r_emitter_outer)
            self.solids['tec_gap_solid'] = HeatConduction2D(
                name=f"{self.name}_TECGapMesh", mesh=mesh_tec, material=cfg_tec.material, initial_temp=1200.0)

        # [C] 氦气隙
        cfg_he = self.gap_configs['he_gap']
        if cfg_he.mode == 'meshed':
            d_he = self.geom.r_inner_clad_inner - self.geom.r_collector_outer
            mesh_he = Mesh2D(x_dim=d_he, n_x=cfg_he.n_radial_nodes, y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                             geometry_type='cylindrical', inner_radius=self.geom.r_collector_outer)
            self.solids['he_gap_solid'] = HeatConduction2D(
                name=f"{self.name}_HeGapMesh", mesh=mesh_he, material=cfg_he.material, initial_temp=850.0)

        # [D] CO2气隙
        cfg_co2 = self.gap_configs['co2_gap']
        if not self.strict_adiabatic_single_tfe and cfg_co2.mode == 'meshed':
            d_co2 = self.geom.r_moderator_inner - self.geom.r_outer_clad_outer
            mesh_co2 = Mesh2D(x_dim=d_co2, n_x=cfg_co2.n_radial_nodes, y_dim=H, n_y=Ny, y_faces=self.common_y_faces,
                              geometry_type='cylindrical', inner_radius=self.geom.r_outer_clad_outer)
            self.solids['co2_gap_solid'] = HeatConduction2D(
                name=f"{self.name}_CO2GapMesh", mesh=mesh_co2, material=cfg_co2.material, initial_temp=750.0)

    def _build_couplers(self):
        """2. 建立固体-固体，固体-流体耦合"""
        # 在这里建立各种 Couple，并装配到实体的 Boundary 上

        # 提取各个实体，简化后续代码引用
        s_pellet = self.solids['pellet']
        s_emitter = self.solids['emitter']
        s_collector = self.solids['collector']
        s_iclad = self.solids['inner_clad']
        s_oclad = self.solids['outer_clad']

        # =========================================================
        # [Interface A] 芯块 -> 发射极 (裂变气隙)
        # =========================================================
        cfg_fg = self.gap_configs['fission_gas']
        w_fg = self.geom.r_fission_gas_outer - self.geom.r_pellet_outer

        if cfg_fg.mode == 'simplified':
            self.couplers['pellet_emitter_gap'] = GapCouple2D(
                obj1=s_pellet, obj2=s_emitter, direction='right', gap_width=w_fg,
                gas_conductivity=cfg_fg.h_eq * w_fg,  # 计算等效导热系数
                emissivity1=cfg_fg.emissivity_inner, emissivity2=cfg_fg.emissivity_outer
            )
        else:  # 'meshed'
            s_fg = self.solids['fission_gas_solid']
            # 1. 固体网格传导
            self.couplers['pellet_fg_cond'] = SolidSolidCouple2D(obj1=s_pellet, obj2=s_fg, direction='right')
            self.couplers['fg_emitter_cond'] = SolidSolidCouple2D(obj1=s_fg, obj2=s_emitter, direction='right')
            # 2. 纯辐射并联通道 (屏蔽气体导热)
            self.couplers['pellet_emitter_rad'] = GapCouple2D(
                obj1=s_pellet, obj2=s_emitter, direction='right', gap_width=w_fg,
                gas_conductivity=0.0,
                emissivity1=cfg_fg.emissivity_inner, emissivity2=cfg_fg.emissivity_outer
            )

        # =========================================================
        # [Interface B] 发射极 -> 接收极 (TEC 极间隙)
        # =========================================================
        # 无论哪种模式，TECCouple2D 是必须存在的，因为它要接收等离子体热流！
        cfg_tec = self.gap_configs['tec_gap']
        w_tec = self.geom.r_collector_inner - self.geom.r_emitter_outer

        if cfg_tec.mode == 'simplified':
            self.couplers['tec_couple'] = TECCouple2D(
                obj1=s_emitter, obj2=s_collector, direction='right', gap_width=w_tec,
                gas_conductivity=cfg_tec.h_eq * w_tec,
                emissivity1=cfg_tec.emissivity_inner, emissivity2=cfg_tec.emissivity_outer
            )
        else:  # 'meshed'
            s_tec = self.solids['tec_gap_solid']
            # 1. 固体网格传导
            self.couplers['emitter_tec_cond'] = SolidSolidCouple2D(obj1=s_emitter, obj2=s_tec, direction='right')
            self.couplers['tec_collector_cond'] = SolidSolidCouple2D(obj1=s_tec, obj2=s_collector,
                                                                     direction='right')
            # 2. 核心电热耦合通道 (屏蔽气体导热，只保留辐射与电子传热接口)
            self.couplers['tec_couple'] = TECCouple2D(
                obj1=s_emitter, obj2=s_collector, direction='right', gap_width=w_tec,
                gas_conductivity=0.0,
                emissivity1=cfg_tec.emissivity_inner, emissivity2=cfg_tec.emissivity_outer
            )

        # =========================================================
        # [Interface C] 接收极 -> 内套管 (氦气隙)
        # =========================================================
        cfg_he = self.gap_configs['he_gap']
        w_he = self.geom.r_inner_clad_inner - self.geom.r_collector_outer

        if cfg_he.mode == 'simplified':
            self.couplers['collector_iclad_gap'] = GapCouple2D(
                obj1=s_collector, obj2=s_iclad, direction='right', gap_width=w_he,
                gas_conductivity=cfg_he.h_eq * w_he,
                emissivity1=cfg_he.emissivity_inner, emissivity2=cfg_he.emissivity_outer
            )
        else:  # 'meshed'
            s_he = self.solids['he_gap_solid']
            self.couplers['collector_he_cond'] = SolidSolidCouple2D(obj1=s_collector, obj2=s_he, direction='right')
            self.couplers['he_iclad_cond'] = SolidSolidCouple2D(obj1=s_he, obj2=s_iclad, direction='right')
            self.couplers['collector_iclad_rad'] = GapCouple2D(
                obj1=s_collector, obj2=s_iclad, direction='right', gap_width=w_he,
                gas_conductivity=0.0,
                emissivity1=cfg_he.emissivity_inner, emissivity2=cfg_he.emissivity_outer
            )

        # =========================================================
        # [Interface D] 冷却剂流道流固耦合
        # =========================================================
        if self.coolant is not None:
            r_fluid_in = self.geom.r_inner_clad_outer
            r_fluid_out = self.geom.r_coolant_outer

            def nu_ringpipe_adapter(Re: float, Pr: float, pd_dummy: float = 1.1) -> float:
                return nu_ringpipe(R_out=r_fluid_out, R_in=r_fluid_in, Re=Re, Pr=Pr)

            # 强制初始化套管的固体状态以获取热容
            s_iclad.initialize_state()
            s_oclad.initialize_state()
            cap_inner = s_iclad.get_boundary_node_capacitance('right')
            cap_outer = s_oclad.get_boundary_node_capacitance('left')

            self.couplers['iclad_coolant_fsc'] = FluidSolidCouple(
                name=f"{self.name}_iclad_coolant_couple", fluid=self.coolant,
                solid_boundary_region=s_iclad.boundaries['right'],
                heated_perimeter=2.0 * math.pi * r_fluid_in,
                correlation_func=nu_ringpipe_adapter, solid_node_capacitance=cap_inner
            )

            self.couplers['oclad_coolant_fsc'] = FluidSolidCouple(
                name=f"{self.name}_oclad_coolant_couple", fluid=self.coolant,
                solid_boundary_region=s_oclad.boundaries['left'],
                heated_perimeter=2.0 * math.pi * r_fluid_out,
                correlation_func=nu_ringpipe_adapter, solid_node_capacitance=cap_outer
            )

        if self.strict_adiabatic_single_tfe:
            s_oclad.boundaries['right'].clear_conditions()
            return

        # =========================================================
        # [Interface E] 外套管 -> 慢化剂 (CO2 气隙)
        # =========================================================
        s_mod = self.solids['moderator']
        cfg_co2 = self.gap_configs['co2_gap']
        w_co2 = self.geom.r_moderator_inner - self.geom.r_outer_clad_outer

        if cfg_co2.mode == 'simplified':
            self.couplers['oclad_mod_gap'] = GapCouple2D(
                obj1=s_oclad, obj2=s_mod, direction='right', gap_width=w_co2,
                gas_conductivity=cfg_co2.h_eq * w_co2,
                emissivity1=cfg_co2.emissivity_inner, emissivity2=cfg_co2.emissivity_outer
            )
        else:  # 'meshed'
            s_co2 = self.solids['co2_gap_solid']
            self.couplers['oclad_co2_cond'] = SolidSolidCouple2D(obj1=s_oclad, obj2=s_co2, direction='right')
            self.couplers['co2_mod_cond'] = SolidSolidCouple2D(obj1=s_co2, obj2=s_mod, direction='right')
            self.couplers['oclad_mod_rad'] = GapCouple2D(
                obj1=s_oclad, obj2=s_mod, direction='right', gap_width=w_co2,
                gas_conductivity=0.0,
                emissivity1=cfg_co2.emissivity_inner, emissivity2=cfg_co2.emissivity_outer
            )

        # =========================================================
        # [Interface F] 虚拟慢化剂外边界 (与全局慢化剂的纽带)
        # =========================================================
        s_mod = self.solids['moderator']
        # 【物理修正】：绝对不能用 1e-8！这会导致严重的 ODE 刚性崩溃（Stiff ODE）。
        # 我们使用一个合理的等效接触换热系数 (比如 h=500 W/m^2K) 来“软化”这个边界。
        A_surf = s_mod.boundaries['right'].area
        R_soft = 1.0 / (500.0 * A_surf)  # 换算为绝对热阻 [K/W]

        self.mod_outer_bc = s_mod.boundaries['right'].add_resistance_condition(
            T_ext=self.boundary_data.moderator_temperature,
            R_ext=R_soft
        )

    # ==========================================
    # 专用参数更新接口
    # ==========================================
    def update_neutronic_power(self, p_total: float, p_fiss: float = 0.0, p_decay: float = 0.0, alpha: float = 1.0):
        """
        接收外部点堆传入的反应堆总功率，计算后映射至芯块网格。

        :param p_total: 反应堆当前总功率 [W]
        :param p_fiss: 裂变功率分量 [W] (预留给深度燃耗或多物理耦合)
        :param p_decay: 衰变功率分量 [W] (预留)
        :param alpha: 亚松弛因子，防止功率阶跃冲击数值求解器
        """
        # 1. 对标量总功率进行亚松弛计算
        p_total_new = alpha * p_total + (1.0 - alpha) * self.neutronic_data._total_power_old

        # 2. 更新内部状态缓存
        self.neutronic_data._total_power_old = p_total_new
        self.neutronic_data.total_power = p_total_new

        # 3. 调用 Fuel 类自带的高效接口下发功率
        pellet = self.solids['pellet']
        pellet.set_nuclear_power(p_fiss=p_fiss, p_decay=p_decay, p_total=p_total_new)

        # 4. TEASA 强制保护机制 (Double Check)
        # 虽然 Fuel 类中重写了 Q_source，但为了防止被底层 SystemManager 的源项清理逻辑覆盖，
        # 必须显式激活外源缓冲标记。(参考 Electord.py 的安全逻辑)
        pellet.use_external_source_buffer = True

    def update_electric_fields(self,
                               dU_emit: np.ndarray, rho_emit: np.ndarray,
                               dU_coll: np.ndarray, rho_coll: np.ndarray,
                               alpha: float = 1.0):
        """
        更新电场分布并计算焦耳热映射 (严格参考 TECPair 的设值逻辑)
        :param dU_emit, dU_coll: 发射极/接收极的轴向网格电压降 [V] (长度 n_axial)
        :param rho_emit, rho_coll: 发射极/接收极的电阻率 [Ohm*m] (长度 n_axial)
        :param alpha: 亚松弛因子，提高松耦合电热迭代的稳定性
        """
        node_lengths = np.diff(
            np.asarray(getattr(self, 'common_y_faces', self.solids['emitter'].mesh.y_faces), dtype=float)
        )
        self.update_electric_field_sources(
            E_emit=np.asarray(dU_emit, dtype=float) / node_lengths,
            rho_emit=rho_emit,
            E_coll=np.asarray(dU_coll, dtype=float) / node_lengths,
            rho_coll=rho_coll,
            alpha=alpha,
        )

    def update_electric_field_sources(self,
                                      E_emit: np.ndarray, rho_emit: np.ndarray,
                                      E_coll: np.ndarray, rho_coll: np.ndarray,
                                      alpha: float = 1.0):
        """
        Apply axial electric fields [V/m] as Joule sources.

        The field must already be computed on the axial node centers. This keeps
        non-uniform-grid differencing outside the source mapping step.
        """
        emitter = self.solids['emitter']
        collector = self.solids['collector']

        E_emit = np.asarray(E_emit, dtype=float)
        E_coll = np.asarray(E_coll, dtype=float)
        rho_emit = np.asarray(rho_emit, dtype=float)
        rho_coll = np.asarray(rho_coll, dtype=float)

        self.electric_data.emitter_voltage[:] = E_emit
        self.electric_data.emitter_resistivity[:] = rho_emit
        self.electric_data.collector_voltage[:] = E_coll
        self.electric_data.collector_resistivity[:] = rho_coll

        q_watts_e_flat, _ = joule_power_from_electric_field(
            E_emit,
            rho_emit,
            emitter.vols_flat,
            emitter.mesh.shape_nodes,
        )
        q_watts_c_flat, _ = joule_power_from_electric_field(
            E_coll,
            rho_coll,
            collector.vols_flat,
            collector.mesh.shape_nodes,
        )

        old_e = self.electric_data.emitter_joule_heat
        old_c = self.electric_data.collector_joule_heat
        if old_e.shape != q_watts_e_flat.shape:
            old_e = np.zeros_like(q_watts_e_flat)
        if old_c.shape != q_watts_c_flat.shape:
            old_c = np.zeros_like(q_watts_c_flat)

        q_watts_e_new = alpha * q_watts_e_flat + (1.0 - alpha) * old_e
        q_watts_c_new = alpha * q_watts_c_flat + (1.0 - alpha) * old_c

        self.electric_data.emitter_joule_heat = q_watts_e_new.copy()
        self.electric_data.collector_joule_heat = q_watts_c_new.copy()

        emitter.set_joule_heating(q_watts_e_new)
        collector.set_joule_heating(q_watts_c_new)

    def update_plasma_flux(self, q_e_flux: np.ndarray, q_c_flux: np.ndarray, alpha: float = 1.0):
        """
        施加等离子体放电产生的面热流密度 (严格参考 TECPair 逻辑)

        :param q_e_flux: 发射极外表面电子冷却热流 [W/m^2] (通常带走热量，需注意外部正负号定义)
        :param q_c_flux: 接收极内表面电子加热热流 [W/m^2]
        :param alpha: 亚松弛因子
        """
        # 1. 结合内部缓存对热流密度进行亚松弛
        q_e_new = alpha * q_e_flux + (1.0 - alpha) * self.plasma_data.electron_cooling_flux
        q_c_new = alpha * q_c_flux + (1.0 - alpha) * self.plasma_data.electron_heating_flux

        # 刷新内部数据状态
        self.plasma_data.electron_cooling_flux[:] = q_e_new
        self.plasma_data.electron_heating_flux[:] = q_c_new

        # 2. 获取网格面的面积
        # (因为外部电路计算面热流 q 往往是基于 Emitter 外表面积的)
        A_emit_surf = self.solids['emitter'].boundaries['right'].area
        A_coll_surf = self.solids['collector'].boundaries['left'].area

        # 3. 转换为绝对传输功率 [W]
        # (保证能量守恒：W/m^2 * m^2 = Watts)
        Q_e_watts = q_e_new * A_emit_surf
        Q_c_watts = q_c_new * A_emit_surf

        # 4. 严格调用 TECCouple2D 提供的专属设值接口
        self.plasma_area_diagnostics = {
            "basis": "ThermoCalc current density uses emitter area; collector electron power is converted on the same basis.",
            "sideAreaE_m2": np.asarray(A_emit_surf, dtype=float).copy(),
            "sideAreaC_m2": np.asarray(A_coll_surf, dtype=float).copy(),
            "emitter_power_emitter_area_w": np.asarray(Q_e_watts, dtype=float).copy(),
            "collector_power_emitter_area_w": np.asarray(Q_c_watts, dtype=float).copy(),
            "collector_power_collector_area_w": np.asarray(q_c_new * A_coll_surf, dtype=float).copy(),
            "collector_area_power_delta_w": np.asarray(q_c_new * (A_coll_surf - A_emit_surf), dtype=float).copy(),
        }
        self.couplers['tec_couple'].set_tec_sources(Q_emitter=Q_e_watts, Q_collector=Q_c_watts)

    # ==========================================
    # 系统管理与生命周期重载 (对接 SystemManager)
    # ==========================================

    def pre_step(self, dt: float, current_time: float):
        """[生命周期钩子] 在求解器推进前更新动态边界"""
        if hasattr(self, 'mod_outer_bc'):
            self.mod_outer_bc.update_params(T_ext=self.boundary_data.moderator_temperature)

    @staticmethod
    def _capture_solid_source_state(solid):
        state = {}
        if hasattr(solid, 'Q_source'):
            state['Q_source_ref'] = solid.Q_source
            state['Q_source_value'] = np.asarray(solid.Q_source, dtype=float).copy()
        if hasattr(solid, 'Q_vol'):
            state['Q_vol_value'] = np.asarray(solid.Q_vol, dtype=float).copy()
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
        if 'Q_vol_value' in state and hasattr(solid, 'Q_vol'):
            solid.Q_vol[:] = state['Q_vol_value']
        if 'use_external_source_buffer' in state:
            solid.use_external_source_buffer = state['use_external_source_buffer']
        if 'source_callback' in state:
            solid.source_callback = state['source_callback']

    @staticmethod
    def _capture_coupler_source_state(coupler):
        state = {}
        for attr in ('Q_source_1', 'Q_source_2', '_current_Q_source'):
            if hasattr(coupler, attr):
                state[attr] = np.asarray(getattr(coupler, attr), dtype=float).copy()
        return state

    @staticmethod
    def _restore_coupler_source_state(coupler, state):
        for attr, value in state.items():
            getattr(coupler, attr)[:] = value

    def save_step_state(self):
        electric_attrs = [
            'emitter_voltage', 'emitter_resistivity', 'collector_voltage', 'collector_resistivity',
            'current_density', 'emitter_joule_heat', 'collector_joule_heat'
        ]
        plasma_attrs = [
            'emitter_work_function', 'collector_work_function', 'barrier_voltage_drop',
            'emitter_temperature', 'electron_cooling_flux', 'electron_heating_flux'
        ]

        return {
            'neutronic_total_power': self.neutronic_data.total_power,
            'neutronic_total_power_old': self.neutronic_data._total_power_old,
            'electric': {
                attr: getattr(self.electric_data, attr).copy()
                for attr in electric_attrs
            },
            'plasma': {
                attr: getattr(self.plasma_data, attr).copy()
                for attr in plasma_attrs
            },
            'moderator_temperature': self.boundary_data.moderator_temperature.copy(),
            'solid_sources': [
                (solid, self._capture_solid_source_state(solid))
                for solid in self.solids.values()
            ],
            'coupler_sources': [
                (coupler, self._capture_coupler_source_state(coupler))
                for coupler in self.couplers.values()
            ],
        }

    def load_step_state(self, state):
        if state is None:
            return

        self.neutronic_data.total_power = state['neutronic_total_power']
        self.neutronic_data._total_power_old = state['neutronic_total_power_old']

        for attr, value in state['electric'].items():
            target = getattr(self.electric_data, attr)
            value = np.asarray(value, dtype=float)
            if target.shape == value.shape:
                target[:] = value
            else:
                setattr(self.electric_data, attr, value.copy())
        for attr, value in state['plasma'].items():
            getattr(self.plasma_data, attr)[:] = value

        self.boundary_data.moderator_temperature[:] = state['moderator_temperature']

        for solid, source_state in state['solid_sources']:
            self._restore_solid_source_state(solid, source_state)

        for coupler, source_state in state['coupler_sources']:
            self._restore_coupler_source_state(coupler, source_state)

    def get_solids(self) -> list:
        """
        [接口] 返回组件内部包含的所有底层固体导热求解器对象。
        SystemManager 会在初始化时收集它们以构建全局计算矩阵。
        """
        # 将字典中的值转化为列表返回
        return list(self.solids.values())

    def get_couplers(self) -> list:
        """
        [接口] 返回组件内部包含的所有耦合器。
        SystemManager 会遍历它们来施加网格边界条件和跨域源项。
        """
        # 将字典中的值转化为列表返回
        return list(self.couplers.values())

    # ==========================================
    # 断点续算与状态持久化 (Restart & Persistence)
    # ==========================================

    def get_state_dict(self, prefix: str) -> dict:
        """
        [持久化接口] 将 TFE 宏观组件的缓存状态序列化为扁平字典

        注意：底层的固体网格（如 pellet、emitter）已经由 SystemManager
        通过 get_solids() 扁平化提取并单独保存。这里只保存 TFE 独有的
        电学、等离子体和中子学宏观缓存数据。
        """
        state = {
            # 1. 维度校验指纹 (确保加载时的网格轴向节点数完全一致)
            f"{prefix}/n_axial": np.array([self.mesh.n_axial]),

            # 2. 中子学数据 (标量转为 1D 数组存储)
            f"{prefix}/neutronic/total_power": np.array([self.neutronic_data.total_power]),
            f"{prefix}/neutronic/_total_power_old": np.array([self.neutronic_data._total_power_old]),
        }

        # 3. 电学数据 (1D 数组批量提取)
        electric_attrs = [
            'emitter_voltage', 'emitter_resistivity', 'collector_voltage', 'collector_resistivity',
            'current_density', 'emitter_joule_heat', 'collector_joule_heat'
        ]
        for attr in electric_attrs:
            state[f"{prefix}/electric/{attr}"] = getattr(self.electric_data, attr)

        # 4. 等离子体数据 (1D 数组批量提取)
        plasma_attrs = [
            'emitter_work_function', 'collector_work_function', 'barrier_voltage_drop',
            'emitter_temperature', 'electron_cooling_flux', 'electron_heating_flux'
        ]
        for attr in plasma_attrs:
            state[f"{prefix}/plasma/{attr}"] = getattr(self.plasma_data, attr)

        # 5. 边界条件缓存
        state[f"{prefix}/boundary/moderator_temperature"] = self.boundary_data.moderator_temperature

        return state

    def load_state_dict(self, data: dict, prefix: str):
        """
        [持久化接口] 从扁平字典中恢复 TFE 宏观组件的缓存状态
        """
        # --- 1. 致命级维度校验 ---
        shape_key = f"{prefix}/n_axial"
        if shape_key not in data:
            raise KeyError(f"Missing fingerprint '{shape_key}' for TFEUnit '{self.name}'.")

        saved_n_axial = int(data[shape_key][0])
        if saved_n_axial != self.mesh.n_axial:
            raise ValueError(f"CRITICAL: Axial nodes mismatch for TFEUnit '{self.name}'. "
                             f"File has {saved_n_axial}, current memory has {self.mesh.n_axial}.")

        # --- 2. 恢复中子学数据 ---
        if f"{prefix}/neutronic/total_power" in data:
            self.neutronic_data.total_power = float(data[f"{prefix}/neutronic/total_power"][0])
        if f"{prefix}/neutronic/_total_power_old" in data:
            self.neutronic_data._total_power_old = float(data[f"{prefix}/neutronic/_total_power_old"][0])

        # --- 3. 恢复电学数据 (使用 [:] 原地更新以保护内存地址) ---
        electric_attrs = [
            'emitter_voltage', 'emitter_resistivity', 'collector_voltage', 'collector_resistivity',
            'current_density', 'emitter_joule_heat', 'collector_joule_heat'
        ]
        for attr in electric_attrs:
            key = f"{prefix}/electric/{attr}"
            if key in data:
                target = getattr(self.electric_data, attr)
                value = np.asarray(data[key], dtype=float)
                if target.shape == value.shape:
                    target[:] = value
                else:
                    setattr(self.electric_data, attr, value.copy())

        # --- 4. 恢复等离子体数据 ---
        plasma_attrs = [
            'emitter_work_function', 'collector_work_function', 'barrier_voltage_drop',
            'emitter_temperature', 'electron_cooling_flux', 'electron_heating_flux'
        ]
        for attr in plasma_attrs:
            key = f"{prefix}/plasma/{attr}"
            if key in data:
                getattr(self.plasma_data, attr)[:] = data[key]

        # --- 5. 恢复边界条件缓存 ---
        mod_temp_key = f"{prefix}/boundary/moderator_temperature"
        if mod_temp_key in data:
            self.boundary_data.moderator_temperature[:] = data[mod_temp_key]
