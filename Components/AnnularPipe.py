from Components.BaseComponent import BaseComponent
from Solvers.Couplers import FluidSolidCouple


class AnnularPipe(BaseComponent):
    """
    TASTIN 宏观组件 (Layer 3) - 环形管道 (Annular Pipe)

    本质：由内管壁、外管壁以及它们之间的环形流体通道组成的装配体。
    作用：自动建立流体与两侧固体壁面的双向对流换热耦合，向系统暴露 2 个固体和 2 个耦合器。
    """

    def __init__(self,
                 name: str,
                 solid_inner,
                 solid_outer,
                 fluid_channel,
                 heated_perimeter_inner: float,
                 heated_perimeter_outer: float,
                 correlation_func,
                 boundary_inner_solid: str = 'right',
                 boundary_outer_solid: str = 'left'):
        """
        初始化环形管道组件

        :param name: 组件名称
        :param solid_inner: 内部固体管壁对象 (通常受内侧热源加热)
        :param solid_outer: 外部固体管壁对象 (通常向外散热或绝热)
        :param fluid_channel: 环形间隙中的流体支路对象/节点列表
        :param heated_perimeter_inner: 流体与内管壁接触的湿周长 [m]
        :param heated_perimeter_outer: 流体与外管壁接触的湿周长 [m]
        :param correlation_func: 对流换热关联式 (内外壁默认使用同一种流动状态下的关联式)
        :param boundary_inner_solid: 内管壁与流体接触的边界名称 (默认为 'right' 或 'outer')
        :param boundary_outer_solid: 外管壁与流体接触的边界名称 (默认为 'left' 或 'inner')
        """
        super().__init__(name)

        # 1. 挂载物理实体
        self.solid_inner = solid_inner
        self.solid_outer = solid_outer
        self.fluid_channel = fluid_channel

        # 2. 强制初始化固体状态，提取热容数据以备耦合稳定性分析
        self.solid_inner.initialize_state()
        self.solid_outer.initialize_state()

        cap_inner = self.solid_inner.get_boundary_node_capacitance(boundary_inner_solid)
        cap_outer = self.solid_outer.get_boundary_node_capacitance(boundary_outer_solid)

        # 3. 核心组装 A：内管壁-流体 耦合器
        self.coupler_inner = FluidSolidCouple(
            name=f"{name}_coupler_inner",
            fluid=self.fluid_channel,
            solid_boundary_region=self.solid_inner.boundaries[boundary_inner_solid],
            heated_perimeter=heated_perimeter_inner,
            correlation_func=correlation_func,
            solid_node_capacitance=cap_inner
        )

        # 4. 核心组装 B：外管壁-流体 耦合器
        self.coupler_outer = FluidSolidCouple(
            name=f"{name}_coupler_outer",
            fluid=self.fluid_channel,
            solid_boundary_region=self.solid_outer.boundaries[boundary_outer_solid],
            heated_perimeter=heated_perimeter_outer,
            correlation_func=correlation_func,
            solid_node_capacitance=cap_outer
        )

    def get_solids(self) -> list:
        """[接口实现] 向系统上缴两个管壁对象"""
        return [self.solid_inner, self.solid_outer]

    def get_couplers(self) -> list:
        """[接口实现] 向系统上缴两组流固耦合器"""
        return [self.coupler_inner, self.coupler_outer]
