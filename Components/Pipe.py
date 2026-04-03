from Components.BaseComponent import BaseComponent
from Solvers.Couplers import FluidSolidCouple


class Pipe(BaseComponent):
    """
    TASTIN 宏观组件 (Layer 3) - 管道 (Pipe)

    本质：一个包含固体管壁和内部流体的耦合装配体。
    作用：自动建立固体内壁与流体之间的对流换热耦合关系，并向系统管理器暴露底层打工仔。
    """

    def __init__(self,
                 name: str,
                 solid_wall,
                 fluid_channel,
                 heated_perimeter: float,
                 correlation_func,
                 coupled_boundary_name: str = 'left'):
        """
        初始化管道组件

        :param name: 组件名称
        :param solid_wall: 固体管壁对象 (如 HeatConduction2D)
        :param fluid_channel: 流体支路对象 (如 FluidChannel，或者 Volumes 列表，具体取决于您的 Coupler 定义)
        :param heated_perimeter: 加热周长 [m]，用于计算对流换热面积
        :param correlation_func: 对流换热关联式函数 (如 Lyon-Martinelli)
        :param coupled_boundary_name: 固体上与流体接触的边界名称 (默认为内壁 'left')
        """
        super().__init__(name)

        # 1. 挂载物理实体
        self.solid = solid_wall
        self.fluid_channel = fluid_channel

        # 2. 强制初始化固体状态（防呆设计）
        # 确保固体内部网格和物性已经就位，以便提取出准确的热容数据
        self.solid.initialize_state()

        # 3. 提取固体接触面的节点热容 (用于计算自适应步长或耦合稳定性)
        cap = self.solid.get_boundary_node_capacitance(coupled_boundary_name)

        # 4. 核心组装：实例化专属的流固耦合器
        self.coupler = FluidSolidCouple(
            name=f"{name}_coupler",
            fluid=self.fluid_channel,
            solid_boundary_region=self.solid.boundaries[coupled_boundary_name],
            heated_perimeter=heated_perimeter,
            correlation_func=correlation_func,
            solid_node_capacitance=cap
        )

    def get_solids(self) -> list:
        """[接口实现] 交出管壁"""
        return [self.solid]

    def get_couplers(self) -> list:
        """[接口实现] 交出刚组装好的流固耦合器"""
        return [self.coupler]

    # pre_step 和 post_step 目前不需要额外逻辑，继承基类的 pass 即可
