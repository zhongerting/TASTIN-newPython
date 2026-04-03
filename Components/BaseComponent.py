from abc import ABC


class BaseComponent(ABC):
    """
    TASTIN 宏观组件基类 (Layer 3)

    所有复杂的组装体（如 Pipe, TFE_Unit）都必须继承此基类。
    它的核心作用是向 SystemManager 暴露底层物理实体。
    """

    def __init__(self, name: str):
        self.name = name

    def get_solids(self) -> list:
        """
        [接口] 返回组件内部包含的所有底层固体导热求解器对象。
        默认返回空列表，子类需重写。
        """
        return []

    def get_couplers(self) -> list:
        """
        [接口] 返回组件内部包含的所有耦合器 (如 SolidFluidCoupler, SolidSolidCoupler)。
        默认返回空列表，子类需重写。
        """
        return []

    def pre_step(self, dt: float, current_time: float):
        """
        [生命周期钩子] 在每个时间步开始，SystemManager 求解方程之前调用。
        用于处理组件特有的预处理逻辑（例如更新轴向功率分布等）。
        """
        pass

    def post_step(self, dt: float, current_time: float):
        """
        [生命周期钩子] 在每个时间步结束时调用。
        用于组件内部状态的收尾或数据记录。
        """
        pass
