import numpy as np
import sys


# ==============================================================================
# 1. 模拟环境 (Mock Objects)
#    用于隔离测试，不依赖 TASTIN 庞大的系统库
# ==============================================================================

class MockResistanceBC:
    """模拟边界条件对象"""

    def __init__(self):
        self.T_ext = 0.0
        self.R_ext = 0.0

    def update_params(self, T_ext, R_ext):
        self.T_ext = float(T_ext)
        self.R_ext = float(R_ext)


class MockBoundary:
    """模拟边界区域"""

    def __init__(self, name, T_node, R_internal):
        self.name = name
        self.T_node = float(T_node)
        self.R_internal = float(R_internal)
        self.shape = (1,)  # 模拟单节点
        self.bc = MockResistanceBC()

    def get_coupling_snapshot(self):
        """返回 (T_node, R_internal) 快照"""
        return np.array([self.T_node]), np.array([self.R_internal])


class MockGapCoupleBase:
    """
    模拟基类 (替代 GapCouple2D)
    强制设定固定的间隙热阻，以便与手算结果对比
    """

    def __init__(self, obj1_bound, obj2_bound, R_gap_fixed=0.8):
        self.bound1 = obj1_bound
        self.bound2 = obj2_bound
        self.bc1 = obj1_bound.bc
        self.bc2 = obj2_bound.bc
        self.R_gap_total = R_gap_fixed

    def sync(self):
        """
        基类行为：只计算基础的串联热阻，不加源项
        """
        # 获取快照
        T1_node, R1_int = self.bound1.get_coupling_snapshot()
        T2_node, R2_int = self.bound2.get_coupling_snapshot()

        # 基础计算 (无源)
        # Obj1 看到: R_gap + R2_int, 对侧温度 T2_node
        self.bc1.update_params(T_ext=T2_node, R_ext=self.R_gap_total + R2_int)

        # Obj2 看到: R_gap + R1_int, 对侧温度 T1_node
        self.bc2.update_params(T_ext=T1_node, R_ext=self.R_gap_total + R1_int)


# ==============================================================================
# 2. 待测类定义 (TECCouple2D)
#    (将之前生成的代码逻辑嵌入此处，继承自 MockGapCoupleBase 进行测试)
# ==============================================================================

class TECCouple2D(MockGapCoupleBase):
    """
    待验证的 TEC 耦合器逻辑
    """

    def __init__(self, obj1_bound, obj2_bound, R_gap_fixed=0.8):
        super().__init__(obj1_bound, obj2_bound, R_gap_fixed)
        self.shape = (1,)
        self.Q_source_1 = np.zeros(self.shape)
        self.Q_source_2 = np.zeros(self.shape)

    def set_tec_sources(self, Q_emitter, Q_collector):
        self.Q_source_1[:] = Q_emitter
        self.Q_source_2[:] = Q_collector

    def sync(self):
        # 1. 执行基类基础计算
        super().sync()

        # 2. 获取基类计算好的参数
        R_ext_1 = self.bc1.R_ext  # R_gap + R_int2
        R_ext_2 = self.bc2.R_ext  # R_gap + R_int1

        _, R_int_1 = self.bound1.get_coupling_snapshot()
        _, R_int_2 = self.bound2.get_coupling_snapshot()

        T_node_2 = self.bc1.T_ext
        T_node_1 = self.bc2.T_ext

        # 3. 戴维南等效修正 (核心验证逻辑)
        # T_ext_new = T_base + Q_self * R_ext_self + Q_other * R_int_other
        T_ext_1_new = T_node_2 + (self.Q_source_1 * R_ext_1) + (self.Q_source_2 * R_int_2)
        T_ext_2_new = T_node_1 + (self.Q_source_2 * R_ext_2) + (self.Q_source_1 * R_int_1)

        # 4. 更新边界
        self.bc1.update_params(T_ext=T_ext_1_new, R_ext=R_ext_1)
        self.bc2.update_params(T_ext=T_ext_2_new, R_ext=R_ext_2)


# ==============================================================================
# 3. 验证主程序
# ==============================================================================

def run_verification():
    print("=" * 60)
    print("TECCouple2D 逻辑验证程序 (基于双节点手算模型)")
    print("=" * 60)

    # --- 设定参数 (与手算一致) ---
    T_emit = 2000.0  # Obj1 Node Temp
    T_coll = 1000.0  # Obj2 Node Temp
    R_int1 = 0.1
    R_int2 = 0.1
    R_gap = 0.8
    # R_loop = 0.1 + 0.8 + 0.1 = 1.0

    print(f"参数设定:")
    print(f"  T_emitter (Node)   = {T_emit} K")
    print(f"  T_collector (Node) = {T_coll} K")
    print(f"  R_int1 = {R_int1}, R_int2 = {R_int2}, R_gap = {R_gap} K/W")
    print(f"  理论回路总热阻 R_loop = 1.0 K/W")
    print("-" * 60)

    # --- 初始化对象 ---
    bound1 = MockBoundary("Emitter", T_emit, R_int1)
    bound2 = MockBoundary("Collector", T_coll, R_int2)

    # 实例化耦合器
    coupler = TECCouple2D(bound1, bound2, R_gap_fixed=R_gap)

    # ==========================================================================
    # 工况 1: 纯电子冷却 (Emitter Cooling Only)
    # Q1 = -100 W, Q2 = 0
    # ==========================================================================
    print("\n【工况 1: 纯电子冷却测试】")
    Q1_val = -100.0
    Q2_val = 0.0
    print(f"  输入: Q_emit = {Q1_val} W, Q_coll = {Q2_val} W")

    # 1. 设置源项并同步
    coupler.set_tec_sources(Q1_val, Q2_val)
    coupler.sync()

    # 2. 获取计算出的 T_ext
    T_ext_1_calc = bound1.bc.T_ext

    # 3. 计算求解器将会看到的热流 (Flux = (T_node - T_ext) / (R_int + R_ext))
    #    注意：R_ext 在 sync 中已被设定为 R_gap + R_int2 = 0.9
    #    分母总和 = R_int1 + R_ext = 0.1 + 0.9 = 1.0
    Flux_1_solver = (T_emit - T_ext_1_calc) / (R_int1 + bound1.bc.R_ext)

    # 4. 手算真值
    # 基础热流 = (2000-1000)/1.0 = 1000 W
    # 源项贡献 = 100 W 全部抽出。
    # 发射极承担流出增量 = 100 * (R_ext / R_loop) = 100 * 0.9 = 90 W
    # 总流出 = 1090 W
    Flux_1_truth = 1090.0
    T_ext_1_truth = 910.0  # 2000 - 1090*1.0 + 0.1*1090? No, T_ext = T_node - Flux * (R_tot) + ?
    # T_ext_truth check: Flux = (2000 - T_ext)/1.0 => 1090 = 2000 - T_ext => T_ext = 910. Correct.

    print(f"  结果对比 (Emitter Side):")
    print(f"    T_ext (Solver) : {T_ext_1_calc:.4f} K  [预期: {T_ext_1_truth:.4f} K]")
    print(f"    Flux  (Solver) : {Flux_1_solver:.4f} W  [预期: {Flux_1_truth:.4f} W]")

    if abs(Flux_1_solver - Flux_1_truth) < 1e-5:
        print("  >> 状态: PASS ✅")
    else:
        print("  >> 状态: FAILED ❌")

    # ==========================================================================
    # 工况 2: 混合模式 (Cooling + Heating)
    # Q1 = -100 W, Q2 = +50 W
    # ==========================================================================
    print("\n【工况 2: 混合模式测试】")
    Q1_val = -100.0
    Q2_val = 50.0
    print(f"  输入: Q_emit = {Q1_val} W, Q_coll = {Q2_val} W")

    coupler.set_tec_sources(Q1_val, Q2_val)
    coupler.sync()

    T_ext_1_calc = bound1.bc.T_ext
    Flux_1_solver = (T_emit - T_ext_1_calc) / (R_int1 + bound1.bc.R_ext)

    # 手算真值
    # 基础: 1000
    # Q1贡献: +90 (流出)
    # Q2贡献: -5 (流出减少，因为对面加热阻碍了散热) -> 50 * (0.1 / 1.0) = 5
    # 总流出: 1000 + 90 - 5 = 1085 W
    Flux_1_truth = 1085.0
    T_ext_1_truth = 915.0  # 2000 - 1085 = 915. Correct.

    print(f"  结果对比 (Emitter Side):")
    print(f"    T_ext (Solver) : {T_ext_1_calc:.4f} K  [预期: {T_ext_1_truth:.4f} K]")
    print(f"    Flux  (Solver) : {Flux_1_solver:.4f} W  [预期: {Flux_1_truth:.4f} W]")

    if abs(Flux_1_solver - Flux_1_truth) < 1e-5:
        print("  >> 状态: PASS ✅")
    else:
        print("  >> 状态: FAILED ❌")


if __name__ == "__main__":
    run_verification()
