import numpy as np
import matplotlib.pyplot as plt
from Solvers.Neutronics.PointReactor import PointReactor  # 确保你的 PointReactor.py 与此文件在同级目录


def run_test():
    reactor = PointReactor()

    # 1. 设定初始状态
    P0 = 115000.0  # 115 kW
    reactor.initialize_steady_state(P0)

    # 2. 瞬态模拟设置
    t_end = 20.0  # 模拟总时长 50 秒
    dt = 0.1  # 步长 0.1 秒
    times = np.arange(0, t_end + dt, dt)

    # 用于记录数据的列表
    powers = []
    fission_powers = []
    decay_powers = []
    reactivities = []

    print(f"开始测试计算，初始总功率: {reactor.total_power / 1000.0:.2f} kW")

    for t in times:
        # 记录当前收敛状态（本步的初始状态）
        powers.append(reactor.total_power)
        fission_powers.append(reactor.fission_power)
        decay_powers.append(reactor.decay_power)

        # 3. 施加边界条件
        r_ctrl = 0.007  # 恒定的控制反应性
        r_fb = -0.007 * (t / t_end)  # 线性降低的反馈反应性

        reactivities.append(r_ctrl + r_fb)

        # 4. 步进到下一个时间点
        if t < t_end:
            reactor.step(dt, r_ctrl, r_fb)
            reactor.commit()  # 验证通过，直接固化状态

    print(f"计算完成。功率峰值: {max(powers) / 1e6:.2f} MW")

    # 5. 绘制验证图表
    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(times, np.array(powers) / 1000.0, 'b-', linewidth=2, label='Total Power')
    ax1.plot(times, np.array(fission_powers) / 1000.0, 'b--', alpha=0.7, label='Fission Power')
    ax1.plot(times, np.array(decay_powers) / 1000.0, 'b:', alpha=0.7, label='Decay Power')

    ax1.set_xlabel('Time (s)', fontsize=12)
    ax1.set_ylabel('Power (kW)', color='b', fontsize=12)
    ax1.tick_params('y', colors='b')
    ax1.set_yscale('log')  # 对于这种大范围瞬态，使用对数坐标看得更清晰
    ax1.legend(loc='upper left')
    ax1.grid(True, linestyle='--', alpha=0.6)

    ax2 = ax1.twinx()
    ax2.plot(times, reactivities, 'r-', linewidth=2, label='Total Reactivity')
    ax2.set_ylabel('Reactivity', color='r', fontsize=12)
    ax2.tick_params('y', colors='r')
    ax2.legend(loc='upper right')

    plt.title('Point Reactor Kinetics Verification')
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_test()
