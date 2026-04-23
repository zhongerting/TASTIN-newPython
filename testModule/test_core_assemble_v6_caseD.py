import os
import sys
from typing import Any, Dict, Optional

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from profiler import TEASAProfiler
from test_core_assemble_v6 import run_test_v6_case_a as _run_case_core


def run_test_v6_case_d(
    run_duration_s: float = 20.0,
    total_power_w: float = 115000.0,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: float = 0.0351,
    restart_file: Optional[str] = "test_core_assemble_v5_restart_t10000.npz",
    save_interval: float = 0.0,
    enable_plot: bool = False,
    max_dt: float = 1.0,
    safety_factor: float = 20.0,
    enable_tec_coupled: bool = True,
    tec_update_interval_s: float = 1.0,
    tec_target_voltage_v: float = 27.2,
    tec_initial_current_a: float = 220.0,
) -> Dict[str, Any]:
    """
    Case D:
    1. 从 v5 的长时间稳态场重启动
    2. 标定当前温度场为反馈基准点
    3. 开启 TEC，保持 rho_control = 0
    4. 观察点堆 + 热工 + 电计算三场联动是否稳定

    当前默认值:
    - restart_file = test_core_assemble_v5_restart_t10000.npz
    - run_duration_s = 20.0 s
    - enable_tec_coupled = True
    - tec_update_interval_s = 1.0 s
    - tec_target_voltage_v = 27.2 V
    - tec_initial_current_a = 220.0 A
    """
    return _run_case_core(
        run_duration_s=run_duration_s,
        total_power_w=total_power_w,
        inlet_temperature_k=inlet_temperature_k,
        channel_inlet_flow_kg_s=channel_inlet_flow_kg_s,
        restart_file=restart_file,
        save_interval=save_interval,
        enable_plot=enable_plot,
        max_dt=max_dt,
        safety_factor=safety_factor,
        enable_tec_coupled=enable_tec_coupled,
        tec_update_interval_s=tec_update_interval_s,
        tec_target_voltage_v=tec_target_voltage_v,
        tec_initial_current_a=tec_initial_current_a,
        case_name="D",
        rho_step=0.0,
        step_time_s=5.0,
    )


if __name__ == "__main__":
    run_test_v6_case_d()
    TEASAProfiler.report()
