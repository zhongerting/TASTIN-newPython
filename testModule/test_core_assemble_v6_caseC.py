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


def run_test_v6_case_c(
    run_duration_s: float = 50.0,
    total_power_w: float = 115000.0,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: float = 0.0351,
    restart_file: Optional[str] = "test_core_assemble_v5_restart_t5000.npz",
    save_interval: float = 0.0,
    enable_plot: bool = False,
    max_dt: float = 1.0,
    safety_factor: float = 20.0,
    enable_tec_coupled: bool = False,
    tec_update_interval_s: float = 1.0,
    tec_target_voltage_v: float = 27.2,
    tec_initial_current_a: float = 220.0,
    rho_step: float = -5.0e-5,
    step_time_s: float = 5.0,
) -> Dict[str, Any]:
    """
    Case C:
    1. 从 v5 热工稳态场重启动
    2. 用当前温度场标定反馈基准点
    3. 在 step_time_s 时刻施加一个很小的负阶跃反应性 rho_step

    当前默认值:
    - rho_step = -5.0e-5
    - step_time_s = 5.0 s
    - run_duration_s = 50.0 s
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
        case_name="C",
        rho_step=rho_step,
        step_time_s=step_time_s,
    )


if __name__ == "__main__":
    run_test_v6_case_c()
    TEASAProfiler.report()
