import numpy as np
from scipy.integrate import solve_ivp
from typing import Callable, Optional, Tuple, Dict, Any
import logging

# 配置日志
logger = logging.getLogger(__name__)


class NuclearODESolver:
    """
    用于核反应堆物理与热工水力耦合计算的常微分方程组(ODE)求解器。

    该类封装了 scipy.integrate.solve_ivp，专门针对核工程中的刚性方程组
    (Stiff Systems) 进行了预设配置。它旨在替代传统的 GEAR.FOR 算法。
    """

    def __init__(self, method: str = 'BDF', rtol: float = 1e-6, atol: float = 1e-8):
        """
        初始化求解器。

        Parameters
        ----------
        method : str, optional
            求解算法。推荐用于刚性系统的算法：
            - 'BDF': 后向差分公式 (默认)，稳定性好，适合大规模刚性方程。
            - 'Radau': Radau IIA 族隐式方法，精度高，适合高刚性。
            - 'LSODA': 自动在 Adams (非刚性) 和 BDF (刚性) 间切换 (类似原始 ODEPACK)。
        rtol : float
            相对误差容限 (Relative Tolerance)。
        atol : float
            绝对误差容限 (Absolute Tolerance)。
        """
        self.method = method
        self.rtol = rtol
        self.atol = atol
        logger.info(f"Solver initialized with method={method}, rtol={rtol}, atol={atol}")

    def solve(self,
              fun: Callable[[float, np.ndarray], np.ndarray],
              t_span: Tuple[float, float],
              y0: np.ndarray,
              t_eval: Optional[np.ndarray] = None,
              args: Tuple = (),
              jac: Optional[Callable] = None,
              max_step: float = np.inf,  # 建议默认改为 np.inf，以免无意中限制了步长，原代码是 1.0
              first_step: Optional[float] = None,  # [新增] 支持强制设定第一步长
              **kwargs  # [新增] 支持透传其他参数
              ) -> Dict[str, Any]:
        """
        求解 ODE 方程组 dy/dt = f(t, y)。

        Parameters
        ----------
        ... (原有注释保持不变) ...
        first_step : float, optional
            强制求解器使用的第一个时间步长。对于初始时刻有剧烈突变的刚性问题非常重要。
        **kwargs :
            传递给 scipy.integrate.solve_ivp 的其他关键字参数。
            :param first_step:
            :param max_step:
            :param jac:
            :param args:
            :param t_eval:
            :param y0:
            :param t_span:
            :param fun:
        """

        logger.debug(f"Starting integration from {t_span[0]} to {t_span[1]}")

        try:
            # 调用 SciPy 的核心求解器
            sol = solve_ivp(
                fun=fun,
                t_span=t_span,
                y0=y0,
                method=self.method,
                t_eval=t_eval,
                args=args,
                jac=jac,
                rtol=self.rtol,
                atol=self.atol,
                max_step=max_step,
                first_step=first_step,  # [新增] 传递第一步长
                vectorized=False,
                **kwargs  # [新增] 传递其他参数
            )

            if sol.success:
                logger.info("Integration successful.")
                return {
                    'success': True,
                    't': sol.t,
                    'y': sol.y,
                    'message': sol.message
                }
            else:
                logger.warning(f"Integration failed: {sol.message}")
                return {
                    'success': False,
                    't': sol.t,
                    'y': sol.y,
                    'message': sol.message
                }

        except Exception as e:
            logger.error(f"Solver encountered a critical error: {str(e)}")
            raise e


# ==========================================
# 使用示例 (对应架构设计中的验证阶段)
# ==========================================

if __name__ == "__main__":
    # 配置简单的日志输出
    logging.basicConfig(level=logging.INFO)

    # 1. 定义物理模型 (模拟 Fortran 中的 FCN 子程序)
    # 这是一个简化的点堆动力学模型示例 (Point Kinetics)
    def point_kinetics_rhs(t, y, rho_ext, beta, Lambda):
        """
        y[0]: 中子密度 n
        y[1]: 先驱核浓度 c (简化为1组)
        """
        n = y[0]
        c = y[1]

        # 反应性 rho(t) 可以是时间的函数
        rho = rho_ext

        dn_dt = (rho - beta) / Lambda * n + 0.1 * c  # 假设 lambda_precursor = 0.1
        dc_dt = beta / Lambda * n - 0.1 * c

        return np.array([dn_dt, dc_dt])


    # 2. 准备初始条件和参数
    y_initial = np.array([1.0, 10.0])  # 初始 n, c
    t_start, t_end = 0.0, 100.0
    parameters = (0.002, 0.0065, 0.0001)  # rho, beta, Lambda

    step_size = 0.05
    eval_points = np.arange(t_start, t_end, step_size)

    # 3. 实例化求解器 (使用 BDF 处理刚性)
    solver = NuclearODESolver(method='BDF', rtol=1e-5, atol=1e-8)

    # 4. 执行求解
    result = solver.solve(fun=point_kinetics_rhs, t_span=(t_start, t_end), y0=y_initial, t_eval=eval_points,
                          args=parameters, jac=None)

    print(result['y'])

    # 5. 输出结果
    if result['success']:
        print(f"Final Neutron Density: {result['y'][0, -1]:.4f}")
    else:
        print("Calculation failed.")
