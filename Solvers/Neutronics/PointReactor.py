import numpy as np
from scipy.integrate import solve_ivp
from numba import njit

# =============================================================================
# 底层纯数学/物理核心 (使用 Numba JIT 编译以绕过 Python 循环开销)
# =============================================================================

@njit(cache=True)
def _prke_rhs(t, y, rho, Lambda, beta_total, beta_fra, lambda_c, gamma_fra, lambda_d):
    """
    点堆动力学方程的右端项 (RHS)
    y[0]    : 裂变功率 P_fiss
    y[1:7]  : 缓发中子先驱核等效功率 C_i (6组)
    y[7:11] : 衰变热等效能量 W_j (4组)
    """
    dydt = np.zeros(11)
    P_fiss = y[0]

    # 计算缓发中子先驱核衰变总和
    sum_lambda_C = 0.0
    for i in range(6):
        sum_lambda_C += lambda_c[i] * y[i + 1]
        # 缓发中子方程 dCi/dt = (beta_i / Lambda) * P_fiss - lambda_i * C_i
        dydt[i + 1] = (beta_fra[i] / Lambda) * P_fiss - lambda_c[i] * y[i + 1]

    # 裂变中子方程 dP/dt = ((rho - beta_total) / Lambda) * P_fiss + sum(lambda_i * C_i)
    dydt[0] = (rho - beta_total) / Lambda * P_fiss + sum_lambda_C

    # 衰变热方程 dWj/dt = gamma_j * P_fiss - lambda_dj * W_j
    for j in range(4):
        dydt[7 + j] = gamma_fra[j] * P_fiss - lambda_d[j] * y[7 + j]

    return dydt


@njit(cache=True)
def _prke_jac(t, y, rho, Lambda, beta_total, beta_fra, lambda_c, gamma_fra, lambda_d):
    """
    点堆动力学方程的解析雅可比矩阵 (Jacobian)
    极大地提升隐式求解器 (BDF/Radau) 处理刚性方程的效率
    """
    jac = np.zeros((11, 11))

    # d(dP/dt) / dP
    jac[0, 0] = (rho - beta_total) / Lambda

    for i in range(6):
        # d(dP/dt) / dCi
        jac[0, i + 1] = lambda_c[i]
        # d(dCi/dt) / dP
        jac[i + 1, 0] = beta_fra[i] / Lambda
        # d(dCi/dt) / dCi
        jac[i + 1, i + 1] = -lambda_c[i]

    for j in range(4):
        # d(dWj/dt) / dP
        jac[7 + j, 0] = gamma_fra[j]
        # d(dWj/dt) / dWj
        jac[7 + j, 7 + j] = -lambda_d[j]

    return jac

# =============================================================================
# 点堆动力学管理器 (TASTIN Layer 2 Solver)
# =============================================================================

class PointReactor:
    def __init__(self):
        # --- 物理参数定义 (根据原 C++ 代码设定) ---
        self.Lambda = 0.2e-4  # 中子代时间 neulife

        # 6组缓发中子参数 (betaFra, betaCon)
        self.beta_fra = np.array([2.618e-4, 1.737e-3, 1.555e-3, 3.133e-3, 9.123e-4, 3.330e-4])
        self.lambda_c = np.array([1.271e-2, 3.174e-2, 1.160e-1, 3.110e-2, 1.397e+0, 3.872e+0])
        self.beta_total = np.sum(self.beta_fra)  # = 0.0079321

        # 4组衰变热参数 (decayFra, decayCon)
        self.gamma_fra = np.array([0.01728, 0.01365, 0.01024, 0.02304])
        self.lambda_d = np.array([1.00e-5, 3.00e-3, 1.00e-2, 4.00e-2])
        self.gamma_total = np.sum(self.gamma_fra)  # = 0.06421

        # --- 状态变量 (11维向量) ---
        # _y_committed: 上一个收敛的时间步的真实状态 (用于回滚和新一步的起点)
        # _y_trial: 当前 Picard 迭代中的试探状态
        self._y_committed = np.zeros(11)
        self._y_trial = np.zeros(11)

        # 数值求解器配置
        self.integrator_method = 'Radau'  # 推荐用于高刚性系统，也可选 'BDF'
        self.rtol = 1e-6
        self.atol = 1e-8

    def initialize_steady_state(self, total_power_initial: float):
        """
        根据给定的初始总功率，解析计算并初始化绝对稳态
        """
        # 在稳态下：总功率 = 裂变功率 + 衰变热功率
        # 根据衰变热方程 dWj/dt = 0 -> Wj * lambda_dj = gamma_j * P_fiss
        # 你的 C++ 代码中衰变热功率定义为 P_decay = sum(W_j * lambda_dj)
        # 结合以上两者，稳态衰变热功率 = sum(gamma_j) * P_fiss = gamma_total * P_fiss
        # 因此: P_tot_0 = P_fiss_0 + gamma_total * P_fiss_0
        p_fiss_0 = total_power_initial / (1.0 + self.gamma_total)

        self._y_committed[0] = p_fiss_0

        # 稳态缓发中子 dCi/dt = 0 -> Ci = (beta_i / Lambda / lambda_i) * P_fiss_0
        for i in range(6):
            self._y_committed[i + 1] = (self.beta_fra[i] / self.Lambda / self.lambda_c[i]) * p_fiss_0

        # 稳态衰变热 dWj/dt = 0 -> Wj = (gamma_j / lambda_dj) * p_fiss_0
        for j in range(4):
            self._y_committed[7 + j] = (self.gamma_fra[j] / self.lambda_d[j]) * p_fiss_0

        # 同步试探状态
        self._y_trial = self._y_committed.copy()

    def step(self, dt: float, reactivity_control: float, reactivity_feedback: float):
        """
        执行一个时间步的积分（用于全局 Picard 迭代中）。
        传入本迭代步最新的 feedback，从上一次 commit 的状态重新积分。
        """
        total_reactivity = reactivity_control + reactivity_feedback

        # 封装传递给 SciPy 的 RHS 和 Jac
        args = (total_reactivity, self.Lambda, self.beta_total,
                self.beta_fra, self.lambda_c, self.gamma_fra, self.lambda_d)

        # 求解从 t=0 到 t=dt 的瞬态
        sol = solve_ivp(
            fun=_prke_rhs,
            t_span=(0.0, dt),
            y0=self._y_committed,
            method=self.integrator_method,
            jac=_prke_jac,
            args=args,
            rtol=self.rtol,
            atol=self.atol
        )

        if not sol.success:
            raise RuntimeError(f"PointReactor ODE integration failed: {sol.message}")

        # 提取当前步的最终状态，作为试探状态存入
        self._y_trial = sol.y[:, -1]

    def commit(self):
        """
        当全局 Picard 迭代收敛时，由 SystemManager 调用，固化当前时间步状态。
        """
        self._y_committed = self._y_trial.copy()

    # --- 属性接口：方便其他模块获取功率状态 ---

    @property
    def fission_power(self) -> float:
        return self._y_trial[0]

    @property
    def decay_power(self) -> float:
        """
        按照原 C++ 代码逻辑，衰变热功率 = sum(W_j * lambda_dj)
        """
        d_power = 0.0
        for j in range(4):
            d_power += self._y_trial[7 + j] * self.lambda_d[j]
        return d_power

    @property
    def total_power(self) -> float:
        return self.fission_power + self.decay_power

    @property
    def current_state(self) -> np.ndarray:
        return self._y_trial
