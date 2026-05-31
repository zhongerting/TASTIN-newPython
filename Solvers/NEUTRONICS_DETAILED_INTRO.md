# Neutronics 中子动力学模块详细介绍

## 1. 模块定位

`Neutronics/` 是 TASTIN `Solvers` 层中的中子动力学求解子系统。当前目录中只有一个核心文件：

```text
Neutronics/
  PointReactor.py    # 点堆中子动力学 + 衰变热求解器
```

该模块实现了一个 11 维点堆动力学模型，用于计算反应堆总功率中的裂变功率和衰变热功率，并为 `SystemManager`、`ReactorCore` 等上层组件提供瞬态功率源。

在全局多物理场流程中，Neutronics 的典型职责是：

- 根据控制棒/外部控制给定的反应性计算裂变功率变化；
- 根据温度反馈或组件反馈修正总反应性；
- 跟踪 6 组缓发中子先驱核；
- 跟踪 4 组衰变热能量；
- 输出裂变功率、衰变热功率和总功率；
- 在 Picard 内迭代中使用双缓冲状态，支持试探积分、提交和回滚。

## 2. 物理模型

当前 `PointReactor` 使用经典点堆动力学方程，状态向量长度为 11：

```text
y[0]     = P_fiss     裂变功率
y[1:7]   = C_i        6 组缓发中子先驱核等效功率变量
y[7:11]  = W_j        4 组衰变热等效能量变量
```

### 2.1 裂变功率方程

裂变功率满足：

```text
dP_fiss/dt =
  ((rho - beta_total) / Lambda) * P_fiss
  + sum(lambda_i * C_i)
```

其中：

- `rho` 是总反应性；
- `beta_total` 是总缓发中子份额；
- `Lambda` 是中子代时间；
- `C_i` 是第 `i` 组缓发中子先驱核变量；
- `lambda_i` 是第 `i` 组先驱核衰变常数。

### 2.2 缓发中子方程

6 组缓发中子先驱核方程为：

```text
dC_i/dt = (beta_i / Lambda) * P_fiss - lambda_i * C_i
```

在稳态下：

```text
C_i = (beta_i / Lambda / lambda_i) * P_fiss
```

### 2.3 衰变热方程

4 组衰变热变量满足：

```text
dW_j/dt = gamma_j * P_fiss - lambda_dj * W_j
```

衰变热功率按原 C++ 代码逻辑定义为：

```text
P_decay = sum(W_j * lambda_dj)
```

总功率为：

```text
P_total = P_fiss + P_decay
```

## 3. 物理参数

`PointReactor.__init__()` 中内置了所有点堆参数。

### 3.1 中子代时间

```python
Lambda = 0.2e-4
```

即：

```text
Lambda = 2.0e-5 s
```

### 3.2 缓发中子参数

缓发中子采用 6 组模型：

```python
beta_fra = [
    2.618e-4,
    1.737e-3,
    1.555e-3,
    3.133e-3,
    9.123e-4,
    3.330e-4,
]

lambda_c = [
    1.271e-2,
    3.174e-2,
    1.160e-1,
    3.110e-2,
    1.397e+0,
    3.872e+0,
]
```

总缓发中子份额：

```text
beta_total = sum(beta_fra) = 0.0079321
```

### 3.3 衰变热参数

衰变热采用 4 组模型：

```python
gamma_fra = [0.01728, 0.01365, 0.01024, 0.02304]
lambda_d  = [1.00e-5, 3.00e-3, 1.00e-2, 4.00e-2]
```

总衰变热份额：

```text
gamma_total = sum(gamma_fra) = 0.06421
```

这意味着在稳态下：

```text
P_decay = gamma_total * P_fiss
P_total = (1 + gamma_total) * P_fiss
```

## 4. 数值实现

### 4.1 RHS 函数 `_prke_rhs`

`_prke_rhs()` 是点堆动力学 ODE 的右端项函数：

```python
_prke_rhs(t, y, rho, Lambda, beta_total,
          beta_fra, lambda_c, gamma_fra, lambda_d)
```

它返回 11 维 `dydt`：

```text
dydt[0]      裂变功率导数
dydt[1:7]    缓发中子先驱核导数
dydt[7:11]   衰变热变量导数
```

该函数使用 `@njit(cache=True)` 编译，避免 SciPy ODE 求解器高频调用 RHS 时产生 Python 循环开销。

若环境中没有安装 Numba，源码会提供纯 Python fallback：

```python
try:
    from numba import njit
except ImportError:
    def njit(*args, **kwargs):
        ...
```

因此模块可以在没有 Numba 的环境中运行，只是性能会下降。

### 4.2 解析雅可比 `_prke_jac`

`_prke_jac()` 返回 11x11 解析雅可比矩阵，用于加速隐式 ODE 求解器：

```python
_prke_jac(t, y, rho, Lambda, beta_total,
          beta_fra, lambda_c, gamma_fra, lambda_d)
```

主要非零项包括：

```text
d(dP/dt)/dP      = (rho - beta_total) / Lambda
d(dP/dt)/dC_i    = lambda_i
d(dC_i/dt)/dP    = beta_i / Lambda
d(dC_i/dt)/dC_i  = -lambda_i
d(dW_j/dt)/dP    = gamma_j
d(dW_j/dt)/dW_j  = -lambda_dj
```

点堆动力学方程通常具有较强刚性，尤其是中子代时间很小而衰变热时间尺度很长。解析雅可比能显著提升 `Radau` 或 `BDF` 的稳定性和效率。

### 4.3 积分器配置

`PointReactor` 默认使用：

```python
integrator_method = "Radau"
rtol = 1e-6
atol = 1e-8
```

`Radau` 是隐式 Runge-Kutta 方法，适合刚性系统。源码注释也说明可选 `BDF`。

每次 `step()` 内部调用：

```python
solve_ivp(
    fun=_prke_rhs,
    t_span=(0.0, dt),
    y0=self._y_committed,
    method=self.integrator_method,
    jac=_prke_jac,
    args=args,
    rtol=self.rtol,
    atol=self.atol,
)
```

注意：这里的积分区间总是从 `0.0` 到 `dt`，因为点堆方程只依赖当前步内给定的总反应性，不依赖绝对时间。

## 5. `PointReactor` 类

### 5.1 状态双缓冲设计

`PointReactor` 内部维护两个状态：

```python
self._y_committed
self._y_trial
```

含义：

| 状态 | 含义 |
| --- | --- |
| `_y_committed` | 上一个已经收敛并提交的真实时间步状态 |
| `_y_trial` | 当前 Picard 内迭代中的试探积分状态 |

`step()` 总是从 `_y_committed` 重新积分到 `dt` 末端，并将结果写入 `_y_trial`。

只有当全局 Picard 迭代收敛后，`SystemManager` 才调用：

```python
point_reactor.commit()
```

把 `_y_trial` 固化为新的 `_y_committed`。

这种设计的优点是：如果当前 Picard 迭代的热工反馈或耦合结果未收敛，可以反复从同一个已提交起点重新积分，而不会污染真实历史状态。

### 5.2 初始化稳态

`initialize_steady_state(total_power_initial)` 根据给定总功率解析初始化稳态。

由于稳态衰变热功率满足：

```text
P_decay = gamma_total * P_fiss
```

所以：

```text
P_fiss_0 = P_total_0 / (1 + gamma_total)
```

随后按稳态条件初始化先驱核和衰变热变量：

```text
C_i = (beta_i / Lambda / lambda_i) * P_fiss_0
W_j = (gamma_j / lambda_dj) * P_fiss_0
```

初始化完成后，`_y_trial` 与 `_y_committed` 保持一致。

### 5.3 时间步推进

`step(dt, reactivity_control, reactivity_feedback)` 接收两个反应性输入：

```text
total_reactivity = reactivity_control + reactivity_feedback
```

其中：

- `reactivity_control` 通常表示外部控制输入，例如控制棒、设定扰动；
- `reactivity_feedback` 通常来自热工反馈，例如燃料温度、冷却剂温度或组件层计算。

积分成功后：

```python
self._y_trial = sol.y[:, -1]
return True
```

若 SciPy 积分失败，会抛出：

```python
RuntimeError("PointReactor ODE integration failed: ...")
```

### 5.4 提交状态

`commit()` 在全局时间步收敛后调用：

```python
self._y_committed = self._y_trial.copy()
```

这一步代表中子动力学状态正式进入下一时间步。

### 5.5 步内保存与恢复

`save_step_state()` 和 `load_step_state(state)` 用于保存和恢复当前点堆对象状态：

```python
{
    "y_committed": ...,
    "y_trial": ...,
}
```

这通常服务于 `SystemManager` 的异常回滚或 Picard 迭代回滚机制。

## 6. 输出接口

`PointReactor` 通过属性接口向其他模块输出功率：

### 6.1 裂变功率

```python
point_reactor.fission_power
```

返回：

```text
_y_trial[0]
```

### 6.2 衰变热功率

```python
point_reactor.decay_power
```

计算：

```text
sum(_y_trial[7 + j] * lambda_d[j])
```

### 6.3 总功率

```python
point_reactor.total_power
```

计算：

```text
fission_power + decay_power
```

### 6.4 当前完整状态

```python
point_reactor.current_state
```

返回当前试探状态 `_y_trial`。

## 7. 与 SystemManager 和组件层的关系

根据 `SOLVERS_ANALYSIS.md` 和 `AI_AGENT_SOLVERS_ANALYSIS.md` 的说明，Neutronics 通常在 `SystemManager.step()` 的 Picard 内迭代中被调用。

典型流程：

```text
SystemManager.step(dt)
  1. 保存入口状态
  2. 组件 pre_step 更新功率分配、TEC 等
  3. Picard 内迭代:
     3.1 运行耦合器，刷新热工边界
     3.2 推进点堆中子动力学:
         point_reactor.step(dt, reactivity_control, reactivity_feedback)
     3.3 将裂变功率和衰变热分配到固体热源
     3.4 推进水力网络
     3.5 推进固体导热
     3.6 检查流体/固体温度收敛
  4. 收敛后 point_reactor.commit()
  5. global_time += dt
```

因此，点堆求解器本身只负责给定反应性下的功率时间演化。温度反馈、功率空间分配、热源映射和全局收敛控制由上层组件和 `SystemManager` 完成。

## 8. 断点续算接口

`PointReactor` 支持扁平字典形式的状态保存：

```python
get_state_dict(prefix)
```

保存内容：

| 键 | 含义 |
| --- | --- |
| `{prefix}/y_committed` | 已提交状态 |
| `{prefix}/y_trial` | 当前试探状态 |
| `{prefix}/Lambda` | 中子代时间 |
| `{prefix}/rtol` | ODE 相对误差 |
| `{prefix}/atol` | ODE 绝对误差 |

恢复接口：

```python
load_state_dict(data, prefix)
```

若只有 `y_committed` 而没有 `y_trial`，恢复时会自动令：

```python
_y_trial = _y_committed.copy()
```

## 9. 典型使用流程

### 9.1 初始化点堆对象

```python
from Solvers.Neutronics.PointReactor import PointReactor

reactor = PointReactor()
reactor.initialize_steady_state(total_power_initial=1.0e6)
```

### 9.2 执行一个试探时间步

```python
reactor.step(
    dt=0.01,
    reactivity_control=0.0,
    reactivity_feedback=-2.0e-4,
)

print(reactor.fission_power)
print(reactor.decay_power)
print(reactor.total_power)
```

### 9.3 全局收敛后提交

```python
reactor.commit()
```

### 9.4 保存和恢复

```python
state = reactor.get_state_dict("neutronics")

reactor2 = PointReactor()
reactor2.load_state_dict(state, "neutronics")
```

## 10. 重要实现细节与注意事项

### 10.1 反应性单位

`reactivity_control` 和 `reactivity_feedback` 在代码中直接作为点堆方程中的 `rho` 使用，应保持无量纲反应性单位。

### 10.2 `step()` 不自动提交

调用 `step()` 只更新 `_y_trial`，不会改变 `_y_committed`。必须在全局时间步确认收敛后调用 `commit()`。

### 10.3 时间步内反应性保持常值

当前 `step()` 接收的是一个标量总反应性，并在 `[0, dt]` 积分期间保持不变。若未来需要时间相关反应性，需要修改 RHS 参数传递方式或在 RHS 中引入时间函数。

### 10.4 功率空间分布不在本模块中处理

`PointReactor` 只输出总裂变功率和总衰变热功率。功率如何分配到燃料、包壳、冷却剂或热电组件，由 `ReactorCore`、组件层或其他源项管理逻辑负责。

### 10.5 Numba 是性能优化，不是硬依赖

没有 Numba 时模块仍可运行，但 RHS 和 Jacobian 将退化为普通 Python 函数。对于频繁调用的瞬态计算，建议安装 Numba。

### 10.6 刚性求解器更适合该模型

点堆动力学同时包含 prompt neutron、delayed neutron 和 decay heat 等不同时间尺度。默认 `Radau` 是合理选择；若切换为显式方法，可能需要非常小的时间步。

## 11. 模块优点

当前 Neutronics 实现具有以下特点：

- 模型紧凑，只包含点堆动力学和衰变热核心；
- 11 维状态清晰，便于调试和持久化；
- 使用解析稳态初始化，避免初始先驱核和衰变热不平衡；
- 使用 Numba 编译 RHS 和雅可比，提高隐式积分性能；
- 使用双缓冲状态，天然适配 Picard 内迭代；
- 输出接口简单，直接提供裂变功率、衰变热功率和总功率；
- 与 `SystemManager` 的回滚和断点续算机制兼容。

## 12. 推荐阅读顺序

维护或扩展 Neutronics 时，建议按以下顺序阅读：

1. `Neutronics/PointReactor.py`
   - 先看 `_prke_rhs()`，理解 11 维方程；
   - 再看 `_prke_jac()`，理解解析雅可比；
   - 最后看 `PointReactor.step()`、`commit()` 和状态接口；
2. `SystemManager.py`
   - 查看点堆求解器在 Picard 内迭代中的调用位置；
3. 组件层中的 ReactorCore 或功率分配逻辑
   - 理解点堆总功率如何映射为固体热源。

## 13. 源码校验信息

最后按源码校验日期：2026-05-31。

涉及文件：

- `Neutronics/PointReactor.py`
- `SystemManager.py`
