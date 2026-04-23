# TASTIN 源码 Bug 反馈报告

> **日期**：2026-04-15
> **背景**：在搭建 TOPAZ-II 主冷却剂回路仿真（V3）过程中发现的源码 bug。
> 修复方案均为**运行时 monkey-patch**，未改动源码。

---

## Bug 1：`HeatPipe2D.__init__` 调用 `super().__init__(name=self.name)` 但 `self.name` 未赋值

### 文件
`Components/basicComponents/HeatPipe2D.py` 第 83 行

### 现象
```python
class HeatPipe2D(HeatConduction2D):
    def __init__(self, ..., name="Unnamed_Solid", ...):
        # ... 一堆赋值, 但缺 self.name = name
        super().__init__(mesh, self.wall_mat, name=self.name, initial_temp=initial_temp)
        #                                          ^^^^^^^^^
        # AttributeError: 'HeatPipe2D' object has no attribute 'name'
```

### 报错信息
```
AttributeError: 'HeatPipe2D' object has no attribute 'name'
```

### 触发场景
任何实例化 `HeatPipe2D`（直接或通过 `HPwithFin` / `RingHP`）的代码。

### 建议修复
在 `super().__init__` 调用之前加一行：
```python
self.name = name   # 必须在 super() 调用之前赋值
super().__init__(mesh, self.wall_mat, name=name, initial_temp=initial_temp)
```

---

## Bug 2：`HeatPipe2D._compute_fluxes` 使用旧 API `R_x_inner` / `R_y_inner`，但父类已重构为 `G_*_inner`

### 文件
`Components/basicComponents/HeatPipe2D.py` 第 274 行附近

### 现象
父类 `HeatConduction2D` 在某次重构中把"热阻 R"改为"热导 G"（数值更稳定），并提供 `G_x_inner` / `G_y_inner` 属性。但 `HeatPipe2D` 重写的 `_compute_fluxes` 方法**没跟着改**：

```python
def _compute_fluxes(self, t):
    # ...
    flux_x = (T_2d[:-1, :] - T_2d[1:, :]) / self.R_x_inner   # ❌ R_x_inner 不存在
    flux_y = (T_2d[:, :-1] - T_2d[:, 1:]) / self.R_y_inner   # ❌ R_y_inner 不存在
```

### 报错信息
```
AttributeError: 'HeatPipe2D' object has no attribute 'R_x_inner'.
Did you mean: 'G_x_inner'?
```

Python 提示了正确属性名 `G_x_inner` 

### 建议修复
把 `/ R` 改为 `* G`（数学等价，G = 1/R）：
```python
flux_x = (T_2d[:-1, :] - T_2d[1:, :]) * self.G_x_inner   # ✓
flux_y = (T_2d[:, :-1] - T_2d[:, 1:]) * self.G_y_inner   # ✓
```

---

## Bug 3：`HPwithFin` 调用 `HeatPipe2D` 时未传递 `name` 参数，导致所有热管同名 → SystemManager 注册冲突

### 文件
`Components/HPwithFin.py` 第 90 行附近

### 现象
```python
class HPwithFin(BaseComponent):
    def __init__(self, name, ...):
        super().__init__(name)
        # ...
        self.hp = HeatPipe2D(
            mesh=self.hp_mesh,
            solid1=wall_mat, solid2=fluid_mat, solid3=wick_struct_mat,
            n_wick=n_wick, porosity=porosity,
            n_eva=n_eva, n_aba=n_aba, n_con=n_con,
            emissivity=emissivity, ...
            # ❌ 没有传 name 参数!
        )
```

结果：所有内部 HeatPipe2D 都默认叫 `"Unnamed_Solid"`。

### 报错信息
```
ValueError: Solid component with name 'Unnamed_Solid' already exists in SystemManager!
```

### 触发场景
当一个 `RingHP` 包含 ≥2 个节点的热管时（即 `len(hp_multipliers) ≥ 2`），第二根开始注册时就抛错。

### 建议修复
在 `HPwithFin.__init__` 调用 `HeatPipe2D` 时显式传 `name`：
```python
self.hp = HeatPipe2D(
    name=f"{self.name}_HP_inner",   # 派生唯一名字
    mesh=self.hp_mesh,
    ...
)
```

---

## ~~Bug 4：`HPwithFin` 冷凝段散热双重计算~~ ❌ 撤回 (误诊)

### 撤回说明

> **此 Bug 经过 V4.1 测试已被证伪，原 HPwithFin 散热计算是正确的。**

### 误诊起因

在 V3 / V3.5 / V4 中观察到 `Q_rej / Q_loop ≈ 2.00x`，曾误判为
"冷凝段裸壁辐射与翅片辐射双重计算"。

### 真相

V4.1 通过 `print(return_leg.volumes[-1].T)` 诊断发现：

```
[V4.1 修正后, T_in=968K]
  T_outlet (旧, 读 core_inlet_boundary.T):  862.95 K  ← 假
  T_outlet (新, 读 return_leg.volumes[-1].T): 768.44 K  ← 真实流体温度

  Q_loop (旧, 用假 T_out): 49.93 kW   ← 严重低估
  Q_loop (新, 用真 T_out): 94.86 kW   ← 真实带走的热量
  Q_rej:                  99.95 kW   ← HP 算的辐射散热

  Q_rej / Q_loop = 1.05x  ← 5% 内能量守恒成立 ✓
```

也就是说，**HPwithFin 的散热计算从一开始就是对的**，约 5% 的偏差来自
我用了固定 Cp(T_init) 而不是 Cp(T_avg) 的近似。

### 误诊根源

我们把 `IncompressibleBoundaryVolume.T` 当真实出口流体温度读取，但它默认
是"无限热池" (T 锁定为初始值)，详见下方 **Bug 4 (新)：使用陷阱**。

### 教训

下次遇到"严格整数比"的能量平衡偏差时，先查"读数位置是否正确"，
而不是急着论证"代码计算是否错误"。

---

## Bug 4 (新)：`IncompressibleBoundaryVolume` 默认作"无限热池"，T 锁定为初始值，导致用户从出口边界读温度时拿到错误的值

### 文件
`Solvers/Hydrodynamics/BoundaryVolume.py`

### 现象

`IncompressibleBoundaryVolume` 用作压力出口边界 (`is_pressure_boundary=True`) 时，
其 `T` 属性默认锁定为 `__init__` 时给定的初始值，**不响应**实际流入的流体温度。

```python
# 用户代码:
core_inlet_boundary = IncompressibleBoundaryVolume(
    name="CoreInlet_Outlet", material=nak,
    P=161e3, T=863.0,    # 初始 T = 863
)
core_inlet_boundary.is_pressure_boundary = True

# 仿真跑完后:
print(core_inlet_boundary.T)   # 仍然是 863.0! 即使流入流体温度是 768K
```

代码内部似乎是 `mixing_enabled` 默认为 `False`，导致跳过了流入流体的混合计算。

### 影响

V3 / V3.5 / V4 都用 `core_inlet_boundary.T` 作为"回路出口温度"输出，得到
错误的 ΔT_loop 和 Q_loop。在我们的工况下:
- 表观 T_outlet = 863 K (假)
- 真实 T_outlet = 768 K (从 `return_leg.volumes[-1].T` 读)
- 真实 Q_loop = 95 kW，被表观 Q_loop = 50 kW 严重低估

最坑的地方：`T_init = 863 K` 凑巧接近真实稳态出口温度，所以 V3.5 (T_in=968K) 看起来
"误差很小、像是物理对了"，掩盖了 bug。直到 T_in=1000K 参数扫描才暴露出
"无论入口温度怎么变，出口温度都是 863K"的反常现象。

### 建议修复

两种方案：

**方案 A**：改默认 `mixing_enabled=True`，让边界容器正常计算流入流体的热混合。

**方案 B**：保留默认行为（无限热池在某些场景有用），但在文档/注释里
**显著提示**用户："要读真实出口流体温度，请用 `from_vol.T`（即上游通道末端节点），
不要用 `boundary_volume.T`"。

V4.1 中我用的临时方案 (workaround)：
```python
# ⚠️ 不要读 core_inlet_boundary.T (无限热池)
# T_outlet_real = core_inlet_boundary.T   ← 错!

# ✅ 读上游通道末端节点的真实 T
T_outlet_real = return_leg.volumes[-1].T
```

---

## Bug 5（设计提示）：`SystemManager.step()` 的 Step A 会清空流体 `Q_wall`，导致用户主动设置的流体源项失效

### 文件
`Solvers/SystemManager.py` 的 `step()` 方法

### 现象
在 SystemManager 调度下（V3 必须用），任何在主循环里设置 `vol.Q_wall = Q` 的代码都不生效，因为 `step()` 内部会先清空 `Q_wall` 再让所有 coupler 写入。

但**这一行为没在文档/注释里说明**，导致用户必须做以下三件事之一：
1. 写一个继承 `BaseComponent` 的"伪组件"，在 `pre_step` 里设 Q_wall —— **失败**（pre_step 在 Step A 之前，仍被清掉）
2. 改用 `add_coupling_source(expl, impl)` —— **失败**（同样在 Step A 被清）
3. **唯一可行**：写一个**伪 coupler**（只实现 `execute()` 方法）注册到 `sys_mgr.couplers`，让 SystemManager 在 Step A 里调用它来写 Q_wall

V3 用了方案 3：
```python
class CoreHeatSourceCoupler:
    def execute(self):
        for vol in self.core_channel.volumes:
            vol.Q_wall = self.q_per_node

sys_mgr.couplers.append(core_heater)   # 直接 append 而不是 add_coupler
```

### 建议
两种改进方向，二选一：

**方向 A**：提供一个**官方接口** `sys_mgr.add_persistent_heat_source(volumes, q)`，封装方案 3 的逻辑，避免用户摸索一晚上。

**方向 B**：在 `step()` 的源项清空逻辑里**保留**用户在 `pre_step` 设置的 `Q_wall`（区分 user-set vs coupler-set）。

---

## Bug 6（次要）：`IncompressibleBoundaryVolume.is_pressure_boundary` 默认为 `False`，类名与默认行为有歧义

### 文件
`Solvers/Hydrodynamics/BoundaryVolume.py`

### 现象
```python
class IncompressibleBoundaryVolume:
    def __init__(self, ...):
        # ...
        self.is_pressure_boundary = False   # ⚠️ 类名暗示是边界, 但默认不锚定压力
```

如果用户实例化后忘了手动 `vol.is_pressure_boundary = True`，整个流网会**没有压力锚点**，求解器矩阵会奇异（或出现非物理压力）。

V1→V2 修复就栽在这里 —— V1 把压力锚"侧挂"通过高 K 阻力连接，结果主回路压力飙到 3.17 TPa。

### 建议
两种改法：
1. **改默认值**：`self.is_pressure_boundary = True`（更符合直觉）
2. **强制构造参数**：`__init__(self, ..., is_pressure_boundary)` 不给默认值

---

## 总结表

| # | 严重性 | 文件 | 问题 |
|---|---|---|---|
| 1 | 🔴 高 | HeatPipe2D.py | `__init__` 用 `self.name` 但未赋值 → AttributeError |
| 2 | 🔴 高 | HeatPipe2D.py | `_compute_fluxes` 用废弃的 `R_*_inner` → AttributeError |
| 3 | 🔴 高 | HPwithFin.py | 调用 HeatPipe2D 没传 `name` → SystemManager 注册冲突 |
| 4 | 🟡 中 | BoundaryVolume.py | IncompressibleBoundaryVolume.T 默认锁死 (无限热池) → 出口温度读取陷阱 |
| 5 | 🟡 中 | SystemManager.py | Q_wall 清零行为未文档化 → 用户难以注入流体热源 |
| 6 | 🟢 低 | BoundaryVolume.py | `is_pressure_boundary` 默认 False, 类名有歧义 |
| ~~Bug 4 (旧)~~ | — | HPwithFin.py | ~~冷凝段散热双重计算~~ → ❌ 已撤回 (误诊, 实为 Bug 4 新) |

---

## 附：V3 仿真用到的 monkey-patch 实现

为了让 V3 跑通而不修改源码，我们写了运行时 patch（详见 `test_coolant_loop_v3_final.py` 开头）：

```python
from Components.HPwithFin import HeatPipe2D as _HeatPipe2D_cls

_orig_hp2d_init = _HeatPipe2D_cls.__init__
_hp2d_counter = [0]

def _patched_hp2d_init(self, *args, **kwargs):
    """同时修复 Bug 1 + Bug 3"""
    name = kwargs.get('name', 'Unnamed_Solid')
    if name == 'Unnamed_Solid' or not name:
        _hp2d_counter[0] += 1
        name = f"HeatPipe2D_auto_{_hp2d_counter[0]}"
        kwargs['name'] = name
    self.name = name
    _orig_hp2d_init(self, *args, **kwargs)

def _patched_compute_fluxes(self, t):
    """修复 Bug 2: R → G"""
    # ... (用 self.G_x_inner / self.G_y_inner 替代)
    pass

_HeatPipe2D_cls.__init__ = _patched_hp2d_init
_HeatPipe2D_cls._compute_fluxes = _patched_compute_fluxes
```

这些 patch 是临时方案，正式版本应该回到源码里修。

---

