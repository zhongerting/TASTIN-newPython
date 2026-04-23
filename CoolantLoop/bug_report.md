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

## Bug 4：`HPwithFin` 冷凝段散热**双重计算**，导致总散热量偏高 ~2 倍

### 文件
`Components/HPwithFin.py` 第 2.2 节注释附近

### 现象
冷凝段（外径裸壁 + 翅片）的散热路径有两条：
1. **裸壁直视辐射**：`q_bare = ε σ A_con (T⁴ - T_space⁴)`
2. **翅片排热**（准稳态）：`q_fin = (T_con - T_space) / R_fin_eq`

代码里两条路径**互相独立计算并相加**，但**裸壁辐射用了 100% 面积**（`bare_area_array=area_con`），没有扣除被翅片包覆的部分：

```python
# 2.2 冷凝段支路 A: 热管外壁面直视辐射 (100% 面积参与辐射)
self.bc_rad_con = ... add_dynamic_radiation_condition(
    bare_area_array=area_con,    # ❌ 应该是 area_con × (1 - wrap_ratio)
    ...
)
```

物理上：被翅片包覆的那部分壁面**不该再算裸壁辐射**（它的散热已经通过翅片路径计算了）。

### 影响
我们 V3 仿真显示 `Q_rej / Q_input ≈ 1.92`（应该接近 1.0），刚好对应"裸壁面积全算 + 翅片再算一遍"的情形。这导致：
- 集流环平衡温度被人为压低（散热被高估）
- 长瞬态会出现"散热持续超过加热"的非物理现象（300s 仿真未达稳态）

### 建议修复
裸壁散热扣除翅片包覆面积：
```python
self.bc_rad_con = ... add_dynamic_radiation_condition(
    bare_area_array=area_con * (1 - self.fin_wrap_ratio),   # ✓ 扣除翅片占用部分
    ...
)
```

或者另一种思路：让翅片路径和裸壁路径**严格按面积比例划分**，互不重叠。

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
| 4 | 🟡 中 | HPwithFin.py | 冷凝段散热双重计算 → Q_rej 偏高 ~2x |
| 5 | 🟡 中 | SystemManager.py | Q_wall 清零行为未文档化 → 用户难以注入流体热源 |
| 6 | 🟢 低 | BoundaryVolume.py | `is_pressure_boundary` 默认 False, 类名有歧义 |

---

## 附：V3 仿真用到的 monkey-patch 实现

为了让 V3 跑通而不修改源码，我们写了运行时 patch（详见 `test_coolant_loop_v3 (3)_final.py` 开头）：

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

## 求确认

> **最终只仿真冷却剂回路（不含堆芯）**，把堆芯出口当作流量入口边界、堆芯入口当作压力出口边界。命名为 V3.5（开环架构）。


