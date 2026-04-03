# TASTIN-Python

TASTIN-Python 是一个基于 Python 的航天器核热耦合系统仿真平台，用于模拟空间反应堆的瞬态热工水力行为。

## 项目简介

TASTIN-Python 是 TASTIN 程序的 Python 重构版本，采用了模块化和面向对象的设计思想，支持多物理场耦合计算。系统集成了热传导、流体网络、中子动力学和热电转换等核心模块。

## 主要功能

### 核心组件 (Components)

| 组件                    | 说明                                       |
| --------------------- | ---------------------------------------- |
| `ReactorCore`         | 堆芯宏观容器，集成点堆动力学和热电耦合电路                    |
| `TFEUnit`             | 热离子能量转换单元 (Thermionic Energy Conversion) |
| `RingHP`              | 环形热管组件                                   |
| `HPwithFin`           | 带鳍片的热管组件                                 |
| `AnnularPipe`         | 环形管道                                     |
| `Pipe`                | 管道组件                                     |
| `TECCircuitManager`   | 热电耦合电路管理器                                |
| `ExternalHeatSources` | 外部热源（太阳直射、地球反照、地球红外）                     |

### 求解器 (Solvers)

- **热传导求解器** (`HeatConduction/`): 2D 瞬态热传导计算
- **流体网络求解器** (`Hydrodynamics/`): 液压网络系统仿真
- **中子动力学求解器** (`Neutronics/`): 点堆动力学模型
- **系统管理器** (`SystemManager`): 多物理场协调调度与数据同步

### 材料库 (Materials)

#### 流体材料

- `Sodium`: 钠 (Na) - 快堆冷却剂
- `NaK78`: 钠钾合金 (NaK78)
- `Potassium`: 钾 (K)

#### 固体材料

- `UO2`: 二氧化铀 - 燃料
- `StainlessSteel`: 不锈钢 - 结构材料
- `Molybdenum`: 钼 - 结构材料
- `BerylliumOxide`: 氧化铍 - 慢化剂
- `ZrH`: 氢化锆 - 慢化剂
- `NaHP`: 钠热管材料
- `WickMaterial`: 热管毛细材料
- `WallMaterial`: 壁面材料

### 相关公式 (Correlations)

热工计算中使用的经验关联式，包括传热、压降、沸腾换热等。

### 热电耦合计算 (ThermoCalc)

C++ 实现的 热离子能量转换 (Thermionic Energy Conversion) 计算模块，支持:

- 电路 TEC 建模
- 非线性求解器

## 系统架构

```
TASTIN-Python
├── Components/          # 宏观组件层
│   ├── BaseComponent.py  # 组件基类
│   ├── ReactorCore.py     # 堆芯组件
│   ├── TFEUnit.py         # 热电转换单元
│   ├── ExternalHeatSources/  # 外部热源
│   └── basicComponents/   # 基础组件
├── Solvers/              # 求解器层
│   ├── HeatConduction/    # 热传导求解器
│   ├── Hydrodynamics/     # 流体动力学求解器
│   ├── Neutronics/        # 中子动力学求解器
│   ├── SystemManager.py   # 系统管理器
│   └── Couplers.py       # 耦合器
├── Materials/            # 材料属性层
│   ├── Fluids/           # 流体材料
│   └── Solids/           # 固体材料
├── Correlations/         # 工程关联式
├── MathSolvers/          # 数学求解器
└── ThermoCalc/           # C++ 热电耦合计算
```

## 技术特点

1. **模块化设计**: 基于抽象基类的插件式架构
2. **多物理场耦合**: 热传导、流体流动、中子动力学、TEC 电计算
3. **瞬态仿真**: 支持时间步进仿真
4. **轨道热分析**: 内置轨道环境热流模型（太阳辐射、地球反照、红外热辐射）
5. **Fortran 兼容性**: 参考 Fortran 版本算法设计

## 依赖环境

- Python 3.8+
- NumPy
- SciPy

## 使用方法

### 基本仿真流程

```python
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager

# 创建流体网络
fluid_network = HydraulicNetwork(...)

# 创建系统管理器
system = SystemManager(fluid_network)

# 添加组件
system.add_solid_component(reactor_core)

# 运行瞬态仿真
system.run_transient(end_time=3600.0, dt=1.0)
```

## 目录说明

- `TASTIN-python.7z`: 完整项目压缩包
- `.idea/`: PyCharm IDE 配置
- `__pycache__/`: Python 字节码缓存（自动生成）

## 致谢

<br />

## 许可证

本项目仅供研究使用。
