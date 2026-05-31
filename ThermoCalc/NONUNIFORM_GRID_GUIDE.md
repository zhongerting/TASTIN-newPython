# Thermionic Energy Conversion Solver - Non-Uniform Grid Support Guide

## 1. 修改总结 (Summary of Changes)

本次重构使求解器从原有的“均匀轴向网格”限制中解放出来，现在支持任意比例的非均匀节点分布（如局部加密、带有端部反射层的非等长划分）。

### C++ 核心层修改
- **数据结构升级**：在 `singleThermionicEnergyConversion.h` 中，将 `sideAreaE` 和 `sideAreaC` 从标量修改为 `vector<double>`，支持每个节点拥有独立的侧面积。
- **积分逻辑重构**：修改了 `singleThermionicEnergyConversion.cpp` 中所有涉及电流和功率的积分循环。现在使用局部的 `sideAreaE[i]` 计算电流贡献，并在计算截面电流（`IEsec`/`ICsec`）时采用了修正的加权梯形公式，确保非均匀网格下的电荷守恒。
- **电势求解修正 (FVM)**：更新了 `VcalcFVM()` 函数。节点间的界面电导（Face Conductance）计算现在使用精确的物理中心距 `0.5 * (dl[i] + dl[i+1])`，从而在非均匀几何下维持静电场方程的计算精度。

### 接口与脚本修改
- **Python 绑定**：更新 `bindings.cpp`，使得 `InputData` 能够接受二维 NumPy 数组作为侧面积输入（$N_{elements} \times n_{nodes}$）。
- **自动化适配**：批量更新了所有测试脚本，确保输入数据的维度与重构后的接口匹配。

---

## 2. 迁移与使用说明 (Migration & Usage Instructions)

如果您要将此模型迁移到新环境，请遵循以下配置规范：

### 2.1 输入数据准备 (Python 端)
在定义 `InputData` 时，必须确保以下四个分布参数的维度完全一致且物理匹配：

1.  **节点步长 (`dlE`, `dlC`)**：定义每个小段的轴向长度。
2.  **侧面积 (`sideAreaE`, `sideAreaC`)**：定义该小段对应的发射极/接收极表面积。
    - *注意*：对于等截面圆柱电极，如果某一段长度为 $dl_i$，则其面积应设为 $A_{total} \cdot (dl_i / L_{total})$。

**示例代码：**
```python
# 假设 6个反射层(0.0108m) + 25个活性区(0.015m) + 6个反射层(0.0108m)
dl_array = np.array([0.0108]*6 + [0.015]*25 + [0.0108]*6)
input_data.dlE = np.tile(dl_array, (N_elem, 1))

# 面积也要按比例分配
sideE_mid = 0.000928  # 活性区单节点面积参考值
sideE_ref = sideE_mid * (0.0108 / 0.015)
sideE_array = np.array([sideE_ref]*6 + [sideE_mid]*25 + [sideE_ref]*6)
input_data.sideAreaE = np.tile(sideE_array, (N_elem, 1))
```

### 2.2 物理边界对齐要求
- **节点对应**：目前的物理模型（如 $J$ 的计算）假定第 $i$ 个发射极节点与第 $i$ 个接收极节点在空间上是正对着的。因此，虽然支持 `dlE` 和 `dlC` 不同，但强烈建议保持两者分布一致，以确保径向间隙 $d_{gap}$ 定义的准确性。
- **累积长度**：确保 $\sum dl$ 等于您设计的总物理长度，否则会导致计算出的总欧姆电阻偏离预期。

### 2.3 编译要求
- 必须使用支持 **C++17** 的编译器。
- 确保链接了最新的 `bindings.cpp` 和 `singleThermionicEnergyConversion.cpp` 源码。
- 建议使用 CMake 进行构建，CMakeLists.txt 已配置为自动处理 `pybind11`。

---

## 3. 性能与验证建议

- **基线测试**：迁移后请优先运行 `test_real_case_v2.py`。如果在均匀网格下结果与历史数据（如 `Iout=701.57A`）一致，则证明基础环境配置正确。
- **非均匀验证**：运行 `test_nonuniform_v2_37rings_corrected.py`。该脚本涵盖了真实的非均匀物理映射，是验证“几何-温度-电阻”耦合逻辑的最佳案例。
- **缓存建议**：目前的缓存机制通过完整数组比对，对非均匀网格是安全的。只要不动态在求解中途改变 `dl` 划分，不需要额外维护。

---
*文档生成日期：2026年5月19日*
