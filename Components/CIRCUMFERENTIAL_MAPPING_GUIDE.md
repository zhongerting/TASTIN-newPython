# 周向分区映射指南

## 1. 用途

`Components/circumferential_mapping.py` 用于把一个完整圆周上的均匀源分区映射到另一套均匀目标分区。首个目标场景是把 V14 的 18 个代表热管温度区映射为遮热罩 `fortran_shield2` 所需的 12 个周向辐射器区。

当前工具仅提供映射计算，尚未替换 `RadiatorThermalShield` 的现有分区逻辑。

## 2. 分配矩阵

设源区数量为 `Ns`，目标区数量为 `Nt`：

```text
source_width = 360 / Ns
target_width = 360 / Nt
F[j,i] = overlap(target_j, source_i) / source_width
```

`F` 的形状为 `(Nt, Ns)`。每个源区会被完整分配，因此每列之和为 1：

```python
mapping = build_uniform_circumferential_mapping(18, 12)
assert np.allclose(mapping.sum(axis=0), 1.0)
```

源区和目标区可以分别指定起始角：

```python
mapping = build_uniform_circumferential_mapping(
    18,
    12,
    source_offset_deg=0.0,
    target_offset_deg=10.0,
)
```

跨越 `360°/0°` 的区间由函数自动处理。

## 3. 总量映射

热功率、面积和热管数量属于总量，直接乘分配矩阵：

```python
target_power = mapping @ source_power
```

由于矩阵列守恒：

```text
sum(target_power) = sum(source_power)
```

## 4. 强度量映射

`T^4` 和热流密度属于强度量，需要按源区代表的物理权重归一化：

```python
target_t4 = map_circumferential_intensive(
    source_t4,
    mapping,
    source_weights=hp_multipliers,
)
```

计算关系为：

```text
numerator[j]   = sum_i(F[j,i] * weight[i] * value[i])
denominator[j] = sum_i(F[j,i] * weight[i])
target[j]      = numerator[j] / denominator[j]
```

调用方必须传入已经计算好的 `mean(T_cells ** 4)`，工具不会把温度自动变成四次方。

## 5. V14 的 18→12 映射

18 个源区每区宽 `20°`，12 个目标区每区宽 `30°`。零偏置时：

```text
target[2m]   = source[3m] 全部 + source[3m+1] 一半
target[2m+1] = source[3m+1] 一半 + source[3m+2] 全部
m = 0...5
```

例如前三个 V14 倍率为 `[5, 6, 6]`：

```text
Tc4[0] = (5*T4[0] + 3*T4[1]) / 8
Tc4[1] = (3*T4[1] + 6*T4[2]) / 9
```

V14 的 `symmetric_ring_multiplier=2` 对归一化后的 `T^4` 不产生影响，但后续计算总辐射面积和总热功率时必须保留。

## 6. 批量数据

`source_values` 的最后一维必须是源区维度，因此支持批量映射：

```python
source_t4_history.shape == (n_times, 18)
target_t4_history = map_circumferential_intensive(
    source_t4_history,
    mapping,
    source_weights=hp_multipliers,
)
assert target_t4_history.shape == (n_times, 12)
```

## 7. 输入约束

以下情况会抛出 `ValueError`：

- 分区数量不是正整数；
- 角度偏置不是有限数；
- 数值或权重包含非有限值；
- 权重为负；
- 源数据最后一维与矩阵不匹配；
- 某个目标区没有获得任何正物理权重。

## 8. 后续接入

V14 遮热罩接入还需单独完成：

1. 从 18 个代表热管区提取 `mean(T^4)`；
2. 使用本工具聚合为 12 个遮热罩输入区；
3. 将遮热罩返回的 12 区背景 `T^4` 按重叠关系回写到 18 区；
4. 用 200 根真实热管的面积倍率检查能量闭合；
5. 校核原 78 根辐射管角系数对 200 根热管结构的适用性。
