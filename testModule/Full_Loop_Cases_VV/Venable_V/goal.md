# Venable_V 验证算例目标与执行约束

本文用于说明 `Venable_V` 后续验证算例计算的目标、允许操作范围和过程记录要求。后续所有 Venable Table 7-1 单 TFE 最大输出功率验证工作，应优先遵守本文，再参考 `PARAMETER_ADJUSTMENT_GUIDE.md` 进行参数调整。

## 1. 最终目标

最终目标是完成 Venable 1995 单电池 TFE 电特性试验 Table 7-1 数据的 V&V 验证工作。

验证对象是文献给出的 14 个 active-zone power 工况下的最大输出功率和效率。计算侧应在每个工况下建立对应的单 TFE 电输出模型，通过负载扫描、I-V 曲线扫描或等价优化方式得到计算最大输出功率，再与文献数据对比。

核心输出应至少包括：

- 每个工况的 `Q_az`、文献 `P_out,max`、计算 `P_out,max`、文献效率、计算效率；
- 每个工况的绝对误差和相对误差；
- 14 个工况的整体误差指标；
- 当前输入参数、热边界假设、Cs 压力/温度映射和求解设置；
- 是否存在未收敛、非有限值、边界不合理或趋势异常。

## 2. 代码修改边界

如果计算过程中遇到报错，首先应从算例本身找原因，包括但不限于：

- 输入参数单位是否错误；
- 工况边界是否超出模型可用范围；
- 温度、压力、电压、电流、负载扫描范围是否不合理；
- `Pcs -> Tcs` 映射是否给出了不物理的 Cs 储罐温度；
- 初值是否导致 ThermoCalc 或电路求解进入退化状态；
- 单 TFE 几何、面积、节点数、热边界数组形状是否与接口要求不一致；
- 输出后处理是否错误解释了模型结果。

除非确认是恶性 bug，否则不得修改原本程序的 Python 或 C++ 源码文件。这里的“原本程序源码”包括仓库主体模块中的 `.py`、`.cpp`、`.h`、pybind11 绑定和 ThermoCalc 后端实现。

允许优先修改或新增的内容包括：

- `Venable_V` 算例输入文件；
- `Venable_V` 下的工况构建脚本；
- `Venable_V` 下的后处理脚本；
- `Venable_V` 下的参数配置、结果表和过程记录文档；
- 只服务于该验证算例的 wrapper、runner 或 smoke test。

如果怀疑主体代码存在恶性 bug，必须先记录复现条件、最小输入、报错栈、影响范围和为什么无法通过算例输入规避。只有在确认不是工况设置问题后，才讨论是否修改主体代码。

## 3. 验证偏差较大时的调参原则

如果测试验证效果较差，即目标参数结果偏差较大，应按照 `PARAMETER_ADJUSTMENT_GUIDE.md` 的边界进行输入参数调整，尽量使误差最小化。

调参优先级为：

1. 先排除数值问题：扫描范围、扫描分辨率、网格无关性、求解容差、非有限值。
2. 再修正明确的输入映射问题，例如 `Pcs -> Tcs`、单位换算、数组形状、面积口径。
3. 再替换 placeholder 物理闭合，例如 emitter/collector 温度闭合、轴向温度分布、热损失分配。
4. 再根据文献或几何推导修正面积、active length、gap 等结构参数。
5. 最后才考虑有界的全局表面参数或发射参数校准。

不得通过逐工况任意调参来强行压低误差。正式验证模型应尽量使用同一套全局参数解释 14 个工况的整体趋势。

## 4. 实事求是原则

如果在合理调参后仍无法让误差继续降低，应如实说明结果，而不是继续增加缺乏依据的自由参数。

最终报告应区分：

- 已严格对应文献的实验锚点；
- 已由文献或几何推导确认的模型参数；
- 仍为 placeholder 或工程假设的输入；
- 做过但被拒绝的敏感性调整；
- 误差无法继续降低的可能原因。

可接受的结论包括但不限于：

- 当前公开文献缺少某些关键边界，导致只能完成趋势验证或半定量验证；
- 现有 ThermoCalc 模型的物理假设与 Venable 单电池试验台存在差异；
- 温度闭合或 Cs 压力映射不确定性主导误差；
- 单 TFE 试验台外部电路、表面状态或热损失路径无法从当前资料中完全复原。

不能把无法验证成功简单归因于“程序错误”，也不能把调参后凑出的结果描述成无条件验证成功。

## 5. 过程跟踪要求

在算例计算、验证和参数修改过程中，必须创建并维护过程跟踪 `.md` 文档，保证后续可以追溯每一步调整，防止参数反复横跳。

建议在 `Venable_V` 下维护以下记录文件：

- `validation_process_log.md`：记录每次运行、报错、修正、结果摘要和下一步判断；
- `adjustment_log.md`：记录每次参数调整的原因、数值、范围、误差变化和是否保留；
- `validation_result_summary.md`：在阶段性验证完成后整理最终模型、误差表和结论。

每次参数调整建议使用以下格式：

```text
Date:
Run ID:
Parameter:
Previous value:
New value:
Scope: global / case group / single diagnostic case
Source status: directly sourced / derived / sensitivity / fitted
Rationale:
Affected files:
Baseline error summary:
New error summary:
Decision: keep / reject / diagnostic only
Next action:
```

每次计算运行建议记录：

```text
Date:
Run ID:
Command:
Input files:
Output files:
Code version / git status:
ThermoCalc mode:
Case count:
Convergence status:
Main result:
Observed problem:
Decision:
Next action:
```

## 6. 阶段性完成标准

一个阶段的 Venable_V 验证工作只有在满足以下条件时，才可以认为完成：

- 14 个 Table 7-1 工况均完成计算或明确记录无法计算的原因；
- 输出表中包含文献值、计算值和误差；
- 所有输入参数来源和调整历史可追溯；
- 若误差较大，已经按 `PARAMETER_ADJUSTMENT_GUIDE.md` 做过合理排查；
- 若仍无法降低误差，已经给出清晰、具体、可复查的原因分析；
- 没有在未确认恶性 bug 的情况下修改主体 Python/C++ 源码。

## 7. 启动 goal 模式前的补充约束

后续如果使用 goal 模式执行本文任务，应按阶段推进，不要一开始直接进入全量调参。

建议执行阶段如下：

1. `precheck`：检查解释器、ThermoCalc 导入、当前输入文件、输出目录、已有日志和 git 状态。
2. `single_case_smoke`：先选 1 个代表性工况做最小计算，确认输入形状、求解入口和后处理链条能跑通。
3. `baseline_14_cases`：使用当前 baseline 参数完成 14 个 Table 7-1 工况，不做调参。
4. `error_diagnosis`：分析误差随 `Q_az` 和 `Pcs` 的趋势，区分数值问题、输入映射问题和物理闭合问题。
5. `parameter_adjustment_loop`：按 `PARAMETER_ADJUSTMENT_GUIDE.md` 逐类调整输入参数，每轮只改变一类主导参数。
6. `stage_summary`：整理当前模型、误差、保留/拒绝的参数调整和下一阶段建议。

### 7.1 恶性 bug 判定标准

“恶性 bug”不是指验证结果偏差大，也不是指某个工况不收敛。只有满足以下条件之一时，才应考虑修改主体 Python/C++ 源码：

- 用最小 Venable_V 输入仍能稳定复现崩溃、非有限数组或接口异常；
- 报错由主体接口实现和公开调用约定明显矛盾导致，且无法通过算例输入规避；
- 同一固定输入在无随机源的情况下出现非确定性结果；
- 错误影响主体模块的通用行为，不只影响 Venable_V 的某个临时输入；
- 后处理确认输出字段、单位或数组含义与主体实现存在确定性错误。

即使怀疑恶性 bug，也必须先在过程日志中记录最小复现、报错栈、输入文件、影响范围和为什么不能通过算例侧修正。

### 7.2 运行环境要求

执行 Venable_V 算例脚本时，默认使用项目指定解释器：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" <script>
```

进入正式计算前应做 precheck：

- 解释器路径是否为项目指定 Conda 环境；
- `numpy` 等基础依赖是否可导入；
- ThermoCalc/Python 绑定是否能在当前环境中导入；
- `Venable_V` 输入文件是否完整；
- 当前 `PARAMETER_ADJUSTMENT_GUIDE.md` 和 `goal.md` 是否与执行意图一致；
- 是否存在未记录的旧输出会被覆盖。

如果 ThermoCalc ABI 或环境导入失败，应先作为环境/入口问题记录，不应直接修改主体源码。

### 7.3 run_id 和输出目录规则

后续每次计算运行必须使用唯一 `run_id`，建议格式：

```text
YYYYMMDD_HHMMSS_<stage>_<short_note>
```

例如：

```text
20260701_231500_single_case_smoke_qaz3162
20260701_234000_baseline_14_cases
```

建议输出目录结构：

```text
Venable_V/
  runs/
    <run_id>/
      input_snapshot/
      results/
      logs/
      plots/
      run_summary.json
```

要求：

- 不覆盖已有 `run_id` 目录；
- 每轮运行的 CSV/JSON/日志都放入对应 `runs/<run_id>/`；
- `validation_process_log.md` 必须引用对应 `run_id`；
- 参数调整前后结果必须能追溯到具体 `run_id`；
- 如果需要临时诊断输出，应标记为 diagnostic，不得混入正式结果。

### 7.4 参数调整循环规则

进入参数调整阶段后，每轮只改变一类主要参数，例如只改 `Pcs -> Tcs`，或只改温度闭合，或只改面积口径。不要在同一轮同时改变多个强影响参数，否则无法判断误差改善来自哪里。

每轮调整后至少记录：

- 调整前 baseline 或上一轮 `run_id`；
- 调整后新 `run_id`；
- 改动参数和取值范围；
- 14 个工况误差变化；
- 趋势是否变得更物理；
- 该调整是保留、拒绝还是仅作诊断。

正式模型只能吸收全局参数调整或有文献/几何依据的调整，不吸收逐点拟合参数。

### 7.5 停止准则

当完成以下工作后，如果误差仍无法继续降低，应停止继续增加自由参数，转入原因分析和阶段总结：

- 已排除扫描范围、扫描分辨率、网格、容差和后处理错误；
- 已检查并修正明确的单位、数组形状和输入映射问题；
- 已替换或评估 `Pcs -> Tcs` 映射；
- 已替换或评估 emitter/collector 温度闭合；
- 已检查几何面积、active length 和 gap 的来源；
- 已完成有限的全局表面/发射参数敏感性分析；
- 继续调参只能依赖逐工况自由参数或缺乏依据的经验倍率。

达到停止准则后，阶段性结论应如实写明：当前误差水平、最可能的限制因素、哪些参数仍缺少公开文献约束，以及后续若要继续降低误差需要补充哪些实验信息。
