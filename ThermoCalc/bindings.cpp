#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace std;

#include "circuitTECs.h"
#include "emissionLookup.h"
#include "singleThermionicEnergyConversion.h"
#include "thermionicEmission.h"

// -------------------------------------------------------------------------
// 1. 枚举与辅助结构体定义
// -------------------------------------------------------------------------

// 电路计算模式枚举
enum class CalculationMode {
    FixedCurrent,
    FixedVoltage, // 固定电压模式
    FixedResistance, // 固定电阻模式
    ParallelFixedVoltage,
    ParallelFixedCurrent,
    ParallelLoadCurve
};

struct InputData {
    // ---------------------------------------------------------
    // 1. 维度控制
    // ---------------------------------------------------------
    int N_elements; // 串联元件的总数量 (N)
    int n_axi;      // 轴向离散节点数量 (用于校验数据尺寸)

    // ---------------------------------------------------------
    // 2. 物理场分布数据 (二维数组: N x n_axi)
    // 对应构造函数中的 vector<double> 类型的分布参数
    // ---------------------------------------------------------
    py::array_t<double> Temitter;   // [0] 发射极温度分布
    py::array_t<double> Tcollector; // [1] 接收极温度分布
    py::array_t<double> dlE;        // [2] 发射极单元长度分布 (通常是网格步长)
    py::array_t<double> dlC;        // [3] 接收极单元长度分布
    py::array_t<double> Tcs;        // [9] 铯蒸汽温度分布
    py::array_t<double> V_init;     // [10] 极板电势差初值分布 (guess values)

    // ---------------------------------------------------------
    // 3. 几何与标量参数 (一维数组: N)
    // 对应构造函数中 input[4], input[5]... 等标量参数
    // 每个元件有一个独立的值
    // ---------------------------------------------------------
    py::array_t<double> crossAreaE; // [4][0] 发射极横截面积
    py::array_t<double> crossAreaC; // [4][1] 接收极横截面积
    py::array_t<double> sideAreaE;  // [5] 发射极逐节点侧面积
    py::array_t<double> sideAreaC;  // [6] 接收极逐节点侧面积
    py::array_t<double> U_init;     // [8][0] 总电压初值
    py::array_t<double> d_gap;      // [8][1] 电极间距
    py::array_t<double> Itarget;    // [11][0] 目标电流 (虽然是 output 初始猜测，但也需要输入)

    // ---------------------------------------------------------
    // 4. 边界与导线参数 (二维数组: N x 4)
    // 对应构造函数中 input[7] 和 input[12]
    // 假设每个元件有 4 个导线电阻值/电压边界值
    // ---------------------------------------------------------
    py::array_t<double> resistanceWire; // [7] 导线电阻 [r_wire0, r_wire1, r_wire2, r_wire3]
    py::array_t<double> wireU;          // [12] 导线电压初值 [u0, u1, u2, u3]

    // ---------------------------------------------------------
    // 5. 全局电路控制参数 (直接赋值给 circuitTECs)
    // ---------------------------------------------------------
    CalculationMode mode;   // 计算模式 (固定电压/电阻)
    double target_val;      // 目标值 (Utarget 或 Rload)
    double I_total_init;    // 整个电路的初始电流猜测 Iout
    double R_load_init;     // 初始负载电阻 (如果是定电压模式可能用不到，但留着备用)
    py::array_t<double> loadCurveCurrent; // 并联负载曲线电流轴 [A]
    py::array_t<double> loadCurveVoltage; // 并联负载曲线电压轴 [V]
};

// -------------------------------------------------------------------------
// 2. 辅助函数：NumPy 数据提取
// -------------------------------------------------------------------------

// 从 2D 数组 (N x M) 中提取第 i 行，转换为 std::vector<double>
std::vector<double> get_row_vector(const py::array_t<double>& arr, int i) {
    auto r = arr.unchecked<2>(); // 不检查边界，速度快 (2D)
    int n_cols = r.shape(1);
    std::vector<double> vec(n_cols);
    for (int j = 0; j < n_cols; ++j) {
        vec[j] = r(i, j);
    }
    return vec;
}

// 从 1D 数组 (N) 中提取第 i 个元素
double get_scalar(const py::array_t<double>& arr, int i) {
    auto r = arr.unchecked<1>();
    return r(i);
}

std::vector<double> get_vector_1d(const py::array_t<double>& arr) {
    auto r = arr.unchecked<1>();
    std::vector<double> vec(static_cast<std::size_t>(r.shape(0)));
    for (py::ssize_t i = 0; i < r.shape(0); ++i) {
        vec[static_cast<std::size_t>(i)] = r(i);
    }
    return vec;
}

void require_1d(const char* name, const py::array_t<double>& arr, int n_rows) {
    if (arr.ndim() != 1 || arr.shape(0) != n_rows) {
        throw py::value_error(
            string(name) + " must have shape (" + to_string(n_rows) + ",)."
        );
    }
}

void require_2d(const char* name, const py::array_t<double>& arr, int n_rows, int n_cols) {
    if (arr.ndim() != 2 || arr.shape(0) != n_rows || arr.shape(1) != n_cols) {
        throw py::value_error(
            string(name) + " must have shape (" + to_string(n_rows) + ", "
            + to_string(n_cols) + ")."
        );
    }
}

void validate_input(const InputData& data) {
    if (data.N_elements <= 0) {
        throw py::value_error("N_elements must be positive.");
    }
    if (data.n_axi <= 0) {
        throw py::value_error("n_axi must be positive.");
    }
    require_2d("Temitter", data.Temitter, data.N_elements, data.n_axi);
    require_2d("Tcollector", data.Tcollector, data.N_elements, data.n_axi);
    require_2d("dlE", data.dlE, data.N_elements, data.n_axi);
    require_2d("dlC", data.dlC, data.N_elements, data.n_axi);
    require_2d("sideAreaE", data.sideAreaE, data.N_elements, data.n_axi);
    require_2d("sideAreaC", data.sideAreaC, data.N_elements, data.n_axi);
    require_2d("Tcs", data.Tcs, data.N_elements, data.n_axi);
    require_2d("V_init", data.V_init, data.N_elements, data.n_axi);
    require_2d("resistanceWire", data.resistanceWire, data.N_elements, 4);
    require_2d("wireU", data.wireU, data.N_elements, 4);
    require_1d("crossAreaE", data.crossAreaE, data.N_elements);
    require_1d("crossAreaC", data.crossAreaC, data.N_elements);
    require_1d("U_init", data.U_init, data.N_elements);
    require_1d("d_gap", data.d_gap, data.N_elements);
    require_1d("Itarget", data.Itarget, data.N_elements);
    if (data.mode == CalculationMode::ParallelLoadCurve) {
        if (data.loadCurveCurrent.ndim() != 1 || data.loadCurveVoltage.ndim() != 1 ||
            data.loadCurveCurrent.shape(0) != data.loadCurveVoltage.shape(0) ||
            data.loadCurveCurrent.shape(0) < 2) {
            throw py::value_error("ParallelLoadCurve requires loadCurveCurrent/loadCurveVoltage with matching shape (n>=2,).");
        }
    }
}

// -------------------------------------------------------------------------
// 3. 工厂函数：构建电路
// -------------------------------------------------------------------------

// 这是一个独立的 C++ 函数，负责“翻译” InputData 并创建对象树
// 返回 circuitTECs 指针，所有权移交给 Python
std::unique_ptr<circuitTECs> create_circuit(const InputData& data) {
    validate_input(data);
    auto circuit = std::make_unique<circuitTECs>();

    // 预先分配内存，防止 realloc
    circuit->TECs.reserve(data.N_elements);

    // 循环创建每一个 TEC 单元
    for (int i = 0; i < data.N_elements; ++i) {
        // 准备构造函数需要的 input 向量 (共13个槽位)
        std::vector<std::vector<double>> input(13);

        // [0-3] 向量分布参数
        input[0] = get_row_vector(data.Temitter, i);
        input[1] = get_row_vector(data.Tcollector, i);
        input[2] = get_row_vector(data.dlE, i);
        input[3] = get_row_vector(data.dlC, i);

        // [4] 截面积 {E, C}
        input[4] = { get_scalar(data.crossAreaE, i), get_scalar(data.crossAreaC, i) };

        // [5-6] 逐节点侧面积
        input[5] = get_row_vector(data.sideAreaE, i);
        input[6] = get_row_vector(data.sideAreaC, i);

        // [7] 导线电阻 (向量)
        input[7] = get_row_vector(data.resistanceWire, i);

        // [8] {U, d}
        input[8] = { get_scalar(data.U_init, i), get_scalar(data.d_gap, i) };

        // [9-10] Tcs, V
        input[9] = get_row_vector(data.Tcs, i);
        input[10] = get_row_vector(data.V_init, i);

        // [11] {Itarget}
        input[11] = { get_scalar(data.Itarget, i) };

        // [12] wireU
        input[12] = get_row_vector(data.wireU, i);

        // --- 核心创建逻辑 ---
        // 1. new 对象 (堆上分配)
        singleThermionicEnergyConversion* tec = new singleThermionicEnergyConversion(input);
        
        // 2. 初始化 (分配内部内存，构建 thermionicUnits 池)
        tec->initial();

        // 3. 将指针移交给电路对象
        // 注意：根据您原本的代码，circuitTECs 使用 raw pointers。
        // Python 销毁 circuitTECs 时，如果 circuitTECs 析构函数没写 delete，这里会内存泄漏。
        // 但为了保持与您 C++ 逻辑一致，我们这里只负责传递。
        circuit->TECs.push_back(tec);
    }

    // --- 设置全局电路参数 ---
    circuit->nTECs = data.N_elements;
    circuit->Iout = data.I_total_init;
    circuit->isFixedU = false;
    circuit->isFixedI = false;
    circuit->isFixedR = false;
    circuit->isParallelFixedU = false;
    circuit->isParallelFixedI = false;
    circuit->isParallelLoadCurve = false;
    circuit->Utarget = data.target_val;
    circuit->Itarget = data.target_val;
    circuit->Rload = data.target_val;

    if (data.mode == CalculationMode::FixedCurrent) {
        circuit->isFixedI = true;
        circuit->Itarget = data.target_val;
        circuit->Iout = data.target_val;
        circuit->Uout = 0.0;
        for (int i = 0; i < data.N_elements; ++i) {
            circuit->Uout += get_scalar(data.U_init, i);
        }
    }
    else if (data.mode == CalculationMode::FixedVoltage) {
        circuit->isFixedU = true;
        circuit->Utarget = data.target_val;
    } 
    else if (data.mode == CalculationMode::FixedResistance) {
        circuit->isFixedR = true;
        circuit->Rload = data.target_val;
    }
    else if (data.mode == CalculationMode::ParallelFixedVoltage) {
        circuit->isParallelFixedU = true;
        circuit->Utarget = data.target_val;
        circuit->Uout = data.target_val;
    }
    else if (data.mode == CalculationMode::ParallelFixedCurrent) {
        circuit->isParallelFixedI = true;
        circuit->Itarget = data.target_val;
        circuit->Iout = data.I_total_init;
    }
    else if (data.mode == CalculationMode::ParallelLoadCurve) {
        circuit->isParallelLoadCurve = true;
        circuit->Utarget = data.target_val;
        circuit->Uout = data.target_val;
        circuit->setLoadCurve(get_vector_1d(data.loadCurveCurrent), get_vector_1d(data.loadCurveVoltage));
    }

    return circuit;
}

py::dict calc_emission_point(double TE, double TC, double Vo, double Tcs, double d_gap) {
    if (!std::isfinite(TE) || !std::isfinite(TC) || !std::isfinite(Vo) ||
        !std::isfinite(Tcs) || !std::isfinite(d_gap)) {
        throw py::value_error("TE, TC, Vo, Tcs, and d_gap must be finite.");
    }
    if (TE <= 0.0 || TC <= 0.0 || Tcs <= 0.0 || d_gap <= 0.0) {
        throw py::value_error("TE, TC, Tcs, and d_gap must be positive.");
    }

    std::vector<double> input = { TE, TC, Tcs, d_gap, Vo, -1.0, -1.0, -1.0 };
    thermionicEmission unit(input);
    ThermionicEmissionDiagnosticResult result = unit.calcDiagnostics(true);

    py::dict out;
    out["J"] = result.J;
    out["Vd"] = result.Vd;
    out["delta_V"] = result.delta_V;
    out["phiE"] = result.phiE;
    out["phiC"] = result.phiC;
    out["regime"] = result.regime;
    out["converged"] = result.converged;
    out["finite_flag"] = result.finite;
    out["iteration_count"] = result.iteration_count;
    out["obstructed_iterations"] = result.obstructed_iterations;
    out["transition_iterations"] = result.transition_iterations;
    out["saturation_iterations"] = result.saturation_iterations;
    out["obstructed_residual"] = result.obstructed_residual;
    out["transition_residual"] = result.transition_residual;
    out["saturation_residual"] = result.saturation_residual;
    return out;
}

py::dict calc_emission_point_production(double TE, double TC, double Vo, double Tcs, double d_gap) {
    if (!std::isfinite(TE) || !std::isfinite(TC) || !std::isfinite(Vo) ||
        !std::isfinite(Tcs) || !std::isfinite(d_gap)) {
        throw py::value_error("TE, TC, Vo, Tcs, and d_gap must be finite.");
    }
    if (TE <= 0.0 || TC <= 0.0 || Tcs <= 0.0 || d_gap <= 0.0) {
        throw py::value_error("TE, TC, Tcs, and d_gap must be positive.");
    }
    std::vector<double> input = { TE, TC, Tcs, d_gap, Vo, -1.0, -1.0, -1.0 };
    thermionicEmission unit(input);
    double J = unit.calc();
    py::dict out;
    out["J"] = J;
    out["Vd"] = unit.Vd;
    out["delta_V"] = unit.delta_V;
    out["phiE"] = unit.phiE;
    out["phiC"] = unit.phiC;
    return out;
}

std::vector<double> array_to_vector_1d(const py::array_t<double>& arr, const char* name) {
    if (arr.ndim() != 1) {
        throw py::value_error(string(name) + " must be one-dimensional.");
    }
    auto r = arr.unchecked<1>();
    std::vector<double> values(static_cast<std::size_t>(r.shape(0)));
    for (py::ssize_t i = 0; i < r.shape(0); ++i) {
        values[static_cast<std::size_t>(i)] = r(i);
    }
    return values;
}

std::vector<float> array_to_vector_flat_float64(const py::array_t<double>& arr, const char* name, std::size_t expected) {
    if (static_cast<std::size_t>(arr.size()) != expected) {
        throw py::value_error(string(name) + " size does not match axis product.");
    }
    py::array_t<double, py::array::c_style | py::array::forcecast> c_arr(arr);
    py::buffer_info info = c_arr.request();
    const double* ptr = static_cast<const double*>(info.ptr);
    std::vector<float> values(expected);
    for (std::size_t i = 0; i < expected; ++i) {
        values[i] = static_cast<float>(ptr[i]);
    }
    return values;
}

std::vector<float> array_to_vector_flat_float32(const py::array_t<float>& arr, const char* name, std::size_t expected) {
    if (static_cast<std::size_t>(arr.size()) != expected) {
        throw py::value_error(string(name) + " size does not match axis product.");
    }
    py::array_t<float, py::array::c_style | py::array::forcecast> c_arr(arr);
    py::buffer_info info = c_arr.request();
    const float* ptr = static_cast<const float*>(info.ptr);
    std::vector<float> values(expected);
    for (std::size_t i = 0; i < expected; ++i) {
        values[i] = ptr[i];
    }
    return values;
}

std::vector<uint8_t> array_to_vector_flat_u8(const py::array_t<uint8_t>& arr, const char* name, std::size_t expected) {
    if (static_cast<std::size_t>(arr.size()) != expected) {
        throw py::value_error(string(name) + " size does not match axis product.");
    }
    py::array_t<uint8_t, py::array::c_style | py::array::forcecast> c_arr(arr);
    py::buffer_info info = c_arr.request();
    const uint8_t* ptr = static_cast<const uint8_t*>(info.ptr);
    std::vector<uint8_t> values(expected);
    for (std::size_t i = 0; i < expected; ++i) {
        values[i] = ptr[i];
    }
    return values;
}

void add_emission_lookup_block(
    const std::string& name,
    int priority,
    const py::array_t<double>& TE_axis,
    const py::array_t<double>& TC_axis,
    const py::array_t<double>& Vo_axis,
    const py::array_t<double>& Tcs_axis,
    const py::array_t<double>& J,
    const py::array_t<double>& Vd,
    const py::array_t<double>& delta_V,
    const py::array_t<double>& phiE,
    const py::array_t<double>& phiC,
    const py::array_t<uint8_t>& lookup_safe)
{
    EmissionLookupBlock block;
    block.name = name;
    block.priority = priority;
    block.TE_axis = array_to_vector_1d(TE_axis, "TE_axis");
    block.TC_axis = array_to_vector_1d(TC_axis, "TC_axis");
    block.Vo_axis = array_to_vector_1d(Vo_axis, "Vo_axis");
    block.Tcs_axis = array_to_vector_1d(Tcs_axis, "Tcs_axis");
    const std::size_t n = block.TE_axis.size() * block.TC_axis.size() * block.Vo_axis.size() * block.Tcs_axis.size();
    block.J = array_to_vector_flat_float64(J, "J", n);
    block.Vd = array_to_vector_flat_float64(Vd, "Vd", n);
    block.delta_V = array_to_vector_flat_float64(delta_V, "delta_V", n);
    block.phiE = array_to_vector_flat_float64(phiE, "phiE", n);
    block.phiC = array_to_vector_flat_float64(phiC, "phiC", n);
    block.lookup_safe = array_to_vector_flat_u8(lookup_safe, "lookup_safe", n);
    addEmissionLookupBlock(block);
}

void add_emission_runtime_block(
    const std::string& name,
    int priority,
    int region_id,
    const py::array_t<double>& TE_axis,
    const py::array_t<double>& TC_axis,
    const py::array_t<double>& Vo_axis,
    const py::array_t<double>& Tcs_axis,
    const py::array_t<float>& J,
    const py::array_t<float>& Vd,
    const py::array_t<float>& delta_V,
    const py::array_t<float>& phiE,
    const py::array_t<float>& phiC,
    const py::array_t<uint8_t>& lookup_safe,
    const py::array_t<uint8_t>& zero_mask)
{
    EmissionLookupBlock block;
    block.name = name;
    block.priority = priority;
    block.region_id = region_id;
    block.TE_axis = array_to_vector_1d(TE_axis, "TE_axis");
    block.TC_axis = array_to_vector_1d(TC_axis, "TC_axis");
    block.Vo_axis = array_to_vector_1d(Vo_axis, "Vo_axis");
    block.Tcs_axis = array_to_vector_1d(Tcs_axis, "Tcs_axis");
    const std::size_t n = block.TE_axis.size() * block.TC_axis.size() * block.Vo_axis.size() * block.Tcs_axis.size();
    block.J = array_to_vector_flat_float32(J, "J", n);
    block.Vd = array_to_vector_flat_float32(Vd, "Vd", n);
    block.delta_V = array_to_vector_flat_float32(delta_V, "delta_V", n);
    block.phiE = array_to_vector_flat_float32(phiE, "phiE", n);
    block.phiC = array_to_vector_flat_float32(phiC, "phiC", n);
    block.lookup_safe = array_to_vector_flat_u8(lookup_safe, "lookup_safe", n);
    block.zero_mask = array_to_vector_flat_u8(zero_mask, "zero_mask", n);
    addEmissionLookupBlock(block);
}

void add_emission_dense_region(
    const std::string& name,
    int priority,
    int region_id,
    double d_gap,
    const py::array_t<double>& TE_axis,
    const py::array_t<double>& TC_axis,
    const py::array_t<double>& Vo_axis,
    const py::array_t<double>& Tcs_axis,
    const py::array_t<float>& J,
    const py::array_t<float>& Vd,
    const py::array_t<float>& delta_V,
    const py::array_t<float>& phiE,
    const py::array_t<float>& phiC,
    const py::array_t<uint8_t>& lookup_safe_bits,
    const py::array_t<uint8_t>& zero_mask_bits,
    std::size_t point_count)
{
    DenseEmissionLookupRegion region;
    region.name = name;
    region.priority = priority;
    region.region_id = region_id;
    region.d_gap = d_gap;
    region.TE_axis = array_to_vector_1d(TE_axis, "TE_axis");
    region.TC_axis = array_to_vector_1d(TC_axis, "TC_axis");
    region.Vo_axis = array_to_vector_1d(Vo_axis, "Vo_axis");
    region.Tcs_axis = array_to_vector_1d(Tcs_axis, "Tcs_axis");
    const std::size_t n = region.TE_axis.size() * region.TC_axis.size() * region.Vo_axis.size() * region.Tcs_axis.size();
    if (point_count != n) {
        throw py::value_error("point_count does not match dense axis product.");
    }
    const std::size_t bit_bytes = (n + 7u) / 8u;
    region.point_count = n;
    region.J = array_to_vector_flat_float32(J, "J", n);
    region.Vd = array_to_vector_flat_float32(Vd, "Vd", n);
    region.delta_V = array_to_vector_flat_float32(delta_V, "delta_V", n);
    region.phiE = array_to_vector_flat_float32(phiE, "phiE", n);
    region.phiC = array_to_vector_flat_float32(phiC, "phiC", n);
    region.lookup_safe_bits = array_to_vector_flat_u8(lookup_safe_bits, "lookup_safe_bits", bit_bytes);
    region.zero_mask_bits = array_to_vector_flat_u8(zero_mask_bits, "zero_mask_bits", bit_bytes);
    addEmissionDenseRegion(region);
}

py::dict lookup_emission_point(double TE, double TC, double Vo, double Tcs, double d_gap) {
    EmissionLookupQueryResult result = queryEmissionLookup(TE, TC, Vo, Tcs, d_gap);
    py::dict out;
    out["found"] = result.found;
    out["source"] = result.source;
    out["J"] = result.J;
    out["Vd"] = result.Vd;
    out["delta_V"] = result.delta_V;
    out["phiE"] = result.phiE;
    out["phiC"] = result.phiC;
    return out;
}

py::dict lookup_emission_points(
    const py::array_t<double>& TE,
    const py::array_t<double>& TC,
    const py::array_t<double>& Vo,
    const py::array_t<double>& Tcs,
    double d_gap)
{
    if (TE.ndim() != 1 || TC.ndim() != 1 || Vo.ndim() != 1 || Tcs.ndim() != 1) {
        throw py::value_error("TE, TC, Vo, and Tcs must be one-dimensional arrays.");
    }
    const py::ssize_t n = TE.shape(0);
    if (TC.shape(0) != n || Vo.shape(0) != n || Tcs.shape(0) != n) {
        throw py::value_error("TE, TC, Vo, and Tcs must have the same length.");
    }
    auto te = TE.unchecked<1>();
    auto tc = TC.unchecked<1>();
    auto vo = Vo.unchecked<1>();
    auto tcs = Tcs.unchecked<1>();
    py::array_t<double> J({ n }, { static_cast<py::ssize_t>(sizeof(double)) });
    py::array_t<double> Vd({ n }, { static_cast<py::ssize_t>(sizeof(double)) });
    py::array_t<double> delta_V({ n }, { static_cast<py::ssize_t>(sizeof(double)) });
    py::array_t<double> phiE({ n }, { static_cast<py::ssize_t>(sizeof(double)) });
    py::array_t<double> phiC({ n }, { static_cast<py::ssize_t>(sizeof(double)) });
    py::array_t<uint8_t> found({ n }, { static_cast<py::ssize_t>(sizeof(uint8_t)) });
    auto j_out = J.mutable_unchecked<1>();
    auto vd_out = Vd.mutable_unchecked<1>();
    auto dv_out = delta_V.mutable_unchecked<1>();
    auto pe_out = phiE.mutable_unchecked<1>();
    auto pc_out = phiC.mutable_unchecked<1>();
    auto found_out = found.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < n; ++i) {
        EmissionLookupQueryResult result = queryEmissionLookup(te(i), tc(i), vo(i), tcs(i), d_gap);
        found_out(i) = result.found ? 1 : 0;
        j_out(i) = result.J;
        vd_out(i) = result.Vd;
        dv_out(i) = result.delta_V;
        pe_out(i) = result.phiE;
        pc_out(i) = result.phiC;
    }
    py::dict out;
    out["found"] = found;
    out["J"] = J;
    out["Vd"] = Vd;
    out["delta_V"] = delta_V;
    out["phiE"] = phiE;
    out["phiC"] = phiC;
    return out;
}

// -------------------------------------------------------------------------
// 4. Pybind11 模块定义
// -------------------------------------------------------------------------

PYBIND11_MODULE(te_solver, m) {
    m.doc() = "Thermionic Energy Conversion Solver Binding";

    // 1. 绑定枚举
    py::enum_<CalculationMode>(m, "CalculationMode")
        .value("FixedCurrent", CalculationMode::FixedCurrent)
        .value("FixedVoltage", CalculationMode::FixedVoltage)
        .value("FixedResistance", CalculationMode::FixedResistance)
        .value("ParallelFixedVoltage", CalculationMode::ParallelFixedVoltage)
        .value("ParallelFixedCurrent", CalculationMode::ParallelFixedCurrent)
        .value("ParallelLoadCurve", CalculationMode::ParallelLoadCurve)
        .export_values();

    // 2. 绑定 InputData (只读写属性即可)
    py::class_<InputData>(m, "InputData")
        .def(py::init<>())
        .def_readwrite("N_elements", &InputData::N_elements)
        .def_readwrite("n_axi", &InputData::n_axi)
        .def_readwrite("Temitter", &InputData::Temitter)
        .def_readwrite("Tcollector", &InputData::Tcollector)
        .def_readwrite("dlE", &InputData::dlE)
        .def_readwrite("dlC", &InputData::dlC)
        .def_readwrite("Tcs", &InputData::Tcs)
        .def_readwrite("V_init", &InputData::V_init)
        .def_readwrite("crossAreaE", &InputData::crossAreaE)
        .def_readwrite("crossAreaC", &InputData::crossAreaC)
        .def_readwrite("sideAreaE", &InputData::sideAreaE)
        .def_readwrite("sideAreaC", &InputData::sideAreaC)
        .def_readwrite("U_init", &InputData::U_init)
        .def_readwrite("d_gap", &InputData::d_gap)
        .def_readwrite("Itarget", &InputData::Itarget)
        .def_readwrite("resistanceWire", &InputData::resistanceWire)
        .def_readwrite("wireU", &InputData::wireU)
        .def_readwrite("mode", &InputData::mode)
        .def_readwrite("target_val", &InputData::target_val)
        .def_readwrite("I_total_init", &InputData::I_total_init)
        .def_readwrite("loadCurveCurrent", &InputData::loadCurveCurrent)
        .def_readwrite("loadCurveVoltage", &InputData::loadCurveVoltage);

    // 3. 绑定 singleThermionicEnergyConversion 类
    // 重点：暴露核心物理向量，以便 Python 端可以直接修改
    py::class_<singleThermionicEnergyConversion>(m, "SingleTEC")
        // 不需要暴露构造函数，因为我们用 create_circuit
        // 关键物理场 (读写)
        .def_readwrite("Temitter", &singleThermionicEnergyConversion::Temitter)
        .def_readwrite("Tcollector", &singleThermionicEnergyConversion::Tcollector)
        .def_readwrite("Tcs", &singleThermionicEnergyConversion::Tcs)
        // 结果场 (读写，通常是只读，但为了调试方便给读写)
        .def_readwrite("J", &singleThermionicEnergyConversion::J)
        .def_readwrite("V", &singleThermionicEnergyConversion::V)
        .def_readwrite("UE", &singleThermionicEnergyConversion::UE)
        .def_readwrite("UC", &singleThermionicEnergyConversion::UC)
        .def_readwrite("rhoE", &singleThermionicEnergyConversion::rhoE)
        .def_readwrite("rhoC", &singleThermionicEnergyConversion::rhoC)
        .def_readwrite("phiE", &singleThermionicEnergyConversion::phiE)
        .def_readwrite("phiC", &singleThermionicEnergyConversion::phiC)
        .def_readwrite("Vd", &singleThermionicEnergyConversion::Vd)
        .def_readwrite("joulePowerE", &singleThermionicEnergyConversion::joulePowerE)
        .def_readwrite("joulePowerC", &singleThermionicEnergyConversion::joulePowerC)
        .def_readwrite("terminalPointUE1", &singleThermionicEnergyConversion::terminalPointUE1)
        .def_readwrite("terminalPointUE2", &singleThermionicEnergyConversion::terminalPointUE2)
        .def_readwrite("terminalPointUC1", &singleThermionicEnergyConversion::terminalPointUC1)
        .def_readwrite("terminalPointUC2", &singleThermionicEnergyConversion::terminalPointUC2)
        // >>>>> 新增：内部截面电流分布 <<<<<
        .def_readwrite("IEsecSingle", &singleThermionicEnergyConversion::IEsecSingle)
        .def_readwrite("ICsecSingle", &singleThermionicEnergyConversion::ICsecSingle)
        // 标量结果
        .def_readwrite("I", &singleThermionicEnergyConversion::I)
        .def_readwrite("U", &singleThermionicEnergyConversion::U)
        // 计算方法 (如果需要单独调试单个元件)
        .def("calc_current", &singleThermionicEnergyConversion::Icalc);

    // 4. 绑定 circuitTECs 类
    py::class_<circuitTECs>(m, "CircuitTECs")
        // 暴露 TECs 列表
        // Pybind11 会自动将其转换为 Python list (包含 SingleTEC 对象的引用)
        .def_readwrite("TECs", &circuitTECs::TECs)
        // 全局控制参数
        .def_readwrite("Utarget", &circuitTECs::Utarget)
        .def_readwrite("Rload", &circuitTECs::Rload)
        .def_readwrite("Iout", &circuitTECs::Iout)
        .def_readwrite("Uout", &circuitTECs::Uout)
        .def_readwrite("Itarget", &circuitTECs::Itarget)
        .def_readwrite("isFixedU", &circuitTECs::isFixedU)
        .def_readwrite("isFixedI", &circuitTECs::isFixedI)
        .def_readwrite("isFixedR", &circuitTECs::isFixedR)
        .def_readwrite("isParallelFixedU", &circuitTECs::isParallelFixedU)
        .def_readwrite("isParallelFixedI", &circuitTECs::isParallelFixedI)
        .def_readwrite("isParallelLoadCurve", &circuitTECs::isParallelLoadCurve)
        .def_readwrite("converged", &circuitTECs::converged)
        .def_readwrite("iterationCount", &circuitTECs::iterationCount)
        .def_readwrite("branchCurrents", &circuitTECs::branchCurrents)
        .def_readwrite("branchVoltages", &circuitTECs::branchVoltages)
        .def("set_tcs", [](circuitTECs& circuit, const py::array_t<double>& values) {
            int n_elements = static_cast<int>(circuit.TECs.size());
            int n_axi = n_elements == 0 ? 0 : static_cast<int>(circuit.TECs[0]->Tcs.size());
            require_2d("Tcs", values, n_elements, n_axi);
            vector<vector<double>> rows;
            rows.reserve(n_elements);
            for (int i = 0; i < n_elements; ++i) {
                rows.push_back(get_row_vector(values, i));
            }
            circuit.setTcs(rows);
        })
        .def("set_load_curve", [](circuitTECs& circuit, const py::array_t<double>& current, const py::array_t<double>& voltage) {
            if (current.ndim() != 1 || voltage.ndim() != 1 || current.shape(0) != voltage.shape(0)) {
                throw py::value_error("current and voltage must be one-dimensional arrays with the same length.");
            }
            circuit.setLoadCurve(get_vector_1d(current), get_vector_1d(voltage));
        })
        // 核心计算函数
        .def("calc", &circuitTECs::circuitTECsCalc);

    // 5. 绑定工厂函数
    m.def("create_circuit", &create_circuit, "Create and initialize the circuit from InputData");
    m.def("clear_emission_lookup", &clearEmissionLookup);
    m.def("set_emission_lookup_enabled", &setEmissionLookupEnabled, py::arg("enabled"));
    m.def("is_emission_lookup_enabled", &isEmissionLookupEnabled);
    m.def("emission_lookup_block_count", &emissionLookupBlockCount);
    m.def("emission_lookup_region_count", &emissionLookupRegionCount);
    m.def("emission_lookup_dense_region_count", &emissionLookupDenseRegionCount);
    m.def(
        "add_emission_lookup_block",
        &add_emission_lookup_block,
        py::arg("name"),
        py::arg("priority"),
        py::arg("TE_axis"),
        py::arg("TC_axis"),
        py::arg("Vo_axis"),
        py::arg("Tcs_axis"),
        py::arg("J"),
        py::arg("Vd"),
        py::arg("delta_V"),
        py::arg("phiE"),
        py::arg("phiC"),
        py::arg("lookup_safe")
    );
    m.def(
        "add_emission_runtime_block",
        &add_emission_runtime_block,
        py::arg("name"),
        py::arg("priority"),
        py::arg("region_id"),
        py::arg("TE_axis"),
        py::arg("TC_axis"),
        py::arg("Vo_axis"),
        py::arg("Tcs_axis"),
        py::arg("J"),
        py::arg("Vd"),
        py::arg("delta_V"),
        py::arg("phiE"),
        py::arg("phiC"),
        py::arg("lookup_safe"),
        py::arg("zero_mask")
    );
    m.def(
        "add_emission_dense_region",
        &add_emission_dense_region,
        py::arg("name"),
        py::arg("priority"),
        py::arg("region_id"),
        py::arg("d_gap"),
        py::arg("TE_axis"),
        py::arg("TC_axis"),
        py::arg("Vo_axis"),
        py::arg("Tcs_axis"),
        py::arg("J"),
        py::arg("Vd"),
        py::arg("delta_V"),
        py::arg("phiE"),
        py::arg("phiC"),
        py::arg("lookup_safe_bits"),
        py::arg("zero_mask_bits"),
        py::arg("point_count")
    );
    m.def(
        "load_emission_dense_file",
        &loadEmissionDenseFile,
        py::arg("path")
    );
    m.def(
        "lookup_emission_point",
        &lookup_emission_point,
        py::arg("TE"),
        py::arg("TC"),
        py::arg("Vo"),
        py::arg("Tcs"),
        py::arg("d_gap") = 0.5
    );
    m.def(
        "lookup_emission_points",
        &lookup_emission_points,
        py::arg("TE"),
        py::arg("TC"),
        py::arg("Vo"),
        py::arg("Tcs"),
        py::arg("d_gap") = 0.5
    );
    m.def(
        "calc_emission_point",
        &calc_emission_point,
        py::arg("TE"),
        py::arg("TC"),
        py::arg("Vo"),
        py::arg("Tcs"),
        py::arg("d_gap") = 0.5,
        "Evaluate one local thermionic-emission point with diagnostic metadata."
    );
    m.def(
        "calc_emission_point_production",
        &calc_emission_point_production,
        py::arg("TE"),
        py::arg("TC"),
        py::arg("Vo"),
        py::arg("Tcs"),
        py::arg("d_gap") = 0.5,
        "Evaluate one local thermionic-emission point through production calc()."
    );
}

