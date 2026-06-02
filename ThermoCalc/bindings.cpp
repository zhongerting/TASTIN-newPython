#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace std;

#include "circuitTECs.h"
#include "singleThermionicEnergyConversion.h"

// -------------------------------------------------------------------------
// 1. 枚举与辅助结构体定义
// -------------------------------------------------------------------------

// 电路计算模式枚举
enum class CalculationMode {
    FixedVoltage, // 固定电压模式
    FixedResistance // 固定电阻模式
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
    if (data.mode == CalculationMode::FixedVoltage) {
        circuit->isFixedU = true;
        circuit->isFixedR = false;
        circuit->Utarget = data.target_val;
    } 
    else {
        circuit->isFixedU = false;
        circuit->isFixedR = true;
        circuit->Rload = data.target_val;
    }

    return circuit;
}

// -------------------------------------------------------------------------
// 4. Pybind11 模块定义
// -------------------------------------------------------------------------

PYBIND11_MODULE(te_solver, m) {
    m.doc() = "Thermionic Energy Conversion Solver Binding";

    // 1. 绑定枚举
    py::enum_<CalculationMode>(m, "CalculationMode")
        .value("FixedVoltage", CalculationMode::FixedVoltage)
        .value("FixedResistance", CalculationMode::FixedResistance)
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
        .def_readwrite("I_total_init", &InputData::I_total_init);

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
        .def_readwrite("isFixedU", &circuitTECs::isFixedU)
        .def_readwrite("isFixedR", &circuitTECs::isFixedR)
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
        // 核心计算函数
        .def("calc", &circuitTECs::circuitTECsCalc);

    // 5. 绑定工厂函数
    m.def("create_circuit", &create_circuit, "Create and initialize the circuit from InputData");
}

