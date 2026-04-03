cmake_minimum_required(VERSION 3.10)

# 项目名称（建议与 Python 模块名保持一致）
project(te_solver)

# 设置 C++ 标准 (推荐 C++14 或 C++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------------------------------------
# 1. 自动下载并配置 pybind11 (使用 FetchContent)
#    这样您不需要手动安装 pybind11，CMake 会联网自动搞定
# -------------------------------------------------------------
include(FetchContent)
message(STATUS "Fetching pybind11...")
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://gitee.com/mirrors/pybind11.git
  GIT_TAG        v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# -------------------------------------------------------------
# 2. 定义 Python 扩展模块
# -------------------------------------------------------------
# pybind11_add_module 会自动处理文件后缀 (.pyd on Windows, .so on Linux)
# 这里的 "te_solver" 必须与 bindings.cpp 中 PYBIND11_MODULE(te_solver, m) 的名字一致
pybind11_add_module(te_solver
    bindings.cpp
    circuitTECs.cpp
    singleThermionicEnergyConversion.cpp
    thermionicEmission.cpp
    NonLinerSolver.cpp
)

# -------------------------------------------------------------
# 3. 编译选项优化 (Windows)
# -------------------------------------------------------------
if(MSVC)
    # 添加 /utf-8 选项解决中文乱码警告
    target_compile_options(te_solver PRIVATE /O2 /fp:fast /utf-8)
endif()

message(STATUS "Build setup complete for te_solver.")