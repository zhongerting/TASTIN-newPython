# TASTIN-Python

## Project Overview

TASTIN-Python is a spacecraft nuclear thermal-hydraulic coupled system simulation platform based on Python. It simulates the transient thermal-hydraulic behavior of space reactors. It is a Python refactoring of the TASTIN program, designed with modularity and object-oriented principles, supporting multi-physics coupled calculations. 

**Key Technologies:**
- **Language:** Python 3.8+, with a C++ extension for specific calculations.
- **Dependencies:** NumPy, SciPy, Matplotlib (for visualization).
- **C++ Extension:** Uses `pybind11` (automatically fetched via CMake) for the Thermionic Energy Conversion (TEC) module (`ThermoCalc`).

**Architecture:**
- `Components/`: Macro-level components like `ReactorCore`, `TFEUnit`, `RingHP` (Heat Pipes), and `ExternalHeatSources`.
- `Solvers/`: Physics solvers for Heat Conduction (2D transient), Hydrodynamics (hydraulic network), Neutronics (point reactor kinetics), and a `SystemManager` for coordinating multi-physics fields.
- `Materials/`: Material property definitions for fluids (e.g., Sodium, NaK78) and solids (e.g., UO2, StainlessSteel).
- `Correlations/`: Engineering correlations for heat transfer, pressure drop, etc.
- `ThermoCalc/`: A C++ module for Thermionic Energy Conversion calculation.
- `testModule/`: Contains numerous test cases and verification scripts.
- `*优化重构/`: Various directories (e.g., `HeatConduction优化重构`, `SystemManager优化重构`) indicating active refactoring and performance optimization efforts, containing TODO lists, benchmark scripts, and optimization logs.

## Building and Running

### Prerequisites
- Python 3.8 or higher.
- `pip install numpy scipy matplotlib`

### Building the C++ Extension (`ThermoCalc`)
The project uses CMake and `pybind11` to build the `te_solver` C++ module.
```bash
# Navigate to ThermoCalc directory
cd ThermoCalc

# Create a build directory and compile
mkdir build
cd build
cmake ..
cmake --build . --config Release
# Copy the generated shared library (.so or .pyd) to the appropriate Python import path (usually back to the root or ThermoCalc dir).
```

### Running Tests
The project contains numerous test cases in the `testModule/` directory and benchmark scripts in the refactoring directories. These are implemented as standard Python scripts rather than using a framework like `pytest`.
```bash
# Run a specific test case from the root directory
python testModule/test_case_3_1.py
```
Test results often include matplotlib plots saved as PNG files (e.g., `Case_3_1_Result.png`) and logging outputs directly to the console.

## Development Conventions
- **Modular Design:** The project heavily uses Object-Oriented Programming (OOP) with abstract base classes to allow a plug-and-play architecture for components.
- **Coupling & Solvers:** Physics solvers are separated from the components. `Coupler` classes (e.g., `SolidSolidCouple2D`, `FluidSolidCouple`) are used to manage data exchange and synchronization between different physical solver objects. The `SystemManager` acts as the unified orchestrator for multi-physics coordination.
- **Logging:** Python's standard `logging` module is used for progress and result output (e.g., `logger.info()`, `logger.error()`).
- **Data Handling:** `numpy` arrays are extensively used for state vectors, meshes, and calculations. The codebase relies heavily on ODE solvers like `scipy.integrate.solve_ivp` (e.g. via `NuclearODESolver` wrapping BDF method).
- **Optimization and Refactoring:** The project is actively undergoing performance optimization. Changes, especially to core solvers, must follow a strict "one idea per round" rule, be verified against frozen benchmark scripts, and pass smoke/stability tests before integration to ensure physical calculation accuracy is not compromised for speed. Optimization logs and TODOs are maintained in their respective `*优化重构/` directories.
