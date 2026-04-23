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
mkdir build && cd build
cmake ..
cmake --build .
# TODO: Copy the generated shared library (.so or .pyd) to the appropriate Python import path (usually back to the root or ThermoCalc dir).
```

### Running Tests
The project contains numerous test cases in the `testModule/` directory. These are implemented as standard Python scripts rather than using a framework like `pytest`.
```bash
# Run a specific test case from the root directory
python testModule/test_case_3_1.py
```
Test results often include matplotlib plots saved as PNG files (e.g., `Case_3_1_Result.png`) and logging outputs directly to the console.

## Development Conventions
- **Modular Design:** The project heavily uses Object-Oriented Programming (OOP) with abstract base classes to allow a plug-and-play architecture for components.
- **Coupling & Solvers:** Physics solvers are separated from the components. For example, `CoupledModelWrapper` and `SolidSolidCouple2D` are used to manage the data exchange between different solver objects.
- **Logging:** Python's standard `logging` module is used for progress and result output (e.g., `logger.info()`, `logger.error()`).
- **Data Handling:** `numpy` arrays are extensively used for state vectors, meshes, and calculations.
