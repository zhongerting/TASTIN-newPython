# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

TASTIN-Python is a multiphysics transient simulator for a thermionic space nuclear-power reactor. Python assembles systems and runs the transient thermal-hydraulic, solid heat-conduction, neutronics, heat-pipe/radiator, and TEC coupling models. The TEC electrical backend is a C++17/pybind11 extension in `ThermoCalc/`.

There is no production entry point in the repository root: `main.py` is empty. Run actual cases from `testModule/` or `CoolantLoop/`.

## Repository guidance and reading order

- Read `AGENTS.md` first. It is the project-wide handoff and module-routing guide; the current source takes precedence over historical notes.
- For the affected subsystem, read its first-entry document before opening broad sets of source files:
  - `Components/COMPONENTS_DETAILED_INTRO.md`
  - `Solvers/AI_AGENT_SOLVERS_ANALYSIS.md`
  - `Materials/AI_AGENT_MATERIALS_ANALYSIS.md`
  - `Correlations/AI_AGENT_CORRELATIONS_ANALYSIS.md`
  - `ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`
  - `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`
  - `CoolantLoop/AI_AGENT_COOLANTLOOP_ANALYSIS.md`
- `ARCHITECTURE.md` contains the longer dependency map. `README.md` contains the current top-level navigation and smoke examples.
- When a public interface, physical formula, unit convention, timestep order, restart format, grid mapping, or cross-module contract changes, update the relevant module handoff document as well as the code.

## Python environment and commands

Use the repository's Conda Python 3.12 environment for normal execution, especially anything importing ThermoCalc:

```powershell
$Py = 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe'
& $Py --version
```

Run commands from the repository root. There is no project-level packaging configuration or unified install command. The code imports directly from the checkout.

### Syntax checks

Compile only the files affected by a change first:

```powershell
& $Py -m py_compile Components\ReactorCore.py Solvers\SystemManager.py
```

### Unit tests

Tests are primarily standard-library `unittest` modules and executable verification scripts; there is no repository-level pytest configuration. Run a module, a class, or one method as follows:

```powershell
# One test module
& $Py -m unittest testModule.test_system_manager_lifecycle

# One test class
& $Py -m unittest testModule.test_system_manager_lifecycle.SystemManagerLifecycleTests

# One test method
& $Py -m unittest testModule.test_system_manager_lifecycle.SystemManagerLifecycleTests.test_restart_round_trip
```

Use the impact map in `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md` to choose tests. Do not start full discovery or long/overnight cases as a first check; many scripts construct sizeable models, write artifacts, or perform long transient runs.

### Common smoke commands

```powershell
# V9 topology and no-TEC open-loop smoke
& $Py -m py_compile testModule\test_core_assemble_v9_caseA.py testModule\run_v9_caseA_open_loop.py testModule\test_v9_caseA_topology.py
& $Py -m unittest testModule.test_v9_caseA_topology
& $Py testModule\run_v9_caseA_open_loop.py --duration 0.1 --record-interval 0.1 --restart-interval 0.1 --max-dt 0.01 --disable-tec-coupled

# Focused ThermoCalc checks (requires the matching .pyd)
& $Py testModule\test_thermocalc_interface.py
& $Py testModule\test_thermocalc_parallel.py
& $Py testModule\test_thermocalc_lookup.py
```

For a requested full-loop or long transient, read the corresponding runner README and verify the input restart, output directory, timestep, and runtime cost before launching it.

## ThermoCalc C++ extension

`ThermoCalc/CMakeLists.txt` builds the `te_solver` module with C++17 and fetches pybind11 through CMake `FetchContent`. Configure and build an ABI-compatible extension with:

```powershell
cmake -S ThermoCalc -B ThermoCalc\build_cp312
cmake --build ThermoCalc\build_cp312 --config Release
```

The build may need network access to fetch pybind11. The resulting `.pyd` is normally under `ThermoCalc\build_cp312\Release`. Test a rebuilt extension without replacing the checkout copy:

```powershell
$env:THERMOCALC_PYD_DIR = (Resolve-Path ThermoCalc\build_cp312\Release).Path
& $Py testModule\test_thermocalc_interface.py
Remove-Item Env:THERMOCALC_PYD_DIR
```

`THERMOCALC_PYD_DIR` is inserted at the front of `sys.path` by `ThermoCalc/ThermoCalcWrapper.py`. Do not assume that the default `python` command has the ABI required by `te_solver.cp312-win_amd64.pyd`; check the interpreter first.

Optional emission lookup testing is controlled explicitly with `THERMOCALC_ENABLE_LOOKUP`, `THERMOCALC_LOOKUP_DB`, and `THERMOCALC_LOOKUP_REGIONS`, or with the corresponding explicit `ReactorCore`/`ThermoCalcModel` arguments. Generated emission databases and build products are runtime data, not source-of-truth files.

## Architecture

The main dependency direction is:

```text
case runners (testModule/, CoolantLoop/)
  -> macro components (Components/)
  -> global orchestration and physics solvers (Solvers/)
  -> materials, correlations, and numerical helpers

ReactorCore / TECCircuitManager
  -> ThermoCalcWrapper.py
  -> te_solver (pybind11)
  -> circuitTECs -> singleThermionicEnergyConversion -> thermionicEmission
```

### Solvers and timestep lifecycle

`Solvers/SystemManager.py` is the global coordinator. Components expose their registered solid solvers and couplers through `get_solids()` and `get_couplers()`, and may implement `pre_step()`/`post_step()` plus state save/load hooks. A normal global step coordinates component preprocessing, Picard coupling iterations, solid heat conduction, hydraulic-network advancement, optional point-reactor advancement/commit, convergence, and postprocessing. `Solvers/Couplers.py` exchanges boundary temperatures, resistances, and heat/source terms between fluid, solid, gap, and TEC domains.

The solver layer contains:

- `Solvers/Hydrodynamics/`: control volumes, channels, junctions/pumps, pressure/enthalpy network assembly, and sparse hydraulic solves.
- `Solvers/HeatConduction/`: 1D/2D finite-volume meshes, boundary conditions, radiation, and solid ODE/integrator paths.
- `Solvers/Neutronics/PointReactor.py`: point kinetics and decay heat.
- `Solvers/Couplers.py`: fluid-solid, solid-solid, gap, active-gap, and TEC coupling.

### Macro components

- `ReactorCore` assembles representative TFE units, global moderator/structural regions, point-reactor feedback, power/multiplier mappings, and one or more TEC circuits.
- `TFEUnit` assembles fuel, electrodes, cladding, moderator, gaps, coolant coupling, neutronic power, TEC heat flux, and electrical/Joule-heat state.
- `RingHP` assembles a collector-ring fluid/header solid with representative `HPwithFin` heat pipes and optional external orbital heat sources.
- `HPwithFin` uses a resolved `HeatPipe2D` plus a reduced-order quasi-steady fin branch; it does not expose the fin as an independent solid solver.
- `RadiatorPipeWithFin` is a separate explicit TOPAZ-II NaK tube-and-copper-fin radiator model; do not confuse it with the heat-pipe `HPwithFin` path.

The principal case families are versioned and topology-specific. V7/V8/V9/V10/V11/V12/V13 runners in `testModule/` and the collector-ring/TOPAZ-II runners in `CoolantLoop/` are not interchangeable by filename alone.

## High-value contracts and pitfalls

- Temperatures are K. Discrete node/source powers are W. External heat fluxes and TEC plasma surface fluxes are W/m²; `ExternalHeatFluxBC` performs the area multiplication, so callers must not multiply those fluxes by area again.
- `ReactorCore` and `RingHP` representative multipliers affect power, heat sources, flow, heat rejection, moderator mapping, and diagnostics. Check the multiplier and physical-vs-representative flow convention before changing topology or statistics.
- `HydraulicNetwork` topology and its pressure-boundary/reference set are fixed after construction. Rebuild the network when changing nodes, junctions, or pressure references. A passive pressure reference is not the same as a fixed-pressure thermal boundary.
- Restart files are topology-, geometry-, material-, multiplier-, and solver-configuration-specific. A successful `.npz` load does not prove physical compatibility. V8/V9/V10/V11 migrations must use their documented migration/injection runners rather than loading an arbitrary older restart directly.
- After `load_global_state()`, check the documented time synchronization and boundary/coupler refresh sequence before reading time-dependent pumps, sources, or diagnostics.
- ThermoCalc electrical results require the matching Python ABI. Production Joule heating uses C++ `joulePowerE/joulePowerC` mapped by `Components/tec_electric.py`; do not recreate production Joule heat from diagnostic electric-field gradients.
- `Components/basicComponents/Electord.py` is intentionally misspelled and existing imports depend on that filename; do not rename it casually.
- The main `SodiumPotassium78` implementation is in `Materials/Fluids/SodiumPotassium78.py`; do not infer the active material solely from the similarly named `NaK78.py` file.
- A zero TEC multiplier means the representative TFE must keep passive gap behavior while its active TEC sources and electrical diagnostics are cleared by the component contract.

## Outputs and repository hygiene

Simulation outputs such as `.npz`, `.csv`, `.png`, `.log`, `.err`, `.out`, restart directories, emission databases, and CMake build directories are ignored or treated as local artifacts. Keep source/configuration/documentation changes separate from generated results unless the user explicitly requests the artifacts. Before adding a new generated file, check `.gitignore` and the target directory.

There is no configured repository-wide formatter or linter. Match the surrounding Python/C++ style and use targeted `py_compile`, `unittest`, topology smoke, energy audit, and (when relevant) ABI-matched ThermoCalc checks as the verification ladder.
