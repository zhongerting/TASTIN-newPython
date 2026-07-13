# ThermoCalc Series Fixed-Current Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a production series `fixed_i` TEC circuit mode with open-circuit fallback for non-generating states.

**Architecture:** Add one C++ mode that reuses `circuitCalc(I)`. A requested current is accepted only when the series solve converges with finite positive terminal voltage; otherwise the circuit is recomputed at zero current and reported as unconverged. Python only maps the public mode and reports existing diagnostics.

**Tech Stack:** C++17, pybind11, Python 3.12, NumPy, unittest-style assertions.

## Global Constraints

- Do not overwrite, rename, or delete `ThermoCalc/te_solver.cp312-win_amd64.pyd`.
- Build the test extension in a new independent directory and load it with `THERMOCALC_PYD_DIR`.
- Use `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe` for Python commands.
- Preserve every existing circuit mode and default.

---

### Task 1: Define failing public-interface tests

**Files:**
- Create: `testModule/test_thermocalc_series_fixed_current.py`

**Interfaces:**
- Consumes: `ThermoCalcModel.setup_circuit_mode(mode_str, target_value, I_guess)`
- Produces: executable checks for zero-current open circuit, feasible fixed current, and rejected current fallback

- [ ] Write a test helper that builds one or two hot series TECs.
- [ ] Assert `setup_circuit_mode("fixed_i", ...)` is accepted.
- [ ] Assert zero target current returns finite electrical and thermal arrays.
- [ ] Assert a feasible target reports the requested current and positive voltage.
- [ ] Assert an excessive target returns zero current and `converged == false`.
- [ ] Run against the current production extension and confirm failure because `FixedCurrent` is absent.

### Task 2: Add the minimal fixed-current mode

**Files:**
- Modify: `ThermoCalc/circuitTECs.h`
- Modify: `ThermoCalc/circuitTECs.cpp`
- Modify: `ThermoCalc/bindings.cpp`
- Modify: `ThermoCalc/ThermoCalcWrapper.py`

**Interfaces:**
- Produces: `CalculationMode.FixedCurrent`, `circuitTECs::isFixedI`, and `circuitTECs::iFixedCircuitCalc()`

- [ ] Add and initialize `isFixedI`.
- [ ] Implement `iFixedCircuitCalc()` by calling `circuitCalc(Itarget)`.
- [ ] Accept only a converged, finite, positive generated voltage.
- [ ] On rejection call `circuitCalc(0.0)`; preserve finite open-circuit voltage, otherwise set voltage to zero.
- [ ] Always return zero current and unconverged status after rejection.
- [ ] Add the enum and factory mapping in pybind11.
- [ ] Map Python `fixed_i` to the new enum and initialize the C++ mode flags.
- [ ] Keep the low-temperature zero-emission guard behavior finite.

### Task 3: Isolated build, regression, and documentation

**Files:**
- Modify: `ThermoCalc/AI_AGENT_THERMOCALC_ANALYSIS.md`

- [ ] Record the production pyd SHA-256 before building.
- [ ] Configure and compile into a new `ThermoCalc/build_series_fixed_i_test` directory.
- [ ] Run the new tests using `THERMOCALC_PYD_DIR` pointing at that build.
- [ ] Run `test_thermocalc_interface.py` and `test_thermocalc_parallel.py` against the same test extension.
- [ ] Run Python syntax checks.
- [ ] Record the new mode and fallback semantics in the ThermoCalc module guide.
- [ ] Verify the production pyd SHA-256 is unchanged.