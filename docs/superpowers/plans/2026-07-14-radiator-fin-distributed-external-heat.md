# Radiator Fin Distributed External Heat Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in distributed fin absorption to V15 without changing default results.

**Architecture:** Keep `lumped_root_area` as default. Reuse the `HPwithFin` equation: wall absorption enters one boundary BC; fin absorption enters the quasi-steady fin RHS.

**Tech Stack:** Python 3.12, NumPy, existing external-heat and boundary classes, unittest.

## Global Constraints

- Use `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`.
- Keep all existing constructor calls valid.
- Use single-sided fin projected area for absorption and double-sided area for radiation.
- Apply each tube's N78 density uniformly over axial slices, each with its own area.
- Preserve unrelated worktree changes.

---

## Task 1: Add Focused Failing Tests

**Status:** completed

**Files:**
- Modify: `testModule/test_v15_caseA_topology.py`
- Create: `testModule/test_radiator_pipe_fin_external_heat.py`

- [ ] Build one reduced V15 pipe-fin unit in both modes.
- [ ] Assert the default remains `lumped_root_area` and its wall-BC area is unchanged.
- [ ] Assert the legacy mode keeps its old effective area and the new mode uses geometric single-sided fin area.
- [ ] Assert distributed fin absorption equals flux times single-sided projected fin area per axial slice.
- [ ] Run the focused tests and record the expected pre-implementation failure.

## Task 2: Add Fin Absorption to RadiatorPipeWithFin

**Status:** completed

**Files:**
- Modify: `Components/RadiatorPipeWithFin.py`
- Test: `testModule/test_radiator_pipe_fin_external_heat.py`

- [ ] Add `set_fin_external_heat_source(source, illuminated_area_scale=1.0)`.
- [ ] Apply one tube heat-flux density uniformly to all axial fin slices.
- [ ] Add `q_flux * fin_strip_width * dx * scale` to each quasi-steady fin cell.
- [ ] Store `last_fin_absorption_distribution` and preserve `Q_net = Q_rad - Q_abs`.
- [ ] Add tube/fin/total external-heat accounting diagnostics.
- [ ] Run the focused test and confirm the fin energy identity closes.

## Task 3: Wire the Opt-In V15 Mode

**Status:** completed

**Files:**
- Modify: `testModule/Full_Loop_Cases/v15_pipefin_radiator.py`
- Modify: `testModule/test_v15_caseA_topology.py`

- [ ] Add `external_heat_fin_loading_mode = lumped_root_area` to the V15 config.
- [ ] Reject mode names other than `lumped_root_area` and `distributed_fin`.
- [ ] Preserve the legacy lumped wall BC exactly.
- [ ] In distributed mode, put only tube projected area on the wall BC and route fin loading through `set_fin_external_heat_source`.
- [ ] Configure accounting in both modes while preserving `external_heat_source` and `external_heat_bc` aliases.
- [ ] Run focused and V15 topology tests.

## Task 4: Regression, Smoke Test, and Documentation

**Status:** completed

**Files:**
- Modify: `Components/RADIATORPIPEWITHFIN_DETAILED_INTRO.md`
- Modify: `Components/EXTERNALHEATSOURCES_DETAILED_INTRO.md`
- Modify: `testModule/Full_Loop_Cases/AI_AGENT_FULL_LOOP_CASES_ANALYSIS.md`
- Temporary: `E:\tmp\v15_distributed_external_heat_smoke.py`

- [ ] Document mode semantics, areas, units, defaults, and the energy identity.
- [ ] Run a temporary short V15 transient with distributed mode and verify finite temperatures and heat rates.
- [ ] Run `py_compile` for changed Python files.
- [ ] Run the focused external-heat, V14 topology, V15 topology, and matrix-source test modules.
- [ ] Review the final diff for default compatibility and unrelated changes.
