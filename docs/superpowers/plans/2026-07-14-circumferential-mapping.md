# Circumferential Mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reusable, conservative periodic circumferential mapping utility and document its V14 18-zone to shield 12-zone use.

**Architecture:** A small NumPy utility builds a column-conservative source-allocation matrix from exact angular overlaps. A second function maps intensive values using optional physical source weights. The current shield implementation remains unchanged.

**Tech Stack:** Python 3.12, NumPy, unittest.

## Global Constraints

- Use `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`.
- Add no dependency.
- Map `T^4` supplied by the caller; do not transform temperature internally.
- Preserve extensive totals through matrix column conservation.
- Do not connect this utility to `RadiatorThermalShield` in this task.

---

### Task 1: Periodic Circumferential Mapping Utility

**Files:**
- Create: `Components/circumferential_mapping.py`
- Create: `testModule/test_circumferential_mapping.py`

**Interfaces:**
- Produces: `build_uniform_circumferential_mapping(n_source, n_target, source_offset_deg=0.0, target_offset_deg=0.0) -> np.ndarray`
- Produces: `map_circumferential_intensive(source_values, mapping, source_weights=None) -> np.ndarray`

- [ ] **Step 1: Write failing tests**

Add unittest cases that assert:

```python
mapping = build_uniform_circumferential_mapping(18, 12)
np.testing.assert_allclose(mapping.sum(axis=0), 1.0)
np.testing.assert_allclose(mapping[0, :3], [1.0, 0.5, 0.0])
np.testing.assert_allclose(mapping[1, :3], [0.0, 0.5, 1.0])
```

Also test offset wrap, extensive conservation using `mapping @ totals`, weighted V14 aggregation with `[5, 6, 6]`, constant-field preservation, batched final-axis input, and invalid counts/shapes/non-finite/negative/zero-effective-weight inputs.

- [ ] **Step 2: Run tests and confirm import failure**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_circumferential_mapping
```

Expected: FAIL because `Components.circumferential_mapping` does not exist.

- [ ] **Step 3: Implement the minimum utility**

Use exact wrapped interval intersections. Matrix shape is `(n_target, n_source)`; each entry is overlap divided by source width. Validate finite offsets and positive integer counts.

For intensive mapping:

```python
weighted = source_values * source_weights
numerator = weighted @ mapping.T
denominator = source_weights @ mapping.T
return numerator / denominator
```

Broadcast weights across leading value dimensions and reject zero target denominators.

- [ ] **Step 4: Run focused tests**

Run the Task 1 unittest command. Expected: all tests pass.

- [ ] **Step 5: Run neighboring shield tests**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.test_radiator_thermal_shield testModule.test_v14_caseA_topology testModule.test_v15_caseA_topology
```

Expected: pass with no shield behavior change.

### Task 2: Reuse Documentation

**Files:**
- Create: `Components/CIRCUMFERENTIAL_MAPPING_GUIDE.md`
- Modify: `Components/COMPONENTS_DETAILED_INTRO.md`

**Interfaces:**
- Documents both Task 1 functions and the V14 `hp_multipliers` call pattern.

- [ ] **Step 1: Write the guide**

Document `F[j,i]=overlap/source_width`, extensive and intensive formulas, the aligned 18-to-12 pattern, V14 `5/6` multiplier example, offsets, validation, and deferred shield integration.

- [ ] **Step 2: Add one navigation link**

Add `CIRCUMFERENTIAL_MAPPING_GUIDE.md` to the Components documentation navigation near the thermal-shield/radiator material.

- [ ] **Step 3: Verify documentation and source**

Run:

```powershell
rg -n "circumferential_mapping|CIRCUMFERENTIAL_MAPPING_GUIDE|build_uniform_circumferential_mapping" Components
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile Components\circumferential_mapping.py testModule\test_circumferential_mapping.py
```

Expected: references are present and compilation succeeds.

- [ ] **Step 4: Final regression**

Run Task 1 focused tests and neighboring shield tests again. Expected: all pass.
