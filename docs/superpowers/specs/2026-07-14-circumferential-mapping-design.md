# Circumferential 18-to-12 Mapping Design

Date: 2026-07-14

## Goal

Add a reusable Components-level method for mapping uniformly partitioned periodic circumferential data from any source segment count to any target segment count. The first consumer will be the V14 18-zone heat-pipe representation mapped to the 12 circumferential radiator zones required by `RadiatorThermalShield.fortran_shield2`.

This change only adds and validates the reusable mapping method. It does not yet replace `RadiatorThermalShield._shield2_sector_indices()`, because a complete V14 integration must also define the reverse 12-zone shield-background application to the 18 source zones.

## Geometry

Source and target zones cover one 360-degree period with uniform widths:

```text
source_width = 360 / n_source
target_width = 360 / n_target
```

Both grids accept independent angular offsets. Wrapped intervals crossing 360/0 degrees are split before intersection.

The allocation matrix is:

```text
F[j,i] = overlap_angle(target_j, source_i) / source_width
```

Columns represent source-zone allocation and must satisfy:

```text
sum_j F[j,i] = 1
```

within floating-point tolerance.

## Data Mapping

### Extensive quantities

For source totals such as power, area, or represented heat-pipe count:

```python
target_totals = F @ source_totals
```

Column conservation guarantees preservation of the global total.

### Intensive quantities

For `T^4`, heat-flux density, or other intensive values with source physical weights `a`:

```python
numerator = F @ (a * source_values)
denominator = F @ a
target_values = numerator / denominator
```

A target zone with zero denominator is invalid and raises `ValueError`.

For thermal-shield radiation, callers pass `mean(T_fin_cells ** 4)`; the mapping utility does not raise temperatures to the fourth power internally.

## V14 18-to-12 Case

With aligned zero angles:

```text
target[2m]   receives all source[3m] and half source[3m+1]
target[2m+1] receives half source[3m+1] and all source[3m+2]
m = 0..5
```

For V14, `source_weights` are the 18 `hp_multipliers` values. Equal per-pipe radiation area cancels from the normalized intensive mapping. The global symmetric-ring multiplier also cancels for `T^4`, but must remain in later total-area and total-power accounting.

## Public API

Create `Components/circumferential_mapping.py` with:

```python
build_uniform_circumferential_mapping(
    n_source: int,
    n_target: int,
    source_offset_deg: float = 0.0,
    target_offset_deg: float = 0.0,
) -> np.ndarray
```

and:

```python
map_circumferential_intensive(
    source_values,
    mapping,
    source_weights=None,
) -> np.ndarray
```

The mapper accepts the source dimension on the final axis so scalar zone arrays and batched arrays can share the same matrix.

No class, registry, or configuration object is introduced.

## Validation

Reject:

- non-positive segment counts;
- non-finite offsets;
- mapping/source shape mismatch;
- non-finite source values or weights;
- negative source weights;
- target zones with zero total physical weight.

Tests cover:

- analytic aligned 18-to-12 coefficients;
- periodic wrap with angular offsets;
- matrix column sums;
- extensive total conservation;
- V14 5/6 multiplier-weighted `T^4` aggregation;
- constant-field preservation;
- invalid input handling.

## Documentation

Add `Components/CIRCUMFERENTIAL_MAPPING_GUIDE.md` with formulas and the V14 example. Add a short navigation entry to `Components/COMPONENTS_DETAILED_INTRO.md`.

## Deferred Integration

A later V14 shield task will:

1. extract the 18 representative heat-pipe `mean(T^4)` values;
2. map them to 12 shield input zones using this utility;
3. map the 12 shield background `T^4` values back across the 18 source zones using angular overlap;
4. include the symmetric physical radiator area in energy closure;
5. reassess whether the existing 78-pipe-derived shield view coefficients are valid for the 200-heat-pipe geometry.
