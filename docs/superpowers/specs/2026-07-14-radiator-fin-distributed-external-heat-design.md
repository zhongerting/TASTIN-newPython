# Radiator Fin Distributed External Heat Design

## Goal

Add an opt-in high-fidelity external-heat mode to `RadiatorPipeWithFin`. Tube-wall absorption remains a wall boundary load, while fin absorption enters the quasi-steady fin equation using the single-sided projected illuminated area. Existing cases and numerical results remain unchanged unless the new mode is explicitly selected.

## Compatibility

- Keep the current lumped loading mode as the default.
- Preserve the existing `external_heat_bc` and `external_heat_source` attributes.
- A `RadiatorPipeWithFin` without a fin external source must follow its current equations exactly.
- Existing V15 constructors and runners require no new arguments.
- Reject unsupported loading-mode names during model construction.

## Modes

### `lumped_root_area` (default)

Keep the current behavior:

```text
Q_wall_boundary = q_external * (A_tube_projected + A_fin_projected)
Q_fin_absorption = 0
```

One `ExternalHeatFluxBC` applies the combined absorbed power to the tube-wall outer boundary.

### `distributed_fin_absorption`

Split the same total illuminated area:

```text
Q_wall_boundary = q_external * A_tube_projected
Q_fin_absorption = q_external * A_fin_projected
```

The tube-wall load uses one `ExternalHeatFluxBC`. The fin load is passed to the quasi-steady fin solver and is not added to the wall BC.

## Fin Equation

Each radiator tube has one circumferential table value at a given time. That heat-flux density is applied uniformly to all axial slices. Each slice uses only its own projected area.

For one axial slice, the two physical half-fins are represented by the existing doubled-width strip:

```text
A_fin_projected = fin_strip_width * fin_height
Q_fin_absorption = q_external * A_fin_projected
```

For one fin-height cell:

```text
dx = fin_height / n_fin_width
Q_abs_cell = q_external * fin_strip_width * dx
```

`Q_abs_cell` enters the tridiagonal fin equation right-hand side. Fin radiation continues to use the existing double-sided radiating perimeter. The root heat flow becomes:

```text
Q_fin_net_from_root = Q_fin_radiation - Q_fin_absorption
```

This matches the physical and numerical treatment already used by `HPwithFin.distributed_fin_absorption`.

## API

`RadiatorPipeWithFin` gains:

- `set_fin_external_heat_source(source, illuminated_area_scale=1.0)`
- `configure_external_heat_accounting(source, tube_area_array, fin_area_array)`
- `get_external_heat_absorption_distribution(current_time)`
- `last_fin_absorption_distribution`

`V15PipeFinRadiatorConfig` gains:

- `external_heat_fin_loading_mode: str = lumped_root_area`

The V15 adapter selects the existing lumped path or the new distributed path. Accounting is diagnostic only and must not create another heat source.

## Diagnostics

Expose per-axial-slice values in watts:

- `tube_absorption`
- `fin_absorption`
- `total_absorption`
- existing `fin_radiation`
- corrected `fin_net_from_root`

Repeated source evaluation is allowed, but only the wall BC and fin equation apply energy. Diagnostics must remain side-effect free.

## Validation

1. Existing V15 topology test passes with the default lumped mode.
2. Existing external-heat and radiator tests remain unchanged.
3. A focused distributed-mode test verifies:
   - exactly one wall `ExternalHeatFluxBC`;
   - wall BC area excludes fin area;
   - fin absorption equals `q_external * single-sided projected fin area`;
   - axial repetition does not multiply total power beyond summed segment area;
   - `Q_fin_net_from_root = Q_fin_radiation - Q_fin_absorption`.
4. Build two reduced V15 systems, one per mode, and verify equal prescribed total external power at the same table time.
5. Run a short temporary V15 smoke in distributed mode and require finite temperatures and finite heat diagnostics.

## Documentation

Update:

- `Components/RADIATORPIPEWITHFIN_DETAILED_INTRO.md`
- `Components/EXTERNALHEATSOURCES_DETAILED_INTRO.md`
- `testModule/Full_Loop_Cases/AI_AGENT_FULL_LOOP_CASES_ANALYSIS.md`

Record the default compatibility mode, the opt-in distributed mode, area definitions, energy identity, and validation commands.
