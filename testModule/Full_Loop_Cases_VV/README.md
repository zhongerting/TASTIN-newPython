# Full Loop Cases V&V

This directory is reserved for full-loop validation and verification cases built from public TOPAZ-II experimental references.

## Reference Material

- `TOPAZII_VV_public_experimental_data.md`
- `TOPAZII_VV_public_figures/`

The public TOPAZ-II references point to a layered V&V strategy:

- V-71/TSET unit tests: primary full-loop thermal-hydraulic reference for NaK flow, core inlet/outlet temperature, pressure, radiator/header behavior, and 106 kW local coolant/saturation/subcooling data.
- Venable 1995 single-cell TFE data: single-TFE electrical reference for I-V curves, maximum power, efficiency, and cesium pressure trends.
- Benke 1994 single-cell TFE test stand: single-TFE thermal resistance and insulation reference, especially He gap and collector sleeve temperatures.

## Modeling Policy

Do not assume an existing case such as V12 or V13 is the validation structure.

Each V&V case should be modeled from the reference experiment being validated. For the first full-loop case, use the V-71/TSET documentation to define a standalone test model with:

- 37 electrically heated TFE coolant channels.
- NaK-78 primary loop geometry from the V-71 coolant-loop tables.
- EM pump / pressure-loss behavior from the V-71 pump and flow figures.
- Volume compensator pressure boundary from the documented gas/NaK volumes.
- 78 radiator tubes plus upper/lower headers under the ground TSET boundary, not a space-only radiator assumption.

Existing repository components and runners may be reused only after checking that their geometry, topology, boundary conditions, and diagnostics match the target experiment.
