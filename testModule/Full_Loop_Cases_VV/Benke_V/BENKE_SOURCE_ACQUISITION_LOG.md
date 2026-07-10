# Benke source acquisition log

## Current local source state

Local workspace inspection found no Benke thesis PDF, no Benke-specific sleeve-temperature figure crop, and no Benke water-balance figure crop. The available local evidence is limited to the curated summary in `TOPAZII_VV_public_experimental_data.md`, which provides geometry, governing equations, reference ranges, and the list of recommended validation outputs.

## Searches performed on 2026-07-02

Local searches:

```text
rg --files | rg -i "benke|single-cell|single_cell|thermionic fuel element|operational testing|thermal modeling|\.pdf$|\.docx$"
Get-ChildItem -Recurse -File testModule\Full_Loop_Cases_VV | Select-String "Benke|benke|thermal|sleeve|water"
```

Web searches attempted:

```text
"Operational Testing and Thermal Modeling" "TOPAZ-II Single-Cell" Benke PDF
"Benke" "TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand"
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand"
site:calhoun.nps.edu Benke TOPAZ thermionic fuel element
site:apps.dtic.mil Benke TOPAZ II single cell thermionic
site:archive.org Benke TOPAZ-II thermionic fuel element
TOPAZ-II single cell thermionic fuel element test stand Benke Naval Postgraduate School thesis 1994
```

Result: no directly usable Benke PDF or digitizable Benke thermal validation figure was located in the current local workspace or through the attempted public search. The search result quality was poor and mostly returned generic TOPAZ references.

## Consequence for validation status

The Benke thermal model can currently be checked against literature ranges:

- active-zone power near 3003 W for the typical case;
- regulated He effective conductivity range 0.073-0.087 W/(m K);
- water-side heat-transfer coefficient range 528-1012 W/(m2 K);
- exact water-side energy balance closure.

It cannot yet be called a completed curve-level Benke validation because these required digitized outputs are missing:

- 12-point sleeve thermocouple temperatures;
- cooling-water outlet temperature or water-side delta-T;
- axial temperature distribution;
- any original Benke figure/page traceability for those values.

## Data admission rule

Do not fill `experimental_data/*_digitized.csv` from Venable output power or from model output. Those files must contain digitized or otherwise traceable Benke experimental data only.

## Additional title-level search on 2026-07-02

A second public title-level search was attempted for:

```text
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand" Benke PDF
"S. M. Benke" "TOPAZ-II" "Single-Cell"
"Operational Testing" "Thermal Modeling" "TOPAZ-II" "Benke"
"TOPAZ-II Single Cell Thermionic Fuel Element Test Stand" "Benke"
"single cell thermionic fuel element test stand" "thermal modeling"
"Benke" "Venable" "AIP Conference Proceedings" "TOPAZ-II"
```

No directly usable Benke PDF, sleeve thermocouple curve, water-balance curve, or digitizable figure page was found. The validation framework therefore remains data-ready but not curve-level complete.

## Manifest and domain-targeted search update on 2026-07-02

Local figure manifest check:

```text
testModule/Full_Loop_Cases_VV/TOPAZII_VV_public_figures/00_manifest.csv
```

Result: the current figure package contains Paramonov/El-Genk V-71/TSET figures and Venable 1995 single-cell electrical figures only. It does not contain Benke sleeve thermocouple figures, Benke water-balance figures, or Benke thesis page crops.

Additional web searches attempted with exact title and domain targeting:

```text
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand" site:calhoun.nps.edu
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand" site:apps.dtic.mil
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand" site:osti.gov
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand" site:worldcat.org
Benke TOPAZ-II Naval Postgraduate School thesis thermionic
"Benke" "TOPAZ-II" "Naval Postgraduate School"
"Benke" "Thermionic Fuel Element" "thesis"
"Operational Testing" "Thermionic Fuel Element Test Stand"
"Single-Cell Thermionic Fuel Element Test Stand"
"single-cell TFE test stand" Benke
"TISA" "Benke" "TOPAZ"
"regulated He gap" "TOPAZ" "Benke"
"Operational Testing and Thermal Modeling" "AIP Conference Proceedings"
"Operational Testing and Thermal Modeling of a TOPAZ-II Single Cell"
"TOPAZ-II Single Cell Thermionic Fuel Element" "AIP"
"AIP Conference Proceedings" "TOPAZ-II" "Benke" "Venable"
```

Result: no directly usable Benke PDF, figure image, curve data, or table data was located. Current Benke_V validation status therefore remains data-ready but blocked from completed curve-level validation by missing traceable Benke experimental measurements.
