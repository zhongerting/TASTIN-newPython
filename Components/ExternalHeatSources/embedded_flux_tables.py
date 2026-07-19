from __future__ import annotations

"""
Embedded orbital heat-flux tables derived from the legacy Fortran radiator input.

The original workflow read ``RadiatorInput.txt`` at runtime and selected a
partition heat-flux table by index.  For the Python refactor we keep the same
lookup idea, but store the table data directly in code so the model no longer
depends on external input-card files.
"""

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from .w0_8p12_sum_data import load_w0_8p12_sum_matrices

W0_8P12_ORBIT_PERIOD_S = 6552.0


def _parse_numeric_row(text: str) -> np.ndarray:
    compact = " ".join(text.strip().split())
    return np.fromstring(compact, sep=" ", dtype=float)


def _parse_numeric_matrix(text: str) -> np.ndarray:
    rows = []
    for line in text.strip().splitlines():
        stripped = line.strip()
        if stripped:
            rows.append(_parse_numeric_row(stripped))
    return np.vstack(rows)


_ZONE_TIME_TEXT = """
0.0 193.1 386.3 579.4 772.5 965.6 1158.8 1351.9 1545.0 1738.1 1931.3 2124.4
2317.5 2510.6 2703.8 2896.9 3090.0 3283.1 3476.3 3669.4 3862.5 4055.6 4248.8
4441.9 4635.0 4828.1 5021.3 5214.4 5407.5 5600.6 5793.8 5986.9 6180.0
"""

_ZONE_TABLES_TEXT = """
0.0 0.0 194.2 455.4 696.7 907.3 904.2 1076.0 1200.8 1272.9 1289.4 1249.4 1154.8 1009.6 820.2 594.9 343.5 77.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
0.1 0.1 199.6 428.9 639.3 821.7 819.9 967.1 1071.9 1129.8 1138.2 1096.8 1007.4 873.9 702.1 499.6 275.1 38.7 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
6.3 38.1 217.2 388.2 542.7 673.7 673.6 775.8 844.4 876.3 870.0 826.0 746.0 633.6 493.8 332.5 156.9 8.3 6.7 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3 6.3
47.0 151.2 251.5 346.3 428.6 493.1 494.7 539.0 560.7 558.9 533.7 486.2 418.4 333.3 234.6 126.6 40.0 32.4 26.0 24.2 24.1 24.1 24.1 24.1 24.1 24.1 24.1 24.1 24.1 24.1 24.1 24.1 47.0
265.2 283.2 293.4 302.2 303.2 293.3 296.3 275.6 244.9 205.5 159.2 124.0 121.0 114.7 105.5 93.7 79.9 64.7 51.6 47.6 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 265.2
497.5 429.3 350.2 276.5 201.1 163.0 162.6 184.4 201.8 214.1 220.7 221.3 215.9 204.8 188.3 167.3 142.7 115.8 91.8 84.6 84.2 84.2 84.2 84.2 84.2 84.2 84.2 84.2 84.2 84.2 84.2 84.2 497.5
706.9 563.0 406.3 262.1 200.0 236.2 235.6 266.9 291.9 309.4 318.7 319.3 311.3 294.9 271.0 240.5 204.8 165.9 131.7 121.9 195.2 121.4 121.4 121.4 121.4 121.4 121.4 121.4 121.4 121.4 121.4 121.4 706.9
870.2 670.0 454.3 257.6 254.5 301.7 300.9 342.0 374.9 398.3 411.1 412.7 403.2 382.8 352.6 313.7 268.0 217.9 172.0 157.4 336.5 156.6 156.6 156.6 156.6 156.6 156.6 156.6 156.6 156.6 156.6 156.6 870.2
956.7 726.5 480.6 245.8 283.8 336.3 335.4 380.9 417.4 443.3 457.4 459.1 448.3 425.5 391.8 348.5 297.6 241.8 191.3 175.4 411.0 174.3 174.3 174.3 174.3 174.3 174.3 174.3 174.3 174.3 174.3 174.3 956.7
958.0 728.0 482.8 252.8 285.4 337.8 336.9 382.3 418.7 444.4 458.3 459.8 448.7 425.7 391.7 348.2 297.1 241.3 191.3 175.7 412.0 174.6 174.6 174.6 174.6 174.6 174.6 174.6 174.6 174.6 174.6 174.6 958.0
872.7 672.5 457.3 261.0 254.7 301.7 300.8 341.6 374.3 397.5 410.1 411.5 401.9 381.4 351.1 312.3 266.6 216.6 171.1 157.0 338.0 156.2 156.2 156.2 156.2 156.2 156.2 156.2 156.2 156.2 156.2 156.2 872.7
709.5 565.6 408.3 263.2 195.9 231.8 231.1 262.2 287.1 304.6 314.1 315.1 307.5 291.6 268.3 238.4 203.4 165.0 130.5 120.2 196.3 119.7 119.7 119.7 119.7 119.7 119.7 119.7 119.7 119.7 119.7 119.7 709.5
501.9 433.9 355.1 281.1 205.5 161.2 160.7 182.3 199.5 211.7 218.2 218.8 213.5 202.5 186.3 165.6 141.3 114.7 91.1 83.7 83.3 83.3 83.3 83.3 83.3 83.3 83.3 83.3 83.3 83.3 83.3 83.3 501.9
270.5 288.7 299.2 308.1 309.2 299.3 296.8 276.2 245.6 206.3 160.0 123.0 119.9 113.6 104.4 92.7 78.9 64.0 51.0 47.2 46.9 46.9 46.9 46.9 46.9 46.9 46.9 46.9 46.9 46.9 46.9 46.9 270.5
49.4 153.8 253.8 347.9 429.4 493.3 489.7 533.5 554.9 553.0 527.8 480.5 413.1 328.7 230.8 123.8 35.6 28.9 23.1 21.3 21.2 21.2 21.2 21.2 21.2 21.2 21.2 21.2 21.2 21.2 21.2 21.2 49.4
4.8 41.0 220.1 390.9 545.1 675.7 671.4 773.4 841.8 873.6 867.5 823.5 743.8 631.8 492.4 331.7 156.7 6.3 5.1 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8 4.8
0.0 0.0 202.5 431.8 642.4 824.8 820.3 967.5 1072.4 1130.4 1138.9 1097.6 1008.3 874.9 703.2 500.8 276.5 40.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 195.3 456.5 697.8 908.5 904.4 1076.3 1201.0 1273.2 1289.7 1249.8 1155.2 1010.0 820.6 595.4 344.0 77.7 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 32.7 300.3 554.7 784.8 781.2 977.7 1131.4 1235.6 1285.7 1279.5 1217.3 1101.9 938.2 733.4 496.5 237.9 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 38.1 273.2 496.3 697.6 695.5 866.9 1000.3 1090.0 1131.9 1124.3 1067.4 963.9 818.1 636.4 426.9 198.7 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
5.4 5.4 54.7 230.2 396.2 545.0 544.5 670.1 766.6 829.8 856.8 846.6 799.5 717.6 604.5 465.1 305.6 132.9 5.9 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4 5.4
19.8 19.9 85.4 182.2 273.6 353.9 355.1 421.0 469.3 497.8 505.5 491.8 457.5 403.9 333.5 249.3 155.0 54.8 22.1 20.0 19.8 19.8 19.8 19.8 19.8 19.8 19.8 19.8 19.8 19.8 19.8 19.8 19.8
120.8 127.1 132.1 142.4 152.3 158.0 160.7 161.5 157.2 148.1 134.5 123.3 120.5 114.4 105.4 93.9 80.3 65.4 52.0 47.5 47.1 47.1 47.1 47.1 47.1 47.1 47.1 47.1 47.1 47.1 47.1 47.1 120.8
352.2 270.8 184.3 108.6 128.1 152.2 151.8 172.8 189.8 202.0 208.8 210.0 205.5 195.5 180.4 161.0 138.0 112.6 89.1 80.6 163.1 79.9 79.9 79.9 79.9 79.9 79.9 79.9 79.9 79.9 79.9 79.9 352.2
567.3 408.9 243.0 155.0 193.8 230.1 229.5 261.0 286.4 304.5 314.6 316.1 309.0 293.7 270.8 241.3 206.4 168.1 132.8 174.9 348.2 120.0 120.0 120.0 120.0 120.0 120.0 120.0 120.0 120.0 120.0 120.0 567.3
735.4 520.0 295.1 203.8 254.7 302.6 301.8 343.5 377.0 401.0 414.4 416.5 407.3 387.2 357.1 318.2 272.4 221.9 175.0 262.3 494.0 158.0 158.0 158.0 158.0 158.0 158.0 158.0 158.0 158.0 158.0 158.0 735.4
827.2 581.3 325.5 231.9 290.0 344.4 343.5 390.9 429.2 456.5 471.7 474.2 463.8 440.9 406.6 362.4 310.2 252.8 199.6 310.8 573.7 179.9 179.9 179.9 179.9 179.9 179.9 179.9 179.9 179.9 179.9 179.9 827.2
826.0 580.2 324.2 228.9 286.2 340.1 339.1 386.1 424.0 451.0 466.2 468.7 458.5 436.0 402.2 358.5 307.0 250.3 197.6 309.3 572.2 177.8 177.8 177.8 177.8 177.8 177.8 177.8 177.8 177.8 177.8 177.8 826.0
738.0 522.6 297.8 203.8 254.8 302.5 301.6 343.2 376.6 400.4 413.7 415.7 406.4 386.2 356.1 317.2 271.4 221.0 174.2 263.6 495.6 157.7 157.7 157.7 157.7 157.7 157.7 157.7 157.7 157.7 157.7 157.7 738.0
575.1 416.9 251.5 160.3 200.2 237.4 236.8 269.1 295.2 313.7 323.9 325.3 317.9 301.9 278.2 247.7 211.8 172.3 136.1 181.0 354.5 123.5 123.5 123.5 123.5 123.5 123.5 123.5 123.5 123.5 123.5 123.5 575.1
360.8 279.6 193.8 119.5 134.1 159.1 158.7 180.4 198.0 210.5 217.4 218.5 213.6 203.0 187.2 166.8 142.8 116.3 92.2 83.8 169.8 83.1 83.1 83.1 83.1 83.1 83.1 83.1 83.1 83.1 83.1 83.1 360.8
126.8 133.3 138.3 148.7 159.0 165.0 162.1 163.2 159.2 150.4 137.1 124.1 121.3 115.3 106.3 94.8 81.1 66.1 52.5 47.8 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 47.3 126.8
21.1 21.2 92.5 189.6 281.7 362.5 358.5 424.9 473.7 502.6 510.6 497.1 462.9 409.4 338.9 254.6 160.1 59.7 23.5 21.2 21.1 21.1 21.1 21.1 21.1 21.1 21.1 21.1 21.1 21.1 21.1 21.1 21.1
4.5 4.5 58.4 233.9 399.7 548.4 543.7 669.2 765.6 828.8 855.9 845.8 798.9 717.3 604.5 465.5 306.4 134.1 4.9 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5
0.0 0.0 41.1 276.4 499.6 700.9 696.0 867.5 1001.0 1090.8 1132.8 1125.3 1068.5 965.0 819.3 637.8 428.4 200.3 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
0.0 0.0 33.8 301.4 555.9 785.9 781.4 977.9 1131.7 1235.9 1286.0 1279.9 1217.7 1102.3 938.6 733.9 497.0 238.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
"""

_TOP_TIME_TEXT = """
0 91 182 273 364 455 546 637 728 819 910 1001 1092 1183 1274 1365 1456
1547 1638 1729 1820 1911 2002 2093 2184 2275 2366 2457 2548 2639 2730 2821
2912 3003 3094 3185 3276 3367 3458 3549 3640 3731 3822 3913 4004 4095 4186
4277 4368 4459 4550 4641 4732 4823 4914 5005 5096 5187 5278 5369 5460 5551
5642 5733 5824 5915 6006 6097 6188 6279 6370 6461 6552
"""

_TOP_TABLE_TEXT = """
103.1478863 1.40E+03 1392.669119 1383.365086 1374.880855 1366.908096
1360.31961 1356.850887 1350.108202 1343.499361 1337.081829 1331.321806
1325.218118 1319.491524 1315.385833 1311.756238 1309.901345 1308.502945
1308.041249 1308 1308 1308 1308 1308 1308 1308 1308 1308 1308 1308 1308
1308.016247 1308.17201 1308.831047 1311.157312 1314.404023 1318.882478
1324.319785 1329.709523 1335.891041 1341.853556 1348.684463 1355.26429
1362.598509 1369.137999 1376.868237 1385.071621 85.57576334 9.42E+01
103.3864644 112.9141958 122.411311 131.5243958 140.2847084 148.6145775
155.9226914 161.7054077 166.8411959 170.6186944 172.5733271 173.9686754
173.6007124 172.0703389 169.625775 165.3989577 159.3052414 153.2061211
145.6357688 137.5827657 128.9103664 120.2558878 111.7340314 102.7921365
"""

FORTRAN_RADIATOR_ZONE_TIME = _parse_numeric_row(_ZONE_TIME_TEXT)
FORTRAN_RADIATOR_ZONE_TABLES = _parse_numeric_matrix(_ZONE_TABLES_TEXT)
FORTRAN_TOP_SURFACE_TIME = _parse_numeric_row(_TOP_TIME_TEXT)
FORTRAN_TOP_SURFACE_TABLE = _parse_numeric_row(_TOP_TABLE_TEXT)


@dataclass(frozen=True)
class EmbeddedFluxTable:
    table_id: int
    name: str
    time: np.ndarray
    values: np.ndarray
    periodic: bool = True


class EmbeddedFluxTableLibrary:
    """Lookup container for embedded orbital heat-flux tables."""

    def __init__(self, tables: Dict[int, EmbeddedFluxTable]):
        self._tables = dict(tables)

    def get_table(self, table_id: int) -> EmbeddedFluxTable:
        key = int(table_id)
        if key not in self._tables:
            available = ", ".join(str(item) for item in self.available_ids())
            raise KeyError(f"Unknown orbital heat table id {key}. Available ids: {available}")
        return self._tables[key]

    def available_ids(self) -> Tuple[int, ...]:
        return tuple(sorted(self._tables))

    def has_table(self, table_id: int) -> bool:
        return int(table_id) in self._tables


def _scale_sample_time_to_period(sample_time: np.ndarray, orbit_period_s: float = None) -> np.ndarray:
    time = np.asarray(sample_time, dtype=float)
    if orbit_period_s is None:
        return time.copy()
    period = float(orbit_period_s)
    if not np.isfinite(period) or period <= 0.0:
        raise ValueError("orbit_period_s must be finite and positive")
    return time * (period / float(time[-1]))


@lru_cache(maxsize=None)
def load_csv_flux_table_library(
    path: str,
    orbit_period_s: float = None,
) -> EmbeddedFluxTableLibrary:
    """Load sample-indexed theta heat-flux columns as scalar lookup tables."""
    csv_path = Path(path).resolve()
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1, dtype=float)
    if data.ndim != 2 or data.shape[0] < 2 or data.shape[1] < 2:
        raise ValueError(f"External heat CSV '{csv_path}' must contain time plus flux columns")
    if not np.all(np.isfinite(data)) or not np.all(np.diff(data[:, 0]) > 0.0):
        raise ValueError(f"External heat CSV '{csv_path}' must contain finite values and increasing time")
    sample_time = _scale_sample_time_to_period(data[:, 0], orbit_period_s)
    return EmbeddedFluxTableLibrary({
        index: EmbeddedFluxTable(
            table_id=index,
            name=f"{csv_path.stem}_theta_{index:03d}",
            time=sample_time.copy(),
            values=data[:, index + 1].copy(),
            periodic=True,
        )
        for index in range(data.shape[1] - 1)
    })


def _build_fortran_orbital_heat_tables() -> Dict[int, EmbeddedFluxTable]:
    tables: Dict[int, EmbeddedFluxTable] = {}

    for table_id, values in enumerate(FORTRAN_RADIATOR_ZONE_TABLES, start=1):
        tables[table_id] = EmbeddedFluxTable(
            table_id=table_id,
            name=f"fortran_rad_zone_{table_id:02d}",
            time=FORTRAN_RADIATOR_ZONE_TIME.copy(),
            values=np.array(values, dtype=float, copy=True),
            periodic=True,
        )

    tables[1001] = EmbeddedFluxTable(
        table_id=1001,
        name="fortran_top_surface_reference",
        time=FORTRAN_TOP_SURFACE_TIME.copy(),
        values=FORTRAN_TOP_SURFACE_TABLE.copy(),
        periodic=True,
    )

    return tables

FORTRAN_ORBITAL_HEAT_TABLES = _build_fortran_orbital_heat_tables()
FORTRAN_ORBITAL_HEAT_TABLE_LIBRARY = EmbeddedFluxTableLibrary(FORTRAN_ORBITAL_HEAT_TABLES)


@dataclass(frozen=True)
class EmbeddedFluxMatrix:
    """Embedded orbital heat-flux matrix with one time column and N flux columns."""

    key: str
    time: np.ndarray
    values: np.ndarray
    source: str = ""
    periodic: bool = True


class EmbeddedFluxMatrixLibrary:
    """Lookup helper for embedded orbital heat-flux matrices."""

    def __init__(self, matrices: Dict[str, EmbeddedFluxMatrix]):
        self.matrices = dict(matrices)

    def get_matrix(self, key: str) -> EmbeddedFluxMatrix:
        try:
            return self.matrices[key]
        except KeyError as exc:
            available = ", ".join(sorted(self.matrices)) or "<none>"
            raise KeyError(f"Unknown embedded flux matrix '{key}'. Available: {available}") from exc

    def available_keys(self) -> Tuple[str, ...]:
        return tuple(sorted(self.matrices.keys()))


def _build_w0_8p12_orbital_heat_matrices() -> Dict[str, EmbeddedFluxMatrix]:
    """Build fixed sum-flux matrices for w0=8.12deg and i_s=58.5deg."""
    matrices = {}
    for key, data in load_w0_8p12_sum_matrices().items():
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(f"Embedded flux matrix '{key}' must contain time plus flux columns")
        matrices[key] = EmbeddedFluxMatrix(
            key=key,
            time=_scale_sample_time_to_period(data[:, 0], W0_8P12_ORBIT_PERIOD_S),
            values=data[:, 1:].astype(float),
            source="Components.ExternalHeatSources.w0_8p12_sum_data",
            periodic=True,
        )
    n18_path = Path(__file__).with_name("is58p5_w0_8p12_N18_sum.csv")
    n18_data = np.loadtxt(n18_path, delimiter=",", skiprows=1, dtype=float)
    matrices["is58p5_w0_8p12_N18_sum"] = EmbeddedFluxMatrix(
        key="is58p5_w0_8p12_N18_sum",
        time=_scale_sample_time_to_period(n18_data[:, 0], W0_8P12_ORBIT_PERIOD_S),
        values=n18_data[:, 1:].copy(),
        source=str(n18_path),
        periodic=True,
    )
    return matrices


W0_8P12_ORBITAL_HEAT_MATRICES = _build_w0_8p12_orbital_heat_matrices()
W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY = EmbeddedFluxMatrixLibrary(W0_8P12_ORBITAL_HEAT_MATRICES)

