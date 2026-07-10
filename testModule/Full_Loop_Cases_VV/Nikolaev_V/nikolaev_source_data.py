from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple


SOURCE_TITLE = 'A single-cell TFE mock-up of the thermionic nuclear power system "Space-R"'
SOURCE_CITATION = "AIP Conference Proceedings 324, 815 (1995)"
SOURCE_DOI = "10.1063/1.47120"
LOCAL_PDF = r"e:\文献阅读\nikolaev1995.pdf"


@dataclass(frozen=True)
class Table1Characteristic:
    parameter: str
    unit: str
    topaz_ii: float
    space_r: float
    note: str = ""


@dataclass(frozen=True)
class Table2OperatingPoint:
    voltage_v: float
    current_a: float
    thermal_power_kw: float
    emitter_temperature_k: float
    efficiency_percent: float

    @property
    def electric_power_w(self) -> float:
        return self.voltage_v * self.current_a


@dataclass(frozen=True)
class FuelTemperaturePoint:
    free_volume_percent: float
    radial_factor: float
    max_fuel_temperature_k: float


@dataclass(frozen=True)
class CapillaryLimitPoint:
    free_volume_percent: float
    radial_factor: float
    max_capillary_diameter_mm: float


TABLE1_CHARACTERISTICS: Tuple[Table1Characteristic, ...] = (
    Table1Characteristic("output_power", "W", 132.0, 300.0),
    Table1Characteristic("efficiency", "percent", 5.0, 7.3),
    Table1Characteristic("effective_height", "cm", 46.0, 70.0, "OCR row label appears as H~n (era)."),
    Table1Characteristic(
        "emitter_cladding_thickness",
        "mm",
        1.15,
        2.3,
        "The extracted table label says emitter diameter, but the text describes thickening the emitter cladding to 2.3 mm.",
    ),
    Table1Characteristic("emitter_temperature", "K", 1820.0, 1880.0),
    Table1Characteristic("collector_temperature", "K", 870.0, 870.0),
)


TABLE2_OPERATING_POINTS: Tuple[Table2OperatingPoint, ...] = (
    Table2OperatingPoint(0.7, 429.0, 4.2, 1880.0, 7.1),
    Table2OperatingPoint(0.8, 375.0, 4.1, 1890.0, 7.3),
    Table2OperatingPoint(0.9, 333.0, 4.1, 1910.0, 7.3),
)


TABLE3_FUEL_TEMPERATURES: Tuple[FuelTemperaturePoint, ...] = (
    FuelTemperaturePoint(20.0, 1.0, 2110.0),
    FuelTemperaturePoint(20.0, 1.15, 2220.0),
    FuelTemperaturePoint(30.0, 1.0, 2070.0),
    FuelTemperaturePoint(30.0, 1.15, 2175.0),
    FuelTemperaturePoint(40.0, 1.0, 2040.0),
    FuelTemperaturePoint(40.0, 1.15, 2140.0),
)


TABLE4_CAPILLARY_LIMITS: Tuple[CapillaryLimitPoint, ...] = (
    CapillaryLimitPoint(20.0, 1.0, 0.26),
    CapillaryLimitPoint(20.0, 1.15, 0.15),
    CapillaryLimitPoint(30.0, 1.0, 0.35),
    CapillaryLimitPoint(30.0, 1.15, 0.20),
    CapillaryLimitPoint(40.0, 1.0, 0.40),
    CapillaryLimitPoint(40.0, 1.15, 0.23),
)


MOCKUP_FACTS = {
    "prototype_name": "SC-320",
    "emitter_outer_diameter_mm": 17.3,
    "emitter_thickness_mm": 2.3,
    "heater_tungsten_length_mm": 350.0,
    "interelectrode_gap_mm": 0.5,
    "collector_temperature_for_vac_k": 870.0,
    "emitter_work_function_ev": 4.85,
    "vac_thermal_power_kw_values": [2.5, 3.0, 3.5, 4.0, 4.5, 5.0],
}
