import os
import sys
import numpy as np
import time
import logging
from pathlib import Path

from profiler import TEASAProfiler

logger = logging.getLogger(__name__)

# =========================================================
# [核心修复] 动态环境变量注入 (Dynamic Path Injection)
# 让 Wrapper 自己找到并加载同一目录下的 te_solver.pyd
# =========================================================
# 获取当前 ThermoCalcWrapper.py 所在的绝对路径
current_dir = os.path.dirname(os.path.abspath(__file__))

# 将该目录置于系统搜索路径的最高优先级 (Index 0)
_pyd_dir = os.environ.get("THERMOCALC_PYD_DIR")
if _pyd_dir:
    _pyd_dir = os.path.abspath(_pyd_dir)
    if _pyd_dir not in sys.path:
        sys.path.insert(0, _pyd_dir)

if current_dir not in sys.path:
    sys.path.insert(1 if _pyd_dir else 0, current_dir)

# =========================================================
# 宽容导入机制 (Graceful Import)
# =========================================================
try:
    import te_solver  # 现在 Python 会在 sys.path[0] 也就是 current_dir 里找到它
    HAS_TE_SOLVER = True
except ImportError as e:
    HAS_TE_SOLVER = False
    logger.warning(f"⚠️ Module 'te_solver' (C++ PyBind11) not found in {current_dir}! "
                   f"ThermoCalcModel instantiation will fail if TEC coupling is enabled. "
                   f"Detail: {e}")


_LOOKUP_LOADED_DB = None
_LOOKUP_LOADED_REGIONS = None


def _env_flag(name: str) -> bool:
    value = os.environ.get(name, "")
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _env_float(name: str, default: float) -> float:
    value = os.environ.get(name, "").strip()
    if not value:
        return float(default)
    try:
        return float(value)
    except ValueError:
        logger.warning("Invalid %s=%r; using %.6g", name, value, default)
        return float(default)

def _normalize_lookup_regions(regions=None):
    if regions is None:
        text = os.environ.get("THERMOCALC_LOOKUP_REGIONS", "").strip()
        if text:
            regions = [part.strip() for part in text.split(",") if part.strip()]
    if regions is None:
        return ("core",)
    if isinstance(regions, str):
        regions = [part.strip() for part in regions.split(",") if part.strip()]
    result = tuple(str(region).strip() for region in regions if str(region).strip())
    return result or ("core",)


def _find_default_lookup_database():
    runtime_root = Path(current_dir) / "emission_runtime_db_v2"
    for name in ("pcs_0p02_5torr_tc1500", "pcs_0p02_5torr"):
        candidate = runtime_root / name
        if (candidate / "runtime_dense_manifest.json").exists():
            return candidate
    return None


def _load_runtime_lookup_database(db_path: Path, manifest: dict, regions: tuple[str, ...]) -> int:
    if not hasattr(te_solver, "add_emission_runtime_block"):
        raise RuntimeError("te_solver does not expose runtime lookup API.")
    wanted = set(regions)
    for chunk in manifest["chunks"]:
        region = str(chunk["region"])
        if region not in wanted:
            continue
        data_path = db_path / chunk["output"]
        with np.load(data_path, allow_pickle=False) as data:
            te_solver.add_emission_runtime_block(
                str(chunk["chunk_id"]),
                int(chunk["priority"]),
                int(chunk["region_id"]),
                np.asarray(data["TE_axis"], dtype=np.float64),
                np.asarray(data["TC_axis"], dtype=np.float64),
                np.asarray(data["Vo_axis"], dtype=np.float64),
                np.asarray(data["Tcs_axis"], dtype=np.float64),
                np.asarray(data["J"], dtype=np.float32),
                np.asarray(data["Vd"], dtype=np.float32),
                np.asarray(data["delta_V"], dtype=np.float32),
                np.asarray(data["phiE"], dtype=np.float32),
                np.asarray(data["phiC"], dtype=np.float32),
                np.asarray(data["lookup_safe"], dtype=np.uint8),
                np.asarray(data["zero_mask"], dtype=np.uint8),
            )
    return int(te_solver.emission_lookup_block_count())


def _load_dense_runtime_lookup_database(db_path: Path, manifest: dict, regions: tuple[str, ...]) -> int:
    if not hasattr(te_solver, "add_emission_dense_region"):
        raise RuntimeError("te_solver does not expose dense runtime lookup API.")
    wanted = set(regions)
    for region_name, meta in manifest["regions"].items():
        region = str(region_name)
        if region not in wanted:
            continue
        outputs = dict(meta.get("outputs", {}))
        binary_output = outputs.get("binary")
        if binary_output and hasattr(te_solver, "load_emission_dense_file"):
            te_solver.load_emission_dense_file(str(db_path / binary_output))
            continue
        npz_output = outputs.get("npz")
        if not npz_output:
            raise FileNotFoundError(f"Dense runtime region {region} has no loadable npz/binary output.")
        data_path = db_path / npz_output
        with np.load(data_path, allow_pickle=False) as data:
            point_count = int(np.prod(data["J"].shape))
            te_solver.add_emission_dense_region(
                region,
                int(meta["priority"]),
                int(meta["region_id"]),
                float(meta.get("d_gap", manifest.get("d_gap", 0.5))),
                np.asarray(data["TE_axis"], dtype=np.float64),
                np.asarray(data["TC_axis"], dtype=np.float64),
                np.asarray(data["Vo_axis"], dtype=np.float64),
                np.asarray(data["Tcs_axis"], dtype=np.float64),
                np.asarray(data["J"], dtype=np.float32),
                np.asarray(data["Vd"], dtype=np.float32),
                np.asarray(data["delta_V"], dtype=np.float32),
                np.asarray(data["phiE"], dtype=np.float32),
                np.asarray(data["phiC"], dtype=np.float32),
                np.asarray(data["lookup_safe_bits"], dtype=np.uint8),
                np.asarray(data["zero_mask_bits"], dtype=np.uint8),
                point_count,
            )
    if hasattr(te_solver, "emission_lookup_dense_region_count"):
        return int(te_solver.emission_lookup_dense_region_count())
    return int(te_solver.emission_lookup_region_count())


def _load_full_lookup_database(db_path: Path, manifest: dict, plan: dict, regions: tuple[str, ...]) -> int:
    wanted = set(regions)
    for chunk in plan["chunks"]:
        region = str(chunk["region"])
        if region not in wanted:
            continue
        priority = int(manifest["regions"][region]["priority"])
        raw_path = db_path / chunk["output"]
        data_path = raw_path.with_suffix(".optimized.npz")
        if not data_path.exists():
            data_path = raw_path
        with np.load(data_path, allow_pickle=False) as data:
            if "lookup_safe_flag" in data.files:
                safe = np.asarray(data["lookup_safe_flag"], dtype=np.uint8)
            else:
                safe = np.asarray(data["done"] & data["finite_flag"] & data["converged"], dtype=np.uint8)
            te_solver.add_emission_lookup_block(
                str(chunk["chunk_id"]),
                priority,
                np.asarray(data["TE_axis"], dtype=np.float64),
                np.asarray(data["TC_axis"], dtype=np.float64),
                np.asarray(data["Vo_axis"], dtype=np.float64),
                np.asarray(data["Tcs_axis"], dtype=np.float64),
                np.asarray(data["J"], dtype=np.float64),
                np.asarray(data["Vd"], dtype=np.float64),
                np.asarray(data["delta_V"], dtype=np.float64),
                np.asarray(data["phiE"], dtype=np.float64),
                np.asarray(data["phiC"], dtype=np.float64),
                safe,
            )
    return int(te_solver.emission_lookup_block_count())


def load_emission_lookup_database(db_dir: str, *, enable: bool = True, force: bool = False, regions=None) -> int:
    """Load the emission lookup database into the C++ te_solver singleton."""
    global _LOOKUP_LOADED_DB, _LOOKUP_LOADED_REGIONS
    if not HAS_TE_SOLVER:
        raise RuntimeError("te_solver is unavailable; cannot load emission lookup database.")
    required = ("clear_emission_lookup", "add_emission_lookup_block", "set_emission_lookup_enabled")
    missing = [name for name in required if not hasattr(te_solver, name)]
    if missing:
        raise RuntimeError(f"te_solver does not expose lookup API: {missing}")

    db_path = Path(db_dir).resolve()
    region_tuple = _normalize_lookup_regions(regions)
    if _LOOKUP_LOADED_DB == str(db_path) and _LOOKUP_LOADED_REGIONS == region_tuple and not force:
        te_solver.set_emission_lookup_enabled(bool(enable))
        return int(te_solver.emission_lookup_block_count())

    dense_manifest_path = db_path / "runtime_dense_manifest.json"
    runtime_manifest_path = db_path / "runtime_manifest.json"
    manifest_path = db_path / "manifest.json"
    plan_path = db_path / "chunk_plan.json"
    if not dense_manifest_path.exists() and not runtime_manifest_path.exists() and (not manifest_path.exists() or not plan_path.exists()):
        raise FileNotFoundError(f"Missing emission lookup manifest/chunk_plan under {db_path}")

    import json

    te_solver.clear_emission_lookup()
    if dense_manifest_path.exists():
        with dense_manifest_path.open("r", encoding="utf-8") as f:
            dense_manifest = json.load(f)
        loaded_blocks = _load_dense_runtime_lookup_database(db_path, dense_manifest, region_tuple)
    elif runtime_manifest_path.exists():
        with runtime_manifest_path.open("r", encoding="utf-8") as f:
            runtime_manifest = json.load(f)
        loaded_blocks = _load_runtime_lookup_database(db_path, runtime_manifest, region_tuple)
    else:
        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)
        with plan_path.open("r", encoding="utf-8") as f:
            plan = json.load(f)
        loaded_blocks = _load_full_lookup_database(db_path, manifest, plan, region_tuple)

    te_solver.set_emission_lookup_enabled(bool(enable))
    _LOOKUP_LOADED_DB = str(db_path)
    _LOOKUP_LOADED_REGIONS = region_tuple
    return int(loaded_blocks)


class ThermoCalcModel:
    """
    热离子能量转换器 (TEC) 的纯 Python 封装外壳
    负责处理几何参数分配、状态初始化、底层 C++ 对象通信以及结果提取。
    """

    def __init__(self,
                 n_elements: int,
                 n_nodes: int,
                 lookup_db: str = None,
                 enable_lookup: bool = None,
                 lookup_regions=None):
        self.N_elem = n_elements
        self.n_node = n_nodes

        if enable_lookup is None:
            enable_lookup = (
                _env_flag("THERMOCALC_ENABLE_LOOKUP")
                if "THERMOCALC_ENABLE_LOOKUP" in os.environ
                else True
            )
        default_lookup_db = _find_default_lookup_database()
        selected_lookup_db = (
            lookup_db
            if lookup_db is not None
            else os.environ.get("THERMOCALC_LOOKUP_DB")
            or (
                str(default_lookup_db)
                if default_lookup_db is not None
                else None
            )
        )
        self.lookup_db = selected_lookup_db
        self.lookup_enabled = bool(enable_lookup and selected_lookup_db)
        self.lookup_regions = lookup_regions
        self.lookup_loaded_blocks = 0
        if selected_lookup_db and self.lookup_enabled:
            self.lookup_loaded_blocks = load_emission_lookup_database(
                selected_lookup_db,
                enable=True,
                regions=lookup_regions,
            )
            logger.info("Loaded ThermoCalc emission lookup database: %s blocks", self.lookup_loaded_blocks)
        elif HAS_TE_SOLVER and hasattr(te_solver, "set_emission_lookup_enabled"):
            te_solver.set_emission_lookup_enabled(False)
        # 保存底层的 C++ 计算核心对象
        self._circuit = None
        self._input_data = te_solver.InputData()

        # 1. 初始化基础维度参数
        self._input_data.N_elements = self.N_elem
        self._input_data.n_axi = self.n_node

        # 2. 提供默认几何与物理参数矩阵 (用户可通过后续方法覆盖)
        self._input_data.Tcs = np.full((self.N_elem, self.n_node), 600.0)
        self._input_data.V_init = np.full((self.N_elem, self.n_node), 0.2)

        # 默认几何分配 (以单根 TFE 的经验值为基准，后续提供自定义接口)
        dl_val = 0.507 / float(self.n_node)
        self._input_data.dlE = np.full((self.N_elem, self.n_node), dl_val)
        self._input_data.dlC = np.full((self.N_elem, self.n_node), dl_val)

        self._input_data.crossAreaE = np.full(self.N_elem, 6.667e-5)
        self._input_data.crossAreaC = np.full(self.N_elem, 1.0786e-4)

        sideE_val = 0.00092855424159680002 * 25.0 / float(self.n_node)
        self._input_data.sideAreaE = np.full((self.N_elem, self.n_node), sideE_val)

        sideC_val = 0.00097592945800480005 * 25.0 / float(self.n_node)
        self._input_data.sideAreaC = np.full((self.N_elem, self.n_node), sideC_val)

        self._input_data.U_init = np.full(self.N_elem, 1.6)
        self._input_data.d_gap = np.full(self.N_elem, 0.5)
        self._input_data.Itarget = np.full(self.N_elem, 200.0)

        # 默认导线连接电阻/电压
        self._input_data.resistanceWire = np.zeros((self.N_elem, 4))
        wireU_single = np.array([0.8, 0.8, 0.0, 0.0])
        self._input_data.wireU = np.tile(wireU_single, (self.N_elem, 1))

        # 电路模式默认值 (初始以定电阻模式为主)
        # self._input_data.mode = te_solver.CalculationMode.FixedResistance
        self._input_data.mode = te_solver.CalculationMode.FixedVoltage
        # self._input_data.target_val = 0.0044
        self._input_data.target_val = 0.89 * 34
        self._input_data.I_total_init = 284.0
        self._mode_str = "fixed_u"
        self._load_curve_current = None
        self._load_curve_voltage = None
        if hasattr(self._input_data, "loadCurveCurrent"):
            self._input_data.loadCurveCurrent = np.array([0.0, 1000.0], dtype=float)
        if hasattr(self._input_data, "loadCurveVoltage"):
            self._input_data.loadCurveVoltage = np.array([0.0, 100.0], dtype=float)

        # --- 温度场占位符 ---
        # 初始默认处于 600K 均匀状态
        self._T_emitter = np.full((self.N_elem, self.n_node), 600.0)
        self._T_collector = np.full((self.N_elem, self.n_node), 600.0)
        self._zero_emission_skipped = False
        self._zero_emission_reason = None
    def set_temperatures(self, T_em: np.ndarray, T_co: np.ndarray):
        """
        更新发射极与接收极温度场
        :param T_em: 形状为 (N_elem, n_node) 的 numpy 数组
        :param T_co: 形状为 (N_elem, n_node) 的 numpy 数组
        """
        if T_em.shape != (self.N_elem, self.n_node) or T_co.shape != (self.N_elem, self.n_node):
            raise ValueError(f"温度场维度错误。预期 ({self.N_elem}, {self.n_node})")

        self._T_emitter = T_em.copy()
        self._T_collector = T_co.copy()

        # 如果底层电路已经创建，直接通过 Pybind11 更新属性
        if self._circuit is not None:
            for i in range(self.N_elem):
                self._circuit.TECs[i].Temitter = self._T_emitter[i, :]
                self._circuit.TECs[i].Tcollector = self._T_collector[i, :]

    def setup_circuit_mode(self, mode_str: str, target_value: float, I_guess: float = 150.0):
        """
        设置全局电路求解模式
        :param mode_str: 'fixed_R' (定电阻) 或 'fixed_U' (定电压) 或 'fixed_I' (定电流)
        :param target_value: 目标阻值(Ohm) 或 电压(V) 或 电流(A)
        :param I_guess: 迭代猜测初始电流 (A)
        """
        mode_key = mode_str.lower()
        self._mode_str = mode_key
        if mode_key == 'fixed_r':
            self._input_data.mode = te_solver.CalculationMode.FixedResistance
        elif mode_key == 'fixed_u':
            self._input_data.mode = te_solver.CalculationMode.FixedVoltage
        elif mode_key == 'fixed_i':
            if not hasattr(te_solver.CalculationMode, "FixedCurrent"):
                raise ValueError("fixed_i mode is not exposed by the ThermoCalc C++ binding.")
            if not np.isfinite(target_value) or target_value < 0.0:
                raise ValueError("fixed_i target current must be finite and non-negative.")
            self._input_data.mode = te_solver.CalculationMode.FixedCurrent
        elif mode_key == 'parallel_fixed_u':
            if not hasattr(te_solver.CalculationMode, "ParallelFixedVoltage"):
                raise ValueError("parallel_fixed_u mode is not exposed by the ThermoCalc C++ binding.")
            self._input_data.mode = te_solver.CalculationMode.ParallelFixedVoltage
        elif mode_key == 'parallel_fixed_i':
            if not hasattr(te_solver.CalculationMode, "ParallelFixedCurrent"):
                raise ValueError("parallel_fixed_i mode is not exposed by the ThermoCalc C++ binding.")
            self._input_data.mode = te_solver.CalculationMode.ParallelFixedCurrent
        elif mode_key == 'parallel_load_curve':
            if not hasattr(te_solver.CalculationMode, "ParallelLoadCurve"):
                raise ValueError("parallel_load_curve mode is not exposed by the ThermoCalc C++ binding.")
            self._input_data.mode = te_solver.CalculationMode.ParallelLoadCurve
            if self._load_curve_current is None or self._load_curve_voltage is None:
                max_i = max(float(I_guess) * 10.0, 1000.0)
                self.set_load_curve(
                    np.array([0.0, max_i], dtype=float),
                    np.array([0.0, float(target_value) * max_i], dtype=float),
                )
        else:
            raise ValueError(f"不支持的模式: {mode_str}")

        self._input_data.target_val = target_value
        self._input_data.I_total_init = I_guess
        if self._circuit is not None:
            self.build()

    def set_load_curve(self, current_a, voltage_v):
        """
        设置并联外部负载 U-I 曲线。

        :param current_a: 严格递增的总电流数组 [A]
        :param voltage_v: 对应外部负载端电压数组 [V]
        """
        current = np.asarray(current_a, dtype=float)
        voltage = np.asarray(voltage_v, dtype=float)
        if current.ndim != 1 or voltage.ndim != 1 or current.shape != voltage.shape:
            raise ValueError("load curve current_a and voltage_v must be 1D arrays with the same shape.")
        if current.size < 2:
            raise ValueError("load curve must contain at least two points.")
        if not np.all(np.isfinite(current)) or not np.all(np.isfinite(voltage)):
            raise ValueError("load curve values must be finite.")
        if np.any(np.diff(current) <= 0.0):
            raise ValueError("load curve current axis must be strictly increasing.")
        self._load_curve_current = current.copy()
        self._load_curve_voltage = voltage.copy()
        if hasattr(self._input_data, "loadCurveCurrent"):
            self._input_data.loadCurveCurrent = self._load_curve_current
        if hasattr(self._input_data, "loadCurveVoltage"):
            self._input_data.loadCurveVoltage = self._load_curve_voltage
        if self._circuit is not None:
            if not hasattr(self._circuit, "set_load_curve"):
                raise AttributeError("[ThermoCalcWrapper] C++ circuit lacks set_load_curve().")
            self._circuit.set_load_curve(self._load_curve_current, self._load_curve_voltage)

    def build(self):
        """
        执行 C++ 层的实例化，准备开始计算。
        仅在首次运行或几何、连接方式发生改变时调用。
        """
        # 将最新的温度场推入 InputData
        self._input_data.Temitter = self._T_emitter
        self._input_data.Tcollector = self._T_collector

        # 构建电路
        self._circuit = te_solver.create_circuit(self._input_data)

        # 根据所选模式，初始化底层 circuit 的运行标记
        if self._input_data.mode == te_solver.CalculationMode.FixedResistance:
            self._circuit.isFixedR = True
            self._circuit.isFixedU = False
            self._circuit.Rload = self._input_data.target_val
            # 初始猜测
            self._circuit.Iout = self._input_data.I_total_init
            self._circuit.Uout = self._circuit.Iout * self._circuit.Rload
            self._circuit.Utarget = self._circuit.Uout

        elif self._input_data.mode == te_solver.CalculationMode.FixedVoltage:
            self._circuit.isFixedR = False
            self._circuit.isFixedU = True
            self._circuit.Utarget = self._input_data.target_val
            self._circuit.Uout = self._input_data.target_val
            self._circuit.Iout = self._input_data.I_total_init
        elif hasattr(te_solver.CalculationMode, "FixedCurrent") and self._input_data.mode == te_solver.CalculationMode.FixedCurrent:
            self._circuit.isFixedR = False
            self._circuit.isFixedU = False
            self._circuit.isFixedI = True
            self._circuit.isParallelFixedU = False
            self._circuit.isParallelFixedI = False
            self._circuit.isParallelLoadCurve = False
            self._circuit.Itarget = self._input_data.target_val
            self._circuit.Iout = self._input_data.target_val
        elif hasattr(te_solver.CalculationMode, "ParallelFixedVoltage") and self._input_data.mode == te_solver.CalculationMode.ParallelFixedVoltage:
            self._circuit.isFixedR = False
            self._circuit.isFixedU = False
            self._circuit.isParallelFixedU = True
            self._circuit.isParallelFixedI = False
            self._circuit.isParallelLoadCurve = False
            self._circuit.Utarget = self._input_data.target_val
            self._circuit.Uout = self._input_data.target_val
            self._circuit.Iout = self._input_data.I_total_init
        elif hasattr(te_solver.CalculationMode, "ParallelFixedCurrent") and self._input_data.mode == te_solver.CalculationMode.ParallelFixedCurrent:
            self._circuit.isFixedR = False
            self._circuit.isFixedU = False
            self._circuit.isParallelFixedU = False
            self._circuit.isParallelFixedI = True
            self._circuit.isParallelLoadCurve = False
            self._circuit.Itarget = self._input_data.target_val
            self._circuit.Iout = self._input_data.I_total_init
        elif hasattr(te_solver.CalculationMode, "ParallelLoadCurve") and self._input_data.mode == te_solver.CalculationMode.ParallelLoadCurve:
            self._circuit.isFixedR = False
            self._circuit.isFixedU = False
            self._circuit.isParallelFixedU = False
            self._circuit.isParallelFixedI = False
            self._circuit.isParallelLoadCurve = True
            self._circuit.Utarget = self._input_data.target_val
            self._circuit.Uout = self._input_data.target_val
            self._circuit.Iout = self._input_data.I_total_init
            if self._load_curve_current is not None and hasattr(self._circuit, "set_load_curve"):
                self._circuit.set_load_curve(self._load_curve_current, self._load_curve_voltage)

    def _should_skip_zero_emission(self) -> bool:
        if _env_flag("THERMOCALC_DISABLE_ZERO_EMISSION_GUARD"):
            return False
        if self._T_emitter.size == 0 or self._T_collector.size == 0:
            return False
        if not np.all(np.isfinite(self._T_emitter)) or not np.all(np.isfinite(self._T_collector)):
            return False
        te_max = float(np.max(self._T_emitter))
        cutoff = _env_float("THERMOCALC_ZERO_EMISSION_TE_MAX_K", 1000.0)
        if te_max >= cutoff:
            return False
        self._zero_emission_reason = (
            f"max emitter temperature {te_max:.3f} K below "
            f"zero-emission cutoff {cutoff:.3f} K"
        )
        return True

    def _apply_zero_emission_result(self):
        if self._circuit is None:
            return
        zeros = np.zeros(self.n_node, dtype=float)
        for tec in self._circuit.TECs:
            tec.I = 0.0
            tec.U = 0.0
            tec.J = zeros.copy()
            tec.V = zeros.copy()
            tec.UE = zeros.copy()
            tec.UC = zeros.copy()
            tec.IEsecSingle = zeros.copy()
            tec.ICsecSingle = zeros.copy()
            tec.phiE = zeros.copy()
            tec.phiC = zeros.copy()
            tec.Vd = zeros.copy()
            tec.joulePowerE = zeros.copy()
            tec.joulePowerC = zeros.copy()
            tec.terminalPointUE1 = 0.0
            tec.terminalPointUE2 = 0.0
            tec.terminalPointUC1 = 0.0
            tec.terminalPointUC2 = 0.0
        self._circuit.Iout = 0.0
        if getattr(self._circuit, "isFixedU", False) or getattr(self._circuit, "isParallelFixedU", False):
            self._circuit.Uout = float(self._input_data.target_val)
        else:
            self._circuit.Uout = 0.0
        fixed_i_unavailable = getattr(self._circuit, "isFixedI", False) and float(self._input_data.target_val) > 0.0
        self._circuit.converged = not fixed_i_unavailable
        self._circuit.iterationCount = 1 if fixed_i_unavailable else 0
        if getattr(self._circuit, "isParallelFixedU", False):
            self._circuit.branchCurrents = np.zeros(self.N_elem, dtype=float)
            self._circuit.branchVoltages = np.full(self.N_elem, float(self._input_data.target_val), dtype=float)
        elif getattr(self._circuit, "isParallelFixedI", False) or getattr(self._circuit, "isParallelLoadCurve", False):
            self._circuit.branchCurrents = np.zeros(self.N_elem, dtype=float)
            self._circuit.branchVoltages = np.zeros(self.N_elem, dtype=float)
        self._zero_emission_skipped = True
    @TEASAProfiler.profile
    def calculate(self, verbose: bool = False) -> float:
        """
        触发底层计算 (瞬态单步更新时调用)
        :return: 计算耗时 (ms)
        """
        if self._circuit is None:
            self.build()

        t0 = time.time()
        if self._should_skip_zero_emission():
            self._apply_zero_emission_result()
        else:
            self._zero_emission_skipped = False
            self._zero_emission_reason = None
            self._circuit.calc()
        t1 = time.time()

        dt_ms = (t1 - t0) * 1000.0
        if verbose:
            print(f"[ThermoCalc] TEC System converged in {dt_ms:.2f} ms")
        return dt_ms

    def get_global_results(self):
        """ 获取系统级电路宏观结果 """
        if self._circuit is None:
            return None
        return {
            "Iout": self._circuit.Iout,
            "Uout": self._circuit.Uout,
            "Rload": self._circuit.Rload,
            "mode": self._mode_str,
            "converged": bool(getattr(self._circuit, "converged", True)),
            "iteration_count": int(getattr(self._circuit, "iterationCount", 0)),
            "branch_currents": np.array(getattr(self._circuit, "branchCurrents", []), dtype=float),
            "branch_voltages": np.array(getattr(self._circuit, "branchVoltages", []), dtype=float),
            "effective_rload": self._circuit.Rload,
            "zero_emission_skipped": bool(self._zero_emission_skipped),
            "zero_emission_reason": self._zero_emission_reason,
        }

    def get_tec_results(self, idx: int):
        """
        提取特定 TFE 元件的详细物理场结果
        :param idx: 元件索引 (0 到 N_elem-1)
        :return: 包含各项数据的字典
        """
        if self._circuit is None or idx < 0 or idx >= self.N_elem:
            return None

        tec = self._circuit.TECs[idx]
        return {
            "I": tec.I,
            "U": tec.U,
            "J": np.array(tec.J),
            "V": np.array(tec.V),
            "UE": np.array(tec.UE),
            "UC": np.array(tec.UC),
            "rhoE": np.array(tec.rhoE),
            "rhoC": np.array(tec.rhoC),
            "IEsecSingle": np.array(tec.IEsecSingle),
            "ICsecSingle": np.array(tec.ICsecSingle),
            "phiE": np.array(tec.phiE),
            "phiC": np.array(tec.phiC),
            "Vd": np.array(tec.Vd),
            "joulePowerE": np.array(tec.joulePowerE),
            "joulePowerC": np.array(tec.joulePowerC),
            "terminalPointUE1": tec.terminalPointUE1,
            "terminalPointUE2": tec.terminalPointUE2,
            "terminalPointUC1": tec.terminalPointUC1,
            "terminalPointUC2": tec.terminalPointUC2,
            "TE": np.array(tec.Temitter),
            "TC": np.array(tec.Tcollector)
        }

    def set_tcs(self, tcs_val):
        """
        [通信接口 -> C++] 动态设置/更新铯池温度 (Tcs)

        :param tcs_val: 可以是全局统一的标量温度 [K]，也可以是与网格匹配的 (N_elem, n_node) 的 NumPy 数组
        """
        # --- 1. 数据形状安全对齐 (Safety Alignment) ---
        if np.isscalar(tcs_val):
            # 如果传入的是标量，直接利用 NumPy 的广播机制填满整个数组
            self._input_data.Tcs[:] = float(tcs_val)
        else:
            # 如果传入的是数组，必须进行严格的维度校验，防止 C++ 底层越界崩溃
            tcs_arr = np.asarray(tcs_val, dtype=float)
            if tcs_arr.shape != (self.N_elem, self.n_node):
                raise ValueError(
                    f"[ThermoCalcWrapper] Tcs array shape mismatch! Expected {(self.N_elem, self.n_node)}, got {tcs_arr.shape}")
            self._input_data.Tcs[:] = tcs_arr

        # --- 2. C++ 运行时热更新 (Hot Update) ---
        if self._circuit is not None:
            if not hasattr(self._circuit, 'set_tcs'):
                raise AttributeError(
                    "[ThermoCalcWrapper] The underlying C++ object lacks 'set_tcs'. Cannot update Tcs at runtime!")
            self._circuit.set_tcs(self._input_data.Tcs)

    def set_rload(self, rload_val: float):
        """
        [通信接口 -> C++] 动态设置/更新系统的外部负载电阻

        :param rload_val: 外部负载电阻 [Ohm]
        """
        # --- 1. 更新初始化结构体 (兜底保护) ---
        # 兼容你 _input_data 中可能存在的各种命名习惯
        r_val = float(rload_val)
        if hasattr(self._input_data, 'Rload'):
            self._input_data.Rload = r_val
        elif hasattr(self._input_data, 'R_load'):
            self._input_data.R_load = r_val

        # --- 2. C++ 运行时热更新 (Hot Update) ---
        if self._circuit is not None:
            # 同样进行防御性探测
            if hasattr(self._circuit, 'set_rload'):
                self._circuit.set_rload(r_val)
            elif hasattr(self._circuit, 'Rload'):
                # 根据您代码中 get_global_results() 的实现，_circuit 应该暴露了 Rload 属性
                self._circuit.Rload = r_val
            else:
                raise AttributeError(
                    "[ThermoCalcWrapper] The underlying C++ object lacks 'set_rload' method or 'Rload' property. Cannot update load resistance at runtime!")
