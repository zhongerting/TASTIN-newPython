"""
TASTIN 外热流模块 (External Heat Sources)
用于轨道空间核电源系统的外部热流边界条件抽象与实现。

典型外热流来源:
1. 太阳辐射热流 (Solar Heat Flux) - 基于轨道几何计算
2. 地球反照率热流 (Albedo Heat Flux) - 地球反射的太阳辐射
3. 地球红外热流 (Earth IR Flux) - 地球大气层红外辐射

设计原则:
- ExternalHeatSource: 物理模型接口，返回热流密度 [W/m²]
- ExternalHeatFluxBC: 边界条件封装类，继承 FluxBC，
  在 update_state() 中将 [W/m²] 转换为 [W]（乘以面积），
  并利用 BoundaryRegion 的自动更新机制实现时间相关的热流更新
- 与 BoundaryRegion 的 DynamicRadiationFluxBC 架构保持一致
"""

from abc import ABC, abstractmethod
from typing import Tuple, Union, Optional
import numpy as np

from .embedded_flux_tables import (
    FORTRAN_ORBITAL_HEAT_TABLE_LIBRARY,
    W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY,
    W0_8P12_ORBIT_PERIOD_S,
    EmbeddedFluxMatrixLibrary,
    EmbeddedFluxTableLibrary,
    load_csv_flux_table_library,
)


class BaseExternalHeatSource(ABC):
    """
    外热流抽象基类
    定义时间相关和位置相关的热流计算接口
    返回热流密度 [W/m²]
    """

    def __init__(self, shape: Tuple[int, ...]):
        """
        :param shape: 热流数据的维度形状，对应边界节点数
        """
        self.shape = shape
        self.bc_type = "external_heat_source"

    def _broadcast_flux(self, flux: Union[float, np.ndarray]) -> np.ndarray:
        """将标量或数组热流安全转换为与边界一致的数组形状。"""
        flux_array = np.asarray(flux, dtype=float)
        if flux_array.shape == ():
            return np.full(self.shape, float(flux_array), dtype=float)
        if flux_array.shape != self.shape:
            flux_array = np.broadcast_to(flux_array, self.shape)
        return np.array(flux_array, dtype=float, copy=True)

    @abstractmethod
    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        根据当前时间计算热流密度分布
        :param time: 当前仿真时间 [s]
        :return: 热流密度数组 [W/m²]，正值表示热量流入边界
        """
        pass

    @abstractmethod
    def update_params(self, **kwargs):
        """
        更新物理参数 (如轨道参数、姿态角等)
        """
        pass


class ExternalHeatFluxBC:
    """
    外热流边界条件封装类

    功能：
    1. 持有 ExternalHeatSource 引用
    2. 在 update_state() 中将热流密度 [W/m²] 转换为热流 [W]
    3. 利用 BoundaryRegion.update_internal_state() 的自动调用机制实现时间更新

    设计原理：
    - 继承自 FluxBC 的接口模式
    - update_state() 会被 BoundaryRegion.update_internal_state() 在每次迭代时自动调用
    - current_time 属性由 BoundaryRegion 设置
    - 这与 DynamicRadiationFluxBC 的工作原理完全一致
    """

    def __init__(self,
                 heat_source: BaseExternalHeatSource,
                 area_array: np.ndarray):
        """
        :param heat_source: 外热源对象
        :param area_array: 边界面积数组 [m²]，用于将热流密度转换为热流
        """
        self.heat_source = heat_source
        self.area_array = np.array(area_array, dtype=float)

        if self.area_array.shape != self.heat_source.shape:
            raise ValueError(f"Area array shape {self.area_array.shape} mismatch with heat source shape {self.heat_source.shape}")

        self.q_flux = np.zeros(self.heat_source.shape, dtype=float)
        self.bc_type = "flux"
        self.shape = self.heat_source.shape
        self.current_time = 0.0

    def update_state(self, T_node: np.ndarray):
        """
        [关键方法] 更新热流值

        此方法会被 BoundaryRegion.update_internal_state() 在每次迭代时自动调用
        实现了时间相关的热流自动更新

        :param T_node: 节点温度（未使用，但保留接口一致性）
        """
        q_density = self.heat_source.get_heat_flux(self.current_time)
        self.q_flux[:] = q_density * self.area_array

    def compute_flux_from_node(self, T_node: np.ndarray, R_internal: np.ndarray) -> np.ndarray:
        """
        返回当前热流 [W]
        :param T_node: 节点温度
        :param R_internal: 内部热阻
        :return: 热流数组 [W]
        """
        return self.q_flux

    def update_params(self, **kwargs):
        """转发参数更新给热源"""
        self.heat_source.update_params(**kwargs)


class OrbitalHeatSource(BaseExternalHeatSource):
    """
    轨道太阳辐射热源

    物理模型 (参考 Fortran sun_heat.f90):
    - 考虑轨道高度、轨道倾角、太阳方位角
    - 计算任意时刻飞行器表面的太阳入射因子
    - 考虑地球阴影遮挡 (Umbra/Penumbra)

    轨道参数:
    - 轨道高度 h [km]
    - 轨道倾角 i_s [deg] (太阳同步轨道通常为 0° 或 180°)
    - 轨道周期 TT [s]
    - 太阳常数 S [W/m²] (通常 1361 W/m²)

    表面参数:
    - 表面微元法向与太阳入射方向的夹角余弦 X = cos(iz_s)*sin(w0) + sin(iz_s)*cos(w0)*cos(l0-l_s)
    - 其中 iz_s 是太阳天顶角，w0/l0 是微元法向的球坐标
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 solar_constant: float = 1361.0,
                 orbit_height: float = 800.0,
                 orbit_period: float = 7644.0,
                 orbit_inclination: float = 0.0,
                 surface_normal_angles: Union[float, np.ndarray] = None):
        """
        :param shape: 边界形状
        :param solar_constant: 太阳常数 [W/m²]，默认 1361 W/m²
        :param orbit_height: 轨道高度 [km]，默认 800 km (LEO)
        :param orbit_period: 轨道周期 [s]，默认 7644 s (约 127 分钟)
        :param orbit_inclination: 轨道倾角 [deg]，0° 为太阳同步轨道
        :param surface_normal_angles: 表面法向球坐标 (w0, l0) 或标量
                                     w0: 极角 [rad]，从 z 轴测量
                                     l0: 方位角 [rad]，从 x 轴逆时针
        """
        super().__init__(shape)

        self.S = solar_constant
        self.h = orbit_height
        self.TT = orbit_period
        self.i_s = orbit_inclination

        self.Re = 6371.0
        self.pi = np.pi

        if surface_normal_angles is None:
            self.w0 = np.full(shape, self.pi / 2)
            self.l0 = np.zeros(shape)
        elif np.isscalar(surface_normal_angles):
            self.w0 = np.full(shape, surface_normal_angles)
            self.l0 = np.zeros(shape)
        else:
            if len(surface_normal_angles) != 2:
                raise ValueError("surface_normal_angles must be (w0, l0) tuple or scalar")
            w0_arr, l0_arr = surface_normal_angles
            if np.isscalar(w0_arr):
                self.w0 = np.full(shape, w0_arr)
            else:
                self.w0 = np.array(w0_arr, dtype=float)
            if np.isscalar(l0_arr):
                self.l0 = np.full(shape, l0_arr)
            else:
                self.l0 = np.array(l0_arr, dtype=float)

    @staticmethod
    def _clip_unit_interval(value: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """抑制浮点误差带来的 arccos/asin 定义域越界。"""
        return np.clip(value, -1.0, 1.0)

    def _compute_solar_angle(self, time: float) -> float:
        """
        计算任意时刻的太阳入射角

        :param time: 仿真时间 [s]
        :return: 太阳入射因子 X (0~1)
        """
        i_j = time / self.TT * 2 * self.pi
        while i_j > 2 * self.pi:
            i_j -= 2 * self.pi

        i_e0 = -self.i_s * self.pi / 180.0
        ie_s = np.arccos(self._clip_unit_interval(np.cos(i_e0) * np.cos(i_j)))

        w_iz = 0.0
        l_iz = 0.0

        if abs(i_e0) < 1e-10:
            w_iz = 0.0
            l_iz = self.pi + i_j
        else:
            tan_i_e0 = np.tan(i_e0)
            if abs(tan_i_e0) < 1e-15:
                w_iz = 0.0
                l_iz = self.pi + i_j
            else:
                w_iz = np.arcsin(np.sin(i_e0) * np.sin(self.pi + i_j))
                tempa = np.arcsin(self._clip_unit_interval(np.tan(w_iz) / tan_i_e0))
                tempb = np.arccos(
                    self._clip_unit_interval(np.cos(self.pi + i_j) / np.cos(w_iz))
                )
                if tempb * np.cos(tempa) > 0:
                    l_iz = tempa
                else:
                    l_iz = self.pi - tempa

                if l_iz < 0:
                    l_iz += 2.0 * self.pi

        iz_s = np.arccos(
            self._clip_unit_interval(np.cos(w_iz) * np.cos(l_iz - self.pi / 2.0))
        )
        sin_iz_s = np.sin(iz_s)
        with np.errstate(divide='ignore', invalid='ignore'):
            l_s = np.arccos(self._clip_unit_interval(-np.cos(ie_s) / sin_iz_s))
        if self.i_s < 0:
            l_s = abs(l_s)
        else:
            l_s = -abs(l_s)
        if (abs(self.i_s) < 1e-10 or abs(self.i_s - 180.0) < 1e-10) and np.isnan(l_s):
            if i_j >= self.pi / 2 and i_j < self.pi:
                l_s = 0.0
            else:
                l_s = self.pi

        X = np.cos(iz_s) * np.sin(self.w0) + np.sin(iz_s) * np.cos(self.w0) * np.cos(self.l0 - l_s)

        Sha = 0.0
        sin_i_s = np.sin(self.i_s * self.pi / 180.0)
        cos_i_s = np.cos(self.i_s * self.pi / 180.0)

        if abs(sin_i_s) < 1e-15:
            Sha = 0.0
        elif abs(cos_i_s) < 1e-15:
            Sha = 0.0
        elif (self.Re + self.h) >= (self.Re / abs(sin_i_s)):
            Sha = 0.0
        else:
            Sha = np.arccos(
                self._clip_unit_interval(
                    np.sqrt(((1 + self.h / self.Re) ** 2 - 1) /
                            ((1 + self.h / self.Re) ** 2 * cos_i_s ** 2))
                )
            )
            if i_j > (self.pi - Sha) and i_j < (self.pi + Sha):
                X = 0.0

        X = np.where(X < 0, 0.0, X)

        return X

    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        获取当前时刻的热流密度分布

        :param time: 当前仿真时间 [s]
        :return: 热流密度 [W/m²]
        """
        X = self._compute_solar_angle(time)
        flux = self.S * X
        return self._broadcast_flux(flux)

    def update_params(self,
                      orbit_height: float = None,
                      orbit_period: float = None,
                      orbit_inclination: float = None,
                      solar_constant: float = None,
                      surface_normal_angles: Tuple[float, float] = None):
        """
        更新轨道参数
        """
        if orbit_height is not None:
            self.h = orbit_height
        if orbit_period is not None:
            self.TT = orbit_period
        if orbit_inclination is not None:
            self.i_s = orbit_inclination
        if solar_constant is not None:
            self.S = solar_constant
        if surface_normal_angles is not None:
            w0, l0 = surface_normal_angles
            if np.isscalar(w0):
                self.w0.fill(w0)
            else:
                self.w0[:] = w0
            if np.isscalar(l0):
                self.l0.fill(l0)
            else:
                self.l0[:] = l0


class OrbitalIntegralHeatSource(BaseExternalHeatSource):
    """
    Reserved interface for strict orbital integral models.

    This class keeps the extension point for the full analytical/integral
    treatment described in the theory notes.  The solar direct model is
    already implemented in ``OrbitalHeatSource``.  Earth albedo and Earth IR
    can later be filled in by subclassing this interface without changing the
    BoundaryRegion integration path.
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 model_name: str,
                 orbit_height: float = 800.0,
                 orbit_period: float = 7644.0,
                 orbit_inclination: float = 0.0,
                 surface_normal_angles: Union[float, np.ndarray] = None):
        super().__init__(shape)
        self.model_name = model_name
        self.orbit_height = orbit_height
        self.orbit_period = orbit_period
        self.orbit_inclination = orbit_inclination

        if surface_normal_angles is None:
            self.surface_normal_angles = None
        else:
            self.surface_normal_angles = surface_normal_angles

    def get_heat_flux(self, time: float) -> np.ndarray:
        raise NotImplementedError(
            f"{self.model_name} is reserved for the strict orbital integral model "
            "and has not been implemented yet."
        )

    def update_params(self, **kwargs):
        if "orbit_height" in kwargs and kwargs["orbit_height"] is not None:
            self.orbit_height = kwargs["orbit_height"]
        if "orbit_period" in kwargs and kwargs["orbit_period"] is not None:
            self.orbit_period = kwargs["orbit_period"]
        if "orbit_inclination" in kwargs and kwargs["orbit_inclination"] is not None:
            self.orbit_inclination = kwargs["orbit_inclination"]
        if "surface_normal_angles" in kwargs and kwargs["surface_normal_angles"] is not None:
            self.surface_normal_angles = kwargs["surface_normal_angles"]


class IntegralAlbedoHeatSource(OrbitalIntegralHeatSource):
    """Reserved placeholder for the strict Earth albedo integral model."""

    def __init__(self, shape: Tuple[int, ...], **kwargs):
        super().__init__(shape=shape, model_name="IntegralAlbedoHeatSource", **kwargs)


class IntegralEarthIRHeatSource(OrbitalIntegralHeatSource):
    """Reserved placeholder for the strict Earth IR integral model."""

    def __init__(self, shape: Tuple[int, ...], **kwargs):
        super().__init__(shape=shape, model_name="IntegralEarthIRHeatSource", **kwargs)


class OrbitalTableHeatSource(BaseExternalHeatSource):
    """
    Legacy-Fortran-style orbital heat source driven by embedded lookup tables.

    The default table library reproduces the old workflow:
    1. choose a heat-table id for each boundary segment
    2. wrap time by the table period
    3. linearly interpolate the current heat flux

    The table data are stored directly inside the Python package, so runtime
    file reading is no longer required.
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 table_ids: Union[int, np.ndarray, list, tuple],
                 table_library: EmbeddedFluxTableLibrary = FORTRAN_ORBITAL_HEAT_TABLE_LIBRARY,
                 scale_factor: float = 1.0,
                 offset: float = 0.0,
                 periodic: bool = True):
        super().__init__(shape)
        self.table_library = table_library
        self.scale_factor = float(scale_factor)
        self.offset = float(offset)
        self.periodic = bool(periodic)
        self.table_ids = self._normalize_table_ids(table_ids)

    def _normalize_table_ids(self, table_ids: Union[int, np.ndarray, list, tuple]) -> np.ndarray:
        if np.isscalar(table_ids):
            table_array = np.full(self.shape, int(table_ids), dtype=int)
        else:
            table_array = np.array(table_ids, dtype=int)
            if table_array.shape != self.shape:
                table_array = np.broadcast_to(table_array, self.shape)
                table_array = np.array(table_array, dtype=int, copy=True)

        for table_id in np.unique(table_array):
            self.table_library.get_table(int(table_id))
        return table_array

    @staticmethod
    def _wrap_time(time: float, period: float) -> float:
        wrapped = float(time)
        if period <= 0.0:
            return wrapped
        if wrapped > period:
            wrapped = np.fmod(wrapped, period)
            if np.isclose(wrapped, 0.0):
                wrapped = period
        elif wrapped < 0.0:
            wrapped = np.fmod(wrapped, period)
            if wrapped < 0.0:
                wrapped += period
        return wrapped

    def _interpolate_table(self, table_id: int, time: float) -> float:
        table = self.table_library.get_table(int(table_id))
        sample_time = float(time)
        if self.periodic and table.periodic:
            sample_time = self._wrap_time(sample_time, float(table.time[-1]))
        return float(np.interp(sample_time, table.time, table.values))

    def get_heat_flux(self, time: float) -> np.ndarray:
        flux = np.empty(self.shape, dtype=float)
        flux_flat = flux.reshape(-1)
        for index, table_id in enumerate(self.table_ids.reshape(-1)):
            flux_flat[index] = self._interpolate_table(int(table_id), time)
        flux *= self.scale_factor
        if self.offset != 0.0:
            flux += self.offset
        return flux

    def update_params(self,
                      table_ids: Union[int, np.ndarray, list, tuple] = None,
                      scale_factor: float = None,
                      offset: float = None,
                      periodic: bool = None,
                      table_library: EmbeddedFluxTableLibrary = None):
        if table_library is not None:
            self.table_library = table_library
            self.table_ids = self._normalize_table_ids(self.table_ids)
        if table_ids is not None:
            self.table_ids = self._normalize_table_ids(table_ids)
        if scale_factor is not None:
            self.scale_factor = float(scale_factor)
        if offset is not None:
            self.offset = float(offset)
        if periodic is not None:
            self.periodic = bool(periodic)



class OrbitalMatrixHeatSource(BaseExternalHeatSource):
    """Matrix lookup orbital heat source driven by packaged w0=8.12deg sum data."""

    def __init__(self,
                 shape: Tuple[int, ...],
                 matrix_key: str,
                 matrix_library: EmbeddedFluxMatrixLibrary = W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY,
                 scale_factor: float = 1.0,
                 offset: float = 0.0,
                 periodic: bool = True):
        super().__init__(shape)
        self.matrix_library = matrix_library
        self.matrix_key = str(matrix_key)
        self.matrix = None
        self.scale_factor = float(scale_factor)
        self.offset = float(offset)
        self.periodic = bool(periodic)
        self._load_and_validate_matrix()

    def _load_and_validate_matrix(self):
        self.matrix = self.matrix_library.get_matrix(self.matrix_key)
        expected = int(np.prod(self.shape))
        actual = int(self.matrix.values.shape[1])
        if actual != expected:
            raise ValueError(
                f"Orbital heat matrix {self.matrix_key!r} column count {actual} does not match heat source shape {self.shape} ({expected} nodes)."
            )

    @staticmethod
    def _wrap_time_to_matrix(time: float, sample_time: np.ndarray) -> float:
        start = float(sample_time[0])
        end = float(sample_time[-1])
        step = float(sample_time[1] - sample_time[0]) if sample_time.size > 1 else 0.0
        period = (end - start) + step
        if period <= 0.0:
            return float(time)
        return ((float(time) - start) % period) + start

    def get_heat_flux(self, time: float) -> np.ndarray:
        sample_time = float(time)
        if self.periodic and self.matrix.periodic:
            sample_time = self._wrap_time_to_matrix(sample_time, self.matrix.time)
        values = np.empty(self.matrix.values.shape[1], dtype=float)
        for index in range(self.matrix.values.shape[1]):
            values[index] = np.interp(sample_time, self.matrix.time, self.matrix.values[:, index])
        values *= self.scale_factor
        if self.offset != 0.0:
            values += self.offset
        return values.reshape(self.shape)

    def update_params(self,
                      matrix_key: str = None,
                      scale_factor: float = None,
                      offset: float = None,
                      periodic: bool = None,
                      matrix_library: EmbeddedFluxMatrixLibrary = None):
        reload_matrix = False
        if matrix_library is not None:
            self.matrix_library = matrix_library
            reload_matrix = True
        if matrix_key is not None:
            self.matrix_key = str(matrix_key)
            reload_matrix = True
        if reload_matrix:
            self._load_and_validate_matrix()
        if scale_factor is not None:
            self.scale_factor = float(scale_factor)
        if offset is not None:
            self.offset = float(offset)
        if periodic is not None:
            self.periodic = bool(periodic)


class AlbedoHeatSource(BaseExternalHeatSource):
    """
    地球反照率热源

    物理模型:
    - 地球反射的太阳辐射，通常为太阳常数的 30% 左右
    - 反照率因子通常取 0.3 (30%)

    用于模拟近地轨道飞行器受到地球反射太阳辐射加热的情况。
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 albedo_factor: float = 0.3,
                 solar_constant: float = 1361.0,
                 earth_temperature: float = 255.0):
        """
        :param shape: 边界形状
        :param albedo_factor: 地球反照率 (0~1)，默认 0.3
        :param solar_constant: 太阳常数 [W/m²]
        :param earth_temperature: 地球等效辐射温度 [K]
        """
        super().__init__(shape)

        self.albedo_factor = albedo_factor
        self.S = solar_constant
        self.T_earth = earth_temperature
        self.sigma = 5.670374419e-8

        self.q_albedo = self.albedo_factor * self.S

    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        反照率热流通常不随时间快速变化（除非轨道阴影变化）
        这里简化为恒定热流

        :param time: 当前仿真时间 [s] (未使用)
        :return: 热流密度 [W/m²]
        """
        return self._broadcast_flux(self.q_albedo)

    def update_params(self, albedo_factor: float = None, solar_constant: float = None):
        """
        更新反照率参数
        """
        if albedo_factor is not None:
            self.albedo_factor = albedo_factor
            self.q_albedo = self.albedo_factor * self.S
        if solar_constant is not None:
            self.S = solar_constant
            self.q_albedo = self.albedo_factor * self.S


class EarthIRHeatSource(BaseExternalHeatSource):
    """
    地球红外热源

    物理模型:
    - 地球大气层向太空辐射红外热量
    - 近地轨道飞行器会接收到这部分热流
    - 通常取值为 230 ~ 260 W/m²

    这是一个次要热源，通常与太阳辐射和反照率一起考虑。
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 earth_ir_flux: float = 237.0,
                 emissivity: float = 1.0):
        """
        :param shape: 边界形状
        :param earth_ir_flux: 地球红外热流密度 [W/m²]，默认 237 W/m²
        :param emissivity: 表面发射率
        """
        super().__init__(shape)

        self.q_ir = earth_ir_flux
        self.emissivity = emissivity

    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        地球红外热流近似为恒定值

        :param time: 当前仿真时间 [s] (未使用)
        :return: 热流密度 [W/m²]
        """
        return self._broadcast_flux(self.q_ir)

    def update_params(self, earth_ir_flux: float = None, emissivity: float = None):
        """
        更新红外热流参数
        """
        if earth_ir_flux is not None:
            self.q_ir = earth_ir_flux
        if emissivity is not None:
            self.emissivity = emissivity


class TimeVaryingHeatSource(BaseExternalHeatSource):
    """
    时间变化热源

    通用包装器，允许将任意时间函数包装为热源
    用于模拟:
    - 瞬态功率变化
    - 控制信号驱动的加热/冷却
    - 用户自定义的时间曲线
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 time_func=None,
                 amplitude: float = 1.0):
        """
        :param shape: 边界形状
        :param time_func: 时间函数，签名为 func(time) -> float 或 np.ndarray
        :param amplitude: 振幅缩放因子
        """
        super().__init__(shape)

        self.time_func = time_func
        self.amplitude = amplitude

    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        根据时间函数计算热流

        :param time: 当前仿真时间 [s]
        :return: 热流密度 [W/m²]
        """
        if self.time_func is None:
            return np.zeros(self.shape, dtype=float)
        flux = self.time_func(time)
        return self._broadcast_flux(flux * self.amplitude)

    def update_params(self, time_func=None, amplitude=None):
        """
        更新时间函数或振幅
        """
        if time_func is not None:
            self.time_func = time_func
        if amplitude is not None:
            self.amplitude = amplitude


class FixedPowerHeatSource(BaseExternalHeatSource):
    """
    固定功率热源

    提供恒定的热流边界条件，
    类似于 ReactorCore 中的敷设功率 (Wall Heat Flux) 概念。
    """

    def __init__(self,
                 shape: Tuple[int, ...],
                 q_flux: Union[float, np.ndarray]):
        """
        :param shape: 边界形状
        :param q_flux: 热流密度 [W/m²]，可以是标量或数组
        """
        super().__init__(shape)

        if np.isscalar(q_flux):
            self.q_fixed = np.full(shape, q_flux, dtype=float)
        else:
            self.q_fixed = np.array(q_flux, dtype=float)

    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        :param time: 当前仿真时间 [s] (未使用)
        :return: 恒定热流密度 [W/m²]
        """
        return self._broadcast_flux(self.q_fixed)

    def update_params(self, q_flux: Union[float, np.ndarray] = None):
        """
        更新热流值
        """
        if q_flux is not None:
            if np.isscalar(q_flux):
                self.q_fixed.fill(q_flux)
            else:
                self.q_fixed[:] = q_flux


class CompositeHeatSource(BaseExternalHeatSource):
    """
    组合热源

    将多个外热源组合成一个，支持叠加多个热流效果
    例如：太阳辐射 + 地球反照率 + 地球红外
    """

    def __init__(self, shape: Tuple[int, ...]):
        super().__init__(shape)
        self.sources = []

    def add_source(self, source: BaseExternalHeatSource):
        """添加一个热源"""
        if source.shape != self.shape:
            raise ValueError(f"Source shape {source.shape} mismatch with CompositeHeatSource shape {self.shape}")
        self.sources.append(source)

    def get_heat_flux(self, time: float) -> np.ndarray:
        """
        叠加所有热源的热流密度

        :param time: 当前仿真时间 [s]
        :return: 总热流密度 [W/m²]
        """
        if not self.sources:
            return np.zeros(self.shape, dtype=float)
        total = self.sources[0].get_heat_flux(time).copy()
        for source in self.sources[1:]:
            total += source.get_heat_flux(time)
        return total

    def update_params(self, **kwargs):
        """转发参数更新给所有热源"""
        for source in self.sources:
            source.update_params(**kwargs)
