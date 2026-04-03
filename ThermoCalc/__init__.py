# 在项目根目录或任何测试文件中
from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

# 如果遇到 te_solver 找不到的问题，有时需要在 __init__.py 中做相对路径处理，
# 或者在运行 TASTIN 时，将 TASTIN 根目录加入 PYTHONPATH。
