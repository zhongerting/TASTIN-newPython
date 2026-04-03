import time
import functools

class TEASAProfiler:
    """
    TEASA 轻量级性能剖析器
    用于统计函数的调用次数和累积执行时间。
    """
    stats = {}

    @classmethod
    def profile(cls, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.perf_counter()  # 记录开始时间 (高精度)
            result = func(*args, **kwargs)    # 执行原函数
            end_time = time.perf_counter()    # 记录结束时间

            elapsed = end_time - start_time
            name = func.__qualname__          # 获取类名.方法名 (如 HydraulicNetwork._calc_momentum_coeffs)

            if name not in cls.stats:
                cls.stats[name] = {'count': 0, 'time': 0.0}

            cls.stats[name]['count'] += 1
            cls.stats[name]['time'] += elapsed

            return result
        return wrapper

    @classmethod
    def report(cls):
        """打印最终的性能统计报告"""
        print("\n" + "="*60)
        print(f"{'Function Name':<40} | {'Calls':<8} | {'Total Time (s)':<12}")
        print("-" * 60)
        # 按总耗时降序排列
        sorted_stats = sorted(cls.stats.items(), key=lambda x: x[1]['time'], reverse=True)
        for name, data in sorted_stats:
            print(f"{name:<40} | {data['count']:<8} | {data['time']:<12.6f}")
        print("="*60 + "\n")