from functools import wraps
import logging
import time


main_logger = logging.getLogger("time_check")  # logger
main_logger.setLevel(logging.INFO)
main_logger.propagate = False

def configure_logging(name: str):
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(message)s")
    handler.setFormatter(formatter)
    main_logger.addHandler(handler)
    main_logger.addHandler(logging.FileHandler(f"experiments/results/{name}.log"))
    main_logger.info(f"Logging started at {time.strftime('%Y-%m-%d %H:%M:%S')}")



def time_check(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        main_logger.info(f"{func.__qualname__:60} in {end-start:>7.3f}" + "s")
        return result

    return wrapper
