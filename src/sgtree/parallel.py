from __future__ import annotations

import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, TypeVar


T = TypeVar("T")
R = TypeVar("R")


def _bounded_workers(workers: int, n_tasks: int) -> int:
    if n_tasks <= 0:
        return 0
    return max(1, min(int(workers), int(n_tasks)))


def map_threaded(func: Callable[[T], R], args: list[T], workers: int) -> list[R]:
    if not args:
        return []
    n_workers = _bounded_workers(workers, len(args))
    if n_workers <= 1:
        return [func(item) for item in args]
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        return list(executor.map(func, args))


def map_processed(func: Callable[[T], R], args: list[T], workers: int) -> list[R]:
    if not args:
        return []
    n_workers = _bounded_workers(workers, len(args))
    if n_workers <= 1:
        return [func(item) for item in args]
    try:
        with mp.Pool(n_workers) as pool:
            return pool.map(func, args)
    except (PermissionError, OSError) as exc:
        print(f"warning: multiprocessing unavailable ({exc}); falling back to serial execution")
        return [func(item) for item in args]
