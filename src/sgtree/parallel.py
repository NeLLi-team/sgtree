"""Bounded thread and process mapping for independent pipeline tasks."""

from __future__ import annotations

import logging
import multiprocessing as mp
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import TypeVar

logger = logging.getLogger(__name__)

T = TypeVar("T")
R = TypeVar("R")


def _bounded_workers(workers: int, n_tasks: int) -> int:
    if n_tasks <= 0:
        return 0
    return max(1, min(int(workers), int(n_tasks)))


def map_threaded(func: Callable[[T], R], args: list[T], workers: int) -> list[R]:
    """Map with no more threads than tasks, using serial execution for one worker."""
    if not args:
        return []
    n_workers = _bounded_workers(workers, len(args))
    if n_workers <= 1:
        return [func(item) for item in args]
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        return list(executor.map(func, args))


def map_processed(func: Callable[[T], R], args: list[T], workers: int) -> list[R]:
    """Map with bounded processes, falling back if a pool cannot be created."""
    if not args:
        return []
    n_workers = _bounded_workers(workers, len(args))
    if n_workers <= 1:
        return [func(item) for item in args]
    try:
        pool = mp.Pool(n_workers)
    except (PermissionError, OSError) as exc:
        logger.warning(
            "multiprocessing unavailable (%s); falling back to serial execution", exc
        )
        return [func(item) for item in args]
    with pool:
        return pool.map(func, args)
