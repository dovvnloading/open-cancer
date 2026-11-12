"""Bounded rechecking of published releases during Pages propagation/outages.

Only a complete, successful integrity check counts as healthy. Every failed attempt
is retried from the beginning, so mixed-version asset responses cannot pass.
"""
from __future__ import annotations

from collections.abc import Callable
import time


def wait_for_healthy_release(
    check: Callable[[], list[str]], *, wait_seconds: float = 0,
    interval_seconds: float = 5, clock: Callable[[], float] = time.monotonic,
    sleep: Callable[[float], None] = time.sleep,
    on_retry: Callable[[list[str], float], None] | None = None,
) -> list[str]:
    """Retry a *complete* verification until it passes or the deadline expires."""
    if wait_seconds < 0 or interval_seconds <= 0:
        raise ValueError("Expected nonnegative wait and positive polling interval")
    deadline = clock() + wait_seconds
    while True:
        errors = check()
        if not errors:
            return []
        remaining = deadline - clock()
        if remaining <= 0:
            return errors
        delay = min(interval_seconds, remaining)
        if on_retry is not None:
            on_retry(errors, delay)
        sleep(delay)
