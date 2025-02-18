"""Conditionally import numba based on the TCPYPI_DISABLE_NUMBA environment variable.

Set the TCPYPI_DISABLE_NUMBA environment variable to 1 to disable numba.

This is useful for coverage testing.
"""

import os
from functools import wraps


def noop_njit(*args, **kwargs):
    """No-op decorator to replace @njit when numba is disabled."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapper

    # Handle both @njit and @njit() syntax
    if len(args) == 1 and callable(args[0]):
        return decorator(args[0])
    return decorator


if os.getenv("TCPYPI_DISABLE_NUMBA") == "1":
    njit = noop_njit
else:
    import numba as nb

    njit = nb.njit
