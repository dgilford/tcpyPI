"""Conditionally import numba based on the TCPYPI_DISABLE_NUMBA environment variable.

Set the TCPYPI_DISABLE_NUMBA environment variable to 1 to disable numba.

This is useful for coverage testing.
"""

import os
from functools import wraps

import numpy as np


def noop_njit(*args, **kwargs):
    """No-op decorator to replace @njit when numba is disabled.

    Note that this decorator isn't entirely a noop, because numba does some casting
    of NumPy types to Python types. Since NumPy v2, NumPy types like np.float64 are
    printed explicitly, so the doctests would fail without casting scalars back
    to Python types.

    See <https://github.com/numba/numba/issues/704> for more details about Numba's
    behavior.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            if isinstance(result, np.floating):
                result = float(result)
            elif isinstance(result, tuple):
                # Cast components to float
                result = tuple(float(x) if isinstance(x, np.floating) else x for x in result)
            return result

        return wrapper

    # Handle both @njit and @njit() syntax
    if len(args) == 1 and callable(args[0]):
        return decorator(args[0])
    return decorator


if os.getenv("TCPYPI_DISABLE_NUMBA") == "1":
    _numba_available = False
else:
    try:
        import numba as nb

        _numba_available = True
    except ImportError:
        # Graceful degradation: tcpyPI remains importable and correct without
        # Numba, just slow (roadmap issue #75, "make numba a lazy dependency").
        import warnings

        warnings.warn(
            "Numba is not installed; tcpyPI will run in pure-Python mode "
            "(orders of magnitude slower). Install it with `pip install numba` "
            "to enable the compiled kernels.",
            RuntimeWarning,
            stacklevel=2,
        )
        _numba_available = False

if not _numba_available:
    njit = noop_njit
    # No pure-Python generalized-ufunc equivalent; consumers (pi_field) fall back
    # to a per-column loop when this is None.
    guvectorize = None
else:
    guvectorize = nb.guvectorize

    def njit(*args, **kwargs):
        """``numba.njit`` with on-disk caching enabled by default.

        ``cache=True`` avoids recompiling the PI/CAPE kernels in every fresh Python
        process (first-call JIT is otherwise repeated on each run). Override per call
        with ``@njit(cache=False)`` if a specific function is not cacheable.

        Do NOT enable ``fastmath`` here or on any tcpyPI kernel: the algorithm relies
        on IEEE NaN semantics (e.g. ``NaN > 0`` is False, and ``np.min``/``max`` NaN
        propagation) for its missing-value handling, which ``fastmath`` would silently
        break.
        """
        kwargs.setdefault("cache", True)
        if len(args) == 1 and callable(args[0]) and not isinstance(args[0], str):
            # Bare ``@njit`` usage: apply options via a second call.
            return nb.njit(**kwargs)(args[0])
        return nb.njit(*args, **kwargs)
