"""Conditionally import numba based on the DISABLE_NUMBA environment variable.

Set the DISABLE_NUMBA environment variable to 1 to disable numba.

This is useful for coverage testing.
"""

import os
from functools import wraps

if os.getenv('DISABLE_NUMBA') == '1':
    # Create a no-op decorator to use when numba is disabled
    def njit(*args, **kwargs):
        def decorator(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return wrapper
        # Handle both @njit and @njit() syntax
        if len(args) == 1 and callable(args[0]):
            return decorator(args[0])
        return decorator
else:
    import numba as nb
    njit = nb.njit
