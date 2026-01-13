"""Analysis helpers for tcpyPI.

These functions provide a small, user-facing API for decomposing potential intensity
into logarithmic components following Wing et al. (2015), Eq. 2.

Notes
-----
- The decomposition implemented here matches `tcpyPI.utilities.decompose_pi`:
  `lnpi` is `ln(V^2)` (i.e., `2 * ln(V)`), not `ln(V)`.
- Inputs can be scalars or NumPy arrays; invalid physical states yield `nan`.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple, Union

import numpy as np


_SSTUnits = str


def _maybe_scalar(value: Union[np.ndarray, float]) -> Any:
    return float(value) if np.ndim(value) == 0 else value


def _sst_to_kelvin(sst: Any, sst_units: _SSTUnits) -> Any:
    units = sst_units.upper()
    if units == "K":
        return sst
    if units == "C":
        return np.asarray(sst) + 273.15
    raise ValueError(f"Unsupported sst_units={sst_units!r}; expected 'C' or 'K'.")


def log_decompose_pi(
    pi: Any,
    sst: Any,
    t0: Any,
    CKCD: float = 0.9,
    *,
    sst_units: _SSTUnits = "K",
) -> Tuple[Any, Any, Any, float]:
    """Log-decompose potential intensity into efficiency, disequilibrium, and Ck/Cd.

    Parameters
    ----------
    pi
        Potential intensity wind speed (m/s).
    sst
        Sea surface temperature, in units given by `sst_units`.
    t0
        Outflow temperature (K).
    CKCD
        Ratio of exchange coefficients (Ck/Cd).
    sst_units
        Units for `sst`: `'K'` (default) or `'C'`.

    Returns
    -------
    (lnpi, lneff, lndiseq, lnCKCD)
        `lnpi` is `ln(pi^2)`.

    Notes
    -----
    This matches the behavior of `tcpyPI.utilities.decompose_pi`, including edge cases:
    - If efficiency <= 0: all terms are `nan` except `lnCKCD`.
    - If efficiency > 0 and pi <= 0: `lneff` is returned, but `lnpi` and `lndiseq` are `nan`.
    """

    pi_arr = np.asarray(pi, dtype=float)
    t0_arr = np.asarray(t0, dtype=float)
    sst_k = np.asarray(_sst_to_kelvin(sst, sst_units), dtype=float)

    # Fast path: preserve the original Numba-compiled scalar behavior.
    if pi_arr.ndim == 0 and t0_arr.ndim == 0 and sst_k.ndim == 0:
        from .utilities import decompose_pi as _decompose_pi

        return _decompose_pi(float(pi_arr), float(sst_k), float(t0_arr), CKCD=CKCD)

    pi_arr, sst_k, t0_arr = np.broadcast_arrays(pi_arr, sst_k, t0_arr)

    lnCKCD = float(np.log(CKCD))

    efficiency = (sst_k - t0_arr) / t0_arr
    valid_eff = efficiency > 0
    valid_pi = pi_arr > 0

    lneff = np.full(efficiency.shape, np.nan, dtype=float)
    lneff[valid_eff] = np.log(efficiency[valid_eff])

    lnpi = np.full(pi_arr.shape, np.nan, dtype=float)
    lnpi[valid_pi] = 2.0 * np.log(pi_arr[valid_pi])

    lndiseq = np.full(pi_arr.shape, np.nan, dtype=float)
    valid = valid_eff & valid_pi
    lndiseq[valid] = lnpi[valid] - lneff[valid] - lnCKCD

    return (
        _maybe_scalar(lnpi),
        _maybe_scalar(lneff),
        _maybe_scalar(lndiseq),
        lnCKCD,
    )


def pi_log_decomposition(
    SSTC: Any,
    MSL: Any,
    P: Any,
    TC: Any,
    R: Any,
    CKCD: float = 0.9,
    *,
    ascent_flag: int = 0,
    diss_flag: int = 1,
    V_reduc: float = 0.8,
    ptop: float = 50,
    miss_handle: int = 1,
) -> Dict[str, Any]:
    """Run `tcpyPI.pi` and return outputs plus the Wing et al. (2015) log decomposition.

    Parameters match `tcpyPI.pi` (SST in Celsius).

    Returns
    -------
    dict
        Keys: `vmax`, `pmin`, `ifl`, `t0`, `otl`, `lnpi`, `lneff`, `lndiseq`, `lnCKCD`.
    """

    from .pi import pi as _pi

    vmax, pmin, ifl, t0, otl = _pi(
        SSTC,
        MSL,
        P,
        TC,
        R,
        CKCD=CKCD,
        ascent_flag=ascent_flag,
        diss_flag=diss_flag,
        V_reduc=V_reduc,
        ptop=ptop,
        miss_handle=miss_handle,
    )

    lnpi, lneff, lndiseq, lnCKCD = log_decompose_pi(
        vmax,
        SSTC,
        t0,
        CKCD=CKCD,
        sst_units="C",
    )

    return {
        "vmax": vmax,
        "pmin": pmin,
        "ifl": ifl,
        "t0": t0,
        "otl": otl,
        "lnpi": lnpi,
        "lneff": lneff,
        "lndiseq": lndiseq,
        "lnCKCD": lnCKCD,
    }
