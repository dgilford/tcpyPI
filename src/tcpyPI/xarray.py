"""Whole-field potential intensity on xarray Datasets.

This module packages the gridded-workflow entry point that previously lived in
the repository's ``run_sample.py`` script: one call computes PI over an entire
``xarray.Dataset`` through :func:`tcpyPI.pi_field` (a compiled, multithreaded
generalized ufunc — no per-column Python loop).

Requires the optional xarray extra (``pip install tcpypi[xarray]``). The module
itself imports cleanly without xarray installed (mirroring the graceful numba
fallback in ``tcpyPI.numba``) so that test collection works in minimal
environments; calling :func:`pi_dataset` then raises an informative error.
"""

from __future__ import annotations

try:
    import xarray as xr
except ImportError:  # pragma: no cover - exercised only in minimal environments
    xr = None

from .pi import pi_field

_REQUIRED_VARS = ("sst", "msl", "t", "r")


def pi_dataset(
    ds,
    dim="p",
    *,
    CKCD=0.9,
    ascent_flag=0,
    diss_flag=1,
    V_reduc=0.8,
    ptop=50,
    miss_handle=1,
    outflow_source="cape_star",
):
    """Calculate potential intensity over a gridded dataset in one whole-field call.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset containing ``sst`` (sea surface temperature, degrees C),
        ``msl`` (mean sea level pressure, hPa), ``t`` (temperature profiles,
        degrees C), and ``r`` (mixing ratio profiles, g/kg), with the vertical
        pressure coordinate named by ``dim`` (hPa). Same units contract as
        :func:`tcpyPI.pi`.
    dim : str, default="p"
        Name of the vertical pressure coordinate in ``ds``.
    CKCD : float, default=0.9
        Ratio of exchange coefficients (Ck/Cd).
    ascent_flag : int, default=0
        0 = reversible ascent, 1 = pseudo-adiabatic ascent.
    diss_flag : int, default=1
        1 = include dissipative heating, 0 = exclude it.
    V_reduc : float, default=0.8
        Gradient-to-surface wind reduction factor applied to VMAX.
    ptop : float, default=50
        Pressure (hPa) above which the profile is ignored.
    miss_handle : int, default=1
        Missing-value handling in the CAPE calculation (see :func:`tcpyPI.pi`).
    outflow_source : {"cape_star", "cape_env"}, default="cape_star"
        Which CAPE calculation supplies the outflow temperature and level.

    Returns
    -------
    xarray.Dataset
        PI outputs ``vmax`` (m/s), ``pmin`` (hPa), ``ifl`` (status flag),
        ``t0`` (K), and ``otl`` (hPa) on the input horizontal grid, with
        ``standard_name``/``units`` attributes set.

        Flags tripped by the environmental-CAPE call persist, so ``ifl`` values
        of 0 or 3 can accompany finite outputs — ``ifl == 1`` is the
        trustworthy-output gate.
    """
    if xr is None:
        raise ImportError(
            "pi_dataset requires xarray; install the optional extra with "
            "`pip install tcpypi[xarray]`"
        )

    missing = [v for v in _REQUIRED_VARS if v not in ds] + ([] if dim in ds else [dim])
    if missing:
        raise ValueError(
            f"input dataset is missing required variable(s) {missing}; "
            f"expected variables {list(_REQUIRED_VARS)} and vertical coordinate {dim!r}"
        )

    vmax, pmin, ifl, t0, otl = xr.apply_ufunc(
        pi_field,
        ds["sst"],
        ds["msl"],
        ds[dim],
        ds["t"],
        ds["r"],
        kwargs=dict(
            CKCD=CKCD,
            ascent_flag=ascent_flag,
            diss_flag=diss_flag,
            V_reduc=V_reduc,
            ptop=ptop,
            miss_handle=miss_handle,
            outflow_source=outflow_source,
        ),
        input_core_dims=[[], [], [dim], [dim], [dim]],
        output_core_dims=[[], [], [], [], []],
    )

    out_ds = xr.Dataset({"vmax": vmax, "pmin": pmin, "ifl": ifl, "t0": t0, "otl": otl})

    out_ds.vmax.attrs["standard_name"], out_ds.vmax.attrs["units"] = (
        "Maximum Potential Intensity",
        "m/s",
    )
    out_ds.pmin.attrs["standard_name"], out_ds.pmin.attrs["units"] = (
        "Minimum Central Pressure",
        "hPa",
    )
    out_ds.ifl.attrs["standard_name"] = "pyPI Flag"
    out_ds.t0.attrs["standard_name"], out_ds.t0.attrs["units"] = "Outflow Temperature", "K"
    out_ds.otl.attrs["standard_name"], out_ds.otl.attrs["units"] = (
        "Outflow Temperature Level",
        "hPa",
    )

    return out_ds
