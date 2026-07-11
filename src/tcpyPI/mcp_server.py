"""tcpyPI MCP server: conversational access to validated potential-intensity code.

Exposes five tools over the Model Context Protocol (stdio transport):

1. ``compute_pi`` — single sounding -> VMAX, PMIN, IFL, TO, OTL.
2. ``compute_pi_grid`` — netCDF path in -> netCDF path out + JSON summary.
3. ``decompose_pi`` — Wing et al. (2015) log decomposition of PI.
4. ``ventilation_index`` — Tang & Emanuel (2012) ventilation index.
5. ``genesis_potential_index`` — GPI (Emanuel & Nolan 2004 or Emanuel 2010).

Design contract (see AGENTS.md): the server is I/O and validation only — every
number comes from the tested `tcpyPI` functions. Units are explicit and never
silently converted; gridded inputs must carry matching netCDF ``units``
attributes or the call fails with a structured error. Every result carries
provenance (package version, DOI, citations, knob settings) and the IFL status
contract: flags tripped by the environmental-CAPE call persist, so ``ifl`` of
0 or 3 can accompany finite outputs — ``ifl == 1`` is the only
trustworthy-output gate.

Run with ``tcpypi-mcp`` (requires ``pip install tcpypi[mcp]``). The module
imports without fastmcp/xarray installed; the server entry point and the grid
tool then raise informative errors.
"""

from __future__ import annotations

import math
import os
from importlib.metadata import PackageNotFoundError, version

import numpy as np

from . import numba as _numba_shim
from .gpi import genesis_potential_index as _genesis_potential_index
from .pi import pi as _pi
from .pi import pi_log_decomposition as _pi_log_decomposition
from .vi import ventilation_index as _ventilation_index

_DOI = "10.5281/zenodo.3756005"  # concept DOI (all versions)
_CITATIONS = [
    "Bister, M. and Emanuel, K. A. (2002): doi:10.1029/2001JD000776",
    "Gilford, D. M. (2021), Geosci. Model Dev. 14, 2351-2369: doi:10.5194/gmd-14-2351-2021",
]

_IFL_MEANING = {
    0: "improper sounding or parcel (e.g. missing SST); outputs not trustworthy",
    1: "converged; trustworthy output",
    2: "CAPE routine failed to converge; outputs not trustworthy",
    3: "missing data in the profile or pressures; outputs not trustworthy",
}

# Inline soundings above this size are refused and directed to compute_pi_grid.
_MAX_INLINE_LEVELS = 5000

# Accepted netCDF `units` attribute spellings (normalized to lowercase). The
# grid tool REJECTS anything else rather than converting: silent attr-driven
# unit conversion is how plausible-but-wrong numbers happen.
_UNIT_SYNONYMS = {
    "degC": {
        "degc",
        "deg c",
        "deg_c",
        "degree_c",
        "degrees c",
        "degrees_c",
        "degree celsius",
        "degrees celsius",
        "celsius",
        "c",
        "°c",
    },
    "hPa": {
        "hpa",
        "hectopascal",
        "hectopascals",
        "mb",
        "mbar",
        "millibar",
        "millibars",
    },
    "g/kg": {
        "g/kg",
        "g kg-1",
        "g kg**-1",
        "g kg^-1",
        "g.kg-1",
        "grams per kilogram",
        "gram/kilogram",
        "grams/kilogram",
    },
}


def _tcpypi_version():
    try:
        return version("tcpypi")
    except PackageNotFoundError:  # pragma: no cover - only for odd installs
        return "unknown"


def _provenance(**settings):
    """Provenance block attached to every tool result."""
    return {
        "tcpypi_version": _tcpypi_version(),
        "doi": _DOI,
        "citations": _CITATIONS,
        "numba": _numba_shim.guvectorize is not None,
        "settings": settings,
    }


def _error(message, **extra):
    return {"status": "error", "message": message, **extra}


def _num(x):
    """One float, JSON-safe: NaN/inf become None."""
    x = float(x)
    return x if math.isfinite(x) else None


def _nums(x):
    """Scalar-or-array -> JSON-safe float or nested list (NaN/inf -> None)."""
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 0:
        return _num(arr)
    return [_nums(v) for v in arr]


def _scalar_in(x):
    """JSON number-or-null -> float (null means missing; JSON cannot carry NaN)."""
    return np.nan if x is None else float(x)


def _array_in(seq):
    """JSON array of numbers-or-nulls -> float array (nulls become NaN)."""
    return np.asarray([np.nan if v is None else float(v) for v in seq], dtype=float)


def compute_pi(
    sst_c: float | None,
    msl_hpa: float | None,
    p_hpa: list[float | None],
    t_c: list[float | None],
    r_gkg: list[float | None],
    CKCD: float = 0.9,
    ascent_flag: int = 0,
    diss_flag: int = 1,
    V_reduc: float = 0.8,
    ptop: float = 50.0,
    miss_handle: int = 1,
    outflow_source: str = "cape_star",
) -> dict:
    """Bister & Emanuel (2002) tropical cyclone potential intensity for one sounding.

    Units (fixed, never converted): sst_c and t_c in degrees Celsius, msl_hpa
    and p_hpa in hPa, r_gkg (water vapor mixing ratio) in g/kg. The profile
    arrays must be the same length and may be in any vertical order. Missing
    mixing ratios may be zeros; the temperature profile should reach at least
    the tropopause. Missing values are passed as JSON null (JSON cannot carry
    NaN); their handling follows miss_handle and is reported through ifl.

    Returns vmax_ms (m/s, includes the V_reduc surface-wind reduction),
    pmin_hpa (hPa), ifl (status flag with ifl_meaning and trustworthy fields —
    only ifl == 1 output is trustworthy), t0_k (outflow temperature, K),
    otl_hpa (outflow level, hPa), and provenance. For gridded fields use
    compute_pi_grid instead of inlining large arrays.
    """
    n = len(p_hpa)
    if not (len(t_c) == n and len(r_gkg) == n):
        return _error(
            f"profile arrays must have equal lengths; got p_hpa={n}, "
            f"t_c={len(t_c)}, r_gkg={len(r_gkg)}"
        )
    if n > _MAX_INLINE_LEVELS:
        return _error(
            f"profile has {n} levels (> {_MAX_INLINE_LEVELS}); "
            "use compute_pi_grid for gridded/large inputs"
        )
    try:
        vmax, pmin, ifl, t0, otl = _pi(
            _scalar_in(sst_c),
            _scalar_in(msl_hpa),
            _array_in(p_hpa),
            _array_in(t_c),
            _array_in(r_gkg),
            CKCD=CKCD,
            ascent_flag=ascent_flag,
            diss_flag=diss_flag,
            V_reduc=V_reduc,
            ptop=ptop,
            miss_handle=miss_handle,
            outflow_source=outflow_source,
        )
    except Exception as exc:
        return _error(f"pi() failed: {exc!r}")
    ifl = int(ifl)
    return {
        "status": "success",
        "vmax_ms": _num(vmax),
        "pmin_hpa": _num(pmin),
        "ifl": ifl,
        "ifl_meaning": _IFL_MEANING.get(ifl, "unknown flag"),
        "trustworthy": ifl == 1,
        "t0_k": _num(t0),
        "otl_hpa": _num(otl),
        "provenance": _provenance(
            CKCD=CKCD,
            ascent_flag=ascent_flag,
            diss_flag=diss_flag,
            V_reduc=V_reduc,
            ptop=ptop,
            miss_handle=miss_handle,
            outflow_source=outflow_source,
        ),
    }


def _check_units(da, name, expected):
    """Return an error message unless `da` carries an accepted `units` attr."""
    units = da.attrs.get("units")
    if units is None:
        return (
            f"variable {name!r} has no `units` attribute; expected {expected}. "
            "Set the attribute explicitly — this server never guesses units."
        )
    if str(units).strip().lower() not in _UNIT_SYNONYMS[expected]:
        return (
            f"variable {name!r} has units {units!r} but {expected} is required. "
            "Convert the data (and update the attribute) before calling; "
            "this server never converts units silently."
        )
    return None


def compute_pi_grid(
    input_nc: str,
    output_nc: str,
    dim: str = "p",
    CKCD: float = 0.9,
    ascent_flag: int = 0,
    diss_flag: int = 1,
    V_reduc: float = 0.8,
    ptop: float = 50.0,
    miss_handle: int = 1,
    outflow_source: str = "cape_star",
) -> dict:
    """Compute potential intensity over a gridded netCDF file (whole-field, compiled).

    Reads input_nc, computes PI for every column in one compiled multithreaded
    call (tcpyPI.xarray.pi_dataset / pi_field), writes vmax/pmin/ifl/t0/otl to
    output_nc, and returns a JSON summary — file paths cross the protocol, not
    arrays. The input must contain variables sst (degC), msl (hPa), t (degC),
    r (mixing ratio, g/kg) with vertical pressure coordinate `dim` (hPa), and
    each must carry a matching netCDF `units` attribute: mismatched or missing
    units are rejected with an error, never converted silently.

    The summary reports the IFL flag histogram; only ifl == 1 columns are
    trustworthy (flags from the environmental-CAPE call persist, so other
    flags can accompany finite-looking numbers).
    """
    try:
        import xarray as xr

        from .xarray import pi_dataset
    except ImportError:
        return _error(
            "compute_pi_grid requires xarray/h5netcdf; install with `pip install tcpypi[mcp]`"
        )

    if not os.path.exists(input_nc):
        return _error(f"input file not found: {input_nc}", parameter="input_nc")
    if os.path.abspath(input_nc) == os.path.abspath(output_nc):
        return _error("output_nc must differ from input_nc", parameter="output_nc")

    try:
        ds = xr.open_dataset(input_nc)
    except Exception as exc:
        return _error(f"could not open {input_nc}: {exc!r}", parameter="input_nc")

    try:
        missing = [v for v in ("sst", "msl", "t", "r") if v not in ds] + (
            [] if dim in ds else [dim]
        )
        if missing:
            return _error(
                f"input dataset is missing required variable(s) {missing}",
                parameter="input_nc",
            )
        for name, expected in (
            ("sst", "degC"),
            ("msl", "hPa"),
            ("t", "degC"),
            ("r", "g/kg"),
            (dim, "hPa"),
        ):
            problem = _check_units(ds[name], name, expected)
            if problem:
                return _error(problem, variable=name)

        knobs = dict(
            CKCD=CKCD,
            ascent_flag=ascent_flag,
            diss_flag=diss_flag,
            V_reduc=V_reduc,
            ptop=ptop,
            miss_handle=miss_handle,
            outflow_source=outflow_source,
        )
        try:
            out = pi_dataset(ds, dim=dim, **knobs)
        except Exception as exc:
            return _error(f"PI calculation failed: {exc!r}")

        try:
            out.to_netcdf(output_nc)
        except Exception as exc:
            return _error(f"could not write {output_nc}: {exc!r}", parameter="output_nc")

        ifl = np.asarray(out.ifl)
        ifl_hist = {str(k): int((ifl == k).sum()) for k in np.unique(ifl[np.isfinite(ifl)])}
        ok = ifl == 1
        vmax_ok = np.asarray(out.vmax)[ok]
        return {
            "status": "success",
            "output_file": output_nc,
            "output_variables": ["vmax", "pmin", "ifl", "t0", "otl"],
            "n_columns": int(ifl.size),
            "ifl_histogram": ifl_hist,
            "ifl_note": "only ifl == 1 columns are trustworthy",
            "n_trustworthy": int(ok.sum()),
            "vmax_ms_over_trustworthy": {
                "min": _num(vmax_ok.min()) if ok.any() else None,
                "mean": _num(vmax_ok.mean()) if ok.any() else None,
                "max": _num(vmax_ok.max()) if ok.any() else None,
            },
            "provenance": _provenance(dim=dim, **knobs),
        }
    finally:
        ds.close()


def decompose_pi(
    vmax_ms: float | list[float],
    sst_c: float | list[float],
    t0_k: float | list[float],
    CKCD: float = 0.9,
) -> dict:
    """Wing et al. (2015) log decomposition of potential intensity.

    Separates ln(vmax^2) additively into thermodynamic efficiency, air-sea
    disequilibrium, and the Ck/Cd constant: lnpi = lneff + lndiseq + lnCKCD.
    Inputs: vmax_ms (m/s, e.g. from compute_pi), sst_c (degrees Celsius),
    t0_k (outflow temperature, K); scalars or equal-shaped arrays. Invalid
    physical states yield null entries. Related scalar helpers (efficiency,
    disequilibrium residual) live in tcpyPI.utilities for Python users.
    """
    try:
        lnpi, lneff, lndiseq, lnckcd = _pi_log_decomposition(
            np.asarray(vmax_ms, dtype=float),
            np.asarray(sst_c, dtype=float),
            np.asarray(t0_k, dtype=float),
            CKCD=CKCD,
            sst_units="C",
        )
    except Exception as exc:
        return _error(f"pi_log_decomposition() failed: {exc!r}")
    return {
        "status": "success",
        "lnpi": _nums(lnpi),
        "lneff": _nums(lneff),
        "lndiseq": _nums(lndiseq),
        "lnCKCD": _nums(lnckcd),
        "note": "additive in log space: lnpi = lneff + lndiseq + lnCKCD",
        "provenance": _provenance(CKCD=CKCD, sst_units="C"),
    }


def ventilation_index(
    v_shear_ms: float | list[float],
    v_pot_ms: float | list[float],
    chi_m: float | list[float],
    formulation: str = "te12",
) -> dict:
    """Tang & Emanuel (2012) ventilation index: VI = v_shear * chi_m / v_pot.

    Inputs: v_shear_ms (850-200 hPa environmental wind shear, m/s), v_pot_ms
    (potential intensity, m/s, e.g. from compute_pi), chi_m (nondimensional
    midlevel entropy deficit; compute it from soundings with
    tcpyPI.entropy_deficit_te12_from_profile in Python). Scalars or
    equal-shaped arrays. Higher VI is less favorable for genesis and
    intensification.
    """
    try:
        vi = _ventilation_index(
            np.asarray(v_shear_ms, dtype=float),
            np.asarray(v_pot_ms, dtype=float),
            np.asarray(chi_m, dtype=float),
            formulation=formulation,
        )
    except Exception as exc:
        return _error(f"ventilation_index() failed: {exc!r}")
    return {
        "status": "success",
        "ventilation_index": _nums(vi),
        "provenance": _provenance(formulation=formulation),
    }


def genesis_potential_index(
    abs_vort_s1: float | list[float],
    rh_mid_pct: float | list[float] | None = None,
    v_shear_ms: float | list[float] | None = None,
    v_pot_ms: float | list[float] | None = None,
    formulation: str = "en04",
    chi: float | list[float] | None = None,
) -> dict:
    """Genesis Potential Index (Emanuel & Nolan 2004 "en04" or Emanuel 2010 "e10").

    Inputs: abs_vort_s1 (low-level absolute vorticity, 1/s; sign is ignored),
    v_shear_ms (850-200 hPa wind shear, m/s), v_pot_ms (potential intensity,
    m/s, e.g. from compute_pi). en04 additionally requires rh_mid_pct
    (midlevel relative humidity, %); e10 instead requires chi (nondimensional
    entropy deficit). Scalars or equal-shaped arrays.
    """
    if formulation == "en04" and rh_mid_pct is None:
        return _error("formulation 'en04' requires rh_mid_pct", parameter="rh_mid_pct")
    if formulation == "e10" and chi is None:
        return _error("formulation 'e10' requires chi", parameter="chi")
    try:
        gpi = _genesis_potential_index(
            np.asarray(abs_vort_s1, dtype=float),
            None if rh_mid_pct is None else np.asarray(rh_mid_pct, dtype=float),
            None if v_shear_ms is None else np.asarray(v_shear_ms, dtype=float),
            None if v_pot_ms is None else np.asarray(v_pot_ms, dtype=float),
            formulation=formulation,
            chi=None if chi is None else np.asarray(chi, dtype=float),
        )
    except Exception as exc:
        return _error(f"genesis_potential_index() failed: {exc!r}")
    return {
        "status": "success",
        "gpi": _nums(gpi),
        "provenance": _provenance(formulation=formulation),
    }


_TOOLS = (
    compute_pi,
    compute_pi_grid,
    decompose_pi,
    ventilation_index,
    genesis_potential_index,
)


def create_server():
    """Build the FastMCP server with the five tcpyPI tools registered."""
    from fastmcp import FastMCP

    mcp = FastMCP(
        "tcpyPI",
        instructions=(
            "Validated tropical cyclone potential intensity (Bister & Emanuel 2002) "
            "and genesis diagnostics from the tcpyPI package. Units are fixed and "
            "never converted; only results with ifl == 1 are trustworthy."
        ),
    )
    for tool in _TOOLS:
        mcp.tool(tool)
    return mcp


def main():  # pragma: no cover - thin CLI shim around create_server().run()
    try:
        server = create_server()
    except ImportError as exc:
        raise SystemExit(
            f"tcpypi-mcp requires the MCP extra: pip install 'tcpypi[mcp]' ({exc})"
        ) from exc
    server.run()  # stdio transport


if __name__ == "__main__":  # pragma: no cover
    main()
