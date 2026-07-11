"""Tests for the whole-field fast path pi_field() (see GitHub issue #86).

pi_field wraps the same compiled per-column core as pi() in a parallel
generalized ufunc, so for identical float64 inputs the outputs must be
EXACTLY equal (bit-for-bit) — that is the acceptance criterion here.
"""

import numpy as np

import tcpyPI

# Compact ocean sounding (from the pi() docstring example, thinned).
P = np.array([1000.0, 975, 950, 925, 900, 850, 800, 750, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50])
TC = np.array([28.0, 25, 24, 23, 22, 19, 16, 13, 11, 5, -2, -11, -27, -37, -49, -65, -79, -73, -64])
R = np.array([18.0, 18, 16, 13, 12, 10, 9, 7, 4, 1.7, 1.7, 0.2, 0.1, 0.05, 0.014, 0.003, 0.002, 0.002, 0.002])


def _field_inputs(nlat=4, nlon=5):
    """Build a small synthetic field with varied SSTs and per-cell profiles."""
    rng = np.random.default_rng(42)
    sst = 24.0 + 8.0 * rng.random((nlat, nlon))
    msl = 1005.0 + 10.0 * rng.random((nlat, nlon))
    tc = TC + rng.uniform(-1.0, 1.0, (nlat, nlon, P.size))
    r = R * rng.uniform(0.8, 1.2, (nlat, nlon, P.size))
    return sst, msl, tc, r


def _loop_pi(sst, msl, tc, r, **kw):
    outs = [np.empty(sst.shape) for _ in range(5)]
    for idx in np.ndindex(sst.shape):
        vals = tcpyPI.pi(sst[idx], msl[idx], P, tc[idx], r[idx], **kw)
        for o, v in zip(outs, vals):
            o[idx] = v
    return outs


def _assert_exact(field_outs, loop_outs):
    for got, ref in zip(field_outs, loop_outs):
        np.testing.assert_array_equal(np.asarray(got, dtype=float), ref)


def test_pi_field_exact_equality_with_pi():
    sst, msl, tc, r = _field_inputs()
    field = tcpyPI.pi_field(sst, msl, P, tc, r, CKCD=0.9)
    loop = _loop_pi(sst, msl, tc, r, CKCD=0.9)
    _assert_exact(field, loop)
    assert field[2].dtype == np.int64  # IFL is an integer array


def test_pi_field_nan_cells_and_flags():
    sst, msl, tc, r = _field_inputs()
    sst[0, 0] = np.nan   # land-like cell -> IFL=0
    msl[1, 1] = np.nan   # missing MSL cell -> IFL=3
    vmax, pmin, ifl, to, otl = tcpyPI.pi_field(sst, msl, P, tc, r)
    assert ifl[0, 0] == 0 and np.isnan(vmax[0, 0])
    assert ifl[1, 1] == 3 and np.isnan(vmax[1, 1])
    # other cells unaffected — still exactly equal to the per-column path
    loop = _loop_pi(sst, msl, tc, r)
    _assert_exact((vmax, pmin, ifl, to, otl), loop)


def test_pi_field_order_agnostic_profiles():
    sst, msl, tc, r = _field_inputs(2, 2)
    ref = tcpyPI.pi_field(sst, msl, P, tc, r)
    rev = tcpyPI.pi_field(sst, msl, P[::-1], tc[..., ::-1], r[..., ::-1])
    _assert_exact(rev, [np.asarray(x, dtype=float) for x in ref])


def test_pi_field_nan_pressure_returns_missing_field():
    sst, msl, tc, r = _field_inputs(2, 2)
    Pbad = P.copy()
    Pbad[3] = np.nan
    vmax, pmin, ifl, to, otl = tcpyPI.pi_field(sst, msl, Pbad, tc, r)
    assert np.all(ifl == 3) and np.all(np.isnan(vmax))


def test_pi_field_broadcasting_scalars_and_shared_profile():
    # scalar MSL + shared 1-D profile broadcast across an SST grid
    sst = np.array([[28.0, 30.0]])
    out = tcpyPI.pi_field(sst, 1010.0, P, TC, R)
    assert out[0].shape == (1, 2)
    ref0 = tcpyPI.pi(28.0, 1010.0, P, TC, R)
    ref1 = tcpyPI.pi(30.0, 1010.0, P, TC, R)
    assert out[0][0, 0] == ref0[0] and out[0][0, 1] == ref1[0]


def test_pi_field_input_validation():
    sst, msl, tc, r = _field_inputs(2, 2)
    with np.testing.assert_raises(ValueError):
        tcpyPI.pi_field(sst, msl, np.stack([P, P]), tc, r)  # 2-D pressure array
    with np.testing.assert_raises(ValueError):
        tcpyPI.pi_field(sst, msl, P[:-1], tc, r)  # level-dim mismatch
    with np.testing.assert_raises(ValueError):
        tcpyPI.pi_field(sst, msl, P, tc, r, outflow_source="nope")
