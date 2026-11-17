"""Regression tests for the JMIN definition in ``cape()`` (GitHub issue #89).

JMIN must be the lowest profile level at or above the parcel level (Luke Davis
1/7/2017 fix, pcmin.m since 7/26/2017). Profile levels beneath the parcel
(P > PP) must not influence CAPE, the LNB, or the outflow temperature.
"""

import numpy as np
import pytest

from tcpyPI import utilities
from tcpyPI.pi import cape

from test_outflow_source import _example_profile


def _profile_si():
    """Docstring example profile in cape() units: K and g/g."""
    SSTC, MSL, P, TC, R = _example_profile()
    return SSTC, MSL, P.astype(float), utilities.T_Ctok(TC), R * 0.001


def _parcel(T, R, PP):
    # environmental parcel as lifted in the PI iteration: lowest-level temperature
    # and mixing ratio, displaced to the (sub-sounding) pressure PP
    return T[0], R[0], PP


@pytest.mark.parametrize("PP", [990.0, 962.5, 950.0, 930.0, 905.0])
def test_levels_below_parcel_do_not_affect_cape(PP):
    _, _, P, T, R = _profile_si()
    TP, RP, PP = _parcel(T, R, PP)

    full = cape(TP, RP, PP, T, R, P)

    # same sounding with every level beneath the parcel removed
    above = P <= PP
    trimmed = cape(TP, RP, PP, T[above], R[above], P[above])
    np.testing.assert_array_equal(full, trimmed)

    # same sounding with the sub-parcel levels perturbed: result must not change
    T_pert, R_pert = T.copy(), R.copy()
    T_pert[~above] += 5.0
    R_pert[~above] *= 0.5
    perturbed = cape(TP, RP, PP, T_pert, R_pert, P)
    np.testing.assert_array_equal(full, perturbed)


def test_parcel_at_lowest_level_unchanged():
    # PP == P[0] (e.g. environmental CAPE) gives JMIN = 0: must reproduce the
    # pre-fix (v1.4) values exactly
    _, _, P, T, R = _profile_si()
    out = cape(T[0], R[0], P[0], T, R, P)
    assert out[3] == 1
    np.testing.assert_allclose(
        out[:3], (54.22338721577129, 277.94160236005564, 596.5267060009274), rtol=1e-12
    )


def _saturated_parcel_column():
    # Column from data/sample_data.nc (MERRA-2 2004, March, 55S 160W), rounded:
    # the saturated parcel is buoyant only at the first level above it when PP
    # sits just below 1000 hPa.
    SSTC = 7.54
    level_data = np.array(
        [
            [1000.0, 6.92, 5.5282], [975.0, 6.3, 5.6225], [950.0, 5.1, 5.5152],
            [925.0, 3.78, 5.2711], [900.0, 2.55, 4.9444], [875.0, 1.83, 4.5559],
            [850.0, 1.05, 4.1817], [825.0, 0.2, 3.8001], [800.0, -0.56, 3.4001],
            [775.0, -1.4, 3.0324], [750.0, -2.4, 2.7354], [725.0, -3.63, 2.4676],
            [700.0, -5.04, 2.1945], [650.0, -8.02, 1.6923], [600.0, -11.62, 1.2356],
            [550.0, -15.93, 0.9181], [500.0, -20.96, 0.6881], [450.0, -26.57, 0.487],
            [400.0, -32.8, 0.3162], [350.0, -39.69, 0.1848], [300.0, -46.55, 0.0938],
            [250.0, -51.33, 0.0409], [200.0, -51.83, 0.0172], [150.0, -50.88, 0.0055],
            [100.0, -53.09, 0.003], [70.0, -53.18, 0.0028], [50.0, -52.91, 0.003],
            [40.0, -53.0, 0.0031], [30.0, -52.62, 0.0032], [20.0, -50.93, 0.0033],
            [10.0, -46.35, 0.0036],
        ]
    )
    P = level_data[:, 0]
    return SSTC, P, utilities.T_Ctok(level_data[:, 1]), level_data[:, 2] * 0.001


def test_zero_cape_unless_buoyant_strictly_above_jmin():
    # pcmin.m (2017) returns zero CAPE when no level strictly above JMIN is
    # buoyant, even if the parcel is buoyant at JMIN itself. tcpyPI v1.x keeps
    # this for parity; it makes CAPE discontinuous as PP crosses a level.
    SSTC, P, T, R = _saturated_parcel_column()
    TP = utilities.T_Ctok(SSTC)
    ES0 = utilities.es_cc(SSTC)

    def cape_star(PP):
        return cape(TP, utilities.rv(ES0, PP), PP, T, R, P)

    assert cape_star(1000.0)[0] > 0  # JMIN = 0: buoyant layer above P[0]
    assert cape_star(990.0)[0] > 0
    CAPED, TOB, LNB, IFLAG = cape_star(999.9)  # JMIN = 1, buoyant only at JMIN
    assert (CAPED, LNB, IFLAG) == (0, 0, 1)


def test_no_level_above_parcel_returns_zero_cape():
    # parcel above every retained level (all levels have P > PP): no CAPE, no LNB
    _, _, P, T, R = _profile_si()
    CAPED, TOB, LNB, IFLAG = cape(T[0], R[0], 45.0, T, R, P, ptop=50)
    assert CAPED == 0
    assert LNB == 0
    assert TOB == T[0]
    assert IFLAG == 1
