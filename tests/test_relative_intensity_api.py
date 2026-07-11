import numpy as np

import tcpyPI


def test_relative_intensity_scalar():
    out = tcpyPI.relative_intensity(20.0, 40.0)
    np.testing.assert_allclose(out, 0.5, rtol=0, atol=0)


def test_relative_intensity_masks_invalid_pi():
    out = tcpyPI.relative_intensity(np.array([10.0, 10.0]), np.array([0.0, -1.0]))
    assert np.all(np.isnan(out))


def test_relative_intensity_clip():
    out = tcpyPI.relative_intensity(
        np.array([10.0, 200.0]), np.array([20.0, 20.0]), clip=(0.0, 2.0)
    )
    np.testing.assert_allclose(out, np.array([0.5, 2.0]), rtol=0, atol=0)
