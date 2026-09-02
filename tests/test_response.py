import healpy as hp
import numpy as np

from leakage.response import (
    HD_NORMALIZATION,
    antenna_patterns,
    hellings_downs,
    pair_indices,
    pixel_response_matrix,
    radec_to_unit,
)


PULSARS = radec_to_unit(
    np.deg2rad([5.0, 43.0, 119.0, 181.0, 251.0, 319.0]),
    np.deg2rad([-61.0, -22.0, 11.0, 47.0, -37.0, 68.0]),
)


def test_radec_vectors_are_unit_normalized():
    np.testing.assert_allclose(np.linalg.norm(PULSARS, axis=1), 1.0, atol=1e-15)


def test_pair_order_and_response_shape():
    left, right = pair_indices(len(PULSARS))
    response = pixel_response_matrix(PULSARS, nside=8)
    assert len(left) == len(right) == 15
    assert response.shape == (15, hp.nside2npix(8))


def test_antenna_patterns_are_finite():
    fplus, fcross = antenna_patterns(PULSARS, nside=8)
    assert np.all(np.isfinite(fplus))
    assert np.all(np.isfinite(fcross))


def test_isotropic_pixel_integral_recovers_hellings_downs():
    response = pixel_response_matrix(PULSARS, nside=32)
    recovered = response @ np.ones(response.shape[1])
    expected = hellings_downs(PULSARS)
    np.testing.assert_allclose(recovered, expected, atol=8e-4, rtol=0.0)


def test_single_pixel_matches_direct_pair_kernel():
    nside = 8
    pixel = 137
    fplus, fcross = antenna_patterns(PULSARS, nside)
    left, right = pair_indices(len(PULSARS))
    direct = (
        HD_NORMALIZATION
        * hp.nside2pixarea(nside)
        * (fplus[left, pixel] * fplus[right, pixel]
           + fcross[left, pixel] * fcross[right, pixel])
    )
    response = pixel_response_matrix(PULSARS, nside)
    np.testing.assert_allclose(response[:, pixel], direct, atol=1e-15, rtol=0.0)


def test_nside8_is_sufficient_for_percent_level_hd_quadrature():
    response = pixel_response_matrix(PULSARS, nside=8)
    recovered = response @ np.ones(response.shape[1])
    expected = hellings_downs(PULSARS)
    np.testing.assert_allclose(recovered, expected, atol=1.2e-2, rtol=0.0)
