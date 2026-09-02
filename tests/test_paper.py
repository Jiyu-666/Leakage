import healpy as hp
import numpy as np
import pytest

from leakage.paper import (
    PAPER_ANISOTROPY_LEVELS,
    PAPER_PIXEL_DETECTION_RATE,
    PAPER_REPORTED_FREQUENCIES_NHZ,
    PAPER_SNR_DETECTION_RATE,
    angular_power_ratio,
    candidate_uniform_pixels,
    correlations_from_sky,
    fourier_frequencies_nhz,
    isotropic_plus_point_map,
    load_ng15_positions,
)
from leakage.response import hellings_downs, pixel_response_matrix


def test_public_analysis_array_has_67_pulsars_and_2211_pairs():
    names, vectors = load_ng15_positions()
    assert len(names) == 67
    assert len(names) * (len(names) - 1) // 2 == 2211
    assert len(np.unique(names)) == 67
    np.testing.assert_allclose(np.linalg.norm(vectors, axis=1), 1.0, atol=2e-15)


def test_literal_paper_branch_has_68_pulsars_and_2278_pairs():
    names, vectors = load_ng15_positions(include_short_baseline=True)
    assert len(names) == 68
    assert len(names) * (len(names) - 1) // 2 == 2278
    assert names[-1] == "J0614-3329"
    np.testing.assert_allclose(np.linalg.norm(vectors, axis=1), 1.0, atol=2e-15)


def test_paper_frequency_rounding_discrepancy_is_preserved():
    exact_grid = fourier_frequencies_nhz()
    # The listed values are multiples of the rounded 1.98 nHz, whereas k/T is
    # not.  Preserve both branches instead of silently weakening the check.
    np.testing.assert_allclose(
        PAPER_REPORTED_FREQUENCIES_NHZ,
        np.arange(1, 6) * 1.98,
        atol=1e-14,
        rtol=0.0,
    )
    assert round(exact_grid[0], 2) == 1.98
    assert 0.01 < np.max(np.abs(exact_grid - PAPER_REPORTED_FREQUENCIES_NHZ)) < 0.02


def test_reported_rates_are_multiples_of_one_over_48_to_rounding():
    rates = np.concatenate((PAPER_SNR_DETECTION_RATE.ravel(), PAPER_PIXEL_DETECTION_RATE.ravel()))
    counts = np.rint(48 * rates)
    np.testing.assert_allclose(rates, counts / 48.0, atol=5.1e-5, rtol=0.0)


def test_candidate_direction_set_has_48_distinct_nside8_pixels():
    pixels = candidate_uniform_pixels()
    assert len(pixels) == 48
    assert len(np.unique(pixels)) == 48
    assert np.all((pixels >= 0) & (pixels < hp.nside2npix(8)))


@pytest.mark.parametrize("ratio", PAPER_ANISOTROPY_LEVELS)
def test_point_injection_has_requested_c1_over_c0(ratio):
    sky = isotropic_plus_point_map(pixel=137, c1_over_c0=ratio)
    assert np.all(sky >= 0.0)
    np.testing.assert_allclose(np.mean(sky), 1.0, atol=2e-15)
    np.testing.assert_allclose(angular_power_ratio(sky), ratio, atol=2e-11)


def test_theoretical_isotropic_correlations_close_hd_at_nside8():
    _, vectors = load_ng15_positions()
    response = pixel_response_matrix(vectors, nside=8)
    rho = correlations_from_sky(response, np.ones(hp.nside2npix(8)))
    np.testing.assert_allclose(rho, hellings_downs(vectors), atol=1.2e-2, rtol=0.0)
