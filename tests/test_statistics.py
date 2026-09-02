import numpy as np

from leakage.paper import correlations_from_sky, isotropic_plus_point_map, load_ng15_positions
from leakage.response import pixel_response_matrix
from leakage.statistics import (
    anisotropy_likelihood_ratio,
    fit_isotropic_amplitude,
    radiometer_map,
    real_spherical_harmonic_design,
)


def _paper_array_problem():
    _, vectors = load_ng15_positions(include_short_baseline=True)
    response = pixel_response_matrix(vectors, nside=8)
    sigma = np.full(response.shape[0], 0.1)
    return response, sigma


def test_real_spherical_harmonics_are_orthonormal_under_pixel_quadrature():
    design, modes = real_spherical_harmonic_design(nside=16, lmax=3)
    gram = design.T @ design * (4.0 * np.pi / len(design))
    assert len(modes) == 16
    np.testing.assert_allclose(gram, np.eye(16), atol=1.5e-3, rtol=0.0)


def test_isotropic_amplitude_fit_recovers_known_scale():
    response, sigma = _paper_array_problem()
    correlations = 1.7 * correlations_from_sky(response, np.ones(response.shape[1]))
    amplitude, chi2 = fit_isotropic_amplitude(response, correlations, sigma)
    np.testing.assert_allclose(amplitude, 1.7, atol=2e-15)
    assert chi2 < 1e-26


def test_isotropic_likelihood_ratio_is_zero():
    response, sigma = _paper_array_problem()
    correlations = correlations_from_sky(response, np.ones(response.shape[1]))
    result = anisotropy_likelihood_ratio(
        response,
        correlations,
        sigma,
        n_starts=2,
        seed=11,
    )
    assert result.anisotropic_fit.success
    assert result.snr < 1e-8


def test_known_point_map_is_better_fit_than_isotropy():
    response, sigma = _paper_array_problem()
    sky = isotropic_plus_point_map(pixel=137, c1_over_c0=0.3)
    correlations = correlations_from_sky(response, sky)
    result = anisotropy_likelihood_ratio(
        response,
        correlations,
        sigma,
        n_starts=3,
        seed=12,
    )
    assert result.anisotropic_fit.success
    assert result.anisotropic_fit.chi2 < result.isotropic_chi2
    assert result.snr > 1.0


def test_radiometer_outputs_finite_pixel_estimates():
    response, sigma = _paper_array_problem()
    correlations = correlations_from_sky(response, np.ones(response.shape[1]))
    estimate, uncertainty = radiometer_map(response, correlations, sigma)
    assert estimate.shape == uncertainty.shape == (response.shape[1],)
    assert np.all(np.isfinite(estimate))
    assert np.all(uncertainty > 0.0)
