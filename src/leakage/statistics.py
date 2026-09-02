"""Correlation-space estimators for the paper reproduction.

The likelihood is the diagonal Gaussian of Eq. (20).  This module deliberately
does not manufacture the unpublished per-pair PFOS uncertainties: callers must
supply them explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass

import healpy as hp
import numpy as np
from scipy.optimize import least_squares
from scipy.special import sph_harm_y


@dataclass(frozen=True)
class SqrtSphericalFit:
    amplitude: float
    coefficients: np.ndarray
    sky_map: np.ndarray
    chi2: float
    success: bool


@dataclass(frozen=True)
class AnisotropyLikelihoodRatio:
    isotropic_amplitude: float
    isotropic_chi2: float
    anisotropic_fit: SqrtSphericalFit
    snr_squared: float

    @property
    def snr(self) -> float:
        return float(np.sqrt(max(self.snr_squared, 0.0)))


def _validate_data(
    response: np.ndarray,
    correlations: np.ndarray,
    sigma: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = np.asarray(response, dtype=float)
    rho = np.asarray(correlations, dtype=float)
    errors = np.asarray(sigma, dtype=float)
    if matrix.ndim != 2 or rho.shape != (matrix.shape[0],):
        raise ValueError("correlations must match the rows of response")
    if errors.shape != rho.shape or np.any(errors <= 0.0):
        raise ValueError("sigma must be positive and match correlations")
    return matrix, rho, errors


def fit_isotropic_amplitude(
    response: np.ndarray,
    correlations: np.ndarray,
    sigma: np.ndarray,
    amplitude_bounds: tuple[float, float] = (1.0e-2, 1.0e2),
) -> tuple[float, float]:
    """Fit the positive isotropic amplitude and return ``(amplitude, chi2)``."""
    matrix, rho, errors = _validate_data(response, correlations, sigma)
    template = matrix @ np.ones(matrix.shape[1])
    weights = errors**-2
    amplitude = np.sum(weights * template * rho) / np.sum(weights * template**2)
    amplitude = float(np.clip(amplitude, *amplitude_bounds))
    chi2 = float(np.sum(((rho - amplitude * template) / errors) ** 2))
    return amplitude, chi2


def real_spherical_harmonic_design(nside: int, lmax: int) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Return an orthonormal real-Ylm design in ``l, m=-l..l`` order."""
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=False)
    columns: list[np.ndarray] = []
    modes: list[tuple[int, int]] = []
    for ell in range(lmax + 1):
        for order in range(-ell, ell + 1):
            if order < 0:
                harmonic = (
                    np.sqrt(2.0)
                    * (-1.0) ** order
                    * sph_harm_y(ell, -order, theta, phi).imag
                )
            elif order == 0:
                harmonic = sph_harm_y(ell, 0, theta, phi).real
            else:
                harmonic = (
                    np.sqrt(2.0)
                    * (-1.0) ** order
                    * sph_harm_y(ell, order, theta, phi).real
                )
            columns.append(harmonic)
            modes.append((ell, order))
    return np.column_stack(columns), modes


def fit_sqrt_spherical_power(
    response: np.ndarray,
    correlations: np.ndarray,
    sigma: np.ndarray,
    *,
    power_lmax: int = 6,
    n_starts: int = 8,
    seed: int = 0,
    amplitude_bounds: tuple[float, float] = (1.0e-2, 1.0e2),
) -> SqrtSphericalFit:
    """Fit a non-negative square-root spherical-harmonic power map.

    Following the paper/MAPS convention, ``power_lmax=6`` uses a square-root
    field through ``l=3``.  Its b00 coefficient is fixed and a separate positive
    amplitude is fitted.  Multiple starts make optimizer instability visible and
    reproducible; one start is always the isotropic solution.
    """
    matrix, rho, errors = _validate_data(response, correlations, sigma)
    if power_lmax < 0 or power_lmax % 2:
        raise ValueError("power_lmax must be a non-negative even integer")
    if n_starts < 1:
        raise ValueError("n_starts must be positive")
    nside = hp.npix2nside(matrix.shape[1])
    design, _ = real_spherical_harmonic_design(nside, power_lmax // 2)
    fixed_b00 = design[:, 0]
    free_design = design[:, 1:]

    iso_amp, _ = fit_isotropic_amplitude(
        matrix,
        rho,
        errors,
        amplitude_bounds=amplitude_bounds,
    )
    lower_log, upper_log = np.log10(amplitude_bounds)

    def model(parameters: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        field = fixed_b00 + free_design @ parameters[1:]
        power = field**2
        sky = power / np.mean(power)
        return 10.0 ** parameters[0] * (matrix @ sky), sky

    def residual(parameters: np.ndarray) -> np.ndarray:
        prediction, _ = model(parameters)
        return (rho - prediction) / errors

    rng = np.random.default_rng(seed)
    starts = [np.r_[np.log10(iso_amp), np.zeros(free_design.shape[1])]]
    for _ in range(n_starts - 1):
        starts.append(
            np.r_[
                rng.uniform(lower_log, upper_log),
                rng.uniform(-1.0, 1.0, free_design.shape[1]),
            ]
        )

    lower = np.r_[lower_log, np.full(free_design.shape[1], -np.inf)]
    upper = np.r_[upper_log, np.full(free_design.shape[1], np.inf)]
    solutions = [
        least_squares(residual, start, bounds=(lower, upper), method="trf")
        for start in starts
    ]
    best = min(solutions, key=lambda result: np.sum(result.fun**2))
    _, sky = model(best.x)
    return SqrtSphericalFit(
        amplitude=float(10.0 ** best.x[0]),
        coefficients=best.x[1:].copy(),
        sky_map=sky,
        chi2=float(np.sum(best.fun**2)),
        success=bool(best.success),
    )


def anisotropy_likelihood_ratio(
    response: np.ndarray,
    correlations: np.ndarray,
    sigma: np.ndarray,
    **fit_options,
) -> AnisotropyLikelihoodRatio:
    """Evaluate Eq. (25): SNR^2 = 2 log(L_anis/L_iso)."""
    iso_amp, iso_chi2 = fit_isotropic_amplitude(response, correlations, sigma)
    anisotropic = fit_sqrt_spherical_power(
        response,
        correlations,
        sigma,
        **fit_options,
    )
    return AnisotropyLikelihoodRatio(
        isotropic_amplitude=iso_amp,
        isotropic_chi2=iso_chi2,
        anisotropic_fit=anisotropic,
        snr_squared=max(iso_chi2 - anisotropic.chi2, 0.0),
    )


def radiometer_map(
    response: np.ndarray,
    correlations: np.ndarray,
    sigma: np.ndarray,
    *,
    normalize_mean: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the per-pixel radiometer estimate and its null uncertainty.

    Each pixel is fitted alone, matching the pixel-upper-limit construction used
    by MAPS.  ``normalize_mean`` is explicit because the paper does not state
    whether its plotted maps received MAPS' post-fit all-sky normalization.
    """
    matrix, rho, errors = _validate_data(response, correlations, sigma)
    weights = errors**-2
    dirty = matrix.T @ (weights * rho)
    fisher_diagonal = np.sum(weights[:, None] * matrix**2, axis=0)
    estimate = dirty / fisher_diagonal
    uncertainty = 1.0 / np.sqrt(fisher_diagonal)
    if normalize_mean:
        normalization = 1.0 / np.mean(estimate)
        estimate = estimate * normalization
        uncertainty = uncertainty * abs(normalization)
    return estimate, uncertainty
