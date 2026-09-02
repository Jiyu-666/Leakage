"""Earth-term angular response used by the paper-reproduction pipeline.

Conventions
-----------
``omega`` points along gravitational-wave propagation; the physical source is
at ``-omega``.  The antenna denominator is therefore ``1 - omega . p``, as in
Eq. (4) of the target paper.  HEALPix maps use RING ordering.

The pixel response is normalized so that a map filled with one has unit sky
mean and reproduces the distinct-pulsar Hellings--Downs curve:

    Gamma_ab = 3/(8*pi) integral dOmega P(Omega)
               sum_A F_a^A(Omega) F_b^A(Omega).
"""

from __future__ import annotations

import healpy as hp
import numpy as np

HEALPIX_RING = True
HD_NORMALIZATION = 3.0 / (8.0 * np.pi)


def radec_to_unit(ra_rad, dec_rad) -> np.ndarray:
    """Convert right ascension and declination in radians to unit vectors."""
    ra = np.asarray(ra_rad, dtype=float)
    dec = np.asarray(dec_rad, dtype=float)
    return np.column_stack(
        (
            np.cos(dec) * np.cos(ra),
            np.cos(dec) * np.sin(ra),
            np.sin(dec),
        )
    )


def pair_indices(n_pulsar: int) -> tuple[np.ndarray, np.ndarray]:
    """Return distinct pulsar pairs in NumPy upper-triangle order."""
    if n_pulsar < 2:
        raise ValueError("at least two pulsars are required")
    return np.triu_indices(n_pulsar, k=1)


def _healpix_geometry(nside: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pixels = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, pixels, nest=not HEALPIX_RING)
    omega = np.column_stack(
        (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        )
    )
    # Right-handed transverse basis with m x n = omega.
    m = np.column_stack((np.sin(phi), -np.cos(phi), np.zeros_like(phi)))
    n = np.column_stack(
        (
            -np.cos(theta) * np.cos(phi),
            -np.cos(theta) * np.sin(phi),
            np.sin(theta),
        )
    )
    return omega, m, n


def antenna_patterns(pulsar_vectors: np.ndarray, nside: int) -> tuple[np.ndarray, np.ndarray]:
    """Return Earth-term plus/cross patterns with shape ``(npsr, npix)``."""
    pulsars = np.asarray(pulsar_vectors, dtype=float)
    if pulsars.ndim != 2 or pulsars.shape[1] != 3:
        raise ValueError("pulsar_vectors must have shape (npsr, 3)")
    norms = np.linalg.norm(pulsars, axis=1)
    if not np.allclose(norms, 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("pulsar vectors must be unit normalized")

    omega, m, n = _healpix_geometry(nside)
    pm = pulsars @ m.T
    pn = pulsars @ n.T
    denominator = 1.0 - pulsars @ omega.T

    numerator_plus = 0.5 * (pm**2 - pn**2)
    numerator_cross = pm * pn
    fplus = np.divide(
        numerator_plus,
        denominator,
        out=np.zeros_like(numerator_plus),
        where=np.abs(denominator) > 1e-14,
    )
    fcross = np.divide(
        numerator_cross,
        denominator,
        out=np.zeros_like(numerator_cross),
        where=np.abs(denominator) > 1e-14,
    )
    return fplus, fcross


def pixel_response_matrix(pulsar_vectors: np.ndarray, nside: int) -> np.ndarray:
    """Map relative sky power pixels to normalized distinct-pair ORFs.

    A vector of ones represents an isotropic sky with unit mean.  The returned
    matrix has shape ``(npsr * (npsr - 1) / 2, 12 * nside**2)``.
    """
    fplus, fcross = antenna_patterns(pulsar_vectors, nside)
    left, right = pair_indices(len(fplus))
    pixel_area = hp.nside2pixarea(nside)
    return (
        HD_NORMALIZATION
        * pixel_area
        * (fplus[left] * fplus[right] + fcross[left] * fcross[right])
    )


def hellings_downs(pulsar_vectors: np.ndarray) -> np.ndarray:
    """Return the normalized HD curve for all distinct pulsar pairs."""
    pulsars = np.asarray(pulsar_vectors, dtype=float)
    left, right = pair_indices(len(pulsars))
    cos_zeta = np.clip(np.einsum("ij,ij->i", pulsars[left], pulsars[right]), -1.0, 1.0)
    x = 0.5 * (1.0 - cos_zeta)
    result = np.full_like(x, 0.5)
    nonzero = x > 0.0
    result[nonzero] += 1.5 * x[nonzero] * np.log(x[nonzero]) - 0.25 * x[nonzero]
    return result

