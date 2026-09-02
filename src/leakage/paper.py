"""Locked numerical inputs and targets from the paper's point-source study.

Only values stated in the paper are labelled ``PAPER_*``.  The 48-direction
construction is explicitly labelled as a candidate because the paper does not
publish its pixel indices.
"""

from __future__ import annotations

from pathlib import Path

import healpy as hp
import numpy as np
from scipy.optimize import brentq

PAPER_NSIDE = 8
PAPER_LMAX = 6
PAPER_OBSERVATION_SPAN_YR = 16.03
PAPER_REPORTED_FREQUENCIES_NHZ = np.array([1.98, 3.96, 5.94, 7.92, 9.90])
PAPER_ANISOTROPY_LEVELS = np.array([0.3, 0.4, 0.5])
PAPER_MEDIAN_PAIR_SIGMA = np.array([11.76, 11.95, 15.61, 20.83, 30.78])

# Rows are C1/C0 = 0.3, 0.4, 0.5; columns are the first five frequency bins.
PAPER_SNR_DETECTION_RATE = np.array(
    [
        [0.0625, 0.0417, 0.0, 0.0, 0.0],
        [0.1875, 0.0417, 0.0, 0.0, 0.0],
        [0.3125, 0.1875, 0.0208, 0.0, 0.0],
    ]
)
PAPER_PIXEL_DETECTION_RATE = np.array(
    [
        [0.4583, 0.4583, 0.2500, 0.2292, 0.2083],
        [0.4583, 0.5208, 0.3542, 0.3125, 0.1667],
        [0.5833, 0.3958, 0.3125, 0.2708, 0.2500],
    ]
)

PULSAR_TABLE = Path(__file__).parents[2] / "data" / "nanograv15_67_positions.csv"
EXCLUDED_PULSAR_TABLE = (
    Path(__file__).parents[2] / "data" / "nanograv15_excluded_position.csv"
)


def fourier_frequencies_nhz(
    span_yr: float = PAPER_OBSERVATION_SPAN_YR,
    nbin: int = 5,
) -> np.ndarray:
    """Return the positive Fourier frequencies k/T in nHz."""
    seconds_per_year = 365.25 * 24.0 * 3600.0
    return np.arange(1, nbin + 1) / (span_yr * seconds_per_year) * 1.0e9


def load_ng15_positions(
    path: str | Path = PULSAR_TABLE,
    *,
    include_short_baseline: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Load the public NG15 sky positions as names and unit vectors.

    The default is the standard 67-pulsar GWB array.  Set
    ``include_short_baseline=True`` for the literal 68-pulsar statement in the
    target paper; this appends J0614-3329, which standard NG15 analyses omit.
    """
    table = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    names = np.asarray(table["name"])
    theta = np.asarray(table["theta_rad"], dtype=float)
    phi = np.asarray(table["phi_rad"], dtype=float)
    if include_short_baseline:
        excluded = np.genfromtxt(
            EXCLUDED_PULSAR_TABLE,
            delimiter=",",
            names=True,
            dtype=None,
            encoding="utf-8",
        )
        names = np.append(names, excluded["name"].item())
        theta = np.append(theta, float(excluded["theta_rad"]))
        phi = np.append(phi, float(excluded["phi_rad"]))
    vectors = np.column_stack(
        (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        )
    )
    return names, vectors


def candidate_uniform_pixels(nside: int = PAPER_NSIDE) -> np.ndarray:
    """Candidate reconstruction of the paper's unpublished 48 directions.

    The centers of all 48 Nside=2 RING pixels are projected onto the target
    Nside grid.  This is a reproducible sensitivity-analysis branch, not a
    claim about the authors' exact pixel list.
    """
    theta, phi = hp.pix2ang(2, np.arange(hp.nside2npix(2)), nest=False)
    return hp.ang2pix(nside, theta, phi, nest=False)


def angular_power_ratio(sky_map: np.ndarray) -> float:
    """Compute C1/C0 using the paper's C_l=(2l+1)^-1 sum_m |c_lm|^2."""
    values = np.asarray(sky_map, dtype=float)
    hp.npix2nside(len(values))
    cl = hp.anafast(values, lmax=1, iter=3)
    return float(cl[1] / cl[0])


def isotropic_plus_point_map(
    pixel: int,
    c1_over_c0: float,
    nside: int = PAPER_NSIDE,
) -> np.ndarray:
    """Make a unit-mean isotropic sky plus one non-negative point pixel.

    The point amplitude is solved on the actual HEALPix grid, then the whole
    map is rescaled to unit mean.  The rescaling preserves C1/C0 and matches
    the amplitude normalization removed by an optimal-statistic estimate.
    """
    npix = hp.nside2npix(nside)
    if pixel < 0 or pixel >= npix:
        raise ValueError("pixel is outside the HEALPix map")
    if not 0.0 < c1_over_c0 < 1.0:
        raise ValueError("c1_over_c0 must lie strictly between zero and one")

    base = np.ones(npix)
    point = np.zeros(npix)
    point[pixel] = 1.0

    def mismatch(amplitude: float) -> float:
        return angular_power_ratio(base + amplitude * point) - c1_over_c0

    upper = float(npix)
    while mismatch(upper) < 0.0:
        upper *= 2.0
    amplitude = brentq(mismatch, 0.0, upper, xtol=1e-12, rtol=1e-12)
    result = base + amplitude * point
    return result / np.mean(result)


def correlations_from_sky(response: np.ndarray, sky_map: np.ndarray) -> np.ndarray:
    """Generate the paper's unsampled theoretical pair correlations."""
    matrix = np.asarray(response, dtype=float)
    values = np.asarray(sky_map, dtype=float)
    if matrix.ndim != 2 or values.ndim != 1 or matrix.shape[1] != len(values):
        raise ValueError("response and sky_map shapes are incompatible")
    return matrix @ values
