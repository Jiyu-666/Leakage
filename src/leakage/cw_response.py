"""Coherent Earth-term CW frequency-response helpers.

The functions in this module intentionally contain no stochastic-background
or radiometer assumptions.  They implement the four-amplitude, coherent
Earth-term frequency response. They are not the MeerKAT mapping pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

DAY_S = 86400.0


@dataclass(frozen=True)
class ProjectionDiagnostics:
    names: tuple[tuple[str, ...], ...]
    ranks: np.ndarray
    idempotence_error: float
    orthogonality_error: float


def source_basis(theta: float, phi: float):
    """Return the (m, n, Omega) convention used by pta_replicator.add_cgw."""
    ct, st = np.cos(theta), np.sin(theta)
    cp, sp = np.cos(phi), np.sin(phi)
    m = np.array([sp, -cp, 0.0])
    n = np.array([-ct * cp, -ct * sp, st])
    omega = np.array([-st * cp, -st * sp, -ct])
    return m, n, omega


def pulsar_vectors(psrs):
    ra = np.array([p.model.RAJ.quantity.to_value("rad") for p in psrs])
    dec = np.array([p.model.DECJ.quantity.to_value("rad") for p in psrs])
    return np.column_stack(
        (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
    )


def antenna_patterns(psrs, theta: float, phi: float):
    """Earth-term plus/cross antenna patterns for every pulsar."""
    m, n, omega = source_basis(theta, phi)
    p = pulsar_vectors(psrs)
    mp, np_, op = p @ m, p @ n, p @ omega
    denom = 1.0 + op
    safe = np.abs(denom) > 1e-12
    fp = np.zeros(len(psrs))
    fc = np.zeros(len(psrs))
    fp[safe] = 0.5 * (mp[safe] ** 2 - np_[safe] ** 2) / denom[safe]
    fc[safe] = mp[safe] * np_[safe] / denom[safe]
    return fp, fc


def sample_times(psrs, tref_mjd: float):
    return [
        np.asarray((p.toas.get_mjds().value - tref_mjd) * DAY_S, dtype=float)
        for p in psrs
    ]


def four_amplitude_blocks(psrs, times, frequency: float, theta: float, phi: float):
    """Four columns spanning every fixed-polarization monochromatic Earth term."""
    fp, fc = antenna_patterns(psrs, theta, phi)
    blocks = []
    for fplus, fcross, t in zip(fp, fc, times):
        phase = 2.0 * np.pi * frequency * t
        s, c = np.sin(phase), np.cos(phase)
        blocks.append(
            np.column_stack((-fplus * s, -fplus * c, -fcross * s, -fcross * c))
        )
    return blocks


def timing_projectors(psrs, sigma_s: float):
    """Weighted linear timing-model projectors from the active PINT design columns."""
    projectors, names, ranks = [], [], []
    max_idempotence = 0.0
    max_orthogonality = 0.0
    for psr in psrs:
        design, params, _ = psr.model.designmatrix(psr.toas)
        design = np.asarray(design, dtype=float)
        weighted = design / sigma_s
        # QR without rank truncation projects out extra directions for dependent
        # design columns. SVD retains only the actual timing-model column space.
        left, singular, _ = np.linalg.svd(weighted, full_matrices=False)
        rank = np.linalg.matrix_rank(weighted)
        q = left[:, :rank]
        projector = np.eye(len(design)) - q @ q.T
        projectors.append(projector)
        names.append(tuple(params))
        ranks.append(np.linalg.matrix_rank(weighted))
        max_idempotence = max(
            max_idempotence,
            np.linalg.norm(projector @ projector - projector, ord=2),
        )
        max_orthogonality = max(
            max_orthogonality,
            np.linalg.norm(q.T @ projector, ord=2),
        )
    return projectors, ProjectionDiagnostics(
        names=tuple(names),
        ranks=np.asarray(ranks),
        idempotence_error=max_idempotence,
        orthogonality_error=max_orthogonality,
    )


def apply_blocks(blocks, operators):
    return [operator @ block for operator, block in zip(operators, blocks)]


def coherent_statistic(data_blocks, template_blocks, sigma_s: float):
    """Return 2F=Delta chi2, normal-matrix condition number, and amplitudes."""
    normal = np.zeros((4, 4))
    rhs = np.zeros(4)
    inv_var = sigma_s**-2
    for data, template in zip(data_blocks, template_blocks):
        normal += inv_var * template.T @ template
        rhs += inv_var * template.T @ data
    amplitudes = np.linalg.pinv(normal, rcond=1e-12) @ rhs
    two_f = float(rhs @ amplitudes)
    condition = float(np.linalg.cond(normal))
    return two_f, condition, amplitudes


def coherent_profile(
    psrs,
    times,
    data_blocks,
    frequencies,
    theta,
    phi,
    sigma_s,
    operators=None,
):
    values = np.empty(len(frequencies))
    conditions = np.empty(len(frequencies))
    for i, frequency in enumerate(frequencies):
        templates = four_amplitude_blocks(psrs, times, frequency, theta, phi)
        if operators is not None:
            templates = apply_blocks(templates, operators)
        values[i], conditions[i], _ = coherent_statistic(
            data_blocks, templates, sigma_s
        )
    return values, conditions






def exact_window_projection(signal, times, frequencies, injection_frequency):
    """Direct DFT and exact positive+negative-frequency window prediction."""
    basis = np.column_stack(
        (
            np.cos(2 * np.pi * injection_frequency * times),
            np.sin(2 * np.pi * injection_frequency * times),
        )
    )
    coeffs, *_ = np.linalg.lstsq(basis, signal, rcond=None)
    reconstructed = basis @ coeffs
    mono_error = np.linalg.norm(signal - reconstructed) / np.linalg.norm(signal)
    complex_amplitude = coeffs[0] - 1j * coeffs[1]

    def window(delta):
        return np.exp(-2j * np.pi * np.outer(delta, times)).sum(axis=1)

    predicted = (
        0.5 * complex_amplitude * window(frequencies - injection_frequency)
        + 0.5 * np.conjugate(complex_amplitude)
        * window(frequencies + injection_frequency)
    )
    direct = np.exp(-2j * np.pi * np.outer(frequencies, times)) @ signal
    scale = max(np.max(np.abs(direct) ** 2), np.finfo(float).tiny)
    window_error = np.max(np.abs(np.abs(direct) ** 2 - np.abs(predicted) ** 2)) / scale
    return np.abs(direct) ** 2, np.abs(predicted) ** 2, mono_error, window_error


def fejer_power(n_toa: int, cadence_s: float, delta_frequency):
    """Unnormalised finite-N Fejer power with a stable zero-offset limit."""
    x = np.pi * np.asarray(delta_frequency) * cadence_s
    denominator = np.sin(x)
    ratio = np.empty_like(x, dtype=float)
    near_zero = np.abs(denominator) < 1e-14
    ratio[near_zero] = n_toa
    ratio[~near_zero] = np.sin(n_toa * x[~near_zero]) / denominator[~near_zero]
    return ratio**2


def angular_separation(theta_a, phi_a, theta_b, phi_b):
    va = np.array(
        [np.sin(theta_a) * np.cos(phi_a), np.sin(theta_a) * np.sin(phi_a), np.cos(theta_a)]
    )
    vb = np.array(
        [np.sin(theta_b) * np.cos(phi_b), np.sin(theta_b) * np.sin(phi_b), np.cos(theta_b)]
    )
    return np.arccos(np.clip(va @ vb, -1.0, 1.0))


def batch_coherent_profiles(
    psrs,
    times,
    data_sets,
    frequencies,
    theta,
    phi,
    sigma_s,
    operators=None,
):
    """Evaluate many injections against one precomputed trial-frequency bank."""
    n_sample = sum(len(t) for t in times)
    bank = np.empty((len(frequencies), n_sample, 4))
    inverse_normal = np.empty((len(frequencies), 4, 4))
    conditions = np.empty(len(frequencies))
    inv_var = sigma_s**-2
    for i, frequency in enumerate(frequencies):
        blocks = four_amplitude_blocks(psrs, times, frequency, theta, phi)
        if operators is not None:
            blocks = apply_blocks(blocks, operators)
        template = np.concatenate(blocks)
        normal = inv_var * template.T @ template
        bank[i] = template
        inverse_normal[i] = np.linalg.pinv(normal, rcond=1e-12)
        conditions[i] = np.linalg.cond(normal)

    data_matrix = np.column_stack([np.concatenate(blocks) for blocks in data_sets])
    rhs = inv_var * bank.transpose(0, 2, 1).reshape(-1, n_sample) @ data_matrix
    rhs = rhs.reshape(len(frequencies), 4, len(data_sets))
    values = np.einsum("fai,fab,fbi->fi", rhs, inverse_normal, rhs)
    return values.T, conditions
