"""Power-law core-mass-function estimation utilities."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def valid_mass_mask(masses: Sequence[float]) -> np.ndarray:
    """Return the mask for finite, strictly positive masses."""
    values = np.asarray(masses, dtype=float)
    return np.isfinite(values) & (values > 0)


def powerlaw_alpha_mle(
    masses: Sequence[float],
    minimum_mass: float,
) -> float:
    """Estimate alpha for ``dN/dM proportional to M**(-alpha)``."""
    values = np.asarray(masses, dtype=float)
    minimum = float(minimum_mass)
    if not np.isfinite(minimum) or minimum <= 0:
        raise ValueError("minimum_mass must be finite and positive")
    tail = values[
        valid_mass_mask(values) & (values >= minimum)
    ]
    if tail.size < 2:
        raise ValueError("At least two masses are required above minimum_mass")
    denominator = np.sum(np.log(tail / minimum))
    if denominator <= 0:
        raise ValueError("Masses must extend above minimum_mass")
    return float(1.0 + tail.size / denominator)


def powerlaw_pdf(
    masses: Sequence[float],
    alpha: float,
    minimum_mass: float,
) -> np.ndarray:
    """Evaluate the normalized lower-truncated power-law density."""
    values = np.asarray(masses, dtype=float)
    exponent = float(alpha)
    minimum = float(minimum_mass)
    if exponent <= 1 or minimum <= 0:
        raise ValueError("alpha must exceed one and minimum_mass be positive")
    density = np.zeros_like(values, dtype=float)
    mask = np.isfinite(values) & (values >= minimum)
    density[mask] = (
        (exponent - 1.0)
        / minimum
        * (values[mask] / minimum) ** (-exponent)
    )
    return density


def powerlaw_cdf(
    masses: Sequence[float],
    alpha: float,
    minimum_mass: float,
) -> np.ndarray:
    """Evaluate the lower-truncated power-law cumulative distribution."""
    values = np.asarray(masses, dtype=float)
    exponent = float(alpha)
    minimum = float(minimum_mass)
    if exponent <= 1 or minimum <= 0:
        raise ValueError("alpha must exceed one and minimum_mass be positive")
    cumulative = np.zeros_like(values, dtype=float)
    above = np.isfinite(values) & (values >= minimum)
    cumulative[above] = 1.0 - (
        values[above] / minimum
    ) ** (1.0 - exponent)
    cumulative[values == np.inf] = 1.0
    cumulative[~np.isfinite(values) & (values != np.inf)] = np.nan
    return cumulative


def powerlaw_ks_distance(
    masses: Sequence[float],
    alpha: float,
    minimum_mass: float,
) -> float:
    """Return the two-sided empirical/model CDF distance."""
    values = np.sort(
        np.asarray(masses, dtype=float)[
            valid_mass_mask(masses)
            & (np.asarray(masses, dtype=float) >= minimum_mass)
        ]
    )
    if values.size < 2:
        raise ValueError("At least two tail masses are required")
    model = powerlaw_cdf(values, alpha, minimum_mass)
    empirical_upper = np.arange(1, values.size + 1) / values.size
    empirical_lower = np.arange(values.size) / values.size
    return float(
        max(
            np.max(np.abs(empirical_upper - model)),
            np.max(np.abs(empirical_lower - model)),
        )
    )


def diagnose_minimum_mass(
    masses: Sequence[float],
    *,
    minimum_count: int = 30,
    completeness_mass: float | None = None,
) -> dict[str, np.ndarray | float | int]:
    """Choose the lower cutoff that minimizes the power-law KS distance."""
    values = np.asarray(masses, dtype=float)
    values = values[valid_mass_mask(values)]
    if minimum_count < 2:
        raise ValueError("minimum_count must be at least two")
    if values.size < minimum_count:
        raise ValueError(
            f"Only {values.size} valid masses; minimum_count={minimum_count}"
        )
    candidates = np.unique(values)
    if completeness_mass is not None:
        completeness = float(completeness_mass)
        if not np.isfinite(completeness) or completeness <= 0:
            raise ValueError(
                "completeness_mass must be finite and positive"
            )
        candidates = candidates[candidates >= completeness]

    retained_minima: list[float] = []
    counts: list[int] = []
    exponents: list[float] = []
    distances: list[float] = []
    for minimum in candidates:
        count = int(np.count_nonzero(values >= minimum))
        if count < minimum_count:
            continue
        alpha = powerlaw_alpha_mle(values, float(minimum))
        distance = powerlaw_ks_distance(values, alpha, float(minimum))
        retained_minima.append(float(minimum))
        counts.append(count)
        exponents.append(alpha)
        distances.append(distance)
    if not retained_minima:
        raise ValueError("No cutoff candidate satisfies minimum_count")

    distances_array = np.asarray(distances)
    best = int(np.argmin(distances_array))
    return {
        "minimum_masses": np.asarray(retained_minima),
        "counts": np.asarray(counts, dtype=int),
        "alphas": np.asarray(exponents),
        "ks_distances": distances_array,
        "best_index": best,
        "minimum_mass": retained_minima[best],
        "count": counts[best],
        "alpha": exponents[best],
        "ks_distance": distances[best],
    }


def log_alpha_posterior(
    alpha: float,
    masses: Sequence[float],
    minimum_mass: float,
) -> float:
    """Return the log posterior under a flat prior for alpha > 1."""
    exponent = float(alpha)
    if not np.isfinite(exponent) or exponent <= 1:
        return -np.inf
    values = np.asarray(masses, dtype=float)
    tail = values[
        valid_mass_mask(values) & (values >= minimum_mass)
    ]
    if tail.size < 2:
        raise ValueError("At least two tail masses are required")
    return float(
        tail.size * np.log(exponent - 1.0)
        - exponent * np.sum(np.log(tail / minimum_mass))
    )


def sample_alpha_posterior(
    masses: Sequence[float],
    minimum_mass: float,
    *,
    steps: int = 60_000,
    burn_in: int = 10_000,
    proposal_width: float = 0.03,
    seed: int = 42,
) -> dict[str, np.ndarray | float]:
    """Sample the one-parameter power-law posterior by Metropolis steps."""
    if steps <= 1 or burn_in < 0 or burn_in >= steps:
        raise ValueError("Require steps > burn_in >= 0")
    if proposal_width <= 0:
        raise ValueError("proposal_width must be positive")
    values = np.asarray(masses, dtype=float)
    alpha = powerlaw_alpha_mle(values, minimum_mass)
    log_probability = log_alpha_posterior(alpha, values, minimum_mass)
    rng = np.random.default_rng(seed)
    chain = np.empty(steps, dtype=float)
    accepted = 0
    for index in range(steps):
        proposal = alpha + rng.normal(0.0, proposal_width)
        proposal_probability = log_alpha_posterior(
            proposal,
            values,
            minimum_mass,
        )
        if np.log(rng.random()) < proposal_probability - log_probability:
            alpha = proposal
            log_probability = proposal_probability
            accepted += 1
        chain[index] = alpha
    posterior = chain[burn_in:]
    return {
        "chain": chain,
        "posterior": posterior,
        "acceptance_fraction": accepted / steps,
        "median": float(np.median(posterior)),
        "p16": float(np.percentile(posterior, 16)),
        "p84": float(np.percentile(posterior, 84)),
    }


def bootstrap_powerlaw_alpha(
    masses: Sequence[float],
    minimum_mass: float,
    *,
    samples: int = 3_000,
    seed: int = 123,
) -> np.ndarray:
    """Bootstrap alpha while holding the selected lower cutoff fixed."""
    if samples <= 0:
        raise ValueError("samples must be positive")
    values = np.asarray(masses, dtype=float)
    tail = values[
        valid_mass_mask(values) & (values >= minimum_mass)
    ]
    if tail.size < 2:
        raise ValueError("At least two tail masses are required")
    rng = np.random.default_rng(seed)
    output = np.empty(samples, dtype=float)
    for index in range(samples):
        resampled = rng.choice(tail, size=tail.size, replace=True)
        output[index] = powerlaw_alpha_mle(resampled, minimum_mass)
    return output


def differential_mass_function(
    masses: Sequence[float],
    bin_edges: Sequence[float],
) -> dict[str, np.ndarray]:
    """Return counts and ``dN/dM`` on supplied mass-bin edges."""
    values = np.asarray(masses, dtype=float)
    values = values[valid_mass_mask(values)]
    edges = np.asarray(bin_edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("bin_edges must be a one-dimensional edge array")
    if np.any(~np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError("bin_edges must be finite and strictly increasing")
    counts, _ = np.histogram(values, bins=edges)
    widths = np.diff(edges)
    return {
        "edges": edges,
        "centers": np.sqrt(edges[:-1] * edges[1:]),
        "counts": counts,
        "widths": widths,
        "density": counts / widths,
        "density_error": np.sqrt(counts) / widths,
    }
