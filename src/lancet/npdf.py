"""Column-density probability-distribution utilities."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np


def positive_column_density_excess(
    column_density: np.ndarray,
    background: float,
) -> np.ndarray:
    """Return finite positive values after subtracting a fixed background."""
    values = np.asarray(column_density, dtype=float) - float(background)
    return values[np.isfinite(values) & (values > 0)]


def shared_excess_mean(
    regions: Mapping[str, np.ndarray],
    background: float,
) -> float:
    """Return the mean positive background-subtracted value across regions."""
    samples = [
        positive_column_density_excess(values, background)
        for values in regions.values()
    ]
    if not samples or any(values.size == 0 for values in samples):
        raise ValueError(
            "Every region must contain positive background-subtracted pixels"
        )
    mean = float(np.mean(np.concatenate(samples)))
    if not np.isfinite(mean) or mean <= 0:
        raise ValueError("The shared column-density mean must be positive")
    return mean


def logarithmic_column_density_contrast(
    column_density: np.ndarray,
    *,
    background: float,
    normalization: float,
) -> np.ndarray:
    """Return eta = ln((N-background)/normalization) for valid pixels."""
    scale = float(normalization)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("normalization must be finite and positive")
    excess = positive_column_density_excess(column_density, background)
    return np.log(excess / scale)


def npdf_histogram(
    eta: Sequence[float],
    bin_edges: Sequence[float],
) -> dict[str, np.ndarray]:
    """Return a normalized eta histogram with Poisson counting errors."""
    values = np.asarray(eta, dtype=float)
    values = values[np.isfinite(values)]
    edges = np.asarray(bin_edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError("bin_edges must be a one-dimensional edge array")
    widths = np.diff(edges)
    if np.any(~np.isfinite(edges)) or np.any(widths <= 0):
        raise ValueError("bin_edges must be finite and strictly increasing")
    counts, _ = np.histogram(values, bins=edges)
    total = int(np.sum(counts))
    if total == 0:
        raise ValueError("No eta values fall inside the configured bins")
    density = counts / (total * widths)
    error = np.sqrt(counts) / (total * widths)
    return {
        "edges": edges,
        "centers": 0.5 * (edges[:-1] + edges[1:]),
        "widths": widths,
        "counts": counts,
        "density": density,
        "density_error": error,
    }
