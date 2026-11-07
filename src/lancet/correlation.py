"""Pairwise correlations and power-law regression utilities."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from scipy import odr, optimize, stats
from .physical_properties import resolved_core_mask


def log10_positive(values: np.ndarray | Sequence[float]) -> np.ndarray:
    """Return log10 values, using NaN for non-positive or non-finite entries."""
    array = np.asarray(values, dtype=float)
    transformed = np.full(array.shape, np.nan, dtype=float)
    valid = np.isfinite(array) & (array > 0)
    transformed[valid] = np.log10(array[valid])
    return transformed


def prepare_correlation_data(
    table: Any,
    columns: Sequence[str],
    *,
    log_columns: Sequence[str] = (),
    resolved_only_columns: Sequence[str] = (),
) -> dict[str, np.ndarray]:
    """Extract configured columns and apply declared log10 transformations."""
    requested = [str(column) for column in columns]
    missing = set(requested) - set(table.colnames)
    if missing:
        raise KeyError(
            "Correlation catalogue is missing columns: "
            + ", ".join(sorted(missing))
        )
    logarithmic = {str(column) for column in log_columns}
    undeclared = logarithmic - set(requested)
    if undeclared:
        raise ValueError(
            "Log columns are not present in the configured order: "
            + ", ".join(sorted(undeclared))
        )

    data: dict[str, np.ndarray] = {}
    for column in requested:
        values = np.asarray(table[column], dtype=float)
        data[column] = (
            log10_positive(values) if column in logarithmic else values.copy()
        )
        if column in resolved_only_columns:
            data[column][~resolved_core_mask(table)] = np.nan
    return data


def pairwise_correlations(
    data: Mapping[str, np.ndarray],
    *,
    minimum_points: int = 3,
) -> dict[str, np.ndarray | list[str]]:
    """Calculate Pearson and Spearman statistics with pairwise masking."""
    if minimum_points < 2:
        raise ValueError("minimum_points must be at least 2")
    keys = list(data)
    if not keys:
        raise ValueError("At least one correlation variable is required")
    lengths = {len(np.asarray(data[key])) for key in keys}
    if len(lengths) != 1:
        raise ValueError("All correlation arrays must have the same length")

    shape = (len(keys), len(keys))
    pearson = np.full(shape, np.nan, dtype=float)
    pearson_p = np.full(shape, np.nan, dtype=float)
    spearman = np.full(shape, np.nan, dtype=float)
    spearman_p = np.full(shape, np.nan, dtype=float)
    n_points = np.zeros(shape, dtype=np.int32)

    for row, y_key in enumerate(keys):
        y = np.asarray(data[y_key], dtype=float)
        for column, x_key in enumerate(keys):
            x = np.asarray(data[x_key], dtype=float)
            valid = np.isfinite(x) & np.isfinite(y)
            x_valid = x[valid]
            y_valid = y[valid]
            n_points[row, column] = len(x_valid)
            if len(x_valid) < minimum_points:
                continue
            if np.ptp(x_valid) == 0 or np.ptp(y_valid) == 0:
                continue
            pearson_result = stats.pearsonr(x_valid, y_valid)
            spearman_result = stats.spearmanr(x_valid, y_valid)
            pearson[row, column] = pearson_result.statistic
            pearson_p[row, column] = pearson_result.pvalue
            spearman[row, column] = spearman_result.statistic
            spearman_p[row, column] = spearman_result.pvalue

    return {
        "keys": keys,
        "pearson": pearson,
        "pearson_p": pearson_p,
        "spearman": spearman,
        "spearman_p": spearman_p,
        "n_points": n_points,
    }


def prepare_powerlaw_data(
    table: Any,
    *,
    x_column: str,
    y_column: str,
    y_error_column: str,
    x_fractional_error: float,
    minimum_points: int = 5,
    resolved_only_columns: Sequence[str] = (),
) -> dict[str, np.ndarray | int]:
    """Prepare positive finite measurements for a log-log power-law fit."""
    required = {x_column, y_column, y_error_column}
    missing = required - set(table.colnames)
    if missing:
        raise KeyError(
            "Regression catalogue is missing columns: "
            + ", ".join(sorted(missing))
        )
    fraction = float(x_fractional_error)
    if not np.isfinite(fraction) or fraction < 0:
        raise ValueError("x_fractional_error must be finite and non-negative")

    x_all = np.asarray(table[x_column], dtype=float)
    y_all = np.asarray(table[y_column], dtype=float)
    y_error_all = np.asarray(table[y_error_column], dtype=float)
    x_error_all = fraction * x_all
    valid = (
        np.isfinite(x_all)
        & np.isfinite(y_all)
        & np.isfinite(x_error_all)
        & np.isfinite(y_error_all)
        & (x_all > 0)
        & (y_all > 0)
        & (x_error_all >= 0)
        & (y_error_all >= 0)
    )
    if {x_column, y_column} & set(resolved_only_columns):
        valid &= resolved_core_mask(table)
    x = x_all[valid]
    y = y_all[valid]
    x_error = x_error_all[valid]
    y_error = y_error_all[valid]
    if len(x) < minimum_points:
        raise ValueError(
            f"Need at least {minimum_points} valid positive regression points"
        )

    return {
        "x": x,
        "y": y,
        "x_error": x_error,
        "y_error": y_error,
        "log_x": np.log10(x),
        "log_y": np.log10(y),
        "log_x_error": x_error / (x * np.log(10.0)),
        "log_y_error": y_error / (y * np.log(10.0)),
        "valid_mask": valid,
        "n_total": len(table),
        "n_valid": len(x),
    }


def powerlaw_negative_log_likelihood(
    parameters: np.ndarray | Sequence[float],
    log_x: np.ndarray,
    log_y: np.ndarray,
    log_x_error: np.ndarray,
    log_y_error: np.ndarray,
) -> float:
    """Return the configured power-law negative log likelihood."""
    slope, intercept, log_intrinsic_scatter = np.asarray(
        parameters,
        dtype=float,
    )
    variance = (
        np.asarray(log_y_error, dtype=float) ** 2
        + (slope * np.asarray(log_x_error, dtype=float)) ** 2
        + np.exp(2.0 * log_intrinsic_scatter)
    )
    residual = (
        np.asarray(log_y, dtype=float)
        - slope * np.asarray(log_x, dtype=float)
        - intercept
    )
    if np.any(~np.isfinite(variance)) or np.any(variance <= 0):
        return np.inf
    return float(
        0.5
        * np.sum(
            residual**2 / variance + np.log(2.0 * np.pi * variance)
        )
    )


def fit_powerlaw_mle(
    log_x: np.ndarray,
    log_y: np.ndarray,
    log_x_error: np.ndarray,
    log_y_error: np.ndarray,
) -> optimize.OptimizeResult:
    """Fit the power-law likelihood inside the priors used for MCMC."""
    slope, intercept = np.polyfit(log_x, log_y, 1)
    scatter = max(
        float(np.std(log_y - slope * log_x - intercept)),
        np.exp(-9.0),
    )
    result = optimize.minimize(
        powerlaw_negative_log_likelihood,
        np.asarray([slope, intercept, np.log(scatter)]),
        args=(log_x, log_y, log_x_error, log_y_error),
        method="L-BFGS-B",
        bounds=((-20.0, 20.0), (-50.0, 50.0), (-10.0, 5.0)),
        options={"maxiter": 100_000, "ftol": 1e-12, "gtol": 1e-9},
    )
    return result


def fit_effective_chi2(
    log_x: np.ndarray,
    log_y: np.ndarray,
    log_x_error: np.ndarray,
    log_y_error: np.ndarray,
) -> optimize.OptimizeResult:
    """Fit an ODR-style effective chi-squared sensitivity model."""

    def objective(parameters: np.ndarray) -> float:
        slope, intercept = parameters
        variance = log_y_error**2 + (slope * log_x_error) ** 2
        if np.any(variance <= 0) or np.any(~np.isfinite(variance)):
            return np.inf
        residual = log_y - slope * log_x - intercept
        return float(np.sum(residual**2 / variance))

    return optimize.minimize(
        objective,
        np.polyfit(log_x, log_y, 1),
        method="Nelder-Mead",
        options={"maxiter": 100_000, "xatol": 1e-11, "fatol": 1e-10},
    )


def run_powerlaw_mcmc(
    data: Mapping[str, np.ndarray | int],
    *,
    n_walkers: int,
    n_steps: int,
    burn_in: int,
    thin: int,
    seed: int,
) -> dict[str, np.ndarray | float | int]:
    """Sample the configured power-law model deterministically with emcee."""
    import emcee

    if n_walkers < 8:
        raise ValueError("n_walkers must be at least 8 for a 3-D model")
    if n_steps <= 0 or burn_in < 0 or burn_in >= n_steps or thin <= 0:
        raise ValueError("Require n_steps > burn_in >= 0 and thin > 0")

    log_x = np.asarray(data["log_x"], dtype=float)
    log_y = np.asarray(data["log_y"], dtype=float)
    log_x_error = np.asarray(data["log_x_error"], dtype=float)
    log_y_error = np.asarray(data["log_y_error"], dtype=float)
    mle = fit_powerlaw_mle(
        log_x,
        log_y,
        log_x_error,
        log_y_error,
    )

    def log_probability(parameters: np.ndarray) -> float:
        slope, intercept, log_scatter = parameters
        if not (
            -20.0 < slope < 20.0
            and -50.0 < intercept < 50.0
            and -10.0 < log_scatter < 5.0
        ):
            return -np.inf
        return -powerlaw_negative_log_likelihood(
            parameters,
            log_x,
            log_y,
            log_x_error,
            log_y_error,
        )

    np.random.seed(seed)
    random = np.random.default_rng(seed)
    initial = np.array(mle.x, dtype=float, copy=True)
    initial[2] = np.clip(initial[2], -9.0, 4.0)
    positions = initial + 1.0e-4 * random.standard_normal((n_walkers, 3))
    sampler = emcee.EnsembleSampler(n_walkers, 3, log_probability)
    sampler.run_mcmc(positions, n_steps, progress=False)
    samples = sampler.get_chain(
        discard=burn_in,
        thin=thin,
        flat=True,
    )
    if len(samples) == 0:
        raise RuntimeError("No MCMC samples remain after burn-in and thinning")
    quantiles = np.percentile(samples, [2.5, 16, 50, 84, 97.5], axis=0)
    scatter_quantiles = np.percentile(
        np.exp(samples[:, 2]),
        [2.5, 16, 50, 84, 97.5],
    )
    try:
        autocorrelation_time = sampler.get_autocorr_time(quiet=True)
    except emcee.autocorr.AutocorrError:
        autocorrelation_time = np.full(3, np.nan)

    return {
        "samples": samples,
        "quantiles": quantiles,
        "scatter_quantiles": scatter_quantiles,
        "acceptance_fraction": float(
            np.mean(sampler.acceptance_fraction)
        ),
        "autocorrelation_time": autocorrelation_time,
        "burn_in": burn_in,
        "thin": thin,
        "seed": seed,
        "mle": np.asarray(mle.x, dtype=float),
    }


def compare_powerlaw_methods(
    data: Mapping[str, np.ndarray | int],
    *,
    sensitivity_fractions: Sequence[float],
) -> dict[str, Any]:
    """Compare common log-log regressions and uncertainty prescriptions."""
    log_x = np.asarray(data["log_x"], dtype=float)
    log_y = np.asarray(data["log_y"], dtype=float)
    log_x_error = np.asarray(data["log_x_error"], dtype=float)
    log_y_error = np.asarray(data["log_y_error"], dtype=float)
    pearson = stats.pearsonr(log_x, log_y)
    spearman = stats.spearmanr(log_x, log_y)
    ols = stats.linregress(log_x, log_y)
    theil = stats.theilslopes(log_y, log_x, alpha=0.95)

    model = odr.Model(
        lambda coefficients, values: (
            coefficients[0] * values + coefficients[1]
        )
    )
    odr_data = odr.RealData(
        log_x,
        log_y,
        sx=log_x_error,
        sy=log_y_error,
    )
    odr_fit = odr.ODR(
        odr_data,
        model,
        beta0=np.polyfit(log_x, log_y, 1),
        maxit=10_000,
    ).run()

    fractions = np.asarray(sensitivity_fractions, dtype=float)
    if (
        fractions.ndim != 1
        or len(fractions) == 0
        or np.any(~np.isfinite(fractions))
        or np.any(fractions < 0)
    ):
        raise ValueError(
            "sensitivity_fractions must be finite and non-negative"
        )
    likelihood_slopes = np.empty(len(fractions), dtype=float)
    effective_chi2_slopes = np.empty(len(fractions), dtype=float)
    for index, fraction in enumerate(fractions):
        # For xerr = f*x, sigma(log10 x) = f/ln(10).
        trial_log_x_error = np.full_like(
            log_x,
            fraction / np.log(10.0),
        )
        likelihood_slopes[index] = fit_powerlaw_mle(
            log_x,
            log_y,
            trial_log_x_error,
            log_y_error,
        ).x[0]
        effective_chi2_slopes[index] = fit_effective_chi2(
            log_x,
            log_y,
            trial_log_x_error,
            log_y_error,
        ).x[0]

    return {
        "pearson": pearson,
        "spearman": spearman,
        "ols": ols,
        "theil": theil,
        "odr": odr_fit,
        "sensitivity_fractions": fractions,
        "likelihood_slopes": likelihood_slopes,
        "effective_chi2_slopes": effective_chi2_slopes,
    }
