"""Tests for the core–environment correlation workflow."""

from pathlib import Path
import subprocess
import sys

from astropy.table import Table
import numpy as np
import pytest

from lancet.correlation import (
    pairwise_correlations,
    powerlaw_negative_log_likelihood,
    prepare_correlation_data,
    prepare_powerlaw_data,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = (
    PROJECT_ROOT
    / "results"
    / "environment_measurements"
    / "tables"
    / "core_environment_catalogue.ecsv"
)


def test_correlation_preparation_masks_invalid_log_values() -> None:
    table = Table(
        {
            "linear": [1.0, np.nan, 3.0, 4.0],
            "logged": [1.0, 10.0, 0.0, -1.0],
        }
    )
    data = prepare_correlation_data(
        table,
        ["linear", "logged"],
        log_columns=["logged"],
    )

    np.testing.assert_allclose(
        data["logged"][:2],
        [0.0, 1.0],
    )
    assert np.isnan(data["logged"][2:]).all()


def test_pairwise_correlations_use_pairwise_finite_masks() -> None:
    data = {
        "x": np.array([1.0, 2.0, 3.0, np.nan]),
        "y": np.array([2.0, 4.0, 6.0, 8.0]),
        "z": np.array([3.0, np.nan, 1.0, 0.0]),
    }
    result = pairwise_correlations(data)

    assert result["n_points"][1, 0] == 3
    assert result["pearson"][1, 0] == pytest.approx(1.0)
    assert result["spearman"][1, 0] == pytest.approx(1.0)
    assert result["n_points"][2, 0] == 2
    assert np.isnan(result["pearson"][2, 0])


def test_powerlaw_data_uses_fraction_over_ln10() -> None:
    table = Table(
        {
            "x": [1.0, 10.0, 100.0, 1000.0, 1.0e4],
            "y": [2.0, 20.0, 200.0, 2000.0, 2.0e4],
            "yerr": [0.2, 2.0, 20.0, 200.0, 2000.0],
        }
    )
    data = prepare_powerlaw_data(
        table,
        x_column="x",
        y_column="y",
        y_error_column="yerr",
        x_fractional_error=0.33,
    )

    np.testing.assert_allclose(
        data["log_x_error"],
        0.33 / np.log(10.0),
    )
    np.testing.assert_allclose(
        data["log_y_error"],
        0.10 / np.log(10.0),
    )


def test_powerlaw_likelihood_prefers_generating_relation() -> None:
    log_x = np.linspace(-1.0, 1.0, 20)
    log_y = 1.5 * log_x + 0.4
    error = np.full_like(log_x, 0.05)

    expected = powerlaw_negative_log_likelihood(
        [1.5, 0.4, -5.0],
        log_x,
        log_y,
        error,
        error,
    )
    incorrect = powerlaw_negative_log_likelihood(
        [0.5, -0.4, -5.0],
        log_x,
        log_y,
        error,
        error,
    )

    assert expected < incorrect


def test_current_catalogue_preserves_paper_correlation_slope() -> None:
    if not CATALOGUE.is_file():
        pytest.skip("Generated environment catalogue is unavailable")
    table = Table.read(CATALOGUE, format="ascii.ecsv")
    data = prepare_correlation_data(
        table,
        ["Surface_intercorebkg", "Interco_env_median"],
        log_columns=["Surface_intercorebkg", "Interco_env_median"],
    )
    result = pairwise_correlations(data)

    assert result["n_points"][0, 1] == 228
    assert result["pearson"][0, 1] == pytest.approx(
        0.7200655933,
        rel=1e-8,
    )
    assert result["spearman"][0, 1] == pytest.approx(
        0.7215132073,
        rel=1e-8,
    )


def test_correlation_task_inspection_reports_configured_profile() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/correlation.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "rows=231" in completed.stdout
    assert "x fractional uncertainty=0.33" in completed.stdout
    assert "valid rows=228/231" in completed.stdout
    assert "Pearson=0.720066" in completed.stdout


def test_publication_correlation_figure_is_available() -> None:
    figure_root = PROJECT_ROOT / "results" / "correlation" / "figures"
    outputs = [
        figure_root / "Correlation_Surface_env-Surface_core.pdf",
        figure_root / "Correlation_Surface_env-Surface_core.png",
    ]
    if not all(path.is_file() for path in outputs):
        pytest.skip("Generated publication correlation figure is unavailable")
    for path in outputs:
        assert path.stat().st_size > 10_000
