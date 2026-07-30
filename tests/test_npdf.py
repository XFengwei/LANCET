"""Tests for the ongoing regional N-PDF workflow."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from lancet.npdf import (
    logarithmic_column_density_contrast,
    npdf_histogram,
    positive_column_density_excess,
    shared_excess_mean,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_positive_excess_masks_invalid_and_nonpositive_values() -> None:
    result = positive_column_density_excess(
        np.array([np.nan, 1.0, 3.0, 5.0]),
        2.0,
    )
    np.testing.assert_allclose(result, [1.0, 3.0])


def test_shared_mean_uses_positive_excess_domain() -> None:
    regions = {
        "a": np.array([1.0, 3.0, 5.0]),
        "b": np.array([2.0, 4.0]),
    }
    result = shared_excess_mean(regions, 2.0)
    assert result == pytest.approx(np.mean([1.0, 3.0, 2.0]))
    eta = logarithmic_column_density_contrast(
        regions["a"],
        background=2.0,
        normalization=result,
    )
    np.testing.assert_allclose(eta, np.log(np.array([1.0, 3.0]) / result))


def test_npdf_histogram_integrates_to_one() -> None:
    edges = np.linspace(-2.0, 2.0, 9)
    result = npdf_histogram(
        [-1.5, -0.5, 0.0, 0.5, 1.5],
        edges,
    )
    integral = np.sum(result["density"] * result["widths"])
    assert integral == pytest.approx(1.0)
    assert np.sum(result["counts"]) == 5


def test_npdf_task_inspection_marks_assumptions_ongoing() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/n_pdf.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ONGOING / NOT PAPER-READY" in completed.stdout
    assert "Young: shape=(240, 240)" in completed.stdout
    assert "shared positive-excess mean: 3.97753e+22 cm^-2" in (
        completed.stdout
    )
