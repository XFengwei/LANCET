"""Tests for elliptical-Gaussian source deconvolution."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from lancet.deconvolution import (
    UnresolvedSourceError,
    convolve_gaussians,
    deconvolve_gaussian,
)


def _pa_difference(first: float, second: float) -> float:
    return abs((first - second + 90.0) % 180.0 - 90.0)


def test_reference_example_matches_expected_deconvolution() -> None:
    result = deconvolve_gaussian(
        0.7747,
        0.6559,
        63.38,
        0.712,
        0.606,
        47.34,
    )

    assert result[0] == pytest.approx(0.351105983847, rel=1e-10)
    assert result[1] == pytest.approx(0.181409724402, rel=1e-10)
    assert _pa_difference(result[2], 90.9774502434) < 1e-9


def test_aligned_ellipses_reduce_to_quadrature_subtraction() -> None:
    result = deconvolve_gaussian(5.0, 3.0, 30.0, 4.0, 2.0, 30.0)

    assert result[0] == pytest.approx(3.0)
    assert result[1] == pytest.approx(np.sqrt(5.0))
    assert _pa_difference(result[2], 30.0) < 1e-10


def test_reconvolution_recovers_the_observed_ellipse() -> None:
    observed = (0.7747, 0.6559, 63.38)
    beam = (0.712, 0.606, 47.34)
    source = deconvolve_gaussian(*observed, *beam)
    reconstructed = convolve_gaussians(*source, *beam)

    assert reconstructed[:2] == pytest.approx(observed[:2], rel=1e-10)
    assert _pa_difference(reconstructed[2], observed[2]) < 1e-9


def test_reversed_axes_are_normalized_with_a_pa_rotation() -> None:
    normal = deconvolve_gaussian(5.0, 3.0, 30.0, 4.0, 2.0, 30.0)
    reversed_axes = deconvolve_gaussian(
        3.0, 5.0, 120.0, 2.0, 4.0, 120.0
    )

    assert reversed_axes[:2] == pytest.approx(normal[:2])
    assert _pa_difference(reversed_axes[2], normal[2]) < 1e-10


def test_unresolved_source_can_return_nan_or_raise() -> None:
    result = deconvolve_gaussian(1.0, 1.0, 0.0, 2.0, 2.0, 0.0)
    assert all(np.isnan(value) for value in result)

    with pytest.raises(UnresolvedSourceError):
        deconvolve_gaussian(
            1.0,
            1.0,
            0.0,
            2.0,
            2.0,
            0.0,
            on_unresolved="raise",
        )


def test_invalid_axes_are_rejected() -> None:
    with pytest.raises(ValueError, match="positive"):
        deconvolve_gaussian(0.0, 1.0, 0.0, 0.5, 0.5, 0.0)


def test_diagnostic_script_runs() -> None:
    project_root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/source_deconvolution.py",
        ],
        cwd=project_root,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Deconvolved:" in completed.stdout
    assert "Reconstruction residual:" in completed.stdout
