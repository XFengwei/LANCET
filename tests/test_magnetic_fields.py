"""Tests for the ongoing magnetic-field diagnostics."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from lancet.magnetic_fields import (
    debiased_polarized_intensity,
    magnetic_field_angle_degrees,
    planck_polarization_products,
    polarization_angle_degrees,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_debiased_polarized_intensity_is_finite_and_below_raw() -> None:
    q = np.array([0.0, 3.0, 4.0])
    u = np.array([0.0, 4.0, 3.0])
    result = debiased_polarized_intensity(q, u, 1.0)
    raw = np.hypot(q, u)
    assert result[0] == 0
    assert np.all(result >= 0)
    assert np.all(result[1:] < raw[1:])


def test_polarization_and_field_angles_use_ninety_degree_rotation() -> None:
    q = np.array([1.0])
    u = np.array([0.0])
    polarization = polarization_angle_degrees(q, u)
    magnetic = magnetic_field_angle_degrees(q, u)
    assert polarization[0] == pytest.approx(0.0)
    assert magnetic[0] == pytest.approx(90.0)


def test_planck_products_apply_independent_signal_masks() -> None:
    intensity = np.array([[0.01, 0.04]])
    q = np.array([[0.0, 0.003]])
    u = np.array([[0.0, 0.004]])
    products = planck_polarization_products(
        intensity,
        q,
        u,
        rms_i=0.01,
        rms_qu=0.001,
    )
    np.testing.assert_array_equal(
        products["intensity_snr3_mask"],
        [[False, True]],
    )
    assert not np.isfinite(products["magnetic_field_angle"][0, 0])
    assert np.isfinite(products["magnetic_field_angle"][0, 1])
    assert np.isfinite(
        products["polarization_fraction_debiased"][0, 1]
    )


def test_magnetic_field_task_inspection_marks_workflow_ongoing() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/magnetic_fields.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ONGOING / NOT PAPER-READY" in completed.stdout
    assert "common shape: (120, 120)" in completed.stdout
    assert "no field-strength estimate" in completed.stdout
