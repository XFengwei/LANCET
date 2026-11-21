"""Tests for the Planck polarization task."""

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


@pytest.mark.local_inputs("data")
def test_planck_polarization_inspection_reports_inputs() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "tasks/planck_polarization/run.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ongoing" not in completed.stdout.lower()
    assert "provisional" not in completed.stdout.lower()
    assert "common shape: (120, 120)" in completed.stdout
    assert "no field-strength estimate" in completed.stdout


def _load_task():
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location(
        "planck_polarization_script",
        PROJECT_ROOT / "tasks" / "planck_polarization" / "run.py",
    )
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_vector_step_is_half_beam_in_whole_pixels() -> None:
    task = _load_task()

    assert task.vector_step_pixels(4.82, 0.5, 1.0) == 2
    assert task.vector_step_pixels(4.82, 1.0, 1.0) == 5
    assert task.vector_step_pixels(1.0, 0.5, 1.0) == 1


def test_view_limits_give_a_centred_box() -> None:
    from astropy.wcs import WCS

    task = _load_task()
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]
    wcs.wcs.crval = [316.7686, -0.003]
    wcs.wcs.crpix = [60.5, 60.5]
    wcs.wcs.cdelt = [-1 / 60, 1 / 60]

    (x0, x1), (y0, y1) = task.view_limits(wcs, 316.7686, -0.003, 30.0)

    assert x1 - x0 == pytest.approx(30.0, abs=0.05)
    assert y1 - y0 == pytest.approx(30.0, abs=0.05)
    assert (x0 + x1) / 2 == pytest.approx(59.5, abs=0.05)


def test_plot_zooms_to_thirty_arcmin_with_half_beam_vectors() -> None:
    import tomllib

    plot = tomllib.loads(
        (PROJECT_ROOT / "tasks" / "planck_polarization" / "task.toml").read_text()
    )["plot"]

    assert plot["view_size_arcmin"] == 30.0
    assert plot["vector_spacing_beam"] == 0.5
    assert plot["beam_fwhm_arcmin"] == 4.82


def test_task_has_no_ongoing_or_provisional_wording() -> None:
    folder = PROJECT_ROOT / "tasks" / "planck_polarization"
    for path in folder.iterdir():
        text = path.read_text().lower()
        assert "ongoing" not in text, path.name
        assert "provisional" not in text, path.name
        assert "paper_ready" not in text, path.name
