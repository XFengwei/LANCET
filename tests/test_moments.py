"""Tests for C18O moment-map utilities."""

from pathlib import Path
import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results" / "c18o_moments"
RETAINED = (
    PROJECT_ROOT / "data" / "delivered" / "total_power" / "SEDIGISM"
)


def test_task_has_no_reprojection() -> None:
    import tomllib

    config = tomllib.loads(
        (PROJECT_ROOT / "tasks" / "c18o_moments" / "task.toml").read_text()
    )

    assert set(config["inputs"]) == {"c18o_cube"}
    assert "reprojection" not in config


def test_moment_task_inspection_is_read_only() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "tasks/c18o_moments/run.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "c18o_2-1_apex-sedigism.fits" in completed.stdout
    assert "snr_hi=9.0, snr_lo=2.0" in completed.stdout
    assert "results/c18o_moments" in completed.stdout


def test_reconstructed_cube_derivatives_follow_raw_cube() -> None:
    velocity = RESULTS / "cubes" / "c18o_jybeam_velocity_kms.fits"
    temperature = (
        RESULTS / "cubes" / "c18o_brightness_temperature.fits"
    )
    if not velocity.is_file() or not temperature.is_file():
        pytest.skip("C18O cube derivatives have not been generated")

    raw_cube = RETAINED / "c18o_2-1_apex-sedigism.fits"
    actual_velocity = fits.getdata(velocity, memmap=True)
    raw_data = fits.getdata(
        raw_cube,
        memmap=True,
    )
    np.testing.assert_array_equal(actual_velocity, raw_data)
    assert fits.getheader(velocity)["BITPIX"] == -32
    assert fits.getheader(velocity)["CUNIT3"] == "km s-1"

    assert fits.getdata(temperature, memmap=True).shape == raw_data.shape
    assert fits.getheader(temperature)["BITPIX"] == -32
    assert fits.getheader(temperature)["BUNIT"] == "K"
    assert fits.getheader(temperature)["CUNIT3"] == "m/s"


def test_reconstructed_mask_and_moments_are_float32() -> None:
    generated_root = RESULTS / "moments"
    if not (generated_root / "c18o.mom2.fits.gz").is_file():
        pytest.skip("C18O moment products have not been generated")
    retained_root = RETAINED / "Moments"

    products = ("mask", "mom0", "mom1", "mom2", "emom0", "emom1", "emom2")
    for product in products:
        actual_path = generated_root / f"c18o.{product}.fits.gz"
        assert actual_path.is_file()
        assert fits.getheader(actual_path)["BITPIX"] == -32

    retained_mask = retained_root / "TP_C18O.mask.fits.gz"
    if retained_mask.is_file():
        np.testing.assert_array_equal(
            fits.getdata(generated_root / "c18o.mask.fits.gz"),
            fits.getdata(retained_mask),
        )
