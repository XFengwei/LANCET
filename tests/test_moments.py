"""Tests for C18O moment-map utilities."""

from pathlib import Path
import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest

from lancet.moments import reproject_image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results" / "moment_maps"
RETAINED = (
    PROJECT_ROOT / "data" / "delivered" / "total_power" / "SEDIGISM"
)


def _header(shape: tuple[int, int], pixel_arcsec: float) -> fits.Header:
    return fits.Header(
        {
            "NAXIS": 2,
            "NAXIS1": shape[1],
            "NAXIS2": shape[0],
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRPIX1": (shape[1] + 1) / 2,
            "CRPIX2": (shape[0] + 1) / 2,
            "CRVAL1": 221.0,
            "CRVAL2": -59.0,
            "CDELT1": -pixel_arcsec / 3600.0,
            "CDELT2": pixel_arcsec / 3600.0,
            "BUNIT": "km s-1",
            "BMAJ": 30.0 / 3600.0,
            "BMIN": 30.0 / 3600.0,
            "BPA": 0.0,
        }
    )


def test_reprojection_uses_target_shape_and_preserves_unit() -> None:
    source = np.ones((5, 5), dtype=np.float32)
    output, footprint, header = reproject_image(
        source,
        _header((5, 5), 2.0),
        _header((9, 9), 1.0),
    )

    assert output.shape == (9, 9)
    assert footprint.shape == (9, 9)
    assert header["BUNIT"] == "km s-1"
    assert np.isfinite(output[4, 4])


def test_moment_task_inspection_is_read_only() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/c18o_moments.py",
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
    assert "results/moment_maps" in completed.stdout


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


def test_reprojected_moment2_uses_alma_shape_and_float32() -> None:
    actual_path = (
        RESULTS / "reprojected" / "c18o_moment2_on_alma_grid.fits"
    )
    if not actual_path.is_file():
        pytest.skip("C18O reprojected product has not been generated")
    actual_header = fits.getheader(actual_path)
    alma_header = fits.getheader(
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "alma"
        / "TP7MTM2TM1"
        / "cont"
        / "allscale_combined_v3_dust_ff.MJysr.rp.fits"
    )
    assert actual_header["BITPIX"] == -32
    assert actual_header["NAXIS1"] == alma_header["NAXIS1"]
    assert actual_header["NAXIS2"] == alma_header["NAXIS2"]
    assert actual_header["BUNIT"] == "km s-1"
