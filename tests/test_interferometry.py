"""Tests for interferometric coverage and image-comparison utilities."""

import os
from pathlib import Path
import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest

from lancet.interferometry import (
    beam_match_and_reproject,
    deduplicate_uv_coordinates,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _image_header(beam_arcsec: float) -> fits.Header:
    return fits.Header(
        {
            "NAXIS": 2,
            "NAXIS1": 51,
            "NAXIS2": 51,
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRPIX1": 26.0,
            "CRPIX2": 26.0,
            "CRVAL1": 221.0,
            "CRVAL2": -59.0,
            "CDELT1": -0.2 / 3600.0,
            "CDELT2": 0.2 / 3600.0,
            "BMAJ": beam_arcsec / 3600.0,
            "BMIN": beam_arcsec / 3600.0,
            "BPA": 0.0,
            "BUNIT": "MJy/sr",
        }
    )


def test_uv_deduplication_keeps_first_coordinate_in_each_bin() -> None:
    u_values = np.array([1.01, 1.04, 2.0, 2.04])
    v_values = np.array([3.01, 3.04, 4.0, 4.04])

    u_unique, v_unique = deduplicate_uv_coordinates(
        u_values,
        v_values,
        decimals=1,
    )

    np.testing.assert_array_equal(u_unique, [1.01, 2.0])
    np.testing.assert_array_equal(v_unique, [3.01, 4.0])


def test_beam_match_uses_reference_shape_wcs_and_beam() -> None:
    image = np.zeros((51, 51), dtype=float)
    image[25, 25] = 1.0
    reference = np.zeros((51, 51), dtype=float)
    image_header = _image_header(1.0)
    reference_header = _image_header(2.0)

    matched, output_header = beam_match_and_reproject(
        image,
        image_header,
        reference,
        reference_header,
    )

    assert matched.shape == reference.shape
    assert np.isfinite(matched[25, 25])
    assert output_header["BMAJ"] * 3600.0 == pytest.approx(2.0)
    assert output_header["BMIN"] * 3600.0 == pytest.approx(2.0)
    assert output_header["BUNIT"] == "MJy/sr"


def test_array_comparison_inspection_is_read_only() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/combination_diagnostic.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "IRDC316.8.TM1TM2.bin.combselfcal" in completed.stdout
    assert "combined_pbcut: shape=(5127, 7096), BITPIX=-32" in (
        completed.stdout
    )


def test_uv_reader_does_not_create_casa_log(
    tmp_path: Path,
) -> None:
    measurement_set = tmp_path / "empty.ms"
    measurement_set.mkdir()
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    code = (
        "import sys\n"
        "from lancet.interferometry import read_unique_uv_coordinates\n"
        "try:\n"
        "    read_unique_uv_coordinates('empty.ms')\n"
        "except Exception:\n"
        "    pass\n"
        "assert 'casatools' in sys.modules\n"
    )
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )

    assert list(tmp_path.glob("casa-*.log")) == []


def test_generated_beam_match_agrees_with_retained_product() -> None:
    generated = (
        PROJECT_ROOT
        / "results"
        / "combination_diagnostic"
        / "beam_match"
        / "combined_convolved_to_ACA.fits"
    )
    if not generated.is_file():
        pytest.skip("Beam-match task product has not been generated")
    retained = (
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "alma"
        / "7MTM2TM1"
        / "cont"
        / "miriad"
        / "IRDC316.8.7MTM2TM1.robust0.5.clean.drogdeg.pbcut.MJysr.convolve2ACA.fits"
    )
    actual = fits.getdata(generated, memmap=True)
    expected = fits.getdata(retained, memmap=True)

    assert fits.getheader(generated)["BITPIX"] == -32
    np.testing.assert_array_equal(np.isfinite(actual), np.isfinite(expected))
    np.testing.assert_allclose(
        actual,
        expected,
        rtol=2e-5,
        atol=3e-5,
        equal_nan=True,
    )
