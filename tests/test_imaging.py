"""Tests for continuum-image unit conversion and primary-beam masking."""

from pathlib import Path

from astropy.io import fits
import numpy as np
import pytest

from lancet.imaging import (
    apply_primary_beam_mask,
    frequency_label,
    gaussian_beam_solid_angle_sr,
    jybeam_to_mjysr,
    jybeam_to_mjysr_from_header,
    read_frequency_ghz,
    read_mjysr,
)


def test_gaussian_beam_area_and_unit_conversion() -> None:
    major = 2.0 / 3600.0
    minor = 1.0 / 3600.0
    area = gaussian_beam_solid_angle_sr(major, minor)
    expected_area = (
        np.pi
        * np.deg2rad(major)
        * np.deg2rad(minor)
        / (4.0 * np.log(2.0))
    )

    assert area == pytest.approx(expected_area)
    converted = jybeam_to_mjysr(
        np.array([1.0]),
        major_fwhm_deg=major,
        minor_fwhm_deg=minor,
    )
    assert converted[0] == pytest.approx(1e-6 / expected_area)


def test_primary_beam_mask_validates_shape_and_threshold() -> None:
    image = np.arange(4, dtype=float).reshape(2, 2)
    primary_beam = np.array([[0.1, 0.2], [0.5, np.nan]])
    masked = apply_primary_beam_mask(
        image,
        primary_beam,
        minimum_response=0.2,
    )

    assert np.isnan(masked[0, 0])
    assert masked[0, 1] == 1.0
    assert masked[1, 0] == 2.0
    assert np.isnan(masked[1, 1])
    with pytest.raises(ValueError, match="shape"):
        apply_primary_beam_mask(image, np.ones((3, 3)))


def test_unit_conversion_from_header_uses_the_restoring_beam() -> None:
    major = 0.712 / 3600.0
    minor = 0.606 / 3600.0
    header = fits.Header({"BMAJ": major, "BMIN": minor, "BUNIT": "Jy/beam"})
    image = np.array([[1.0, -2.0], [np.nan, 0.5]])

    converted = jybeam_to_mjysr_from_header(image, header)

    factor = 1e-6 / gaussian_beam_solid_angle_sr(major, minor)
    np.testing.assert_allclose(converted, image * factor, equal_nan=True)


def test_read_mjysr_converts_jy_per_beam_images(tmp_path: Path) -> None:
    header = fits.Header({"BMAJ": 2.0 / 3600, "BMIN": 1.0 / 3600, "BUNIT": "Jy/beam"})
    path = tmp_path / "image.fits"
    fits.writeto(path, np.full((1, 1, 2, 3), 2.0, dtype=np.float32), header)

    data, out_header = read_mjysr(path)

    expected = jybeam_to_mjysr_from_header(np.full((2, 3), 2.0), header)
    assert data.shape == (2, 3)
    assert data.dtype == np.float32
    np.testing.assert_allclose(data, expected, rtol=1e-6)
    assert out_header["BUNIT"] == "MJy/sr"


def test_read_mjysr_passes_mjysr_and_rejects_other_units(tmp_path: Path) -> None:
    path = tmp_path / "mjysr.fits"
    fits.writeto(path, np.ones((2, 2), dtype=np.float32), fits.Header({"BUNIT": "MJy/sr"}))
    data, _ = read_mjysr(path)
    assert np.all(data == 1.0)

    other = tmp_path / "kelvin.fits"
    fits.writeto(other, np.ones((2, 2), dtype=np.float32), fits.Header({"BUNIT": "K"}))
    with pytest.raises(ValueError, match="BUNIT"):
        read_mjysr(other)


def test_read_frequency_ghz(tmp_path):
    path = tmp_path / "tp.fits"
    fits.PrimaryHDU(np.zeros((2, 2), np.float32)).writeto(path)
    with pytest.raises(ValueError, match="FREQGHZ"):
        read_frequency_ghz(path)
    fits.setval(path, "FREQGHZ", value=223.8)
    assert read_frequency_ghz(path) == 223.8
    assert frequency_label(223.8) == "224GHz"
