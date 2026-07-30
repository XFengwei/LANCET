"""Tests for reusable continuum-image preparation."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest

from lancet.imaging import (
    apply_primary_beam_mask,
    gaussian_beam_solid_angle_sr,
    jybeam_to_mjysr,
    jybeam_to_mjysr_from_header,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ALMA_7MTM2 = (
    PROJECT_ROOT / "data" / "delivered" / "alma" / "7MTM2" / "cont"
)
BASE_IMAGE = (
    ALMA_7MTM2
    / "miriad"
    / "IRDC316.8.7MTM2.robust0.5.clean.dropdeg.fits"
)
MASKED_IMAGE = (
    ALMA_7MTM2
    / "miriad"
    / "IRDC316.8.7MTM2.robust0.5.clean.dropdeg.pbcut.fits"
)
CONVERTED_IMAGE = (
    ALMA_7MTM2
    / "miriad"
    / "IRDC316.8.7MTM2.robust0.5.clean.dropdeg.pbcut.MJysr.fits"
)
PRIMARY_BEAM = ALMA_7MTM2 / "casa" / "primary_beam_rp2miriad.fits"
SCRIPT = PROJECT_ROOT / "scripts" / "prepare_continuum_images.py"
SPEC = spec_from_file_location("prepare_continuum_images_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
PREPARATION_SCRIPT = module_from_spec(SPEC)
sys.modules[SPEC.name] = PREPARATION_SCRIPT
SPEC.loader.exec_module(PREPARATION_SCRIPT)


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


def test_unit_conversion_exactly_reproduces_delivered_product() -> None:
    source, header = fits.getdata(MASKED_IMAGE, header=True, memmap=True)
    expected = fits.getdata(CONVERTED_IMAGE, memmap=True)
    converted = jybeam_to_mjysr_from_header(source, header)

    assert np.array_equal(converted, expected, equal_nan=True)


def test_primary_beam_mask_exactly_reproduces_delivered_product() -> None:
    source = fits.getdata(BASE_IMAGE, memmap=True)
    primary_beam = fits.getdata(PRIMARY_BEAM, memmap=True)
    expected = fits.getdata(MASKED_IMAGE, memmap=True)
    masked = apply_primary_beam_mask(source, primary_beam)

    assert np.array_equal(
        masked.astype(np.float32),
        expected,
        equal_nan=True,
    )


def test_unit_conversion_writer_uses_float32_storage(tmp_path: Path) -> None:
    source = tmp_path / "source.fits"
    output = tmp_path / "converted.fits"
    header = fits.Header(
        {
            "BMAJ": 2.0 / 3600.0,
            "BMIN": 1.0 / 3600.0,
            "BUNIT": "Jy/beam",
        }
    )
    fits.writeto(
        source,
        np.array([[1.0, np.nan]], dtype=np.float64),
        header,
    )

    PREPARATION_SCRIPT.run_unit_conversions(
        [{"input": source, "output": output}],
    )

    stored, stored_header = fits.getdata(output, header=True)
    expected = jybeam_to_mjysr_from_header(
        fits.getdata(source),
        fits.getheader(source),
    )
    assert stored_header["BITPIX"] == -32
    assert stored.dtype.itemsize == 4
    np.testing.assert_array_equal(
        stored,
        expected.astype(np.float32),
    )


def test_image_preparation_task_inspection_is_read_only() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/prepare_continuum_images.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Primary-beam masks: 4 at response >= 0.2" in completed.stdout
    assert "Jy/beam to MJy/sr conversions: 11" in completed.stdout
    assert "7MTM2_pbcut: task_output" in completed.stdout


def test_generated_floating_point_fits_use_float32_storage() -> None:
    violations: list[str] = []
    for path in (PROJECT_ROOT / "results").rglob("*.fits"):
        bitpix = int(fits.getheader(path)["BITPIX"])
        if bitpix < 0 and bitpix != -32:
            violations.append(
                f"{path.relative_to(PROJECT_ROOT)}: BITPIX={bitpix}"
            )

    assert not violations, "Float64 task outputs found:\n" + "\n".join(
        violations
    )
