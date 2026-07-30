"""Tests for the Paper II free-free model task."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

from astropy.io import fits
import numpy as np
import pytest

from lancet.freefree import (
    beam_solid_angle,
    fit_freefree_map,
    fit_freefree_pixel,
    fractional_flux_difference,
    freefree_flux_density,
    meerkat_subband_frequencies,
    prepare_sed_cube,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FREEFREE_RESULTS = PROJECT_ROOT / "results" / "construct_freefree_model"
MEERKAT_CUBE = (
    PROJECT_ROOT
    / "data"
    / "delivered"
    / "archive"
    / "radio"
    / "MeerKAT_subcube_reprojected.fits"
)
SCRIPT = PROJECT_ROOT / "scripts" / "construct_freefree_model.py"
SPEC = spec_from_file_location("construct_freefree_model_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
FREEFREE_SCRIPT = module_from_spec(SPEC)
sys.modules[SPEC.name] = FREEFREE_SCRIPT
SPEC.loader.exec_module(FREEFREE_SCRIPT)


def test_beam_solid_angle_matches_gaussian_definition() -> None:
    omega = beam_solid_angle(8.0, 8.0)
    arcsec_to_radian = np.pi / (180.0 * 3600.0)
    expected = (
        np.pi * (8.0 * arcsec_to_radian) ** 2 / (4.0 * np.log(2.0))
    )
    assert omega == pytest.approx(expected)

    with pytest.raises(ValueError, match="positive"):
        beam_solid_angle(0.0, 8.0)


def test_exact_meerkat_frequencies_are_read_from_delivered_header() -> None:
    header = fits.getheader(MEERKAT_CUBE)
    frequency = meerkat_subband_frequencies(header)

    assert frequency.shape == (14,)
    assert frequency[0] == pytest.approx(0.9081416015625)
    assert frequency[-1] == pytest.approx(1.656723632812)
    assert np.all(np.diff(frequency) > 0)


def test_pixel_fit_recovers_synthetic_emission_measure() -> None:
    frequency = np.array([0.9, 1.1, 1.4, 5.5])
    omega = beam_solid_angle(8.0, 8.0)
    expected_emission_measure = 2.5e6
    flux = freefree_flux_density(
        frequency,
        expected_emission_measure,
        electron_temperature_k=8000.0,
        beam_solid_angle_sr=omega,
    )
    uncertainty = np.full(4, 1e-4)

    emission_measure, _, chi_squared, _, model = fit_freefree_pixel(
        frequency,
        flux,
        uncertainty,
        electron_temperature_k=8000.0,
        beam_solid_angle_sr=omega,
    )

    assert emission_measure == pytest.approx(
        expected_emission_measure,
        rel=1e-8,
    )
    assert chi_squared == pytest.approx(0.0, abs=1e-15)
    assert model == pytest.approx(flux)


def test_fractional_flux_difference_masks_invalid_reference_pixels() -> None:
    comparison = np.array([[3.0, 1.0, 4.0], [1.0, np.nan, 2.0]])
    reference = np.array([[2.0, 0.0, np.nan], [2.0, 1.0, 2.0]])

    signed = fractional_flux_difference(comparison, reference)
    absolute = fractional_flux_difference(
        comparison,
        reference,
        absolute=True,
    )

    assert signed[0, 0] == 0.5
    assert signed[1, 0] == -0.5
    assert absolute[1, 0] == 0.5
    assert np.isnan(signed[0, 1])
    assert np.isnan(signed[0, 2])
    assert np.isnan(signed[1, 1])


def test_prepare_sed_cube_preserves_band_order() -> None:
    cube = np.zeros((1, 16, 2, 3), dtype=float)
    cube[0, 0] = 100.0
    cube[0, 1] = 200.0
    for index in range(14):
        cube[0, index + 2] = index + 1
    cornish = np.full((2, 3), 99.0)

    mfs, alpha, sed = prepare_sed_cube(cube, cornish)

    assert np.all(mfs == 100.0)
    assert np.all(alpha == 200.0)
    assert sed.shape == (15, 2, 3)
    assert np.all(sed[:14, 0, 0] == np.arange(1, 15))
    assert np.all(sed[14] == 99.0)


def test_map_fit_respects_mask_and_minimum_band_count() -> None:
    frequency = np.array([0.9, 1.1, 1.4])
    rms = np.full(3, 1e-5)
    omega = beam_solid_angle(8.0, 8.0)
    spectrum = freefree_flux_density(
        frequency,
        1e6,
        beam_solid_angle_sr=omega,
    )
    sed = np.zeros((3, 1, 2))
    sed[:, 0, 0] = spectrum
    sed[:, 0, 1] = spectrum

    results = fit_freefree_map(
        sed,
        frequency,
        rms,
        beam_solid_angle_sr=omega,
        fit_mask=np.array([[True, False]]),
        minimum_bands=3,
    )

    assert results["EM_map"][0, 0] == pytest.approx(1e6, rel=1e-5)
    assert np.isnan(results["EM_map"][0, 1])
    assert results["nband_map"].tolist() == [[3, 0]]


def test_paper2_profile_reproduces_reference_13_band_pixel() -> None:
    with (
        PROJECT_ROOT
        / "tasks"
        / "construct_freefree_model"
        / "task.toml"
    ).open("rb") as handle:
        fit = tomllib.load(handle)["fit"]

    y_pixel, x_pixel = 497, 383
    cube = fits.getdata(MEERKAT_CUBE, memmap=True)[0]
    cornish = fits.getdata(
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "archive"
        / "radio"
        / "CORNISH_5p5GHz_8arcsec_on_MeerKATgrid.fits",
        memmap=True,
    )
    expected_emission_measure = 1_454_721.875

    spectrum = np.concatenate(
        [
            np.asarray(cube[2:16, y_pixel, x_pixel], dtype=float),
            [float(cornish[y_pixel, x_pixel])],
        ]
    )
    frequency = np.concatenate(
        [np.asarray(fit["meerkat_frequencies_ghz"], dtype=float), [5.5]]
    )
    rms = np.concatenate(
        [
            np.asarray(fit["rms_meerkat_jy_per_beam"], dtype=float),
            [float(fit["rms_cornish_jy_per_beam"])],
        ]
    )
    valid = (
        np.isfinite(spectrum)
        & (spectrum > 0)
        & (spectrum > float(fit["snr_threshold"]) * rms)
    )
    emission_measure = fit_freefree_pixel(
        frequency[valid],
        spectrum[valid],
        rms[valid],
        electron_temperature_k=8000.0,
        beam_solid_angle_sr=beam_solid_angle(8.0, 8.0),
    )[0]

    assert np.count_nonzero(valid) == 13
    assert emission_measure == pytest.approx(
        expected_emission_measure,
        rel=1e-7,
    )


def test_task_inspection_reports_paper_profile_without_writing() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/construct_freefree_model.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "0.908141602" in completed.stdout
    assert "1.656723633" in completed.stdout
    assert "Frequencies selected by task profile (configured)" in completed.stdout
    assert "0.908000000" in completed.stdout
    assert "Fit readiness: ready (15 positive per-band RMS values configured)" in (
        completed.stdout
    )


def test_temperature_comparison_plot_writes_pdf_and_png(
    tmp_path: Path,
) -> None:
    for temperature in (6000, 8000, 10000):
        directory = tmp_path / f"Te_{temperature}K"
        directory.mkdir()
        shutil.copyfile(
            FREEFREE_RESULTS
            / f"Te_{temperature}K"
            / "ff230_map.fits",
            directory / "ff230_map.fits",
        )

    with (
        PROJECT_ROOT
        / "tasks"
        / "construct_freefree_model"
        / "task.toml"
    ).open("rb") as handle:
        config = tomllib.load(handle)
    outputs = FREEFREE_SCRIPT.plot_temperature_flux_comparison(
        {"output_root": tmp_path},
        config,
    )

    assert len(outputs) == 2
    assert {path.suffix for path in outputs} == {".pdf", ".png"}
    for path in outputs:
        assert path.is_file()
        assert path.stat().st_size > 10_000


def test_fit_result_writer_enforces_compact_fits_storage(
    tmp_path: Path,
) -> None:
    shape = (2, 2)
    results = {
        "EM_map": np.ones(shape, dtype=np.float64),
        "EMerr_map": np.ones(shape, dtype=np.float64),
        "ff230_map": np.ones(shape, dtype=np.float64),
        "chi2_map": np.ones(shape, dtype=np.float64),
        "rchi2_map": np.ones(shape, dtype=np.float64),
        "nband_map": np.ones(shape, dtype=np.float64),
        "fitmask_map": np.ones(shape, dtype=np.float64),
        "model_cube": np.ones((15, *shape), dtype=np.float64),
    }
    written = FREEFREE_SCRIPT.save_fit_results(
        results,
        fits.Header(),
        np.linspace(1.0, 5.5, 15),
        tmp_path,
        electron_temperature_k=8000.0,
    )

    for path in written:
        bitpix = int(fits.getheader(path)["BITPIX"])
        if path.name in {"nband_map.fits", "fitmask_map.fits"}:
            assert bitpix == 16
        else:
            assert bitpix == -32
