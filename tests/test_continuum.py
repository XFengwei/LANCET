"""Tests for continuum synthesis and SIMBA validation utilities."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import shutil
import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest

from lancet.continuum import (
    calibrate_simba,
    scale_power_law,
    smooth_circular_jybeam,
    synthesize_continuum,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LEGACY_SEDFITTER = (
    PROJECT_ROOT / "data" / "delivered" / "sedfitter"
)
FREEFREE_RESULTS = PROJECT_ROOT / "results" / "construct_freefree_model"
DUST_230 = (
    PROJECT_ROOT
    / "data"
    / "delivered"
    / "paper1"
    / "TP_SEDsynthetic230GHz_cont.fits"
)
FREEFREE_230 = (
    FREEFREE_RESULTS
    / "Te_8000K"
    / "ff230_map_smoothed_18arcsec_regrid2dust.fits"
)
CONTINUUM_230 = (
    LEGACY_SEDFITTER / "SEDsynthetic_Dust_FF_230GHz.fits"
)
CONTINUUM_250_SMOOTHED = (
    LEGACY_SEDFITTER
    / "SEDsynthetic_Dust_FF_250GHz_smooth24arcsec.fits"
)
SIMBA_REPROJECTED = (
    PROJECT_ROOT
    / "data"
    / "delivered"
    / "archive"
    / "infrared"
    / "SIMBA_1.2mm_newhdr_rp2sed.fits"
)
SCRIPT = PROJECT_ROOT / "scripts" / "synthesize_continuum.py"
SPEC = spec_from_file_location("synthesize_continuum_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SYNTHESIS_SCRIPT = module_from_spec(SPEC)
sys.modules[SPEC.name] = SYNTHESIS_SCRIPT
SPEC.loader.exec_module(SYNTHESIS_SCRIPT)


def test_power_law_scaling() -> None:
    image = np.array([[1.0, 2.0]])
    scaled = scale_power_law(
        image,
        reference_frequency_ghz=100.0,
        target_frequency_ghz=200.0,
        spectral_index=2.0,
    )

    np.testing.assert_allclose(scaled, [[4.0, 8.0]])


def test_synthesis_preserves_dust_mask_and_fills_missing_freefree() -> None:
    dust = np.array([[1.0, np.nan], [2.0, 3.0]])
    freefree = np.array([[0.5, 4.0], [np.nan, 1.0]])

    dust_scaled, freefree_scaled, total = synthesize_continuum(
        dust,
        freefree,
    )

    assert dust_scaled[0, 0] == 1.0
    assert np.isnan(dust_scaled[0, 1])
    assert freefree_scaled[1, 0] == 0.0
    assert total[1, 0] == 2.0
    assert np.isnan(total[0, 1])


def test_paper2_synthesis_reproduces_reference_maps() -> None:
    dust = fits.getdata(DUST_230, memmap=True)
    freefree = fits.getdata(FREEFREE_230, memmap=True)

    _, _, model_230 = synthesize_continuum(
        dust,
        freefree,
        target_frequency_ghz=230.0,
    )
    _, _, model_250 = synthesize_continuum(
        dust,
        freefree,
        target_frequency_ghz=250.0,
    )

    expected_230 = fits.getdata(
        LEGACY_SEDFITTER / "SEDsynthetic_Dust_FF_230GHz.fits",
        memmap=True,
    )
    expected_250 = fits.getdata(
        LEGACY_SEDFITTER / "SEDsynthetic_Dust_FF_250GHz.fits",
        memmap=True,
    )
    np.testing.assert_allclose(
        model_230,
        expected_230,
        rtol=2e-7,
        atol=3e-7,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        model_250,
        expected_250,
        rtol=2e-7,
        atol=3e-7,
        equal_nan=True,
    )


def test_250_ghz_smoothing_exactly_reproduces_legacy_map() -> None:
    model_250, header = fits.getdata(
        LEGACY_SEDFITTER / "SEDsynthetic_Dust_FF_250GHz.fits",
        header=True,
        memmap=True,
    )
    smoothed = smooth_circular_jybeam(
        model_250,
        pixel_scale_arcsec=abs(float(header["CDELT1"])) * 3600.0,
        original_beam_arcsec=18.0,
        target_beam_arcsec=24.0,
    )
    expected = fits.getdata(
        LEGACY_SEDFITTER
        / "SEDsynthetic_Dust_FF_250GHz_smooth24arcsec.fits",
        memmap=True,
    )

    assert np.array_equal(
        smoothed.astype(np.float32),
        expected,
        equal_nan=True,
    )


def test_simba_calibration_exactly_reproduces_legacy_array() -> None:
    raw = fits.getdata(
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "archive"
        / "infrared"
        / "SIMBA_1.2mm.fits",
        memmap=True,
    )
    calibrated = calibrate_simba(raw)
    expected = fits.getdata(
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "archive"
        / "infrared"
        / "SIMBA_1.2mm_newhdr.fits",
        memmap=True,
    )

    assert np.array_equal(
        calibrated.astype(np.float32),
        expected,
        equal_nan=True,
    )


def test_continuum_task_inspects_declared_inputs() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/synthesize_continuum.py",
            "--stage",
            "inspect",
            "--freefree",
            str(FREEFREE_230),
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Model readiness: ready" in completed.stdout
    assert "dust=3.8, free-free=-0.1" in completed.stdout


def test_plot_stage_writes_pdf_and_png_products(tmp_path: Path) -> None:
    model_directory = tmp_path / "model"
    validation_directory = tmp_path / "validation"
    model_directory.mkdir()
    validation_directory.mkdir()

    for source, destination in (
        (DUST_230, model_directory / "dust_230GHz.fits"),
        (FREEFREE_230, model_directory / "freefree_230GHz.fits"),
        (CONTINUUM_230, model_directory / "continuum_230GHz.fits"),
        (
            CONTINUUM_250_SMOOTHED,
            validation_directory
            / "continuum_250GHz_smoothed_24arcsec.fits",
        ),
        (
            SIMBA_REPROJECTED,
            validation_directory / "simba_250GHz_on_model_grid.fits",
        ),
    ):
        shutil.copyfile(source, destination)

    config = SYNTHESIS_SCRIPT.load_task_config(
        PROJECT_ROOT / "tasks" / "synthesize_continuum" / "task.toml"
    )
    outputs = SYNTHESIS_SCRIPT.plot_products(
        {"output_root": tmp_path},
        config,
    )

    assert len(outputs) == 6
    assert {path.suffix for path in outputs} == {".pdf", ".png"}
    for path in outputs:
        assert path.is_file()
        assert path.stat().st_size > 10_000
