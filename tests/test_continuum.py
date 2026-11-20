"""Tests for continuum synthesis and SIMBA validation utilities."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest

from lancet.column_density import modified_blackbody_mjysr
from lancet.continuum import (
    calibrate_simba,
    dust_model_jybeam,
    dust_opacity_power_law,
    fill_missing,
    scale_power_law,
    smooth_circular_jybeam,
    synthesize_continuum,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DUST = (
    PROJECT_ROOT
    / "data"
    / "delivered"
    / "sedfitter"
    / "Dustmodel_223.8GHz.fits"
)
FREEFREE = (
    PROJECT_ROOT
    / "results"
    / "construct_freefree_model"
    / "Te_8000K"
    / "ff224_map_smoothed_18arcsec_regrid2dust.fits"
)
SIMBA = (
    PROJECT_ROOT
    / "data"
    / "delivered"
    / "archive"
    / "infrared"
    / "SIMBA_1.2mm.fits"
)
SCRIPT = PROJECT_ROOT / "tasks" / "synthesize_continuum" / "run.py"
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
        reference_frequency_ghz=223.8,
        target_frequency_ghz=223.8,
    )

    assert dust_scaled[0, 0] == 1.0
    assert np.isnan(dust_scaled[0, 1])
    assert freefree_scaled[1, 0] == 0.0
    assert total[1, 0] == 2.0
    assert np.isnan(total[0, 1])


def test_synthesis_requires_explicit_frequencies() -> None:
    with pytest.raises(TypeError, match="frequency_ghz"):
        synthesize_continuum(np.ones((2, 2)), np.ones((2, 2)))


def test_task_model_frequency_is_the_alma_continuum_frequency() -> None:
    config = SYNTHESIS_SCRIPT.load_task_config(
        PROJECT_ROOT / "tasks" / "synthesize_continuum" / "task.toml"
    )

    assert config["model"]["reference_frequency_ghz"] == 223.8
    assert config["model"]["frequencies_ghz"] == [223.8, 250.0]
    assert SYNTHESIS_SCRIPT.frequency_label(223.8) == "224GHz"
    assert SYNTHESIS_SCRIPT.frequency_label(250.0) == "250GHz"


def test_task_model_adds_unscaled_components_at_model_frequency(
    tmp_path: Path,
) -> None:
    config = SYNTHESIS_SCRIPT.load_task_config(
        PROJECT_ROOT / "tasks" / "synthesize_continuum" / "task.toml"
    )
    written = SYNTHESIS_SCRIPT.build_models(
        {"dust": DUST, "freefree": FREEFREE, "output_root": tmp_path},
        config,
    )

    dust = fits.getdata(DUST).astype(float)
    freefree = np.nan_to_num(fits.getdata(FREEFREE).astype(float))
    model, header = fits.getdata(written[223.8]["continuum"], header=True)
    assert written[223.8]["continuum"].name == "continuum_224GHz.fits"
    assert header["FREQGHZ"] == 223.8
    assert np.array_equal(
        model, (dust + freefree).astype(np.float32), equal_nan=True
    )

    model_250 = fits.getdata(written[250.0]["continuum"]).astype(float)
    expected_250 = (
        dust * (250.0 / 223.8) ** 3.8
        + freefree * (250.0 / 223.8) ** -0.1
    )
    np.testing.assert_allclose(
        model_250, expected_250, rtol=1e-6, equal_nan=True
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
            "tasks/synthesize_continuum/run.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Model readiness: ready" in completed.stdout
    assert "dust=3.8, free-free=-0.1" in completed.stdout
    assert "Configured frequencies [GHz]: 223.8, 250.0" in completed.stdout
    assert "ff224_map_smoothed_18arcsec_regrid2dust.fits" in completed.stdout


def test_plot_stage_writes_pdf_and_png_products(tmp_path: Path) -> None:
    config = SYNTHESIS_SCRIPT.load_task_config(
        PROJECT_ROOT / "tasks" / "synthesize_continuum" / "task.toml"
    )
    paths = {
        "dust": DUST,
        "freefree": FREEFREE,
        "simba_250ghz": SIMBA,
        "output_root": tmp_path,
    }
    SYNTHESIS_SCRIPT.build_models(paths, config)
    SYNTHESIS_SCRIPT.build_simba_validation(paths, config)
    outputs = SYNTHESIS_SCRIPT.plot_products(paths, config)

    assert len(outputs) == 6
    assert {path.suffix for path in outputs} == {".pdf", ".png"}
    assert {path.stem for path in outputs} == {
        "continuum_components_224GHz",
        "simba_comparison_250GHz",
        "sed_demonstration_224GHz",
    }
    for path in outputs:
        assert path.is_file()
        assert path.stat().st_size > 10_000


def test_dust_opacity_power_law_gives_gas_opacity() -> None:
    kappa = dust_opacity_power_law(
        223.8,
        reference_opacity_cm2_per_g=10.0,
        reference_frequency_ghz=1000.0,
        beta=2.0,
        gas_to_dust=100.0,
    )

    assert kappa == pytest.approx(0.1 * 0.2238**2, rel=1e-12)


def test_dust_model_jybeam_uses_beam_solid_angle() -> None:
    column = np.array([1.0e22, 5.0e22])
    temperature = np.array([15.0, 20.0])
    beam_sr = np.pi / (4 * np.log(2)) * np.deg2rad(18.0 / 3600.0) ** 2

    model = dust_model_jybeam(
        column,
        temperature,
        frequency_ghz=223.8,
        gas_opacity_cm2_per_g=0.005,
        mean_molecular_weight=2.8,
        beam_fwhm_arcsec=18.0,
    )

    expected = modified_blackbody_mjysr(
        column, temperature, frequency_ghz=223.8, opacity_cm2_per_g=0.005
    ) * 1.0e6 * beam_sr
    np.testing.assert_allclose(model, expected, rtol=1e-12)


def test_dust_model_jybeam_blanks_nonpositive_temperature() -> None:
    model = dust_model_jybeam(
        np.array([1.0e22, 1.0e22]),
        np.array([0.0, 15.0]),
        frequency_ghz=223.8,
        gas_opacity_cm2_per_g=0.005,
        mean_molecular_weight=2.8,
        beam_fwhm_arcsec=18.0,
    )

    assert np.isnan(model[0]) and np.isfinite(model[1])


def test_fill_missing_keeps_primary_pixels() -> None:
    primary = np.array([1.0, np.nan, np.nan])
    filler = np.array([9.0, 2.0, np.nan])

    filled = fill_missing(primary, filler)

    np.testing.assert_array_equal(filled, [1.0, 2.0, np.nan])


def test_fill_stage_fills_only_outside_the_dust_model(tmp_path: Path) -> None:
    config = SYNTHESIS_SCRIPT.load_task_config(
        PROJECT_ROOT / "tasks" / "synthesize_continuum" / "task.toml"
    )
    inputs = config["inputs"]
    paths = {
        "dust": DUST,
        "freefree": FREEFREE,
        "fill_column_density": PROJECT_ROOT / inputs["fill_column_density"],
        "fill_temperature": PROJECT_ROOT / inputs["fill_temperature"],
        "output_root": tmp_path,
    }

    written = SYNTHESIS_SCRIPT.build_filled_models(paths, config)

    assert written["dust"].name == "dust_224GHz_fill36.fits"
    assert written["continuum"].name == "continuum_224GHz_fill36.fits"
    assert "freefree" not in written
    dust = np.squeeze(fits.getdata(DUST)).astype(float)
    filled = fits.getdata(written["dust"]).astype(float)
    inside = np.isfinite(dust)
    np.testing.assert_array_equal(filled[inside], dust[inside].astype(np.float32))
    assert np.isfinite(filled).sum() > 4 * inside.sum()
    report = written["report"].read_text()
    assert "overlap_ratio_median" in report
