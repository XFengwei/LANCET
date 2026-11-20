"""Tests for the combine_continuum diagnostic (phase 2) stages and their inputs."""

from importlib.util import module_from_spec, spec_from_file_location
import inspect
from pathlib import Path
import sys

from astropy.io import fits
import numpy as np
import pytest

from lancet import ProjectPaths

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "tasks" / "combine_continuum" / "run.py"
CONFIG = PROJECT_ROOT / "tasks" / "combine_continuum" / "task.toml"


def _task():
    name = "combine_continuum_run"
    if name not in sys.modules:
        spec = spec_from_file_location(name, SCRIPT)
        module = module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


TASK = _task()
MISSING_PRODUCT = r"Run stage '(jcomb|feather)' first: missing .*/products/"


def _config():
    return TASK.load_task_config(CONFIG)


def _paths():
    return TASK.resolve_task_paths(ProjectPaths.discover(PROJECT_ROOT), _config())


def _paths_without_products(tmp_path):
    """Resolved inputs, with outputs and products in an empty temporary folder."""
    paths = TASK.diagnostic_paths(_paths())
    paths["output_root"] = tmp_path / "diagnostics"
    paths["products"] = TASK.product_paths(tmp_path)
    return paths


def _forbid(*_args, **_kwargs):
    raise AssertionError("read before the product check")


def test_diagnostic_paths_use_the_diagnostics_root():
    paths = _paths()
    diagnostic = TASK.diagnostic_paths(paths)
    assert diagnostic["output_root"] == paths["output_root"] / "diagnostics"
    assert paths["output_root"] == diagnostic["output_root"].parent   # the original is not changed
    assert diagnostic["products"] == paths["products"]


def test_diagnostic_stages_are_registered():
    assert set(TASK.DIAGNOSTIC_STAGES) <= set(TASK.STAGES)
    from lancet.combination_tests import require_inputs
    assert TASK.require_outputs is require_inputs


def test_diagnostic_wrapper_passes_diagnostic_paths(tmp_path):
    calls = []

    def stage(paths, config, *, save_figure, overwrite=False):
        calls.append((paths["output_root"], save_figure, overwrite))

    paths = {"output_root": tmp_path, "diagnostics_root": tmp_path / "diagnostics"}
    TASK._diagnostic(stage)(paths, {}, True)
    assert calls[0][0] == tmp_path / "diagnostics" and calls[0][2] is True
    assert callable(calls[0][1])


def test_require_inputs_names_the_stage(tmp_path):
    from lancet.combination_tests import require_inputs
    with pytest.raises(FileNotFoundError, match="Run stage 'jcomb' first: missing"):
        require_inputs([tmp_path / "x.fits"], "jcomb")


def test_require_products_names_the_stage_that_writes_each_product(tmp_path):
    from lancet.combination_tests import require_products
    paths = {"products": TASK.product_paths(tmp_path)}
    paths["products"]["jcomb_continuum"].parent.mkdir(parents=True)
    paths["products"]["jcomb_continuum"].write_text("x")
    require_products(paths, ["jcomb_continuum"])
    with pytest.raises(FileNotFoundError, match="Run stage 'feather' first: missing .*feather_continuum"):
        require_products(paths, ["jcomb_continuum", "feather_continuum"])
    assert set(TASK.STAGES) >= {name.split("_")[0] for name in paths["products"]}


def test_overlap_stops_without_jcomb_products(tmp_path, monkeypatch):
    from lancet import combination_tests
    monkeypatch.setattr(combination_tests, "read_test_rows", _forbid)
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        combination_tests.run_overlap_tests(_paths_without_products(tmp_path), _config(),
                                            save_figure=None, overwrite=False)


def test_full_sampling_stops_without_jcomb_products(tmp_path, monkeypatch):
    from lancet import full_sampling
    monkeypatch.setattr(full_sampling, "read_test_rows", _forbid)
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        full_sampling.build_full_sampling_profiles(_paths_without_products(tmp_path), _config())


def test_aca_forward_test_stops_without_products(tmp_path):
    from lancet.aca_forward_test import run_aca_forward_test
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        run_aca_forward_test(_paths_without_products(tmp_path), _config(), save_figure=None)


def test_amplitude_alignment_stops_without_products(tmp_path):
    from lancet.amplitude_alignment import run_amplitude_alignment
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        run_amplitude_alignment(_paths_without_products(tmp_path), _config(), save_figure=None)


def test_figure2_requires_products(tmp_path):
    from lancet.combination_figure import make_combination_figure
    with pytest.raises(FileNotFoundError, match="Run stage 'jcomb' first: missing .*/products/"):
        make_combination_figure(_paths_without_products(tmp_path), _config(), save_figure=None)


def test_figure_has_no_three_panel_switch():
    from lancet.combination_figure import make_combination_figure
    assert "three_panel" not in inspect.signature(make_combination_figure).parameters


def test_image_curve_limit_is_the_lowest_of_beam_nyquist_and_interpolation():
    from lancet.amplitude_alignment import image_curve_limit, interpolation_limit, response_limit
    pixel_rad = 0.5 * np.pi / 648000
    half_nyquist = 0.5 / (2 * pixel_rad) / 1000                       # 0.5 x 206.26 klambda
    # Half-pixel linear interpolation response |cos(pi u p)| = 0.98 at u = arccos(0.98) / (pi p).
    analytic = np.arccos(0.98) / (np.pi * pixel_rad) / 1000
    assert abs(analytic - 26.306) < 1e-3
    assert abs(interpolation_limit(0.5, 0.02) - analytic) < 1e-9
    assert np.isclose(np.cos(np.pi * interpolation_limit(0.5, 0.02) * 1000 * pixel_rad), 0.98)
    assert response_limit(0.66, 0.1) > 240
    assert abs(image_curve_limit(0.66, 0.1, 0.5, 0.02) - analytic) < 1e-9
    assert abs(image_curve_limit(0.66, 0.1, 0.5, 1.0) - half_nyquist) < 1e-3   # loss 1: Nyquist itself
    assert image_curve_limit(18.1, 0.1, 0.5, 0.02) == response_limit(18.1, 0.1)


def test_pixel_binning_transfer_is_the_box_response():
    from lancet.amplitude_alignment import pixel_binning_transfer
    pixel_rad = 0.5 * np.pi / 648000
    assert pixel_binning_transfer(0.0, 0.5) == 1.0
    assert np.isclose(pixel_binning_transfer(1 / pixel_rad / 1000, 0.5), 0.0, atol=1e-12)


def _gaussian_image(path, *, source_fwhm, beam_fwhm, size=2048, pixel=0.1, centre=(221.25, -59.8)):
    """Restored Jy/beam image of a unit-flux circular Gaussian at the reference pixel.

    CRPIX = 1023 makes the reference pixel the centre of a 5-pixel block, so the 0.5" binned grid
    and the 0.5" sampling patch share pixel centres and the patch interpolation is exact.
    """
    restored = np.hypot(source_fwhm, beam_fwhm)
    offset = (np.arange(size) - 1022) * pixel                          # 0-based index 1022 = CRPIX 1023
    r2 = offset[None, :]**2 + offset[:, None]**2
    data = (beam_fwhm / restored)**2 * np.exp(-4 * np.log(2) * r2 / restored**2)
    header = fits.Header({"BUNIT": "JY/BEAM", "BMAJ": beam_fwhm / 3600, "BMIN": beam_fwhm / 3600, "BPA": 0.0,
                          "CTYPE1": "RA---SIN", "CTYPE2": "DEC--SIN", "CRVAL1": centre[0], "CRVAL2": centre[1],
                          "CDELT1": -pixel / 3600, "CDELT2": pixel / 3600, "CRPIX1": 1023.0, "CRPIX2": 1023.0})
    fits.PrimaryHDU(data.astype(np.float32), header).writeto(path)
    return restored


IMAGE_CURVE_TASK = {"tp_minimum_response": 0.1, "continuous_min_klambda": 1.0, "image_curve_points": 60,
                    "image_curve_max_interpolation_loss": 0.02}


def _gaussian_curve(tmp_path, shift_arcsec):
    """Corrected curve of a unit-flux 2" Gaussian (beam 0.7"), the pointing shifted in RA and Dec.

    One pointing pair (0), one band (216.74 GHz), one orientation sector; the PB is flat over the
    source. Returns the grid, the corrected curve and the analytic truth: continuous_rms returns
    mJy of the image times the Gaussian PB, so a unit-flux source keeps the PB-weighted flux
    pb^2 / (pb^2 + restored^2).
    """
    from lancet.amplitude_alignment import _image_curve, gaussian_transfer
    path = tmp_path / "gaussian.fits"
    restored = _gaussian_image(path, source_fwhm=2.0, beam_fwhm=0.7)
    settings = {**_config()["overlap_tests"], "pb_12m_fwhm_arcsec": 2000.0}
    ra, dec = 221.25, -59.8
    centre = {0: (ra + shift_arcsec / 3600 / np.cos(np.deg2rad(dec)), dec + shift_arcsec / 3600)}
    grid, rms, info = _image_curve(path, centre, np.array([0]), 1, 8, settings, IMAGE_CURVE_TASK)
    assert np.isclose(info["beam_fwhm_arcsec"], 0.7)
    pb = 2000.0 * float(settings["pb_reference_frequency_ghz"]) / 216.74
    return grid, rms, 1000 * pb**2 / (pb**2 + restored**2) * gaussian_transfer(grid, 2.0)


def test_beam_and_binning_corrected_curve_recovers_a_gaussian(tmp_path):
    from lancet.amplitude_alignment import gaussian_transfer, image_curve_limit, pixel_binning_transfer
    grid, rms, expected = _gaussian_curve(tmp_path, 0.0)          # pointing on a binned pixel centre
    limit = image_curve_limit(0.7, 0.1, 0.5, 0.02)
    assert np.isclose(grid.max(), limit)
    assert np.all(np.abs(rms / expected - 1) < 0.02), np.max(np.abs(rms / expected - 1))
    # Both corrections are applied: at the limit the beam is ~2.8 %, the binning ~0.7 %.
    assert rms[-1] * gaussian_transfer(limit, 0.7) / expected[-1] < 0.98
    assert rms[-1] * pixel_binning_transfer(limit, 0.5) / expected[-1] < 0.995


def test_image_curve_stops_before_half_pixel_interpolation_loss(tmp_path):
    # Pointing half a binned pixel (0.25") off the grid on both axes: the patch interpolation
    # (|cos(pi u p)| per axis) is not corrected, so the curve is ~2 % low at the limit and no lower.
    grid, rms, expected = _gaussian_curve(tmp_path, 0.25)
    error = rms / expected - 1
    assert np.all(np.abs(error) < 0.02), np.max(np.abs(error))
    assert error[-1] < -0.015


def test_read_products_returns_image_curves(tmp_path):
    from lancet.amplitude_alignment import read_products
    (tmp_path / "array_profiles.csv").write_text("array,uv_klambda,rms_mjy\naca_7m,5.0,100.0\n")
    (tmp_path / "overlap.csv").write_text("uv_klambda,aca_native_mjy\n5.0,90.0\n")
    (tmp_path / "tp_continuous.csv").write_text("uv_klambda,tp_raw_common_mjy\n1.0,1000.0\n2.0,800.0\n")
    (tmp_path / "manifest.json").write_text('{"sky_model": "x"}\n')
    (tmp_path / "image_curves.csv").write_text(
        "name,uv_klambda,rms_mjy\ninterferometer,1.0,50.0\ninterferometer,2.0,60.0\njcomb,1.0,900.0\n")
    profiles, overlap, continuous, manifest, curves = read_products(tmp_path)
    assert manifest == {"sky_model": "x"} and list(curves) == ["interferometer", "jcomb"]
    assert np.array_equal(curves["interferometer"][0], [1.0, 2.0])
    assert np.array_equal(curves["interferometer"][1], [50.0, 60.0])
    assert np.array_equal(curves["jcomb"][1], [900.0])


def test_figure2_names_feather_when_only_jcomb_exists(tmp_path):
    from lancet.combination_figure import make_combination_figure
    paths = _paths_without_products(tmp_path)
    paths["products"]["jcomb_continuum"].parent.mkdir(parents=True)
    paths["products"]["jcomb_continuum"].write_text("x")
    with pytest.raises(FileNotFoundError, match="Run stage 'feather' first"):
        make_combination_figure(paths, _config(), save_figure=None)


def test_figure2_refuses_rerun_before_computing(tmp_path, monkeypatch):
    from lancet import full_sampling
    from lancet.combination_figure import make_combination_figure
    monkeypatch.setattr(full_sampling, "build_full_sampling_profiles", _forbid)
    paths = _paths_without_products(tmp_path)
    for product in paths["products"].values():
        product.parent.mkdir(parents=True, exist_ok=True)
        product.write_text("x")
    config = _config()
    figures = paths["output_root"] / "figures"
    figures.mkdir(parents=True)
    for suffix in config["figure2"]["output_formats"]:
        (figures / f"figure2_combination_diagnostics.{suffix}").write_text("x")
    with pytest.raises(FileExistsError):
        make_combination_figure(paths, config, save_figure=None)


def test_binned_sky_model_accepts_4d_jybeam(tmp_path):
    from lancet.combination_prep import jybeam_to_mjysr
    from lancet.combination_tests import binned_sky_model

    header = fits.Header({"BUNIT": "JY/BEAM", "BMAJ": 1/3600, "BMIN": 1/3600, "BPA": 0.0,
                          "CTYPE1": "RA---SIN", "CTYPE2": "DEC--SIN", "CRVAL1": 0.0, "CRVAL2": 0.0,
                          "CDELT1": -0.1/3600, "CDELT2": 0.1/3600, "CRPIX1": 10.5, "CRPIX2": 10.5,
                          "CTYPE3": "FREQ", "CRVAL3": 2.24e11, "CDELT3": 1e9, "CRPIX3": 1.0,
                          "CTYPE4": "STOKES", "CRVAL4": 1.0, "CDELT4": 1.0, "CRPIX4": 1.0})
    path = tmp_path / "image.fits"
    fits.PrimaryHDU(np.ones((1, 1, 20, 20), np.float32), header).writeto(path)
    data, out = binned_sky_model(path, 0.5)
    assert data.shape == (4, 4)
    assert np.allclose(data, jybeam_to_mjysr(1, 1, 1), rtol=1e-6, atol=0)
    assert np.isclose(abs(out["CDELT2"]) * 3600, 0.5)


def test_joint_and_dirty_beam_take_paths():
    from lancet.dirty_beam import run_dirty_beam
    from lancet.joint_image_comparison import run_joint_image_comparison
    for function in (run_joint_image_comparison, run_dirty_beam):
        assert next(iter(inspect.signature(function).parameters)) == "paths"


def test_config_has_no_path_keys_in_stage_tables():
    config = _config()
    removed = {"casa", "miriad", "casa_pb", "psf_fits", "imfit_log", "paper_image", "gain", "pb",
               "frequency_ghz"}
    for table in ("joint_image_comparison", "dirty_beam", "pointing_coverage"):
        assert not removed & set(config[table]), table


def test_sky_model_sources_resolve_to_path_groups():
    config, paths = _config(), _paths()
    sources = [spec["source"] for spec in config["aca_forward_test"]["sky_models"].values()]
    sources.append(config["amplitude_alignment"]["sky_model"])
    curves = config["amplitude_alignment"]["image_curves"]
    assert curves == {"interferometer": "inputs.interferometer", "jcomb": "products.jcomb_continuum",
                      "feather": "products.feather_continuum"}
    sources.extend(curves.values())
    for source in sources:
        group, name = source.split(".")
        assert name in paths[group], source
    models = config["aca_forward_test"]["sky_models"]
    assert set(models) == {"interferometer", "jcomb", "feather", "aca_only"}
    assert models["jcomb"]["source"] == "products.jcomb_continuum"
    assert not any("tp_scale" in spec for spec in models.values())     # no TP scale is derived
    assert config["amplitude_alignment"]["sky_model"] == "inputs.interferometer"
