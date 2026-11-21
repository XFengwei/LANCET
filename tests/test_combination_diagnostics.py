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


@pytest.mark.local_inputs("data", "results")
def test_diagnostic_paths_use_the_diagnostics_root():
    paths = _paths()
    diagnostic = TASK.diagnostic_paths(paths)
    assert diagnostic["output_root"] == paths["output_root"] / "diagnostics"
    assert paths["output_root"] == diagnostic["output_root"].parent   # the original is not changed
    assert diagnostic["products"] == paths["products"]


def test_diagnostic_stages_are_registered():
    assert set(TASK.DIAGNOSTIC_STAGES) <= set(TASK.STAGES)
    from lancet.combination.overlap_tests import require_inputs
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
    from lancet.combination.overlap_tests import require_inputs
    with pytest.raises(FileNotFoundError, match="Run stage 'jcomb' first: missing"):
        require_inputs([tmp_path / "x.fits"], "jcomb")


def test_require_products_names_the_stage_that_writes_each_product(tmp_path):
    from lancet.combination.overlap_tests import require_products
    paths = {"products": TASK.product_paths(tmp_path)}
    paths["products"]["jcomb_continuum"].parent.mkdir(parents=True)
    paths["products"]["jcomb_continuum"].write_text("x")
    require_products(paths, ["jcomb_continuum"])
    with pytest.raises(FileNotFoundError, match="Run stage 'feather' first: missing .*feather_continuum"):
        require_products(paths, ["jcomb_continuum", "feather_continuum"])
    assert set(TASK.STAGES) >= {name.split("_")[0] for name in paths["products"]}


@pytest.mark.local_inputs("data", "results")
def test_overlap_stops_without_jcomb_products(tmp_path, monkeypatch):
    from lancet.combination import overlap_tests as combination_tests
    monkeypatch.setattr(combination_tests, "read_test_rows", _forbid)
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        combination_tests.run_overlap_tests(_paths_without_products(tmp_path), _config(),
                                            save_figure=None, overwrite=False)


@pytest.mark.local_inputs("data", "results")
def test_full_sampling_stops_without_jcomb_products(tmp_path, monkeypatch):
    from lancet.combination import full_sampling
    monkeypatch.setattr(full_sampling, "read_test_rows", _forbid)
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        full_sampling.build_full_sampling_profiles(_paths_without_products(tmp_path), _config())


@pytest.mark.local_inputs("data", "results")
def test_aca_forward_test_stops_without_products(tmp_path):
    from lancet.combination.aca_forward_test import run_aca_forward_test
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        run_aca_forward_test(_paths_without_products(tmp_path), _config(), save_figure=None)


@pytest.mark.local_inputs("data", "results")
def test_amplitude_alignment_stops_without_products(tmp_path):
    from lancet.combination.amplitude_alignment import run_amplitude_alignment
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        run_amplitude_alignment(_paths_without_products(tmp_path), _config(), save_figure=None)


@pytest.mark.local_inputs("data", "results")
def test_figure2_requires_products(tmp_path):
    from lancet.combination.coverage_figure import make_combination_figure
    with pytest.raises(FileNotFoundError, match="Run stage 'jcomb' first: missing .*/products/"):
        make_combination_figure(_paths_without_products(tmp_path), _config(), save_figure=None)


def test_figure_has_no_three_panel_switch():
    from lancet.combination.coverage_figure import make_combination_figure
    assert "three_panel" not in inspect.signature(make_combination_figure).parameters


def test_image_curve_limit_is_the_lowest_of_beam_nyquist_and_interpolation():
    from lancet.combination.amplitude_alignment import image_curve_limit, interpolation_limit, response_limit
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
    from lancet.combination.amplitude_alignment import pixel_binning_transfer
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
    from lancet.combination.amplitude_alignment import _image_curve, gaussian_transfer
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
    from lancet.combination.amplitude_alignment import gaussian_transfer, image_curve_limit, pixel_binning_transfer
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
    from lancet.combination.amplitude_alignment import read_products
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


@pytest.mark.local_inputs("data", "results")
def test_figure2_names_feather_when_only_jcomb_exists(tmp_path):
    from lancet.combination.coverage_figure import make_combination_figure
    paths = _paths_without_products(tmp_path)
    paths["products"]["jcomb_continuum"].parent.mkdir(parents=True)
    paths["products"]["jcomb_continuum"].write_text("x")
    with pytest.raises(FileNotFoundError, match="Run stage 'feather' first"):
        make_combination_figure(paths, _config(), save_figure=None)


@pytest.mark.local_inputs("data", "results")
def test_figure2_refuses_rerun_before_computing(tmp_path, monkeypatch):
    from lancet.combination import full_sampling
    from lancet.combination.coverage_figure import make_combination_figure
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
    from lancet.combination.prep import jybeam_to_mjysr
    from lancet.combination.overlap_tests import binned_sky_model

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
    from lancet.combination.dirty_beam import run_dirty_beam
    from lancet.combination.joint_image_comparison import run_joint_image_comparison
    for function in (run_joint_image_comparison, run_dirty_beam):
        assert next(iter(inspect.signature(function).parameters)) == "paths"


def test_config_has_no_path_keys_in_stage_tables():
    config = _config()
    removed = {"casa", "miriad", "casa_pb", "psf_fits", "imfit_log", "paper_image", "gain", "pb",
               "frequency_ghz"}
    for table in ("joint_image_comparison", "dirty_beam", "pointing_coverage"):
        assert not removed & set(config[table]), table


@pytest.mark.local_inputs("data", "results")
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


# ---------------------------------------------------------------- stage dependencies (spec section 2)

OVERLAP_FILES = ("summary.json", "aca_7m_rows.npz", "aca_7m_metadata.json", "configuration_c5_metadata.json")
FULL_SAMPLING_FILES = ("summary.json", "amplitudes.csv", "cell_statistics.npz", "configuration_c2_rows.npz",
                       "configuration_c5_rows.npz")


def _paths_with(tmp_path, products=None, **stage_files):
    """Diagnostic paths with placeholder products and placeholder files of earlier stages."""
    paths = _paths_without_products(tmp_path)
    for name in products if products is not None else paths["products"]:
        paths["products"][name].parent.mkdir(parents=True, exist_ok=True)
        paths["products"][name].write_text("x")
    for folder, names in stage_files.items():
        (paths["output_root"] / folder).mkdir(parents=True, exist_ok=True)
        for name in names:
            (paths["output_root"] / folder / name).write_text("x")
    return paths


@pytest.mark.local_inputs("data", "results")
def test_full_sampling_names_overlap_before_reading(tmp_path, monkeypatch):
    from lancet.combination import full_sampling
    monkeypatch.setattr(full_sampling, "read_test_rows", _forbid)
    paths = _paths_with(tmp_path)
    with pytest.raises(FileNotFoundError, match="Run stage 'overlap' first"):
        full_sampling.build_full_sampling_profiles(paths, _config())
    assert not (paths["output_root"] / "full_sampling").exists()


@pytest.mark.local_inputs("data", "results")
def test_aca_forward_test_names_overlap_and_full_sampling(tmp_path):
    from lancet.combination.aca_forward_test import run_aca_forward_test
    with pytest.raises(FileNotFoundError, match="Run stage 'overlap' first"):
        run_aca_forward_test(_paths_with(tmp_path / "a"), _config(), save_figure=None)
    paths = _paths_with(tmp_path / "b", overlap_tests=OVERLAP_FILES)
    with pytest.raises(FileNotFoundError, match="Run stage 'full_sampling' first"):
        run_aca_forward_test(paths, _config(), save_figure=None)
    assert not (paths["output_root"] / "aca_forward_test").exists()


@pytest.mark.local_inputs("data", "results")
def test_amplitude_alignment_names_overlap_and_full_sampling(tmp_path):
    from lancet.combination.amplitude_alignment import run_amplitude_alignment
    with pytest.raises(FileNotFoundError, match="Run stage 'overlap' first"):
        run_amplitude_alignment(_paths_with(tmp_path / "a"), _config(), save_figure=None)
    paths = _paths_with(tmp_path / "b", overlap_tests=OVERLAP_FILES)
    with pytest.raises(FileNotFoundError, match="Run stage 'full_sampling' first"):
        run_amplitude_alignment(paths, _config(), save_figure=None)
    assert not (paths["output_root"] / "amplitude_alignment").exists()


@pytest.mark.local_inputs("data", "results")
def test_amplitude_alignment_requires_only_the_configured_image_curve_products(tmp_path):
    from lancet.combination.amplitude_alignment import run_amplitude_alignment
    config = _config()
    config["amplitude_alignment"]["image_curves"] = {"interferometer": "inputs.interferometer",
                                                     "jcomb": "products.jcomb_continuum"}
    paths = _paths_with(tmp_path, products=["jcomb_continuum"])   # no feather product
    with pytest.raises(FileNotFoundError, match="Run stage 'overlap' first"):
        run_amplitude_alignment(paths, config, save_figure=None)


@pytest.mark.local_inputs("data", "results")
def test_figure2_names_each_missing_stage(tmp_path, monkeypatch):
    from lancet.combination import full_sampling
    from lancet.combination.coverage_figure import make_combination_figure
    monkeypatch.setattr(full_sampling, "build_full_sampling_profiles", _forbid)
    config = _config()
    assert config["figure2"]["amplitude_source"] == "amplitude_alignment"
    uv = [f"{name}.npz" for name in ("aca_7m", "configuration_c2", "configuration_c5")]
    cases = [({}, "uv"), ({"uv_coordinates": uv}, "overlap"),
             ({"uv_coordinates": uv, "overlap_tests": OVERLAP_FILES}, "full_sampling"),
             ({"uv_coordinates": uv, "overlap_tests": OVERLAP_FILES, "full_sampling": FULL_SAMPLING_FILES},
              "amplitude_alignment")]
    for i, (files, stage) in enumerate(cases):
        paths = _paths_with(tmp_path / str(i), **files)
        with pytest.raises(FileNotFoundError, match=f"Run stage '{stage}' first"):
            make_combination_figure(paths, config, save_figure=None)
        assert not (paths["output_root"] / "figures").exists()


def test_figure2_dependencies_include_the_alignment_products():
    from lancet.combination.amplitude_alignment import ALIGNMENT_OUTPUTS
    from lancet.combination.coverage_figure import figure_dependencies
    stages = figure_dependencies("amplitude_alignment")
    assert list(stages) == ["uv", "overlap", "full_sampling", "amplitude_alignment"]
    assert {"image_curves.csv", "array_profiles.csv", "tp_continuous.csv"} <= set(ALIGNMENT_OUTPUTS)
    assert stages["amplitude_alignment"] == ("amplitude_alignment", ALIGNMENT_OUTPUTS)
    assert "amplitude_alignment" not in figure_dependencies("full_sampling")


ALIGNMENT_RESULTS = PROJECT_ROOT / "results" / "combine_continuum" / "diagnostics" / "amplitude_alignment"


@pytest.mark.skipif(not (ALIGNMENT_RESULTS / "manifest.json").is_file(), reason="needs an amplitude_alignment run")
def test_alignment_figure_has_no_old_tp_factor_reference_line():
    from matplotlib import pyplot as plt
    from lancet.combination.amplitude_alignment import plot_alignment, read_products
    profiles, overlap, continuous, manifest, _ = read_products(ALIGNMENT_RESULTS)   # read only
    figures = []
    plot_alignment(profiles, overlap, continuous, _config()["amplitude_alignment"],
                   manifest["tp_response_limit_klambda"], ALIGNMENT_RESULTS,
                   lambda fig, outputs, **_: figures.append(fig))
    ratio = figures[0].axes[2]
    assert not any("2.7" in t.get_text() for t in ratio.texts)
    assert not any(np.allclose(line.get_ydata(), 1 / 2.7) for line in ratio.lines)
    plt.close("all")


def test_diagnostic_manifests_carry_no_fixed_date():
    import re
    from lancet.combination import aca_forward_test, amplitude_alignment, coverage_figure as combination_figure
    for module in (aca_forward_test, amplitude_alignment, combination_figure):
        assert not re.search(r'"date":\s*"\d{4}-', inspect.getsource(module)), module.__name__


def test_pb_uv_footprint_fwhm():
    from lancet.combination.amplitude_alignment import pb_uv_fwhm_klambda, response_limit
    # Gaussian PB of FWHM theta: uv-plane Gaussian of FWHM 4 ln2 / (pi theta).
    for theta in (26.0 * 225.0 / 216.74, 26.0 * 225.0 / 230.95):
        assert np.isclose(pb_uv_fwhm_klambda(theta), 2 * response_limit(theta, 0.5))
    assert abs(pb_uv_fwhm_klambda(26.0 * 225.0 / 216.74) - 6.744) < 1e-3
    assert abs(pb_uv_fwhm_klambda(26.0 * 225.0 / 230.95) - 7.187) < 1e-3


def test_image_curve_method_states_the_pb_smoothing():
    from lancet.combination.amplitude_alignment import image_curve_method, pb_smoothing_note
    config = _config()
    note = pb_smoothing_note(config["overlap_tests"])
    for value in ("26.99", "25.33", "6.74", "7.19", "4 ln2 / (pi theta)"):
        assert value in note, value
    assert note in image_curve_method(config["amplitude_alignment"], config["overlap_tests"])


def test_image_curves_follow_configured_set_and_labels():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from lancet.combination.coverage_figure import _plot_image_curves
    u = np.array([1.0, 2.0, 3.0])
    curves = {name: (u, u * k) for name, k in [("interferometer", 1), ("jcomb", 2), ("feather", 3)]}
    figure, axis = plt.subplots()
    _plot_image_curves(axis, curves, ["interferometer", "jcomb"],
                       {"interferometer": "12m+7m", "jcomb": "combine"})
    assert [line.get_label() for line in axis.get_lines()] == ["12m+7m", "combine"]
    assert [t.get_text() for t in axis.get_legend().get_texts()] == ["12m+7m", "combine"]
    assert not any(np.allclose(line.get_ydata(), 3 * u) for line in axis.get_lines())
    plt.close(figure)


def test_image_curve_settings_are_validated():
    from lancet.combination.coverage_figure import check_image_curve_settings
    have = {"interferometer", "jcomb"}
    check_image_curve_settings(["interferometer", "jcomb"], {"jcomb": "combine"}, have)
    with pytest.raises(ValueError, match="bogus"):
        check_image_curve_settings(["bogus"], {}, have)
    with pytest.raises(ValueError, match="feather"):
        check_image_curve_settings(["feather"], {}, have)
    with pytest.raises(ValueError, match="stray"):
        check_image_curve_settings(["jcomb"], {"stray": "x"}, have)


def test_shown_curve_values_need_no_curves_without_alignment():
    from lancet.combination.coverage_figure import shown_curve_values
    assert shown_curve_values({}, ["interferometer", "jcomb"]) == []
    curves = {"jcomb": (np.array([1.0]), np.array([2.0, 3.0])), "feather": (np.array([1.0]), np.array([9.0]))}
    assert shown_curve_values(curves, ["jcomb"]) == [2.0, 3.0]


def test_dark_theme_colours_every_array_and_curve_of_the_light_theme():
    from lancet.combination.coverage_figure import THEMES, image_curve_styles

    light, dark = THEMES["light"], THEMES["dark"]
    assert set(light) == set(dark)
    assert set(dark["arrays"]) == set(light["arrays"])
    assert set(dark["curves"]) == set(light["curves"])
    styles = image_curve_styles(["jcomb"], {}, dark["curves"])
    assert styles["jcomb"][2] == dark["curves"]["jcomb"]
    assert image_curve_styles(["jcomb"], {})["jcomb"][2] == light["curves"]["jcomb"]
