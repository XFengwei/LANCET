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


def test_figure2_stops_without_products(tmp_path):
    from lancet.combination_figure import make_combination_figure
    with pytest.raises(FileNotFoundError, match=MISSING_PRODUCT):
        make_combination_figure(_paths_without_products(tmp_path), _config(), save_figure=None)


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
        for stem in ("Figure1_combination_diagnostics", "figure2_combination_diagnostics"):
            (figures / f"{stem}.{suffix}").write_text("x")
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
    for source in sources:
        group, name = source.split(".")
        assert name in paths[group], source
    models = config["aca_forward_test"]["sky_models"]
    assert set(models) == {"interferometer", "jcomb", "feather", "aca_only"}
    assert models["jcomb"]["source"] == "products.jcomb_continuum"
    assert config["amplitude_alignment"]["sky_model"] == "inputs.interferometer"
