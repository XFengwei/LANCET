import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from lancet import ProjectPaths
from lancet.paths import DataOutputError

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "linear_combine_continuum.py"
SPEC = spec_from_file_location("linear_combine_continuum_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "linear_combine_continuum" / "task.toml"


def _config(changes):
    config = TASK.load_task_config(CONFIG)
    for (table, key), value in changes.items():
        config[table][key] = value
    return config


def test_output_root_inside_data_rejected():
    config = _config({("outputs", "root"): "data/delivered/x"})
    with pytest.raises(DataOutputError):
        TASK.resolve_task_paths(ProjectPaths.discover(PROJECT_ROOT), config)


def test_tp_inputs_must_be_delivered():
    config = _config({("inputs", "tp_dust"): "results/synthesize_continuum/model/dust_230GHz.fits"})
    with pytest.raises(ValueError, match="tp_dust"):
        TASK.resolve_task_paths(ProjectPaths.discover(PROJECT_ROOT), config)


def test_guard_outputs_refuses_existing(tmp_path):
    existing = tmp_path / "product.fits"
    existing.write_text("x")
    with pytest.raises(FileExistsError, match="product.fits"):
        TASK.guard_outputs([existing, tmp_path / "new.fits"], overwrite=False)
    TASK.guard_outputs([existing], overwrite=True)


def _header(n, pixel_arcsec, beam_arcsec):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "RA---SIN", "DEC--SIN"
    h["CRVAL1"], h["CRVAL2"] = 221.25, -59.8
    h["CRPIX1"] = h["CRPIX2"] = (n + 1) / 2
    h["CDELT1"], h["CDELT2"] = -pixel_arcsec / 3600, pixel_arcsec / 3600
    h["BUNIT"], h["BMAJ"], h["BMIN"], h["BPA"] = "Jy/beam", beam_arcsec / 3600, beam_arcsec / 3600, 0.0
    return h


def test_prepare_fails_on_small_margin(tmp_path):
    clean = np.full((60, 60), np.nan, dtype=np.float32)
    clean[10:50, 10:50] = 1.0                     # 40" mosaic at 1"/px
    paths = {}
    for name, data in (("interferometer", clean), ("interferometer_attenuated", 0.8 * clean)):
        paths[name] = tmp_path / f"{name}.fits"
        fits.PrimaryHDU(data, _header(60, 1.0, 1.0)).writeto(paths[name])
    tp = np.full((21, 21), np.nan, dtype=np.float32)
    tp[6:15, 6:15] = 1.0                          # 54" footprint: ~5" beyond the mosaic
    for name in ("tp_dust", "tp_dust_ff"):
        paths[name] = tmp_path / f"{name}.fits"
        fits.PrimaryHDU(tp, _header(21, 6.0, 18.1)).writeto(paths[name])
    paths["output_root"] = tmp_path / "out"
    with pytest.raises(RuntimeError, match="margin"):
        TASK.run_prepare(paths, _config({("prepare", "margin_decimation"): 1}))
    assert not (tmp_path / "out").exists()      # nothing is written when an input check fails


def test_product_blank_fills_tp_only_region():
    gain = np.zeros((6, 8), dtype=np.float32)
    gain[2:4, 2:5] = 0.8                      # mosaic
    footprint = np.ones((6, 8), dtype=bool)
    footprint[:, 7] = False                   # TP model missing in the last column
    filled = TASK.product_blank(gain, footprint, fill_outside_mosaic=True)
    np.testing.assert_array_equal(filled, ~footprint)
    cut = TASK.product_blank(gain, footprint, fill_outside_mosaic=False)
    np.testing.assert_array_equal(cut, ~footprint | (gain == 0))
    assert TASK.load_task_config(CONFIG)["combine"]["fill_outside_mosaic"] is True


def test_tp_footprint_on_grid_follows_native_nans():
    from astropy.wcs import WCS

    def wcs(pix_deg, n):
        w = WCS(naxis=2)
        w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
        w.wcs.crval = [221.0, -59.8]
        w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
        w.wcs.cdelt = [-pix_deg, pix_deg]
        return w

    tp = np.ones((21, 21), dtype=np.float32)
    tp[:, :10] = np.nan                       # east half of the TP missing
    footprint = TASK.tp_footprint_on_grid(tp, wcs(6 / 3600, 21), wcs(1 / 3600, 60), (60, 60))
    assert footprint.dtype == bool
    assert footprint[:, -5:].all() and not footprint[:, :5].any()


def test_interferometer_outside_defined_gain_is_not_used():
    clean = np.full((4, 5), 2.0, dtype=np.float32)
    gain = np.zeros((4, 5), dtype=np.float32)
    gain[1:3, 1:4] = 0.6
    used = TASK.interferometer_within_gain(clean, gain)
    np.testing.assert_array_equal(used, np.where(gain > 0, 2.0, 0.0))
    assert used.dtype == np.float32


def test_append_manifest_keeps_earlier_entries(tmp_path):
    import json

    TASK.append_manifest(tmp_path, {"stage": "prepare", "outputs": ["a.fits"]})
    TASK.append_manifest(tmp_path, {"stage": "combine", "outputs": [tmp_path / "b.fits"]})
    entries = json.loads((tmp_path / "manifest.json").read_text())
    assert [e["stage"] for e in entries] == ["prepare", "combine"]
    assert all("utc" in e for e in entries)
    assert entries[1]["outputs"] == [str(tmp_path / "b.fits")]


def test_ratio_flag_and_output_paths(tmp_path):
    assert TASK.ratio_check(1.04, 0.05) == {"value": 1.04, "tolerance": 0.05, "flag": False}
    assert TASK.ratio_check(0.78, 0.10)["flag"] is True
    found = TASK.output_paths({"a": tmp_path / "x.fits", "b": 3, "c": {"d": tmp_path / "y.csv"}})
    assert found == [str(tmp_path / "x.fits"), str(tmp_path / "y.csv")]
    assert TASK.output_paths(tmp_path / "z.json") == [str(tmp_path / "z.json")]


def test_product_metadata_records_padding():
    config = TASK.load_task_config(CONFIG)
    settings = TASK.combine_settings(config, "jcomb", 1.4)
    meta = TASK.product_metadata(config, "jcomb", "dust_ff", settings)
    assert meta["padding_arcsec"] == config["prepare"]["padding_arcsec"]
    assert "reflect" in meta["padding_mode"]
    assert meta["k1"] == config["jcomb"]["k1_klambda"] and meta["fill_outside_mosaic"] is True
    assert meta["tp_factor"] == 1.4 and meta["absolute_scale"] == "pending"
    assert TASK.combine_settings(config, "feather", 1)["tp_factor"] == 1.0
    assert config["tp_factor"]["value"] == 1.0 and config["tp_factor"]["determined_from"] == "dust_ff"


def test_stages_and_prepared_layout(tmp_path):
    assert list(TASK.STAGES) == ["prepare", "tp_factor", "jcomb", "feather", "box_spectrum", "planck_calibration"]
    names = sorted(p.name for p in TASK.prepared_paths(tmp_path).values())
    assert names == ["gain.fits", "input_check.json", "interferometer_mjysr.fits", "tp_dust_ff_footprint.fits",
                     "tp_dust_ff_mjysr.fits", "tp_dust_footprint.fits", "tp_dust_mjysr.fits"]
    with pytest.raises(FileNotFoundError, match="prepare"):
        TASK.load_prepared({"output_root": tmp_path}, "dust")


def test_catalogue_positions_skip_pixels_rounding_off_the_grid():
    gain = np.ones((10, 12), dtype=np.float32)
    gain[:, :3] = 0.2
    pos = TASK.catalogue_positions(np.array([1.0, 6.0, 11.6, -0.6]), np.array([5.0, 5.0, 5.0, 5.0]), gain, 0.5)
    assert pos == [(6.0, 5.0)]
