"""Tests for the combine_continuum phase 1 stages (src/lancet/combination_stages.py)."""

import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

from astropy.io import fits
import numpy as np
import pytest
from scipy import ndimage

from lancet import combination_stages as stages

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUN = PROJECT_ROOT / "tasks" / "combine_continuum" / "run.py"
SPEC = spec_from_file_location("combine_continuum_run_stages", RUN)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "combine_continuum" / "task.toml"


def _config(changes):
    config = TASK.load_task_config(CONFIG)
    for (table, key), value in changes.items():
        config[table][key] = value
    return config


def _header(n, pixel_arcsec, beam_arcsec):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "RA---SIN", "DEC--SIN"
    h["CRVAL1"], h["CRVAL2"] = 221.25, -59.8
    h["CRPIX1"] = h["CRPIX2"] = (n + 1) / 2
    h["CDELT1"], h["CDELT2"] = -pixel_arcsec / 3600, pixel_arcsec / 3600
    h["BUNIT"], h["BMAJ"], h["BMIN"], h["BPA"] = "Jy/beam", beam_arcsec / 3600, beam_arcsec / 3600, 0.0
    return h


def _truesky_header(n, pixel_arcsec, beam_arcsec):
    """4-D header (RA, Dec, FREQ, STOKES) as in the delivered MIRIAD true-sky image."""
    h = _header(n, pixel_arcsec, beam_arcsec)
    h["BUNIT"] = "JY/BEAM"
    h["CTYPE3"], h["CRVAL3"], h["CRPIX3"], h["CDELT3"] = "FREQ", 2.24e11, 1.0, 1.0e9
    h["CTYPE4"], h["CRVAL4"], h["CRPIX4"], h["CDELT4"] = "STOKES", 1.0, 1.0, 1.0
    return h


def _gain_header(n):
    """Header of the delivered gain map: 4-D, BUNIT 'GAIN', no beam."""
    h = _truesky_header(n, 0.5, 1.0)
    h["BUNIT"] = "GAIN"
    for key in ("BMAJ", "BMIN", "BPA"):
        del h[key]
    return h


def _smooth_field(n=61, border=5, seed=3):
    """Positive field smoothed to about the TP beam, NaN in a border (as outside a real TP footprint)."""
    field = ndimage.gaussian_filter(np.random.default_rng(seed).normal(size=(n, n)), 1.3, mode="wrap")
    field = (1.0 + field / np.abs(field).max()).astype(np.float32)
    field[:border], field[-border:], field[:, :border], field[:, -border:] = np.nan, np.nan, np.nan, np.nan
    return field


def _prepare_config():
    # The synthetic 61 x 61 TP is far coarser than a real model, so the regrid check gets a looser bound.
    return _config({("prepare", "margin_decimation"): 1, ("prepare", "regrid_power_tolerance"): 0.1})


def _inputs(tmp_path, tp_frequency=223.8, truesky_edge=100):
    """Synthetic inputs: 400 x 400 px at 0.5" (200" field), TP 61 x 61 px at 6"."""
    n = 400
    truesky = np.full((n, n), np.nan, dtype=np.float32)
    truesky[truesky_edge:n - truesky_edge, truesky_edge:n - truesky_edge] = 1.0
    paths = {"interferometer": tmp_path / "interferometer.fits", "gain": tmp_path / "gain.fits"}
    fits.PrimaryHDU(truesky[None, None], _truesky_header(n, 0.5, 1.0)).writeto(paths["interferometer"])
    fits.PrimaryHDU(np.full((n, n), 0.9, dtype=np.float32)[None, None], _gain_header(n)).writeto(paths["gain"])
    tp = _smooth_field(61)
    for name in ("tp_dust", "tp_dust_ff"):
        header = _header(61, 6.0, 18.1)
        if tp_frequency is not None:
            header["FREQGHZ"] = tp_frequency
        paths[name] = tmp_path / f"{name}.fits"
        fits.PrimaryHDU(tp, header).writeto(paths[name])
    paths["output_root"] = tmp_path / "out"
    return paths


def test_mosaic_gain_is_zero_where_the_interferometer_is_blank():
    interferometer = np.ones((4, 4), dtype=np.float32)
    interferometer[0] = np.nan
    gain = stages.mosaic_gain(np.full((4, 4), 0.8), interferometer)
    assert gain.dtype == np.float32
    np.testing.assert_array_equal(gain[0], 0.0)
    np.testing.assert_allclose(gain[1:], 0.8, rtol=1e-6)


def test_prepare_reads_4d_truesky_and_tp_frequency(tmp_path):
    paths = _inputs(tmp_path)
    out = stages.run_prepare(paths, _prepare_config())
    gain = fits.getdata(out["gain"])
    assert gain[200, 200] == pytest.approx(0.9) and gain[50, 50] == 0 and gain[200, 350] == 0
    np.testing.assert_allclose(gain[100:300, 100:300], 0.9, rtol=1e-6)
    assert gain[:100].max() == 0 and gain[300:].max() == 0
    for name in ("tp_dust", "tp_dust_ff"):
        assert fits.getheader(out[name])["HIERARCH COMB_TP_FREQUENCY_GHZ"] == pytest.approx(223.8)
    interferometer = fits.getdata(out["interferometer"])
    assert fits.getheader(out["interferometer"])["BITPIX"] == -32 and interferometer.ndim == 2
    assert interferometer.shape == (400, 400) and np.isfinite(interferometer).all()
    assert interferometer[200, 200] > 0 and interferometer[0, 0] == 0


def test_prepare_stops_on_tp_without_freqghz(tmp_path):
    paths = _inputs(tmp_path)
    data = fits.getdata(paths["tp_dust"])
    fits.PrimaryHDU(data, _header(61, 6.0, 18.1)).writeto(paths["tp_dust"], overwrite=True)
    with pytest.raises(ValueError, match="FREQGHZ"):
        stages.run_prepare(paths, _prepare_config())
    prepared = paths["output_root"] / "prepared"
    assert not prepared.exists() or not any(prepared.iterdir())    # no partial outputs block a rerun


def test_tp_frequency_requires_the_same_freqghz_in_both_models(tmp_path):
    paths = _inputs(tmp_path)
    assert stages.tp_frequency_ghz(paths) == pytest.approx(223.8)
    data, header = fits.getdata(paths["tp_dust"], header=True)
    header["FREQGHZ"] = 250.0
    fits.PrimaryHDU(data, header).writeto(paths["tp_dust"], overwrite=True)
    with pytest.raises(ValueError, match="tp_dust .*tp_dust_ff"):
        stages.tp_frequency_ghz(paths)


def test_tp_factor_names_prepare_when_prepared_inputs_are_missing(tmp_path, monkeypatch):
    paths = _inputs(tmp_path)
    monkeypatch.setattr(stages, "read_plane", _forbid)
    with pytest.raises(FileNotFoundError, match="Run stage 'prepare' first"):
        stages.run_tp_factor(paths, _prepare_config())


def _forbid(*_args, **_kwargs):
    raise AssertionError("input read before the prepared outputs were checked")


def test_prepare_fails_on_small_margin(tmp_path):
    paths = _inputs(tmp_path, truesky_edge=0)       # 200" mosaic finite to the edge of the grid
    for name in ("tp_dust", "tp_dust_ff"):          # TP footprint 198" wide: it ends inside the mosaic
        header = _header(61, 6.0, 18.1)
        header["FREQGHZ"] = 223.8
        fits.PrimaryHDU(_smooth_field(border=14), header).writeto(paths[name], overwrite=True)
    with pytest.raises(RuntimeError, match="margin"):
        stages.run_prepare(paths, _prepare_config())
    assert not (tmp_path / "out").exists()          # nothing is written when an input check fails


# ---------------------------------------------------------------- phase 1 checks and products

def test_product_metadata_takes_tp_frequency():
    meta = stages.product_metadata(_config({}), "jcomb", "dust_ff", {}, 223.8)
    assert meta["tp_frequency_ghz"] == 223.8


def test_product_blank_outside_mosaic_when_not_filling():
    blank = stages.product_blank(np.array([[0.0, 0.5]]), np.ones((1, 2), dtype=bool), False)
    np.testing.assert_array_equal(blank, [[True, False]])


def test_tp_factor_has_only_the_prepared_variant():
    assert set(stages.TP_FACTOR_VARIANTS) == {"prepared"}


def test_tp_factor_stage_runs_without_psd_matched(tmp_path):
    paths = _inputs(tmp_path)
    config = _prepare_config()
    config["tp_factor"] = {k: config["tp_factor"][k] for k in (
        "value", "determined_from", "matching_range_arcsec", "window_taper_arcsec", "bins_klambda")}
    stages.run_prepare(paths, config)
    out = stages.run_tp_factor(paths, config)
    summary_path = paths["output_root"] / "checks" / "tp_factor" / "tp_factor.json"
    assert out["summary"] == summary_path and summary_path.is_file()
    text = summary_path.read_text()
    assert "psd_matched" not in text and "as_determined" not in json.loads(text)["dust_ff"]
    assert set(json.loads(text)["dust_ff"]) == {"prepared"}


def test_box_spectrum_requires_the_jcomb_product(tmp_path):
    paths = _method_inputs(tmp_path)
    with pytest.raises(FileNotFoundError, match="jcomb"):
        stages.run_box_spectrum(paths, _config({}))
    assert not (paths["output_root"] / "checks").exists()


def _method_inputs(tmp_path, n=32):
    paths = _inputs(tmp_path)
    paths["products"] = TASK.product_paths(paths["output_root"], "224GHz")
    paths["pb"] = tmp_path / "pb.fits"
    fits.PrimaryHDU(np.full((n, n), 0.9, dtype=np.float32)[None, None], _gain_header(n)).writeto(paths["pb"])
    return paths


def _patch_method_stages(monkeypatch, calls):
    n = 32
    header = _header(n, 1.0, 2.0)
    ctx = {"low": np.ones((n, n), np.float32), "high": np.ones((n, n), np.float32),
           "gain": np.ones((n, n), np.float32), "header": header, "tp_header": header, "pix": 1.0,
           "high_beam": (2.0, 2.0, 0.0), "low_beam": (6.0, 6.0, 0.0), "clean": np.ones((n, n)),
           "tp_native": np.ones((n, n)), "blank": np.zeros((n, n), bool),
           "footprint": np.ones((n, n), bool), "sigma": 0.1, "centres": np.array([1.0, 2.0])}

    def combine(*args, **kwargs):
        calls.append(1)
        return np.ones((n, n))

    check = {"mosaic": {"flux_ratio": 1.0, "n": 4}}
    measured = {"checks": {"large_scale_ratio": {"value": 1.0, "tolerance": 0.1, "flag": False}},
                "power": {"interferometer": [1.0], "scaled_tp": [1.0], "combined": [1.0]},
                "amplitude_ratio_to_interferometer": [1.0], "amplitude_ratio_to_scaled_tp_at_high_beam": [1.0]}
    monkeypatch.setattr(stages, "measurement_context", lambda *a, **k: dict(ctx))
    monkeypatch.setattr(stages, "combine_prepared", combine)
    monkeypatch.setattr(stages, "measure_product", lambda *a, **k: {**measured, "checks": dict(measured["checks"])})
    monkeypatch.setattr(stages, "zero_spacing", lambda *a, **k: (dict(check), np.ones((n, n)), np.ones((n, n)),
                                                                 None, None))
    monkeypatch.setattr(stages, "synthetic_gate", lambda *a, **k: {})
    for name in ("plot_comparison", "plot_zero_spacing", "plot_power_spectrum"):
        monkeypatch.setattr(stages, name, lambda *a, **k: None)


def test_method_output_layout(tmp_path, monkeypatch):
    paths, calls = _method_inputs(tmp_path), []
    _patch_method_stages(monkeypatch, calls)
    config, root = _config({}), paths["output_root"]
    stages.run_method("jcomb", paths, config)
    for kind in ("continuum", "dust"):
        product = root / "products" / f"jcomb_{kind}_224GHz.fits"
        header = fits.getheader(product)
        assert header["BITPIX"] == -32 and header["HIERARCH COMB_TP_FREQUENCY_GHZ"] == 223.8
        assert (root / "checks" / "zero_spacing" / f"jcomb_{kind}_224GHz_at_tp_beam.fits").is_file()
        assert (root / "checks" / "zero_spacing" / f"zero_spacing_ratio_jcomb_{kind}_224GHz.fits").is_file()
    assert (root / "checks" / "jcomb" / "measurements.json").is_file()
    assert not (root / "jcomb").exists()
    done = len(calls)
    with pytest.raises(FileExistsError):
        stages.run_method("jcomb", paths, config)
    assert len(calls) == done


def test_pb_taper():
    pb = np.array([np.nan, 0.1, 0.2, 0.25, 0.3, 0.9], dtype=np.float32)
    weight = stages.pb_taper(pb, 0.2, 0.3)
    assert weight.dtype == np.float32
    np.testing.assert_allclose(weight, [0, 0, 0, 0.5, 1, 1], atol=1e-6)
    for bad in ((0.3, 0.2), (0.0, 0.3), (0.2, 1.1)):
        with pytest.raises(ValueError):
            stages.pb_taper(pb, *bad)


def _extended_setup(monkeypatch, n=128, pix=2.0):
    """Real combination and zero-spacing check on a synthetic mosaic (pb 0.2 at its edge) inside a wider TP."""
    from astropy.wcs import WCS
    yy, xx = np.mgrid[:n, :n]
    radius = np.hypot(yy - n / 2 + 0.5, xx - n / 2 + 0.5)
    pb = np.clip(1.0 - 0.8 * radius / 24.0, 0.0, 1.0).astype(np.float32)   # 0.2 at 24 px
    pb[radius > 24] = 0.1
    mosaic = pb >= 0.2
    gain = np.where(mosaic, np.minimum(1.0, pb / 0.508), 0.0).astype(np.float32)
    high = np.where(mosaic, 50.0 * np.exp(-radius ** 2 / (2 * 5.0 ** 2)), 0.0).astype(np.float32)
    header, tp_header = _header(n, pix, 4.0), _header(n // 4, 4 * pix, 18.1)
    ctx = {"low": np.ones((n, n), np.float32), "high": high, "gain": gain, "header": header,
           "tp_header": tp_header, "tp_wcs": WCS(tp_header), "pix": pix, "high_beam": (4.0, 4.0, 0.0),
           "low_beam": (18.1, 18.1, 0.0), "clean": high, "tp_native": np.ones((n // 4, n // 4)),
           "footprint": np.ones((n, n), bool), "centres": np.array([1.0, 2.0]), "sigma": 1.0e-5,
           "blank": ~mosaic}
    measured = {"checks": {"large_scale_ratio": {"value": 1.0, "tolerance": 0.1, "flag": False}},
                "power": {"interferometer": [1.0], "scaled_tp": [1.0], "combined": [1.0]},
                "amplitude_ratio_to_interferometer": [1.0], "amplitude_ratio_to_scaled_tp_at_high_beam": [1.0]}
    monkeypatch.setattr(stages, "measurement_context", lambda *a, **k: dict(ctx))
    monkeypatch.setattr(stages, "measure_product", lambda *a, **k: {**measured, "checks": dict(measured["checks"])})
    monkeypatch.setattr(stages, "synthetic_gate", lambda *a, **k: {})
    for name in ("plot_comparison", "plot_zero_spacing", "plot_power_spectrum"):
        monkeypatch.setattr(stages, name, lambda *a, **k: None)
    return pb


def test_method_writes_extended_products(tmp_path, monkeypatch):
    n = 128
    pb = _extended_setup(monkeypatch, n)
    paths = _method_inputs(tmp_path, n)
    fits.PrimaryHDU(pb[None, None], _gain_header(n)).writeto(paths["pb"], overwrite=True)
    config, root = _config({}), paths["output_root"]
    stages.run_method("jcomb", paths, config)
    for kind in ("continuum", "dust"):
        standard = fits.getdata(root / "products" / f"jcomb_{kind}_224GHz.fits")
        path = root / "products" / f"jcomb_{kind}_224GHz_extended.fits"
        header, extended = fits.getheader(path), fits.getdata(path)
        assert header["BITPIX"] == -32 and header["HIERARCH COMB_EXTENDED"] is True
        assert header["HIERARCH COMB_TAPER_PB0"] == 0.2 and header["HIERARCH COMB_TAPER_PB1"] == 0.3
        assert header["HIERARCH COMB_FILL_OUTSIDE_MOSAIC"] is True
        assert fits.getheader(root / "products" / f"jcomb_{kind}_224GHz.fits")["HIERARCH COMB_FILL_OUTSIDE_MOSAIC"] is False
        outside = pb < 0.2
        assert np.isfinite(extended[outside]).all() and np.isnan(standard[outside]).all()
        core = (pb >= 0.3) & np.isfinite(standard)
        np.testing.assert_allclose(extended[core], standard[core], rtol=1e-3)
        zero = root / "checks" / "zero_spacing"
        assert (zero / f"jcomb_{kind}_224GHz_extended_at_tp_beam.fits").is_file()
        assert (zero / f"zero_spacing_ratio_jcomb_{kind}_224GHz_extended.fits").is_file()
    measured = json.loads((root / "checks" / "jcomb" / "measurements.json").read_text())
    extended = measured["dust_ff"]["extended"]
    assert extended["zero_spacing"]["tp_only"]["n"] > 0
    assert extended["tp_only_ratio_to_scaled_tp"]["n"] > 0
    assert 0.99 < extended["ratio_to_standard_pb_ge_0.35"] < 1.01


def test_extended_outputs_guarded(tmp_path, monkeypatch):
    paths, calls = _method_inputs(tmp_path), []
    _patch_method_stages(monkeypatch, calls)
    config = _config({})
    stages.run_method("jcomb", paths, config)
    assert len(calls) == 4                      # standard and extended, for both flavours
    extended = paths["products"]["jcomb_dust_extended"]
    extended.unlink()                           # only an extended output is left to collide
    (paths["output_root"] / "products" / "jcomb_dust_224GHz.fits").unlink()
    paths["products"]["jcomb_continuum_extended"].write_bytes(b"x")
    paths["products"]["jcomb_continuum"].unlink()
    with pytest.raises(FileExistsError, match="extended"):
        stages.run_method("jcomb", paths, config)
    assert len(calls) == 4


def test_prepare_rejects_gain_of_another_shape(tmp_path):
    paths = _inputs(tmp_path)
    fits.PrimaryHDU(np.full((300, 300), 0.9, dtype=np.float32)[None, None], _gain_header(300)).writeto(
        paths["gain"], overwrite=True)
    with pytest.raises(ValueError, match="shape"):
        stages.run_prepare(paths, _prepare_config())


def test_measurements_run_on_a_product_blank_outside_the_mosaic():
    from astropy.wcs import WCS
    n, pix = 64, 2.0
    config = _config({})
    header, tp_header = _header(n, pix, 4.0), _header(16, 4 * pix, 18.1)
    gain = np.zeros((n, n), np.float32)
    gain[16:48, 16:48] = 1.0
    field = np.nan_to_num(_smooth_field(n, border=1), nan=1.0)
    product = np.where(gain > 0, field, np.nan)       # blank outside the mosaic
    edges = np.geomspace(1.0, 30.0, 7)
    ctx = {"high": 1.0e3 * field, "low": np.ones((n, n), np.float32), "tp_high": np.ones((n, n)),
           "tp_native": np.ones((16, 16)), "tp_wcs": WCS(tp_header), "header": header, "gain": gain,
           "blank": gain == 0, "mask": gain > 0, "window": (gain > 0).astype(float), "pix": pix,
           "high_beam": (4.0, 4.0, 0.0), "clean": np.nan_to_num(field), "positions": [], "sigma": 0.1,
           "edges": edges, "centres": np.sqrt(edges[:-1] * edges[1:])}
    measured = stages.measure_product(product, ctx, config, 1.0)
    assert measured["tp_only_ratio_to_scaled_tp"]["n"] == 0
    assert measured["tp_only_ratio_to_scaled_tp"]["median"] is None
    check = stages.zero_spacing(product, ctx, config, 1.0)[0]
    assert check["tp_only"] == {"n": 0}
