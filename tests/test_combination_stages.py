"""Tests for the combine_continuum phase 1 stages (src/lancet/combination_stages.py)."""

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
    fits.PrimaryHDU(np.full((n, n), 0.9, dtype=np.float32)[None, None],
                    _truesky_header(n, 0.5, 1.0)).writeto(paths["gain"])
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


def test_prepare_fails_on_small_margin(tmp_path):
    paths = _inputs(tmp_path, truesky_edge=0)       # 200" mosaic finite to the edge of the grid
    for name in ("tp_dust", "tp_dust_ff"):          # TP footprint 198" wide: it ends inside the mosaic
        header = _header(61, 6.0, 18.1)
        header["FREQGHZ"] = 223.8
        fits.PrimaryHDU(_smooth_field(border=14), header).writeto(paths[name], overwrite=True)
    with pytest.raises(RuntimeError, match="margin"):
        stages.run_prepare(paths, _prepare_config())
    assert not (tmp_path / "out").exists()          # nothing is written when an input check fails
