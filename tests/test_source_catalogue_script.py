"""Tests for the source_catalogue task runner."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

from astropy.io import fits
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "tasks" / "source_catalogue" / "run.py"
SPEC = spec_from_file_location("source_catalogue_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "source_catalogue" / "task.toml"
DELIVERED = "data/delivered/alma/continuum/"
GETSF = "data/extraction/getsf/triple/results/"


def test_task_reads_extraction_from_data_and_new_combination() -> None:
    inputs = TASK.load_task_config(CONFIG)["inputs"]

    assert inputs["getsf_catalogue"] == GETSF + "G316.sw.sources.ok.cat"
    assert inputs["intercore_background"] == (
        GETSF + "+components/G316.sw.1302.obs.sbackground.fits"
    )
    assert inputs["source_model"] == (
        GETSF + "+components/G316.sw.1302.obs.sources.fits"
    )
    assert inputs["combined_continuum"] == (
        DELIVERED + "IRDC316.8.TP7MTM2TM1.robust0.5.cont.image.pbcor.fits"
    )
    assert inputs["gain"] == (
        DELIVERED + "IRDC316.8.7MTM2TM1.robust0.5.cont.gain.fits"
    )
    assert "core_background" not in inputs


def _header() -> fits.Header:
    header = fits.Header()
    for key, value in {
        "CTYPE1": "RA---SIN", "CTYPE2": "DEC--SIN",
        "CRVAL1": 221.2, "CRVAL2": -59.8, "CRPIX1": 2.0, "CRPIX2": 2.0,
        "CDELT1": -2.7777778e-5, "CDELT2": 2.7777778e-5,
        "CUNIT1": "deg", "CUNIT2": "deg",
        "BMAJ": 1.97858e-4, "BMIN": 1.68404e-4, "BPA": 47.4,
    }.items():
        header[key] = value
    return header


def test_background_stage_writes_float32_true_sky_background(
    tmp_path: Path,
) -> None:
    header = _header()
    combined = header.copy()
    combined["BUNIT"] = "Jy/beam"
    fits.writeto(tmp_path / "c.fits", np.full((3, 3), 1.0e-3, np.float32), combined)
    model = header.copy()
    model["BUNIT"] = "MJy/sr"
    sources = np.zeros((3, 3), np.float32)
    sources[1, 1] = 10.0
    fits.writeto(tmp_path / "s.fits", sources, model)
    gain = np.ones((1, 1, 3, 3), np.float32)
    gain[0, 0, 1, 1] = 0.5
    fits.writeto(tmp_path / "g.fits", gain, header)
    paths = {
        "combined_continuum": tmp_path / "c.fits",
        "source_model": tmp_path / "s.fits",
        "gain": tmp_path / "g.fits",
        "output_root": tmp_path / "out",
    }

    output = TASK.build_background(paths)

    data, written = fits.getdata(output, header=True)
    assert data.dtype == np.dtype(">f4")
    assert written["BUNIT"] == "MJy/sr"
    assert data[0, 0] == pytest.approx(data[1, 1] + 20.0, rel=1e-5)
    with pytest.raises(FileExistsError):
        TASK.build_background(paths)
