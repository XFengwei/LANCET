"""Tests for the plot_allscale_continuum task runner."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

from astropy.io import fits
import numpy as np

from lancet.imaging import jybeam_to_mjysr_from_header


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "tasks" / "plot_allscale_continuum" / "run.py"
SPEC = spec_from_file_location("plot_allscale_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "plot_allscale_continuum" / "task.toml"
DELIVERED = "data/delivered/alma/continuum/"


def _write_image(path: Path, data: np.ndarray, unit: str) -> None:
    header = fits.Header()
    for key, value in {
        "CTYPE1": "RA---SIN", "CTYPE2": "DEC--SIN",
        "CRVAL1": 221.2, "CRVAL2": -59.8, "CRPIX1": 2.0, "CRPIX2": 2.0,
        "CDELT1": -1.0e-4, "CDELT2": 1.0e-4,
        "BMAJ": 3.0e-4, "BMIN": 2.0e-4, "BPA": 0.0, "BUNIT": unit,
    }.items():
        header[key] = value
    fits.writeto(path, data.astype(np.float32), header)


def test_jybeam_image_is_plotted_in_mjysr(tmp_path, monkeypatch) -> None:
    mjysr = np.arange(1.0, 10.0).reshape(3, 3)
    header = fits.Header({"BMAJ": 3.0e-4, "BMIN": 2.0e-4})
    jybeam = mjysr / jybeam_to_mjysr_from_header(np.ones((3, 3)), header)
    path = tmp_path / "image.fits"
    _write_image(path, jybeam, "Jy/beam")
    captured = {}

    def fake_figure(data, header, **kwargs):
        captured["data"] = np.asarray(data)
        captured["unit"] = header["BUNIT"]
        return None

    monkeypatch.setattr(TASK, "create_dual_colour_figure", fake_figure)
    config = TASK.load_task_config(CONFIG)
    TASK._create_figure(
        path, {"primary_beam": None}, config,
        style_name="allscale", annotated=False,
    )

    np.testing.assert_allclose(captured["data"], mjysr, rtol=1e-5)
    assert captured["unit"] == "MJy/sr"


def test_task_reads_delivered_jcomb_and_miriad_pb() -> None:
    config = TASK.load_task_config(CONFIG)
    inputs = config["inputs"]

    assert inputs["combined_continuum"] == (
        DELIVERED
        + "IRDC316.8.TP7MTM2TM1.robust0.5.cont.image.pbcor.extended.fits"
    )
    assert inputs["dust_continuum"] == (
        DELIVERED
        + "IRDC316.8.TP7MTM2TM1.robust0.5.cont_dust.image.pbcor.extended.fits"
    )
    assert inputs["primary_beam"] == (
        DELIVERED + "IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits"
    )
    assert {"pdf", "png"} <= set(config["figures"]["annotated_formats"])
