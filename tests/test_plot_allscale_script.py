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


def test_sqrt_style_is_single_inferno_scale() -> None:
    style = TASK.load_task_config(CONFIG)["styles"]["sqrt"]

    assert style["colour_map"] == "inferno"
    assert style["stretch"] == "sqrt"
    assert style["vmin"] < style["vmax"]


def test_sqrt_stage_writes_annotated_and_clean_views(tmp_path, monkeypatch) -> None:
    calls = []

    class FakeFigure:
        def savefig(self, output, **kwargs):
            Path(output).write_bytes(b"x")

    def fake_figure(data, header, **kwargs):
        calls.append(kwargs)
        return FakeFigure()

    image = tmp_path / "image.fits"
    _write_image(image, np.ones((3, 3)), "MJy/sr")
    monkeypatch.setattr(TASK, "create_single_scale_figure", fake_figure)
    monkeypatch.setattr(TASK, "_close", lambda figure: None)
    config = TASK.load_task_config(CONFIG)
    paths = {
        "combined_continuum": image,
        "primary_beam": None,
        "output_root": tmp_path / "out",
    }

    written = TASK.plot_sqrt(paths, config)

    names = sorted(path.name for path in written)
    assert names == [
        "allscale_continuum_dust_freefree_sqrt.pdf",
        "allscale_continuum_dust_freefree_sqrt.png",
        "allscale_continuum_dust_freefree_sqrt_clean.png",
    ]
    assert [call["annotated"] for call in calls] == [True, False]
    assert all(call["colour_map"] == "inferno" for call in calls)
