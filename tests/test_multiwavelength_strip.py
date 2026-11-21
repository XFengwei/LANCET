"""Tests for the multiwavelength_strip task (Paper II Fig. 1)."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import tomllib

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from lancet.imaging import jybeam_to_mjysr_from_header


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "multiwavelength_strip_script",
    PROJECT_ROOT / "tasks" / "multiwavelength_strip" / "run.py",
)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "multiwavelength_strip" / "task.toml"


def _header() -> fits.Header:
    header = fits.Header()
    for key, value in {
        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
        "CRVAL1": 221.2, "CRVAL2": -59.8, "CRPIX1": 10.5, "CRPIX2": 10.5,
        "CDELT1": -1.0e-4, "CDELT2": 1.0e-4,
        "BMAJ": 2.0e-4, "BMIN": 2.0e-4, "BPA": 0.0, "BUNIT": "Jy/beam",
    }.items():
        header[key] = value
    return header


def test_panel_image_can_be_converted_from_jybeam(tmp_path: Path) -> None:
    header = _header()
    data = np.full((20, 20), 2.0e-3, np.float32)
    path = tmp_path / "jybeam.fits"
    fits.writeto(path, data, header)
    wcs_out = WCS(header).celestial

    image = TASK.reproject_image(path, wcs_out, (20, 20), to_mjysr=True)

    expected = jybeam_to_mjysr_from_header(np.array(2.0e-3), header)
    np.testing.assert_allclose(image[10, 10], expected, rtol=1e-5)


def test_mm_panel_reads_the_delivered_jcomb_product() -> None:
    config = tomllib.loads(CONFIG.read_text())
    panel = next(p for p in config["panels"] if p["name"] == "mm")

    assert panel["input"] == (
        "data/delivered/alma/continuum/"
        "IRDC316.8.TP7MTM2TM1.robust0.5.cont.image.pbcor.extended36.fits"
    )
    assert panel["to_mjysr"] is True
    assert "fill_input" not in panel


def test_fill_gaps_takes_model_only_where_image_is_blank() -> None:
    image = np.array([[1.0, np.nan], [np.nan, 4.0]], np.float32)
    model = np.array([[9.0, 2.0], [3.0, 9.0]], np.float32)

    filled = TASK.fill_gaps(image, model)

    np.testing.assert_array_equal(filled, [[1.0, 2.0], [3.0, 4.0]])


def test_mm_panel_fills_gaps_from_the_tp_model() -> None:
    config = tomllib.loads(CONFIG.read_text())
    panel = next(p for p in config["panels"] if p["name"] == "mm")

    assert panel["gap_fill_input"] == (
        "results/synthesize_continuum/model/continuum_224GHz_fill36.fits"
    )


def test_pb_blend_fades_from_model_to_image() -> None:
    image = np.full(5, 10.0, np.float32)
    model = np.full(5, 2.0, np.float32)
    pb = np.array([np.nan, 0.1, 0.35, 0.5, 0.9])

    blended = TASK.pb_blend(image, model, pb, (0.2, 0.5))

    np.testing.assert_allclose(blended, [2.0, 2.0, 6.0, 10.0, 10.0])


def test_mm_panel_blends_across_the_mosaic_edge() -> None:
    config = tomllib.loads(CONFIG.read_text())
    panel = next(p for p in config["panels"] if p["name"] == "mm")

    assert panel["blend_pb_input"] == (
        "data/delivered/alma/continuum/"
        "IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits"
    )
    assert panel["blend_pb_range"] == [0.2, 0.5]


def test_sqrt_scale_maps_absolute_range() -> None:
    image = np.array([-5.0, 0.0, 300.0, 1200.0, 5000.0, np.nan])

    scaled = TASK.sqrt_scale(image, 0.0, 1200.0)

    np.testing.assert_allclose(scaled[:5], [0.0, 0.0, 0.5, 1.0, 1.0])
    assert np.isnan(scaled[5])


def test_mm_panel_uses_the_fig3_inferno_view() -> None:
    config = tomllib.loads(CONFIG.read_text())
    panel = next(p for p in config["panels"] if p["name"] == "mm")

    assert panel["style"] == "sqrt"
    assert panel["cmap"] == "inferno"
    assert (panel["vmin"], panel["vmax"]) == (0.0, 1200.0)
