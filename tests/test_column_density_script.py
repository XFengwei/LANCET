"""Tests for the column_density task runner."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pytest

from lancet import ProjectPaths
from lancet.imaging import jybeam_to_mjysr_from_header


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "tasks" / "column_density" / "run.py"
SPEC = spec_from_file_location("column_density_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
TASK = module_from_spec(SPEC)
sys.modules[SPEC.name] = TASK
SPEC.loader.exec_module(TASK)

CONFIG = PROJECT_ROOT / "tasks" / "column_density" / "task.toml"


def _celestial_header(ny: int, nx: int) -> fits.Header:
    header = fits.Header()
    header["CTYPE1"] = "RA---SIN"
    header["CTYPE2"] = "DEC--SIN"
    header["CRVAL1"] = 221.2
    header["CRVAL2"] = -59.8
    header["CRPIX1"] = nx / 2
    header["CRPIX2"] = ny / 2
    header["CDELT1"] = -1.0e-4
    header["CDELT2"] = 1.0e-4
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["BMAJ"] = 3.0e-4
    header["BMIN"] = 2.0e-4
    header["BPA"] = 40.0
    return header


def test_grid_header_drops_degenerate_axes(tmp_path: Path) -> None:
    header = _celestial_header(3, 4)
    header["CTYPE3"] = "FREQ"
    header["CRVAL3"] = 2.238e11
    header["CDELT3"] = 1.0e9
    header["CRPIX3"] = 1.0
    header["CTYPE4"] = "STOKES"
    header["CRVAL4"] = 1.0
    header["CDELT4"] = 1.0
    header["CRPIX4"] = 1.0
    header["EPOCH"] = 2000.0
    path = tmp_path / "cube.fits"
    fits.writeto(path, np.zeros((1, 1, 3, 4), dtype=np.float32), header)

    grid = TASK.grid_header(path)

    assert grid["NAXIS"] == 2
    assert (grid["NAXIS1"], grid["NAXIS2"]) == (4, 3)
    assert WCS(grid).naxis == 2
    assert grid["CRPIX1"] == pytest.approx(2.0)
    assert grid["CRPIX2"] == pytest.approx(1.5)
    assert WCS(grid).wcs.equinox == 2000.0


def _write_map_inputs(root: Path, allscale_unit: str) -> dict[str, Path]:
    header = _celestial_header(3, 4)
    temperature = root / "out" / "prepared" / "temperature_allscale.fits"
    temperature.parent.mkdir(parents=True)
    fits.writeto(temperature, np.full((3, 4), 20.0, np.float32), header)

    mjysr = np.linspace(1.0, 50.0, 12, dtype=np.float32).reshape(3, 4)
    allscale_header = header.copy()
    allscale_header["BUNIT"] = allscale_unit
    if allscale_unit == "Jy/beam":
        scale = jybeam_to_mjysr_from_header(np.ones((3, 4)), header)
        allscale = (mjysr / scale).astype(np.float32)
    else:
        allscale = mjysr
    allscale_path = root / f"allscale_{allscale_unit.replace('/', '')}.fits"
    fits.writeto(allscale_path, allscale, allscale_header)

    intercore_header = header.copy()
    intercore_header["BUNIT"] = "MJy/sr"
    intercore_path = root / "intercore.fits"
    fits.writeto(intercore_path, mjysr / 2, intercore_header)
    return {
        "dust_continuum_allscale": allscale_path,
        "intercore_background": intercore_path,
        "output_root": root / "out",
    }


def test_allscale_jybeam_input_is_converted_to_mjysr(tmp_path: Path) -> None:
    config = TASK.load_task_config(CONFIG)
    results = {}
    for unit in ("MJy/sr", "Jy/beam"):
        root = tmp_path / unit.replace("/", "")
        root.mkdir()
        paths = _write_map_inputs(root, unit)
        products = TASK.build_maps(paths, config)
        results[unit] = fits.getdata(products["column_density_allscale"])

    assert np.all(np.isfinite(results["MJy/sr"]))
    np.testing.assert_allclose(results["Jy/beam"], results["MJy/sr"], rtol=1e-5)


def _fake_project(tmp_path: Path, allscale: str) -> tuple[ProjectPaths, dict]:
    config = TASK.load_task_config(CONFIG)
    inputs = dict(config["inputs"])
    inputs["dust_continuum_allscale"] = allscale
    for value in inputs.values():
        path = tmp_path / value
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"")
    config = {**config, "inputs": inputs}
    return ProjectPaths(tmp_path), config


def test_allscale_may_come_from_combine_continuum(tmp_path: Path) -> None:
    project, config = _fake_project(
        tmp_path,
        "results/combine_continuum/products/jcomb_dust_224GHz_extended.fits",
    )

    paths = TASK.resolve_task_paths(project, config)

    assert paths["dust_continuum_allscale"].name == (
        "jcomb_dust_224GHz_extended.fits"
    )


def test_allscale_from_another_task_is_rejected(tmp_path: Path) -> None:
    project, config = _fake_project(
        tmp_path,
        "results/synthesize_continuum/model/dust_224GHz.fits",
    )

    with pytest.raises(ValueError, match="combine_continuum"):
        TASK.resolve_task_paths(project, config)


def test_task_reads_new_combination_and_miriad_products() -> None:
    config = TASK.load_task_config(CONFIG)
    inputs = config["inputs"]

    assert inputs["dust_continuum_allscale"] == (
        "results/combine_continuum/products/jcomb_dust_224GHz_extended.fits"
    )
    assert inputs["grid_reference"] == (
        "data/delivered/alma/continuum/"
        "IRDC316.8.7MTM2TM1.robust0.5.cont.clean.fits"
    )
    assert inputs["primary_beam"] == (
        "data/delivered/alma/continuum/"
        "IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits"
    )
    assert {"pdf", "png"} <= set(config["figures"]["formats"])


def test_task_floors_nonpositive_intensity() -> None:
    config = TASK.load_task_config(CONFIG)

    assert config["model"]["zero_nonpositive_intensity"] is True
