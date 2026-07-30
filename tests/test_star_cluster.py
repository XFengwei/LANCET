"""Tests for the standalone star-cluster orbit script."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

from astropy import units as u
from astropy.table import Table
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "star_cluster.py"
SPEC = spec_from_file_location("star_cluster_script", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
STAR_CLUSTER = module_from_spec(SPEC)
sys.modules[SPEC.name] = STAR_CLUSTER
SPEC.loader.exec_module(STAR_CLUSTER)


def test_angle_wrap_uses_minus_180_to_180_interval() -> None:
    wrapped = STAR_CLUSTER.angle_wrap([0, 180, 270, 360])
    assert wrapped == pytest.approx([0, -180, -90, 0])


def test_orbit_selection_uses_one_id_and_past_time_window() -> None:
    catalogue = Table(
        {
            "ID": [1.0, 1.0, 1.0, 2.0],
            "time": [-5.0, -4.9, 0.0, -1.0],
            "l": np.zeros(4),
            "b": np.zeros(4),
            "ra": np.zeros(4),
            "dec": np.zeros(4),
        }
    )
    orbit = STAR_CLUSTER.select_orbit(catalogue, 1.0, 5.0)

    assert len(orbit) == 2
    assert list(orbit["time"]) == [-4.9, 0.0]


def test_delivered_catalogue_matches_reference_selection() -> None:
    catalogue = STAR_CLUSTER.load_catalogue(
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "starcluster"
        / "orbit_onlyback_OC100Myr_step01.fits"
    )
    orbit = STAR_CLUSTER.select_orbit(catalogue, 3944.0, 5.0)
    center = STAR_CLUSTER.filament_center()
    nearby = STAR_CLUSTER.close_approaches(orbit, center, 10 * u.arcmin)

    assert len(catalogue) == 581650
    assert len(orbit) == 50
    assert len(nearby) == 3
    assert nearby["sep_arcmin"][0] == pytest.approx(4.49653719)


def test_plot_functions_write_both_diagnostics(tmp_path: Path) -> None:
    catalogue = STAR_CLUSTER.load_catalogue(
        PROJECT_ROOT
        / "data"
        / "delivered"
        / "starcluster"
        / "orbit_onlyback_OC100Myr_step01.fits"
    )
    orbit = STAR_CLUSTER.select_orbit(catalogue, 3944.0, 5.0)
    center = STAR_CLUSTER.filament_center()
    whole_sky = tmp_path / "whole_sky.png"
    zoom = tmp_path / "zoom.png"

    STAR_CLUSTER.plot_whole_sky(
        orbit, center, whole_sky, cluster_id=3944.0
    )
    STAR_CLUSTER.plot_zoom(
        orbit,
        center,
        zoom,
        cluster_id=3944.0,
        radius_arcmin=10.0,
        zoom_deg=3.0,
    )

    assert whole_sky.is_file() and whole_sky.stat().st_size > 0
    assert zoom.is_file() and zoom.stat().st_size > 0
