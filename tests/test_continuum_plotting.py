"""Tests for all-scale continuum plotting utilities."""

from pathlib import Path
import subprocess
import sys

from astropy.io import fits
import matplotlib
import numpy as np
import pytest

from lancet.continuum_plotting import (
    add_primary_beam_contour,
    create_dual_colour_figure,
    highlighted_emission,
    image_plane,
    world_crop_limits,
)


matplotlib.use("Agg")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _header() -> fits.Header:
    return fits.Header(
        {
            "NAXIS": 2,
            "NAXIS1": 41,
            "NAXIS2": 31,
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRPIX1": 21.0,
            "CRPIX2": 16.0,
            "CRVAL1": 221.2346,
            "CRVAL2": -59.7738,
            "CDELT1": -1.0 / 3600.0,
            "CDELT2": 1.0 / 3600.0,
            "BMAJ": 1.0 / 3600.0,
            "BMIN": 0.8 / 3600.0,
            "BPA": 0.0,
            "BUNIT": "MJy/sr",
        }
    )


def test_image_plane_squeezes_only_degenerate_axes() -> None:
    cube = np.zeros((1, 3, 4, 1), dtype=np.float32)

    image = image_plane(cube)

    assert image.shape == (3, 4)
    assert image.dtype == np.float32
    with pytest.raises(ValueError, match="two-dimensional"):
        image_plane(np.zeros((2, 3, 4)))


def test_highlight_mask_includes_transition_and_nonfinite_values() -> None:
    data = np.array([[0.0, 2.0], [3.0, np.nan]])

    highlighted = highlighted_emission(data, transition=2.0)

    np.testing.assert_array_equal(
        np.ma.getmaskarray(highlighted),
        [[True, True], [False, True]],
    )


def test_world_crop_is_centred_on_reference_pixel() -> None:
    x_limits, y_limits = world_crop_limits(
        _header(),
        center_ra_deg=221.2346,
        center_dec_deg=-59.7738,
        width_deg=20.0 / 3600.0,
        height_deg=10.0 / 3600.0,
    )

    assert sum(x_limits) / 2.0 == pytest.approx(20.0)
    assert sum(y_limits) / 2.0 == pytest.approx(15.0)
    assert x_limits[1] - x_limits[0] == pytest.approx(20.0)
    assert y_limits[1] - y_limits[0] == pytest.approx(10.0)


def test_clean_figure_can_be_rendered(tmp_path: Path) -> None:
    data = np.linspace(-1.0, 10.0, 31 * 41).reshape(31, 41)
    figure = create_dual_colour_figure(
        data,
        _header(),
        center_ra_deg=221.2346,
        center_dec_deg=-59.7738,
        width_deg=20.0 / 3600.0,
        height_deg=10.0 / 3600.0,
        background_vmin=-1.0,
        transition=3.0,
        highlight_vmax=10.0,
        annotated=False,
    )
    output = tmp_path / "clean.png"

    figure.savefig(output, dpi=50)

    assert output.is_file()
    assert output.stat().st_size > 0


def test_primary_beam_contour_uses_its_own_wcs(tmp_path: Path) -> None:
    from astropy.wcs import WCS
    from matplotlib import pyplot as plt

    image_header = _header()
    beam_header = _header()
    beam_header["CRPIX1"] = 18.0
    response = np.tile(np.linspace(0.0, 1.0, 41), (31, 1))
    path = tmp_path / "primary_beam.fits"
    fits.writeto(path, response, beam_header)
    figure = plt.figure()
    axis = figure.add_subplot(projection=WCS(image_header).celestial)

    contour = add_primary_beam_contour(axis, path, minimum_response=0.2)

    # The 0.2 contour is at native x=8, shifted by three image pixels.
    np.testing.assert_allclose(contour.levels, [0.2])
    np.testing.assert_allclose(
        contour.get_paths()[0].vertices[:, 0], 11.0, atol=1e-6
    )
    plt.close(figure)


def test_allscale_plot_inspection_is_read_only() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "scripts/plot_allscale_continuum.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "allscale_combined_v3_dust_ff.MJysr.fits" in completed.stdout
    assert "synthetic_continuum" in completed.stdout
    assert "Output root:" in completed.stdout
