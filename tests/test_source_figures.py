"""Tests for the source-map and decomposition figure geometry."""

from astropy.wcs import WCS
import numpy as np
import pytest

from lancet.catalogues import remove_source_model, source_model_true_sky
from lancet.source_figures import box_window, ellipse_geometry


def _north_up_wcs(pixel_arcsec: float = 0.1) -> WCS:
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    wcs.wcs.crval = [221.25, -59.8]
    wcs.wcs.crpix = [501, 501]
    wcs.wcs.cdelt = [-pixel_arcsec / 3600, pixel_arcsec / 3600]
    return wcs


def test_ellipse_angle_follows_position_angle_east_of_north() -> None:
    wcs = _north_up_wcs()

    north = ellipse_geometry(221.25, -59.8, 1.0, 0.5, 0.0, wcs)
    east = ellipse_geometry(221.25, -59.8, 1.0, 0.5, 90.0, wcs)

    assert north.angle_deg == pytest.approx(90.0)
    assert east.angle_deg == pytest.approx(180.0)
    assert north.x == pytest.approx(500.0) and north.y == pytest.approx(500.0)


def test_ellipse_size_is_fwhm_in_pixels_times_scale() -> None:
    ellipse = ellipse_geometry(221.25, -59.8, 1.2, 0.6, 30.0, _north_up_wcs(), scale=1.5)

    assert ellipse.width_px == pytest.approx(18.0, rel=1e-3)
    assert ellipse.height_px == pytest.approx(9.0, rel=1e-3)


def test_ellipse_refuses_a_rotated_grid() -> None:
    wcs = _north_up_wcs()
    wcs.wcs.pc = [[0.0, -1.0], [1.0, 0.0]]

    with pytest.raises(ValueError, match="north-up"):
        ellipse_geometry(221.25, -59.8, 1.0, 0.5, 0.0, wcs)


def test_box_window_is_square_and_clipped_to_the_image() -> None:
    wcs = _north_up_wcs(pixel_arcsec=0.5)

    rows, cols = box_window(wcs, 221.25, -59.8, 2.0 / 60, shape=(1000, 1000))
    edge_rows, edge_cols = box_window(wcs, 221.25, -59.8 + 0.07, 2.0 / 60, shape=(1000, 1000))

    assert rows.stop - rows.start == cols.stop - cols.start == 240
    assert edge_rows.stop == 1000 and edge_rows.start > 0


def test_source_model_true_sky_is_what_remove_source_model_subtracts() -> None:
    image = np.full((2, 3), 10.0)
    sources = np.array([[0.0, 2.0, 4.0], [1.0, np.nan, 0.0]])
    gain = np.array([[0.5, 0.5, 1.0], [0.25, 1.0, 1.0]])

    true_sky = source_model_true_sky(sources, gain)

    np.testing.assert_allclose(true_sky, [[0.0, 4.0, 4.0], [4.0, 0.0, 0.0]])
    np.testing.assert_allclose(remove_source_model(image, sources, gain), image - true_sky)
