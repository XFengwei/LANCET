import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from scipy import fft

from lancet.combination.uvplane import gaussian_transfer, uv_distance_klambda
from lancet.combination.prep import (
    azimuthal_power,
    derive_gain,
    footprint_margin_arcsec,
    jybeam_to_mjysr,
    mjysr_to_jybeam,
    pad,
    read_plane,
    regrid_to,
    taper_edge,
)


def _celestial_header(nx, ny, pixel_arcsec):
    h = fits.Header()
    h["CTYPE1"], h["CTYPE2"] = "RA---SIN", "DEC--SIN"
    h["CRVAL1"], h["CRVAL2"] = 221.25, -59.8
    h["CRPIX1"], h["CRPIX2"] = (nx + 1) / 2, (ny + 1) / 2
    h["CDELT1"], h["CDELT2"] = -pixel_arcsec / 3600, pixel_arcsec / 3600
    h["CUNIT1"], h["CUNIT2"] = "deg", "deg"
    return h


def test_read_plane_miriad_layout(tmp_path):
    h = _celestial_header(30, 20, 0.1)
    h["CTYPE3"], h["CRVAL3"], h["CDELT3"], h["CRPIX3"] = "VELO-LSR", 0.0, 1.0, 1.0
    h["CTYPE4"], h["CRVAL4"], h["CDELT4"], h["CRPIX4"] = "STOKES", 1.0, 1.0, 1.0
    h["BUNIT"], h["BMAJ"], h["BMIN"], h["BPA"] = "JY/BEAM", 0.7 / 3600, 0.6 / 3600, 47.0
    path = tmp_path / "miriad.fits"
    fits.PrimaryHDU(np.ones((1, 1, 20, 30), dtype=">f8"), h).writeto(path)
    data, header, wcs = read_plane(path)
    assert data.shape == (20, 30)
    assert data.dtype == np.float32
    assert wcs.naxis == 2
    h["BUNIT"] = "K"
    fits.PrimaryHDU(np.ones((1, 1, 20, 30)), h).writeto(tmp_path / "k.fits")
    with pytest.raises(ValueError):
        read_plane(tmp_path / "k.fits")


def test_unit_round_trip():
    x = np.linspace(-1, 3, 50)
    np.testing.assert_allclose(mjysr_to_jybeam(jybeam_to_mjysr(x, 0.712, 0.606), 0.712, 0.606), x,
                               rtol=1e-6)
    np.testing.assert_allclose(jybeam_to_mjysr(1.0, 18.1, 18.1), 114.6116, rtol=1e-4)


def test_derive_gain():
    g = np.linspace(0.1, 1.0, 400).reshape(20, 20)
    clean = np.full((20, 20), 2.0)
    clean[:5, :5] = np.nan
    gain = derive_gain(clean, 2.0 * g)
    np.testing.assert_allclose(gain[5:, 5:], g[5:, 5:], atol=1e-6)
    assert np.all(gain[:5, :5] == 0)
    with pytest.raises(ValueError):
        derive_gain(clean, 2.2 * clean)


def _gaussian_image(shape, pixel_arcsec, fwhm_arcsec, centre):
    y, x = np.indices(shape, dtype=float)
    r2 = ((x - centre[0]) ** 2 + (y - centre[1]) ** 2) * pixel_arcsec**2
    return np.exp(-4 * np.log(2) * r2 / fwhm_arcsec**2)


def test_regrid_preserves_power():
    coarse = _gaussian_image((61, 61), 6.0, 18.1, (30, 30))
    hc = _celestial_header(61, 61, 6.0)
    hf = _celestial_header(300, 300, 1.0)
    fine = regrid_to(coarse, WCS(hc), WCS(hf), (300, 300))
    truth = _gaussian_image((300, 300), 1.0, 18.1, (149.5, 149.5))
    edges = np.arange(0.0, 10.0, 1.0)
    p_fine = azimuthal_power(np.nan_to_num(fine), 1.0, edges)
    p_true = azimuthal_power(truth, 1.0, edges)
    np.testing.assert_allclose(p_fine, p_true, rtol=0.02)


def test_margin_negative_when_mosaic_outside():
    y, x = np.indices((200, 200))
    footprint = (x - 100) ** 2 + (y - 100) ** 2 <= 50**2
    inside = np.zeros_like(footprint)
    inside[80:120, 80:120] = True
    assert footprint_margin_arcsec(footprint, inside, 0.5) > 0
    beyond = np.zeros_like(footprint)
    beyond[100, 90:154] = True  # the row through the disc tip, 3 px past x = 150
    assert footprint_margin_arcsec(footprint, beyond, 0.5) == pytest.approx(-3 * 0.5, abs=0.5)


def test_taper_edge_shape():
    y, x = np.indices((101, 101))
    footprint = x < 50
    data = np.where(footprint, 2.0, np.nan)
    out = taper_edge(data, footprint, width_arcsec=10.0, pixel_arcsec=1.0)
    np.testing.assert_array_equal(out[footprint], 2.0)
    assert np.all(out[:, 60:] == 0)
    row = out[50, 49:62]
    assert np.all(np.diff(row) <= 0)
    assert np.all(np.isfinite(out))


def test_pad_no_wraparound():
    img = np.zeros((50, 50))
    img[0, 0] = 1.0
    padded, cut = pad(img, 9)
    u, v, _ = uv_distance_klambda(padded.shape, 1.0)
    fwhm = 3.0 * np.sqrt(8 * np.log(2))  # sigma = 3 px at 1"/px
    out = fft.irfft2(fft.rfft2(padded) * gaussian_transfer(u, v, fwhm, fwhm, 0.0), s=padded.shape)[cut]
    assert out[0, 49] < 1e-6 * out.max()
    assert out[49, 0] < 1e-6 * out.max()


def test_write_product_header(tmp_path):
    from lancet.combination.prep import write_product

    h = _celestial_header(30, 20, 0.1)
    h["CTYPE3"], h["CRVAL3"], h["CDELT3"], h["CRPIX3"] = "VELO-LSR", 0.0, 1.0, 1.0
    h["NAXIS"], h["NAXIS1"], h["NAXIS2"], h["NAXIS3"] = 3, 30, 20, 1
    h["BUNIT"], h["BMAJ"], h["BMIN"], h["BPA"] = "JY/BEAM", 0.7 / 3600, 0.6 / 3600, 47.0
    path = tmp_path / "product.fits"
    write_product(path, np.ones((20, 30)), h, {"method": "jcomb", "k1": 5.0}, {"interferometer": "abc123"})
    with fits.open(path) as hdul:
        out = hdul[0].header
        assert out["BITPIX"] == -32
        assert out["NAXIS"] == 2
        assert out["BUNIT"] == "Jy/beam"
        assert out["BMAJ"] == pytest.approx(0.7 / 3600)
        assert out["COMB_METHOD"] == "jcomb"
        assert out["COMB_K1"] == 5.0
        assert any("abc123" in str(line) for line in out["HISTORY"])


def test_combine_prepared_end_to_end_small():
    from lancet.combination.prep import combine_prepared

    shape, pix = (200, 240), 0.5
    rng = np.random.default_rng(3)
    sky = rng.random(shape) + 0.5
    u, v, _ = uv_distance_klambda(shape, pix)
    low_beam, high_beam = (18.1, 18.1, 0.0), (2.0, 1.5, 20.0)

    def conv(beam):
        return fft.irfft2(fft.rfft2(sky) * gaussian_transfer(u, v, *beam), s=shape)

    settings = {"k1": 5.0, "k2": 9.0, "order": 5, "tp_factor": 1.0, "minimum_tp_response": 0.1,
                "padding_arcsec": 0.0}
    blank = np.zeros(shape, dtype=bool)
    blank[:, :10] = True
    out = combine_prepared(conv(low_beam), conv(high_beam), pix, high_beam, low_beam, "jcomb", settings,
                           blank=blank)
    expected = mjysr_to_jybeam(conv(high_beam), 2.0, 1.5)
    assert np.all(np.isnan(out[blank]))
    np.testing.assert_allclose(out[~blank], expected[~blank], rtol=0, atol=1e-3 * expected.max())


def test_jcomb_tp_alone_has_no_grid_edge_ringing():
    """A smooth TP cut by the grid edge must come back unchanged where there is no
    interferometer: zero-padding the TP leaves a step that J-comb's G_high/G_low rings on."""
    from lancet.combination.prep import combine_prepared

    shape, pix = (300, 400), 0.5
    yy, xx = np.indices(shape)
    low = 20.0 + 30.0 * np.exp(-0.5 * (((xx - 500) * pix / 120.0) ** 2 + ((yy - 120) * pix / 90.0) ** 2))
    high = np.zeros(shape)
    high_beam, low_beam = (0.71, 0.61, 47.4), (18.1, 18.1, 0.0)
    settings = {"k1": 5.0, "k2": 9.0, "order": 5, "tp_factor": 1.0, "padding_arcsec": 120.0}
    out = combine_prepared(low, high, pix, high_beam, low_beam, "jcomb", settings)
    error = np.abs(out / mjysr_to_jybeam(low, 0.71, 0.61) - 1)
    assert error.max() < 0.03                     # zero-padding the TP gives 0.72
    assert error[40:-40, 40:-40].max() < 0.015    # more than 20" from the grid edge


def test_derive_gain_keeps_exact_ratio_at_weak_pixels():
    g = np.linspace(0.1, 1.0, 400).reshape(20, 20)
    clean = np.full((20, 20), 2.0)
    clean[10, 10] = 2e-4                       # 1e-4 of the peak, but exactly clean x gain
    gain = derive_gain(clean, clean * g)
    assert abs(gain[10, 10] - g[10, 10]) < 1e-6
