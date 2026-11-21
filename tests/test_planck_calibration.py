import numpy as np
import pytest

from lancet.combination import planck_calibration as pc


def test_dbdt_cmb_matches_numerical_derivative():
    nu = np.array([100e9, 217e9, 353e9])
    step = 1e-3
    numeric = (pc.planck_function(nu, pc.TCMB + step) - pc.planck_function(nu, pc.TCMB - step)) / (2 * step)
    assert np.allclose(pc.dbdt_cmb(nu), numeric, rtol=1e-6)


def test_kcmb_conversion_for_a_narrow_band_is_dbdt_times_the_colour_ratio():
    nu = np.linspace(216.9e9, 217.1e9, 201)
    tau = np.ones_like(nu)
    same = pc.kcmb_to_mjysr(nu, tau, 217e9, 1.8, 20.0)
    assert same == pytest.approx(pc.dbdt_cmb(217e9) * 1e20, rel=1e-4)
    shape = lambda f: f**1.8 * pc.planck_function(f, 20.0)
    assert pc.kcmb_to_mjysr(nu, tau, 224e9, 1.8, 20.0) / same == pytest.approx(shape(224e9) / shape(217e9), rel=1e-4)
    assert pc.kcmb_to_mjysr_iras(nu, tau, 217e9) == pytest.approx(same, rel=1e-4)


def test_smooth_with_coverage_tracks_the_footprint():
    data = np.full((80, 80), np.nan)
    data[:, :40] = 2.0
    smooth, coverage = pc.smooth_with_coverage(data, 10.0, 1.0, kernel_pixels=41)
    assert coverage[40, 20] == pytest.approx(1.0, abs=1e-4) and coverage[40, 70] == pytest.approx(0.0, abs=1e-4)
    assert coverage[40, 40] == pytest.approx(0.5, abs=0.06)
    assert smooth[40, 20] == pytest.approx(2.0, abs=1e-3)
    assert smooth[40, 40] / coverage[40, 40] == pytest.approx(2.0, rel=1e-6)   # edge value is diluted by coverage


def test_linear_fit_and_ratio_stats_recover_gain_offset():
    rng = np.random.default_rng(0)
    planck = rng.uniform(10, 30, (20, 20))
    model = 1.15 * planck - 4.7
    sel = np.ones(planck.shape, dtype=bool)
    fit = pc.linear_fit(model, planck, sel, n_boot=50)
    assert fit["gain"] == pytest.approx(1.15) and fit["offset_mjysr"] == pytest.approx(-4.7)
    assert fit["background_mjysr"] == pytest.approx(4.7 / 1.15)
    stats = pc.ratio_stats(model, planck - fit["background_mjysr"], sel)
    assert stats["median"] == pytest.approx(1.15) and stats["n"] == 400


def test_bright_region_and_its_three_ratios():
    yy, xx = np.indices((41, 41))
    blob = np.exp(-((yy - 20) ** 2 + (xx - 25) ** 2) / (2 * 6.0**2))
    planck = 20.0 * blob + 5.0                 # Planck: source plus a uniform background of 5
    model = 0.9 * 20.0 * blob                  # model: 10% low, no background
    allowed = np.ones(model.shape, dtype=bool)
    allowed[:, 30:] = False
    peak, region = pc.bright_region(model, allowed, 3.0)
    assert peak == (20, 25) and region[20, 25] and not region[20, 30] and not region[20, 21]
    reference = allowed & (planck <= np.quantile(planck[allowed], 0.2))
    out = pc.bright_region_calibration(model, planck, peak, region, reference, background_mjysr=5.0)
    r = out["ratio_model_over_planck"]["peak_pixel"]
    assert r["direct"] == pytest.approx(18.0 / 25.0)
    assert r["background_removed"] == pytest.approx(0.9)
    assert r["differential"] == pytest.approx(0.9)
    assert out["region_mean"]["n"] == int(region.sum())
