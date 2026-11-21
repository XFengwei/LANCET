import numpy as np
from astropy.io import fits

from lancet.combination import aca_forward_test
from lancet.combination.aca_forward_test import binned_rms, pair_weights, sky_model


def test_bootstrap_weights_repeat_pairs():
    assert np.array_equal(pair_weights(3, np.array([0, 0, 2])), [2, 0, 1])
    power = np.array([1.0, 4.0, 9.0])
    use = np.ones(3, bool)
    assert np.isclose(binned_rms(power, use, np.arange(3), pair_weights(3)), np.sqrt(14/3))
    assert np.isclose(binned_rms(power, use, np.arange(3), np.array([2.0, 0, 1])), np.sqrt(11/3))


def test_jy_per_beam_is_converted_to_mjy_per_sr(tmp_path):
    header = fits.Header({"BUNIT": "Jy/beam", "BMAJ": 1/3600, "BMIN": 1/3600, "CDELT1": -0.1/3600,
                          "CDELT2": 0.1/3600, "CRPIX1": 5.5, "CRPIX2": 5.5, "CTYPE1": "RA---SIN",
                          "CTYPE2": "DEC--SIN", "CRVAL1": 0.0, "CRVAL2": 0.0})
    path = tmp_path/"image.fits"
    fits.PrimaryHDU(np.ones((10, 10), np.float32), header).writeto(path)
    data, out = sky_model(path, 0.5)
    area = np.pi/(4*np.log(2))*(np.pi/648000)**2
    assert data.shape == (2, 2) and out["BUNIT"] == "MJy/sr"
    assert np.allclose(data, 1e-6/area, rtol=1e-6) and np.isclose(out["CDELT2"]*3600, 0.5)


def test_gaussian_transfer_and_limit_are_consistent():
    from lancet.combination.amplitude_alignment import gaussian_transfer, response_limit
    limit = response_limit(18.1, 0.1)
    assert np.isclose(gaussian_transfer(limit, 18.1), 0.1)
    half = response_limit(18.1, 0.5)
    assert np.isclose(half, 5.0287, rtol=1e-3)  # matches the TP summary half-response radius


def test_no_tp_scale_is_derived_from_the_forward_test():
    # The 12m+7m image is not the f = 0 member of J-comb or feather, so no TP scale is
    # interpolated between them; the stage reports only the predicted / measured ratios.
    assert not hasattr(aca_forward_test, "tp_scales_reproducing_aca")
    assert not hasattr(aca_forward_test, "scale_for_unit_ratio")
    import inspect
    assert "tp_scale" not in inspect.getsource(aca_forward_test).replace("tp_scale_given", "")


def test_signal_panel_draws_the_flagged_models(tmp_path):
    from matplotlib import pyplot as plt
    models = {"interferometer": {"label": "12m+7m", "signal_panel": True},
              "jcomb": {"label": "J-comb", "signal_panel": True}, "aca_only": {"label": "ACA-only image"}}
    test = {"sky_models": models, "figure_selection": "interior_pairs", "output_formats": ["png"]}
    profiles = []
    for uv in (6.0, 8.0):
        record = {"selection": "interior_pairs", "sideband": "both", "uv_klambda": uv, "aca_rms_mjy": 10.0,
                  "tp_rms_mjy": 9.0, "tp_to_aca": 0.9}
        for model in models:
            record.update({model+"_rms_mjy": 8.0, model+"_to_aca": 0.8, model+"_to_aca_p16": 0.7,
                           model+"_to_aca_p84": 0.9})
        profiles.append(record)
    figures = []
    aca_forward_test._plot(profiles, test, tmp_path, lambda fig, outputs, **_: figures.append(fig), plt)
    signal = figures[0].axes[1]
    assert len(signal.lines) == 2 + 2          # measured ACA, TP model, and the two flagged models
