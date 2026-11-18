import numpy as np


def _written(stem):
    return all(stem.with_suffix(f".{ext}").stat().st_size > 0 for ext in ("pdf", "png"))


def test_plot_tp_factor_writes_both_formats(tmp_path):
    from lancet.combination_figures import plot_tp_factor

    centres = np.geomspace(0.6, 28.0, 20)
    p_int = centres**-3.0
    stem = tmp_path / "tp_factor" / "tp_factor_dust_ff"
    plot_tp_factor(stem, centres, p_int, p_int / 1.96, [4.1, 10.3], 1.4, 1.41)
    assert _written(stem)


def test_plot_comparison_writes_both_formats(tmp_path):
    from lancet.combination_figures import plot_comparison

    image = np.random.default_rng(2).random((60, 80)) * 1e-3
    coverage = np.zeros(image.shape, dtype=bool)
    coverage[20:40, 10:70] = True
    stem = tmp_path / "jcomb" / "comparison_dust_ff"
    plot_comparison(stem, image, 1.4 * image, 2.4 * image, coverage, "J-comb", "TP")
    assert _written(stem)


def test_plot_power_spectrum_accepts_missing_tp_ratios(tmp_path):
    from lancet.combination_figures import plot_power_spectrum

    centres = np.array([1.4, 2.4, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11.0, 13.4, 17.3, 24.5, 38.7])
    p = centres**-3.0
    to_tp = [1.0 if c < 9 else None for c in centres]           # the TP term ends at k2
    stem = tmp_path / "feather" / "power_spectrum_dust"
    plot_power_spectrum(stem, centres, {"interferometer": p, "tp": 2 * p, "combined": 1.5 * p},
                        {"interferometer": np.ones_like(centres), "tp": to_tp}, 9.0, (), "Feather", "TP × 1.4")
    assert _written(stem)





def test_plot_zero_spacing_writes_both_formats(tmp_path):
    from lancet.combination_figures import plot_zero_spacing

    tp = np.random.default_rng(5).random((40, 50)) * 50 + 10
    mosaic = np.zeros(tp.shape, dtype=bool)
    mosaic[15:25, 10:40] = True
    tp_only = ~mosaic
    tp_only[:3] = False
    stats = {"mosaic": {"flux_ratio": 1.001}, "tp_only": {"flux_ratio": 0.999}}
    stem = tmp_path / "jcomb" / "zero_spacing_dust_ff"
    plot_zero_spacing(stem, tp, 1.001 * tp, mosaic, tp_only, "J-comb", "TP dust+free-free", stats)
    assert _written(stem)
