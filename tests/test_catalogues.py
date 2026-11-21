"""Tests for source-catalogue utilities and task integration."""

from pathlib import Path
import json
import subprocess
import sys

from astropy.table import Table
import numpy as np
import pytest

from lancet.catalogues import (
    effective_source_sizes,
    gain_correct,
    getsf_peak_indices,
    peak_to_total_flux_ratio,
    relative_flux_difference,
    remove_source_model,
    select_science_sources,
    uniform_brightness_flux_ratio,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_getsf_peak_index_rule_preserves_central_axis_offset() -> None:
    x, y = getsf_peak_indices(
        np.array([2.9, 7.1]),
        np.array([1.2, 8.2]),
        shape=(10, 10),
    )

    np.testing.assert_array_equal(x, [2, 8])
    np.testing.assert_array_equal(y, [1, 9])


def test_effective_sizes_apply_resolved_and_unresolved_rules() -> None:
    result = effective_source_sizes(
        np.array([2.0, 0.5]),
        np.array([1.0, 0.5]),
        beam_major_arcsec=0.71,
        beam_minor_arcsec=0.60,
        distance_pc=2800.0,
    )

    assert result["resolved"].tolist() == ["Y", "N"]
    assert result["major_arcsec"][0] == pytest.approx(
        np.sqrt(2.0 - 0.71 * 0.60)
    )
    assert result["major_arcsec"][1] == pytest.approx(
        0.5 * np.sqrt(0.71 * 0.60)
    )


def test_science_selection_uses_either_measurement_scale() -> None:
    table = Table(
        {
            "FXT_BST01": [1.0, -1.0],
            "FXT_BST02": [-1.0, 1.0],
            "FXP_BST01": [1.0, -1.0],
            "FXP_BST02": [-1.0, 1.0],
            "SIGNM01": [5.0, 1.0],
            "SIGNM02": [1.0, 5.0],
            "AFWHM01": [2.0, 4.0],
            "BFWHM01": [2.0, 4.0],
            "AFWHM02": [4.0, 2.0],
            "BFWHM02": [4.0, 2.0],
        }
    )

    assert len(select_science_sources(table)) == 2


def test_peak_to_total_flux_ratio_propagates_both_errors() -> None:
    ratio, uncertainty = peak_to_total_flux_ratio(
        np.array([2.0]),
        np.array([0.2]),
        np.array([4.0]),
        np.array([0.4]),
    )

    np.testing.assert_allclose(ratio, [0.5])
    np.testing.assert_allclose(
        uncertainty,
        [np.sqrt((0.2 / 4.0) ** 2 + (0.4 * 2.0 / 4.0**2) ** 2)],
    )


def test_relative_flux_difference_uses_12m7m_as_reference() -> None:
    difference = relative_flux_difference(
        np.array([10.0, 20.0]),
        np.array([12.0, 18.0]),
    )

    np.testing.assert_allclose(difference, [0.2, -0.1])


def test_uniform_brightness_ratio_is_one_at_beam_scale() -> None:
    ratio = uniform_brightness_flux_ratio(
        np.array([0.65, 1.30]),
        beam_fwhm_arcsec=0.65,
    )

    np.testing.assert_allclose(ratio, [1.0, 0.25])


def test_gain_correct_divides_only_below_unit_gain() -> None:
    values = np.array([2.0, 2.0, 2.0, 2.0, 2.0])
    gain = np.array([1.0, 0.5, 0.8, 1.0000001, np.nan])

    corrected = gain_correct(values, gain)

    np.testing.assert_allclose(corrected[:4], [2.0, 4.0, 2.5, 2.0])
    assert np.isnan(corrected[4])


def test_gain_correct_rejects_nonpositive_gain() -> None:
    with pytest.raises(ValueError, match="gain"):
        gain_correct(np.array([1.0]), np.array([0.0]))


def test_remove_source_model_uses_true_sky_source_model() -> None:
    image = np.array([[10.0, 10.0], [10.0, np.nan]])
    sources = np.array([[2.0, 0.0], [1.0, 0.0]])
    gain = np.array([[1.0, 0.5], [0.5, np.nan]])

    background = remove_source_model(image, sources, gain)

    np.testing.assert_allclose(background[0], [8.0, 10.0])
    assert background[1, 0] == pytest.approx(8.0)
    assert np.isnan(background[1, 1])


def test_remove_source_model_needs_gain_where_sources_exist() -> None:
    background = remove_source_model(
        np.array([5.0, 5.0]),
        np.array([1.0, 0.0]),
        np.array([np.nan, np.nan]),
    )

    assert np.isnan(background[0])
    assert background[1] == 5.0


@pytest.mark.local_inputs("data", "results")
def test_source_catalogue_inspection_is_read_only() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "tasks/source_catalogue/run.py",
            "--stage",
            "inspect",
        ],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "rows=234, rows after science cuts=231" in completed.stdout
    assert "results/source_catalogue" in completed.stdout


def test_generated_source_catalogue_schema_and_regression_values() -> None:
    generated = (
        PROJECT_ROOT
        / "results"
        / "source_catalogue"
        / "tables"
        / "source_catalogue.ecsv"
    )
    if not generated.is_file():
        pytest.skip("Generated source catalogue is unavailable")
    actual = Table.read(generated, format="ascii.ecsv")

    assert actual.colnames == [
        "ID",
        "NO",
        "XCO_P",
        "YCO_P",
        "WCS_ACOOR",
        "WCS_DCOOR",
        "RA(hms)",
        "DEC(dms)",
        "SIGN",
        "GOOD",
        "SIGNM01",
        "GOODM01",
        "SIGNM02",
        "GOODM02",
        "FXP_BST01",
        "FXP_ERR01",
        "FXT_BST01",
        "FXT_ERR01",
        "FXP_BST02",
        "FXP_ERR02",
        "FXT_BST02",
        "FXT_ERR02",
        "AFWHM01",
        "BFWHM01",
        "THETA01",
        "AFWHM02",
        "BFWHM02",
        "THETA02",
        "Ibkg1",
        "Ibkg2",
        "Tenv",
        "Gain",
        "FWHM_OBS_str",
        "Deconv_maj",
        "Deconv_min",
        "Deconv_theta",
        "Deconv_FWHM_angle",
        "Deconv_FWHM_linear",
        "Resolvability",
        "FWHM_SRC_str",
        "FXP_WITH_ERROR",
        "FXT_WITH_ERROR",
        "FXP_BST01_RAW",
        "FXP_ERR01_RAW",
        "FXT_BST01_RAW",
        "FXT_ERR01_RAW",
        "Ibkg1_RAW",
        "Subregion",
    ]
    gain = np.asarray(actual["Gain"], dtype=float)
    assert np.count_nonzero(gain < 1) == 19
    np.testing.assert_allclose(
        actual["FXT_BST01"],
        actual["FXT_BST01_RAW"] / np.minimum(gain, 1.0),
    )
    assert len(actual) == 231
    indices = [0, 115, 230]
    np.testing.assert_array_equal(actual["ID"][indices], [1, 116, 231])
    np.testing.assert_array_equal(actual["NO"][indices], [1, 103, 213])
    np.testing.assert_allclose(
        actual["FXT_BST01"][indices],
        [0.1211, 0.003502, 0.00208],
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        actual["Deconv_FWHM_linear"][indices],
        [1335.7144605041904, 2826.786036473224, 2564.636457668026],
        rtol=0,
        atol=1e-10,
    )
    np.testing.assert_array_equal(
        np.asarray(actual["Subregion"][indices], dtype=str),
        ["E", "I", "I"],
    )


def test_generated_flux_diagnostics_are_complete() -> None:
    root = PROJECT_ROOT / "results" / "source_catalogue"
    table_path = root / "tables" / "flux_diagnostics.ecsv"
    summary_path = root / "tables" / "flux_diagnostics_summary.json"
    figures = [
        root / "figures" / "GetsfResult_FluxMeasure.pdf",
        root / "figures" / "GetsfResult_FluxMeasure.png",
        root / "figures" / "GetsfResult_FluxRatio.pdf",
        root / "figures" / "GetsfResult_FluxRatio.png",
    ]
    if not table_path.is_file() or not summary_path.is_file():
        pytest.skip("Flux-diagnostic products are unavailable")

    table = Table.read(table_path, format="ascii.ecsv")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    assert len(table) == summary["science_source_count"] == 231
    assert summary["relative_difference_source_count"] == 169
    assert summary["central_source_count"] == 154
    assert summary["display_sigma_percent"] == 10
    assert {
        "Peak_to_integrated_12m7m",
        "Peak_to_integrated_12m7mTP",
        "Relative_integrated_flux_difference",
        "Peak_SNR_both",
        "Central_relative_difference",
    } <= set(table.colnames)
    for figure in figures:
        assert figure.is_file()
        assert figure.stat().st_size > 0
