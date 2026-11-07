"""Full-range amplitude profiles with a controlled overlap comparison."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .combination_tests import (
    _broad_ddids, _concat_rows, _write_rows_csv, matching_cells,
    predict_primary_beams, read_test_rows, row_values, subset_rows,
)


def _cell_statistics(rows, keys, valid, size, kappa, models=None):
    value, weight, radius = row_values(rows)
    histogram = lambda x: np.bincount(keys[valid], weights=x[valid], minlength=size)
    total = histogram(weight)
    mean = lambda x: np.divide(histogram(weight*x), total, out=np.full(size, np.nan), where=total > 0)
    result = {"count": histogram(np.ones(len(value))), "radius": mean(radius),
              "power": mean(np.abs(value)**2 - 2*kappa/weight),
              "cross_power": mean(np.real(rows["x"]*rows["y"].conj()))}
    power = np.real(rows["x"]*rows["y"].conj())
    w2, w2p, w2p2 = [histogram(weight**2*x) for x in (np.ones(len(value)), power, power**2)]
    variance_sum = w2p2 - 2*result["cross_power"]*w2p + result["cross_power"]**2*w2
    # Variance of the weighted cell mean, with effective-count correction.
    denominator = total**2-w2
    result["cross_variance"] = np.divide(np.maximum(variance_sum, 0), denominator,
        out=np.full(size, np.nan), where=denominator > 0)
    for key, model in (models or {}).items():
        result[key] = mean(np.abs(model)**2)
    return result


def combine_sampling_cells(aca, alma, n_pairs, n_bins, sectors):
    """Match occupied cells where both arrays sample; retain unique outer modes.

    Natural weights have already been used inside cells. Cell means receive
    equal weight, as in the completed overlap test. No PB factor is fitted to
    measured amplitudes. Signed noise-subtracted powers are not clipped.
    """
    size = n_pairs*2*n_bins*sectors
    radial = (np.arange(size)//sectors) % n_bins
    result = []
    sqrt_positive = lambda x: float(np.sqrt(x)) if x > 0 else np.nan
    for i in range(n_bins):
        a = (radial == i) & (aca["count"] >= 4)
        b = (radial == i) & (alma["count"] >= 4)
        common = a & b
        matched = bool(common.any())
        if matched:
            a = b = common
        record = {"bin": i, "matched_overlap": matched, "common_cells": int(common.sum())}
        for label, data, use in [("aca", aca, a), ("alma", alma, b)]:
            record[label+"_cells"] = int(use.sum())
            record[label+"_rows"] = int(data["count"][use].sum())
            record[label+"_uv_klambda"] = float(np.mean(data["radius"][use])) if use.any() else np.nan
            auto = float(np.mean(data["power"][use])) if use.any() else np.nan
            cross = float(np.mean(data["cross_power"][use])) if use.any() else np.nan
            variance = data.get("cross_variance", np.zeros_like(data["count"]))
            error = float(np.sqrt(np.sum(variance[use]))/use.sum()) if use.any() else np.nan
            detected = bool(np.isfinite(cross) and cross > 2*error)
            record[label+"_signal_power_jy2"] = cross
            record[label+"_auto_noise_subtracted_power_jy2"] = auto
            record[label+"_rms_mjy"] = 1000*sqrt_positive(cross)
            record[label+"_cross_rms_mjy"] = 1000*sqrt_positive(cross)
            record[label+"_power_standard_error_jy2"] = error
            record[label+"_detected"] = detected
            record[label+"_upper_rms_mjy"] = 1000*sqrt_positive(max(cross, 0)+2*error) if use.any() else np.nan
        for model in ("allscale", "twelve_m"):
            factor = sqrt_positive(np.mean(aca[model+"_target"][a])/np.mean(aca[model+"_native"][a])) if a.any() else np.nan
            record[model+"_aca_to_12m_factor"] = factor
            record[model+"_aca_12m_pb_mjy"] = record["aca_rms_mjy"]*factor
        result.append(record)
    return result


def separate_array_profiles(statistics, n_bins, sectors):
    """Compare arrays on common cells where possible, retaining unique modes.

    Common cells must have at least four rows in every array present in that
    radial bin. If there are no common cells, each array retains its own cells
    and the CSV explicitly records that the bin is not matched.
    """
    size = len(next(iter(statistics.values()))["count"])
    radial = (np.arange(size)//sectors) % n_bins
    records = []
    for i in range(n_bins):
        occupied = {name: (radial == i) & (s["count"] >= 4) for name, s in statistics.items()}
        present = [use for use in occupied.values() if use.any()]
        common = np.logical_and.reduce(present) if present else np.zeros(size, dtype=bool)
        matched = len(present) > 1 and common.any()
        for name, data in statistics.items():
            use = common if matched and occupied[name].any() else occupied[name]
            if not use.any():
                continue
            power = float(np.mean(data["cross_power"][use]))
            error = float(np.sqrt(np.sum(data["cross_variance"][use]))/use.sum())
            factor = 1.0
            if name == "aca_7m":
                factor = float(np.sqrt(np.mean(data["allscale_target"][use])/np.mean(data["allscale_native"][use])))
            records.append({"array": name, "bin": i, "matched_arrays": bool(matched),
                "cells": int(use.sum()), "rows": int(data["count"][use].sum()),
                "uv_klambda": float(np.mean(data["radius"][use])),
                "signal_power_jy2": power, "power_standard_error_jy2": error,
                "pb_amplitude_factor": factor, "detected": bool(power > 2*error),
                "rms_mjy": 1000*np.sqrt(power)*factor if power > 0 else np.nan,
                "upper_rms_mjy": 1000*np.sqrt(max(power, 0)+2*error)*factor})
    return records


def build_separate_array_profiles(paths, coverage):
    """Split C2/C5 from existing full-range caches; do not reread the MSs."""
    root = Path(paths["output_root"])
    out = root/"full_sampling"
    summary = json.loads((root/"overlap_tests"/"summary.json").read_text())
    edges = np.array(coverage["uv_edges_klambda"])
    sectors = int(coverage["settings"]["orientation_sectors"])
    nb = len(edges)-1
    size = len(coverage["target_pointings"])*2*nb*sectors
    with np.load(out/"cell_statistics.npz") as z:
        statistics = {"aca_7m": {k.removeprefix("aca_"): z[k] for k in z.files if k.startswith("aca_")}}
    for name in ("configuration_c2", "configuration_c5"):
        with np.load(out/f"{name}_rows.npz") as z:
            rows = {k: z[k] for k in z.files}
        key, good = matching_cells(rows, coverage["target_pointings"], "alma", edges, sectors)
        kappa = np.array([summary["noise"][name][str(int(d))]["kappa_time"] for d in rows["ddid"]])
        statistics[name] = _cell_statistics(rows, key, good, size, kappa)
    records = separate_array_profiles(statistics, nb, sectors)
    _write_rows_csv(out/"separate_array_amplitudes.csv", records)
    return records


def build_full_sampling_profiles(paths, config, *, overwrite=False):
    """Keep ACA's shorter modes and the entire observed 12 m baseline range."""
    root = Path(paths["output_root"])
    out = root/"full_sampling"
    out.mkdir(parents=True, exist_ok=True)
    if (out/"amplitudes.csv").exists() and not overwrite:
        raise FileExistsError("Full sampling profiles exist; use --overwrite")
    test = root/"overlap_tests"
    summary = json.loads((test/"summary.json").read_text())
    settings = config["overlap_tests"]
    edges = np.array(config["figure1"]["amplitude_uv_edges_klambda"])
    fields = sorted(p["alma_field"] for p in summary["matched_pointings"])
    rows, meta = {}, {}
    for name in ("configuration_c5", "configuration_c2"):
        cache = out/f"{name}_rows.npz"
        metapath = out/f"{name}_metadata.json"
        if cache.is_file() and metapath.is_file():
            meta[name] = json.loads(metapath.read_text())
            if meta[name].get("test_field_ids") != fields or meta[name]["test_max_klambda"] != config["visibility_profile"]["uv_max_klambda"]:
                raise ValueError("Full-range row cache has a different pointing selection")
            with np.load(cache) as z:
                rows[name] = {k: z[k] for k in z.files}
        else:
            rows[name], meta[name] = read_test_rows(Path(paths["measurement_sets"][name]),
                float(config["visibility_profile"]["uv_max_klambda"]), field_ids=fields)
            np.savez_compressed(cache, **rows[name])
            metapath.write_text(json.dumps(meta[name], indent=2)+"\n")
    with np.load(test/"aca_7m_rows.npz") as z:
        aca = {k: z[k] for k in z.files}
    meta["aca_7m"] = json.loads((test/"aca_7m_metadata.json").read_text())
    use = np.isin(aca["ddid"], _broad_ddids(meta["aca_7m"])) & np.isin(
        aca["field"], [p["aca_field"] for p in summary["matched_pointings"]])
    aca = subset_rows(aca, use)
    alma = _concat_rows([rows["configuration_c5"], rows["configuration_c2"]])
    noise = lambda name, values: np.array([summary["noise"][name][str(int(d))]["kappa_time"] for d in values["ddid"]])
    ka = noise("aca_7m", aca)
    kb = np.r_[noise("configuration_c5", rows["configuration_c5"]), noise("configuration_c2", rows["configuration_c2"])]
    # Use the actual matched 12 m centre as the target aperture, not a midpoint.
    targets = []
    for p in summary["matched_pointings"]:
        ra, dec = np.rad2deg(np.asarray(meta["configuration_c5"]["subtables"]["FIELD"][p["alma_field"]]["PHASE_DIR"]).ravel())
        targets.append({**p, "ra_deg": float(ra % 360), "dec_deg": float(dec)})
    models = {}
    for label, image in [("allscale", "allscale_sky"), ("twelve_m", "twelve_m_pbcor")]:
        native, target, _ = predict_primary_beams(aca, meta["aca_7m"], targets,
            "aca", Path(paths["images"][image]), settings)
        models[label+"_native"], models[label+"_target"] = native, target
    sectors = int(settings["orientation_sectors"])
    ca, ga = matching_cells(aca, targets, "aca", edges, sectors)
    cb, gb = matching_cells(alma, targets, "alma", edges, sectors)
    nb = len(edges)-1
    size = len(targets)*2*nb*sectors
    astat = _cell_statistics(aca, ca, ga, size, ka, models)
    bstat = _cell_statistics(alma, cb, gb, size, kb)
    records = combine_sampling_cells(astat, bstat, len(targets), nb, sectors)
    for i, r in enumerate(records):
        r["uv_low_klambda"], r["uv_high_klambda"] = float(edges[i]), float(edges[i+1])
    _write_rows_csv(out/"amplitudes.csv", records)
    ranges = {}
    for name, data in [("aca", aca), ("alma", alma)]:
        _, _, radius = row_values(data)
        ranges[name] = [float(radius.min()), float(radius.max())]
        if np.any((radius < edges[0]) | (radius > edges[-1])):
            raise ValueError("Figure-1 radial bins omit usable full-range data")
    report = {"ranges_klambda": ranges,
              "overlap_klambda": [max(ranges["aca"][0], ranges["alma"][0]), min(ranges["aca"][1], ranges["alma"][1])],
              "target_pointings": targets, "settings": dict(settings), "uv_edges_klambda": edges.tolist(),
              "rows": {"aca": len(aca["x"]), "alma": len(alma["x"])},
              "interpretation": "Matched occupied pointing/sideband/orientation cells within overlap; own occupied cells elsewhere. RMS estimated from independent-hand cross-power throughout, with auto-power noise subtraction retained as a check. Same 19 pointing pairs throughout. 7m model translation targets the native matched 12m pointing PB. No extrapolated correction factors or fitted flux normalization.",
              "long_baseline_limits": "Two standard errors from weighted within-cell cross-power scatter, propagated over the fixed cell population. For non-detections, conservative amplitude limit sqrt(max(mean_power,0)+2*SE). Does not include calibration/systematic errors or full visibility covariance."}
    (out/"summary.json").write_text(json.dumps(report, indent=2)+"\n")
    np.savez_compressed(out/"cell_statistics.npz", **{"aca_"+k:v for k,v in astat.items()}, **{"alma_"+k:v for k,v in bstat.items()})
    return report
