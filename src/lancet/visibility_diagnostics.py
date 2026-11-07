"""Read-only continuum visibility profiles with explicit flags and averaging."""

from __future__ import annotations

from pathlib import Path

import numpy as np


C_METRES_PER_SECOND = 299792458.0
PARALLEL_HANDS = {1, 5, 8, 9, 12}  # I, RR, LL, XX, YY (MS correlation codes)


def average_visibility_rows(data, flags, weights, frequencies_hz):
    """Vector-average selected correlations/channels within each MS row.

    Inputs have CASA ordering (correlation, channel, row). Weights are either
    per-channel or broadcast from WEIGHT. A row is retained only if at least
    one unflagged, finite datum has positive finite weight. No averaging is
    performed across fields, times, baselines, or spectral windows here.
    """
    data = np.asarray(data)
    flags = np.asarray(flags, dtype=bool)
    weights = np.broadcast_to(weights, data.shape)
    frequency = np.asarray(frequencies_hz)[None, :, None]
    valid = (~flags & np.isfinite(data) & np.isfinite(weights)
             & (weights > 0) & np.isfinite(frequency) & (frequency > 0))
    w = np.where(valid, weights, 0.0).astype(float)
    total = w.sum(axis=(0, 1))
    visibility = np.full(total.shape, np.nan + 1j * np.nan)
    effective_frequency = np.full(total.shape, np.nan)
    np.divide((np.where(valid, data, 0.0) * w).sum(axis=(0, 1)),
              total, out=visibility, where=total > 0)
    np.divide((w * frequency).sum(axis=(0, 1)), total,
              out=effective_frequency, where=total > 0)
    return visibility, total, effective_frequency, valid.sum(axis=(0, 1))


def radial_visibility_sums(uv_klambda, visibility, weight, edges):
    """Accumulate additive row statistics; never average amplitudes coherently."""
    radius = np.asarray(uv_klambda)
    values = np.asarray(visibility)
    weight = np.asarray(weight)
    good = (np.isfinite(radius) & np.isfinite(values) & np.isfinite(weight)
            & (weight > 0))
    amplitude = np.abs(values[good])
    r, w = radius[good], weight[good]
    histogram = lambda x: np.histogram(r, bins=edges, weights=x)[0]
    return np.stack([
        histogram(np.ones(r.size)), histogram(w), histogram(w * amplitude),
        histogram(w * amplitude**2), histogram(w * r),
    ])


def _table_tool():
    # Disable update/network activity and file logs before importing CASA.
    from casaconfig import config

    config.nologfile = True
    config.logfile = None
    config.data_auto_update = False
    config.measures_auto_update = False
    from casatools import table

    return table()


def measurement_set_metadata(path: Path):
    """Read data layout, pointing coordinates and channel setup, without writes."""
    tool = _table_tool()
    metadata = {"path": str(path), "subtables": {}}
    tool.open(str(path), nomodify=True)
    try:
        metadata["rows"] = int(tool.nrows())
        metadata["columns"] = list(tool.colnames())
        metadata["data_unit"] = tool.getcolkeywords("DATA").get("QuantumUnits", [])
    finally:
        tool.close()
    for name, columns in {
        "SPECTRAL_WINDOW": ["CHAN_FREQ", "CHAN_WIDTH", "NUM_CHAN"],
        "DATA_DESCRIPTION": ["SPECTRAL_WINDOW_ID", "POLARIZATION_ID"],
        "POLARIZATION": ["CORR_TYPE"],
        "FIELD": ["NAME", "PHASE_DIR"],
        "ANTENNA": ["DISH_DIAMETER"],
    }.items():
        tool.open(str(path / name), nomodify=True)
        try:
            metadata["subtables"][name] = [
                {c: np.asarray(tool.getcell(c, row)).tolist() for c in columns}
                for row in range(tool.nrows())
            ]
        finally:
            tool.close()
    return metadata


def read_visibility_profiles(path: Path, edges, *, chunk_rows=20000):
    """Stream all DDIDs in bounded chunks, respecting variable channel shapes.

    DATA is used deliberately: these delivered self-calibrated MSs have no
    CORRECTED_DATA column. Prefer WEIGHT_SPECTRUM when defined, otherwise
    broadcast WEIGHT across channels, as specified by the MS convention.
    """
    metadata = measurement_set_metadata(path)
    subtables = metadata["subtables"]
    nfields = len(subtables["FIELD"])
    sums = np.zeros((5, len(edges) - 1))
    field_sums = np.zeros((nfields, *sums.shape))
    ddid_sums = []
    counts = {"total_rows": metadata["rows"], "autocorrelation_rows": 0,
              "flagged_cross_rows": 0, "eligible_cross_rows": 0,
              "eligible_parallel_samples": 0, "usable_parallel_samples": 0,
              "usable_averaged_rows": 0, "rows_outside_bins": 0}
    extrema = {"uv_metres": [np.inf, -np.inf],
               "uv_klambda": [np.inf, -np.inf], "frequency_hz": [np.inf, -np.inf]}
    tool = _table_tool()
    tool.open(str(path), nomodify=True)
    weight_columns = set()
    try:
        counts["autocorrelation_rows"] = _query_count(tool, "ANTENNA1 == ANTENNA2")
        counts["flagged_cross_rows"] = _query_count(
            tool, "ANTENNA1 != ANTENNA2 && FLAG_ROW")
        for ddid, desc in enumerate(subtables["DATA_DESCRIPTION"]):
            spw = subtables["SPECTRAL_WINDOW"][desc["SPECTRAL_WINDOW_ID"]]
            freq = np.array(spw["CHAN_FREQ"])
            corr = subtables["POLARIZATION"][desc["POLARIZATION_ID"]]["CORR_TYPE"]
            selected = np.flatnonzero([c in PARALLEL_HANDS for c in corr])
            if not selected.size:
                raise ValueError(f"No parallel-hand correlations in DDID {ddid}")
            group = tool.query(
                f"DATA_DESC_ID == {ddid} && ANTENNA1 != ANTENNA2 && !FLAG_ROW")
            dd_sums = np.zeros_like(sums)
            try:
                weight_column = "WEIGHT"
                if ("WEIGHT_SPECTRUM" in metadata["columns"] and group.nrows()
                        and group.iscelldefined("WEIGHT_SPECTRUM", 0)):
                    weight_column = "WEIGHT_SPECTRUM"
                weight_columns.add(weight_column)
                counts["eligible_cross_rows"] += int(group.nrows())
                for start in range(0, group.nrows(), chunk_rows):
                    count = min(chunk_rows, group.nrows() - start)
                    get = lambda name: np.asarray(group.getcol(name, startrow=start, nrow=count))
                    data = get("DATA")[selected]
                    flag = get("FLAG")[selected]
                    weight = get(weight_column)[selected]
                    if weight_column == "WEIGHT":
                        weight = weight[:, None, :]
                    vis, w, f, valid_samples = average_visibility_rows(data, flag, weight, freq)
                    counts["eligible_parallel_samples"] += int(data.size)
                    counts["usable_parallel_samples"] += int(valid_samples.sum())
                    uvw = get("UVW")
                    metres = np.hypot(uvw[0], uvw[1])
                    radius = metres * f / C_METRES_PER_SECOND / 1000
                    good = (w > 0) & np.isfinite(radius) & (radius > 0)
                    counts["usable_averaged_rows"] += int(good.sum())
                    counts["rows_outside_bins"] += int(
                        (good & ((radius < edges[0]) | (radius > edges[-1]))).sum())
                    if good.any():
                        for key, value in [("uv_metres", metres), ("uv_klambda", radius),
                                           ("frequency_hz", f)]:
                            extrema[key][0] = min(extrema[key][0], float(value[good].min()))
                            extrema[key][1] = max(extrema[key][1], float(value[good].max()))
                    radial = radial_visibility_sums(radius[good], vis[good], w[good], edges)
                    sums += radial
                    dd_sums += radial
                    fields = get("FIELD_ID")
                    for field in np.unique(fields[good]):
                        take = good & (fields == field)
                        field_sums[field] += radial_visibility_sums(
                            radius[take], vis[take], w[take], edges)
            finally:
                group.close()
            ddid_sums.append(dd_sums)
            print(f"  {path.name}: DDID {ddid} read", flush=True)
    finally:
        tool.close()
    metadata.update(counts=counts, ranges=extrema, weight_columns=sorted(weight_columns),
                    data_column="DATA", averaging="Weighted complex mean over parallel hands and channels per MS row; then scalar amplitude mean in radial bins")
    return sums, field_sums, np.array(ddid_sums), metadata


def _query_count(tool, expression):
    group = tool.query(expression)
    try:
        return int(group.nrows())
    finally:
        group.close()


def profile_columns(sums, edges):
    """Return radial centers, scalar amplitudes, spread, counts and weights."""
    count, weight, amp_sum, amp2_sum, radius_sum = sums
    divide = lambda a, b: np.divide(a, b, out=np.full_like(a, np.nan), where=b > 0)
    amplitude = divide(amp_sum, weight)
    return {
        "uv_low_klambda": np.asarray(edges[:-1]),
        "uv_high_klambda": np.asarray(edges[1:]),
        "uv_mean_klambda": divide(radius_sum, weight),
        "mean_amplitude_jy": amplitude,
        "amplitude_std_jy": np.sqrt(np.maximum(0, divide(amp2_sum, weight) - amplitude**2)),
        "averaged_row_count": count,
        "sum_weight": weight,
        "mean_weight": divide(weight, count),
    }
