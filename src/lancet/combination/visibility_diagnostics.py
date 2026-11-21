"""Read-only measurement-set access and per-row visibility averaging with explicit flags."""

from __future__ import annotations

from pathlib import Path

import numpy as np


C_METRES_PER_SECOND = 299792458.0


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
