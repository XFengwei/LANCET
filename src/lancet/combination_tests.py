"""Controlled noise, pointing and primary-beam tests for continuum overlap."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .visibility_diagnostics import (
    C_METRES_PER_SECOND, _table_tool, average_visibility_rows,
    measurement_set_metadata,
)


ARRAYS = {
    "configuration_c5": ("12 m C5", "#7254A3"),
    "configuration_c2": ("12 m C2", "#007A9E"),
    "aca_7m": ("ACA 7 m", "#D87821"),
}


def read_test_rows(path, maximum_klambda, chunk_rows=20000, field_ids=None):
    """Read parallel hands separately, preserving row identifiers for time pairs."""
    meta = measurement_set_metadata(path)
    parts = {k: [] for k in ["x", "y", "wx", "wy", "frequency", "uvw",
                              "field", "ddid", "time", "scan", "observation",
                              "antenna1", "antenna2"]}
    t = _table_tool()
    t.open(str(path), nomodify=True)
    selected_rows = 0
    try:
        for ddid, desc in enumerate(meta["subtables"]["DATA_DESCRIPTION"]):
            corr = meta["subtables"]["POLARIZATION"][desc["POLARIZATION_ID"]]["CORR_TYPE"]
            if 9 not in corr or 12 not in corr:
                raise ValueError("XX and YY are required for the independent-hand control")
            frequencies = np.asarray(meta["subtables"]["SPECTRAL_WINDOW"][desc["SPECTRAL_WINDOW_ID"]]["CHAN_FREQ"])
            selection = f"DATA_DESC_ID == {ddid} && ANTENNA1 != ANTENNA2 && !FLAG_ROW"
            if field_ids is not None:
                selection += " && FIELD_ID IN [" + ",".join(str(int(f)) for f in field_ids) + "]"
            g = t.query(selection)
            try:
                weight_column = "WEIGHT_SPECTRUM" if (
                    "WEIGHT_SPECTRUM" in meta["columns"] and g.nrows()
                    and g.iscelldefined("WEIGHT_SPECTRUM", 0)) else "WEIGHT"
                for start in range(0, g.nrows(), chunk_rows):
                    n = min(chunk_rows, g.nrows() - start)
                    get = lambda name: np.asarray(g.getcol(name, startrow=start, nrow=n))
                    data, flags, weights = get("DATA"), get("FLAG"), get(weight_column)
                    if weight_column == "WEIGHT":
                        weights = weights[:, None, :]
                    hands = []
                    for c in (9, 12):
                        i = corr.index(c)
                        hands.append(average_visibility_rows(
                            data[i:i+1], flags[i:i+1], weights[i:i+1], frequencies))
                    x, wx, fx, _ = hands[0]
                    y, wy, fy, _ = hands[1]
                    frequency = np.divide(wx * fx + wy * fy, wx + wy,
                                          out=np.full(n, np.nan), where=(wx + wy) > 0)
                    uvw = get("UVW").T
                    radius = np.hypot(uvw[:, 0], uvw[:, 1]) * frequency / C_METRES_PER_SECOND / 1000
                    use = ((wx > 0) & (wy > 0) & np.isfinite(x) & np.isfinite(y)
                           & np.isfinite(radius) & (radius > 0) & (radius <= maximum_klambda))
                    selected_rows += int(use.sum())
                    values = {"x": x, "y": y, "wx": wx, "wy": wy,
                              "frequency": frequency, "uvw": uvw,
                              "ddid": np.full(n, ddid, dtype=np.int16)}
                    for key, column in [("field", "FIELD_ID"), ("time", "TIME"),
                                        ("scan", "SCAN_NUMBER"), ("observation", "OBSERVATION_ID"),
                                        ("antenna1", "ANTENNA1"), ("antenna2", "ANTENNA2")]:
                        values[key] = get(column)
                    for key in parts:
                        parts[key].append(values[key][use])
            finally:
                g.close()
            print(f"  Cached {path.name} DDID {ddid}", flush=True)
    finally:
        t.close()
    meta["test_rows_with_both_hands"] = selected_rows
    meta["test_max_klambda"] = maximum_klambda
    if field_ids is not None:
        meta["test_field_ids"] = sorted(int(f) for f in field_ids)
    return {k: np.concatenate(v) for k, v in parts.items()}, meta


def row_values(rows):
    weight = rows["wx"] + rows["wy"]
    value = (rows["wx"] * rows["x"] + rows["wy"] * rows["y"]) / weight
    radius = np.hypot(rows["uvw"][:, 0], rows["uvw"][:, 1]) * rows["frequency"] / C_METRES_PER_SECOND / 1000
    return value, weight, radius


def adjacent_time_pairs(rows, maximum_gap):
    """Non-overlapping adjacent pairs within observation/scan/field/DDID/baseline."""
    group_keys = ["observation", "scan", "field", "ddid", "antenna1", "antenna2"]
    order = np.lexsort(tuple([rows["time"]] + [rows[k] for k in reversed(group_keys)]))
    change = np.ones(len(order), dtype=bool)
    change[1:] = False
    for key in group_keys:
        change[1:] |= rows[key][order[1:]] != rows[key][order[:-1]]
    position = np.arange(len(order))
    starts = np.maximum.accumulate(np.where(change, position, 0))
    left = position[:-1][((position[:-1] - starts[:-1]) % 2 == 0) & ~change[1:]]
    i, j = order[left], order[left + 1]
    gap = rows["time"][j] - rows["time"][i]
    keep = (gap > 0) & (gap <= maximum_gap)
    return i[keep], j[keep]


def empirical_noise(rows, maximum_gap):
    """Estimate component-variance/nominal-variance from short time differences.

    For z=(Vi-Vj)/sqrt(1/wi+1/wj), a calibrated inverse-component-variance
    weight gives unit real/imaginary variance. Robust MAD limits sky drift and
    outlier sensitivity. XX-YY gives a second check, potentially containing Q.
    """
    value, weight, _ = row_values(rows)
    i, j = adjacent_time_pairs(rows, maximum_gap)
    result = {}
    kappa = np.full(len(value), np.nan)
    robust_variance = lambda z: float(np.mean([
        (np.median(np.abs(c - np.median(c))) / 0.6744897501960817)**2
        for c in (z.real, z.imag)]))
    for ddid in np.unique(rows["ddid"]):
        use = rows["ddid"][i] == ddid
        if use.sum() < 30:
            raise ValueError(f"Too few time pairs to estimate noise for DDID {ddid}")
        z = (value[i[use]] - value[j[use]]) / np.sqrt(1 / weight[i[use]] + 1 / weight[j[use]])
        same = rows["ddid"] == ddid
        pol = (rows["x"][same] - rows["y"][same]) / np.sqrt(1 / rows["wx"][same] + 1 / rows["wy"][same])
        estimate = robust_variance(z)
        result[str(int(ddid))] = {
            "time_pair_count": int(use.sum()), "kappa_time": estimate,
            "kappa_hand_difference": robust_variance(pol),
            "median_pair_gap_seconds": float(np.median(rows["time"][j[use]] - rows["time"][i[use]])),
        }
        kappa[same] = estimate
    return kappa, result, (i, j)


def radial_test_profile(rows, edges, kappa):
    """Keep signed debiased powers until after averaging; RMS is a new statistic."""
    value, weight, radius = row_values(rows)
    hist = lambda x: np.histogram(radius, bins=edges, weights=x)[0]
    total = hist(weight)
    divide = lambda x: np.divide(hist(weight * x), total,
                                 out=np.full_like(total, np.nan), where=total > 0)
    auto = divide(np.abs(value)**2)
    noise = divide(2 * kappa / weight)
    cross = divide(np.real(rows["x"] * rows["y"].conj()))
    return {"uv_klambda": divide(radius), "count": hist(np.ones(len(value))),
            "mean_amplitude_jy": divide(np.abs(value)), "raw_power_jy2": auto,
            "noise_power_jy2": noise, "debiased_power_jy2": auto - noise,
            "cross_power_jy2": cross, "sum_weight": total}


def subset_rows(rows, use):
    return {k: v[use] for k, v in rows.items()}


def average_time_pairs(rows, indices):
    i, j = indices
    out = {k: v[i].copy() for k, v in rows.items()}
    for name, weight_name in [("x", "wx"), ("y", "wy")]:
        w1, w2 = rows[weight_name][i], rows[weight_name][j]
        out[name] = (w1 * rows[name][i] + w2 * rows[name][j]) / (w1 + w2)
        out[weight_name] = w1 + w2
    for key in ("uvw", "frequency", "time"):
        out[key] = (rows[key][i] + rows[key][j]) / 2
    return out


def time_pair_smearing_bound(rows, indices, fwhm_arcsec, reference_ghz):
    """Worst two-sample fringe loss for emission within the PB=0.01 radius."""
    i, j = indices
    frequency = (rows["frequency"][i] + rows["frequency"][j]) / 2
    delta = (rows["uvw"][i] - rows["uvw"][j]) * frequency[:, None] / C_METRES_PER_SECOND
    theta = (fwhm_arcsec * reference_ghz * 1e9 / frequency * np.pi / 648000
             * np.sqrt(np.log(100) / (4*np.log(2))))
    phase = 2*np.pi*(np.hypot(delta[:, 0], delta[:, 1])*theta + np.abs(delta[:, 2])*theta**2/2)
    return {"maximum_phase_change_rad": float(phase.max()),
            "maximum_two_sample_amplitude_loss": float((1-np.cos(phase/2)).max())}


def match_pointings(metadata_aca, metadata_12m, maximum_arcsec):
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    def coordinates(meta):
        pos = np.array([f["PHASE_DIR"] for f in meta["subtables"]["FIELD"]]).squeeze()
        return SkyCoord(pos[:, 0] * u.rad, pos[:, 1] * u.rad)
    aca, alma = coordinates(metadata_aca), coordinates(metadata_12m)
    ids, separation, _ = aca.match_to_catalog_sky(alma)
    candidates = sorted((float(separation[a].arcsec), a, int(ids[a]))
                        for a in range(len(aca)) if separation[a].arcsec <= maximum_arcsec)
    # Deterministic nearest-first one-to-one matching, so no field is duplicated.
    used_12m = set()
    pairs = []
    for sep, a, b in candidates:
        if b in used_12m:
            continue
        used_12m.add(b)
        midpoint = aca[a].directional_offset_by(aca[a].position_angle(alma[b]), aca[a].separation(alma[b]) / 2)
        pairs.append({"pair": len(pairs), "aca_field": a, "alma_field": b,
                      "separation_arcsec": sep, "ra_deg": float(midpoint.ra.deg),
                      "dec_deg": float(midpoint.dec.deg)})
    return pairs


def sample_image_visibilities(image, header, ra_deg, dec_deg, u_lambda, v_lambda,
                              *, fwhm_arcsec, pixel_arcsec=0.5,
                              patch_arcsec=256, padding=2, pb_center=None):
    """FFT a common restored sky through a Gaussian PB, about one phase centre.

    Positive FFT x is east; input astronomical images usually have negative
    CDELT1 and are sampled through WCS. The SIN tangent-plane mapping is exact;
    the w term is neglected and bounded separately. Output units are Jy.
    """
    from astropy.wcs import WCS
    from scipy.ndimage import map_coordinates

    pixel = pixel_arcsec * np.pi / 648000
    n = int(round(patch_arcsec / pixel_arcsec))
    offset = (np.arange(n) - n // 2) * pixel
    l, m = np.meshgrid(offset, offset)
    normal = np.sqrt(1 - l*l - m*m)
    ra0, dec0 = np.deg2rad([ra_deg, dec_deg])
    dec = np.arcsin(m * np.cos(dec0) + normal * np.sin(dec0))
    ra = ra0 + np.arctan2(l, normal * np.cos(dec0) - m * np.sin(dec0))
    x, y = WCS(header, naxis=2, fix=False).celestial.world_to_pixel_values(np.rad2deg(ra) % 360, np.rad2deg(dec))
    sampled = map_coordinates(image, [y, x], order=1, mode="constant", cval=0)
    if pb_center is None:
        distance2 = l*l + m*m
    else:
        center_ra, center_dec = np.deg2rad(pb_center)
        cos_sep = (np.sin(dec) * np.sin(center_dec)
                   + np.cos(dec) * np.cos(center_dec) * np.cos(ra - center_ra))
        distance2 = np.arccos(np.clip(cos_sep, -1, 1))**2
    fwhm = fwhm_arcsec * np.pi / 648000
    patch = sampled * np.exp(-4 * np.log(2) * distance2 / fwhm**2)
    patch *= 1e6 * pixel**2 / normal  # MJy/sr -> Jy/pixel
    size = n * padding
    padded = np.zeros((size, size))
    start = size // 2 - n // 2
    padded[start:start+n, start:start+n] = patch
    transform = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(padded)))
    delta = 1 / (size * pixel)
    coordinates = [np.asarray(v_lambda) / delta + size // 2,
                   np.asarray(u_lambda) / delta + size // 2]
    return (map_coordinates(transform.real, coordinates, order=3, mode="constant", cval=np.nan)
            + 1j * map_coordinates(transform.imag, coordinates, order=3, mode="constant", cval=np.nan))


def binned_sky_model(path, pixel_arcsec):
    """Area-average restored MJy/sr pixels before sampling the sky model."""
    from astropy.io import fits
    from astropy.nddata import block_reduce

    data, header = fits.getdata(path, header=True, memmap=True)
    data = np.squeeze(data)
    if header.get("BUNIT", "").strip() != "MJy/sr":
        raise ValueError("Sky model must be in MJy/sr")
    factor = int(round(pixel_arcsec / (abs(header["CDELT2"]) * 3600)))
    if factor < 1:
        raise ValueError("Model pixels cannot be smaller than the native grid")
    data = block_reduce(np.where(np.isfinite(data), data, 0), factor, func=np.mean)
    for i in (1, 2):
        header[f"CRPIX{i}"] = (header[f"CRPIX{i}"] - 0.5) / factor + 0.5
        header[f"CDELT{i}"] *= factor
        header[f"NAXIS{i}"] = data.shape[2-i]
    return np.asarray(data, dtype=float), header


def _write_rows_csv(path, records):
    if not records:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def _broad_ddids(meta):
    return [i for i, d in enumerate(meta["subtables"]["DATA_DESCRIPTION"])
            if np.sum(np.abs(meta["subtables"]["SPECTRAL_WINDOW"][d["SPECTRAL_WINDOW_ID"]]["CHAN_WIDTH"])) >= 1e9]


def _concat_rows(items):
    return {key: np.concatenate([r[key] for r in items]) for key in items[0]}


def matching_cells(rows, pairs, array, edges, sectors):
    """Pair field, sideband, radial bin and orientation modulo Hermitian symmetry."""
    mapping = {p["aca_field" if array == "aca" else "alma_field"]: p["pair"] for p in pairs}
    pair = np.array([mapping.get(int(f), -1) for f in rows["field"]], dtype=int)
    _, _, radius = row_values(rows)
    radial = np.searchsorted(edges, radius, side="right") - 1
    band = (rows["frequency"] > 225e9).astype(int)
    angle = np.arctan2(rows["uvw"][:, 1], rows["uvw"][:, 0]) % np.pi
    sector = np.minimum((angle / np.pi * sectors).astype(int), sectors - 1)
    good = (pair >= 0) & (radial >= 0) & (radial < len(edges) - 1)
    key = (((pair * 2 + band) * (len(edges) - 1) + radial) * sectors + sector)
    return key, good


def _cell_means(rows, key, good, size, kappa, models=None):
    value, w, radius = row_values(rows)
    histogram = lambda x: np.bincount(key[good], weights=x[good], minlength=size)
    total = histogram(w)
    mean = lambda x: np.divide(histogram(w * x), total, out=np.full(size, np.nan), where=total > 0)
    out = {"count": histogram(np.ones(len(w))), "radius": mean(radius),
           "power": mean(np.abs(value)**2 - 2 * kappa / w),
           "raw_power": mean(np.abs(value)**2)}
    if models:
        for name, values in models.items():
            out[name] = mean(np.abs(values)**2)
    return out


def common_cell_ratios(aca, alma, key_aca, good_aca, key_alma, good_alma,
                       kappa_aca, kappa_alma, pairs, edges, settings,
                       models_aca=None, models_alma=None):
    nb, ns = len(edges) - 1, int(settings["orientation_sectors"])
    size = len(pairs) * 2 * nb * ns
    a = _cell_means(aca, key_aca, good_aca, size, kappa_aca, models_aca)
    b = _cell_means(alma, key_alma, good_alma, size, kappa_alma, models_alma)
    use = (a["count"] >= 4) & (b["count"] >= 4)
    all_keys = np.arange(size)
    radial = (all_keys // ns) % nb
    pair = all_keys // (ns * nb * 2)
    rng = np.random.default_rng(int(settings["random_seed"]))
    multiplicities = np.array([np.bincount(rng.integers(0, len(pairs), len(pairs)), minlength=len(pairs))
                               for _ in range(int(settings["bootstrap_repeats"]))])
    records = []
    for index in range(nb):
        selected = use & (radial == index)
        record = {"uv_low_klambda": float(edges[index]), "uv_high_klambda": float(edges[index+1]),
                  "common_cells": int(selected.sum()), "pointing_pairs": int(len(np.unique(pair[selected])))}
        record["uv_klambda"] = float(np.mean((a["radius"][selected] + b["radius"][selected]) / 2)) if selected.any() else np.nan
        for metric in ["raw_power", "power", *(models_aca or {})]:
            av, bv = a[metric][selected], b[metric][selected]
            numerator, denominator = np.mean(av), np.mean(bv)
            record[metric + "_aca_rms_jy"] = float(np.sqrt(numerator)) if numerator > 0 else np.nan
            record[metric + "_alma_rms_jy"] = float(np.sqrt(denominator)) if denominator > 0 else np.nan
            record[metric + "_ratio"] = float(np.sqrt(numerator / denominator)) if numerator > 0 and denominator > 0 else np.nan
        if selected.any():
            mult = multiplicities[:, pair[selected]]
            av = np.sum(mult * a["power"][selected], axis=1)
            bv = np.sum(mult * b["power"][selected], axis=1)
            valid = (av > 0) & (bv > 0)
            draws = np.sqrt(av[valid] / bv[valid])
            record["bootstrap_positive_fraction"] = float(valid.mean())
            bounds = np.percentile(draws, [16, 84]) if len(draws) >= 0.8 * len(mult) else [np.nan, np.nan]
            record["power_ratio_low"], record["power_ratio_high"] = [float(x) for x in bounds]
        else:
            record.update(bootstrap_positive_fraction=0., power_ratio_low=np.nan, power_ratio_high=np.nan)
        records.append(record)
    return records, {"common_cells": int(use.sum()), "eligible_aca_rows": int(good_aca.sum()),
                     "eligible_alma_rows": int(good_alma.sum()),
                     "retained_aca_rows": int(a["count"][use].sum()),
                     "retained_alma_rows": int(b["count"][use].sum()),
                     "minimum_rows_per_array_per_cell": 4}


def predict_primary_beams(rows, meta, pairs, array, model_path, settings):
    image, header = binned_sky_model(model_path, float(settings["model_pixel_arcsec"]))
    native, common = np.full(len(rows["x"]), np.nan+0j), np.full(len(rows["x"]), np.nan+0j)
    nominal = float(settings["pb_7m_fwhm_arcsec" if array == "aca" else "pb_12m_fwhm_arcsec"])
    max_phase = 0.
    for pair in pairs:
        field = pair["aca_field" if array == "aca" else "alma_field"]
        ra, dec = np.rad2deg(np.asarray(meta["subtables"]["FIELD"][field]["PHASE_DIR"]).ravel())
        for band in (0, 1):
            use = (rows["field"] == field) & ((rows["frequency"] > 225e9) == bool(band))
            if not use.any():
                continue
            frequency = float(np.median(rows["frequency"][use]))
            fwhm = nominal * float(settings["pb_reference_frequency_ghz"]) * 1e9 / frequency
            uv = rows["uvw"][use] * rows["frequency"][use, None] / C_METRES_PER_SECOND
            kwargs = dict(pixel_arcsec=float(settings["model_pixel_arcsec"]),
                          patch_arcsec=float(settings["model_patch_arcsec"]),
                          padding=int(settings["model_padding_factor"]))
            native[use] = sample_image_visibilities(image, header, ra % 360, dec, uv[:, 0], uv[:, 1], fwhm_arcsec=fwhm, **kwargs)
            common[use] = sample_image_visibilities(image, header, ra % 360, dec, uv[:, 0], uv[:, 1],
                fwhm_arcsec=float(settings["pb_12m_fwhm_arcsec"]) * float(settings["pb_reference_frequency_ghz"]) * 1e9 / frequency,
                pb_center=(pair["ra_deg"], pair["dec_deg"]), **kwargs)
            # At PB=0.01, the neglected w phase is small; outer emission is strongly attenuated.
            theta = fwhm * np.pi / 648000 * np.sqrt(np.log(100) / (4 * np.log(2)))
            max_phase = max(max_phase, float(np.max(np.abs(uv[:, 2])) * np.pi * theta**2))
        print(f"  PB model {model_path.name}: {array} pair {pair['pair']+1}/{len(pairs)}", flush=True)
    return native, common, max_phase


def check_model_sampling(rows, meta, pair, array, model_path, settings, native, common):
    """Halve model-grid sampling for the nearest matched pointing, all its rows."""
    fine_pixel = float(settings["model_pixel_arcsec"]) / 2
    image, header = binned_sky_model(model_path, fine_pixel)
    field = pair["aca_field" if array == "aca" else "alma_field"]
    ra, dec = np.rad2deg(np.asarray(meta["subtables"]["FIELD"][field]["PHASE_DIR"]).ravel())
    nominal = settings["pb_7m_fwhm_arcsec" if array == "aca" else "pb_12m_fwhm_arcsec"]
    checks = []
    for band in (0, 1):
        use = (rows["field"] == field) & ((rows["frequency"] > 225e9) == bool(band))
        if not use.any():
            continue
        frequency = float(np.median(rows["frequency"][use]))
        uv = rows["uvw"][use] * rows["frequency"][use, None] / C_METRES_PER_SECOND
        weights = rows["wx"][use] + rows["wy"][use]
        for mode, width, center, reference in [
            ("native", nominal, None, native[use]),
            ("common", settings["pb_12m_fwhm_arcsec"], (pair["ra_deg"], pair["dec_deg"]), common[use])]:
            fine = sample_image_visibilities(image, header, ra % 360, dec, uv[:, 0], uv[:, 1],
                fwhm_arcsec=width*settings["pb_reference_frequency_ghz"]*1e9/frequency,
                pixel_arcsec=fine_pixel, patch_arcsec=settings["model_patch_arcsec"],
                padding=settings["model_padding_factor"], pb_center=center)
            error = np.sqrt(np.sum(weights*np.abs(fine-reference)**2)/np.sum(weights*np.abs(fine)**2))
            checks.append({"array": array, "band": band, "pb": mode,
                           "rows": int(use.sum()), "relative_complex_l2_change": float(error)})
    return checks


def run_overlap_tests(paths, config, *, save_figure, overwrite=False):
    """Run and retain all three requested factor tests in the owning task."""
    settings = config["overlap_tests"]
    edges = np.array(settings["uv_edges_klambda"], dtype=float)
    root = Path(paths["output_root"]) / "overlap_tests"
    root.mkdir(parents=True, exist_ok=True)
    if (root / "summary.json").exists() and not overwrite:
        raise FileExistsError("Overlap tests already exist; use --overwrite")
    rows, metadata, kappa, noise, profiles, smearing = {}, {}, {}, {}, {}, {}
    for name in ARRAYS:
        cache = root / f"{name}_rows.npz"
        meta_path = root / f"{name}_metadata.json"
        if cache.is_file() and meta_path.is_file():
            with np.load(cache) as z:
                rows[name] = {k: z[k] for k in z.files}
            metadata[name] = json.loads(meta_path.read_text())
            if metadata[name]["test_max_klambda"] != settings["read_max_klambda"]:
                raise ValueError("Row-cache range differs from configuration; remove this cache before rerunning")
        else:
            rows[name], metadata[name] = read_test_rows(Path(paths["measurement_sets"][name]), float(settings["read_max_klambda"]))
            np.savez_compressed(cache, **rows[name])
            meta_path.write_text(json.dumps(metadata[name], indent=2) + "\n")
        kappa[name], noise[name], pairs_time = empirical_noise(rows[name], float(settings["maximum_time_pair_gap_seconds"]))
        broad = np.isin(rows[name]["ddid"], _broad_ddids(metadata[name]))
        profiles[name + "_all"] = radial_test_profile(rows[name], edges, kappa[name])
        profiles[name + "_broad"] = radial_test_profile(subset_rows(rows[name], broad), edges, kappa[name][broad])
        i, j = pairs_time
        use = broad[i] & broad[j]
        i, j = i[use], j[use]
        smearing[name] = time_pair_smearing_bound(rows[name], (i, j),
            float(settings["pb_7m_fwhm_arcsec" if name == "aca_7m" else "pb_12m_fwhm_arcsec"]),
            float(settings["pb_reference_frequency_ghz"]))
        before = np.r_[i, j]
        profiles[name + "_pair_before"] = radial_test_profile(subset_rows(rows[name], before), edges, kappa[name][before])
        averaged = average_time_pairs(rows[name], (i, j))
        profiles[name + "_pair_after"] = radial_test_profile(averaged, edges, kappa[name][i])
        print(name, "noise factors", noise[name], flush=True)
    pairs = match_pointings(metadata["aca_7m"], metadata["configuration_c5"], float(settings["maximum_pair_separation_arcsec"]))
    _write_rows_csv(root / "matched_pointings.csv", pairs)
    aca_broad = np.isin(rows["aca_7m"]["ddid"], _broad_ddids(metadata["aca_7m"]))
    aca = subset_rows(rows["aca_7m"], aca_broad)
    ka = kappa["aca_7m"][aca_broad]
    alma = _concat_rows([rows["configuration_c5"], rows["configuration_c2"]])
    kb = np.r_[kappa["configuration_c5"], kappa["configuration_c2"]]
    ca, ga = matching_cells(aca, pairs, "aca", edges, int(settings["orientation_sectors"]))
    cb, gb = matching_cells(alma, pairs, "alma", edges, int(settings["orientation_sectors"]))
    # Forward-model only rows in paired fields and the tested uv interval.
    aca, ka, ca = subset_rows(aca, ga), ka[ga], ca[ga]
    alma, kb, cb = subset_rows(alma, gb), kb[gb], cb[gb]
    ga, gb = np.ones(len(ca), bool), np.ones(len(cb), bool)
    model_a, model_b, model_notes, sampling_checks = {}, {}, {}, []
    for label, image_name in [("allscale", "allscale_sky"), ("twelve_m", "twelve_m_pbcor")]:
        path = Path(paths["images"][image_name])
        a, ac, wa = predict_primary_beams(aca, metadata["aca_7m"], pairs, "aca", path, settings)
        b, bc, wb = predict_primary_beams(alma, metadata["configuration_c5"], pairs, "alma", path, settings)
        model_a[label + "_native"], model_a[label + "_common"] = a, ac
        model_b[label + "_native"], model_b[label + "_common"] = b, bc
        model_notes[label] = {"input": str(path), "maximum_w_phase_at_pb_001_rad": max(wa, wb)}
        if label == "allscale":
            sampling_checks.extend(check_model_sampling(aca, metadata["aca_7m"], pairs[0], "aca", path, settings, a, ac))
            sampling_checks.extend(check_model_sampling(alma, metadata["configuration_c5"], pairs[0], "alma", path, settings, b, bc))
    ratios, cell_counts = common_cell_ratios(aca, alma, ca, ga, cb, gb, ka, kb,
        pairs, edges, settings, model_a, model_b)
    strict = [{**p, "pair": i} for i, p in enumerate(
        p for p in pairs if p["separation_arcsec"] <= settings["strict_pair_separation_arcsec"])]
    sca, sga = matching_cells(aca, strict, "aca", edges, int(settings["orientation_sectors"]))
    scb, sgb = matching_cells(alma, strict, "alma", edges, int(settings["orientation_sectors"]))
    strict_ratios, strict_counts = common_cell_ratios(aca, alma, sca, sga, scb, sgb, ka, kb,
        strict, edges, settings, model_a, model_b)
    for collection in (ratios, strict_ratios):
        for record in collection:
            for label in model_notes:
                factor = record[label + "_common_ratio"] / record[label + "_native_ratio"]
                record[label + "_pb_adjusted_ratio"] = record["power_ratio"] * factor
    _write_rows_csv(root / "common_cell_ratios.csv", ratios)
    _write_rows_csv(root / "strict_common_cell_ratios.csv", strict_ratios)
    profile_records = []
    for name, profile in profiles.items():
        for i in range(len(edges) - 1):
            profile_records.append({"selection": name, **{k: float(v[i]) for k, v in profile.items()}})
    _write_rows_csv(root / "noise_profiles.csv", profile_records)
    np.savez_compressed(root / "matched_model_predictions.npz", aca_keys=ca, alma_keys=cb,
                        **{"aca_" + k: v for k, v in model_a.items()},
                        **{"alma_" + k: v for k, v in model_b.items()})
    summary = {"noise": noise, "matched_pointings": pairs, "common_cells": cell_counts,
               "models": model_notes, "settings": dict(settings),
               "common_cell_ratios": ratios, "strict_common_cell_ratios": strict_ratios,
               "strict_common_cells": strict_counts, "numerical_sampling_checks": sampling_checks,
               "time_pair_smearing_bound": smearing}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _plot_overlap_tests(root, profiles, pairs, metadata, ratios, save_figure)
    _plot_pb_amplitude_curve(root, ratios, save_figure)
    return summary


def run_pb_amplitude_curve(paths, config, *, save_figure, overwrite=False):
    """Reconstruct absolute matched-cell amplitudes using the existing test caches."""
    root = Path(paths["output_root"]) / "overlap_tests"
    if (root / "aca_pb_amplitude_curve.pdf").exists() and not overwrite:
        raise FileExistsError("PB amplitude curve exists; use --overwrite")
    summary = json.loads((root / "summary.json").read_text())
    settings = config["overlap_tests"]
    if dict(settings) != summary["settings"]:
        raise ValueError("Settings changed; rerun overlap-tests before plotting cached predictions")
    edges = np.asarray(settings["uv_edges_klambda"])
    pairs = summary["matched_pointings"]
    rows, kappa, metadata = {}, {}, {}
    for name in ARRAYS:
        with np.load(root / f"{name}_rows.npz") as z:
            rows[name] = {k: z[k] for k in z.files}
        metadata[name] = json.loads((root / f"{name}_metadata.json").read_text())
        kappa[name] = np.array([summary["noise"][name][str(int(d))]["kappa_time"]
                                for d in rows[name]["ddid"]])
    take = np.isin(rows["aca_7m"]["ddid"], _broad_ddids(metadata["aca_7m"]))
    aca, ka = subset_rows(rows["aca_7m"], take), kappa["aca_7m"][take]
    alma = _concat_rows([rows["configuration_c5"], rows["configuration_c2"]])
    kb = np.r_[kappa["configuration_c5"], kappa["configuration_c2"]]
    ca, ga = matching_cells(aca, pairs, "aca", edges, int(settings["orientation_sectors"]))
    cb, gb = matching_cells(alma, pairs, "alma", edges, int(settings["orientation_sectors"]))
    aca, ka, ca = subset_rows(aca, ga), ka[ga], ca[ga]
    alma, kb, cb = subset_rows(alma, gb), kb[gb], cb[gb]
    with np.load(root / "matched_model_predictions.npz") as z:
        if not np.array_equal(ca, z["aca_keys"]) or not np.array_equal(cb, z["alma_keys"]):
            raise ValueError("Cached row selection differs from the saved PB predictions")
        model_a = {k[4:]: z[k] for k in z.files if k.startswith("aca_") and k != "aca_keys"}
        model_b = {k[5:]: z[k] for k in z.files if k.startswith("alma_") and k != "alma_keys"}
    records, _ = common_cell_ratios(aca, alma, ca, np.ones(len(ca), bool),
        cb, np.ones(len(cb), bool), ka, kb, pairs, edges, settings, model_a, model_b)
    return _plot_pb_amplitude_curve(root, records, save_figure)


def _plot_pb_amplitude_curve(root, records, save_figure):
    """Translate each array's ensemble signal RMS to the same modelled PB."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    output = []
    for r in records:
        row = {"uv_klambda": r["uv_klambda"], "common_cells": r["common_cells"],
               "aca_native_mjy": 1000*r["power_aca_rms_jy"],
               "alma_native_mjy": 1000*r["power_alma_rms_jy"]}
        for model in ("allscale", "twelve_m"):
            for array in ("aca", "alma"):
                factor = r[f"{model}_common_{array}_rms_jy"] / r[f"{model}_native_{array}_rms_jy"]
                row[f"{array}_{model}_common_pb_mjy"] = 1000*r[f"power_{array}_rms_jy"]*factor
                row[f"{array}_{model}_pb_factor"] = factor
        output.append(row)
    _write_rows_csv(root / "aca_pb_amplitude_curve.csv", output)
    columns = {k: np.array([r[k] for r in output]) for k in output[0]}
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.16, top=0.74)
    x = columns["uv_klambda"]
    ax.plot(x, columns["aca_native_mjy"], color="#D87821", linestyle="--", linewidth=1.8,
            marker="o", markersize=4, label="ACA 7 m: native PB")
    ax.plot(x, columns["aca_allscale_common_pb_mjy"], color="#D87821", linewidth=2,
            marker="o", markersize=4, label="ACA 7 m: common PB (model)")
    ax.plot(x, columns["alma_allscale_common_pb_mjy"], color="#007A9E", linewidth=1.8,
            marker="s", markersize=3.5, label="12 m: common PB (model)")
    first, second = columns["aca_allscale_common_pb_mjy"], columns["aca_twelve_m_common_pb_mjy"]
    ax.fill_between(x, np.minimum(first, second), np.maximum(first, second), color="#D87821", alpha=0.15,
                    label="Range from two sky models")
    ax.set(xlabel=r"uv distance (k$\lambda$)", ylabel="Noise-corrected RMS amplitude (mJy)", xlim=(8, 32))
    fig.legend(*ax.get_legend_handles_labels(), fontsize=8.5, ncol=2,
               loc="upper center", bbox_to_anchor=(0.58, 0.93))
    ax.tick_params(direction="out")
    ax.grid(axis="y", alpha=0.12, linewidth=0.5)
    fig.text(0.16, 0.98, "Matched sky and uv sampling", fontsize=10, ha="left", va="top")
    save_figure(fig, [root / f"aca_pb_amplitude_curve.{s}" for s in ("pdf", "svg", "png")], dpi=300)
    plt.close(fig)
    return output


def _positive_sqrt(values):
    values = np.asarray(values)
    return np.sqrt(np.where(values > 0, values, np.nan))


def _plot_overlap_tests(root, profiles, pairs, metadata, ratios, save_figure):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    def save(fig, name):
        save_figure(fig, [root / f"{name}.{suffix}" for suffix in ["pdf", "svg", "png"]], dpi=300)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.4))
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.11, top=0.88, wspace=0.38, hspace=0.48)
    a = profiles["aca_7m_all"]
    axes[0, 0].plot(a["uv_klambda"], 1000*a["mean_amplitude_jy"], color="0.55", linestyle="--", label="ACA: all windows")
    for name, (label, color) in ARRAYS.items():
        p = profiles[name + "_broad"]
        axes[0, 0].plot(p["uv_klambda"], 1000*p["mean_amplitude_jy"], color=color, label=label+": broad")
        for stage, style in [("before", ":"), ("after", "-")]:
            p = profiles[name + "_pair_" + stage]
            axes[0, 1].plot(p["uv_klambda"], 1000*p["mean_amplitude_jy"], color=color, linestyle=style)
    axes[0, 0].legend(fontsize=7.5, loc="upper right")
    axes[0, 0].set_title("a  Spectral-window selection", loc="left", fontsize=10)
    axes[0, 1].set_title("b  Same time pairs", loc="left", fontsize=10)
    axes[0, 1].legend(handles=[Line2D([], [], color="0.3", linestyle=":", label="Before averaging"),
                             Line2D([], [], color="0.3", label="Two integrations averaged")], fontsize=7.5)
    for axis in axes[0]:
        axis.set(xlabel=r"uv distance (k$\lambda$)", ylabel="Mean amplitude (mJy)", xlim=(8, 32))
    # This panel compares two independent ways of measuring the same noise correction.
    for name, (label, color) in ARRAYS.items():
        p = profiles[name + "_broad"]
        axes[1, 0].plot(p["uv_klambda"], 1000*_positive_sqrt(p["debiased_power_jy2"]), color=color, label=label)
        axes[1, 0].plot(p["uv_klambda"], 1000*_positive_sqrt(p["cross_power_jy2"]), color=color, linestyle=":")
        p1, p2 = profiles[name + "_pair_before"], profiles[name + "_pair_after"]
        axes[1, 1].plot(p1["uv_klambda"], p2["mean_amplitude_jy"] / p1["mean_amplitude_jy"], color=color)
    axes[1, 0].set_title("c  Signal RMS: two estimators", loc="left", fontsize=10)
    axes[1, 0].legend(handles=[Line2D([], [], color="0.3", label="Auto-power minus noise"),
                             Line2D([], [], color="0.3", linestyle=":", label="XX/YY cross-power")], fontsize=7.5)
    axes[1, 0].set(xlabel=r"uv distance (k$\lambda$)", ylabel="Signal RMS (mJy)", xlim=(8, 32))
    axes[1, 1].axhline(1, color="0.6", linewidth=0.8, linestyle="--")
    axes[1, 1].set_title("d  Coherent-averaging effect", loc="left", fontsize=10)
    axes[1, 1].set(xlabel=r"uv distance (k$\lambda$)", ylabel="Mean amplitude: after / before", xlim=(8, 32))
    fig.legend(handles=[Line2D([], [], color=c, label=l) for l, c in ARRAYS.values()],
               ncol=3, loc="upper center", bbox_to_anchor=(0.55, 0.99), fontsize=8)
    save(fig, "noise_bias_tests")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.8))
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.18, top=0.84, wspace=0.42)
    ra0, dec0 = 221.25, -59.8
    for name, color, marker, label in [("aca_7m", "#D87821", "o", "ACA pointings"),
                                       ("configuration_c5", "#007A9E", "+", "12 m pointings")]:
        pos = np.array([x["PHASE_DIR"] for x in metadata[name]["subtables"]["FIELD"]]).squeeze()
        ra, dec = np.rad2deg(pos[:, 0]) % 360, np.rad2deg(pos[:, 1])
        axes[0].scatter((ra-ra0)*np.cos(np.deg2rad(dec0))*60, (dec-dec0)*60,
                        color=color, marker=marker, s=10, linewidths=0.6, alpha=0.65, label=label)
    x = [(p["ra_deg"]-ra0)*np.cos(np.deg2rad(dec0))*60 for p in pairs]
    y = [(p["dec_deg"]-dec0)*60 for p in pairs]
    axes[0].scatter(x, y, marker="s", s=24, facecolors="none", edgecolors="0.15", linewidths=0.8, label="Matched pairs")
    axes[0].set(xlabel="RA offset (arcmin)", ylabel="Dec offset (arcmin)")
    axes[0].set_aspect("equal", adjustable="datalim")
    axes[0].invert_xaxis()
    axes[0].set_title("a  Pointing correspondence", loc="left", fontsize=10)
    axes[0].legend(fontsize=7, loc="upper left")
    aca = profiles["aca_7m_broad"]
    b1, b2 = profiles["configuration_c5_broad"], profiles["configuration_c2_broad"]
    combined_weight = b1["sum_weight"] + b2["sum_weight"]
    for metric, label, color, style in [("raw_power_jy2", "All sky: raw RMS", "0.6", ":"),
                                        ("debiased_power_jy2", "All sky: noise corrected", "#D87821", "--")]:
        alma = (b1[metric]*b1["sum_weight"] + b2[metric]*b2["sum_weight"]) / combined_weight
        axes[1].plot(aca["uv_klambda"], _positive_sqrt(aca[metric] / alma), color=color, linestyle=style, label=label)
    x = np.array([r["uv_klambda"] for r in ratios])
    ratio = np.array([r["power_ratio"] for r in ratios])
    low, high = [np.array([r[k] for r in ratios]) for k in ["power_ratio_low", "power_ratio_high"]]
    axes[1].plot(x, ratio, color="0.15", marker="o", markersize=3, label="Matched sky, band and uv cells")
    axes[1].fill_between(x, low, high, color="0.3", alpha=0.12)
    axes[1].axhline(1, color="0.6", linewidth=0.8)
    axes[1].set(xlabel=r"uv distance (k$\lambda$)", ylabel="ACA / 12 m RMS ratio", xlim=(8, 32))
    axes[1].set_title("b  Sky and sampling selection", loc="left", fontsize=10)
    fig.legend(*axes[1].get_legend_handles_labels(), fontsize=7, ncol=3,
               loc="upper center", bbox_to_anchor=(0.54, 1.02))
    save(fig, "sky_match_tests")

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 4.1))
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.17, top=0.76, wspace=0.42)
    for model, color, label in [("allscale", "#3B7D4A", "All-scale sky model"),
                                 ("twelve_m", "#777777", "12 m sky model")]:
        native, common = [np.array([r[model + k + "_ratio"] for r in ratios]) for k in ["_native", "_common"]]
        axes[0].plot(x, native, color=color, label=label)
        axes[0].plot(x, common, color=color, linestyle=":")
        adjusted = np.array([r[model + "_pb_adjusted_ratio"] for r in ratios])
        axes[1].plot(x, adjusted, color=color, label=label+": PB adjusted")
    axes[1].plot(x, ratio, color="#D87821", marker="o", markersize=3, label="Measured: noise and sky matched")
    axes[1].fill_between(x, low, high, color="#D87821", alpha=0.10)
    for axis in axes:
        axis.axhline(1, color="0.6", linewidth=0.8)
        axis.set(xlabel=r"uv distance (k$\lambda$)", ylabel="ACA / 12 m RMS ratio", xlim=(8, 32))
    axes[0].set_title("a  Common-sky PB experiment", loc="left", fontsize=10)
    axes[1].set_title("b  Conditional PB adjustment", loc="left", fontsize=10)
    handles = [Line2D([], [], color="#3B7D4A", label="All-scale sky model"),
               Line2D([], [], color="#777777", label="12 m sky model"),
               Line2D([], [], color="#D87821", marker="o", markersize=3, label="Measured, matched sky")]
    fig.legend(handles=handles, ncol=3, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.54, 1.01))
    fig.text(0.30, 0.90, "Panel a: solid = native beams; dotted = common beam", ha="center", va="top", fontsize=7.5)
    save(fig, "primary_beam_tests")
