"""Phase 1 stages of the combine_continuum task: prepare, TP factor, J-comb and Feather, box spectra, Planck calibration."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
import json
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np
from scipy import ndimage

from lancet.combination.uvplane import gaussian_transfer, uv_distance_klambda
from lancet.combination import planck_calibration as pc
from lancet.combination.figures import (
    ARCSEC_PER_KLAMBDA,
    plot_comparison,
    plot_box_spectrum,
    plot_planck_maps,
    plot_planck_scatter,
    plot_power_spectrum,
    plot_tp_factor,
    plot_zero_spacing,
)
from lancet.combination.overlap_tests import require_inputs
from lancet.imaging import frequency_label, read_frequency_ghz
from lancet.paths import guard_outputs
from lancet.combination.prep import (
    azimuthal_power,
    beam_parameters,
    combine_prepared,
    footprint_margin_arcsec,
    jybeam_to_mjysr,
    mjysr_to_jybeam,
    read_plane,
    regrid_to,
    taper_edge,
    write_product,
)
from lancet.combination.validation import (
    cross_spectrum_slopes,
    remove_window_mean,
    beam_volume_outside,
    core_contrast_ratio,
    large_scale_ratio,
    mosaic_window,
    negative_fraction,
    psd_match_factor,
    smooth_to,
    synthetic_gate,
    tapered_window,
    tp_at_high_resolution,
    tp_only_ratio,
    windowed_amplitude_ratio,
    windowed_power,
    zero_spacing_check,
)

TP_FLAVOURS = ("dust", "dust_ff")
TP_LABELS = {"dust": "TP dust", "dust_ff": "TP dust+free-free"}
METHOD_LABELS = {"jcomb": "J-comb", "feather": "Feather"}
FLAVOUR_KIND = {"dust_ff": "continuum", "dust": "dust"}
PADDING_MODE = "TP: mirror reflection with a cosine roll-off; 12m+7m: zeros"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def input_hashes(paths: Mapping[str, Path], names: Sequence[str]) -> dict[str, str]:
    return {name: sha256(paths[name]) for name in names}


def write_json(path: Path, payload: Mapping | Sequence) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1, default=str) + "\n")


def ratio_check(value: float, tolerance: float) -> dict:
    """A ratio expected to be 1, its tolerance and whether it is exceeded."""
    return {"value": float(value), "tolerance": float(tolerance), "flag": bool(abs(value - 1.0) > tolerance)}


# ---------------------------------------------------------------- prepare

def _pixel_arcsec(wcs) -> float:
    return float(proj_plane_pixel_scales(wcs)[0] * 3600.0)


def _finer_wcs(wcs, shape, pixel_arcsec):
    """Same extent as ``wcs``/``shape`` with a finer pixel."""
    factor = _pixel_arcsec(wcs) / pixel_arcsec
    fine = wcs.deepcopy()
    fine.wcs.cdelt = wcs.wcs.cdelt / factor
    fine.wcs.crpix = (wcs.wcs.crpix - 0.5) * factor + 0.5
    return fine, (int(round(shape[0] * factor)), int(round(shape[1] * factor)))


def regrid_power_deviation(tp, tp_wcs, settings) -> float:
    """Max |power ratio - 1| between the tapered native TP and its spline regrid."""
    pix = _pixel_arcsec(tp_wcs)
    tapered = taper_edge(tp, np.isfinite(tp), settings["edge_taper_arcsec"], pix)
    fine_wcs, fine_shape = _finer_wcs(tp_wcs, tp.shape, settings["regrid_check_pixel_arcsec"])
    fine = regrid_to(tapered, tp_wcs, fine_wcs, fine_shape, order=settings["regrid_order"])
    edges = settings["regrid_check_edges_klambda"]
    native = azimuthal_power(tapered, pix, edges)
    regridded = azimuthal_power(np.nan_to_num(fine), settings["regrid_check_pixel_arcsec"], edges)
    return float(np.max(np.abs(regridded / native - 1.0)))


def tp_footprint_on_grid(tp, tp_wcs, target_wcs, shape) -> np.ndarray:
    """Pixels of the target grid where the native TP model is defined."""
    marker = np.where(np.isfinite(tp), 1.0, np.nan).astype(np.float32)
    return np.isfinite(regrid_to(marker, tp_wcs, target_wcs, shape, order=0))


def read_gain(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    """Read the delivered gain map as a float32 2-D plane of ``shape``; no unit or beam is expected."""
    data = np.squeeze(np.asarray(fits.getdata(path)))
    if data.shape != tuple(shape):
        raise ValueError(f"{path}: gain shape {data.shape} differs from the interferometer's {tuple(shape)}")
    return data.astype(np.float32)


def mosaic_gain(gain: np.ndarray, interferometer: np.ndarray) -> np.ndarray:
    """The delivered gain where the interferometer image is finite, 0 elsewhere (float32)."""
    return np.where(np.isfinite(interferometer), gain, 0.0).astype(np.float32)


def check_tp(tp, tp_wcs, gain, wcs, settings) -> dict:
    """Footprint margin around the mosaic and the regridding power deviation of one TP model."""
    step = int(settings["margin_decimation"])
    coarse_wcs, mosaic = wcs[::step, ::step], gain[::step, ::step] > 0
    footprint = tp_footprint_on_grid(tp, tp_wcs, coarse_wcs, mosaic.shape)
    margin = footprint_margin_arcsec(footprint, mosaic, _pixel_arcsec(wcs) * step)
    deviation = regrid_power_deviation(tp, tp_wcs, settings)
    problems = []
    if margin < settings["minimum_margin_arcsec"]:
        problems.append(f"footprint margin {margin:.1f} arcsec < {settings['minimum_margin_arcsec']} arcsec")
    if deviation > settings["regrid_power_tolerance"]:
        problems.append(f"regrid power deviation {deviation:.3f} > {settings['regrid_power_tolerance']}")
    return {"margin_arcsec": margin, "regrid_power_deviation": deviation, "problems": problems}


def tp_frequency_ghz(paths: Mapping[str, Path]) -> float:
    """``FREQGHZ`` of the TP models; ValueError when a model lacks it or the models differ."""
    frequencies = {f"tp_{flavour}": read_frequency_ghz(paths[f"tp_{flavour}"]) for flavour in TP_FLAVOURS}
    if len(set(frequencies.values())) > 1:
        raise ValueError("TP models differ in FREQGHZ: "
                         + ", ".join(f"{name} {value} GHz" for name, value in frequencies.items()))
    return next(iter(frequencies.values()))


def prepared_paths(root: Path) -> dict[str, Path]:
    prepared = Path(root) / "prepared"
    out = {"interferometer": prepared / "interferometer_mjysr.fits", "gain": prepared / "gain.fits",
           "report": prepared / "input_check.json"}
    for flavour in TP_FLAVOURS:
        out[f"tp_{flavour}"] = prepared / f"tp_{flavour}_mjysr.fits"
        out[f"tp_{flavour}_footprint"] = prepared / f"tp_{flavour}_footprint.fits"
    return out


def extended_prepared_paths(root: Path) -> dict[str, Path]:
    """Prepared planes of the optional extended TP models (``tp_<flavour>_extended`` inputs)."""
    prepared = Path(root) / "prepared"
    out = {}
    for flavour in TP_FLAVOURS:
        out[f"tp_{flavour}"] = prepared / f"tp_{flavour}_extended_mjysr.fits"
        out[f"tp_{flavour}_footprint"] = prepared / f"tp_{flavour}_extended_footprint.fits"
    return out


def _has_extended_tp(paths: Mapping[str, Path]) -> bool:
    return all(f"tp_{flavour}_extended" in paths for flavour in TP_FLAVOURS)


def extended_context(ctx: Mapping, paths: Mapping[str, Path], flavour: str) -> Mapping:
    """``ctx`` with the TP planes of the extended TP model, when one is declared; else ``ctx``."""
    if f"tp_{flavour}_extended" not in paths:
        return ctx
    prep = extended_prepared_paths(paths["output_root"])
    require_inputs([prep[f"tp_{flavour}"], prep[f"tp_{flavour}_footprint"]], "prepare")
    tp_native, tp_header, tp_wcs = read_plane(paths[f"tp_{flavour}_extended"])
    tp_native = jybeam_to_mjysr(tp_native.astype(float), *beam_parameters(tp_header)[:2])
    return {**ctx, "low": fits.getdata(prep[f"tp_{flavour}"]).astype(np.float32),
            "footprint": fits.getdata(prep[f"tp_{flavour}_footprint"]) > 0.5,
            "tp_native": tp_native, "tp_header": tp_header, "tp_wcs": tp_wcs}


def interferometer_min_pb(config: Mapping) -> float:
    """``[prepare] interferometer_min_pb``: 0.2 (the delivered cut) up to, not including, 1."""
    value = float(config["prepare"]["interferometer_min_pb"])
    if not 0.2 <= value < 1.0:
        raise ValueError(f"[prepare] interferometer_min_pb = {value} must satisfy 0.2 <= value < 1 "
                         "(the delivered true-sky image is already blank below pb 0.2)")
    return value


def cut_at_pb(truesky: np.ndarray, pb: np.ndarray, minimum_pb: float) -> np.ndarray:
    """The true-sky image with NaN wherever the primary beam is below ``minimum_pb``."""
    return np.where(pb < minimum_pb, np.nan, truesky).astype(np.float32)


def run_prepare(paths: Mapping[str, Path], config: Mapping, overwrite: bool = False) -> dict[str, Path]:
    """Check the inputs and write them on the interferometer grid in MJy/sr, plus the mosaic gain."""
    settings = config["prepare"]
    out = prepared_paths(paths["output_root"])
    extended_out = extended_prepared_paths(paths["output_root"]) if _has_extended_tp(paths) else {}
    guard_outputs(list(out.values()) + list(extended_out.values()), overwrite)
    minimum_pb = interferometer_min_pb(config)
    frequency = tp_frequency_ghz(paths)                 # before anything is read or written
    for flavour in TP_FLAVOURS if extended_out else ():
        extended_frequency = read_frequency_ghz(paths[f"tp_{flavour}_extended"])
        if extended_frequency != frequency:
            raise ValueError(f"tp_{flavour}_extended FREQGHZ {extended_frequency} differs from the TP "
                             f"models' {frequency} GHz")
    truesky, header, wcs = read_plane(paths["interferometer"])
    truesky = cut_at_pb(truesky, read_gain(paths["pb"], truesky.shape), minimum_pb)
    delivered_gain = read_gain(paths["gain"], truesky.shape)
    high_beam = beam_parameters(header)
    gain = mosaic_gain(delivered_gain, truesky)
    del delivered_gain
    report = {"interferometer": {"path": str(paths["interferometer"]), "shape": list(truesky.shape),
                                 "beam_arcsec_deg": high_beam, "pixel_arcsec": _pixel_arcsec(wcs),
                                 "minimum_pb": minimum_pb},
              "gain": {"min_nonzero": float(gain[gain > 0].min()), "max": float(gain.max())}}
    tp_planes, problems = {}, []
    for flavour in TP_FLAVOURS:
        tp, tp_header, tp_wcs = read_plane(paths[f"tp_{flavour}"])
        check = check_tp(tp, tp_wcs, gain, wcs, settings)
        problems += [f"tp_{flavour}: {p}" for p in check.pop("problems")]
        report[f"tp_{flavour}"] = {"path": str(paths[f"tp_{flavour}"]), "shape": list(tp.shape),
                                   "beam_arcsec_deg": beam_parameters(tp_header), **check}
        tp_planes[flavour] = (tp, tp_header, tp_wcs)
    print(json.dumps(report, indent=1, default=str), flush=True)
    if problems:
        raise RuntimeError("; ".join(problems))

    hashes = input_hashes(paths, ["interferometer", "gain", "pb"])
    write_product(out["gain"], gain, header,
                  {"stage": "prepare", "interferometer_min_pb": minimum_pb,
                   "quantity": f"12m+7m mosaic gain; 0 where the image is blank (pb < {minimum_pb})"}, hashes, bunit="")
    write_product(out["interferometer"], jybeam_to_mjysr(np.nan_to_num(truesky), *high_beam[:2]),
                  header, {"stage": "prepare", "quantity": "12m+7m true-sky image; 0 outside the mosaic"}, hashes,
                  bunit="MJy/sr")
    planes = [(out, f"tp_{flavour}", flavour, plane) for flavour, plane in tp_planes.items()]
    planes += [(extended_out, f"tp_{flavour}_extended", flavour, read_plane(paths[f"tp_{flavour}_extended"]))
               for flavour in TP_FLAVOURS if extended_out]
    for target, name, flavour, (tp, tp_header, tp_wcs) in planes:
        tp_beam = beam_parameters(tp_header)
        tapered = taper_edge(tp, np.isfinite(tp), settings["edge_taper_arcsec"], _pixel_arcsec(tp_wcs))
        grid = np.nan_to_num(regrid_to(tapered, tp_wcs, wcs, truesky.shape, order=settings["regrid_order"]))
        tp_hash = input_hashes(paths, [name])
        meta = {"stage": "prepare", "tp_frequency_ghz": frequency,
                "edge_taper_arcsec": settings["edge_taper_arcsec"]}
        write_product(target[f"tp_{flavour}"], jybeam_to_mjysr(grid, *tp_beam[:2]), header, meta, tp_hash,
                      bunit="MJy/sr", beam=tp_beam)
        footprint = tp_footprint_on_grid(tp, tp_wcs, wcs, truesky.shape).astype(np.float32)
        write_product(target[f"tp_{flavour}_footprint"], footprint, header,
                      {"stage": "prepare", "quantity": "TP model defined (1) or not (0)"}, tp_hash, bunit="")
    write_json(out["report"], report)
    return out


def require_prepared(paths: Mapping[str, Path]) -> dict[str, Path]:
    """The prepared images; fails naming ``prepare`` when one is missing."""
    prep = prepared_paths(paths["output_root"])
    require_inputs([p for name, p in prep.items() if name != "report"], "prepare")
    return prep


def load_prepared(paths: Mapping[str, Path], flavour: str):
    """(tp MJy/sr, interferometer MJy/sr, gain, header, pixel", high beam, low beam)."""
    prep = require_prepared(paths)
    high = fits.getdata(prep["interferometer"]).astype(np.float32)
    header = fits.getheader(prep["interferometer"])
    low = fits.getdata(prep[f"tp_{flavour}"]).astype(np.float32)
    low_beam = beam_parameters(fits.getheader(prep[f"tp_{flavour}"]))
    gain = fits.getdata(prep["gain"]).astype(np.float32)
    return low, high, gain, header, _pixel_arcsec(WCS(header)), beam_parameters(header), low_beam


# ---------------------------------------------------------------- TP factor

def measure_tp_factor(tp, interferometer_at_tp, mask, pix, settings) -> dict:
    """TP factor that lifts the TP power spectrum onto that of the 12m+7m image at the TP resolution.

    Both spectra are taken inside ``mask`` tapered at its edge; the factor is the 12m+7m / TP
    amplitude ratio at the TP's lowest point within the matching range.
    """
    window = tapered_window(mask, pix, settings["window_taper_arcsec"])
    lo, hi, n = settings["bins_klambda"]
    edges = np.geomspace(lo, hi, int(n) + 1)
    centres = np.sqrt(edges[:-1] * edges[1:])
    p_high = windowed_power(interferometer_at_tp, window, pix, edges)
    p_tp = windowed_power(tp, window, pix, edges)
    band = sorted(float(x) for x in ARCSEC_PER_KLAMBDA / np.asarray(settings["matching_range_arcsec"], dtype=float))
    match = psd_match_factor(p_high, p_tp, centres, band)
    return {"factor": match["factor"], "at_klambda": match["at_klambda"], "matching_range_klambda": band,
            "window_area_arcsec2": float(np.sum(mask) * pix**2), "bin_centres_klambda": centres.tolist(),
            "power_interferometer_at_tp_resolution": p_high.tolist(), "power_tp": p_tp.tolist(),
            "power_unit": "MJy2 sr-2 arcsec2"}


TP_FACTOR_VARIANTS = {
    "prepared": "Power-spectrum measurement on the prepared inputs used for the combination: 12m+7m only where "
                "the gain is defined, smoothing normalised at the edge, TP regridded with the task's spline.",
}


def run_tp_factor(paths: Mapping[str, Path], config: Mapping, overwrite: bool = False) -> dict[str, Path]:
    """Power-spectrum-matched TP factor measured on the prepared inputs; documents the applied factor."""
    settings = config["tp_factor"]
    root = paths["output_root"] / "checks" / "tp_factor"
    out = {"summary": root / "tp_factor.json",
           **{f"figure_{f}_{v}_{ext}": root / f"tp_factor_{f}_{v}.{ext}"
              for f in [settings["determined_from"]] for v in TP_FACTOR_VARIANTS for ext in ("pdf", "png")}}
    guard_outputs(list(out.values()), overwrite)
    if settings["determined_from"] not in TP_FLAVOURS:
        raise ValueError(f"determined_from must be one of {TP_FLAVOURS}")
    require_prepared(paths)
    tp_fwhm = config["prepare"]["tp_fwhm_arcsec"]
    clean, header, wcs = read_plane(paths["interferometer"])
    high_beam, pix = beam_parameters(header), _pixel_arcsec(wcs)
    clean = np.nan_to_num(jybeam_to_mjysr(clean, *high_beam[:2])).astype(np.float32)
    defined = fits.getdata(prepared_paths(paths["output_root"])["gain"]) > 0
    smoothed = np.nan_to_num(smooth_to(np.where(defined, clean, np.nan), pix, high_beam,
                                       tp_fwhm)).astype(np.float32)
    del clean
    summary = {"value_combined": settings["value"], "settings": dict(settings), "variants": TP_FACTOR_VARIANTS,
               "note": "The factor in value is applied to every TP model. This stage measures the "
                       "power-spectrum-matched factor on the prepared inputs; it is not applied."}
    flavour = settings["determined_from"]
    summary[flavour] = {}
    for variant in TP_FACTOR_VARIANTS:
        m = measure_tp_factor(load_prepared(paths, flavour)[0], smoothed, defined, pix, settings)
        summary[flavour][variant] = m
        plot_tp_factor(root / f"tp_factor_{flavour}_{variant}", m["bin_centres_klambda"],
                       m["power_interferometer_at_tp_resolution"], m["power_tp"], m["matching_range_klambda"],
                       settings["value"], m["factor"])
        print(flavour, variant, "factor", round(m["factor"], 3), "at", round(m["at_klambda"], 2), "klambda",
              flush=True)
    write_json(out["summary"], summary)
    return out


# ---------------------------------------------------------------- combination and its measurements

def combine_settings(config: Mapping, method: str, tp_factor: float) -> dict:
    settings = {"tp_factor": float(tp_factor), "padding_arcsec": config["prepare"]["padding_arcsec"]}
    if method == "jcomb":
        section = config["jcomb"]
        settings.update(k1=section["k1_klambda"], k2=section["k2_klambda"], order=section["butterworth_order"],
                        minimum_tp_response=section["minimum_tp_response"])
    return settings


def pb_taper(pb: np.ndarray, pb0: float, pb1: float) -> np.ndarray:
    """Cosine weight in the primary beam: 0 for pb <= pb0 (and NaN), 1 for pb >= pb1, float32."""
    if not 0 < pb0 < pb1 <= 1:
        raise ValueError(f"The taper needs 0 < pb0 < pb1 <= 1, got ({pb0}, {pb1})")
    pb = np.asarray(pb, dtype=np.float32)
    t = np.clip((pb - np.float32(pb0)) / np.float32(pb1 - pb0), 0.0, 1.0)
    weight = np.float32(0.5) * (1.0 - np.cos(np.float32(np.pi) * t))
    return np.where(np.isfinite(pb), weight, 0.0).astype(np.float32)


def product_metadata(config: Mapping, method: str, flavour: str, settings: Mapping,
                     tp_frequency_ghz: float, fill_outside_mosaic: bool | None = None) -> dict:
    """Header settings (COMB_*) of one combined product; ``fill_outside_mosaic`` overrides the configured value."""
    fill = config["combine"]["fill_outside_mosaic"] if fill_outside_mosaic is None else fill_outside_mosaic
    return {"method": method, "tp_flavour": flavour, "tp_frequency_ghz": tp_frequency_ghz,
            "working_unit": "MJy/sr", "fill_outside_mosaic": bool(fill),
            "padding_mode": PADDING_MODE, "absolute_scale": config["absolute_scale"]["status"], **settings}


def product_blank(gain, tp_footprint, fill_outside_mosaic: bool) -> np.ndarray:
    """Pixels left blank in a product.

    Always outside the TP model. With ``fill_outside_mosaic`` the product keeps
    the TP-only region; otherwise it is also blank where the mosaic gain is undefined.
    """
    blank = ~np.asarray(tp_footprint, dtype=bool)
    if not fill_outside_mosaic:
        blank |= np.asarray(gain) == 0
    return blank


def catalogue_positions(px, py, gain, min_gain) -> list[tuple[float, float]]:
    """Catalogue pixel positions whose nearest pixel is on the grid with gain >= ``min_gain``."""
    ny, nx = np.shape(gain)
    out = []
    for x, y in zip(np.asarray(px, dtype=float), np.asarray(py, dtype=float)):
        ix, iy = int(round(x)), int(round(y))
        if 0 <= ix < nx and 0 <= iy < ny and gain[iy, ix] >= min_gain:
            out.append((float(x), float(y)))
    return out


def clean_reference(truesky: np.ndarray, gain: np.ndarray) -> np.ndarray:
    """The delivered true-sky image with 0 outside the prepared mosaic (gain 0), so the pb cut applies."""
    return np.where(gain > 0, np.nan_to_num(truesky), 0.0)


def measurement_context(paths: Mapping[str, Path], config: Mapping, flavour: str) -> dict:
    """Arrays and reference quantities shared by the measurements of one TP flavour."""
    from astropy.table import Table
    low, high, gain, header, pix, high_beam, low_beam = load_prepared(paths, flavour)
    m = config["measure"]
    clean = clean_reference(read_plane(paths["interferometer"])[0], gain)
    x0, y0, x1, y1 = m["noise_box_pix"]
    box = clean[y0:y1, x0:x1]
    sigma = float(1.4826 * np.median(np.abs(box - np.median(box))))
    table = Table.read(paths["core_catalogue"])
    px, py = WCS(header).celestial.world_to_pixel_values(np.asarray(table["WCS_ACOOR"]),
                                                         np.asarray(table["WCS_DCOOR"]))
    edges = np.asarray(m["profile_edges_klambda"], dtype=float)
    footprint = fits.getdata(prepared_paths(paths["output_root"])[f"tp_{flavour}_footprint"]) > 0.5
    tp_native, tp_header, tp_wcs = read_plane(paths[f"tp_{flavour}"])
    tp_native = jybeam_to_mjysr(tp_native.astype(float), *beam_parameters(tp_header)[:2])
    return {"tp_native": tp_native, "tp_header": tp_header, "tp_wcs": tp_wcs, "low": low, "high": high, "gain": gain, "header": header, "pix": pix, "high_beam": high_beam,
            "low_beam": low_beam, "clean": clean, "sigma": sigma,
            "positions": catalogue_positions(px, py, gain, m["min_gain"]), "mask": gain >= m["min_gain"],
            "edges": edges, "centres": np.sqrt(edges[:-1] * edges[1:]),
            "window": mosaic_window(gain, m["min_gain"], pix, m["window_smooth_arcsec"]),
            "tp_high": tp_at_high_resolution(low, pix, low_beam, high_beam, config["jcomb"]["k2_klambda"],
                                             config["prepare"]["padding_arcsec"]),
            "footprint": footprint, "blank": product_blank(gain, footprint, bool(config["combine"]["fill_outside_mosaic"]))}


def measure_product(product: np.ndarray, ctx: Mapping, config: Mapping, tp_factor: float) -> dict:
    """Measurements of one combined image (Jy/beam) against its two inputs; see the task README."""
    m = config["measure"]
    surface = jybeam_to_mjysr(np.nan_to_num(product), ctx["high_beam"][0], ctx["high_beam"][1])
    to_interferometer = windowed_amplitude_ratio(surface, ctx["high"], ctx["window"], ctx["pix"], ctx["edges"])
    k2 = config["jcomb"]["k2_klambda"]
    with np.errstate(invalid="ignore", divide="ignore"):
        to_tp = windowed_amplitude_ratio(surface, tp_factor * ctx["tp_high"], ctx["window"], ctx["pix"], ctx["edges"])
    to_tp = [float(r) if hi <= k2 else None for r, hi in zip(to_tp, ctx["edges"][1:])]   # TP term ends at k2
    join = (ctx["edges"][:-1] >= k2) & (ctx["edges"][1:] <= m["high_k_join_max_klambda"])
    large = large_scale_ratio(product, tp_factor * ctx["low"], ctx["pix"], ctx["high_beam"],
                              config["prepare"]["tp_fwhm_arcsec"], ctx["mask"])
    cores = core_contrast_ratio(product, ctx["clean"], ctx["positions"], ctx["pix"], m["core_annulus_arcsec"])
    high_k = float(np.median(to_interferometer[join]))
    return {
        "tp_factor": tp_factor,
        "large_scale_ratio_to_scaled_tp": large,
        "tp_only_ratio_to_scaled_tp": tp_only_ratio(product, tp_factor * ctx["low"], ctx["high_beam"],
                                                    ~ctx["blank"] & (ctx["gain"] == 0)),
        "high_k_join": high_k,
        "core_contrast": cores,
        "negative_fraction": negative_fraction(product, ctx["sigma"], ctx["mask"], m["negative_sigma"]),
        "negative_fraction_interferometer": negative_fraction(ctx["clean"], ctx["sigma"], ctx["mask"],
                                                              m["negative_sigma"]),
        "sigma_jybeam": ctx["sigma"],
        "uv_bins_klambda": ctx["centres"].tolist(),
        "amplitude_ratio_to_interferometer": to_interferometer.tolist(),
        "amplitude_ratio_to_scaled_tp_at_high_beam": to_tp,
        "power": {"interferometer": windowed_power(ctx["high"], ctx["window"], ctx["pix"], ctx["edges"]).tolist(),
                  "scaled_tp": windowed_power(tp_factor * ctx["low"], ctx["window"], ctx["pix"],
                                              ctx["edges"]).tolist(),
                  "combined": windowed_power(surface, ctx["window"], ctx["pix"], ctx["edges"]).tolist(),
                  "unit": "MJy2 sr-2 arcsec2"},
        "checks": {"large_scale_ratio": ratio_check(large["median"], m["large_scale_tolerance"]),
                   "high_k_join": ratio_check(high_k, m["high_k_join_tolerance"]),
                   "core_contrast": ratio_check(cores["median"], m["core_contrast_tolerance"])},
    }


def _feather_beam_outside(pix, high_beam, low_beam, size=2048) -> float:
    u, v, _ = uv_distance_klambda((size, size), pix)
    g_low, g_high = gaussian_transfer(u, v, *low_beam), gaussian_transfer(u, v, *high_beam)
    return beam_volume_outside(g_low + (1.0 - g_low) * g_high, g_high, (size, size))


def zero_spacing(product, ctx: Mapping, config: Mapping, tp_factor: float):
    """Smooth a combined image to the TP beam and compare it with the TP on the TP's own grid.

    Returns (summary, combined image at the TP beam on the TP grid, ratio map, mosaic mask,
    TP-only mask). Only positions with a full TP beam inside the product are compared.
    """
    from reproject import reproject_interp
    tp_fwhm, pix, wcs = config["prepare"]["tp_fwhm_arcsec"], ctx["pix"], WCS(ctx["header"]).celestial
    surface = jybeam_to_mjysr(np.where(np.isfinite(product), product, np.nan), *ctx["high_beam"][:2])
    smooth = smooth_to(surface, pix, ctx["high_beam"], tp_fwhm)
    size = 2 * int(round(tp_fwhm / pix)) + 1
    inside = ndimage.uniform_filter((~ctx["blank"]).astype(np.float32), size=size, mode="constant") > 0.999
    shape = ctx["tp_native"].shape
    on_tp = lambda data, order: reproject_interp((data, wcs), ctx["tp_wcs"], shape_out=shape, order=order)[0]
    combined_tp = on_tp(np.where(inside, smooth, np.nan), "bilinear")
    gain_tp = np.nan_to_num(on_tp(ctx["gain"].astype(float), "nearest-neighbor"))
    inside_tp = np.nan_to_num(on_tp(inside.astype(float), "nearest-neighbor")) > 0.5
    tp = tp_factor * ctx["tp_native"]
    valid = inside_tp & np.isfinite(combined_tp) & np.isfinite(tp)
    mosaic, tp_only = valid & (gain_tp >= config["measure"]["min_gain"]), valid & (gain_tp == 0)
    summary = zero_spacing_check(combined_tp, tp, {"mosaic": mosaic, "tp_only": tp_only, "all": valid})
    summary["note"] = ("Combined image smoothed to the TP beam, sampled on the TP grid, against the TP times the "
                       "TP factor. flux_ratio is the ratio of summed brightness; 1 means the zero-spacing "
                       "(total-flux) information of the TP is kept.")
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(valid, combined_tp / tp, np.nan)
    return summary, np.where(valid, combined_tp, np.nan), ratio, mosaic, tp_only


def _run_extended(method, flavour, ctx, standard, paths, config, settings, out, frequency) -> dict:
    """Write the extended product of one flavour and its zero-spacing check; return its measurements.

    The 12m+7m plane is weighted by the primary-beam taper, then combined with the same TP and
    settings; only the TP footprint blanks it. Runs after the standard product's work, so that the
    full-size planes of the two combinations are never held together.
    """
    factor, (pb0, pb1) = settings["tp_factor"], config["combine"]["extended_taper_pb"]
    tp_fwhm = config["prepare"]["tp_fwhm_arcsec"]
    ctx = extended_context(ctx, paths, flavour)     # the extended TP model, when one is declared
    tp_name = f"tp_{flavour}_extended" if f"tp_{flavour}_extended" in paths else f"tp_{flavour}"
    hashes = input_hashes(paths, ["interferometer", "gain", "pb", tp_name])
    pb = read_gain(paths["pb"], ctx["high"].shape)
    high = ctx["high"] * pb_taper(pb, pb0, pb1)
    blank = product_blank(ctx["gain"], ctx["footprint"], True)
    product = combine_prepared(ctx["low"], high, ctx["pix"], ctx["high_beam"], ctx["low_beam"], method, settings,
                               blank=blank)
    del high
    meta = product_metadata(config, method, flavour, {**settings, "extended": True, "taper_pb0": pb0, "taper_pb1": pb1},
                            frequency, fill_outside_mosaic=True)
    write_product(out[f"image_extended_{flavour}"], product, ctx["header"], meta, hashes, beam=ctx["high_beam"])
    ctx = {**ctx, "blank": blank}
    check, at_tp_beam, ratio, _, _ = zero_spacing(product, ctx, config, factor)
    info = {"taper_pb": [pb0, pb1], "zero_spacing": check,
            "tp_only_ratio_to_scaled_tp": tp_only_ratio(product, factor * ctx["low"], ctx["high_beam"],
                                                        ~blank & (ctx["gain"] == 0))}
    threshold = max(0.35, pb1)
    keep = (pb >= threshold) & np.isfinite(standard) & (standard > 3.0 * ctx["sigma"])
    info["ratio_to_standard_threshold_pb"] = threshold
    info["ratio_to_standard_pb_ge_0.35"] = float(np.median(product[keep] / standard[keep])) if keep.any() else None
    tp_meta = {"stage": method, "tp_flavour": flavour, "tp_factor": factor, "extended": True}
    write_product(out[f"at_tp_beam_extended_{flavour}"], at_tp_beam, ctx["tp_header"],
                  {**tp_meta, "quantity": "extended image smoothed to the TP beam, on the TP grid"}, hashes,
                  bunit="MJy/sr", beam=(tp_fwhm, tp_fwhm, 0.0))
    write_product(out[f"ratio_extended_{flavour}"], ratio, ctx["tp_header"],
                  {**tp_meta, "quantity": "extended image at the TP beam / (TP factor x TP)"}, hashes, bunit="",
                  beam=(tp_fwhm, tp_fwhm, 0.0))
    return info


def run_method(method: str, paths: Mapping[str, Path], config: Mapping, overwrite: bool = False) -> dict[str, Path]:
    """Combined images of one method for both TP models, each with the zero-spacing check
    (combined image at the TP beam against the TP), measurements and figures."""
    output_root = paths["output_root"]
    minimum_pb, taper_pb0 = interferometer_min_pb(config), config["combine"]["extended_taper_pb"][0]
    if taper_pb0 < minimum_pb:
        raise ValueError(f"[combine] extended_taper_pb starts at {taper_pb0}, below [prepare] interferometer_min_pb "
                         f"{minimum_pb}: the taper would act on blank data")
    root, zero_root = output_root / "checks" / method, output_root / "checks" / "zero_spacing"
    frequency = tp_frequency_ghz(paths)
    label = frequency_label(frequency)
    out = {"measurements": root / "measurements.json"}
    for flavour in TP_FLAVOURS:
        kind = FLAVOUR_KIND[flavour]
        out[f"image_{flavour}"] = paths["products"][f"{method}_{kind}"]
        out[f"at_tp_beam_{flavour}"] = zero_root / f"{method}_{kind}_{label}_at_tp_beam.fits"
        out[f"ratio_{flavour}"] = zero_root / f"zero_spacing_ratio_{method}_{kind}_{label}.fits"
        out[f"image_extended_{flavour}"] = paths["products"][f"{method}_{kind}_extended"]
        out[f"at_tp_beam_extended_{flavour}"] = zero_root / f"{method}_{kind}_{label}_extended_at_tp_beam.fits"
        out[f"ratio_extended_{flavour}"] = zero_root / f"zero_spacing_ratio_{method}_{kind}_{label}_extended.fits"
        for name in ("comparison", "zero_spacing", "power_spectrum"):
            for ext in ("pdf", "png"):
                out[f"{name}_{flavour}_{ext}"] = root / f"{name}_{flavour}.{ext}"
    guard_outputs(list(out.values()), overwrite)
    factor = float(config["tp_factor"]["value"])
    settings = combine_settings(config, method, factor)
    measurements = {"method": method, "settings": settings, "absolute_scale": config["absolute_scale"]["status"]}
    tp_fwhm = config["prepare"]["tp_fwhm_arcsec"]
    k1, k2 = config["jcomb"]["k1_klambda"], config["jcomb"]["k2_klambda"]
    for flavour in TP_FLAVOURS:
        ctx = measurement_context(paths, config, flavour)
        hashes = input_hashes(paths, ["interferometer", "gain", f"tp_{flavour}"])
        product = combine_prepared(ctx["low"], ctx["high"], ctx["pix"], ctx["high_beam"], ctx["low_beam"], method,
                                   settings, blank=ctx["blank"])
        write_product(out[f"image_{flavour}"], product, ctx["header"],
                      product_metadata(config, method, flavour, settings, frequency), hashes,
                      beam=ctx["high_beam"])
        m = measure_product(product, ctx, config, factor)
        check, at_tp_beam, ratio, mosaic, tp_only = zero_spacing(product, ctx, config, factor)
        m["zero_spacing"] = check
        m["checks"]["zero_spacing_flux_in_mosaic"] = ratio_check(check["mosaic"]["flux_ratio"],
                                                                 config["measure"]["large_scale_tolerance"])
        if method == "feather":
            m["beam_volume_outside_high_beam"] = _feather_beam_outside(ctx["pix"], ctx["high_beam"], ctx["low_beam"])
        measurements[flavour] = m
        meta = {"stage": method, "tp_flavour": flavour, "tp_factor": factor}
        write_product(out[f"at_tp_beam_{flavour}"], at_tp_beam, ctx["tp_header"],
                      {**meta, "quantity": "combined image smoothed to the TP beam, on the TP grid"}, hashes,
                      bunit="MJy/sr", beam=(tp_fwhm, tp_fwhm, 0.0))
        write_product(out[f"ratio_{flavour}"], ratio, ctx["tp_header"],
                      {**meta, "quantity": "combined image at the TP beam / (TP factor x TP)"}, hashes, bunit="",
                      beam=(tp_fwhm, tp_fwhm, 0.0))
        tp_label = TP_LABELS[flavour] + ("" if factor == 1.0 else f" × {factor:g}")
        plot_comparison(root / f"comparison_{flavour}", ctx["clean"],
                        mjysr_to_jybeam(factor * ctx["low"], *ctx["high_beam"][:2]), product, ctx["gain"] > 0,
                        METHOD_LABELS[method], tp_label)
        plot_zero_spacing(root / f"zero_spacing_{flavour}", factor * ctx["tp_native"], at_tp_beam, mosaic, tp_only,
                          METHOD_LABELS[method], tp_label, check)
        plot_power_spectrum(root / f"power_spectrum_{flavour}", ctx["centres"],
                            {"interferometer": m["power"]["interferometer"], "tp": m["power"]["scaled_tp"],
                             "combined": m["power"]["combined"]},
                            {"interferometer": m["amplitude_ratio_to_interferometer"],
                             "tp": m["amplitude_ratio_to_scaled_tp_at_high_beam"]},
                            k2, (k1, k2) if method == "jcomb" else (), METHOD_LABELS[method], tp_label)
        print(method, flavour, {k: round(c["value"], 4) for k, c in m["checks"].items()},
              "zero spacing", {r: round(v["flux_ratio"], 4) for r, v in check.items() if isinstance(v, dict) and v.get("n")},
              flush=True)
        if method == "jcomb" and flavour == TP_FLAVOURS[-1]:
            v = config["measure"]
            measurements["mock_sky_test"] = synthetic_gate(
                ctx["gain"], ctx["pix"], ctx["high_beam"], ctx["low_beam"], combine_settings(config, method, 1.0),
                v["aca_shortest_klambda"], v["min_gain"], flux_tolerance=v["gate_flux_tolerance"],
                rms_tolerance=v["gate_rms_tolerance"])
            print("mock sky test", measurements["mock_sky_test"], flush=True)
        measurements[flavour]["extended"] = _run_extended(method, flavour, ctx, product, paths, config, settings,
                                                          out, frequency)
        del ctx, product
    measurements["flags"] = sorted(f"{flavour}/{name}" for flavour in TP_FLAVOURS
                                   for name, c in measurements[flavour]["checks"].items() if c["flag"])
    write_json(out["measurements"], measurements)
    return out


def run_jcomb(paths, config, overwrite=False):
    return run_method("jcomb", paths, config, overwrite)


def run_feather(paths, config, overwrite=False):
    return run_method("feather", paths, config, overwrite)


# ---------------------------------------------------------------- power spectra in one box

def box_estimates(tp, interferometer, defined, pix, high_beam, low_beam, tp_fwhm, core, settings, k2, padding) -> dict:
    """TP factors that align the TP with the 12m+7m image inside the ``core`` slice of a cutout.

    The cutout is larger than the box so that smoothing and the TP beam correction are
    not affected by the box edge. Three estimates per uv bin (bin width = 1 / box size):
    ``psd_tp_res``: sqrt(P_12m+7m / P_TP) with the 12m+7m image smoothed to the TP beam;
    ``amplitude_high_res`` and ``cross_high_res``: the same amplitude ratio, and the
    cross-spectrum (shared signal only), with the TP brought to the 12m+7m beam instead.
    """
    high_tp = np.nan_to_num(smooth_to(np.where(defined, interferometer, np.nan), pix, high_beam, tp_fwhm))[core]
    tp_high = tp_at_high_resolution(tp, pix, low_beam, high_beam, k2, padding)[core]
    tp_box, int_box, mask = tp[core], interferometer[core], defined[core]
    window = tapered_window(mask, pix, settings["window_taper_arcsec"])
    step = ARCSEC_PER_KLAMBDA / (min(tp_box.shape) * pix)      # bin width set by the shorter side
    edges = step * (0.5 + np.arange(0, int(20.0 / step) + 1))
    centres = 0.5 * (edges[:-1] + edges[1:])
    spectra = {"interferometer_tp_res": windowed_power(high_tp, window, pix, edges),
               "tp": windowed_power(tp_box, window, pix, edges),
               "interferometer": windowed_power(int_box, window, pix, edges),
               "tp_high_res": windowed_power(tp_high, window, pix, edges)}
    band = sorted(float(x) for x in ARCSEC_PER_KLAMBDA / np.asarray(settings["matching_range_arcsec"], dtype=float))
    in_band = (centres >= band[0]) & (centres <= band[1])
    a, b = remove_window_mean(tp_high, window), remove_window_mean(int_box, window)
    s_ab, s_aa = cross_spectrum_slopes(a, b, window, pix, edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        psd = np.sqrt(spectra["interferometer_tp_res"] / spectra["tp"])
        amplitude = np.sqrt(spectra["interferometer"] / spectra["tp_high_res"])
        cross = 1.0 / s_aa
    usable = centres <= k2                      # the TP beam correction ends at k2
    band_ab, band_aa = (float(x[0]) for x in cross_spectrum_slopes(a, b, window, pix, band))
    clean = lambda v, keep: [float(x) if ok and np.isfinite(x) else None for x, ok in zip(v, keep)]
    return {"bin_centres_klambda": centres.tolist(), "bin_width_klambda": float(step),
            "matching_range_klambda": band, "window": window, "interferometer_tp_res": high_tp,
            "spectra": {k: v.tolist() for k, v in spectra.items()},
            "factor_per_bin": {"psd_tp_res": clean(psd, np.ones(len(psd), bool)),
                               "amplitude_high_res": clean(amplitude, usable), "cross_high_res": clean(cross, usable)},
            "factor_in_matching_range": {
                "psd_tp_res_lowest_tp_point": float(np.max(psd[in_band])),
                "psd_tp_res_median": float(np.median(psd[in_band])),
                "amplitude_high_res_median": float(np.median(amplitude[in_band & usable])),
                "cross_high_res": 1.0 / band_aa,
                "correlation": float(np.sqrt(band_ab / band_aa))},
            "bins_in_matching_range": int(np.sum(in_band)),
            "fraction_of_box_with_gain_defined": float(np.mean(mask))}


def run_box_spectrum(paths: Mapping[str, Path], config: Mapping, overwrite: bool = False) -> dict[str, Path]:
    """Power spectra of the TP and the 12m+7m image inside each box, and the factor that aligns them."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from scipy import fft as _fft
    s = config["box_spectrum"]
    root = paths["output_root"] / "checks" / "box_spectrum"
    out = {f"{box['name']}_{ext}": root / f"box_{box['name']}.{ext}" for box in s["boxes"]
           for ext in ("json", "pdf", "png")}
    guard_outputs(list(out.values()), overwrite)
    flavour = config["tp_factor"]["determined_from"]
    truth_path = paths["products"][f"jcomb_{FLAVOUR_KIND[flavour]}"]
    require_inputs([truth_path], "jcomb")
    low, high, gain, header, pix, high_beam, low_beam = load_prepared(paths, flavour)
    tp_fwhm, k2 = config["prepare"]["tp_fwhm_arcsec"], config["jcomb"]["k2_klambda"]
    margin = int(np.ceil(3.0 * tp_fwhm / pix))
    truth_full = fits.getdata(truth_path)
    for box in s["boxes"]:
        centre = SkyCoord(box["centre"], unit=(u.hourangle, u.deg))
        x, y = (int(round(float(v))) for v in WCS(header).celestial.world_to_pixel(centre))
        half_x, half_y = (int(round(side / pix / 2)) for side in box["size_arcsec"])
        if y - half_y < 0 or x - half_x < 0 or y + half_y > gain.shape[0] or x + half_x > gain.shape[1]:
            raise ValueError(f"Box {box['name']} does not fit on the grid")
        # Cutout = box plus a three-beam margin for the smoothing, clipped at the grid edge.
        y0, x0 = max(y - half_y - margin, 0), max(x - half_x - margin, 0)
        y1, x1 = min(y + half_y + margin, gain.shape[0]), min(x + half_x + margin, gain.shape[1])
        cut = np.s_[y0:y1, x0:x1]
        core = np.s_[y - half_y - y0:y + half_y - y0, x - half_x - x0:x + half_x - x0]
        defined = gain[cut] > 0
        args = (pix, high_beam, low_beam, tp_fwhm, core, s, k2, config["prepare"]["padding_arcsec"])
        data = box_estimates(low[cut], high[cut], defined, *args)
        window = data.pop("window")
        interferometer_tp_res = data.pop("interferometer_tp_res")
        summary = {"name": box["name"], "tp_model": flavour, "centre": box["centre"],
                   "centre_deg": [float(centre.ra.deg), float(centre.dec.deg)], "size_arcsec": box["size_arcsec"],
                   "distance_to_grid_edge_arcsec": float(min(y - half_y, x - half_x, gain.shape[0] - y - half_y,
                                                             gain.shape[1] - x - half_x) * pix),
                   **data}
        # Known-answer check: truth = the J-comb image (unit factor); mock TP = truth at the TP beam / F.
        f_true = float(s["mock_true_factor"])
        truth = np.nan_to_num(jybeam_to_mjysr(truth_full[cut], *high_beam[:2]))
        mock_tp = np.nan_to_num(smooth_to(truth, pix, high_beam, tp_fwhm)) / f_true
        _, _, k = uv_distance_klambda(truth.shape, pix)
        cutoff = float(s["mock_missing_below_klambda"])
        filtered = _fft.irfft2(_fft.rfft2(truth) * (1.0 - 1.0 / (1.0 + (k / cutoff) ** 10)), s=truth.shape)
        summary["mock_check"] = {"true_factor": f_true}
        for name, mock_int in (("interferometer_complete", truth),
                               ("interferometer_without_large_scales", filtered)):
            mock = box_estimates(mock_tp, mock_int * defined, defined, *args)
            summary["mock_check"][name] = mock["factor_in_matching_range"]
        summary["mock_check"]["note"] = (
            "Mock TP = the J-comb image at the TP beam divided by the true factor. Mock 12m+7m = the J-comb "
            f"image inside the mosaic, as is and with scales below {cutoff:g} klambda removed.")
        write_json(out[f"{box['name']}_json"], summary)
        aligned = data["factor_in_matching_range"]["psd_tp_res_lowest_tp_point"]
        plot_box_spectrum(root / f"box_{box['name']}", np.where(defined[core], high[cut][core], np.nan),
                          low[cut][core], np.where(defined[core], interferometer_tp_res, np.nan), window,
                          (half_x * pix, half_y * pix), data["bin_centres_klambda"], data["spectra"],
                          {**data["factor_per_bin"], "aligned": aligned}, data["matching_range_klambda"], k2)
        print(box["name"], json.dumps({k: summary[k] for k in ("factor_in_matching_range",
                                                                "fraction_of_box_with_gain_defined")}
                                      | {"mock_check": summary.get("mock_check")}, indent=1), flush=True)
    return out


# ---------------------------------------------------------------- Planck calibration

def run_planck_calibration(paths: Mapping[str, Path], config: Mapping, overwrite: bool = False) -> dict[str, Path]:
    """Compare the TP model with Planck at the TP frequency: over the well-covered area (linear
    fit with a uniform background) and in the brightest Planck beam alone."""
    from astropy.wcs import WCS as _WCS
    from reproject import reproject_exact, reproject_interp
    s = config["planck_calibration"]
    root = paths["output_root"] / "checks" / "planck_calibration"
    out = {"summary": root / "planck_calibration.json", "planck": root / "planck_at_tp_frequency.fits",
           "model": root / "tp_at_planck_beam.fits", "coverage": root / "tp_coverage_of_planck_beam.fits",
           **{f"{name}_{ext}": root / f"{name}.{ext}" for name in ("planck_maps", "planck_scatter")
              for ext in ("pdf", "png")}}
    guard_outputs(list(out.values()), overwrite)
    flavour, band, pix = s["tp_model"], s["band"], s["pixel_arcsec"]
    nu_tp = read_frequency_ghz(paths[f"tp_{flavour}"]) * 1e9
    tp_fwhm = config["prepare"]["tp_fwhm_arcsec"]

    # Planck: K_CMB -> CO(2-1) removed -> MJy/sr at the TP frequency, on the working grid.
    tp, tp_header, tp_wcs = read_plane(paths[f"tp_{flavour}"])
    tp = jybeam_to_mjysr(tp.astype(float), *beam_parameters(tp_header)[:2])
    ra, dec = tp_wcs.wcs_pix2world([[tp.shape[1] / 2 - 0.5, tp.shape[0] / 2 - 0.5]], 0)[0]
    grid = pc.grid_header(float(ra), float(dec), pix, s["width_arcmin"])
    record = json.loads(paths["planck_fetch_record"].read_text())
    nside = int(record["nside"])
    co_factor = float(record["maps"]["co_type1"]["header"]["CNV2-1"])
    rimo = fits.getdata(paths["planck_rimo"], "MAP_PARAMS")
    planck_fwhm = {str(r["FREQUENCY"][0]).strip(): float(r["FWHM"]) * 60 for r in rimo}[band]
    nu, tau = pc.read_bandpass(paths["planck_rimo"], band)
    factor = pc.kcmb_to_mjysr(nu, tau, nu_tp, s["dust_beta"], s["dust_temperature_k"])
    z, zco = np.load(paths["planck_217"]), np.load(paths["planck_co"])
    t_map = pc.healpix_on_grid(z["pixel"], z["I_STOKES"], nside, grid)
    co = pc.healpix_on_grid(zco["pixel"], zco["INTEN21"], nside, grid) * co_factor
    planck, planck_with_co = (t_map - co) * factor, t_map * factor

    # TP model at the Planck beam, with the fraction of each beam inside its footprint.
    shape = (grid["NAXIS2"], grid["NAXIS1"])
    kernel_fwhm = float(np.sqrt(planck_fwhm**2 - tp_fwhm**2))
    tp_grid, _ = reproject_exact((tp, tp_wcs), grid, shape_out=shape)
    model, coverage = pc.smooth_with_coverage(tp_grid, kernel_fwhm, pix)
    # Variant: outside the footprint, the 870 um map scaled by the median TP / 870 um ratio.
    atl, atl_header = fits.getdata(paths["atlasgal_870"], header=True)
    atl = jybeam_to_mjysr(np.squeeze(atl).astype(float), atl_header["BMAJ"] * 3600, atl_header["BMIN"] * 3600)
    atl_on_tp, _ = reproject_interp((atl, _WCS(atl_header).celestial), tp_wcs, shape_out=tp.shape)
    fill_scale = float(np.nanmedian(tp / atl_on_tp))
    filled_grid, _ = reproject_exact((np.where(np.isfinite(tp), tp, atl_on_tp * fill_scale), tp_wcs), grid,
                                     shape_out=shape)
    filled, coverage_filled = pc.smooth_with_coverage(filled_grid, kernel_fwhm, pix)

    finite = np.isfinite(planck)
    selected = (coverage >= s["coverage_min"]) & finite
    beam_pixels = 1.1331 * planck_fwhm**2 / pix**2
    selections = {f"coverage_ge_{s['coverage_min']:g}": (model, planck, selected),
                  "coverage_ge_0.8": (model, planck, (coverage >= 0.8) & finite),
                  "coverage_ge_0.7": (model, planck, (coverage >= 0.7) & finite),
                  "filled_with_870um": (filled, planck, (coverage_filled >= 0.99) & finite & (coverage >= 0.5)),
                  f"no_co_removal_coverage_ge_{s['coverage_min']:g}": (model, planck_with_co, selected)}
    fits_ = {}
    for name, (m, p, sel) in selections.items():
        fits_[name] = {**pc.linear_fit(m, p, sel), "independent_beams": float(sel.sum() / beam_pixels),
                       "ratio_direct": pc.ratio_stats(m, p, sel)}
        fits_[name]["ratio_background_removed"] = pc.ratio_stats(m, p - fits_[name]["background_mjysr"], sel)
    main_fit = fits_[f"coverage_ge_{s['coverage_min']:g}"]

    # Brightest beam (southern and middle part of the filament), on its own.
    radius = s["bright_radius_arcmin"] * 60.0 / pix
    peak, region = pc.bright_region(model, selected, radius)
    threshold = np.quantile(planck[selected], s["reference_fraction"])
    reference = selected & (planck <= threshold)
    bright = pc.bright_region_calibration(model, planck, peak, region, reference, main_fit["background_mjysr"])
    bright_with_co = pc.bright_region_calibration(
        model, planck_with_co, peak, region, reference,
        fits_[f"no_co_removal_coverage_ge_{s['coverage_min']:g}"]["background_mjysr"])
    peak_ra, peak_dec = _WCS(grid).wcs_pix2world([[peak[1], peak[0]]], 0)[0]
    bright.update({"peak_ra_dec_deg": [float(peak_ra), float(peak_dec)], "radius_arcmin": s["bright_radius_arcmin"],
                   "coverage_at_peak": float(coverage[peak]), "coverage_min_in_region": float(coverage[region].min()),
                   "independent_beams": float(region.sum() / beam_pixels),
                   "co_fraction_of_planck_at_peak": float((planck_with_co[peak] - planck[peak]) / planck_with_co[peak]),
                   "without_co_removal": bright_with_co["ratio_model_over_planck"]})
    ratios = bright["ratio_model_over_planck"]["region_mean"]
    summary = {
        "tp_model": flavour, "tp_frequency_ghz": nu_tp / 1e9, "planck_band": band, "planck_fwhm_arcsec": planck_fwhm,
        "kcmb_to_mjysr_at_tp_frequency": factor,
        "kcmb_to_mjysr_iras_convention_check": pc.kcmb_to_mjysr_iras(nu, tau, float(band) * 1e9),
        "co21_kcmb_per_k_kms": co_factor, "fill_scale_tp_over_870um": fill_scale,
        "linear_fits_tp_vs_planck": fits_,
        "ratio_by_coverage": {f"{lo:.1f}-{hi:.1f}": pc.ratio_stats(model, planck - main_fit["background_mjysr"],
                                                                  (coverage >= lo) & (coverage < hi) & finite)
                              for lo, hi in ((0.5, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01))},
        "brightest_beam": bright,
        "implied_tp_factor": {
            "linear_fit_with_background (1/gain)": 1.0 / main_fit["gain"],
            "brightest_beam_direct": 1.0 / ratios["direct"],
            "brightest_beam_background_removed": 1.0 / ratios["background_removed"],
            "brightest_beam_differential": 1.0 / ratios["differential"]},
        "note": "The factor is what the TP model would be multiplied by to match Planck. Nothing is applied.",
    }
    write_json(out["summary"], summary)
    hashes = input_hashes(paths, ["planck_217", "planck_co", "planck_rimo", f"tp_{flavour}"])
    for key, data, quantity, beam in (
            ("planck", planck, f"Planck {band} GHz, CO(2-1) removed, at {nu_tp / 1e9:g} GHz", planck_fwhm),
            ("model", np.where(coverage >= 0.5, model, np.nan), f"TP {flavour} at the Planck beam", planck_fwhm),
            ("coverage", coverage, "fraction of the Planck beam inside the TP model", planck_fwhm)):
        header = grid.copy()
        header["BUNIT"] = "" if key == "coverage" else "MJy/sr"
        header["BMAJ"] = header["BMIN"] = beam / 3600.0
        header["BPA"] = 0.0
        write_product(out[key], data, header, {"stage": "planck_calibration", "quantity": quantity}, hashes,
                      bunit=header["BUNIT"], beam=(beam, beam, 0.0))
    label = TP_LABELS[flavour]
    plot_planck_maps(root / "planck_maps", planck, model, coverage, s["coverage_min"], peak, radius, pix, label)
    plot_planck_scatter(root / "planck_scatter", planck, model, selected, region, reference, main_fit, bright, label)
    print(json.dumps({"fits": {k: (round(v["gain"], 3), round(v["offset_mjysr"], 2)) for k, v in fits_.items()},
                      "brightest_beam": bright, "implied_tp_factor": summary["implied_tp_factor"]}, indent=1), flush=True)
    return out
