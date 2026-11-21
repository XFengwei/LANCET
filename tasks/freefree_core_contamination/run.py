#!/usr/bin/env python
"""Free-free (CORNISH) contamination of the cores and fragments (Paper II R7).

Stages: ``contamination`` (smooth radio emission in each core footprint),
``plot`` (overlay and flux comparison), ``mass_bias`` (fragment and core
mass-ratio scenarios), ``all``.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping
import json
from pathlib import Path
import sys
import tomllib

import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from astropy.wcs import WCS

from lancet import ProjectPaths


DEFAULT_CONFIG = Path("tasks/freefree_core_contamination/task.toml")
PRODUCTS = {
    "contamination": "summary.json",
    "plot": "cornish_core_comparison.png",
    "mass_bias": "mass_efficiency_bias.json",
}


def load_task_config(path: Path) -> dict[str, object]:
    with Path(path).open("rb") as handle:
        return tomllib.load(handle)


def resolve_task_paths(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> dict[str, Path]:
    paths = {
        name: (project.root / str(value)).resolve()
        for name, value in config["inputs"].items()
    }
    paths["output_root"] = project.validate_output(str(config["outputs"]["root"]))
    for name, path in paths.items():
        if name != "output_root" and not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    return paths


def read_footprints(path: Path) -> dict[int, object]:
    """getsf footprint ellipses keyed by getsf source number."""
    from regions import Regions

    footprints = {}
    for line in Path(path).read_text().splitlines():
        line = line.lstrip("# ")
        if line.startswith("ellipse("):
            shape, number = line.split("#")
            footprints[int(number.strip())] = Regions.parse(
                "fk5\n" + shape, format="ds9"
            )[0]
    return footprints


def spectral_factor(settings: Mapping[str, object]) -> float:
    """Optically thin free-free scaling from the radio to the mm frequency."""
    return (
        float(settings["frequency_ghz"]) / float(settings["radio_frequency_ghz"])
    ) ** float(settings["spectral_index"])


def _check_output(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )


def build_contamination(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, object]:
    """Smooth CORNISH free-free emission within the getsf core footprints."""
    from astropy.stats import SigmaClip
    from astropy.wcs.utils import proj_plane_pixel_area
    from photutils.background import Background2D, StdBackgroundRMS
    from scipy.ndimage import map_coordinates, uniform_filter

    out = paths["output_root"]
    _check_output(out / PRODUCTS["contamination"], overwrite)
    out.mkdir(parents=True, exist_ok=True)
    settings = config["core_contamination"]
    radio, rh = fits.getdata(paths["cornish_map"], header=True)
    assert radio.ndim == 2 and rh["BUNIT"].upper() == "JY/BEAM"
    noise_box = int(round(float(settings["noise_box_beams"]) * rh["BMAJ"] / abs(rh["CDELT1"])))
    noise_box += 1 - noise_box % 2
    background = Background2D(
        radio, box_size=noise_box, filter_size=3,
        sigma_clip=SigmaClip(sigma=3, maxiters=10),
        bkgrms_estimator=StdBackgroundRMS(),
    )
    thermal_rms = background.background_rms
    # Negative excursions retain broad bowls/sidelobes that symmetric clipping
    # could reject as outliers. For symmetric noise their conditional RMS is sigma.
    negative = radio < 0
    negative_power = uniform_filter(np.where(negative, radio**2, 0.0), size=noise_box, mode="reflect")
    negative_fraction = uniform_filter(negative.astype(float), size=noise_box, mode="reflect")
    negative_rms = np.sqrt(np.maximum(negative_power, 0) / np.maximum(negative_fraction, 1e-12))
    noise = np.maximum(thermal_rms, negative_rms)
    assert np.all(np.isfinite(noise)) and np.all(noise > 0)
    rms_mjy = float(settings["adopted_radio_rms_mjy_per_beam"])
    threshold = rms_mjy * float(settings["detection_sigma"]) / 1000
    limit_brightness = rms_mjy * float(settings["upper_limit_sigma"]) / 1000
    robust = np.isfinite(radio) & (radio > threshold)
    retained_radio = np.where(robust, radio, 0.0)
    for name, values in [
        ("radio_noise_jybeam", noise),
        ("radio_clipped_rms_jybeam", thermal_rms),
        ("radio_negative_rms_jybeam", negative_rms),
        ("radio_retained_uniform3sigma_jybeam", retained_radio),
    ]:
        fits.writeto(out / f"{name}.fits", values.astype("float32"), rh, overwrite=True)
    mask_header = rh.copy()
    mask_header["BUNIT"] = "1"
    fits.writeto(out / "radio_positive_uniform3sigma_mask.fits", robust.astype("uint8"), mask_header, overwrite=True)

    rw = WCS(rh).celestial
    beam_area = np.pi * rh["BMAJ"] * rh["BMIN"] / (4 * np.log(2))  # deg2
    catalogue = Table.read(paths["source_catalogue"])
    alma, ah = fits.getdata(paths["continuum_image"], header=True)
    alma = np.squeeze(alma)
    assert ah["BUNIT"].lower() == "jy/beam"
    aw = WCS(ah).celestial
    native_alma_beam_area = np.pi * ah["BMAJ"] * ah["BMIN"] / (4 * np.log(2))
    pixel_area = abs(proj_plane_pixel_area(aw))
    apertures = read_footprints(paths["getsf_footprints"])
    factor = spectral_factor(settings)
    alma_beam_area = np.pi * (float(settings["alma_beam_arcsec"]) / 3600) ** 2 / (4 * np.log(2))
    minimum_coverage = float(settings["minimum_coverage"])
    rows = []
    for source in catalogue:
        aperture = apertures[int(source["NO"])].to_pixel(aw)
        weights = aperture.to_mask(mode="subpixels", subpixels=8)
        assert abs(weights.data.sum() / aperture.area - 1) < 0.01
        box = weights.bbox
        y, x = np.mgrid[box.iymin:box.iymax, box.ixmin:box.ixmax]
        sky = aw.pixel_to_world(x, y)
        rx, ry = rw.world_to_pixel(sky)
        values = map_coordinates(radio, [ry, rx], order=1, mode="constant", cval=np.nan, prefilter=False)
        retained_values = map_coordinates(retained_radio, [ry, rx], order=1, mode="constant", cval=np.nan, prefilter=False)
        retained_mask = map_coordinates(robust.astype(float), [ry, rx], order=0, mode="constant", cval=0, prefilter=False) > 0.5
        retained_values = np.where(retained_mask, retained_values, 0.0)
        keep = weights.data > 0
        finite = np.isfinite(values) & keep
        coverage = float(weights.data[finite].sum() / weights.data[keep].sum())
        total = (
            float(np.sum(retained_values[finite] * weights.data[finite]) * pixel_area / beam_area * factor)
            if coverage > minimum_coverage else np.nan
        )
        retained_area_fraction = float(weights.data[keep & retained_mask].sum() / weights.data[keep].sum())
        upper_limit = (
            limit_brightness * weights.data.sum() * pixel_area / beam_area * factor
            if coverage > minimum_coverage else np.nan
        )
        detected = bool(np.isfinite(total) and total > 0)
        native_values = alma[box.iymin:box.iymax, box.ixmin:box.ixmax]
        native_finite = np.isfinite(native_values) & keep
        alma_coverage = float(weights.data[native_finite].sum() / weights.data[keep].sum())
        observed = (
            float(np.sum(native_values[native_finite] * weights.data[native_finite]) * pixel_area / native_alma_beam_area)
            if alma_coverage > minimum_coverage else np.nan
        )
        native_peak = float(map_coordinates(alma, [[aperture.center.y], [aperture.center.x]], order=1, mode="constant", cval=np.nan, prefilter=False)[0])
        cx, cy = rw.world_to_pixel(aperture.to_sky(aw).center)
        radio_center = float(map_coordinates(retained_radio, [[cy], [cx]], order=1, mode="constant", cval=np.nan, prefilter=False)[0])
        if map_coordinates(robust.astype(float), [[cy], [cx]], order=0, mode="constant", cval=0)[0] < 0.5:
            radio_center = 0.0
        predicted_peak = radio_center * factor * alma_beam_area / beam_area
        rows.append({
            "ID": int(source["ID"]), "NO": int(source["NO"]),
            "Subregion": str(source["Subregion"]), "coverage": coverage,
            "alma_coverage": alma_coverage,
            "retained_area_fraction": retained_area_fraction,
            "radio_detected": detected,
            "freefree_upper_limit_jy": upper_limit if not detected else np.nan,
            "observed_12m7m_flux_jy": observed,
            "catalogue_background_subtracted_flux_jy": float(source["FXT_BST01"]),
            "smooth_freefree_flux_jy": total,
            "smooth_freefree_fraction": total / observed,
            "freefree_fraction_upper_limit": upper_limit / observed if not detected else np.nan,
            "smooth_freefree_peak_jy_per_065beam": predicted_peak,
            "smooth_freefree_peak_fraction": predicted_peak * native_alma_beam_area / alma_beam_area / native_peak,
            "aperture_area_arcsec2": float(weights.data.sum() * pixel_area * 3600**2),
        })
    table = Table(rows=rows)
    table.write(out / "smooth_core_contamination.ecsv", format="ascii.ecsv", overwrite=True)
    table.write(out / "smooth_core_contamination.csv", format="ascii.csv", overwrite=True)
    root = Path.cwd().resolve()

    def rel(path: Path) -> str:
        return str(path.relative_to(root)) if path.is_relative_to(root) else str(path)

    summary = {
        "radio_map": rel(paths["cornish_map"]), "footprints": rel(paths["getsf_footprints"]),
        "continuum_image": rel(paths["continuum_image"]),
        "adopted_radio_rms_mjy_per_beam": rms_mjy,
        "detection_sigma": settings["detection_sigma"],
        "upper_limit_sigma": settings["upper_limit_sigma"],
        "minimum_radio_mjy_per_beam": threshold * 1000,
        "native_alma_beam_arcsec": [ah["BMAJ"] * 3600, ah["BMIN"] * 3600],
        "radio_beam_arcsec": [rh["BMAJ"] * 3600, rh["BMIN"] * 3600],
        "frequency_ghz": settings["frequency_ghz"], "spectral_index": settings["spectral_index"],
        "spectral_factor": factor,
        "model": "Smooth brightness within the native CORNISH beam; optically thin extrapolation. Bilinear interpolation, no additional convolution or deconvolution.",
        "noise_box_pixels": noise_box, "noise_box_arcsec": noise_box * abs(rh["CDELT1"]) * 3600,
        "noise_method": "Maximum of sigma-clipped mesh RMS and local RMS of all negative excursions, relative to zero. No radio brightness background is subtracted.",
        "noise_range_mjybeam": [float(noise.min() * 1000), float(noise.max() * 1000)],
        "retained_native_pixels": int(robust.sum()), "native_pixels": int(radio.size),
        "measurement": "Native pixels >3 times the adopted uniform rms retained. 12m7m true-sky aperture flux denominator includes local background. For covered nondetections, nominal 1-sigma brightness limit is integrated over the full footprint under the smooth-emission assumption; not a sum of independent pixel errors. Spatial noise maps remain diagnostic.",
        "regions": {},
    }
    for code in ["Y", "I", "E"]:
        select = (table["Subregion"] == code) & np.isfinite(table["smooth_freefree_fraction"])
        t = table[select]
        summary["regions"][code] = {
            "catalogue_count": int(np.sum(table["Subregion"] == code)),
            "covered_count": len(t),
            "median_fraction": float(np.median(t["smooth_freefree_fraction"])) if len(t) else None,
            "with_retained_radio": int(sum(t["smooth_freefree_flux_jy"] > 0)),
            "max_fraction": float(max(t["smooth_freefree_fraction"])) if len(t) else None,
            "flux_weighted_fraction": float(sum(t["smooth_freefree_flux_jy"]) / sum(t["observed_12m7m_flux_jy"])) if len(t) else None,
            "count_above_10percent": int(sum(t["smooth_freefree_fraction"] > 0.1)),
            "count_above_100percent": int(sum(t["smooth_freefree_fraction"] > 1)),
        }
        limits = np.asarray(t["freefree_fraction_upper_limit"])
        limits = limits[np.isfinite(limits)]
        summary["regions"][code]["nondetections"] = {
            "count": len(limits),
            "median_fraction_limit": float(np.median(limits)) if len(limits) else None,
            "limits_above_10percent": int(sum(limits > 0.1)),
            "limits_above_50percent": int(sum(limits > 0.5)),
            "limits_above_100percent": int(sum(limits > 1)),
        }
    (out / PRODUCTS["contamination"]).write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {out / PRODUCTS['contamination']}")
    return summary


def plot_contamination(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
    qa_scripts: Path | None = None,
) -> None:
    """Overlay (a) and aperture flux comparison (b) for one subregion."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import AsinhNorm, LogNorm
    from matplotlib.lines import Line2D
    import seaborn as sns
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.nddata import Cutout2D

    out = paths["output_root"]
    _check_output(out / PRODUCTS["plot"], overwrite)
    settings = config["plot"]
    metadata = json.loads((out / PRODUCTS["contamination"]).read_text())
    table = Table.read(out / "smooth_core_contamination.ecsv")
    source = Table.read(paths["source_catalogue"])
    radio, rh = fits.getdata(paths["cornish_map"], header=True)
    noise = fits.getdata(out / "radio_noise_jybeam.fits")
    data, h = fits.getdata(paths["continuum_image"], header=True)
    pb = np.squeeze(fits.getdata(paths["primary_beam"]))
    data = np.where(np.isfinite(pb) & (pb >= float(settings["minimum_pb"])), np.squeeze(data), np.nan)
    assert h["BUNIT"].lower() == "jy/beam"
    center = SkyCoord(float(settings["center_ra_deg"]) * u.deg, float(settings["center_dec_deg"]) * u.deg)
    size = float(settings["size_arcsec"]) * u.arcsec
    cut = Cutout2D(data * 1000, center, (size, size), wcs=WCS(h).celestial, mode="partial", fill_value=np.nan)
    code = str(settings["subregion"])
    selected = table[(table["Subregion"] == code) & np.isfinite(table["smooth_freefree_fraction"])]
    nondetections = selected[~selected["radio_detected"]]
    selected = selected[selected["radio_detected"]]

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(8.2, 4.3))
    height = 0.34 * 8.2 / 4.3
    left = fig.add_axes([0.085, 0.18, 0.34, height], projection=cut.wcs)
    right = fig.add_axes([0.635, 0.18, 0.34, height])
    im = left.imshow(cut.data, origin="lower", cmap=sns.color_palette("mako", as_cmap=True),
                     norm=AsinhNorm(linear_width=0.3, vmin=-0.2, vmax=12), rasterized=True)
    left.set_facecolor(".12")
    sigma = metadata["adopted_radio_rms_mjy_per_beam"]
    levels = np.asarray(settings["contour_sigma"], dtype=float) * sigma
    contours = left.contour(radio * 1000, levels=levels, transform=left.get_transform(WCS(rh).celestial),
                            colors="#f5a442", linewidths=0.65)
    contours.set_rasterized(True)
    footprints = read_footprints(paths["getsf_footprints"])
    for row in source[source["Subregion"] == code]:
        artist = footprints[int(row["NO"])].to_pixel(cut.wcs).as_artist(facecolor="none", edgecolor="white", lw=0.5)
        artist.set_rasterized(True)
        left.add_patch(artist)
    left.set_xlim(-0.5, cut.data.shape[1] - 0.5)
    left.set_ylim(-0.5, cut.data.shape[0] - 0.5)
    left.coords[0].set_ticks(number=3)
    left.coords[1].set_ticks(number=3)
    left.set_xlabel("RA (ICRS)")
    left.set_ylabel("Dec (ICRS)")
    left.text(0.025, 0.975, "(a)", color="white", transform=left.transAxes, va="top", fontsize=10)
    contour_text = ", ".join(f"{value:g}" for value in settings["contour_sigma"])
    left.legend(handles=[Line2D([], [], color="#f5a442", lw=1, label=f"CORNISH: {contour_text} sigma")],
                loc="lower left", bbox_to_anchor=(0.02, 1.015), frameon=False, handlelength=1.3)
    cax = left.inset_axes([1.035, 0, 0.035, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"1.3 mm intensity (mJy beam$^{-1}$)", fontsize=8)
    cb.ax.tick_params(labelsize=8)
    cb.set_ticks([0, 0.5, 2, 6, 12])
    cb.set_ticklabels(["0", "0.5", "2", "6", "12"])

    x = np.asarray(selected["observed_12m7m_flux_jy"]) * 1000
    y = np.asarray(selected["smooth_freefree_flux_jy"]) * 1000
    assert np.all(x > 0) and np.all(y > 0)
    right.scatter(x, y, s=18, c="#007b96", alpha=0.8, edgecolors="none", zorder=3, label="SNR > 3")
    limit_x = np.asarray(nondetections["observed_12m7m_flux_jy"]) * 1000
    limit_y = np.asarray(nondetections["freefree_upper_limit_jy"]) * 1000
    assert np.all(limit_x > 0) and np.all(limit_y > 0)
    right.errorbar(limit_x, limit_y, yerr=0.35 * limit_y, uplims=True, fmt="o", ms=4,
                   mfc="none", mec="#007b96", ecolor="#007b96", elinewidth=0.6, capsize=1.5,
                   alpha=0.65, label="1-sigma upper limit", zorder=2)
    line = np.geomspace(0.25, max(160, x.max() * 1.2), 100)
    right.plot(line, line, "--", c=".3", lw=0.9, label="1:1")
    right.plot(line, 0.1 * line, ":", c=".5", lw=0.9, label="10%")
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlim(0.25, max(160, x.max() * 1.2))
    right.set_ylim(min(0.001, y.min() * 0.7), 80)
    right.set_xlabel("1.3 mm flux within core footprint (mJy)")
    right.set_ylabel(f"Expected free–free flux at {metadata['frequency_ghz']:g} GHz (mJy)")
    right.text(0.025, 0.975, "(b)", transform=right.transAxes, va="top", fontsize=10)
    handles, labels = right.get_legend_handles_labels()
    legend_handles = dict(zip(labels, handles))
    legend_order = ["SNR > 3", "1-sigma upper limit", "1:1", "10%"]
    right.legend([legend_handles[label] for label in legend_order], legend_order, loc="lower left", frameon=False)
    for row in selected[selected["smooth_freefree_fraction"] > 1]:
        right.annotate(f"ID {int(row['ID'])}",
                       (float(row["observed_12m7m_flux_jy"]) * 1000, float(row["smooth_freefree_flux_jy"]) * 1000),
                       xytext=(5, 3), textcoords="offset points", fontsize=8)
    if qa_scripts is not None:
        sys.path.insert(0, str(qa_scripts))
        from audit_panel_alignment import require_matplotlib_panel_alignment

        require_matplotlib_panel_alignment(
            fig, axes=[left, right], panel_ids=["a", "b"], row_groups=[["a", "b"]],
            json_out=out / "cornish_core_comparison.alignment.json", strict=True,
        )
    for fmt in settings["formats"]:
        fig.savefig(out / f"cornish_core_comparison.{fmt}", dpi=int(settings["dpi"]),
                    bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    (out / "figure_caption.txt").write_text(
        f"Left: 12m+7m true-sky 1.3 mm continuum, displayed within PB >= {settings['minimum_pb']:g}. "
        f"CORNISH contours are {contour_text} times the adopted rms of {sigma:g} mJy per native "
        "2.5-arcsec beam. White ellipses mark core footprints. Right: true-sky aperture flux, without "
        f"background subtraction, versus optically thin {metadata['frequency_ghz']:g} GHz free-free flux. "
        "Filled circles use native radio pixels above 3 sigma. Hollow circles with downward arrows show "
        f"nominal 1-sigma limits for covered nondetections, converting the {sigma:g} mJy/beam brightness "
        "limit into integrated flux over each full footprint under the smooth-emission assumption. "
        "Dashed and dotted lines show equality and 10 percent.\n"
    )

    fig, ax = plt.subplots(figsize=(4.5, 4))
    im = ax.imshow(noise * 1000, origin="lower", cmap="viridis", norm=LogNorm(vmin=0.1, vmax=12), rasterized=True)
    ax.contour(radio * 1000, levels=[3 * sigma], colors="white", linewidths=0.5)
    ax.set_xlabel("Native radio pixel x")
    ax.set_ylabel("Native radio pixel y")
    fig.colorbar(im, ax=ax, label="Local noise estimate (mJy/beam)")
    for fmt in ("pdf", "png"):
        fig.savefig(out / f"radio_noise_map.{fmt}", dpi=int(settings["dpi"]), bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out / PRODUCTS['plot']}")


def build_mass_bias(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> dict[str, object]:
    """Propagate the smooth-radio model into fragment/core mass-ratio scenarios."""
    from astropy.convolution import Gaussian2DKernel, convolve_fft
    from reproject import reproject_interp

    out = paths["output_root"]
    _check_output(out / PRODUCTS["mass_bias"], overwrite)
    metadata = json.loads((out / PRODUCTS["contamination"]).read_text())
    radio, rh = fits.getdata(paths["cornish_map"], header=True)
    threshold = metadata["minimum_radio_mjy_per_beam"] / 1000
    template = np.where(radio > threshold, radio, 0.0)
    factor = metadata["spectral_factor"]
    aca, ah = fits.getdata(paths["aca_image"], header=True)
    pb = fits.getdata(paths["aca_pb"])
    aca, pb = np.squeeze(aca), np.squeeze(pb)
    labels = fits.getdata(paths["fragment_label_mask"])
    assert aca.shape == pb.shape == labels.shape
    cat = Table.read(paths["fragment_catalogue"], format="ascii")
    selection = Table.read(paths["fragment_selection"])
    if "flux" in selection.colnames:
        selection.remove_column("flux")
    fragments = join(selection, cat[["ID", "_idx", "flux", "area_exact"]], keys="ID")
    sigma_scale = 1 / (2 * np.sqrt(2 * np.log(2)))
    pa = np.deg2rad(ah["BPA"])
    v = np.array([np.sin(pa), np.cos(pa)])
    target_cov = ((ah["BMAJ"] * 3600) ** 2 * np.outer(v, v) + (ah["BMIN"] * 3600) ** 2 * (np.eye(2) - np.outer(v, v))) * sigma_scale**2
    radio_cov = np.eye(2) * (rh["BMAJ"] * 3600 * sigma_scale) ** 2
    jac = WCS(rh).celestial.pixel_scale_matrix * 3600
    cov = np.linalg.inv(jac) @ (target_cov - radio_cov) @ np.linalg.inv(jac).T
    ev, evec = np.linalg.eigh(cov)
    assert np.all(ev > 0)
    kernel = Gaussian2DKernel(np.sqrt(ev[1]), np.sqrt(ev[0]), theta=np.arctan2(evec[1, 1], evec[0, 1]))
    kernel.normalize()
    radio_area = np.pi * (rh["BMAJ"] * 3600) * (rh["BMIN"] * 3600) / (4 * np.log(2))
    smooth = convolve_fft(template / radio_area * factor, kernel, boundary="fill", fill_value=0)
    support = convolve_fft(np.ones(radio.shape), kernel, boundary="fill", fill_value=0)
    aca_wcs = WCS(ah).celestial
    model, foot = reproject_interp((smooth, WCS(rh).celestial), aca_wcs, shape_out=aca.shape)
    coverage, _ = reproject_interp((support, WCS(rh).celestial), aca_wcs, shape_out=aca.shape)
    pixel_area = abs(ah["CDELT1"] * ah["CDELT2"]) * 3600**2
    aca_area = np.pi * (ah["BMAJ"] * 3600) * (ah["BMIN"] * 3600) / (4 * np.log(2))
    scale = float(config["fragments"]["mass_opacity_scale"])
    rows = []
    for row in fragments:
        mask = labels == row["_idx"]
        assert mask.sum() == round(row["area_exact"])
        native_flux = float(np.sum(aca[mask]) * pixel_area / aca_area)
        assert np.isclose(native_flux, row["flux"], rtol=2e-5)
        complete = bool(np.all(np.isfinite(model[mask])) and np.all(coverage[mask] > 0.99999) and np.all(foot[mask] > 0))
        ff = float(np.sum(model[mask] * pb[mask]) * pixel_area) if complete else np.nan
        rows.append({"ID": row["ID"], "Subregion": row["Subregion"], "selected": row["selected"],
                     "radio_covered": complete, "aca_flux_jy": native_flux, "freefree_flux_jy": ff,
                     "fragment_freefree_fraction": ff / native_flux,
                     "fragment_mass_msun": float(row["mass_dense"]) * scale})
    table = Table(rows=rows)
    table.write(out / "fragment_freefree_contamination.ecsv", format="ascii.ecsv", overwrite=True)
    table.write(out / "fragment_freefree_contamination.csv", format="ascii.csv", overwrite=True)
    cores = join(Table.read(paths["physical_catalogue"]), Table.read(out / "smooth_core_contamination.ecsv"),
                 keys="ID", table_names=["mass", "ff"])
    efficiency = Table.read(paths["core_formation_summary"])
    result = {
        "radio_rms_mjybeam": metadata["adopted_radio_rms_mjy_per_beam"],
        "fragment_method": "Convolve thresholded 2.5 arcsec radio surface brightness to the ACA beam, reproject, apply ACA PB and integrate true leaf masks. Catalogue fluxes independently reproduced.",
        "core_method": "Sensitivity scenario: multiply each adopted core mass by 1 minus its true-sky aperture free-free fraction, holding temperature, geometry and optical-depth factor fixed. This is not a calibrated correction of background-subtracted core photometry.",
        "regions": {},
    }
    for r in efficiency:
        f = table[(table["selected"]) & (table["Subregion"] == r["label"])]
        c = cores[cores["Subregion_mass"] == r["subregion_code"]]
        cf = np.asarray(c["smooth_freefree_fraction"])
        cf = np.where(np.isfinite(cf), cf, 0)
        ff = np.asarray(f["fragment_freefree_fraction"])
        ff = np.where(np.isfinite(ff), ff, 0)
        assert np.all((cf >= 0) & (cf < 1)) and np.all((ff >= 0) & (ff < 1))
        mc = float(sum(c["Mass_core"]))
        mf = float(sum(f["fragment_mass_msun"]))
        cm = float(sum(c["Mass_core"] * (1 - cf)))
        fm = float(sum(f["fragment_mass_msun"] * (1 - ff)))
        result["regions"][str(r["region"])] = {
            "core_mass": mc, "fragment_mass": mf,
            "radio_covered_fragments": int(sum(f["radio_covered"])), "selected_fragments": len(f),
            "core_mass_scenario": cm, "fragment_mass_scenario": fm,
            "fractional_core_mass_reduction": 1 - cm / mc, "fractional_fragment_mass_reduction": 1 - fm / mf,
            "core_cloud_before": mc / float(r["cloud_mass"]), "core_cloud_scenario": cm / float(r["cloud_mass"]),
            "fragment_cloud_before": mf / float(r["cloud_mass"]), "fragment_cloud_scenario": fm / float(r["cloud_mass"]),
            "core_fragment_before": mc / mf, "core_fragment_scenario": cm / fm,
        }
    (out / PRODUCTS["mass_bias"]).write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {out / PRODUCTS['mass_bias']}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--stage", choices=("contamination", "plot", "mass_bias", "all"), default="all")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--qa-scripts", type=Path,
                        help="Optional folder with audit_panel_alignment.py for the figure alignment check")
    args = parser.parse_args()
    project = ProjectPaths.discover()
    config = load_task_config(project.root / args.config)
    paths = resolve_task_paths(project, config)
    if args.stage in ("contamination", "all"):
        build_contamination(paths, config, overwrite=args.overwrite)
    if args.stage in ("plot", "all"):
        plot_contamination(paths, config, overwrite=args.overwrite, qa_scripts=args.qa_scripts)
    if args.stage in ("mass_bias", "all"):
        build_mass_bias(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
