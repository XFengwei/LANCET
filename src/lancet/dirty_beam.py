"""MIRIAD 12m+7m joint-imaging dirty beam (Paper II referee comment 4b)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FWHM_PER_SIGMA = np.sqrt(8*np.log(2))


def gaussian_beam(shape, pixel_arcsec, bmaj_arcsec, bmin_arcsec, bpa_deg, center):
    """Unit-peak elliptical Gaussian on the PSF grid (PA east of north, RA increasing left)."""
    y, x = np.indices(shape, dtype=float)
    dx = -(x-center[0])*pixel_arcsec
    dy = (y-center[1])*pixel_arcsec
    angle = np.deg2rad(bpa_deg)
    major = dx*np.sin(angle)+dy*np.cos(angle)
    minor = dx*np.cos(angle)-dy*np.sin(angle)
    sig_maj, sig_min = bmaj_arcsec/FWHM_PER_SIGMA, bmin_arcsec/FWHM_PER_SIGMA
    return np.exp(-0.5*((major/sig_maj)**2+(minor/sig_min)**2))


def radial_profile(image, center, pixel_arcsec, edges_arcsec):
    """Azimuthal mean and minimum in annuli around the PSF peak."""
    y, x = np.indices(image.shape, dtype=float)
    radius = np.hypot(x-center[0], y-center[1])*pixel_arcsec
    index = np.digitize(radius.ravel(), edges_arcsec)-1
    values = image.ravel()
    mean = np.full(len(edges_arcsec)-1, np.nan)
    minimum = np.full_like(mean, np.nan)
    for k in range(len(mean)):
        use = index == k
        if use.any():
            mean[k] = values[use].mean()
            minimum[k] = values[use].min()
    return mean, minimum


def cumulative_volume(image, center, pixel_arcsec, radii_arcsec):
    """Integral of the unit-peak PSF within each radius, in beam-centre-normalized arcsec^2."""
    y, x = np.indices(image.shape, dtype=float)
    radius = np.hypot(x-center[0], y-center[1])*pixel_arcsec
    order = np.argsort(radius.ravel())
    sorted_radius = radius.ravel()[order]
    running = np.cumsum(image.ravel()[order])*pixel_arcsec**2
    position = np.searchsorted(sorted_radius, radii_arcsec, side="right")-1
    return running[np.clip(position, 0, None)]


def _read_imfit(log_path):
    """Parse MIRIAD imfit object=beam output into FWHM major/minor and PA."""
    text = Path(log_path).read_text()
    result = {}
    for line in text.splitlines():
        if "Major axis (arcsec)" in line:
            result["bmaj_arcsec"] = float(line.split(":")[1].split()[0])
        elif "Minor axis (arcsec)" in line:
            result["bmin_arcsec"] = float(line.split(":")[1].split()[0])
        elif "Position angle (degrees)" in line:
            result["bpa_deg"] = float(line.split(":")[1].split()[0])
    if len(result) != 3:
        raise ValueError(f"Could not parse beam fit from {log_path}")
    return result


def run_dirty_beam(project_root, config, *, save_figure, overwrite=False):
    """Plot the effective mosaic PSF from MIRIAD mospsf with its Gaussian fit."""
    from astropy.io import fits
    from matplotlib import pyplot as plt
    from matplotlib.patches import Ellipse

    settings = config["dirty_beam"]
    root = Path(project_root)
    out = root/config["outputs"]["root"]/"dirty_beam"
    out.mkdir(parents=True, exist_ok=True)
    summary_path = out/"summary.json"
    if summary_path.exists() and not overwrite:
        raise FileExistsError("Dirty-beam outputs exist; use --overwrite")

    psf_path = root/settings["psf_fits"]
    with fits.open(psf_path) as hdul:
        psf = np.squeeze(hdul[0].data).astype(float)
        header = hdul[0].header
    pixel = abs(header["CDELT1"])*3600
    if not np.isclose(pixel, abs(header["CDELT2"])*3600):
        raise ValueError("PSF pixels must be square")
    peak_y, peak_x = np.unravel_index(np.nanargmax(psf), psf.shape)
    psf = psf/psf[peak_y, peak_x]
    center = (peak_x, peak_y)

    imfit = _read_imfit(root/settings["imfit_log"])
    paper = fits.getheader(root/settings["paper_image"])
    paper_beam = {"bmaj_arcsec": paper["BMAJ"]*3600, "bmin_arcsec": paper["BMIN"]*3600,
                  "bpa_deg": float(paper["BPA"])}
    # The restoring beam of the paper image is the reference clean beam.
    fit = paper_beam
    gauss = gaussian_beam(psf.shape, pixel, fit["bmaj_arcsec"], fit["bmin_arcsec"], fit["bpa_deg"], center)

    edges = np.arange(0, settings["profile_max_arcsec"]+pixel, pixel)
    mids = 0.5*(edges[1:]+edges[:-1])
    psf_mean, psf_min = radial_profile(psf, center, pixel, edges)
    gauss_mean, _ = radial_profile(gauss, center, pixel, edges)
    radii = mids
    psf_volume = cumulative_volume(psf, center, pixel, radii)
    gauss_volume = np.pi*fit["bmaj_arcsec"]*fit["bmin_arcsec"]/(4*np.log(2))
    first_null = float(mids[np.argmax(psf_mean < 0)]) if np.any(psf_mean < 0) else None
    outside = np.hypot(*(np.indices(psf.shape)[::-1]-np.array(center)[:, None, None]))*pixel > 2*fit["bmaj_arcsec"]
    sidelobe = {"max_positive": float(psf[outside].max()), "min_negative": float(psf[outside].min())}
    ratio_at = {f"{r:g}": float(cumulative_volume(psf, center, pixel, np.array([r]))[0]/gauss_volume)
                for r in settings["volume_ratio_radii_arcsec"]}

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
        "axes.linewidth": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    fig, axes = plt.subplots(2, 2, figsize=(settings["width_inches"], settings["width_inches"]*0.86))
    fig.subplots_adjust(left=0.09, right=0.97, bottom=0.08, top=0.97, wspace=0.32, hspace=0.30)
    a, b, c, d = axes.flat
    for axis in axes.flat:
        axis.tick_params(direction="in", top=True, right=True, length=3, width=0.7)

    def offsets(half):
        n = int(round(half/pixel))
        cut = psf[peak_y-n:peak_y+n+1, peak_x-n:peak_x+n+1]
        return cut, [n*pixel+pixel/2, -n*pixel-pixel/2, -n*pixel-pixel/2, n*pixel+pixel/2]

    zoom, extent = offsets(settings["zoom_half_arcsec"])
    image = a.imshow(zoom, origin="lower", extent=extent, cmap="RdBu_r", vmin=-0.2, vmax=1.0)
    a.contour(zoom, levels=[0.1, 0.3, 0.5, 0.7, 0.9], extent=extent, colors="k", linewidths=0.5)
    a.contour(zoom, levels=[-0.05, -0.02], extent=extent, colors="k", linewidths=0.5, linestyles="--")
    a.add_patch(Ellipse((0, 0), fit["bmin_arcsec"], fit["bmaj_arcsec"], angle=-fit["bpa_deg"],
                        fill=False, edgecolor="gold", linewidth=1.2))
    a.set(xlabel=r"$\Delta\alpha$ (arcsec)", ylabel=r"$\Delta\delta$ (arcsec)")
    fig.colorbar(image, ax=a, fraction=0.046, pad=0.03)

    wide, extent = offsets(settings["wide_half_arcsec"])
    lim = settings["wide_color_limit"]
    image = b.imshow(wide, origin="lower", extent=extent, cmap="RdBu_r", vmin=-lim, vmax=lim)
    b.set(xlabel=r"$\Delta\alpha$ (arcsec)", ylabel=r"$\Delta\delta$ (arcsec)")
    fig.colorbar(image, ax=b, fraction=0.046, pad=0.03, extend="both")

    c.plot(mids, psf_mean, color="k", linewidth=1.0, label="Dirty beam, azimuthal mean")
    c.plot(mids, psf_min, color="0.55", linewidth=0.8, linestyle=":", label="Dirty beam, annulus minimum")
    c.plot(mids, gauss_mean, color="goldenrod", linewidth=1.0, linestyle="--", label="Restoring beam")
    c.axhline(0, color="0.6", linewidth=0.5)
    c.set(xscale="log", xlim=(pixel, settings["profile_max_arcsec"]), ylim=(-0.12, 1.05),
          xlabel="Radius (arcsec)", ylabel="Normalized response")
    c.legend(loc="upper right", fontsize=7.5)

    d.plot(radii, psf_volume/gauss_volume, color="k", linewidth=1.0)
    d.axhline(1, color="goldenrod", linewidth=0.8, linestyle="--")
    d.set(xscale="log", xlim=(pixel, settings["profile_max_arcsec"]),
          xlabel="Radius (arcsec)", ylabel=r"Dirty / clean beam volume within $r$")

    for axis, label in zip(axes.flat, "abcd"):
        axis.text(0.04, 0.95, f"({label})", transform=axis.transAxes, ha="left", va="top",
                  fontsize=9, color="k", bbox={"facecolor": "white", "edgecolor": "none", "pad": 1, "alpha": 0.8})

    formats = settings["output_formats"]
    save_figure(fig, [out/f"miriad_dirty_beam.{fmt}" for fmt in formats], dpi=settings["dpi"])
    plt.close(fig)

    np.savetxt(out/"radial_profile.csv",
               np.column_stack([mids, psf_mean, psf_min, gauss_mean, psf_volume/gauss_volume]),
               delimiter=",", comments="",
               header="radius_arcsec,psf_mean,psf_min,gaussian_mean,cumulative_volume_ratio")
    summary = {"psf_fits": str(settings["psf_fits"]), "pixel_arcsec": pixel,
               "peak_pixel_xy": [int(peak_x), int(peak_y)], "imfit_beam_6arcsec_box": imfit,
               "paper_restoring_beam": paper_beam, "rerun_restor_beams": settings["rerun_restor_beams"], "first_null_arcsec": first_null,
               "sidelobes_beyond_2bmaj": sidelobe, "gaussian_volume_arcsec2": gauss_volume,
               "dirty_to_clean_volume_ratio": ratio_at}
    summary_path.write_text(json.dumps(summary, indent=2)+"\n")
    return summary
