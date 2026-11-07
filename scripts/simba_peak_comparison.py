"""Compare bright SIMBA and model peaks on the existing common grid.

Run from the project root in the py312 environment. Inputs are read-only.
"""
from pathlib import Path
import argparse
import csv
import json
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.ndimage import maximum_filter

root = Path.cwd()
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--overwrite", action="store_true", help="Replace this task's existing outputs")
args = parser.parse_args()
config = tomllib.loads((root / "tasks/synthesize_continuum/task.toml").read_text())
out = root / config["outputs"]["root"] / "figures"
out.mkdir(parents=True, exist_ok=True)
validation = root / config["outputs"]["root"] / "validation"
outputs = [validation / "simba_peak_measurements.csv", validation / "simba_peak_measurements.json"]
outputs += [out / f"simba_bright_peak_comparison.{ext}" for ext in ("pdf", "svg", "png", "tiff")]
existing = [p for p in outputs if p.exists()]
if existing and not args.overwrite:
    raise FileExistsError(f"Output exists: {existing[0]}; use --overwrite to regenerate")
simba, header = fits.getdata(validation / "simba_250GHz_on_model_grid.fits", header=True)
model, model_header = fits.getdata(validation / "continuum_250GHz_smoothed_24arcsec.fits", header=True)
assert simba.shape == model.shape
assert all(np.isclose(header[k], model_header[k]) for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "BMAJ", "BMIN"))
assert np.isclose(abs(header["CDELT1"]) * 3600, 6.0)
assert np.isclose(header["BMAJ"] * 3600, config["beam"]["simba_fwhm_arcsec"])
assert header["BUNIT"] == model_header["BUNIT"] == "Jy/beam"
raw, raw_header = fits.getdata(root / config["inputs"]["simba_250ghz"], header=True)
valid_raw = (np.isfinite(raw) & (raw > config["simba"]["invalid_floor"])).astype(float)
valid_map, footprint = reproject_interp((valid_raw, WCS(raw_header).celestial),
                                        WCS(header).celestial, shape_out=simba.shape,
                                        order="nearest-neighbor")
valid = np.isfinite(simba) & np.isfinite(model) & (footprint > 0) & (valid_map > 0.5)
roi = np.zeros(simba.shape, dtype=bool)
lo, hi = map(int, config["plot"]["crop_pixels"])
roi[lo:hi, lo:hi] = True  # Same displayed region as the existing comparison.
peaks = (valid & roi & (simba >= 1.0) &
         (simba == maximum_filter(np.where(valid, simba, -np.inf), size=9)))
ys, xs = np.where(peaks)
order = np.argsort(simba[ys, xs])[::-1]
yy, xx = np.indices(simba.shape)
wcs = WCS(header).celestial
rows = []
for rank, idx in enumerate(order, 1):
    y, x = int(ys[idx]), int(xs[idx])
    radius = np.hypot(xx - x, yy - y)
    annulus = (radius >= 6) & (radius < 10) & valid
    if np.count_nonzero(annulus) < 0.8 * np.count_nonzero((radius >= 6) & (radius < 10)):
        raise ValueError(f"Insufficient valid annulus coverage for peak {rank}")
    ra, dec = wcs.all_pix2world(x, y, 0)
    s, m = float(simba[y, x]), float(model[y, x])
    sb, mb = float(np.median(simba[annulus])), float(np.median(model[annulus]))
    rows.append(dict(id=rank, x=x, y=y, ra_deg=float(ra), dec_deg=float(dec),
                     simba_peak=s, model_at_simba_peak=m, simba_background=sb,
                     model_background=mb, simba_contrast=s-sb, model_contrast=m-mb,
                     background_pixels=int(annulus.sum())))
assert len(rows) == 4, f"Expected four peaks from inspected selection, found {len(rows)}"
with (validation / "simba_peak_measurements.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader(); writer.writerows(rows)
metadata = {"frequency_ghz": 250, "common_beam_fwhm_arcsec": config["beam"]["simba_fwhm_arcsec"],
            "pixel_arcsec": 6, "simba_calibration_factor": config["simba"]["flux_scale"],
            "selection": f"SIMBA peaks >= 1 Jy/beam, 9x9 pixel maximum filter, inside x/y [{lo},{hi}), valid original-map footprint",
            "calibration_provenance": "0.20 is the observational SIMBA flux calibration, confirmed by the author on 2026-10-31.",
            "background": "Median in 36--60 arcsec annulus at each SIMBA-selected peak; same positions and masks for both maps",
            "model": "dust and free-free extrapolation, smoothed from 18 to 24 arcsec",
            "limits": "Four bright peaks; no independent source or absolute calibration uncertainty estimated; local annuli can include filtered emission, negative bowls and neighbouring sources.",
            "rows": rows}
(validation / "simba_peak_measurements.json").write_text(json.dumps(metadata, indent=2) + "\n")

plt.rcParams.update({"font.family":"sans-serif", "pdf.fonttype":42, "svg.fonttype":"none", "font.size":8,
                     "axes.labelsize":8, "xtick.labelsize":7, "ytick.labelsize":7,
                     "legend.fontsize":7, "axes.linewidth":0.7})
fig, ax = plt.subplots(figsize=(3.50, 3.25))
fig.subplots_adjust(left=.17, right=.97, bottom=.16, top=.95)
ax.plot([0,6], [0,6], color="0.55", lw=.8, ls="--", zorder=0, label="1:1")
for p in rows:
    ax.plot([p['simba_peak'], p['simba_contrast']],
            [p['model_at_simba_peak'], p['model_contrast']], color="0.72", lw=.8, zorder=1)
ax.scatter([p['simba_peak'] for p in rows], [p['model_at_simba_peak'] for p in rows],
           s=33, color="#007b96", zorder=3, label="Peak brightness")
ax.scatter([p['simba_contrast'] for p in rows], [p['model_contrast'] for p in rows],
           s=32, facecolors="white", edgecolors="#d17925", marker="s", linewidths=1,
           zorder=3, label="Peak minus local background")
ax.set(xlim=(0,6), ylim=(0,6), xlabel=r"SIMBA (Jy beam$^{-1}$)",
       ylabel=r"Reconstruction (Jy beam$^{-1}$)")
ax.set_aspect("equal", adjustable="box")
ax.legend(loc="upper left", frameon=False, handlelength=1.3, borderpad=.25)
fig.savefig(out / "simba_bright_peak_comparison.pdf")
fig.savefig(out / "simba_bright_peak_comparison.svg")
fig.savefig(out / "simba_bright_peak_comparison.png", dpi=300)
fig.savefig(out / "simba_bright_peak_comparison.tiff", dpi=300)
plt.close(fig)
print(json.dumps(metadata, indent=2))
