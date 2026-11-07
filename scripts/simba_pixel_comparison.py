"""Plot positive paired pixels in the common SIMBA/model comparison field."""
from pathlib import Path
import json
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp

root = Path(__file__).resolve().parent.parent
cfg = tomllib.loads((root / "tasks/synthesize_continuum/task.toml").read_text())
out = root / cfg["outputs"]["root"]
simba, sh = fits.getdata(out / "validation/simba_250GHz_on_model_grid.fits", header=True)
model, mh = fits.getdata(out / "validation/continuum_250GHz_smoothed_24arcsec.fits", header=True)
assert simba.shape == model.shape
assert all(np.isclose(sh[k], mh[k]) for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "BMAJ", "BMIN"))
assert sh["BUNIT"] == mh["BUNIT"] == "Jy/beam"
raw, rh = fits.getdata(root / cfg["inputs"]["simba_250ghz"], header=True)
raw_valid = (np.isfinite(raw) & (raw > cfg["simba"]["invalid_floor"])).astype(float)
mask, footprint = reproject_interp((raw_valid, WCS(rh).celestial), WCS(sh).celestial,
                                   shape_out=simba.shape, order="nearest-neighbor")
lo, hi = map(int, cfg["plot"]["crop_pixels"])
cut = np.s_[lo:hi, lo:hi]
s, m, mask, footprint = simba[cut], model[cut], mask[cut], footprint[cut]
valid = np.isfinite(s) & np.isfinite(m) & (footprint > 0) & (mask > .5)
positive = valid & (s > .25) & (m > 0)
x, y = s[positive], m[positive]
assert len(x) > 0
# Color encodes the reconstructed intensity at each pixel.
brightness = y
limit_min, limit_max = .1, 10.0

plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
                     "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(4.3, 3.5))
fig.subplots_adjust(left=.16, right=.79, bottom=.17, top=.97)
ax.set_box_aspect(1)
points = ax.scatter(x, y, s=8, c=brightness, cmap="viridis", norm=LogNorm(vmin=.1, vmax=10),
                    alpha=.55, linewidths=0, rasterized=True)
colorbar_axis = ax.inset_axes([1.04, 0, .06, 1])
colorbar = fig.colorbar(points, cax=colorbar_axis)
colorbar.set_label(r"Reconstructed intensity (Jy beam$^{-1}$)")
ax.plot([limit_min, limit_max], [limit_min, limit_max], color="0.35", ls="--", lw=.8, label="1:1")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(limit_min, limit_max); ax.set_ylim(limit_min, limit_max)
ax.set_xlabel(r"SIMBA (Jy beam$^{-1}$)")
ax.set_ylabel(r"Reconstruction (Jy beam$^{-1}$)")
ax.legend(frameon=False, loc="upper left")
for ext in ("pdf", "svg", "png"):
    fig.savefig(out / f"figures/simba_pixel_comparison.{ext}", dpi=300,
                bbox_inches="tight", pad_inches=.02)
plt.close(fig)
info = {"common_valid_pixels": int(valid.sum()), "simba_above_threshold_pixels": int((valid & (s > .25)).sum()),
        "selected_positive_pairs": int(positive.sum()), "simba_threshold_jy_per_beam": .25, "crop_pixels": [lo, hi],
        "inside_plot_limits": int(((x >= limit_min) & (x <= limit_max) & (y >= limit_min) & (y <= limit_max)).sum()),
        "axis_limits_jy_per_beam": [limit_min, limit_max],
        "beam_arcsec": 24, "pixel_arcsec": 6,
        "marker_size": "constant", "marker_color": "logarithmic reconstructed intensity, clipped to color scale 0.1--10 Jy/beam",
        "note": "Pixels are correlated within the 24 arcsec beam; positive selection is not a significance threshold."}
(out / "validation/simba_pixel_comparison.json").write_text(json.dumps(info, indent=2) + "\n")
print(json.dumps(info, indent=2))
