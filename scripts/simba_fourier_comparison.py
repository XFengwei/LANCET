"""Compare radial Fourier amplitudes of the beam-matched SIMBA and model maps.

Run from the repository root after the synthesize_continuum validation stage.
The Fourier samples are image-derived diagnostics, not measured visibilities.
"""
from pathlib import Path
import argparse
import csv
import json
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.ndimage import distance_transform_edt
from scipy.signal.windows import tukey

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--overwrite", action="store_true")
args = parser.parse_args()
root = Path.cwd()
cfg = tomllib.loads((root / "tasks/synthesize_continuum/task.toml").read_text())
out = root / cfg["outputs"]["root"]
figdir, valdir = out / "figures", out / "validation"
figdir.mkdir(parents=True, exist_ok=True)
targets = [valdir / "simba_fourier_comparison.csv", valdir / "simba_fourier_comparison.json"]
targets += [figdir / f"simba_fourier_comparison.{ext}" for ext in ("pdf", "svg", "png", "tiff")]
if not args.overwrite and any(p.exists() for p in targets):
    raise FileExistsError("Fourier comparison output exists; use --overwrite to regenerate")

simba, header = fits.getdata(valdir / "simba_250GHz_on_model_grid.fits", header=True)
model, mh = fits.getdata(valdir / "continuum_250GHz_smoothed_24arcsec.fits", header=True)
assert simba.shape == model.shape
assert all(np.isclose(header[k], mh[k]) for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "BMAJ", "BMIN"))
assert header["BUNIT"] == mh["BUNIT"] == "Jy/beam"
assert np.isclose(abs(header["CDELT1"])*3600, 6.0)
assert np.isclose(header["BMAJ"]*3600, 24.0) and np.isclose(header["BMIN"]*3600, 24.0)
raw, rh = fits.getdata(root / cfg["inputs"]["simba_250ghz"], header=True)
mask_input = (np.isfinite(raw) & (raw > cfg["simba"]["invalid_floor"])).astype(float)
mask, footprint = reproject_interp((mask_input, WCS(rh).celestial), WCS(header).celestial,
                                   shape_out=simba.shape, order="nearest-neighbor")
valid = np.isfinite(simba) & np.isfinite(model) & (footprint > 0) & (mask > .5)
lo, hi = map(int, cfg["plot"]["crop_pixels"])
cut = np.s_[lo:hi, lo:hi]
s, m, valid = simba[cut], model[cut], valid[cut]
n = hi - lo
distance = distance_transform_edt(valid)
mask_taper = np.sin(.5*np.pi*np.clip(distance/4, 0, 1))**2
window = np.outer(tukey(n, .25), tukey(n, .25)) * mask_taper * valid
assert window.sum() > 0
beam_pixels = np.pi * 24**2 / (4*np.log(2)*6**2)
fourier_s = np.fft.fftshift(np.fft.fft2(np.where(valid, s, 0)*window)) / beam_pixels
fourier_m = np.fft.fftshift(np.fft.fft2(np.where(valid, m, 0)*window)) / beam_pixels
freq = np.fft.fftshift(np.fft.fftfreq(n, d=6/206265)) / 1000
u, v = np.meshgrid(freq, freq)
distance_uv = np.hypot(u, v)
edges = np.geomspace(206.265/(n*6), 7.0, 15)
rows = []
for low, high in zip(edges[:-1], edges[1:]):
    select = (distance_uv >= low) & (distance_uv < high)
    if select.sum() < 4:
        continue
    fs, fm = fourier_s[select], fourier_m[select]
    a, b = float(np.mean(np.abs(fs))), float(np.mean(np.abs(fm)))
    cross = np.sum(fs*np.conjugate(fm))
    rows.append(dict(uv_klambda=float(np.sqrt(low*high)), lo_klambda=float(low),
                     hi_klambda=float(high), angular_scale_arcsec=float(206.265/np.sqrt(low*high)),
                     simba_mean_amplitude_jy=a, model_mean_amplitude_jy=b,
                     amplitude_ratio=a/b, cross_transfer=float(cross.real/np.sum(abs(fm)**2)),
                     coherence=float(abs(cross)**2/(np.sum(abs(fs)**2)*np.sum(abs(fm)**2))),
                     fourier_grid_cells=int(select.sum())))
assert all(r["uv_klambda"] > 0 and r["simba_mean_amplitude_jy"] > 0 and r["model_mean_amplitude_jy"] > 0 for r in rows)
with (valdir / "simba_fourier_comparison.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys()); writer.writeheader(); writer.writerows(rows)
beam_half = 206.265 * (4*np.log(2)**2)**.5 / (np.pi*24)
metadata = {"crop_pixels": [lo, hi], "pixel_arcsec": 6, "beam_arcsec": 24,
            "simba_flux_calibration": cfg["simba"]["flux_scale"],
            "window": "Identical 25% Tukey edge window and four-pixel sine-squared taper of common valid footprint; no mean subtraction or rescaling.",
            "valid_fraction": float(valid.mean()), "frequency_definition": "radial Fourier frequency in k-lambda",
            "amplitude": "arithmetic mean modulus of 2D FFT after Jy/beam to Jy/pixel conversion; no noise-bias correction or beam deconvolution",
            "cross_transfer": "Re(sum(F_SIMBA * conj(F_model))) / sum(|F_model|^2), descriptive and window-dependent",
            "limitations": "Finite window couples modes; +/- modes are conjugate; SIMBA noise adds positive modulus bias. Beam and map filtering are not deconvolved; no uncertainties estimated.",
            "beam_50_percent_klambda": float(beam_half), "rows": rows}
(valdir / "simba_fourier_comparison.json").write_text(json.dumps(metadata, indent=2)+"\n")

plt.rcParams.update({"font.family":"sans-serif", "pdf.fonttype":42, "svg.fonttype":"none",
                     "font.size":8, "axes.labelsize":8, "xtick.labelsize":8,
                     "ytick.labelsize":8, "legend.fontsize":8, "axes.linewidth":.7})
fig, ax = plt.subplots(figsize=(3.50, 3.15))
fig.subplots_adjust(left=.18, right=.97, bottom=.16, top=.78)
ks = np.array([r['uv_klambda'] for r in rows])
ax.loglog(ks, [r['simba_mean_amplitude_jy'] for r in rows], color="#d67924", marker="o", markersize=3.2, lw=1.1, label="SIMBA")
ax.loglog(ks, [r['model_mean_amplitude_jy'] for r in rows], color="#007b96", marker="s", markersize=3.0, lw=1.1, label="Reconstruction")
ax.axvline(beam_half, color="0.5", lw=.8, ls=":", label="Beam: 50% response")
ax.set_xlim(edges[0], edges[-1]); ax.set_xlabel(r"Spatial frequency ($k\lambda$)")
ax.set_ylabel("Mean Fourier amplitude (Jy)")
ax.xaxis.set_major_locator(FixedLocator([.3, .5, 1, 2, 4, 6]))
ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.legend(loc="lower left", frameon=False, handlelength=1.5)
top = ax.twiny(); top.set_xscale("log"); top.set_xlim(ax.get_xlim())
scales = np.array([600, 300, 150, 75, 40])
top.set_xticks(206.265/scales, [str(x) for x in scales])
top.xaxis.set_minor_locator(NullLocator())
top.set_xlabel(r"Angular scale (arcsec)", labelpad=5)
fig.savefig(figdir / "simba_fourier_comparison.pdf")
fig.savefig(figdir / "simba_fourier_comparison.svg")
fig.savefig(figdir / "simba_fourier_comparison.png", dpi=300)
fig.savefig(figdir / "simba_fourier_comparison.tiff", dpi=300)
plt.close(fig)
print(json.dumps(metadata, indent=2))
