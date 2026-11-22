"""Render the ALMA pointing-footprint figure for the dark Observations page.

Run from the repository root with the project's py312 Python environment.
Circles mark the approximate primary-beam FWHM at 225 GHz (43" for 7 m, 26" for
12 m) around each actual pointing centre; they are not sensitivity contours.

Dark-page rules (website/README.md, "Design essentials"): page background
#080b10, text #ecece8, axes #a6adb6; each array keeps its Fig. 2 dark-mode
colour (lancet.combination.coverage_figure.THEMES["dark"]): 7 m orange,
12 m teal (C2 and C5 share it, as in Fig. 2 panel a). The maps use a dark
grey stretch so the circles read on top of them.
"""
from pathlib import Path
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from astropy.visualization.wcsaxes import SphericalCircle
import astropy.units as u
from reproject import reproject_interp
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from lancet.combination.coverage_figure import THEMES

DARK = THEMES['dark']
BACKGROUND, TEXT, AXIS = DARK['background'], DARK['text'], DARK['axis']
ACA, TWELVE = DARK['arrays']['aca_7m'], DARK['arrays']['configuration_c2']
DIAGNOSTICS = ROOT/'results/combine_continuum/diagnostics'
ARRAYS = [
    (DIAGNOSTICS/'overlap_tests/aca_7m_metadata.json', 43, ACA),
    (DIAGNOSTICS/'full_sampling/configuration_c2_metadata.json', 26, TWELVE),
    (DIAGNOSTICS/'full_sampling/configuration_c5_metadata.json', 26, TWELVE),
]
MAPS = [('image_HHLI-SPIRE250.fits', 'Herschel / 250 μm'), ('image_GLIMPSE-IRAC3.fits', 'Spitzer / 5.8 μm')]
# Grey map shown between black and this level, so the filament stays mid-grey.
MAP_CEILING = 0.72


def grid_wcs():
    target = WCS(naxis=2)
    target.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    target.wcs.crval = [221.25, -59.8]
    target.wcs.crpix = [301, 301]
    target.wcs.cdelt = [-2.0/3600, 2.0/3600]
    return target


def fill_blanks(image):
    """Fill blank (saturated) pixels from their nearest valid neighbour."""
    blank = ~np.isfinite(image)
    if not blank.any():
        return image
    _, (rows, cols) = ndimage.distance_transform_edt(blank, return_indices=True)
    return image[rows, cols]


def pointings(path):
    for field in json.loads(path.read_text())['subtables']['FIELD']:
        ra, dec = np.rad2deg(np.asarray(field['PHASE_DIR']).flatten())
        yield ra % 360, dec


plt.rcParams.update({'text.color': TEXT, 'axes.labelcolor': TEXT, 'xtick.color': AXIS,
                     'ytick.color': AXIS, 'axes.edgecolor': AXIS, 'font.size': 10})
fig = plt.figure(figsize=(13, 5.4), facecolor=BACKGROUND)
wcs = grid_wcs()
for index, (name, title) in enumerate(MAPS):
    with fits.open(ROOT/'data/delivered/archive/infrared'/name) as hdus:
        hdu = next(x for x in hdus if x.data is not None and x.data.ndim >= 2)
        data, _ = reproject_interp((np.squeeze(hdu.data), WCS(hdu.header).celestial), wcs, shape_out=(600, 600))
    data = fill_blanks(data)
    norm = ImageNormalize(data, interval=PercentileInterval(99.5), stretch=AsinhStretch())
    ax = fig.add_subplot(1, 2, index + 1, projection=wcs)
    ax.set_facecolor(BACKGROUND)
    ax.imshow(MAP_CEILING*norm(data), origin='lower', cmap='gray', vmin=0, vmax=1)
    xy = []
    for path, beam_arcsec, colour in ARRAYS:
        for ra, dec in pointings(path):
            xy.append(wcs.world_to_pixel_values(ra, dec))
            ax.add_patch(SphericalCircle((ra*u.deg, dec*u.deg), beam_arcsec/7200*u.deg, edgecolor=colour,
                                         facecolor='none', lw=.7, alpha=.9, transform=ax.get_transform('world')))
    xy = np.array(xy)
    ax.set_xlim(xy[:, 0].min() - 35, xy[:, 0].max() + 35)
    ax.set_ylim(xy[:, 1].min() - 35, xy[:, 1].max() + 35)
    ax.set_title(title, pad=10, color=TEXT)
    ax.set_xlabel('Right ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.coords.grid(alpha=.10, color='white')
    for coord in ax.coords:
        coord.set_ticklabel(color=TEXT)
fig.legend([Line2D([], [], color=ACA, lw=1.6), Line2D([], [], color=TWELVE, lw=1.6)],
           ['ACA 7 m', 'ALMA 12 m (C2 + C5)'], loc='lower center', ncol=2, frameon=False)
fig.subplots_adjust(left=.08, right=.98, top=.93, bottom=.18, wspace=.30)
output = ROOT/'website/assets/footprints.jpg'
fig.savefig(output, dpi=180, facecolor=BACKGROUND, pil_kwargs={'quality': 90})
print(f'{output.name}: {fig.get_size_inches()*180}')
