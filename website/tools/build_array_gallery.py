"""Five continuum products on the multiwavelength-strip field, for the Observations page.

ACA 7 m, ALMA 12 m and 12 m + 7 m are the delivered images; the total-power
panel is the synthesize_continuum TP model (factor 1.0, as used in the
combination); the fully combined panel is the strip's millimetre grid
(delivered J-comb product, TP model beyond the mosaic). All five use the
strip's millimetre display: its colour map and square-root stretch, in MJy/sr.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import tomllib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import PowerNorm
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from reproject import reproject_interp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from lancet.imaging import read_mjysr

spec = spec_from_file_location('strip', ROOT/'tasks/multiwavelength_strip/run.py')
strip = module_from_spec(spec)
spec.loader.exec_module(strip)
strip_config = tomllib.loads((ROOT/'tasks/multiwavelength_strip/task.toml').read_text())
combine_inputs = tomllib.loads((ROOT/'tasks/combine_continuum/task.toml').read_text())['inputs']
mm = next(panel for panel in strip_config['panels'] if panel['name'] == 'mm')
assert mm['style'] == 'sqrt'
wcs, shape = strip.strip_wcs(strip_config['grid'])
BACKGROUND = '#080b10'

PRODUCTS = [
    ('aca', ROOT/combine_inputs['aca_pbcor'], 'ACA 7 m'),
    ('12m', ROOT/combine_inputs['twelve_m_pbcor'], 'ALMA 12 m'),
    ('combined', ROOT/combine_inputs['interferometer'], 'ALMA 12 m + ACA 7 m'),
    ('tp', ROOT/'results/synthesize_continuum/model/continuum_224GHz_fill36.fits', 'Total-power model'),
    ('full', ROOT/'results/multiwavelength_strip/grid/mm_i.fits', 'Fully combined 12 m + 7 m + TP'),
]
TICKS = [0, 100, 300, 600, 1200]


def on_strip_grid(path):
    data, header = read_mjysr(path)
    if header.get('BUNIT', '').strip().lower() != 'mjy/sr':
        raise ValueError(f'{path}: expected MJy/sr after conversion')
    data = np.squeeze(data)
    if data.shape == shape and WCS(header).celestial.wcs.compare(wcs.wcs):
        return data
    image, _ = reproject_interp((data, WCS(header).celestial), wcs, shape_out=shape)
    return image


def render(name, data, label):
    fig = plt.figure(figsize=(14, 4), dpi=160, facecolor=BACKGROUND)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(BACKGROUND)
    cmap = matplotlib.colormaps[mm['cmap']].with_extremes(bad=BACKGROUND)
    # PowerNorm(0.5) is the strip's sqrt_scale: sqrt((I - vmin) / (vmax - vmin)), clipped.
    norm = PowerNorm(0.5, vmin=mm['vmin'], vmax=mm['vmax'], clip=True)
    plotted = ax.imshow(np.ma.masked_invalid(data), origin='lower', cmap=cmap, norm=norm,
                        interpolation='nearest')
    ax.set_axis_off()
    outline = [pe.withStroke(linewidth=2, foreground='black')]
    ax.text(.02, .95, label + ' · 1.3 mm', transform=ax.transAxes, color='white', va='top',
            fontsize=14, path_effects=outline)
    origin = np.array([shape[1]*.94, shape[0]*.72])
    for direction, vector in strip.compass_vectors(wcs, shape).items():
        ax.annotate('', xy=origin + 100*vector, xytext=origin,
                    arrowprops=dict(arrowstyle='->', color='white', lw=1.5))
        ax.text(*(origin + 130*vector), direction, color='white', ha='center', va='center',
                path_effects=outline)
    ax.plot(*origin, marker='o', markersize=3, color='white', markeredgewidth=0, zorder=5)
    bar_ax = fig.add_axes([.03, .15, .30, .025])
    bar = fig.colorbar(plotted, cax=bar_ax, orientation='horizontal')
    bar.set_ticks(TICKS, labels=[str(t) for t in TICKS])
    bar.ax.tick_params(colors='white', labelsize=8)
    bar.set_label('Surface brightness (MJy/sr)', color='white', fontsize=9)
    bar.outline.set_edgecolor('white')
    figure = strip_config['figure']
    length = figure['scalebar_length_deg']*3600/strip_config['grid']['pixel_scale_arcsec']
    x, y = shape[1]*.96, shape[0]*.09
    ax.plot([x - length, x], [y, y], color='white', lw=2)
    ax.text(x - length/2, y + 25, figure['scalebar_label'], ha='center', color='white',
            path_effects=outline)
    output = ROOT/f'website/assets/array-{name}.jpg'
    fig.savefig(output, dpi=160, facecolor=BACKGROUND, pil_kwargs={'quality': 88})
    plt.close(fig)
    finite = data[np.isfinite(data)]
    print(f'{output.name}: coverage {np.isfinite(data).mean():.3f}, '
          f'above {mm["vmax"]:.0f} MJy/sr {np.mean(finite > mm["vmax"]):.4f}')


for name, path, label in PRODUCTS:
    render(name, on_strip_grid(path), label)
