"""Render five registered print cards from a vertically expanded sky field."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from matplotlib.patches import Polygon as PolygonPatch
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'src')]
import importlib.util
_spec = importlib.util.spec_from_file_location(
    'multiwavelength_strip', ROOT / 'tasks' / 'multiwavelength_strip' / 'run.py')
strip = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(strip)
from lancet.infrared_mosaic import fill_saturated_cores

OUT = Path(__file__).resolve().parent
BASE = ROOT / 'results/multiwavelength_strip'
WIDTH, HEIGHT, BORDER, DPI = 85.6, 54.0, 5.0, 1200


def stretched(data, reference, low, high, soft):
    lo, hi = np.nanpercentile(reference, [low, high])
    value = np.clip((data - lo) / (hi - lo), 0, 1)
    return np.arcsinh(value / soft) / np.arcsinh(1 / soft)


def rgb(panel, project):
    name = panel['name']
    def read(key, base=None):
        base = OUT if base is None else base
        return fits.getdata(strip.grid_path(base, name, key))
    if panel['kind'] == 'vvv':
        settings = strip.load_toml(ROOT / panel['vvv_config'])['rgb']
        channels = []
        for colour, weight in zip(('red', 'green', 'blue'), settings['weights']):
            key = settings[colour]
            reference = read(key, BASE)
            sample = reference[np.isfinite(reference)]
            sample = sample[::max(1, sample.size // 5_000_000)]
            _, sky, noise = sigma_clipped_stats(sample, sigma=3, maxiters=5)
            black = sky - settings['black_sigma'] * noise
            white = np.percentile(sample, settings['high_percentile'])
            data = fill_saturated_cores(np.nan_to_num(read(key), nan=sky), white)
            value = np.clip((data - black) / (white - black), 0, 1)
            soft = settings['asinh_softening']
            value = np.arcsinh(value / soft) / np.arcsinh(1 / soft)
            channels.append(np.clip(value * weight, 0, 1))
        return strip.boost_saturation(np.stack(channels, axis=-1), settings['saturation'])
    if panel['kind'] == 'rgb':
        channels = [np.clip(stretched(read(key), read(key, BASE),
                    panel['low_percentile'][i], panel['high_percentile'][i],
                    panel['asinh_softening'][i]) * panel['weights'][i], 0, 1)
                    for i, key in enumerate(('r', 'g', 'b'))]
        return np.nan_to_num(np.stack(channels, axis=-1), nan=0)
    data = read('i')
    if panel.get('style') == 'sqrt':
        value = strip.sqrt_scale(data, panel['vmin'], panel['vmax'])
        return matplotlib.colormaps[panel['cmap']].with_extremes(bad='black')(
            np.ma.masked_invalid(value))[..., :3]
    if panel.get('style') == 'inferno':
        norm = matplotlib.colors.AsinhNorm(linear_width=100., vmin=20., vmax=1500., clip=True)
        return matplotlib.colormaps['inferno'](norm(data))[..., :3]
    if panel.get('style') == 'dual_colour':
        return strip.dual_colour_rgb(data, panel['background_vmin'],
                                    panel['transition'], panel['highlight_vmax'])
    value = stretched(data, read('i', BASE), panel['low_percentile'],
                      panel['high_percentile'], panel['asinh_softening'])
    return matplotlib.colormaps[panel['cmap']].with_extremes(bad='black')(
        np.ma.masked_invalid(value))[..., :3]


def render(config):
    project = strip.ProjectPaths(ROOT)
    wcs, shape = strip.strip_wcs(config['grid'])
    font = config['figure']['label_fontsize'] * WIDTH / (6.7 * 25.4)
    outline = [pe.withStroke(linewidth=1.25, foreground='black')]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    report = {'image_mm': [WIDTH, HEIGHT], 'border_mm': BORDER,
              'total_mm': [WIDTH + 2 * BORDER, HEIGHT + 2 * BORDER],
              'dpi': DPI, 'grid_shape': shape, 'panels': {}}
    previews = []
    for i, panel in enumerate(config['panels'], 1):
        name = panel['name']
        image = rgb(panel, project)
        fig = plt.figure(figsize=((WIDTH+2*BORDER)/25.4, (HEIGHT+2*BORDER)/25.4),
                         facecolor='black')
        ax = fig.add_axes([BORDER/(WIDTH+2*BORDER), BORDER/(HEIGHT+2*BORDER),
                           WIDTH/(WIDTH+2*BORDER), HEIGHT/(HEIGHT+2*BORDER)])
        ax.imshow(image, origin='lower', interpolation='none', aspect='auto')
        ax.set_axis_off()
        ax.set_autoscale_on(False)
        colour = panel.get('label_colour', 'white')
        ax.text(.008, 1 - .07 * (4/14*WIDTH)/HEIGHT, panel['label'],
                transform=ax.transAxes, color=colour, fontsize=font*1.5, va='top',
                path_effects=[pe.withStroke(linewidth=1.25,
                    foreground='white' if colour == 'black' else 'black')])
        if panel.get('subtitle'):
            ax.text(.008, .905, panel['subtitle'], transform=ax.transAxes,
                    color='white', fontsize=5.2, va='top', path_effects=outline)
        strip.draw_annotations(ax, panel.get('annotations', []), shape,
                               config['grid']['pixel_scale_arcsec'], font)
        for text in ax.texts[1:]:
            text.set_path_effects(outline)
        # Preserve the reference line weights in proportion to the image width.
        for line in ax.lines:
            line.set_linewidth(line.get_linewidth() * WIDTH/(6.7*25.4))
        for patch in ax.patches:
            patch.set_linewidth(patch.get_linewidth() * WIDTH/(6.7*25.4))
        if name == 'nir' or panel.get('compass', False):
            origin = np.array([shape[1]-170., shape[0]-240.])
            pixels_per_point = shape[1] / (WIDTH/25.4*72)
            shaft_half = .3*pixels_per_point
            head_length, head_half = 1.75*pixels_per_point, .85*pixels_per_point
            parts = [Point(*origin).buffer(.6*pixels_per_point)]
            for label, vector in strip.compass_vectors(wcs, shape).items():
                side = np.array([-vector[1], vector[0]])
                tip = origin+100*vector
                base = tip-head_length*vector
                parts.append(Polygon([origin+shaft_half*side, base+shaft_half*side,
                    base+head_half*side, tip, base-head_half*side,
                    base-shaft_half*side, origin-shaft_half*side]))
                ax.text(*(origin+135*vector), label, color='white', fontsize=font,
                        ha='center', va='center', path_effects=outline)
            compass = unary_union(parts)
            assert compass.geom_type == 'Polygon'
            ax.add_patch(PolygonPatch(np.asarray(compass.exterior.coords), facecolor='white',
                edgecolor='white', linewidth=0, zorder=5,
                path_effects=[pe.withStroke(linewidth=.65, foreground='black')]))
            ax.plot(*origin, marker='o', markersize=1.2, markeredgewidth=0,
                    color='white', zorder=6)
        ax.text(.982, .070, 'https://xfengwei.github.io/LANCET', transform=ax.transAxes,
                color='white', fontsize=5.5, ha='right', va='bottom', path_effects=outline)
        ax.text(.982, .025, 'fengweilookuper@gmail.com', transform=ax.transAxes,
                color='white', fontsize=5.5, ha='right', va='bottom', path_effects=outline)
        if name == 'radio' or panel.get('scalebar', False):
            length = config['figure']['scalebar_length_deg']*3600/config['grid']['pixel_scale_arcsec']
            x1, y = length+60., 70.
            ax.plot([x1-length, x1], [y, y], color='white', lw=1.25)
            ax.text(x1-length/2, y+20, config['figure']['scalebar_label'],
                    color='white', fontsize=font, ha='center', va='bottom',
                    path_effects=outline)
        prefix = OUT / f'{i:02d}_{name}_card'
        for fmt in ('png', 'tiff', 'pdf'):
            kwargs = {'pil_kwargs': {'compression': 'tiff_lzw'}} if fmt == 'tiff' else {}
            fig.savefig(prefix.with_suffix('.'+fmt), dpi=DPI, facecolor='black', **kwargs)
        fig.canvas.draw()
        bounds = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        report['panels'][name] = {'image_bounds_mm': [v*25.4 for v in bounds.bounds],
                                 'output_prefix': prefix.name}
        fig.savefig(OUT/f'{i:02d}_{name}_preview.png', dpi=200, facecolor='black')
        previews.append((panel['label'], OUT/f'{i:02d}_{name}_preview.png'))
        plt.close(fig)
        print(f'Exported {prefix.name}', flush=True)
    if len(previews) <= 2:
        rows, cols = 1, len(previews)
    elif len(previews) <= 4:
        rows, cols = 2, 2
    else:
        rows, cols = 3, 2
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 5*rows),
                             squeeze=False, facecolor='#eeeeee')
    for ax, (_, path) in zip(axes.flat, previews):
        ax.imshow(plt.imread(path))
        ax.set_axis_off()
    for ax in list(axes.flat)[len(previews):]:
        ax.set_axis_off()
    fig.subplots_adjust(left=.02, right=.98, top=.98, bottom=.02, hspace=.04, wspace=.04)
    fig.savefig(OUT/'contact_sheet.png', dpi=180)
    plt.close(fig)
    (OUT/'geometry.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    config = strip.load_toml(ROOT / strip.DEFAULT_CONFIG)
    config['grid'] = json.loads((OUT / 'grid_settings.json').read_text())
    config['outputs']['root'] = str(OUT.relative_to(ROOT))
    mm = next(p for p in config['panels'] if p['name'] == 'mm')
    # Same inferno view as Fig. 1 and Fig. 3: square root from 0 to 1200 MJy/sr.
    mm.update(style='sqrt', cmap='inferno', vmin=0.0, vmax=1200.0, label_colour='white')
    next(p for p in config['panels'] if p['name'] == 'radio')['label'] = 'MeerKAT  21.7 cm'
    next(p for p in config['panels'] if p['name'] == 'mir')['annotations'] = json.loads((OUT / 'annotations.json').read_text())
    header = fits.getheader(OUT / 'grid/mm_i.fits')
    cached_wcs = strip.WCS(header).celestial
    cached_shape = (header['NAXIS2'], header['NAXIS1'])
    strip.strip_wcs = lambda grid: (cached_wcs, cached_shape)
    render(config)
