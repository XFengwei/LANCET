"""Render five registered print cards from a vertically expanded sky field."""
from pathlib import Path
import argparse
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
import multiwavelength_strip as strip
from lancet.infrared_mosaic import fill_saturated_cores

OUT = ROOT / 'results/multiwavelength_cards'
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
                         facecolor='white')
        ax = fig.add_axes([BORDER/(WIDTH+2*BORDER), BORDER/(HEIGHT+2*BORDER),
                           WIDTH/(WIDTH+2*BORDER), HEIGHT/(HEIGHT+2*BORDER)])
        ax.imshow(image, origin='lower', interpolation='none', aspect='auto')
        ax.set_axis_off()
        ax.set_autoscale_on(False)
        colour = panel.get('label_colour', 'white')
        ax.text(.008, 1 - .07 * (4/14*WIDTH)/HEIGHT, panel['label'],
                transform=ax.transAxes, color=colour, fontsize=font, va='top',
                path_effects=[pe.withStroke(linewidth=1.25,
                    foreground='white' if colour == 'black' else 'black')])
        strip.draw_annotations(ax, panel.get('annotations', []), shape,
                               config['grid']['pixel_scale_arcsec'], font)
        for text in ax.texts[1:]:
            text.set_path_effects(outline)
        # Preserve the reference line weights in proportion to the image width.
        for line in ax.lines:
            line.set_linewidth(line.get_linewidth() * WIDTH/(6.7*25.4))
        for patch in ax.patches:
            patch.set_linewidth(patch.get_linewidth() * WIDTH/(6.7*25.4))
        if name == 'nir':
            origin = np.array([shape[1]-170., 150.])
            for label, vector in strip.compass_vectors(wcs, shape).items():
                ax.annotate('', xy=origin+100*vector, xytext=origin,
                            arrowprops=dict(arrowstyle='-|>,head_length=0.35,head_width=0.17',
                                            mutation_scale=5, color='white', lw=.6,
                                            shrinkA=0, shrinkB=0,
                                            joinstyle='miter', capstyle='butt'))
                ax.text(*(origin+135*vector), label, color='white', fontsize=font,
                        ha='center', va='center', path_effects=outline)
            ax.plot(*origin, marker='o', markersize=2, color='white')
        if name == 'radio':
            length = config['figure']['scalebar_length_deg']*3600/config['grid']['pixel_scale_arcsec']
            x1, y = shape[1]-60., 70.
            ax.plot([x1-length, x1], [y, y], color='white', lw=1.25)
            ax.text(x1-length/2, y+20, config['figure']['scalebar_label'],
                    color='white', fontsize=font, ha='center', va='bottom',
                    path_effects=outline)
        prefix = OUT / f'{i:02d}_{name}_card'
        for fmt in ('png', 'tiff', 'pdf'):
            kwargs = {'pil_kwargs': {'compression': 'tiff_lzw'}} if fmt == 'tiff' else {}
            fig.savefig(prefix.with_suffix('.'+fmt), dpi=DPI, facecolor='white', **kwargs)
        fig.canvas.draw()
        bounds = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        report['panels'][name] = {'image_bounds_mm': [v*25.4 for v in bounds.bounds],
                                 'output_prefix': prefix.name}
        fig.savefig(OUT/f'{i:02d}_{name}_preview.png', dpi=200, facecolor='white')
        previews.append((panel['label'], OUT/f'{i:02d}_{name}_preview.png'))
        plt.close(fig)
        print(f'Exported {prefix.name}', flush=True)
    fig, axes = plt.subplots(3, 2, figsize=(10, 10), facecolor='#eeeeee')
    for ax, (_, path) in zip(axes.flat, previews):
        ax.imshow(plt.imread(path))
        ax.set_axis_off()
    axes.flat[-1].set_axis_off()
    fig.subplots_adjust(left=.02, right=.98, top=.98, bottom=.02, hspace=.04, wspace=.04)
    fig.savefig(OUT/'contact_sheet.png', dpi=180)
    plt.close(fig)
    (OUT/'geometry.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot-only', action='store_true')
    parser.add_argument('--refresh-mm', action='store_true')
    args = parser.parse_args()
    config = strip.load_toml(ROOT / strip.DEFAULT_CONFIG)
    next(p for p in config['panels'] if p['name'] == 'radio')['label'] = 'MeerKAT  21.7 cm'
    config['grid']['height_arcmin'] = 14 * HEIGHT / WIDTH
    config['outputs']['root'] = str(OUT.relative_to(ROOT))
    mm_panel = next(p for p in config['panels'] if p['name'] == 'mm')
    OUT.mkdir(parents=True, exist_ok=True)
    if args.refresh_mm:
        strip.run_reproject(strip.ProjectPaths(ROOT), {**config, 'panels': [mm_panel]})
    elif not args.plot_only:
        strip.run_reproject(strip.ProjectPaths(ROOT), config)
    (OUT/'grid_settings.json').write_text(json.dumps(config['grid'], indent=2)+'\n')
    render(config)
