"""Render unannotated wavelength panels from registered scientific FITS grids."""
from pathlib import Path
import importlib.util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('cards', ROOT/'results/multiwavelength_cards/render_cards.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.OUT = m.BASE
config = m.strip.load_toml(ROOT/m.strip.DEFAULT_CONFIG)
wcs, shape = m.strip.strip_wcs(config['grid'])
for panel in config['panels']:
    if panel['name'] == 'mm':
        panel = dict(panel, style='inferno')
    image = m.rgb(panel, m.strip.ProjectPaths(ROOT))
    fig = plt.figure(figsize=(14, 4), dpi=160)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(image, origin='lower', aspect='auto')
    ax.set_axis_off()
    outline = [pe.withStroke(linewidth=2, foreground='black')]
    ax.text(.02, .94, panel['label'], transform=ax.transAxes, color='white', fontsize=13, va='top', path_effects=outline)
    origin = np.array([shape[1]*.94, shape[0]*.72])
    for label, vector in m.strip.compass_vectors(wcs, shape).items():
        tip = origin + 100*vector
        ax.annotate('', xy=tip, xytext=origin, arrowprops=dict(arrowstyle='->', color='white', lw=1.5))
        ax.text(*(origin+130*vector), label, color='white', ha='center', va='center', path_effects=outline)
    ax.plot(*origin, marker='o', markersize=3, color='white', markeredgewidth=0, zorder=5)
    length = config['figure']['scalebar_length_deg']*3600/config['grid']['pixel_scale_arcsec']
    x, y = shape[1]*.96, shape[0]*.10
    ax.plot([x-length, x], [y,y], color='white', lw=2)
    ax.text(x-length/2, y+25, config['figure']['scalebar_label'], color='white', ha='center', path_effects=outline)
    fig.savefig(ROOT/f"website/assets/gallery-{panel['name']}.jpg", dpi=160)
    plt.close(fig)
