"""Render the three continuum cards with the multiwavelength frame and styling."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
MULTI = ROOT/'design/multiwavelength_cards'
parser = argparse.ArgumentParser()
parser.add_argument('--refresh', action='store_true', help='Reproject the FITS inputs first')
args = parser.parse_args()
if args.refresh:
    subprocess.run([sys.executable, str(OUT/'prepare_inputs.py')], check=True)

spec = importlib.util.spec_from_file_location('multiwavelength_renderer', MULTI/'render_cards.py')
shared = importlib.util.module_from_spec(spec)
spec.loader.exec_module(shared)
# One brightness mapping for all three cards: the ALMA scale of the
# multiwavelength cards (inferno, square root, 0-1200 MJy/sr).
settings = json.loads((MULTI/'display_settings.json').read_text())
scale = settings['alma']
assert scale['stretch'] == 'sqrt'
cmap = shared.matplotlib.colormaps[scale['cmap']]
manifest = json.loads((OUT/'manifest.json').read_text())
header = fits.getheader(MULTI/'grid/mm_i.fits')
wcs = WCS(header).celestial
shape = (header['NAXIS2'], header['NAXIS1'])


def continuum_rgb(panel, project):
    data = fits.getdata(OUT/'grid'/f'{panel["name"]}.fits')
    display = np.nan_to_num(data, nan=0., posinf=0., neginf=0.)
    value = shared.strip.sqrt_scale(display, scale['vmin_mjysr'], scale['vmax_mjysr'])
    return cmap(value)[..., :3]


config = shared.strip.load_toml(ROOT/shared.strip.DEFAULT_CONFIG)
config['grid'] = json.loads((MULTI/'grid_settings.json').read_text())
config['outputs']['root'] = str(OUT.relative_to(ROOT))
config['panels'] = []
keys = ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2']
for item in manifest['inputs']:
    data, h = fits.getdata(OUT/'grid'/f'{item["key"]}.fits', header=True)
    assert data.shape == shape
    assert all(h[k] == header[k] for k in keys)
    config['panels'].append(dict(name=item['key'], label=item['label'],
        label_colour='white', compass=item['key'] == 'tp',
        scalebar=item['key'] == '12m7mtp'))
    item['finite_fraction'] = float(np.isfinite(data).mean())
    item['display_zero_replacements'] = int(np.count_nonzero(~np.isfinite(data)))
    finite = data[np.isfinite(data)]
    item['clipped_low_fraction'] = float((finite < scale['vmin_mjysr']).mean())
    item['clipped_high_fraction'] = float((finite > scale['vmax_mjysr']).mean())
manifest.update(shape=list(shape),
    frame_reference='design/multiwavelength_cards/grid/mm_i.fits',
    field=config['grid'], scale={k: scale[k] for k in ('cmap', 'stretch', 'vmin_mjysr', 'vmax_mjysr')},
    tp_display_factor=1.0, colorbars_shown=False, nan_display_value=0.,
    image_mm=[85.6, 54.], border_mm=5., border_colour='black',
    card_panels=[dict(key=p['name'], label=p['label'], compass=p['compass'],
                      scalebar=p['scalebar']) for p in config['panels']])
(OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
(OUT/'display_settings.json').write_text(json.dumps(settings, indent=2)+'\n')
(OUT/'grid_settings.json').write_text(json.dumps(config['grid'], indent=2)+'\n')
shared.OUT = OUT
shared.rgb = continuum_rgb
shared.strip.strip_wcs = lambda grid: (wcs, shape)
for pattern in ('*_card.*', '*_preview.png'):
    for path in OUT.glob(pattern):
        path.unlink()
shared.render(config)
