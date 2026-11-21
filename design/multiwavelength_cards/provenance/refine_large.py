"""Keep the large north-up field; tidy its unsupported right-hand extension."""
import json
import os
import argparse

import make_cards as cards
import numpy as np
from astropy.io import fits
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve

parser = argparse.ArgumentParser()
parser.add_argument('--zoom', action='store_true')
args = parser.parse_args()
source = cards.ROOT / 'results/multiwavelength_cards_north_up'
cards.OUT = cards.ROOT / 'results/multiwavelength_cards_north_up_refined'
if args.zoom:
    cards.OUT = cards.ROOT / 'results/multiwavelength_cards_north_up_refined_zoom'
out = cards.OUT
(out/'grid').mkdir(parents=True, exist_ok=True)
config = cards.strip.load_toml(cards.ROOT / cards.strip.DEFAULT_CONFIG)
old_wcs, old_shape = cards.strip.strip_wcs(config['grid'])
config['grid'] = json.loads((source/'grid_settings.json').read_text())
config['outputs']['root'] = str(out.relative_to(cards.ROOT))
next(p for p in config['panels'] if p['name'] == 'radio')['label'] = 'MeerKAT  21.7 cm'
wcs, shape = cards.strip.strip_wcs(config['grid'])
for path in (source/'grid').glob('*.fits'):
    if path.name != 'mm_i.fits' and not (out/'grid'/path.name).exists():
        os.symlink(path, out/'grid'/path.name)

# Limit the display correction to the circled right-hand background patch.
data, header = fits.getdata(source/'grid/mm_i.fits', header=True)
step = 6
ys, xs = np.mgrid[0:shape[0]:step, 0:shape[1]:step]
small = data[::step, ::step].astype(float)
protected = ndimage.binary_dilation(small > 80, iterations=10)
ellipse = ((xs/shape[1]-.86)/.17)**2 + ((ys/shape[0]-.35)/.19)**2 < 1
mask = ellipse & ~protected

# Harmonic continuation uses the surrounding values instead of nearest-pixel
# Voronoi cells. This is a local display correction, not a scientific map update.
ids = np.full(mask.shape, -1, dtype=int)
ids[mask] = np.arange(mask.sum())
rows, cols, vals = [], [], []
rhs = np.zeros(mask.sum())
iy, ix = np.where(mask)
degree = np.zeros(mask.sum())
for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
    jy, jx = iy+dy, ix+dx
    inside = (jy >= 0) & (jy < mask.shape[0]) & (jx >= 0) & (jx < mask.shape[1])
    row = ids[iy[inside], ix[inside]]
    jy, jx = jy[inside], jx[inside]
    degree[row] += 1
    unknown = mask[jy, jx]
    rows.extend(row[unknown]); cols.extend(ids[jy[unknown], jx[unknown]])
    vals.extend(-np.ones(unknown.sum()))
    rhs[row[~unknown]] += small[jy[~unknown], jx[~unknown]]
rows.extend(np.arange(mask.sum())); cols.extend(np.arange(mask.sum())); vals.extend(degree)
matrix = sparse.csr_matrix((vals, (rows, cols)), shape=(mask.sum(), mask.sum()))
solution = spsolve(matrix, rhs)
surface = small.copy()
surface[mask] = solution
result = data.copy()
affected = np.zeros(shape, dtype=bool)
for y0 in range(0, shape[0], 256):
    yy, xx = np.mgrid[y0:min(y0+256, shape[0]), :shape[1]]
    coords = [yy/step, xx/step]
    allowed = ndimage.map_coordinates(mask.astype(float), coords, order=1,
                                      mode='nearest') > .999
    smooth = ndimage.map_coordinates(surface, coords, order=3, mode='nearest')
    block = result[y0:y0+len(yy)]
    block[allowed] = smooth[allowed]
    affected[y0:y0+len(yy)] = allowed
fits.writeto(out/'grid/mm_i.fits', result, header, overwrite=True)
assert np.array_equal(result[~affected], data[~affected])
assert np.array_equal(result[data > 233], data[data > 233])
(out/'refinement.json').write_text(json.dumps({
    'field_preserved': config['grid'], 'affected_pixel_fraction': float(affected.mean()),
    'all_other_mm_pixels_identical': True, 'all_bright_mm_pixels_identical': True,
    'method': 'harmonic background continuation inside the circled lower-right ellipse',
    'bright_filament_and_its_surroundings_protected': True,
    'source_scientific_maps_unchanged': True,
    'ellipse_fractional_center_and_radii': [.86, .35, .17, .19],
    'coarse_linear_system_max_residual': float(np.max(np.abs(matrix@solution-rhs))),
}, indent=2)+'\n')
np.save(out/'affected_extrapolation_mask.npy', affected)

# Keep sky anchor positions; draw the marker shapes upright in RA/Dec axes.
cx, cy = (old_shape[1]-1)/2, (old_shape[0]-1)/2
ncx, ncy = (shape[1]-1)/2, (shape[0]-1)/2
per_arcmin = 60/config['grid']['pixel_scale_arcsec']
mir = next(p for p in config['panels'] if p['name'] == 'mir')
for item in mir['annotations']:
    if item['kind'] == 'box':
        x, y = item['center_arcmin']
        px, py = wcs.world_to_pixel(old_wcs.pixel_to_world(cx+x*per_arcmin, cy+y*per_arcmin))
        item['center_arcmin'] = [float((px-ncx)/per_arcmin), float((py-ncy)/per_arcmin)]
        item['leader_arcmin'] = [.35, .90]
    elif item['kind'] == 'bracket':
        x0, x1 = item['x_range_arcmin']
        y = cy+item['y_arcmin']*per_arcmin
        px, py = wcs.world_to_pixel(old_wcs.pixel_to_world(
            [cx+x0*per_arcmin, cx+x1*per_arcmin], [y, y]))
        item['x_range_arcmin'] = [float(v) for v in (px-ncx)/per_arcmin]
        item['y_arcmin'] = float((min(py)-ncy)/per_arcmin - .20)
if args.zoom:
    # Follow the red outline: retain the lower-left and remove top/right sky.
    nx, ny = int(round(shape[1]*.92)), int(round(shape[0]*.92))
    cropped_wcs = wcs.slice((slice(0, ny), slice(0, nx)))
    centre = cropped_wcs.pixel_to_world((nx-1)/2, (ny-1)/2)
    dx, dy = (shape[1]-nx)/2/per_arcmin, (shape[0]-ny)/2/per_arcmin
    for item in mir['annotations']:
        if item['kind'] == 'box':
            item['center_arcmin'][0] += dx
            item['center_arcmin'][1] += dy
        else:
            item['x_range_arcmin'] = [x+dx for x in item['x_range_arcmin']]
            item['y_arcmin'] += dy
    for source_path in (source/'grid').glob('*.fits'):
        target = out/'grid'/source_path.name
        original = target if source_path.name == 'mm_i.fits' else source_path
        image = fits.getdata(original)[:ny, :nx].copy()
        if target.is_symlink():
            target.unlink()
        fits.writeto(target, image, cropped_wcs.to_header(), overwrite=True)
    config['grid'].update(center_ra_deg=float(centre.ra.deg),
                          center_dec_deg=float(centre.dec.deg),
                          width_arcmin=nx*.3/60, height_arcmin=ny*.3/60)
    np.save(out/'affected_extrapolation_mask.npy', affected[:ny, :nx])
    (out/'zoom.json').write_text(json.dumps({
        'source_shape_yx': shape, 'crop_x0_y0_x1_y1': [0, 0, nx, ny],
        'output_grid': config['grid'], 'same_crop_for_all_panels': True,
        'data_reinterpolated': False,
    }, indent=2)+'\n')
    # Use the exact sliced WCS for annotation directions and the scale bar.
    cards.strip.strip_wcs = lambda grid: (cropped_wcs, (ny, nx))
    # Positions follow the author's red sketch, as fractions of the image.
    for item in mir['annotations']:
        if item['kind'] != 'bracket':
            continue
        x0, x1, y = ((.23, .455, .15) if item['label'] == 'H II region'
                     else (.62, .865, .405))
        item['x_range_arcmin'] = [
            (x-.5)*config['grid']['width_arcmin'] for x in (x0, x1)]
        item['y_arcmin'] = (y-.5)*config['grid']['height_arcmin']
(out/'annotations.json').write_text(json.dumps(mir['annotations'], indent=2)+'\n')
(out/'grid_settings.json').write_text(json.dumps(config['grid'], indent=2)+'\n')
cards.render(config)
