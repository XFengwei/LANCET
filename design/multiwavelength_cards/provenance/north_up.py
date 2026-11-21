"""Create a separate RA/Dec-aligned set of five wavelength cards."""
import argparse
import json
from pathlib import Path

import make_cards as cards
import numpy as np
from astropy.io import fits
from matplotlib.transforms import Affine2D

parser = argparse.ArgumentParser()
parser.add_argument('--plot-only', action='store_true')
parser.add_argument('--mm-bounded', action='store_true')
args = parser.parse_args()
cards.OUT = cards.ROOT / 'results/multiwavelength_cards_north_up'
if args.mm_bounded:
    cards.OUT = cards.ROOT / 'results/multiwavelength_cards_north_up_mm_fov'
cards.OUT.mkdir(parents=True, exist_ok=True)
config = cards.strip.load_toml(cards.ROOT / cards.strip.DEFAULT_CONFIG)
old_wcs, old_shape = cards.strip.strip_wcs(config['grid'])
config['grid'].update(long_axis_pa_deg=270., width_arcmin=18.,
                      height_arcmin=18 * cards.HEIGHT / cards.WIDTH)
mm_panel = next(p for p in config['panels'] if p['name'] == 'mm')
if args.mm_bounded:
    _, native_wcs = cards.strip.read_image(cards.ROOT / mm_panel['input'])
    centre = native_wcs.pixel_to_world(3547.5, 2700.)
    config['grid'].update(center_ra_deg=float(centre.ra.deg),
                          center_dec_deg=float(centre.dec.deg),
                          width_arcmin=11.78,
                          height_arcmin=11.78 * cards.HEIGHT / cards.WIDTH)
    mm_panel.pop('fill_input', None)
    # Keep the region markers inside the tighter observed field.
    mir = next(p for p in config['panels'] if p['name'] == 'mir')
    for item in mir.get('annotations', []):
        if item['label'] == 'H II region':
            item.update(x_range_arcmin=[-5.1, -2.5], y_arcmin=1.0)
        elif item['label'] == 'infrared dark cloud':
            item.update(x_range_arcmin=[2.2, 5.3], y_arcmin=-.45)
config['outputs']['root'] = str(cards.OUT.relative_to(cards.ROOT))
next(p for p in config['panels'] if p['name'] == 'radio')['label'] = 'MeerKAT  21.7 cm'
new_wcs, shape = cards.strip.strip_wcs(config['grid'])

# The tangent points coincide in the original version. For the tighter view,
# check the small affine approximation error against the exact sky transform.
old_points = np.array([[0., 0.], [1., 0.], [0., 1.]])
x, y = new_wcs.world_to_pixel(old_wcs.pixel_to_world(old_points[:, 0], old_points[:, 1]))
mapped = np.column_stack([x, y])
origin, vx, vy = mapped[0], mapped[1]-mapped[0], mapped[2]-mapped[0]
transform = Affine2D.from_values(vx[0], vx[1], vy[0], vy[1], origin[0], origin[1])
original_annotations = cards.strip.draw_annotations


def sky_annotations(ax, annotations, _shape, scale, fontsize):
    before = {id(a) for a in ax.get_children()}
    original_annotations(ax, annotations, old_shape, scale, fontsize)
    for artist in ax.get_children():
        if id(artist) not in before:
            artist.set_transform(transform + ax.transData)
    # Horizontal text needs clearance below the whole diagonal bracket.
    cx, cy = (old_shape[1]-1)/2, (old_shape[0]-1)/2
    per_arcmin = 60 / scale
    for item in annotations:
        if item['kind'] != 'bracket':
            continue
        xs = [cx + v*per_arcmin for v in item['x_range_arcmin']]
        y0 = cy + item['y_arcmin']*per_arcmin
        points = transform.transform([[xs[0], y0], [xs[1], y0]])
        text = next(t for t in ax.texts if t.get_text() == item['label'])
        text.set_transform(ax.transData)
        text.set_position((points[:, 0].mean(), points[:, 1].min() - .08*per_arcmin))
        if args.mm_bounded and item['label'] == 'infrared dark cloud':
            text.set_position((points[:, 0].mean(), points[:, 1].max() + .08*per_arcmin))
            text.set_va('bottom')


cards.strip.draw_annotations = sky_annotations
if not args.plot_only:
    cards.strip.run_reproject(cards.strip.ProjectPaths(cards.ROOT), config)
(cards.OUT/'grid_settings.json').write_text(json.dumps(config['grid'], indent=2)+'\n')
# Verify sky orientation and that the original 14 x 4 arcmin strip fits inside.
corners = np.array([[-.5, -.5], [old_shape[1]-.5, -.5],
                    [old_shape[1]-.5, old_shape[0]-.5], [-.5, old_shape[0]-.5]])
px, py = new_wcs.world_to_pixel(old_wcs.pixel_to_world(corners[:, 0], corners[:, 1]))
old_strip_fits = bool(np.all((px >= -.5) & (px <= shape[1]-.5)) and
                     np.all((py >= -.5) & (py <= shape[0]-.5)))
if args.mm_bounded:
    assert np.max(np.linalg.norm(transform.transform(corners) - np.column_stack([px, py]), axis=1)) < 1
    # Sample the whole frame perimeter to verify containment in the native map.
    xx = np.linspace(-.5, shape[1]-.5, 200)
    yy = np.linspace(-.5, shape[0]-.5, 200)
    frame_x = np.r_[xx, xx, np.full(200, -.5), np.full(200, shape[1]-.5)]
    frame_y = np.r_[np.full(200, -.5), np.full(200, shape[0]-.5), yy, yy]
    nx, ny = native_wcs.world_to_pixel(new_wcs.pixel_to_world(frame_x, frame_y))
    assert np.all((nx > 0) & (nx < 7095) & (ny > 0) & (ny < 7095))
else:
    assert old_strip_fits
vectors = cards.strip.compass_vectors(new_wcs, shape)
assert np.allclose(vectors['N'], [0, 1], atol=1e-5)
assert np.allclose(vectors['E'], [-1, 0], atol=1e-5)
assert np.isclose(new_wcs.wcs.pc[0, 1], 0)
(cards.OUT/'orientation.json').write_text(json.dumps({
    'north_up': True, 'east_left': True, 'old_strip_corners_fit': old_strip_fits,
    'bounded_by_native_mm_map': args.mm_bounded,
    'mm_model_fill_enabled': 'fill_input' in mm_panel,
    'old_strip_corners_xy': np.column_stack([px, py]).tolist(),
    'annotation_transform': transform.get_matrix().tolist(),
}, indent=2)+'\n')
cards.render(config)
