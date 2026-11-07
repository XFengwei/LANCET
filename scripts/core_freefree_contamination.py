"""Estimate smooth CORNISH free-free emission within saved getsf footprints."""
from pathlib import Path
import json
import tomllib
import argparse

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from regions import Regions
from scipy.ndimage import map_coordinates, uniform_filter
from astropy.stats import SigmaClip
from photutils.background import Background2D, StdBackgroundRMS

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--radio-rms-mjybeam', type=float)
parser.add_argument('--output-root', type=Path, default=Path('results/construct_freefree_model/core_contamination'))
args = parser.parse_args()
cfg = tomllib.loads((root / 'tasks/construct_freefree_model/task.toml').read_text())
radio_path = root / cfg['inputs']['cornish_map']
radio, rh = fits.getdata(radio_path, header=True)
assert radio.ndim == 2 and rh['BUNIT'].upper() == 'JY/BEAM'
out = root / args.output_root
out.mkdir(parents=True, exist_ok=True)
settings = cfg['core_contamination']
if args.radio_rms_mjybeam is not None:
    if args.radio_rms_mjybeam <= 0:
        raise ValueError('Radio rms must be positive')
    settings['adopted_radio_rms_mjy_per_beam'] = args.radio_rms_mjybeam
noise_box = int(round(settings['noise_box_beams']*rh['BMAJ']/abs(rh['CDELT1'])))
noise_box += 1 - noise_box % 2
background = Background2D(radio, box_size=noise_box, filter_size=3,
                         sigma_clip=SigmaClip(sigma=3, maxiters=10),
                         bkgrms_estimator=StdBackgroundRMS())
thermal_rms = background.background_rms
# Negative excursions retain broad bowls/sidelobes that symmetric clipping
# could reject as outliers. For symmetric noise their conditional RMS is sigma.
negative = radio < 0
negative_power = uniform_filter(np.where(negative, radio**2, 0.), size=noise_box, mode='reflect')
negative_fraction = uniform_filter(negative.astype(float), size=noise_box, mode='reflect')
negative_rms = np.sqrt(np.maximum(negative_power, 0)/np.maximum(negative_fraction, 1e-12))
noise = np.maximum(thermal_rms, negative_rms)
assert np.all(np.isfinite(noise)) and np.all(noise > 0)
threshold = settings['adopted_radio_rms_mjy_per_beam']*settings['detection_sigma']/1000
limit_brightness = settings['adopted_radio_rms_mjy_per_beam']*settings['upper_limit_sigma']/1000
robust = np.isfinite(radio) & (radio > threshold)
retained_radio = np.where(robust, radio, 0.)
for name, values in [('radio_noise_jybeam',noise),('radio_clipped_rms_jybeam',thermal_rms),
                     ('radio_negative_rms_jybeam',negative_rms),('radio_retained_uniform3sigma_jybeam',retained_radio)]:
    fits.writeto(out/(name+'.fits'),values.astype('float32'),rh,overwrite=True)
mask_header=rh.copy();mask_header['BUNIT']='1'
fits.writeto(out/'radio_positive_uniform3sigma_mask.fits',robust.astype('uint8'),mask_header,overwrite=True)
rw = WCS(rh).celestial
beam_area = np.pi * rh['BMAJ'] * rh['BMIN'] / (4*np.log(2))  # deg2
catalogue = Table.read(root / 'results/source_catalogue/tables/source_catalogue.ecsv')
image_path = root / settings['continuum_image']
alma, ah = fits.getdata(image_path, header=True)
alma = np.squeeze(alma)
assert ah['BUNIT'].lower() == 'jy/beam'
aw = WCS(ah).celestial
native_alma_beam_area = np.pi*ah['BMAJ']*ah['BMIN']/(4*np.log(2))
pixel_area = abs(proj_plane_pixel_area(aw))
aperture_path = root / 'extraction/getsf/triple/results/+sources/+catalogs/G316.sw.sources.ok.cat.1302.foo.fk5.reg'
apertures = {}
for line in aperture_path.read_text().splitlines():
    line = line.lstrip('# ')
    if line.startswith('ellipse('):
        shape, number = line.split('#')
        apertures[int(number.strip())] = Regions.parse('fk5\n'+shape, format='ds9')[0]
frequency = tomllib.loads((root / 'tasks/physical_properties/task.toml').read_text())['physics']['frequency_ghz']
spectral_factor = (frequency/5.5)**(-.1)
alma_beam_area = np.pi*(.65/3600)**2/(4*np.log(2))
rows = []
for source in catalogue:
    aperture = apertures[int(source['NO'])].to_pixel(aw)
    weights = aperture.to_mask(mode='subpixels', subpixels=8)
    assert abs(weights.data.sum()/aperture.area - 1) < .01
    box = weights.bbox
    y, x = np.mgrid[box.iymin:box.iymax, box.ixmin:box.ixmax]
    sky = aw.pixel_to_world(x, y)
    rx, ry = rw.world_to_pixel(sky)
    values = map_coordinates(radio, [ry, rx], order=1, mode='constant', cval=np.nan, prefilter=False)
    retained_values = map_coordinates(retained_radio, [ry, rx], order=1, mode='constant', cval=np.nan, prefilter=False)
    retained_mask = map_coordinates(robust.astype(float), [ry, rx], order=0, mode='constant', cval=0, prefilter=False) > .5
    retained_values = np.where(retained_mask, retained_values, 0.)
    keep = weights.data > 0
    finite = np.isfinite(values) & keep
    coverage = float(weights.data[finite].sum()/weights.data[keep].sum())
    total = float(np.sum(retained_values[finite]*weights.data[finite])*pixel_area/beam_area*spectral_factor) if coverage > .999999 else np.nan
    retained_area_fraction = float(weights.data[keep & retained_mask].sum()/weights.data[keep].sum())
    upper_limit = limit_brightness*weights.data.sum()*pixel_area/beam_area*spectral_factor if coverage > .999999 else np.nan
    detected = bool(np.isfinite(total) and total > 0)
    native_values = alma[box.iymin:box.iymax,box.ixmin:box.ixmax]
    native_finite = np.isfinite(native_values) & keep
    alma_coverage = float(weights.data[native_finite].sum()/weights.data[keep].sum())
    observed = float(np.sum(native_values[native_finite]*weights.data[native_finite])*pixel_area/native_alma_beam_area) if alma_coverage > .999999 else np.nan
    native_peak = float(map_coordinates(alma,[[aperture.center.y],[aperture.center.x]],order=1,mode='constant',cval=np.nan,prefilter=False)[0])
    cx, cy = rw.world_to_pixel(aperture.to_sky(aw).center)
    radio_center = float(map_coordinates(retained_radio, [[cy], [cx]], order=1, mode='constant', cval=np.nan, prefilter=False)[0])
    if map_coordinates(robust.astype(float),[[cy],[cx]],order=0,mode='constant',cval=0)[0] < .5:
        radio_center = 0.
    predicted_peak = radio_center*spectral_factor*alma_beam_area/beam_area
    rows.append({'ID':int(source['ID']), 'NO':int(source['NO']),
                 'Subregion':str(source['Subregion']), 'coverage':coverage,
                 'alma_coverage':alma_coverage,
                 'retained_area_fraction':retained_area_fraction,
                 'radio_detected':detected,
                 'freefree_upper_limit_jy':upper_limit if not detected else np.nan,
                 'observed_12m7m_flux_jy':observed,
                 'catalogue_background_subtracted_flux_jy':float(source['FXT_BST01']),
                 'smooth_freefree_flux_jy':total,
                 'smooth_freefree_fraction':total/observed,
                 'freefree_fraction_upper_limit':upper_limit/observed if not detected else np.nan,
                 'smooth_freefree_peak_jy_per_065beam':predicted_peak,
                 'smooth_freefree_peak_fraction':predicted_peak*native_alma_beam_area/alma_beam_area/native_peak,
                 'aperture_area_arcsec2':float(weights.data.sum()*pixel_area*3600**2)})
table = Table(rows=rows)
table.write(out/'smooth_core_contamination.ecsv', format='ascii.ecsv', overwrite=True)
table.write(out/'smooth_core_contamination.csv', format='ascii.csv', overwrite=True)
summary = {'radio_map':str(radio_path.relative_to(root)), 'footprints':str(aperture_path.relative_to(root)),
           'continuum_image':str(image_path.relative_to(root)),
           'adopted_radio_rms_mjy_per_beam':settings['adopted_radio_rms_mjy_per_beam'],
           'detection_sigma':settings['detection_sigma'],
           'upper_limit_sigma':settings['upper_limit_sigma'],
           'minimum_radio_mjy_per_beam':threshold*1000,
           'native_alma_beam_arcsec':[ah['BMAJ']*3600,ah['BMIN']*3600],
           'radio_beam_arcsec':[rh['BMAJ']*3600, rh['BMIN']*3600],
           'frequency_ghz':frequency, 'spectral_index':-.1, 'spectral_factor':spectral_factor,
           'model':'Smooth brightness within the native CORNISH beam; optically thin extrapolation. Bilinear interpolation, no additional convolution or deconvolution.',
           'noise_box_pixels':noise_box, 'noise_box_arcsec':noise_box*abs(rh['CDELT1'])*3600,
           'noise_method':'Maximum of sigma-clipped mesh RMS and local RMS of all negative excursions, relative to zero. No radio brightness background is subtracted.',
           'noise_range_mjybeam':[float(noise.min()*1000),float(noise.max()*1000)],
           'retained_native_pixels':int(robust.sum()), 'native_pixels':int(radio.size),
           'measurement':'Native pixels >3 times the adopted uniform rms retained. Original-image 12m7m aperture flux denominator includes local background. For covered nondetections, nominal 1-sigma brightness limit is integrated over the full footprint under the smooth-emission assumption; not a sum of independent pixel errors. Spatial noise maps remain diagnostic.',
           'regions':{}}
for code in ['Y','I','E']:
    select = (table['Subregion']==code) & np.isfinite(table['smooth_freefree_fraction'])
    t=table[select]
    summary['regions'][code]={'catalogue_count':int(np.sum(table['Subregion']==code)),
        'covered_count':len(t), 'median_fraction':float(np.median(t['smooth_freefree_fraction'])) if len(t) else None,
        'with_retained_radio':int(sum(t['smooth_freefree_flux_jy']>0)),
        'max_fraction':float(max(t['smooth_freefree_fraction'])) if len(t) else None,
        'flux_weighted_fraction':float(sum(t['smooth_freefree_flux_jy'])/sum(t['observed_12m7m_flux_jy'])) if len(t) else None,
        'count_above_10percent':int(sum(t['smooth_freefree_fraction']>.1)),
        'count_above_100percent':int(sum(t['smooth_freefree_fraction']>1))}
    limits = np.asarray(t['freefree_fraction_upper_limit'])
    limits = limits[np.isfinite(limits)]
    summary['regions'][code]['nondetections'] = {
        'count':len(limits),
        'median_fraction_limit':float(np.median(limits)) if len(limits) else None,
        'limits_above_10percent':int(sum(limits>.1)),
        'limits_above_50percent':int(sum(limits>.5)),
        'limits_above_100percent':int(sum(limits>1))}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
