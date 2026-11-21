"""Three 1.3 mm continuum views on the multiwavelength-card sky frame, in MJy/sr."""
import json
from pathlib import Path
import sys
import warnings

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
import numpy as np
from reproject import reproject_adaptive

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from lancet.imaging import read_mjysr

OUT = Path(__file__).resolve().parent
MULTI = ROOT/'design/multiwavelength_cards'
FRAME = MULTI/'grid/mm_i.fits'
# Reprojected here (Jy/beam converted to MJy/sr with each image's own beam).
INPUTS = [
    ('tp', 'Total-power',
     'results/synthesize_continuum/model/continuum_224GHz_fill36.fits'),
    ('12m7m', 'ALMA Interferometry',
     'data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.clean.truesky.pb0.2.fits'),
]
# The combined card uses the grid of the ALMA multiwavelength card unchanged:
# delivered J-comb extended36 product, TP model beyond the 12m+7m grid, and
# the display cross-fade over mosaic pb 0.2-0.5.
COMBINED = ('12m7mtp', 'Total-power + Interferometry')


def prepare():
    frame_header = fits.getheader(FRAME)
    wcs = WCS(frame_header).celestial
    shape = (frame_header['NAXIS2'], frame_header['NAXIS1'])
    (OUT/'grid').mkdir(exist_ok=True)
    records = []
    for key, label, relative in INPUTS:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)
            data, header = read_mjysr(ROOT/relative)
            native_wcs = WCS(header).celestial
        image, footprint = reproject_adaptive(
            (data, native_wcs), wcs, shape_out=shape, conserve_flux=False,
            kernel='gaussian', boundary_mode='strict')
        image = image.astype(np.float32)
        image[footprint <= 0] = np.nan
        output_header = wcs.to_header()
        output_header['BUNIT'] = 'MJy/sr'
        for name in ('BMAJ', 'BMIN', 'BPA'):
            output_header[name] = float(header.get(name, 0))
        output_header.add_history('Beam-unit conversion and reprojection only.')
        fits.writeto(OUT/'grid'/f'{key}.fits', image, output_header, overwrite=True)
        records.append(dict(key=key, label=label, source=relative,
            beam_arcsec=[3600*float(header['BMAJ']), 3600*float(header['BMIN'])],
            finite_fraction=float(np.isfinite(image).mean())))
        print(f'{label}: coverage {records[-1]["finite_fraction"]:.3f}', flush=True)
    key, label = COMBINED
    combined, combined_header = fits.getdata(FRAME, header=True)
    combined_header['BUNIT'] = 'MJy/sr'
    fits.writeto(OUT/'grid'/f'{key}.fits', combined, combined_header, overwrite=True)
    records.append(dict(key=key, label=label,
        source='design/multiwavelength_cards/grid/mm_i.fits',
        finite_fraction=float(np.isfinite(combined).mean())))
    print(f'{label}: coverage {records[-1]["finite_fraction"]:.3f}', flush=True)
    manifest = dict(inputs=records, frame_reference=str(FRAME.relative_to(ROOT)),
        shape=list(shape), tp_display_factor=1.0)
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    prepare()
