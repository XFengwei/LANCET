"""Verify common brightness presentation and package the three continuum views."""
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np
from astropy.io import fits
from PIL import Image
import pymupdf

out = Path(__file__).resolve().parent
manifest = json.loads((out/'manifest.json').read_text())
geometry = json.loads((out/'geometry.json').read_text())['panels']
records = {}
reference = None
multi = out.parent/'multiwavelength_cards'
frame_header = fits.getheader(multi/'grid/mm_i.fits')
wcs_keys = ('CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2',
            'PC1_1','PC1_2','PC2_1','PC2_2')
expected_wcs = tuple(frame_header.get(k) for k in wcs_keys)
multi_scale = json.loads((multi/'display_settings.json').read_text())['alma']
assert manifest['scale']['cmap'] == multi_scale['cmap'] == 'inferno'
assert manifest['scale']['stretch'] == multi_scale['stretch'] == 'sqrt'
assert manifest['scale']['vmin_mjysr'] == multi_scale['vmin_mjysr']
assert manifest['scale']['vmax_mjysr'] == multi_scale['vmax_mjysr']
assert manifest['tp_display_factor'] == 1.0
for number, item in enumerate(manifest['inputs'], 1):
    key = item['key']
    data, h = fits.getdata(out/'grid'/f'{key}.fits', header=True)
    assert h['BUNIT'] == 'MJy/sr'
    assert list(data.shape) == manifest['shape']
    signature = tuple(h.get(k) for k in ('CRPIX1','CRPIX2','CRVAL1','CRVAL2',
        'CDELT1','CDELT2','PC1_1','PC1_2','PC2_1','PC2_2'))
    if reference is None:
        reference = signature
    assert signature == reference
    assert signature == expected_wcs
    assert data.shape == (frame_header['NAXIS2'],frame_header['NAXIS1'])
    prefix = out/f'{number:02d}_{key}_card'
    with Image.open(prefix.with_suffix('.png')) as im:
        a = np.asarray(im.convert('RGB'))
        dpi = im.info['dpi']
        assert np.allclose(dpi, [1200,1200], atol=.01)
        assert abs(im.width - 95.6/25.4*1200) <= 1
        assert abs(im.height - 64/25.4*1200) <= 1
        margin = int(5/25.4*1200)-1
        assert np.all(a[:margin] == 0) and np.all(a[-margin:] == 0)
        assert np.all(a[:, :margin] == 0) and np.all(a[:, -margin:] == 0)
        pixels = list(im.size)
    with Image.open(prefix.with_suffix('.tiff')) as im:
        assert list(im.size) == pixels
    with pymupdf.open(prefix.with_suffix('.pdf')) as doc:
        mm = [doc[0].rect.width*25.4/72, doc[0].rect.height*25.4/72]
        assert np.allclose(mm, [95.6,64], atol=.001)
    records[key] = dict(pixels=pixels, pdf_mm=mm, finite_fraction=float(np.isfinite(data).mean()),
                        colorbars_shown=False)
images = [r['image_bounds_mm'] for r in geometry.values()]
assert np.allclose(images, images[0], atol=1e-6)
assert np.allclose(images, [5.,5.,85.6,54.], atol=1e-6)
assert manifest['colorbars_shown'] is False
assert list(records) == ['tp', '12m7m', '12m7mtp']
assert len(list(out.glob('*_card.png'))) == 3
(out/'verification.json').write_text(json.dumps(dict(all_passed=True, cards=records), indent=2)+'\n')
with ZipFile(out/'continuum_comparison_cards.zip', 'w', ZIP_DEFLATED) as bundle:
    for pattern in ('*_card.png','*_card.tiff','*_card.pdf','README.md','manifest.json'):
        for path in sorted(out.glob(pattern)):
            bundle.write(path,path.name)
print('All three cards passed: 85.6 x 54 mm images, 5 mm black borders, no colourbars.')
