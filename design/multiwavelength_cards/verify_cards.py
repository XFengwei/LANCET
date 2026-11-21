"""Check exported card dimensions, black margins, coverage and registration."""
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import fitz
import numpy as np
from astropy.io import fits
from PIL import Image

out = Path(__file__).resolve().parent
geometry = json.loads((out / 'geometry.json').read_text())
report = {'cards': {}, 'grids': {}, 'all_passed': True}
reference_header = None
for path in sorted((out / 'grid').glob('*.fits')):
    data, header = fits.getdata(path, header=True)
    keys = ('CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2',
            'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2')
    signature = {k: header.get(k) for k in keys}
    if reference_header is None:
        reference_header = signature
    assert signature == reference_header, path
    assert list(data.shape) == geometry['grid_shape'], path
    report['grids'][path.name] = {'finite_fraction': float(np.isfinite(data).mean())}
for path in sorted(out.glob('*_card.png')):
    with Image.open(path) as im:
        a = np.asarray(im.convert('RGB'))
        dpi = im.info['dpi']
        assert abs(dpi[0] - 1200) < .1
        assert abs(im.width - 95.6/25.4*1200) <= 1
        assert abs(im.height - 64/25.4*1200) <= 1
        margin = int(5/25.4*1200) - 1
        assert np.all(a[:margin] == 0) and np.all(a[-margin:] == 0)
        assert np.all(a[:, :margin] == 0) and np.all(a[:, -margin:] == 0)
        size = list(im.size)
    with Image.open(path.with_suffix('.tiff')) as im:
        assert list(im.size) == size
    with fitz.open(path.with_suffix('.pdf')) as doc:
        mm = [doc[0].rect.width * 25.4/72, doc[0].rect.height * 25.4/72]
        assert np.allclose(mm, [95.6, 64], atol=.001)
    report['cards'][path.stem] = {'pixels': size, 'dpi': list(dpi),
                                'pdf_mm': mm, 'black_border_passed': True, 'border_mm': 5.0}
assert len(report['cards']) == 5
bounds = [p['image_bounds_mm'] for p in geometry['panels'].values()]
assert np.allclose(bounds, [5, 5, 85.6, 54], atol=1e-6)
(out / 'verification.json').write_text(json.dumps(report, indent=2)+'\n')
with ZipFile(out / 'five_wavelength_cards.zip', 'w', ZIP_DEFLATED) as bundle:
    for pattern in ('*_card.png', '*_card.tiff', '*_card.pdf', 'README.md'):
        for path in sorted(out.glob(pattern)):
            bundle.write(path, path.name)
print(json.dumps(report, indent=2))
