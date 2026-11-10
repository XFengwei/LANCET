"""Export full-field array comparisons directly from the continuum FITS images."""
from pathlib import Path
import sys
import tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'scripts'), str(ROOT/'src')]
import multiwavelength_strip as strip
config = tomllib.loads((ROOT/'tasks/combination_diagnostic/task.toml').read_text())
view = tomllib.loads((ROOT/'tasks/multiwavelength_strip/task.toml').read_text())['grid']
# Include the full filament and margin around the mosaic, on a shared rotated grid.
view.update(width_arcmin=18., height_arcmin=8., pixel_scale_arcsec=.5)
wcs, shape = strip.strip_wcs(view)
for name, key, label in [('aca','aca_pbcor','ACA 7 m'), ('12m','twelve_m_pbcor','ALMA 12 m'), ('combined','combined_pbcut','ALMA 12 m + ACA 7 m')]:
    with fits.open(ROOT/config['images'][key], memmap=True) as hdul:
        data, footprint = reproject_interp((np.squeeze(hdul[0].data), WCS(hdul[0].header).celestial), wcs, shape_out=shape)
    fig = plt.figure(figsize=(14, 14*8/18), dpi=160, facecolor='#080b10')
    ax = fig.add_axes([0,0,1,1]); ax.set_facecolor('#080b10')
    ax.imshow(data, origin='lower', cmap='inferno', norm=AsinhNorm(linear_width=20, vmin=-20, vmax=500), interpolation='nearest')
    ax.set_axis_off()
    ax.text(.02,.95,label+' · 1.3 mm',transform=ax.transAxes,color='white',va='top',fontsize=14)
    ax.text(.02,.025,'Shared colour scale · −20 to 500 MJy/sr',transform=ax.transAxes,color='white',fontsize=9)
    length=.02*3600/view['pixel_scale_arcsec']; x=shape[1]*.96; y=shape[0]*.09
    ax.plot([x-length,x],[y,y],color='white',lw=2)
    ax.text(x-length/2,y+25,'1 pc',ha='center',color='white')
    fig.savefig(ROOT/f'website/assets/array-{name}.jpg',dpi=160,facecolor=fig.get_facecolor());plt.close(fig)
