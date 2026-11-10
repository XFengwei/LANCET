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
# Use the exact field, rotation, and sampling of the multiwavelength gallery.
wcs, shape = strip.strip_wcs(view)
inputs = [(name, ROOT/config['images'][key], label) for name, key, label in [('aca','aca_pbcor','ACA 7 m'), ('12m','twelve_m_pbcor','ALMA 12 m'), ('combined','combined_pbcut','ALMA 12 m + ACA 7 m')]]
inputs += [('tp', ROOT/'results/synthesize_continuum/model/continuum_230GHz.fits', 'Total-power continuum model'), ('full', ROOT/'results/multiwavelength_strip/grid/mm_i.fits', 'Fully combined 12 m + 7 m + TP')]
for name, path, label in inputs:
    with fits.open(path, memmap=True) as hdul:
        values = np.squeeze(hdul[0].data)
        if hdul[0].header.get('BUNIT', '').lower() == 'jy/beam':
            # Convert Jy/beam to MJy/sr before applying the common display scale.
            beam_sr = np.pi/(4*np.log(2))*np.deg2rad(hdul[0].header['BMAJ'])*np.deg2rad(hdul[0].header['BMIN'])
            values = values / beam_sr / 1e6
        data, footprint = reproject_interp((values, WCS(hdul[0].header).celestial), wcs, shape_out=shape)
    fig = plt.figure(figsize=(14, 4), dpi=160, facecolor='black')
    ax = fig.add_axes([0,0,1,1]); ax.set_facecolor('black')
    plotted = ax.imshow(data, origin='lower', cmap='inferno', norm=AsinhNorm(linear_width=100, vmin=20, vmax=1500, clip=True), interpolation='nearest')
    ax.set_axis_off()
    ax.text(.02,.95,label+' · 1.3 mm',transform=ax.transAxes,color='white',va='top',fontsize=14)
    bar_ax = fig.add_axes([.03,.15,.30,.025])
    bar = fig.colorbar(plotted, cax=bar_ax, orientation='horizontal', ticks=[20,100,500,1500])
    bar.set_ticks([20,100,500,1500], labels=['20','100','500','1500'])
    bar.ax.tick_params(colors='white', labelsize=8)
    bar.set_label('Surface brightness (MJy/sr)', color='white', fontsize=9)
    bar.outline.set_edgecolor('white')
    length=.02*3600/view['pixel_scale_arcsec']; x=shape[1]*.96; y=shape[0]*.09
    ax.plot([x-length,x],[y,y],color='white',lw=2)
    ax.text(x-length/2,y+25,'1 pc',ha='center',color='white')
    fig.savefig(ROOT/f'website/assets/array-{name}.jpg',dpi=160,facecolor=fig.get_facecolor());plt.close(fig)
