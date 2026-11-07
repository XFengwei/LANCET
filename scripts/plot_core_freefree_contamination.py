"""Display spatial association and smooth-template aperture flux comparison."""
from pathlib import Path
import sys
import tomllib
import argparse
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm, LogNorm
from matplotlib.lines import Line2D
import seaborn as sns
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy import units as u
from regions import Regions

sys.path.insert(0, '/Users/fexu/.codex/skills/nature-figure/scripts')
from audit_panel_alignment import require_matplotlib_panel_alignment

root=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-root',type=Path,default=Path('results/construct_freefree_model/core_contamination'))
args=parser.parse_args()
config=tomllib.loads((root/'tasks/construct_freefree_model/task.toml').read_text())
out=root/args.output_root
metadata=json.loads((out/'summary.json').read_text())
table=Table.read(out/'smooth_core_contamination.ecsv')
source=Table.read(root/'results/source_catalogue/tables/source_catalogue.ecsv')
radio,rh=fits.getdata(root/'data/delivered/archive/radio/CORNISH_5p5GHz.fits',header=True)
noise=fits.getdata(out/'radio_noise_jybeam.fits')
snr=radio/noise
data,h=fits.getdata(root/config['core_contamination']['continuum_image'],header=True)
gain=fits.getdata(root/'data/delivered/alma/7MTM2TM1/cont/casa/primary_beam_rp2miriad.fits')
data=np.where(np.isfinite(gain) & (gain>=.2),np.squeeze(data),np.nan)
assert h['BUNIT'].lower()=='jy/beam'
center=SkyCoord(221.3468*u.deg,-59.8230*u.deg)
cut=Cutout2D(np.squeeze(data)*1000,center,(150*u.arcsec,150*u.arcsec),wcs=WCS(h).celestial,mode='partial',fill_value=np.nan)
selected=table[(table['Subregion']=='E') & np.isfinite(table['smooth_freefree_fraction'])]
assert len(selected)==94
covered_count=len(selected)
nondetections=selected[~selected['radio_detected']]
selected=selected[selected['radio_detected']]

plt.rcParams.update({'font.family':'sans-serif','font.size':9,'axes.labelsize':9,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(8.2,4.3))
height=.34*8.2/4.3
left=fig.add_axes([.085,.18,.34,height],projection=cut.wcs)
right=fig.add_axes([.635,.18,.34,height])
im=left.imshow(cut.data,origin='lower',cmap=sns.color_palette('mako',as_cmap=True),norm=AsinhNorm(linear_width=.3,vmin=-.2,vmax=12),rasterized=True)
left.set_facecolor('.12')
sigma=metadata['adopted_radio_rms_mjy_per_beam']
levels=np.array([3,6,12,24,48])*sigma
contours=left.contour(radio*1000,levels=levels,transform=left.get_transform(WCS(rh).celestial),colors='#f5a442',linewidths=.65)
contours.set_rasterized(True)
footpath=root/'extraction/getsf/triple/results/+sources/+catalogs/G316.sw.sources.ok.cat.1302.foo.fk5.reg'
footprints={}
for line in footpath.read_text().splitlines():
    line=line.lstrip('# ')
    if line.startswith('ellipse('):
        shape,num=line.split('#')
        footprints[int(num.strip())]=Regions.parse('fk5\n'+shape,format='ds9')[0]
for row in source[source['Subregion']=='E']:
    artist=footprints[int(row['NO'])].to_pixel(cut.wcs).as_artist(facecolor='none',edgecolor='white',lw=.5)
    artist.set_rasterized(True)
    left.add_patch(artist)
left.set_xlim(-.5,cut.data.shape[1]-.5);left.set_ylim(-.5,cut.data.shape[0]-.5)
left.coords[0].set_ticks(number=3);left.coords[1].set_ticks(number=3)
left.set_xlabel('RA (ICRS)');left.set_ylabel('Dec (ICRS)')
left.text(.025,.975,'(a)',color='white',transform=left.transAxes,va='top',fontsize=10)
left.legend(handles=[Line2D([],[],color='#f5a442',lw=1,label='CORNISH: 3, 6, 12, 24, 48 sigma')],loc='lower left',bbox_to_anchor=(.02,1.015),frameon=False,handlelength=1.3)
cax=left.inset_axes([1.035,0,.035,1]);cb=fig.colorbar(im,cax=cax)
cb.set_label(r'1.3 mm intensity (mJy beam$^{-1}$)',fontsize=8)
cb.ax.tick_params(labelsize=8)
cb.set_ticks([0,.5,2,6,12]);cb.set_ticklabels(['0','0.5','2','6','12'])

x=np.asarray(selected['observed_12m7m_flux_jy'])*1000
y=np.asarray(selected['smooth_freefree_flux_jy'])*1000
assert np.all(x>0)
right.scatter(x,y,s=18,c='#007b96',alpha=.8,edgecolors='none',zorder=3,label='SNR > 3')
limit_x=np.asarray(nondetections['observed_12m7m_flux_jy'])*1000
limit_y=np.asarray(nondetections['freefree_upper_limit_jy'])*1000
assert np.all(limit_x>0) and np.all(limit_y>0)
right.errorbar(limit_x,limit_y,yerr=.35*limit_y,uplims=True,fmt='o',ms=4,
 mfc='none',mec='#007b96',ecolor='#007b96',elinewidth=.6,capsize=1.5,
 alpha=.65,label='1-sigma upper limit',zorder=2)
line=np.geomspace(.25,max(160,x.max()*1.2),100)
right.plot(line,line,'--',c='.3',lw=.9,label='1:1')
right.plot(line,.1*line,':',c='.5',lw=.9,label='10%')
assert np.all(y>0)
right.set_xscale('log');right.set_yscale('log')
right.set_xlim(.25,max(160,x.max()*1.2));right.set_ylim(min(.001,y.min()*.7),80)
right.set_xlabel('1.3 mm flux within core footprint (mJy)')
right.set_ylabel('Expected free–free flux at 223.8 GHz (mJy)')
right.text(.025,.975,'(b)',transform=right.transAxes,va='top',fontsize=10)
handles,labels=right.get_legend_handles_labels()
legend_handles=dict(zip(labels,handles))
legend_order=['SNR > 3','1-sigma upper limit','1:1','10%']
right.legend([legend_handles[label] for label in legend_order],legend_order,
 loc='lower left',frameon=False)
for row in selected[selected['smooth_freefree_fraction']>1]:
    xx=float(row['observed_12m7m_flux_jy'])*1000;yy=float(row['smooth_freefree_flux_jy'])*1000
    if int(row['ID'])==213:
        right.annotate('ID 213',(xx,yy),xytext=(5,3),textcoords='offset points',fontsize=8)
require_matplotlib_panel_alignment(fig,axes=[left,right],panel_ids=['a','b'],
 row_groups=[['a','b']],json_out=out/'cornish_core_comparison.alignment.json',strict=True)
fig.savefig(out/'cornish_core_comparison.pdf',dpi=300,bbox_inches='tight',pad_inches=.04)
fig.savefig(out/'cornish_core_comparison.svg',dpi=300,bbox_inches='tight',pad_inches=.04)
fig.savefig(out/'cornish_core_comparison.png',dpi=300,bbox_inches='tight',pad_inches=.04)
fig.savefig(out/'cornish_core_comparison.tiff',dpi=300,bbox_inches='tight',pad_inches=.04)
plt.close(fig)
(out/'figure_caption.txt').write_text(f'Left: original 12m+7m 1.3 mm continuum, displayed within PB >= 0.2. CORNISH contours are 3, 6, 12, 24 and 48 times the adopted rms of {sigma:g} mJy per native 2.5-arcsec beam. White ellipses mark core footprints. Right: original-image aperture flux, without background subtraction, versus optically thin 223.8 GHz free-free flux. Filled circles use native radio pixels above 3 sigma. Hollow circles with downward arrows show nominal 1-sigma limits for covered nondetections, converting the {sigma:g} mJy/beam brightness limit into integrated flux over each full footprint under the smooth-emission assumption. Dashed and dotted lines show equality and 10 percent.\n')

fig,ax=plt.subplots(figsize=(4.5,4))
im=ax.imshow(noise*1000,origin='lower',cmap='viridis',norm=LogNorm(vmin=.1,vmax=12),rasterized=True)
ax.contour(radio*1000,levels=[3*sigma],colors='white',linewidths=.5)
ax.set_xlabel('Native radio pixel x');ax.set_ylabel('Native radio pixel y')
fig.colorbar(im,ax=ax,label='Local noise estimate (mJy/beam)')
fig.savefig(out/'radio_noise_map.pdf',dpi=300,bbox_inches='tight')
fig.savefig(out/'radio_noise_map.png',dpi=300,bbox_inches='tight')
plt.close(fig)
