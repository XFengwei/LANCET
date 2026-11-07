"""Render web footprint maps from existing FITS images and pointing metadata.
Run from the repository root with the project's py312 Python environment.
Circles indicate approximate primary-beam FWHM at 225 GHz, not sensitivity contours.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from reproject import reproject_interp
from astropy.visualization import simple_norm
from astropy.visualization.wcsaxes import SphericalCircle
import astropy.units as u
root=Path(__file__).resolve().parents[2]
arrays=[('ACA 7 m','results/combination_diagnostic/overlap_tests/aca_7m_metadata.json',43,'#ffbc69'),('ALMA 12 m C2','results/combination_diagnostic/full_sampling/configuration_c2_metadata.json',26,'#66cfff'),('ALMA 12 m C5','results/combination_diagnostic/full_sampling/configuration_c5_metadata.json',26,'#66cfff')]
plt.rcParams.update({'text.color':'#eeeeee','axes.labelcolor':'#eeeeee','xtick.color':'#eeeeee','ytick.color':'#eeeeee','font.size':10})
fig=plt.figure(figsize=(13,7),facecolor='#080b10')
for i,(file,title) in enumerate([('image_HHLI-SPIRE250.fits','Herschel / 250 μm'),('image_GLIMPSE-IRAC3.fits','Spitzer / 5.8 μm')]):
 hdus=fits.open(root/'data/delivered/archive/infrared'/file); h=next(x for x in hdus if x.data is not None and x.data.ndim>=2); w=WCS(h.header).celestial; d=np.squeeze(h.data)
 target=WCS(naxis=2); target.wcs.ctype=['RA---TAN','DEC--TAN']; target.wcs.crval=[221.25,-59.8]; target.wcs.crpix=[301,301]; target.wcs.cdelt=[-2.0/3600,2.0/3600]; d,_=reproject_interp((d,w),target,shape_out=(600,600)); w=target
 ax=fig.add_subplot(1,2,i+1,projection=w);ax.set_facecolor('#080b10')
 ax.imshow(d,origin='lower',cmap='gray',norm=simple_norm(d,'asinh',percent=99.5))
 coords=[]
 for name,path,beam,color in arrays:
  fields=json.loads((root/path).read_text())['subtables']['FIELD']
  for field in fields:
   ra,dec=np.rad2deg(np.asarray(field['PHASE_DIR']).flatten());ra%=360
   coords.append(w.world_to_pixel_values(ra,dec))
   ax.add_patch(SphericalCircle((ra*u.deg,dec*u.deg),beam/7200*u.deg,edgecolor=color,facecolor='none',lw=.55,alpha=.75,transform=ax.get_transform('world')))
 xy=np.array(coords); ax.set_xlim(xy[:,0].min()-35,xy[:,0].max()+35);ax.set_ylim(xy[:,1].min()-35,xy[:,1].max()+35)
 ax.set_title(title,pad=18);ax.set_xlabel('Right ascension (J2000)');ax.set_ylabel('Declination (J2000)');ax.coords.grid(alpha=.12,color='white')
fig.legend([Line2D([],[],color='#ffbc69'),Line2D([],[],color='#66cfff')],['ACA 7 m','ALMA 12 m (C2 + C5)'],loc='lower center',ncol=2,frameon=False)
fig.subplots_adjust(left=.08,right=.98,top=.91,bottom=.15,wspace=.32)
fig.savefig(root/'website/assets/footprints.jpg',dpi=180,facecolor=fig.get_facecolor())
