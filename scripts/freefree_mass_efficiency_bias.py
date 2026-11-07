"""Propagate the smooth-radio model into fragment/core mass-ratio scenarios."""
from pathlib import Path
import json
import argparse

import numpy as np
from astropy.io import fits
from astropy.table import Table, join
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel, convolve_fft
from reproject import reproject_interp

root=Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--radio-output-root',type=Path,default=Path('results/construct_freefree_model/core_contamination'))
args=parser.parse_args()
out=root/args.radio_output_root
metadata=json.loads((out/'summary.json').read_text())
radio,rh=fits.getdata(root/metadata['radio_map'],header=True)
threshold=metadata['minimum_radio_mjy_per_beam']/1000
template=np.where(radio>threshold,radio,0.)
factor=metadata['spectral_factor']
aca,ah=fits.getdata(root/'data/delivered/alma/7M/cont/IRDC316.8.ACA.p1.cont.image.fits',header=True)
pb=fits.getdata(root/'data/delivered/alma/7M/cont/IRDC316.8.ACA.p1.cont.pb.fits')
labels,lh=fits.getdata(root/'extraction/astrodendro/my_mask_idx.fits',header=True)
assert aca.shape==pb.shape==labels.shape
cat=Table.read(root/'extraction/astrodendro/my_dendrogram.dat',format='ascii')
selection=Table.read(root/'results/core_formation/tables/fragment_pb_selection.ecsv')
selection.remove_column('flux')
fragments=join(selection,cat[['ID','_idx','flux','area_exact']],keys='ID')
sigma_scale=1/(2*np.sqrt(2*np.log(2)))
pa=np.deg2rad(ah['BPA']);v=np.array([np.sin(pa),np.cos(pa)])
target_cov=((ah['BMAJ']*3600)**2*np.outer(v,v)+(ah['BMIN']*3600)**2*(np.eye(2)-np.outer(v,v)))*sigma_scale**2
radio_cov=np.eye(2)*(rh['BMAJ']*3600*sigma_scale)**2
jac=WCS(rh).celestial.pixel_scale_matrix*3600
cov=np.linalg.inv(jac)@(target_cov-radio_cov)@np.linalg.inv(jac).T
ev,evec=np.linalg.eigh(cov);assert np.all(ev>0)
kernel=Gaussian2DKernel(np.sqrt(ev[1]),np.sqrt(ev[0]),theta=np.arctan2(evec[1,1],evec[0,1]))
kernel.normalize()
radio_area=np.pi*(rh['BMAJ']*3600)*(rh['BMIN']*3600)/(4*np.log(2))
smooth=convolve_fft(template/radio_area*factor,kernel,boundary='fill',fill_value=0)
support=convolve_fft(np.ones(radio.shape),kernel,boundary='fill',fill_value=0)
model,foot=reproject_interp((smooth,WCS(rh).celestial),WCS(ah).celestial,shape_out=aca.shape)
coverage,_=reproject_interp((support,WCS(rh).celestial),WCS(ah).celestial,shape_out=aca.shape)
pixel_area=abs(ah['CDELT1']*ah['CDELT2'])*3600**2
aca_area=np.pi*(ah['BMAJ']*3600)*(ah['BMIN']*3600)/(4*np.log(2))
rows=[]
for row in fragments:
    mask=labels==row['_idx'];assert mask.sum()==round(row['area_exact'])
    native_flux=float(np.sum(aca[mask])*pixel_area/aca_area)
    assert np.isclose(native_flux,row['flux'],rtol=2e-5)
    complete=bool(np.all(np.isfinite(model[mask])) and np.all(coverage[mask]>.99999) and np.all(foot[mask]>0))
    ff=float(np.sum(model[mask]*pb[mask])*pixel_area) if complete else np.nan
    rows.append({'ID':row['ID'],'Subregion':row['Subregion'],'selected':row['selected'],
      'radio_covered':complete,'aca_flux_jy':native_flux,'freefree_flux_jy':ff,
      'fragment_freefree_fraction':ff/native_flux,'fragment_mass_msun':float(row['mass_dense'])*1.8})
table=Table(rows=rows)
table.write(out/'fragment_freefree_contamination.ecsv',format='ascii.ecsv',overwrite=True)
table.write(out/'fragment_freefree_contamination.csv',format='ascii.csv',overwrite=True)
cores=join(Table.read(root/'results/physical_properties/tables/core_physical_properties.ecsv'),Table.read(out/'smooth_core_contamination.ecsv'),keys='ID',table_names=['mass','ff'])
eff=Table.read(root/'results/core_formation/tables/core_formation_summary.ecsv')
result={'radio_rms_mjybeam':metadata['adopted_radio_rms_mjy_per_beam'],
 'fragment_method':'Convolve thresholded 2.5 arcsec radio surface brightness to the ACA beam, reproject, apply ACA PB and integrate true leaf masks. Catalogue fluxes independently reproduced.',
 'core_method':'Sensitivity scenario: multiply each adopted core mass by 1 minus its original-image aperture free-free fraction, holding temperature, geometry and optical-depth factor fixed. This is not a calibrated correction of background-subtracted core photometry.',
 'regions':{}}
for r in eff:
    f=table[(table['selected'])&(table['Subregion']==r['label'])]
    c=cores[cores['Subregion_mass']==r['subregion_code']]
    cf=np.asarray(c['smooth_freefree_fraction']);cf=np.where(np.isfinite(cf),cf,0)
    ff=np.asarray(f['fragment_freefree_fraction']);ff=np.where(np.isfinite(ff),ff,0)
    assert np.all((cf>=0)&(cf<1)) and np.all((ff>=0)&(ff<1))
    mc=float(sum(c['Mass_core']));mf=float(sum(f['fragment_mass_msun']))
    cm=float(sum(c['Mass_core']*(1-cf)));fm=float(sum(f['fragment_mass_msun']*(1-ff)))
    result['regions'][str(r['region'])]={'core_mass':mc,'fragment_mass':mf,
     'radio_covered_fragments':int(sum(f['radio_covered'])),'selected_fragments':len(f),
     'core_mass_scenario':cm,'fragment_mass_scenario':fm,
     'fractional_core_mass_reduction':1-cm/mc,'fractional_fragment_mass_reduction':1-fm/mf,
     'core_cloud_before':mc/float(r['cloud_mass']),'core_cloud_scenario':cm/float(r['cloud_mass']),
     'fragment_cloud_before':mf/float(r['cloud_mass']),'fragment_cloud_scenario':fm/float(r['cloud_mass']),
     'core_fragment_before':mc/mf,'core_fragment_scenario':cm/fm}
(out/'mass_efficiency_bias.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
print(table[table['Subregion']=='Evolved'])
