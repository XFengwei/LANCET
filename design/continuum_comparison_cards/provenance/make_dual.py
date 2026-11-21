"""Crop to the 12m+7m frame and render the adopted dual brightness scale."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts')]
from multiwavelength_strip import dual_colour_rgb

OUT = Path(__file__).resolve().parent
SOURCE = ROOT/'results/continuum_comparison_cards'
(OUT/'grid').mkdir(exist_ok=True)
manifest = json.loads((SOURCE/'manifest.json').read_text())
tp_factor = 2.7
next(item for item in manifest['inputs'] if item['key'] == 'tp')['label'] = 'Total-power'
manifest['tp_display_factor'] = tp_factor
manifest['tp_display_factor_scope'] = 'TP image display only; native FITS values unchanged'
benchmark, h = fits.getdata(SOURCE/'grid/12m7m.fits', header=True)
y, x = np.where(np.isfinite(benchmark))
x0, y0, x1, y1 = int(x.min()), int(y.min()), int(x.max()+1), int(y.max()+1)
# Match the physical card aspect ratio without stretching the sky.
width = int(round((y1-y0)*85.6/54.))
centre_x = (x0+x1)/2
x0 = int(round(centre_x-width/2)); x1 = x0+width
assert 0 <= x0 < x1 <= benchmark.shape[1]
wcs = WCS(h).celestial.slice((slice(y0,y1), slice(x0,x1)))
shape = (y1-y0, x1-x0)
for item in manifest['inputs']:
    data, header = fits.getdata(SOURCE/'grid'/f'{item["key"]}.fits', header=True)
    cropped = data[y0:y1, x0:x1]
    out_header = wcs.to_header()
    for key in ('BUNIT','BMAJ','BMIN','BPA'):
        out_header[key] = header[key]
    fits.writeto(OUT/'grid'/f'{item["key"]}.fits',cropped,out_header,overwrite=True)
    item['finite_fraction'] = float(np.isfinite(cropped).mean())
    item['display_zero_replacements'] = int(np.count_nonzero(~np.isfinite(cropped)))
manifest.update(shape=list(shape), frame_reference='12m+7m height; width expanded to card aspect ratio',
    crop_xyxy=[x0,y0,x1,y1], nan_display_value=0.,
    colorbars_shown=False, image_mm=[85.6,54.], border_mm=5.,
    scale=dict(style='dual_colour', background_vmin=-21.,transition=233.,highlight_vmax=1500.,
               background_cmap='Greys',background_stretch='linear',highlight_cmap='truncated mako',
               highlight_stretch='square root',unit='MJy/sr'))
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':5.2,'pdf.fonttype':42,'svg.fonttype':'none'})
page_w,page_h=95.6,64.
sky_w,sky_h=85.6,54.
left,bottom=5.,5.
geometry={};previews=[]
for number,item in enumerate(manifest['inputs'],1):
    data=fits.getdata(OUT/'grid'/f'{item["key"]}.fits')
    display=np.nan_to_num(data,nan=0.,posinf=0.,neginf=0.)
    if item['key'] == 'tp':
        display = display * tp_factor
    rgb=dual_colour_rgb(display,-21.,233.,1500.)
    fig=plt.figure(figsize=(page_w/25.4,page_h/25.4),facecolor='white')
    ax=fig.add_axes([left/page_w,bottom/page_h,sky_w/page_w,sky_h/page_h])
    ax.imshow(rgb,origin='lower',interpolation='none',aspect='auto');ax.set_axis_off();ax.set_autoscale_on(False)
    stroke=[pe.withStroke(linewidth=1.1,foreground='white')]
    ax.text(.012,.975,item['label']+'  ·  1.3 mm',transform=ax.transAxes,color='black',
            va='top',fontsize=6.2,path_effects=stroke)
    bmaj,bmin=item['beam_arcsec']
    beam_label = f'Beam {bmaj:.2f}″ × {bmin:.2f}″'
    if item['key'] == 'tp':
        beam_label += f'  ·  display ×{tp_factor:g}'
    ax.text(.012,.905,beam_label,transform=ax.transAxes,
            color='black',va='top',fontsize=5.2,path_effects=stroke)
    length=.02*3600/(abs(wcs.wcs.cdelt[0])*3600)
    end,y=shape[1]-65,70
    ax.plot([end-length,end],[y,y],color='black',lw=1.1,path_effects=stroke)
    ax.text(end-length/2,y+25,'1 pc',ha='center',va='bottom',color='black',
            fontsize=5.2,path_effects=stroke)
    prefix=OUT/f'{number:02d}_{item["key"]}_card'
    for fmt in ('png','tiff','pdf'):
        kwargs={'pil_kwargs':{'compression':'tiff_lzw'}} if fmt=='tiff' else {}
        fig.savefig(prefix.with_suffix('.'+fmt),dpi=1200,facecolor='white',**kwargs)
    preview=OUT/f'{number:02d}_{item["key"]}_preview.png'
    fig.savefig(preview,dpi=200,facecolor='white');previews.append(preview)
    fig.canvas.draw()
    bounds=ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted()).bounds
    assert np.allclose(np.array(bounds)*25.4,[5.,5.,85.6,54.])
    geometry[item['key']]=dict(image_bounds_mm=(np.array(bounds)*25.4).tolist(),colorbars_shown=False)
    plt.close(fig);print(f'Exported {prefix.name}',flush=True)
(OUT/'geometry.json').write_text(json.dumps(geometry,indent=2)+'\n')
fig,axes=plt.subplots(3,2,figsize=(10,10),facecolor='#eeeeee')
for ax,path in zip(axes.flat,previews):
    ax.imshow(plt.imread(path));ax.set_axis_off()
axes.flat[-1].set_axis_off()
fig.subplots_adjust(left=.02,right=.98,top=.98,bottom=.02,hspace=.04,wspace=.04)
fig.savefig(OUT/'contact_sheet.png',dpi=180);plt.close(fig)
