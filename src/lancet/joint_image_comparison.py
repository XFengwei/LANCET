"""Beam-matched CASA/MIRIAD cutouts for the referee response."""
from pathlib import Path
import hashlib
import json
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from radio_beam import Beam
from reproject import reproject_interp
from .interferometry import beam_match_and_reproject


def run_joint_image_comparison(root, config, *, save_figure, overwrite=False):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from astropy.visualization import ImageNormalize, SqrtStretch
    import seaborn as sns
    from matplotlib.patches import Ellipse

    settings = config['joint_image_comparison']
    out = root / config['outputs']['root'] / 'joint_image_comparison'
    out.mkdir(exist_ok=True, parents=True)
    stem = out / 'CASA_MIRIAD_central_south'
    if stem.with_suffix('.pdf').exists() and not overwrite:
        raise FileExistsError(stem)
    paths = {name: root/settings[name] for name in ['casa', 'miriad', 'casa_pb']}
    arrays, headers = {}, {}
    for name, path in paths.items():
        arrays[name], headers[name] = fits.getdata(path, header=True)
        arrays[name] = np.squeeze(arrays[name])
    for name in ['casa', 'miriad']:
        assert headers[name]['BUNIT'].strip().lower() == 'mjy/sr'
    beam = Beam.from_fits_header(headers['miriad'])
    kernel = beam.deconvolve(Beam.from_fits_header(headers['casa']))
    plt.rcParams.update({'font.size': 10, 'font.family': 'sans-serif',
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 2, figsize=(settings['width_inches'], 6.6))
    fig.subplots_adjust(left=.10, right=.82, bottom=.09, top=.98, wspace=.12, hspace=.20)
    scale = config['image_comparison_plot']
    vmin, vmax = float(scale['vmin_mjysr']), float(scale['vmax_mjysr'])
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())
    colormap = sns.color_palette(scale['colormap'], n_colors=256, as_cmap=True)
    stats = []
    for row, (name, centre) in enumerate(zip(settings['region_names'], settings['region_centres_deg'])):
        position = SkyCoord(*centre, unit='deg', frame='fk5')
        size = settings['cutout_arcsec']*u.arcsec
        ref = Cutout2D(arrays['miriad'], position, size, wcs=WCS(headers['miriad']).celestial)
        source = Cutout2D(arrays['casa'], position, size+20*u.arcsec, wcs=WCS(headers['casa']).celestial)
        rh, sh = ref.wcs.to_header(), source.wcs.to_header()
        for key in ['BMAJ', 'BMIN', 'BPA', 'BUNIT']:
            rh[key], sh[key] = headers['miriad'][key], headers['casa'][key]
        matched, _ = beam_match_and_reproject(source.data, sh, ref.data, rh)
        pbcut = Cutout2D(arrays['casa_pb'], position, size+20*u.arcsec, wcs=WCS(headers['casa_pb']).celestial)
        pb, _ = reproject_interp((pbcut.data, pbcut.wcs), ref.wcs, shape_out=ref.shape)
        common = np.isfinite(matched) & np.isfinite(ref.data) & (pb >= settings['minimum_pb'])
        left, right = np.where(common, matched, np.nan), np.where(common, ref.data, np.nan)
        pix = float(proj_plane_pixel_scales(ref.wcs)[0]*3600)
        extent = [ref.shape[1]*pix/2, -ref.shape[1]*pix/2, -ref.shape[0]*pix/2, ref.shape[0]*pix/2]
        for col, (label, data) in enumerate([('CASA', left), ('MIRIAD', right)]):
            ax = axes[row, col]
            im = ax.imshow(data, origin='lower', extent=extent, cmap=colormap, norm=norm, rasterized=True)
            ax.text(.04, .95, f'({"abcd"[row*2+col]}) {name} · {label}', transform=ax.transAxes,
                    va='top', color='black', fontsize=10,
                    bbox={'facecolor':'white', 'edgecolor':'none', 'alpha':.85, 'pad':2})
            half = settings['cutout_arcsec']/2
            ax.set(xlim=(half,-half), ylim=(-half,half), xticks=[40,0,-40], yticks=[-40,0,40])
            ax.tick_params(direction='in', top=True, right=True)
            if row == 1:
                ax.set_xlabel('RA offset (arcsec)')
            if col == 0:
                ax.set_ylabel('Dec offset (arcsec)')
            else:
                ax.tick_params(labelleft=False)
            ax.add_patch(Ellipse((40,-25), beam.minor.to_value(u.arcsec), beam.major.to_value(u.arcsec),
                                 angle=-beam.pa.to_value(u.deg), color='white'))
            ax.plot([-20, -20-14.7332], [-30,-30], color='white', lw=2)
            ax.text(-27.4,-25,'0.2 pc',color='white',ha='center',fontsize=9)
            h = rh.copy()
            h.add_history('Common PB>=0.2 and finite-data mask; CASA convolved to MIRIAD beam in surface-brightness units.')
            fits.writeto(out/f'{name.lower()}_{label.lower()}_matched.fits',data,h,overwrite=True)
        bounds=axes[row,1].get_position()
        cax=fig.add_axes([.85,bounds.y0,.018,bounds.height])
        ticks = [0,100,200,300,400,500]
        cb = fig.colorbar(im,cax=cax,label='MJy/sr', ticks=ticks)
        cb.set_ticklabels([str(v) for v in ticks])
        delta = left[common]-right[common]
        stats.append({'region':name,'centre_fk5_deg':centre,'common_pixels':int(common.sum()),
                      'total_pixels':int(common.size),'display_vmin':vmin,'display_vmax':vmax,'display_stretch':'sqrt','colormap':scale['colormap'],
                      'median_casa_minus_miriad_mjysr':float(np.median(delta)),
                      'rms_difference_mjysr':float(np.sqrt(np.mean(delta**2)))})
    save_figure(fig,[stem.with_suffix(e) for e in ['.pdf','.svg','.png']],dpi=400,bbox_inches=None,
                exclude_axes=fig.axes[4:])
    plt.close(fig)
    report={'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},
            'target_beam_arcsec':[beam.major.to_value(u.arcsec),beam.minor.to_value(u.arcsec)],
            'target_pa_deg':beam.pa.to_value(u.deg),
            'convolution_kernel_arcsec':[kernel.major.to_value(u.arcsec),kernel.minor.to_value(u.arcsec)],
            'processing':'PB-corrected MJy/sr; CASA convolved before reprojection. Identical common finite/PB>=0.2 mask and per-row colour scale; no flux rescaling.',
            'scope':'Comparison of existing reductions; imaging settings can differ. Dirty PSF is assessed independently.',
            'regions':stats}
    stem.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
