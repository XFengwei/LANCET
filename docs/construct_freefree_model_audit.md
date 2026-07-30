# Free-free model data and provenance audit

Audit date: 2026-07-28

The historical `Archive/Radio/` paths now resolve under
`data/delivered/archive/radio/`.

## FITS classification

| FITS file | Classification | Role or transformation |
| --- | --- | --- |
| `Tile724_CORNISH_5GHz.fits` | Upstream CORNISH-South tile | Parent survey tile |
| `CORNISH_5p5GHz.fits` | Post-processed upstream cutout | Read-only input |
| `CORNISH_5p5GHz_8arcsec_native.fits` | Legacy derived product | Smooth 2.5 to 8 arcsec and rescale Jy/beam |
| `CORNISH_5p5GHz_8arcsec_on_MeerKATgrid.fits` | Legacy derived product | Reproject smoothed CORNISH map to 1200 x 1200 MeerKAT grid |
| `MeerKAT_subcube_reprojected.fits` | Post-processed upstream MeerKAT cube | Contains repaired 8-arcsec beam cards |
| `MeerKAT_subimage_1.38GHz_reprojected.fits` | Post-processed upstream MeerKAT image | Contains repaired beam cards and a two-dimensional header |

The 301 x 301 `CORNISH_5p5GHz.fits` image was verified to be an exact
zero-difference cutout of `Tile724_CORNISH_5GHz.fits`, corresponding to
zero-based tile pixels `x=1174:1475` and `y=341:642`.

The new task's preparation stage independently reproduced both legacy CORNISH
derivatives. The native 301 x 301 and reprojected 1200 x 1200 arrays have
identical finite-pixel masks and a maximum absolute difference of exactly zero
from the files currently stored under `data/delivered/archive/radio/`.

The two MeerKAT files already have the repaired 8-arcsec beam metadata. Their
pre-reprojection source data and their earlier headers
are not present here, so the current files are retained as stable delivered
inputs. New task runs must not overwrite them.

## Historical task outputs

The historical workflow wrote eight fit products for each of 6000, 8000, and
10000 K:

- `EM_map.fits`
- `EMerr_map.fits`
- `ff230_map.fits`
- `chi2_map.fits`
- `rchi2_map.fits`
- `nband_map.fits`
- `fitmask_map.fits`
- `model_subbands_cube.fits`

For the 8000 K model it additionally writes:

- `ff230_map_smoothed_18arcsec.fits`
- `ff230_map_smoothed_18arcsec_regrid2dust.fits`

These 26 products were migrated intact from `sedfitter/` on 2026-07-30 and
remain recorded in `data/path_ledger.tsv`. The delivered
copies were subsequently removed after all corresponding products were
verified under `results/construct_freefree_model/`.

## Reproducibility findings

- Exact MeerKAT subband frequencies are stored in FITS cards `FREQ0001` to
  `FREQ0014`. For strict Paper II reproduction, the task profile uses the
  rounded Paper II list; the FITS-header frequencies remain available
  through `frequency_source = "fits_header"`.
- The transferred standalone prototype contained an older example
  configuration: `1e-5 Jy/beam` for MeerKAT, `2e-5 Jy/beam` for CORNISH, no
  S/N cut, and a different approximate frequency list. Its verified fitting
  functions now live in `src/lancet/freefree.py`, but those example settings
  did not generate the stored Paper II maps.
- The stored Paper II products use `3e-4 Jy/beam` for every nominal MeerKAT
  subband, `1e-2 Jy/beam` for CORNISH, and a 5-sigma per-band cut.
- The stored fit mask is reproduced exactly by requiring finite MFS and alpha
  values, `MFS > 5e-3 Jy/beam`, and `-0.5 <= alpha <= 2.0`. This confirms that
  the alpha cut was active when the legacy products were generated.
- The Paper II map fit uses a 5-sigma per-band cut.
- The model cube has 15 irregularly spaced bands, including CORNISH at
  5.5 GHz; it should be saved with explicit band-frequency metadata rather
  than pretending it has a linear FITS frequency axis.

## RMS provenance check

The numerical RMS values are not measured by the fitting implementation; both
the old example configuration and the Paper II profile use fixed constants.

Forensic checks made before removal of the historical 8000 K products gave:

| Check | Paper II profile | Older script profile |
| --- | ---: | ---: |
| Fit-mask pixel mismatch | 0 | 232,181 |
| Band-count-map pixel mismatch | 0 | 235,729 |
| Median absolute EM difference in 200 sampled fitted pixels | 0 | 2.4% |
| Maximum absolute EM difference in the sample | 0 | 4.5% |

The delivered MeerKAT cube has two entirely NaN subband planes. The stored
band-count map therefore contains 23,766 pixels fitted with 12 MeerKAT bands
and 304 pixels fitted with those 12 bands plus CORNISH.

The fixed MeerKAT value of 0.3 mJy/beam is comparable to robust local
100 x 100-pixel background estimates in several subbands, although the values
vary by band and location. The fixed CORNISH value of 10 mJy/beam is not its
thermal noise: CORNISH-South reports about 0.11 mJy/beam at 2.5-arcsec
resolution, and the delivered cutout has a local robust scatter around
0.10-0.20 mJy/beam before smoothing. It functions as a conservative effective
uncertainty and bright-central-region selection threshold.
