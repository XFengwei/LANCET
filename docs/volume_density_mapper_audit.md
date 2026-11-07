# Volume-density mapper execution audit

Audit date: 2026-09-03

## Scope and software

The tested input is
`results/column_density/maps/surface_density_allscale.fits`, a 5127 by 7096
float32 surface-density image in g cm-2 on 0.1-arcsec pixels. The adopted
distance is 2.8 kpc.

The test used `gxli/volume-density-mapper` commit
`0495302cb3eb229892deb7b6871e651d125ca206` (package version 0.1.4) with
`constrained-diffusion` 1.2.5, NumPy 2.5.1, SciPy 1.18.0, and Astropy 8.0.1.
Delivered inputs were read only.

The map does not have one spatial resolution. The finite mask of the
primary-beam-cut 12m+7m+TP continuum image defines the 0.65-arcsec ALMA
footprint. Pixels outside that footprint use the retained 18-arcsec
surface-density map. The transition is flagged out to one 18-arcsec beam on
both sides.

## Central 1500 by 1500 pilot

The exact requested centre, ICRS `14:44:58.79 -59:48:01.94`, maps to the
zero-based cutout bounds `[y: 886:2386, x: 2536:4036]`. The cutout covers 150
arcsec per side and is 60.8323% inside the ALMA footprint.

The native 2-D mapper run used no artificial zero padding and completed in
79.33 s with about 1.00 GiB maximum resident memory. Twenty-one negative
surface-density pixels were clipped to zero for the physical calculation;
the original values remain in `surface_density_input.fits`.

The unmodified upstream mean-density output is not safe to use directly:
65,804 pixels (2.925%) received zero characteristic width. Its call to
`numpy.nan_to_num` consequently replaces division-by-zero values with the
largest float64 number, which then overflows when written as float32. The
LANCET wrapper retains the raw zero widths for diagnosis, floors the effective
width at the local 0.65- or 18-arcsec beam, recomputes density, and records both
conditions in a bit mask.

Away from the one-beam resolution-transition zone, the pilot medians are:

| Region | Effective width | H2 number density |
| --- | ---: | ---: |
| 0.65-arcsec ALMA footprint | 0.236 pc | 3.81e5 cm-3 |
| 18-arcsec outer region | 0.321 pc | 5.11e4 cm-3 |

These are method-dependent mean-density estimates, not independent
line-of-sight depth measurements.

## 3-D closure gate

A 250 by 250 by 250 smoke cube was reconstructed from an exact 6 by 6 block
average of the pilot. The cube occupies about 62.5 MB and completed in 0.70 s.
After projection along depth, the recovered-to-input column ratio has
percentiles 0.000, 0.000, 0.389, 0.805, and 0.963 at 0, 1, 50, 99, and 100%.

The median missing column is therefore about 61%. Inspection of the upstream
implementation confirms the cause: constrained diffusion returns decomposed
layers plus a final residual, but `density_reconstruction_3d` projects only
the layers. The residual is not placed into the cube.

The declared gate requires at least 95% projected-column recovery at the
first percentile. It fails, so no native 3-D pilot or full-map cube was
generated. This also avoids impractical products: float32 storage alone would
be 13.5 GB for the native pilot and 1.033 TB for the full map.

## Resolution-aware full map

The full 2-D calculation uses two independent branches:

1. The 18-arcsec retained surface-density map is analyzed on its native
   6-arcsec pixels and reprojected to the target grid.
2. The all-scale map over the ALMA footprint is analyzed at 0.3-arcsec pixels,
   a factor-three block average that still samples the 0.65-arcsec beam with
   2.17 pixels per FWHM.

The branches are combined only after decomposition. This prevents the mapper,
which accepts only one scalar pixel size and no resolution map, from treating
the mixed-resolution image as a uniform observation. The low- and
high-resolution mapper runs completed in 0.78 s and 158.11 s, respectively.

For pixels outside the transition zone and with surface density above
1e-3 g cm-2, the full-map diagnostics are:

| Region | Pixels | Width 5/50/95% | n(H2) 5/50/95% |
| --- | ---: | ---: | ---: |
| 0.65-arcsec branch | 1,367,944 | 0.179 / 0.365 / 0.591 pc | 8.80e4 / 2.32e5 / 7.54e5 cm-3 |
| 18-arcsec branch | 24,051,621 | 1.625 / 4.130 / 7.288 pc | 4.89e2 / 3.06e3 / 1.73e4 cm-3 |

The 1e-3 g cm-2 threshold is used only for these diagnostic summaries and
figure contrast. It is not applied to the output FITS maps.

## Pilot-to-full consistency check

Within 631,312 high-resolution pixels that avoid both transition masks, the
full-map and pilot log-width maps have Spearman correlation 0.696 and the
log-density maps have correlation 0.852. The full/pilot width ratio has
16/50/84% values 1.362/1.629/2.272; the density ratio is
0.440/0.614/0.734.

This systematic shift is expected from the full run's larger spatial context
and factor-three rebinning, but it is scientifically material. Absolute width
and density values should therefore be treated as preliminary until a
rebin/context convergence study is completed. Morphology is more stable than
normalization.

## Products and quality mask

Pilot products are under `results/volume_density_mapper/pilot/`; full-map
products are under `results/volume_density_mapper/full/`. Both directories
contain raw and beam-floored width maps, mass and H2 number-density maps, a
local-resolution map, a summary JSON file, and a diagnostic figure. The pilot
also contains the downsampled smoke cube.

The uint8 quality mask uses bit 0 for the ALMA footprint, bit 1 for a finite
positive upstream width, bit 2 when the beam floor is applied, and bit 3 for
the resolution-transition exclusion zone. Analyses comparing the two
resolution regimes should exclude pixels with bit 3 set.

## Primary-beam masking update

Update date: 2026-09-04

Visual inspection of the two-branch full-map result showed the expected
strong discontinuity where the 0.65-arcsec and 18-arcsec decompositions met.
At the user's request, the final science products now use the numerical
primary-beam image
`data/delivered/alma/7MTM2TM1/cont/casa/primary_beam_rp2miriad.fits` and retain
only pixels with response greater than or equal to 0.20. This is the same
cutoff configured by `prepare_continuum_images`.

The numerical PB >= 0.20 selection and the independently generated finite
PB-cut continuum footprint agree exactly: both contain 4,546,473 pixels, with
zero pixels in their exclusive-or difference. The full mask retains 12.4968%
of the 5127 by 7096 grid. The central pilot retains 1,368,726 of 2,250,000
pixels, or 60.8323%.

Recommended outputs are under
`results/volume_density_mapper/full_pbmasked/` and
`results/volume_density_mapper/pilot_pbmasked/`. For every floating-point
science product, all rejected pixels were verified to be NaN; the uint8
quality mask is zero outside the PB selection. The previous unmasked
two-branch products are retained for method diagnostics but should not be used
as the default science maps.
