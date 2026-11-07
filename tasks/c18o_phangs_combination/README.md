# C18O PHANGS combination pilot

This task tests the current PHANGS-ALMA pipeline on the existing continuum-
subtracted IRDC 316.8 C18O(2-1) inputs. It jointly images the 12m and 7m
measurement sets, then feathers the LSRK TP cube. All outputs are isolated
under `data/reduction/alma/TP7MTM2/cube/phangs_test`; the legacy MIRIAD tree is
never used as an output.

The checked-in smoke-test settings use a 0.34 km/s window centered at -40 km/s
and request 0.085 km/s channels. The TP 95th-percentile spectrum peaks near
-39.83 km/s. This keeps the test much smaller than the 470-channel MIRIAD
product. For a full-cube run, change `velocity_width_kms` to `40.0` only after
the smoke test passes.

The spatial image grid follows the legacy MIRIAD `combined.clean.fits`: a SIN
projection with 0.4 arcsec pixels. CASA's mosaic FFT cannot process MIRIAD's
odd 1716 x 1225 shape, so the CASA grid is 1716 x 1226 with one padding row.
Its reference pixel is `J2000 14h44m55.306994s -59d46m28.813092s`, chosen so
the original MIRIAD pixel centers remain aligned. PHANGS dynamic sizing is
disabled so it cannot replace this grid.

The staged MS is already on the requested LSRK spectral grid, so tclean uses
nearest-neighbor channel selection. CASA's linear interpolation otherwise
blanks the final channel because an additional sample beyond that edge is not
available.

The generated clean recipe also writes the configured channel count
explicitly. `analysisUtils.ValueMapping` assumes that every four-channel
spectral window is WVR data, so automatic inference would discard this valid
four-channel science pilot before the PHANGS clean loop starts.

For the single-scale signal mask, the runner uses PHANGS's direct CASA image
write path. Its large-cube FITS workaround mixes radian-valued `imhead` WCS
numbers with `CUNIT='deg'`; this small pilot mask fits comfortably in memory,
so avoiding that round trip preserves an image-aligned mask WCS.

The pilot clean policy is configured separately under `[cleaning]`. It builds
a positive-only mask with 3.5-sigma seeds grown through connected 1.5-sigma
emission, clips that deep-clean mask at PB > 0.3, and initially cleans to
1.5 sigma with a one-percent flux-convergence criterion. The image itself
retains the broader `pblimit = 0.2` footprint. Spectral rolling is disabled for
the four-channel pilot because `numpy.roll` would otherwise make its first and
last channels artificial neighbours. These controls apply to PHANGS's signal
mask; the auto-multithresh-looking values in the generated tclean recipe do
not define the final single-scale mask because that stage uses `usemask=user`.

The independently configured PSF center is
`J2000 14h44m56.8141s -59d47m20.1360s`, the position of 12m field 67 and only
6.90 arcsec from ACA field 86. The earlier PSF center
(`14h45m00.0s -59d45m00.0s`) lay outside the convex hull of the 12m pointing
centers and was 142 arcsec from its nearest 12m field. In this asymmetric
filament mosaic that caused CASA to place the PSF peak outside the inner image
quarter and return invalid restoring beams.

## Requirements

PHANGS master currently documents Python 3.12, CASA >= 6.7.3, and
`analysisUtils`. Install the CASA extras described by the
[PHANGS repository](https://github.com/PhangsTeam/phangs_imaging_scripts).
The inspection and key-generation stages do not perform imaging; run the
imaging stages only in an environment with both CASA and `analysisUtils`.

## Safe test sequence

Inspect inputs without writing anything:

```bash
conda run -n py312 python scripts/c18o_phangs_combination.py --stage inspect
```

Generate the absolute-path PHANGS key set:

```bash
conda run -n py312 python scripts/c18o_phangs_combination.py --stage keys
```

Validate those keys against a local checkout of the PHANGS repository:

```bash
conda run -n py312 python scripts/c18o_phangs_combination.py \
  --stage validate \
  --pipeline-root /path/to/phangs_imaging_scripts
```

In the CASA-enabled PHANGS environment, run one stage at a time:

```bash
python scripts/c18o_phangs_combination.py \
  --analysis-utils /path/to/analysis_scripts --stage stage
python scripts/c18o_phangs_combination.py \
  --analysis-utils /path/to/analysis_scripts --stage dirty
python scripts/c18o_phangs_combination.py \
  --analysis-utils /path/to/analysis_scripts --stage image
python scripts/c18o_phangs_combination.py \
  --analysis-utils /path/to/analysis_scripts --stage postprocess
```

`stage` deliberately calls PHANGS with `do_contsub=False` because both input
MSs are already continuum-subtracted. It also sets `statwt_line=False` so the
short emission-containing pilot window does not replace the calibrated 12m/7m
relative weights. Staged inputs are retained for diagnosis.

To test revised single-scale mask settings without repeating the dirty and
multiscale stages, revert the current imaging stack to its retained
`_multiscale` checkpoint and clean from there:

```bash
python scripts/c18o_phangs_combination.py \
  --analysis-utils /path/to/analysis_scripts --stage reclean
```

This replaces the current unsuffixed, `_singlescale`, and `_prev` stacks but
leaves `_dirty` and `_multiscale` intact. Use `--stage image` instead when the
multiscale parameters or visibility weighting have changed.

Line extraction forces `clear_pointing=True` before the separately extracted
spectral windows are concatenated. The clean recipe uses `usepointing=False`,
so CASA obtains the mosaic centers from the `FIELD` table and the large
channel-independent `POINTING` table is unnecessary. Without this policy the
eight ACA spectral windows each contribute a copy of 7,597,135 pointing rows,
inflating a four-channel combined MS to about 12 GB.

To replace only an existing oversized combined MS while reusing the retained
12m and 7m staging files, run:

```bash
python scripts/c18o_phangs_combination.py \
  --analysis-utils /path/to/analysis_scripts \
  --stage extract --overwrite
```

The `dirty` stage is the first real combination smoke test. Inspect its PSF,
beam, channel RMS, and mosaic footprint before launching `image`. The latter
runs PHANGS multiscale plus single-scale cleaning and can consume substantial
time and disk space. `postprocess` performs primary-beam correction, prepares
the TP cube on the interferometric grid, feathers without apodization, converts
units, trims, and exports the final FITS products.
