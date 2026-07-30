# Prepare continuum images

This task provides the common image-preparation operations used by the
multi-array continuum analysis.

## Operations

- `mask` applies a primary-beam response cutoff of 0.20. It does not divide
  the image by the primary beam and therefore does not perform a flux
  correction.
- `convert` transforms Jy/beam images to MJy/sr using the elliptical Gaussian
  beam area recorded by `BMAJ` and `BMIN`.
- `all` runs masking before unit conversion so converted products can consume
  freshly generated masked images.

Every delivered input is read-only. Generated FITS products are written under
`results/prepare_continuum_images/`.

The complete configured output is approximately 1.8 GiB. Floating-point FITS
products are stored as float32 (`BITPIX = -32`); the conversion is calculated
at float64 precision before this storage cast. Across the current images, the
maximum relative float32 rounding difference is below `6e-8`. Use
`--only NAME` to test or regenerate selected operations without creating every
product.

The two delivered all-scale `.MJysr.rp.fits` images are transitional stable
inputs: their earlier reprojection step is not part of this task. The mask is
recreated under `results/`; the delivered versions are never overwritten.

Run:

```bash
conda run -n py312 python scripts/prepare_continuum_images.py --stage inspect
conda run -n py312 python scripts/prepare_continuum_images.py \
  --stage mask --only 7MTM2
conda run -n py312 python scripts/prepare_continuum_images.py \
  --stage convert --only 7MTM2_pbcut
```

After the selected operations are verified, omit `--only` to run a complete
stage or use `--stage all` for the entire configured workflow.
