# Continuum image workflow audit

## Retained scientific workflow

The continuum-image sequence separates into three reproducible tasks:

1. `prepare_continuum_images` applies primary-beam cutoff masks and converts
   Jy/beam images to MJy/sr.
2. `compare_interferometer_arrays` reads declared measurement sets for UV
   coverage, compares ACA and 12 m continuum images, and performs the
   beam-matched consistency check.
3. `plot_allscale_continuum` produces the interferometer plus total-power dust
   and dust-plus-free-free figures.

This ordering follows the actual data dependencies rather than the order in
which the calculations were explored.

## Implemented preparation task

`prepare_continuum_images` declares four primary-beam masks and eleven
brightness-unit conversions. Its reusable functions live in
`src/lancet/imaging.py`.

The unit conversion is

```text
I [MJy/sr] = S [Jy/beam] / Omega_beam [sr] * 1e-6
Omega_beam = pi / (4 ln 2) * BMAJ * BMIN
```

The primary-beam operation retains pixels with response greater than or equal
to 0.20. It is a cutoff mask, not a primary-beam flux correction.

Regression checks reproduce the retained masked arrays exactly. Unit-converted
arrays are calculated identically and then stored as float32 to reduce disk
use. Across the current conversion set, the maximum relative difference from
the earlier float64 storage is below `6e-8`.

## Multi-array comparison

The `compare_interferometer_arrays` task is implemented with `inspect`, `uv`,
`compare-images`, and `beam-match` stages.

The retained products are:

- the composed UV-coverage manuscript figure;
- the vertical ACA, 12 m, and 12 m plus 7 m image comparison;
- convolution of the combined interferometer image to the ACA beam;
- reprojection to the ACA grid; and
- a residual diagnostic against the ACA image.

The four measurement sets under
`data/reduction/alma/7MTM2TM1/cont/ms/` are declared read-only intermediate
inputs used only to illustrate UV coverage. They are not scientific delivered
products and will never be modified by the plotting task.

The final UV manuscript figure is a composed two-panel layout. Its component
plots and the saved manuscript PDF establish the intended geometry, including
the zoom connectors and the short-spacing overlap annotations.

UV coordinates are read in bounded chunks, rounded to 0.1 m for duplicate
removal, and cached as compact float32 arrays. Dense marker layers are
rasterized in the PDF while annotations, axes, overlap circles, and connectors
remain vector graphics.

## All-scale continuum display

The `plot_allscale_continuum` task is implemented with `inspect`, `combined`,
`synthetic`, `dust`, and `all` stages. Its reusable dual-colour plotting
functions live in `src/lancet/continuum_plotting.py`.

The retained outputs are:

- the dust-plus-free-free dual-color manuscript PDF and clean PNG;
- the dust-only dual-color manuscript PDF and clean PNG; and
- the clean 230 GHz synthetic-continuum image where it remains useful.

The earlier primary-beam reprojection calculation is unused by the displayed
image and will not be retained. The displayed primary-beam boundary comes from
the declared DS9 region.

The delivered 7096-by-7096 all-scale maps are read-only float64 products. The
plotting task reads them with memory mapping and writes only PDF and PNG
figures, so it does not duplicate their FITS storage.

## Transitional inputs

The delivered all-scale `.MJysr.rp.fits` files predate the task framework.
Their primary-beam masks are now reproducible, but their earlier reprojection
step is not present in this source sequence. They remain declared stable inputs
until that upstream transformation is recovered or replaced.

## Removed duplication

- The vertical interferometer comparison is the manuscript product; the
  alternative 2-by-2 layout is a redundant presentation.
- One nominal low-resolution and high-resolution panel used the same FITS
  path, so it is not a valid four-way comparison.
- Standalone interactive image displays without saved quantitative output are
  replaced by explicit residual statistics and optional diagnostic figures.
- Repeated plotting helpers, style blocks, unused imports, and unused distance
  calculations will be consolidated.
