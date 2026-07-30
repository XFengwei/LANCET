# SED-fitter `combine/` provenance audit

Audit date: 2026-07-30

The cleanup approved after this audit was completed on 2026-07-30. Nine FITS
derivatives and one `.DS_Store` file were removed, recovering 890,021,188
bytes (about 849 MiB). The historical migration manifest
was verified as fully merged into `data/path_ledger.tsv`, which remains the
record of every file's former root-level and migration-time path.

## Retain as delivered scientific products

These four files are adopted SED-fitting or multiresolution-combination
products. No project code was found that writes them, and current tasks use
three of them as read-only scientific inputs.

| File | Role |
| --- | --- |
| `NH2_18arcsec+36arcsec.fits` | Adopted low-resolution combined column density; active input |
| `NH2_6arcsec+18arcsec.fits` | Adopted intermediate-resolution combined column density |
| `NH2_6arcsec+18arcsec+36arcsec.fits` | Adopted multiscale column density; active input |
| `Tdust_18arcsec+36arcsec.fits` | Adopted low-resolution combined temperature; active input |

The 18+36 arcsec maps are not copies of the `version2/fitresult/` maps; their
pixel values and finite masks differ. They therefore contain scientific
combination information that is not replaced by retaining `version2/` alone.

## Removed derivatives with current replacements

The following files were written by the historical column-density analysis
and were removed on 2026-07-30. Current tasks regenerate them under
`results/column_density/`.

| Historical file | Current replacement | Verification |
| --- | --- | --- |
| `Tdust_18arcsec+36arcsec_regrid_allscale.fits` | `prepared/temperature_allscale.fits` | Pixel-identical, including NaN mask |
| `NH2_allscale.fits` | `maps/column_density_allscale.fits` | Same NaN mask; relative difference below 1.2e-7 for NH2 >= 1e20 cm^-2 |
| `NH2_intercore.fits` | `maps/column_density_intercore.fits` | Same NaN mask; relative difference below 1.2e-7 for NH2 >= 1e20 cm^-2 |
| `Sigma_18arcsec+36arcsec.fits` | `prepared/surface_density_lowres.fits` | Pixel-identical |
| `Sigma_6arcsec+18arcsec+36arcsec.fits` | `prepared/surface_density_multiscale.fits` | Pixel-identical |
| `Sigma_allscale.fits` | `maps/surface_density_allscale.fits` | Same NaN mask; float32 rounding only at scientifically meaningful values |
| `Sigma_intercore.fits` | `maps/surface_density_intercore.fits` | Same NaN mask; float32 rounding only |

The very small differences in the newly calculated all-scale column-density
maps come from the numerically stable `log1p` implementation and float32
rounding. Pixels above 1e20 cm^-2 agree to better than 1.2e-7 in relative
value.

## Other removed derivatives

| File | Classification | Recommendation |
| --- | --- | --- |
| `NH2_18arcsec+36arcsec_regrid_allscale.fits` | Historical reprojection | Removed. It is reproduced exactly by `hcongrid` from the retained low-resolution NH2 map and the delivered all-scale grid. |
| `Tdust_18arcsec+36arcsec_regrid_ACA7m.fits` | Historical reprojection | Removed by user decision. It was not an original SED-fit result and no active task read it. Its precise legacy interpolation recipe was not found, and its inherited CASA header incorrectly reported `Jy/beam` for temperature data. |
| `.DS_Store` | macOS metadata | Removed. |

The four retained scientific products occupy about 25 MiB.

## Free-free cleanup

The removed `data/delivered/sedfitter/FFmodel/` directory contained 26
historical task outputs: eight maps/cubes at each of 6000, 8000, and 10000 K,
plus two 8000 K post-processing products. All 26 corresponding products are
present under `results/construct_freefree_model/`. Active task configuration
and regression fixtures now use the `results/` products.
