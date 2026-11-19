# Free-free contamination of the cores

Paper II referee item R7 (closed). Estimates the smooth radio (free-free)
emission inside each core footprint and its effect on core and fragment
masses. Split out of `construct_freefree_model` on 2026-11-20, because it
depends on later tasks: `source_catalogue`, `physical_properties` and
`core_formation`.

## Status: not reviewed, does not run

Its ALMA inputs were removed or renamed in `data/delivered/alma/` on
2026-11-19. They will be repointed when this task is reviewed, after
`core_formation`:

| Script | Outdated input | Replacement to check |
|---|---|---|
| `core_freefree_contamination.py`, `plot_core_freefree_contamination.py` | `continuum_image` (`…7MTM2TM1.robust0.5.clean.drogdeg.fits`) | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.clean.fits` |
| `plot_core_freefree_contamination.py` | `7MTM2TM1/cont/casa/primary_beam_rp2miriad.fits` (CASA pb) | MIRIAD `IRDC316.8.7MTM2TM1.robust0.5.cont.pb.fits` |
| `freefree_mass_efficiency_bias.py` | `7M/cont/IRDC316.8.ACA.p1.cont.{image,pb}.fits` | `data/delivered/alma/continuum/IRDC316.8.7M.robust0.5.cont.{image,pb}.fits` |

Paths to `results/`, `extraction/` and the other task configurations are
still written in the scripts, not in `task.toml`. That will also be fixed at
the review, together with merging the three scripts into one `run.py`
with stages.

## Method (adopted profile)

1. `core_freefree_contamination.py` extrapolates the native CORNISH 5.5 GHz
   map to the `physical_properties` frequency (223.8 GHz). It integrates it
   within the saved getsf 12m+7m footprint ellipses, with a uniform CORNISH
   rms of 1.0 mJy/beam. Detections need native pixels > 3σ; nondetections get
   1σ aperture limits. The 12m+7m image, including its local background,
   integrated over the same footprint, is the denominator. Settings are in
   `[core_contamination]`.
2. `plot_core_freefree_contamination.py` draws the spatial overlay and the
   flux comparison (PDF/SVG, PNG/TIFF, caption, alignment and collision QA).
3. `freefree_mass_efficiency_bias.py` matches the radio template to the ACA
   beam and measures it in the retained fragment leaf masks. It then compares
   the mass-fraction sensitivity of the evolved-region fragments and cores.
   The adopted catalogues are not changed.

## Outputs

`results/freefree_core_contamination/` (moved from
`results/construct_freefree_model/core_contamination/`; some QA files still
name the old path). `rms_1mjy/` and `rms_1p5mjy/` hold earlier profiles;
several files at the top level (`*_5mjy_*`, `*_3sigma_*`) come from profiles
that were superseded. They will be sorted at the review.
