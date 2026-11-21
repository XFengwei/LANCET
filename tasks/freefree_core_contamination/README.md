# Free-free contamination of the cores

Paper II referee item R7 (closed). Estimates the smooth radio (free-free)
emission inside each core footprint and its effect on core and fragment
masses. Split out of `construct_freefree_model` on 2026-11-20 because it
depends on later tasks (`source_catalogue`, `core_physical_properties`,
`core_formation_efficiency`). All inputs and settings are in `task.toml`.

## Inputs

| Key | File |
|---|---|
| `cornish_map` | `data/delivered/archive/radio/CORNISH_5p5GHz.fits` (native 2.5″) |
| `source_catalogue` | `results/source_catalogue/tables/source_catalogue.ecsv` |
| `getsf_footprints` | getsf 12m+7m footprint ellipses (`…G316.sw.sources.ok.cat.1302.foo.fk5.reg`) |
| `continuum_image` | `data/delivered/alma/continuum/IRDC316.8.7MTM2TM1.robust0.5.cont.clean.truesky.fits` |
| `primary_beam` | `…7MTM2TM1.robust0.5.cont.pb.fits` (display cut, pb ≥ 0.2) |
| `aca_image`, `aca_pb` | `…IRDC316.8.7M.robust0.5.cont.{image,pb}.fits` (astrodendro fragment image) |
| `fragment_label_mask`, `fragment_catalogue` | `data/extraction/astrodendro/my_mask_idx.fits`, `my_dendrogram.dat` |
| `fragment_selection`, `core_formation_summary` | `results/core_formation_efficiency/tables/` |
| `physical_catalogue` | `results/core_physical_properties/tables/core_physical_properties.ecsv` |

## Method (adopted profile)

1. `contamination`: the native CORNISH 5.5 GHz map is extrapolated to
   223.8 GHz (optically thin, α = −0.1) and integrated within the getsf
   12m+7m footprints. Native pixels above 3σ of a uniform 1.0 mJy/beam rms
   are retained; covered nondetections get 1σ aperture limits. The 12m+7m
   true-sky aperture flux, including its local background, is the
   denominator. Writes `smooth_core_contamination.{ecsv,csv}`, `summary.json`
   and the radio noise and retained maps.
2. `plot`: spatial overlay of the evolved region and the flux comparison
   (`cornish_core_comparison.{pdf,svg,png,tiff}`, `figure_caption.txt`,
   `radio_noise_map`). `--qa-scripts <folder>` runs the optional
   panel-alignment audit (`audit_panel_alignment.py`).
3. `mass_bias`: the radio template is matched to the 7m beam and measured in
   the retained fragment leaf masks; the fragment and core mass ratios are
   recomputed as a sensitivity scenario (`fragment_freefree_contamination`,
   `mass_efficiency_bias.json`). The adopted catalogues are not changed.

```bash
PYTHONPATH=src conda run -n py312 python tasks/freefree_core_contamination/run.py --stage all --overwrite
```

Stages: `contamination`, `plot`, `mass_bias`, `all`. Existing products are
kept unless `--overwrite` is given.

## Review (2026-11-21)

- The three scripts were merged into `run.py`; every path is now in
  `task.toml` (the frequency is set here instead of being read from another
  task's configuration), and the alignment audit no longer needs a local
  skill path.
- Inputs: the removed `drogdeg` 12m+7m image is replaced by the MIRIAD
  true-sky image, matching the gain-corrected core masses; the CASA beam by
  the MIRIAD `pb`; the removed `ACA.p1` 7m image by the delivered
  `7M.robust0.5.cont.image`, which reproduces the fragment catalogue fluxes
  to 1e-5.
- Check on the Paper II catalogue: every radio-side column (free-free flux,
  limits, coverage, retained area) and both radio maps are identical; only
  the 12m+7m denominator changes.
- Rerun on the current catalogues: evolved region 94 covered cores,
  32 with radio, 13 above 10 %, none above 100 %; flux-weighted fraction
  0.0517 → 0.0516. Scenario mass reductions unchanged: cores 3.3 %,
  fragments 9.9 % (evolved); the core/fragment ratios follow the new core
  masses (Young 0.534 → 0.556, Intermediate 0.384 → 0.392, Evolved 0.187).
- Earlier noise profiles (`rms_1mjy/`, `rms_1p5mjy/`, the `*_3sigma_*` and
  `*_5mjy_*` maps) and the old alignment audits are kept in
  `results/freefree_core_contamination/superseded/`.
