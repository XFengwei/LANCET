# Core formation efficiency

This task measures how much cloud mass is contained in Paper I-scale
fragments and in the current compact-core catalogue. It also displays the
compact cores and larger fragment footprints on the continuum maps.

The cloud masses, dust opacities, physical scales, and plot settings are
declared in `task.toml`. Fragments are retained only when their entire dendrogram leaf footprints
lie within the MIRIAD 12m+7m mosaic primary beam (`cont.pb`) of at least 0.2; their full `mass_dense` values
from `data/extraction/astrodendro/my_catalog.dat` are summed by Paper I subregion.
The fixed Paper I regional totals in the configuration are retained as the
pre-selection reference, not used for the selected sums. Paper I used
an opacity of 0.9 cm2/g for the fragment masses. The task multiplies those
masses by 0.9/0.5 so that they use the 0.5 cm2/g opacity adopted here. Cloud
masses remain unchanged because their SED fitting used the same opacity law.
Core masses and core counts are calculated from the generated
physical-properties catalogue. The summary reports fragment/cloud,
core/cloud, and direct fragment-to-core mass efficiencies. The uncertainty
beside the direct ratio propagates the quadrature core-mass error only because
no uncertainty is assigned to the adopted fragment masses.

The delivered Paper I catalogue is retained as an optional provenance
reference. The working catalogue preserves the decimal-place masses needed to
reproduce the Paper I regional totals, unlike the rounded delivered table.
Selection tests the leaf masks in `my_mask_idx.fits`, using the `_idx` mapping
in `my_dendrogram.dat`. Every source pixel, including its boundary, is sampled
at the PB-map resolution. Fragment masks are not clipped. Cloud masses remain
full-subregion values and are not matched to the primary-beam footprint.

The regional summary is written under `results/core_formation_efficiency/tables/`.
The same directory contains the per-fragment selection and sampled PB values.
Publication PDFs and matching PNG previews are written under
`results/core_formation_efficiency/figures/`.

## Temperature sampling (Paper II Fig. 9)

Stage `temperature_sampling` (moved from the R8d referee diagnostic on
2026-11-20) draws explicit core-temperature priors between Tenv and the
regional upper temperatures 30/40/50 K (Young/Intermediate/Evolved): a
uniform prior and a peak-brightness beta prior, with independent core draws
and a shared regional percentile. Each draw recomputes the optical depth
with the background kept at Tenv; 100 000 draws, fixed seed, about 30 s.
Settings are in `[temperature_sampling]`; products under
`results/core_formation_efficiency/temperature_sampling/`. These are
conditional temperature-prior ranges, not data-constrained uncertainties;
the adopted catalogues keep Tenv.

Stage `efficiency` draws Fig. 9 (`figures/multiscale_mass_efficiency`): the
adopted Tenv hierarchy with counting errors, and the brightness-prior
regional masses with their 16–84 % intervals. It needs the sampling made
from the current regional summary.

```bash
PYTHONPATH=src conda run -n py312 python tasks/core_formation_efficiency/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/core_formation_efficiency/run.py --stage all --overwrite
```

Stages: `inspect`, `summary`, `temperature_sampling`, `efficiency`,
`overlays`, `all` (in that order).

## Change from Paper II (2026-11-20)

- Selection pb: the removed CASA primary beam is replaced by the MIRIAD
  `cont.pb`. With the archived CASA pb the code reproduces the Paper II
  summary exactly; with the MIRIAD pb the same 49 of 65 fragments are kept,
  so the fragment masses are unchanged.
- With the reviewed `core_physical_properties`: core mass Young 41.9 → 43.6,
  Intermediate 255.9 → 261.0, Evolved 247.9 → 248.0 M☉; fragment-to-core
  efficiency 53.4 → 55.6 %, 38.4 → 39.2 %, 18.7 → 18.7 %; core/cloud
  0.38 → 0.40 %, 1.51 → 1.54 %, 3.10 → 3.10 %.
- Temperature sampling: the moved code reproduces the R8d products exactly
  from the Paper II inputs (summary and draws identical). On the new
  catalogue the brightness-prior independent medians are Young 22.0 → 22.8,
  Intermediate 143.0 → 145.6, Evolved 137.8 → 138.1 M☉; the ordering
  fractions change by at most 0.02. The task was renamed from
  `core_formation` (author).
