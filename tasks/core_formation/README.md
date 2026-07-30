# Core formation

This task measures how much cloud mass is contained in Paper I-scale
fragments and in the current compact-core catalogue. It also displays the
compact cores and larger fragment footprints on the continuum maps.

The cloud masses, published fragment masses, fragment counts, dust opacities,
physical scales, and plot settings are declared in `task.toml`. Paper I used
an opacity of 0.9 cm2/g for the fragment masses. The task multiplies those
masses by 0.9/0.5 so that they use the 0.5 cm2/g opacity adopted here. Cloud
masses remain unchanged because their SED fitting used the same opacity law.
Core masses and core counts are calculated from the generated
physical-properties catalogue. The summary reports fragment/cloud,
core/cloud, and direct fragment-to-core mass efficiencies. The uncertainty
beside the direct ratio propagates the quadrature core-mass error only because
no uncertainty is assigned to the adopted fragment masses.

The delivered Paper I catalogue is retained as an optional provenance
reference and is not read by the calculation.

```bash
conda run -n py312 python scripts/core_formation.py --stage inspect
conda run -n py312 python scripts/core_formation.py --stage all
```

To update only one product group, use `--stage summary`, `--stage efficiency`,
or `--stage overlays`. Add `--overwrite` when replacing existing task
products.

The regional summary is written under `results/core_formation/tables/`.
Publication PDFs and matching PNG previews are written under
`results/core_formation/figures/`.
