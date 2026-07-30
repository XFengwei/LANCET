# Magnetic-field analysis — ongoing

Status: **ongoing idea; not paper-ready**.

This task reads the delivered Planck 353 GHz Stokes I, Q, and U maps and
creates signal-to-noise masks, modified-asymptotic debiased polarized
intensity, polarization fraction, polarization angle and error, and the
plane-of-sky magnetic-field orientation. It also creates a diagnostic vector
map. All derivatives are written under `results/magnetic_fields/`; the
delivered Planck maps remain read-only.

The configured I and Q/U RMS values are exploratory values inherited from the
current analysis. Before this becomes a publication workflow, the project
must establish a reproducible Planck noise estimate, confirm the Q/U and
position-angle convention on the delivered Galactic grid, decide how beam
correlation and covariance are handled, and select and validate a
magnetic-field-strength method. The present task estimates orientation only;
it does not calculate field strength.

```bash
conda run -n py312 python scripts/magnetic_fields.py --stage inspect
conda run -n py312 python scripts/magnetic_fields.py --stage all
```

Use `--stage build` or `--stage plot` for individual product groups. Add
`--overwrite` when intentionally replacing an earlier ongoing iteration.
