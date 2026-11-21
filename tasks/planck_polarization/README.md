# Planck polarization

This task reads the delivered Planck 353 GHz Stokes I, Q, and U maps and
creates signal-to-noise masks, modified-asymptotic debiased polarized
intensity, polarization fraction, polarization angle and error, and the
plane-of-sky magnetic-field orientation. It also creates a vector map of the
central 30′ × 30′. All products are written under
`results/planck_polarization/`; the delivered Planck maps remain read-only.

The rms of the Stokes I and Q/U maps is set in `[noise]`, the signal-to-noise
thresholds in `[thresholds]`. The task gives field orientations; it does not
estimate the field strength.

```bash
PYTHONPATH=src conda run -n py312 python tasks/planck_polarization/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/planck_polarization/run.py --stage all --overwrite
```

Use `--stage build` or `--stage plot` for individual product groups. Existing
products are kept unless `--overwrite` is given.

## Review (2026-11-21)

Renamed from `magnetic_fields` (author): the task makes Planck polarization
products and only the field orientation, no field strength. Rerun: all nine
maps and `polarization_summary.json` are identical to the previous products;
results moved from `results/magnetic_fields/` to `results/planck_polarization/`.
The code module `lancet.magnetic_fields` and the product file names
(`magnetic_field_angle.fits`, `planck_magnetic_field_diagnostic.*`) are kept.

The status label box on the figure was removed (author); the top colour-bar
ticks (0.40, 0.45) still overlap.

Figure (author, 2026-11-21): the diagnostic now shows a 30′ × 30′ box centred
on the field (`view_*` in `[plot]`) with a vector every half beam. The beam is
4.82′ (Planck 2018 353 GHz effective FWHM, confirmed by the author; the
delivered cut-outs carry no beam keyword), which on the 1′ grid gives a vector every
2 pixels (2′, 0.41 beam). The maps themselves are not cropped.
