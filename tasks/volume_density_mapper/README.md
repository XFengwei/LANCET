# Volume-density mapper test

This task applies Guang-Xing Li's `volume-density-mapper` to
`surface_density_allscale.fits` without changing the input. The pilot is the
exact 1500 by 1500 native-pixel cutout centred at ICRS
`14:44:58.79 -59:48:01.94`.

The map has a 0.1-arcsec pixel grid but not a uniform beam. Pixels inside the
finite ALMA primary-beam-cut footprint are assigned 0.65 arcsec resolution;
the remainder are assigned 18 arcsec. Inferred widths are never allowed below
the local beam. The output quality bit mask records:

- bit 0 (1): inside the high-resolution ALMA footprint;
- bit 1 (2): upstream characteristic width is finite and positive;
- bit 2 (4): the local resolution floor was applied;
- bit 3 (8): within one 18-arcsec beam of the resolution transition.

Install the mapper in the active environment or pass a source checkout:

```bash
python -m pip install volume-density-mapper
conda run -n py312 python scripts/volume_density_mapper.py inspect
conda run -n py312 python scripts/volume_density_mapper.py pilot
conda run -n py312 python scripts/volume_density_mapper.py full
conda run -n py312 python scripts/volume_density_mapper.py mask
```

For a checkout that has not been installed, append
`--mapper-source /path/to/volume-density-mapper/src` to the `pilot` or `full`
command.

The upstream 3-D function is tested only at 250 cubed in the pilot. A native
1500-cubed pilot would occupy 13.5 GB, and a native full-map cube would occupy
about 1.0 TB before FITS overhead. More importantly, the current upstream
implementation omits the final constrained-diffusion residual from the cube.
The task therefore does not promote a full 3-D reconstruction unless column
closure is fixed and revalidated.

The `mask` stage applies the numerical 12m+7m primary-beam response map with
the same response cutoff, 0.20, used by `prepare_continuum_images`. It writes
the recommended science products under `pilot_pbmasked/` and
`full_pbmasked/`. Every floating-point value outside the retained primary-beam
region is NaN, and the quality mask is zero there. The earlier unmasked
products remain available only for diagnosing the resolution boundary.

Measured runtimes, numerical checks, pilot-to-full differences, and the
interpretation boundary are recorded in
[`docs/volume_density_mapper_audit.md`](../../docs/volume_density_mapper_audit.md).
