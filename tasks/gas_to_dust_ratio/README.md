# Gas-to-dust ratio

Gas-to-dust mass ratio at the filament's Galactocentric radius from the
Giannetti et al. (2017, A&A 606, L12) Eq. 2 relation,
log10 R_gd = 0.087 R_gal/kpc + 1.44, evaluated by `lancet.gas_to_dust`.
This is the R_gd quoted in Paper II (Sect. 3, surface density and core
masses).

## Input

`[galaxy]` in `task.toml`: R_gal = 6.5 ± 0.2 kpc (Watkins et al. 2019). No
data files.

## Output

`results/gas_to_dust_ratio/gas_to_dust_ratio.json`: the ratio and its
statistical (coefficient and radius errors), systematic (paired coefficient
branches) and intrinsic-scatter (0.2 dex) intervals, kept separate.

```bash
PYTHONPATH=src conda run -n py312 python tasks/gas_to_dust_ratio/run.py --stage all --overwrite
```

## Result (2026-11-21)

R_gd = 101.3; statistical −12.5 / +14.3; systematic −30.8 / +11.7;
intrinsic scatter −37.4 / +59.2. Paper II quotes 101 [+12 −30] ± 13.

## Readers

None read this file yet. `core_physical_properties` uses
`gas_to_dust_ratio = 100.0 ± 24` and `column_density` the total-gas opacity
0.005 cm² g⁻¹ (= 0.5 / 100), i.e. R_gd = 100 rounded; `synthesize_continuum`
also uses 100.
