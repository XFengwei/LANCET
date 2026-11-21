# Background-emission effects

This supporting task quantifies how foreground and background dust emission
changes the inferred optical depth of a compact core. It has no file inputs;
all physical and plotting assumptions are declared in `task.toml`.

For an observed environmental intensity `I_bg`, dust temperature `T`, and
Planck function `B_nu(T)`, the total optical depth of a symmetric
foreground-plus-background layer is

```text
tau_env = -ln(1 - I_bg / B_nu)
```

The environmental emission in front of the core is

```text
I_fg = I_bg / (sqrt(1 - I_bg / B_nu) + 1)
```

and the corrected core optical depth is

```text
tau_core = -ln[
    1 - I_core / (B_nu - I_fg) / sqrt(1 - I_bg / B_nu)
]
```

The default calculation uses 223.8 GHz, 20 K, a 0.71 by 0.61 arcsec Gaussian
beam, and core brightnesses of 1, 10, and 100 mJy/beam. The figure shows the
core optical depth relative to the zero-background result, along with the
total environmental optical depth.

```bash
PYTHONPATH=src conda run -n py312 python tasks/background_effects/run.py --stage inspect
PYTHONPATH=src conda run -n py312 python tasks/background_effects/run.py --stage all
```

The machine-readable table and JSON summary are written under
`results/background_effects/tables/`. PDF and PNG figures are written under
`results/background_effects/figures/`. Add `--overwrite` only when replacing
existing products.

## Review (2026-11-21)

Rerun: table and summary byte-identical; the model has no data inputs, so
the new combination does not change it. Where the cores sit on the curve:
the reviewed `source_catalogue` core background `Ibkg2` has median 148 and
maximum 758 MJy/sr (Paper II: 309 and 961). The correction τ_core(I_bg)/τ_core(0)
is 1.006–1.008 at the new median and
1.034–1.043 at the new maximum (1–100 mJy/beam), against
1.013–1.017 and 1.043–1.056 in Paper II.
