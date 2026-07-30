# N-PDF analysis — ongoing

Status: **ongoing idea; not paper-ready**.

This task makes three 240-arcsec cutouts from the retained
`NH2_6arcsec+18arcsec+36arcsec.fits` map and calculates regional
column-density PDFs. It subtracts the configured fixed background, discards
non-positive excess values before taking a logarithm, and normalizes all
regions by the shared mean of their positive background-subtracted pixels.
Cutouts, histogram tables, assumptions, and figures are written under
`results/n_pdf/`; the delivered column-density map remains read-only.

The current map replaces the removed historical `NH2_allscale.fits` filename,
but its suitability as the final N-PDF input has not been established. Before
publication use, the project must decide which angular-scale combination is
scientifically appropriate, replace or justify the fixed background,
establish masks and completeness, account for correlated pixels and
resolution, and implement and validate any lognormal or power-law-tail
fitting. Changing these choices in `task.toml` creates a reproducible new
iteration without changing the script.

```bash
conda run -n py312 python scripts/n_pdf.py --stage inspect
conda run -n py312 python scripts/n_pdf.py --stage all
```

Use `--stage build` or `--stage plot` for individual product groups. Add
`--overwrite` when intentionally replacing an earlier ongoing iteration.
