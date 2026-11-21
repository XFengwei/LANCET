# LANCET data

The `data/` directory contains the inputs of the LANCET project. It has two strictly
different zones:

- `reduction/` is a mutable workspace for active ALMA and total-power data
  reduction.
- `delivered/` contains the stable products that project analyses are allowed
  to read.

## Project rule

Analysis tasks and scientific scripts:

1. should use `data/delivered/` for stable scientific inputs;
2. may read `data/reduction/` read-only when explicitly documenting or
   visualizing an intermediate reduction step;
3. must not modify reduction data as part of an analysis run; and
4. must write all intermediate products, tables, FITS files, figures, and logs
   outside `data/`.

A reduction input used by an analysis workflow must be declared as a
`read_only_intermediate` dependency. Such an input is allowed for reduction
diagnostics and method illustrations, but should not silently become the source
of a stable scientific result. Stable downstream dependencies should still be
promoted to `data/delivered/`.

Reduction-specific CASA, MIRIAD, or export scripts may operate inside
`data/reduction/`. Once a reduction product is accepted for analysis, it should
be copied or promoted into `data/delivered/` and then treated as immutable.

Some current files predate this rule. Their locations are documented in
`STATUS.md`; they will be corrected when the corresponding workflows are
refactored.

## Tracking files

- `MANIFEST.tsv` lists every file under `delivered/` and the top two levels
  of `reduction/` and `extraction/`: path, size, origin (`observed`,
  `archive`, `collaborator`, `reduced`, `derived`), whether it can be
  reproduced (`download`, `reduction`, `task:<name>`, `no`), the tasks that
  read it, and a note. It replaces the old path ledgers
  (`path_aliases.yaml`, `path_ledger*.tsv`), which were removed. Regenerate it
  from the repository root with `python docs/make_data_manifest.py`.
- `STATUS.md` summarizes the current audit and known transition issues.
- Git tracks `README.md` and `MANIFEST.tsv` only, so a fresh clone states
  what data it needs.

## Delivered ancillary catalogues

`delivered/starcluster/` (the open-cluster orbit catalogue) was removed on
2026-11-21 together with the `star_cluster` script (author).

## Delivered SED-fitting products

`delivered/sedfitter/` began as the read-only migration of the former
root-level `sedfitter/` tree. Its migration-time paths remain recorded even
when a verified task-output duplicate is later removed. The directory mixes
external SED-fitting results with historical task products; its README and
the combine audit classify those products and point active workflows to their
`results/` replacements.

## Delivered archive products

Historical paths under `Archive/` now resolve to
`delivered/archive/`. These are read-only inputs even when their filenames
indicate that they were previously cropped, reprojected, convolved, or had
their headers repaired. The free-free task audit identifies the radio
derivatives used by the Paper II free-free workflow in
[`docs/construct_freefree_model_audit.md`](../docs/construct_freefree_model_audit.md).
