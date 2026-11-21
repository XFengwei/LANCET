# Modestly zoomed north-up wavelength cards

This version follows the red top/right outline on the preferred larger cards.
The lower-left image boundary is retained, and approximately 8% is removed from
the top and right of every wavelength view. The field is approximately
16.56 × 10.445 arcmin; exact dimensions and centre are in `grid_settings.json`.
`zoom.json` records the exact pixel crop. All nine component maps receive the
same crop and sliced WCS. Cropping does not interpolate the sky data again.

The north-up/east-left orientation, display colour scaling (with ALMA updated as described below), upright
Spitzer brackets and box, horizontal callout arm, small VVV compass heads,
and MeerKAT 21.7 cm label are retained. Compass and 1 pc scale bar are drawn
for the cropped field, and annotation sky anchors are preserved.

The bracket positions follow the subsequent red sketch: the infrared-dark-cloud
bracket is moved downward and left, while the H II bracket is shifted left and
slightly upward and widened. Both labels remain centred below their brackets.
The exact annotation settings are recorded in `annotations.json`.

The local harmonic background cleanup previously accepted for the right-hand
millimetre extension is retained. It is a presentation adjustment to that
background patch, not added observational coverage or a new scientific map.
The original scientific FITS inputs are preserved; obsolete card versions have been removed.

Each image is still 85.6 × 54.0 mm with a 5 mm black border: 95.6 × 64.0 mm total.
PNG and LZW TIFF exports use 1200 dpi, and PDFs have exact physical page size.
Print PDFs at 100% scale. These are five individual views for lens-specific
interlacing, rather than already interlaced artwork.

The full cards and changed panels are visually checked for cropping and label
placement. `verification.json` records image dimensions, DPI, borders, shared
WCS and coverage. Outlined PDF text cannot be assessed by the automatic
text/collision auditors; their reports are retained without claiming a pass.

Rerun from the repository root:

```sh
PYTHONPATH=src conda run --no-capture-output -n py312 python design/multiwavelength_cards/render_cards.py
```

Run this folder's `verify_cards.py` to validate the exports and rebuild the ZIP.

Only this multiwavelength card version is retained. Its finalized grids and
annotation settings are self-contained for rendering. Historical generation
scripts are preserved under `provenance/`; those scripts reference removed
intermediate versions and are retained as a processing record. The runnable
`render_cards.py` uses the retained grids directly.

The outer border is black on all five cards and is exactly 5 mm on each side.
The ALMA 1.3 mm panel uses the single inferno palette with an asinh stretch,
linear width 100 MJy/sr and fixed limits 20 to 1500 MJy/sr. Its label is white
for contrast. There is no separate colourbar. Only its display mapping changes;
the retained scientific/display FITS grid and other wavelength palettes remain
unchanged. `verification.json` records measured image bounds and black-border
pixel checks.

Upper-left wavelength labels are enlarged by 1.5×. The VVV compass is at the
upper right, with two lower-right contact lines on every image: https://xfengwei.github.io/LANCET
and fengweilookuper@gmail.com. Contacts
stay inside the artwork, leaving the 5 mm black cutting margin free of text.
The ALMA intensity floor and softer faint-background stretch are tuned to the
near-black background appearance of current Paper II Figure 6; the inferno
palette is retained. The manuscript and its figures are unchanged.

The contact details are right-aligned at the lower right on all five cards.
The VVV compass is at the upper right. The MeerKAT 1 pc scale bar is moved
to the lower left to avoid overlapping the contacts.

The compass arrows and central hub form one united silhouette with a continuous
black outline. The white hub is reduced to 1.2 pt diameter and overlays the
arrow junction. Shapely is used by the renderer to form the united outline.

## Location (2026-11-21)

Moved from `results/multiwavelength_cards/` to `design/` (author): these are
print and website artwork, not task results. Only the scripts, this README and
the JSON settings are tracked in git; the rendered PDF/PNG/TIFF/ZIP files, the
`grid/` FITS cache and logs stay local (gitignored). `render_cards.py` loads
the strip code from `tasks/multiwavelength_strip/run.py`.
`provenance/` keeps the scripts of earlier card versions as written.

## ALMA card updated (2026-11-21)

The millimetre card now shows the new combination, as Paper II Fig. 1 does:
`grid/mm_i.fits` was rebuilt on the cached card geometry with the
`multiwavelength_strip` mm settings (delivered J-comb `extended36` product,
beyond the 12m+7m grid the `synthesize_continuum` fill36 TP model, display
cross-fade over mosaic pb 0.2–0.5). This replaces the v3 map and its local
harmonic background cleanup described above. Because the new map is about
3 × fainter in extended gas (TP factor 1.0 instead of the v3 ~3.1), the card
stretch changed from asinh 20–1500 MJy/sr to inferno square root 0–1200
MJy/sr, the scale of Fig. 1 and Fig. 3 (`display_settings.json`). The other
four cards are unchanged.
