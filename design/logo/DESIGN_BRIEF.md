# LANCET logo design handoff

This brief carries the established logo context into a new design discussion. It records the current direction and assets; it does not prescribe a redesign.

## Identity

- **Name:** LANCET.
- **Expansion:** Linear Filament and Nested Cluster Evolution Tomography.
- **Scientific context:** LANCET studies the Galactic filament G316.8 and the evolution and assembly of stars and clusters across nested scales. The logo represents the research programme; it should not imply that G316.8 is the only future target unless that scope is intentionally chosen.
- The wordmark spells **LANCET** in uppercase. The expansion is used in the manuscript and website, but is not part of the compact logo artwork.

## Established visual direction

The documented compact mark combines:

- a diagonal, sweeping **lancet / filament motif**;
- a **gold four-point star** near its tip;
- a **radio-telescope dish**;
- the uppercase **LANCET** wordmark.

The motif and telescope connect filamentary structure with astronomical observation; the star evokes star formation. These are readings of the existing artwork, not extra constraints on future exploration. The README identifies this compact layout as the approved direction. Keep its recognizable composition and symbols when preparing variants unless the user asks to revisit the concept.

The existing palette pairs white or charcoal lettering and symbols with warm gold for the star. Both dark- and light-background treatments exist. Halo treatments also exist, but should be chosen for a specific use rather than assumed to be the default.

## Current assets

All files below are 1536 × 1024 PNGs. All six include the mathematically symmetric antenna supports, with axial height reduced by a factor of 1.2, applied on 2026-11-10. The original dish and remaining artwork are preserved.

| File | Background and treatment | Notes |
|---|---|---|
| `lancet-white-gold-on-black-compact.png` | Black, white mark, gold star | Plain dark-background version |
| `lancet-black-gold-on-white-compact.png` | White, charcoal mark, gold star | Plain light-background version |
| `lancet-white-gold-on-black-compact-halo.png` | Black, white/ gold mark with glow | Halo version |
| `lancet-black-gold-on-white-compact-halo.png` | White, charcoal/ gold mark with glow | Halo version |
| `lancet-white-gold-compact-transparent-halo.png` | Transparent outer canvas; white/gold mark with glow | RGBA; glow remains visible around the mark |
| `lancet-black-gold-compact-transparent-halo.png` | Transparent outer canvas; charcoal/gold mark with glow | RGBA; inspect over a light background |

The transparent files are **not** clean, edge-only cutouts: the soft halo remains. These are generated raster images, not editable vector masters. Earlier loose-lettering versions and an unsuccessful transparent draft were removed; generation history is not stored in this folder.

`website/assets/lancet-logo.png` is the website logo source. Website CSS crops and scales that image for the header/footer, so its appearance there is not a separate logo lockup. Check small-size legibility and background behavior in the intended placement before selecting an asset for another use.

## Practical guidance for the next design chat

- Start from the existing compact artwork as the visual reference. Preserve the wordmark spelling, the lancet/filament sweep, gold star, and radio dish unless explicitly asked to explore a different concept.
- Choose the dark, light, or transparent treatment based on where it will appear; a halo changes the appearance of the transparent variants.
- Because the available sources are PNGs, request a clean vector master (for example SVG/PDF) only if the intended use needs scalable print or precise layout control. Do not treat an automatic trace as an approved master without checking its curves, lettering, and fine telescope details.
- Useful open decisions for future use: primary version for general branding, minimum readable size, whether the expansion should ever accompany the wordmark, and whether a genuinely halo-free transparent lockup or vector master is needed.

## Source notes

The asset list and status above are based on this folder's existing `README.md` and files. The project expansion and scientific context are reflected in the Paper II title and LANCET website. No color values, font family, trademark rules, or formal brand manual are recorded here; do not infer exact specifications from the raster artwork.
