# LANCET compact logo materials

For a self-contained design handoff and asset guidance, see
[`DESIGN_BRIEF.md`](DESIGN_BRIEF.md).

Approved compact wordmark, with diagonal lancet/filament motif, gold star, and radio telescope.

- `lancet-white-gold-on-black-compact.png`: plain dark-background edition.
- `lancet-black-gold-on-white-compact.png`: plain light-background edition.
- `lancet-white-gold-on-black-compact-halo.png`: white glow and gold star halo on black.
- `lancet-black-gold-on-white-compact-halo.png`: charcoal halo and gold star halo on white.

Generated using the built-in image-generation tool. Halo-edit brief: preserve the approved compact layout, typography, and emblem; add a soft white halo on black or a soft charcoal halo on white, with a warm gold star halo. These are raster PNGs with opaque backgrounds. Superseded loose-lettering editions and the unsuccessful transparent draft were removed from this folder. Generated source history remains outside the repository.

Transparent halo editions (RGBA PNG, with transparent outer backgrounds):
- `lancet-white-gold-compact-transparent-halo.png`
- `lancet-black-gold-compact-transparent-halo.png`

The image tool retained soft halos in these transparent exports; they are not plain edge-only cutouts. The black logo should be viewed over a light background.

## Antenna correction (2026-11-10)

All six original PNG filenames include the approved symmetric, shortened upper
antenna assembly. The original dish bowl, pedestal, filament, star, wordmark,
background treatment, and transparency are retained. Test editions, construction
previews, and obsolete iteration scripts have been removed.

Geometry uses the dish-opening midpoint (1078, 552) and the perpendicular to the
rim chord from (1004, 591) to (1152, 513). The two support attachment points are
reflections across this axis. Support height from their baseline is 85 px,
compressed from 102 px by a factor of 1.2. The receiver radius is 12 px and the
support width is 6 px. The same geometry is applied to every style, with local
colour/alpha compositing and halo blending at the original curved rim.

Validation confirmed the six original filenames, 1536 × 1024 dimensions, RGB or
RGBA mode as appropriate, exact geometric reflection, and unchanged pixels
outside the upper antenna and joining region. Website asset copies are separate
and were not part of this six-file update.
