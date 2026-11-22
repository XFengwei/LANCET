# LANCET website

Four static HTML pages (plus `data.html`, a placeholder for the coming data release), shared CSS, and optimized project imagery. Preview from the repository root with `python3 -m http.server 8765 --directory website`.

The homepage animation cycles through registered VVV, Spitzer + Herschel 70 μm, Herschel 250 μm, and MeerKAT views. The hero frames share the same sky grid and exclude the 1.3 mm image. Regenerate them with `MPLCONFIGDIR=/tmp/lancet-mpl conda run -n py312 python website/tools/build_hero.py` in the project environment. The script writes the WebP frames directly (1800 px, quality 82).

`website/tools/build_brand.py` exports the header/footer logo with a true alpha channel (the black canvas of `design/logo/lancet-white-gold-on-black-compact.png` converted to transparency, so no dark box shows over the hero), the 180 px touch icon, and the 1200 × 630 social preview `assets/og-image.jpg`. `assets/favicon.svg` is the gold four-point star, drawn by hand. The social metadata uses absolute URLs at the expected GitHub Pages address.

`website/tools/build_footprints.py` regenerates the WCS-aligned pointing overlays from local FITS maps and cached measurement-set FIELD metadata. The six narrow spectral windows listed on the Observations page are the transitions targeted in Paper I; their window centres are not rest frequencies. This list does not claim to inventory every serendipitous line.

`website/tools/build_array_gallery.py` renders the five continuum products of the Observations page (ACA 7 m, ALMA 12 m, 12 m + 7 m, total-power model, fully combined) on the multiwavelength-strip field with the strip's millimetre display. `website/tools/build_coverage.py` copies the website version of Paper II Fig. 2, drawn by the `combine_continuum` stage `figure2_dark`, to `assets/coverage.jpg`. `website/tools/build_gallery.py` renders the five wavelength frames from the `multiwavelength_strip` grids. Run them in the project environment after the tasks they read have run.

To publish, select **GitHub Actions** under repository Settings → Pages. The workflow deploys the `website/` directory from pushes to `main` or a manual workflow run. Expected address: https://xfengwei.github.io/LANCET/ .

Project imagery and diagnostics are optimized copies of materials in `presentation/` and `results/`. The team page uses the Paper II manuscript collaborator list, with PI and data-reduction roles and affiliations supplied by the project lead. Publication details and research projects follow the current project information.

## Logo animation

Open `logo-animation.html` through the preview server for the standalone
12-second motion study. It uses local SVG, CSS, and JavaScript without runtime
libraries or network dependencies.

- 0–3.4 s: an irregular cloud contracts into one curved filamentary spine.
- 3.4–6 s: the broad knife and telescope appear together. The white middle slash
  is progressively exposed behind the blade along its track as the filament opens.
- 5.8–8 s: the primary star appears, followed by smaller cluster members.
- 7.5–9.2 s: the cluster dissolves into the emblem and the ribbon morph completes.
- 9.2–12 s: the exact non-halo emblem settles in and LANCET fades in last.

The former 10 s, 14 s, and 18 s landmarks fall at 7.5 s, 9 s, and 12 s.
Monotone cubic timing blends speed through these landmarks. Quintic easing
softens phase starts and finishes; overlapping morph phases avoid an intermediate
pause. The knife and incision boundary interpolate between contour samples. The slash is part of the cutting action,
with no separate delayed appearance.

`assets/lancet-animation-final.png` is a byte-identical copy of
`design/logo/lancet-white-gold-on-black-compact.png`, including the corrected,
shortened symmetric antenna. The logo source editions are not modified by the
animation. SVG target contours are embedded in the page for smooth interpolation;
the original PNG supplies the exact final appearance.

The timing is deterministic: seeking backward or replaying reconstructs the same
scene. Playback runs once and holds, with pause, replay and a seek slider.
Reduced-motion preferences show the completed logo immediately; explicitly
starting playback remains possible. Switching away pauses playback. The star
sequence is the project's visual narrative, not a general chronological model
for every stellar cluster. The preview remains separate from the homepage.

## Design essentials

The site is dark: page `#080b10`, text `#ecece8`, muted text and axes
`#a6adb6`, accent gold `#e4bc87`, hairlines `#ffffff20` (`style.css`). Every
image sits on that page without a visible box.

**Science figures are replotted for the page, not recoloured.** A paper figure
on white is never inverted, filtered or recoloured as a picture: inversion
changes every hue, and the numbers would then come from somewhere other than
the paper. The task that owns the figure draws a second version from the same
data and layout with a dark theme (for Fig. 2: `combine_continuum` stage
`figure2_dark`, colours in `lancet.combination.coverage_figure.THEMES`). The
tools here only copy or convert task outputs, so the website figure and the
paper figure cannot drift apart. Rules for a dark theme:

- Figure, axes and inset backgrounds are the page colour `#080b10`; label boxes
  and legend frames use it too, never white.
- Text is `#ecece8`; spines and ticks `#a6adb6`. The paper's blacks and greys
  become these neutrals or darker greys; a neutral never takes a series hue.
- Each array, curve or category keeps its hue from the paper figure, so a reader
  can move between the two. Only its lightness changes: dark-mode colours sit in
  OKLCH L 0.48–0.67 with at least 3:1 contrast on `#080b10`.
- Check categorical colours with a palette validator (OKLab ΔE ×100, normal
  vision and protan/deutan/tritan simulation). Series that overlap in a scatter
  must pass all pairs (Fig. 2 arrays: worst colour-blind ΔE 19.6). Where a pair
  falls to ΔE 6–8, a second cue is required: markers, line style or a direct
  label (Fig. 2 image curves: dashed 12 m + 7 m, solid combined).
- Faint elements that the white version showed through contrast (the TP
  response discs in the Fig. 2 inset) get lighter fills and edges on dark, so
  they read at the same strength.
- The light (paper) version must stay pixel-identical when the theme code
  changes.

**Sky maps share one display.** Comparable views use one field, rotation, colour
map and stretch: for 1.3 mm continuum the Paper II Fig. 3 view (inferno,
square-root stretch, 0–1200 MJy/sr) on the multiwavelength-strip grid. No
display factor brightens one product to look like another; blank sky is the
page colour. Labels are white with a black outline, with a compass and a 1 pc
bar.

**Pages reveal their sections while scrolling** (all four pages). Each
block below the hero or page title carries `class="reveal"`;
`reveal.js` fades it in and lifts it 40 px when its top passes 80 % of the
window height, and it then stays. Revealing a block also reveals every block
above it, so a fast scroll or an anchor jump never leaves a hidden section
behind. A flag set in `<head>` (`html.js-reveal`) hides the blocks only when the
script will run: without JavaScript, with reduced motion, or in print, the
whole page shows at once. A new page needs the flag script in `<head>`,
`class="reveal"` on its sections below the opening and `reveal.js`.

**Galleries rotate, and readers control them.** Image galleries and the
telescope cards change every 5 s and hold 8 s after a click. Galleries of links
(the telescope rows) pause while the pointer or keyboard focus is on them, and
are stepped with ‹ and › arrows; hidden rows are not clickable or focusable.
On narrow screens and with reduced motion the telescope rows are listed
instead. Units in upper-case labels are wrapped in `<span class="unit">` so
they keep their case (μm, mm, GHz).

**Every image has a source.** Each asset is regenerated by a script in
`tools/` from task results or delivered data, never edited by hand. When an
asset changes, bump its `?v=` query in the HTML so browsers reload it.

