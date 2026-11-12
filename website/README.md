# LANCET website

Four static HTML pages, shared CSS, and optimized project imagery. Preview from the repository root with `python3 -m http.server 8765 --directory website`.

The homepage animation cycles through registered VVV, Spitzer + Herschel 70 μm, Herschel 250 μm, and MeerKAT views. The hero frames share the same sky grid and exclude the 1.3 mm image. Regenerate them with `MPLCONFIGDIR=/tmp/lancet-mpl /Users/fexu/miniforge/envs/py312/bin/python website/tools/build_hero.py` in the project environment. The script writes the WebP frames directly (1800 px, quality 82).

`website/tools/build_brand.py` exports the header/footer logo with a true alpha channel (the black canvas of `design/logo/lancet-white-gold-on-black-compact.png` converted to transparency, so no dark box shows over the hero), the 180 px touch icon, and the 1200 × 630 social preview `assets/og-image.jpg`. `assets/favicon.svg` is the gold four-point star, drawn by hand. The social metadata uses absolute URLs at the expected GitHub Pages address.

`website/tools/build_footprints.py` regenerates the WCS-aligned pointing overlays from local FITS maps and cached measurement-set FIELD metadata. The six narrow spectral windows listed on the Observations page are the transitions targeted in Paper I; their window centres are not rest frequencies. This list does not claim to inventory every serendipitous line.

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
