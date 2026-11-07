# LANCET website

Four static HTML pages, shared CSS, and optimized project imagery. Preview from the repository root with `python3 -m http.server 8765 --directory website`.

The homepage animation cycles through registered VVV, Spitzer + Herschel 70 μm, Herschel 250 μm, and MeerKAT views. The hero frames share the same sky grid and exclude the 1.3 mm image. Regenerate them with `MPLCONFIGDIR=/tmp/lancet-mpl /Users/fexu/miniforge/envs/py312/bin/python website/tools/build_hero.py` in the project environment.

`website/tools/build_footprints.py` regenerates the WCS-aligned pointing overlays from local FITS maps and cached measurement-set FIELD metadata. The six narrow spectral windows listed on the Observations page are the transitions targeted in Paper I; their window centres are not rest frequencies. This list does not claim to inventory every serendipitous line.

To publish, select **GitHub Actions** under repository Settings → Pages. The workflow deploys the `website/` directory from pushes to `main` or a manual workflow run. Expected address: https://xfengwei.github.io/LANCET/ .

Project imagery and diagnostics are optimized copies of materials in `presentation/` and `results/`. The team page uses the Paper II manuscript collaborator list, with PI and data-reduction roles and affiliations supplied by the project lead. Publication details and research projects follow the current project information.
