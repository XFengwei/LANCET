"""Copy the website version of Paper II Fig. 2 into assets/coverage.jpg.

The figure is drawn by the combine_continuum task (stage figure2_dark: the same
data and layout as Fig. 2, on the page background with dark-mode colour steps).
This script only converts the verified task PNG to JPEG; it does not replot.
"""
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'results/combine_continuum/diagnostics/figures/figure2_combination_diagnostics_dark.png'
OUTPUT = ROOT/'website/assets/coverage.jpg'

image = Image.open(SOURCE).convert('RGB')
image.save(OUTPUT, quality=90, optimize=True)
print(f'{OUTPUT.name}: {image.size[0]} x {image.size[1]} from {SOURCE.relative_to(ROOT)}')
