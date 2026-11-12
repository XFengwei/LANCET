"""Export web brand assets from the approved dark-background logo.

Outputs (all in website/assets/):
- lancet-logo-web.webp: header/footer logo with a true alpha channel. The black
  canvas is converted to transparency ("unscreened"), so the mark composites
  over the starfield exactly as it looks over black, without a dark box.
- apple-touch-icon.png: 180 px gold four-point star on the site background.
- og-image.jpg: 1200 x 630 social preview from the near-infrared hero frame.

The logo source in design/logo/ is read only and never modified.
"""
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "website/assets"
SOURCE = ROOT / "design/logo/lancet-white-gold-on-black-compact.png"
BACKGROUND = (8, 11, 16)  # --bg in style.css
GOLD = (228, 188, 135)  # --gold in style.css
BLACK_POINT = 3 / 255  # generated canvas carries faint noise just above zero
LOGO_BOX = (307, 77, 1229, 891)  # visible crop used by .logo-brand in the CSS


def unscreen(image: Image.Image) -> Image.Image:
    """Return RGBA whose normal composite over black equals the RGB input."""
    rgb = np.asarray(image.convert("RGB"), dtype=np.float64) / 255
    peak = rgb.max(axis=2)
    alpha = np.clip((peak - BLACK_POINT) / (1 - BLACK_POINT), 0, 1)
    colour = np.clip(rgb / np.maximum(alpha, 1e-6)[..., None], 0, 1)
    rgba = np.dstack([colour, alpha])
    return Image.fromarray((rgba * 255 + 0.5).astype("uint8"), "RGBA")


def star(size: int, scale: float, centre: tuple[float, float]) -> list[tuple[float, float]]:
    """Concave four-point star matching favicon.svg, sampled as a polygon."""
    def bezier(p0, p1, p2, p3, n=24):
        t = np.linspace(0, 1, n)[:, None]
        return ((1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1
                + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3)
    # One quarter in a 32-unit box (top point to right point), rotated four times.
    quarter = bezier(np.array([16, 0]), np.array([17, 10]), np.array([22, 15]), np.array([32, 16]))
    points = []
    for k in range(4):
        angle = k * np.pi / 2
        rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        points.extend((quarter - 16) @ rot.T)
    cx, cy = centre
    return [(cx + x * scale, cy + y * scale) for x, y in points]


def main() -> None:
    source = Image.open(SOURCE)

    logo = unscreen(source.resize((600, 400), Image.LANCZOS))
    logo.save(ASSETS / "lancet-logo-web.webp", quality=92, alpha_quality=100, method=6)

    icon = Image.new("RGB", (720, 720), BACKGROUND)
    ImageDraw.Draw(icon).polygon(star(720, 13.5, (360, 360)), fill=GOLD)
    icon.resize((180, 180), Image.LANCZOS).save(ASSETS / "apple-touch-icon.png", optimize=True)

    sky = Image.open(ASSETS / "hero-nir.webp").convert("RGB")
    sky = sky.resize((1200, round(sky.height * 1200 / sky.width)), Image.LANCZOS)
    top = round((sky.height - 630) * 0.65)  # matches the hero's object-position
    card = np.asarray(sky.crop((0, top, 1200, top + 630)), dtype=np.float64)
    fade = np.clip((0.78 - np.linspace(0, 1, 1200)) / 0.45, 0, 1) ** 1.2 * 0.82
    card = card * (1 - fade[None, :, None]) + np.array(BACKGROUND) * fade[None, :, None]
    card = Image.fromarray(card.astype("uint8"))
    mark = unscreen(source.crop(LOGO_BOX))
    mark = mark.resize((round(mark.width * 400 / mark.height), 400), Image.LANCZOS)
    card.paste(mark, (90, 115), mark)
    card.save(ASSETS / "og-image.jpg", quality=86, optimize=True, progressive=True)

    print("wrote lancet-logo-web.webp, apple-touch-icon.png, og-image.jpg")


if __name__ == "__main__":
    main()
