"""Seamless mosaicking and colour stretching of near-infrared survey tiles."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import warnings

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import FITSFixedWarning, WCS
import numpy as np
from reproject import reproject_interp
from scipy import ndimage


def load_tile(
    path: Path, reference_zero_point: float, edge_trim: int = 0
) -> tuple[np.ndarray, WCS]:
    """Read one tile, remove its sky pedestal and rescale to a common zero point.

    VVV tiles carry a padded, under-exposed border; *edge_trim* pixels on each
    side are blanked so they never enter the mosaic.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        with fits.open(path, memmap=False) as hdul:
            header = hdul[0].header
            data = hdul[0].data.astype(np.float32)
            wcs = WCS(header)
    _, sky, _ = sigma_clipped_stats(data[::4, ::4], sigma=3.0, maxiters=5)
    scale = 10.0 ** (0.4 * (reference_zero_point - float(header["MAGZPT"])))
    data = (data - np.float32(sky)) * np.float32(scale)
    if edge_trim > 0:
        data[:edge_trim, :] = np.nan
        data[-edge_trim:, :] = np.nan
        data[:, :edge_trim] = np.nan
        data[:, -edge_trim:] = np.nan
    return data, wcs


def block_median(image: np.ndarray, block: int) -> np.ndarray:
    """Median of each block x block cell; NaN where a cell is not fully covered."""
    ny, nx = (-(-size // block) for size in image.shape)
    padded = np.full((ny * block, nx * block), np.nan, dtype=np.float32)
    padded[: image.shape[0], : image.shape[1]] = image
    cells = padded.reshape(ny, block, nx, block).transpose(0, 2, 1, 3)
    return np.median(cells.reshape(ny, nx, -1), axis=-1)


def block_any(mask: np.ndarray, block: int) -> np.ndarray:
    """True for each block x block cell containing at least one True pixel."""
    ny, nx = (-(-size // block) for size in mask.shape)
    padded = np.zeros((ny * block, nx * block), dtype=bool)
    padded[: mask.shape[0], : mask.shape[1]] = mask
    return padded.reshape(ny, block, nx, block).any(axis=(1, 3))


def solve_offset_surfaces(
    backgrounds: Sequence[np.ndarray],
    domains: Sequence[np.ndarray],
    smoothness: float,
    decay_cells: float,
) -> list[np.ndarray]:
    """Find smooth per-tile offsets that make overlapping tiles agree.

    For every coarse cell covered by two tiles, the offset difference should
    equal their background difference. Neighbouring cells of one tile are tied
    together by *smoothness*, and a weak damping term pulls every offset to
    zero, so corrections decay into tile interiors over ~*decay_cells*.
    Only tile differences enter, so sky structure common to the tiles (the
    nebula) is left untouched.
    """
    from scipy import sparse
    from scipy.sparse.linalg import spsolve

    index = []
    count = 0
    for domain in domains:
        ids = np.full(domain.shape, -1, dtype=np.int64)
        ids[domain] = np.arange(count, count + domain.sum())
        count += int(domain.sum())
        index.append(ids)

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    rhs: list[np.ndarray] = []
    n_rows = 0

    def add(columns: Sequence[np.ndarray], weights: Sequence[float], target: np.ndarray):
        nonlocal n_rows
        n = target.size
        row = np.arange(n_rows, n_rows + n)
        for column, weight in zip(columns, weights):
            rows.append(row)
            cols.append(column)
            vals.append(np.full(n, weight))
        rhs.append(target)
        n_rows += n

    for i in range(len(backgrounds)):
        for j in range(i + 1, len(backgrounds)):
            difference = backgrounds[i] - backgrounds[j]
            both = np.isfinite(difference) & domains[i] & domains[j]
            if both.any():
                add((index[i][both], index[j][both]), (1.0, -1.0), difference[both])

    root_smooth = np.sqrt(smoothness)
    root_damp = np.sqrt(smoothness) / decay_cells
    for ids in index:
        for a, b in ((ids[:, :-1], ids[:, 1:]), (ids[:-1, :], ids[1:, :])):
            pair = (a >= 0) & (b >= 0)
            add((a[pair], b[pair]), (root_smooth, -root_smooth), np.zeros(pair.sum()))
        valid = ids[ids >= 0]
        add((valid,), (root_damp,), np.zeros(valid.size))

    design = sparse.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_rows, count),
    )
    solution = spsolve(
        (design.T @ design).tocsc(), design.T @ np.concatenate(rhs)
    )

    surfaces = []
    for ids in index:
        surface = np.full(ids.shape, np.nan)
        surface[ids >= 0] = solution[ids[ids >= 0]]
        surfaces.append(surface)
    return surfaces


def upsample_surface(surface: np.ndarray, block: int, shape: tuple[int, int]) -> np.ndarray:
    """Bilinearly expand a coarse surface to full resolution, extending its edges."""
    missing = ~np.isfinite(surface)
    if missing.any():
        nearest = ndimage.distance_transform_edt(
            missing, return_distances=False, return_indices=True
        )
        surface = surface[tuple(nearest)]
    full = ndimage.zoom(surface, block, order=1, mode="nearest", grid_mode=True)
    return full[: shape[0], : shape[1]].astype(np.float32)


def feather_weight(footprint: np.ndarray, width: float) -> np.ndarray:
    """Weight rising smoothly from 0 at the footprint edge to 1 at *width* pixels."""
    padded = np.pad(footprint, 1, constant_values=False)
    distance = ndimage.distance_transform_edt(padded)[1:-1, 1:-1]
    ramp = np.clip(distance / width, 0.0, 1.0)
    return (ramp * ramp * (3.0 - 2.0 * ramp)).astype(np.float32)



def coadd_tiles(
    tiles: Sequence[tuple[np.ndarray, WCS]],
    wcs_out: WCS,
    shape_out: tuple[int, int],
    *,
    block: int,
    smoothness: float,
    decay_pixels: float,
    feather_pixels: float,
) -> np.ndarray:
    """Reproject tiles onto one grid, flatten their seams and feather them."""
    images = []
    footprints = []
    for tile in tiles:
        image, footprint = reproject_interp(tile, wcs_out, shape_out=shape_out)
        footprint = (footprint > 0.99) & np.isfinite(image)
        if footprint.any():
            images.append(np.where(footprint, image, np.nan).astype(np.float32))
            footprints.append(footprint)

    backgrounds = [block_median(image, block) for image in images]
    # A tile's offset surface covers every cell it touches, but only fully
    # covered cells (finite block medians) constrain it.
    domains = [block_any(footprint, block) for footprint in footprints]
    surfaces = solve_offset_surfaces(
        backgrounds, domains, smoothness, decay_pixels / block
    )

    total = np.zeros(shape_out, dtype=np.float32)
    weights = np.zeros(shape_out, dtype=np.float32)
    for image, footprint, surface in zip(images, footprints, surfaces):
        weight = feather_weight(footprint, feather_pixels)
        corrected = image - upsample_surface(surface, block, shape_out)
        total += np.where(footprint, corrected * weight, 0.0)
        weights += weight
    mosaic = np.full(shape_out, np.nan, dtype=np.float32)
    covered = weights > 0
    mosaic[covered] = total[covered] / weights[covered]
    return mosaic


def fill_saturated_cores(image: np.ndarray, white: float) -> np.ndarray:
    """Fill the low-valued centres that VISTA leaves in saturated stars."""
    saturated = image >= white
    holes = ndimage.binary_fill_holes(saturated) & ~saturated
    filled = image.copy()
    filled[holes] = white
    return filled


def stretch_channel(
    image: np.ndarray,
    high_percentile: float,
    softening: float,
    black_sigma: float,
) -> np.ndarray:
    """Map sky - *black_sigma* noise to 0 and *high_percentile* to 1 (asinh)."""
    finite = image[np.isfinite(image)]
    sample = finite[:: max(1, finite.size // 5_000_000)]
    _, sky, noise = sigma_clipped_stats(sample, sigma=3.0, maxiters=5)
    black = sky - black_sigma * noise
    white = np.percentile(sample, high_percentile)
    image = fill_saturated_cores(np.nan_to_num(image, nan=sky), white)
    scaled = np.clip((image - black) / (white - black), 0.0, 1.0)
    scaled = np.arcsinh(scaled / softening) / np.arcsinh(1.0 / softening)
    return scaled.astype(np.float32)


def boost_saturation(rgb: np.ndarray, factor: float) -> np.ndarray:
    """Scale each pixel's colour away from its grey level."""
    grey = rgb.mean(axis=-1, keepdims=True)
    return np.clip(grey + factor * (rgb - grey), 0.0, 1.0)
