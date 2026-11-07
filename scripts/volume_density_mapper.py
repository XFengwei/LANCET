#!/usr/bin/env python3
"""Run a resolution-aware test of gxli/volume-density-mapper on G316."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import tomllib

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "lancet-matplotlib"))

import astropy.units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy import ndimage

from lancet.volume_density import (
    apply_resolution_floor,
    cube_storage_bytes,
    density_from_width,
    number_density_h2,
)


DEFAULT_CONFIG = PROJECT_ROOT / "tasks/volume_density_mapper/task.toml"


def _load_config(path: Path) -> tuple[dict[str, object], dict[str, Path]]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    inputs = config["inputs"]
    outputs = config["outputs"]
    paths = {
        name: (PROJECT_ROOT / str(value)).resolve()
        for name, value in inputs.items()
    }
    paths["output_root"] = (PROJECT_ROOT / str(outputs["root"])).resolve()
    return config, paths


def _import_mapper(source: Path | None):
    if source is not None:
        source = source.expanduser().resolve()
        if not source.is_dir():
            raise FileNotFoundError(f"Mapper source directory not found: {source}")
        sys.path.insert(0, str(source))
    try:
        from volume_density_mapper import (  # type: ignore[import-not-found]
            compute_mean_density_width,
            density_reconstruction_3d,
        )
    except ImportError as exc:
        raise RuntimeError(
            "volume-density-mapper is unavailable; install it or pass "
            "--mapper-source /path/to/volume-density-mapper/src"
        ) from exc
    return compute_mean_density_width, density_reconstruction_3d


def _squeezed(path: Path) -> tuple[np.ndarray, fits.Header]:
    data, header = fits.getdata(path, header=True, memmap=True)
    data = np.squeeze(data)
    if data.ndim != 2:
        raise ValueError(f"Expected a 2D image in {path}, found {data.shape}")
    return data, header


def _pixel_arcsec(header: fits.Header) -> float:
    scales = np.abs(
        [
            value.to_value(u.arcsec)
            for value in WCS(header).celestial.proj_plane_pixel_scales()
        ]
    )
    if scales.size != 2 or np.any(~np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("Could not determine a finite celestial pixel scale")
    if not np.isclose(scales[0], scales[1], rtol=1.0e-4):
        raise ValueError(f"Mapper requires square pixels; found {scales}")
    return float(np.mean(scales))


def _pixel_cm(header: fits.Header, distance_pc: float) -> float:
    return _pixel_arcsec(header) * float(distance_pc) * const.au.cgs.value


def _celestial_grids_aligned(
    first: fits.Header,
    second: fits.Header,
    shape: tuple[int, int],
) -> bool:
    y = np.array([0.0, 0.0, shape[0] - 1.0, shape[0] - 1.0, (shape[0] - 1.0) / 2.0])
    x = np.array([0.0, shape[1] - 1.0, 0.0, shape[1] - 1.0, (shape[1] - 1.0) / 2.0])
    first_sky = WCS(first).celestial.pixel_to_world(x, y)
    second_sky = WCS(second).celestial.pixel_to_world(x, y)
    return bool(np.all(first_sky.separation(second_sky).to_value(u.arcsec) < 1.0e-4))


def _cutout_bounds(
    header: fits.Header,
    shape: tuple[int, int],
    *,
    center_ra: str,
    center_dec: str,
    size: int,
) -> tuple[int, int, int, int, SkyCoord]:
    center = SkyCoord(center_ra, center_dec, unit=(u.hourangle, u.deg), frame="icrs")
    x, y = WCS(header).celestial.world_to_pixel(center)
    xc, yc = int(np.rint(x)), int(np.rint(y))
    x0, y0 = xc - size // 2, yc - size // 2
    x1, y1 = x0 + size, y0 + size
    if not (0 <= x0 < x1 <= shape[1] and 0 <= y0 < y1 <= shape[0]):
        raise ValueError(f"Requested cutout {(y0, y1, x0, x1)} exceeds {shape}")
    return y0, y1, x0, x1, center


def _crop_header(header: fits.Header, y0: int, x0: int) -> fits.Header:
    cropped = header.copy()
    cropped["CRPIX1"] = float(cropped["CRPIX1"]) - x0
    cropped["CRPIX2"] = float(cropped["CRPIX2"]) - y0
    return cropped


def _product_header(header: fits.Header, unit: str, history: str) -> fits.Header:
    output = header.copy()
    output["BUNIT"] = unit
    for key in ("BMAJ", "BMIN", "BPA"):
        output.pop(key, None)
    output.add_history(history)
    output.add_history("Local angular resolution is stored in resolution_arcsec.fits")
    return output


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; pass --overwrite to replace it: {path}")


def _write_fits(
    path: Path,
    data: np.ndarray,
    header: fits.Header,
    *,
    overwrite: bool,
) -> None:
    _prepare_output(path, overwrite=overwrite)
    fits.writeto(path, data, header, overwrite=overwrite)
    print(f"Wrote {path}")


def _write_json(path: Path, payload: dict[str, object], *, overwrite: bool) -> None:
    _prepare_output(path, overwrite=overwrite)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {path}")


def _percentiles(data: np.ndarray, mask: np.ndarray) -> dict[str, float | None]:
    values = np.asarray(data)[mask & np.isfinite(data)]
    if values.size == 0:
        return {str(value): None for value in (1, 5, 50, 95, 99)}
    return {
        str(value): float(np.percentile(values, value))
        for value in (1, 5, 50, 95, 99)
    }


def _transition_mask(
    high_resolution_mask: np.ndarray,
    buffer_pixels: int,
) -> np.ndarray:
    if buffer_pixels <= 0:
        return np.zeros(high_resolution_mask.shape, dtype=bool)
    inside = ndimage.distance_transform_edt(high_resolution_mask)
    outside = ndimage.distance_transform_edt(~high_resolution_mask)
    return (inside + outside) <= float(buffer_pixels)


def _quality_mask(
    high_resolution: np.ndarray,
    upstream_defined: np.ndarray,
    floor_applied: np.ndarray,
    transition: np.ndarray,
) -> np.ndarray:
    quality = np.zeros(high_resolution.shape, dtype=np.uint8)
    quality |= high_resolution.astype(np.uint8)
    quality |= upstream_defined.astype(np.uint8) << 1
    quality |= floor_applied.astype(np.uint8) << 2
    quality |= transition.astype(np.uint8) << 3
    return quality


def _run_width_mapper(
    surface_density: np.ndarray,
    dx_cm: float,
    compute_mean_density_width,
    *,
    padding: bool,
) -> tuple[np.ndarray, float]:
    clean = np.clip(np.nan_to_num(surface_density, nan=0.0), 0.0, None).astype(np.float32)
    start = time.perf_counter()
    _, width = compute_mean_density_width(
        clean,
        dx_cm,
        padding=padding,
        verbose=False,
    )
    return np.asarray(width, dtype=np.float64), time.perf_counter() - start


def _block_mean(data: np.ndarray, factor: int) -> np.ndarray:
    ny = data.shape[0] // factor * factor
    nx = data.shape[1] // factor * factor
    trimmed = np.asarray(data[:ny, :nx], dtype=np.float32)
    return trimmed.reshape(ny // factor, factor, nx // factor, factor).mean(axis=(1, 3))


def _summary_regions(
    effective_width_pc: np.ndarray,
    number_density: np.ndarray,
    high_mask: np.ndarray,
    transition: np.ndarray,
    upstream_defined: np.ndarray,
    floor_applied: np.ndarray,
) -> dict[str, object]:
    result: dict[str, object] = {}
    for name, selection in (
        ("high_resolution", high_mask & ~transition),
        ("low_resolution", ~high_mask & ~transition),
        ("transition", transition),
        ("all", np.ones(high_mask.shape, dtype=bool)),
    ):
        count = int(np.count_nonzero(selection))
        positive_density = selection & np.isfinite(number_density) & (number_density > 0)
        result[name] = {
            "pixels": count,
            "positive_density_fraction": (
                float(np.count_nonzero(positive_density) / count) if count else None
            ),
            "upstream_width_defined_fraction": (
                float(np.count_nonzero(upstream_defined & selection) / count) if count else None
            ),
            "resolution_floor_fraction": (
                float(np.count_nonzero(floor_applied & selection) / count) if count else None
            ),
            "effective_width_pc_percentiles": _percentiles(effective_width_pc, selection),
            "number_density_h2_cm-3_percentiles": _percentiles(
                number_density,
                positive_density,
            ),
        }
    return result


def inspect(config: dict[str, object], paths: dict[str, Path]) -> None:
    surface, header = _squeezed(paths["surface_density"])
    footprint, footprint_header = _squeezed(paths["alma_footprint"])
    primary_beam, primary_beam_header = _squeezed(paths["primary_beam"])
    low, low_header = _squeezed(paths["low_resolution_surface_density"])
    if footprint.shape != surface.shape:
        raise ValueError("ALMA footprint and surface-density shapes differ")
    if primary_beam.shape != surface.shape:
        raise ValueError("Primary beam and surface-density shapes differ")
    geometry = config["geometry"]
    size = int(geometry["pilot_size_pixels"])
    y0, y1, x0, x1, center = _cutout_bounds(
        header,
        surface.shape,
        center_ra=str(geometry["center_ra"]),
        center_dec=str(geometry["center_dec"]),
        size=size,
    )
    high_mask = np.isfinite(footprint)
    minimum_response = float(config["primary_beam"]["minimum_response"])
    numerical_pb_mask = np.isfinite(primary_beam) & (primary_beam >= minimum_response)
    pilot_mask = high_mask[y0:y1, x0:x1]
    payload = {
        "surface_density": str(paths["surface_density"]),
        "surface_shape_yx": list(surface.shape),
        "surface_unit": header.get("BUNIT"),
        "native_pixel_arcsec": _pixel_arcsec(header),
        "low_resolution_shape_yx": list(low.shape),
        "low_resolution_pixel_arcsec": _pixel_arcsec(low_header),
        "footprint_celestial_grid_aligned": _celestial_grids_aligned(
            header,
            footprint_header,
            surface.shape,
        ),
        "primary_beam_celestial_grid_aligned": _celestial_grids_aligned(
            header,
            primary_beam_header,
            surface.shape,
        ),
        "primary_beam_minimum_response": minimum_response,
        "primary_beam_retained_pixels": int(np.count_nonzero(numerical_pb_mask)),
        "primary_beam_matches_existing_pbcut": bool(
            np.array_equal(numerical_pb_mask, high_mask)
        ),
        "pilot_center_icrs": center.to_string("hmsdms", precision=2),
        "pilot_bounds_y0_y1_x0_x1": [y0, y1, x0, x1],
        "pilot_shape_yx": [size, size],
        "pilot_high_resolution_fraction": float(pilot_mask.mean()),
        "pilot_native_cube_bytes": cube_storage_bytes((size, size)),
        "full_native_cube_bytes": cube_storage_bytes(surface.shape),
        "warning": (
            "The mapper assumes one dx and has no spatial beam model; use the "
            "resolution map and quality flags produced by this wrapper."
        ),
    }
    print(json.dumps(payload, indent=2))


def _pilot_plot(
    root: Path,
    surface: np.ndarray,
    high_mask: np.ndarray,
    width_pc: np.ndarray,
    number_density: np.ndarray,
    quality: np.ndarray,
    *,
    overwrite: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    path = root / "figures/pilot_diagnostics.png"
    _prepare_output(path, overwrite=overwrite)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    finite_surface = surface[np.isfinite(surface) & (surface > 0)]
    lo, hi = np.percentile(finite_surface, [1, 99.7])
    panels = (
        (np.log10(np.clip(surface, lo, None)), "magma", "log surface density"),
        (high_mask.astype(float), "gray", "ALMA high-resolution footprint"),
        (np.log10(width_pc), "viridis", "log effective width [pc]"),
        (
            np.log10(np.where(number_density > 0, number_density, np.nan)),
            "plasma",
            "log n(H2) [cm$^{-3}$]",
        ),
        (((quality & 4) > 0).astype(float), "gray", "resolution floor applied"),
        (((quality & 8) > 0).astype(float), "gray", "transition exclusion zone"),
    )
    for axis, (image, cmap, title) in zip(axes.flat, panels):
        shown = axis.imshow(image, origin="lower", cmap=cmap)
        axis.set_title(title)
        axis.set_xlabel("x pixel")
        axis.set_ylabel("y pixel")
        if image.dtype.kind == "f" and np.unique(image[np.isfinite(image)]).size > 2:
            fig.colorbar(shown, ax=axis, shrink=0.8)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Wrote {path}")


def _full_plot(
    root: Path,
    surface: np.ndarray,
    width_pc: np.ndarray,
    number_density: np.ndarray,
    resolution_arcsec: np.ndarray,
    quality: np.ndarray,
    *,
    overwrite: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    path = root / "figures/full_diagnostics.png"
    _prepare_output(path, overwrite=overwrite)
    step = max(1, int(np.ceil(max(surface.shape) / 1600)))
    view = np.s_[::step, ::step]
    positive_surface = surface[np.isfinite(surface) & (surface > 1.0e-3)]
    lo, hi = np.percentile(positive_surface, [1, 99.7])
    log_surface = np.log10(np.clip(surface[view], lo, None))
    log_width = np.log10(width_pc[view])
    log_density = np.log10(
        np.where(number_density[view] > 0, number_density[view], np.nan)
    )
    panels = (
        (
            log_surface,
            "magma",
            "log surface density",
            np.log10(lo),
            np.log10(hi),
        ),
        (
            log_width,
            "viridis",
            "log effective width [pc]",
            *np.nanpercentile(log_width, [1, 99]),
        ),
        (
            log_density,
            "plasma",
            "log n(H2) [cm$^{-3}$]",
            *np.nanpercentile(log_density, [1, 99.7]),
        ),
        (
            resolution_arcsec[view],
            "cividis",
            "effective resolution [arcsec]",
            None,
            None,
        ),
        (
            ((quality[view] & 4) > 0).astype(float),
            "gray",
            "resolution floor applied",
            0.0,
            1.0,
        ),
        (
            ((quality[view] & 8) > 0).astype(float),
            "gray",
            "transition exclusion zone",
            0.0,
            1.0,
        ),
    )
    fig, axes = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
    for axis, (image, cmap, title, vmin, vmax) in zip(axes.flat, panels):
        shown = axis.imshow(image, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        axis.set_title(title)
        axis.set_xlabel(f"x / {step} native pixels")
        axis.set_ylabel(f"y / {step} native pixels")
        if image.dtype.kind == "f" and np.unique(image[np.isfinite(image)]).size > 2:
            fig.colorbar(shown, ax=axis, shrink=0.8)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Wrote {path}")


def pilot(config: dict[str, object], paths: dict[str, Path], args: argparse.Namespace) -> None:
    compute_mean_density_width, density_reconstruction_3d = _import_mapper(args.mapper_source)
    surface_full, header = _squeezed(paths["surface_density"])
    footprint, _ = _squeezed(paths["alma_footprint"])
    geometry = config["geometry"]
    resolution = config["resolution"]
    model = config["model"]
    size = int(geometry["pilot_size_pixels"])
    distance_pc = float(geometry["distance_pc"])
    y0, y1, x0, x1, center = _cutout_bounds(
        header,
        surface_full.shape,
        center_ra=str(geometry["center_ra"]),
        center_dec=str(geometry["center_dec"]),
        size=size,
    )
    surface = np.asarray(surface_full[y0:y1, x0:x1], dtype=np.float32)
    surface_clean = np.clip(np.nan_to_num(surface, nan=0.0), 0.0, None)
    high_mask = np.isfinite(footprint[y0:y1, x0:x1])
    crop_header = _crop_header(header, y0, x0)
    dx_cm = _pixel_cm(crop_header, distance_pc)

    if args.cached_native is not None:
        cached = np.load(args.cached_native)
        width_raw = np.asarray(cached["width"], dtype=np.float64)
        if width_raw.shape != surface.shape:
            raise ValueError(f"Cached width shape {width_raw.shape} does not match {surface.shape}")
        mapper_seconds = None
    else:
        width_raw, mapper_seconds = _run_width_mapper(
            surface_clean,
            dx_cm,
            compute_mean_density_width,
            padding=bool(model["padding"]),
        )

    resolution_arcsec = np.where(
        high_mask,
        float(resolution["high_arcsec"]),
        float(resolution["low_arcsec"]),
    ).astype(np.float32)
    effective_width, upstream_defined, floor_applied = apply_resolution_floor(
        width_raw,
        resolution_arcsec,
        distance_pc=distance_pc,
    )
    effective_density = density_from_width(surface_clean, effective_width)
    number_density = number_density_h2(
        effective_density,
        mean_molecular_weight=float(model["mean_molecular_weight"]),
    )
    buffer_pixels = int(
        np.ceil(float(resolution["transition_buffer_arcsec"]) / _pixel_arcsec(crop_header))
    )
    transition = _transition_mask(high_mask, buffer_pixels)
    quality = _quality_mask(high_mask, upstream_defined, floor_applied, transition)
    root = paths["output_root"] / "pilot"
    history = "LANCET resolution-aware gxli volume-density-mapper pilot"

    _write_fits(root / "surface_density_input.fits", surface, crop_header, overwrite=args.overwrite)
    _write_fits(
        root / "resolution_arcsec.fits",
        resolution_arcsec,
        _product_header(crop_header, "arcsec", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "width_characteristic_raw_pc.fits",
        (width_raw / const.pc.cgs.value).astype(np.float32),
        _product_header(crop_header, "pc", history),
        overwrite=args.overwrite,
    )
    effective_width_pc = effective_width / const.pc.cgs.value
    _write_fits(
        root / "width_effective_pc.fits",
        effective_width_pc.astype(np.float32),
        _product_header(crop_header, "pc", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "mass_density_g_cm3.fits",
        effective_density.astype(np.float32),
        _product_header(crop_header, "g cm^-3", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "number_density_h2_cm3.fits",
        number_density.astype(np.float32),
        _product_header(crop_header, "cm^-3", history),
        overwrite=args.overwrite,
    )
    quality_header = _product_header(crop_header, "bitmask", history)
    quality_header["QBIT0"] = "inside high-resolution ALMA footprint"
    quality_header["QBIT1"] = "upstream width finite and positive"
    quality_header["QBIT2"] = "local resolution floor applied"
    quality_header["QBIT3"] = "within transition exclusion zone"
    _write_fits(root / "quality_mask.fits", quality, quality_header, overwrite=args.overwrite)

    cube_summary: dict[str, object] | None = None
    if not args.skip_cube_smoke:
        factor = 6
        reduced = _block_mean(surface_clean, factor)
        reduced_dx = dx_cm * factor
        cube_start = time.perf_counter()
        if args.cached_cube is not None:
            cube = np.load(args.cached_cube)
            if cube.shape != (reduced.shape[0], *reduced.shape):
                raise ValueError(f"Cached cube shape {cube.shape} does not match reduced pilot")
            cube_seconds = None
        else:
            cube = density_reconstruction_3d(
                reduced,
                reduced_dx,
                padding=bool(model["padding"]),
                verbose=False,
            )
            cube_seconds = time.perf_counter() - cube_start
        recovered = cube.sum(axis=0, dtype=np.float64) * reduced_dx
        good = reduced > 0
        ratio = recovered[good] / reduced[good]
        cube_header = crop_header.copy()
        cube_header["CDELT1"] = float(cube_header["CDELT1"]) * factor
        cube_header["CDELT2"] = float(cube_header["CDELT2"]) * factor
        cube_header["CRPIX1"] = (float(cube_header["CRPIX1"]) - 0.5) / factor + 0.5
        cube_header["CRPIX2"] = (float(cube_header["CRPIX2"]) - 0.5) / factor + 0.5
        cube_header["CTYPE3"] = "DEPTH"
        cube_header["CUNIT3"] = "cm"
        cube_header["CDELT3"] = reduced_dx
        cube_header["CRPIX3"] = cube.shape[0] // 2 + 1
        cube_header["CRVAL3"] = 0.0
        cube_header["BUNIT"] = "g cm^-3"
        cube_header.add_history("Upstream 250^3 smoke cube; diffuse residual is not restored")
        _write_fits(
            root / "density_cube_smoke_250.fits",
            np.asarray(cube, dtype=np.float32),
            cube_header,
            overwrite=args.overwrite,
        )
        cube_summary = {
            "shape": list(cube.shape),
            "runtime_seconds": cube_seconds,
            "projected_column_recovery_ratio_percentiles": {
                str(value): float(np.percentile(ratio, value))
                for value in (0, 1, 50, 99, 100)
            },
            "scientific_gate_passed": bool(np.percentile(ratio, 1) >= 0.95),
            "gate": "at least 95% projected-column recovery at the first percentile",
        }

    summary = {
        "upstream_repository": "https://github.com/gxli/volume-density-mapper",
        "center_icrs": center.to_string("hmsdms", precision=2),
        "bounds_y0_y1_x0_x1": [y0, y1, x0, x1],
        "shape_yx": list(surface.shape),
        "pixel_arcsec": _pixel_arcsec(crop_header),
        "pixel_cm": dx_cm,
        "distance_pc": distance_pc,
        "mapper_padding": bool(model["padding"]),
        "mapper_runtime_seconds": mapper_seconds,
        "nonpositive_input_pixels_clipped": int(np.count_nonzero(~np.isfinite(surface) | (surface < 0))),
        "high_resolution_fraction": float(high_mask.mean()),
        "transition_buffer_pixels": buffer_pixels,
        "regions": _summary_regions(
            effective_width_pc,
            number_density,
            high_mask,
            transition,
            upstream_defined,
            floor_applied,
        ),
        "native_3d_storage_bytes": cube_storage_bytes(surface.shape),
        "cube_smoke": cube_summary,
        "full_3d_authorized": False,
        "full_3d_blocker": (
            "Current upstream cube omits the constrained-diffusion residual; "
            "the pilot projection-closure gate must pass before full 3-D work."
        ),
    }
    _write_json(root / "summary.json", summary, overwrite=args.overwrite)
    _pilot_plot(
        root,
        surface,
        high_mask,
        effective_width_pc,
        number_density,
        quality,
        overwrite=args.overwrite,
    )
    print(json.dumps(summary, indent=2))


def full(config: dict[str, object], paths: dict[str, Path], args: argparse.Namespace) -> None:
    from reproject import reproject_interp

    compute_mean_density_width, _ = _import_mapper(args.mapper_source)
    surface, header = _squeezed(paths["surface_density"])
    footprint, _ = _squeezed(paths["alma_footprint"])
    low_surface, low_header = _squeezed(paths["low_resolution_surface_density"])
    geometry = config["geometry"]
    resolution = config["resolution"]
    model = config["model"]
    full_config = config["full_map"]
    distance_pc = float(geometry["distance_pc"])
    high_mask = np.isfinite(footprint)
    yy, xx = np.where(high_mask)
    if yy.size == 0:
        raise ValueError("ALMA footprint contains no finite pixels")
    factor = int(full_config["high_resolution_rebin"])
    y0 = int(yy.min()) // factor * factor
    x0 = int(xx.min()) // factor * factor
    y1 = min(surface.shape[0], int(np.ceil((yy.max() + 1) / factor)) * factor)
    x1 = min(surface.shape[1], int(np.ceil((xx.max() + 1) / factor)) * factor)
    y1 = y0 + (y1 - y0) // factor * factor
    x1 = x0 + (x1 - x0) // factor * factor

    clean_low = np.clip(np.nan_to_num(low_surface, nan=0.0), 0.0, None).astype(np.float32)
    if args.cached_low_width is not None:
        low_width = np.load(args.cached_low_width)
        low_seconds = None
    else:
        low_width, low_seconds = _run_width_mapper(
            clean_low,
            _pixel_cm(low_header, distance_pc),
            compute_mean_density_width,
            padding=bool(model["padding"]),
        )
    low_width_full = reproject_interp(
        (low_width, WCS(low_header).celestial),
        WCS(header).celestial,
        shape_out=surface.shape,
        order="bilinear",
        return_footprint=False,
    )

    clean_high = np.clip(
        np.nan_to_num(surface[y0:y1, x0:x1], nan=0.0),
        0.0,
        None,
    ).astype(np.float32)
    reduced_high = _block_mean(clean_high, factor)
    if args.cached_high_width is not None:
        high_width_reduced = np.load(args.cached_high_width)
        high_seconds = None
    else:
        high_width_reduced, high_seconds = _run_width_mapper(
            reduced_high,
            _pixel_cm(header, distance_pc) * factor,
            compute_mean_density_width,
            padding=bool(model["padding"]),
        )
    if high_width_reduced.shape != reduced_high.shape:
        raise ValueError("Cached or calculated high-resolution width has the wrong shape")
    high_width = ndimage.zoom(
        high_width_reduced,
        factor,
        order=1,
        mode="nearest",
        prefilter=False,
    )[: y1 - y0, : x1 - x0]

    raw_width = np.asarray(low_width_full, dtype=np.float64)
    local_mask = high_mask[y0:y1, x0:x1]
    raw_width[y0:y1, x0:x1][local_mask] = high_width[local_mask]
    resolution_arcsec = np.where(
        high_mask,
        float(resolution["high_arcsec"]),
        float(resolution["low_arcsec"]),
    ).astype(np.float32)
    effective_width, upstream_defined, floor_applied = apply_resolution_floor(
        raw_width,
        resolution_arcsec,
        distance_pc=distance_pc,
    )
    clean_surface = np.clip(np.nan_to_num(surface, nan=0.0), 0.0, None)
    mass_density = density_from_width(clean_surface, effective_width)
    number_density = number_density_h2(
        mass_density,
        mean_molecular_weight=float(model["mean_molecular_weight"]),
    )

    buffer_pixels = int(
        np.ceil(float(resolution["transition_buffer_arcsec"]) / _pixel_arcsec(header))
    )
    margin = buffer_pixels + 2
    by0, by1 = max(0, y0 - margin), min(surface.shape[0], y1 + margin)
    bx0, bx1 = max(0, x0 - margin), min(surface.shape[1], x1 + margin)
    transition = np.zeros(surface.shape, dtype=bool)
    transition[by0:by1, bx0:bx1] = _transition_mask(
        high_mask[by0:by1, bx0:bx1],
        buffer_pixels,
    )
    quality = _quality_mask(high_mask, upstream_defined, floor_applied, transition)
    effective_width_pc = effective_width / const.pc.cgs.value
    root = paths["output_root"] / "full"
    history = (
        "LANCET resolution-aware mapper: 0.65arcsec branch rebinned 3x; "
        "18arcsec branch from low-resolution map"
    )
    _write_fits(
        root / "resolution_arcsec.fits",
        resolution_arcsec,
        _product_header(header, "arcsec", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "width_characteristic_raw_pc.fits",
        (raw_width / const.pc.cgs.value).astype(np.float32),
        _product_header(header, "pc", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "width_effective_pc.fits",
        effective_width_pc.astype(np.float32),
        _product_header(header, "pc", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "mass_density_g_cm3.fits",
        mass_density.astype(np.float32),
        _product_header(header, "g cm^-3", history),
        overwrite=args.overwrite,
    )
    _write_fits(
        root / "number_density_h2_cm3.fits",
        number_density.astype(np.float32),
        _product_header(header, "cm^-3", history),
        overwrite=args.overwrite,
    )
    quality_header = _product_header(header, "bitmask", history)
    quality_header["QBIT0"] = "inside high-resolution ALMA footprint"
    quality_header["QBIT1"] = "upstream width finite and positive"
    quality_header["QBIT2"] = "local resolution floor applied"
    quality_header["QBIT3"] = "within transition exclusion zone"
    _write_fits(root / "quality_mask.fits", quality, quality_header, overwrite=args.overwrite)
    summary = {
        "shape_yx": list(surface.shape),
        "distance_pc": distance_pc,
        "native_pixel_arcsec": _pixel_arcsec(header),
        "high_resolution_rebin": factor,
        "high_resolution_analysis_pixel_arcsec": _pixel_arcsec(header) * factor,
        "high_resolution_crop_y0_y1_x0_x1": [y0, y1, x0, x1],
        "low_resolution_mapper_runtime_seconds": low_seconds,
        "high_resolution_mapper_runtime_seconds": high_seconds,
        "transition_buffer_pixels": buffer_pixels,
        "regions": _summary_regions(
            effective_width_pc,
            number_density,
            high_mask,
            transition,
            upstream_defined,
            floor_applied,
        ),
        "full_native_3d_storage_bytes": cube_storage_bytes(surface.shape),
        "full_3d_written": False,
        "full_3d_blocker": "Pilot projection closure failed because upstream omits residual.",
    }
    _write_json(root / "summary.json", summary, overwrite=args.overwrite)
    _full_plot(
        root,
        clean_surface,
        effective_width_pc,
        number_density,
        resolution_arcsec,
        quality,
        overwrite=args.overwrite,
    )
    print(json.dumps(summary, indent=2))


def _pbmasked_plot(
    root: Path,
    surface: np.ndarray,
    width_pc: np.ndarray,
    number_density: np.ndarray,
    primary_beam: np.ndarray,
    retained: np.ndarray,
    *,
    overwrite: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    path = root / "figures/primary_beam_masked_diagnostics.png"
    _prepare_output(path, overwrite=overwrite)
    step = max(1, int(np.ceil(max(surface.shape) / 1600)))
    view = np.s_[::step, ::step]
    display_surface = np.where(retained, surface, np.nan)[view]
    display_width = np.where(retained, width_pc, np.nan)[view]
    display_density = np.where(retained, number_density, np.nan)[view]
    display_pb = np.where(retained, primary_beam, np.nan)[view]
    positive_surface = display_surface[np.isfinite(display_surface) & (display_surface > 0)]
    positive_density = display_density[np.isfinite(display_density) & (display_density > 0)]
    width_values = display_width[np.isfinite(display_width) & (display_width > 0)]
    surface_limits = np.log10(np.percentile(positive_surface, [1, 99.7]))
    density_limits = np.log10(np.percentile(positive_density, [1, 99.7]))
    width_limits = np.log10(np.percentile(width_values, [1, 99]))
    panels = (
        (
            np.log10(np.where(display_surface > 0, display_surface, np.nan)),
            "magma",
            "PB-masked log surface density",
            *surface_limits,
        ),
        (
            np.log10(display_width),
            "viridis",
            "PB-masked log effective width [pc]",
            *width_limits,
        ),
        (
            np.log10(np.where(display_density > 0, display_density, np.nan)),
            "plasma",
            "PB-masked log n(H2) [cm$^{-3}$]",
            *density_limits,
        ),
        (display_pb, "cividis", "Retained primary-beam response", None, None),
    )
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    for axis, (image, cmap, title, vmin, vmax) in zip(axes.flat, panels):
        shown = axis.imshow(image, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        axis.set_title(title)
        axis.set_xlabel(f"x / {step} native pixels")
        axis.set_ylabel(f"y / {step} native pixels")
        fig.colorbar(shown, ax=axis, shrink=0.82)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Wrote {path}")


def primary_beam_mask(
    config: dict[str, object],
    paths: dict[str, Path],
    args: argparse.Namespace,
) -> None:
    """Create authoritative PB-cut copies of the pilot and full 2-D products."""
    surface, surface_header = _squeezed(paths["surface_density"])
    primary_beam, primary_beam_header = _squeezed(paths["primary_beam"])
    if primary_beam.shape != surface.shape or not _celestial_grids_aligned(
        surface_header,
        primary_beam_header,
        surface.shape,
    ):
        raise ValueError("Primary beam is not pixel-aligned with the surface-density map")
    threshold = float(config["primary_beam"]["minimum_response"])
    if not 0.0 < threshold <= 1.0:
        raise ValueError("Primary-beam minimum response must lie in (0, 1]")
    full_mask = np.isfinite(primary_beam) & (primary_beam >= threshold)
    geometry = config["geometry"]
    size = int(geometry["pilot_size_pixels"])
    y0, y1, x0, x1, _ = _cutout_bounds(
        surface_header,
        surface.shape,
        center_ra=str(geometry["center_ra"]),
        center_dec=str(geometry["center_dec"]),
        size=size,
    )
    product_names = (
        "resolution_arcsec.fits",
        "width_characteristic_raw_pc.fits",
        "width_effective_pc.fits",
        "mass_density_g_cm3.fits",
        "number_density_h2_cm3.fits",
    )
    summaries: dict[str, object] = {}
    for name, selection in (
        ("pilot", (slice(y0, y1), slice(x0, x1))),
        ("full", (slice(None), slice(None))),
    ):
        source_root = paths["output_root"] / name
        output_root = paths["output_root"] / f"{name}_pbmasked"
        region_mask = full_mask[selection]
        region_pb = np.asarray(primary_beam[selection], dtype=np.float32)
        region_surface = np.asarray(surface[selection], dtype=np.float32)
        source_reference = source_root / "width_effective_pc.fits"
        if not source_reference.is_file():
            raise FileNotFoundError(f"Run the {name} stage before PB masking: {source_reference}")
        reference_header = fits.getheader(source_reference)
        history = f"LANCET primary-beam response mask: retained PB >= {threshold:.2f}"
        masked_surface = np.where(region_mask, region_surface, np.nan).astype(np.float32)
        _write_fits(
            output_root / "surface_density_input_pbmasked.fits",
            masked_surface,
            _product_header(reference_header, "g cm^-2", history),
            overwrite=args.overwrite,
        )
        response_header = _product_header(reference_header, "1", history)
        _write_fits(
            output_root / "primary_beam_response.fits",
            region_pb,
            response_header,
            overwrite=args.overwrite,
        )
        mask_header = _product_header(reference_header, "mask", history)
        mask_header["PBMIN"] = (threshold, "Minimum retained primary-beam response")
        _write_fits(
            output_root / "primary_beam_mask.fits",
            region_mask.astype(np.uint8),
            mask_header,
            overwrite=args.overwrite,
        )
        finite_counts: dict[str, int] = {}
        masked_products: dict[str, np.ndarray] = {}
        for product_name in product_names:
            data, product_header = fits.getdata(
                source_root / product_name,
                header=True,
                memmap=True,
            )
            if data.shape != region_mask.shape:
                raise ValueError(
                    f"{source_root / product_name} shape {data.shape} does not match "
                    f"PB mask {region_mask.shape}"
                )
            masked = np.where(region_mask, data, np.nan).astype(np.float32)
            product_header.add_history(history)
            _write_fits(
                output_root / product_name,
                masked,
                product_header,
                overwrite=args.overwrite,
            )
            if np.any(np.isfinite(masked[~region_mask])):
                raise RuntimeError(f"PB masking failed outside retained region for {product_name}")
            finite_counts[product_name] = int(np.count_nonzero(np.isfinite(masked)))
            masked_products[product_name] = masked
        quality, quality_header = fits.getdata(
            source_root / "quality_mask.fits",
            header=True,
            memmap=True,
        )
        masked_quality = np.where(region_mask, quality, 0).astype(np.uint8)
        quality_header.add_history(history + "; zero marks rejected pixels")
        _write_fits(
            output_root / "quality_mask.fits",
            masked_quality,
            quality_header,
            overwrite=args.overwrite,
        )
        summary = {
            "source_directory": str(source_root),
            "primary_beam": str(paths["primary_beam"]),
            "minimum_response": threshold,
            "shape_yx": list(region_mask.shape),
            "retained_pixels": int(np.count_nonzero(region_mask)),
            "rejected_pixels": int(region_mask.size - np.count_nonzero(region_mask)),
            "retained_fraction": float(region_mask.mean()),
            "outside_mask_is_nan_for_all_float_products": True,
            "finite_pixels_by_product": finite_counts,
            "quality_zero_outside_mask": bool(np.all(masked_quality[~region_mask] == 0)),
        }
        _write_json(output_root / "summary.json", summary, overwrite=args.overwrite)
        _pbmasked_plot(
            output_root,
            region_surface,
            masked_products["width_effective_pc.fits"],
            masked_products["number_density_h2_cm3.fits"],
            region_pb,
            region_mask,
            overwrite=args.overwrite,
        )
        summaries[name] = summary
    print(json.dumps(summaries, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("inspect", "pilot", "full", "mask"))
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--mapper-source", type=Path)
    parser.add_argument("--cached-native", type=Path)
    parser.add_argument("--cached-cube", type=Path)
    parser.add_argument("--cached-low-width", type=Path)
    parser.add_argument("--cached-high-width", type=Path)
    parser.add_argument("--skip-cube-smoke", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    task_config, task_paths = _load_config(arguments.config)
    if arguments.stage == "inspect":
        inspect(task_config, task_paths)
    elif arguments.stage == "pilot":
        pilot(task_config, task_paths, arguments)
    elif arguments.stage == "full":
        full(task_config, task_paths, arguments)
    else:
        primary_beam_mask(task_config, task_paths, arguments)
