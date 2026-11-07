#!/usr/bin/env python3
"""Prepare and summarize a non-destructive SA-JEPA continuum pilot.

The training itself is deliberately run with Guang-Xing Li's published
SA-JEPA checkout.  This helper only bridges a LANCET FITS product to the
package's scalar-field input and turns the saved latent field into compact,
discussion-ready diagnostic figures.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / (
    "data/delivered/alma/TP7MTM2TM1/cont/"
    "allscale_combined_v3_dust.MJysr.rp.pbcut.fits"
)
DEFAULT_ROOT = PROJECT_ROOT / "results/sa_jepa_continuum_pilot"


def _finite_bbox(data: np.ndarray) -> tuple[int, int, int, int]:
    yy, xx = np.where(np.isfinite(data))
    if yy.size == 0:
        raise ValueError("Input contains no finite pixels")
    return int(yy.min()), int(yy.max()) + 1, int(xx.min()), int(xx.max()) + 1


def _weighted_resize(data: np.ndarray, out_shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(data)
    values = np.where(valid, data, 0.0).astype(np.float32)
    weights = valid.astype(np.float32)
    zoom = (out_shape[0] / data.shape[0], out_shape[1] / data.shape[1])
    numerator = ndimage.zoom(values, zoom, order=1, mode="nearest", prefilter=False)
    denominator = ndimage.zoom(weights, zoom, order=1, mode="nearest", prefilter=False)
    resized = np.full(out_shape, np.nan, dtype=np.float32)
    good = denominator > 1.0e-4
    resized[good] = numerator[good] / denominator[good]
    return resized, denominator.astype(np.float32)


def prepare(args: argparse.Namespace) -> None:
    out_root = args.output_root.resolve()
    data_dir = out_root / "data"
    figure_dir = out_root / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    with fits.open(args.input, memmap=True) as hdul:
        raw = np.squeeze(np.asarray(hdul[0].data))
        header = hdul[0].header.copy()
    if raw.ndim != 2:
        raise ValueError(f"Expected a 2D FITS image after squeeze, got {raw.shape}")

    requested_center = None
    if args.center_ra is not None or args.center_dec is not None or args.native_crop_size is not None:
        if args.center_ra is None or args.center_dec is None or args.native_crop_size is None:
            raise ValueError("--center-ra, --center-dec, and --native-crop-size must be supplied together")
        requested_center = SkyCoord(
            args.center_ra,
            args.center_dec,
            unit=(u.hourangle, u.deg),
            frame="icrs",
        )
        center_x, center_y = map(float, WCS(header).celestial.world_to_pixel(requested_center))
        crop_size = int(args.native_crop_size)
        x0 = int(round(center_x)) - crop_size // 2
        y0 = int(round(center_y)) - crop_size // 2
        x1 = x0 + crop_size
        y1 = y0 + crop_size
        if not (0 <= x0 < x1 <= raw.shape[1] and 0 <= y0 < y1 <= raw.shape[0]):
            raise ValueError(
                f"Requested native crop {(y0, y1, x0, x1)} is outside source shape {raw.shape}"
            )
    else:
        y0, y1, x0, x1 = _finite_bbox(raw)
    cropped = raw[y0:y1, x0:x1]
    scale = min(1.0, args.max_dimension / max(cropped.shape))
    out_shape = tuple(max(32, int(round(size * scale))) for size in cropped.shape)
    field, coverage = _weighted_resize(cropped, out_shape)
    full_prepared_shape = field.shape
    prepared_crop = None
    if args.prepared_crop is not None:
        py0, py1, px0, px1 = args.prepared_crop
        if not (0 <= py0 < py1 <= field.shape[0] and 0 <= px0 < px1 <= field.shape[1]):
            raise ValueError(
                f"Prepared crop {(py0, py1, px0, px1)} is outside field shape {field.shape}"
            )
        field = field[py0:py1, px0:px1]
        coverage = coverage[py0:py1, px0:px1]
        prepared_crop = [py0, py1, px0, px1]
    data_mask = coverage >= args.minimum_coverage
    target_mask = data_mask.copy()
    snr_threshold_value = None
    snr_background_value = None
    if args.snr_min is not None:
        if args.noise_rms is None or args.noise_rms <= 0:
            raise ValueError("--noise-rms must be positive when --snr-min is used")
        if not (0 <= args.background_percentile <= 100):
            raise ValueError("--background-percentile must lie in [0, 100]")
        snr_background_value = float(np.percentile(field[data_mask], args.background_percentile))
        snr_threshold_value = snr_background_value + float(args.snr_min) * float(args.noise_rms)
        target_mask &= field >= snr_threshold_value
    if args.target_mask_erosion > 0:
        target_mask = ndimage.binary_erosion(
            target_mask,
            iterations=args.target_mask_erosion,
        )
    # Erosion limits eligible target centres only. Keep every reliably
    # resampled datum available to the model as surrounding context.
    field[~data_mask] = np.nan

    np.save(data_dir / "input.npy", field.astype(np.float32))
    np.save(data_dir / "target_mask.npy", target_mask.astype(np.float32))

    if requested_center is not None:
        crop_header = header.copy()
        crop_header["CRPIX1"] = float(crop_header["CRPIX1"]) - x0
        crop_header["CRPIX2"] = float(crop_header["CRPIX2"]) - y0
        crop_header["HISTORY"] = "LANCET SA-JEPA exact native cutout; delivered input unchanged"
        fits.PrimaryHDU(np.asarray(cropped, dtype=np.float32), header=crop_header).writeto(
            data_dir / "native_crop.fits",
            overwrite=False,
        )

    cdelt1 = abs(float(header.get("CDELT1", np.nan))) * 3600.0
    cdelt2 = abs(float(header.get("CDELT2", np.nan))) * 3600.0
    resize_y = cropped.shape[0] / out_shape[0]
    resize_x = cropped.shape[1] / out_shape[1]
    pixel_arcsec = float(np.nanmean([cdelt1 * resize_x, cdelt2 * resize_y]))
    pixel_pc = args.distance_kpc * 1000.0 * np.deg2rad(pixel_arcsec / 3600.0)
    valid_values = field[np.isfinite(field)]
    metadata = {
        "source_fits": str(args.input.resolve()),
        "source_bunit": header.get("BUNIT"),
        "source_shape_yx": list(raw.shape),
        "source_crop_y0_y1_x0_x1": [y0, y1, x0, x1],
        "crop_selection": "sky_center" if requested_center is not None else "finite_bbox",
        "requested_center_icrs": (
            {
                "ra": requested_center.ra.to_string(unit=u.hourangle, sep=":", precision=2),
                "dec": requested_center.dec.to_string(unit=u.deg, sep=":", precision=2, alwayssign=True),
            }
            if requested_center is not None
            else None
        ),
        "native_crop_shape_yx": list(cropped.shape),
        "full_prepared_shape_yx": list(full_prepared_shape),
        "prepared_crop_y0_y1_x0_x1": prepared_crop,
        "prepared_shape_yx": list(field.shape),
        "prepared_valid_fraction": float(np.isfinite(field).mean()),
        "target_mask_fraction": float(target_mask.mean()),
        "minimum_resampling_coverage": args.minimum_coverage,
        "target_mask_erosion_pixels": args.target_mask_erosion,
        "snr_mask": {
            "snr_min": args.snr_min,
            "noise_rms_mjysr": args.noise_rms,
            "background_percentile": args.background_percentile if args.snr_min is not None else None,
            "background_mjysr": snr_background_value,
            "threshold_mjysr": snr_threshold_value,
        },
        "pixel_scale_arcsec_approx": pixel_arcsec,
        "pixel_scale_pc_approx": pixel_pc,
        "assumed_distance_kpc": args.distance_kpc,
        "cdd_scales_pixels": [2, 4, 8, 16, 32],
        "cdd_scales_pc_approx": [float(pixel_pc * s) for s in [2, 4, 8, 16, 32]],
        "finite_percentiles": {
            str(p): float(np.percentile(valid_values, p))
            for p in [0, 1, 5, 50, 95, 99, 100]
        },
    }
    (out_root / "preparation_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    fig, ax = plt.subplots(figsize=(11, 7), constrained_layout=True)
    lo, hi = np.nanpercentile(field, [1, 99.7])
    image = ax.imshow(field, origin="lower", cmap="magma", vmin=lo, vmax=hi)
    ax.contour(target_mask, levels=[0.5], colors="cyan", linewidths=0.5, alpha=0.7)
    ax.set(title="SA-JEPA pilot input: G316 all-scale dust continuum", xlabel="resampled x pixel", ylabel="resampled y pixel")
    fig.colorbar(image, ax=ax, label=header.get("BUNIT", "input units"))
    fig.savefig(figure_dir / "prepared_input.png", dpi=180)
    plt.close(fig)
    print(json.dumps(metadata, indent=2))


def _latent_from_outputs(outputs: dict) -> tuple[str, np.ndarray]:
    # SA-JEPA calls the dense EMA-target embedding ``gt_map``.  Its
    # ``target_map`` key is a one-channel sampling diagnostic, not a latent.
    for key in ("gt_map", "pred_map", "context_map"):
        value = outputs.get(key)
        if value is None:
            continue
        array = value.detach().cpu().numpy() if hasattr(value, "detach") else np.asarray(value)
        if array.ndim == 4:
            array = array[0]
        if array.ndim == 3 and array.shape[0] >= 3:
            return key, np.asarray(array, dtype=np.float32)
    available = ", ".join(sorted(outputs))
    raise KeyError(f"No dense latent map found; available keys: {available}")


def summarize(args: argparse.Namespace) -> None:
    import torch

    out_root = args.output_root.resolve()
    session_dir = args.session_dir.resolve()
    figure_dir = out_root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.summary_name == "pilot" else f"_{args.summary_name}"
    outputs = torch.load(session_dir / "inference_outputs.pt", map_location="cpu", weights_only=False)
    latent_key, latent = _latent_from_outputs(outputs)
    _, height, width = latent.shape
    input_field = np.load(out_root / "data/input.npy")
    target_mask = np.load(out_root / "data/target_mask.npy") > 0
    if args.valid_erosion > 0:
        target_mask = ndimage.binary_erosion(target_mask, iterations=args.valid_erosion)
    if input_field.shape != (height, width):
        input_field = ndimage.zoom(input_field, (height / input_field.shape[0], width / input_field.shape[1]), order=1)
        target_mask = ndimage.zoom(target_mask.astype(float), (height / target_mask.shape[0], width / target_mask.shape[1]), order=0) > 0

    tokens = latent.transpose(1, 2, 0).reshape(-1, latent.shape[0])
    valid = target_mask.reshape(-1) & np.isfinite(tokens).all(axis=1) & np.isfinite(input_field).reshape(-1)
    z = tokens[valid]
    if z.shape[0] < 100:
        raise RuntimeError(f"Only {z.shape[0]} valid latent tokens found")

    pca = PCA(n_components=3, random_state=42)
    xyz_valid = pca.fit_transform(z)
    xyz = np.full((tokens.shape[0], 3), np.nan, dtype=np.float32)
    xyz[valid] = xyz_valid
    xyz_map = xyz.reshape(height, width, 3)
    rgb = np.zeros_like(xyz_map)
    for channel in range(3):
        values = xyz_valid[:, channel]
        lo, hi = np.percentile(values, [1, 99])
        rgb[..., channel] = np.clip((xyz_map[..., channel] - lo) / max(hi - lo, 1e-12), 0, 1)
    rgb[~valid.reshape(height, width)] = np.nan

    umap_xyz_valid = None
    umap_xyz_map = None
    umap_rgb = None
    if args.umap:
        import torchdr

        umap_input = np.asarray(z, dtype=np.float32)
        if args.umap_standardize:
            umap_mean = umap_input.mean(axis=0, keepdims=True)
            umap_std = umap_input.std(axis=0, keepdims=True)
            umap_input = (umap_input - umap_mean) / np.maximum(umap_std, 1.0e-6)
        reducer = torchdr.UMAP(
            n_components=3,
            n_neighbors=args.umap_neighbors,
            min_dist=args.umap_min_dist,
            metric="sqeuclidean",
            random_state=42,
            max_iter=args.umap_max_iter,
            device="cpu",
            backend=None,
        )
        with torch.no_grad():
            reduced = reducer.fit_transform(torch.from_numpy(umap_input))
        if isinstance(reduced, torch.Tensor):
            reduced = reduced.detach().cpu().numpy()
        umap_xyz_valid = np.asarray(reduced, dtype=np.float32)
        umap_xyz = np.full((tokens.shape[0], 3), np.nan, dtype=np.float32)
        umap_xyz[valid] = umap_xyz_valid
        umap_xyz_map = umap_xyz.reshape(height, width, 3)
        umap_rgb = np.zeros_like(umap_xyz_map)
        for channel in range(3):
            values = umap_xyz_valid[:, channel]
            lo_umap, hi_umap = np.percentile(values, [1, 99])
            umap_rgb[..., channel] = np.clip(
                (umap_xyz_map[..., channel] - lo_umap) / max(hi_umap - lo_umap, 1e-12),
                0,
                1,
            )
        umap_rgb[~valid.reshape(height, width)] = np.nan

    max_fit = min(50000, z.shape[0])
    rng = np.random.default_rng(42)
    fit_idx = rng.choice(z.shape[0], size=max_fit, replace=False) if z.shape[0] > max_fit else np.arange(z.shape[0])
    clusters = MiniBatchKMeans(n_clusters=args.clusters, random_state=42, n_init=10, batch_size=4096)
    clusters.fit(z[fit_idx])
    labels_valid = clusters.predict(z)
    labels = np.full(tokens.shape[0], np.nan, dtype=np.float32)
    labels[valid] = labels_valid
    label_map = labels.reshape(height, width)

    log_intensity = np.log1p(np.clip(input_field.reshape(-1)[valid], 0, None))
    correlations = [float(np.corrcoef(log_intensity, xyz_valid[:, i])[0, 1]) for i in range(3)]
    centered = z - z.mean(axis=0, keepdims=True)
    singular = np.linalg.svd(centered, full_matrices=False, compute_uv=False)
    variance = singular**2
    probability = variance / variance.sum()
    effective_rank = float(np.exp(-np.sum(probability * np.log(np.clip(probability, 1e-30, None)))))
    valid_y, valid_x = np.where(valid.reshape(height, width))
    valid_intensity = input_field.reshape(-1)[valid]
    cluster_stats = []
    for label in range(args.clusters):
        selected = labels_valid == label
        intensities = valid_intensity[selected]
        cluster_stats.append({
            "label": label,
            "count": int(selected.sum()),
            "fraction": float(selected.mean()),
            "median_intensity": float(np.median(intensities)),
            "intensity_p16_p84": [float(x) for x in np.percentile(intensities, [16, 84])],
            "median_x_y": [float(np.median(valid_x[selected])), float(np.median(valid_y[selected]))],
        })

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    lo, hi = np.nanpercentile(input_field, [1, 99.7])
    base = axes[0, 0].imshow(input_field, origin="lower", cmap="magma", vmin=lo, vmax=hi)
    axes[0, 0].set_title("Input continuum")
    fig.colorbar(base, ax=axes[0, 0], fraction=0.046)
    axes[0, 1].imshow(rgb, origin="lower")
    axes[0, 1].set_title("Dense latent atlas (PCA RGB)")
    clustered = axes[1, 0].imshow(label_map, origin="lower", cmap="tab10", vmin=-0.5, vmax=args.clusters - 0.5)
    axes[1, 0].set_title(f"Latent morphology groups (k={args.clusters})")
    fig.colorbar(clustered, ax=axes[1, 0], ticks=range(args.clusters), fraction=0.046)
    scatter_idx = rng.choice(z.shape[0], size=min(30000, z.shape[0]), replace=False)
    axes[1, 1].scatter(xyz_valid[scatter_idx, 0], xyz_valid[scatter_idx, 1], c=labels_valid[scatter_idx], cmap="tab10", s=1, alpha=0.35)
    axes[1, 1].set(title="Latent topology", xlabel="PC1", ylabel="PC2")
    for ax in axes.flat[:3]:
        ax.set(xlabel="x pixel", ylabel="y pixel")
    if args.focus_valid_bbox:
        pad = 3
        xmin = max(int(valid_x.min()) - pad, 0)
        xmax = min(int(valid_x.max()) + pad + 1, width)
        ymin = max(int(valid_y.min()) - pad, 0)
        ymax = min(int(valid_y.max()) + pad + 1, height)
        for ax in axes.flat[:3]:
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
    fig.suptitle("G316 SA-JEPA continuum pilot", fontsize=16)
    figure_name = "latent_atlas_pilot.png" if not suffix else f"latent_atlas{suffix}.png"
    fig.savefig(figure_dir / figure_name, dpi=180)
    plt.close(fig)

    if umap_xyz_valid is not None and umap_rgb is not None:
        umap_fig, umap_axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
        umap_base = umap_axes[0, 0].imshow(input_field, origin="lower", cmap="magma", vmin=lo, vmax=hi)
        umap_axes[0, 0].set_title("Input continuum")
        umap_fig.colorbar(umap_base, ax=umap_axes[0, 0], fraction=0.046)
        umap_axes[0, 1].imshow(umap_rgb, origin="lower")
        umap_axes[0, 1].set_title("Dense latent atlas (UMAP RGB)")
        umap_axes[1, 0].imshow(rgb, origin="lower")
        umap_axes[1, 0].set_title("Dense latent atlas (PCA RGB)")
        umap_scatter_idx = rng.choice(z.shape[0], size=min(30000, z.shape[0]), replace=False)
        umap_points = umap_axes[1, 1].scatter(
            umap_xyz_valid[umap_scatter_idx, 0],
            umap_xyz_valid[umap_scatter_idx, 1],
            c=log_intensity[umap_scatter_idx],
            cmap="magma",
            s=1,
            alpha=0.35,
        )
        umap_axes[1, 1].set(title="UMAP topology (colour = log intensity)", xlabel="UMAP-1", ylabel="UMAP-2")
        umap_fig.colorbar(umap_points, ax=umap_axes[1, 1], fraction=0.046)
        for ax in (umap_axes[0, 0], umap_axes[0, 1], umap_axes[1, 0]):
            ax.set(xlabel="x pixel", ylabel="y pixel")
        if args.focus_valid_bbox:
            for ax in (umap_axes[0, 0], umap_axes[0, 1], umap_axes[1, 0]):
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)
        umap_fig.suptitle("G316 SA-JEPA continuum pilot: UMAP test", fontsize=16)
        umap_fig.savefig(figure_dir / f"umap_atlas{suffix}.png", dpi=180)
        plt.close(umap_fig)

    np.save(out_root / f"pca_xyz{suffix}.npy", xyz_map.astype(np.float32))
    if umap_xyz_map is not None:
        np.save(out_root / f"umap_xyz{suffix}.npy", umap_xyz_map.astype(np.float32))
    np.save(out_root / f"cluster_labels{suffix}.npy", label_map.astype(np.float32))
    summary = {
        "session_dir": str(session_dir),
        "latent_source_key": latent_key,
        "latent_shape_chw": list(latent.shape),
        "valid_tokens": int(z.shape[0]),
        "effective_rank": effective_rank,
        "pca_explained_variance_ratio": [float(x) for x in pca.explained_variance_ratio_],
        "correlation_log_intensity_with_pc1_pc2_pc3": correlations,
        "umap": (
            {
                "method": "TorchDR UMAP",
                "n_neighbors": args.umap_neighbors,
                "min_dist": args.umap_min_dist,
                "max_iter": args.umap_max_iter,
                "standardize": args.umap_standardize,
                "correlation_log_intensity_with_umap1_umap2_umap3": [
                    float(np.corrcoef(log_intensity, umap_xyz_valid[:, i])[0, 1])
                    for i in range(3)
                ],
            }
            if umap_xyz_valid is not None
            else None
        ),
        "clusters": args.clusters,
        "cluster_stats": cluster_stats,
        "valid_mask_erosion_pixels": args.valid_erosion,
        "figure_focus_valid_bbox": args.focus_valid_bbox,
        "interpretation_warning": "Pilot clusters are unsupervised coordinates, not physical classes or a catalogue.",
    }
    (out_root / f"pilot_summary{suffix}.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "summarize"))
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--session-dir", type=Path, default=DEFAULT_ROOT / "sessions/lancet_continuum_cpu_pilot")
    parser.add_argument("--max-dimension", type=int, default=768)
    parser.add_argument("--minimum-coverage", type=float, default=0.95)
    parser.add_argument("--distance-kpc", type=float, default=2.8)
    parser.add_argument("--center-ra", help="ICRS right ascension, interpreted as hour angle")
    parser.add_argument("--center-dec", help="ICRS declination in degrees")
    parser.add_argument("--native-crop-size", type=int, help="Square native-FITS crop size in pixels")
    parser.add_argument("--prepared-crop", type=int, nargs=4, metavar=("Y0", "Y1", "X0", "X1"))
    parser.add_argument("--target-mask-erosion", type=int, default=0)
    parser.add_argument("--snr-min", type=float, help="Minimum excess S/N for eligible target pixels")
    parser.add_argument("--noise-rms", type=float, help="Adopted RMS in the FITS brightness units")
    parser.add_argument("--background-percentile", type=float, default=10.0)
    parser.add_argument("--clusters", type=int, default=6)
    parser.add_argument("--summary-name", default="pilot")
    parser.add_argument("--valid-erosion", type=int, default=0)
    parser.add_argument("--focus-valid-bbox", action="store_true")
    parser.add_argument("--umap", action="store_true")
    parser.add_argument("--umap-neighbors", type=int, default=50)
    parser.add_argument("--umap-min-dist", type=float, default=0.2)
    parser.add_argument("--umap-max-iter", type=int, default=500)
    parser.add_argument("--umap-standardize", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


if __name__ == "__main__":
    namespace = parse_args()
    if namespace.stage == "prepare":
        prepare(namespace)
    else:
        summarize(namespace)
