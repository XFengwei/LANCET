#!/usr/bin/env python
"""Inspect a stellar-cluster orbit near the LANCET filament."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from lancet import ProjectPaths


DEFAULT_CATALOGUE = (
    "data/delivered/starcluster/orbit_onlyback_OC100Myr_step01.fits"
)
DEFAULT_CENTER_RA = "14:45:13.199"
DEFAULT_CENTER_DEC = "-59:48:48.89"
REQUIRED_COLUMNS = {"ID", "time", "l", "b", "ra", "dec"}


def angle_wrap(angle_deg: np.ndarray | float) -> np.ndarray:
    """Wrap angles in degrees to the interval ``[-180, 180)``."""
    return (np.asarray(angle_deg, dtype=float) + 180.0) % 360.0 - 180.0


def load_catalogue(path: str | Path) -> Table:
    """Read and validate the collaborator's orbit catalogue."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Star-cluster catalogue not found: {path}")
    catalogue = Table.read(path, memmap=True)
    missing = REQUIRED_COLUMNS - set(catalogue.colnames)
    if missing:
        raise ValueError(
            "Star-cluster catalogue is missing required columns: "
            + ", ".join(sorted(missing))
        )
    return catalogue


def select_orbit(
    catalogue: Table,
    cluster_id: float,
    lookback_myr: float,
) -> Table:
    """Select one cluster over the requested interval before the present."""
    if not np.isfinite(cluster_id):
        raise ValueError("cluster_id must be finite")
    if not np.isfinite(lookback_myr) or lookback_myr <= 0:
        raise ValueError("lookback_myr must be positive and finite")

    time = np.asarray(catalogue["time"], dtype=float)
    mask = (
        np.asarray(catalogue["ID"], dtype=float) == float(cluster_id)
    ) & (time > -float(lookback_myr)) & (time <= 0.0)
    orbit = catalogue[mask].copy()
    if not len(orbit):
        raise ValueError(
            f"No rows found for cluster ID {cluster_id:g} over the last "
            f"{lookback_myr:g} Myr"
        )
    orbit.sort("time")
    return orbit


def filament_center(
    ra: str = DEFAULT_CENTER_RA,
    dec: str = DEFAULT_CENTER_DEC,
) -> SkyCoord:
    """Return the LANCET center from sexagesimal ICRS coordinates."""
    return SkyCoord(
        ra,
        dec,
        unit=(u.hourangle, u.deg),
        frame="icrs",
    )


def orbit_separations(orbit: Table, center: SkyCoord) -> u.Quantity:
    """Calculate angular separation of every orbit sample from the center."""
    coordinates = SkyCoord(
        ra=np.asarray(orbit["ra"], dtype=float) * u.deg,
        dec=np.asarray(orbit["dec"], dtype=float) * u.deg,
        frame="icrs",
    )
    return center.separation(coordinates)


def close_approaches(
    orbit: Table,
    center: SkyCoord,
    radius: u.Quantity,
) -> Table:
    """Return orbit samples inside a cone, with separation in arcminutes."""
    if radius <= 0 * u.arcmin:
        raise ValueError("search radius must be positive")
    separations = orbit_separations(orbit, center)
    selected = orbit[separations < radius].copy()
    selected["sep_arcmin"] = separations[separations < radius].to_value(
        u.arcmin
    )
    if len(selected):
        selected.sort("sep_arcmin")
    return selected


def _prepare_output(output: Path, overwrite: bool) -> None:
    if output.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {output}; pass --overwrite to replace it"
        )
    output.parent.mkdir(parents=True, exist_ok=True)


def plot_whole_sky(
    orbit: Table,
    center: SkyCoord,
    output: Path,
    *,
    cluster_id: float,
    overwrite: bool = False,
) -> Path:
    """Plot the selected orbit in Galactic coordinates."""
    _prepare_output(output, overwrite)
    longitude_plot = -angle_wrap(np.asarray(orbit["l"], dtype=float))
    latitude = np.asarray(orbit["b"], dtype=float)
    center_galactic = center.galactic
    center_longitude_plot = -angle_wrap(center_galactic.l.to_value(u.deg))

    fig = plt.figure(figsize=(10, 5))
    axis = fig.add_subplot(111, projection="mollweide")
    axis.scatter(
        np.deg2rad(longitude_plot),
        np.deg2rad(latitude),
        s=10,
        alpha=0.6,
        label=f"Cluster {cluster_id:g} orbit",
    )
    axis.scatter(
        np.deg2rad(center_longitude_plot),
        np.deg2rad(center_galactic.b.to_value(u.deg)),
        s=80,
        marker="*",
        color="tab:red",
        label="LANCET filament center",
    )
    longitude_ticks = np.arange(-150, 180, 30)
    longitude_labels = (-longitude_ticks) % 360
    axis.set_xticks(np.deg2rad(longitude_ticks))
    axis.set_xticklabels(
        [rf"${value:g}^\circ$" for value in longitude_labels]
    )
    axis.grid(True, alpha=0.4)
    axis.set_xlabel("Galactic longitude $l$")
    axis.set_ylabel("Galactic latitude $b$")
    axis.set_title("Selected cluster orbit in Galactic coordinates")
    axis.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_zoom(
    orbit: Table,
    center: SkyCoord,
    output: Path,
    *,
    cluster_id: float,
    radius_arcmin: float,
    zoom_deg: float,
    overwrite: bool = False,
) -> Path:
    """Plot the selected orbit around the LANCET center."""
    _prepare_output(output, overwrite)
    center_galactic = center.galactic
    center_l = center_galactic.l.to_value(u.deg)
    center_b = center_galactic.b.to_value(u.deg)
    lookback = -np.asarray(orbit["time"], dtype=float)

    fig, axis = plt.subplots(figsize=(7, 6))
    points = axis.scatter(
        np.asarray(orbit["l"], dtype=float),
        np.asarray(orbit["b"], dtype=float),
        marker="o",
        s=16,
        alpha=0.9,
        c=lookback,
        cmap="inferno",
        label=f"Cluster {cluster_id:g} orbit",
    )
    colorbar = fig.colorbar(points, ax=axis)
    colorbar.set_label("Lookback time [Myr]")
    axis.scatter(
        center_l,
        center_b,
        s=120,
        marker="*",
        color="red",
        label="LANCET filament center",
    )
    axis.add_patch(
        Circle(
            (center_l, center_b),
            radius=radius_arcmin / 60.0,
            edgecolor="red",
            facecolor="none",
            linewidth=2,
            linestyle="--",
            label=f"{radius_arcmin:g} arcmin radius",
        )
    )
    axis.set_xlim(center_l + zoom_deg, center_l - zoom_deg)
    axis.set_ylim(center_b - zoom_deg, center_b + zoom_deg)
    axis.set_xlabel("Galactic longitude $l$ [deg]")
    axis.set_ylabel("Galactic latitude $b$ [deg]")
    axis.set_title("Selected cluster orbit around the LANCET filament")
    axis.grid(True, alpha=0.4)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=Path(DEFAULT_CATALOGUE))
    parser.add_argument("--cluster-id", type=float, default=3944.0)
    parser.add_argument("--lookback-myr", type=float, default=5.0)
    parser.add_argument("--center-ra", default=DEFAULT_CENTER_RA)
    parser.add_argument("--center-dec", default=DEFAULT_CENTER_DEC)
    parser.add_argument("--search-radius-arcmin", type=float, default=10.0)
    parser.add_argument("--zoom-deg", type=float, default=3.0)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/starcluster"),
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Print the selection summary without creating figures",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    project = ProjectPaths.discover()
    catalogue_path = args.catalogue.expanduser()
    if not catalogue_path.is_absolute():
        catalogue_path = project.root / catalogue_path

    try:
        catalogue = load_catalogue(catalogue_path)
        orbit = select_orbit(
            catalogue,
            cluster_id=args.cluster_id,
            lookback_myr=args.lookback_myr,
        )
        center = filament_center(args.center_ra, args.center_dec)
        radius = args.search_radius_arcmin * u.arcmin
        nearby = close_approaches(orbit, center, radius)
    except (FileNotFoundError, ValueError) as error:
        print(error)
        return 2

    separations = orbit_separations(orbit, center)
    print(f"Catalogue rows: {len(catalogue)}")
    print(f"Selected cluster ID: {args.cluster_id:g}")
    print(f"Selected orbit samples: {len(orbit)}")
    print(
        "Time range: "
        f"{np.min(orbit['time']):.1f} to {np.max(orbit['time']):.1f} Myr"
    )
    print(
        "LANCET center: "
        f"l={center.galactic.l.deg:.6f} deg, "
        f"b={center.galactic.b.deg:.6f} deg"
    )
    print(f"Closest approach: {separations.min().arcmin:.3f} arcmin")
    print(
        f"Samples within {args.search_radius_arcmin:g} arcmin: {len(nearby)}"
    )

    if args.summary_only:
        return 0

    try:
        output_dir = project.validate_output(args.output_dir)
        whole_sky = output_dir / "Cluster_whole_sky_galactic.png"
        zoom = output_dir / "Cluster_LANCET_zoomin.png"
        plot_whole_sky(
            orbit,
            center,
            whole_sky,
            cluster_id=args.cluster_id,
            overwrite=args.overwrite,
        )
        plot_zoom(
            orbit,
            center,
            zoom,
            cluster_id=args.cluster_id,
            radius_arcmin=args.search_radius_arcmin,
            zoom_deg=args.zoom_deg,
            overwrite=args.overwrite,
        )
    except (FileExistsError, ValueError) as error:
        print(error)
        return 2

    print(f"Wrote {whole_sky.relative_to(project.root)}")
    print(f"Wrote {zoom.relative_to(project.root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
