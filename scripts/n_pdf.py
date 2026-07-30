#!/usr/bin/env python
"""Build provisional regional column-density PDFs."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import json
from pathlib import Path
import tomllib

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table, vstack
from astropy.wcs import WCS
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from lancet import ProjectPaths
from lancet.npdf import (
    logarithmic_column_density_contrast,
    npdf_histogram,
    positive_column_density_excess,
    shared_excess_mean,
)


DEFAULT_CONFIG = Path("tasks/n_pdf/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
    """Load task settings from TOML."""
    if not path.is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _mapping(
    config: Mapping[str, object],
    name: str,
) -> Mapping[str, object]:
    value = config[name]
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a TOML table")
    return value


def _project_path(project: ProjectPaths, value: object) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = project.root / path
    return path.resolve()


def resolve_task_paths(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> dict[str, Path]:
    """Resolve the delivered column-density input and safe result root."""
    inputs = _mapping(config, "inputs")
    outputs = _mapping(config, "outputs")
    paths = {
        name: _project_path(project, value)
        for name, value in inputs.items()
    }
    paths["output_root"] = project.validate_output(str(outputs["root"]))
    for name in inputs:
        if not paths[name].is_file():
            raise FileNotFoundError(f"{name} not found: {paths[name]}")
        if not paths[name].is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input")
    return paths


def _load_image(path: Path) -> tuple[np.ndarray, fits.Header, WCS]:
    data, header = fits.getdata(path, header=True, memmap=True)
    array = np.asarray(data)
    if array.ndim < 2:
        raise ValueError("Column-density FITS data must be at least 2D")
    if array.ndim > 2:
        array = array[
            tuple(0 for _ in range(array.ndim - 2))
            + (slice(None), slice(None))
        ]
    array = np.asarray(array, dtype=float)
    if array.ndim != 2:
        raise ValueError("Column-density input could not be reduced to 2D")
    return array, header, WCS(header).celestial


def _region_definitions(
    config: Mapping[str, object],
) -> list[tuple[str, Mapping[str, object]]]:
    regions = _mapping(config, "regions")
    ordering = _mapping(config, "ordering")
    definitions: list[tuple[str, Mapping[str, object]]] = []
    for key in ordering["region_order"]:
        value = regions[str(key)]
        if not isinstance(value, Mapping):
            raise TypeError(f"regions.{key} must be a TOML table")
        definitions.append((str(key), value))
    return definitions


def _make_cutouts(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> tuple[
    dict[str, np.ndarray],
    dict[str, fits.Header],
    fits.Header,
]:
    image, source_header, wcs = _load_image(paths["column_density"])
    cutouts: dict[str, np.ndarray] = {}
    headers: dict[str, fits.Header] = {}
    for key, region in _region_definitions(config):
        center = SkyCoord(
            float(region["center_ra_deg"]) * u.deg,
            float(region["center_dec_deg"]) * u.deg,
            frame="icrs",
        )
        cutout = Cutout2D(
            image,
            center,
            size=(
                float(region["height_arcsec"]) * u.arcsec,
                float(region["width_arcsec"]) * u.arcsec,
            ),
            wcs=wcs,
            mode="partial",
            fill_value=np.nan,
        )
        cutouts[key] = np.asarray(cutout.data, dtype=float)
        header = cutout.wcs.to_header()
        for name in (
            "BUNIT",
            "BTYPE",
            "OBJECT",
            "TELESCOP",
            "INSTRUME",
            "DATE-OBS",
        ):
            if name in source_header:
                header[name] = source_header[name]
        header["BTYPE"] = "Column density cutout"
        header.add_history(
            "LANCET n_pdf ongoing workflow; assumptions are provisional"
        )
        headers[key] = header
    return cutouts, headers, source_header


def _analysis_settings(
    config: Mapping[str, object],
) -> tuple[float, np.ndarray]:
    analysis = _mapping(config, "analysis")
    background = float(analysis["background_cm2"])
    if not np.isfinite(background) or background < 0:
        raise ValueError("analysis.background_cm2 must be non-negative")
    edges = np.linspace(
        float(analysis["eta_min"]),
        float(analysis["eta_max"]),
        int(analysis["eta_bins"]) + 1,
    )
    return background, edges


def _task_status(config: Mapping[str, object]) -> str:
    return str(_mapping(config, "task")["status"])


def _cutout_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> dict[str, Path]:
    return {
        key: (
            paths["output_root"]
            / "cutouts"
            / f"{key}_column_density.fits"
        )
        for key, _ in _region_definitions(config)
    }


def _table_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "regional_npdf.ecsv"


def _summary_path(paths: Mapping[str, Path]) -> Path:
    return paths["output_root"] / "tables" / "regional_npdf_summary.json"


def _figure_paths(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> list[Path]:
    plot = _mapping(config, "plot")
    formats = [str(value).lower() for value in plot["formats"]]
    unsupported = set(formats) - {"pdf", "png"}
    if unsupported:
        raise ValueError(
            "Unsupported figure formats: " + ", ".join(sorted(unsupported))
        )
    return [
        paths["output_root"]
        / "figures"
        / f"regional_npdf.{extension}"
        for extension in formats
    ]


def _prepare_outputs(
    outputs: Sequence[Path],
    *,
    overwrite: bool,
) -> None:
    for output in outputs:
        if output.exists() and not overwrite:
            raise FileExistsError(
                f"Output already exists: {output}; "
                "pass --overwrite to replace it"
            )
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)


def _analyse(
    cutouts: Mapping[str, np.ndarray],
    config: Mapping[str, object],
) -> tuple[Table, dict[str, object]]:
    background, edges = _analysis_settings(config)
    normalization = shared_excess_mean(cutouts, background)
    tables: list[Table] = []
    region_summaries: dict[str, object] = {}
    definitions = dict(_region_definitions(config))
    for key, values in cutouts.items():
        region = definitions[key]
        eta = logarithmic_column_density_contrast(
            values,
            background=background,
            normalization=normalization,
        )
        histogram = npdf_histogram(eta, edges)
        count = len(np.asarray(histogram["centers"]))
        tables.append(
            Table(
                {
                    "region": np.full(count, key),
                    "label": np.full(count, str(region["label"])),
                    "eta_center": histogram["centers"],
                    "eta_left": histogram["edges"][:-1],
                    "eta_right": histogram["edges"][1:],
                    "pixel_count": histogram["counts"],
                    "pdf": histogram["density"],
                    "pdf_error": histogram["density_error"],
                }
            )
        )
        positive = positive_column_density_excess(values, background)
        region_summaries[key] = {
            "label": str(region["label"]),
            "shape": list(values.shape),
            "finite_pixel_count": int(np.isfinite(values).sum()),
            "positive_excess_pixel_count": int(positive.size),
            "histogram_pixel_count": int(
                np.sum(histogram["counts"])
            ),
            "mean_positive_excess_cm2": float(np.mean(positive)),
            "median_positive_excess_cm2": float(np.median(positive)),
        }
    table = vstack(tables, metadata_conflicts="silent")
    table.meta["status"] = _task_status(config)
    table.meta["paper_ready"] = False
    table.meta["background_cm2"] = background
    table.meta["normalization_cm2"] = normalization
    table.meta["eta_definition"] = (
        "ln((N_H2 - background) / shared positive-excess mean)"
    )
    table.meta["normalization_domain"] = (
        "finite positive background-subtracted pixels in all regions"
    )
    summary: dict[str, object] = {
        "status": _task_status(config),
        "paper_ready": False,
        "input_role": (
            "provisional retained 6+18+36 arcsec all-scale column-density map"
        ),
        "background_cm2": background,
        "normalization_cm2": normalization,
        "normalization_domain": (
            "finite positive background-subtracted pixels across all regions"
        ),
        "regions": region_summaries,
        "continuity_note": (
            "Revisit map selection, background model, masks, correlated "
            "pixels, completeness, and distribution fitting before "
            "scientific interpretation."
        ),
    }
    return table, summary


def inspect_inputs(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
) -> None:
    """Report cutout coverage and provisional N-PDF assumptions."""
    cutouts, _, source_header = _make_cutouts(paths, config)
    _, summary = _analyse(cutouts, config)
    print("N-PDF workflow:")
    print(f"  status: {_task_status(config).upper()} / NOT PAPER-READY")
    print(f"  column_density: {paths['column_density']}")
    print(
        f"  BUNIT={source_header.get('BUNIT', 'not recorded')}; "
        "adopted as the provisional retained 6+18+36 arcsec all-scale map"
    )
    regions = summary["regions"]
    assert isinstance(regions, Mapping)
    for key, region in _region_definitions(config):
        values = regions[key]
        assert isinstance(values, Mapping)
        print(
            f"  {region['label']}: shape={tuple(values['shape'])}, "
            f"positive pixels after background="
            f"{values['positive_excess_pixel_count']}"
        )
    print(
        f"Background: {summary['background_cm2']:.6g} cm^-2; "
        f"shared positive-excess mean: "
        f"{summary['normalization_cm2']:.6g} cm^-2"
    )
    print(
        "Continuity note: keep this task ongoing until map choice, "
        "background/mask policy, correlated pixels, completeness, and "
        "distribution fitting are established."
    )
    print(f"Output root: {paths['output_root']}")


def build_products(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> tuple[Table, dict[str, object]]:
    """Write regional cutouts, normalized histograms, and assumptions."""
    cutouts, headers, _ = _make_cutouts(paths, config)
    table, summary = _analyse(cutouts, config)
    cutout_outputs = _cutout_paths(paths, config)
    outputs = [
        *cutout_outputs.values(),
        _table_path(paths),
        _summary_path(paths),
    ]
    _prepare_outputs(outputs, overwrite=overwrite)
    for key, output in cutout_outputs.items():
        fits.writeto(
            output,
            np.asarray(cutouts[key], dtype=np.float32),
            headers[key],
            overwrite=True,
        )
        print(f"Wrote {output}")
    table.write(
        _table_path(paths),
        format="ascii.ecsv",
        overwrite=True,
    )
    _summary_path(paths).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {_table_path(paths)}")
    print(f"Wrote {_summary_path(paths)}")
    return table, summary


def plot_npdf(
    paths: Mapping[str, Path],
    config: Mapping[str, object],
    *,
    overwrite: bool,
) -> list[Path]:
    """Plot the current regional normalized column-density PDFs."""
    table_path = _table_path(paths)
    summary_path = _summary_path(paths)
    if not table_path.is_file() or not summary_path.is_file():
        raise FileNotFoundError(
            "N-PDF tables are unavailable; run --stage build first"
        )
    table = Table.read(table_path, format="ascii.ecsv")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    outputs = _figure_paths(paths, config)
    _prepare_outputs(outputs, overwrite=overwrite)
    plot = _mapping(config, "plot")
    figure, axis = plt.subplots(
        figsize=(
            float(plot["width_inches"]),
            float(plot["height_inches"]),
        )
    )
    for key, region in _region_definitions(config):
        mask = np.asarray(table["region"]).astype(str) == key
        subset = table[mask]
        positive = np.asarray(subset["pixel_count"], dtype=int) > 0
        centers = np.asarray(subset["eta_center"], dtype=float)
        density = np.asarray(subset["pdf"], dtype=float)
        error = np.asarray(subset["pdf_error"], dtype=float)
        displayed_density = np.where(positive, density, np.nan)
        axis.step(
            centers,
            displayed_density,
            where="mid",
            color=str(region["color"]),
            linewidth=1.5,
            label=str(region["label"]),
        )
        axis.errorbar(
            centers[positive],
            density[positive],
            yerr=error[positive],
            fmt=".",
            markersize=4,
            capsize=2,
            color=str(region["color"]),
        )
    axis.text(
        0.02,
        0.98,
        (
            "Ongoing diagnostic\n"
            r"$N_0=\langle N_{\rm H_2}-N_{\rm bkg}\rangle_+=$"
            f"{float(summary['normalization_cm2']):.2e} "
            r"cm$^{-2}$"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )
    axis.set_xlabel(
        r"$\eta=\ln[(N_{\rm H_2}-N_{\rm bkg})/N_0]$"
    )
    axis.set_ylabel("Probability density")
    axis.set_xlim(float(plot["eta_min"]), float(plot["eta_max"]))
    axis.set_yscale("log")
    axis.set_ylim(float(plot["pdf_min"]), float(plot["pdf_max"]))
    axis.grid(alpha=0.15, which="both")
    axis.legend(frameon=False, loc=str(plot["legend_location"]))
    figure.tight_layout()
    for output in outputs:
        figure.savefig(
            output,
            dpi=int(plot["dpi"]),
            bbox_inches="tight",
        )
        print(f"Wrote {output}")
    plt.close(figure)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "build", "plot", "all"),
        default="inspect",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    project = ProjectPaths.discover()
    config = load_task_config(_project_path(project, args.config))
    paths = resolve_task_paths(project, config)
    if args.stage == "inspect":
        inspect_inputs(paths, config)
    elif args.stage == "build":
        build_products(paths, config, overwrite=args.overwrite)
    elif args.stage == "plot":
        plot_npdf(paths, config, overwrite=args.overwrite)
    else:
        build_products(paths, config, overwrite=args.overwrite)
        plot_npdf(paths, config, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
