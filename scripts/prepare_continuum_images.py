#!/usr/bin/env python
"""Mask primary-beam edges and convert continuum-image brightness units."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
import tomllib

from astropy.io import fits
import numpy as np

from lancet import ProjectPaths
from lancet.imaging import (
    apply_primary_beam_mask,
    jybeam_to_mjysr_from_header,
)


DEFAULT_CONFIG = Path("tasks/prepare_continuum_images/task.toml")


def load_task_config(path: Path) -> dict[str, object]:
    """Load task settings from TOML."""
    if not path.is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _project_path(project: ProjectPaths, value: object) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = project.root / path
    return path.resolve()


def _output_path(output_root: Path, value: object) -> Path:
    relative = Path(str(value))
    if relative.is_absolute():
        raise ValueError(f"Task output must be relative: {relative}")
    path = (output_root / relative).resolve()
    if not path.is_relative_to(output_root):
        raise ValueError(f"Task output escapes its root: {path}")
    return path


def resolve_task(
    project: ProjectPaths,
    config: Mapping[str, object],
) -> tuple[Path, list[dict[str, object]], list[dict[str, object]]]:
    """Resolve and validate task inputs and outputs."""
    outputs = config["outputs"]
    if not isinstance(outputs, Mapping):
        raise TypeError("Task outputs must be a TOML table")
    output_root = project.validate_output(str(outputs["root"]))

    masks: list[dict[str, object]] = []
    for raw_entry in config["primary_beam_masks"]:
        if not isinstance(raw_entry, Mapping):
            raise TypeError("Primary-beam mask entries must be TOML tables")
        entry = dict(raw_entry)
        for key in ("image", "primary_beam"):
            path = _project_path(project, entry[key])
            if not path.is_relative_to(project.delivered.resolve()):
                raise ValueError(f"{entry['name']} {key} must be delivered")
            if not path.is_file():
                raise FileNotFoundError(f"{entry['name']} {key} missing: {path}")
            entry[key] = path
        entry["output"] = _output_path(output_root, entry["output"])
        masks.append(entry)

    conversions: list[dict[str, object]] = []
    for raw_entry in config["unit_conversions"]:
        if not isinstance(raw_entry, Mapping):
            raise TypeError("Unit-conversion entries must be TOML tables")
        entry = dict(raw_entry)
        source_kind = str(entry["source_kind"])
        if source_kind == "delivered":
            input_path = _project_path(project, entry["input"])
            if not input_path.is_relative_to(project.delivered.resolve()):
                raise ValueError(
                    f"{entry['name']} input must be under delivered data"
                )
            if not input_path.is_file():
                raise FileNotFoundError(
                    f"{entry['name']} input missing: {input_path}"
                )
        elif source_kind == "task_output":
            input_path = _output_path(output_root, entry["input"])
        else:
            raise ValueError(
                f"{entry['name']} has unsupported source_kind={source_kind}"
            )
        entry["input"] = input_path
        entry["output"] = _output_path(output_root, entry["output"])
        conversions.append(entry)
    return output_root, masks, conversions


def _prepare_output(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}; pass --overwrite to replace it"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def inspect_task(
    output_root: Path,
    masks: Sequence[Mapping[str, object]],
    conversions: Sequence[Mapping[str, object]],
    config: Mapping[str, object],
) -> None:
    """Report declared operations and dependency readiness."""
    mask_settings = config["primary_beam_mask"]
    if not isinstance(mask_settings, Mapping):
        raise TypeError("Primary-beam settings must be a TOML table")
    print(f"Output root: {output_root}")
    print(
        "Primary-beam masks: "
        f"{len(masks)} at response >= "
        f"{float(mask_settings['minimum_response']):g}"
    )
    for entry in masks:
        print(f"  {entry['name']}: {entry['image']}")
        print(f"    primary beam: {entry['primary_beam']}")
        print(f"    output: {entry['output']}")

    print(f"Jy/beam to MJy/sr conversions: {len(conversions)}")
    for entry in conversions:
        status = "ready" if Path(entry["input"]).is_file() else "needs mask stage"
        print(
            f"  {entry['name']}: {entry['source_kind']} ({status})"
            f"\n    input: {entry['input']}"
            f"\n    output: {entry['output']}"
        )


def run_primary_beam_masks(
    masks: Sequence[Mapping[str, object]],
    config: Mapping[str, object],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Apply the configured primary-beam cutoff masks."""
    settings = config["primary_beam_mask"]
    if not isinstance(settings, Mapping):
        raise TypeError("Primary-beam settings must be a TOML table")
    threshold = float(settings["minimum_response"])
    written: list[Path] = []
    for entry in masks:
        output = Path(entry["output"])
        _prepare_output(output, overwrite=overwrite)
        image, header = fits.getdata(
            Path(entry["image"]),
            header=True,
        )
        primary_beam = fits.getdata(Path(entry["primary_beam"]))
        masked = apply_primary_beam_mask(
            image,
            primary_beam,
            minimum_response=threshold,
        )
        if "BLANK" in header and int(header["BITPIX"]) < 0:
            del header["BLANK"]
        header["PBTHRESH"] = (threshold, "Minimum retained PB response")
        header.add_history(
            "LANCET prepare_continuum_images: applied primary-beam mask"
        )
        fits.writeto(
            output,
            masked.astype(np.float32),
            header,
            overwrite=overwrite,
        )
        written.append(output)
        print(f"Wrote primary-beam masked image: {output}")
    return written


def run_unit_conversions(
    conversions: Sequence[Mapping[str, object]],
    *,
    overwrite: bool = False,
) -> list[Path]:
    """Convert configured images from Jy/beam to MJy/sr."""
    written: list[Path] = []
    for entry in conversions:
        input_path = Path(entry["input"])
        if not input_path.is_file():
            raise FileNotFoundError(
                f"Conversion input missing: {input_path}; run --stage mask first"
            )
        output = Path(entry["output"])
        _prepare_output(output, overwrite=overwrite)
        image, header = fits.getdata(input_path, header=True)
        converted = jybeam_to_mjysr_from_header(image, header)
        if "BLANK" in header and int(header["BITPIX"]) < 0:
            del header["BLANK"]
        header["BUNIT"] = "MJy/sr"
        header.add_history(
            "LANCET prepare_continuum_images: converted Jy/beam to MJy/sr"
        )
        fits.writeto(
            output,
            np.asarray(converted, dtype=np.float32),
            header,
            overwrite=overwrite,
        )
        written.append(output)
        print(f"Wrote brightness-unit conversion: {output}")
    return written


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("inspect", "mask", "convert", "all"),
        default="inspect",
    )
    parser.add_argument(
        "--only",
        action="append",
        metavar="NAME",
        help="Run only the named operation; repeat to select multiple names",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    project = ProjectPaths.discover()
    config_path = args.config
    if not config_path.is_absolute():
        config_path = project.root / config_path
    config = load_task_config(config_path.resolve())
    output_root, masks, conversions = resolve_task(project, config)
    if args.only:
        selected = set(args.only)
        available = {
            str(entry["name"]) for entry in (*masks, *conversions)
        }
        unknown = selected - available
        if unknown:
            raise ValueError(
                "Unknown operation names: " + ", ".join(sorted(unknown))
            )
        masks = [
            entry for entry in masks if str(entry["name"]) in selected
        ]
        conversions = [
            entry
            for entry in conversions
            if str(entry["name"]) in selected
        ]

    if args.stage in {"inspect", "all"}:
        inspect_task(output_root, masks, conversions, config)
    if args.stage in {"mask", "all"}:
        run_primary_beam_masks(
            masks,
            config,
            overwrite=args.overwrite,
        )
    if args.stage in {"convert", "all"}:
        run_unit_conversions(
            conversions,
            overwrite=args.overwrite,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
