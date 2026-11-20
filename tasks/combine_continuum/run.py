#!/usr/bin/env python
"""Combine the 12m+7m continuum image with the TP models, and test the combination.

Combination stages: ``inspect``, ``prepare``, ``tp_factor``, ``box_spectrum``,
``planck_calibration``, ``jcomb``, ``feather``. Diagnostic stages: ``uv``, ``overlap``,
``full_sampling``, ``amplitude_alignment``, ``figure2``, ``compare_images``, ``dirty_beam``,
``aca_forward_test``, ``joint_images``, ``pointing_coverage``. ``combine``, ``diagnostics`` and
``all`` select groups of stages. ``inspect`` writes nothing. Diagnostic stages write under
``diagnostics/`` in the output root.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
import json
from pathlib import Path
import tomllib

from astropy.io import fits

from lancet import ProjectPaths
from lancet import combination_diagnostics as diagnostics
from lancet import combination_stages
from lancet.aca_forward_test import run_aca_forward_test
from lancet.combination_figure import make_combination_figure
from lancet.combination_stages import input_hashes, write_json
from lancet.combination_tests import require_inputs, run_overlap_tests
from lancet.dirty_beam import run_dirty_beam
from lancet.full_sampling import build_full_sampling_profiles
from lancet.imaging import frequency_label, read_frequency_ghz
from lancet.joint_image_comparison import run_joint_image_comparison

DEFAULT_CONFIG = Path("tasks/combine_continuum/task.toml")
DELIVERED_INPUTS = ("interferometer", "gain", "psf", "psf_imfit_log", "aca_pbcor", "twelve_m_pbcor",
                    "casa_image", "casa_pb", "planck_217", "planck_co", "planck_rimo",
                    "planck_fetch_record", "atlasgal_870")
# TP models are the products of the synthesize_continuum task.
SYNTHESIS_INPUTS = ("tp_dust", "tp_dust_ff")
SYNTHESIS_RESULTS = Path("results/synthesize_continuum")
# box_spectrum reads the J-comb product as its mock truth image, so it runs after jcomb.
COMBINE_STAGES = ("inspect", "prepare", "tp_factor", "jcomb", "feather", "box_spectrum",
                  "planck_calibration")
DIAGNOSTIC_STAGES = ("uv", "overlap", "full_sampling", "amplitude_alignment", "figure2",
                     "compare_images", "dirty_beam", "aca_forward_test", "joint_images",
                     "pointing_coverage")
DEFAULT_LABEL = "224GHz"
IMAGE_INPUTS = ("interferometer", "gain", "psf", "aca_pbcor", "twelve_m_pbcor", "casa_image", "casa_pb",
                "atlasgal_870", "tp_dust", "tp_dust_ff")
TP_INPUTS = ("tp_dust", "tp_dust_ff")


# ---------------------------------------------------------------- configuration and bookkeeping

def load_task_config(path: Path) -> dict:
    if not Path(path).is_file():
        raise FileNotFoundError(f"Task configuration not found: {path}")
    with open(path, "rb") as handle:
        return tomllib.load(handle)


def _project_path(project: ProjectPaths, value: object) -> Path:
    path = Path(str(value)).expanduser()
    return (path if path.is_absolute() else project.root / path).resolve()


def product_paths(output_root: Path, label: str = DEFAULT_LABEL) -> dict[str, Path]:
    """Combined images, named by method, TP model and the TP frequency label."""
    return {f"{method}_{kind}": Path(output_root) / "products" / f"{method}_{kind}_{label}.fits"
            for method in ("jcomb", "feather") for kind in ("continuum", "dust")}


def resolve_task_paths(project: ProjectPaths, config: Mapping) -> dict:
    """Resolve inputs, measurement sets, regions and output folders, checking each source rule."""
    inputs = {name: _project_path(project, value) for name, value in config["inputs"].items()}
    for name in DELIVERED_INPUTS:
        if not inputs[name].is_relative_to(project.delivered.resolve()):
            raise ValueError(f"{name} must be a delivered input: {inputs[name]}")
    synthesis = (project.root / SYNTHESIS_RESULTS).resolve()
    for name in SYNTHESIS_INPUTS:
        if not inputs[name].is_relative_to(synthesis):
            raise ValueError(f"{name} must be a synthesize_continuum product: {inputs[name]}")
    for name, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found: {path}")
    sets = {name: _project_path(project, value) for name, value in config["measurement_sets"].items()}
    for name, path in sets.items():
        if not path.is_relative_to(project.reduction.resolve()):
            raise ValueError(f"measurement set {name} must be a reduction input: {path}")
        if not path.is_dir():
            raise FileNotFoundError(f"measurement set {name} not found: {path}")
    regions = {name: _project_path(project, value) for name, value in config["regions"].items()}
    output_root = project.validate_output(str(config["outputs"]["root"]))
    label = (frequency_label(read_frequency_ghz(inputs["tp_dust_ff"])) if inputs["tp_dust_ff"].is_file()
             else DEFAULT_LABEL)
    return {**inputs, "inputs": dict(inputs), "measurement_sets": sets, "regions": regions,
            "output_root": output_root, "diagnostics_root": output_root / "diagnostics",
            "products": product_paths(output_root, label)}


def stages_for(selection: str) -> list[str]:
    """Stage names for ``combine``, ``diagnostics``, ``all`` or one stage."""
    if selection == "combine":
        return list(COMBINE_STAGES)
    if selection == "diagnostics":
        return list(DIAGNOSTIC_STAGES)
    if selection == "all":
        return [*COMBINE_STAGES, *DIAGNOSTIC_STAGES]
    if selection in (*COMBINE_STAGES, *DIAGNOSTIC_STAGES):
        return [selection]
    raise ValueError(f"Unknown stage: {selection}")


# Fails naming the stage when one of its products is missing.
require_outputs = require_inputs


def diagnostic_paths(paths: Mapping) -> dict:
    """The task paths with the output root moved to the diagnostics folder."""
    return {**paths, "output_root": paths["diagnostics_root"]}


def output_paths(result) -> list[str]:
    """Paths found in a stage's return value (a path, or a mapping or list of paths)."""
    if isinstance(result, Path):
        return [str(result)]
    if isinstance(result, Mapping):
        return [p for value in result.values() for p in output_paths(value)]
    if isinstance(result, (list, tuple)):
        return [p for value in result for p in output_paths(value)]
    return []


def append_manifest(root: Path, entry: Mapping) -> None:
    """Append one stage record to ``manifest.json`` in the output root."""
    path = Path(root) / "manifest.json"
    entries = json.loads(path.read_text()) if path.exists() else []
    record = {"stage": entry["stage"], "utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    record.update({k: ([str(x) for x in v] if k == "outputs" else v) for k, v in entry.items() if k != "stage"})
    entries.append(record)
    write_json(path, entries)


# ---------------------------------------------------------------- inspect

def _describe_image(name: str, path: Path) -> str:
    header = fits.getheader(path)
    shape = tuple(header[f"NAXIS{i}"] for i in range(header["NAXIS"], 0, -1))
    beam = tuple(header.get(key) for key in ("BMAJ", "BMIN", "BPA"))
    return f"{name}: shape {shape}, BUNIT {header.get('BUNIT')}, beam (BMAJ, BMIN, BPA) {beam}\n  {path}"


def run_inspect(paths: dict, config: dict, overwrite: bool) -> None:
    """Print the inputs, the TP frequencies and the output root. Writes nothing."""
    for name in IMAGE_INPUTS:
        print(_describe_image(name, paths[name]))
        if name in TP_INPUTS:
            print(f"  TP frequency: {read_frequency_ghz(paths[name])} GHz")
    for name, path in paths["measurement_sets"].items():
        print(f"measurement set {name}: {path}")
    print(f"output root: {paths['output_root']}")


# ---------------------------------------------------------------- diagnostics

def _diagnostic(function: Callable, *, figure: bool = True) -> Callable[[dict, dict, bool], object]:
    """A stage that runs ``function`` on the diagnostic paths."""
    def stage(paths: dict, config: dict, overwrite: bool) -> object:
        extra = {"save_figure": diagnostics._save_figure} if figure else {}
        return function(diagnostic_paths(paths), config, overwrite=overwrite, **extra)
    return stage


STAGES: dict[str, Callable[[dict, dict, bool], object]] = {
    "inspect": run_inspect, "prepare": combination_stages.run_prepare,
    "tp_factor": combination_stages.run_tp_factor, "box_spectrum": combination_stages.run_box_spectrum,
    "planck_calibration": combination_stages.run_planck_calibration, "jcomb": combination_stages.run_jcomb,
    "feather": combination_stages.run_feather,
    "uv": _diagnostic(diagnostics.plot_uv_coverage, figure=False),
    "overlap": _diagnostic(run_overlap_tests),
    "full_sampling": _diagnostic(build_full_sampling_profiles, figure=False),
    "amplitude_alignment": _diagnostic(diagnostics.run_amplitude_alignment),
    "figure2": _diagnostic(make_combination_figure),
    "compare_images": _diagnostic(diagnostics.plot_image_comparison, figure=False),
    "dirty_beam": _diagnostic(run_dirty_beam),
    "aca_forward_test": _diagnostic(run_aca_forward_test),
    "joint_images": _diagnostic(run_joint_image_comparison),
    "pointing_coverage": _diagnostic(diagnostics.run_pointing_coverage, figure=False)}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--stage", choices=[*COMBINE_STAGES, *DIAGNOSTIC_STAGES, "combine", "diagnostics", "all"],
                        default="all")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    project = ProjectPaths.discover()
    config = load_task_config(project.root / args.config if not args.config.is_absolute() else args.config)
    paths = resolve_task_paths(project, config)
    for stage in stages_for(args.stage):
        result = STAGES[stage](paths, config, args.overwrite)
        if stage == "inspect":
            continue
        append_manifest(paths["output_root"], {
            "stage": stage, "config": str(args.config),
            "settings": {k: v for k, v in config.items() if k not in ("inputs", "outputs")},
            "inputs_sha256": input_hashes(paths["inputs"], list(config["inputs"])),
            "outputs": output_paths(result)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
