#!/usr/bin/env python
"""Gas-to-dust mass ratio at the filament's Galactocentric radius.

Giannetti et al. (2017), A&A 606, L12, Eq. 2 (``lancet.gas_to_dust``).
Writes ``results/gas_to_dust_ratio/gas_to_dust_ratio.json``.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping
import json
from pathlib import Path
import tomllib

from lancet import ProjectPaths
from lancet.gas_to_dust import estimate_gas_to_dust_ratio


DEFAULT_CONFIG = Path("tasks/gas_to_dust_ratio/task.toml")
OUTPUT_NAME = "gas_to_dust_ratio.json"
RELATION = "Giannetti et al. (2017), A&A 606, L12, Eq. 2"


def load_task_config(path: Path) -> dict[str, object]:
    with Path(path).open("rb") as handle:
        return tomllib.load(handle)


def build_estimate(
    config: Mapping[str, object],
    output_root: Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Evaluate the relation and write the estimate as JSON."""
    output = Path(output_root) / OUTPUT_NAME
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output}; pass --overwrite to replace it")
    galaxy = config["galaxy"]
    estimate = estimate_gas_to_dust_ratio(
        float(galaxy["galactocentric_radius_kpc"]),
        float(galaxy["galactocentric_radius_error_kpc"]),
    )
    result = {"relation": RELATION, **estimate.to_dict()}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"R_gd = {estimate.ratio:.1f} "
        f"[stat -{estimate.statistical_ratio_minus:.1f} / +{estimate.statistical_ratio_plus:.1f}; "
        f"syst -{estimate.systematic_ratio_minus:.1f} / +{estimate.systematic_ratio_plus:.1f}; "
        f"intrinsic -{estimate.intrinsic_ratio_minus:.1f} / +{estimate.intrinsic_ratio_plus:.1f}] "
        f"at R_gal = {estimate.radius_kpc:g} ± {estimate.radius_error_kpc:g} kpc"
    )
    print(f"Wrote {output}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--stage", choices=("estimate", "all"), default="all")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    project = ProjectPaths.discover()
    config = load_task_config(project.root / args.config)
    output_root = project.validate_output(str(config["outputs"]["root"]))
    build_estimate(config, output_root, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
