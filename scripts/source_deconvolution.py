#!/usr/bin/env python
"""Diagnostic calculation and plot for elliptical-Gaussian deconvolution."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

from lancet import (
    ProjectPaths,
    UnresolvedSourceError,
    convolve_gaussians,
    deconvolve_gaussian,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed-major", type=float, default=0.7747)
    parser.add_argument("--observed-minor", type=float, default=0.6559)
    parser.add_argument("--observed-pa", type=float, default=63.38)
    parser.add_argument("--beam-major", type=float, default=0.712)
    parser.add_argument("--beam-minor", type=float, default=0.606)
    parser.add_argument("--beam-pa", type=float, default=47.34)
    parser.add_argument(
        "--unit",
        default="arcsec",
        help="Axis unit used only for printed and plotted labels",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional PNG or PDF diagnostic figure",
    )
    return parser


def _add_ellipse(
    axis: plt.Axes,
    major: float,
    minor: float,
    pa_deg: float,
    **kwargs: object,
) -> None:
    # In an East/North tangent plane, astronomical PA maps to a Matplotlib
    # angle of 90 - PA, measured counterclockwise from +x (East).
    axis.add_patch(
        Ellipse(
            (0.0, 0.0),
            width=major,
            height=minor,
            angle=90.0 - pa_deg,
            fill=False,
            **kwargs,
        )
    )


def _position_angle_difference(first: float, second: float) -> float:
    """Smallest absolute difference for a 180-degree-periodic angle."""
    return abs((first - second + 90.0) % 180.0 - 90.0)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        source = deconvolve_gaussian(
            args.observed_major,
            args.observed_minor,
            args.observed_pa,
            args.beam_major,
            args.beam_minor,
            args.beam_pa,
            on_unresolved="raise",
        )
    except UnresolvedSourceError as error:
        print(f"Unresolved: {error}")
        return 2

    reconstructed = convolve_gaussians(
        *source,
        args.beam_major,
        args.beam_minor,
        args.beam_pa,
    )
    axis_residual = max(
        abs(reconstructed[0] - args.observed_major),
        abs(reconstructed[1] - args.observed_minor),
    )
    pa_residual = _position_angle_difference(
        reconstructed[2], args.observed_pa
    )

    print(
        "Observed:    "
        f"major={args.observed_major:.6g} {args.unit}, "
        f"minor={args.observed_minor:.6g} {args.unit}, "
        f"PA={args.observed_pa % 180:.6g} deg"
    )
    print(
        "Beam:        "
        f"major={args.beam_major:.6g} {args.unit}, "
        f"minor={args.beam_minor:.6g} {args.unit}, "
        f"PA={args.beam_pa % 180:.6g} deg"
    )
    print(
        "Deconvolved: "
        f"major={source[0]:.6g} {args.unit}, "
        f"minor={source[1]:.6g} {args.unit}, "
        f"PA={source[2]:.6g} deg"
    )
    print(
        "Reconstruction residual: "
        f"axes={axis_residual:.3e} {args.unit}, PA={pa_residual:.3e} deg"
    )

    if args.output is not None:
        project = ProjectPaths.discover()
        output = project.validate_output(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)

        fig, axis = plt.subplots(figsize=(6, 6), dpi=120)
        _add_ellipse(
            axis,
            args.observed_major,
            args.observed_minor,
            args.observed_pa,
            color="black",
            linewidth=2,
            label="Observed",
        )
        _add_ellipse(
            axis,
            args.beam_major,
            args.beam_minor,
            args.beam_pa,
            color="tab:red",
            linewidth=2,
            label="Beam",
        )
        _add_ellipse(
            axis,
            *source,
            color="tab:blue",
            linewidth=2,
            linestyle="--",
            label="Deconvolved source",
        )
        radius = 0.75 * max(
            args.observed_major,
            args.observed_minor,
            args.beam_major,
            args.beam_minor,
            source[0],
            source[1],
        )
        axis.set_xlim(-radius, radius)
        axis.set_ylim(-radius, radius)
        axis.set_aspect("equal")
        axis.set_xlabel(f"East offset ({args.unit})")
        axis.set_ylabel(f"North offset ({args.unit})")
        axis.legend(frameon=False)
        axis.set_title("Elliptical Gaussian deconvolution")
        fig.savefig(output, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {output.relative_to(project.root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
