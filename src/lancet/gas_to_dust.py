"""Gas-to-dust ratio from the Giannetti et al. (2017) relation.

The implementation follows Eq. 2 of:

Giannetti, A., Leurini, S., König, C., et al. 2017,
Astronomy & Astrophysics, 606, L12.
https://doi.org/10.1051/0004-6361/201731728
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import argparse
import json
import math
from typing import Sequence


@dataclass(frozen=True)
class Giannetti2017Relation:
    """Coefficients and uncertainties reported for Eq. 2."""

    slope: float = 0.087
    slope_statistical: float = 0.007
    slope_systematic_plus: float = 0.045
    slope_systematic_minus: float = 0.025
    intercept: float = 1.44
    intercept_statistical: float = 0.03
    intercept_systematic_plus: float = 0.21
    intercept_systematic_minus: float = 0.45
    intrinsic_scatter_dex: float = 0.2
    minimum_radius_kpc: float = 2.0
    maximum_radius_kpc: float = 20.0


@dataclass(frozen=True)
class GasToDustEstimate:
    """Central estimate with separate statistical and systematic intervals."""

    radius_kpc: float
    radius_error_kpc: float
    log10_ratio: float
    ratio: float
    statistical_log10_sigma: float
    statistical_ratio_minus: float
    statistical_ratio_plus: float
    systematic_log10_minus: float
    systematic_log10_plus: float
    systematic_ratio_minus: float
    systematic_ratio_plus: float
    intrinsic_scatter_dex: float
    intrinsic_ratio_minus: float
    intrinsic_ratio_plus: float

    def to_dict(self) -> dict[str, float]:
        """Return a JSON-serializable representation."""
        return asdict(self)


DEFAULT_RELATION = Giannetti2017Relation()


def estimate_gas_to_dust_ratio(
    radius_kpc: float,
    radius_error_kpc: float = 0.0,
    *,
    relation: Giannetti2017Relation = DEFAULT_RELATION,
    allow_extrapolation: bool = False,
) -> GasToDustEstimate:
    """Estimate the gas-to-dust mass ratio at a Galactocentric radius.

    Statistical uncertainty is propagated from the reported statistical
    coefficient errors and the radius uncertainty, treating them as
    independent because the coefficient covariance is not published in Eq. 2.
    The two systematic coefficient branches are evaluated as paired,
    anti-correlated shifts, as reported by Giannetti et al. (2017).
    """
    radius_kpc = float(radius_kpc)
    radius_error_kpc = float(radius_error_kpc)

    if not math.isfinite(radius_kpc):
        raise ValueError("radius_kpc must be finite")
    if not math.isfinite(radius_error_kpc) or radius_error_kpc < 0:
        raise ValueError("radius_error_kpc must be finite and non-negative")
    if (
        not allow_extrapolation
        and not relation.minimum_radius_kpc
        <= radius_kpc
        <= relation.maximum_radius_kpc
    ):
        raise ValueError(
            "Giannetti et al. (2017) Eq. 2 is calibrated for approximately "
            f"{relation.minimum_radius_kpc:g}–{relation.maximum_radius_kpc:g} kpc; "
            "pass allow_extrapolation=True to evaluate outside this range"
        )

    log10_ratio = relation.slope * radius_kpc + relation.intercept
    ratio = 10.0**log10_ratio

    statistical_log10_sigma = math.sqrt(
        (radius_kpc * relation.slope_statistical) ** 2
        + (relation.slope * radius_error_kpc) ** 2
        + relation.intercept_statistical**2
    )
    statistical_lower = 10.0 ** (log10_ratio - statistical_log10_sigma)
    statistical_upper = 10.0 ** (log10_ratio + statistical_log10_sigma)

    # The systematic brackets in Eq. 2 are paired: a steeper slope is
    # accompanied by a lower intercept, and vice versa.
    systematic_branch_1 = (
        relation.slope + relation.slope_systematic_plus
    ) * radius_kpc + (
        relation.intercept - relation.intercept_systematic_minus
    )
    systematic_branch_2 = (
        relation.slope - relation.slope_systematic_minus
    ) * radius_kpc + (
        relation.intercept + relation.intercept_systematic_plus
    )
    systematic_lower_log = min(systematic_branch_1, systematic_branch_2)
    systematic_upper_log = max(systematic_branch_1, systematic_branch_2)
    systematic_lower = 10.0**systematic_lower_log
    systematic_upper = 10.0**systematic_upper_log

    intrinsic_lower = 10.0 ** (
        log10_ratio - relation.intrinsic_scatter_dex
    )
    intrinsic_upper = 10.0 ** (
        log10_ratio + relation.intrinsic_scatter_dex
    )

    return GasToDustEstimate(
        radius_kpc=radius_kpc,
        radius_error_kpc=radius_error_kpc,
        log10_ratio=log10_ratio,
        ratio=ratio,
        statistical_log10_sigma=statistical_log10_sigma,
        statistical_ratio_minus=ratio - statistical_lower,
        statistical_ratio_plus=statistical_upper - ratio,
        systematic_log10_minus=log10_ratio - systematic_lower_log,
        systematic_log10_plus=systematic_upper_log - log10_ratio,
        systematic_ratio_minus=ratio - systematic_lower,
        systematic_ratio_plus=systematic_upper - ratio,
        intrinsic_scatter_dex=relation.intrinsic_scatter_dex,
        intrinsic_ratio_minus=ratio - intrinsic_lower,
        intrinsic_ratio_plus=intrinsic_upper - ratio,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate the gas-to-dust mass ratio with Giannetti et al. "
            "(2017), A&A, 606, L12, Eq. 2."
        )
    )
    parser.add_argument(
        "--radius-kpc",
        type=float,
        required=True,
        help="Galactocentric radius in kpc",
    )
    parser.add_argument(
        "--radius-error-kpc",
        type=float,
        default=0.0,
        help="One-sigma Galactocentric-radius uncertainty in kpc",
    )
    parser.add_argument(
        "--allow-extrapolation",
        action="store_true",
        help="Allow evaluation outside the approximately 2–20 kpc fit range",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print machine-readable JSON",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    args = _build_parser().parse_args(argv)
    estimate = estimate_gas_to_dust_ratio(
        args.radius_kpc,
        args.radius_error_kpc,
        allow_extrapolation=args.allow_extrapolation,
    )

    if args.json:
        print(json.dumps(estimate.to_dict(), indent=2, sort_keys=True))
        return 0

    radius = f"{estimate.radius_kpc:g} kpc"
    if estimate.radius_error_kpc:
        radius += f" ± {estimate.radius_error_kpc:g} kpc"

    print("Giannetti et al. (2017), A&A, 606, L12, Eq. 2")
    print(f"R_GC: {radius}")
    print(f"log10(gas-to-dust ratio): {estimate.log10_ratio:.4f}")
    print(f"gas-to-dust ratio: {estimate.ratio:.2f}")
    print(
        "statistical: "
        f"-{estimate.statistical_ratio_minus:.2f} "
        f"/ +{estimate.statistical_ratio_plus:.2f}"
    )
    print(
        "systematic:  "
        f"-{estimate.systematic_ratio_minus:.2f} "
        f"/ +{estimate.systematic_ratio_plus:.2f}"
    )
    print(
        f"intrinsic scatter ({estimate.intrinsic_scatter_dex:.1f} dex): "
        f"-{estimate.intrinsic_ratio_minus:.2f} "
        f"/ +{estimate.intrinsic_ratio_plus:.2f}"
    )
    return 0
