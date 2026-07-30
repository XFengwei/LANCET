"""Core physical-property calculations with quantity-aware uncertainty."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

from astropy import constants as const
from astropy import units as u
import numpy as np


def blackbody_cgs(nu: u.Quantity, temperature: u.Quantity) -> u.Quantity:
    """Planck function in cgs surface-brightness units."""
    return (
        2
        * const.h
        * nu**3
        / const.c**2
        / np.expm1(const.h * nu / const.k_B / temperature)
        / u.sr
    ).cgs


def brightness_to_column_density(
    *,
    nu: u.Quantity,
    kappa: u.Quantity,
    Iv: u.Quantity,
    temp: u.Quantity,
    mean_molecular_weight: float = 2.8,
) -> u.Quantity:
    """Invert brightness into H2 column density."""
    tau = -np.log(1 - Iv / blackbody_cgs(nu, temp))
    return (
        tau / (mean_molecular_weight * const.m_p * kappa)
    ).to(u.cm**-2)


def core_optical_depth(
    *,
    nu: u.Quantity,
    temp: u.Quantity,
    Isr: u.Quantity,
    Ibg: u.Quantity,
) -> u.Quantity:
    """Calculate core opacity including foreground/background emission."""
    blackbody = blackbody_cgs(nu, temp)
    foreground = Ibg / (np.sqrt(1 - Ibg / blackbody) + 1)
    return -np.log(
        1
        - Isr
        / (blackbody - foreground)
        * (1 - Ibg / blackbody) ** -0.5
    )


def symmetric_environment_optical_depth(
    *,
    nu: u.Quantity,
    temp: u.Quantity,
    Ibg: u.Quantity,
) -> u.Quantity:
    """Return total foreground-plus-background opacity for equal halves."""
    blackbody = blackbody_cgs(nu, temp)
    return -np.log(1 - Ibg / blackbody)


def thick_mass(
    *,
    nu: u.Quantity,
    kappa: u.Quantity,
    ftot: u.Quantity,
    temp: u.Quantity,
    tau_c: np.ndarray,
    dist: u.Quantity,
) -> u.Quantity:
    """Calculate opacity-corrected core mass."""
    thin_mass = (
        ftot * dist**2 / kappa / blackbody_cgs(nu, temp) / u.sr
    )
    correction = tau_c / (1 - np.exp(-tau_c))
    return (thin_mass * correction).to(u.M_sun)


def core_surface_density(
    *,
    nu: u.Quantity,
    kappa: u.Quantity,
    Isr: u.Quantity,
    temp: u.Quantity,
    Ibg: u.Quantity,
) -> u.Quantity:
    """Calculate core gas surface density."""
    return (
        core_optical_depth(
            nu=nu,
            temp=temp,
            Isr=Isr,
            Ibg=Ibg,
        )
        / kappa
    ).to(u.g / u.cm**2)


def core_volume_density(
    *,
    nu: u.Quantity,
    kappa: u.Quantity,
    ftot: u.Quantity,
    temp: u.Quantity,
    tau_c: np.ndarray,
    dist: u.Quantity,
    angle: u.Quantity,
) -> u.Quantity:
    """Calculate mean mass density inside the adopted FWHM radius."""
    mass = thick_mass(
        nu=nu,
        kappa=kappa,
        ftot=ftot,
        temp=temp,
        tau_c=tau_c,
        dist=dist,
    )
    linear = angle.value * dist.to_value(u.pc) * u.au
    density = 3 * mass / (4 * np.pi * linear**3)
    return density.to(1e-16 * u.g / u.cm**3)


def propagate_uncertainty(
    function: Callable[..., Any],
    parameters: Mapping[str, Mapping[str, Any]],
    *,
    sample_count: int = 10_000,
    random_seed: int | None = 42,
) -> dict[str, Any]:
    """Propagate normal, uniform, and fixed parameters by Monte Carlo."""
    if sample_count <= 0:
        raise ValueError("sample_count must be positive")
    rng = np.random.default_rng(random_seed)
    sampled: dict[str, Any] = {}
    for name, information in parameters.items():
        distribution = str(information["dist"]).lower()
        if distribution == "fixed":
            value = information["value"]
            sampled[name] = (
                np.full(sample_count, value.value) * value.unit
                if isinstance(value, u.Quantity)
                else np.full(sample_count, value)
            )
        elif distribution == "normal":
            mean = information["mean"]
            standard_deviation = information["std"]
            if isinstance(mean, u.Quantity):
                standard_deviation = standard_deviation.to(mean.unit)
                sampled[name] = rng.normal(
                    mean.value,
                    standard_deviation.value,
                    sample_count,
                ) * mean.unit
            else:
                sampled[name] = rng.normal(
                    mean,
                    standard_deviation,
                    sample_count,
                )
        elif distribution == "uniform":
            low = information["low"]
            high = information["high"]
            if isinstance(low, u.Quantity):
                high = high.to(low.unit)
                sampled[name] = rng.uniform(
                    low.value,
                    high.value,
                    sample_count,
                ) * low.unit
            else:
                sampled[name] = rng.uniform(low, high, sample_count)
        else:
            raise ValueError(f"Unsupported distribution: {distribution}")

    result = function(**sampled)
    if isinstance(result, u.Quantity):
        unit = result.unit
        values = np.asarray(result.value)
    else:
        unit = None
        values = np.asarray(result)
    summary = {
        "mean": np.mean(values),
        "std": np.std(values, ddof=1),
        "median": np.percentile(values, 50),
        "p16": np.percentile(values, 16),
        "p84": np.percentile(values, 84),
    }
    if unit is not None:
        summary = {name: value * unit for name, value in summary.items()}
    return summary
