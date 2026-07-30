"""Thermal Bonnor--Ebert stability calculations."""

from __future__ import annotations

from astropy import constants as const
from astropy import units as u
import numpy as np


def _require_positive(quantity: u.Quantity, name: str) -> u.Quantity:
    values = np.asarray(quantity.value, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError(f"{name} must contain finite positive values")
    return quantity


def sound_speed(
    temperature: u.Quantity,
    *,
    mean_particle_weight: float = 2.33,
) -> u.Quantity:
    """Return the isothermal sound speed."""
    temperature = _require_positive(
        u.Quantity(temperature).to(u.K),
        "temperature",
    )
    if not np.isfinite(mean_particle_weight) or mean_particle_weight <= 0:
        raise ValueError("mean_particle_weight must be finite and positive")
    return np.sqrt(
        const.k_B * temperature / (mean_particle_weight * const.m_p)
    ).to(u.km / u.s)


def nonthermal_velocity_dispersion(
    fwhm: u.Quantity,
    *,
    temperature: u.Quantity | None = None,
    tracer_molecular_weight: float | None = None,
) -> u.Quantity:
    """Convert a line FWHM to non-thermal one-dimensional dispersion.

    When ``temperature`` and ``tracer_molecular_weight`` are supplied, the
    tracer's thermal dispersion is removed in quadrature.
    """
    fwhm = _require_positive(
        u.Quantity(fwhm).to(u.km / u.s),
        "fwhm",
    )
    sigma_squared = (fwhm / np.sqrt(8.0 * np.log(2.0))) ** 2
    supplied = (temperature is not None, tracer_molecular_weight is not None)
    if any(supplied) and not all(supplied):
        raise ValueError(
            "temperature and tracer_molecular_weight must be supplied together"
        )
    if all(supplied):
        if (
            tracer_molecular_weight is None
            or not np.isfinite(tracer_molecular_weight)
            or tracer_molecular_weight <= 0
        ):
            raise ValueError(
                "tracer_molecular_weight must be finite and positive"
            )
        tracer_thermal = sound_speed(
            u.Quantity(temperature),
            mean_particle_weight=float(tracer_molecular_weight),
        )
        sigma_squared = sigma_squared - tracer_thermal**2
        if np.any(np.asarray(sigma_squared.value) < 0):
            raise ValueError(
                "The observed linewidth is narrower than the tracer thermal "
                "linewidth"
            )
    return np.sqrt(sigma_squared).to(u.km / u.s)


def mass_density_from_h2_number_density(
    number_density_h2: u.Quantity,
    *,
    mean_mass_per_h2: float = 2.8,
) -> u.Quantity:
    """Convert H2 number density to gas mass density."""
    number_density_h2 = _require_positive(
        u.Quantity(number_density_h2).to(u.cm**-3),
        "number_density_h2",
    )
    if not np.isfinite(mean_mass_per_h2) or mean_mass_per_h2 <= 0:
        raise ValueError("mean_mass_per_h2 must be finite and positive")
    return (
        mean_mass_per_h2 * const.m_p * number_density_h2
    ).to(u.g / u.cm**3)


def mass_density_from_surface_density(
    surface_density: u.Quantity,
    path_length: u.Quantity,
) -> u.Quantity:
    """Convert gas surface density to line-of-sight average mass density."""
    surface_density = _require_positive(
        u.Quantity(surface_density).to(u.g / u.cm**2),
        "surface_density",
    )
    path_length = _require_positive(
        u.Quantity(path_length).to(u.cm),
        "path_length",
    )
    return (surface_density / path_length).to(u.g / u.cm**3)


def effective_support_speed(
    temperature: u.Quantity,
    *,
    nonthermal_dispersion: u.Quantity | None = None,
    mean_particle_weight: float = 2.33,
) -> u.Quantity:
    """Return thermal or thermal-plus-nonthermal isotropic support speed."""
    thermal = sound_speed(
        temperature,
        mean_particle_weight=mean_particle_weight,
    )
    if nonthermal_dispersion is None:
        return thermal
    nonthermal = u.Quantity(nonthermal_dispersion).to(u.km / u.s)
    values = np.asarray(nonthermal.value, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise ValueError(
            "nonthermal_dispersion must contain finite non-negative values"
        )
    return np.sqrt(thermal**2 + nonthermal**2).to(u.km / u.s)


def external_pressure_from_density(
    *,
    temperature: u.Quantity,
    number_density_h2: u.Quantity | None = None,
    mass_density: u.Quantity | None = None,
    nonthermal_dispersion: u.Quantity | None = None,
    total_velocity_dispersion: u.Quantity | None = None,
    mean_particle_weight: float = 2.33,
    mean_mass_per_h2: float = 2.8,
) -> u.Quantity:
    """Calculate external pressure from density and velocity support.

    Supply a non-thermal dispersion to add gas thermal pressure in
    quadrature, or a total dispersion that already includes every component.
    """
    if (number_density_h2 is None) == (mass_density is None):
        raise ValueError(
            "Supply exactly one of number_density_h2 or mass_density"
        )
    if (
        nonthermal_dispersion is not None
        and total_velocity_dispersion is not None
    ):
        raise ValueError(
            "Supply nonthermal_dispersion or total_velocity_dispersion, "
            "not both"
        )
    if mass_density is None:
        density = mass_density_from_h2_number_density(
            u.Quantity(number_density_h2),
            mean_mass_per_h2=mean_mass_per_h2,
        )
    else:
        density = _require_positive(
            u.Quantity(mass_density).to(u.g / u.cm**3),
            "mass_density",
        )
    if total_velocity_dispersion is not None:
        velocity = _require_positive(
            u.Quantity(total_velocity_dispersion).to(u.km / u.s),
            "total_velocity_dispersion",
        )
    else:
        velocity = effective_support_speed(
            temperature,
            nonthermal_dispersion=nonthermal_dispersion,
            mean_particle_weight=mean_particle_weight,
        )
    return (density * velocity**2).to(u.dyne / u.cm**2)


def bonnor_ebert_mass_external_pressure(
    support_speed: u.Quantity,
    external_pressure: u.Quantity,
) -> u.Quantity:
    """Return the critical BE mass from support speed and measured pressure."""
    support = _require_positive(
        u.Quantity(support_speed).to(u.cm / u.s),
        "support_speed",
    )
    pressure = _require_positive(
        u.Quantity(external_pressure).to(u.dyne / u.cm**2),
        "external_pressure",
    )
    return (
        1.18 * support**4 / (const.G ** 1.5 * np.sqrt(pressure))
    ).to(u.M_sun)


def bonnor_ebert_mass_pressure(
    temperature: u.Quantity,
    *,
    number_density_h2: u.Quantity | None = None,
    mass_density: u.Quantity | None = None,
    nonthermal_dispersion: u.Quantity | None = None,
    mean_particle_weight: float = 2.33,
    mean_mass_per_h2: float = 2.8,
) -> u.Quantity:
    """Return the pressure-confined critical Bonnor--Ebert mass.

    The numerator uses the thermal sound speed. A supplied non-thermal
    dispersion contributes to the confining external pressure.
    """
    thermal_support = sound_speed(
        temperature,
        mean_particle_weight=mean_particle_weight,
    )
    pressure = external_pressure_from_density(
        temperature=temperature,
        number_density_h2=number_density_h2,
        mass_density=mass_density,
        nonthermal_dispersion=nonthermal_dispersion,
        mean_particle_weight=mean_particle_weight,
        mean_mass_per_h2=mean_mass_per_h2,
    )
    return bonnor_ebert_mass_external_pressure(
        thermal_support,
        pressure,
    )


def bonnor_ebert_mass_radius(
    temperature: u.Quantity,
    radius: u.Quantity,
    *,
    nonthermal_dispersion: u.Quantity | None = None,
    mean_particle_weight: float = 2.33,
) -> u.Quantity:
    """Return the radius-based critical BE mass with chosen support."""
    radius = _require_positive(u.Quantity(radius).to(u.pc), "radius")
    support = effective_support_speed(
        temperature,
        nonthermal_dispersion=nonthermal_dispersion,
        mean_particle_weight=mean_particle_weight,
    )
    return (2.4 * radius * support**2 / const.G).to(u.M_sun)


def implied_critical_boundary(
    temperature: u.Quantity,
    radius: u.Quantity,
    *,
    nonthermal_dispersion: u.Quantity | None = None,
    mean_particle_weight: float = 2.33,
    mean_mass_per_h2: float = 2.8,
) -> dict[str, u.Quantity]:
    """Return pressure and edge density implied by the critical radius form."""
    radius = _require_positive(u.Quantity(radius).to(u.cm), "radius")
    support = effective_support_speed(
        temperature,
        nonthermal_dispersion=nonthermal_dispersion,
        mean_particle_weight=mean_particle_weight,
    ).to(u.cm / u.s)
    pressure = (0.49**2 * support**4 / (const.G * radius**2)).to(
        u.dyne / u.cm**2
    )
    density = (pressure / support**2).to(u.g / u.cm**3)
    number_density = (
        density / (mean_mass_per_h2 * const.m_p)
    ).to(u.cm**-3)
    return {
        "external_pressure": pressure,
        "edge_mass_density": density,
        "edge_number_density_h2": number_density,
    }


def bonnor_ebert_ratio(
    core_mass: u.Quantity,
    critical_mass: u.Quantity,
) -> np.ndarray:
    """Return ``M_core / M_BE``; values above one are thermally unstable."""
    core_mass = _require_positive(
        u.Quantity(core_mass).to(u.M_sun),
        "core_mass",
    )
    critical_mass = _require_positive(
        u.Quantity(critical_mass).to(u.M_sun),
        "critical_mass",
    )
    return np.asarray(
        (core_mass / critical_mass).to_value(u.dimensionless_unscaled),
        dtype=float,
    )
