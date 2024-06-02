""" Particle Vapor Equilibrium, condensation and evaporation
based on partial pressures to get dm/dt or other forms of
particle growth and decay.

From Seinfeld and Pandas: The Condensation (chapter 13) Equation 13.3
This function calculates the rate of change of the mass of an aerosol
particle with diameter Dp.

The rate of change of mass, dm, is given by the formula:
    dm/dt = 4π * radius * Di * Mi * f(Kn, alpha) * delta_pi / RT
where:
    radius is the radius of the particle,
    Di is the diffusion coefficient of species i,
    Mi is the molar mass of species i,
    f(Kn, alpha) is a transition function of the Knudsen number (Kn) and the
    mass accommodation coefficient (alpha),
    delta_pi is the partial pressure of species i in the gas phase vs the
    particle phase.
    R is the gas constant,
    T is the temperature.

    An additional denominator term is added to acount for temperature changes,
    which is need for cloud droplets, but can be used in general too.

This is also in Aerosol Modeling Chapter 2 Equation 2.40

Seinfeld, J. H., & Pandis, S. N. (2016). Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change. In Wiley (3rd ed.).
John Wiley & Sons, Inc.

Topping, D., & Bane, M. (2022). Introduction to Aerosol Modelling
(D. Topping & M. Bane, Eds.). Wiley. https://doi.org/10.1002/9781119625728

Units are all Base SI units.
"""

import numpy as np

from funticula.funx.mass_transport import (
    mass_transfer_rate, first_order_mass_transport_k)
from funticula.funx.mean_free_path import (
    molecule_mean_free_path)
from funticula.funx.partial_pressure import (
    partial_pressure_delta)
from funticula.funx.vapor_correction import (
    vapor_transition_correction)
from funticula.funx.knudsen_number import (
    calculate_knudsen_number)
from funticula.funx.dynamic_viscosity import (
    get_dynamic_viscosity)
from funticula.funx.kelvin_effect import (
    kelvin_term, kelvin_radius)


def mass_condensation_rate(
    dynamic_viscosity_ref, temperature, dynamic_viscosity_temp_ref, sutherland_constant,
    air_dynamic_viscosity, pressure, molar_mass_air, pi, gas_constant, radius,
    mass_accommodation_coefficient, effective_surface_tension, molar_mass, effective_density,
    partial_pressure_gas, partial_pressure_particle, diffusion_coefficient_of_vapor, 
):
    """ get the mass cond rate in kg/s """
    air_dynamic_viscosity = get_dynamic_viscosity(
        dynamic_viscosity_ref, temperature, dynamic_viscosity_temp_ref, SUTHERLAND_CONSTANT
    )
    air_mean_free_path = molecule_mean_free_path(
        air_dynamic_viscosity, pressure, molar_mass_air, PI, GAS_CONSTANT, temperature,
    )
    knudsen_number = calculate_knudsen_number(
        air_mean_free_path, radius
    )
    vapor_transition = vapor_transition_correction(
            knudsen_number, mass_accommodation_coefficient
    )
    kelvin_value = kelvin_term(
        effective_surface_tension, molar_mass, GAS_CONSTANT, temperature, effective_density, radius
    )
    pressure_delta = partial_pressure_delta(
        partial_pressure_gas,
        partial_pressure_particle,
        kelvin_value
    )
    mass_condensation_rate = mass_transfer_rate(
        PI, radius, diffusion_coefficient_of_vapor, vapor_transition,
        pressure_delta,
        GAS_CONSTANT, 
        temperature,
        molar_mass
    )
    return mass_condensation_rate
