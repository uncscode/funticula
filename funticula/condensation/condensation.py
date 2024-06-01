"""Particle Vapor Equilibrium, condensation and evaporation
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

# Functions
from funticula.condensation.functions.mass_transport import (
    mass_transfer_rate, first_order_mass_transport_k)
from funticula.condensation.functions.mean_free_path import (
    molecule_mean_free_path)
from funticula.condensation.functions.partial_pressure_module import (
    partial_pressure_delta
)
from funticula.condensation.functions.vapor_correction_module import (
    vapor_transition_correction)
from funticula.condensation.functions.knudsen_number_module import (
    calculate_knudsen_number)
from funticula.condensation.functions.dynamic_viscosity import (
    get_dynamic_viscosity)
from funticula.condensation.functions.kelvin_effect_module import (
    kelvin_term, kelvin_radius)


# input parameters
temperature = 298.15  # K
pressure = 101325  # Pa
dynamic_viscosity_ref = 1.8e-5  # Pa*s
dynamic_viscosity_temp_ref = 273.15  # K
radius = 1e-8  # m
molar_mass_air = 0.029  # kg/mol
mass_accommodation_coefficient = 0.1  # unitless
diffusion_coefficient_of_vapor = 2.0e-5  # m^2/s
# condensation species properties
effective_surface_tension = 0.072  # N/m
effective_density = 1000  # kg/m^3
molar_mass = 0.018  # kg/mol

# assume condensation partial pressure are given by another function
partial_pressure_gas = 1e-6  # Pa, based on vapor concentration
partial_pressure_particle = 1e-3  # Pa, based on activity and pure vap pressure
# How it is done in particula.next,
# Calculate the partial pressure
# partial_pressure_particle = particle.activity.partial_pressure(
#     pure_vapor_pressure=gas_species.get_pure_vapor_pressure(
#         temperature),
#     mass_concentration=particle.get_mass()
# )
# partial_pressure_gas = gas_species.get_partial_pressure(temperature)


# ********Start of the calculation********

# Dynamic viscosity
air_dynamic_viscosity = get_dynamic_viscosity(
    temperature=temperature,
    reference_viscosity=dynamic_viscosity_ref,
    reference_temperature=dynamic_viscosity_temp_ref
)
# mass transfer abstract class
air_mean_free_path = molecule_mean_free_path(
            molar_mass=molar_mass_air,
            temperature=temperature,
            pressure=pressure,
            dynamic_viscosity=air_dynamic_viscosity
)
# Calculate the Knudsen number of the particle
knudsen_number = calculate_knudsen_number(
    mean_free_path=air_mean_free_path,
    particle_radius=radius
)
# Calculate the vapor transition correction factor
vapor_transition = vapor_transition_correction(
        knudsen_number=knudsen_number,
        mass_accommodation=mass_accommodation_coefficient
    )
# Calculate the first-order mass transport coefficient
first_order_transport_k = first_order_mass_transport_k(
        radius=radius,
        vapor_transition=vapor_transition,
        diffusion_coefficient=diffusion_coefficient_of_vapor,
    )
# calculate the kelvin radius
kevlvin_radius_for_species = kelvin_radius(
    effective_surface_tension=effective_surface_tension,
    effective_density=effective_density,
    molar_mass=molar_mass,
    temperature=temperature
)
# calculate the kelvin term
kelvin_value = kelvin_term(
    radius=radius,
    kelvin_radius_value=kevlvin_radius_for_species,
)
# calculate the pressure delta accounting for the kelvin term
pressure_delta = partial_pressure_delta(
    partial_pressure_gas=partial_pressure_gas,
    partial_pressure_particle=partial_pressure_particle,
    kelvin_term=kelvin_value
)

# Calculate the mass transfer rate per particle
mass_condensation_rate = mass_transfer_rate(
    pressure_delta=pressure_delta,
    first_order_mass_transport=first_order_transport_k,
    temperature=temperature,
    molar_mass=molar_mass
)

print(f"Mass transfer rate: {mass_condensation_rate} kg/s")
