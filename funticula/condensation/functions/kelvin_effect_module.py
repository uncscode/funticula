""" kelvin effect: radius and term """

import numpy as np


def kelvin_radius(
    effective_surface_tension, molar_mass, gas_constant,
    temperature, effective_density
):
    """ https://en.wikipedia.org/wiki/Kelvin_equation """
    return (
        (2 * effective_surface_tension * molar_mass) /
        (gas_constant * temperature * effective_density)
    )


def kelvin_term(
    effective_surface_tension, molar_mass, gas_constant,
    temperature, effective_density, radius
):
    """ https://en.wikipedia.org/wiki/Kelvin_equation """
    kelvin_radius_value = kelvin_radius(
        effective_surface_tension, molar_mass, gas_constant,
        temperature, effective_density
    )
    return np.exp(kelvin_radius_value / radius)
