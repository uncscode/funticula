""" mean free path of air """


def molecule_mean_free_path(
    dynamic_viscosity, pressure, molar_mass, pi, gas_constant, temperature
):
    """ https://en.wikipedia.org/wiki/Mean_free_path """
    return (
        (2 * dynamic_viscosity / pressure)
        / (8 * molar_mass / (pi * gas_constant * temperature))**0.5,
    )
