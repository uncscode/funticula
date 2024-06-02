""" partial pressure """


def partial_pressure_delta(
    partial_pressure_gas, partial_pressure_particle, kelvin_term,
):
    """ ??? """
    return partial_pressure_gas - partial_pressure_particle * kelvin_term
