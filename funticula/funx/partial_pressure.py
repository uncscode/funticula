""" partial pressure """


def partial_pressure_delta(
    partial_pressure_gas, partial_pressure_particle, kelvin_term,
):
    """ Partial pressure difference gas phase vs curved liquid surface """
    return partial_pressure_gas - partial_pressure_particle * kelvin_term
