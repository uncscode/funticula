""" dynamic viscosity via sutherland formulation """


def get_dynamic_viscosity(
    reference_viscosity, temperature, reference_temperature, sutherland_constant
):
    """ https://resources.wolframcloud.com/FormulaRepository/resources/Sutherlands-Formula """
    return (
        reference_viscosity * (temperature / reference_temperature)**1.5 *
        (reference_temperature + sutherland_constant) /
        (temperature + sutherland_constant)
    )
