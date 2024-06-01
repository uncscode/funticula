""" mass_transport """

def first_order_mass_transport_k(
    pi, radius, diffusion_coefficient, vapor_transition, 
):
    """ https://en.wikipedia.org/wiki/Mass_diffusivity """
    return 4 * pi * radius * diffusion_coefficient * vapor_transition


def mass_transfer_rate(
        first_order_mass_transport,
        pressure_delta,
        gas_constant,
        temperature,
        molar_mass
):
    """
    - Aerosol Modeling, eq 2.41 (excluding particle number)
    - Seinfeld and Pandis equation 13.3
    """
    return (
        first_order_mass_transport * pressure_delta /
        (gas_constant / molar_mass * temperature)
    )
