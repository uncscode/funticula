""" vapor transition correction """


def vapor_transition_correction(knudsen_number, mass_accommodation):
    """ Seinfeld and Pandis equation 12.43 """
    return (
        (0.75 * mass_accommodation * (1 + knudsen_number)) /
        ((knudsen_number**2 + knudsen_number)
         + 0.283 * mass_accommodation * knudsen_number
         + 0.75 * mass_accommodation)
    )
