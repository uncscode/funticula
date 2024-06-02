""" knudsen number """


def calculate_knudsen_number(mean_free_path, particle_radius):
    """ https://en.wikipedia.org/wiki/Knudsen_number """
    return mean_free_path / particle_radius
