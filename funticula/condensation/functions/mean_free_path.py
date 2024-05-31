""" calculating the mean free path of air

    The mean free path is the average distance
    traveled by a molecule between collisions
    with other molecules present in a medium (air).

    The expeected mean free path of air is approx.
    65 nm at 298 K and 101325 Pa.

"""

from typing import Union
from numpy.typing import NDArray
import numpy as np


GAS_CONSTANT = 8.31446261815324  # J/(mol K)


def molecule_mean_free_path(
    molar_mass: Union[float, NDArray[np.float_]],
    temperature: float,
    pressure: float,
    dynamic_viscosity: float
) -> Union[float, NDArray[np.float_]]:
    """
    Calculate the mean free path of a gas molecule in air based on the
    temperature, pressure, and molar mass of the gas. The mean free path
    is the average distance traveled by a molecule between collisions with
    other molecules present in a medium (air).

    Args:
    -----
    - molar_mass (Union[float, NDArray[np.float_]]): The molar mass
    of the gas molecule [kg/mol]. Default is the molecular weight of air.
    - temperature (float): The temperature of the gas [K]. Default is 298.15 K.
    - pressure (float): The pressure of the gas [Pa]. Default is 101325 Pa.
    - dynamic_viscosity (Optional[float]): The dynamic viscosity of the gas
    [Pa*s]. If not provided, it will be calculated based on the temperature.

    Returns:
    --------
    - Union[float, NDArray[np.float_]]: The mean free path of the gas molecule
    in meters (m).

    References:
    ----------
    - https://en.wikipedia.org/wiki/Mean_free_path
    """

    return (
        (2 * dynamic_viscosity / pressure)
        / (8 * molar_mass / (np.pi * GAS_CONSTANT * temperature))**0.5,
    )
