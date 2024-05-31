"""Mass Transport Functions
"""
from typing import Union
from numpy.typing import NDArray
import numpy as np

GAS_CONSTANT = 8.31446261815324  # J/(mol K)


def first_order_mass_transport_k(
        radius: Union[float, NDArray[np.float_]],
        vapor_transition: Union[float, NDArray[np.float_]],
        diffusion_coefficient: Union[float, NDArray[np.float_]]
) -> Union[float, NDArray[np.float_]]:
    """
    Calculate the first-order mass transport coefficient, K, for a given radius
    diffusion coefficient, and vapor transition correction factor. For a
    single particle.

    Args:
    -----
    - radius (Union[float, NDArray[np.float_]]): The radius of the particle
    [m].
    - diffusion_coefficient (Union[float, NDArray[np.float_]]): The diffusion
    coefficient of the vapor [m^2/s]
    - vapor_transition (Union[float, NDArray[np.float_]]): The vapor transition
    correction factor. [unitless]

    Returns:
    --------
    - Union[float, NDArray[np.float_]]: The first-order mass transport
    coefficient per particle (m^3/s).

    References:
    ----------
    - Aerosol Modeling, Chapter 2, Equation 2.49 (excluding particle number)
    - https://en.wikipedia.org/wiki/Mass_diffusivity
    """
    return 4 * np.pi * radius * diffusion_coefficient * vapor_transition


def mass_transfer_rate(
        pressure_delta: Union[float, NDArray[np.float_]],
        first_order_mass_transport: Union[float, NDArray[np.float_]],
        temperature: Union[float, NDArray[np.float_]],
        molar_mass: Union[float, NDArray[np.float_]]
) -> Union[float, NDArray[np.float_]]:
    """
    Calculate the mass transfer rate based on the difference in partial
    pressure and the first-order mass transport coefficient.

    Args:
    -----
    - pressure_delta (Union[float, NDArray[np.float_]]): The difference in
    partial pressure between the gas phase and the particle phase.
    - first_order_mass_transport (Union[float, NDArray[np.float_]]): The
    first-order mass transport coefficient per particle.
    - temperature (Union[float, NDArray[np.float_]]): The temperature at which
    the mass transfer rate is to be calculated.

    Returns:
    --------
    - Union[float, NDArray[np.float_]]: The mass transfer rate for the particle
    [kg/s].

    References:
    ----------
    - Aerosol Modeling, Chapter 2, Equation 2.41 (excluding particle number)
    - Seinfeld and Pandis, "Atmospheric Chemistry and Physics", Equation 13.3
    """
    return np.array(
        first_order_mass_transport * pressure_delta
        / (GAS_CONSTANT / molar_mass * temperature),
        dtype=np.float_
    )
