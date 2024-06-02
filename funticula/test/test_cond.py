""" test """

from funticula.cond import mass_condensation_rate


def test_funticula():
    """ test funticula """
    dynamic_viscosity_ref = 1.8e-5  # Pa*s
    temperature = 298.15  # K
    dynamic_viscosity_temp_ref = 273.15  # K
    SUTHERLAND_CONSTANT = 110.4
    pressure = 101325  # Pa
    molar_mass_air = 0.029  # kg/mol
    PI = np.pi    
    GAS_CONSTANT = 8.314  # J/(mol*K)
    radius = 1e-8  # m
    mass_accommodation_coefficient = 0.1  # unitless
    effective_surface_tension = 0.072  # N/m
    molar_mass = 0.018  # kg/mol
    effective_density = 1000  # kg/m^3
    partial_pressure_gas = 1e-6  # Pa, based on vapor concentration
    partial_pressure_particle = 1e-3  # Pa, based on activity and pure vap pressure
    diffusion_coefficient_of_vapor = 2.0e-5  # m^2/s
    mcrv = mass_condensation_rate(
        dynamic_viscosity_ref, temperature, dynamic_viscosity_temp_ref, sutherland_constant,
        air_dynamic_viscosity, pressure, molar_mass_air, pi, gas_constant, radius,
        mass_accommodation_coefficient, effective_surface_tension, molar_mass, effective_density,
        partial_pressure_gas, partial_pressure_particle, diffusion_coefficient_of_vapor, 
    )
