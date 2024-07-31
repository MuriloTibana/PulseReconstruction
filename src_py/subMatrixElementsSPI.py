import numpy as np

def subMatrixElementsSPI(dipole, initial_energy, final_energy, 
                            final_state, laser_parameters, position):
    # Laser Parameters
    frequency = laser_parameters[0]
    period = 1 / np.sqrt(1 / laser_parameters[1]**2 + 1j * laser_parameters[5])
    amplitude = laser_parameters[2] + 1j * laser_parameters[3]

    # Matrix Element
    prefactor = -1j * np.sqrt(np.pi / 2) * amplitude * laser_parameters[1]

    delta_fi = final_energy - initial_energy - frequency
    laser_factor = -1j * delta_fi * position - 0.5 * (delta_fi**2 * period**2)
    laser = np.exp(laser_factor)

    S = prefactor * dipole[final_state] * laser

    return S
