import numpy as np
from errorFunction import errorFunction
def subMatrixElementsTPI(
    initial_energy, N_states, dipole_middle, middle_energies,
    final_state, dipole_final, final_energy,
    laser_parameters_one, laser_parameters_two,
    position_one, position_two
):
    # Laser One Paramters
    frequency_one = laser_parameters_one[0]
    period_one = 1 / np.sqrt(1 / laser_parameters_one[1]**2 + 1j * laser_parameters_one[5])
    amplitude_one = laser_parameters_one[2] + 1j * laser_parameters_one[3]

    # Laser Two Parameters
    frequency_two = laser_parameters_two[0]
    period_two = 1 / np.sqrt(1 / laser_parameters_two[1]**2 + 1j * laser_parameters_two[5])
    amplitude_two = laser_parameters_two[2] + 1j * laser_parameters_two[3]

    # Matrix Element
    S = 0
    for mid_state in range(N_states):
        # Energies and Energy Gaps
        middle_energy = middle_energies[mid_state]
        delta_fm = final_energy - middle_energy - frequency_two
        delta_mi = middle_energy - initial_energy - frequency_one

        prefactor = -0.25 * np.pi * amplitude_two * amplitude_one * period_two * period_one
        laser_factor_real = -0.5 * (period_one**2 * delta_mi**2 + period_two**2 * delta_fm**2)
        laser_factor_imag = 1j * (position_two * delta_fm + position_one * delta_mi)
        erf_factor = 1 / np.sqrt(2) * 1 / np.sqrt(period_one**2 + period_two**2) * (
            position_two - position_one + 1j * (delta_fm * period_two**2 - delta_mi * period_one**2))
    
        laser1 = np.exp(laser_factor_imag)
        laser2 = np.exp(laser_factor_real)
        laser3 = errorFunction(erf_factor, laser_factor_real)
        laser = laser1 * (laser2 + laser3)

        S += prefactor * dipole_final[mid_state, final_state] * dipole_middle[mid_state] * laser

    return S