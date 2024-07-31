import numpy as np
from matrixElementsSPI import matrixElementsSPI
from matrixElementsTPI import matrixElementsTPI

def matrixElementsCalculation(initial_energy, 
                              N_free_states_l0, two_photon_dipoles_l0, l0_free_energies, 
                              N_free_states_l2, two_photon_dipoles_l2, l2_free_energies, 
                              N_bound_states_l1, two_photon_dipoles_l1, l1_bound_energies, 
                              N_free_states_l1, one_photon_dipoles_l1, l1_free_energies, 
                              laser_list, correlation_delay):

    # Setup Storage Arrays
    max_states = max(N_free_states_l0, N_free_states_l1, N_free_states_l2)
    psi = np.zeros((len(correlation_delay), max_states, 3), dtype=complex)

    # Setup Padding Arrays
    l0_padding = np.zeros((len(correlation_delay), max_states - N_free_states_l0), dtype=complex)
    l1_padding = np.zeros((len(correlation_delay), max_states - N_free_states_l1), dtype=complex)
    l2_padding = np.zeros((len(correlation_delay), max_states - N_free_states_l2), dtype=complex)

    # First Laser Loop
    for i in range(laser_list.shape[0]):
        laser_one = laser_list[i, :]
        # Second Laser Loop
        for j in range(laser_list.shape[0]):
            laser_two = laser_list[j, :]
            # L=0 Final State Two-Photon
            T0 = matrixElementsTPI(initial_energy, 
                                   N_bound_states_l1, two_photon_dipoles_l1, l1_bound_energies, 
                                   N_free_states_l0, two_photon_dipoles_l0, l0_free_energies, 
                                   laser_one, laser_two, correlation_delay)
            psi[:, :, 0] += np.column_stack((np.squeeze(T0[1, :, :]).T, l0_padding))

            # L=2 Final State Two-Photon
            T2 = matrixElementsTPI(initial_energy, 
                                   N_bound_states_l1, two_photon_dipoles_l1, l1_bound_energies, 
                                   N_free_states_l2, two_photon_dipoles_l2, l2_free_energies, 
                                   laser_one, laser_two, correlation_delay)
            psi[:, :, 2] += np.column_stack((np.squeeze(T2[1, :, :]).T, l2_padding))

        # L=1 Final State One-Photon
        O = matrixElementsSPI(initial_energy, 
                              N_free_states_l1, one_photon_dipoles_l1, l1_free_energies, 
                              laser_one, correlation_delay)
        psi[:, :, 1] += np.column_stack((np.squeeze(O[1, :, :]).T, l1_padding))
    
    return psi
