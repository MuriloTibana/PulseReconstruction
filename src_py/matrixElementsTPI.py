from subMatrixElementsTPI import *

def matrixElementsTPI(initial_energy, 
                        N_states_middle, dipole_middle, middle_energies, 
                        N_states_final, dipole_final, final_energies, 
                        laser_parameters_one, laser_parameters_two, correlation_delay):
    # Prepare Variables
    A = np.zeros((2, N_states_final, len(correlation_delay)), dtype=complex)
    position_one = laser_parameters_one[4]
    position_two = laser_parameters_two[4]
    omega_one = laser_parameters_one[0]
    omega_two = laser_parameters_two[0]
    # Final State Summation
    for end_state in range(N_states_final):
        final_energy = final_energies[end_state]
        if final_energy < 3 * max(omega_one, omega_two) + initial_energy:
            E_fi = final_energy - initial_energy
            A[0, end_state, 0] = final_energy
            # Both Photons from Same Laser
            term_one = subMatrixElementsTPI(initial_energy, 
                                              N_states_middle, dipole_middle, middle_energies, 
                                              end_state, dipole_final, final_energy, 
                                              laser_parameters_one, laser_parameters_two, 
                                              position_one, position_two) * (1 + np.exp(1j * E_fi * correlation_delay))

            # First Photon from Delayed
            term_two = subMatrixElementsTPI(initial_energy, 
                                              N_states_middle, dipole_middle, middle_energies, 
                                              end_state, dipole_final, final_energy, 
                                              laser_parameters_one, laser_parameters_two, 
                                              position_one + correlation_delay, position_two) * np.exp(1j * omega_one * correlation_delay)

            # Second Photon from Delayed
            term_three = subMatrixElementsTPI(initial_energy, 
                                                N_states_middle, dipole_middle, middle_energies, 
                                                end_state, dipole_final, final_energy, 
                                                laser_parameters_one, laser_parameters_two, 
                                                position_one, position_two + correlation_delay) * np.exp(1j * omega_two * correlation_delay)

            # Sum Contributions
            A[1, end_state, :] = term_one + term_two + term_three

    return A