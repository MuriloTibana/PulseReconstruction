from subMatrixElementsSPI import * 

def matrixElementsSPI(initial_energy, 
                        N_states, dipole, final_energies, 
                        laser_parameters, correlation_delay):
    # Prepare Variables
    A = np.zeros((2, N_states, len(correlation_delay)), dtype=complex)
    position = laser_parameters[4]
    
    # Final State Summation
    for end_state in range(N_states):
        final_energy = final_energies[end_state]
        E_fi = final_energy - initial_energy
        A[0, end_state, 0] = final_energy

        # Calculate Contribution of State
        A[1, end_state, :] = subMatrixElementsSPI(dipole, initial_energy, 
                                                     final_energy, end_state, 
                                                     laser_parameters, position) * (
                                                     1 + np.exp(1j * E_fi * correlation_delay))

    return A