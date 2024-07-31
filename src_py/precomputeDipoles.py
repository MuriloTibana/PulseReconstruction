import h5py
import numpy as np
from Wigner3j import Wigner3j

def precomputeDipoles(data_dir):
    # Initialize storage
    energies = []
    states = []

    for l in range(3):
        with h5py.File(f'{data_dir}He.h5', 'r') as f:
            # ===== Read Energies and States ======
            data_energy = f[f'/Energy_l{l}']
            data_states = f[f'/Psi_r_l{l}']

            data_energy = np.array(data_energy)
            data_states = np.array(data_states)

            data_energy = data_energy.T
            data_states = data_states.T
            
            if l == 0:
                nmax = data_energy.shape[0]
                rmax = data_states.shape[0]
                energies = np.zeros((nmax, 3), dtype=complex)
                states = np.zeros((rmax, nmax, 3), dtype=complex)

            energies[:nmax-l, l] = data_energy[:, 0] + 1j * data_energy[:, 1]
            states[:, :nmax-l, l] = data_states[:, :, 0] + 1j * data_states[:, :, 1]

    with h5py.File(f'{data_dir}parameters.h5', 'r') as f: 
        grid = f[f'/EPS/r']
        grid = np.array(grid)

    # Organize Atomic States 
    angular_momentum = 0 
    l0_free_energies = np.real(energies[np.real(energies[:, angular_momentum]) >= 0, angular_momentum])
    l0_free_states = states[:, np.real(energies[:, angular_momentum]) >= 0, angular_momentum]
    N_free_states_l0 = len(l0_free_states[1])

    angular_momentum = 1
    l1_bound_energies = np.real(energies[np.real(energies[:, angular_momentum]) < 0, angular_momentum])
    l1_bound_states = states[:, np.real(energies[:, angular_momentum]) < 0, angular_momentum]
    N_bound_states_l1 = len(l1_bound_energies)
    l1_free_energies = np.real(energies[np.real(energies[:, angular_momentum]) >= 0, angular_momentum])
    l1_free_states = states[:, np.real(energies[:, angular_momentum]) >= 0, angular_momentum]
    N_free_states_l1 = len(l1_free_energies)

    angular_momentum = 2
    l2_free_energies = np.real(energies[np.real(energies[:, angular_momentum]) >= 0, angular_momentum])
    l2_free_states = states[:, np.real(energies[:, angular_momentum]) >= 0, angular_momentum]
    N_free_states_l2 = len(l2_free_energies)
    initial_energy = np.real(energies[0, 0])
    initial_state = states[:, 0, 0]

    # Calculate Dipoles
    one_photon_dipoles_l1 = np.zeros(len(l1_free_energies), dtype=complex)
    for state in range(len(l1_free_energies)):
        one_photon_dipoles_l1[state] = np.sum(l1_free_states[:, state] * grid * initial_state) * np.sqrt(3) * Wigner3j(1, 1, 0, 0, 0, 0)**2

    two_photon_dipoles_l1 = np.zeros(len(l1_bound_energies), dtype=complex)
    two_photon_dipoles_l0 = np.zeros((len(l1_bound_energies), len(l0_free_energies)), dtype=complex)
    two_photon_dipoles_l2 = np.zeros((len(l1_bound_energies), len(l2_free_energies)), dtype=complex)

    for state_one in range(len(l1_bound_energies)):
        for state_two in range(len(l0_free_energies)):
            two_photon_dipoles_l0[state_one, state_two] = np.sum(l1_bound_states[:, state_one] * grid * l0_free_states[:, state_two]) * np.sqrt(3) * Wigner3j(0, 1, 1, 0, 0, 0)**2

        for state_two in range(len(l2_free_energies)):
            two_photon_dipoles_l2[state_one, state_two] = np.sum(l1_bound_states[:, state_one] * grid * l2_free_states[:, state_two]) * np.sqrt(15) * Wigner3j(2, 1, 1, 0, 0, 0)**2

        two_photon_dipoles_l1[state_one] = np.sum(l1_bound_states[:, state_one] * grid * initial_state) * np.sqrt(3) * Wigner3j(1, 1, 0, 0, 0, 0)**2

    return (l0_free_energies, l0_free_states, N_free_states_l0, 
            l1_bound_energies, l1_bound_states, N_bound_states_l1, 
            l1_free_energies, l1_free_states, N_free_states_l1, 
            l2_free_energies, l2_free_states, N_free_states_l2, 
            initial_energy, initial_state, 
            one_photon_dipoles_l1, 
            two_photon_dipoles_l1, 
            two_photon_dipoles_l0, 
            two_photon_dipoles_l2)
