import numpy as np
import time
import h5py
import matplotlib.pyplot as plt

import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
src_py_path = os.path.join(current_dir, '../', 'src_py')
sys.path.insert(0, src_py_path)

from laser import Laser
from matrixElementsCalculation import matrixElementsCalculation
from precomputeDipoles import precomputeDipoles

def generate_random_lasers(num_gaussians):
    amplitudes = np.random.uniform(0, 5e-3, num_gaussians)
    photon_energies = np.random.uniform(0.2, 0.8, num_gaussians)
    phases = np.random.uniform(0, 2 * np.pi, num_gaussians)
    positions = np.random.uniform(-800, 800, num_gaussians)
    laser_chirps = np.random.uniform(0.5e-4, 5e-4, num_gaussians)
    FWHMs = Laser.SI2au_duration(np.random.uniform(5, 200, num_gaussians))
    
    lasers = [
        Laser(amplitude, photon_energy, FWHM, laser_chirp, position, phase)
        for amplitude, photon_energy, FWHM, laser_chirp, position, phase 
        in zip(amplitudes, photon_energies, FWHMs, laser_chirps, positions, phases)
    ]
    return lasers

start = time.time()
data_sets = 1

data_dir = '../data/'
(l0_free_energies, l0_free_states, N_free_states_l0, l1_bound_energies, 
 l1_bound_states, N_bound_states_l1, l1_free_energies, l1_free_states, 
 N_free_states_l1, l2_free_energies, l2_free_states, N_free_states_l2, 
 initial_energy, initial_state, one_photon_dipoles_l1, two_photon_dipoles_l1, 
 two_photon_dipoles_l0, two_photon_dipoles_l2) = precomputeDipoles(data_dir)

# Open the HDF5 file for writing all data sets
with h5py.File('test.h5', 'w') as hf:
    for data in range(data_sets):
        gaussians = 1
        guess_lasers = generate_random_lasers(gaussians)

        tau_max = 1000
        dtau = 0.2
        correlation_delay = np.linspace(-tau_max, tau_max, int(tau_max / dtau) + 1)

        guess_parameters = np.array([laser.params() for laser in guess_lasers])
        print(f'Parameters Data Set {data}:  {guess_parameters}')

        psi = matrixElementsCalculation(
            initial_energy,
            N_free_states_l0, two_photon_dipoles_l0, l0_free_energies,
            N_free_states_l2, two_photon_dipoles_l2, l2_free_energies,
            N_bound_states_l1, two_photon_dipoles_l1, l1_bound_energies,
            N_free_states_l1, one_photon_dipoles_l1, l1_free_energies,
            guess_parameters, correlation_delay)

        values = np.squeeze(np.sum(np.sum(np.abs(psi) ** 2, axis=2), axis=1))
        combined_data = np.vstack((correlation_delay, values)).T
        print(f'Combined Data Set {data}: {combined_data.shape}')
        # Create a group for each data set
        group = hf.create_group(f'data_set_{data}')
        group.create_dataset('guess_parameters', data=guess_parameters)
        group.create_dataset('combined_data', data=combined_data)

end = time.time()
print(f"Time Taken: {(end-start):.03f}s")
