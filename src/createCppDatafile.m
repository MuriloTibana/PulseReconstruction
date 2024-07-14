harm7 = gaussian_train_7.params();
temp = harm7(:,2);
harm7(:,2) = harm7(:,4);
harm7(:,4) = temp;

harm9 = gaussian_train_9.params();
temp = harm9(:,2);
harm9(:,2) = harm9(:,4);
harm9(:,4) = temp;

harm11 = gaussian_train_11.params();
temp = harm11(:,2);
harm11(:,2) = harm11(:,4);
harm11(:,4) = temp;

harm13 = gaussian_train_13.params();
temp = harm13(:,2);
harm13(:,2) = harm13(:,4);
harm13(:,4) = temp;

h5create("CppData.h5",'/LASER/HARM7',[6 7]);
h5write("CppData.h5",'/LASER/HARM7',harm7');

h5create("CppData.h5",'/LASER/HARM9',[6 7]);
h5write("CppData.h5",'/LASER/HARM9',harm9');

h5create("CppData.h5",'/LASER/HARM11',[6 7]);
h5write("CppData.h5",'/LASER/HARM11',harm11');

h5create("CppData.h5",'/LASER/HARM13',[6 7]);
h5write("CppData.h5",'/LASER/HARM13',harm13');

h5create("CppData.h5",'/ATOM/N_FREE_L0',1);
h5write("Cppdata.h5",'/ATOM/N_FREE_L0',N_free_states_l0);

h5create("CppData.h5",'/ATOM/N_FREE_L1',1);
h5write("Cppdata.h5",'/ATOM/N_FREE_L1',N_free_states_l1);

h5create("CppData.h5",'/ATOM/N_FREE_L2',1);
h5write("Cppdata.h5",'/ATOM/N_FREE_L2',N_free_states_l2);

h5create("CppData.h5",'/ATOM/N_BOUND_L1',1);
h5write("Cppdata.h5",'/ATOM/N_BOUND_L1',N_bound_states_l1);

h5create("CppData.h5",'/ATOM/TWO_PHOTON_DIPOLE_L0',[N_free_states_l0 N_bound_states_l1]);
h5write("CppData.h5",'/ATOM/TWO_PHOTON_DIPOLE_L0',real(two_photon_dipoles_l0'));

h5create("CppData.h5",'/ATOM/TWO_PHOTON_DIPOLE_L2',[N_free_states_l2 N_bound_states_l1]);
h5write("CppData.h5",'/ATOM/TWO_PHOTON_DIPOLE_L2',real(two_photon_dipoles_l2'));

h5create("CppData.h5",'/ATOM/TWO_PHOTON_DIPOLE_L1',N_bound_states_l1);
h5write("CppData.h5",'/ATOM/TWO_PHOTON_DIPOLE_L1',real(two_photon_dipoles_l1'));

h5create("CppData.h5",'/ATOM/ONE_PHOTON_DIPOLE_L1',N_free_states_l1);
h5write("CppData.h5",'/ATOM/ONE_PHOTON_DIPOLE_L1',real(one_photon_dipoles_l1'));

h5create("CppData.h5",'/ATOM/FREE_ENERGIES_L0',N_free_states_l0);
h5write("CppData.h5",'/ATOM/FREE_ENERGIES_L0',l0_free_energies');

h5create("CppData.h5",'/ATOM/FREE_ENERGIES_L2',N_free_states_l2);
h5write("CppData.h5",'/ATOM/FREE_ENERGIES_L2',l2_free_energies');

h5create("CppData.h5",'/ATOM/FREE_ENERGIES_L1',N_free_states_l1);
h5write("CppData.h5",'/ATOM/FREE_ENERGIES_L1',l1_free_energies');

h5create("CppData.h5",'/ATOM/BOUND_ENERGIES_L1',N_bound_states_l1);
h5write("CppData.h5",'/ATOM/BOUND_ENERGIES_L1',l1_bound_energies');

h5create("CppData.h5",'/ATOM/INITIAL_ENERGY',1);
h5write("CppData.h5",'/ATOM/INITIAL_ENERGY',initial_energy);
