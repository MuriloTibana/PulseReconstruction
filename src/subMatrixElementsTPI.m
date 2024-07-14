% Return contribution of given state and photon energy for two photons
function S = subMatrixElementsTPI(initial_energy, ...
    N_states,dipole_middle,middle_energies, ...
    final_state,dipole_final,final_energy, ...
    laser_parameters_one,laser_parameters_two, ...
    position_one,position_two)

%===== Laser One Parameters =================
frequency_one = laser_parameters_one(1);
period_one = 1/sqrt(1/laser_parameters_one(2)^2 + 1i*laser_parameters_one(6));
amplitude_one = laser_parameters_one(3) + 1i*laser_parameters_one(4);

%===== Laser Two Parameters =================
frequency_two = laser_parameters_two(1);
period_two = 1/sqrt(1/laser_parameters_two(2)^2 + 1i*laser_parameters_two(6));
amplitude_two = laser_parameters_two(3) + 1i*laser_parameters_two(4);

%===== Matrix Element =======================
S=0;
for mid_state = 1:N_states
    %===== Energies and Energy Gaps =========
    middle_energy = middle_energies(mid_state);
    delta_fm = final_energy - middle_energy - frequency_two;
    delta_mi = middle_energy - initial_energy - frequency_one;

    prefactor = -0.25 * pi * amplitude_two * amplitude_one * period_two * period_one;
    laser_factor_real = -0.5.*(period_one.^2 .* delta_mi.^2 + period_two.^2 .* delta_fm.^2);
    laser_factor_imag = 1i.*(position_two.*delta_fm + position_one.*delta_mi);
    erf_factor = 1./sqrt(2) .* 1./sqrt(period_one.^2 + period_two^2) .* (position_two - position_one + 1i.*(delta_fm.*period_two^2 - delta_mi.*period_one^2));
    laser1 = exp(laser_factor_imag);
    laser2 = exp(laser_factor_real);
    laser3 = error_function(erf_factor,laser_factor_real);
    laser = laser1 .* (laser2 + laser3);

    S = S + (prefactor .* dipole_final(mid_state,final_state) .* dipole_middle(mid_state) .* laser);
end % Loop over middle states
end % Function end