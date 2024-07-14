% Return contribution of given state and photon energy for one photon
function S = subMatrixElementsSPI(dipole,initial_energy,final_energy, ...
    final_state,laser_parameters,position)

%===== Laser Parameters =====================
frequency = laser_parameters(1);
period = 1/sqrt(1/laser_parameters(2)^2 + 1i*laser_parameters(6));
amplitude = laser_parameters(3) + 1i*laser_parameters(4);

%===== Matrix Element =======================
prefactor = -1i .* sqrt(pi/2) .* amplitude .* laser_parameters(2);

delta_fi = final_energy - initial_energy - frequency;
laser_factor = -1i .* delta_fi .* position - 0.5 .* (delta_fi.^2 .* period.^2);
laser = exp(laser_factor);

S = prefactor .* dipole(final_state,1,1) .* laser;

end % Function end