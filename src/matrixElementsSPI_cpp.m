function A = matrixElementsSPI_cpp(initial_energy, ...
    N_states,dipole,final_energies, ...
    laser_parameters,correlation_delay)

%===== Final State Summation ================
A = zeros(2,N_states,length(correlation_delay));
position = laser_parameters(5);
for end_state = 1:N_states
    final_energy = final_energies(end_state);
    E_fi = final_energy - initial_energy;
    A(1,end_state,1) = final_energy;

    A(2,end_state,:) = clib.multiphoton.OnePhotonAmplitude(initial_energy, ...
        dipole(end_state,1,1),final_energy,laser_parameters,position) ...
        .* (1 + exp(1i.*E_fi.*correlation_delay));
end
end