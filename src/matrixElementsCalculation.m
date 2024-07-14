% Returns the total photoionization autocorrelation wavefunction
function psi = matrixElementsCalculation(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    laser_list,correlation_delay,omega_list,chirp,N_gaussians,...
    del_pos,first_color,second_color,gaussian_variation)

%===== Add Removed Parameters to Laser ======
for pos = flip(del_pos)
    laser_list = [laser_list(1:pos-1) 0 laser_list(pos:end)];
end % Loop over deleted positions

if gaussian_variation
    % Handle different numbers of Gaussians for each color
    total_gaussians = first_color + second_color;
else
    total_gaussians = N_gaussians;
end

laser_list = reshape(laser_list, total_gaussians, []);

%===== Pad Laser Array with Zeros ===========
if and(~chirp, size(laser_list, 2) < 6)
    laser_list = [laser_list zeros(size(laser_list, 1), 1)];
end

%===== Add Frequency Parameter to Laser =====
color_override = [];
if gaussian_variation
    % Assign colors based on the number of Gaussians for each color
    for i = 1:first_color
        color_override = [color_override; omega_list(1)];
    end
    for i = 1:second_color
        color_override = [color_override; omega_list(2)];
    end
else
    % Default behavior
    for color = omega_list
        for pulse = 1:size(laser_list, 1) / size(omega_list, 2)
            color_override = [color_override; color];
        end % Loop over number of gaussians of each color
    end % Loop over each color
end

if size(laser_list, 2) < 6
    laser_list = [color_override laser_list(:,:)];
end

%===== Setup Storage Arrays =================
max_states = max([N_free_states_l0,N_free_states_l1,N_free_states_l2]);
psi = zeros(size(correlation_delay,2), ...
    max_states,3);

%===== Setup Padding Arrays =================
l0_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l0);
l1_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l1);
l2_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l2);

%===== First Laser Loop =====================
for i = 1:size(laser_list,1)
    laser_one = laser_list(i,:);
    %===== Second Laser Loop ================
    for j = 1:size(laser_list,1)
        laser_two = laser_list(j,:);
        %===== L=0 Final State Two-Photon ===
        T0 = matrixElementsTPI(initial_energy, ...
                N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
                N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
                laser_one,laser_two,correlation_delay);
        psi(:,:,1) = psi(:,:,1) + [squeeze(T0(2,:,:))' l0_padding];
        %===== L=2 Final State Two-Photon ===
        T2 = matrixElementsTPI(initial_energy, ...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
            laser_one,laser_two,correlation_delay);
        psi(:,:,3) = psi(:,:,3) + [squeeze(T2(2,:,:))' l2_padding];
    end % Loop over second laser
    %===== L=1 Final State One-Photon =======
    O = matrixElementsSPI(initial_energy, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        laser_one,correlation_delay);
    psi(:,:,2) = psi(:,:,2) + [squeeze(O(2,:,:))' l1_padding];
end % Loop over first laser
end % Function end