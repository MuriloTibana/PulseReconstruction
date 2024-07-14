% Returns the total photoionization cross-correlation wavefunction
function psi = matrixElementsCalculation_xcorr(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    unknown_laser_list,known_laser_list,correlation_delay,omega_list,chirp,...
    N_gaussians,del_pos)

%===== Add Removed Parameters to Laser ======
for pos = flip(del_pos)
    unknown_laser_list = [unknown_laser_list(1:pos-1) 0 unknown_laser_list(pos:end)];
end % Loop over deleted positions
unknown_laser_list = reshape(unknown_laser_list,N_gaussians,[]);


%===== Pad Laser Array with Zeros ===========
if(and(~chirp,size(unknown_laser_list,2)<6))
    unknown_laser_list = [unknown_laser_list zeros(size(unknown_laser_list,1),1)];
end

%===== Add Frequency Parameter to Laser =====
color_override = [];
for color = omega_list
    for pulse = 1:size(unknown_laser_list,1)/size(omega_list,2)
        color_override = [color_override; color];
    end % Loop over number of gaussians of each color
end % Loop over each color
if(size(unknown_laser_list,2) < 6)
    unknown_laser_list = [color_override unknown_laser_list(:,:)];
end

%===== Setup Storage Arrays =================
max_states = max([N_free_states_l0,N_free_states_l1,N_free_states_l2]);
psi = zeros(size(correlation_delay,2), ...
    max_states,3);

%===== Setup Padding Arrays =================
l0_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l0);
l1_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l1);
l2_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l2);

%===== L=1 Final State One-Photon ===========
%===== Photon from Unknown Laser ============
for i = 1:size(unknown_laser_list,1)
    unknown_laser = unknown_laser_list(i,:);
    O_unknown = matrixElementsSPI_xcorr(initial_energy,...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies,...
        unknown_laser,zeros(size(correlation_delay)));
    psi(:,:,2) = psi(:,:,2) + [squeeze(O_unknown(2,:,:))' l1_padding];
end % Loop over unknown laser
%===== Photon from Known Laser ==============
for i = 1:size(known_laser_list,1)
    known_laser = known_laser_list(i,:);
    O_known = matrixElementsSPI_xcorr(initial_energy,...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies,...
        known_laser,correlation_delay);
    %===== Phase Accumulation ===============
    O_known(2,:,:) = squeeze(O_known(2,:,:)) .* exp(1i.*ones(size(O_known,2),1)*known_laser(1)*correlation_delay);
    psi(:,:,2) = psi(:,:,2) + [squeeze(O_known(2,:,:))' l1_padding];
end % Loop over known laser

%===== First Photon from Unknown Laser ======
for i = 1:size(unknown_laser_list,1)
    laser_one = unknown_laser_list(i,:);
    %===== Second Photon from Unknown Laser =
    for j = 1:size(unknown_laser_list,1)
        laser_two = unknown_laser_list(j,:);
        %===== L=0 Final State Two-Photon ===
        T0_ff = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),zeros(size(correlation_delay)));
        psi(:,:,1) = psi(:,:,1) + [squeeze(T0_ff(2,:,:))' l0_padding];
        %===== L=2 Final State Two-Photon ===
        T2_ff = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),zeros(size(correlation_delay)));
        psi(:,:,3) = psi(:,:,3) + [squeeze(T2_ff(2,:,:))' l2_padding];
    end % Loop over second laser
end % Loop over first laser

%===== First Photon from Unknown Laser ======
for i = 1:size(unknown_laser_list,1)
    laser_one = unknown_laser_list(i,:);
    %===== Second Photon from Known Laser ===
    for j = 1:size(known_laser_list,1)
        laser_two = known_laser_list(j,:);
        %===== L=0 Final State Two-Photon ===
        T0_fg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),correlation_delay);
        %===== Phase Accumulation ===========
        T0_fg(2,:,:) = squeeze(T0_fg(2,:,:)) .* exp(1i.*ones(size(T0_fg,2),1)*laser_two(1)*correlation_delay);
        psi(:,:,1) = psi(:,:,1) + [squeeze(T0_fg(2,:,:))' l0_padding];
        %===== L=0 Final State Two-Photon ===
        T2_fg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),correlation_delay);
        %===== Phase Accumulation ===========
        T2_fg(2,:,:) = squeeze(T2_fg(2,:,:)) .* exp(1i.*ones(size(T2_fg,2),1)*laser_two(1)*correlation_delay);
        psi(:,:,3) = psi(:,:,3) + [squeeze(T2_fg(2,:,:))' l2_padding];
    end % Loop over second laser
end % Loop over first laser

%===== First Photon from Known Laser ========
for i = 1:size(known_laser_list,1)
    laser_one = known_laser_list(i,:);
    %===== Second Photon from Unknown Laser =
    for j = 1:size(unknown_laser_list,1)
        laser_two = unknown_laser_list(j,:);
        %===== L=0 Final State Two-Photon ===
        T0_gf = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            correlation_delay,zeros(size(correlation_delay)));
        %===== Phase Accumulation ===========
        T0_gf(2,:,:) = squeeze(T0_gf(2,:,:)) .* exp(1i.*ones(size(T0_gf,2),1)*laser_one(1)*correlation_delay);
        psi(:,:,1) = psi(:,:,1) + [squeeze(T0_gf(2,:,:))' l0_padding];
        %===== L=2 Final State Two-Photon ===
        T2_gf = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            correlation_delay,zeros(size(correlation_delay)));
        %===== Phase Accumulation ===========
        T2_gf(2,:,:) = squeeze(T2_gf(2,:,:)) .* exp(1i.*ones(size(T2_gf,2),1)*laser_one(1)*correlation_delay);
        psi(:,:,3) = psi(:,:,3) + [squeeze(T2_gf(2,:,:))' l2_padding];
    end % Loop over second laser
end % Loop over first laser

%===== First Photon from Known Laser ========
for i = 1:size(known_laser_list,1)
    laser_one = known_laser_list(i,:);
    %===== Second Photon from Known Laser ===
    for j = 1:size(known_laser_list,1)
        laser_two = known_laser_list(j,:);
        %===== L=0 Final State Two-Photon ===
        T0_gg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            correlation_delay,correlation_delay);
        %===== Phase Accumulation ===========
        T0_gg(2,:,:) = squeeze(T0_gg(2,:,:)) .* exp(1i.*ones(size(T0_gg,2),1)*(laser_one(1)+laser_two(1))*correlation_delay);
        psi(:,:,1) = psi(:,:,1) + [squeeze(T0_gg(2,:,:))' l0_padding];
        %===== L=2 Final State Two-Photon ===
        T2_gg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            correlation_delay,correlation_delay);
        %===== Phase Accumulation ===========
        T2_gg(2,:,:) = squeeze(T2_gg(2,:,:)) .* exp(1i.*ones(size(T2_gg,2),1)*(laser_one(1)+laser_two(1))*correlation_delay);
        psi(:,:,3) = psi(:,:,3) + [squeeze(T2_gg(2,:,:))' l2_padding];
    end % Loop over second laser
end % Loop over first laser
end % Function end