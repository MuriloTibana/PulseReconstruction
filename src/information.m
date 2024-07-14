%===== Clear Previous Data ==================
close all; clear all;
%%
%===== Runtime Parameters ===================
%========== File and Directory Names ========
data_dir = './data/';
%========== Setup Delay Sweep ===============
tau_max = 1000; dtau = 0.2;
correlation_delay = linspace(-tau_max,tau_max,2*tau_max/dtau+1);

%%
%===== Read Data ============================
N_states = 500;
for l=0:2
    %===== Read Energies ====================
    data = h5read([data_dir 'He.h5'],char("/Energy_l"+l));
    if l==0
        nmax = length(data(:,1));
        energies = zeros(nmax,3);
    end
    energies(1:nmax-l,l+1) = data(:,1) + 1i*data(:,2);
    %===== Read States ======================
    data = h5read([data_dir 'He.h5'],char("/Psi_r_l"+l));
    if l==0
        rmax = length(data(:,1,1));
        states = zeros(rmax,nmax,3);
    end
    states(:,1:nmax-l,l+1) = data(:,:,1) + 1i*data(:,:,2);
end % Loop of final L
%========== Read Grid Parameters ============
grid = h5read([data_dir 'parameters.h5'],char("/EPS/r"));
dr = h5read([data_dir 'parameters.h5'],char("/EPS/delta_x"));

%%
%===== Organize Atomic States ===============
%========== L=0 Block =======================
angular_momentum = 0;
l0_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l0_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l0 = length(l0_free_energies);
%========== L=1 Block =======================
angular_momentum = 1;
l1_bound_energies = real(energies(real(energies(:,angular_momentum+1)) < 0,angular_momentum+1));
l1_bound_states = states(:,real(energies(:,angular_momentum+1)) < 0,angular_momentum+1);
N_bound_states_l1 = length(l1_bound_energies);
l1_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l1_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l1 = length(l1_free_energies);
%========== L=2 Block =======================
angular_momentum = 2;
l2_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l2_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l2 = length(l2_free_energies);
%========== Initial State ===================
initial_energy = real(energies(1,1));
initial_state = states(:,1,1);

%%
%===== Precompute Dipoles ===================
one_photon_dipoles_l1 = zeros(length(l1_free_energies),1);
for state = 1:length(l1_free_energies)
    one_photon_dipoles_l1(state) = sum(l1_free_states(:,state).*grid.*initial_state).*sqrt(3).*Wigner3j(1,1,0,0,0,0).^2;
end % Loop over free L=1 states
two_photon_dipoles_l1 = zeros(length(l1_bound_energies),1);
two_photon_dipoles_l0 = zeros(length(l1_bound_energies),length(l0_free_energies));
two_photon_dipoles_l2 = zeros(length(l1_bound_energies),length(l2_free_energies));
for state_one = 1:length(l1_bound_energies)
    for state_two = 1:length(l0_free_energies)
        two_photon_dipoles_l0(state_one,state_two) = sum(l1_bound_states(:,state_one).*grid.*l0_free_states(:,state_two)).*sqrt(3).*Wigner3j(0,1,1,0,0,0).^2;
    end % Loop over free L=0 states
    for state_two = 1:length(l2_free_energies)
        two_photon_dipoles_l2(state_one,state_two) = sum(l1_bound_states(:,state_one).*grid.*l2_free_states(:,state_two)).*sqrt(15).*Wigner3j(2,1,1,0,0,0).^2;
    end % Loop over free L=2 states
    two_photon_dipoles_l1(state_one) = sum(l1_bound_states(:,state_one).*grid.*initial_state).*sqrt(3).*Wigner3j(1,1,0,0,0,0).^2;
end % Loop over bound L=1 states

%%
%===== Load Laser Data ======================
load([data_dir 'helium_experiment_16g.mat']);

%%
%===== Setup Correlation Functions ==========
calc_cross = @(laser1, laser2, delay) squeeze(sum(sum(abs( ...
            matrixElementsCalculation_xcorr(initial_energy, ...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
            N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
            laser1.params(),laser2.params(),delay,[],true,size(laser1,1),[])).^2,2),3));

calc_auto = @(laser,delay) squeeze(sum(sum(abs( ...
            matrixElementsCalculation(initial_energy, ...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
            N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
            laser.params(),delay,[],true,size(laser,1),[])).^2,2),3));

%%
%===== 9th Harmonic Autocorrelation =========
harm9 = calc_auto(gaussian_train_9,correlation_delay);
info_harm9 = info(harm9,correlation_delay);

%%
%===== 11th Harmonic Autocorrelation ========
harm11 = calc_auto(gaussian_train_11,correlation_delay);
info_harm11 = info(harm11,correlation_delay);

%%
%===== 13th Harmonic Autocorrelation ========
harm13 = calc_auto(gaussian_train_13,correlation_delay);
info_harm13 = info(harm13,correlation_delay);

%%
%===== 11th x 9th Harmonic Autocorrelation ==
harm9x11 = calc_cross(gaussian_train_9, gaussian_train_11, correlation_delay);
info_harm9x11 = info(harm9x11,correlation_delay);

%%
%===== 13th x 9th Harmonic Autocorrelation ==
harm9x13 = calc_cross(gaussian_train_9, gaussian_train_13, correlation_delay);
info_harm9x13 = info(harm9x13,correlation_delay);

%%
%===== 11th x 13th Harmonic Autocorrelation ==
harm11x13 = calc_cross(gaussian_train_11, gaussian_train_13, correlation_delay);
info_harm11x13 = info(harm11x13,correlation_delay);

%%
%===== Setup Information Function ===========
function I = info(spectrum, delay)
    entropy = log2(spectrum);
    entropy(isnan(entropy)) = 0;
    I = trapz(delay, -1 * spectrum .* entropy);
end
