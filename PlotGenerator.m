%% Plotting Pulse and Autocorrelation graphs 
clear; clc; close all;
addpath('src');

%% ====== Data File ====== %%
% setting up directories and file name
name_file = '9h_11h_5gx1g_c_windowing_blur.mat';
has11h = contains(name_file, '11h');
has10h = contains(name_file, '10h');
save_everything = true;
if has11h 
    plot_dir = 'figures\multi_color\';
    error_dir = 'errors\multi_color\';
else
    plot_dir = 'figures\single_color\';
    error_dir = 'errors\single_color\';
end

% load data
load(fullfile('results\', name_file));

% parsing file name
dot_index = strfind(name_file, 'c');
name_plot = name_file(1:dot_index);

%% ====== Reconstruction Parameters ====== %%
% defining parameters
time_reverse = false; initial_guess = [0 200];
tau_max = 1000; dtau = 0.2;
correlation_delay = linspace(-tau_max,tau_max,tau_max/dtau+1);

time = linspace(-1000,1000,10001);
dt = abs(time(2) - time(1));
% defining calc function to calculate autocorrelation graph of the estimated laser
calc = @(basis,delay) squeeze(sum(sum(abs( ...
        matrixElementsCalculation(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,delay,[],true,size(basis,1),[])).^2,2),3));

%% ====== Calculating Field Error ====== %%
% exact pulse: 
load('data/helium_experiment_16g.mat')
harmonics = [9, 11, 13];

experiment_components = cell(numel(harmonics), 1);
experiment_lasers = cell(numel(harmonics), 1);
experiment_parameters = cell(numel(harmonics), 1);
known_autocorrelation = cell(numel(harmonics, 1));

for i = 1:numel(harmonics)
    harmonic = harmonics(i);
    experiment_components{i} = ['gaussian_train_', num2str(harmonic)];
    experiment_lasers{i} = eval([experiment_components{i}, ';']);
end

harmonic9_laser = experiment_lasers{1};
harmonic11_laser = experiment_lasers{2};

if has11h || has10h
    experiment_laser = [harmonic9_laser; harmonic11_laser];
else
    experiment_laser = harmonic9_laser;
end
experiment_vals = calculate(experiment_laser, time);

% fit error: 
% define fit error function 
if time_reverse
    fit_guess = @(params) abs(exp(1i*params(1)) * flip(conj(Laser.generate(guess,true).calculate(time - params(2)))) ...
    - experiment_vals) / (max(abs(experiment_vals)));
else
    fit_guess = @(params) abs(exp(1i*params(1)) * Laser.generate(guess,true).calculate(time - params(2)) ...
    - experiment_vals) / (max(abs(experiment_vals)));
end

% minimizing fit error function 
options = optimoptions(@lsqnonlin, ...
    'FunctionTolerance', 1e-18, ...            
    'StepTolerance', 1e-18, ...                
    'OptimalityTolerance', 1e-15, ...           
    'MaxFunctionEvaluations', 1e7, ...        
    'MaxIterations', 500, ...                 
    'FiniteDifferenceType', 'forward', ...     
    'UseParallel', true, ...                   
    'Display', 'iter-detailed');               

parameters = lsqnonlin(fit_guess, initial_guess, [], [], options);

if time_reverse
    estimated_vals = exp(1i*parameters(1)) * flip(conj(Laser.generate(guess,true).calculate(time - parameters(2) + 10)));
    estimated_vals_9th = exp(1i*parameters(1)) * flip(conj(Laser.generate(guess,true).calculate(time - parameters(2) + 10)));
    estimated_vals_11th = exp(1i*parameters(1)) * flip(conj(Laser.generate(guess,true).calculate(time - parameters(2) + 10)));
    center_graph = parameters(2);
else
    estimated_vals = exp(1i*parameters(1)) * Laser.generate(guess,true).calculate(time - parameters(2));
    estimated_vals_9th = exp(1i*parameters(1)) * Laser.generate(guess,true).calculate(time - parameters(2));
    estimated_vals_11th = exp(1i*parameters(1)) * Laser.generate(guess,true).calculate(time - parameters(2));
    center_graph = -parameters(2);
end

% calculating
field_error = sqrt(dtau * sum(abs(experiment_vals - estimated_vals).^2));

estimated_phase = unwrap(angle(estimated_vals));
experiment_phase = unwrap(angle(experiment_vals));
estimated_omega = zeros(size(estimated_phase));
experiment_omega = zeros(size(experiment_phase));

for t = 2:length(time)-1
    estimated_omega(t) = -(estimated_phase(t+1) - estimated_phase(t-1)) / (2.0 *dt);
    experiment_omega(t) = -(experiment_phase(t+1) - experiment_phase(t-1)) / (2.0 *dt);
end

%% ====== Calculating Ionization Error ====== %%
obj = abs(estimated_correlation - known)./max(abs(known));
ion_error = sqrt(dtau * sum(obj.^2));

%% ====== Plotting ====== %% 
property_label = {'FontSize', 14, 'FontName', 'Times New Roman'};
property_title = {'FontSize', 12, 'FontName', 'Times New Roman'};
center_graph = -250;
% pulse comparison
fig1 = figure(1);
yyaxis left
plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0.3, 1, 0.3], 'LineWidth', 1.5); hold off; grid on;
ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex', property_label{:});

yyaxis right 
plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); hold on;
plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); grid on; 
ylabel('\omega(t) [a.u.]', property_label{:});

legend({'$|\tilde{f}_{exact}(t)|$', '$|\tilde{f}_{estimated}(t)|$','$\omega_{exact}(t)$', '$\omega_{estimated}(t)$'}, 'Location', 'best', 'Interpreter', 'latex'); 
xlabel('Time [a.u.]', 'Interpreter', 'latex', property_label{:});

title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none', property_title{:});
if has11h 
    xlim([-600 600]);
else
    xlim([-500 500]);
end
grid on;

% autocorrelation comparison
fig2 = figure(2);
plot(correlation_delay, known, '-','Color', [0.1, 0.4, 0.6], 'LineWidth', 2); hold on;
plot(correlation_delay, estimated_correlation, ':', 'Color', [0.5, 0.9, 0.5], 'LineWidth', 2); hold off; grid on;

legend('Exact', 'Estimated', 'Location', 'best');
xlabel('Correlation Delay [a.u.]', property_label{:});
ylabel('Ionization Probability', property_label{:});
title(sprintf('Autocorrelation Comparison: %s\nIonization Error: %.10f', name_plot, ion_error), 'Interpreter', 'none', property_title{:});
xlim([-600 600]);
grid on;

%% ======= Saving Everything ======= % 
if save_everything
    file_name_pulse = strcat(name_plot,'_PulseComparison');
    saveas(fig1,fullfile(plot_dir, file_name_pulse),'jpeg');
    
    file_name_autocorrelation = strcat(name_plot,'_AutocorrelationComparison');
    saveas(fig2,fullfile(plot_dir, file_name_autocorrelation),'jpeg');
    
    file_name_errors = strcat(name_plot, '_erros');
    save(fullfile(error_dir, file_name_errors), "field_error", "ion_error")

    save(fullfile('results\', name_file),"field_error", "ion_error", "experiment_vals", "estimated_vals", ...
        "experiment_omega", "estimated_omega", "center_graph", "name_plot" ,"-append")
end