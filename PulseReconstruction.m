%% Pulse Reconstruction 
% Reconstruct a single-color or multicolor pulse using single-color or multicolor Gaussians. 
% Choose the reconstruction between: 
%   - fitting chirp
%   - fitting frequency 
%   - using different number of gaussians for each color (only valid for 6 gaussians)
%   - overload from a previous calculation
% Choose the procedure between:
%   - raw calculation
%   - windowing 
%   - windowing with gaussian blur
%   - ultimate (mixing windowing and blur, choosing which one is best)
close all; clear; clc;
addpath('src');

%% --------- Run Setup --------- %% 
% Reconstruction Parameters: 
harmonicOne = 9; 
harmonicTwo = 11; 
harmonicThree = 10;
gaussians = 3; 
harmonics = 3;
reconstruct_single_color = false;
    % true = reconstruct a single-color pulse (9th)
    % false = reconstruct a multicolor pulse (9th + 11th)
pulse_single_color = false; 
    % true = use single-color gaussians
    % false = use multicolor gaussians
overload = true; 
    % true = load precalculated guess parameters
    % false = create guess parameters
    overload_data = '../3gaussian_nochirp/9h_10h_11h_3g_nc_windowing_blur.mat';
chirp = true;
    % true = don't remove chirp 
    % false = remove chirp
remove_frequency = true;
    % true = remove frequency 
    % false = don't remove frequency
gaussian_variation = false; % (only valid for 6 gaussians multicolor)
    % true = each color have different number of gaussians
    % false = each color have the same number of gaussians
    first_color = 4; 
    second_color = 6 - first_color;

N_iterations = 500; 
tau_max = 1000; dtau = 0.2;
correlation_delay = linspace(-tau_max, tau_max, tau_max / dtau + 1);

% Procedure Parameters: 
windowing = true; 
    checkpoint = false;
    N_windows = 10; 
    percentages = [5 15 25 50 75 100]; 
gaussian_blur = true; 
    blur_percentage = 0.1;
ultimate = false;

if ultimate 
    windowing = false; 
    gaussian_blur = false;
end

% Create data name 
if pulse_single_color
    data_name = sprintf('%dh_%dg', harmonicOne, gaussians);
else
    if gaussian_variation
        data_name = sprintf('%dh_%dh_%dgx%dg', harmonicOne, harmonicTwo, first_color, second_color);
    else
        if harmonics == 2
            data_name = sprintf('%dh_%dh_%dg', harmonicOne, harmonicTwo, gaussians);
        elseif harmonics == 3
            data_name = sprintf('%dh_%dh_%dh_%dg', harmonicOne, harmonicThree, harmonicTwo, gaussians);
        end
    end
end

data_name = [data_name '_' conditionalString(chirp, 'c', 'nc')];
data_name = conditionalAppend(data_name, overload, '_overload');
data_name = conditionalAppend(data_name, windowing, '_windowing');
data_name = conditionalAppend(data_name, gaussian_blur, '_blur');
data_name = conditionalAppend(data_name, ultimate, '_ultimate');
data_name = conditionalAppend(data_name, ~remove_frequency, '_f');
disp(['Data Name: ' data_name])

if overload 
    index_start = strfind(overload_data, '/');
    index_end = strfind(overload_data, '.');
    overload_name = overload_data(index_start(2)+1:index_end(3)-1);
    disp(['Overload Data: ' overload_name])
end
data_dir = 'data/';

%% --------- Precompute Dipoles for calc --------- %% 
[initial_energy, initial_state, l0_free_energies, l0_free_states, N_free_states_l0, ...
 l1_bound_energies, l1_bound_states, N_bound_states_l1, l1_free_energies, l1_free_states, ...
 N_free_states_l1, l2_free_energies, l2_free_states, N_free_states_l2, ...
 one_photon_dipoles_l1, two_photon_dipoles_l0, two_photon_dipoles_l1, two_photon_dipoles_l2] = precomputeDipoles(data_dir);


%% --------- Create Guess Parameters --------- %% 
if ~overload % create guess parameters (single-color or multicolor gaussians)
    amplitude = 2.4*10^(-3);
    FWHM = Laser.SI2au_duration(1e4);
    laser_chirp = 1e-4;
    position = zeros(1,gaussians);
    phase = pi;

    if pulse_single_color % single-color gaussians
        guess_laser = [];
        photon_energy = Laser.SI2au_wavelength(800 / harmonicOne);
  
        for i = 1:gaussians
            if mod(i,2) == 0   
                position(i) = (i/2)*50;
            else               
                position(i) = -((i+1)/2 - 1)*50;  
            end
            guess_laser = [guess_laser; Laser(amplitude, photon_energy, FWHM, laser_chirp, position(i), phase)];
        end
        omega = photon_energy;
    else  % multicolor gaussians
        if harmonics == 2
            guess_laserOne = []; guess_laserTwo = [];
            p_energyOne = Laser.SI2au_wavelength(800 / harmonicOne); 
            p_energyTwo = Laser.SI2au_wavelength(800 / harmonicTwo); 
            if gaussian_variation
                for i = 1:first_color
                    if mod(i,2) == 0   
                        position(i) = (i/2)*50;
                    else               
                        position(i) = -((i+1)/2 - 1)*50;  
                    end
                    guess_laserOne = [guess_laserOne; Laser(amplitude, p_energyOne, FWHM, laser_chirp, position(i), phase)];
                end
                for j = first_color+1:6
                    if mod(j,2) == 0   
                        position(j-first_color) = ((j)/2)*50;
                    else               
                        position(j-first_color) = -((j+1)/2 - 1)*50;  
                    end
                    guess_laserTwo = [guess_laserTwo; Laser(amplitude, p_energyTwo, FWHM, laser_chirp, position(j-first_color), phase)];
                end                
                guess_laser = [guess_laserOne; guess_laserTwo];
                omega = [p_energyOne p_energyTwo];
            else
                for i = 1:gaussians
                    if mod(i,2) == 0    
                        position = (i/2)*50;
                        guess_laserTwo = [guess_laserTwo; Laser(amplitude, p_energyTwo, FWHM, laser_chirp, position, phase)];
                    else                
                        position = -((i+1)/2 - 1)*50;  
                        guess_laserOne = [guess_laserOne; Laser(amplitude, p_energyOne, FWHM, laser_chirp, position, phase)];
                    end
                end
                guess_laser = [guess_laserOne; guess_laserTwo];
                omega = [p_energyOne p_energyTwo];
            end
        elseif harmonics == 3
            guess_laserOne = []; guess_laserTwo = []; guess_laserThree = [];
            p_energyOne = Laser.SI2au_wavelength(800 / harmonicOne); 
            p_energyTwo = Laser.SI2au_wavelength(800 / harmonicThree);
            p_energyThree = Laser.SI2au_wavelength(800 / harmonicTwo);  
            
            for i = 1:gaussians
                if mod(i,2) == 0    
                    position(i) = (i/2)*50;
                else               
                    position(i) = -((i+1)/2 - 1)*50;  
                end

                if mod(i, 3) == 1
                    guess_laserOne = [guess_laserOne; Laser(amplitude, p_energyOne, FWHM, laser_chirp, position(i), phase)];
                elseif mod(i, 3) == 2
                    guess_laserTwo = [guess_laserTwo; Laser(amplitude, p_energyTwo, FWHM, laser_chirp, position(i), phase)];
                else
                    guess_laserThree = [guess_laserThree; Laser(amplitude, p_energyThree, FWHM, laser_chirp, position(i), phase)];
                end
            end
        
            guess_laser = [guess_laserOne; guess_laserTwo; guess_laserThree];
            omega = [p_energyOne p_energyTwo p_energyThree];
        end
    end 
else % load precalculated guess parameters
    temp_chirp = chirp;
    temp_f = remove_frequency;
    temp_dataName = data_name;
    load(overload_data)
    guess_laser = estimated_laser;
    chirp = temp_chirp;
    remove_frequency = temp_f;
    lower_bound = []; upper_bound = []; 
    data_name = temp_dataName;
end
% guess_laser: [omega, T, fr, fi, phi, alpha]
% reshape parameters (remove frequency/remove chirp/remove imaginary amplitude and phase of the first gaussian)
guess = guess_laser.params(remove_frequency, chirp);

if remove_frequency % remove frequency
    lower_bound = ones(gaussians,1) * Laser(-1 - 100i,0.2,1,-0.1,-tau_max).params(remove_frequency,chirp);
    upper_bound = ones(gaussians,1) * Laser(1 + 100i,1.5,1000,0.1,tau_max).params(remove_frequency,chirp);
    delete_position = [4 3] * size(guess,1) - size(guess,1) + 1;

    calc = @(basis,delay,gaussians,del_pos)squeeze(sum(sum(abs( ...
    matrixElementsCalculation(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    basis,delay,omega,chirp,gaussians,del_pos,first_color,second_color,gaussian_variation)).^2,2),3));

else % don't remove frequency
    lower_bound(1:gaussians/2,:) = ones(gaussians/2,1) * Laser(-1 - 100i,0.4556,1,-0.1,-tau_max).params(remove_frequency,chirp);
    lower_bound(gaussians/2 + 1:gaussians,:) = ones(gaussians/2,1) * Laser(-1 - 100i,0.5695,1,-0.1,-tau_max).params(remove_frequency,chirp);
    upper_bound(1:gaussians/2,:) = ones(gaussians/2,1) * Laser(1 + 100i,0.5695,1000,0.1,tau_max).params(remove_frequency,chirp);
    upper_bound(gaussians/2 + 1:gaussians,:) = ones(gaussians/2,1) * Laser(1 + 100i,0.6834,1000,0.1,tau_max).params(remove_frequency,chirp);
    delete_position = [5 4] * size(guess,1) - size(guess,1) + 1;
    
    calc = @(basis,delay,gaussians,del_pos)squeeze(sum(sum(abs( ...
    matrixElementsCalculation(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    basis,delay,[],chirp,gaussians,del_pos,first_color,second_color,gaussian_variation)).^2,2),3));
end

disp('Parameters Before Reshape'); 
disp('Guess Parameters:')
for i = 1:size(guess, 1)
    disp(num2str(guess(i, :)))
end
disp('Lower Bounds:')
for j = 1:size(lower_bound,1)
    disp(num2str(lower_bound(i,:)))
end
disp('Upper Bounds:')
for j = 1:size(upper_bound,1)
    disp(num2str(upper_bound(i,:)))
end

guess = reshape(guess,1,[]);
lower_bound = reshape(lower_bound,1,[]);
upper_bound = reshape(upper_bound,1,[]);

for position = delete_position
    guess(position) = [];
    lower_bound(position) = [];
    upper_bound(position) = [];
end
disp('Parameters After Reshape'); 
disp(['Guess Parameters:' num2str(guess)])
disp(['Lower Bounds:' num2str(lower_bound)])
disp(['Upper Bounds:' num2str(upper_bound)])

%% --------- Load Autocorrelation to Reconstruct --------- %%
if reconstruct_single_color 
    % get known data for pulse made of 9th harmonic
    load([data_dir 'autocorrelation_singlecolor.mat'])
    known_autocorrelations = cell2mat(known_autocorrelation);
    known = known_autocorrelations(:, 1);
else
     % get known data for pulse train made of 9th + 11th harmonic
    load([data_dir 'autocorrelation_multicolor.mat'])
end

%% --------- Fitting Procedure --------- %%
options = optimoptions(@lsqnonlin, 'FunctionTolerance', 1e-16, ...
            'StepTolerance', 1e-16, 'OptimalityTolerance', 1e-16, ...
            'MaxFunctionEvaluations', N_iterations * 20, ...
            'MaxIterations', N_iterations, 'FiniteDifferenceType', 'forward', ...
            'UseParallel', true, 'Display', 'iter');

if isempty(gcp('nocreate'))
    parpool('local', 8);
end
disp(gcp('nocreate'));
if ~ultimate 
    if windowing 
        filter = @(time) interp1(correlation_delay,known,time);
        window_width = max(correlation_delay)/N_windows;
        max_omega = max(omega) * 4;
        max_frequency = max_omega / (2*pi);
        sample_step = 0.5 / max_frequency;

        for percent = percentages
            disp(['Percent of Data: ' num2str(percent) '%']);

            sample_tau = 0:sample_step:percent * window_width/100;
            tau = sample_tau;
            for window = 1:N_windows-1
                tau = [tau sample_tau + window * window_width];
            end  
            sample_data = filter(tau)';
            err = @(basis) abs(calc(basis,tau,gaussians,delete_position) - sample_data)./max(abs(sample_data));

            if gaussian_blur % run windowing with gaussian blur
                guess_idx = size(guess, 2);
                disp(['Guess parameters before gaussian blur: ' num2str(guess)]);
                for i = 1:guess_idx
                    guess(1, i) = guess(1,i) + blur_percentage * abs(guess(1,i)) * normrnd(0, 1/3);
                end
                disp(['Guess parameters after gaussian blur: ' num2str(guess)]);
                new_guess = lsqnonlin(err, guess, lower_bound, upper_bound, options);
                guess = new_guess; 

            else % run windowing without gaussian blur
                new_guess = lsqnonlin(err, guess, lower_bound, upper_bound, options);
                guess = new_guess;
            end

            if checkpoint 
                temp_guess = new_guess;
                save(['checkpoint_' data_name '_' num2str(percent) 'percent.mat'],'temp_guess');
            end 
        end 
    else % run raw lsqnonlin
        err = @(basis) abs(calc(basis,correlation_delay,gaussians,delete_position) - known)./max(abs(known));
        new_guess = lsqnonlin(err, guess, lower_bound, upper_bound, options);
        guess = new_guess;
    end
else % start the ultimate procedure
    % windowing with 5% of the data 
    filter = @(time) interp1(correlation_delay,known,time);
    window_width = max(correlation_delay)/N_windows;
    max_omega = max(omega) * 4;
    max_frequency = max_omega / (2*pi);
    sample_step = 0.5 / max_frequency;
    percent = 5;
    disp(['Percent of Data: ' num2str(percent) '%']);
    sample_tau = 0:sample_step:percent * window_width/100;
    tau = sample_tau;
    for window = 1:N_windows-1
        tau = [tau sample_tau + window * window_width];
    end  
    sample_data = filter(tau)';
    err = @(basis) abs(calc(basis,tau,gaussians,delete_position) - sample_data)./max(abs(sample_data));
    new_guess = lsqnonlin(err, guess, lower_bound, upper_bound, options);
    guess = new_guess;

    % Begin comparison between windowing and gaussian blur %  
    for percent = percentages(2:end)
        disp(['Percent of Data: ' num2str(percent) '%']);
        sample_tau = 0:sample_step:percent * window_width/100;
        tau = sample_tau;
        for window = 1:N_windows-1
            tau = [tau sample_tau + window * window_width];
        end  
        sample_data = filter(tau)';
        err = @(basis) abs(calc(basis,tau,gaussians,delete_position) - sample_data)./max(abs(sample_data));
        [guess_windowing, resnorm_windowing] = lsqnonlin(err, guess, lower_bound, upper_bound, options);

        guess_idx = size(guess, 2);
        disp(['Guess parameters before gaussian blur: ' num2str(guess)]);
        for i = 1:guess_idx
            guess(1,i) = guess(1,i) + blur_percentage * abs(guess(1,i)) * normrnd(0, 1/3);
        end
        disp(['Guess parameters after gaussian blur: ' num2str(guess)]);
        [guess_blur,resnorm_blur] = lsqnonlin(err, guess, lower_bound, upper_bound, options);

        disp(['Resnorm for Windowing: ' num2str(resnorm_windowing, '%.10f')])
        disp(['Resnorm for Gaussian Blur: ' num2str(resnorm_blur, '%.10f')])

        if resnorm_blur < resnorm_windowing
            disp('Gaussian Blur is better')
            guess = guess_blur;
        else
            disp('Windowing is better')
            guess = guess_windowing;
        end
    end 
end
guess = repairParameters(guess,delete_position,gaussians,omega,chirp,first_color,second_color,gaussian_variation);
disp('Final Guess')
disp('Guess Parameters:')
for i = 1:size(guess, 1)
    disp(num2str(guess(i, :)))
end

estimated_laser = Laser.generate(guess,[],[]);
calc = @(basis, delay) squeeze(sum(sum(abs( ...
    matrixElementsCalculation(initial_energy, ...
    N_free_states_l0, two_photon_dipoles_l0, l0_free_energies, ...
    N_free_states_l2, two_photon_dipoles_l2, l2_free_energies, ...
    N_bound_states_l1, two_photon_dipoles_l1, l1_bound_energies, ...
    N_free_states_l1, one_photon_dipoles_l1, l1_free_energies, ...
    basis, delay, [], true, size(basis, 1), [], [],[],[])).^2, 2), 3));

estimated_correlation = calc(estimated_laser.params(), correlation_delay);
save(data_name)

