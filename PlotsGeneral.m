clc; clear; close all;
addpath('src');

%% ====== Reconstruction Parameters ====== %%
% defining parameters
tau_max = 1000; dtau = 0.2;
correlation_delay = linspace(-tau_max,tau_max,tau_max/dtau+1);
directory = 'results/';

time = linspace(-1000,1000,10001);
dt = abs(time(2) - time(1));

datas = dir(fullfile(directory, '*.mat'));
filename = datas.name;
for idx = 1:numel(datas)
    filename = datas(idx).name;
    load(fullfile(directory, filename))
    
    figure(idx);
    yyaxis left
    plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
    plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0.3, 1, 0.3], 'LineWidth', 1.5); 
    ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
    
    yyaxis right 
    plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
    plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
    ylabel('\omega(t) [a.u.]');
    
    legend({'$|\tilde{f}_{exact}(t)|$', '$|\tilde{f}_{estimated}(t)|$','$\omega_{exact}(t)$', '$\omega_{estimated}(t)$'}, 'Location', 'best', 'Interpreter', 'latex'); 
    xlabel('Time [a.u.]', 'Interpreter', 'latex');
    
    title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
    xlim([-600 600]); grid on;
end