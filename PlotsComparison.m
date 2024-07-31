clc; clear; close all;
addpath('src');

%% ====== Plotting Everything ====== %%
% defining parameters
directory = 'results/';
time = linspace(-1000,1000,10001);
% datas = dir(fullfile(directory, '*.mat'));
% filename = datas.name;
% for idx = 1:numel(datas)
%     filename = datas(idx).name;
%     load(fullfile(directory, filename))
% 
%     figure(idx);
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0.3, 1, 0.3], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
% 
%     legend({'$|\tilde{f}_{exact}(t)|$', '$|\tilde{f}_{estimated}(t)|$','$\omega_{exact}(t)$', '$\omega_{estimated}(t)$'}, 'Location', 'best', 'Interpreter', 'latex'); 
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

%% ======== All Plots for Reconstructing the 9th + 11th harmonics with 9th + 11th harmonic gaussians ======== %% 
% figureNames = {'9h_11h_2g_c'  '9h_11h_4g_c_windowing_blur' '9h_11h_6g_c_windowing_blur'...
%                '9h_11h_2g_nc' '9h_11h_4g_nc_windowing' '9h_11h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(2,3);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

%% ======== All Plots for Reconstructing the 9th + 11th harmonic with 10th harmonic gaussians ======== %% 
% figureNames = {'10h_2g_c' '10h_4g_c' '10h_6g_c_windowing_blur'...
%                '10h_2g_nc' '10h_4g_nc' '10h_6g_nc_windowing_blur'};
% figure('units','normalized','outerposition',[0 0 1 1])
% t = tiledlayout(2,3);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

%% ======== All Plots for Reconstructing the 9th + 11th harmonics with 9th + 10th + 11th harmonic gaussians ======== %% 
figureNames = {'9h_10h_11h_3g_c_windowing_blur'  '9h_10h_11h_6g_c_windowing_blur'...
               '9h_10h_11h_3g_nc_windowing_blur' '9h_10h_11h_6g_nc_windowing_blur'};
figure()
t = tiledlayout(2,2);
for idx = 1:numel(figureNames)
    load(fullfile(directory, figureNames{idx}))

    nexttile
    yyaxis left
    plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
    plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
    ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');

    yyaxis right 
    plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
    plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
    ylabel('\omega(t) [a.u.]');
    xlabel('Time [a.u.]', 'Interpreter', 'latex');

    title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
    xlim([-600 600]); grid on;
end

%% ====== Comparing Chirp vs Non-chirp ======
% figureNames = {'9h_11h_6g_c_windowing_blur' '9h_11h_6g_nc_windowing_blur' '9h_10h_11h_6g_c_windowing_blur'...
%     '9h_10h_11h_6g_nc_windowing_blur' '10h_6g_c_windowing_blur' '10h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(3,2);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', figureNames{idx}, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

% figureNames = {'9h_11h_6g_c_windowing_blur' '9h_11h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(2,1);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', figureNames{idx}, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

%% ====== Comparing 9th vs 9th and 11th vs 9th, 10th, and 11th ======
% figureNames = {'10h_6g_c_windowing_blur' '9h_11h_6g_c_windowing_blur'  '9h_10h_11h_6g_c_windowing_blur'...
%                '10h_6g_nc_windowing_blur' '9h_11h_6g_nc_windowing_blur' '9h_10h_11h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(2,3);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

%% ====== Different Gaussians per color ====== %%
% figureNames = {'9h_6g_nc_windowing_blur' '9h_11h_5gx1g_nc_windowing_blur' '9h_11h_4gx2g_nc_windowing_blur'...
%                '9h_11h_2gx4g_nc_windowing_blur' '9h_11h_1gx5g_nc_windowing_blur' '11h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(2,3);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

% figureNames = {'9h_6g_c_windowing_blur' '9h_11h_5gx1g_c_windowing_blur' '9h_11h_4gx2g_c_windowing_blur'...
%                '9h_11h_2gx4g_c_windowing_blur' '9h_11h_1gx5g_c_windowing_blur' '11h_6g_c_windowing_blur'};
% figure()
% t = tiledlayout(2,3);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

% figureNames = {'9h_6g_c_windowing_blur' '10h_6g_c_windowing_blur' '11h_6g_c_windowing_blur'...
%                '9h_6g_nc_windowing_blur' '10h_6g_nc_windowing_blur' '11h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(2,3);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end

%% ===== Best Results ===== 
% figureNames = {'9h_11h_2gx4g_c_windowing_blur' '9h_11h_4gx2g_nc_windowing_blur' '9h_11h_6g_nc_windowing_blur'};
% figure()
% t = tiledlayout(3,1);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', name_plot, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end


%% ====== Comparing Methods ====== 
% figureNames = {'9h_11h_2g_c_overload_windowing_blur_f' '9h_11h_2g_nc'};
% figure()
% t = tiledlayout(2,1);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', figureNames{idx}, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end
% 
% figureNames = {'9h_11h_6g_c_overload_4g_c_f' '9h_11h_4g_nc_windowing'};
% figure()
% t = tiledlayout(2,1);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', figureNames{idx}, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end
% 
% figureNames = {'9h_11h_6g_nc_overload_windowing_blur_f' '9h_11h_6g_c_windowing_blur'};
% figure()
% t = tiledlayout(2,1);
% for idx = 1:numel(figureNames)
%     load(fullfile(directory, figureNames{idx}))
% 
%     nexttile
%     yyaxis left
%     plot(time + center_graph, abs(experiment_vals), '-', 'Color', [0.4, 0.6, 1], 'LineWidth', 2); hold on;
%     plot(time + center_graph, abs(estimated_vals), '-', 'Color', [0, 0.8, 0], 'LineWidth', 1.5); 
%     ylabel('$|\tilde{f}(t)|$ [a.u.]', 'Interpreter', 'latex');
% 
%     yyaxis right 
%     plot(time + center_graph, experiment_omega,'--', 'Color', [1, 0.6, 0.2], 'LineWidth', 1); 
%     plot(time + center_graph, estimated_omega, '-','Color', [1, 0.3, 0], 'LineWidth', 1); 
%     ylabel('\omega(t) [a.u.]');
%     xlabel('Time [a.u.]', 'Interpreter', 'latex');
% 
%     title(sprintf('Pulse Comparison: %s\nField Error: %.10f', figureNames{idx}, field_error), 'Interpreter', 'none');
%     xlim([-600 600]); grid on;
% end





