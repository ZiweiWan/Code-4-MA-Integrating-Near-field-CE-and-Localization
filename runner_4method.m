%% 1. 初始化环境
clear; clc; close all;
warning off; 

rng(2024);
addpath(genpath(pwd));

%% 2. 仿真参数设置
SNR_vec = -15:5:15;
Niter = 2000; 
antenna_measurement_ratio       =0.75;
subcarrier_measurement_ratio    =0.5;

channel_params_file = 'pregenerated_channel_params_L6.mat';
file_template = 'RE_All_5_Methods_L6_%g_%g.mat'; % 
output_results_file = sprintf(file_template, antenna_measurement_ratio, subcarrier_measurement_ratio);


%% 3. 加载预生成的信道数据
if ~exist(channel_params_file, 'file')
    fprintf('错误: 找不到信道文件 %s\n', channel_params_file);
    fprintf('请先运行 generate_channels.m (确保 L=6) 来生成信道数据。\n');
    return;
end
load(channel_params_file);

%% 4. 预分配结果存储
num_snr_points = length(SNR_vec);

% --- 方法一 (OMP-2D+Geo) ---
nmse_initial_all = zeros(num_snr_points, Niter);
angle_rmse_initial_all = zeros(num_snr_points, Niter);
angle_mae_initial_all  = zeros(num_snr_points, Niter);
r_rmse_initial_all     = zeros(num_snr_points, Niter);
r_mae_initial_all      = zeros(num_snr_points, Niter);
ospa_initial_all       = zeros(num_snr_points, Niter);

% --- 方法二 (OMP-3D) ---
nmse_omp_all     = zeros(num_snr_points, Niter);
angle_rmse_omp_all = zeros(num_snr_points, Niter);
angle_mae_omp_all  = zeros(num_snr_points, Niter);
r_rmse_omp_all     = zeros(num_snr_points, Niter);
r_mae_omp_all      = zeros(num_snr_points, Niter);
ospa_omp_all       = zeros(num_snr_points, Niter);

% --- 方法三 (NOMP+Geo) ---
nmse_isac_all    = zeros(num_snr_points, Niter);
angle_rmse_isac_all = zeros(num_snr_points, Niter);
angle_mae_isac_all  = zeros(num_snr_points, Niter);
r_rmse_isac_all     = zeros(num_snr_points, Niter);
r_mae_isac_all      = zeros(num_snr_points, Niter);
ospa_isac_all      = zeros(num_snr_points, Niter);

% --- 方法四 (OMP-2D @ rmax) ---
nmse_omp_rmax_all = zeros(num_snr_points, Niter);
angle_rmse_omp_rmax_all = zeros(num_snr_points, Niter);
angle_mae_omp_rmax_all  = zeros(num_snr_points, Niter);
r_rmse_omp_rmax_all     = zeros(num_snr_points, Niter);
r_mae_omp_rmax_all      = zeros(num_snr_points, Niter);
ospa_omp_rmax_all       = zeros(num_snr_points, Niter);


%% 5. 主仿真循环
global_tic = tic;

for i_SNR = 1:num_snr_points
    snr_dB = SNR_vec(i_SNR);
    snr_tic = tic;
    
    for iter = 1:Niter
        current_channel_params = channel_params_dataset{iter};
        
        [nmse_i, nmse_o, nmse_is, nmse_or,  ...
         metrics_i, metrics_o, metrics_is, metrics_or] = ...
            main_comp(current_channel_params, snr_dB, ant_pos_full, 2, 2,antenna_measurement_ratio,subcarrier_measurement_ratio);

        % --- 存储方法一 (OMP-2D+Geo) ---
        nmse_initial_all(i_SNR, iter) = nmse_i;
        angle_rmse_initial_all(i_SNR, iter) = metrics_i.angle_rmse_deg;
        angle_mae_initial_all(i_SNR, iter)  = metrics_i.angle_mae_deg;
        r_rmse_initial_all(i_SNR, iter)     = metrics_i.r_rmse;
        r_mae_initial_all(i_SNR, iter)      = metrics_i.r_mae;
        ospa_initial_all(i_SNR, iter)       = metrics_i.ospa;

        % --- 存储方法二 (OMP-3D) ---
        nmse_omp_all(i_SNR, iter)     = nmse_o;
        angle_rmse_omp_all(i_SNR, iter) = metrics_o.angle_rmse_deg;
        angle_mae_omp_all(i_SNR, iter)  = metrics_o.angle_mae_deg;
        r_rmse_omp_all(i_SNR, iter)     = metrics_o.r_rmse;
        r_mae_omp_all(i_SNR, iter)      = metrics_o.r_mae;
        ospa_omp_all(i_SNR, iter)       = metrics_o.ospa;

        % --- 存储方法三 (NOMP+Geo) ---
        nmse_isac_all(i_SNR, iter)    = nmse_is;
        angle_rmse_isac_all(i_SNR, iter) = metrics_is.angle_rmse_deg;
        angle_mae_isac_all(i_SNR, iter)  = metrics_is.angle_mae_deg;
        r_rmse_isac_all(i_SNR, iter)     = metrics_is.r_rmse;
        r_mae_isac_all(i_SNR, iter)      = metrics_is.r_mae;
        ospa_isac_all(i_SNR, iter)      = metrics_is.ospa;

        % --- 存储方法四 (OMP-2D @ rmax) ---
        nmse_omp_rmax_all(i_SNR, iter) = nmse_or;
        angle_rmse_omp_rmax_all(i_SNR, iter) = metrics_or.angle_rmse_deg;
        angle_mae_omp_rmax_all(i_SNR, iter)  = metrics_or.angle_mae_deg;
        r_rmse_omp_rmax_all(i_SNR, iter)     = metrics_or.r_rmse;
        r_mae_omp_rmax_all(i_SNR, iter)      = metrics_or.r_mae;
        ospa_omp_rmax_all(i_SNR, iter)       = metrics_or.ospa;
        
    end
    
    snr_time = toc(snr_tic);
    fprintf('  SNR = %3d dB 完成 | 耗时: %.1f 秒\n', snr_dB, snr_time);
end

total_time = toc(global_tic);
fprintf('\n--- 所有仿真完成，总耗时 %.2f 分钟 ---\n\n', total_time/60);

%% 6. 计算最终平均结果

all_results_matrices = {
    'nmse_initial_all', 'nmse_omp_all', 'nmse_isac_all', 'nmse_omp_rmax_all', ... 
    'angle_rmse_initial_all', 'angle_mae_initial_all', 'r_rmse_initial_all', 'r_mae_initial_all', 'ospa_initial_all', ...
    'angle_rmse_omp_all', 'angle_mae_omp_all', 'r_rmse_omp_all', 'r_mae_omp_all', 'ospa_omp_all', ...
    'angle_rmse_isac_all', 'angle_mae_isac_all', 'r_rmse_isac_all', 'r_mae_isac_all', 'ospa_isac_all', ...
    'angle_rmse_omp_rmax_all', 'angle_mae_omp_rmax_all', 'r_rmse_omp_rmax_all', 'r_mae_omp_rmax_all', 'ospa_omp_rmax_all' ...
};

for i = 1:length(all_results_matrices)
    mat_name = all_results_matrices{i};
    mat_data = eval(mat_name);
    mat_data(isinf(mat_data)) = NaN;
    assignin('base', mat_name, mat_data); 
end

% --- 通信 ---
avg_nmse_initial_db = 10*log10(mean(nmse_initial_all, 2, 'omitnan'));
avg_nmse_omp_db     = 10*log10(mean(nmse_omp_all, 2, 'omitnan'));
avg_nmse_isac_db    = 10*log10(mean(nmse_isac_all, 2, 'omitnan'));
avg_nmse_omp_rmax_db = 10*log10(mean(nmse_omp_rmax_all, 2, 'omitnan'));

% --- 感知 ---
% 方法一
avg_angle_rmse_initial = mean(angle_rmse_initial_all, 2, 'omitnan');
avg_angle_mae_initial  = mean(angle_mae_initial_all, 2, 'omitnan');
avg_r_rmse_initial     = mean(r_rmse_initial_all, 2, 'omitnan');
avg_r_mae_initial      = mean(r_mae_initial_all, 2, 'omitnan');
avg_ospa_initial       = mean(ospa_initial_all, 2, 'omitnan');
% 方法二
avg_angle_rmse_omp = mean(angle_rmse_omp_all, 2, 'omitnan');
avg_angle_mae_omp  = mean(angle_mae_omp_all, 2, 'omitnan');
avg_r_rmse_omp     = mean(r_rmse_omp_all, 2, 'omitnan');
avg_r_mae_omp      = mean(r_mae_omp_all, 2, 'omitnan');
avg_ospa_omp       = mean(ospa_omp_all, 2, 'omitnan');
% 方法三
avg_angle_rmse_isac = mean(angle_rmse_isac_all, 2, 'omitnan');
avg_angle_mae_isac  = mean(angle_mae_isac_all, 2, 'omitnan');
avg_r_rmse_isac     = mean(r_rmse_isac_all, 2, 'omitnan');
avg_r_mae_isac      = mean(r_mae_isac_all, 2, 'omitnan');
avg_ospa_isac      = mean(ospa_isac_all, 2, 'omitnan');
% 方法四
avg_angle_rmse_omp_rmax = mean(angle_rmse_omp_rmax_all, 2, 'omitnan');
avg_angle_mae_omp_rmax  = mean(angle_mae_omp_rmax_all, 2, 'omitnan');
avg_r_rmse_omp_rmax     = mean(r_rmse_omp_rmax_all, 2, 'omitnan');
avg_r_mae_omp_rmax      = mean(r_mae_omp_rmax_all, 2, 'omitnan');
avg_ospa_omp_rmax       = mean(ospa_omp_rmax_all, 2, 'omitnan');


fprintf('--- 正在保存所有结果到: %s ---\n', output_results_file);
save(output_results_file);