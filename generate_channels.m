%% 1. 初始化环境
clear; clc; close all;
warning off;
addpath(genpath(pwd));

%% 2. 参数定义
Niter = 5000;
L = 6;
output_file = 'pregenerated_channel_params_L6.mat';

N_full_side = 32;
fc = 10e9; c = 3e8; lambda_c = c / fc;
ant_area_side = N_full_side * 0.5 * lambda_c; Z = 2 * ant_area_side^2 / lambda_c;
rmin = 2; rmax = 0.8 * Z;

%% 
[x_idx, z_idx] = meshgrid(0:N_full_side-1, 0:N_full_side-1);
ant_pos_full = [x_idx(:)'; zeros(1, size(x_idx(:),1)); z_idx(:)'] * 0.5 * lambda_c;
ant_pos_full = ant_pos_full - mean(ant_pos_full, 2);

%% 
channel_params_dataset = cell(Niter, 1);
rng(2024); % 设置固定的随机种子以保证可复现性

fprintf('--- 开始生成 %d 个信道参数实例 (固定 L=%d) ---\n', Niter, L);
tic;
for iter = 1:Niter
    % 生成角度参数，确保路径之间有最小间隔
    theta = [];
    phi = [];
    while length(theta) < L
        th_cand = 2*pi/3 * rand(1, 1) + pi/6;
        ph_cand = 2*pi/3 * rand(1, 1) + pi/6;
        if isempty(theta) || (all(abs(theta - th_cand) > deg2rad(5)) && all(abs(phi - ph_cand) > deg2rad(5)))
            theta(end+1) = th_cand;
            phi(end+1) = ph_cand;
        end
    end
    
    % 生成其他参数
    r = (rmax - rmin) * rand(1, L) + rmin;
    tau = ((rmax+10) / c) * rand(1, L) + (r/c);
    beta = (randn(1, L) + 1j * randn(1, L)) / sqrt(2);
    
    % 将所有参数存入一个结构体，并显式包含L
    channel_inst.L = L;
    channel_inst.r = r;
    channel_inst.theta = theta;
    channel_inst.phi = phi;
    channel_inst.tau = tau;
    channel_inst.beta = beta;
    
    channel_params_dataset{iter} = channel_inst;
end
generation_time = toc;

%% 5. 保存到文件
fprintf('--- 正在保存数据到 %s ---\n', output_file);
save(output_file, 'channel_params_dataset', 'ant_pos_full', '-v7.3');
fprintf('--- 数据保存成功 ---\n');

