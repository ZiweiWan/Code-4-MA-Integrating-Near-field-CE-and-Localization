function [nmse_initial, nmse_omp, nmse_isac, nmse_omp_rmax, metrics_initial, metrics_omp, metrics_isac, metrics_omp_rmax] = ...
    main_comp(channel_params, snr_dB, ant_pos_full, num_div_x, num_div_z,antenna_measurement_ratio,subcarrier_measurement_ratio)
%% --- 1. 参数设置与默认值 ---
G_r_omp = 30;

N_full_side = sqrt(size(ant_pos_full, 2));
N_total = N_full_side^2;
fc = 10e9; c = 3e8; lambda_c = c / fc; k_c = 2 * pi / lambda_c;
K_full = 64; subcarrier_spacing = 200e3;
ant_area_side = N_full_side * 0.5 * lambda_c; Z = 2 * ant_area_side^2 / lambda_c;
rmin = 2; rmax = 0.8 * Z;
subcarrier_freqs_full = (-(K_full-1)/2:(K_full-1)/2)' * subcarrier_spacing + fc;
L = length(channel_params.r);
L_est_sensing = 10;L_est_initial = 10;L_est_omp_full = 10;
refine_max_iters = 10;
refine_conv_tol = 1e-5;
% antenna_measurement_ratio = 0.5;
% subcarrier_measurement_ratio = 0.5;
gain_filter_threshold_ratio = 0.1;
angle_matching_threshold_deg = 10.0; 

ospa_c = 3.0;

regularization_param = 1e-3;

%% 步骤 1: 根据参数生成信道并接收信号
H_full = zeros(N_total, K_full);
for l = 1:L
    s_l = [channel_params.r(l) * sin(channel_params.theta(l)) * cos(channel_params.phi(l));
           channel_params.r(l) * sin(channel_params.theta(l)) * sin(channel_params.phi(l));
           channel_params.r(l) * cos(channel_params.theta(l))];
    dist_vec = sqrt(sum((ant_pos_full - s_l).^2, 1))';
    a_l = exp(-1j * k_c * (dist_vec - channel_params.r(l)));
    f_l = exp(-1j * 2 * pi * subcarrier_freqs_full * channel_params.tau(l));
    H_full = H_full + channel_params.beta(l) * (a_l * f_l.');
end
signal_power = norm(H_full, 'fro')^2 / (N_total * K_full);
noise_power = signal_power / (10^(snr_dB/10));
noise_matrix = sqrt(noise_power/2) * (randn(N_total, K_full) + 1j * randn(N_total, K_full));
Y_full = H_full + noise_matrix;
num_ant_measurements_global = round(N_total * antenna_measurement_ratio);
sampled_ant_indices_global = sort(randperm(N_total, num_ant_measurements_global));
num_sc_measurements = round(K_full * subcarrier_measurement_ratio);
sampled_sc_indices = sort(randperm(K_full, num_sc_measurements));
Y_doubly_sampled_global = Y_full(sampled_ant_indices_global, sampled_sc_indices);
ant_pos_sampled_global = ant_pos_full(:, sampled_ant_indices_global);

%% 步骤 2A: 方法一 (OMP-2D + 几何定位, 无精炼)
num_sub_arrays = num_div_x * num_div_z;
sub_arrays_initial = cell(num_sub_arrays, 1); count = 1;
[x_idx, z_idx] = meshgrid(0:N_full_side-1, 0:N_full_side-1);
N_sub_side_x = N_full_side / num_div_x; N_sub_side_z = N_full_side / num_div_z;
for j = 0:(num_div_z - 1)
    for i = 0:(num_div_x - 1)
        x_start = i * N_sub_side_x; z_start = j * N_sub_side_z;
        mask = (x_idx(:) >= x_start) & (x_idx(:) < x_start + N_sub_side_x) & (z_idx(:) >= z_start) & (z_idx(:) < z_start + N_sub_side_z);
        indices = find(mask);
        ant_pos_sub_global = ant_pos_full(:, indices);
        sub_array_info.indices = indices;
        sub_array_info.ant_pos_global = ant_pos_sub_global;
        sub_array_info.center_pos = mean(ant_pos_sub_global, 2);
        sub_arrays_initial{count} = sub_array_info;
        count = count + 1;
    end
end
G_theta = 60; G_phi = 30;
theta_grid = linspace(pi/6, pi/2 + pi/3, G_theta);
phi_grid = linspace(pi/6, 5*pi/6, G_phi);
for m = 1:num_sub_arrays
    Y_sub_full = Y_full(sub_arrays_initial{m}.indices, :);
    ant_pos_local_full = ant_pos_full(:, sub_arrays_initial{m}.indices) - sub_arrays_initial{m}.center_pos;
    N_sub = size(Y_sub_full, 1);
    num_ant_measurements_sub = round(N_sub * antenna_measurement_ratio);
    sampled_ant_indices_sub = randperm(N_sub, num_ant_measurements_sub);
    Y_sub_doubly_sampled = Y_sub_full(sampled_ant_indices_sub, sampled_sc_indices);
    ant_pos_local_sampled = ant_pos_local_full(:, sampled_ant_indices_sub);
    [D_angle_sub, ~] = build_dictionary_angle(theta_grid, phi_grid, rmax, ant_pos_local_sampled, lambda_c);
    R = Y_sub_doubly_sampled;
    estimated_paths_coarse = cell(L_est_initial, 1);
    D_estimated_sub = zeros(size(D_angle_sub, 1), L_est_initial);
    for l = 1:L_est_initial
        proj_matrix = D_angle_sub' * R;
        gain_sq = sum(abs(proj_matrix).^2, 2);
        [~, best_idx] = max(gain_sq);
        [idx_t, idx_p] = ind2sub([G_theta, G_phi], best_idx);
        path.theta = theta_grid(idx_t);
        path.phi = phi_grid(idx_p);
        estimated_paths_coarse{l} = path;
        D_estimated_sub(:, l) = D_angle_sub(:, best_idx);
        X_temp = D_estimated_sub(:, 1:l) \ Y_sub_doubly_sampled;
        R = Y_sub_doubly_sampled - D_estimated_sub(:, 1:l) * X_temp;
    end
    sub_arrays_initial{m}.est_paths = estimated_paths_coarse;
end
max_geometric_error = 0.1;
clustering_angle_thresh_deg = 10.0;
min_cluster_size = 2;
scatterer_pos_est_initial = [];
for l = 1:L_est_initial
    all_vectors = zeros(3, num_sub_arrays);
    for m = 1:num_sub_arrays
        path = sub_arrays_initial{m}.est_paths{l};
        all_vectors(:, m) = [sin(path.theta)*cos(path.phi); sin(path.theta)*sin(path.phi); cos(path.theta)];
    end
    angular_distances_deg = rad2deg(acos(min(max(all_vectors' * all_vectors, -1), 1)));
    unassigned_indices = 1:num_sub_arrays;
    clusters = {};
    while ~isempty(unassigned_indices)
        seed_idx = unassigned_indices(1);
        new_cluster = seed_idx;
        unassigned_indices(1) = [];
        neighbors = unassigned_indices(angular_distances_deg(seed_idx, unassigned_indices) < clustering_angle_thresh_deg);
        new_cluster = [new_cluster, neighbors];
        unassigned_indices = setdiff(unassigned_indices, neighbors);
        clusters{end+1} = new_cluster;
    end
    for c_idx = 1:length(clusters)
        current_cluster = clusters{c_idx};
        pos_est = [NaN; NaN; NaN];
        if length(current_cluster) >= min_cluster_size
            A = zeros(3, 3);
            b = zeros(3, 1);
            for m_idx = current_cluster
                v_m = all_vectors(:, m_idx);
                c_m_global = sub_arrays_initial{m_idx}.center_pos;
                P_m = eye(3) - v_m * v_m';
                A = A + P_m;
                b = b + P_m * c_m_global;
            end
            if cond(A) > 1e10
                pos_est = [NaN; NaN; NaN]; 
            else
                pos_est = A \ b;
            end
        end
        is_valid = ~any(isnan(pos_est));
        if is_valid
            geometric_error = 0;
            for m_idx = current_cluster
                v_m = all_vectors(:, m_idx);
                c_m_global = sub_arrays_initial{m_idx}.center_pos;
                geometric_error = geometric_error + norm((eye(3) - v_m*v_m') * (pos_est - c_m_global))^2;
            end
            if geometric_error > max_geometric_error || norm(pos_est) < rmin/2 || norm(pos_est) > rmax*1.5 || pos_est(2) < 0
                is_valid = false;
            end
        end
        if is_valid
            scatterer_pos_est_initial = [scatterer_pos_est_initial, pos_est];
        end
    end
end


% 2A.4. 基于定位结果的最终信道估计
L_valid_raw_initial = size(scatterer_pos_est_initial, 2);
H_initial = zeros(N_total, K_full);
tau_grid_fine = linspace(0, 2 * (rmax + 10) / c, 400);
freqs_sampled = subcarrier_freqs_full(sampled_sc_indices);
s_initial_filtered = []; 
if L_valid_raw_initial > 0
    A_est_initial = zeros(num_ant_measurements_global, L_valid_raw_initial);
    for l = 1:L_valid_raw_initial
        s_est_l = scatterer_pos_est_initial(:, l);
        dist_vec = sqrt(sum((ant_pos_sampled_global - s_est_l).^2, 1))';
        r_est_l = norm(s_est_l);
        A_est_initial(:, l) = exp(-1j * k_c * (dist_vec - r_est_l));
    end
    
    if cond(A_est_initial) > 1e3
        % 矩阵病态，使用正则化方法
        lambda = regularization_param;
        I = eye(L_valid_raw_initial);
        A_prime_A = A_est_initial' * A_est_initial;
        A_prime_Y = A_est_initial' * Y_doubly_sampled_global;
        X_est_sampled_initial = (A_prime_A + lambda * I) \ A_prime_Y;
    else
        % 矩阵良态，使用标准最小二乘法
        X_est_sampled_initial = A_est_initial \ Y_doubly_sampled_global;
    end
    
    beta_est_initial = zeros(L_valid_raw_initial, 1);
    tau_est_initial = zeros(L_valid_raw_initial, 1);
    for l = 1:L_valid_raw_initial
        beta_l_sampled = X_est_sampled_initial(l, :).';
        errors = zeros(length(tau_grid_fine), 1);
        beta_candidates = zeros(length(tau_grid_fine), 1);
        for t_idx = 1:length(tau_grid_fine)
            tau_g = tau_grid_fine(t_idx);
            f_g = exp(-1j * 2 * pi * freqs_sampled * tau_g);
            beta_g = f_g \ beta_l_sampled;
            errors(t_idx) = norm(beta_l_sampled - beta_g * f_g)^2;
            beta_candidates(t_idx) = beta_g;
        end
        [~, best_t_idx] = min(errors);
        tau_est_initial(l) = tau_grid_fine(best_t_idx);
        beta_est_initial(l) = beta_candidates(best_t_idx);
    end
    
    beta_mag_initial = abs(beta_est_initial);
    [s_initial_filtered, valid_indices_initial] = prune_paths(beta_mag_initial, gain_filter_threshold_ratio, scatterer_pos_est_initial, L_est_sensing);
    
    for l_idx = 1:length(valid_indices_initial)
        l = valid_indices_initial(l_idx);
        s_est_l = s_initial_filtered(:, l_idx); 
        dist_vec_full = sqrt(sum((ant_pos_full - s_est_l).^2, 1))';
        r_est_l = norm(s_est_l);
        a_l_est_full = exp(-1j * k_c * (dist_vec_full - r_est_l));
        f_l_full = exp(-1j * 2 * pi * subcarrier_freqs_full * tau_est_initial(l));
        H_initial = H_initial + beta_est_initial(l) * (a_l_est_full * f_l_full.');
    end
end
nmse_initial = norm(H_full - H_initial, 'fro')^2 / norm(H_full, 'fro')^2;
metrics_initial = evaluate_performance_metrics(channel_params, s_initial_filtered, angle_matching_threshold_deg,  ospa_c);


%% 步骤 2B: 方法二 (OMP-3D)
r_grid_omp = linspace(rmin, rmax, G_r_omp);
[D_full, dict_params_full] = build_dictionary_full(theta_grid, phi_grid, r_grid_omp, ant_pos_sampled_global, lambda_c);
R_full = Y_doubly_sampled_global;
support_set_full = zeros(1, L_est_omp_full);
for l = 1:L_est_omp_full
    proj_matrix_full = D_full' * R_full;
    gain_sq_full = sum(abs(proj_matrix_full).^2, 2);
    [~, best_idx] = max(gain_sq_full);
    support_set_full(l) = best_idx;
    A_s_full = D_full(:, support_set_full(1:l));
    X_temp_full = A_s_full \ Y_doubly_sampled_global;
    R_full = Y_doubly_sampled_global - A_s_full * X_temp_full;
end
A_omp_sampled = D_full(:, support_set_full);
X_sampled_omp = A_omp_sampled \ Y_doubly_sampled_global;
beta_omp_full = zeros(L_est_omp_full, 1);
tau_omp_full = zeros(L_est_omp_full, 1);
for l = 1:L_est_omp_full
    beta_l_sampled = X_sampled_omp(l, :).';
    errors = zeros(length(tau_grid_fine), 1);
    beta_candidates = zeros(length(tau_grid_fine), 1);
    for t_idx = 1:length(tau_grid_fine)
        tau_g = tau_grid_fine(t_idx);
        f_g = exp(-1j * 2 * pi * freqs_sampled * tau_g);
        beta_g = f_g \ beta_l_sampled;
        errors(t_idx) = norm(beta_l_sampled - beta_g * f_g)^2;
        beta_candidates(t_idx) = beta_g;
    end
    [~, best_t_idx] = min(errors);
    tau_omp_full(l) = tau_grid_fine(best_t_idx);
    beta_omp_full(l) = beta_candidates(best_t_idx);
end
H_omp_full = zeros(N_total, K_full);
s_omp_full_est = zeros(3, L_est_omp_full);
for l = 1:L_est_omp_full
    s_omp_full_est(:, l) = get_scatterer_pos_from_index(support_set_full(l), dict_params_full);
end
beta_mag_omp = abs(beta_omp_full);
[s_omp_filtered, valid_indices_omp] = prune_paths(beta_mag_omp, gain_filter_threshold_ratio, s_omp_full_est, L_est_sensing);
for l_idx = 1:length(valid_indices_omp)
    l = valid_indices_omp(l_idx);
    s_est = s_omp_full_est(:, l);
    r_est = norm(s_est);
    dist_vec_full = sqrt(sum((ant_pos_full - s_est).^2, 1))';
    a_l_full = exp(-1j * k_c * (dist_vec_full - r_est));
    f_l_full = exp(-1j * 2 * pi * subcarrier_freqs_full * tau_omp_full(l));
    H_omp_full = H_omp_full + beta_omp_full(l) * (a_l_full * f_l_full.');
end
nmse_omp = norm(H_full - H_omp_full, 'fro')^2 / norm(H_full, 'fro')^2;
metrics_omp = evaluate_performance_metrics(channel_params, s_omp_filtered, angle_matching_threshold_deg,  ospa_c);

%% 步骤 3 & 4: 方法三 (NOMP + 几何定位, 带精炼)
sub_arrays = cell(num_sub_arrays, 1); count = 1;
for j = 0:(num_div_z - 1)
    for i = 0:(num_div_x - 1)
        x_start = i * N_sub_side_x; z_start = j * N_sub_side_z;
        mask = (x_idx(:) >= x_start) & (x_idx(:) < x_start + N_sub_side_x) & (z_idx(:) >= z_start) & (z_idx(:) < z_start + N_sub_side_z);
        indices = find(mask);
        ant_pos_sub_global = ant_pos_full(:, indices);
        sub_array_info.indices = indices;
        sub_array_info.ant_pos_global = ant_pos_sub_global;
        sub_array_info.center_pos = mean(ant_pos_sub_global, 2);
        sub_arrays{count} = sub_array_info;
        count = count + 1;
    end
end
for m = 1:num_sub_arrays
    Y_sub_full = Y_full(sub_arrays{m}.indices, :);
    ant_pos_local_full = ant_pos_full(:, sub_arrays{m}.indices) - sub_arrays{m}.center_pos;
    N_sub = size(Y_sub_full, 1);
    num_ant_measurements_sub = round(N_sub * antenna_measurement_ratio);
    sampled_ant_indices_sub = randperm(N_sub, num_ant_measurements_sub);
    Y_sub_doubly_sampled = Y_sub_full(sampled_ant_indices_sub, sampled_sc_indices);
    ant_pos_local_sampled = ant_pos_local_full(:, sampled_ant_indices_sub);
    [D_angle_sub, dict_params_sub] = build_dictionary_angle(theta_grid, phi_grid, rmax, ant_pos_local_sampled, lambda_c);
    estimated_paths = MMV_NOMP(Y_sub_doubly_sampled, D_angle_sub, dict_params_sub, rmax, ant_pos_local_sampled, lambda_c, L_est_sensing, refine_max_iters, refine_conv_tol);
    sub_arrays{m}.est_paths = estimated_paths;
end
scatterer_pos_est = [];
for l = 1:L_est_sensing
    all_vectors = zeros(3, num_sub_arrays);
    for m = 1:num_sub_arrays
        path = sub_arrays{m}.est_paths{l};
        all_vectors(:, m) = [sin(path.theta)*cos(path.phi); sin(path.theta)*sin(path.phi); cos(path.theta)];
    end
    angular_distances_deg = rad2deg(acos(min(max(all_vectors' * all_vectors, -1), 1)));
    unassigned_indices = 1:num_sub_arrays;
    clusters = {};
    while ~isempty(unassigned_indices)
        seed_idx = unassigned_indices(1);
        new_cluster = seed_idx;
        unassigned_indices(1) = [];
        neighbors = unassigned_indices(angular_distances_deg(seed_idx, unassigned_indices) < clustering_angle_thresh_deg);
        new_cluster = [new_cluster, neighbors];
        unassigned_indices = setdiff(unassigned_indices, neighbors);
        clusters{end+1} = new_cluster;
    end
    for c_idx = 1:length(clusters)
        current_cluster = clusters{c_idx};
        pos_est = [NaN; NaN; NaN];
        if length(current_cluster) >= min_cluster_size
            A = zeros(3, 3);
            b = zeros(3, 1);
            for m_idx = current_cluster
                v_m = all_vectors(:, m_idx);
                c_m_global = sub_arrays{m_idx}.center_pos;
                P_m = eye(3) - v_m * v_m';
                A = A + P_m;
                b = b + P_m * c_m_global;
            end
            if cond(A) > 1e10
                pos_est = [NaN; NaN; NaN];
            else
                pos_est = A \ b;
            end
        end
        is_valid = ~any(isnan(pos_est));
        if is_valid
            geometric_error = 0;
            for m_idx = current_cluster
                v_m = all_vectors(:, m_idx);
                c_m_global = sub_arrays{m_idx}.center_pos;
                geometric_error = geometric_error + norm((eye(3) - v_m*v_m') * (pos_est - c_m_global))^2;
            end
            if geometric_error > max_geometric_error || norm(pos_est) < rmin/2 || norm(pos_est) > rmax*1.5 || pos_est(2) < 0
                is_valid = false;
            end
        end
        if is_valid
            scatterer_pos_est = [scatterer_pos_est, pos_est];
        end
    end
end

L_valid_raw = size(scatterer_pos_est, 2);
H_est_final = zeros(N_total, K_full);
s_isac_filtered = [];
if L_valid_raw > 0
    A_est = zeros(num_ant_measurements_global, L_valid_raw);
    for l = 1:L_valid_raw
        s_est_l = scatterer_pos_est(:, l);
        dist_vec = sqrt(sum((ant_pos_sampled_global - s_est_l).^2, 1))';
        r_est_l = norm(s_est_l);
        A_est(:, l) = exp(-1j * k_c * (dist_vec - r_est_l));
    end
    
    if cond(A_est) > 1e3 
        % 矩阵病态，使用正则化方法
        lambda = regularization_param;
        I = eye(L_valid_raw);
        A_prime_A = A_est' * A_est;
        A_prime_Y = A_est' * Y_doubly_sampled_global;
        X_est_sampled_final = (A_prime_A + lambda * I) \ A_prime_Y;
    else
        % 矩阵良态，使用标准最小二乘法
        X_est_sampled_final = A_est \ Y_doubly_sampled_global;
    end
        
    beta_est_final = zeros(L_valid_raw, 1);
    tau_est_final = zeros(L_valid_raw, 1);
    for l = 1:L_valid_raw
        beta_l_sampled = X_est_sampled_final(l, :).';
        errors = zeros(length(tau_grid_fine), 1);
        beta_candidates = zeros(length(tau_grid_fine), 1);
        for t_idx = 1:length(tau_grid_fine)
            tau_g = tau_grid_fine(t_idx);
            f_g = exp(-1j * 2 * pi * freqs_sampled * tau_g);
            beta_g = f_g \ beta_l_sampled;
            errors(t_idx) = norm(beta_l_sampled - beta_g * f_g)^2;
            beta_candidates(t_idx) = beta_g;
        end
        [~, best_t_idx] = min(errors);
        tau_est_final(l) = tau_grid_fine(best_t_idx);
        beta_est_final(l) = beta_candidates(best_t_idx);
    end
    
    beta_mag_isac = abs(beta_est_final);
    [s_isac_filtered, valid_indices_isac] = prune_paths(beta_mag_isac, gain_filter_threshold_ratio, scatterer_pos_est, L_est_sensing);

    for l_idx = 1:length(valid_indices_isac)
        l = valid_indices_isac(l_idx);
        s_est_l = s_isac_filtered(:, l_idx); 
        dist_vec_full = sqrt(sum((ant_pos_full - s_est_l).^2, 1))';
        r_est_l = norm(s_est_l);
        a_l_est_full = exp(-1j * k_c * (dist_vec_full - r_est_l));
        f_l_full = exp(-1j * 2 * pi * subcarrier_freqs_full * tau_est_final(l));
        H_est_final = H_est_final + beta_est_final(l) * (a_l_est_full * f_l_full.');
    end
end
nmse_isac = norm(H_full - H_est_final, 'fro')^2 / norm(H_full, 'fro')^2;
metrics_isac = evaluate_performance_metrics(channel_params, s_isac_filtered, angle_matching_threshold_deg,  ospa_c);

%% 步骤 4: 方法4 (OMP-2D @ rmax)
% 1. 构建字典 (仅角度，r固定为rmax)
[D_angle_rmax, dict_params_angle] = build_dictionary_angle(theta_grid, phi_grid, rmax, ant_pos_sampled_global, lambda_c);

% 2. 执行OMP
R_rmax = Y_doubly_sampled_global;
support_set_rmax = zeros(1, L_est_omp_full);
for l = 1:L_est_omp_full
    proj_matrix_rmax = D_angle_rmax' * R_rmax;
    gain_sq_rmax = sum(abs(proj_matrix_rmax).^2, 2);
    [~, best_idx] = max(gain_sq_rmax);
    support_set_rmax(l) = best_idx;
    A_s_rmax = D_angle_rmax(:, support_set_rmax(1:l));
    X_temp_rmax = A_s_rmax \ Y_doubly_sampled_global;
    R_rmax = Y_doubly_sampled_global - A_s_rmax * X_temp_rmax;
end

% 3. 估计增益和时延
A_omp_rmax_sampled = D_angle_rmax(:, support_set_rmax);
X_sampled_omp_rmax = A_omp_rmax_sampled \ Y_doubly_sampled_global;
beta_omp_rmax = zeros(L_est_omp_full, 1);
tau_omp_rmax = zeros(L_est_omp_full, 1);
for l = 1:L_est_omp_full
    beta_l_sampled = X_sampled_omp_rmax(l, :).';
    errors = zeros(length(tau_grid_fine), 1);
    beta_candidates = zeros(length(tau_grid_fine), 1);
    for t_idx = 1:length(tau_grid_fine)
        tau_g = tau_grid_fine(t_idx);
        f_g = exp(-1j * 2 * pi * freqs_sampled * tau_g);
        beta_g = f_g \ beta_l_sampled;
        errors(t_idx) = norm(beta_l_sampled - beta_g * f_g)^2;
        beta_candidates(t_idx) = beta_g;
    end
    [~, best_t_idx] = min(errors);
    tau_omp_rmax(l) = tau_grid_fine(best_t_idx);
    beta_omp_rmax(l) = beta_candidates(best_t_idx);
end

% 4. 重建信道
H_omp_rmax = zeros(N_total, K_full);
s_omp_rmax_est = zeros(3, L_est_omp_full);
for l = 1:L_est_omp_full
    % 从2D字典索引恢复角度
    [idx_t, idx_p] = ind2sub([G_theta, G_phi], support_set_rmax(l));
    theta_est = dict_params_angle.theta_grid(idx_t);
    phi_est = dict_params_angle.phi_grid(idx_p);
    r_est = rmax; 
    s_omp_rmax_est(:, l) = r_est * [sin(theta_est)*cos(phi_est); sin(theta_est)*sin(phi_est); cos(theta_est)];
end

% 5. 剪枝和最终信道估计
beta_mag_omp_rmax = abs(beta_omp_rmax);
[s_omp_rmax_filtered, valid_indices_rmax] = prune_paths(beta_mag_omp_rmax, gain_filter_threshold_ratio, s_omp_rmax_est, L_est_sensing);

for l_idx = 1:length(valid_indices_rmax)
    l = valid_indices_rmax(l_idx);
    s_est_l = s_omp_rmax_est(:, l);
    r_est_l = norm(s_est_l); % 这将是rmax
    dist_vec_full = sqrt(sum((ant_pos_full - s_est_l).^2, 1))';
    a_l_full = exp(-1j * k_c * (dist_vec_full - r_est_l));
    f_l_full = exp(-1j * 2 * pi * subcarrier_freqs_full * tau_omp_rmax(l));
    H_omp_rmax = H_omp_rmax + beta_omp_rmax(l) * (a_l_full * f_l_full.');
end

% 6. 计算指标
nmse_omp_rmax = norm(H_full - H_omp_rmax, 'fro')^2 / norm(H_full, 'fro')^2;
metrics_omp_rmax = evaluate_performance_metrics(channel_params, s_omp_rmax_filtered, angle_matching_threshold_deg,  ospa_c);


end

%% --- 辅助函数 ---
function s_pos = get_scatterer_pos_from_index(dict_idx, dict_params)
    [idx_t, idx_p, idx_r] = ind2sub(dict_params.size, dict_idx);
    theta_est = dict_params.theta_grid(idx_t);
    phi_est = dict_params.phi_grid(idx_p);
    r_est = dict_params.r_grid(idx_r);
    s_pos = r_est * [sin(theta_est)*cos(phi_est); sin(theta_est)*sin(phi_est); cos(theta_est)];
end

function [s_filtered, valid_indices] = prune_paths(beta_mag, threshold_ratio, s_all, max_to_keep)
    s_filtered = [];
    if isempty(beta_mag), valid_indices = []; return; end
    max_gain = max(beta_mag);
    if max_gain == 0, valid_indices = []; return; end
    potential_indices = find(beta_mag >= max_gain * threshold_ratio);
    if isempty(potential_indices), valid_indices = []; return; end
    [~, sort_order] = sort(beta_mag(potential_indices), 'descend');
    sorted_indices = potential_indices(sort_order);
    if nargin == 4 && length(sorted_indices) > max_to_keep
        valid_indices = sorted_indices(1:max_to_keep);
    else
        valid_indices = sorted_indices;
    end
    if nargin >= 3 && ~isempty(s_all)
        s_filtered = s_all(:, valid_indices);
    end
end





