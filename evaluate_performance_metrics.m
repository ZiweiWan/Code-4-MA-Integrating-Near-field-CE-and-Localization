function metrics = evaluate_performance_metrics(true_params, s_est, angle_matching_threshold_deg, ospa_c)
%% --- 1. 初始化 ---
metrics = struct(...
    'angle_rmse_deg', inf, 'angle_mae_deg', inf, ...
    'r_rmse', inf, 'r_mae', inf, ...
    'ospa', inf);

if isempty(s_est) && isempty(true_params.r)
    metrics = struct(...
        'angle_rmse_deg', 0, 'angle_mae_deg', 0, ...
        'r_rmse', 0, 'r_mae', 0, ...
        'ospa', 0);
    return;
end

%% --- 2. 参数转换 ---
s_true = [true_params.r .* sin(true_params.theta) .* cos(true_params.phi);
          true_params.r .* sin(true_params.theta) .* sin(true_params.phi);
          true_params.r .* cos(true_params.theta)];
      
if isempty(s_est) || isempty(s_true)
    metrics.ospa = calculate_ospa_p1(s_true, s_est, ospa_c);
    return; 
end
      
u_true = s_true ./ vecnorm(s_true);
u_est = s_est ./ vecnorm(s_est);
r_est = vecnorm(s_est);
num_true = size(s_true, 2);
num_est = size(s_est, 2);
angle_matching_threshold_rad = deg2rad(angle_matching_threshold_deg);

%% --- 3. RMSE/MAE ---
all_angle_errors_rad = [];
all_range_errors = [];
for i = 1:num_true
    u_true_i = u_true(:, i);
    r_true_i = true_params.r(i);
    
    dot_products = u_true_i' * u_est;
    dot_products = min(max(dot_products, -1), 1);
    angular_dists_rad = acos(dot_products);
    matched_indices = find(angular_dists_rad <= angle_matching_threshold_rad);
    
    if isempty(matched_indices)
        continue;
    end
    
    [min_angle_dist_rad, idx_in_matched] = min(angular_dists_rad(matched_indices));
    closest_idx = matched_indices(idx_in_matched);
    
    angle_error_rad = min_angle_dist_rad;
    all_angle_errors_rad = [all_angle_errors_rad, angle_error_rad];

    r_est_closest = r_est(closest_idx);
    range_error = r_true_i - r_est_closest;
    all_range_errors = [all_range_errors, range_error];
end

if ~isempty(all_angle_errors_rad)
    metrics.angle_rmse_deg = rad2deg(sqrt(mean(all_angle_errors_rad.^2)));
    metrics.angle_mae_deg  = rad2deg(mean(abs(all_angle_errors_rad)));
    metrics.r_rmse = sqrt(mean(all_range_errors.^2));
    metrics.r_mae  = mean(abs(all_range_errors));
end

%% --- 4. 计算OSPA指标 ---
metrics.ospa = calculate_ospa_p1(s_true, s_est, ospa_c);

end

%% --- OSPA(p=1) ---
function ospa_val = calculate_ospa_p1(X, Y, c)
    % 输入 X, Y 是 3xM 和 3xN 的笛卡尔坐标矩阵
    m = size(X, 2);
    n = size(Y, 2);

    if m == 0 || n == 0
        % 如果一个集合为空，误差就是所有未匹配点的惩罚之和，再除以总数
        ospa_val = c;
        return;
    end
    
    % OSPA定义要求 m <= n, 如果不是则交换
    if m > n
        [X, Y, m, n] = deal(Y, X, n, m);
    end

    % 1. 计算截断后的距离/代价矩阵
    dist_matrix = pdist2(X', Y');
    cost_matrix = min(c, dist_matrix);

    % 2. 求解分配问题以获得最优定位代价
    [assignment, ~, ~] = matchpairs(cost_matrix, c * 100);
    

    loc_cost = sum(cost_matrix(sub2ind(size(cost_matrix), assignment(:,1), assignment(:,2))));


    % 3. 计算基数代价
    card_cost = c * (n - m);

    % 4. 合并并归一化
    ospa_val = (loc_cost + card_cost) / n;
end
