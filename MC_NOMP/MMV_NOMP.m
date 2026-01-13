function [estimated_paths,D_estimated,X_est] = MMV_NOMP(Y_observed, D, dict_params, ...
    r_fixed, ant_pos, lambda_c, ...
    L_est, refine_max_iters, refine_conv_tol)


% --- 1. 初始化 ---
Y_residual = Y_observed;
estimated_paths = cell(L_est, 1);
D_estimated = zeros(size(D, 1), L_est);

% --- 2. 主循环 ---
for l = 1:L_est
    proj_matrix = D' * Y_residual;
    
    gain_sq = sum(abs(proj_matrix).^2, 2);
    
    [~, best_idx] = max(gain_sq);
    
    % 提取粗略参数
    [idx_t, idx_p] = ind2sub([length(dict_params.theta_grid), ...
                               length(dict_params.phi_grid)], best_idx);
    theta_coarse = dict_params.theta_grid(idx_t);
    phi_coarse = dict_params.phi_grid(idx_p);
    
    [theta_ref, phi_ref] = refine_theta_phi(...
        Y_residual, theta_coarse, phi_coarse, ...
        r_fixed, ant_pos, lambda_c, ...
        refine_max_iters, refine_conv_tol);
    % theta_ref = theta_coarse;phi_ref = phi_coarse;
        
    % --- 参数边界处理 ---
    phi_ref = mod(phi_ref, 2*pi); 
    if phi_ref > pi, phi_ref = 2*pi - phi_ref; end
        
    s_cart = r_fixed * [sin(theta_ref)*cos(phi_ref); sin(theta_ref)*sin(phi_ref); cos(theta_ref)];
    diff_mat = ant_pos - s_cart;
    dist_vec = sqrt(sum(diff_mat.^2, 1))';
    d_new = exp(-1j * (2*pi/lambda_c) * (dist_vec - r_fixed));
    D_estimated(:, l) = d_new;
    
    path.theta = theta_ref;
    path.phi = phi_ref;
    estimated_paths{l} = path;
    
    X_est = D_estimated(:, 1:l) \ Y_observed; % X_est 是 l x K 矩阵
    
    Y_residual = Y_observed - D_estimated(:, 1:l) * X_est;
end
end
