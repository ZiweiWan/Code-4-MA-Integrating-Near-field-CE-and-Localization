function [theta_refined, phi_refined] = refine_theta_phi(...
    Y_residual, theta_coarse, phi_coarse, ...
    r_fixed, ant_pos, lambda_c, ...
    max_iters, convergence_tol)

% --- 1. 初始化 ---
P = [theta_coarse; phi_coarse];
k_c = 2 * pi / lambda_c;

% --- 2. 迭代精炼循环 ---
for iter = 1:max_iters
    theta_curr = P(1);
    phi_curr   = P(2);

    % --- A. 构建原子及其一阶和二阶导数 ---
    [a, Ja, Ha] = get_spatial_derivatives(theta_curr, phi_curr, r_fixed, ant_pos, k_c);
    
    % --- B. 计算梯度 ---
    x_proj = a' * Y_residual; % 1 x K 投影向量
    g_theta = 2 * real(x_proj * Y_residual' * Ja.theta);
    g_phi   = 2 * real(x_proj * Y_residual' * Ja.phi);
    g = [g_theta; g_phi];
    
    % --- C. 计算完整的海森矩阵 H = H1 + H2 ---
    
    % C1. 计算第一项 (近似项)
    YJ_theta = Y_residual' * Ja.theta;
    YJ_phi   = Y_residual' * Ja.phi;
    
    H1_11 = 2 * real(YJ_theta' * YJ_theta);
    H1_22 = 2 * real(YJ_phi'   * YJ_phi);
    H1_12 = 2 * real(YJ_theta' * YJ_phi);
    H1 = [H1_11, H1_12; H1_12, H1_22];
    
    % C2. 计算第二项 (二阶修正项)
    YH_tt = Y_residual' * Ha.tt;
    YH_pp = Y_residual' * Ha.pp;
    YH_tp = Y_residual' * Ha.tp;
    
    H2_11 = 2 * real(x_proj * YH_tt);
    H2_22 = 2 * real(x_proj * YH_pp);
    H2_12 = 2 * real(x_proj * YH_tp);
    H2 = [H2_11, H2_12; H2_12, H2_22];
    
    % 完整的海森矩阵
    H = H1 + H2;
    
    % --- D. 求解牛顿方向并更新 ---
    % 增加正则化项以提高数值稳定性
    delta_P = -(H + 1e-6 * eye(2)) \ g;
    P = P + delta_P;

    % --- E. 检查收敛 ---
    if norm(delta_P) < convergence_tol
        break;
    end
end

% --- 3. 输出最终结果 ---
theta_refined = P(1);
phi_refined   = P(2);

end

%% 
function [a, Ja, Ha] = get_spatial_derivatives(theta, phi, r, ant_pos, k_c)
    s_cart = [r * sin(theta) * cos(phi); r * sin(theta) * sin(phi); r * cos(theta)];
    diff_mat = ant_pos - s_cart;
    dist_vec = sqrt(sum(diff_mat.^2, 1))';
    
    a = exp(-1j * k_c * (dist_vec - r));
    
    % --- 一阶导数相关 ---
    ds_dtheta = [r*cos(theta)*cos(phi); r*cos(theta)*sin(phi); -r*sin(theta)];
    ds_dphi   = [-r*sin(theta)*sin(phi); r*sin(theta)*cos(phi); 0];
    
    ddist_dtheta = -(diff_mat' * ds_dtheta) ./ dist_vec;
    ddist_dphi   = -(diff_mat' * ds_dphi) ./ dist_vec;
    
    common_factor = -1j * k_c * a;
    Ja.theta = common_factor .* ddist_dtheta;
    Ja.phi   = common_factor .* ddist_dphi;
    
    % --- 二阶导数相关 ---
    d2s_dtheta2 = -s_cart;
    d2s_dphi2   = [-r*sin(theta)*cos(phi); -r*sin(theta)*sin(phi); 0];
    d2s_dthetadphi = [-r*cos(theta)*sin(phi); r*cos(theta)*cos(phi); 0];
    
    d2dist_dtheta2 = (-(ds_dtheta'*ds_dtheta) - (diff_mat'*d2s_dtheta2) + ddist_dtheta.^2) ./ dist_vec;
    d2dist_dphi2   = (-(ds_dphi'*ds_dphi) - (diff_mat'*d2s_dphi2) + ddist_dphi.^2) ./ dist_vec;
    d2dist_dthetadphi = (-(ds_dtheta'*ds_dphi) - (diff_mat'*d2s_dthetadphi) + ddist_dtheta.*ddist_dphi) ./ dist_vec;
    
    Ha.tt = common_factor .* (d2dist_dtheta2 - 1j*k_c*ddist_dtheta.^2);
    Ha.pp = common_factor .* (d2dist_dphi2 - 1j*k_c*ddist_dphi.^2);
    Ha.tp = common_factor .* (d2dist_dthetadphi - 1j*k_c*ddist_dtheta.*ddist_dphi);
end
