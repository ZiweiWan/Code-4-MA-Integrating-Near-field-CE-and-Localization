function [D, dict_params] = build_dictionary_angle(theta_grid, phi_grid, r_fixed, ant_pos, lambda_c)

% --- 1. 初始化参数 ---
G_theta = length(theta_grid);
G_phi = length(phi_grid);
G_total = G_theta * G_phi;
N = size(ant_pos, 2);
k_c = 2 * pi / lambda_c;

dict_params.theta_grid = theta_grid;
dict_params.phi_grid = phi_grid;

% --- 2. 生成所有网格点组合 ---
[T_theta, T_phi] = ndgrid(theta_grid, phi_grid);
theta_vec = T_theta(:);
phi_vec = T_phi(:);

% --- 3. 向量化计算所有空间导向矢量 ---
s_cart_all = [r_fixed * sin(theta_vec) .* cos(phi_vec), ...
              r_fixed * sin(theta_vec) .* sin(phi_vec), ...
              r_fixed * cos(theta_vec)];

ant_pos_sq_sum = sum(ant_pos.^2, 1);
s_cart_sq_sum = sum(s_cart_all.^2, 2);
cross_term = 2 * s_cart_all * ant_pos;

dist_sq_matrix = s_cart_sq_sum + ant_pos_sq_sum - cross_term;
dist_matrix = sqrt(max(dist_sq_matrix, 0));

D = exp(-1j * k_c * (dist_matrix.' - r_fixed));
end
