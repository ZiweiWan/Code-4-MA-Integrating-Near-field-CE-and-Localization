function [D, dict_params] = build_dictionary_full(theta_grid, phi_grid, r_grid, ant_pos, lambda_c)

% --- 1. 初始化参数 ---
G_theta = length(theta_grid);
G_phi = length(phi_grid);
G_r = length(r_grid);
G_total = G_theta * G_phi * G_r;
N = size(ant_pos, 2);
k_c = 2 * pi / lambda_c;

dict_params.theta_grid = theta_grid;
dict_params.phi_grid = phi_grid;
dict_params.r_grid = r_grid;
dict_params.size = [G_theta, G_phi, G_r];

% --- 2. 生成所有网格点组合 ---
[T_theta, T_phi, T_r] = ndgrid(theta_grid, phi_grid, r_grid);
theta_vec = T_theta(:);
phi_vec = T_phi(:);
r_vec = T_r(:);

% --- 3. 向量化计算所有空间导向矢量 ---
% 计算所有网格点对应的散射体三维坐标
s_cart_all = [r_vec .* sin(theta_vec) .* cos(phi_vec), ...
              r_vec .* sin(theta_vec) .* sin(phi_vec), ...
              r_vec .* cos(theta_vec)];

% 向量化计算天线到每个散射体的距离矩阵
ant_pos_sq_sum = sum(ant_pos.^2, 1);      % 尺寸: 1 x N
s_cart_sq_sum = sum(s_cart_all.^2, 2);   % 尺寸: G_total x 1
cross_term = 2 * s_cart_all * ant_pos;   % 尺寸: G_total x N

dist_sq_matrix = s_cart_sq_sum + ant_pos_sq_sum - cross_term;
dist_matrix = sqrt(max(dist_sq_matrix, 0)); % 尺寸: G_total x N

% --- 4. 生成最终字典 ---
D = exp(-1j * k_c * (dist_matrix.' - r_vec.'));

end
