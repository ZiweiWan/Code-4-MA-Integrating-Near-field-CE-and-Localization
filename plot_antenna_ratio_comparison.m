%% 1. 初始化环境
clear; clc; close all;

%% 2. 参数设置
file_list = {
    'RE_All_5_Methods_L6_0.25_0.5.mat', ...
    'RE_All_5_Methods_L6_0.5_0.5.mat', ...
    'RE_All_5_Methods_L6_0.75_0.5.mat'
};

ratio_labels = {
    'port CR = 0.25', ...
    'port CR = 0.5', ...
    'port CR = 0.75'
};
ratio_styles = {
    '--s',  
    '-d',   
    '-.o'   
};

method_names = {
    'OMP-LSRC';
    'OMP 3D dictionary';
    'Proposed NOMP-LSRC';
    'OMP 2D dictionary'
};

method_colors = {
    [0 0.4470 0.7410],    
    [0.8500 0.3250 0.0980], 
    [0.4660 0.6740 0.1880], 
    [0.4940 0.1840 0.5560]  
};

nmse_data_fields = {
    'avg_nmse_initial_db', 
    'avg_nmse_omp_db', 
    'avg_nmse_isac_db', 
    'avg_nmse_omp_rmax_db'
};

ospa_data_fields = {
    'avg_ospa_initial', 
    'avg_ospa_omp', 
    'avg_ospa_isac'
};

font_size = 12;
legend_font_size = 7;
line_width = 2;
marker_size = 6;
num_ratios = length(file_list);

%% 3. 加载数据
all_data = cell(num_ratios, 1);
for r = 1:num_ratios
    file_name = file_list{r};
    if exist(file_name, 'file')
        all_data{r} = load(file_name);
    else
        error('错误: 找不到文件 %s。请确保该文件在当前路径下。', file_name);
    end
end

%% 4. 绘图: NMSE
figure('Name', 'sim1');
hold on; 
set(gca, 'YScale', 'log');

% 循环 4 种方法 (颜色)
for m = 1:4
    method_name = method_names{m};
    method_color = method_colors{m};
    data_field = nmse_data_fields{m};
    
    % 循环 3 种压缩率 (线型)
    for r = 1:num_ratios
        data = all_data{r};
        snr_vec = data.SNR_vec;
        data_to_plot = data.(data_field);
        
        plot_label = sprintf('%s (%s)', method_name, ratio_labels{r});
        
        semilogy(snr_vec, 10.^(data_to_plot/10), ...
             ratio_styles{r}, ...
            'Color', method_color, ...
            'LineWidth', line_width, ...
            'MarkerSize', marker_size, ...
            'DisplayName', plot_label);
    end
end

hold off;
grid on; 
xlabel('SNR (dB)'); 
ylabel('NMSE'); 
ylim([5e-3, 10]); 

lgd = legend('show', 'Location', 'north', 'NumColumns', 2); 
lgd.FontSize = legend_font_size;
set(gca, 'FontSize', font_size);
print(gcf, 'sim1', '-depsc', '-painters');
%% 5. 绘图: OSPA 
figure('Name', 'sim4');
hold on;

% 循环 3 种方法 (颜色)
for m = 1:3
    method_name = method_names{m};
    method_color = method_colors{m};
    data_field = ospa_data_fields{m};
    
    % 循环 3 种压缩率 (线型)
    for r = 1:num_ratios
        data = all_data{r};
        snr_vec = data.SNR_vec;
        data_to_plot = data.(data_field);
        
        plot_label = sprintf('%s (%s)', method_name, ratio_labels{r});
        
        plot(snr_vec, data_to_plot, ...
             ratio_styles{r}, ...
            'Color', method_color, ...
            'LineWidth', line_width, ...
            'MarkerSize', marker_size, ...
            'DisplayName', plot_label);
    end
end

hold off;
grid on; 
xlabel('SNR (dB)'); 
ylabel('OSPA (m)'); 
ylim([0, 4]); 

lgd = legend('show', 'Location', 'best', 'NumColumns', 2); % 2列图例
lgd.FontSize = legend_font_size;
set(gca, 'FontSize', font_size);

print(gcf, 'sim4', '-depsc', '-painters');