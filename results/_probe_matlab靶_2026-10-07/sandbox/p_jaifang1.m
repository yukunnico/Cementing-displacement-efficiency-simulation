clear;clc;   %zzzzzzzzzz2026.6.10 HT1-004更新
%% 初始化数据
% 请确保读取的 .csv 是包含变径点和实测数据的最终版表格
structure_data = readtable('呼1-004井身结构.csv');
c_depth = structure_data.depth_well_logging_m_;
% 节点数
n_segment = height(structure_data);  % 井身分段数量（表格行数）

%% 核心：基于HT1-004固井管柱结构的七段式深度划分
% 定义七个深度区间的阈值（六个分界点）
% 井底钻达深度 (TD): 7660m
depth_threshold1 = 4025.73;   % 第一变径点：149.2mm钻杆 变 127mm钻杆
depth_threshold2 = 5243.21;   % 第二变径点：127mm钻杆 变 168.3mm尾管(悬挂器/尾管顶部位置)
depth_threshold3 = 5578.00;   % 第三变径点：273.1mm套管鞋 变 241.3mm裸眼 【关注点1：套管鞋】
depth_threshold4 = 7378.05;   % 第四变径点：168.3mm尾管 变 139.7mm尾管
depth_threshold5 = 7521.00;   % 第五变径点：241.3mm裸眼 变 215.9mm裸眼
depth_threshold6 = 7660.00;   % 第六变径点：井底TD
length_casing = 7660; % 下入管柱总长度,m

% 根据深度设置不同的直径值
% 预分配
diameter_bit_out = zeros(1, length(c_depth));
diameter_casing_out = zeros(1, length(c_depth));
diameter_casing_in = zeros(1, length(c_depth));

% 使用逻辑索引
idx1 = (c_depth <= depth_threshold1) ;
idx2 = (c_depth > depth_threshold1) & (c_depth <= depth_threshold2);
idx3 = (c_depth > depth_threshold2) & (c_depth <= depth_threshold3);
idx4 = (c_depth > depth_threshold3) & (c_depth <= depth_threshold4);
idx5 = (c_depth > depth_threshold4) & (c_depth <= depth_threshold5);
idx6 = (c_depth > depth_threshold5) & (c_depth <= depth_threshold6);
idx7 = (c_depth > depth_threshold6);

% 分段赋值（统一换算为国际单位: 米 m）
% 段1: 149.2mm 钻杆 在 273.1mm 套管内 (0-4025.73m)
diameter_bit_out(idx1) = 245.37 * 0.001;    % 1.4.2电测井径/环容表：245.37mm
diameter_casing_out(idx1) = 149.2 * 0.001;  % 149.2钻杆外径
diameter_casing_in(idx1) = (149.2 - 9.65 * 2) * 0.001; % 149.2钻杆内径(壁厚9.65mm)
% 段2: 127mm 钻杆 在 273.1mm 套管内 (4025.73-5243.21m)
diameter_bit_out(idx2) = 245.37 * 0.001;
diameter_casing_out(idx2) = 127.0 * 0.001;  % 127钻杆外径
diameter_casing_in(idx2) = (127.0 - 9.65 * 2) * 0.001; % 127钻杆内径(壁厚9.65mm)
% 段3: 168.3mm 尾管 在 273.1mm 套管内 (5243.21-5578m, 重叠段)
diameter_bit_out(idx3) = 245.37 * 0.001;
diameter_casing_out(idx3) = 168.3 * 0.001;  % 168.3尾管外径
diameter_casing_in(idx3) = (168.3 - 15.88 * 2) * 0.001; % 168.3尾管内径(壁厚15.88mm)
% 段4: 168.3mm 尾管 在实测裸眼中 (5578-7378.05m) 【关注点1】
diameter_bit_out(idx4) = structure_data.annulus_radius_array_cm_(idx4) * 0.01; 
diameter_casing_out(idx4) = 168.3 * 0.001;
diameter_casing_in(idx4) = (168.3 - 15.88 * 2) * 0.001; % 壁厚15.88mm
% 段5: 139.7mm 尾管 在实测裸眼中 (7378.05-7521m)
diameter_bit_out(idx5) = structure_data.annulus_radius_array_cm_(idx5) * 0.01;
diameter_casing_out(idx5) = 139.7 * 0.001;
diameter_casing_in(idx5) = (139.7 - 15.88 * 2) * 0.001; % 壁厚15.88mm
% 段6: 139.7mm 尾管 在实测裸眼中 (7521-7660m)
diameter_bit_out(idx6) = structure_data.annulus_radius_array_cm_(idx6) * 0.01;
diameter_casing_out(idx6) = 139.7 * 0.001;
diameter_casing_in(idx6) = (139.7 - 15.88 * 2) * 0.001;
% 段7: 7660m 以下 (尾管鞋以下纯裸眼口袋)
diameter_bit_out(idx7) = structure_data.annulus_radius_array_cm_(idx7) * 0.01;
diameter_casing_out(idx7) = 0;              % 尾管鞋以下，无内管
diameter_casing_in(idx7) = 0;               % 无内管
area_cout = pi * (diameter_bit_out .^ 2-diameter_casing_out .^ 2) / 4;  % 环空截面积（随深度分段变化，数组形式）m2
area_cin = pi * (diameter_casing_in .^ 2) / 4;  % 套管/钻杆内截面积（随深度分段变化，数组形式）m2
out_diam_bole = structure_data.annulus_radius_array_cm_ * 10;  % 井眼直径 (cm转为mm)
pianxin = true;

% 按井段分别计算套管内体积，之后累加求和
volume_in_casing_seg = area_cin .* structure_data.length_segment_array_m_';  % 该段体积 = 本段内径面积 * 本段套管/钻杆长度
volume_in_annual_seg = area_cout .* structure_data.length_segment_array_m_';  % 该段体积 = 本段环空面积 * 本段套管/钻杆长度
idx_liner_top = find(c_depth >= depth_threshold2, 1, 'first');
idx_slurry_top = find(c_depth >= 3381, 1, 'first');
volume_in_fenggu = sum(volume_in_annual_seg(idx_liner_top:end));  % 5243.21m至井底封固段环空体积
volume_in_drilling_casing = sum(volume_in_casing_seg);     % 套管内总体积 = 各段体积之和（m³）
volume_in_annual = sum(volume_in_annual_seg);     % 环空内总体积 = 各段体积之和（m³）
volume_in_drilling_casing_L = volume_in_drilling_casing * 1000;     % 转换为升（L）
volume_in_annual_L = volume_in_annual * 1000;     % 转换为升（L）

volume_in_allfluid = volume_in_drilling_casing + sum(volume_in_annual_seg(idx_slurry_top:end));  % 套管内容积+3381m至井底环空体积


%% ==================== 初始参数值调整模块 ====================
% 依据 HT1-004固井施工设计 7.2/8.11 施工过程模拟划分为变排量体系
% 0-钻井液(井浆) | 1-先导浆 | 2-低失水驱油隔离液1 | 3-低失水驱油隔离液2 | 4-领浆 | 5-尾浆 | 6-压塞液 | 7-替浆钻井液 | 8-保护液 | 9-基液 | 91-替浆钻井液1 | 92-替浆钻井液2 | 93-替浆钻井液3 | 94-替浆钻井液4 | 95-替浆钻井液5

 % 泵排量数据 (单位: L/min) - 来自8.11固井施工过程模拟
 pump_rate1 = 1.4 * 1000;            % 1-先导浆
 pump_rate2 = 1.2 * 1000;            % 2-低失水驱油隔离液1
 pump_rate3 = 1.2 * 1000;            % 3-低失水驱油隔离液2
 pump_rate4 = 1.2 * 1000;            % 4-领浆
 pump_rate5 = 1.25 * 1000;            % 5-尾浆 1.2-1.4
 pump_rate6 = 1.2 * 1000;            % 6-压塞液
 pump_rate7 = 1.5 * 1000;            % 7-替浆钻井液
 pump_rate8 = 1.4 * 1000;            % 8-保护液
 pump_rate9 = 1.4 * 1000;            % 9-基液
 pump_rate91 = 1.2 * 1000;           % 91-替浆钻井液1
 pump_rate92 = 1.0 * 1000;           % 92-替浆钻井液2
 pump_rate93 = 0.9 * 1000;           % 93-替浆钻井液3
 pump_rate94 = 0.8 * 1000;           % 94-替浆钻井液4
 pump_rate95 = 0.7 * 1000;           % 95-替浆钻井液5
 Pump_values = [pump_rate1, pump_rate2, pump_rate3, pump_rate4, pump_rate5, pump_rate6, pump_rate7, pump_rate8, pump_rate9, pump_rate91, pump_rate92, pump_rate93, pump_rate94, pump_rate95];
    
% 控压参数值
pressure_back = 0;                  % 控压钻进过程中的井口回压（MPa）
pressure_back_static = 0;           % 控压静止过程中的井口回压（MPa）

% 固井流体参数：黏度/屈服应力沿用原代码，密度按HT1-004施工工艺表更新

rou0 = 1.90; rou1 = 1.75; rou2 = 1.95; rou3 = 1.75; rou4 = 1.93; rou5 = 1.90; rou6 = 1.70; rou7 = 1.90; rou8 = 1.90; rou9 = 1.02; rou91 = 1.90; rou92 = 1.90; rou93 = 1.90; rou94 = 1.90; rou95 = 1.90; % 密度（g/cm³）
miu0 = 53 ; miu1 = 58; miu2 = 58 ; miu3 = 65; miu4 =  200; miu5 = 180; miu6 = 50; miu7 = 50; miu8 = 50; miu9 = 50; miu91 = 55; miu92 = 55; miu93 = 55; miu94 = 55; miu95 = 55; % 黏度（mPa·s）沿用旧参数
tau0 = 8.5 ; tau1 = 9.8; tau2 = 9.8; tau3 = 10; tau4 = 14; tau5 = 14; tau6 = 9; tau7 = 9.5; tau8 = 9.2; tau9 = 9; tau91 = 9.5; tau92 = 9.5; tau93 = 9.5; tau94 = 9.5; tau95 = 9.5;  % 屈服应力（Pa）沿用旧参数

%% 各流体体积参数 (来自7.1/7.2施工量与8.11施工过程模拟)
v1 = 25 * 1000;   % 1-先导浆体积,L
v2 = 16 * 1000;   % 2-低失水驱油隔离液1体积,L
v3 = 10 * 1000;   % 3-低失水驱油隔离液2体积,L
v4 = 48 * 1000;   % 4-领浆体积,L
v5 = 28 * 1000;   % 5-尾浆体积,L
v6 = 2 * 1000;  % 6-压塞液体积,L (2m³)
v7 = 29 * 1000; % 7-替浆钻井液体积,L (29m³)
v8 = 14 * 1000; % 8-保护液体积,L (14m³)
v9 = 1 * 1000;  % 9-基液体积,L (1m³)
v91 = 14 * 1000;   % 91-替浆钻井液1体积,L
v92 = 10 * 1000;  % 92-替浆钻井液2体积,L
v93 = 10 * 1000;  % 93-替浆钻井液3体积,L
v94 = 10 * 1000; % 94-替浆钻井液4体积,L
v95 = 7.1 * 1000;           % 95-替浆钻井液5体积,L

% 计算总时间步数与分段数量
dt = 1;  % 【高精度优化】：大幅细化时间步长为 min
n_time = floor((v1./Pump_values(1) + v2./Pump_values(2) + v3./Pump_values(3) + v4./Pump_values(4) + v5./Pump_values(5) + v6./Pump_values(6) + v7./Pump_values(7) + v8./Pump_values(8) + v9./Pump_values(9) + v91./Pump_values(10) + v92./Pump_values(11) + v93./Pump_values(12) + v94./Pump_values(13) + v95./Pump_values(14)) / dt);
time = (v1./Pump_values(1) + v2./Pump_values(2) + v3./Pump_values(3) + v4./Pump_values(4) + v5./Pump_values(5) + v6./Pump_values(6) + v7./Pump_values(7) + v8./Pump_values(8) + v9./Pump_values(9) + v91./Pump_values(10) + v92./Pump_values(11) + v93./Pump_values(12) + v94./Pump_values(13) + v95./Pump_values(14)) / dt;


%% 回压施加模块

backpressure_1 = 0;           % 1-先导浆回压，MPa
backpressure_2 = 0;           % 2-低失水驱油隔离液1回压，MPa
backpressure_3 = 0;           % 3-低失水驱油隔离液2回压，MPa
backpressure_4 = 0;           % 4-领浆回压，MPa
backpressure_5 = 0;           % 5-尾浆回压，MPa
backpressure_6 = 0;           % 6-压塞液回压，MPa
backpressure_7 = 0;           % 7-替浆钻井液回压，MPa
backpressure_8 = 0;           % 8-保护液回压，MPa
backpressure_9 = 0;           % 9-基液回压，MPa
backpressure_91 = 0;          % 91-替浆钻井液1回压，MPa
backpressure_92 = 0;          % 92-替浆钻井液2回压，MPa
backpressure_93 = 0;          % 93-替浆钻井液3回压，MPa
backpressure_94 = 0;          % 94-替浆钻井液4回压，MPa
backpressure_95 = 0;          % 95-替浆钻井液5回压，MPa
backpressure = [backpressure_1, backpressure_2, backpressure_3, backpressure_4, backpressure_5, backpressure_6, backpressure_7, backpressure_8, backpressure_9, backpressure_91, backpressure_92, backpressure_93, backpressure_94, backpressure_95];

%% 数组初始化与高精度地热场构建
vertical_length_all_grid = zeros(1, n_segment);
cos_deg = zeros(1, n_segment);
for i = 1:n_segment
    cos_deg(i) = cos(deg2rad(structure_data.deg_for_logging_degree_(i)));
    vertical_length_all_grid(i) = abs(structure_data.length_segment_array_m_(i) * cos_deg(i));
end
TVD_cum = cumsum(vertical_length_all_grid); % 全井累计垂深

% surface_temp = ******; % 地面温度 ℃
% bottom_temp = ******;  % 井底温度 ℃
% Temp_grid = surface_temp + (bottom_temp - surface_temp) * (TVD_cum / max(TVD_cum)); % 地温梯度线性场

% 环空内外径数组初始化
d_i = zeros(1, n_segment);
d_o = zeros(1, n_segment);
for i = 1:n_segment 
    d_i(i) = diameter_casing_out(i); 
    d_o(i) = structure_data.annulus_radius_array_cm_(i) / 100;
end

% 环空注入流体体积数组初始化 (仅前5段流体进入环空)
volume_injected_1_list = zeros(1, n_time);  
volume_injected_2_list = zeros(1, n_time);  
volume_injected_3_list = zeros(1, n_time);  
volume_injected_4_list = zeros(1, n_time);  
volume_injected_5_list = zeros(1, n_time);  

% 套管注入流体体积数组初始化（从井口向下计算，扩容至替浆91/92/93/94）
volume_injected_casing_1_list = zeros(1, n_time);  
volume_injected_casing_2_list = zeros(1, n_time);  
volume_injected_casing_3_list = zeros(1, n_time);  
volume_injected_casing_4_list = zeros(1, n_time);  
volume_injected_casing_5_list = zeros(1, n_time);  
volume_injected_casing_6_list = zeros(1, n_time);  
volume_injected_casing_7_list = zeros(1, n_time);  
volume_injected_casing_8_list = zeros(1, n_time);  
volume_injected_casing_9_list = zeros(1, n_time);  
volume_injected_casing_91_list = zeros(1, n_time);  
volume_injected_casing_92_list = zeros(1, n_time);  
volume_injected_casing_93_list = zeros(1, n_time);  
volume_injected_casing_94_list = zeros(1, n_time);  
volume_injected_casing_95_list = zeros(1, n_time);  

% 环空注入单一流体累计体积数组初始化
volume_into_1_annulus = zeros(1, n_time);  
volume_into_2_annulus = zeros(1, n_time);  
volume_into_3_annulus = zeros(1, n_time);  
volume_into_4_annulus = zeros(1, n_time);  
volume_into_5_annulus = zeros(1, n_time);  

% 环空内界面参数初始化 (追踪 0_1 到 4_5 共 5 个界面)
depth_tag_liquid_0_1_list = zeros(1, n_time); residual_volume_0_1_list = zeros(1, n_time); residual_height_0_1_list = zeros(1, n_time); vertical_residual_height_0_1_list = zeros(1, n_time);
depth_tag_liquid_1_2_list = zeros(1, n_time); residual_volume_1_2_list = zeros(1, n_time); residual_height_1_2_list = zeros(1, n_time); vertical_residual_height_1_2_list = zeros(1, n_time);
depth_tag_liquid_2_3_list = zeros(1, n_time); residual_volume_2_3_list = zeros(1, n_time); residual_height_2_3_list = zeros(1, n_time); vertical_residual_height_2_3_list = zeros(1, n_time);
depth_tag_liquid_3_4_list = zeros(1, n_time); residual_volume_3_4_list = zeros(1, n_time); residual_height_3_4_list = zeros(1, n_time); vertical_residual_height_3_4_list = zeros(1, n_time);
depth_tag_liquid_4_5_list = zeros(1, n_time); residual_volume_4_5_list = zeros(1, n_time); residual_height_4_5_list = zeros(1, n_time); vertical_residual_height_4_5_list = zeros(1, n_time);
tag_interface = zeros(n_time, 5);

% 套管内界面参数初始化 (追踪扩充涵盖 9/91/95/92/93/94 界面)
depth_tag_casing_0_1_list = zeros(1, n_time); residual_volume_casing_0_1_list = zeros(1, n_time); residual_height_casing_0_1_list = zeros(1, n_time);
depth_tag_casing_1_2_list = zeros(1, n_time); residual_volume_casing_1_2_list = zeros(1, n_time); residual_height_casing_1_2_list = zeros(1, n_time);
depth_tag_casing_2_3_list = zeros(1, n_time); residual_volume_casing_2_3_list = zeros(1, n_time); residual_height_casing_2_3_list = zeros(1, n_time);
depth_tag_casing_3_4_list = zeros(1, n_time); residual_volume_casing_3_4_list = zeros(1, n_time); residual_height_casing_3_4_list = zeros(1, n_time);
depth_tag_casing_4_5_list = zeros(1, n_time); residual_volume_casing_4_5_list = zeros(1, n_time); residual_height_casing_4_5_list = zeros(1, n_time);
depth_tag_casing_5_6_list = zeros(1, n_time); residual_volume_casing_5_6_list = zeros(1, n_time); residual_height_casing_5_6_list = zeros(1, n_time);
depth_tag_casing_6_7_list = zeros(1, n_time); residual_volume_casing_6_7_list = zeros(1, n_time); residual_height_casing_6_7_list = zeros(1, n_time);
depth_tag_casing_7_8_list = zeros(1, n_time); residual_volume_casing_7_8_list = zeros(1, n_time); residual_height_casing_7_8_list = zeros(1, n_time);
depth_tag_casing_8_9_list = zeros(1, n_time); residual_volume_casing_8_9_list = zeros(1, n_time); residual_height_casing_8_9_list = zeros(1, n_time);
depth_tag_casing_9_91_list = zeros(1, n_time); residual_volume_casing_9_91_list = zeros(1, n_time); residual_height_casing_9_91_list = zeros(1, n_time);
depth_tag_casing_91_92_list = zeros(1, n_time); residual_volume_casing_91_92_list = zeros(1, n_time); residual_height_casing_91_92_list = zeros(1, n_time);
depth_tag_casing_92_93_list = zeros(1, n_time); residual_volume_casing_92_93_list = zeros(1, n_time); residual_height_casing_92_93_list = zeros(1, n_time);
depth_tag_casing_93_94_list = zeros(1, n_time); residual_volume_casing_93_94_list = zeros(1, n_time); residual_height_casing_93_94_list = zeros(1, n_time);
depth_tag_casing_94_95_list = zeros(1, n_time); residual_volume_casing_94_95_list = zeros(1, n_time); residual_height_casing_94_95_list = zeros(1, n_time);
tag_interface_casing = zeros(n_time, 14);

% 密度、黏度、屈服应力、流速数组初始化
rou_annulus_all_time = zeros(n_time, n_segment);    
miu_annulus_all_time = zeros(n_time, n_segment);    
tau_annulus_all_time = zeros(n_time, n_segment);    
velo_annulus_all_time = zeros(n_time, n_segment);   

% 摩擦阻力与流型数组初始化
Ff_a = zeros(n_time, n_segment);          
flow_pattern_a = zeros(n_time, n_segment); 

% 压力存储数组初始化
pressure_annuli = zeros(n_time, n_segment);          
pressure_annuli_static = zeros(n_time, n_segment);   
pressure_annuli_friction = zeros(n_time, n_segment); 
pressure_casing = zeros(n_time, n_segment);            
pressure_casing_static = zeros(n_time, n_segment);     
pressure_casing_static_MPa = zeros(n_time, n_segment); 
pressure_casing_friction = zeros(n_time, n_segment);   
pressure_casing_friction_MPa = zeros(n_time, n_segment);   
pressure_casing_MPa = zeros(n_time, n_segment);        
pressure_pump_surface = zeros(1, n_time);              
P_bit_drop = 0;                                        

% 泵排量随时间变化数组初始化
Pump_values_time_list = zeros(1, n_time);         
volume_injected_all_list = zeros(1, n_time);   

% 回压随时间变化数组初始化
backpressure_time_list = zeros(1, n_time);

% 从底部计算各分段累加体积和累加长度
volume_all_from_bottom = zeros(1, n_segment);
length_all_from_bottom = zeros(1, n_segment);
for i = 1:n_segment
    volume_all_from_bottom(i) = sum(structure_data.volume_annulus_L_(i:end));  
    length_all_from_bottom(i) = sum(structure_data.length_segment_array_m_(i:end));  
end

% 初始时刻（t=1）参数赋值
Pump_values_time_list(1) = Pump_values(1);
backpressure_time_list(1) = backpressure(1);
volume_injected_all_list(1) = Pump_values_time_list(1) * dt;  

%% 第一步循环：计算各时刻泵排量与累计注入体积
pump_time_node = zeros(14,1);
pump_time_node(1) = (v1/Pump_values(1)) / dt;
pump_time_node(2) = (v1/Pump_values(1)+v2/Pump_values(2)) / dt;
pump_time_node(3) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)) / dt;
pump_time_node(4) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)) / dt;
pump_time_node(5) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)) / dt;
pump_time_node(6) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)) / dt;
pump_time_node(7) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)) / dt;
pump_time_node(8) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)) / dt;
pump_time_node(9) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)+v9/Pump_values(9)) / dt;
pump_time_node(10) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)+v9/Pump_values(9)+v91/Pump_values(10)) / dt;
pump_time_node(11) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)+v9/Pump_values(9)+v91/Pump_values(10)+v92/Pump_values(11)) / dt;
pump_time_node(12) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)+v9/Pump_values(9)+v91/Pump_values(10)+v92/Pump_values(11)+v93/Pump_values(12)) / dt;
pump_time_node(13) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)+v9/Pump_values(9)+v91/Pump_values(10)+v92/Pump_values(11)+v93/Pump_values(12)+v94/Pump_values(13)) / dt;
pump_time_node(14) = (v1/Pump_values(1)+v2/Pump_values(2)+v3/Pump_values(3)+v4/Pump_values(4)+v5/Pump_values(5)+v6/Pump_values(6)+v7/Pump_values(7)+v8/Pump_values(8)+v9/Pump_values(9)+v91/Pump_values(10)+v92/Pump_values(11)+v93/Pump_values(12)+v94/Pump_values(13)+v95/Pump_values(14)) / dt;

for t = 2:n_time
    if t < pump_time_node(1)
        Pump_values_time_list(t) = Pump_values(1);
       backpressure_time_list(t) = backpressure(1);
    elseif t >= pump_time_node(1) && t <= pump_time_node(2)
        Pump_values_time_list(t) = Pump_values(2);
       backpressure_time_list(t) = backpressure(2);
    elseif t >= pump_time_node(2) && t <= pump_time_node(3)
        Pump_values_time_list(t) = Pump_values(3);
       backpressure_time_list(t) = backpressure(3);
    elseif t > pump_time_node(3) && t <= pump_time_node(4)
        Pump_values_time_list(t) = Pump_values(4);
       backpressure_time_list(t) = backpressure(4);
    elseif  t > pump_time_node(4) && t <= pump_time_node(5)
        Pump_values_time_list(t) = Pump_values(5);
       backpressure_time_list(t) = backpressure(5);
    elseif  t > pump_time_node(5) && t <= pump_time_node(6)
        Pump_values_time_list(t) = Pump_values(6);
       backpressure_time_list(t) = backpressure(6);
    elseif  t > pump_time_node(6) && t <= pump_time_node(7)
        Pump_values_time_list(t) = Pump_values(7);
       backpressure_time_list(t) = backpressure(7);
    elseif  t > pump_time_node(7) && t <= pump_time_node(8)
        Pump_values_time_list(t) = Pump_values(8);
       backpressure_time_list(t) = backpressure(8);
    elseif  t > pump_time_node(8) && t <= pump_time_node(9)
        Pump_values_time_list(t) = Pump_values(9);
       backpressure_time_list(t) = backpressure(9);
    elseif  t > pump_time_node(9) && t <= pump_time_node(10)
        Pump_values_time_list(t) = Pump_values(10);
       backpressure_time_list(t) = backpressure(10);
    elseif  t > pump_time_node(10) && t <= pump_time_node(11)
        Pump_values_time_list(t) = Pump_values(11);
       backpressure_time_list(t) = backpressure(11);
    elseif  t > pump_time_node(11) && t <= pump_time_node(12)
        Pump_values_time_list(t) = Pump_values(12);
       backpressure_time_list(t) = backpressure(12);
    elseif  t > pump_time_node(12) && t <= pump_time_node(13)
        Pump_values_time_list(t) = Pump_values(13);
       backpressure_time_list(t) = backpressure(13);
    elseif  t > pump_time_node(13) && t <= pump_time_node(14)
        Pump_values_time_list(t) = Pump_values(14);
       backpressure_time_list(t) = backpressure(14);
    elseif t > pump_time_node(14)
        Pump_values_time_list(t) = 0;
       backpressure_time_list(t) = 0;
    end
    % 【高精度优化】：累计体积累加需要乘以 dt
    volume_injected_all_list(t) =  volume_injected_all_list(t-1) + Pump_values_time_list(t) * dt; 
end
% 泵排量单位转换：L/min → m³/s
Pump_values_time_list_m3_s = Pump_values_time_list ./ 60 ./ 1000;

%% 第二步循环：计算各时刻各流体注入环空的体积
for t = 1:n_time
    if volume_injected_all_list(t) > volume_in_drilling_casing_L           
       volume_injected_1_list(t) = volume_injected_1_list(t-1) + Pump_values_time_list(t) * dt;
    end
    if volume_injected_all_list(t) > volume_in_drilling_casing_L + v1         
       volume_injected_2_list(t) = volume_injected_2_list(t-1) + Pump_values_time_list(t) * dt;
    end
    if volume_injected_all_list(t) > volume_in_drilling_casing_L + v1 + v2
       volume_injected_3_list(t) = volume_injected_3_list(t-1) + Pump_values_time_list(t) * dt;
    end
    if volume_injected_all_list(t) > volume_in_drilling_casing_L + v1 + v2 + v3
       volume_injected_4_list(t) = volume_injected_4_list(t-1) + Pump_values_time_list(t) * dt;
    end
    if volume_injected_all_list(t) > volume_in_drilling_casing_L + v1 + v2 + v3 + v4
       volume_injected_5_list(t) = volume_injected_5_list(t-1) + Pump_values_time_list(t) * dt;
    end
end

%% 第三步主循环：界面追踪、流体属性分配、流速与压力计算
for t = 1:n_time
    for i = 2:n_segment
        % 判断界面深度，钻井液-先导浆
        if volume_injected_1_list(t) >= volume_all_from_bottom(i) && volume_injected_1_list(t) < volume_all_from_bottom(i-1)
            depth_tag_liquid_0_1_list(t) = i - 1; 
            residual_volume_0_1 = volume_injected_1_list(t) - volume_all_from_bottom(i);
            residual_volume_0_1_list(t) = residual_volume_0_1;
            residual_height_0_1 = residual_volume_0_1 / structure_data.square_annulus_dm2_(i-1) / 10; 
            residual_height_0_1_list(t) = residual_height_0_1;
            vertical_residual_height_0_1_list(t) = residual_height_0_1 * cosd(structure_data.deg_for_logging_degree_(i-1));
        elseif volume_injected_1_list(t) > 0 && volume_injected_1_list(t) < volume_all_from_bottom(end) 
            depth_tag_liquid_0_1_list(t) = n_segment;
            residual_volume_0_1 = volume_injected_1_list(t);
            residual_volume_0_1_list(t) = residual_volume_0_1;
            residual_height_0_1 = residual_volume_0_1 / structure_data.square_annulus_dm2_(n_segment) / 10; 
            residual_height_0_1_list(t) = residual_height_0_1;
            vertical_residual_height_0_1_list(t) = residual_height_0_1 * cosd(structure_data.deg_for_logging_degree_(n_segment));
        elseif volume_injected_1_list(t) > volume_all_from_bottom(1)
            depth_tag_liquid_0_1_list(t) = -1;
            residual_volume_0_1_list(t) = -1;
            residual_height_0_1_list(t) = -1;
        elseif volume_injected_1_list(t) <=0    
            depth_tag_liquid_0_1_list(t) = 10000;
            residual_volume_0_1_list(t) = 10000;
            residual_height_0_1_list(t) = 10000;
        end
        tag_interface(t,1) = depth_tag_liquid_0_1_list(t);
        
        % 第二界面,先导浆-低失水驱油隔离液1
        if volume_injected_2_list(t) >= volume_all_from_bottom(i) && volume_injected_2_list(t) < volume_all_from_bottom(i-1)
            depth_tag_liquid_1_2_list(t) = i - 1;  
            residual_volume_1_2 = volume_injected_2_list(t) - volume_all_from_bottom(i);
            residual_volume_1_2_list(t) = residual_volume_1_2;
            residual_height_1_2 = residual_volume_1_2 / structure_data.square_annulus_dm2_(i-1) / 10; 
            residual_height_1_2_list(t) = residual_height_1_2;
            vertical_residual_height_1_2_list(t) = residual_height_1_2 * cosd(structure_data.deg_for_logging_degree_(i-1));
        elseif volume_injected_2_list(t) > 0 && volume_injected_2_list(t) < volume_all_from_bottom(end) 
            depth_tag_liquid_1_2_list(t) = n_segment;
            residual_volume_1_2 = volume_injected_2_list(t);
            residual_volume_1_2_list(t) = residual_volume_1_2;
            residual_height_1_2 = residual_volume_1_2 / structure_data.square_annulus_dm2_(n_segment) / 10; 
            residual_height_1_2_list(t) = residual_height_1_2;
            vertical_residual_height_1_2_list(t) = residual_height_1_2 * cosd(structure_data.deg_for_logging_degree_(n_segment));
        elseif volume_injected_2_list(t) > volume_all_from_bottom(1)
            depth_tag_liquid_1_2_list(t) = -1;
            residual_volume_1_2_list(t) = -1;
            residual_height_1_2_list(t) = -1;
        elseif volume_injected_2_list(t) <=0
            depth_tag_liquid_1_2_list(t) = 10000;
            residual_volume_1_2_list(t) = 10000;
            residual_height_1_2_list(t) = 10000;
        end
        tag_interface(t,2) = depth_tag_liquid_1_2_list(t);
        
        % 第三界面，低失水驱油隔离液1-低失水驱油隔离液2
        if volume_injected_3_list(t) >= volume_all_from_bottom(i) && volume_injected_3_list(t) < volume_all_from_bottom(i-1)
            depth_tag_liquid_2_3_list(t) = i - 1;  
            residual_volume_2_3 = volume_injected_3_list(t) - volume_all_from_bottom(i);
            residual_volume_2_3_list(t) = residual_volume_2_3;
            residual_height_2_3 = residual_volume_2_3 / structure_data.square_annulus_dm2_(i-1) / 10; 
            residual_height_2_3_list(t) = residual_height_2_3;
            vertical_residual_height_2_3_list(t) = residual_height_2_3 * cosd(structure_data.deg_for_logging_degree_(i-1));
        elseif volume_injected_3_list(t) > 0 && volume_injected_3_list(t) < volume_all_from_bottom(end) 
            depth_tag_liquid_2_3_list(t) = n_segment;
            residual_volume_2_3 = volume_injected_3_list(t);
            residual_volume_2_3_list(t) = residual_volume_2_3;
            residual_height_2_3 = residual_volume_2_3 / structure_data.square_annulus_dm2_(n_segment) / 10; 
            residual_height_2_3_list(t) = residual_height_2_3;
            vertical_residual_height_2_3_list(t) = residual_height_2_3 * cosd(structure_data.deg_for_logging_degree_(n_segment));
        elseif volume_injected_3_list(t) > volume_all_from_bottom(1)
            depth_tag_liquid_2_3_list(t) = -1;
            residual_volume_2_3_list(t) = -1;
            residual_height_2_3_list(t) = -1;
        elseif volume_injected_3_list(t) <= 0 
            depth_tag_liquid_2_3_list(t) = 10000;
            residual_volume_2_3_list(t) = 10000;
            residual_height_2_3_list(t) = 10000;
        end
        tag_interface(t,3) = depth_tag_liquid_2_3_list(t);
        
        % 第四界面，低失水驱油隔离液2-领浆
        if volume_injected_4_list(t) >= volume_all_from_bottom(i) && volume_injected_4_list(t) < volume_all_from_bottom(i-1)
            depth_tag_liquid_3_4_list(t) = i - 1;  
            residual_volume_3_4 = volume_injected_4_list(t) - volume_all_from_bottom(i);
            residual_volume_3_4_list(t) = residual_volume_3_4;
            residual_height_3_4 = residual_volume_3_4 / structure_data.square_annulus_dm2_(i-1) / 10; 
            residual_height_3_4_list(t) = residual_height_3_4;
            vertical_residual_height_3_4_list(t) = residual_height_3_4 * cosd(structure_data.deg_for_logging_degree_(i-1));
        elseif volume_injected_4_list(t) > 0 && volume_injected_4_list(t) < volume_all_from_bottom(end) 
            depth_tag_liquid_3_4_list(t) = n_segment;
            residual_volume_3_4 = volume_injected_4_list(t);
            residual_volume_3_4_list(t) = residual_volume_3_4;
            residual_height_3_4 = residual_volume_3_4 / structure_data.square_annulus_dm2_(n_segment) / 10; 
            residual_height_3_4_list(t) = residual_height_3_4;
            vertical_residual_height_3_4_list(t) = residual_height_3_4 * cosd(structure_data.deg_for_logging_degree_(n_segment));
        elseif volume_injected_4_list(t) > volume_all_from_bottom(1)
            depth_tag_liquid_3_4_list(t) = -1;
            residual_volume_3_4_list(t) = -1;
            residual_height_3_4_list(t) = -1;
        elseif volume_injected_4_list(t) <= 0 
            depth_tag_liquid_3_4_list(t) = 10000;
            residual_volume_3_4_list(t) = 10000;
            residual_height_3_4_list(t) = 10000;
        end
        tag_interface(t,4) = depth_tag_liquid_3_4_list(t);

        % 第五界面，领浆-尾浆
        if volume_injected_5_list(t) >= volume_all_from_bottom(i) && volume_injected_5_list(t) < volume_all_from_bottom(i-1)
            depth_tag_liquid_4_5_list(t) = i - 1;  
            residual_volume_4_5 = volume_injected_5_list(t) - volume_all_from_bottom(i);
            residual_volume_4_5_list(t) = residual_volume_4_5;
            residual_height_4_5 = residual_volume_4_5 / structure_data.square_annulus_dm2_(i-1) / 10; 
            residual_height_4_5_list(t) = residual_height_4_5;
            vertical_residual_height_4_5_list(t) = residual_height_4_5 * cosd(structure_data.deg_for_logging_degree_(i-1));
        elseif volume_injected_5_list(t) > 0 && volume_injected_5_list(t) < volume_all_from_bottom(end) 
            depth_tag_liquid_4_5_list(t) = n_segment;
            residual_volume_4_5 = volume_injected_5_list(t);
            residual_volume_4_5_list(t) = residual_volume_4_5;
            residual_height_4_5 = residual_volume_4_5 / structure_data.square_annulus_dm2_(n_segment) / 10; 
            residual_height_4_5_list(t) = residual_height_4_5;
            vertical_residual_height_4_5_list(t) = residual_height_4_5 * cosd(structure_data.deg_for_logging_degree_(n_segment));
        elseif volume_injected_5_list(t) > volume_all_from_bottom(1)
            depth_tag_liquid_4_5_list(t) = -1;
            residual_volume_4_5_list(t) = -1;
            residual_height_4_5_list(t) = -1;
        elseif volume_injected_5_list(t) <= 0 
            depth_tag_liquid_4_5_list(t) = 10000;
            residual_volume_4_5_list(t) = 10000;
            residual_height_4_5_list(t) = 10000;
        end
        tag_interface(t,5) = depth_tag_liquid_4_5_list(t);
    end   
    
% --------- 2. 环空流体属性分配（密度、黏度、屈服应力） --------------------------
    for i = 1:n_segment
        % 保留原有的体积加权(VOF)高精度混合法则
        if i < depth_tag_liquid_0_1_list(t)
            rou_annulus_all_time(t, i) = rou0;
            miu_annulus_all_time(t, i) = miu0;
            tau_annulus_all_time(t, i) = tau0;
        elseif i == tag_interface(t,1)
            proportion_0_1_1 = residual_height_0_1_list(t) / structure_data.length_segment_array_m_(i);
            proportion_0_1_1 = min(max(proportion_0_1_1, 0), 1); % 【高精度安全限制】：防止微小舍入误差导致溢出
            proportion_0_1_0 = 1 - proportion_0_1_1;
            rou_annulus_all_time(t, i) = proportion_0_1_1 * rou1 + proportion_0_1_0 * rou0;
            miu_annulus_all_time(t, i) = proportion_0_1_1 * miu1 + proportion_0_1_0 * miu0;
            tau_annulus_all_time(t, i) = proportion_0_1_1 * tau1 + proportion_0_1_0 * tau0;
        elseif i>tag_interface(t,1) && i<tag_interface(t,2)
            rou_annulus_all_time(t, i) = rou1;
            miu_annulus_all_time(t, i) = miu1;
            tau_annulus_all_time(t, i) = tau1;
        elseif i == tag_interface(t,2)
            proportion_1_2_1 = residual_height_1_2_list(t) / structure_data.length_segment_array_m_(i);
            proportion_1_2_1 = min(max(proportion_1_2_1, 0), 1); 
            proportion_1_2_0 = 1 - proportion_1_2_1;
            rou_annulus_all_time(t, i) = proportion_1_2_1 * rou2 + proportion_1_2_0 * rou1;
            miu_annulus_all_time(t, i) = proportion_1_2_1 * miu2 + proportion_1_2_0 * miu1;
            tau_annulus_all_time(t, i) = proportion_1_2_1 * tau2 + proportion_1_2_0 * tau1;
        elseif i>tag_interface(t,2) && i<tag_interface(t,3)
            rou_annulus_all_time(t, i) = rou2;
            miu_annulus_all_time(t, i) = miu2;
            tau_annulus_all_time(t, i) = tau2;
        elseif i == tag_interface(t,3)
            proportion_2_3_1 = residual_height_2_3_list(t) / structure_data.length_segment_array_m_(i);
            proportion_2_3_1 = min(max(proportion_2_3_1, 0), 1); 
            proportion_2_3_0 = 1 - proportion_2_3_1;
            rou_annulus_all_time(t, i) = proportion_2_3_1 * rou3 + proportion_2_3_0 * rou2;
            miu_annulus_all_time(t, i) = proportion_2_3_1 * miu3 + proportion_2_3_0 * miu2;
            tau_annulus_all_time(t, i) = proportion_2_3_1 * tau3 + proportion_2_3_0 * tau2;
        elseif i>tag_interface(t,3) && i<tag_interface(t,4)
            rou_annulus_all_time(t, i) = rou3;
            miu_annulus_all_time(t, i) = miu3;
            tau_annulus_all_time(t, i) = tau3;
        elseif i == tag_interface(t,4)
            proportion_3_4_1 = residual_height_3_4_list(t) / structure_data.length_segment_array_m_(i);
            proportion_3_4_1 = min(max(proportion_3_4_1, 0), 1); 
            proportion_3_4_0 = 1 - proportion_3_4_1;
            rou_annulus_all_time(t, i) = proportion_3_4_1 * rou4 + proportion_3_4_0 * rou3;
            miu_annulus_all_time(t, i) = proportion_3_4_1 * miu4 + proportion_3_4_0 * miu3;
            tau_annulus_all_time(t, i) = proportion_3_4_1 * tau4 + proportion_3_4_0 * tau3;
        elseif i>tag_interface(t,4) && i<tag_interface(t,5)
            rou_annulus_all_time(t, i) = rou4;
            miu_annulus_all_time(t, i) = miu4;
            tau_annulus_all_time(t, i) = tau4;
        elseif i == tag_interface(t,5)
            proportion_4_5_1 = residual_height_4_5_list(t) / structure_data.length_segment_array_m_(i);
            proportion_4_5_1 = min(max(proportion_4_5_1, 0), 1); 
            proportion_4_5_0 = 1 - proportion_4_5_1;
            rou_annulus_all_time(t, i) = proportion_4_5_1 * rou5 + proportion_4_5_0 * rou4;
            miu_annulus_all_time(t, i) = proportion_4_5_1 * miu5 + proportion_4_5_0 * miu4;
            tau_annulus_all_time(t, i) = proportion_4_5_1 * tau5 + proportion_4_5_0 * tau4;
        elseif i>tag_interface(t,5)
            rou_annulus_all_time(t, i) = rou5;
            miu_annulus_all_time(t, i) = miu5;
            tau_annulus_all_time(t, i) = tau5;
        end
%         
%% --- 【高精度优化核心】：HPHT 温压耦合导致流体密度动态变化 ---
%         % 注意: 原赋值的是常温常压下的基准密度(g/cm3)，超深井必须根据当前深度温度和压力实时修正
%         alpha_t_ann = 3.5e-4; % 综合热膨胀系数(1/℃) (需要化验室提供准确值)
%         beta_p_ann  = 4.0e-4; % 综合压缩系数(1/MPa) 
%         
%         % 提取上一时间步的压力用于显式预测(防死循环)
%         if t == 1
%             P_pred_ann = (9.81 * rou_annulus_all_time(t, i) * 1000 * TVD_cum(i)) / 1e6;
%         else
%             P_pred_ann = pressure_annuli(t-1, i) / 1e6;
%         end
%         
%         % 公式： ρ_real = ρ_std * [1 - α(T - T0) + β(P - P0)]
%         rou_annulus_all_time(t, i) = rou_annulus_all_time(t, i) * (1 - alpha_t_ann * (Temp_grid(i) - 25) + beta_p_ann * (P_pred_ann - 0.1));
%         % -----------------------------------------------------------
        % 3. 环空流速计算
        velo_annulus_all_time(t,i) = Pump_values_time_list(t) / structure_data.square_annulus_dm2_(i) / 10 / 60;
    end
end
%% 第四步循环：环空压力计算（静液柱+循环压耗）
% =========================================================================
% 第一段：深度 <= depth_threshold1 (3321.682m)
segment1 = (c_depth <= depth_threshold1);      
if any(segment1)
    mean_diameter_segment1 = mean(out_diam_bole(segment1));
    mean_casing_segment1 = mean(diameter_casing_out(segment1))*1000;
    PR_segment1 = PR_Friction(mean_diameter_segment1, mean_casing_segment1);
else
    PR_segment1 = 0;
end
% 第二段：depth_threshold1-depth_threshold2 (3321.682-5307.539m)
segment2 = (c_depth > depth_threshold1) & (c_depth <= depth_threshold2);      
if any(segment2)
    mean_diameter_segment2 = mean(out_diam_bole(segment2));
    mean_casing_segment2 = mean(diameter_casing_out(segment2))*1000;
    PR_segment2 = PR_Friction(mean_diameter_segment2, mean_casing_segment2);
else
    PR_segment2 = 0;
end
% 第三段：depth_threshold2-depth_threshold3 (5243.21-5578m)
segment3 = (c_depth > depth_threshold2) & (c_depth <= depth_threshold3); 
if any(segment3)
    mean_diameter_segment3 = mean(out_diam_bole(segment3));
    mean_casing_segment3 = mean(diameter_casing_out(segment3))*1000; 
    PR_segment3 = PR_Friction(mean_diameter_segment3, mean_casing_segment3);
else
    PR_segment3 = 0;
end
% 第四段：depth_threshold3-depth_threshold4 (5578-7378.05m)
segment4 = (c_depth > depth_threshold3) & (c_depth <= depth_threshold4); 
if any(segment4)
    mean_diameter_segment4 = mean(out_diam_bole(segment4));
    mean_casing_segment4 = mean(diameter_casing_out(segment4))*1000;
    PR_segment4 = PR_Friction(mean_diameter_segment4, mean_casing_segment4);
else
    PR_segment4 = 0;
end
% 第五段：depth_threshold4-depth_threshold5 (7089.576-7096m)
segment5 = (c_depth > depth_threshold4) & (c_depth <= depth_threshold5);                          
if any(segment5)
    mean_diameter_segment5 = mean(out_diam_bole(segment5));
    mean_casing_segment5 = mean(diameter_casing_out(segment5))*1000;
    PR_segment5 = PR_Friction(mean_diameter_segment5, mean_casing_segment5);
else
    PR_segment5 = 0;
end
% 第六段：depth_threshold5-depth_threshold6 (7521-7660m)
segment6 = (c_depth > depth_threshold5) & (c_depth <= depth_threshold6);                          
if any(segment6)
    mean_diameter_segment6 = mean(out_diam_bole(segment6));
    mean_casing_segment6 = mean(diameter_casing_out(segment6))*1000;
    PR_segment6 = PR_Friction(mean_diameter_segment6, mean_casing_segment6);
else
    PR_segment6 = 0;
end
% 第七段：深度 > depth_threshold6 (7660m以下口袋)
segment7 = (c_depth > depth_threshold6);                          
if any(segment7)
    mean_diameter_segment7 = mean(out_diam_bole(segment7));
    mean_casing_segment7 = mean(diameter_casing_out(segment7))*1000;
    PR_segment7 = PR_Friction(mean_diameter_segment7, mean_casing_segment7);
else
    PR_segment7 = 0;
end
% 创建分段PR数组（基于7个深度的精确判断）
PR = zeros(1, n_segment);
for i = 1:n_segment
    if (c_depth(i) <= depth_threshold1)
        PR(i) = PR_segment1;
    elseif (c_depth(i) <= depth_threshold2)
        PR(i) = PR_segment2;
    elseif (c_depth(i) <= depth_threshold3)
        PR(i) = PR_segment3;
    elseif (c_depth(i) <= depth_threshold4)
        PR(i) = PR_segment4;
    elseif (c_depth(i) <= depth_threshold5)
        PR(i) = PR_segment5;
    elseif (c_depth(i) <= depth_threshold6)
        PR(i) = PR_segment6;
    else
        PR(i) = PR_segment7;
    end
end
% 环空截面积（dm²→m²）
A_annulus = structure_data.square_annulus_dm2_ / 100;  
rou_annulus_all_time_kg_m3 = rou_annulus_all_time * 1000;  % 密度转换（kg/m³）
for t = 1:n_time
    for i = 1:n_segment
        % 计算环空摩擦阻力
        [Ff_a(t,i), flow_pattern_a(t,i)] = Friction_annulus_bh(...
            rou_annulus_all_time_kg_m3(t,i), ...  
            velo_annulus_all_time(t,i), ...       
            A_annulus(i), ...                   
            miu_annulus_all_time(t,i), ...        
            tau_annulus_all_time(t,i), ...        
            d_o(i), ...                           
            d_i(i), ...                           
            Pump_values_time_list_m3_s(t) ...     
        );
        
        % 应用偏心系数 PR（扩展至7段）
        if pianxin == true
            if c_depth(i) <= depth_threshold1
                current_PR = PR_segment1;
            elseif c_depth(i) <= depth_threshold2
                current_PR = PR_segment2;
            elseif c_depth(i) <= depth_threshold3
                current_PR = PR_segment3;
            elseif c_depth(i) <= depth_threshold4
                current_PR = PR_segment4;
            elseif c_depth(i) <= depth_threshold5
                current_PR = PR_segment5;
            elseif c_depth(i) <= depth_threshold6
                current_PR = PR_segment6;
            else
                current_PR = PR_segment7;
            end
            Ff_a(t,i) = Ff_a(t,i) * current_PR;
        end
        
        if i == 1
            pressure_calculate = backpressure_time_list(t) * 1e6 + Ff_a(t,i) * structure_data.length_segment_array_m_(i) + 9.81 * rou_annulus_all_time_kg_m3(t,i) * vertical_length_all_grid(i); 
            pressure_annuli_static(t,i) = 9.81 * rou_annulus_all_time_kg_m3(t,i) * vertical_length_all_grid(i); 
            pressure_annuli_friction(t,i) = Ff_a(t,i) * structure_data.length_segment_array_m_(i); 
            pressure_annuli(t,i)=pressure_calculate; 
        elseif i>=2
            pressure_calculate = pressure_annuli(t,i-1) + Ff_a(t,i) * structure_data.length_segment_array_m_(i) + 9.81 * rou_annulus_all_time_kg_m3(t,i) * vertical_length_all_grid(i); 
            pressure_annuli_static(t,i) = pressure_annuli_static(t,i-1) + 9.81 * rou_annulus_all_time_kg_m3(t,i) * vertical_length_all_grid(i); 
            pressure_annuli_friction(t,i) = pressure_annuli_friction(t,i-1) + Ff_a(t,i) * structure_data.length_segment_array_m_(i); 
            pressure_annuli(t,i) = pressure_calculate; 
        end
     end
end    
    
pressure_annuli_MPa = pressure_annuli / 1000000;           
pressure_annuli_static_MPa = pressure_annuli_static / 1000000;  
pressure_annuli_friction_MPa = pressure_annuli_friction / 1000000;  
disp("环空压力计算完成")
%% 第五步循环：套管内压力计算（从井底向上反推，静液柱+循环压耗）
rou_casing_all_time = zeros(n_time, n_segment);    
miu_casing_all_time = zeros(n_time, n_segment);    
tau_casing_all_time = zeros(n_time, n_segment);    
velo_casing_all_time = zeros(n_time, n_segment);   
Ff_casing = zeros(n_time, n_segment);          
flow_pattern_casing = zeros(n_time, n_segment); 
volume_all_from_top_casing = zeros(1, n_segment);
length_all_from_top_casing = zeros(1, n_segment);
for i = 1:n_segment
    volume_all_from_top_casing(i) = sum(area_cin(1:i) .* structure_data.length_segment_array_m_(1:i)') * 1000;  
    length_all_from_top_casing(i) = sum(structure_data.length_segment_array_m_(1:i));  
end
   
for t = 2:n_time
    volume_injected_casing_1_list(1) = 1200;
    % 【高精度优化】：套管内注入体积累加也乘上极小的 dt 步长
    if volume_injected_all_list(t) > 0
        volume_injected_casing_1_list(t) = volume_injected_casing_1_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1
        volume_injected_casing_2_list(t) = volume_injected_casing_2_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2
        volume_injected_casing_3_list(t) = volume_injected_casing_3_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3
        volume_injected_casing_4_list(t) = volume_injected_casing_4_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4
        volume_injected_casing_5_list(t) = volume_injected_casing_5_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5
        volume_injected_casing_6_list(t) = volume_injected_casing_6_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6
        volume_injected_casing_7_list(t) = volume_injected_casing_7_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 
        volume_injected_casing_8_list(t) = volume_injected_casing_8_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8
        volume_injected_casing_9_list(t) = volume_injected_casing_9_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9
        volume_injected_casing_91_list(t) = volume_injected_casing_91_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v91
        volume_injected_casing_92_list(t) = volume_injected_casing_92_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v91 + v92
        volume_injected_casing_93_list(t) = volume_injected_casing_93_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v91 + v92 + v93
        volume_injected_casing_94_list(t) = volume_injected_casing_94_list(t-1)+Pump_values_time_list(t)*dt;
    end
    if volume_injected_all_list(t) > v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v91 + v92 + v93 + v94
        volume_injected_casing_95_list(t) = volume_injected_casing_95_list(t-1)+Pump_values_time_list(t)*dt;
    end
end
% 套管内界面追踪（从井口向下）
for t = 1:n_time  
    for i = 2:n_segment
        % 第一界面（钻井液-先导浆） 
        if volume_injected_casing_1_list(t)  >= volume_all_from_top_casing(i-1) && volume_injected_casing_1_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_0_1_list(t) = i;
            residual_volume_casing_0_1 = volume_all_from_top_casing(i) - volume_injected_casing_1_list(t);   
            residual_volume_casing_0_1_list(t) = residual_volume_casing_0_1;
            residual_height_casing_0_1 = residual_volume_casing_0_1 / (area_cin(i) * 1000);  
            residual_height_casing_0_1_list(t) = residual_height_casing_0_1;
        elseif volume_injected_casing_1_list(t) > 0 && volume_injected_casing_1_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_0_1_list(t) = 1;
            residual_volume_casing_0_1 = volume_injected_casing_1_list(t);
            residual_volume_casing_0_1_list(t) = residual_volume_casing_0_1;
            residual_height_casing_0_1 = residual_volume_casing_0_1 / (area_cin(1) * 1000);
            residual_height_casing_0_1_list(t) = residual_height_casing_0_1;
        elseif volume_injected_casing_1_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_0_1_list(t) = 10000;  
            residual_volume_casing_0_1_list(t) = 10000;
            residual_height_casing_0_1_list(t) = 10000;
        elseif volume_injected_casing_1_list(t) <= 0
            depth_tag_casing_0_1_list(t) = -1;  
            residual_volume_casing_0_1_list(t) = -1;
            residual_height_casing_0_1_list(t) = -1;
        end
        tag_interface_casing(t,1) = depth_tag_casing_0_1_list(t);
        
        % 第二界面（先导浆-低失水驱油隔离液1）
        if volume_injected_casing_2_list(t)  >= volume_all_from_top_casing(i-1) && volume_injected_casing_2_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_1_2_list(t) = i;
            residual_volume_casing_1_2 = volume_all_from_top_casing(i) - volume_injected_casing_2_list(t);
            residual_volume_casing_1_2_list(t) = residual_volume_casing_1_2;
            residual_height_casing_1_2 = residual_volume_casing_1_2 / (area_cin(i) * 1000);
            residual_height_casing_1_2_list(t) = residual_height_casing_1_2;
        elseif volume_injected_casing_2_list(t) > 0 && volume_injected_casing_2_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_1_2_list(t) = 1;
            residual_volume_casing_1_2 = volume_injected_casing_2_list(t);
            residual_volume_casing_1_2_list(t) = residual_volume_casing_1_2;
            residual_height_casing_1_2 = residual_volume_casing_1_2 / (area_cin(1) * 1000);
            residual_height_casing_1_2_list(t) = residual_height_casing_1_2;
        elseif volume_injected_casing_2_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_1_2_list(t) = 10000;
            residual_volume_casing_1_2_list(t) = 10000;
            residual_height_casing_1_2_list(t) = 10000;
        elseif volume_injected_casing_2_list(t) <= 0
            depth_tag_casing_1_2_list(t) = -1;
            residual_volume_casing_1_2_list(t) = -1;
            residual_height_casing_1_2_list(t) = -1;
        end
        tag_interface_casing(t,2) = depth_tag_casing_1_2_list(t);
        
        % 第三界面（低失水驱油隔离液1-低失水驱油隔离液2）
        if volume_injected_casing_3_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_3_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_2_3_list(t) = i;
            residual_volume_casing_2_3 = volume_all_from_top_casing(i) - volume_injected_casing_3_list(t);
            residual_volume_casing_2_3_list(t) = residual_volume_casing_2_3;
            residual_height_casing_2_3 = residual_volume_casing_2_3 / (area_cin(i) * 1000);
            residual_height_casing_2_3_list(t) = residual_height_casing_2_3;
        elseif volume_injected_casing_3_list(t) > 0 && volume_injected_casing_3_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_2_3_list(t) = 1;
            residual_volume_casing_2_3 = volume_injected_casing_3_list(t);
            residual_volume_casing_2_3_list(t) = residual_volume_casing_2_3;
            residual_height_casing_2_3 = residual_volume_casing_2_3 / (area_cin(1) * 1000);
            residual_height_casing_2_3_list(t) = residual_height_casing_2_3;
        elseif volume_injected_casing_3_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_2_3_list(t) = 10000;
            residual_volume_casing_2_3_list(t) = 10000;
            residual_height_casing_2_3_list(t) = 10000;
        elseif volume_injected_casing_3_list(t) <= 0
            depth_tag_casing_2_3_list(t) = -1;
            residual_volume_casing_2_3_list(t) = -1;
            residual_height_casing_2_3_list(t) = -1;
        end
        tag_interface_casing(t,3) = depth_tag_casing_2_3_list(t);
        
        % 第四界面（低失水驱油隔离液2-领浆）
        if volume_injected_casing_4_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_4_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_3_4_list(t) = i;
            residual_volume_casing_3_4 = volume_all_from_top_casing(i) - volume_injected_casing_4_list(t);
            residual_volume_casing_3_4_list(t) = residual_volume_casing_3_4;
            residual_height_casing_3_4 = residual_volume_casing_3_4 / (area_cin(i) * 1000);
            residual_height_casing_3_4_list(t) = residual_height_casing_3_4;
        elseif volume_injected_casing_4_list(t) > 0 && volume_injected_casing_4_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_3_4_list(t) = 1;
            residual_volume_casing_3_4 = volume_injected_casing_4_list(t);
            residual_volume_casing_3_4_list(t) = residual_volume_casing_3_4;
            residual_height_casing_3_4 = residual_volume_casing_3_4 / (area_cin(1) * 1000);
            residual_height_casing_3_4_list(t) = residual_height_casing_3_4;
        elseif volume_injected_casing_4_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_3_4_list(t) = 10000;
            residual_volume_casing_3_4_list(t) = 10000;
            residual_height_casing_3_4_list(t) = 10000;
        elseif volume_injected_casing_4_list(t) <= 0
            depth_tag_casing_3_4_list(t) = -1;
            residual_volume_casing_3_4_list(t) = -1;
            residual_height_casing_3_4_list(t) = -1;
        end
        tag_interface_casing(t,4) = depth_tag_casing_3_4_list(t);
        
        % 第五界面（领浆-尾浆）
        if volume_injected_casing_5_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_5_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_4_5_list(t) = i;
            residual_volume_casing_4_5 = volume_all_from_top_casing(i) - volume_injected_casing_5_list(t);
            residual_volume_casing_4_5_list(t) = residual_volume_casing_4_5;
            residual_height_casing_4_5 = residual_volume_casing_4_5 / (area_cin(i) * 1000);
            residual_height_casing_4_5_list(t) = residual_height_casing_4_5;
        elseif volume_injected_casing_5_list(t) > 0 && volume_injected_casing_5_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_4_5_list(t) = 1;
            residual_volume_casing_4_5 = volume_injected_casing_5_list(t);
            residual_volume_casing_4_5_list(t) = residual_volume_casing_4_5;
            residual_height_casing_4_5 = residual_volume_casing_4_5 / (area_cin(1) * 1000);
            residual_height_casing_4_5_list(t) = residual_height_casing_4_5;
        elseif volume_injected_casing_5_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_4_5_list(t) = 10000;
            residual_volume_casing_4_5_list(t) = 10000;
            residual_height_casing_4_5_list(t) = 10000;
        elseif volume_injected_casing_5_list(t) <= 0
            depth_tag_casing_4_5_list(t) = -1;
            residual_volume_casing_4_5_list(t) = -1;
            residual_height_casing_4_5_list(t) = -1;
        end
        tag_interface_casing(t,5) = depth_tag_casing_4_5_list(t);
        
        % 第六界面（尾浆-压塞液）
        if volume_injected_casing_6_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_6_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_5_6_list(t) = i;
            residual_volume_casing_5_6 = volume_all_from_top_casing(i) - volume_injected_casing_6_list(t);
            residual_volume_casing_5_6_list(t) = residual_volume_casing_5_6;
            residual_height_casing_5_6 = residual_volume_casing_5_6 / (area_cin(i) * 1000);
            residual_height_casing_5_6_list(t) = residual_height_casing_5_6;
        elseif volume_injected_casing_6_list(t) > 0 && volume_injected_casing_6_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_5_6_list(t) = 1;
            residual_volume_casing_5_6 = volume_injected_casing_6_list(t);
            residual_volume_casing_5_6_list(t) = residual_volume_casing_5_6;
            residual_height_casing_5_6 = residual_volume_casing_5_6 / (area_cin(1) * 1000);
            residual_height_casing_5_6_list(t) = residual_height_casing_5_6;
        elseif volume_injected_casing_6_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_5_6_list(t) = 10000;
            residual_volume_casing_5_6_list(t) = 10000;
            residual_height_casing_5_6_list(t) = 10000;
        elseif volume_injected_casing_6_list(t) <= 0
            depth_tag_casing_5_6_list(t) = -1;
            residual_volume_casing_5_6_list(t) = -1;
            residual_height_casing_5_6_list(t) = -1;
        end
        tag_interface_casing(t,6) = depth_tag_casing_5_6_list(t);
        
        % 第七界面（压塞液-泥浆）
        if volume_injected_casing_7_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_7_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_6_7_list(t) = i;
            residual_volume_casing_6_7 = volume_all_from_top_casing(i) - volume_injected_casing_7_list(t);
            residual_volume_casing_6_7_list(t) = residual_volume_casing_6_7;
            residual_height_casing_6_7 = residual_volume_casing_6_7 / (area_cin(i) * 1000);
            residual_height_casing_6_7_list(t) = residual_height_casing_6_7;
        elseif volume_injected_casing_7_list(t) > 0 && volume_injected_casing_7_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_6_7_list(t) = 1;
            residual_volume_casing_6_7 = volume_injected_casing_7_list(t);
            residual_volume_casing_6_7_list(t) = residual_volume_casing_6_7;
            residual_height_casing_6_7 = residual_volume_casing_6_7 / (area_cin(1) * 1000);
            residual_height_casing_6_7_list(t) = residual_height_casing_6_7;
        elseif volume_injected_casing_7_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_6_7_list(t) = 10000;
            residual_volume_casing_6_7_list(t) = 10000;
            residual_height_casing_6_7_list(t) = 10000;
        elseif volume_injected_casing_7_list(t) <= 0
            depth_tag_casing_6_7_list(t) = -1;
            residual_volume_casing_6_7_list(t) = -1;
            residual_height_casing_6_7_list(t) = -1;
        end
        tag_interface_casing(t,7) = depth_tag_casing_6_7_list(t);
        
         % 第八界面（替浆钻井液-保护液）
        if volume_injected_casing_8_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_8_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_7_8_list(t) = i;
            residual_volume_casing_7_8 = volume_all_from_top_casing(i) - volume_injected_casing_8_list(t);
            residual_volume_casing_7_8_list(t) = residual_volume_casing_7_8;
            residual_height_casing_7_8 = residual_volume_casing_7_8 / (area_cin(i) * 1000);
            residual_height_casing_7_8_list(t) = residual_height_casing_7_8;
        elseif volume_injected_casing_8_list(t) > 0 && volume_injected_casing_8_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_7_8_list(t) = 1;
            residual_volume_casing_7_8 = volume_injected_casing_8_list(t);
            residual_volume_casing_7_8_list(t) = residual_volume_casing_7_8;
            residual_height_casing_7_8 = residual_volume_casing_7_8 / (area_cin(1) * 1000);
            residual_height_casing_7_8_list(t) = residual_height_casing_7_8;
        elseif volume_injected_casing_8_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_7_8_list(t) = 10000;
            residual_volume_casing_7_8_list(t) = 10000;
            residual_height_casing_7_8_list(t) = 10000;
        elseif volume_injected_casing_8_list(t) <= 0
            depth_tag_casing_7_8_list(t) = -1;
            residual_volume_casing_7_8_list(t) = -1;
            residual_height_casing_7_8_list(t) = -1;
        end
        tag_interface_casing(t,8) = depth_tag_casing_7_8_list(t);

        % 第九界面（保护液-基液）
        if volume_injected_casing_9_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_9_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_8_9_list(t) = i;
            residual_volume_casing_8_9 = volume_all_from_top_casing(i) - volume_injected_casing_9_list(t);
            residual_volume_casing_8_9_list(t) = residual_volume_casing_8_9;
            residual_height_casing_8_9 = residual_volume_casing_8_9 / (area_cin(i) * 1000);
            residual_height_casing_8_9_list(t) = residual_height_casing_8_9;
        elseif volume_injected_casing_9_list(t) > 0 && volume_injected_casing_9_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_8_9_list(t) = 1;
            residual_volume_casing_8_9 = volume_injected_casing_9_list(t);
            residual_volume_casing_8_9_list(t) = residual_volume_casing_8_9;
            residual_height_casing_8_9 = residual_volume_casing_8_9 / (area_cin(1) * 1000);
            residual_height_casing_8_9_list(t) = residual_height_casing_8_9;
        elseif volume_injected_casing_9_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_8_9_list(t) = 10000;
            residual_volume_casing_8_9_list(t) = 10000;
            residual_height_casing_8_9_list(t) = 10000;
        elseif volume_injected_casing_9_list(t) <= 0
            depth_tag_casing_8_9_list(t) = -1;
            residual_volume_casing_8_9_list(t) = -1;
            residual_height_casing_8_9_list(t) = -1;
        end
        tag_interface_casing(t,9) = depth_tag_casing_8_9_list(t);

        % 第十界面（基液-替浆钻井液91）
        if volume_injected_casing_91_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_91_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_9_91_list(t) = i;
            residual_volume_casing_9_91 = volume_all_from_top_casing(i) - volume_injected_casing_91_list(t);
            residual_volume_casing_9_91_list(t) = residual_volume_casing_9_91;
            residual_height_casing_9_91 = residual_volume_casing_9_91 / (area_cin(i) * 1000);
            residual_height_casing_9_91_list(t) = residual_height_casing_9_91;
        elseif volume_injected_casing_91_list(t) > 0 && volume_injected_casing_91_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_9_91_list(t) = 1;
            residual_volume_casing_9_91 = volume_injected_casing_91_list(t);
            residual_volume_casing_9_91_list(t) = residual_volume_casing_9_91;
            residual_height_casing_9_91 = residual_volume_casing_9_91 / (area_cin(1) * 1000);
            residual_height_casing_9_91_list(t) = residual_height_casing_9_91;
        elseif volume_injected_casing_91_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_9_91_list(t) = 10000;
            residual_volume_casing_9_91_list(t) = 10000;
            residual_height_casing_9_91_list(t) = 10000;
        elseif volume_injected_casing_91_list(t) <= 0
            depth_tag_casing_9_91_list(t) = -1;
            residual_volume_casing_9_91_list(t) = -1;
            residual_height_casing_9_91_list(t) = -1;
        end
        tag_interface_casing(t,10) = depth_tag_casing_9_91_list(t);

        % 第十一界面（替浆钻井液91-替浆钻井液92）
        if volume_injected_casing_92_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_92_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_91_92_list(t) = i;
            residual_volume_casing_91_92 = volume_all_from_top_casing(i) - volume_injected_casing_92_list(t);
            residual_volume_casing_91_92_list(t) = residual_volume_casing_91_92;
            residual_height_casing_91_92 = residual_volume_casing_91_92 / (area_cin(i) * 1000);
            residual_height_casing_91_92_list(t) = residual_height_casing_91_92;
        elseif volume_injected_casing_92_list(t) > 0 && volume_injected_casing_92_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_91_92_list(t) = 1;
            residual_volume_casing_91_92 = volume_injected_casing_92_list(t);
            residual_volume_casing_91_92_list(t) = residual_volume_casing_91_92;
            residual_height_casing_91_92 = residual_volume_casing_91_92 / (area_cin(1) * 1000);
            residual_height_casing_91_92_list(t) = residual_height_casing_91_92;
        elseif volume_injected_casing_92_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_91_92_list(t) = 10000;
            residual_volume_casing_91_92_list(t) = 10000;
            residual_height_casing_91_92_list(t) = 10000;
        elseif volume_injected_casing_92_list(t) <= 0
            depth_tag_casing_91_92_list(t) = -1;
            residual_volume_casing_91_92_list(t) = -1;
            residual_height_casing_91_92_list(t) = -1;
        end
        tag_interface_casing(t,11) = depth_tag_casing_91_92_list(t);

        % 第十二界面（替浆钻井液92-替浆钻井液93）
        if volume_injected_casing_93_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_93_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_92_93_list(t) = i;
            residual_volume_casing_92_93 = volume_all_from_top_casing(i) - volume_injected_casing_93_list(t);
            residual_volume_casing_92_93_list(t) = residual_volume_casing_92_93;
            residual_height_casing_92_93 = residual_volume_casing_92_93 / (area_cin(i) * 1000);
            residual_height_casing_92_93_list(t) = residual_height_casing_92_93;
        elseif volume_injected_casing_93_list(t) > 0 && volume_injected_casing_93_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_92_93_list(t) = 1;
            residual_volume_casing_92_93 = volume_injected_casing_93_list(t);
            residual_volume_casing_92_93_list(t) = residual_volume_casing_92_93;
            residual_height_casing_92_93 = residual_volume_casing_92_93 / (area_cin(1) * 1000);
            residual_height_casing_92_93_list(t) = residual_height_casing_92_93;
        elseif volume_injected_casing_93_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_92_93_list(t) = 10000;
            residual_volume_casing_92_93_list(t) = 10000;
            residual_height_casing_92_93_list(t) = 10000;
        elseif volume_injected_casing_93_list(t) <= 0
            depth_tag_casing_92_93_list(t) = -1;
            residual_volume_casing_92_93_list(t) = -1;
            residual_height_casing_92_93_list(t) = -1;
        end
        tag_interface_casing(t,12) = depth_tag_casing_92_93_list(t);

        % 第十三界面（替浆钻井液93-替浆钻井液94）
        if volume_injected_casing_94_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_94_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_93_94_list(t) = i;
            residual_volume_casing_93_94 = volume_all_from_top_casing(i) - volume_injected_casing_94_list(t);
            residual_volume_casing_93_94_list(t) = residual_volume_casing_93_94;
            residual_height_casing_93_94 = residual_volume_casing_93_94 / (area_cin(i) * 1000);
            residual_height_casing_93_94_list(t) = residual_height_casing_93_94;
        elseif volume_injected_casing_94_list(t) > 0 && volume_injected_casing_94_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_93_94_list(t) = 1;
            residual_volume_casing_93_94 = volume_injected_casing_94_list(t);
            residual_volume_casing_93_94_list(t) = residual_volume_casing_93_94;
            residual_height_casing_93_94 = residual_volume_casing_93_94 / (area_cin(1) * 1000);
            residual_height_casing_93_94_list(t) = residual_height_casing_93_94;
        elseif volume_injected_casing_94_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_93_94_list(t) = 10000;
            residual_volume_casing_93_94_list(t) = 10000;
            residual_height_casing_93_94_list(t) = 10000;
        elseif volume_injected_casing_94_list(t) <= 0
            depth_tag_casing_93_94_list(t) = -1;
            residual_volume_casing_93_94_list(t) = -1;
            residual_height_casing_93_94_list(t) = -1;
        end
        tag_interface_casing(t,13) = depth_tag_casing_93_94_list(t);

        % 第十四界面（替浆钻井液94-替浆钻井液95）
        if volume_injected_casing_95_list(t) >= volume_all_from_top_casing(i-1) && volume_injected_casing_95_list(t) < volume_all_from_top_casing(i)
            depth_tag_casing_94_95_list(t) = i;
            residual_volume_casing_94_95 = volume_all_from_top_casing(i) - volume_injected_casing_95_list(t);
            residual_volume_casing_94_95_list(t) = residual_volume_casing_94_95;
            residual_height_casing_94_95 = residual_volume_casing_94_95 / (area_cin(i) * 1000);
            residual_height_casing_94_95_list(t) = residual_height_casing_94_95;
        elseif volume_injected_casing_95_list(t) > 0 && volume_injected_casing_95_list(t) < volume_all_from_top_casing(1)
            depth_tag_casing_94_95_list(t) = 1;
            residual_volume_casing_94_95 = volume_injected_casing_95_list(t);
            residual_volume_casing_94_95_list(t) = residual_volume_casing_94_95;
            residual_height_casing_94_95 = residual_volume_casing_94_95 / (area_cin(1) * 1000);
            residual_height_casing_94_95_list(t) = residual_height_casing_94_95;
        elseif volume_injected_casing_95_list(t) >= volume_all_from_top_casing(end)
            depth_tag_casing_94_95_list(t) = 10000;
            residual_volume_casing_94_95_list(t) = 10000;
            residual_height_casing_94_95_list(t) = 10000;
        elseif volume_injected_casing_95_list(t) <= 0
            depth_tag_casing_94_95_list(t) = -1;
            residual_volume_casing_94_95_list(t) = -1;
            residual_height_casing_94_95_list(t) = -1;
        end
        tag_interface_casing(t,14) = depth_tag_casing_94_95_list(t);
    end
    
% --------- 2. 套管内流体属性分配（密度、黏度、屈服应力） --------------------------
    for i = 1:n_segment
        % 根据 tag_interface_casing 确定扩增后13段流体分布
        if i > depth_tag_casing_0_1_list(t)
            rou_casing_all_time(t, i) = rou0;
            miu_casing_all_time(t, i) = miu0;
            tau_casing_all_time(t, i) = tau0;
        elseif i == tag_interface_casing(t,1)
            proportion_casing_0_1_1 = residual_height_casing_0_1_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_0_1_1 = min(max(proportion_casing_0_1_1, 0), 1); 
            proportion_casing_0_1_0 = 1 - proportion_casing_0_1_1;
            rou_casing_all_time(t, i) = proportion_casing_0_1_1 * rou0 + proportion_casing_0_1_0 * rou1;
            miu_casing_all_time(t, i) = proportion_casing_0_1_1 * miu0 + proportion_casing_0_1_0 * miu1;
            tau_casing_all_time(t, i) = proportion_casing_0_1_1 * tau0 + proportion_casing_0_1_0 * tau1;
        elseif i<tag_interface_casing(t,1) && i>tag_interface_casing(t,2)
            rou_casing_all_time(t, i) = rou1;
            miu_casing_all_time(t, i) = miu1;
            tau_casing_all_time(t, i) = tau1;
        elseif i == tag_interface_casing(t,2)
            proportion_casing_1_2_1 = residual_height_casing_1_2_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_1_2_1 = min(max(proportion_casing_1_2_1, 0), 1);
            proportion_casing_1_2_0 = 1 - proportion_casing_1_2_1;
            rou_casing_all_time(t, i) = proportion_casing_1_2_1 * rou1 + proportion_casing_1_2_0 * rou2;
            miu_casing_all_time(t, i) = proportion_casing_1_2_1 * miu1 + proportion_casing_1_2_0 * miu2;
            tau_casing_all_time(t, i) = proportion_casing_1_2_1 * tau1 + proportion_casing_1_2_0 * tau2;
        elseif i<tag_interface_casing(t,2) && i>tag_interface_casing(t,3)
            rou_casing_all_time(t, i) = rou2;
            miu_casing_all_time(t, i) = miu2;
            tau_casing_all_time(t, i) = tau2;
        elseif i == tag_interface_casing(t,3)
            proportion_casing_2_3_1 = residual_height_casing_2_3_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_2_3_1 = min(max(proportion_casing_2_3_1, 0), 1);
            proportion_casing_2_3_0 = 1 - proportion_casing_2_3_1;
            rou_casing_all_time(t, i) = proportion_casing_2_3_1 * rou2 + proportion_casing_2_3_0 * rou3;
            miu_casing_all_time(t, i) = proportion_casing_2_3_1 * miu2 + proportion_casing_2_3_0 * miu3;
            tau_casing_all_time(t, i) = proportion_casing_2_3_1 * tau1 + proportion_casing_2_3_0 * tau3;
        elseif i<tag_interface_casing(t,3) && i>tag_interface_casing(t,4)
            rou_casing_all_time(t, i) = rou3;
            miu_casing_all_time(t, i) = miu3;
            tau_casing_all_time(t, i) = tau3;
        elseif i == tag_interface_casing(t,4)
            proportion_casing_3_4_1 = residual_height_casing_3_4_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_3_4_1 = min(max(proportion_casing_3_4_1, 0), 1);
            proportion_casing_3_4_0 = 1 - proportion_casing_3_4_1;
            rou_casing_all_time(t, i) = proportion_casing_3_4_1 * rou3 + proportion_casing_3_4_0 * rou4;
            miu_casing_all_time(t, i) = proportion_casing_3_4_1 * miu3 + proportion_casing_3_4_0 * miu4;
            tau_casing_all_time(t, i) = proportion_casing_3_4_1 * tau3 + proportion_casing_3_4_0 * tau4;
        elseif i<tag_interface_casing(t,4) && i>tag_interface_casing(t,5)
            rou_casing_all_time(t, i) = rou4;
            miu_casing_all_time(t, i) = miu4;
            tau_casing_all_time(t, i) = tau4;
        elseif i == tag_interface_casing(t,5)
            proportion_casing_4_5_1 = residual_height_casing_4_5_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_4_5_1 = min(max(proportion_casing_4_5_1, 0), 1);
            proportion_casing_4_5_0 = 1 - proportion_casing_4_5_1;
            rou_casing_all_time(t, i) = proportion_casing_4_5_1 * rou4 + proportion_casing_4_5_0 * rou5;
            miu_casing_all_time(t, i) = proportion_casing_4_5_1 * miu4 + proportion_casing_4_5_0 * miu5;
            tau_casing_all_time(t, i) = proportion_casing_4_5_1 * tau4 + proportion_casing_4_5_0 * tau5;
        elseif i<tag_interface_casing(t,5) && i>tag_interface_casing(t,6)
            rou_casing_all_time(t, i) = rou5;
            miu_casing_all_time(t, i) = miu5;
            tau_casing_all_time(t, i) = tau5;
        elseif i == tag_interface_casing(t,6)
            proportion_casing_5_6_1 = residual_height_casing_5_6_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_5_6_1 = min(max(proportion_casing_5_6_1, 0), 1);
            proportion_casing_5_6_0 = 1 - proportion_casing_5_6_1;
            rou_casing_all_time(t, i) = proportion_casing_5_6_1 * rou5 + proportion_casing_5_6_0 * rou6;
            miu_casing_all_time(t, i) = proportion_casing_5_6_1 * miu5 + proportion_casing_5_6_0 * miu6;
            tau_casing_all_time(t, i) = proportion_casing_5_6_1 * tau5 + proportion_casing_5_6_0 * tau6;
        elseif i<tag_interface_casing(t,6) && i>tag_interface_casing(t,7)
            rou_casing_all_time(t, i) = rou6;
            miu_casing_all_time(t, i) = miu6;
            tau_casing_all_time(t, i) = tau6;
        elseif i == tag_interface_casing(t,7)
            proportion_casing_6_7_1 = residual_height_casing_6_7_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_6_7_1 = min(max(proportion_casing_6_7_1, 0), 1);
            proportion_casing_6_7_0 = 1 - proportion_casing_6_7_1;
            rou_casing_all_time(t, i) = proportion_casing_6_7_1 * rou6 + proportion_casing_6_7_0 * rou7;
            miu_casing_all_time(t, i) = proportion_casing_6_7_1 * miu6 + proportion_casing_6_7_0 * miu7;
            tau_casing_all_time(t, i) = proportion_casing_6_7_1 * tau6 + proportion_casing_6_7_0 * tau7;
        elseif i<tag_interface_casing(t,7) && i>tag_interface_casing(t,8)
            rou_casing_all_time(t, i) = rou7;
            miu_casing_all_time(t, i) = miu7;
            tau_casing_all_time(t, i) = tau7;
        elseif i == tag_interface_casing(t,8)
            proportion_casing_7_8_1 = residual_height_casing_7_8_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_7_8_1 = min(max(proportion_casing_7_8_1, 0), 1);
            proportion_casing_7_8_0 = 1 - proportion_casing_7_8_1;
            rou_casing_all_time(t, i) = proportion_casing_7_8_1 * rou7 + proportion_casing_7_8_0 * rou8;
            miu_casing_all_time(t, i) = proportion_casing_7_8_1 * miu7 + proportion_casing_7_8_0 * miu8;
            tau_casing_all_time(t, i) = proportion_casing_7_8_1 * tau7 + proportion_casing_7_8_0 * tau8;
        elseif i<tag_interface_casing(t,8) && i>tag_interface_casing(t,9)
            rou_casing_all_time(t, i) = rou8;
            miu_casing_all_time(t, i) = miu8;
            tau_casing_all_time(t, i) = tau8;
        elseif i == tag_interface_casing(t,9)
            proportion_casing_8_9_1 = residual_height_casing_8_9_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_8_9_1 = min(max(proportion_casing_8_9_1, 0), 1);
            proportion_casing_8_9_0 = 1 - proportion_casing_8_9_1;
            rou_casing_all_time(t, i) = proportion_casing_8_9_1 * rou8 + proportion_casing_8_9_0 * rou9;
            miu_casing_all_time(t, i) = proportion_casing_8_9_1 * miu8 + proportion_casing_8_9_0 * miu9;
            tau_casing_all_time(t, i) = proportion_casing_8_9_1 * tau8 + proportion_casing_8_9_0 * tau9;
        elseif i<tag_interface_casing(t,9) && i>tag_interface_casing(t,10)
            rou_casing_all_time(t, i) = rou9;
            miu_casing_all_time(t, i) = miu9;
            tau_casing_all_time(t, i) = tau9;
        elseif i == tag_interface_casing(t,10)
            proportion_casing_9_91_1 = residual_height_casing_9_91_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_9_91_1 = min(max(proportion_casing_9_91_1, 0), 1);
            proportion_casing_9_91_0 = 1 - proportion_casing_9_91_1;
            rou_casing_all_time(t, i) = proportion_casing_9_91_1 * rou9 + proportion_casing_9_91_0 * rou91;
            miu_casing_all_time(t, i) = proportion_casing_9_91_1 * miu9 + proportion_casing_9_91_0 * miu91;
            tau_casing_all_time(t, i) = proportion_casing_9_91_1 * tau9 + proportion_casing_9_91_0 * tau91;
        elseif i<tag_interface_casing(t,10) && i>tag_interface_casing(t,11)
            rou_casing_all_time(t, i) = rou91;
            miu_casing_all_time(t, i) = miu91;
            tau_casing_all_time(t, i) = tau91;
        elseif i == tag_interface_casing(t,11)
            proportion_casing_91_92_1 = residual_height_casing_91_92_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_91_92_1 = min(max(proportion_casing_91_92_1, 0), 1);
            proportion_casing_91_92_0 = 1 - proportion_casing_91_92_1;
            rou_casing_all_time(t, i) = proportion_casing_91_92_1 * rou91 + proportion_casing_91_92_0 * rou92;
            miu_casing_all_time(t, i) = proportion_casing_91_92_1 * miu91 + proportion_casing_91_92_0 * miu92;
            tau_casing_all_time(t, i) = proportion_casing_91_92_1 * tau91 + proportion_casing_91_92_0 * tau92;
        elseif i<tag_interface_casing(t,11) && i>tag_interface_casing(t,12)
            rou_casing_all_time(t, i) = rou92;
            miu_casing_all_time(t, i) = miu92;
            tau_casing_all_time(t, i) = tau92;
        elseif i == tag_interface_casing(t,12)
            proportion_casing_92_93_1 = residual_height_casing_92_93_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_92_93_1 = min(max(proportion_casing_92_93_1, 0), 1);
            proportion_casing_92_93_0 = 1 - proportion_casing_92_93_1;
            rou_casing_all_time(t, i) = proportion_casing_92_93_1 * rou92 + proportion_casing_92_93_0 * rou93;
            miu_casing_all_time(t, i) = proportion_casing_92_93_1 * miu92 + proportion_casing_92_93_0 * miu93;
            tau_casing_all_time(t, i) = proportion_casing_92_93_1 * tau92 + proportion_casing_92_93_0 * tau93;
        elseif i<tag_interface_casing(t,12) && i>tag_interface_casing(t,13)
            rou_casing_all_time(t, i) = rou93;
            miu_casing_all_time(t, i) = miu93;
            tau_casing_all_time(t, i) = tau93;
        elseif i == tag_interface_casing(t,13)
            proportion_casing_93_94_1 = residual_height_casing_93_94_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_93_94_1 = min(max(proportion_casing_93_94_1, 0), 1);
            proportion_casing_93_94_0 = 1 - proportion_casing_93_94_1;
            rou_casing_all_time(t, i) = proportion_casing_93_94_1 * rou93 + proportion_casing_93_94_0 * rou94;
            miu_casing_all_time(t, i) = proportion_casing_93_94_1 * miu93 + proportion_casing_93_94_0 * miu94;
            tau_casing_all_time(t, i) = proportion_casing_93_94_1 * tau93 + proportion_casing_93_94_0 * tau94;
        elseif i<tag_interface_casing(t,13) && i>tag_interface_casing(t,14)
            rou_casing_all_time(t, i) = rou94;
            miu_casing_all_time(t, i) = miu94;
            tau_casing_all_time(t, i) = tau94;
        elseif i == tag_interface_casing(t,14)
            proportion_casing_94_95_1 = residual_height_casing_94_95_list(t) / structure_data.length_segment_array_m_(i);
            proportion_casing_94_95_1 = min(max(proportion_casing_94_95_1, 0), 1);
            proportion_casing_94_95_0 = 1 - proportion_casing_94_95_1;
            rou_casing_all_time(t, i) = proportion_casing_94_95_1 * rou94 + proportion_casing_94_95_0 * rou95;
            miu_casing_all_time(t, i) = proportion_casing_94_95_1 * miu94 + proportion_casing_94_95_0 * miu95;
            tau_casing_all_time(t, i) = proportion_casing_94_95_1 * tau94 + proportion_casing_94_95_0 * tau95;
        elseif i<tag_interface_casing(t,14)
            rou_casing_all_time(t, i) = rou95;
            miu_casing_all_time(t, i) = miu95;
            tau_casing_all_time(t, i) = tau95;
        end
        
        % 套管内流速计算
        if Pump_values_time_list(t) > 0
            velo_casing_all_time(t, i) = Pump_values_time_list(t) / 1000 / 60 / area_cin(i);  % m/s
        else
            velo_casing_all_time(t, i) = 0;  % 无流动
        end
    end
end
% 密度单位转换（g/cm³→kg/m³）
rou_casing_all_time_kg_m3 = rou_casing_all_time * 1000;

%% 套管内压力计算循环（从井底向上反推）
%套管静压
for t = 1:n_time
    for i = 1:n_segment
     if i == 1
     pressure_casing_static(t,i) = 9.81 * rou_casing_all_time_kg_m3(t,i) * vertical_length_all_grid(i); % pa
     else
     pressure_casing_static(t,i) = pressure_casing_static(t,i-1) + 9.81 * rou_casing_all_time_kg_m3(t,i) * vertical_length_all_grid(i); % pa
     end
     pressure_casing_static_MPa(t,i) = pressure_casing_static(t,i)/1000000; % MPa
    end
end

for t = 1:n_time
    % 首先计算套管内各段的摩擦阻力
    for i = 1:n_segment
        [Ff_casing(t,i), flow_pattern_casing(t,i)] = Friction_casing_bh(...
            rou_casing_all_time_kg_m3(t,i), ...  
            velo_casing_all_time(t,i), ...       
            area_cin(i), ...                     
            miu_casing_all_time(t,i), ...        
            tau_casing_all_time(t,i), ...        
            diameter_casing_in(i), ...           
            0, ...                               
            Pump_values_time_list_m3_s(t) ...    
        );
    end
    
    % 从井底向上反推套管内压力
    if Pump_values_time_list(t) > 0
        pressure_casing_bottom = pressure_annuli(t, n_segment) + P_bit_drop;  
        pressure_casing(t, n_segment) = pressure_casing_bottom; 
        
        for i = n_segment-1:-1:1    
            pressure_casing(t, i) = pressure_casing(t, i+1) + Ff_casing(t,i+1) * structure_data.length_segment_array_m_(i+1) - (pressure_casing_static(t,i+1) - pressure_casing_static(t,i));    
        end
        
        for i=2:n_segment
            pressure_casing_friction(t,i) = pressure_casing_friction(t,i-1) + Ff_casing(t,i) * structure_data.length_segment_array_m_(i);                 
        end
        
        true_surface_pressure_Pa = pressure_casing(t, 1) + Ff_casing(t, 1) * structure_data.length_segment_array_m_(1) - pressure_casing_static(t, 1);
        pressure_pump_surface(t) = true_surface_pressure_Pa / 1e6;  % 转换为MPa
        pressure_casing_friction_MPa=pressure_casing_friction/1000000; %转换为MPa
        if pressure_pump_surface(t) < 0
           pressure_pump_surface(t) = 0;  % 修正自由下落时的负压
        end
    else
        pressure_casing(t, n_segment) = pressure_annuli(t, n_segment) + P_bit_drop;
        for i = n_segment-1:-1:1
            pressure_casing(t, i) = pressure_casing(t, i+1) - ...
            9.81 * rou_casing_all_time_kg_m3(t, i+1) * vertical_length_all_grid(i+1);
        end
        pressure_pump_surface(t) = pressure_casing(t, 1) / 1e6;
    end
end

pressure_casing_MPa = pressure_casing / 1000000;
fprintf('套管内压力计算完成（从井底向上反推）！\n');

%% 套管内界面深度计算（用于绘图）
depth_every_interface_casing = NaN(size(tag_interface_casing));
tag_interface_casing2 = NaN(size(tag_interface_casing));
for i = 1:numel(depth_every_interface_casing(:,1))
    if tag_interface_casing(i,1) < 10000 && tag_interface_casing(i,1) > -1
        depth_every_interface_casing(i,1) = structure_data.depth_well_logging_m_(tag_interface_casing(i,1)) - residual_height_casing_0_1_list(i);
        tag_interface_casing2(i,1) = tag_interface_casing(i,1);
    end
    if tag_interface_casing(i,2) < 10000 && tag_interface_casing(i,2) > -1
        depth_every_interface_casing(i,2) = structure_data.depth_well_logging_m_(tag_interface_casing(i,2)) - residual_height_casing_1_2_list(i);
        tag_interface_casing2(i,2) = tag_interface_casing(i,2);
    end
    if tag_interface_casing(i,3) < 10000 && tag_interface_casing(i,3) > -1
        depth_every_interface_casing(i,3) = structure_data.depth_well_logging_m_(tag_interface_casing(i,3)) - residual_height_casing_2_3_list(i);
        tag_interface_casing2(i,3) = tag_interface_casing(i,3);
    end
    if tag_interface_casing(i,4) < 10000 && tag_interface_casing(i,4) > -1
        depth_every_interface_casing(i,4) = structure_data.depth_well_logging_m_(tag_interface_casing(i,4)) - residual_height_casing_3_4_list(i);
        tag_interface_casing2(i,4) = tag_interface_casing(i,4);
    end
    if tag_interface_casing(i,5) < 10000 && tag_interface_casing(i,5) > -1
        depth_every_interface_casing(i,5) = structure_data.depth_well_logging_m_(tag_interface_casing(i,5)) - residual_height_casing_4_5_list(i);
        tag_interface_casing2(i,5) = tag_interface_casing(i,5);
    end
    if tag_interface_casing(i,6) < 10000 && tag_interface_casing(i,6) > -1
        depth_every_interface_casing(i,6) = structure_data.depth_well_logging_m_(tag_interface_casing(i,6)) - residual_height_casing_5_6_list(i);
        tag_interface_casing2(i,6) = tag_interface_casing(i,6);
    end
    if tag_interface_casing(i,7) < 10000 && tag_interface_casing(i,7) > -1
        depth_every_interface_casing(i,7) = structure_data.depth_well_logging_m_(tag_interface_casing(i,7)) - residual_height_casing_6_7_list(i);
        tag_interface_casing2(i,7) = tag_interface_casing(i,7);
    end
    if tag_interface_casing(i,8) < 10000 && tag_interface_casing(i,8) > -1
        depth_every_interface_casing(i,8) = structure_data.depth_well_logging_m_(tag_interface_casing(i,8)) - residual_height_casing_7_8_list(i);
        tag_interface_casing2(i,8) = tag_interface_casing(i,8);
    end
    if tag_interface_casing(i,9) < 10000 && tag_interface_casing(i,9) > -1
        depth_every_interface_casing(i,9) = structure_data.depth_well_logging_m_(tag_interface_casing(i,9)) - residual_height_casing_8_9_list(i);
        tag_interface_casing2(i,9) = tag_interface_casing(i,9);
    end
    if tag_interface_casing(i,10) < 10000 && tag_interface_casing(i,10) > -1
        depth_every_interface_casing(i,10) = structure_data.depth_well_logging_m_(tag_interface_casing(i,10)) - residual_height_casing_9_91_list(i);
        tag_interface_casing2(i,10) = tag_interface_casing(i,10);
    end
    if tag_interface_casing(i,11) < 10000 && tag_interface_casing(i,11) > -1
        depth_every_interface_casing(i,11) = structure_data.depth_well_logging_m_(tag_interface_casing(i,11)) - residual_height_casing_91_92_list(i);
        tag_interface_casing2(i,11) = tag_interface_casing(i,11);
    end
    if tag_interface_casing(i,12) < 10000 && tag_interface_casing(i,12) > -1
        depth_every_interface_casing(i,12) = structure_data.depth_well_logging_m_(tag_interface_casing(i,12)) - residual_height_casing_92_93_list(i);
        tag_interface_casing2(i,12) = tag_interface_casing(i,12);
    end
    if tag_interface_casing(i,13) < 10000 && tag_interface_casing(i,13) > -1
        depth_every_interface_casing(i,13) = structure_data.depth_well_logging_m_(tag_interface_casing(i,13)) - residual_height_casing_93_94_list(i);
        tag_interface_casing2(i,13) = tag_interface_casing(i,13);
    end
    if tag_interface_casing(i,14) < 10000 && tag_interface_casing(i,14) > -1
        depth_every_interface_casing(i,14) = structure_data.depth_well_logging_m_(tag_interface_casing(i,14)) - residual_height_casing_94_95_list(i);
        tag_interface_casing2(i,14) = tag_interface_casing(i,14);
    end
end
interface_depth_casing = {};
% 深度数据
for i = 1:14
    x = find(~isnan(tag_interface_casing2(:,i)));
    y = depth_every_interface_casing(:,i);
    interface_depth_casing{i,1} = x * dt;%时间 
    interface_depth_casing{i,2} = y(~isnan(y));%深度
end

%% 环空ECD计算模块
TVD_cum(TVD_cum == 0) = eps;  % 避免除零
ECD_kg_m3 = zeros(n_time, n_segment);      
ECD_g_cm3 = zeros(n_time, n_segment);
ESD_g_cm3 = zeros(n_time, n_segment);

for t = 1:n_time
    for i = 1:n_segment
        if TVD_cum(i) > 0
            ECD_kg_m3(t, i) = pressure_annuli(t, i) / (9.81 * TVD_cum(i));
            ECD_g_cm3(t, i) = ECD_kg_m3(t, i) / 1000;  
            ESD_g_cm3(t, i) = pressure_annuli_static(t, i) / (9.81 * TVD_cum(i))/ 1000;  % kg/m³ → g/cm³

        else
            ECD_kg_m3(t, i) = 0;
            ECD_g_cm3(t, i) = 0;
        end
    end
end
disp("环空ECD计算完成！")

%% 套管内ECD计算模块
ECD_casing_kg_m3 = zeros(n_time, n_segment);
ECD_casing_g_cm3 = zeros(n_time, n_segment);
for t = 1:n_time
    for i = 1:n_segment
        if TVD_cum(i) > 0
            ECD_casing_kg_m3(t, i) = pressure_casing(t, i) / (9.81 * TVD_cum(i));
            ECD_casing_g_cm3(t, i) = ECD_casing_kg_m3(t, i) / 1000; 
        else
            ECD_casing_kg_m3(t, i) = 0;
            ECD_casing_g_cm3(t, i) = 0;
        end
    end
end
disp("套管内ECD模块计算完成！")

%% === 井底压力曲线绘图模块 ===  各关注点安全密度窗口来自HT1-004施工设计
safe_ECD_5568_lower = 1.940; safe_ECD_5568_upper = 1.975; % 5578m套管鞋
safe_ECD_6880_lower = 1.940; safe_ECD_6880_upper = 1.975; % 6600m关注点
safe_ECD_7463_lower = 1.940; safe_ECD_7463_upper = 1.975; % 7498m漏层
safe_ECD_bottom_lower = 2.010; safe_ECD_bottom_upper = 2.050; % 7660m井底TD（统一安全密度窗口2.01-2.05）
safety_margin = 0.003; % 安全余量 g/cm³，回压裁剪时离边界的最小余量
% 四个关键关注点索引：5578m(套管鞋)、6600m(关注点)、7498m(漏层)、井底7660m。
% CSV按30m一个点生成，非整30m关注点采用最近网格点计算。
[~, idx_5568] = min(abs(structure_data.depth_well_logging_m_ - 5578));
[~, idx_6880] = min(abs(structure_data.depth_well_logging_m_ - 6600));
[~, idx_7463] = min(abs(structure_data.depth_well_logging_m_ - 7498));
idx_bottom_critical = n_segment;
if isempty(idx_5568) || isempty(idx_6880) || isempty(idx_7463)
    error('井身结构CSV缺少5578m、6600m或7498m关注点，请检查呼1-004井身结构.csv。');
end
key_indices = [idx_5568, idx_6880, idx_7463, idx_bottom_critical];
key_names = {sprintf('5578m套管鞋(计算点%.0fm)', structure_data.depth_well_logging_m_(idx_5568)), ...
             sprintf('6600m关注点(计算点%.0fm)', structure_data.depth_well_logging_m_(idx_6880)), ...
             sprintf('7498m漏层(计算点%.0fm)', structure_data.depth_well_logging_m_(idx_7463)), ...
             '7660m井底TD'};
TVD_bottom = structure_data.vertical_depth_for_logging_m_(end); 
time_minutes_real = (1:n_time)' * dt;

safe_P_5568_lower_line = ones(n_time, 1) * safe_ECD_5568_lower * 0.00981 * TVD_cum(idx_5568);
safe_P_5568_upper_line = ones(n_time, 1) * safe_ECD_5568_upper * 0.00981 * TVD_cum(idx_5568);
safe_P_6880_lower_line = ones(n_time, 1) * safe_ECD_6880_lower * 0.00981 * TVD_cum(idx_6880);
safe_P_6880_upper_line = ones(n_time, 1) * safe_ECD_6880_upper * 0.00981 * TVD_cum(idx_6880);
safe_P_7463_lower_line = ones(n_time, 1) * safe_ECD_7463_lower * 0.00981 * TVD_cum(idx_7463);
safe_P_7463_upper_line = ones(n_time, 1) * safe_ECD_7463_upper * 0.00981 * TVD_cum(idx_7463);
safe_P_bottom_lower_line = ones(n_time, 1) * safe_ECD_bottom_lower * 0.00981 * TVD_bottom;
safe_P_bottom_upper_line = ones(n_time, 1) * safe_ECD_bottom_upper * 0.00981 * TVD_bottom;

safe_ECD_5568_lower_line = ones(n_time, 1) * safe_ECD_5568_lower;
safe_ECD_5568_upper_line = ones(n_time, 1) * safe_ECD_5568_upper;
safe_ECD_6880_lower_line = ones(n_time, 1) * safe_ECD_6880_lower;
safe_ECD_6880_upper_line = ones(n_time, 1) * safe_ECD_6880_upper;
safe_ECD_7463_lower_line = ones(n_time, 1) * safe_ECD_7463_lower;
safe_ECD_7463_upper_line = ones(n_time, 1) * safe_ECD_7463_upper;
safe_ECD_bottom_lower_line = ones(n_time, 1) * safe_ECD_bottom_lower;
safe_ECD_bottom_upper_line = ones(n_time, 1) * safe_ECD_bottom_upper;

% ---------------- 绘图 1: 井底压力 (MPa) 与安全密度窗口 ----------------
figure('Name', '井底压力与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, pressure_annuli_MPa(:, end), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_P_bottom_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_P_bottom_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('井底压力 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('井底压力随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空井底压力', sprintf('安全下限 (ECD=%.3f)', safe_ECD_bottom_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_bottom_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

% ---------------- 绘图 2: 6600m关注点压力 (MPa) 与安全窗口 ----------------
figure('Name', '6600m关注点压力与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, pressure_annuli_MPa(:, idx_6880), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_P_6880_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_P_6880_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('6600m关注点压力 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('6600m关注点压力随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空压力', sprintf('安全下限 (ECD=%.3f)', safe_ECD_6880_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_6880_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
ylim([120, 150]);

% ---------------- 绘图 3: 井底 ECD (g/cm³) 与安全窗口 ----------------
figure('Name', '井底ECD与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, ECD_g_cm3(:, end), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_ECD_bottom_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_ECD_bottom_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('井底 ECD (g/cm³)', 'FontSize', 12, 'FontWeight', 'bold');
title('井底 ECD 随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空井底 ECD', sprintf('安全下限 (ECD=%.3f)', safe_ECD_bottom_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_bottom_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

% ---------------- 绘图 4: 6600m关注点 ECD (g/cm³) 与安全窗口 ----------------
figure('Name', '6600m关注点ECD与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, ECD_g_cm3(:, idx_6880), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_ECD_6880_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_ECD_6880_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('6600m关注点ECD (g/cm³)', 'FontSize', 12, 'FontWeight', 'bold');
title('6600m关注点 ECD 随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空 ECD', sprintf('安全下限 (ECD=%.3f)', safe_ECD_6880_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_6880_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
% ---------------- 绘图 5: 5578m套管鞋压力 (MPa) 与安全窗口 ----------------
figure('Name', '5578m套管鞋压力与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, pressure_annuli_MPa(:, idx_5568), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_P_5568_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_P_5568_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('5578m套管鞋压力 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('5578m套管鞋压力随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空压力', sprintf('安全下限 (ECD=%.3f)', safe_ECD_5568_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_5568_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
% ---------------- 绘图 6: 5578m套管鞋 ECD (g/cm³) 与安全窗口 ----------------
figure('Name', '5578m套管鞋ECD与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, ECD_g_cm3(:, idx_5568), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_ECD_5568_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_ECD_5568_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('5578m套管鞋ECD (g/cm³)', 'FontSize', 12, 'FontWeight', 'bold');
title('5578m套管鞋 ECD 随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空 ECD', sprintf('安全下限 (ECD=%.3f)', safe_ECD_5568_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_5568_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
% ---------------- 绘图 7: 7498m漏层压力 (MPa) 与安全窗口 ----------------
figure('Name', '7498m漏层压力与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, pressure_annuli_MPa(:, idx_7463), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_P_7463_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_P_7463_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('7498m漏层压力 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('7498m漏层压力随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空压力', sprintf('安全下限 (ECD=%.3f)', safe_ECD_7463_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_7463_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
% ---------------- 绘图 8: 7498m漏层 ECD (g/cm³) 与安全窗口 ----------------
figure('Name', '7498m漏层ECD与安全密度窗口', 'Color', 'w');
plot(time_minutes_real, ECD_g_cm3(:, idx_7463), 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_ECD_7463_lower_line, 'g--', 'LineWidth', 2); 
plot(time_minutes_real, safe_ECD_7463_upper_line, 'r--', 'LineWidth', 2); 
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('7498m漏层ECD (g/cm³)', 'FontSize', 12, 'FontWeight', 'bold');
title('7498m漏层 ECD 随时间变化与安全密度窗口', 'FontSize', 14, 'FontWeight', 'bold');
legend('环空 ECD', sprintf('安全下限 (ECD=%.3f)', safe_ECD_7463_lower), ...
       sprintf('安全上限 (ECD=%.3f)', safe_ECD_7463_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

pressure_pump_surface_original_MPa = pressure_pump_surface(:);

% ---------------- 绘图 9: 泵压曲线 (MPa) ----------------
figure('Name', '泵压曲线', 'Color', 'w');
plot(time_minutes_real(:), pressure_pump_surface_original_MPa, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8); hold on;
[pump_pressure_max, pump_pressure_max_idx] = max(pressure_pump_surface_original_MPa);
plot(time_minutes_real(pump_pressure_max_idx), pump_pressure_max, 'o', ...
     'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('泵压 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('泵压随时间变化曲线', 'FontSize', 14, 'FontWeight', 'bold');
legend('泵压', sprintf('最大泵压 %.2f MPa', pump_pressure_max), 'Location', 'best');
grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.0, 'GridAlpha', 0.25);

%% 结果矩阵模块（含累计注入量 vs 压力）
Cumulative_Volume_m3 = volume_injected_all_list' / 1000;
Result_Matrix_Volume_Pressure = [time_minutes_real, Cumulative_Volume_m3, ECD_casing_g_cm3(:, end), pressure_pump_surface', pressure_annuli_MPa(:, end)];

%% === 8.11固井施工过程模拟：泵压对比模块 ===
design_pump_pressure_file = 'HT1-004_8_11_pump_pressure_reference.csv';
if exist(design_pump_pressure_file, 'file') ~= 2
    error('缺少8.11泵压参考表：%s', design_pump_pressure_file);
end
design_pump_pressure_811 = readtable(design_pump_pressure_file);
design_volume_811 = design_pump_pressure_811.volume_m3;
design_pressure_811 = design_pump_pressure_811.design_pressure_MPa;
[design_volume_unique, design_volume_first_idx] = unique(design_volume_811, 'stable');
design_pressure_unique = design_pressure_811(design_volume_first_idx);
volume_compare_811 = min(max(Cumulative_Volume_m3(:), design_volume_unique(1)), design_volume_unique(end));
design_pressure_interp_MPa = interp1(design_volume_unique, design_pressure_unique, volume_compare_811, 'linear');
pump_pressure_error_MPa = pressure_pump_surface_original_MPa(:) - design_pressure_interp_MPa(:);
Result_Matrix_PumpPressure_Comparison = [time_minutes_real(:), Cumulative_Volume_m3(:), ...
    pressure_pump_surface_original_MPa(:), design_pressure_interp_MPa(:), pump_pressure_error_MPa(:)];

fprintf('\n=== 8.11泵压对比（原始计算泵压 vs 设计表压力） ===\n');
fprintf('  8.11表压力范围(含停泵点): %.3f - %.3f MPa\n', min(design_pressure_811), max(design_pressure_811));
fprintf('  8.11累计量插值压力范围: %.3f - %.3f MPa\n', min(design_pressure_unique), max(design_pressure_unique));
fprintf('  计算泵压范围: %.3f - %.3f MPa\n', min(pressure_pump_surface_original_MPa), max(pressure_pump_surface_original_MPa));
valid_pump_pressure_error = pump_pressure_error_MPa(~isnan(pump_pressure_error_MPa));
fprintf('  泵压误差(计算-设计): 平均 %.3f MPa，最大绝对误差 %.3f MPa\n', ...
        mean(valid_pump_pressure_error), max(abs(valid_pump_pressure_error)));

figure('Name', '泵压对比-8.11施工过程模拟', 'Color', 'w');
plot(design_volume_unique, design_pressure_unique, 'k--', 'LineWidth', 1.8); hold on;
plot(Cumulative_Volume_m3(:), pressure_pump_surface_original_MPa(:), ...
     'Color', [0.00 0.45 0.74], 'LineWidth', 1.5);
xlabel('累计注入量 (m^3)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('泵压 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');
title('计算泵压与8.11固井施工过程模拟压力对比', 'FontSize', 14, 'FontWeight', 'bold');
legend('8.11表压力', '本脚本计算泵压', 'Location', 'best');
grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.0, 'GridAlpha', 0.25);

%% === 8.11固井施工过程模拟：井底ECD对比模块 ===
design_bottom_ecd_file = 'HT1-004_8_11_bottom_ecd_reference.csv';
if exist(design_bottom_ecd_file, 'file') ~= 2
    error('缺少8.11井底ECD参考表：%s', design_bottom_ecd_file);
end
design_bottom_ecd_811 = readtable(design_bottom_ecd_file);
design_bottom_ecd_volume_811 = design_bottom_ecd_811.volume_m3;
design_bottom_ecd_811_values = design_bottom_ecd_811.design_bottom_ECD_g_cm3;
[design_bottom_ecd_volume_unique, design_bottom_ecd_first_idx] = unique(design_bottom_ecd_volume_811, 'stable');
design_bottom_ecd_unique = design_bottom_ecd_811_values(design_bottom_ecd_first_idx);
volume_compare_bottom_ecd_811 = min(max(Cumulative_Volume_m3(:), design_bottom_ecd_volume_unique(1)), design_bottom_ecd_volume_unique(end));
design_bottom_ecd_interp = interp1(design_bottom_ecd_volume_unique, design_bottom_ecd_unique, volume_compare_bottom_ecd_811, 'linear');
calculated_bottom_ECD_g_cm3 = ECD_g_cm3(:, end);
bottom_ecd_error_g_cm3 = calculated_bottom_ECD_g_cm3(:) - design_bottom_ecd_interp(:);
Result_Matrix_BottomECD_Comparison = [time_minutes_real(:), Cumulative_Volume_m3(:), ...
    calculated_bottom_ECD_g_cm3(:), design_bottom_ecd_interp(:), bottom_ecd_error_g_cm3(:)];

fprintf('\n=== 8.11井底ECD对比（原始计算ECD vs 设计表井底ECD） ===\n');
fprintf('  8.11表井底ECD范围: %.3f - %.3f g/cm³\n', min(design_bottom_ecd_811_values), max(design_bottom_ecd_811_values));
fprintf('  8.11累计量插值井底ECD范围: %.3f - %.3f g/cm³\n', min(design_bottom_ecd_unique), max(design_bottom_ecd_unique));
fprintf('  计算井底ECD范围: %.3f - %.3f g/cm³\n', min(calculated_bottom_ECD_g_cm3), max(calculated_bottom_ECD_g_cm3));
valid_bottom_ecd_error = bottom_ecd_error_g_cm3(~isnan(bottom_ecd_error_g_cm3));
fprintf('  井底ECD误差(计算-设计): 平均 %.4f g/cm³，最大绝对误差 %.4f g/cm³\n', ...
        mean(valid_bottom_ecd_error), max(abs(valid_bottom_ecd_error)));

figure('Name', '井底ECD对比-8.11施工过程模拟', 'Color', 'w');
plot(design_bottom_ecd_volume_unique, design_bottom_ecd_unique, 'k--', 'LineWidth', 1.8); hold on;
plot(Cumulative_Volume_m3(:), calculated_bottom_ECD_g_cm3(:), ...
     'Color', [0.00 0.45 0.74], 'LineWidth', 1.5);
xlabel('累计注入量 (m^3)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('井底 ECD (g/cm³)', 'FontSize', 12, 'FontWeight', 'bold');
title('计算井底ECD与8.11固井施工过程模拟井底ECD对比', 'FontSize', 14, 'FontWeight', 'bold');
legend('8.11表井底ECD', '本脚本计算井底ECD', 'Location', 'best');
grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.0, 'GridAlpha', 0.25);

%% === 控压固井：四点联合约束(5578m/6600m/7498m/井底)的ECD控压反算模块 ===
%  约束采用各层位独立安全密度窗口:
%   5578m(套管鞋)：  ECD ∈ [1.940, 1.970]
%   6600m(关注点)：  ECD ∈ [1.940, 1.975]
%   7498m(漏层)：    ECD ∈ [1.940, 1.975]
%   井底TD：         ECD ∈ [2.010, 2.050]
%  逻辑：对每个时间步，计算满足四点约束的回压上下界，取可行中点；若下界>上界则报冲突
fprintf("开始计算基于四点独立窗口约束的所需回压...\n");
fprintf("  5578m窗口: [%.3f, %.3f]；6600m窗口: [%.3f, %.3f]；7498m窗口: [%.3f, %.3f]；井底窗口: [%.3f, %.3f]\n", ...
    safe_ECD_5568_lower, safe_ECD_5568_upper, safe_ECD_6880_lower, safe_ECD_6880_upper, ...
    safe_ECD_7463_lower, safe_ECD_7463_upper, safe_ECD_bottom_lower, safe_ECD_bottom_upper);

required_backpressure_MPa = zeros(n_time, 1);
new_bottom_pressure_MPa = zeros(n_time, 1);
BP_lower_history = zeros(n_time, 1);
BP_upper_history = zeros(n_time, 1);
conflict_flag = zeros(n_time, 1);

% ---- 各点安全压力边界 (MPa) — 不随时间变化的常数 ----
P_safe_l_5568 = safe_ECD_5568_lower * 0.00981 * TVD_cum(idx_5568);
P_safe_u_5568 = safe_ECD_5568_upper * 0.00981 * TVD_cum(idx_5568);
P_safe_l_6880 = safe_ECD_6880_lower * 0.00981 * TVD_cum(idx_6880);
P_safe_u_6880 = safe_ECD_6880_upper * 0.00981 * TVD_cum(idx_6880);
P_safe_l_7463 = safe_ECD_7463_lower * 0.00981 * TVD_cum(idx_7463);
P_safe_u_7463 = safe_ECD_7463_upper * 0.00981 * TVD_cum(idx_7463);
P_safe_l_bottom = safe_ECD_bottom_lower * 0.00981 * TVD_cum(end);
P_safe_u_bottom = safe_ECD_bottom_upper * 0.00981 * TVD_cum(end);

for t = 1:n_time
    % ---- 当前时间步各关键点的基础压力（静液+摩阻，不含回压） ----
    P_base_5568 = pressure_annuli_static_MPa(t, idx_5568) + pressure_annuli_friction_MPa(t, idx_5568);
    P_base_6880 = pressure_annuli_static_MPa(t, idx_6880) + pressure_annuli_friction_MPa(t, idx_6880);
    P_base_7463 = pressure_annuli_static_MPa(t, idx_7463) + pressure_annuli_friction_MPa(t, idx_7463);
    P_base_bottom = pressure_annuli_static_MPa(t, end) + pressure_annuli_friction_MPa(t, end);
    
    % ---- 回压下界（满足各点ECD≥对应层位安全下限）----
    BP_6880_lower = P_safe_l_6880 - P_base_6880;  % 6600m关注点压稳约束
    BP_7463_lower = P_safe_l_7463 - P_base_7463;
    BP_5568_lower = P_safe_l_5568 - P_base_5568;
    BP_bottom_lower = P_safe_l_bottom - P_base_bottom;
    BP_lower = max([BP_6880_lower, BP_7463_lower, BP_5568_lower, BP_bottom_lower, 0]);
    
    % ---- 回压上界（满足各点ECD≤对应层位安全上限）----
    BP_7463_upper = P_safe_u_7463 - P_base_7463;  % 7498m防漏刚性约束
    BP_6880_upper = P_safe_u_6880 - P_base_6880;
    BP_5568_upper = P_safe_u_5568 - P_base_5568;
    BP_bottom_upper = P_safe_u_bottom - P_base_bottom;
    BP_upper = min([BP_7463_upper, BP_6880_upper, BP_5568_upper, BP_bottom_upper]);
    
    BP_lower_history(t) = BP_lower;
    BP_upper_history(t) = BP_upper;
    
    if BP_lower > BP_upper + 1e-6
        conflict_flag(t) = 1;
        required_backpressure_MPa(t) = BP_upper;  % 防漏优先，取上界
    else
        % 在四点综合可行回压走廊内取中点，使采用回压位于上下界之间的中部。
        required_backpressure_MPa(t) = (BP_lower + BP_upper) / 2;
    end
    
    new_bottom_pressure_MPa(t) = P_base_bottom + required_backpressure_MPa(t);
end

% ---- 冲突统计与报警 ----
n_conflict = sum(conflict_flag);
if n_conflict > 0
    fprintf('⚠ 警告：%d/%d 个时间步四点独立窗口约束矛盾（综合回压下界 > 综合回压上界），自动取上界优先防漏！\n', n_conflict, n_time);
    fprintf('  冲突时间段：%.1f - %.1f min\n', ...
            time_minutes_real(find(conflict_flag, 1)), ...
            time_minutes_real(find(conflict_flag, 1, 'last')));
end

%% === 四点约束控压反算结果汇总 ===
fprintf('\n=== 四点联合约束控压反算完成(各点独立窗口，目标余量=%.3f) ===\n', safety_margin);
fprintf('  5578m窗口: [%.3f, %.3f] g/cm³\n', safe_ECD_5568_lower, safe_ECD_5568_upper);
fprintf('  6600m窗口: [%.3f, %.3f] g/cm³\n', safe_ECD_6880_lower, safe_ECD_6880_upper);
fprintf('  7498m窗口: [%.3f, %.3f] g/cm³\n', safe_ECD_7463_lower, safe_ECD_7463_upper);
fprintf('  井底TD窗口: [%.3f, %.3f] g/cm³\n', safe_ECD_bottom_lower, safe_ECD_bottom_upper);
fprintf('  所需动态回压范围: %.3f - %.3f MPa\n', min(required_backpressure_MPa), max(required_backpressure_MPa));
fprintf('  综合回压下界范围: %.3f - %.3f MPa\n', min(BP_lower_history), max(BP_lower_history));
fprintf('  综合回压上界范围: %.3f - %.3f MPa\n', min(BP_upper_history), max(BP_upper_history));
fprintf('  对应新井底压力范围: %.3f - %.3f MPa\n', min(new_bottom_pressure_MPa), max(new_bottom_pressure_MPa));

%% === 控压后四点压力和ECD计算 ===
pressure_5568_new_MPa = pressure_annuli_static_MPa(:, idx_5568) + ...
                        pressure_annuli_friction_MPa(:, idx_5568) + required_backpressure_MPa;
pressure_6880_new_MPa = pressure_annuli_static_MPa(:, idx_6880) + ...
                        pressure_annuli_friction_MPa(:, idx_6880) + required_backpressure_MPa;
pressure_7463_new_MPa = pressure_annuli_static_MPa(:, idx_7463) + ...
                        pressure_annuli_friction_MPa(:, idx_7463) + required_backpressure_MPa;

ECD_5568_new = pressure_5568_new_MPa ./ (0.00981 * TVD_cum(idx_5568));
ECD_6880_new = pressure_6880_new_MPa ./ (0.00981 * TVD_cum(idx_6880));
ECD_7463_new = pressure_7463_new_MPa ./ (0.00981 * TVD_cum(idx_7463));
ECD_bottom_new = new_bottom_pressure_MPa ./ (0.00981 * TVD_cum(end));

%% === 施加四点约束动态回压后的泵压计算 ===
backpressure_base_MPa = backpressure_time_list(:);
pump_backpressure_delta_MPa = required_backpressure_MPa(:) - backpressure_base_MPa;
pressure_pump_surface_with_backpressure_MPa = pressure_pump_surface_original_MPa + pump_backpressure_delta_MPa;
pressure_pump_surface_with_backpressure_MPa = max(pressure_pump_surface_with_backpressure_MPa, 0);
[pump_pressure_with_bp_max, pump_pressure_with_bp_max_idx] = max(pressure_pump_surface_with_backpressure_MPa);
[pump_pressure_with_bp_min, pump_pressure_with_bp_min_idx] = min(pressure_pump_surface_with_backpressure_MPa);

Result_Matrix_Volume_Pressure_WithBP = [time_minutes_real(:), ...
    Cumulative_Volume_m3(:), ...
    ECD_casing_g_cm3(:, end), ...
    pressure_pump_surface_original_MPa, ...
    pressure_pump_surface_with_backpressure_MPa, ...
    required_backpressure_MPa(:), ...
    pressure_annuli_MPa(:, end), ...
    new_bottom_pressure_MPa(:)];

fprintf('  原始泵压范围: %.3f - %.3f MPa\n', min(pressure_pump_surface_original_MPa), max(pressure_pump_surface_original_MPa));
fprintf('  施加回压后泵压范围: %.3f - %.3f MPa\n', min(pressure_pump_surface_with_backpressure_MPa), max(pressure_pump_surface_with_backpressure_MPa));
fprintf('  相对原始工况新增回压范围: %.3f - %.3f MPa\n', min(pump_backpressure_delta_MPa), max(pump_backpressure_delta_MPa));

%% === 绘图: 施加回压后的泵压曲线 ===
figure('Name', '施加回压后的泵压曲线', 'Color', 'w');
yyaxis left;
plot(time_minutes_real(:), pressure_pump_surface_original_MPa, '-', ...
     'Color', [0.55 0.55 0.55], 'LineWidth', 1.2); hold on;
plot(time_minutes_real(:), pressure_pump_surface_with_backpressure_MPa, '-', ...
     'Color', [0.00 0.45 0.74], 'LineWidth', 1.8);
plot(time_minutes_real(pump_pressure_with_bp_max_idx), pump_pressure_with_bp_max, 'o', ...
     'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
plot(time_minutes_real(pump_pressure_with_bp_min_idx), pump_pressure_with_bp_min, 's', ...
     'MarkerSize', 5, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerEdgeColor', 'k');
ylabel('泵压 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');

yyaxis right;
plot(time_minutes_real(:), required_backpressure_MPa(:), '--', ...
     'Color', [0.49 0.18 0.56], 'LineWidth', 1.3);
ylabel('动态回压 (MPa)', 'FontSize', 12, 'FontWeight', 'bold');

xlabel('时间 (min)', 'FontSize', 12, 'FontWeight', 'bold');
title('施加四点约束动态回压后的泵压曲线', 'FontSize', 14, 'FontWeight', 'bold');
legend('原始泵压', '施加回压后泵压', ...
       sprintf('最大 %.2f MPa', pump_pressure_with_bp_max), ...
       sprintf('最小 %.2f MPa', pump_pressure_with_bp_min), ...
       '动态回压', 'Location', 'best');
grid on; set(gca, 'FontSize', 11, 'LineWidth', 1.0, 'GridAlpha', 0.25);

%% === 绘图 9: 四点约束控压回压历史 ===
figure('Name', '四点联合约束控压反算-回压与约束边界', 'Color', 'w');
subplot(2,1,1);
plot(time_minutes_real, required_backpressure_MPa, 'm-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, BP_lower_history, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, BP_upper_history, 'r--', 'LineWidth', 1.2);
if n_conflict > 0
    plot(time_minutes_real(logical(conflict_flag)), required_backpressure_MPa(logical(conflict_flag)), ...
         'ko', 'MarkerSize', 4, 'DisplayName', '约束矛盾点');
end
xlabel('时间 (min)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('井口回压 (MPa)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('四点独立窗口约束下动态回压与上下界(内部含%.3f余量)', safety_margin), 'FontSize', 12, 'FontWeight', 'bold');
legend('采用回压', '综合回压下界', '综合回压上界', 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

subplot(2,1,2);
plot(time_minutes_real, new_bottom_pressure_MPa, 'b-', 'LineWidth', 1.5); hold on;
plot(time_minutes_real, safe_P_bottom_lower_line, 'g--', 'LineWidth', 1.5);
plot(time_minutes_real, safe_P_bottom_upper_line, 'r--', 'LineWidth', 1.5);
xlabel('时间 (min)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('井底总压力 (MPa)', 'FontSize', 11, 'FontWeight', 'bold');
title('施加控压回压后的井底压力 vs 安全密度窗口', 'FontSize', 12, 'FontWeight', 'bold');
legend('井底压力', sprintf('安全下限 (ECD=%.3f)', safe_ECD_bottom_lower), sprintf('安全上限 (ECD=%.3f)', safe_ECD_bottom_upper), 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

%% === 绘图 10: 四点ECD时间序列对比 ===
figure('Name', '四点联合约束-四关注点ECD对比', 'Color', 'w');

subplot(4,1,1);
plot(time_minutes_real, ECD_5568_new, 'b-', 'LineWidth', 1.3); hold on;
plot(time_minutes_real, safe_ECD_5568_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_5568_upper_line, 'r--', 'LineWidth', 1.2);
ylabel('5578m', 'FontSize', 9); title('5578m套管鞋 - 控压后ECD', 'FontSize', 10, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 9);

subplot(4,1,2);
plot(time_minutes_real, ECD_6880_new, 'm-', 'LineWidth', 1.3); hold on;
plot(time_minutes_real, safe_ECD_6880_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_6880_upper_line, 'r--', 'LineWidth', 1.2);
ylabel('6600m', 'FontSize', 9); title(sprintf('6600m关注点(窗口%.3f-%.3f) - 控压后ECD', safe_ECD_6880_lower, safe_ECD_6880_upper), 'FontSize', 10, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 9);

subplot(4,1,3);
plot(time_minutes_real, ECD_7463_new, 'r-', 'LineWidth', 1.3); hold on;
plot(time_minutes_real, safe_ECD_7463_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_7463_upper_line, 'r--', 'LineWidth', 1.2);
ylabel('7498m', 'FontSize', 9); title(sprintf('7498m漏层(窗口%.3f-%.3f) - 控压后ECD', safe_ECD_7463_lower, safe_ECD_7463_upper), 'FontSize', 10, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 9);

subplot(4,1,4);
plot(time_minutes_real, ECD_bottom_new, 'k-', 'LineWidth', 1.3); hold on;
plot(time_minutes_real, safe_ECD_bottom_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_bottom_upper_line, 'r--', 'LineWidth', 1.2);
xlabel('时间 (min)', 'FontSize', 10); ylabel('井底', 'FontSize', 9);
title('井底TD - 控压后ECD', 'FontSize', 10, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 9);

%% === 绘图 11: 6600m与7498m+5578m压力/ECD面板对比 ===
figure('Name', '压稳+防漏+套管鞋-压力与ECD对比', 'Color', 'w');

subplot(2,3,1);
plot(time_minutes_real, pressure_5568_new_MPa, 'b-', 'LineWidth', 1.2); hold on;
yline(safe_ECD_5568_lower * 0.00981 * TVD_cum(idx_5568), 'g--', 'LineWidth', 1.2);
yline(safe_ECD_5568_upper * 0.00981 * TVD_cum(idx_5568), 'r--', 'LineWidth', 1.2);
xlabel('time(min)'); ylabel('P(MPa)');
title('5578m套管鞋-控压后压力'); grid on; set(gca, 'FontSize', 9);

subplot(2,3,2);
plot(time_minutes_real, ECD_5568_new, 'b-', 'LineWidth', 1.2); hold on;
plot(time_minutes_real, safe_ECD_5568_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_5568_upper_line, 'r--', 'LineWidth', 1.2);
xlabel('time(min)'); ylabel('ECD');
title('5578m套管鞋-控压后ECD'); grid on; set(gca, 'FontSize', 9);

subplot(2,3,3);
plot(time_minutes_real, pressure_6880_new_MPa, 'm-', 'LineWidth', 1.2); hold on;
yline(safe_ECD_6880_lower * 0.00981 * TVD_cum(idx_6880), 'g--', 'LineWidth', 1.2);
yline(safe_ECD_6880_upper * 0.00981 * TVD_cum(idx_6880), 'r--', 'LineWidth', 1.2);
xlabel('time(min)'); ylabel('P(MPa)');
title(sprintf('6600m窗口%.3f-%.3f-控压后压力', safe_ECD_6880_lower, safe_ECD_6880_upper)); grid on; set(gca, 'FontSize', 9);

subplot(2,3,4);
plot(time_minutes_real, ECD_6880_new, 'm-', 'LineWidth', 1.2); hold on;
plot(time_minutes_real, safe_ECD_6880_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_6880_upper_line, 'r--', 'LineWidth', 1.2);
xlabel('time(min)'); ylabel('ECD');
title(sprintf('6600m窗口%.3f-%.3f-控压后ECD', safe_ECD_6880_lower, safe_ECD_6880_upper)); grid on; set(gca, 'FontSize', 9);

subplot(2,3,5);
plot(time_minutes_real, pressure_7463_new_MPa, 'r-', 'LineWidth', 1.2); hold on;
yline(safe_ECD_7463_lower * 0.00981 * TVD_cum(idx_7463), 'g--', 'LineWidth', 1.2);
yline(safe_ECD_7463_upper * 0.00981 * TVD_cum(idx_7463), 'r--', 'LineWidth', 1.2);
xlabel('time(min)'); ylabel('P(MPa)');
title(sprintf('7498m窗口%.3f-%.3f-控压后压力', safe_ECD_7463_lower, safe_ECD_7463_upper)); grid on; set(gca, 'FontSize', 9);

subplot(2,3,6);
plot(time_minutes_real, ECD_7463_new, 'r-', 'LineWidth', 1.2); hold on;
plot(time_minutes_real, safe_ECD_7463_lower_line, 'g--', 'LineWidth', 1.2);
plot(time_minutes_real, safe_ECD_7463_upper_line, 'r--', 'LineWidth', 1.2);
xlabel('time(min)'); ylabel('ECD');
title(sprintf('7498m窗口%.3f-%.3f-控压后ECD', safe_ECD_7463_lower, safe_ECD_7463_upper)); grid on; set(gca, 'FontSize', 9);

%% === 四点联合约束控压反算最终汇总 ===
fprintf('\n=== 四点联合约束控压反算最终汇总（各点独立窗口，余量=%.3f） ===\n', safety_margin);
fprintf('  5578m(套管鞋) 控压后: ECD ∈ [%.3f, %.3f] g/cm³，窗口[%.3f, %.3f]\n', min(ECD_5568_new), max(ECD_5568_new), safe_ECD_5568_lower, safe_ECD_5568_upper);
fprintf('  6600m(关注点) 控压后: ECD ∈ [%.3f, %.3f] g/cm³，窗口[%.3f, %.3f]\n', min(ECD_6880_new), max(ECD_6880_new), safe_ECD_6880_lower, safe_ECD_6880_upper);
fprintf('  7498m(漏层) 控压后: ECD ∈ [%.3f, %.3f] g/cm³，窗口[%.3f, %.3f]\n', min(ECD_7463_new), max(ECD_7463_new), safe_ECD_7463_lower, safe_ECD_7463_upper);
fprintf('  井底TD 控压后: ECD ∈ [%.3f, %.3f] g/cm³，窗口[%.3f, %.3f]\n', min(ECD_bottom_new), max(ECD_bottom_new), safe_ECD_bottom_lower, safe_ECD_bottom_upper);
if n_conflict > 0
    fprintf('  ⚠ 存在%d个矛盾时间点(%.0f%%)！各点独立窗口不可同时满足，需人工调整参数(密度/排量/流变/回压策略)。\n', ...
            n_conflict, n_conflict/n_time*100);
else
    fprintf('  ✔ 所有时间步四点独立窗口均满足(含%.3f余量)。\n', safety_margin);
end

%% === 井底安全密度窗口统一[2.010, 2.050]回压施加值反算模块（跌破段持稳+交点后冻结回压保持原形状） ===
%  功能：基于统一后的井底安全密度窗口 [safe_ECD_bottom_lower, safe_ECD_bottom_upper]，
%        逐时间步反算井口回压施加值 BP_apply(t)，复现用户手绘红线的抬升趋势：
%   1) 自然ECD跌破窗口下限的连续段（含停泵谷）：逐时补亏量，将井底ECD持稳在目标持稳密度
%      ECD_setpoint_win（默认窗口下限+0.015≈2.025），跌破段被填平持稳；
%   2) 自然ECD回升与窗口下限相交的时刻（交点）之后：回压冻结为交点时刻的值不再撤回，
%      控压曲线以"交点前加回压后的目标压力值"为起点、保持自然曲线形状平行后移；
%   3) 间隔≤5 min的多次跌破段合并为一段（滤除交点附近的短促回穿扰动，避免控压曲线小突起），
%      交点取合并段的最后回升时刻；
%   4) 尚未进入任何跌破段、且自然ECD在窗口内时：回压为0，保持自然曲线原样；
%   5) 自然ECD > 窗口上限：回压只能抬升不能压低，BP_apply = 0，并统计超压预警；
%   6) 同时给出每时刻回压可行走廊 [BP_min, BP_max]，供现场节流阀调控参考。

% ---- 自然（不含回压）井底压力与ECD ----
P_base_bottom_win  = pressure_annuli_static_MPa(:, end) + pressure_annuli_friction_MPa(:, end);  % MPa
ECD_base_bottom_win = P_base_bottom_win ./ (0.00981 * TVD_cum(end));                              % g/cm³

% ---- 井底窗口压力边界与目标持稳密度 ----
P_safe_l_win   = safe_ECD_bottom_lower * 0.00981 * TVD_cum(end);   % MPa，窗口下限（孔隙压力当量）
P_safe_u_win   = safe_ECD_bottom_upper * 0.00981 * TVD_cum(end);   % MPa，窗口上限（破裂压力当量）
ECD_setpoint_win = safe_ECD_bottom_lower + 0.015;                 % g/cm³，目标持稳密度（默认2.025，可调）

% ---- 回压可行走廊 ----
BP_min_win = max(P_safe_l_win - P_base_bottom_win, 0);  % 保证井底不低于窗口下限的最小回压
BP_max_win = P_safe_u_win - P_base_bottom_win;          % 保证井底不高于窗口上限的最大回压
overpress_flag_win = BP_max_win < 0;                    % 自然ECD已超上限：回压无法压低，仅预警

% ---- 跌破段识别与合并（间隔≤5 min合并为一段，消除短促回穿造成的小突起） ----
below_raw_win = ECD_base_bottom_win < safe_ECD_bottom_lower;
d0_win = diff([0; below_raw_win(:); 0]);
st0_win = find(d0_win == 1);
en0_win = find(d0_win == -1) - 1;
seg_s_win = []; seg_e_win = [];
for k = 1:numel(st0_win)
    if ~isempty(seg_e_win) && (time_minutes_real(st0_win(k)) - time_minutes_real(seg_e_win(end))) <= 5
        seg_e_win(end) = en0_win(k);   % 间隔≤5 min：并入前一段
    else
        seg_s_win(end+1) = st0_win(k); seg_e_win(end+1) = en0_win(k);
    end
end
below_win = false(n_time, 1);
for k = 1:numel(seg_s_win)
    below_win(seg_s_win(k):seg_e_win(k)) = true;
end
seg_count_win = numel(seg_s_win);

% ---- 回压施加值反算：跌破段逐时补亏量持稳；交点后冻结回压、曲线保持原形状 ----
BP_ecd_win = zeros(n_time, 1);                                   % 回压的井底ECD当量
freeze_ecd_win = 0;                                              % 冻结回压（ECD当量），交点时刻更新后保持
for t = 1:n_time
    if overpress_flag_win(t)
        BP_ecd_win(t) = 0;                                       % 自然超上限：回压无法压低
    elseif below_win(t)
        BP_ecd_win(t) = ECD_setpoint_win - ECD_base_bottom_win(t);  % 逐时补亏量持稳至目标密度
        freeze_ecd_win = BP_ecd_win(t);                          % 实时记录，段末(交点)值即为冻结值
    else
        BP_ecd_win(t) = freeze_ecd_win;                          % 交点之后：冻结回压，曲线平行保持原形状
    end
end
BP_apply_win = BP_ecd_win * 0.00981 * TVD_cum(end);              % MPa，井口回压施加值

% ---- 控压后环空关键点压力/ECD（井口回压沿环空不变传递至各深度） ----
P_bottom_ctrl_win = P_base_bottom_win + BP_apply_win;
ECD_bottom_ctrl_win = P_bottom_ctrl_win ./ (0.00981 * TVD_cum(end));
P_shoe_ctrl_win = pressure_annuli_static_MPa(:, idx_5568) + pressure_annuli_friction_MPa(:, idx_5568) + BP_apply_win;
ECD_shoe_ctrl_win = P_shoe_ctrl_win ./ (0.00981 * TVD_cum(idx_5568));
P_wh_ctrl_win = pressure_annuli_static_MPa(:, 1) + pressure_annuli_friction_MPa(:, 1) + BP_apply_win;  % 控压后井口压力≈施加回压

% ---- 汇总输出 ----
fprintf('\n=== 井底安全密度窗口统一回压施加值反算完成（跌破段持稳+交点后冻结） ===\n');
fprintf('  井底安全密度窗口: [%.3f, %.3f] g/cm³，目标持稳密度 %.3f g/cm³\n', ...
        safe_ECD_bottom_lower, safe_ECD_bottom_upper, ECD_setpoint_win);
for k = 1:seg_count_win
    fprintf('  跌破段%d: %.1f - %.1f min（交点 %.1f min），交点处冻结回压 %.3f MPa，交点后控压ECD = 自然ECD + %.4f\n', ...
            k, time_minutes_real(seg_s_win(k)), time_minutes_real(seg_e_win(k)), ...
            time_minutes_real(seg_e_win(k)), BP_apply_win(seg_e_win(k)), BP_ecd_win(seg_e_win(k)));
end
if seg_count_win == 0
    fprintf('  跌破段数: 0（自然井底ECD全程不低于窗口下限，回压为0保持原样）\n');
end
[bp_max_win, bp_max_idx_win] = max(BP_apply_win);
fprintf('  回压施加值范围: %.3f - %.3f MPa，最大回压 %.3f MPa 出现在 %.1f min\n', ...
        min(BP_apply_win), bp_max_win, bp_max_win, time_minutes_real(bp_max_idx_win));
fprintf('  回压可行走廊: BP_min ∈ [%.3f, %.3f] MPa，BP_max ∈ [%.3f, %.3f] MPa\n', ...
        min(BP_min_win), max(BP_min_win), min(BP_max_win), max(BP_max_win));
fprintf('  控压后井底ECD: [%.3f, %.3f] g/cm³（跌破段持稳于%.3f，交点后保持自然形状）\n', ...
        min(ECD_bottom_ctrl_win), max(ECD_bottom_ctrl_win), ECD_setpoint_win);
n_over_win = sum(overpress_flag_win);
if n_over_win > 0
    fprintf('  ⚠ 警告：%d/%d 个时间步自然井底ECD已超窗口上限，回压只能抬升不能压低，该段BP取0\n', n_over_win, n_time);
end
n_ctrl_over_win = sum(ECD_bottom_ctrl_win > safe_ECD_bottom_upper + 1e-6);
if n_ctrl_over_win > 0
    fprintf('  ⚠ 警告：%d/%d 个时间步控压后井底ECD超窗口上限（冻结回压叠加自然上升所致），建议交点后分段撤回回压\n', ...
            n_ctrl_over_win, n_time);
end
fprintf('  控压后5578m套管鞋ECD: [%.3f, %.3f] g/cm³（窗口[%.3f, %.3f]）\n', ...
        min(ECD_shoe_ctrl_win), max(ECD_shoe_ctrl_win), safe_ECD_5568_lower, safe_ECD_5568_upper);

% ---- 绘图：施加回压后的ECD抬升持稳效果（论文图样式） ----
figure('Name', '回压施加井底窗口抬升持稳效果', 'Color', 'w');
yyaxis left;
plot(time_minutes_real, ECD_base_bottom_win, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0); hold on;
plot(time_minutes_real, ECD_bottom_ctrl_win, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.5);
plot(time_minutes_real, ECD_shoe_ctrl_win,  '-', 'Color', [0.85 0.60 0.20], 'LineWidth', 1.2);
plot(time_minutes_real, ones(n_time, 1) * safe_ECD_bottom_lower, '--', 'Color', [0.20 0.63 0.64], 'LineWidth', 1.5);
plot(time_minutes_real, ones(n_time, 1) * safe_ECD_bottom_upper, '--', 'Color', [0.72 0.25 0.30], 'LineWidth', 1.5);
plot(time_minutes_real, ones(n_time, 1) * ECD_setpoint_win, ':', 'Color', [0.85 0.20 0.20], 'LineWidth', 1.2);
ylabel('ECD distribution of annular pressure(g/cm³)', 'FontSize', 12, 'FontWeight', 'bold');
yyaxis right;
plot(time_minutes_real, P_wh_ctrl_win, '-', 'Color', [0.55 0.35 0.75], 'LineWidth', 1.2);
ylabel('Wellhead back pressure(MPa)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time(min)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('统一井底窗口[%.2f, %.2f]回压施加值反算与控压后持稳效果(目标%.3f)', ...
      safe_ECD_bottom_lower, safe_ECD_bottom_upper, ECD_setpoint_win), 'FontSize', 13, 'FontWeight', 'bold');
legend('Bottomhole (no BP)', 'Bottomhole', 'Casing shoe', 'Pore pressure', ...
       'Break-down pressure', 'Setpoint', 'Wellhead', 'Location', 'best');
grid on; set(gca, 'FontSize', 11);
disp("井底窗口回压施加值反算模块执行完成！");
