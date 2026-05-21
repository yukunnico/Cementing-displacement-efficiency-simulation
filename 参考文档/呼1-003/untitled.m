clear;clc;   %zzzzzzzzzz2026.5.21 HT1-003更新
%% 初始化数据
% 请确保读取的 .csv 是包含变径点和实测数据的最终版表格
structure_data = readtable('呼1-003井身结构.csv');
c_depth = structure_data.depth_well_logging_m_;
% 节点数
n_segment = height(structure_data);  % 井身分段数量（表格行数）

%% 核心：基于HT1-003固井管柱结构的七段式深度划分
% 定义七个深度区间的阈值（六个分界点）
% 井底钻达深度 (TD): 7586m
depth_threshold1 = 3291.897;  % 第一变径点：149.2mm钻杆 变 127mm钻杆
depth_threshold2 = 5277.754;  % 第二变径点：127mm钻杆 变 168.3mm尾管(悬挂器位置)
depth_threshold3 = 5568.00;   % 第三变径点：273.1mm套管鞋 变 241.3mm裸眼 【关注点1：尾管悬挂器位置】
depth_threshold4 = 7059.016;  % 第四变径点：168.3mm尾管 变 139.7mm尾管
depth_threshold5 = 7096.00;   % 第五变径点：241.3mm裸眼 变 215.9mm裸眼
depth_threshold6 = 7586.00;   % 第六变径点：井底TD
length_casing = 7586; % 下入管柱总长度,m

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
% 段1: 149.2mm 钻杆 在 273.1mm 套管内 (0-3291.897m)
diameter_bit_out(idx1) = 245.42 * 0.001;    % 273.1套管内径(245.42mm)作为环空外径
diameter_casing_out(idx1) = 149.2 * 0.001;  % 149.2钻杆外径
diameter_casing_in(idx1) = (149.2 - 9.65 * 2) * 0.001; % 149.2钻杆内径(壁厚9.65mm)
% 段2: 127mm 钻杆 在 273.1mm 套管内 (3291.897-5277.754m)
diameter_bit_out(idx2) = 245.42 * 0.001;
diameter_casing_out(idx2) = 127.0 * 0.001;  % 127钻杆外径
diameter_casing_in(idx2) = (127.0 - 9.65 * 2) * 0.001; % 127钻杆内径(壁厚9.65mm)
% 段3: 168.3mm 尾管 在 273.1mm 套管内 (5277.754-5568m, 重叠段)
diameter_bit_out(idx3) = 245.42 * 0.001;
diameter_casing_out(idx3) = 168.3 * 0.001;  % 168.3尾管外径
diameter_casing_in(idx3) = (168.3 - 14.7 * 2) * 0.001; % 168.3尾管内径(壁厚14.7mm)
% 段4: 168.3mm 尾管 在 241.3mm 裸眼中 (5568-7059.016m) 【关注点1】
diameter_bit_out(idx4) = structure_data.annulus_radius_array_cm_(idx4) * 0.01; 
diameter_casing_out(idx4) = 168.3 * 0.001;
diameter_casing_in(idx4) = (168.3 - 14.7 * 2) * 0.001; % 壁厚14.7mm
% 段5: 139.7mm 尾管 在 241.3mm 裸眼中 (7059.016-7096m)
diameter_bit_out(idx5) = structure_data.annulus_radius_array_cm_(idx5) * 0.01;
diameter_casing_out(idx5) = 139.7 * 0.001;
diameter_casing_in(idx5) = (139.7 - 15.88 * 2) * 0.001; % 壁厚15.88mm
% 段6: 139.7mm 尾管 在 215.9mm 裸眼中 (7096-7586m)
diameter_bit_out(idx6) = structure_data.annulus_radius_array_cm_(idx6) * 0.01;
diameter_casing_out(idx6) = 139.7 * 0.001;
diameter_casing_in(idx6) = (139.7 - 15.88 * 2) * 0.001;
% 段7: 7586m 以下 (尾管鞋以下纯裸眼口袋)
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
volume_in_fenggu = sum(volume_in_annual_seg(172:end));  % 该段体积 = 本段环空面积 * 本段套管/钻杆长度
volume_in_drilling_casing = sum(volume_in_casing_seg);     % 套管内总体积 = 各段体积之和（m³）
volume_in_annual = sum(volume_in_annual_seg);     % 环空内总体积 = 各段体积之和（m³）
volume_in_drilling_casing_L = volume_in_drilling_casing * 1000;     % 转换为升（L）
volume_in_annual_L = volume_in_annual * 1000;     % 转换为升（L）

volume_in_allfluid =volume_in_drilling_casing + sum(volume_in_annual_seg(119:end));  % 该段体积 = 套管内总体积 +  隔离液返高至井底环空累加，，，这里目前需要手动改


%% ==================== 初始参数值调整模块 ====================
% 依据 HT1-003固井施工设计 真实注入阶段划分为泥浆变排量体系
% 0-钻井液(井浆) | 1-先导浆 | 2-驱油隔离液 | 3-领浆 | 4-中间浆 | 5-尾浆 | 6-压塞液 | 7-替浆钻井液 | 8-保护液 | 9-基液 | 91-替浆钻井液1 | 92-替浆钻井液2 | 93-替浆钻井液3 | 94-替浆钻井液4 | 95-替浆钻井液5

% 泵排量数据 (单位: L/min) - 来自HT1-003施工设计排量表
pump_rate1 = 1.2 * 1000;            % 1-先导浆 循环排量 1.4m³/min
pump_rate2 = 1.2 * 1000;            % 2-驱油隔离液 循环排量 1.4m³/min
pump_rate3 = 1.2 * 1000;            % 3-领浆 循环排量 1.4m³/min
pump_rate4 = 1.2 * 1000;            % 4-中间浆 循环排量 1.4m³/min
pump_rate5 = 1.2 * 1000;            % 5-尾浆 循环排量 1.4m³/min
pump_rate6 = 1.6 * 1000;            % 6-压塞液 循环排量 1.4m³/min
pump_rate7 = 1.6 * 1000;            % 7-替浆钻井液 循环排量 1.6m³/min
pump_rate8 = 1.4 * 1000;            % 8-保护液 循环排量 1.4m³/min
pump_rate9 = 1.4 * 1000;            % 9-基液 循环排量 1.2m³/min
pump_rate91 = 1.4 * 1000;           % 91-替浆钻井液1 循环排量 1.2m³/min
pump_rate92 = 1.2 * 1000;             % 92-替浆钻井液2 循环排量 0.7m³/min
pump_rate93 = 1 * 1000;          % 93-替浆钻井液3 循环排量 0.7m³/min
pump_rate94 = 0.8 * 1000;           % 94-替浆钻井液4 循环排量 0.7m³/min
pump_rate95 = 0.7 * 1000;           % 95-替浆钻井液5 循环排量 0.7m³/min
Pump_values = [pump_rate1, pump_rate2, pump_rate3, pump_rate4, pump_rate5, pump_rate6, pump_rate7, pump_rate8, pump_rate9, pump_rate91, pump_rate92, pump_rate93, pump_rate94, pump_rate95];

% 控压参数值
pressure_back = 0;                  % 控压钻进过程中的井口回压（MPa）
pressure_back_static = 0;           % 控压静止过程中的井口回压（MPa）

% 固井流体参数（密度、黏度、屈服应力）- 来自HT1-003施工设计
% 后四段替浆钻井液(91,92,93,94)的流体属性在此处与替浆钻井液(7)保持一致
rou0 = 1.95; rou1 = 1.75; rou2 = 2.00; rou3 = 2.05; rou4 = 1.95; rou5 = 1.95; rou6 = 1.97; rou7 = 1.95; rou8 = 1.97; rou9 = 1.02; rou91 = 1.95; rou92 = 1.95; rou93 = 1.95; rou94 = 1.95; rou95 = 1.95; % 密度（g/cm³）
miu0 = 51 ; miu1 = 55; miu2 = 30; miu3 = 50; miu4 = 180; miu5 = 180; miu6 = 30; miu7 = 30; miu8 = 30; miu9 = 30; miu91 = 30; miu92 = 30; miu93 = 30; miu94 = 30; miu95 = 30; % 黏度（mPa·s）沿用旧参数
tau0 = 10 ; tau1 = 9.2; tau2 = 9; tau3 = 11; tau4 = 14; tau5 = 14; tau6 = 8; tau7 = 8.5; tau8 = 8.2; tau9 = 8; tau91 = 8.3; tau92 = 8.3; tau93 = 8.3; tau94 = 8.3; tau95 = 8.3;  % 屈服应力（Pa）沿用旧参数

%% 各流体体积参数 (按框架结构，水泥浆保留深度分段求和模式，行号来自CSV)
tag_interface_top_lead   = 172; % 领浆顶界行号 (对应 5090m)
tag_interface_top_middle = 205; % 中间浆顶界行号 (对应 6000m)
tag_interface_top_tail   = 240; % 尾浆顶界行号 (对应 7000m)

v1 = 28 * 1000; % 1-先导浆体积,L (30m³)
v2 = 26 * 1000; % 2-隔离液体积,L (26m³)
v3 = sum(structure_data.volume_annulus_L_(tag_interface_top_lead : tag_interface_top_middle)); % 3-领浆体积,L
v4 = sum(structure_data.volume_annulus_L_(tag_interface_top_middle+1 : tag_interface_top_tail)); % 4-中间浆体积,L
v5 = sum(structure_data.volume_annulus_L_(tag_interface_top_tail+1 : end))+2 * 1000; % 5-尾浆体积,L，补上2m3
v6 = 2 * 1000;  % 6-压塞液体积,L (2m³)
v7 = 25 * 1000; % 7-替浆钻井液体积,L (25m³)
v8 = 14 * 1000; % 8-保护液体积,L (14m³)
v9 = 1 * 1000;  % 9-基液体积,L (1m³)
v91 = 10 * 1000; % 91-替浆钻井液1体积,L (15m³, 排量1.2)
v92 = 5 * 1000; % 92-替浆钻井液2体积,L (15m³, 排量0.7)
v93 = 15 * 1000; % 93-替浆钻井液3体积,L (10m³, 排量0.7)
v94 = 10 * 1000; % 93-替浆钻井液3体积,L (10m³, 排量0.7)
v95 = volume_in_drilling_casing_L - (v6 + v7 + v8 + v9 + v91 + v92 + v93 + v94); % 95-替浆钻井液4体积,L (确保碰压精确补齐套管容积)，-2* 1000是补的和现场预测的差距，相当于简化连接设备等的体积误差补充

% 计算总时间步数与分段数量
dt = 1;  % 【高精度优化】：大幅细化时间步长为 min
n_time = floor((v1./Pump_values(1) + v2./Pump_values(2) + v3./Pump_values(3) + v4./Pump_values(4) + v5./Pump_values(5) + v6./Pump_values(6) + v7./Pump_values(7) + v8./Pump_values(8) + v9./Pump_values(9) + v91./Pump_values(10) + v92./Pump_values(11) + v93./Pump_values(12) + v94./Pump_values(13) + v95./Pump_values(14)) / dt);
time = (v1./Pump_values(1) + v2./Pump_values(2) + v3./Pump_values(3) + v4./Pump_values(4) + v5./Pump_values(5) + v6./Pump_values(6) + v7./Pump_values(7) + v8./Pump_values(8) + v9./Pump_values(9) + v91./Pump_values(10) + v92./Pump_values(11) + v93./Pump_values(12) + v94./Pump_values(13) + v95./Pump_values(14)) / dt;

%% 回压数值设置
backpressure_1 = 0;           % 1-先导浆回压，MPa
backpressure_2 = 0;           % 2-驱油隔离液回压，MPa
backpressure_3 = 0;           % 3-领浆回压，MPa
backpressure_4 = 0;           % 4-中间浆回压，MPa
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
    elseif t > pump_time_node(13)
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
        
        % 第二界面,先导浆-隔离液
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
        
        % 第三界面，隔离液-领浆
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
        
        % 第四界面，领浆-中间浆
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

        % 第五界面，中间浆-尾浆
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
% 第一段：深度 <= depth_threshold1 (3488.94m)
segment1 = (c_depth <= depth_threshold1);      
if any(segment1)
    mean_diameter_segment1 = mean(out_diam_bole(segment1));
    mean_casing_segment1 = mean(diameter_casing_out(segment1))*1000;
    PR_segment1 = PR_Friction(mean_diameter_segment1, mean_casing_segment1);
else
    PR_segment1 = 0;
end
% 第二段：深度 > depth_threshold1 & <= depth_threshold2 (3742.22m)
segment2 = (c_depth > depth_threshold1) & (c_depth <= depth_threshold2);      
if any(segment2)
    mean_diameter_segment2 = mean(out_diam_bole(segment2));
    mean_casing_segment2 = mean(diameter_casing_out(segment2))*1000;
    PR_segment2 = PR_Friction(mean_diameter_segment2, mean_casing_segment2);
else
    PR_segment2 = 0;
end
% 第三段：深度 > depth_threshold2 & <= depth_threshold3 (5292.52m)
segment3 = (c_depth > depth_threshold2) & (c_depth <= depth_threshold3); 
if any(segment3)
    mean_diameter_segment3 = mean(out_diam_bole(segment3));
    mean_casing_segment3 = mean(diameter_casing_out(segment3))*1000; 
    PR_segment3 = PR_Friction(mean_diameter_segment3, mean_casing_segment3);
else
    PR_segment3 = 0;
end
% 第四段：深度 > depth_threshold3 & <= depth_threshold4 (5630m)
segment4 = (c_depth > depth_threshold3) & (c_depth <= depth_threshold4); 
if any(segment4)
    mean_diameter_segment4 = mean(out_diam_bole(segment4));
    mean_casing_segment4 = mean(diameter_casing_out(segment4))*1000;
    PR_segment4 = PR_Friction(mean_diameter_segment4, mean_casing_segment4);
else
    PR_segment4 = 0;
end
% 第五段：深度 > depth_threshold4 & <= depth_threshold5 (7052.26m)
segment5 = (c_depth > depth_threshold4) & (c_depth <= depth_threshold5);                          
if any(segment5)
    mean_diameter_segment5 = mean(out_diam_bole(segment5));
    mean_casing_segment5 = mean(diameter_casing_out(segment5))*1000;
    PR_segment5 = PR_Friction(mean_diameter_segment5, mean_casing_segment5);
else
    PR_segment5 = 0;
end
% 第六段：深度 > depth_threshold5 & <= depth_threshold6 (7554.00m)
segment6 = (c_depth > depth_threshold5) & (c_depth <= depth_threshold6);                          
if any(segment6)
    mean_diameter_segment6 = mean(out_diam_bole(segment6));
    mean_casing_segment6 = mean(diameter_casing_out(segment6))*1000;
    PR_segment6 = PR_Friction(mean_diameter_segment6, mean_casing_segment6);
else
    PR_segment6 = 0;
end
% 第七段：深度 > depth_threshold6 (7554m - 7559m 口袋)
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