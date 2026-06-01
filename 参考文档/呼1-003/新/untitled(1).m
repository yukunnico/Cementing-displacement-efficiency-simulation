clear;clc;   %zzzzzzzzzz2026.5.21 HT1-003更新
%% 初始化数据
% 请确保读取的 .csv 是包含变径点和实测数据的最终版表格
structure_data = readtable('呼1-003井身结构.csv');
c_depth = structure_data.depth_well_logging_m_;
% 节点数
n_segment = height(structure_data);  % 井身分段数量（表格行数）

%% 核心：基于HT1-003固井管柱结构的七段式深度划分
% 定义七个深度区间的阈值（六个分界点）
% 井底钻达深度 (TD): 7618m
depth_threshold1 = 3321.682;  % 第一变径点：149.2mm钻杆 变 127mm钻杆
depth_threshold2 = 5307.539;  % 第二变径点：127mm钻杆 变 168.3mm尾管(悬挂器位置)
depth_threshold3 = 5568.00;   % 第三变径点：273.1mm套管鞋 变 241.3mm裸眼 【关注点1：套管鞋】
depth_threshold4 = 7089.576;  % 第四变径点：168.3mm尾管 变 139.7mm尾管
depth_threshold5 = 7096.00;   % 第五变径点：241.3mm裸眼 变 215.9mm裸眼
depth_threshold6 = 7618.00;   % 第六变径点：井底TD
length_casing = 7618; % 下入管柱总长度,m

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
% 段1: 149.2mm 钻杆 在 273.1mm 套管内 (0-3321.682m)
diameter_bit_out(idx1) = 245.42 * 0.001;    % 273.1套管内径(245.42mm)作为环空外径
diameter_casing_out(idx1) = 149.2 * 0.001;  % 149.2钻杆外径
diameter_casing_in(idx1) = (149.2 - 9.65 * 2) * 0.001; % 149.2钻杆内径(壁厚9.65mm)
% 段2: 127mm 钻杆 在 273.1mm 套管内 (3321.682-5307.539m)
diameter_bit_out(idx2) = 245.42 * 0.001;
diameter_casing_out(idx2) = 127.0 * 0.001;  % 127钻杆外径
diameter_casing_in(idx2) = (127.0 - 9.65 * 2) * 0.001; % 127钻杆内径(壁厚9.65mm)
% 段3: 168.3mm 尾管 在 273.1mm 套管内 (5307.539-5568m, 重叠段)
diameter_bit_out(idx3) = 245.42 * 0.001;
diameter_casing_out(idx3) = 168.3 * 0.001;  % 168.3尾管外径
diameter_casing_in(idx3) = (168.3 - 14.7 * 2) * 0.001; % 168.3尾管内径(壁厚14.7mm)
% 段4: 168.3mm 尾管 在实测裸眼中 (5568-7089.576m) 【关注点1】
diameter_bit_out(idx4) = structure_data.annulus_radius_array_cm_(idx4) * 0.01; 
diameter_casing_out(idx4) = 168.3 * 0.001;
diameter_casing_in(idx4) = (168.3 - 14.7 * 2) * 0.001; % 壁厚14.7mm
% 段5: 139.7mm 尾管 在实测裸眼中 (7089.576-7096m)
diameter_bit_out(idx5) = structure_data.annulus_radius_array_cm_(idx5) * 0.01;
diameter_casing_out(idx5) = 139.7 * 0.001;
diameter_casing_in(idx5) = (139.7 - 15.88 * 2) * 0.001; % 壁厚15.88mm
% 段6: 139.7mm 尾管 在实测裸眼中 (7096-7618m)
diameter_bit_out(idx6) = structure_data.annulus_radius_array_cm_(idx6) * 0.01;
diameter_casing_out(idx6) = 139.7 * 0.001;
diameter_casing_in(idx6) = (139.7 - 15.88 * 2) * 0.001;
% 段7: 7618m 以下 (尾管鞋以下纯裸眼口袋)
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
idx_slurry_top = find(c_depth >= 5100, 1, 'first');
volume_in_fenggu = sum(volume_in_annual_seg(idx_liner_top:end));  % 5307.539m至井底封固段环空体积
volume_in_drilling_casing = sum(volume_in_casing_seg);     % 套管内总体积 = 各段体积之和（m³）
volume_in_annual = sum(volume_in_annual_seg);     % 环空内总体积 = 各段体积之和（m³）
volume_in_drilling_casing_L = volume_in_drilling_casing * 1000;     % 转换为升（L）
volume_in_annual_L = volume_in_annual * 1000;     % 转换为升（L）

volume_in_allfluid = volume_in_drilling_casing + sum(volume_in_annual_seg(idx_slurry_top:end));  % 套管内容积+5100m至井底环空体积


%% ==================== 初始参数值调整模块 ====================
% 依据 HT1-003固井施工设计 8.11 施工过程模拟划分为变排量体系
% 0-钻井液(井浆) | 1-先导浆 | 2-低失水驱油隔离液1 | 3-低失水驱油隔离液2 | 4-领浆 | 5-尾浆 | 6-压塞液 | 7-替浆钻井液 | 8-保护液 | 9-基液 | 91-替浆钻井液1 | 92-替浆钻井液2 | 93-替浆钻井液3 | 94-替浆钻井液4 | 95-保留段

 % 泵排量数据 (单位: L/min) - 来自8.11固井施工过程模拟
 pump_rate1 = 1.2 * 1000;            % 1-先导浆
 pump_rate2 = 1.2 * 1000;            % 2-低失水驱油隔离液1
 pump_rate3 = 1.2 * 1000;            % 3-低失水驱油隔离液2
 pump_rate4 = 1.2 * 1000;            % 4-领浆
 pump_rate5 = 1.4 * 1000;            % 5-尾浆
 pump_rate6 = 1.6 * 1000;            % 6-压塞液
 pump_rate7 = 1.6 * 1000;            % 7-替浆钻井液
 pump_rate8 = 1.4 * 1000;            % 8-保护液
 pump_rate9 = 1.4 * 1000;            % 9-基液
 pump_rate91 = 1.2 * 1000;           % 91-替浆钻井液1
 pump_rate92 = 1.0 * 1000;           % 92-替浆钻井液2
 pump_rate93 = 0.8 * 1000;           % 93-替浆钻井液3
 pump_rate94 = 0.7 * 1000;           % 94-替浆钻井液4
 pump_rate95 = 0.7 * 1000;           % 95-保留段(体积为0)
 Pump_values = [pump_rate1, pump_rate2, pump_rate3, pump_rate4, pump_rate5, pump_rate6, pump_rate7, pump_rate8, pump_rate9, pump_rate91, pump_rate92, pump_rate93, pump_rate94, pump_rate95];
    
% 控压参数值
pressure_back = 0;                  % 控压钻进过程中的井口回压（MPa）
pressure_back_static = 0;           % 控压静止过程中的井口回压（MPa）

% 固井流体参数：黏度/屈服应力沿用原代码，密度按8.11/施工工艺表更新
rou0 = 1.95; rou1 = 1.75; rou2 = 2.05; rou3 = 1.95; rou4 = 2.05; rou5 = 1.95; rou6 = 1.95; rou7 = 1.95; rou8 = 1.95; rou9 = 1.95; rou91 = 1.95; rou92 = 1.95; rou93 = 1.95; rou94 = 1.95; rou95 = 1.95; % 密度（g/cm³）
miu0 = 51 ; miu1 = 60; miu2 = 65; miu3 = 65; miu4 =  200; miu5 = 180; miu6 = 40; miu7 = 40; miu8 = 40; miu9 = 40; miu91 = 40; miu92 = 40; miu93 = 40; miu94 = 40; miu95 = 40; % 黏度（mPa·s）沿用旧参数
tau0 = 10 ; tau1 = 11; tau2 = 11; tau3 = 11; tau4 = 14; tau5 = 14; tau6 = 9; tau7 = 9.5; tau8 = 9.2; tau9 = 9; tau91 = 9.3; tau92 = 9.3; tau93 = 9.3; tau94 = 9.3; tau95 = 9.3;  % 屈服应力（Pa）沿用旧参数

%% 各流体体积参数 (来自7.1/7.2施工量与8.11施工过程模拟)
v1 = 28 * 1000;   % 1-先导浆体积,L
v2 = 16 * 1000;   % 2-低失水驱油隔离液1体积,L
v3 = 10 * 1000;   % 3-低失水驱油隔离液2体积,L
v4 = 39 * 1000;   % 4-领浆体积,L
v5 = 28 * 1000;   % 5-尾浆体积,L
v6 = 2 * 1000;  % 6-压塞液体积,L (2m³)
v7 = 25 * 1000; % 7-替浆钻井液体积,L (25m³)
v8 = 14 * 1000; % 8-保护液体积,L (14m³)
v9 = 1 * 1000;  % 9-基液体积,L (1m³)
v91 = 8 * 1000;   % 91-替浆钻井液1体积,L
v92 = 14 * 1000;  % 92-替浆钻井液2体积,L
v93 = 16 * 1000;  % 93-替浆钻井液3体积,L
v94 = 12.9 * 1000; % 94-替浆钻井液4体积,L
v95 = 0;           % 95-保留段，当前8.11无第五段替浆体积

% 计算总时间步数与分段数量
dt = 1;  % 【高精度优化】：大幅细化时间步长为 min
n_time = floor((v1./Pump_values(1) + v2./Pump_values(2) + v3./Pump_values(3) + v4./Pump_values(4) + v5./Pump_values(5) + v6./Pump_values(6) + v7./Pump_values(7) + v8./Pump_values(8) + v9./Pump_values(9) + v91./Pump_values(10) + v92./Pump_values(11) + v93./Pump_values(12) + v94./Pump_values(13) + v95./Pump_values(14)) / dt);
time = (v1./Pump_values(1) + v2./Pump_values(2) + v3./Pump_values(3) + v4./Pump_values(4) + v5./Pump_values(5) + v6./Pump_values(6) + v7./Pump_values(7) + v8./Pump_values(8) + v9./Pump_values(9) + v91./Pump_values(10) + v92./Pump_values(11) + v93./Pump_values(12) + v94./Pump_values(13) + v95./Pump_values(14)) / dt;
