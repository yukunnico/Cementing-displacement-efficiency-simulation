% HT1_004_T  HT1-004 固井循环温度场主模型。
% 输入：呼1-004井身结构.csv、施工排量/体积、各流体物性与地层热物性。
% 输出：Temperature_Result 结构体；其中温度场矩阵统一为 [井深分段 x 时间步]，温度单位为 ℃。
% 运行关系：本模型结束后自动写入 T_in.xlsx（管内）和 T_out.xlsx（环空），供 p_jaifang.m 调用。
clear; clc;

%% ==================== 1. 井身结构与施工基础参数 ====================
this_file_dir = fileparts(mfilename('fullpath'));
if isempty(this_file_dir)
    this_file_dir = pwd;
end
structure_file = fullfile(this_file_dir, '呼1-004井身结构.csv');
structure_data = readtable(structure_file);
c_depth = structure_data.depth_well_logging_m_;
n_segment = height(structure_data);

depth_threshold1 = 4025.73;
depth_threshold2 = 5243.21;
depth_threshold3 = 5578.00;
depth_threshold4 = 7378.05;
depth_threshold5 = 7521.00;
depth_threshold6 = 7660.00;
length_casing = 7660;

diameter_bit_out = zeros(1, length(c_depth));
diameter_casing_out = zeros(1, length(c_depth));
diameter_casing_in = zeros(1, length(c_depth));

idx1 = (c_depth <= depth_threshold1);
idx2 = (c_depth > depth_threshold1) & (c_depth <= depth_threshold2);
idx3 = (c_depth > depth_threshold2) & (c_depth <= depth_threshold3);
idx4 = (c_depth > depth_threshold3) & (c_depth <= depth_threshold4);
idx5 = (c_depth > depth_threshold4) & (c_depth <= depth_threshold5);
idx6 = (c_depth > depth_threshold5) & (c_depth <= depth_threshold6);
idx7 = (c_depth > depth_threshold6);

diameter_bit_out(idx1) = 245.37 * 0.001;
diameter_casing_out(idx1) = 149.2 * 0.001;
diameter_casing_in(idx1) = (149.2 - 9.65 * 2) * 0.001;

diameter_bit_out(idx2) = 245.37 * 0.001;
diameter_casing_out(idx2) = 127.0 * 0.001;
diameter_casing_in(idx2) = (127.0 - 9.65 * 2) * 0.001;

diameter_bit_out(idx3) = 245.37 * 0.001;
diameter_casing_out(idx3) = 168.3 * 0.001;
diameter_casing_in(idx3) = (168.3 - 15.88 * 2) * 0.001;

diameter_bit_out(idx4) = structure_data.annulus_radius_array_cm_(idx4) * 0.01;
diameter_casing_out(idx4) = 168.3 * 0.001;
diameter_casing_in(idx4) = (168.3 - 15.88 * 2) * 0.001;

diameter_bit_out(idx5) = structure_data.annulus_radius_array_cm_(idx5) * 0.01;
diameter_casing_out(idx5) = 139.7 * 0.001;
diameter_casing_in(idx5) = (139.7 - 15.88 * 2) * 0.001;

diameter_bit_out(idx6) = structure_data.annulus_radius_array_cm_(idx6) * 0.01;
diameter_casing_out(idx6) = 139.7 * 0.001;
diameter_casing_in(idx6) = (139.7 - 15.88 * 2) * 0.001;

diameter_bit_out(idx7) = structure_data.annulus_radius_array_cm_(idx7) * 0.01;
diameter_casing_out(idx7) = 0;
diameter_casing_in(idx7) = 0;

area_cout = pi * (diameter_bit_out .^ 2 - diameter_casing_out .^ 2) / 4;
area_cin = pi * diameter_casing_in .^ 2 / 4;
out_diam_bole = structure_data.annulus_radius_array_cm_ * 10;
pianxin = true;

volume_in_casing_seg = area_cin .* structure_data.length_segment_array_m_';
volume_in_annual_seg = area_cout .* structure_data.length_segment_array_m_';
idx_liner_top = find(c_depth >= depth_threshold2, 1, 'first');
idx_slurry_top = find(c_depth >= 3381, 1, 'first');
volume_in_fenggu = sum(volume_in_annual_seg(idx_liner_top:end));
volume_in_drilling_casing = sum(volume_in_casing_seg);
volume_in_annual = sum(volume_in_annual_seg);
volume_in_drilling_casing_L = volume_in_drilling_casing * 1000;
volume_in_annual_L = volume_in_annual * 1000;
volume_in_allfluid = volume_in_drilling_casing + sum(volume_in_annual_seg(idx_slurry_top:end));

pump_rate1 = 1.3 * 1000;
pump_rate2 = 1.3 * 1000;
pump_rate3 = 1.3 * 1000;
pump_rate4 = 1.2 * 1000;
pump_rate5 = 1.2 * 1000;
pump_rate6 = 1.0 * 1000;
pump_rate7 = 1.3 * 1000;
pump_rate8 = 1.2 * 1000;
pump_rate9 = 1.2 * 1000;
pump_rate91 = 1.2 * 1000;
pump_rate92 = 1.1 * 1000;
pump_rate93 = 0.9 * 1000;
pump_rate94 = 0.8 * 1000;
pump_rate95 = 0.7 * 1000;
Pump_values = [pump_rate1, pump_rate2, pump_rate3, pump_rate4, pump_rate5, pump_rate6, ...
    pump_rate7, pump_rate8, pump_rate9, pump_rate91, pump_rate92, pump_rate93, pump_rate94, pump_rate95];

pressure_back = 0;
pressure_back_static = 0;
backpressure_1 = 0; backpressure_2 = 0; backpressure_3 = 0; backpressure_4 = 0; backpressure_5 = 0;
backpressure_6 = 0; backpressure_7 = 0; backpressure_8 = 0; backpressure_9 = 0;
backpressure_91 = 0; backpressure_92 = 0; backpressure_93 = 0; backpressure_94 = 0; backpressure_95 = 0;
backpressure = [backpressure_1, backpressure_2, backpressure_3, backpressure_4, backpressure_5, ...
    backpressure_6, backpressure_7, backpressure_8, backpressure_9, backpressure_91, ...
    backpressure_92, backpressure_93, backpressure_94, backpressure_95];

rou0 = 1.90; rou1 = 1.80; rou2 = 1.95; rou3 = 1.75; rou4 = 1.92;
rou5 = 1.90; rou6 = 1.75; rou7 = 1.90; rou8 = 1.75; rou9 = 1.02;
rou91 = 1.90; rou92 = 1.90; rou93 = 1.90; rou94 = 1.90; rou95 = 1.90;

miu0 = 53; miu1 = 58; miu2 = 58; miu3 = 65; miu4 = 200;
miu5 = 180; miu6 = 50; miu7 = 50; miu8 = 50; miu9 = 50;
miu91 = 55; miu92 = 55; miu93 = 55; miu94 = 55; miu95 = 55;

tau0 = 8.5; tau1 = 9.8; tau2 = 9.8; tau3 = 10; tau4 = 14;
tau5 = 14; tau6 = 9; tau7 = 9.5; tau8 = 9.2; tau9 = 9;
tau91 = 9.5; tau92 = 9.5; tau93 = 9.5; tau94 = 9.5; tau95 = 9.5;

tag_interface_top_lead = 171;
tag_interface_top_tail = 253;

vjisuan = sum(structure_data.volume_annulus_L_(tag_interface_top_tail:end));
vjisuan2 = sum(volume_in_casing_seg);
v1 = 33.3 * 1000;
v2 = 16 * 1000;
v3 = 10 * 1000;
v4 = sum(structure_data.volume_annulus_L_(tag_interface_top_lead:tag_interface_top_tail));
v5 = sum(structure_data.volume_annulus_L_(tag_interface_top_tail + 1:end)) + 1 * 1000;
v6 = 2 * 1000;
v7 = 25 * 1000;
v8 = 10 * 1000;
v9 = 2 * 1000;
v91 = 13 * 1000;
v92 = 15 * 1000;
v93 = 10 * 1000;
v94 = 10 * 1000;
v95 = volume_in_drilling_casing_L - (v6 + v7 + v8 + v9 + v91 + v92 + v93 + v94) - 1 * 1000;
if v95 <= 0
    error('v95=%.2f L，替浆钻井液5体积为非正数，请检查前置体积和套管内容积。', v95);
end
stage_volume_L = [v1, v2, v3, v4, v5, v6, v7, v8, v9, v91, v92, v93, v94, v95];

dt = 1;
total_time_min = sum(stage_volume_L ./ Pump_values);
time = total_time_min / dt;
time_min = (0:dt:total_time_min).';
if abs(time_min(end) - total_time_min) > 1e-9
    time_min = [time_min; total_time_min];
end
% time_min 为施工离散时间序列（min）；pump_time_node 记录每种流体施工结束的时间步。
n_time = numel(time_min);
pump_time_node = zeros(14, 1);
pump_time_node(1) = (v1/Pump_values(1)) / dt;
pump_time_node(2) = (v1/Pump_values(1) + v2/Pump_values(2)) / dt;
pump_time_node(3) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3)) / dt;
pump_time_node(4) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4)) / dt;
pump_time_node(5) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5)) / dt;
pump_time_node(6) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6)) / dt;
pump_time_node(7) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7)) / dt;
pump_time_node(8) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8)) / dt;
pump_time_node(9) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8) + v9/Pump_values(9)) / dt;
pump_time_node(10) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8) + v9/Pump_values(9) + v91/Pump_values(10)) / dt;
pump_time_node(11) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8) + v9/Pump_values(9) + v91/Pump_values(10) + v92/Pump_values(11)) / dt;
pump_time_node(12) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8) + v9/Pump_values(9) + v91/Pump_values(10) + v92/Pump_values(11) + v93/Pump_values(12)) / dt;
pump_time_node(13) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8) + v9/Pump_values(9) + v91/Pump_values(10) + v92/Pump_values(11) + v93/Pump_values(12) + v94/Pump_values(13)) / dt;
pump_time_node(14) = (v1/Pump_values(1) + v2/Pump_values(2) + v3/Pump_values(3) + v4/Pump_values(4) + v5/Pump_values(5) + v6/Pump_values(6) + v7/Pump_values(7) + v8/Pump_values(8) + v9/Pump_values(9) + v91/Pump_values(10) + v92/Pump_values(11) + v93/Pump_values(12) + v94/Pump_values(13) + v95/Pump_values(14)) / dt;

%% 数组初始化与高精度地热场构建
% vertical_length_all_grid：每个测深分段对应的垂深长度（m）；
% TVD_cum：从井口到各分段的累计垂深（m），用于地温梯度计算。
vertical_length_all_grid = zeros(1, n_segment);
cos_deg = zeros(1, n_segment);
for i = 1:n_segment
    cos_deg(i) = cos(deg2rad(structure_data.deg_for_logging_degree_(i)));
    vertical_length_all_grid(i) = abs(structure_data.length_segment_array_m_(i) * cos_deg(i));
end
TVD_cum = cumsum(vertical_length_all_grid);

% 环空几何：d_i=管柱外径，d_o=井眼/套管内径，单位 m。
% 后续温度网格中 d/d_o/d_i 分别表示井眼直径、管柱外径和管柱内径。
d_i = zeros(1, n_segment);
d_o = zeros(1, n_segment);
for i = 1:n_segment
    d_i(i) = diameter_casing_out(i);
    d_o(i) = structure_data.annulus_radius_array_cm_(i) / 100;
end

% 累计注入体积，单位 L：环空只需跟踪前五种流体，套管内跟踪全部施工阶段。
volume_injected_1_list = zeros(1, n_time);
volume_injected_2_list = zeros(1, n_time);
volume_injected_3_list = zeros(1, n_time);
volume_injected_4_list = zeros(1, n_time);
volume_injected_5_list = zeros(1, n_time);

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

volume_into_1_annulus = zeros(1, n_time);
volume_into_2_annulus = zeros(1, n_time);
volume_into_3_annulus = zeros(1, n_time);
volume_into_4_annulus = zeros(1, n_time);
volume_into_5_annulus = zeros(1, n_time);

% 环空界面追踪：tag=界面所在分段号；residual_volume/residual_height=
% 该分段内部分充填的体积（L）/沿井深高度（m）；tag_interface 为 [n_time x 5]。
depth_tag_liquid_0_1_list = zeros(1, n_time); residual_volume_0_1_list = zeros(1, n_time); residual_height_0_1_list = zeros(1, n_time); vertical_residual_height_0_1_list = zeros(1, n_time);
depth_tag_liquid_1_2_list = zeros(1, n_time); residual_volume_1_2_list = zeros(1, n_time); residual_height_1_2_list = zeros(1, n_time); vertical_residual_height_1_2_list = zeros(1, n_time);
depth_tag_liquid_2_3_list = zeros(1, n_time); residual_volume_2_3_list = zeros(1, n_time); residual_height_2_3_list = zeros(1, n_time); vertical_residual_height_2_3_list = zeros(1, n_time);
depth_tag_liquid_3_4_list = zeros(1, n_time); residual_volume_3_4_list = zeros(1, n_time); residual_height_3_4_list = zeros(1, n_time); vertical_residual_height_3_4_list = zeros(1, n_time);
depth_tag_liquid_4_5_list = zeros(1, n_time); residual_volume_4_5_list = zeros(1, n_time); residual_height_4_5_list = zeros(1, n_time); vertical_residual_height_4_5_list = zeros(1, n_time);
tag_interface = zeros(n_time, 5);

% 套管内界面追踪，逻辑与环空相同；tag_interface_casing 的 14 列
% 对应 0-1、1-2 ... 94-95 的流体交界面，矩阵维度为 [n_time x 14]。
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

% 压力模型共享物性矩阵，维度均为 [n_time x n_segment]：
% rou= g/cm^3，miu= mPa*s，tau= Pa，velo= m/s。
rou_annulus_all_time = zeros(n_time, n_segment);
miu_annulus_all_time = zeros(n_time, n_segment);
tau_annulus_all_time = zeros(n_time, n_segment);
velo_annulus_all_time = zeros(n_time, n_segment);
rou_casing_all_time = zeros(n_time, n_segment);
miu_casing_all_time = zeros(n_time, n_segment);
tau_casing_all_time = zeros(n_time, n_segment);
velo_casing_all_time = zeros(n_time, n_segment);

Ff_a = zeros(n_time, n_segment);
flow_pattern_a = zeros(n_time, n_segment);
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

Pump_values_time_list = zeros(1, n_time);
volume_injected_all_list = zeros(1, n_time);
backpressure_time_list = zeros(1, n_time);

%% ==================== 2. 温度场物性与计算控制 ====================
% 温度边界与地温梯度：Temp_in 为井口注入流体温度，Temp_surface 为地表地层温度，
% Temp_gradient 为地温梯度（℃/m）。
Temp_in = 16.0;
Temp_surface = 16.0;
Temp_gradient = 0.01760;
% 是否将本次求解的当前温度场写入压力模型读取的 Excel 文件。
export_temperature_xlsx = true;

% 固体热物性，单位依次为 kg/m^3、J/(kg*K)、W/(m*K)。
% 地层在未提供岩性测井/岩心试验时按代表性砂岩处理；砂岩 cp=800、k=2.38
% 取自 298 K 试验数据。套管/钻杆按 20-100 degC 碳钢处理。
% 资料：砂岩 https://html.rhhz.net/bmpg/html/20170308.htm；
% 碳钢 https://www.tdiinternational.com/technical-source-product-info/material-data-sheets-spec-sheets/technical-data-sheet-c-carbon-steel/
rho_f = 2640;
cp_f = 800;
k_f = 2.38;

rho_casing = 7850;
cp_casing = 502;
k_casing = 52;
rho_drill_pipe = 7850;
cp_drill_pipe = 502;
k_drill_pipe = 52;
% 当前模型将套管壁视为准稳态导热层，只有 k_casing 进入 U_ap；
% 钻杆未单独建立几何/温度节点，rho_drill_pipe 与 cp_drill_pipe 仅作为材料台账输出。
k_steel = k_casing;
shoe_pressure_drop = 0;
include_friction_heat = true;

% number_tr 为径向温度节点数：1=管内，2=环空，3=井壁，4~number_tr=地层。
% radial_log_step 控制井壁外地层径向网格的指数扩张间距。
number_tr = 14;
radial_log_step = 0.20;

% stage_ids 为泵注阶段编号；fluid_ids 额外包含初始钻井液 0。
% rho/mu/tau 表与 fluid_ids 一一对应，分别使用 kg/m^3、Pa*s、Pa 参与传热求解。
stage_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 91, 92, 93, 94, 95];
fluid_ids = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 91, 92, 93, 94, 95];
rho_table = [rou0, rou1, rou2, rou3, rou4, rou5, rou6, rou7, rou8, rou9, rou91, rou92, rou93, rou94, rou95] * 1000;
mu_table = [miu0, miu1, miu2, miu3, miu4, miu5, miu6, miu7, miu8, miu9, miu91, miu92, miu93, miu94, miu95] / 1000;
tau_table = [tau0, tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau91, tau92, tau93, tau94, tau95];
% 流体热物性与 rho/mu/tau 使用同一套界面追踪。下列值为公开资料范围内的
% 工程代表值，适用于本井已确认的水基体系；真实设计值应以现场配方实验为准。
% 水基基液(9号)：按 20 degC 水的 cp=4185、k=0.598 并考虑 1.02 g/cm^3 添加剂取值。
% 固井浆：新拌 Class G/H 水泥浆公开范围 cp=1500-2000、k=0.7-1.0；
% 高密度水基钻井液/隔离液按密度和固相含量分别取值。
% 资料：https://www.iapws.org/relguide/ThCond.html；
% https://www.sciencedirect.com/science/article/pii/S2214785321054559
fluid_name = {'Initial water-based drilling mud'; 'Stage 1 water-based spacer'; ...
    'Stage 2 high-density spacer'; 'Stage 3 water-based spacer'; ...
    'Lead cement slurry'; 'Tail cement slurry'; 'Stage 6 plug fluid'; ...
    'Stage 7 displacement mud'; 'Stage 8 protective fluid'; ...
    'Stage 9 water-based base liquid'; 'Displacement mud 1'; ...
    'Displacement mud 2'; 'Displacement mud 3'; 'Displacement mud 4'; ...
    'Displacement mud 5'};
cp_table = [2500, 2700, 2300, 2700, 1800, 1850, 2700, 2500, 2700, 3900, ...
    2500, 2500, 2500, 2500, 2500];
k_table = [0.90, 0.78, 0.92, 0.78, 0.88, 0.86, 0.78, 0.90, 0.78, 0.60, ...
    0.90, 0.90, 0.90, 0.90, 0.90];
if any(~isfinite(cp_table)) || any(cp_table <= 0) || any(~isfinite(k_table)) || any(k_table <= 0)
    error('Fluid thermal property table contains invalid cp or k values.');
end

stage_duration_min = stage_volume_L ./ Pump_values;
stage_end_min = cumsum(stage_duration_min);
stage_start_min = [0, stage_end_min(1:end-1)];
stage_start_volume_L = [0, cumsum(stage_volume_L(1:end-1))];
stage_index_time = zeros(n_time, 1);
fluid_id_time = zeros(n_time, 1);
Q_time = zeros(n_time, 1);
cumulative_volume_L = zeros(n_time, 1);
for t = 1:n_time
    [stage_index_time(t), cumulative_volume_L(t)] = local_stage_and_volume( ...
        time_min(t), stage_start_min, stage_end_min, stage_start_volume_L, Pump_values, stage_volume_L);
    fluid_id_time(t) = stage_ids(stage_index_time(t));
    Q_time(t) = Pump_values(stage_index_time(t)) / 60000;
end
Pump_values_time_list = Q_time.' * 60000;
Pump_values_time_list_m3_s = Q_time.';
backpressure_time_list = backpressure(stage_index_time).';
volume_injected_all_list = cumulative_volume_L.';

volume_all_from_bottom = zeros(1, n_segment);
length_all_from_bottom = zeros(1, n_segment);
volume_all_from_top_casing = zeros(1, n_segment);
length_all_from_top_casing = zeros(1, n_segment);
for i = 1:n_segment
    volume_all_from_bottom(i) = sum(structure_data.volume_annulus_L_(i:end));
    length_all_from_bottom(i) = sum(structure_data.length_segment_array_m_(i:end));
    volume_all_from_top_casing(i) = sum(area_cin(1:i) .* structure_data.length_segment_array_m_(1:i)') * 1000;
    length_all_from_top_casing(i) = sum(structure_data.length_segment_array_m_(1:i));
end

volume_injected_1_list = max(volume_injected_all_list - volume_in_drilling_casing_L, 0);
volume_injected_2_list = max(volume_injected_all_list - volume_in_drilling_casing_L - v1, 0);
volume_injected_3_list = max(volume_injected_all_list - volume_in_drilling_casing_L - v1 - v2, 0);
volume_injected_4_list = max(volume_injected_all_list - volume_in_drilling_casing_L - v1 - v2 - v3, 0);
volume_injected_5_list = max(volume_injected_all_list - volume_in_drilling_casing_L - v1 - v2 - v3 - v4, 0);
volume_into_1_annulus = volume_injected_1_list;
volume_into_2_annulus = volume_injected_2_list;
volume_into_3_annulus = volume_injected_3_list;
volume_into_4_annulus = volume_injected_4_list;
volume_into_5_annulus = volume_injected_5_list;

[depth_tag_liquid_0_1_list, residual_volume_0_1_list, residual_height_0_1_list, vertical_residual_height_0_1_list] = ...
    local_annulus_interface(volume_injected_1_list, volume_all_from_bottom, structure_data.square_annulus_dm2_, structure_data.deg_for_logging_degree_);
[depth_tag_liquid_1_2_list, residual_volume_1_2_list, residual_height_1_2_list, vertical_residual_height_1_2_list] = ...
    local_annulus_interface(volume_injected_2_list, volume_all_from_bottom, structure_data.square_annulus_dm2_, structure_data.deg_for_logging_degree_);
[depth_tag_liquid_2_3_list, residual_volume_2_3_list, residual_height_2_3_list, vertical_residual_height_2_3_list] = ...
    local_annulus_interface(volume_injected_3_list, volume_all_from_bottom, structure_data.square_annulus_dm2_, structure_data.deg_for_logging_degree_);
[depth_tag_liquid_3_4_list, residual_volume_3_4_list, residual_height_3_4_list, vertical_residual_height_3_4_list] = ...
    local_annulus_interface(volume_injected_4_list, volume_all_from_bottom, structure_data.square_annulus_dm2_, structure_data.deg_for_logging_degree_);
[depth_tag_liquid_4_5_list, residual_volume_4_5_list, residual_height_4_5_list, vertical_residual_height_4_5_list] = ...
    local_annulus_interface(volume_injected_5_list, volume_all_from_bottom, structure_data.square_annulus_dm2_, structure_data.deg_for_logging_degree_);
tag_interface = [depth_tag_liquid_0_1_list(:), depth_tag_liquid_1_2_list(:), depth_tag_liquid_2_3_list(:), ...
    depth_tag_liquid_3_4_list(:), depth_tag_liquid_4_5_list(:)];

volume_injected_casing_1_list = max(volume_injected_all_list, 0);
volume_injected_casing_2_list = max(volume_injected_all_list - v1, 0);
volume_injected_casing_3_list = max(volume_injected_all_list - v1 - v2, 0);
volume_injected_casing_4_list = max(volume_injected_all_list - v1 - v2 - v3, 0);
volume_injected_casing_5_list = max(volume_injected_all_list - v1 - v2 - v3 - v4, 0);
volume_injected_casing_6_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5, 0);
volume_injected_casing_7_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6, 0);
volume_injected_casing_8_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7, 0);
volume_injected_casing_9_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7 - v8, 0);
volume_injected_casing_91_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7 - v8 - v9, 0);
volume_injected_casing_92_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7 - v8 - v9 - v91, 0);
volume_injected_casing_93_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7 - v8 - v9 - v91 - v92, 0);
volume_injected_casing_94_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7 - v8 - v9 - v91 - v92 - v93, 0);
volume_injected_casing_95_list = max(volume_injected_all_list - v1 - v2 - v3 - v4 - v5 - v6 - v7 - v8 - v9 - v91 - v92 - v93 - v94, 0);

[depth_tag_casing_0_1_list, residual_volume_casing_0_1_list, residual_height_casing_0_1_list] = ...
    local_casing_interface(volume_injected_casing_1_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_1_2_list, residual_volume_casing_1_2_list, residual_height_casing_1_2_list] = ...
    local_casing_interface(volume_injected_casing_2_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_2_3_list, residual_volume_casing_2_3_list, residual_height_casing_2_3_list] = ...
    local_casing_interface(volume_injected_casing_3_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_3_4_list, residual_volume_casing_3_4_list, residual_height_casing_3_4_list] = ...
    local_casing_interface(volume_injected_casing_4_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_4_5_list, residual_volume_casing_4_5_list, residual_height_casing_4_5_list] = ...
    local_casing_interface(volume_injected_casing_5_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_5_6_list, residual_volume_casing_5_6_list, residual_height_casing_5_6_list] = ...
    local_casing_interface(volume_injected_casing_6_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_6_7_list, residual_volume_casing_6_7_list, residual_height_casing_6_7_list] = ...
    local_casing_interface(volume_injected_casing_7_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_7_8_list, residual_volume_casing_7_8_list, residual_height_casing_7_8_list] = ...
    local_casing_interface(volume_injected_casing_8_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_8_9_list, residual_volume_casing_8_9_list, residual_height_casing_8_9_list] = ...
    local_casing_interface(volume_injected_casing_9_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_9_91_list, residual_volume_casing_9_91_list, residual_height_casing_9_91_list] = ...
    local_casing_interface(volume_injected_casing_91_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_91_92_list, residual_volume_casing_91_92_list, residual_height_casing_91_92_list] = ...
    local_casing_interface(volume_injected_casing_92_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_92_93_list, residual_volume_casing_92_93_list, residual_height_casing_92_93_list] = ...
    local_casing_interface(volume_injected_casing_93_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_93_94_list, residual_volume_casing_93_94_list, residual_height_casing_93_94_list] = ...
    local_casing_interface(volume_injected_casing_94_list, volume_all_from_top_casing, area_cin);

[depth_tag_casing_94_95_list, residual_volume_casing_94_95_list, residual_height_casing_94_95_list] = ...
    local_casing_interface(volume_injected_casing_95_list, volume_all_from_top_casing, area_cin);

tag_interface_casing = [depth_tag_casing_0_1_list(:), depth_tag_casing_1_2_list(:), depth_tag_casing_2_3_list(:), ...
    depth_tag_casing_3_4_list(:), depth_tag_casing_4_5_list(:), depth_tag_casing_5_6_list(:), ...
    depth_tag_casing_6_7_list(:), depth_tag_casing_7_8_list(:), depth_tag_casing_8_9_list(:), ...
    depth_tag_casing_9_91_list(:), depth_tag_casing_91_92_list(:), depth_tag_casing_92_93_list(:), ...
    depth_tag_casing_93_94_list(:), depth_tag_casing_94_95_list(:)];

%% ==================== 3. 温度场网格初始化 ====================
% 深度方向网格：depth/TVD/dz 为 [nz x 1]；几何量 d/d_o/d_i 与面积 A_ann/A_pipe
% 同样按井深分段存储。无内管的纯裸眼段会在 valid_string 筛选中移除。
depth = c_depth(:);
TVD = TVD_cum(:);
dz = structure_data.length_segment_array_m_(:);
d = diameter_bit_out(:);
d_o = diameter_casing_out(:);
d_i = diameter_casing_in(:);
A_ann = area_cout(:);
A_pipe = area_cin(:);

valid_string = d_o > 0 & d_i > 0 & d > d_o & A_ann > 0 & A_pipe > 0;
if ~all(valid_string)
    last_valid = find(valid_string, 1, 'last');
    if any(~valid_string(1:last_valid))
        error('井身结构中间存在无效几何段，请检查井径、管柱外径、管柱内径或环空面积。');
    end
    depth = depth(1:last_valid);
    TVD = TVD(1:last_valid);
    dz = dz(1:last_valid);
    d = d(1:last_valid);
    d_o = d_o(1:last_valid);
    d_i = d_i(1:last_valid);
    A_ann = A_ann(1:last_valid);
    A_pipe = A_pipe(1:last_valid);
end
nz = numel(depth);

pipe_volume_to_node_L = cumsum(A_pipe .* dz) * 1000;
total_pipe_volume_L = pipe_volume_to_node_L(end);
annulus_volume_from_bottom_L = zeros(nz, 1);
for i = 1:nz
    if i < nz
        annulus_volume_from_bottom_L(i) = sum(A_ann(i+1:end) .* dz(i+1:end)) * 1000;
    end
end

% 温度场内部矩阵均采用 [nz x n_time]，与压力模型的 [n_time x n_segment] 相反；
% 回填到 rou_*_all_time 等共享变量时会转置。
fluid_id_pipe = zeros(nz, n_time);
fluid_id_annulus = zeros(nz, n_time);
rho_pipe = zeros(nz, n_time);
rho_annulus = zeros(nz, n_time);
mu_pipe = zeros(nz, n_time);
mu_annulus = zeros(nz, n_time);
tau_pipe = zeros(nz, n_time);
tau_annulus = zeros(nz, n_time);
cp_pipe = zeros(nz, n_time);
cp_annulus = zeros(nz, n_time);
k_pipe = zeros(nz, n_time);
k_annulus = zeros(nz, n_time);

for t = 1:n_time
    for i = 1:nz
        volume_marker_pipe = cumulative_volume_L(t) - pipe_volume_to_node_L(i);
        volume_marker_annulus = cumulative_volume_L(t) - total_pipe_volume_L - annulus_volume_from_bottom_L(i);

        fluid_id_pipe(i, t) = local_fluid_id_from_volume(volume_marker_pipe, stage_volume_L, stage_ids);
        fluid_id_annulus(i, t) = local_fluid_id_from_volume(volume_marker_annulus, stage_volume_L, stage_ids);

        [rho_pipe(i, t), mu_pipe(i, t), tau_pipe(i, t), cp_pipe(i, t), k_pipe(i, t)] = local_fluid_property( ...
            fluid_id_pipe(i, t), fluid_ids, rho_table, mu_table, tau_table, cp_table, k_table);
        [rho_annulus(i, t), mu_annulus(i, t), tau_annulus(i, t), cp_annulus(i, t), k_annulus(i, t)] = local_fluid_property( ...
            fluid_id_annulus(i, t), fluid_ids, rho_table, mu_table, tau_table, cp_table, k_table);
    end
end

rou_annulus_all_time(:, 1:nz) = (rho_annulus.' / 1000);
miu_annulus_all_time(:, 1:nz) = (mu_annulus.' * 1000);
tau_annulus_all_time(:, 1:nz) = tau_annulus.';
rou_casing_all_time(:, 1:nz) = (rho_pipe.' / 1000);
miu_casing_all_time(:, 1:nz) = (mu_pipe.' * 1000);
tau_casing_all_time(:, 1:nz) = tau_pipe.';
for t = 1:n_time
    velo_annulus_all_time(t, 1:nz) = Q_time(t) ./ max(A_ann.', eps);
    velo_casing_all_time(t, 1:nz) = Q_time(t) ./ max(A_pipe.', eps);
end
rou_annulus_all_time_kg_m3 = rou_annulus_all_time * 1000;
rou_casing_all_time_kg_m3 = rou_casing_all_time * 1000;

% T_yuan 为未受施工扰动的地层初始温度场（℃）。
T_yuan = Temp_surface + Temp_gradient .* TVD;
r = zeros(nz, number_tr);
r(:, 1) = d_i ./ 2;
r(:, 2) = d_o ./ 2;
r(:, 3) = d ./ 2;
for j = 4:number_tr
    r(:, j) = r(:, 3) .* exp((j - 3) * radial_log_step);
end

% T 的列依次为管内、环空、井壁及地层径向节点；末列保持为初始地层边界温度。
T = repmat(T_yuan, 1, number_tr);
temp_pi = zeros(nz, n_time);
temp_po = zeros(nz, n_time);
temp_wall = zeros(nz, n_time);
temp_pf = zeros(nz, n_time);
bht = zeros(n_time, 1);
h1_time = zeros(nz, n_time);
h2_time = zeros(nz, n_time);
h3_time = zeros(nz, n_time);
U_ap_time = zeros(nz, n_time);

temp_pi(:, 1) = T(:, 1);
temp_po(:, 1) = T(:, 2);
temp_wall(:, 1) = T(:, 3);
temp_pf(:, 1) = T(:, ceil(number_tr / 2));
bht(1) = T(end, 2);

% 稀疏线性方程组仅求解前 number_tr-1 个径向节点；最外层地层节点作为定温边界。
idx_unknown = @(iz, jr) (iz - 1) * (number_tr - 1) + jr;
unknown_count = nz * (number_tr - 1);

%% ==================== 4. 管内-环空-地层耦合温度场求解 ====================
% 每个时间步先计算管内/环空压降与对流换热系数，再组装能量守恒的稀疏方程组。
% include_friction_heat=true 时，将 dpressure_loss * Q 作为摩擦生热项加入右端项。
for t = 2:n_time
    Q = Q_time(t - 1);
    dt_step_s = (time_min(t) - time_min(t - 1)) * 60;

    dpressure_loss_i = zeros(nz, 1);
    dpressure_loss_o = zeros(nz, 1);
    h1 = zeros(nz, 1);
    h2 = zeros(nz, 1);
    h3 = zeros(nz, 1);
    U_ap = zeros(nz, 1);

    for i = 1:nz
        [dpressure_loss_i(i), h1(i)] = local_caculate_pipe(Q, d_i(i), ...
            mu_pipe(i,t-1), tau_pipe(i,t-1), cp_pipe(i,t-1), k_pipe(i,t-1), rho_pipe(i,t-1), dz(i), A_pipe(i));
        [dpressure_loss_o(i), h2(i), h3(i)] = local_caculate_annuli(Q, d(i), d_o(i), A_ann(i), ...
            mu_annulus(i,t-1), tau_annulus(i,t-1), cp_annulus(i,t-1), k_annulus(i,t-1), rho_annulus(i,t-1), dz(i));
        U_ap(i) = 1 / (1/max(h1(i), eps) + d_i(i)/(d_o(i)*max(h2(i), eps)) + ...
            d_i(i)*log(d_o(i)/d_i(i))/(2*k_steel));
    end
    h1_time(:, t - 1) = h1;
    h2_time(:, t - 1) = h2;
    h3_time(:, t - 1) = h3;
    U_ap_time(:, t - 1) = U_ap;

    % 采用三元组形式累积稀疏矩阵系数，最后一次性构造 A，避免在循环中反复扩容。
    row = zeros(unknown_count * 5, 1);
    col = zeros(unknown_count * 5, 1);
    val = zeros(unknown_count * 5, 1);
    rhs = zeros(unknown_count, 1);
    nnz_count = 0;

    for i = 1:nz
        for j = 1:(number_tr - 1)
            eq = idx_unknown(i, j);

            if j == 1
                if i == 1
                    % 首个管内节点位于 30 m，顶部边界取该实际节点的初始地层温度，
                    % 避免将井口 0 m 的 Temp_in 直接施加到 30 m 节点。
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 1), 1);
                    rhs(eq) = T_yuan(i);
                else
                    rho_cell = rho_pipe(i, t - 1);
                    cp_cell = cp_pipe(i, t - 1);
                    adv = rho_cell * Q * cp_cell / max(dz(i), eps);
                    storage = A_pipe(i) * rho_cell * cp_cell / dt_step_s;
                    exchange = pi * d_i(i) * U_ap(i);

                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 1), -exchange - adv - storage);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i - 1, 1), adv);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 2), exchange);

                    rhs(eq) = -storage * T(i, 1);
                    if include_friction_heat
                        rhs(eq) = rhs(eq) - dpressure_loss_i(i) * Q;
                    end
                end

            elseif j == 2
                if i == nz
                    rho_shoe = rho_annulus(i, t - 1);
                    cp_shoe = cp_annulus(i, t - 1);
                    deltaT_shoe = shoe_pressure_drop / max(rho_shoe * cp_shoe, eps);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 2), 1);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 1), -1);
                    rhs(eq) = deltaT_shoe;
                else
                    rho_cell = rho_annulus(i, t - 1);
                    cp_cell = cp_annulus(i, t - 1);
                    adv = rho_cell * Q * cp_cell / max(dz(i), eps);
                    storage = A_ann(i) * rho_cell * cp_cell / dt_step_s;
                    pipe_exchange = pi * d_i(i) * U_ap(i);
                    wall_exchange = pi * d(i) * h3(i);

                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 2), -pipe_exchange - wall_exchange - adv - storage);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i + 1, 2), adv);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 1), pipe_exchange);
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 3), wall_exchange);

                    rhs(eq) = -storage * T(i, 2);
                    if include_friction_heat
                        rhs(eq) = rhs(eq) - dpressure_loss_o(i) * Q;
                    end
                end

            elseif i == 1
                [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, j), 1);
                rhs(eq) = T_yuan(i);

            elseif j == 3
                dr = r(i, 4) - r(i, 3);
                [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 3), -k_f/dr - h3(i));
                [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 4), k_f/dr);
                [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, 2), h3(i));

            else
                dr_in = r(i, j) - r(i, j - 1);
                dr_out = r(i, j + 1) - r(i, j);
                r_now = r(i, j);
                r_w = 0.5 * (r(i, j - 1) + r(i, j));
                r_e = 0.5 * (r(i, j) + r(i, j + 1));
                dr_cv = 0.5 * (dr_in + dr_out);
                storage_f = rho_f * cp_f / dt_step_s;
                inward = -k_f * r_w / (r_now * dr_cv * dr_in);
                outward = -k_f * r_e / (r_now * dr_cv * dr_out);
                center = storage_f - inward - outward;

                [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, j), center);
                [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, j - 1), inward);
                if j == number_tr - 1
                    rhs(eq) = storage_f * T(i, j) - outward * T_yuan(i);
                else
                    [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, idx_unknown(i, j + 1), outward);
                    rhs(eq) = storage_f * T(i, j);
                end
            end
        end
    end

    A = sparse(row(1:nnz_count), col(1:nnz_count), val(1:nnz_count), unknown_count, unknown_count);
    solution = A \ rhs;
    if any(~isfinite(solution))
        error('温度场在第 %d 个时间步求解失败。', t);
    end

    T(:, 1:number_tr - 1) = reshape(solution, number_tr - 1, nz).';
    T(:, number_tr) = T_yuan;

    temp_pi(:, t) = T(:, 1);
    temp_po(:, t) = T(:, 2);
    temp_wall(:, t) = T(:, 3);
    temp_pf(:, t) = T(:, ceil(number_tr / 2));
    bht(t) = T(end, 2);
end

if n_time > 1
    h1_time(:, n_time) = h1_time(:, n_time - 1);
    h2_time(:, n_time) = h2_time(:, n_time - 1);
    h3_time(:, n_time) = h3_time(:, n_time - 1);
    U_ap_time(:, n_time) = U_ap_time(:, n_time - 1);
end

% depth(1) 是首个计算节点的井深（本井为 30 m），并非井口 0 m。
% 因此将 30 m 环空节点温度沿返流方向延拓至井口，得到真实的环空出口温度。
surface_extension_length_m = depth(1);
temp_po_surface = temp_po(1, :);
surface_wall_temperature_C = Temp_surface + 0.5 * Temp_gradient * TVD(1);
for t = 2:n_time
    heat_capacity_rate = rho_annulus(1, t - 1) * Q_time(t - 1) * cp_annulus(1, t - 1);
    pipe_conductance = pi * d_i(1) * U_ap_time(1, t - 1);
    wall_conductance = pi * d(1) * h3_time(1, t - 1);
    total_conductance = pipe_conductance + wall_conductance;
    if heat_capacity_rate > eps && total_conductance > eps
        equilibrium_temperature = (pipe_conductance * Temp_in + ...
            wall_conductance * surface_wall_temperature_C) / total_conductance;
        attenuation = exp(-total_conductance * surface_extension_length_m / heat_capacity_rate);
        temp_po_surface(t) = equilibrium_temperature + ...
            (temp_po(1, t) - equilibrium_temperature) * attenuation;
    end
end

%% ==================== 5. 结果输出与绘图 ====================
% Temperature_Result 是本脚本的标准输出接口：温度相关矩阵为 [nz x n_time]；
% 施工阶段、界面、流体物性和排量同时保存，便于与压力模型逐时核对。
fprintf('\n=== HT1-004 固井循环温度场计算完成 ===\n');
fprintf('  井深: %.1f m\n', depth(end));
fprintf('  时间步数: %d, 总时间: %.1f min\n', n_time, time_min(end));
fprintf('  初始井底静温: %.2f degC\n', T_yuan(end));
fprintf('  末时刻管内井底温度: %.2f degC\n', temp_pi(end, end));
fprintf('  末时刻环空井底温度: %.2f degC\n', temp_po(end, end));
fprintf('  末时刻环空首节点温度 (%.1f m): %.2f degC\n', depth(1), temp_po(1, end));
fprintf('  末时刻井口环空出口温度 (0 m): %.2f degC\n', temp_po_surface(end));

Temperature_Result.depth_m = depth;
Temperature_Result.TVD_m = TVD;
Temperature_Result.time_min = time_min;
Temperature_Result.initial_formation_C = T_yuan;
Temperature_Result.pipe_C = temp_pi;
Temperature_Result.annulus_C = temp_po;
Temperature_Result.wall_C = temp_wall;
Temperature_Result.formation_mid_C = temp_pf;
Temperature_Result.bottom_annulus_C = bht;
Temperature_Result.radial_position_m = r;
Temperature_Result.stage_index = stage_index_time;
Temperature_Result.fluid_id_time = fluid_id_time;
Temperature_Result.fluid_id_pipe = fluid_id_pipe;
Temperature_Result.fluid_id_annulus = fluid_id_annulus;
Temperature_Result.pump_rate_m3_s = Q_time;
Temperature_Result.cumulative_volume_L = cumulative_volume_L;
Temperature_Result.Pump_values_time_list_L_min = Pump_values_time_list;
Temperature_Result.Pump_values_time_list_m3_s = Pump_values_time_list_m3_s;
Temperature_Result.backpressure_time_list_MPa = backpressure_time_list;
Temperature_Result.volume_injected_all_list_L = volume_injected_all_list;
Temperature_Result.volume_all_from_bottom_L = volume_all_from_bottom;
Temperature_Result.length_all_from_bottom_m = length_all_from_bottom;
Temperature_Result.volume_all_from_top_casing_L = volume_all_from_top_casing;
Temperature_Result.length_all_from_top_casing_m = length_all_from_top_casing;
Temperature_Result.volume_injected_annulus_lists_L = [volume_injected_1_list; volume_injected_2_list; ...
    volume_injected_3_list; volume_injected_4_list; volume_injected_5_list];
Temperature_Result.volume_injected_casing_lists_L = [volume_injected_casing_1_list; volume_injected_casing_2_list; ...
    volume_injected_casing_3_list; volume_injected_casing_4_list; volume_injected_casing_5_list; ...
    volume_injected_casing_6_list; volume_injected_casing_7_list; volume_injected_casing_8_list; ...
    volume_injected_casing_9_list; volume_injected_casing_91_list; volume_injected_casing_92_list; ...
    volume_injected_casing_93_list; volume_injected_casing_94_list; volume_injected_casing_95_list];
Temperature_Result.volume_into_annulus_lists_L = [volume_into_1_annulus; volume_into_2_annulus; ...
    volume_into_3_annulus; volume_into_4_annulus; volume_into_5_annulus];
Temperature_Result.tag_interface_annulus = tag_interface;
Temperature_Result.tag_interface_casing = tag_interface_casing;
Temperature_Result.rou_annulus_all_time_g_cm3 = rou_annulus_all_time;
Temperature_Result.miu_annulus_all_time_mPa_s = miu_annulus_all_time;
Temperature_Result.tau_annulus_all_time_Pa = tau_annulus_all_time;
Temperature_Result.velo_annulus_all_time_m_s = velo_annulus_all_time;
Temperature_Result.rou_casing_all_time_g_cm3 = rou_casing_all_time;
Temperature_Result.miu_casing_all_time_mPa_s = miu_casing_all_time;
Temperature_Result.tau_casing_all_time_Pa = tau_casing_all_time;
Temperature_Result.velo_casing_all_time_m_s = velo_casing_all_time;
Temperature_Result.rho_pipe_kg_m3 = rho_pipe;
Temperature_Result.rho_annulus_kg_m3 = rho_annulus;
Temperature_Result.mu_pipe_Pa_s = mu_pipe;
Temperature_Result.mu_annulus_Pa_s = mu_annulus;
Temperature_Result.tau_pipe_Pa = tau_pipe;
Temperature_Result.tau_annulus_Pa = tau_annulus;
Temperature_Result.cp_pipe_J_kgK = cp_pipe;
Temperature_Result.cp_annulus_J_kgK = cp_annulus;
Temperature_Result.k_pipe_W_mK = k_pipe;
Temperature_Result.k_annulus_W_mK = k_annulus;
Temperature_Result.h_pipe_inner_W_m2K = h1_time;
Temperature_Result.h_pipe_outer_W_m2K = h2_time;
Temperature_Result.h_annulus_wall_W_m2K = h3_time;
Temperature_Result.temperature_field_file_in = fullfile(this_file_dir, 'T_in.xlsx');
Temperature_Result.temperature_field_file_out = fullfile(this_file_dir, 'T_out.xlsx');
if export_temperature_xlsx
    writematrix(temp_pi, Temperature_Result.temperature_field_file_in);
    writematrix(temp_po, Temperature_Result.temperature_field_file_out);
end
% 热物性台账：便于压力/温度耦合计算和后续采用现场实验数据进行替换。
Temperature_Result.fluid_thermal_properties = table(fluid_ids(:), fluid_name(:), rho_table(:), ...
    cp_table(:), k_table(:), 'VariableNames', {'fluid_id', 'fluid_name', 'rho_kg_m3', 'cp_J_kgK', 'k_W_mK'});
Temperature_Result.solid_thermal_properties = struct( ...
    'formation', struct('material', 'Representative sandstone', 'rho_kg_m3', rho_f, ...
        'cp_J_kgK', cp_f, 'k_W_mK', k_f), ...
    'casing', struct('material', 'Carbon steel', 'rho_kg_m3', rho_casing, ...
        'cp_J_kgK', cp_casing, 'k_W_mK', k_casing), ...
    'drill_pipe', struct('material', 'Carbon steel', 'rho_kg_m3', rho_drill_pipe, ...
        'cp_J_kgK', cp_drill_pipe, 'k_W_mK', k_drill_pipe), ...
    'model_note', ['Only casing conductivity is used in the present quasi-steady pipe-wall ', ...
        'resistance. Drill-pipe thermal storage needs an explicit drill-pipe geometry and temperature node.']);

figure('Name', 'HT1-004固井循环温度场', 'Color', 'w');
plot(T_yuan, depth, 'k--', 'LineWidth', 1.3); hold on;
plot_time_min = unique([0, 30, 60, 120, time_min(end)]);
plot_color = lines(numel(plot_time_min));
legend_text = cell(1, 1 + 2*numel(plot_time_min));
legend_text{1} = '原始地温';
for p = 1:numel(plot_time_min)
    [~, t_plot] = min(abs(time_min - plot_time_min(p)));
    plot(temp_pi(:, t_plot), depth, '--', 'Color', plot_color(p, :), 'LineWidth', 1.0);
    plot(temp_po(:, t_plot), depth, '-', 'Color', plot_color(p, :), 'LineWidth', 1.0);
    legend_text{2*p} = sprintf('管内 %.0f min', plot_time_min(p));
    legend_text{2*p + 1} = sprintf('环空 %.0f min', plot_time_min(p));
end
set(gca, 'YDir', 'reverse');
grid on;
xlabel('温度 / degC');
ylabel('井深 / m');
title('HT1-004固井循环温度场');
legend(legend_text, 'Location', 'best');

figure('Name', 'HT1-004环空井底温度', 'Color', 'w');
plot(time_min, bht, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('时间 / min');
ylabel('环空井底温度 / degC');
title('HT1-004环空井底温度随时间变化');

%% ==================== 不同深度位置温度随时间变化 ====================

%% ---------- 用户绘图参数设置区 ----------

% 需要绘制的井深位置，单位：m
% 可以任意增加或删除深度点，点数没有限制
plot_depths_m = [1000, 3000, 5000, 6500, 7200, 7600];

% 选择需要绘制的温度类型
plot_pipe_temperature      = true;    % 是否绘制管内温度
plot_annulus_temperature   = true;    % 是否绘制环空温度
plot_wall_temperature      = false;   % 是否绘制井壁温度
plot_formation_temperature = false;   % 是否绘制地层温度

% 是否将管内和环空温度分别绘制在两个窗口
separate_figure = true;

% 是否在命令行输出实际采用的计算节点
show_actual_depth = true;


%% ---------- 深度数据检查 ----------

% 转换成行向量，避免用户输入列向量时出现问题
plot_depths_m = plot_depths_m(:).';

if isempty(plot_depths_m)
    error('plot_depths_m为空，请至少输入一个绘图深度。');
end

if any(~isfinite(plot_depths_m))
    error('绘图深度中存在NaN或Inf，请检查plot_depths_m。');
end

% 检查是否存在超出计算井深范围的点
invalid_depth = plot_depths_m < min(depth) | ...
                plot_depths_m > max(depth);

if any(invalid_depth)
    fprintf('\n以下绘图深度超出了模型计算范围，将被删除：\n');
    fprintf('  %.2f m\n', plot_depths_m(invalid_depth));

    plot_depths_m(invalid_depth) = [];
end

if isempty(plot_depths_m)
    error('删除无效深度后，没有可以绘制的深度点。');
end


%% ---------- 寻找距离目标深度最近的计算节点 ----------

number_plot_depth = numel(plot_depths_m);

depth_index = zeros(1, number_plot_depth);
actual_depths_m = zeros(1, number_plot_depth);

for k = 1:number_plot_depth

    % 找到距离指定深度最近的网格节点
    [~, depth_index(k)] = min(abs(depth - plot_depths_m(k)));

    % 记录实际采用的计算节点深度
    actual_depths_m(k) = depth(depth_index(k));

end


%% ---------- 删除重复计算节点 ----------

% 当两个设定深度非常接近时，可能对应同一个计算节点
% 此处保留第一次出现的深度点
[depth_index, unique_position] = unique(depth_index, 'stable');

plot_depths_m = plot_depths_m(unique_position);
actual_depths_m = actual_depths_m(unique_position);

number_plot_depth = numel(depth_index);


%% ---------- 输出目标深度与实际节点 ----------

if show_actual_depth

    fprintf('\n=== 温度曲线绘图深度 ===\n');
    fprintf('序号       设定深度/m       实际节点深度/m       深度误差/m\n');

    for k = 1:number_plot_depth

        fprintf('%3d        %10.2f        %14.2f        %10.2f\n', ...
            k, ...
            plot_depths_m(k), ...
            actual_depths_m(k), ...
            actual_depths_m(k) - plot_depths_m(k));

    end

end


%% ---------- 绘图基础设置 ----------

plot_color = lines(number_plot_depth);

legend_depth = cell(1, number_plot_depth);

for k = 1:number_plot_depth
    legend_depth{k} = sprintf('%.1f m', actual_depths_m(k));
end


%% ============================================================
%  方式一：管内和环空分别绘制
% =============================================================

if separate_figure

    %% ---------- 管内温度 ----------

    if plot_pipe_temperature

        figure( ...
            'Name', '不同深度管内温度随时间变化', ...
            'Color', 'w', ...
            'Position', [180, 120, 900, 560]);

        hold on;

        for k = 1:number_plot_depth

            plot( ...
                time_min, ...
                temp_pi(depth_index(k), :).', ...
                'LineWidth', 1.5, ...
                'Color', plot_color(k, :));

        end

        grid on;
        box on;

        xlabel('时间 / min');
        ylabel('管内温度 / ℃');
        title('不同井深位置管内温度随时间变化');

        legend(legend_depth, ...
            'Location', 'best', ...
            'NumColumns', min(3, number_plot_depth));

        xlim([time_min(1), time_min(end)]);

    end


    %% ---------- 环空温度 ----------

    if plot_annulus_temperature

        figure( ...
            'Name', '不同深度环空温度随时间变化', ...
            'Color', 'w', ...
            'Position', [220, 150, 900, 560]);

        hold on;

        for k = 1:number_plot_depth

            plot( ...
                time_min, ...
                temp_po(depth_index(k), :).', ...
                'LineWidth', 1.5, ...
                'Color', plot_color(k, :));

        end

        grid on;
        box on;

        xlabel('时间 / min');
        ylabel('环空温度 / ℃');
        title('不同井深位置环空温度随时间变化');

        legend(legend_depth, ...
            'Location', 'best', ...
            'NumColumns', min(3, number_plot_depth));

        xlim([time_min(1), time_min(end)]);

    end


    %% ---------- 井壁温度 ----------

    if plot_wall_temperature

        figure( ...
            'Name', '不同深度井壁温度随时间变化', ...
            'Color', 'w', ...
            'Position', [260, 180, 900, 560]);

        hold on;

        for k = 1:number_plot_depth

            plot( ...
                time_min, ...
                temp_wall(depth_index(k), :).', ...
                'LineWidth', 1.5, ...
                'Color', plot_color(k, :));

        end

        grid on;
        box on;

        xlabel('时间 / min');
        ylabel('井壁温度 / ℃');
        title('不同井深位置井壁温度随时间变化');

        legend(legend_depth, ...
            'Location', 'best', ...
            'NumColumns', min(3, number_plot_depth));

        xlim([time_min(1), time_min(end)]);

    end


    %% ---------- 地层温度 ----------

    if plot_formation_temperature

        figure( ...
            'Name', '不同深度地层温度随时间变化', ...
            'Color', 'w', ...
            'Position', [300, 210, 900, 560]);

        hold on;

        for k = 1:number_plot_depth

            plot( ...
                time_min, ...
                temp_pf(depth_index(k), :).', ...
                'LineWidth', 1.5, ...
                'Color', plot_color(k, :));

        end

        grid on;
        box on;

        xlabel('时间 / min');
        ylabel('地层温度 / ℃');
        title('不同井深位置地层温度随时间变化');

        legend(legend_depth, ...
            'Location', 'best', ...
            'NumColumns', min(3, number_plot_depth));

        xlim([time_min(1), time_min(end)]);

    end


%% ============================================================
%  方式二：所有温度绘制在同一个窗口中
% =============================================================

else

    number_subplot = ...
        double(plot_pipe_temperature) + ...
        double(plot_annulus_temperature) + ...
        double(plot_wall_temperature) + ...
        double(plot_formation_temperature);

    if number_subplot == 0
        error('至少需要选择一种温度类型进行绘制。');
    end

    figure( ...
        'Name', '不同深度位置温度随时间变化', ...
        'Color', 'w', ...
        'Position', [160, 60, 950, 780]);

    subplot_number = 0;


    % 管内温度
    if plot_pipe_temperature

        subplot_number = subplot_number + 1;
        subplot(number_subplot, 1, subplot_number);
        hold on;

        for k = 1:number_plot_depth
            plot(time_min, temp_pi(depth_index(k), :).', ...
                'LineWidth', 1.4, ...
                'Color', plot_color(k, :));
        end

        grid on;
        box on;
        xlabel('时间 / min');
        ylabel('温度 / ℃');
        title('管内温度');
        legend(legend_depth, 'Location', 'best');
        xlim([time_min(1), time_min(end)]);

    end


    % 环空温度
    if plot_annulus_temperature

        subplot_number = subplot_number + 1;
        subplot(number_subplot, 1, subplot_number);
        hold on;

        for k = 1:number_plot_depth
            plot(time_min, temp_po(depth_index(k), :).', ...
                'LineWidth', 1.4, ...
                'Color', plot_color(k, :));
        end

        grid on;
        box on;
        xlabel('时间 / min');
        ylabel('温度 / ℃');
        title('环空温度');
        legend(legend_depth, 'Location', 'best');
        xlim([time_min(1), time_min(end)]);

    end


    % 井壁温度
    if plot_wall_temperature

        subplot_number = subplot_number + 1;
        subplot(number_subplot, 1, subplot_number);
        hold on;

        for k = 1:number_plot_depth
            plot(time_min, temp_wall(depth_index(k), :).', ...
                'LineWidth', 1.4, ...
                'Color', plot_color(k, :));
        end

        grid on;
        box on;
        xlabel('时间 / min');
        ylabel('温度 / ℃');
        title('井壁温度');
        legend(legend_depth, 'Location', 'best');
        xlim([time_min(1), time_min(end)]);

    end


    % 地层温度
    if plot_formation_temperature

        subplot_number = subplot_number + 1;
        subplot(number_subplot, 1, subplot_number);
        hold on;

        for k = 1:number_plot_depth
            plot(time_min, temp_pf(depth_index(k), :).', ...
                'LineWidth', 1.4, ...
                'Color', plot_color(k, :));
        end

        grid on;
        box on;
        xlabel('时间 / min');
        ylabel('温度 / ℃');
        title('地层温度');
        legend(legend_depth, 'Location', 'best');
        xlim([time_min(1), time_min(end)]);

    end

    sgtitle('不同井深位置温度随时间变化');

end


%% ==================== 本文件局部辅助函数 ====================
% 下列函数仅供本脚本调用：前两项定位流体界面，中间两项计算压降/换热，
% 后三项分别用于稀疏矩阵组装、施工阶段定位和流体物性查表。
function [tag_list, residual_volume_list, residual_height_list, vertical_residual_height_list] = ...
        local_annulus_interface(injected_volume_list, volume_from_bottom, annulus_area_dm2, inclination_deg)
    % 环空界面从井底向上推进；10000 表示尚未进入环空，-1 表示已越过环空顶部。
    n_time_local = numel(injected_volume_list);
    n_segment_local = numel(volume_from_bottom);
    tag_list = zeros(1, n_time_local);
    residual_volume_list = zeros(1, n_time_local);
    residual_height_list = zeros(1, n_time_local);
    vertical_residual_height_list = zeros(1, n_time_local);

    total_annulus_volume = volume_from_bottom(1);
    bottom_segment_volume = volume_from_bottom(end);
    for tt = 1:n_time_local
        injected_volume = injected_volume_list(tt);
        if injected_volume <= 0
            tag_list(tt) = 10000;
            residual_volume_list(tt) = 10000;
            residual_height_list(tt) = 10000;
            vertical_residual_height_list(tt) = 10000;
        elseif injected_volume >= total_annulus_volume
            tag_list(tt) = -1;
            residual_volume_list(tt) = -1;
            residual_height_list(tt) = -1;
            vertical_residual_height_list(tt) = -1;
        elseif injected_volume < bottom_segment_volume
            tag_list(tt) = n_segment_local;
            residual_volume_list(tt) = injected_volume;
            residual_height_list(tt) = injected_volume / max(annulus_area_dm2(n_segment_local), eps) / 10;
            vertical_residual_height_list(tt) = residual_height_list(tt) * cosd(inclination_deg(n_segment_local));
        else
            for ii = 2:n_segment_local
                if injected_volume >= volume_from_bottom(ii) && injected_volume < volume_from_bottom(ii - 1)
                    tag_list(tt) = ii - 1;
                    residual_volume_list(tt) = injected_volume - volume_from_bottom(ii);
                    residual_height_list(tt) = residual_volume_list(tt) / max(annulus_area_dm2(ii - 1), eps) / 10;
                    vertical_residual_height_list(tt) = residual_height_list(tt) * cosd(inclination_deg(ii - 1));
                    break;
                end
            end
        end
    end
end

function [tag_list, residual_volume_list, residual_height_list] = ...
        local_casing_interface(injected_volume_list, volume_from_top, pipe_area_m2)
    % 管内界面从井口向下推进；-1 表示尚未进入管内，10000 表示已到达管柱末端。
    n_time_local = numel(injected_volume_list);
    n_segment_local = numel(volume_from_top);
    tag_list = zeros(1, n_time_local);
    residual_volume_list = zeros(1, n_time_local);
    residual_height_list = zeros(1, n_time_local);

    total_pipe_volume = volume_from_top(end);
    first_segment_volume = volume_from_top(1);
    for tt = 1:n_time_local
        injected_volume = injected_volume_list(tt);
        if injected_volume <= 0
            tag_list(tt) = -1;
            residual_volume_list(tt) = -1;
            residual_height_list(tt) = -1;
        elseif injected_volume >= total_pipe_volume
            tag_list(tt) = 10000;
            residual_volume_list(tt) = 10000;
            residual_height_list(tt) = 10000;
        elseif injected_volume < first_segment_volume
            tag_list(tt) = 1;
            residual_volume_list(tt) = injected_volume;
            residual_height_list(tt) = injected_volume / max(pipe_area_m2(1) * 1000, eps);
        else
            for ii = 2:n_segment_local
                if injected_volume >= volume_from_top(ii - 1) && injected_volume < volume_from_top(ii)
                    tag_list(tt) = ii;
                    residual_volume_list(tt) = volume_from_top(ii) - injected_volume;
                    residual_height_list(tt) = residual_volume_list(tt) / max(pipe_area_m2(ii) * 1000, eps);
                    break;
                end
            end
        end
    end
end

function [dpi, h1] = local_caculate_pipe(Q, d_i, mu_p, tau_y, cp_m, k_m, rho_m, dz, area_pipe)
    % 管内宾汉近似：dpi 为单位长度压降（Pa/m），h1 为管内对流换热系数（W/(m^2*K)）。
    v_i = Q / max(area_pipe, eps);
    shear_rate = max(8 * v_i / max(d_i, eps), 1e-9);
    mu_eff = mu_p + tau_y / shear_rate;
    Pr = mu_eff * cp_m / max(k_m, eps);
    Re = rho_m * v_i * d_i / max(mu_eff, eps);

    if Re < 2300
        Nu = 3.65 + 0.0668 * d_i / max(dz, eps) * Re * Pr / ...
            (1 + 0.04 * (d_i / max(dz, eps) * Re * Pr)^0.67);
        dpi = 32 * mu_p * v_i / max(d_i^2, eps) + 4 * tau_y / max(d_i, eps);
    else
        Nu = 0.023 * Re^0.8 * Pr^0.4;
        f = 0.3164 / max(Re, eps)^0.25;
        dpi = f * rho_m * v_i^2 / (2 * max(d_i, eps));
    end
    h1 = Nu * k_m / max(d_i, eps);
end

function [dpo, h2, h3] = local_caculate_annuli(Q, d_hole, d_pipe, area_ann, mu_p, tau_y, cp_m, k_m, rho_m, dz)
    % 环空宾汉近似：dpo 为单位长度压降（Pa/m）；h2/h3 分别用于管柱侧和井壁侧换热。
    gap = max(d_hole - d_pipe, 1e-6);
    v_o = Q / max(area_ann, eps);
    shear_rate = max(12 * v_o / gap, 1e-9);
    mu_eff = mu_p + tau_y / shear_rate;
    Pr = mu_eff * cp_m / max(k_m, eps);
    Re = rho_m * v_o * gap / max(mu_eff, eps);

    if Re < 2300
        Nu = 3.65 + 0.0668 * gap / max(dz, eps) * Re * Pr / ...
            (1 + 0.04 * (gap / max(dz, eps) * Re * Pr)^0.67);
        dpo = 32 * mu_p * v_o / max(gap^2, eps) + 4 * tau_y / gap;
    else
        Nu = 0.023 * Re^0.8 * Pr^0.4;
        f = 0.3164 / max(Re, eps)^0.25;
        dpo = f * rho_m * v_o^2 / (2 * gap);
    end
    h_ann = Nu * k_m / gap;
    h2 = h_ann;
    h3 = h_ann;
end

function [row, col, val, nnz_count] = add_coeff(row, col, val, nnz_count, eq, variable_index, coeff)
    % 将单个系数写入稀疏矩阵三元组缓存。
    nnz_count = nnz_count + 1;
    row(nnz_count) = eq;
    col(nnz_count) = variable_index;
    val(nnz_count) = coeff;
end

function [stage_idx, cumulative_volume] = local_stage_and_volume(current_time, stage_start, stage_end, ...
        stage_start_volume, pump_values, stage_volume)
    % 根据当前时刻定位施工阶段，并返回从施工开始计的累计泵注体积（L）。
    stage_idx = find(current_time < stage_end - 1e-12, 1, 'first');
    if isempty(stage_idx)
        stage_idx = numel(stage_end);
        cumulative_volume = sum(stage_volume);
        return;
    end
    elapsed_in_stage = max(current_time - stage_start(stage_idx), 0);
    cumulative_volume = stage_start_volume(stage_idx) + elapsed_in_stage * pump_values(stage_idx);
    cumulative_volume = min(cumulative_volume, sum(stage_volume));
end

function fluid_id = local_fluid_id_from_volume(volume_marker, stage_volume, stage_ids)
    % 根据流体前锋累计体积确定当前位置的流体编号；前锋前方保持初始钻井液编号 0。
    if volume_marker <= 0
        fluid_id = 0;
        return;
    end
    cumulative_stage_volume = cumsum(stage_volume);
    stage_idx = find(volume_marker <= cumulative_stage_volume + 1e-9, 1, 'first');
    if isempty(stage_idx)
        stage_idx = numel(stage_ids);
    end
    fluid_id = stage_ids(stage_idx);
end

function [rho_value, mu_value, tau_value, cp_value, k_value] = local_fluid_property( ...
        fluid_id, fluid_ids, rho_table, mu_table, tau_table, cp_table, k_table)
    % 按流体编号返回密度、塑性黏度、屈服应力、定压比热和导热系数。
    fluid_loc = find(fluid_ids == fluid_id, 1, 'first');
    if isempty(fluid_loc)
        error('未找到流体编号 %.0f 对应的物性。', fluid_id);
    end
    rho_value = rho_table(fluid_loc);
    mu_value = mu_table(fluid_loc);
    tau_value = tau_table(fluid_loc);
    cp_value = cp_table(fluid_loc);
    k_value = k_table(fluid_loc);
end
