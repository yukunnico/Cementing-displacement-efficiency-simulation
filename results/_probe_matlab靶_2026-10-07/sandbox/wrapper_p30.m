% Phase 3.0 MATLAB 靶探针包装器（2026-10-07）
% 登记：Phase 3.0 先行闸门（总纲 §3 Phase 3.0 / 续作 §7）——p_jaifang1.m 输出=绘图+fprintf
% 不落 CSV，本包装在 -batch 里补捕获（writematrix/save/diary），**原脚本零改动**。
% 靶复现不了 ⇒ 停下先报（§8-3），产物仍保留供取证。
here = fileparts(mfilename('fullpath'));
cd(here);
diary('matlab_batch_diary.txt');
fprintf('=== wrapper_p30 start %s ===\n', char(datetime('now')));
% ⚠ p_jaifang1.m 开头有 clear;clc ⇒ run 之前的一切 base 工作区变量都会被清掉
% （首跑教训：t0=tic 被 clear 掉导致捕获段未执行）。此处不放任何前置变量。
try
    run('p_jaifang1.m');
catch ME
    fprintf(2, 'P30_RUN_ERROR: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf(2, '  at %s (%s:%d)\n', ME.stack(k).name, ME.stack(k).file, ME.stack(k).line);
    end
    diary off;
    exit(2);
end
fprintf('=== p_jaifang1 done %s ===\n', char(datetime('now')));

% ---- 全工作区 + 关键矩阵落盘（存在才写，防路径差异哑火）----
save('out_workspace_full.mat');
if exist('Result_Matrix_Volume_Pressure', 'var')
    writematrix(Result_Matrix_Volume_Pressure, 'out_result_matrix_volume_pressure.csv');
end
if exist('Result_Matrix_PumpPressure_Comparison', 'var')
    writematrix(Result_Matrix_PumpPressure_Comparison, 'out_pump_pressure_comparison.csv');
end
if exist('Result_Matrix_BottomECD_Comparison', 'var')
    writematrix(Result_Matrix_BottomECD_Comparison, 'out_bottom_ecd_comparison.csv');
end
if exist('required_backpressure_MPa', 'var')
    writematrix([time_minutes_real(:), required_backpressure_MPa(:), ...
                 BP_lower_history(:), BP_upper_history(:), conflict_flag(:)], ...
                'out_four_point_backpressure.csv');
end
if exist('ECD_5568_new', 'var')
    writematrix([time_minutes_real(:), ECD_5568_new(:), ECD_6880_new(:), ...
                 ECD_7463_new(:), ECD_bottom_new(:)], 'out_four_point_ecd_new.csv');
end
if exist('pressure_pump_surface', 'var')
    writematrix([time_minutes_real(:), pressure_pump_surface(:)], 'out_pump_pressure_surface.csv');
end
if exist('pressure_annuli_MPa', 'var')
    writematrix([time_minutes_real(:), pressure_annuli_MPa(:, end)], 'out_annuli_bottom_pressure.csv');
end
if exist('ECD_casing_g_cm3', 'var')
    writematrix(ECD_casing_g_cm3, 'out_ecd_casing_full.csv');
end
if exist('P_wh_ctrl_win', 'var')
    writematrix([time_minutes_real(:), ECD_base_bottom_win(:), ECD_bottom_ctrl_win(:), ...
                 ECD_shoe_ctrl_win(:), P_wh_ctrl_win(:)], 'out_window_backpressure_ctrl.csv');
end

% ---- 标量摘要（Python 侧对账用）----
fid = fopen('out_summary.txt', 'w');
names = {'n_time', 'n_segment', 'TVD_bottom', ...
         'safe_ECD_bottom_lower', 'safe_ECD_bottom_upper', ...
         'safe_ECD_5568_lower', 'safe_ECD_5568_upper', ...
         'safe_ECD_6880_lower', 'safe_ECD_6880_upper', ...
         'safe_ECD_7463_lower', 'safe_ECD_7463_upper'};
for k = 1:numel(names)
    if exist(names{k}, 'var')
        v = eval(names{k});
        if isscalar(v) && isnumeric(v)
            fprintf(fid, '%s = %.10g\n', names{k}, v);
        end
    end
end
if exist('pressure_pump_surface', 'var')
    fprintf(fid, 'max_pump_pressure_MPa = %.10g\n', max(pressure_pump_surface(:)));
    fprintf(fid, 'min_pump_pressure_MPa = %.10g\n', min(pressure_pump_surface(:)));
end
if exist('conflict_flag', 'var')
    fprintf(fid, 'conflict_steps = %d / %d\n', sum(conflict_flag(:) ~= 0), numel(conflict_flag));
end
if exist('pressure_annuli_MPa', 'var')
    fprintf(fid, 'max_annuli_bottom_MPa = %.10g\n', max(pressure_annuli_MPa(:, end)));
end
if exist('ECD_casing_g_cm3', 'var')
    fprintf(fid, 'max_bottom_ECD_g_cm3 = %.10g\n', max(ECD_casing_g_cm3(:, end)));
    fprintf(fid, 'min_bottom_ECD_g_cm3 = %.10g\n', min(ECD_casing_g_cm3(:, end)));
end
fclose(fid);

whos
fprintf('=== wrapper_p30 done ===\n');
diary off;
exit(0);
