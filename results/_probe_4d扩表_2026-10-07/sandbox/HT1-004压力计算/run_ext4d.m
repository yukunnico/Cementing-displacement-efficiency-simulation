% run_ext4d — Phase 4d 扩表运行包装：执行沙箱 HT1_004_T 并落盘关键标量。
% 注意：HT1_004_T 开头有 clear，故其后必须重新取自身路径。
try
    HT1_004_T;
catch err
    here = fileparts(mfilename('fullpath'));
    fid = fopen(fullfile(here, 'ext_scalars.txt'), 'w');
    fprintf(fid, 'RUN_ERROR=%s\n', err.message);
    fclose(fid);
    rethrow(err);
end
here = fileparts(mfilename('fullpath'));
fid = fopen(fullfile(here, 'ext_scalars.txt'), 'w');
fprintf(fid, 'n_depth=%d\n', numel(depth));
fprintf(fid, 'n_time=%d\n', n_time);
fprintf(fid, 'total_time_min=%.12f\n', total_time_min);
fprintf(fid, 'table_end_min=%.12f\n', time_min(end));
fprintf(fid, 'first200_axis_check=%.12f\n', time_min(200));  % 应为 198.8
fprintf(fid, 'T_in_end_end=%.10f\n', temp_pi(end,end));
fprintf(fid, 'T_out_end_end=%.10f\n', temp_po(end,end));
fprintf(fid, 'T_in_end_col200=%.10f\n', temp_pi(end,200));
fprintf(fid, 'T_out_end_col200=%.10f\n', temp_po(end,200));
fprintf(fid, 'Q_final_m3_s=%.6e\n', Q_time(end));
fclose(fid);
fprintf('EXT4D DONE n_time=%d end=%.1f min\n', n_time, time_min(end));
