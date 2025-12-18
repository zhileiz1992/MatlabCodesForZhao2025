function [metrics] = ZZfunc_runXGBoostPython_v1(d_input, pred_name, output_name, group_name, n_fold_split, fd_save_res, fn_py, run_name)
% called python script (fn_py) to run XGBoost 

%% first split d_input into three CSV files: pred.csv, output.csv, group.csv
if ~exist(fd_save_res, 'dir'); mkdir(fd_save_res); end
d_pred = d_input(:, pred_name);
writetable(d_pred, fullfile(fd_save_res, sprintf('%s.predictor.csv', run_name)), 'WriteVariableNames', false, 'WriteRowNames', false);

d_output = d_input(:, output_name);
writetable(d_output, fullfile(fd_save_res, sprintf('%s.output.csv', run_name)), 'WriteVariableNames', false, 'WriteRowNames', false);

d_group = d_input(:, group_name);
writetable(d_group, fullfile(fd_save_res, sprintf('%s.group.csv', run_name)));


%% then replace the input/output in the template python script
lines = fileread(fn_py);
temp = strrep(lines, 'FD_SAVE_RES', fd_save_res);
temp = strrep(temp, 'RUN_NAME', run_name);
temp = strrep(temp, 'N_FOLD_SPLIT', n_fold_split);

% save as a temperory script
fn_py_temp = fullfile(fd_save_res, 'temp.py');
fid = fopen(fn_py_temp, 'w');
fwrite(fid, temp);
fclose(fid); 

% run the script
cmd = sprintf('source activate base && conda activate xgb && python %s', fn_py_temp); 
[status, cmdout] = system(cmd);

% parse the metrics
fn_metrics = fullfile(fd_save_res, sprintf('%s.metrics.csv', run_name));
if exist(fn_metrics, 'file')
  metrics = readtable(fn_metrics);
else
  metrics = [];
end

end

