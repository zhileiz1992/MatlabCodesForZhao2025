function [metrics] = ZZfunc_checkXGBoostOnTest_v1(model_path, fn_test_predictor, fn_test_output, fn_save, fn_py, fn_temp_py)

lines = fileread(fn_py);
temp = strrep(lines, 'MODEL_PATH', model_path);
temp = strrep(temp, 'FN_TEST_PREDICTOR', fn_test_predictor);
temp = strrep(temp, 'FN_TEST_OUTPUT', fn_test_output);
temp = strrep(temp, 'FN_SAVE', fn_save);

% save as a temperory script
fid = fopen(fn_temp_py, 'w');
fwrite(fid, temp);
fclose(fid); 

% run the WhispserSeg script
cmd = sprintf('source activate base && conda activate xgb && python %s', fn_temp_py); 
[status, cmdout] = system(cmd);

% parse the metrics
if exist(fn_save, 'file')
  metrics = readtable(fn_save);
else
  metrics = [];
end

end

