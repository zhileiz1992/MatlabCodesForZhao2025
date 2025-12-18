% remove rhd files that are after 23:15pm, to save disk space
clear; close all;

fd_home = '/mnt/z4/zz367/EphysMONAO/EphysDataNew/2024Ephys1/';
data_dates = {'2024-01-21'};

di = 1;
dd = data_dates{di};
folder_path = fullfile(fd_home, dd);

files = dir(fullfile(folder_path, '*.rhd'));
% cutoff_time = datetime('23:15', 'InputFormat', 'HH:mm');
% Compute the cutoff datetime: 23:15 of the date in file
cutoff_dt = datetime(dd, 'InputFormat', 'yyyy-MM-dd') + hours(23) + minutes(15);

for i = 1:length(files)
  fname = files(i).name;
  tokens = regexp(fname, '_(\d{6})_(\d{6})\.rhd$', 'tokens');
  if isempty(tokens)
    continue; % skip if filename doesn't match pattern
  end
  
  date_str = tokens{1}{1}; % e.g., '240120'
  time_str = tokens{1}{2}; % e.g., '105048'
  
  try
    file_dt = datetime([date_str time_str], 'InputFormat', 'yyMMddHHmmss');
    
    
    if file_dt > cutoff_dt
      fullpath = fullfile(folder_path, fname);
      delete(fullpath);
      fprintf('Deleted: %s\n', fname);
    end
  catch
    warning('Could not parse: %s', fname);
  end
end