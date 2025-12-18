% preprocess the expert segmentation of songbird datasets
% such that segmentation are stored in the same folder as wav files
% XX.time.txt file: syllable onset and offset in unit of data points, seperated by comma
% XX.label.txt file: syllable labels
% Zhilei, 05/30/2025

close all; clear;


%% folder setting
fd_z4 = '/mnt/z4';
fd_data = fullfile(fd_z4, 'zz367', 'WarbleAnalysis', 'Data', 'TweetyNetData');


%% 1. Bengalese finch
% flatten the folders of dates
birdFn = dir(fullfile(fd_data, 'Bengalese', '*'));
for ii=5:size(birdFn,1)
  if strcmp(birdFn(ii).name, '.') || strcmp(birdFn(ii).name, '..')
    continue
  end
  fd = fullfile(birdFn(ii).folder, birdFn(ii).name);
  subfolders = dir(fd);
  subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));
  for k = 1:length(subfolders)
    subfd = fullfile(fd, subfolders(k).name);
    % Get all files in this subfolder
    files = dir(fullfile(subfd, '*'));
    files = files(~[files.isdir]);  % ignore directories
    % Move each file to root folder
    for j = 1:length(files)
      src = fullfile(subfd, files(j).name);
      dest = fullfile(fd, files(j).name);
      movefile(src, dest);
    end
  end
end

% convert expert segmentation to standard format
birdIDs = birdFn([birdFn.isdir] & ~ismember({birdFn.name}, {'.', '..'}));
birdIDs = {birdIDs.name};
for bi=1:size(birdIDs,2)
  fns_wav = dir(fullfile(fd_data, 'Bengalese', birdIDs{bi}, '*.wav'));
  % determine fs
  [signal, fs] = audioread(fullfile(fns_wav(1).folder, fns_wav(1).name));
  disp(fs);
  parfor fi=1:size(fns_wav,1)
    fn = fullfile(fns_wav(fi).folder, [fns_wav(fi).name '.csv']);
    if ~exist(fn, 'file')
      continue;
    end
    % Read the CSV file
    opts = detectImportOptions(fn);
    opts.SelectedVariableNames = {'onset_s', 'offset_s', 'label'};
    opts = setvartype(opts, 'label', 'string');
    T = readtable(fn, opts);
    % Convert times to sample indices
    onset_samples = round(T.onset_s * fs);
    offset_samples = round(T.offset_s * fs);
    % Write fn.time.txt
    time_out = [onset_samples, offset_samples];
    fn_time = strrep(fn, '.wav.csv', '.time.txt');
    writematrix(time_out, fn_time, 'Delimiter', ',');
    % Write fn.label.txt
    labels = T.label;
    fn_label =  strrep(fn, '.wav.csv', '.label.txt');
    fid = fopen(fn_label, 'w', 'n', 'UTF-8');
    for i = 1:numel(labels)
      fprintf(fid, '%s\r\n', labels{i});
    end
    fclose(fid);
  end
end


%% 2. Canary
birdIDs = {'llb3', 'llb11', 'llb16'};
for bi=2:size(birdIDs,2)
  bd = birdIDs{bi};
  % Input paths
  csv_file = fullfile(fd_data, 'Canary', sprintf('%s_annot.csv', bd));
  audio_folder = fullfile(fd_data, 'Canary', sprintf('%s_songs', bd));
  % Read CSV with relevant columns
  opts = detectImportOptions(csv_file, 'NumHeaderLines', 0);
  opts.VariableNamesLine = 1;    % header on first line
  opts.DataLines = [2, Inf];     % data starts from second line
  % Set only needed variables (replace with actual names if needed)
  vars_needed = {'label', 'onset_s', 'offset_s', 'audio_file'};
  opts.SelectedVariableNames = vars_needed;
  opts = setvartype(opts, vars_needed, 'char');
  T = readtable(csv_file, opts);
  
  % Convert label to integer and then to char (0 → 'a', ..., 25 → 'z')
  % Convert numeric label to single-letter char:
  % 0–25 → 'a'–'z', 26–51 → 'A'–'Z'
  label_int = str2double(T.label);  % or use label_str if defined
  if any(label_int < 0 | label_int > 51)
    error('Label values must be between 0 and 51 to map to a-zA-Z.');
  end
  label_char = repmat('?', size(label_int));  % placeholder
  low_mask = label_int <= 25;
  high_mask = label_int > 25;
  label_char(low_mask)  = char(label_int(low_mask)  + double('a'));
  label_char(high_mask) = char(label_int(high_mask) - 26 + double('A'));
  % Get unique audio files
  audio_files = unique(T.audio_file);
  % Process each audio file
  for i = 1:length(audio_files)
    fname = audio_files{i};
    wav_path = fullfile(audio_folder, fname);
    
    if ~isfile(wav_path)
      warning('Audio file not found: %s', wav_path);
      continue;
    end
    
    % Read sample rate
    info = audioinfo(wav_path);
    fs = info.SampleRate;
    
    % Get rows for this file
    idx = strcmp(T.audio_file, fname);
    if isempty(idx)
      continue;
    end
    onset_s = str2double(T.onset_s(idx));
    offset_s = str2double(T.offset_s(idx));
    onset_samples = floor(onset_s * fs);
    offset_samples = floor(offset_s * fs);
    labels = label_char(idx);
    times = [onset_samples, offset_samples];
    
    % Prepare output filenames (strip .wav if present)
    [base, ~] = strtok(fname, '.');
    label_out = fullfile(audio_folder, [base '.label.txt']);
    time_out = fullfile(audio_folder, [base '.time.txt']);
    
    % Write label file
    fid_label = fopen(label_out, 'w', 'n', 'UTF-8');
    for j = 1:numel(labels)
      fprintf(fid_label, '%s\r\n', labels(j));
    end
    fclose(fid_label);
    
    % Write time file
    writematrix(times, time_out, 'Delimiter', ',');
  end
end











