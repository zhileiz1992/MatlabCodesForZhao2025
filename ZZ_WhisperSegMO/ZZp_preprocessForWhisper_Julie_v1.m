% a script to preprocess Julie's mirror dataset annotation 
% to use for train a WhisperSeg model to segment MO ephys data
% Zhilei, 05/08/2025

close all; clear; 

%% 1. Input
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'WarbleAnalysis'); 
% where is Julie's original annotation
fd_data = fullfile(fd_home, 'DataNew', '20240815_JulieMirrorAnnoNew'); 
birdIDs = {'B', 'BG', 'G', 'R'};


%% 2. Check basic information
bi = 1; 
bd = birdIDs{bi};
% what's the original sampling rate
fns_wav = dir(fullfile(fd_data, bd, '*.wav')); 
[signal,fs] = audioread(fullfile(fns_wav(1).folder, fns_wav(1).name));

% how much data in total
fns_wav_all = dir(fullfile(fd_data, '*', '*.wav')); 
fprintf('Total duration: %.2f mins\n', length(fns_wav_all)*15/60);  % each file is 15 seconds

% how many labels per category 
labels_all = {};
for bi=1:length(birdIDs)
  bd = birdIDs{bi};
  fn_dbase = dir(fullfile(fd_data, bd, '*.ZZ.dbase.mat')); 
  load(fullfile(fn_dbase.folder, fn_dbase.name));
  temp = [dbase.SegmentTitles{:}];
  labels_all = [labels_all temp];
end
[chars, ~, idx] = unique(labels_all);
counts = histc(idx, 1:numel(chars));
count_table = table(chars', counts);


%% 3. Downsample audio to a new folder
for bi=1:length(birdIDs)
  bd = birdIDs{bi};
  fd_save = fullfile(fd_home, 'DataNew', '20240815_JulieMirrorAnno20k', bd);
  if ~exist(fd_save)
    mkdir(fd_save);
  end
  fns_wav = dir(fullfile(fd_data, bd, '*.wav')); 
  for fi=1:length(fns_wav)
    [y, fs] = audioread(fullfile(fns_wav(fi).folder, fns_wav(fi).name));  % fs = 50000
    y_ds = resample(y, 2, 5);               % Downsample from 50k to 20k
    fs_ds = 20000;
    fn = fullfile(fd_save, fns_wav(fi).name);
    audiowrite(fn, y_ds, fs_ds);
  end
end


%% 4. Convert the dbase as well
% final dbase birdID.20k.repLabel.dbase.mat
% change labels
label_ori = {'v', 'h', 'e', 'b', 'x', 'f', 'g', 'm', 'p', 's'};
label_new = {'v', 'h', 'e', 'b', 'x', 'x', 'x', 'x', 'x', 'x'};
map = containers.Map(label_ori, label_new);
for bi=1:length(birdIDs)
  bd = birdIDs{bi};
  fn_dbase = dir(fullfile(fd_data, bd, '*.ZZ.dbase.mat')); 
  load(fullfile(fn_dbase.folder, fn_dbase.name));
  dbase2 = dbase;
  dbase2.Fs = 20000;
  fd_wav = fullfile(fd_home, 'DataNew', '20240815_JulieMirrorAnno20k', bd);
  dbase2.PathName = ZZ_linuxPathToWin_v1(fd_wav, 'Y:', '\mnt\z4'); % use Windows path so can be open in eletro_gui
  % change the sounds files from .avi to .wav
  f = dbase2.SoundFiles;
  for si=1:size(f,1)
    f(si).name = strrep(f(si).name, '.avi', '.wav');
  end
  dbase2.SoundFiles = f;
  % change the segment time, save a copy of dbase
  t = dbase2.SegmentTimes;
  t_new = t; 
  for si=1:size(t, 2)
    if ~isempty(t{si})
      t_new{si}= floor(t{si}/2.5);
    end
  end
  dbase2.SegmentTimes = t_new;
  dbase = dbase2; 
  fn_dbase = fullfile(fd_wav, sprintf('%s.20k.oriLabel.dbase.mat', bd));
  save(fn_dbase, 'dbase');
  % convert the labels
  C = dbase.SegmentTitles;
  C_new = cellfun(@(inner) cellfun(@(ch) map(ch), inner, 'UniformOutput', false), C, 'UniformOutput', false);
  dbase.SegmentTitles = C_new;
  fn_dbase = fullfile(fd_wav, sprintf('%s.20k.repLabel.dbase.mat', bd));
  save(fn_dbase, 'dbase');
end


%% 5. Manual inspection and correction
% dupliate birdID.20k.repLabel.dbase.mat to birdID.20k.zzLabel.dbase.mat
% if neccesary modify the labels 


%% 6. How much data after correction?
% how many labels per category 
labels_all = {};
fns_count = 0;
label_inc = {'v', 'b', 'h', 'e', 'x'};
for bi=1:length(birdIDs)
  bd = birdIDs{bi};
  fd_save = fullfile(fd_home, 'DataNew', '20240815_JulieMirrorAnno20k', bd);
  fn_dbase = fullfile(fd_save, sprintf('%s.20k.zzLabel.dbase.mat', bd)); 
  load(fn_dbase);
  temp = [dbase.SegmentTitles{:}];
  labels_all = [labels_all temp];
  fns_count = fns_count + size(dbase.SoundFiles,1);
  % what files have unwanted labels?
  for si=1:size(dbase.SoundFiles,1)
    labels_this = dbase.SegmentTitles{si};
    a = setdiff(labels_this, label_inc);
    if ~isempty(a)
      fprintf('%s.%d.%s\n', bd, si, a{1});
    end
  end
end
fprintf('Total duration: %.2f mins\n', fns_count*15/60);  
[chars, ~, idx] = unique(labels_all);
counts = histc(idx, 1:numel(chars));
count_table = table(chars', counts);
% calculate proportion
count_table.proportion = count_table.counts / sum(count_table.counts);
disp(count_table)













