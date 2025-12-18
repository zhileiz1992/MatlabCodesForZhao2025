% export wav files and segmentations to a temporary folder
% for convenient network training and access
% differ from v1: for birds that don't have many sorted neurons, export unsorted files and run WhisperSeg to segment

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils')); 
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase')); 


%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
fd_data_base = fullfile(fd_home, 'ExtractedIntan');
fd_save_base = fullfile(fd_home, 'vaeWav');
% what bird to extract
birdID = 'pair2RigBCU25';
pairID = 'pair2CU20CU25';
% birdID = 'pair4RigBCU53';
% pairID = 'pair4CU68CU53';
% what WhisperSeg model to use
fd_ckpt = '/mnt/z4/zz367/WarbleAnalysis/Results/ModelBackup/20250514_MOsortedAll/final_checkpoint_ct2';
% where to save file
fd_save_this = fullfile(fd_save_base, birdID, 'Audio', 'unsorted');
disp(fd_save_this);
if exist(fd_save_this, 'dir')
  rmdir(fd_save_this);
end
mkdir(fd_save_this);


%% check the warbleWav folder for each unsorte date, manually pick a few files per date
% save info in a txt, read in here
fn_info = fullfile(fd_save_base, birdID, 'unsorted_select.txt');
opts = detectImportOptions(fn_info, 'Delimiter', '\t', 'TextType', 'char');
opts.VariableTypes(:) = {'char'};
info = readtable(fn_info, opts);


%% Loop through dates, export wav, then whisperseg
fs = 20000; 
for li=1:size(info,1)
  dd = info.Var1{li};
  dd_long = [dd(1:4) '-' dd(5:6) '-' dd(7:8)];
  % what warble episodes to export
  widx = strsplit(info.Var2{li}, ',');
  % locate the folder that stores NC files
  fd_nc = fullfile(fd_data_base, pairID, dd_long, birdID, 'NCFiles', 'warble');
  for wi=1:size(widx,2)
    fns = dir(fullfile(fd_nc, sprintf('warble_%05d_*_chan0.nc', str2num(widx{wi}))));
    for fi=1:size(fns,1)
      fn_nc = fullfile(fns(fi).folder, fns(fi).name);
      signal = ncread(fn_nc, 'data');
      fn_save = fullfile(fd_save_this, strrep(fns(fi).name, '.nc', '.wav'));
      audiowrite(fn_save, signal, fs);
    end
  end
end
  
% how many wav files 
temp = dir(fullfile(fd_save_base, 'pair5RigCCU29', 'Audio', '*', '*.wav'));
temp2 = dir(fullfile(fd_save_base, 'pair2RigBCU25', 'Audio', '*', '*.wav'));


%% Segment with WhipserSeg
seg_res = ZZfunc_runWhisperSeg_v1(fd_save_this, fd_ckpt); 
% save segmentation as txt files
fn_unique = unique(seg_res.audio_path, 'stable');
for fi=1:size(fn_unique,1)
  fn = fn_unique{fi};
  seg_this = seg_res(strcmp(seg_res.audio_path, fn),:);
  seg_t = table2array(seg_this(:,1:2));
  seg_t = floor(seg_t * fs);
  fn_t = fullfile(fd_save_this, strrep(fn, '.wav', '.time.txt'));
  writematrix(seg_t, fn_t, 'Delimiter', ',');
  seg_n = table2array(seg_this(:,3));
  fn_n =fullfile(fd_save_this, strrep(fn, '.wav', '.label.txt'));
  fid = fopen(fn_n, 'w');
  for i = 1:numel(seg_n)
    fprintf(fid, '%s\n', seg_n{i});
  end
  fclose(fid);
end

% how many 'v' syllable?
a = find(strcmp(seg_res.label, 'v'));






  
  
  
  
  
  
  
  
  
  
  
  
  

