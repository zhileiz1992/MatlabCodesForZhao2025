% A pipeline script to automatically segment syllables from sorted dbases 
% Zhilei Zhao, 2025-05-23
% Input: list of dbases with one channel sorted in each dbase
% dbase is constructed by the ZZ_modern4 of electro_gui
% 1. Automatically annotate with WhisperSeg: call (v); chortle (b); click (h); squawk (e); all other (x)
% 2. Output a new dbase with WhisperSeg labels
% differ from v2: input a .mat table that contains sorted neuron information; loop through all sorted birds automatically;
% no division by 10 from nc files to audio

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase')); 
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils')); 

%% 0. Inputs
fd_z4 = '/mnt/z4'; 
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed'); 
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
% what folder has the dbase files
fd_master = fullfile(fd_home, 'DbaseFiles');
% where to save intermediate results
fd_save_itm= fullfile(fd_home, 'Figures');
% what WhisperSeg model to use
fd_ckpt = '/mnt/z4/zz367/WarbleAnalysis/Results/ModelBackup/20250514_MOsortedAll/final_checkpoint_ct2';
% add a suffix to the output dbase
suffix = 'Wsp1'; 


%% 1. Loop through each bird, segment
bi = 4; 
birdID = birdIDs{bi};
pairID = pairIDs{bi};
% load meta info
fn_meta = fullfile(fd_master, pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID)); 
load(fn_meta);

% Loop through each date, segment/annotate the warble with WhispserSeg 
% what date have sorted neurons
date_unique = unique(info.date);
for date_i=1:size(date_unique, 1)
% for date_i=25:29
  date_short = date_unique{date_i};
  data_date = strjoin({date_short(1:4), date_short(5:6), date_short(7:8)}, '-'); % add a dash in between 
  disp(data_date);
  % export sound channel as wav files in a temp folder
  fd_temp = fullfile(fd_home, 'tempWav', birdID, data_date); 
  if exist(fd_temp)
    rmdir(fd_temp, 's');
  end
  mkdir(fd_temp);
  % read in all sorted dbase
  sorted_chs = info{strcmp(info.date, date_short), 'channel'};
  d_all = {};
  for sc=1:length(sorted_chs)
    fn_base = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, sorted_chs{sc});
    fn_d = fullfile(fd_master, pairID, data_date, birdID, 'warble', fn_base);
    load(fn_d);
    d_all{sc} = dbase;
  end
  dbase.Fs = 20000; 
  % loop through each sound file. decide what channel to use as input for WhisperSeg
  fns_sound = dbase.SoundFiles;  % all dbase of the same date share same set of sound files
  p_logic = dbase.Properties;
  parfor si=1:length(fns_sound)
    p_sum = zeros(1, length(dbase.PropertyNames));
    for di=1:length(d_all)
      p_sum = p_sum + d_all{di}.Properties(si,:);
    end
    p_logic(si,:) = p_sum>0;
    if (p_sum(3)>0) || (p_sum(5)>0)  % warble production or call
      ch = 'chan0';
    elseif p_sum(4)>0   % auditory
      ch = 'chan17';
    elseif p_sum(6)>0   % playback copy
      ch = 'chan22';
    else
      continue;
    end
    % save as wav file
    fn_nc = fullfile(fns_sound(si).folder, strrep(fns_sound(si).name, 'chan0', ch));
    signal = ncread(fn_nc, 'data');
    fn_save = fullfile(fd_temp, strrep(fns_sound(si).name, '.nc', '.wav'));
%     audiowrite(fn_save, signal/10, dbase.Fs);  % need to divide the raw values by 10 
    audiowrite(fn_save, signal, dbase.Fs);
  end
  % WhisperSeg annotate
  seg_res = ZZfunc_runWhisperSeg_v1(fd_temp, fd_ckpt); 
  % flush into a new dbase and save
  dbase_seg = dbase;
  dbase_seg.Properties = p_logic;
  replace_str = '_chan0.wav';
  [dbase, dbaseOld] = ZZ_flushDbaseSegmentationFunc_v3(dbase_seg, seg_res, replace_str); 
  fn_temp = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, strrep(data_date,'-',''), suffix); 
  fn_seg = fullfile(fd_master, pairID, data_date, birdID, 'warble', fn_temp);
  save(fn_seg, 'dbase');
end

  
  