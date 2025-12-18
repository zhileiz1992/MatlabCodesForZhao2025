% A pipeline script to automatically segment syllables from sorted dbases 
% Zhilei Zhao, 2025-01-09
% Input: list of dbases with one channel sorted in each dbase
% dbase is constructed by the ZZ_modern4 of electro_gui
% 1. Automatically annotate with WhisperSeg: call (v); chortle (b); click (h); squawk (e); g syllable; all other (x)
% 2. Output a new dbase with WhisperSeg segmentaion labeled 

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase')); 
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils')); 

%% 0. Inputs
fd_z4 = '/mnt/z4'; 
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed'); 
pairID = 'pair5CU29CU55'; 
birdID = 'pair5RigCCU29'; 
data_date = '2024-09-12'; 
% what folder has the dbase files
fd_master = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble');
% what channels have been sorted, read channel and spike shape info from txt files
fn_spike_shape = fullfile(fd_master, sprintf('%s.%s.spikeShape.txt', birdID, strrep(data_date,'-','')));
info = readtable(fn_spike_shape);
sorted_chs = info{:,1}';
spike_shapes = info{:,2}';  % 1 means top spike, 2 means bottom spike
% where to save intermediate results and plots
fd_save_base = fullfile(fd_home, 'Figures', pairID, data_date, birdID);
if ~exist(fd_save_base)
  mkdir(fd_save_base);
end


%% 1. Segment/Annotate the warble with WhispserSeg 
% what WhisperSeg model to use
fd_ckpt = '/mnt/z4/zz367/WarbleAnalysis/Results/ModelBackup/20240907_ZhileiEphysWarble20kFinetune3/final_checkpoint_ct2';
% export sound channel as wav files in a temp folder
fd_temp = fullfile(fd_home, 'tempWav', birdID, data_date); 
if exist(fd_temp)
  rmdir(fd_temp, 's');
end
mkdir(fd_temp);
% read in all the sorted dbase
d_all = {};
for sc=1:length(sorted_chs)
  fn_d = fullfile(fd_master, sprintf('%s.%s.warble.good.ch%d.dbase.mat', birdID, strrep(data_date,'-',''), sorted_chs(sc)));
  load(fn_d);
  d_all{sc} = dbase;
end
% loop through each sound file. decide what channel to use as input for WhisperSeg
fns_sound = dbase.SoundFiles;
p_logic = dbase.Properties;
for si=1:length(fns_sound)
  p_sum = zeros(1, length(dbase.PropertyNames));
  for di=1:length(d_all)
    p_sum = p_sum + d_all{di}.Properties(si,:);
  end
  p_logic(si,:) = p_sum>0;
  if (p_sum(3)>0) || (p_sum(5)>0)  % warble production or call
    ch = 'chan0';
  elseif p_sum(4)>0   % auditory
    ch = 'chan17';
  elseif p_sum(6)>0   % playback copye
    ch = 'chan22';
  else
    continue;
  end
  % save as wav file
  fn_nc = fullfile(fns_sound(si).folder, strrep(fns_sound(si).name, 'chan0', ch));
  signal = ncread(fn_nc, 'data');
  fn_save = fullfile(fd_temp, strrep(fns_sound(si).name, '.nc', '.wav'));
  audiowrite(fn_save, signal/10, dbase.Fs);  % need to divide the raw values by 10 
end
% WhisperSeg annotate
seg_res = ZZfunc_runWhisperSeg_v1(fd_temp, fd_ckpt); 
% flush into a new dbase and save
dbase_seg = dbase;
dbase_seg.Properties = p_logic;
replace_str = '_chan0.wav';
[dbase, dbaseOld] = ZZ_flushDbaseSegmentationFunc_v3(dbase_seg, seg_res, replace_str); 
fn_seg = fullfile(fd_master, sprintf('%s.%s.warble.good.Wsp.dbase.mat', birdID, strrep(data_date,'-','')));
save(fn_seg, 'dbase');

  
  