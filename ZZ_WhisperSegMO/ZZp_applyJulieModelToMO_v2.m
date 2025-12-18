% a script to apply trained WhisperSeg model to MO ephys audio dataset
% Zhilei Zhao, 05/14/2025
% differ from v1: save dbases in dedicated folders

close all; clear;
addpath(genpath('/home/zz367/LabSoftware/ElectroGui'));
addpath(genpath('/home/zz367/LabSoftware/MATLAB-utils'));
cd('/home/zz367/LabSoftware/ElectroGui/source');

%% 1. Input
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'WarbleAnalysis');
% where wav data is stored
fd_wav_base = fullfile(fd_home, 'DataNew', '20250512_MOsorted');
birdIDs = {'pair5RigCCU29', 'pair5Rig0CU55', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
rep_str_all = {'_chan0.wav', '_chan17.wav', '_chan0.wav', '_chan0.wav', '_chan0.wav'};
% what WhisperSeg model to use
fd_ckpt = '/mnt/z4/zz367/WarbleAnalysis/Results/WhisperSeg/20250514_MOsorted2/final_checkpoint_ct2/';
% where to save results
fd_save = fullfile(fd_wav_base, 'dbase3');
if ~exist(fd_save)
  mkdir(fd_save);
end


%% 2. Apply WhisperSeg model
% bi = 1; 
for bi=1:size(birdIDs,2)
  bd = birdIDs{bi};
  fd_wav = fullfile(fd_wav_base, bd);
  % create a dbase for the wav files
  dbase = electro_gui.CreateDbase(defaults_ZZwav2, fd_wav);
  % change to windows path
  dbase.PathName = ZZ_linuxPathToWin_v1(dbase.PathName, 'Y:', '\mnt\z4');
  dbase.Fs = 20000;
  % save dbase
  fn_save_this = fullfile(fd_save, sprintf('%s.original.dbase.mat', bd));
  save(fn_save_this, 'dbase');
  
  % apply WhisperSeg model
  seg_res = ZZfunc_runWhisperSeg_v1(fd_wav, fd_ckpt);
  % flush into a new dbase and save
  replace_str = rep_str_all{bi};
  [dbase, dbaseOld] = ZZ_flushDbaseSegmentationFunc_v3_wav(dbase, seg_res, replace_str, replace_str);
  fn_seg = fullfile(fd_save, sprintf('%s.applyWsp.dbase.mat', bd));
  save(fn_seg, 'dbase');
end