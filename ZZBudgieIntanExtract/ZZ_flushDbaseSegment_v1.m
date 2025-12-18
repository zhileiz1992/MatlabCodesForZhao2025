% A script to flush syllable segmentation into dbase
% Zhilei Zhao, 02/06/2024

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/'));
rmpath(genpath('/mnt/z4/zz367/LabSoftware/ElectroGui-main'));

%% Input settings
% Folder setting
fd_z4 = '/mnt/z4';
fd_base_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
% Dataset setting
expID = 'pair1';
birdID = 'pair1CU21RigB'; % in the order of audio channels
data_dates = {'2023-09-12','2023-11-04', '2023-11-05', '2023-11-08', '2023-11-12'};
% what type of segmentation to flush
seg_type = 'segFlatness';  % segmentation based on flatness produced by Han's script

%% Segment based on flatness, then flush into dbase
for di=[1]
  % read dbase
  fn_dbase = dir(fullfile(fd_base_dbase, expID, data_dates{di}, birdID, 'warble', '*warble.start.dbase.mat'));
  load(fullfile(fn_dbase.folder, fn_dbase.name));
  % segmentation based on spectral flatness
  % loop through each file, find onsets and offsets of syllables
  dbase_path = ZZ_winPathToLinux_v1(dbase.PathName, 'Y:', fd_z4);
  fns = fullfile(dbase_path, {dbase.SoundFiles(:).name});
  % flush the segmentation into dbase
  % parameters to define syllables
  param.maskFrequency = 1000;  %ignore lower freq when calculate amplitude
  param.ampIgnore = -7; %ignore where amplitude is very small
  % threshold to identify peaks
  param.minDuration = 0.02;  %unit is sec, squawks in warble is quite short
  param.maxDuration = 10;   %in case there is long element
  param.minInterval = 0;  % minimal interval between two syllables
  % threshold to identify peaks: this is setup-specific and depends on the
  % quality of recording, use the 'troubleshoot.m' function to plot flatness
  % and check if the threshold makes sense
  param.thresholdFlatness = -0.6;
  param.extendFlatness = -0.75;
  param.gapSize = 5;  % two syllables with gaps smaller than this will be merged
  [onsets, offsets, labels] = ZZ_seg_flatness_NC_v2(fns, param);
  
  [dbase, dbaseOld] = ZZ_flushDbaseSegmentationFunc_v2(dbase, onsets, offsets, labels);
  fn_save_dbase = strrep(fn_dbase.name, 'start', seg_type);
  save(fullfile(fn_dbase.folder, fn_save_dbase), 'dbase');
end

