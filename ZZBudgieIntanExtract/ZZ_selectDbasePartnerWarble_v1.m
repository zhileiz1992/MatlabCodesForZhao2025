% A script to select subsets of files in a dbase
% focus on the warble episodes of the partner bird
% to examine the auditory responeses
% Zhilei Zhao, 02/26/2024

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/'));

%% Input settings
% Folder setting
fd_z4 = '/mnt/z4';
fd_base_dbase = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'DbaseFiles');
% Dataset setting
expID = 'pair1';
birdID = 'pair1CU21RigB'; % in the order of audio channels
data_dates = {'2023-11-04'};

%% loop through each data dates
% what dbase file to load
for di=[1]
  fn_dbase = dir(fullfile(fd_base_dbase, expID, data_dates{di}, birdID, 'warble', sprintf('%s.*.warble.segFlatness.dbase.mat', birdID)));
  load(fullfile(fn_dbase.folder, fn_dbase.name));
  % what warble episodes to select
  fn_select = dir(fullfile(fd_base_dbase, expID, data_dates{di}, birdID, 'warble', sprintf('%s.*.warble.partner.txt', birdID))); 
  lines = readlines(fullfile(fn_select.folder, fn_select.name));
  to_select = [];
  for li=1:length(lines)
    ln_rep = strrep(lines{li}, ',', ' ');
    array = eval(['[', ln_rep, ']']);
    to_select = [to_select array];
  end
  % determine those warble episode corresponds to what dbase files
  fns_all = {dbase.SoundFiles(:).name};
  select_idx = [];
  for fi=1:length(fns_all)
    % get the warble id
    wid = strsplit(fns_all{fi}, '_');
    wid = str2num(wid{2});
    if ismember(wid, to_select)
      select_idx = [select_idx fi];
    end
  end
  % subset the dbase
  [dbase_subset, dbase_old] = ZZ_subsetDbaseFunc_v1(dbase, select_idx);
  % need to swap the sound files: ch0 <-> ch17
  for fi=1:length(dbase_subset.SoundFiles)
    dbase_subset.SoundFiles(fi).name = strrep(dbase_subset.SoundFiles(fi).name, 'chan0', 'chan17'); 
  end
  for fi=1:length(dbase_subset.ChannelFiles{1,17})
    dbase_subset.ChannelFiles{1,17}(fi).name = strrep(dbase_subset.ChannelFiles{1,17}(fi).name, 'chan17', 'chan0');
  end
  % save as a new dbase
  dbase = dbase_subset;
  fn_new = strrep(fn_dbase.name, 'segFlatness', 'partner');
  save(fullfile(fn_dbase.folder, fn_new), 'dbase');
end
  
    
  
  