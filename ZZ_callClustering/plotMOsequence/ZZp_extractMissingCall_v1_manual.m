% extract ephys data for the newly supplemented sorting
% avoid duplication with old sorting by choosing pProd=0
% Zhilei, 07/04/2025
% differ from v1: when not enough renditions, check the previous/next day for more renditions
% need to replace the neuronID at the end

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};


% loop through birds
bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};

% load the info regarding missing neuron-call subtype pairs
% fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_addMissingCall1.txt', birdID));
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_addManualCall1.txt', birdID));
info = readtable(fn_info, 'Delimiter', '\t', 'ReadVariableNames', false);
% load the info regarding old sorting as well, to get spike_shape info
fn_infoOld = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_infoOld); oldInfo=a.info;


%% 1. Loop through call subtype, extract ephys
% how many seconds to extract before and after the annotated syllable, need to consistent with old extracting
pad_pre = 0.1;
pad_post = 0.1;
unique_type = unique(info.Var2, 'stable');
for ci=1:size(unique_type, 1)
  syl = unique_type{ci};
  % loop through neurons in the table
  idx = find(strcmp(info.Var2, syl));
  segments_all = {};
  parfor rii=1:length(idx)
    segments_all{rii} = [];
    ri = idx(rii);
    % get date and channel id
    neuronID = info.Var1{ri};
    temp = strsplit(neuronID, '-');
    dd = temp{1};
    date_long = [dd(1:4) '-' dd(5:6) '-' dd(7:8)];
    ch = temp{2};
    ch_num = regexp(ch, '\d+', 'match');
    ch_num = str2num(ch_num{1});
    % note that segmentaion and ephys are in the same dbase
    fd_dbase = fullfile(fd_home, 'DbaseFiles', pairID, date_long, birdID, 'warble');
    fn_dbase = sprintf('%s.%s.warble.good.%s.addMiss_%s.dbase.mat', birdID, dd, ch, syl);
    fn = fullfile(fd_dbase, fn_dbase);
    if ~exist(fn, 'file')
      continue
    end
    a = load(fullfile(fd_dbase, fn_dbase));
    dbase = a.dbase;
    % get the spike shape information from old info table
    % that neuron doesn't exist in the old table, check what this neuron is supposed to replace
    neuronIDrep = info.Var3{ri};
    spike_shape = str2num(oldInfo{strcmp(oldInfo.neuronID, neuronIDrep), 'spike_shape'}{1});
%     spike_shape = str2num(oldInfo{strcmp(oldInfo.neuronID, neuronID), 'spike_shape'}{1});
    
    % extract ephys data
    segments = ZZ_extractSpikeForMissingFunc_v1(dbase, dbase, {syl}, ch_num, spike_shape, pad_pre, pad_post);
    if ~isempty(segments)
      [segments.data_date] = deal(date_long);  % add extra columns of date and neuron ID
%       [segments.neuronID] = deal(neuronID);
      [segments.neuronID] = deal(neuronIDrep);
      [segments.fn_dbase] = deal(fn_dbase);
    end
    segments_all{rii} = segments;
  end
  
  % Remove empty elements
  nonempty_segments = segments_all(~cellfun(@isempty, segments_all));
  % Concatenate all struct arrays
  segments_all = [nonempty_segments{:}];
  
  % save the results into disk for late use
  fd_save_this = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedAdded');
  if ~exist(fd_save_this, 'dir')
    mkdir(fd_save_this);
  end
  fn_save_seg = fullfile(fd_save_this, sprintf('%s.%s.segments_all.manual.mat', birdID, syl));
  save(fn_save_seg, 'segments_all', '-v7.3');
    
  % plot PSTH to check
%   pad = 0.06;
%   [aligned_spike, aligned_sound, mean_dur] = ZZfunc_linearWarp_v2(segments_all, pad);
%   fig_pos = [10 10 300 800];
%   [fig, axes] = generatePanelGrid_v2(3, 1, [0.2, 0.4, 0.2], [0.05; 0.05], [0.05;0.02], [0.2;0.05], 0.02, [0;0;0], fig_pos);% loop through each cluster
%   % psth_smooth = 200;  % widow size when smooth psth, 200 data point is 10ms
%   bin_size = 0.005;  % size of bin when calcualting psth in unit of sec
%   tick_size = 4; tick_color = 'black';  % control of raster ticks
%   [axes, psth, rt] = ZZfunc_plotRasterWithSpectrogram_v2(axes, aligned_spike, aligned_sound, fs, pad, bin_size, tick_size, tick_color);
end

  
  
  
  
  
  
  
  
  
  
  
  
  
  