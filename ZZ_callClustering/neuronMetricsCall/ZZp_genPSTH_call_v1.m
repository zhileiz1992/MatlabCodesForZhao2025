% script to generate PSTH (peri-stimulus time historgram) plots for all sorted neurons
% Zhilei, 07/07/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_callClustering'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};

% loop through birds
bi = 3;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load information about sorted neurons and call subtypes
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% get unique dates
date_unique = unique(info.date_long);
% where to save intermediate results and plots
fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
% grab what calls have been extracted
fds = dir(fullfile(fd_save_master, birdID, 'extractedPull', '*segments_all.pull.mat'));


%% 2. Generate PSTH for each call type-neuron pair
for li=1:size(fds,1)
  % li = 1;
  % load previous saved data struct
  temp = strsplit(fds(li).name, '.');
  syl_label = temp{2};
  fprintf('Loading data for %s %s...\n', birdID, syl_label);
  load(fullfile(fds(li).folder, fds(li).name));
  fs = seg_selected(1).fs;
  
  % loop through neuron
  for ni=1:size(info,1)
    neuronID = info.neuronID{ni};
    seg_this = seg_selected(strcmp({seg_selected.neuronID}, neuronID));
    if isempty(seg_this)
      continue;
    end
    % time-wrap spike train and sound
    pad = 0.06; % how much pre and post data to include, unit is sec
    [aligned_spike, aligned_sound, mean_dur] = ZZfunc_linearWarp_v2(seg_this, pad);
    
    % calculate IRCC to add to the title
    ircc_sigma = 0.01; % kernel width when calculating IRCC, unit is sec
    ircc_bin_size = 0.01;  % temporal resolution for IRCC, unit is sec
    ircc = ZZfunc_calcIRCC_v1(aligned_spike, fs, ircc_sigma, ircc_bin_size);
    
    % plot raster and PSTH with averaged spectrograms
    close all;
    fig_pos = [10 10 300 800];
    [fig, axes] = generatePanelGrid_v2(3, 1, [0.2, 0.4, 0.2], [0.05; 0.05], [0.05;0.02], [0.2;0.05], 0.02, [0;0;0], fig_pos);% loop through each cluster
    % psth_smooth = 200;  % widow size when smooth psth, 200 data point is 10ms
    bin_size = 0.005;  % size of bin when calcualting psth in unit of sec
    tick_size = 4; tick_color = 'black';  % control of raster ticks
    [axes, psth, rt] = ZZfunc_plotRasterWithSpectrogram_v2(axes, aligned_spike, aligned_sound, fs, pad, bin_size, tick_size, tick_color);
    title(axes(1), sprintf('%s %s %.3f', syl_label, neuronID, ircc), 'FontSize', 10);
    % save as pdf
    fd_save_this = fullfile(fd_save_master, birdID, syl_label, 'PSTH');
    if ~exist(fd_save_this, 'dir')
      mkdir(fd_save_this);
    end
    fn_pdf = fullfile(fd_save_this, sprintf('%s.%s.%s.PSTH.pdf', birdID, syl_label, neuronID));
    print(fig, fn_pdf, '-dpdf', '-painters');
    
    % optional: plot to check IFRs of all renditions
    %   figure;
    %   imagesc(ifrs);
    %   colormap hot;
    %   set(gca, 'YDir', 'Normal');
    
  end
end



















