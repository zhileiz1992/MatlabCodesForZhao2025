% script to calculate IRCC (inter-rendition correlation coeffient) for call subtypes
% Zhilei, 06/26/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_callClustering'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% where call embedding results are stored
fd_embed_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what call subtype dbase to use
suffix = 'Wsp2Call';


% loop through birds
bi = 2;
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
fds = dir(fullfile(fd_save_master, birdID, 'extracted', '*segments_all.mat'));


%% 2. Calculate IRCC for each call subtype and neuron
for li=1:size(fds,1)
  % li = 1;
  % load previous saved data struct
  temp = strsplit(fds(li).name, '.');
  syl_label = temp{2};
  fprintf('Loading data for %s %s...\n', birdID, syl_label);
  load(fullfile(fds(li).folder, fds(li).name));
  % for production choose chan0, for auditory choose chan17
  pt = 'chan0';
  % choose specific dates OR choose all dates
  seg_selected = segments_all(strcmp({segments_all.aud_ch}, pt));
  fs = seg_selected(1).fs;
  
  % loop through neuron
  %   ni = 40;
  ircc_all = zeros(size(info,1),1);
  for ni=1:size(info,1)
    neuronID = info.neuronID{ni};
    seg_this = seg_selected(strcmp({seg_selected.neuronID}, neuronID));
    if isempty(seg_this)
      continue;
    end
    % time-wrap spike train and sound
    pad = 0.06; % how much pre and post data to include, unit is sec
    [aligned_spike, aligned_sound, mean_dur] = ZZfunc_linearWarp_v2(seg_this, pad);
    
    % calculate instantaneous firing rate using Gaussian convolution and then IRCC
    sigma = 0.01; % 10 ms kernel width for each spike, Caleb used 20ms
    binsz = 0.01; % 10 ms bin for temporal resolution
    ifrs = [];
    for ri=1:size(aligned_spike,1)
      ifrs(ri,:) = ZZfunc_IFRwithGaussian_v1(aligned_spike(ri,:), fs, sigma, binsz);
    end
    % calcualte IRCC
    R = corrcoef(ifrs');
    %   flat = reshape(R,[1,size(R,2)^2]);
    %   cc = mean(flat(~isnan(flat)));
    % Get logical index for off-diagonal elements
    off_diag_idx = ~eye(size(R));  % logical matrix: true everywhere except diagonal
    % Extract off-diagonal elements
    off_diag_elements = R(off_diag_idx);
    % Remove NaNs
    valid_elements = off_diag_elements(~isnan(off_diag_elements));
    % Compute the mean
    cc = mean(valid_elements);
    ircc_all(ni) = cc;
    
    % optional: plot raster and PSTH with averaged spectrograms to check
    close all;
    fig_pos = [10 10 300 800];
    [fig, axes] = generatePanelGrid_v2(3, 1, [0.2, 0.4, 0.2], [0.05; 0.05], [0.05;0.02], [0.2;0.05], 0.02, [0;0;0], fig_pos);% loop through each cluster
    % psth_smooth = 200;  % widow size when smooth psth, 200 data point is 10ms
    bin_size = 0.005;  % size of bin when calcualting psth in unit of sec
    tick_size = 4; tick_color = 'black';  % control of raster ticks
    [axes, psth, rt] = ZZfunc_plotRasterWithSpectrogram_v2(axes, aligned_spike, aligned_sound, fs, pad, bin_size, tick_size, tick_color);
    title(axes(1), sprintf('%s %.3f', neuronID, cc), 'FontSize', 12);
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



















