% calculate firing metrics for sparse neurons
% using call subtypes, since spikes can be aligned
% the unit of analysis is a neuronID-call subtype pair
% Zhilei, 08/21/2025
% 1. Sparseness
% 2. IRCC

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
for bi=1:size(birdIDs, 2)
% bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Data inputs
% what neuron to analyze
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% where is the ephys data struct
fd_ephys = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedPull');
fns = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.pull.mat', birdID)));
% where to save metrics results
fd_save = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew', birdID, 'metricsSparse');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 2. Generate PSTH, calculate metrics
% loop through each call subtype, then loop each neuronID
% first time-warp, then generate PSTH, then 
% parameters for time-warping
align_pad = 0.1;  % how much to look before and after the syllable when align renditions, unit is seconds
fs = 20000;
% for calculating PSTH
psth_bin_size = 0.005;  % bin_size for calculating psth, unit is seconds
% for IRCC
ircc_include = [0.05 0];  % how much to include before syllable onset and after syllable offset when calculating IRCC
ircc_sigma = 0.01; % kernel width when calculating IRCC, unit is sec
ircc_bin_size = 0.01;  % temporal resolution for IRCC, unit is sec
% for sparseness
sparse_bin_size = 0.01;

% save results in a big data struct
metrics = []; 
count = 0;
for vi=1:size(fns,1)
  v = strsplit(fns(vi).name, '.'); 
  v = v{2};
  % save into a separate folder
  fd_save_this = fullfile(fd_save, v);
  if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
  % load ephys struct
  fn = fullfile(fns(vi).folder, fns(vi).name);
  load(fn);
  % only look at production now
  seg_selected = seg_selected(strcmp({seg_selected.aud_ch}, 'chan0'));
  
  % loop through neuron
  for ni=1:size(info, 1)
    neuronID = info.neuronID{ni};
    seg_this = seg_selected(strcmp({seg_selected.neuronID}, neuronID));
    if isempty(seg_this)
      continue;
    end
    
    % align the spikes and sounds
    [aligned_spike, aligned_sound, mean_dur] = ZZfunc_linearWarp_v2(seg_this, align_pad);
    
    % plot raster and PSTH with averaged spectrograms
    close all;
    fig_pos = [10 10 300 600];
    [fig, axes] = generatePanelGrid_v2(3, 1, [0.2, 0.4, 0.2], [0.05; 0.05], [0.05;0.02], [0.2;0.05], 0.02, [0;0;0], fig_pos);% loop through each cluster
    tick_size = 2; tick_color = 'black';  % control of raster ticks
    flim = [500 5500];
    [axes, psth, rt, psth_time] = ZZfunc_plotRasterWithSpectrogram_v3(axes, aligned_spike, aligned_sound, fs, align_pad, psth_bin_size, tick_size, tick_color, flim);
    % add number of renditions to the first panel
    text(axes(1), mean(axes(1).XLim), 5000, sprintf('n=%d',size(seg_this,2)), 'FontSize', 10, 'Color', 'white', ...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title(axes(1), sprintf('%s %s', v, neuronID), 'FontSize', 12);
    % save as pdf
    fn_pdf = fullfile(fd_save_this, sprintf('PSTH.%s.%s.%s.pdf', birdID, v, neuronID));
    print(fig, fn_pdf, '-dpdf', '-painters');
    
    % calculate IRCC
    % remove the portion that is unwanted when calculating IRCC
    istart = max([1 floor(fs*(align_pad-ircc_include(1)))]);
    iend = min([size(aligned_spike,2) size(aligned_spike,2)-floor(fs*(align_pad-ircc_include(2)))]);
    aligned_part = aligned_spike(:, istart:iend);
    ircc = ZZfunc_calcIRCC_v1(aligned_part, fs, ircc_sigma, ircc_bin_size);
    
    % calculate sparseness
    [psth_rate, time_vector] = ZZfunc_calcPSTH_v2(aligned_part, fs, sparse_bin_size);
%     figure; plot(time_vector, psth_rate);
    % normalize
    p = psth_rate / sum(psth_rate);
    % calculate sparseness index
    s = p .* log(p);
    sparseness = 1 + nansum(s) / log(size(p,2));
  
    % save results
    count = count+1;
    metrics(count).count = count;
    metrics(count).birdID = birdID;
    metrics(count).call_subtype = v;
    metrics(count).neuronID = neuronID;
    metrics(count).num_rends = size(seg_this,2);
    metrics(count).psth = psth;
    metrics(count).psth_rel_time = psth_time;
    metrics(count).ircc_istart = istart;
    metrics(count).ircc_iend = iend;
    metrics(count).ircc = ircc;
    metrics(count).sparseness = sparseness;
  end
end
  
% save metrics
fn_metrics = fullfile(fd_save, sprintf('%s.SparseNeuron.metrics.mat', birdID));
save(fn_metrics, 'metrics');
  
end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  












