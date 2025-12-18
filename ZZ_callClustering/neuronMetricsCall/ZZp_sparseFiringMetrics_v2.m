% calculate firing metrics for sparse neurons
% using call subtypes, since spikes can be aligned
% the unit of analysis is a neuronID-call subtype pair
% Zhilei, 08/22/2025
% 1. Sparseness
% 2. IRCC
% 3. Total renditions
% 4. Number of renditions that have spikes
% Differ from v1: no plotting of PSTH, add more quantifications

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
  calc_win = [0.05 0];  % how much to include before syllable onset and after syllable offset when calculating metrics
  % for calculating PSTH
  psth_bin_size = 0.005;  % bin_size for calculating psth, unit is seconds
  % for IRCC
  ircc_sigma = 0.01; % kernel width when calculating IRCC, unit is sec
  ircc_bin_size = 0.01;  % temporal resolution for IRCC, unit is sec
  % for sparseness
  sparse_bin_size = 0.01;
  
  % save results in a big data struct
  metrics = [];
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
    neuron_all = unique({seg_selected.neuronID}, 'sorted');
    metrics_this = [];
    % loop through neuron
    for ni=1:size(info, 1)
      neuronID = info.neuronID{ni};
      seg_this = seg_selected(strcmp({seg_selected.neuronID}, neuronID));
      metrics_this(ni).birdID = birdID;
      metrics_this(ni).call_subtype = v;
      metrics_this(ni).neuronID = neuronID;
      metrics_this(ni).num_rends = size(seg_this,2);
      if ~isempty(seg_this)
        % align the spikes and sounds
        [aligned_spike, aligned_sound, mean_dur] = ZZfunc_linearWarp_v2(seg_this, align_pad);
        
        % calculate PSTH
        % remove the portion that is unwanted
        istart = max([1 floor(fs*(align_pad-calc_win(1)))]);
        iend = min([size(aligned_spike,2) size(aligned_spike,2)-floor(fs*(align_pad-calc_win(2)))]);
        aligned_part = aligned_spike(:, istart:iend);
        [psth_rate, time_vector] = ZZfunc_calcPSTH_v2(aligned_part, fs, psth_bin_size);
        
        % calculate IRCC
        ircc = ZZfunc_calcIRCC_v1(aligned_part, fs, ircc_sigma, ircc_bin_size);
        
        % calculate sparseness
        % normalize
        [psth_rate2, time_vector2] = ZZfunc_calcPSTH_v2(aligned_part, fs, sparse_bin_size);
        p = psth_rate2 / sum(psth_rate2);
        % calculate entropy measure of sparseness
        s = p .* log(p);
        sparseness = 1 + nansum(s) / log(size(p,2));
        
        % how many renditions have spikes
        spike_counts = sum(aligned_part, 2);
        spike_rends = sum(spike_counts>0);
        spike_frac = spike_rends / length(spike_counts);
        
        % save results
        metrics_this(ni).psth = psth_rate;
        metrics_this(ni).psth_rel_time = time_vector;
        metrics_this(ni).ircc_istart = istart;
        metrics_this(ni).ircc_iend = iend;
        metrics_this(ni).ircc = ircc;
        metrics_this(ni).sparseness = sparseness;
        metrics_this(ni).spike_counts = spike_counts;
        metrics_this(ni).spike_rends = spike_rends;
        metrics_this(ni).spike_frac = spike_frac;
      end
    end
  
    % concatenate
    metrics = [metrics metrics_this];
  end

% save metrics
fn_metrics = fullfile(fd_save, sprintf('%s.SparseNeuron.metrics.mat', birdID));
save(fn_metrics, 'metrics');
  
end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  












