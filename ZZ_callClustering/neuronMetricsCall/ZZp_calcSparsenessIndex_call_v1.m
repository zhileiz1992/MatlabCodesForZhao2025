% calculate the sparseness index of MO neurons for call syllables
% Zhilei, 07/10/2025
% method adopted from Goldberg & Fee 2010

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
% where to save results
fd_save = fullfile(fd_home, 'Figures', 'CombinedAnalysis', 'Sparseness');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end

% what bin size to use for calculating PSTH
psth_bin_size = 0.01;  % unit is seconds
pad = 0.05;  % how much to look before and after the syllable, unit is seconds

%% 1. loop through birds, calculate sparseness
% save data into a master table
d_all = cell(size(birdIDs,2), 1);
% bi = 1;
for bi=1:size(birdIDs, 2)
  birdID = birdIDs{bi};
  pairID = pairIDs{bi};
  
  
  %% 1. Load data struct, calculate PSTH, then sparseness index
  fd_save_master = fullfile(fd_home, 'Figures', pairID, 'HahnloserNew');
  fds = dir(fullfile(fd_save_master, birdID, 'extractedPull', '*.v*.segments_all.pull.mat'));
  % get the syllable labels
  syls = cellfun(@(s) strsplit(s, '.'), {fds.name}, 'UniformOutput', false);
  syls = cellfun(@(x) x{2}, syls, 'UniformOutput', false);
  syls = sort(syls);
  
  % loop through each call subtype
  d_this = struct();
  count = 0;
  for vi=1:size(fds,1)
    fprintf('%s %d\n', birdID, vi);
    % load the data struct
    fn_d = fullfile(fds(vi).folder, fds(vi).name);
    a = load(fn_d);
    syl_this = strsplit(fds(vi).name, '.');
    syl_this = syl_this{2};
    
    seg_selected = a.seg_selected;
    fs = seg_selected(1).fs;
    
    % load what neurons were selected for this call subtype
    fn_order = fullfile(fd_save_master, birdID, syl_this, sprintf('Hahnloser-%s-chan0.neuron_orderedPlotted5.mat', syl_this));
    b = load(fn_order);
    neuron_ordered = b.neuron_ordered;
    
    % loop through these neurons, calculate PSTH and sparness index
    for ni=1:size(neuron_ordered, 2)
      neuronID = neuron_ordered{ni};
      seg_this = seg_selected(strcmp({seg_selected.neuronID}, neuronID));
      % align the spikes and sounds
      [aligned_spike, ~, ~] = ZZfunc_linearWarp_v2(seg_this, pad);
      % calculate PSTH
      [psth_rate, time_vector] = ZZfunc_calcPSTH_v2(aligned_spike, fs, psth_bin_size);
      % normalize
      p = psth_rate / sum(psth_rate);
      % calculate sparseness index
      s = p .* log(p);
      s = 1 + nansum(s) / log(size(p,2));
      count = count + 1;
      d_this(count).syl = syl_this;
      d_this(count).neuronID = neuronID;
      d_this(count).sparseness = s;
    end
  end
  
  % calculate the mean sparseness of each neuron
  T = struct2table(d_this);
  % save
  fn_t = fullfile(fd_save, sprintf('%s.sparsenessT.mat', birdID));
  save(fn_t, 'T');
  % Compute mean ircc per neuronID
  result = groupsummary(T, 'neuronID', 'mean', 'sparseness');
  % save
  fn_r = fullfile(fd_save, sprintf('%s.sparsenessResult.mat', birdID));
  save(fn_r, 'result');
  
  d_all{bi} = result.mean_sparseness;
end


%% 2. Plot results as violin plot
all_data = [];
group_labels = [];
mean_list = zeros(numel(d_all),1);
for k = 1:numel(d_all)
    this_data = d_all{k}(:); % ensure column vector
    fprintf('%s mean sparseness: %.3f\n', pretty_ids{k}, mean(this_data));
    all_data = [all_data; this_data];
    group_labels = [group_labels; repmat(k, numel(this_data), 1)];
    mean_list(k) = mean(this_data);
end
fprintf('Overall mean sparseness: %.3f\n', mean(all_data));

% Create plot
close all; 
fig = ZZfunc_newFigurePDFsize_v1([50 50 300 400]); 
% boxplot(all_data, group_labels);
violinplot(all_data, group_labels, 'ViolinColor', [0.2 0.2 0.2], 'EdgeColor', [1 1 1], 'MarkerSize', 8, 'MedianMarkerSize', 100);
% add a sample size number to each violin and mean value to each 
for bi=1:size(d_all,1)
  text(bi, 0.15, sprintf('%d', length(d_all{bi})), 'FontSize', 8, 'HorizontalAlignment', 'center');
  text(bi, 0.1, sprintf('%.3f', mean_list(bi)), 'FontSize', 8, 'HorizontalAlignment', 'center');
end
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xlim([0.5 4.5]);
xticks([1 2 3 4]);
xticklabels(pretty_ids);
xlabel('Bird ID');
ylabel('Sparseness index');
title('Sparseness');

% save figure
fn_fig = fullfile(fd_save, 'sparseness.boxplot5.fig');
savefig(fn_fig);
fn_pdf = strrep(fn_fig, '.fig', '.pdf');
print(fig, fn_pdf, '-dpdf', '-painters');












