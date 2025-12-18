% Method 3 of quantifying the difference on neural responses
% alculate difference on the mean firing rate for each neuron, generating two matrix for the two call subtypes. 
% Then compute an absolute difference matrix: x-axis is the relative time, y-axis is the neuronID. 
% Check what time range has the most difference. 
% Zhilei 08/27/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/decodeIdentityCall'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedPull');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.pull.mat', birdID)));
% where is the ephys data struct located,
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedPull');
% what syllable types to analyze
syls = {'v4', 'v5'};
% syls = {'v1', 'v7'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
fd_embed = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID);
suffix_hahn = 'neuron_orderedPlotted5';
suffix_criteria = 'criteria5';
% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_MeanRate']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 1. Read data
% ephys
e_all = cell(size(syls, 2), 1);  % raw ephys struct
n_all = cell(size(syls, 2), 1);  % selected and plotted neurons
c_all = cell(size(syls, 2), 1);  % criteria struct that includes burst info
% syl_i = 1;
for syl_i=1:size(syls, 2)
  ss = syls{syl_i};
  fn_e = fullfile(fd_ephys, sprintf('%s.%s.segments_all.pull.mat', birdID, ss));
  load(fn_e); e_all{syl_i} = seg_selected(strcmp({seg_selected.aud_ch}, 'chan0'));
  fn_n = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_hahn));
  load(fn_n); n_all{syl_i} = neuron_ordered;
  fn_c = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_criteria));
  load(fn_c); c_all{syl_i} = criteria;
end
clear seg_selected;
 
% acoustic MMD
fd_mmd = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD', 'matrixCall', 'v4v5');
sigma = 0.7;
fn_mmd = fullfile(fd_mmd, sprintf('%s.%s%s.%.2f.mmd.mat', birdID, syls{1}, syls{2}, sigma));
load(fn_mmd);
mmd_matrix = mmd.mmd;
% what's the window size for the acoustic MMD analysis
win_size = 32; 
sec_per_frame = 0.001;  % how much sec per frame of data
% plot the diagonal to check if match previous plot
diag_mmd = diag(mmd_matrix);
% time of each sliding window, relative to syllable onset
rel_t_mmd = ((1:length(diag_mmd))-1-win_size) * sec_per_frame;
close all; 
figure; plot(rel_t_mmd, diag_mmd, 'LineWidth', 2);
xlim([rel_t_mmd(1) rel_t_mmd(end)]);


%% 2. Time-warp renditions to each call subtype
% no cross warping
fs = 20000;
pad = 0.032; 
% what's the size of bins when counting spikes, unit is sec
bin = 0.005;
bin_pt = floor(fs*bin);
hop = 0.001; 
hop_pt = floor(fs*hop);

% only look at identified sparse neurons
neuron_comb = unique([n_all{:}], 'sorted');
% sort the neuron by the major burst time
n_sort = [];
c1 = c_all{1}; c2 = c_all{2};
for ni=1:size(neuron_comb, 2)
  n_sort(ni).neuronID = neuron_comb{ni};
  n_sort(ni).t1 = nan; 
  n_sort(ni).t2 = nan; 
  if ismember(neuron_comb{ni}, n_all{1}); n_sort(ni).t1=c1(strcmp({c1.neuronID}, neuron_comb{ni})).psth_max_smooth_t; end
  if ismember(neuron_comb{ni}, n_all{2}); n_sort(ni).t2=c2(strcmp({c2.neuronID}, neuron_comb{ni})).psth_max_smooth_t; end
  n_sort(ni).mean_t = nanmean([n_sort(ni).t1 n_sort(ni).t2]);
end
n_sort = sortrows(struct2table(n_sort), 'mean_t');
neuron_sort = n_sort.neuronID;

% align spikes
seg1 = e_all{1}(ismember({e_all{1}.neuronID}, neuron_sort));
seg2 = e_all{2}(ismember({e_all{2}.neuronID}, neuron_sort));
[aligned_spike1, aligned_sound1, mean_dur1] = ZZfunc_linearWarp_v2(seg1, pad);
[aligned_spike2, aligned_sound2, mean_dur2] = ZZfunc_linearWarp_v2(seg2, pad);

% count spikes in bins
bin_sums1 = ZZfunc_countSpikeBins_v1(aligned_spike1, bin_pt);
bin_sums2 = ZZfunc_countSpikeBins_v1(aligned_spike2, bin_pt);
% calculate relative bin start time
max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
rel_t = bin*(1:max_bin_num) - bin/2 - pad;

% % count spikes in sliding windows
% bin_sums1 = ZZfunc_countSpikeSlidingWin_v1(aligned_spike1, bin_pt, hop_pt);
% bin_sums2 = ZZfunc_countSpikeSlidingWin_v1(aligned_spike2, bin_pt, hop_pt);
% max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
% rel_t = hop*(1:max_bin_num) - hop - pad;



%% 3. Calculate mean firing rate matrix
mean1 = nan(size(neuron_sort,1), size(bin_sums1,2));
mean2 = nan(size(neuron_sort,1), size(bin_sums2,2));
for ni=1:size(neuron_sort,1)
  nID = neuron_sort{ni};
  % get renditions, then average
  idx1 = find(strcmp({seg1.neuronID}, nID));
  mean1(ni,:) = mean(bin_sums1(idx1,:), 1); 
  idx2 = find(strcmp({seg2.neuronID}, nID));
  mean2(ni,:) = mean(bin_sums2(idx2,:), 1); 
end



%% 4. Plot the mean firing rate
% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  
% truncate to the shortest
tmin = min([size(mean1, 2) size(mean2,2)]);
m_diff = abs(mean1(:, 1:tmin) - mean2(:, 1:tmin));
ds = {mean1; mean2; m_diff}; 
durs = [mean_dur1 mean_dur2 mean_dur1];
names = {syls{1}; syls{2}; 'Abs. diff.'};
close all;
[fig, axes] = generatePanelGrid_v2(1, 3, [0.7], [], [0.1;0.05], [0.15;0.05], 0.02, [0], [10 10 1400 500]);
clim_high = max([mean1(:); mean2(:)]);
for di=1:size(ds,1)
  d = ds{di};
  ax = axes(di);
  imagesc(ax, rel_t(1:size(d,2)), 1:size(d,1), d, [0 clim_high]); 
  colormap(custom_map); colorbar(ax);
  xline(ax, 0, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'white'); 
  xline(ax, durs(di)/fs, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'white'); 
  if di==1
    yticks(ax, 1:size(d,1)); 
    yticklabels(ax, neuron_sort);
  else
    yticks(ax, []);
  end
  title(ax, names{di}, 'FontSize', 14);
  xlabel(ax, 'Rel. time (sec)', 'FontSize', 12);
end
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.pseudoMean.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');

% plot a mean curve?
close all;
[fig2, axes] = generatePanelGrid_v2(2, 1, [0.35;0.35], [0.075], [0.1;0.05], [0.1;0.05], 0.05, [0;0], [10 10 600 1000]);
% first plot the MMD curve
plot(axes(1), rel_t_mmd, diag_mmd, 'LineWidth', 2);
xlim(axes(1), [rel_t_mmd(1) rel_t_mmd(end)]);
ylabel(axes(1), 'Acoustic MMD', 'FontSize', 12);
title(axes(1), sprintf('Acoustic distance (MMD): %s vs %s', syls{1}, syls{2}), 'FontSize', 14);
% then plot the neural response
mean_m_diff = mean(m_diff, 1);
plot(axes(2), rel_t(1:length(mean_m_diff)), mean_m_diff, 'LineWidth', 2);
xlim(axes(2), [rel_t(1) 0.15]);
xlabel(axes(2), 'Rel. time (sec)', 'FontSize', 12);
ylabel(axes(2), 'Mean neural difference', 'FontSize', 12);
title(axes(2), sprintf('Neural distance: %s vs %s', syls{1}, syls{2}), 'FontSize', 14);
linkaxes(axes, 'x');
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.AcousticVsNeural.pdf', birdID, syls{1}, syls{2}));
print(fig2, fn_pdf, '-dpdf', '-painters');


% save the mean firing rate matrix
fn_mat = fullfile(fd_save, sprintf('%s.%s_%s.AcousticVsNeural.mat', birdID, syls{1}, syls{2}));
acoustic_neural.diag_mmd = diag_mmd;
acoustic_neural.rel_t_mmd = rel_t_mmd;
acoustic_neural.mean_m_diff = mean_m_diff;
acoustic_neural.rel_t = rel_t; 
save(fn_mat, 'acoustic_neural');

















