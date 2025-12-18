% Method 2 of decoding call identity from neural responses
% pseudo-population analysis. For time-warped data, could treat that all sparse MO neurons (~30) are recorded simultaneously 
% by randomly sampling renditions, then neural responses become a 30-dim vector. We can then ask if the activity at timepoint 
% x of call A is similar to activity at timepoint y of call B. Then we can build a similarity matrix of the neural activity, 
% which can be compared to the acoustics MMD matrix
% Zhilei 08/26/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


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
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePop']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 1. Read data
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

% read the MMD data
fd_mmd = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD', 'matrixCall', 'v4v5');
sigma = 0.7;
fn_mmd = fullfile(fd_mmd, sprintf('%s.%s%s.%.2f.mmd.mat', birdID, syls{1}, syls{2}, sigma));
load(fn_mmd);
mmd_matrix = mmd.mmd;
% what's the window size
win_size = 32; 
pad = win_size;
sec_per_frame = 0.001;  % how much sec per frame of data
% plot the diaganol
diag_mmd = diag(mmd_matrix);
rel_t_mmd = ((1:length(diag_mmd))-win_size/2 - pad) * sec_per_frame;
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

% calculate relative bin time
max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
rel_t = bin*(1:max_bin_num) - bin/2 - pad;


%% 3. Randomly sample renditions to form pseudo-population data
num_runs = 1000;
pop1 = nan(num_runs, size(neuron_comb,2), size(bin_sums1,2));
pop2 = nan(num_runs, size(neuron_comb,2), size(bin_sums2,2));
info = nan(num_runs, size(neuron_comb,2), 2);  % also record what renditions are sampled
rng(1992);
for ii=1:num_runs
  % loop through neurons, one rendition per neuron
  for ni=1:size(neuron_sort,1)
    nID = neuron_sort{ni};
    % for neuron 1
    idx1 = find(strcmp({seg1.neuronID}, nID));
    i1 = randsample(idx1, 1);
    pop1(ii, ni, :) = bin_sums1(i1,:);
    % for neuron 2
    idx2 = find(strcmp({seg2.neuronID}, nID));
    i2 = randsample(idx2, 1);
    pop2(ii, ni, :) = bin_sums2(i2,:);
    info(ii, ni, :) = [i1 i2];
  end
end

% plot some renditions to check
close all;
ds = {pop1; pop2};
for di=1:size(ds, 1)
  n_colors = 256;             % Total number of colors in the colormap
  jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
  custom_map = [0 0 0; jet_map];  % Prepend black to jet
  fig_pos = [10 10 1800 800];
  [fig, axes] = generatePanelGrid_v2(3, 10, [0.25;0.25;0.25], [0.05;0.05], [0.05;0.05], [0.05;0.05], 0.02, [1;1;1], fig_pos);
  ridx = randsample(1:size(d,1), 30);
  for ii=1:length(ridx)
    m = squeeze(d(ridx(ii), :, :)); 
    plot_i = floor((ii-1)/10)+1;
    plot_j = mod(ii-1, 10)+1;
    ax = axes(plot_i, plot_j);
    imagesc(ax, m, [0 5]);
    colormap(custom_map);
    title(ax, sprintf('%s %d', syls{di}, ridx(ii)), 'FontSize', 10);
  end
  fn_pdf = fullfile(fd_save, sprintf('%s.%s.pseudo.pdf', birdID, syls{di}));
  print(fig, fn_pdf, '-dpdf', '-painters');
end


%% 4. Direct comparison on the mean response matrix
mean1 = squeeze(mean(pop1, 1));
mean2 = squeeze(mean(pop2, 1));
% an absolute difference matrix
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


%% 4. Calculate neural distance: conventional metrics
% metric = 'euclidean';
metric = 'cosine';
% metric = 'correlation';
t = min([size(pop1,3) size(pop2,3)]);
dists = nan(num_runs, t);
for ri=1:num_runs
  for ti=1:t
    m1 = squeeze(pop1(ri,:,ti));
    m2 = squeeze(pop2(ri,:,ti));
    x = [m1; m2];
    dist_this = pdist(x, metric);
    dists(ri, ti) = dist_this;
  end
end
dists(isnan(dists)) = 0;
dist_mean = mean(dists, 1);
% similarity = 1./dist_mean;

% calculate relative time of bins, unit is time
bin_start = bin * (0:(length(dist_mean)-1));
bin_end = bin + bin_start;
bin_center = (bin_start+bin_end)/2;
rel_t = bin_center - pad;
% plot results
close all;
fig = ZZfunc_newFigurePDFsize_v1([50 50 600 600]);
% plot(rel_t, similarity);
plot(rel_t, dist_mean, 'LineWidth', 2);
xlabel('Relative time (sec)', 'FontSize', 14);
ylabel('Neural distance', 'FontSize', 14);
title(metric, 'FontSize', 14);
xlim([rel_t(1) 0.14]);


%% 5. Calculate neural distance: subtraction
t = min([size(pop1,3) size(pop2,3)]);
dists = zeros(num_runs, t);
for ri=1:num_runs
  for ti=1:t
    m1 = squeeze(pop1(ri,:,ti));
    m2 = squeeze(pop2(ri,:,ti));
    m1 = m1>0;
    m2 = m2>0;
%     x = [m1; m2];
%     dist_this = pdist(x, metric);
    dist_this = sum(abs(m1-m2));
    dists(ri, ti) = dist_this;
  end
end
% plot results
dist_mean = mean(dists, 1);
close all;
fig = ZZfunc_newFigurePDFsize_v1([50 50 600 600]);
% plot(rel_t, similarity);
plot(rel_t, dist_mean, 'LineWidth', 2);
xlabel('Relative time (sec)', 'FontSize', 14);
ylabel('Neural distance', 'FontSize', 14);
title('Subtraction', 'FontSize', 14);
xlim([rel_t(1) 0.14]);
% xlim([0 0.14]);







