% Method 2 of decoding call identity from neural responses
% pseudo-population analysis. For time-warped data, could treat that all sparse MO neurons (~30) are recorded simultaneously 
% by randomly sampling renditions, then neural responses become a 30-dim vector. We can then ask if the activity at timepoint 
% x of call A is similar to activity at timepoint y of call B. Then we can build a similarity matrix of the neural activity, 
% which can be compared to the acoustics MMD matrix
% Zhilei 08/28/2025
% differ from v1: use SVM to decode

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
sec_per_frame = 0.001;  % how much sec per frame of data
% plot the diaganol
diag_mmd = diag(mmd_matrix);
% relative time of the sliding window start
rel_t_mmd = ((1:length(diag_mmd)) -1 - win_size) * sec_per_frame;
close all; 
figure; plot(rel_t_mmd, diag_mmd, 'LineWidth', 2);
xlim([rel_t_mmd(1) rel_t_mmd(end)]);



%% 2. Time-warp renditions to each call subtype
% no cross warping
fs = 20000;
% pad = 0.032; 
pad = 0.08;
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

% calculate relative time of the bin start
max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
rel_t = bin*(1:(max_bin_num-1)) - pad;


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
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet
close all;
ds = {pop1; pop2};
for di=1:size(ds, 1)
  d = ds{di};
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

% also save the sample renditions for later use
fn_pop = fullfile(fd_save, sprintf('%s.%s.pop.mat', birdID, syls{1}));
save(fn_pop, 'pop1');
fn_pop = fullfile(fd_save, sprintf('%s.%s.pop.mat', birdID, syls{2}));
save(fn_pop, 'pop2');



%% 4. Build SVM models to decode call identity: time aligned (only look at diagonal)
min_bin = min([size(pop1,3) size(pop2,3)]);
% how many train/test shuffling to do
num_runs = 50;
train_ratio = 0.5;
kernel = 'linear';
% loop through time points
accu_all = nan(min_bin, 2, num_runs);
for ti=1:min_bin
  X1 = squeeze(pop1(:, :, ti));
  X2 = squeeze(pop2(:, :, ti));
  Y1 = repmat({syls{1}}, size(X1,1), 1);
  Y2 = repmat({syls{2}}, size(X2,1), 1);
  r_seed = 1992;
  [accu_real, accu_shuffle] = ZZfunc_trainTestSVM_v1(X1, X2, Y1, Y2, train_ratio, kernel, num_runs, r_seed);
  
  accu_all(ti,1,:) = accu_real;
  accu_all(ti,2,:) = accu_shuffle;

 % plot results as violin plots
%   close all;
%   fig = ZZfunc_newFigurePDFsize_v1([10 10 300 600]);
%   d = [accu_real; accu_shuffle];
%   categ = [repmat({'Real'}, length(accu_real), 1); repmat({'Shuffle'}, length(accu_shuffle), 1)];
%   violinplot(d, categ, 'ViolinColor', [0.3 0.68 0.29; 0.2 0.2 0.2], 'EdgeColor', [1 1 1], 'MarkerSize', 12, 'MedianMarkerSize', 200);
%   yline(0.5, 'LineStyle', '--', 'LineWidth', 1);
%   xlim([0.5 2.5]);
%   ylim([0 1]);
end
% save results
fn_accu = fullfile(fd_save, sprintf('%s.%s_%s.accu_all.mat', birdID, syls{1}, syls{2}));
save(fn_accu, 'accu_all');

% load
load(fn_accu);
% calculate then plot mean accuracy
accu_mean = mean(accu_all, 3);
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 800 600]);
ax = gca; hold(ax, 'on');
plot(rel_t(1:min_bin), accu_mean(:,1), 'Color', '#4daf4a', 'LineStyle', '-', 'LineWidth', 1); 
hold on;
plot(rel_t(1:min_bin), accu_mean(:,2), 'Color', '#737373', 'LineStyle', '-', 'LineWidth', 1); 
xline(0, 'LineStyle', '--');
xline(mean_dur1/fs, 'LineStyle', '--');
xlim([rel_t(1) rel_t(end)]);
ylim([0.4 1]);
xlabel('Rel. time (sec)', 'FontSize', 14);
ylabel('SVM decoding accuracy', 'FontSize', 14);
title(sprintf('Decoding accuracy for %s vs %s', syls{1}, syls{2}),'FontSize', 14); 
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.meanAccu.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');


% plot along with the MMD in two panels
close all;
[fig, axes] = generatePanelGrid_v2(2, 1, [0.35;0.35], [0.075], [0.1;0.05], [0.1;0.05], 0.05, [0;0], [10 10 600 1000]);
% determine the x-axis limit
xlim1 = max([rel_t_mmd(1) rel_t(1)]);
xlim2 = min([rel_t_mmd(end) rel_t(end)]);
% first plot the MMD curve
plot(axes(1), rel_t_mmd, diag_mmd, 'LineWidth', 2);
xlim(axes(1), [xlim1 xlim2]);
ylabel(axes(1), 'MMD on VAE latents', 'FontSize', 14);
title(axes(1), sprintf('Acoustic distance (MMD): %s vs %s', syls{1}, syls{2}), 'FontSize', 14);
% then SVM accuracy
ax2 = axes(2);
% plot(axes(2), rel_t(1:min_bin), accu_mean(:,1), 'Color', '#4daf4a', 'LineStyle', '-', 'LineWidth', 1); 
% plot(axes(2), rel_t(1:min_bin), accu_mean(:,2), 'Color', '#737373', 'LineStyle', '-', 'LineWidth', 1); 
plot_mean_with_ci(ax2, rel_t(1:min_bin), squeeze(accu_all(:,1,:))', [0.3 0.7 0.3], 1, [0.3 0.7 0.3], 0.5);
plot_mean_with_ci(ax2, rel_t(1:min_bin), squeeze(accu_all(:,2,:))', [0.5 0.5 0.5], 1, [0.5 0.5 0.5], 0.5);
% xlim(axes(2), [rel_t(1) rel_t(end)]);
xlim(ax2, [xlim1 xlim2]);
ylim(ax2, [0.45 1]);
xlabel(ax2, 'Rel. time (sec)', 'FontSize', 14);
ylabel(ax2, 'SVM decoding accuracy', 'FontSize', 14);
title(ax2, sprintf('Decoding accuracy for %s vs %s', syls{1}, syls{2}),'FontSize', 14); 
linkaxes(axes, 'x');
ax2 = ZZfunc_addSimpleLegend_v2(ax2, {'Real', 'Shuffled'}, {'#4daf4a', '#737373'}, 16, [-0.02 -0.02], [0.95 0.9]);
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.MMDvsSVM.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');


% plot in the same panel
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 650 600]);
ax = gca; hold(ax, 'on');
xlim1 = max([rel_t_mmd(1) rel_t(1)]);
xlim2 = min([rel_t_mmd(end) rel_t(end)]);
yyaxis left
plot_mean_with_ci(ax, rel_t(1:min_bin), squeeze(accu_all(:,1,:))', [0.3 0.7 0.3], 2, [0.3 0.7 0.3], 0.5);
plot_mean_with_ci(ax, rel_t(1:min_bin), squeeze(accu_all(:,2,:))', [0.5 0.5 0.5], 2, [0.5 0.5 0.5], 0.5);
ylabel(ax, 'SVM decoding accuracy', 'FontSize', 14);
ax = ZZfunc_addSimpleLegend_v2(ax, {'Neural: real', 'Neural: shuffled', 'Acoustic'}, {'#4daf4a', '#737373', '#377eb8'}, 16, [-0.02 -0.02 -0.02], [0.95 0.92 0.89]);

yyaxis right
plot(ax, rel_t_mmd, diag_mmd, 'LineWidth', 2, 'Color', '#377eb8');
xlim(ax, [xlim1 xlim2]);
ylabel(ax, 'MMD on VAE latents', 'FontSize', 14);
xlabel(ax, 'Rel. time (sec)', 'FontSize', 14);

ax.YAxis(1).Color = [0.3 0.7 0.3];
ax.YAxis(2).Color = '#377eb8';
title(ax, 'Neural vs acoustic distance', 'FontSize', 16);
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.MMDvsSVM_comb.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 5. Perform cross-correlation
x = accu_mean(:,1);
t_x = rel_t(1:size(x, 1)); 
y = diag_mmd;
t_y = rel_t_mmd;
[corr, time_lags] = cross_corr_plot(x, t_x, y, t_y);
% plot
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 250 600]);
plot(time_lags, corr, 'Color', '#377eb8', 'LineStyle', '-', 'LineWidth', 2);
xline(0, 'Color', '#737373', 'LineStyle', '--', 'LineWidth', 1);
xlabel('Time lags (sec)', 'FontSize', 14);
% ylabel('Normalized cross correlation', 'FontSize', 14);
xlim([-0.07 0.07]);
title('Cross-correlation', 'FontSize', 14);
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.crossCorrelation.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');










