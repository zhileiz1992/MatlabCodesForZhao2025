% Method 2 of decoding call identity from neural responses
% pseudo-population analysis. For time-warped data, could treat that all sparse MO neurons (~30) are recorded simultaneously 
% by randomly sampling renditions, then neural responses become a 30-dim vector. We can then ask if the activity at timepoint 
% x of call A is similar to activity at timepoint y of call B. Then we can build a similarity matrix of the neural activity, 
% which can be compared to the acoustics MMD matrix
% Zhilei 08/28/2025
% also calculate off-diagnonal accuracy, then compare to the MMD matrix
% differ from v1: troubleshoot why bin-sum is not consistent with raw raster plot

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/calcDistance/MMD'));


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
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet



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
fd_mmd = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD', 'matrixCall2', 'v4_v5');
sigma = 0.9;
fn_mmd = fullfile(fd_mmd, sprintf('%s.%s%s.%.2f.mmd.mat', birdID, syls{1}, syls{2}, sigma));
load(fn_mmd);
mmd_matrix = mmd.mmd;
% what's the window size
win_size = 32; 
sec_per_frame = 0.001;  % how much sec per frame of data
% plot the diaganol
diag_mmd = diag(mmd_matrix);
diag_shuff = diag(mmd.mmd11);
diag_shuff = diag_shuff(1:size(diag_mmd,1));
% relative time of the sliding window start
rel_t_mmd = ((1:length(diag_mmd)) -1 - win_size) * sec_per_frame;
close all; 
figure; plot(rel_t_mmd, diag_mmd, 'LineWidth', 2); hold on;
plot(rel_t_mmd, diag_shuff, 'LineWidth', 2, 'Color', '#737373'); hold on;
xlim([rel_t_mmd(1) rel_t_mmd(end)]);



%% 2. Time-warp renditions to each call subtype
% no cross warping
fs = 20000;
% pad = 0.032; 
pad = 0.08;
% what's the size of bins when counting spikes, unit is sec
% bin = 0.005; 
bin = 0.01;
bin_pt = floor(fs*bin);
% hop = 0.001; 
hop = 0.002;
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

% count spikes in non-overlapping bins
% bin_sums1 = ZZfunc_countSpikeBins_v1(aligned_spike1, bin_pt);
% bin_sums2 = ZZfunc_countSpikeBins_v1(aligned_spike2, bin_pt);
% % calculate relative time of the bin start
% max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
% rel_t = bin*(0:(max_bin_num-1)) - pad;

% count spikes in overlapping bins
bin_sums1 = ZZfunc_countSpikeSlidingWin_v1(aligned_spike1, bin_pt, hop_pt);
bin_sums2 = ZZfunc_countSpikeSlidingWin_v1(aligned_spike2, bin_pt, hop_pt);
max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
rel_t = hop*(1:max_bin_num) - hop - pad;

% further smooth?
% smooth_win = 10;
% bin_sums1 = smoothdata(bin_sums1, 2, 'gaussian', smooth_win);
% bin_sums2 = smoothdata(bin_sums2, 2, 'gaussian', smooth_win);

% check a specific neuron
nID = '20240917-ch13';
idx1 = find(strcmp({seg1.neuronID}, nID));
idx2 = find(strcmp({seg2.neuronID}, nID));
% close all;
[fig, axes] = generatePanelGrid_v2(1, 2, [0.8], [], [0.05;0.05], [0.1;0.05], 0.1, [0], [10 10 1600 600]);
imagesc(axes(1), bin_sums1(idx1,:), [0 5]); 
title(axes(1), sprintf('%s %s', syls{1}, nID), 'FontSize', 14); 
colormap(axes(1), custom_map); colorbar(axes(1));

imagesc(axes(2), bin_sums2(idx2,:), [0 5]); 
colormap(axes(2), custom_map); colorbar(axes(2));
title(axes(2), sprintf('%s %s', syls{2}, nID), 'FontSize', 14); 



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


% or load preivously saved pseudo-population trials
fn_pop = fullfile(fd_save, sprintf('%s.%s.pop.mat', birdID, syls{1})); load(fn_pop);
fn_pop = fullfile(fd_save, sprintf('%s.%s.pop.mat', birdID, syls{2})); load(fn_pop);


%% 4. Build SVM models to decode call identity: calculate off-diagonal as well
bin1 = size(pop1, 3); 
bin2 = size(pop2, 3);
% how many train/test shuffling to do
num_runs = 10;
train_ratio = 0.5;
kernel = 'linear';
% loop through time points
accu_all = nan(bin1, bin2, 2, num_runs);
parfor t1=1:bin1
  for t2=1:bin2
    X1 = squeeze(pop1(:, :, t1));
    X2 = squeeze(pop2(:, :, t2));
    Y1 = repmat({syls{1}}, size(X1,1), 1);
    Y2 = repmat({syls{2}}, size(X2,1), 1);
    r_seed = 1992;
    [accu_real, accu_shuffle] = ZZfunc_trainTestSVM_v1(X1, X2, Y1, Y2, train_ratio, kernel, num_runs, r_seed);

    accu_all(t1, t2, :,:) = [accu_real'; accu_shuffle'];
  end
end
% save results
fn_accu = fullfile(fd_save, sprintf('%s.%s_%s.accu_allMatrix.mat', birdID, syls{1}, syls{2}));
save(fn_accu, 'accu_all');

% load
load(fn_accu);
% calculate then plot mean accuracy
accu_mean = mean(accu_all, 4);
% perform rank-1 decomposion
m1 = squeeze(accu_mean(:,:,1));
[B, C, u, v] = mmd_matrix_decompose(m1);
% get time
t1 = rel_t(1:size(m1,1));
t2 = rel_t(1:size(m1,2));

close all;
[fig, axes] = generatePanelGrid_v2(2, 2, [0.35;0.35], [0.07], [0.05;0.05], [0.1;0.05], 0.08, [0;0], [10 10 950 900]);
ds = {m1; squeeze(accu_mean(:,:,2)); B; C};
titles = {'Real', 'Shuffled', 'Rank-1', 'Residual'};
clims = {[0.5 0.9]; [0.5 0.9]; [0.5 0.9]; [-0.5 0]};
% xlims = [-0.05 0.16];
xlims = [-0.032 0.14];
for di=1:size(ds,1)
  plot_i = floor((di-1)/2)+1;
  plot_j = mod(di-1,2)+1;
  ax = axes(plot_i, plot_j);
  imagesc(ax, t1, t2,  ds{di}, clims{di});
  colormap(ax, 'gray');
  colorbar(ax);
  xlim(ax, xlims); 
  ylim(ax, xlims);
  title(ax, titles{di}, 'FontSize', 14);
end
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.matrixSVM.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');

% plot the raw matrix only
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 700 600]);
ax = gca; hold(ax, 'on');
xlims = [-0.032 0.14];
imagesc(ax, t1, t2,  m1, [0.5 0.9]);
colormap(ax, 'gray');
colorbar(ax);
xlim(ax, xlims); 
ylim(ax, xlims);
set(ax, 'YDir', 'reverse');
title(ax, sprintf('Neural SVM: %s vs %s', syls{1}, syls{2}), 'FontSize', 14);
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.matrixSVMraw.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');



%% Troubleshoot the high off-diagonal values
[~, idx1] = min(abs(t1-0.06));
[~, idx2] = min(abs(t2-0.02));
% idx1=29; idx2=21;
% idx1=31; idx2=21;
fprintf('tx=%.3f, ty=%.3f, value=%.3f\n', t1(idx1), t2(idx2), m1(idx1, idx2));
close all; fig = ZZfunc_newFigurePDFsize_v1([10 10 650 600]);
ax = gca; hold(ax, 'on'); set(ax, 'YDir', 'reverse');
imagesc(ax, t1, t2,  m1, [0.5 0.9]);
colormap(ax, 'gray');
colorbar(ax);
plot(ax, t2(idx2), t1(idx1), 'ro'); 
xlim(ax, xlims); 
ylim(ax, xlims);

% check the raw data and trained SVM
X1 = squeeze(pop1(:, :, idx1));
X2 = squeeze(pop2(:, :, idx2));
Y1 = repmat({syls{1}}, size(X1,1), 1);
Y2 = repmat({syls{2}}, size(X2,1), 1);
r_seed = 1992;
train_ratio = 0.5; 
kernel = 'linear'; 
num_runs = 20;
[accu_real, accu_shuffle] = ZZfunc_trainTestSVM_v1(X1, X2, Y1, Y2, train_ratio, kernel, num_runs, r_seed);
disp(mean(accu_real));

% plot sampled vector
idx_rd = randsample(1:size(X1,1), 50);
% idx_rd = 1:size(X1,1);
% close all;
[fig, axes] = generatePanelGrid_v2(1, 2, [0.7], [], [0.05;0.05], [0.1;0.05], 0.1, [0], [10 10 1600 600]);
imagesc(axes(1), X1(idx_rd,:)', [0 5]); 
title(axes(1), sprintf('%s t=%.3f', syls{1}, t1(idx1)), 'FontSize', 14); 
colormap(axes(1), custom_map); colorbar(axes(1));
yticks(axes(1), 1:size(X1,2)); yticklabels(axes(1), neuron_sort);

imagesc(axes(2), X2(idx_rd,:)', [0 5]); title(axes(2), syls{2}, 'FontSize', 14); 
colormap(axes(2), custom_map); colorbar(axes(2));
yticks(axes(2), 1:size(X2,2)); yticklabels(axes(2), neuron_sort);
title(axes(2), sprintf('%s t=%.3f', syls{2}, t2(idx2)), 'FontSize', 14); 



%%  Plot MMD and SVM side-by-side
% plot acoustic MMD and neural SVM side-by-side
close all; 
[fig, axes] = generatePanelGrid_v2(1, 2, [0.8], [], [0.1;0.05], [0.1;0.05], 0.08, [0], [10 10 1200 500]);
C = mmd.C';
rel_t_mmdx = ((1:size(C,1)) -1 - win_size) * sec_per_frame;
rel_t_mmdy = ((1:size(C,2)) -1 - win_size) * sec_per_frame;
imagesc(axes(1), rel_t_mmdx, rel_t_mmdy, C, [-0.2 -0.02]);
colormap(axes(1), 'gray'); colorbar(axes(1));
xlabel(axes(1),sprintf('%s rel. time (sec)', syls{2}), 'FontSize', 14);
ylabel(axes(1),sprintf('%s rel. time (sec)', syls{1}), 'FontSize', 14);
xlim(axes(1), [-0.032 0.14]); ylim(axes(1), [-0.032 0.14]);
title(axes(1), 'Acoustic distance (MMD)', 'FontSize', 14);

imagesc(axes(2), t1, t2,  m1, [0.5 0.9]);
colormap(axes(2), 'gray'); colorbar(axes(2));
xlabel(axes(2),sprintf('%s rel. time (sec)', syls{2}), 'FontSize', 14);
ylabel(axes(2),sprintf('%s rel. time (sec)', syls{1}), 'FontSize', 14);
xlim(axes(2), [-0.032 0.14]); ylim(axes(2), [-0.032 0.14]);
title(axes(2), 'Neural decoding accuracy (SVM)', 'FontSize', 14);

fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.MMDvsSVMmatrix.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');


%% plot the diagonal in the same plot
min_bin = min([size(accu_all,1) size(accu_all,2)]);
diag_svm = zeros(min_bin, size(accu_all,3), size(accu_all,4)); % Size 50 x 2 x 100
% Loop over 3rd and 4th dimensions
for i = 1:size(accu_all,3)
    for j = 1:size(accu_all,4)
        diag_svm(:, i, j) = diag(accu_all(:, :, i, j));
    end
end

close all;
col_list = {'#4db34d', '#f15a29'};
fig = ZZfunc_newFigurePDFsize_v1([10 10 650 600]);
ax = gca; hold(ax, 'on');
xlim1 = max([rel_t_mmd(1) rel_t(1)]);
xlim2 = min([rel_t_mmd(end) rel_t(end)]);
yyaxis left
plot_mean_with_ci(ax, rel_t(1:min_bin), squeeze(diag_svm(:,1,:))', col_list{1}, 2, '-', hexToRGB(col_list{1}), 0.5);
plot_mean_with_ci(ax, rel_t(1:min_bin), squeeze(diag_svm(:,2,:))', col_list{1}, 2, '--', hexToRGB(col_list{1}), 0.1);
ylabel(ax, 'SVM decoding accuracy', 'FontSize', 14);
ylim(ax, [0.45, 1]);
ax = ZZfunc_addSimpleLegend_v2(ax, {'Neural: real', 'Acoustic'}, col_list, 16, [-0.02 -0.02 -0.02], [0.95 0.92 0.89]);

yyaxis right
plot(ax, rel_t_mmd, diag_mmd, 'LineWidth', 2, 'Color', col_list{2});
plot(ax, rel_t_mmd, diag_shuff, 'LineWidth', 2, 'Color', col_list{2});
xlim(ax, [xlim1 xlim2]);
ylim(ax, [0 max(diag_mmd)*1.05]);
ylabel(ax, 'MMD on VAE latents', 'FontSize', 14);
xlabel(ax, 'Rel. time (sec)', 'FontSize', 14);

ax.YAxis(1).Color = col_list{1};
ax.YAxis(2).Color = col_list{2};
title(ax, 'Neural vs acoustic distance', 'FontSize', 16);
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.MMDvsSVM_diag.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');



%% Perform cross-correlation on diagonal
% first on real data
x = mean(squeeze(diag_svm(:,1,:)), 2);
x = x - mean(x); 
t_x = rel_t(1:size(x, 1)); 
y = diag_mmd;
y = y - mean(y);
t_y = rel_t_mmd;
[corr, time_lags] = cross_corr_plot(x, t_x, y, t_y);

% then on shuffled data
x = mean(squeeze(diag_svm(:,2,:)), 2);
x = x - mean(x); 
t_x = rel_t(1:size(x, 1)); 
y = diag_shuff;
y = y - mean(y);
t_y = rel_t_mmd;
[corr2, time_lags2] = cross_corr_plot(x, t_x, y, t_y);

% plot in the same figure
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 400 600]);
ax = gca; hold(ax, 'on');
plot(time_lags, corr, 'Color', '#737373', 'LineStyle', '-', 'LineWidth', 3); 
plot(time_lags2, corr2, 'Color', '#737373', 'LineStyle', '--', 'LineWidth', 2);
xline(0, 'Color', '#737373', 'LineStyle', '--', 'LineWidth', 1);
% annotate the max point
[maxv, maxi] = max(corr);
plot(time_lags(maxi), corr(maxi), 'Marker', '+', 'MarkerSize', 16, 'Color', '#e41a1c', 'LineWidth', 3);
text(time_lags(maxi), corr(maxi)+0.05, sprintf('%d ms', floor(time_lags(maxi)*1000)), 'FontSize', 14, 'Color', '#e41a1c');
xlabel('Time lag (sec)', 'FontSize', 14);
ylabel('Cross correlation', 'FontSize', 14);
% ylabel('Normalized cross correlation', 'FontSize', 14);
xlim([-0.1 0.1]);
ylim([min(corr)-0.1 max(corr)+0.1]);
title('Cross-correlation', 'FontSize', 14);
fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.crossCorrelation.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');







%% Perform cross-correlation on matrix?
% trim the neural SVM matrix to desired range
[A_trimmed, tiA_trimmed, tjA_trimmed] = trim_matrix(m1, t2, t1, [-0.032 0.14], [-0.032 0.14]);

% trim the acoustic MMD matrix to desired range
[C_trimmed, tiC_trimmed, tjC_trimmed] = trim_matrix(C, rel_t_mmdy, rel_t_mmdx, [-0.032 0.14], [-0.032 0.14]);

% interpolate the neural matrix to the same range and resolutions as MMD matrix
[Xi, Yi] = meshgrid(tiC_trimmed, tjC_trimmed);
A_interpolated = interp2(tiA_trimmed, tjA_trimmed, A_trimmed, Xi, Yi, 'linear', 0);

% perform matrix cross-correlation
% first normalize the matrix
A = normalize_matrix(A_interpolated, 0.5, 0.9);
B = normalize_matrix(C_trimmed, -0.2, -0.02);








