% Method 2 of decoding call identity from neural responses
% pseudo-population analysis. For time-warped data, could treat that all sparse MO neurons (~30) are recorded simultaneously 
% by randomly sampling renditions, then neural responses become a 30-dim vector. We can then ask if the activity at timepoint 
% x of call A is similar to activity at timepoint y of call B. Then we can build a similarity matrix of the neural activity, 
% which can be compared to the acoustics MMD matrix
% Zhilei 09/01/2025
% also calculate off-diagnonal accuracy, then compare to the MMD matrix; Run through all pairs
% differ from v3: pre-align all renditions of a subtype to avoid difference in warped duration
% throubleshoot: considering neural coverage, e.g. filtering out missing data portion

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/calcDistance/MMD'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.replaced.mat', birdID)));
% what syllable types to analyze
syls_all = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
fd_embed = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID);
suffix_hahn = 'neuron_orderedPlotted6';
suffix_criteria = 'criteria6';
% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix4'], 'intermediate_ts');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end

% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet

% read the ephys struct into RAM so no need to re-load each time
struct_all = struct();
aligned = struct();
pad = 0.08;
for sii=1:size(syls_all,2)
  ss = syls_all{sii};
  fprintf('Load data for %s\n', ss);
  fn_e = fullfile(fd_ephys, sprintf('%s.%s.segments_all.replaced.mat', birdID, ss));
  a = load(fn_e);
  struct_all.(ss) = a.seg_selected;
  % align all renditions to use later
  [aligned_spike, ~, ~] = ZZfunc_linearWarp_v2(a.seg_selected, pad);
  aligned.(ss) = aligned_spike;
end


for sii=1:size(syls_all,2)
  for sjj=sii:size(syls_all,2)
    syls = {syls_all{sii}, syls_all{sjj}};
    disp(syls);
    run_name = sprintf('%s_%s', syls{1}, syls{2});
    
%% 1. Read ephys data 
e_all = cell(size(syls, 2), 1);  % raw ephys struct
a_all = cell(size(syls, 2), 1);
n_all = cell(size(syls, 2), 1);  % selected and plotted neurons
c_all = cell(size(syls, 2), 1);  % criteria struct that includes burst info
% syl_i = 1;
for syl_i=1:size(syls, 2)
  ss = syls{syl_i};
%   fn_e = fullfile(fd_ephys, sprintf('%s.%s.segments_all.replaced.mat', birdID, ss));
%   load(fn_e); e_all{syl_i} = seg_selected(strcmp({seg_selected.aud_ch}, 'chan0'));
  e_all{syl_i} = struct_all.(ss);
  a_all{syl_i} = aligned.(ss);
  fn_n = fullfile(fd_hahnloser, 'popRaster', ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_hahn));
  load(fn_n); n_all{syl_i} = neuron_ordered;
  fn_c = fullfile(fd_hahnloser, 'popRaster', ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_criteria));
  load(fn_c); c_all{syl_i} = criteria;
end


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
% set a threshold for minimal number of renditions to be include in the analysis
min_rend = 15;
c1 = struct2table(c_all{1}); 
c2 = struct2table(c_all{2}); 
cjoin = outerjoin(c1, c2, 'Keys','neuronID', 'MergeKeys',true);
neuron_comb = unique([n_all{:}], 'sorted');
cjoin = cjoin(ismember(cjoin.neuronID, neuron_comb),:);
cjoin = cjoin((cjoin.num_rends_c1>=min_rend) & (cjoin.num_rends_c2>=min_rend), :);
% compute a mean time
cjoin.mean_t = nanmean([cjoin.psth_max_smooth_t_c1 cjoin.psth_max_smooth_t_c2], 2);
% sort by mean time
n_sort = sortrows(cjoin, 'mean_t');
neuron_sort = n_sort.neuronID;

% retrieve spike data
idx1 = find(ismember({e_all{1}.neuronID}, neuron_sort));
idx2 = find(ismember({e_all{2}.neuronID}, neuron_sort));
seg1 = e_all{1}(idx1);
seg2 = e_all{2}(idx2);
aligned_spike1 = a_all{1}(idx1,:);
aligned_spike2 = a_all{2}(idx2,:);
% [aligned_spike1, aligned_sound1, mean_dur1] = ZZfunc_linearWarp_v2(seg1, pad);
% [aligned_spike2, aligned_sound2, mean_dur2] = ZZfunc_linearWarp_v2(seg2, pad);

% count spikes in overlapping bins
bin_sums1 = ZZfunc_countSpikeSlidingWin_v1(aligned_spike1, bin_pt, hop_pt);
bin_sums2 = ZZfunc_countSpikeSlidingWin_v1(aligned_spike2, bin_pt, hop_pt);
max_bin_num = max([size(bin_sums1,2) size(bin_sums2,2)]);
rel_t = hop*(1:max_bin_num) - hop - pad;


% plot the coverage over time
close all; figure; set(gcf, 'Position', [10 10 1200 800]);
ax1=subplot(2, 1, 1); hold(ax1, 'on');
plot(ax1, rel_t(1:size(bin_sums1,2)), mean(bin_sums1, 1), 'Color', '#737373', 'LineStyle', '-', 'DisplayName', syls{1});
% then add the location of the sparse neurons
col_lines = lines(size(n_sort,1));
for nii=1:size(n_sort, 1); xline(ax1, n_sort.psth_max_smooth_t_c1(nii), 'LineStyle', '--', 'Color', col_lines(nii,:), 'LineWidth', 2); end
ax2=subplot(2, 1, 2); hold(ax2, 'on');
plot(ax2, rel_t(1:size(bin_sums2,2)), mean(bin_sums2, 1), 'Color', '#737373', 'LineStyle', '-', 'DisplayName', syls{2});
% then add the location of the sparse neurons
for nii=1:size(n_sort, 1); xline(ax2, n_sort.psth_max_smooth_t_c2(nii), 'LineStyle', '--', 'Color', col_lines(nii,:), 'LineWidth', 2); end
linkaxes([ax1 ax2], 'x');

% plot the population raster plot according to a common neuron order
neuron_ordered = neuron_sort(end:-1:1);
neuron_ordered = neuron_ordered';
A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#737373'};
N = length(neuron_ordered);
neuron_color = A(mod(0:N-1, numel(A))+1);
to_sample = 20;  % max number of renditions to sample, plot all renditions if -1
pad_grey = true;  % if a neuron doesn't have to_sample number of renditions, pad empty rows
tick_width = 8; % size of spike ticks, default to 8 for old plots
fig_size = [10 10 400 900];  % in unit of pixels, [10 10 500 900] for old plots
r_seed = 1992; % randome seed to ensure reproductivity 
sample_method = 'max';   % how to choose what renditions to plot: random, max; 
sampled_loc = false;  % whether to determine neuron location based on sampled renditions

seg_selected = struct_all.(syls{1});
fn_plot = fullfile(fd_save, sprintf('troubleshoot.popRaster.%s_%s.%s.fig', syls{1}, syls{2}, syls{1}));
[aligned_spike, aligned_sound, neuron_ordered2, fig, sampled_rends_ordered2] =  ZZfunc_alignSpikeCall_v10_mean(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, pad_grey, tick_width, fig_size, r_seed, sample_method, sampled_loc);
% also export as pdf
exportgraphics(fig, strrep(fn_plot, 'fig', 'tiff'), 'Resolution', 600); 

seg_selected = struct_all.(syls{2});
fn_plot = fullfile(fd_save, sprintf('troubleshoot.popRaster.%s_%s.%s.fig', syls{1}, syls{2}, syls{2}));
[aligned_spike, aligned_sound, neuron_ordered2, fig, sampled_rends_ordered2] =  ZZfunc_alignSpikeCall_v10_mean(seg_selected, fn_plot, neuron_ordered, neuron_color, to_sample, pad_grey, tick_width, fig_size, r_seed, sample_method, sampled_loc);
% also export as pdf
exportgraphics(fig, strrep(fn_plot, 'fig', 'tiff'), 'Resolution', 600); 


%% 3. Randomly sample renditions to form pseudo-population data
num_runs = 1000;
pop1 = nan(num_runs, size(neuron_sort,1), size(bin_sums1,2));
pop2 = nan(num_runs, size(neuron_sort,1), size(bin_sums2,2));
info = nan(num_runs, size(neuron_sort,1), 2);  % also record what renditions are sampled
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
  fn_pdf = fullfile(fd_save, sprintf('%s.%s.%s.pseudo.pdf', birdID, run_name, syls{di}));
  print(fig, fn_pdf, '-dpdf', '-painters');
end

% also save the sample renditions for later use
fn_pop = fullfile(fd_save, sprintf('%s.%s.%s.pop.mat', birdID, run_name, syls{1}));
save(fn_pop, 'pop1');
fn_pop = fullfile(fd_save, sprintf('%s.%s.%s.pop.mat', birdID, run_name, syls{2}));
save(fn_pop, 'pop2');



%% 4. Build SVM models to decode call identity: calculate off-diagonal as well
bin1 = size(pop1, 3); 
bin2 = size(pop2, 3);
% how many train/test shuffling to do
num_runs = 3;
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
fn_accu = fullfile(fd_save, sprintf('%s.%s.accu_allMatrix.mat', birdID, run_name));
save(fn_accu, 'accu_all');


%% Mask regions where there is no coverage
width = 0.01;  %given a sparse neuron, how much width from the peak is considered covered, unit is second





% load(fn_accu);
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
fn_pdf = fullfile(fd_save, sprintf('%s.%s.matrixSVM.pdf', birdID, run_name));
print(fig, fn_pdf, '-dpdf', '-painters');

% plot the raw matrix only
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 700 600]);
ax = gca; hold(ax, 'on');
xlims = [-0.032 0.14];
imagesc(ax, t2, t1,  m1, [0.5 0.9]);
colormap(ax, 'gray');
colorbar(ax);
% xlim(ax, xlims); 
% ylim(ax, xlims);
axis(ax, 'equal');
axis(ax, 'tight');
set(ax, 'YDir', 'reverse');
title(ax, sprintf('Neural SVM: %s vs %s', syls{1}, syls{2}), 'FontSize', 14);
fn_pdf = fullfile(fd_save, sprintf('%s.%s.matrixSVMraw.pdf', birdID, run_name));
print(fig, fn_pdf, '-dpdf', '-painters');

% plot defined time range
tpad_n = [0 -0.032];
ty=t1; tx=t2;
tendy = ty(end) - pad + tpad_n(2);
iy = find((ty>=-tpad_n(1)) & (ty<=tendy));
tendx = tx(end) - pad + tpad_n(2);
ix = find((tx>=-tpad_n(1)) & (tx<=tendx));
n_this = m1(iy, ix);
fig = ZZfunc_newFigurePDFsize_v1([10 10 700 600]);
ax = gca; hold(ax, 'on');
imagesc(ax, tx(ix), ty(iy),  n_this, [0.5 0.8]);
colormap(ax, 'gray');
colorbar(ax);
set(ax, 'YDir', 'reverse');
axis(ax, 'equal');
axis(ax, 'tight');


% check a specific time pair
tt1 = 0.05; 
tt2 = 0.014; 
[~, i1] = min(abs(t1-tt1));
[~, i2] = min(abs(t2-tt2));
% get data, build SVM
X1 = squeeze(pop1(:, :, i1));
X2 = squeeze(pop2(:, :, i2));
Y1 = repmat({syls{1}}, size(X1,1), 1);
Y2 = repmat({syls{2}}, size(X2,1), 1);
r_seed = 1992;
[accu_real, accu_shuffle] = ZZfunc_trainTestSVM_v1(X1, X2, Y1, Y2, train_ratio, kernel, num_runs, r_seed);
fprintf('Mean accu: %.3f\n', mean(accu_real));
% plot the raw data to check
rng(1118);
to_plot = 100;
xx1 = X1(randsample(1:size(X1,1), to_plot), :);
xx2 = X2(randsample(1:size(X2,1), to_plot), :);
figure; set(gcf, 'Position', [10 10 1200 600]);
subplot(1,2,1); imagesc(xx1', [0 5]); colormap(custom_map); yticks(1:size(neuron_sort,1)); yticklabels(neuron_sort);
subplot(1,2,2); imagesc(xx2', [0 5]); colormap(custom_map); yticks(1:size(neuron_sort,1)); yticklabels(neuron_sort);


% plot the diagonal
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
plot_mean_with_ci(ax, rel_t(1:min_bin), squeeze(diag_svm(:,1,:))', col_list{1}, 2, '-', hexToRGB(col_list{1}), 0.5);
plot_mean_with_ci(ax, rel_t(1:min_bin), squeeze(diag_svm(:,2,:))', col_list{1}, 2, '--', hexToRGB(col_list{1}), 0.1);
ylabel(ax, 'SVM decoding accuracy', 'FontSize', 14);
ylim(ax, [0.4, 1]);
xlim(ax, [-0.032 0.16]);
fn_pdf = fullfile(fd_save, sprintf('%s.%s.SVMdiagonal.pdf', birdID, run_name));
print(fig, fn_pdf, '-dpdf', '-painters');

  end
end






















