% Method 2 of decoding call identity from neural responses
% pseudo-population analysis. For time-warped data, could treat that all sparse MO neurons (~30) are recorded simultaneously 
% by randomly sampling renditions, then neural responses become a 30-dim vector. We can then ask if the activity at timepoint 
% x of call A is similar to activity at timepoint y of call B. Then we can build a similarity matrix of the neural activity, 
% which can be compared to the acoustics MMD matrix
% Zhilei 09/01/2025
% also calculate off-diagnonal accuracy, then compare to the MMD matrix; Run through all pairs
% differ from v3: pre-align all renditions of a subtype to avoid difference in warped duration

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/calcDistance/MMD'));


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
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix4'], 'intermediate');
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
num_runs = 5;
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
imagesc(ax, t1, t2,  m1, [0.5 0.9]);
colormap(ax, 'gray');
colorbar(ax);
xlim(ax, xlims); 
ylim(ax, xlims);
set(ax, 'YDir', 'reverse');
title(ax, sprintf('Neural SVM: %s vs %s', syls{1}, syls{2}), 'FontSize', 14);
fn_pdf = fullfile(fd_save, sprintf('%s.%s.matrixSVMraw.pdf', birdID, run_name));
print(fig, fn_pdf, '-dpdf', '-painters');

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






















