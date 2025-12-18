% Method 1 of decoding call identity from neural responses
% build a decoding model for each neuron, calculate the decoding accuracy.
% Then plot the accuracy as y-axis, while the burst location as x-axis. Compare the curve to the acoustics MMD curve
% Zhilei 08/25/2025


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
% syls = {'v4', 'v5'};
syls = {'v1', 'v7'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
fd_embed = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID);
suffix_hahn = 'neuron_orderedPlotted5';
suffix_criteria = 'criteria5';
% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decode']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end



%% 1. Read ephys data
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



%% 2. Loop through neurons, build decoding model
% save plots of individual neurons
fd_save_plot = fullfile(fd_save, sprintf('accuracy_plots_%s_%s', syls{1}, syls{2}));
if ~exist(fd_save_plot, 'dir'); mkdir(fd_save_plot); end
% ignore neurons that don't have enough renditions
min_rends = 20;
fs = 20000;
pad = 0.08; % how much pre and post data to include, unit is sec
% what's the size of bins when counting spikes, unit is sec
bin = 0.005;
bin_pt = floor(fs*bin);
num_runs = 100;  % how many runs to perform

neuron_comb = unique([n_all{:}], 'sorted');
accu_all = nan(size(neuron_comb, 2), 2, num_runs);
num_rends = nan(size(neuron_comb,2), 4);
for ni=1:size(neuron_comb, 2)
  % neuronID = '20240912-ch11';
  neuronID = neuron_comb{ni};
  fprintf('Building SVM for %s\n', neuronID);
  seg1 = e_all{1}(strcmp({e_all{1}.neuronID}, neuronID));
  seg2 = e_all{2}(strcmp({e_all{2}.neuronID}, neuronID));
  num_rends(ni,1:2) = [size(seg1,2) size(seg2,2)];
  if (size(seg1,2)>= min_rends) && (size(seg2,2)>= min_rends)
    % warp spikes
    [aligned_spike1, aligned_sound1, mean_dur1] = ZZfunc_linearWarp_v2(seg1, pad);
    [aligned_spike2, aligned_sound2, mean_dur2] = ZZfunc_linearWarp_v2(seg2, pad);
    
    % how many trials have spikes
    total_sum1 = sum(aligned_spike1, 2);
    total_sum2 = sum(aligned_spike2, 2);
    num_rends(ni, 3) = sum(total_sum1>0);
    num_rends(ni, 4) = sum(total_sum2>0);
    
    % count the number of spikes in bins
    bin_sums1 = ZZfunc_countSpikeBins_v1(aligned_spike1, bin_pt);
    bin_sums2 = ZZfunc_countSpikeBins_v1(aligned_spike2, bin_pt);
    % truncate to have the same number of bins
    num_bins = min([size(bin_sums1,2) size(bin_sums2,2)]);
    X1 = bin_sums1(:, 1:num_bins);
    X2 = bin_sums2(:, 1:num_bins);
    
    % train a SVM to decode identity
    Y1 = repmat({syls{1}}, size(X1,1), 1);
    Y2 = repmat({syls{2}}, size(X2,1), 1);
    train_ratio = 0.7;
    kernel = 'linear';
    r_seed = 1992;
%     select_method = 'max';
    select_method = 'random';
%     [accu_real, accu_shuffle] = ZZfunc_trainTestSVM_v1(X1, X2, Y1, Y2, train_ratio, kernel, num_runs, r_seed);
    [accu_real, accu_shuffle] = ZZfunc_trainTestSVM_v2(X1, X2, Y1, Y2, train_ratio, kernel, num_runs, r_seed, select_method);
    
    accu_all(ni, 1, :) = accu_real;
    accu_all(ni, 2, :) = accu_shuffle;
    
    % plot results as violin plots
    close all;
    fig = ZZfunc_newFigurePDFsize_v1([10 10 300 600]);
    d = [accu_real; accu_shuffle];
    categ = [repmat({'Real'}, length(accu_real), 1); repmat({'Shuffle'}, length(accu_shuffle), 1)];
    violinplot(d, categ, 'ViolinColor', [0.3 0.68 0.29; 0.2 0.2 0.2], 'EdgeColor', [1 1 1], 'MarkerSize', 12, 'MedianMarkerSize', 200);
    yline(0.5, 'LineStyle', '--', 'LineWidth', 1);
    xlim([0.5 2.5]);
    ylim([0 1]);
    ax = gca; ax.XAxis.FontSize = 16;
    title(sprintf('%s: %s-%s', neuronID, syls{1}, syls{2}), 'FontSize', 12);
    fn_pdf = fullfile(fd_save_plot, sprintf('%s.%s-%s.%s.SVM.pdf', birdID, syls{1}, syls{2}, neuronID));
    print(fig, fn_pdf, '-dpdf', '-painters');
  end
end



%% 3. Plot accuracy against burst time
accu_mean = nanmean(accu_all, 3);
accu_time = [];
c1 = c_all{1}; c2 = c_all{2};
for ni=1:size(neuron_comb, 2)
  neuronID = neuron_comb{ni};
  accu_time(ni).neuronID = neuronID;
  accu_time(ni).accu_real = accu_mean(ni, 1);
  accu_time(ni).accu_shuffle = accu_mean(ni, 2);
  accu_time(ni).rends1 = num_rends(ni, 1);
  accu_time(ni).rends2 = num_rends(ni, 2);
  accu_time(ni).fire_rends1 = num_rends(ni, 3);
  accu_time(ni).fire_rends2 = num_rends(ni, 4);
  accu_time(ni).prop_fire1 = num_rends(ni,3)/num_rends(ni,1);
  accu_time(ni).prop_fire2 = num_rends(ni,4)/num_rends(ni,2);
  % find the time of burst for the neuron
  if ismember(neuronID, n_all{1})
    accu_time(ni).t1 = c1(strcmp({c1.neuronID}, neuronID)).psth_max_smooth_t;
  else
    accu_time(ni).t1 = nan;
  end
  if ismember(neuronID, n_all{2})
    accu_time(ni).t2 = c2(strcmp({c2.neuronID}, neuronID)).psth_max_smooth_t;
  else
    accu_time(ni).t2 = nan;
  end
end
% calculate mean time
accu_time = struct2table(accu_time);
accu_time.mean_t = nanmean(accu_time{:, {'t1', 't2'}}, 2);

% sort by burst time
sorted_table = sortrows(accu_time, 'mean_t');
% save results
fn_table = fullfile(fd_save_plot, 'sorted_table.mat');
save(fn_table, 'sorted_table');

% plot
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 700 500]);
x = sorted_table.mean_t;
y = sorted_table.accu_real;
y2 = sorted_table.accu_shuffle;
labels = sorted_table.neuronID; 
prop_thre = 0; 
% idx_nan = find((~isnan(x)) & (~isnan(y)));
idx_nan = find((~isnan(x)) & (~isnan(y)) & ((sorted_table.prop_fire1>=prop_thre) | (sorted_table.prop_fire2>=prop_thre)));
plot(x(idx_nan), y(idx_nan), 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', '#4daf4a', ...
  'MarkerEdgeColor', 'none', 'LineStyle', '-', 'LineWidth', 1, 'Color', '#4daf4a'); hold on;
plot(x(idx_nan), y2(idx_nan), 'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', '#737373', ...
  'MarkerEdgeColor', 'none', 'LineStyle', '-', 'LineWidth', 1, 'Color', '#737373'); hold on;
% add a small text for easy troubleshoot
ylim([0.4 1]);
xlabel('Neuron burst time (sec)', 'FontSize', 14);
ylabel('SVM eecoding accuracy', 'FontSize', 14);
ax = gca;
ax = ZZfunc_addSimpleLegend_v2(ax, {'Real', 'Shuffle'}, {'#4daf4a', '#737373'}, 16, [-0.02 -0.02], [0.95 0.9]);
title(ax, sprintf('%s: %s vs %s', birdID, syls{1}, syls{2}), 'FontSize', 14)

text(x(idx_nan), y(idx_nan), labels(idx_nan), 'FontSize', 8);

fn_pdf = fullfile(fd_save_plot, sprintf('%s.%s-%s.mean.SVM.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');








