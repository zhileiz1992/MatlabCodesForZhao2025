% quantity acoustic distance on VAE latents with different time lags to the spike onset
% Zhilei, 10/16/2025
% differ from v2_winCenter: compare between calls and compounds 

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/overlaySpikes'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/compoundTraj'));


%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
% what's the range to include spikes: default only consider spikes within 1-window of syllable boundaries
epad_frame = win_frame;

for bi=2:size(birdIDs,2)
% bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
% where is VAE data located
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
  
% information regarding spikes and sliding windows (zero lags)
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField']);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
load(fn_spike);

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end



%% Loop through neurons, calculate acoustic distance with different lags to spike
lags = -100:5:200;
dist_mean_all = nan(size(info_neu,1), length(lags));
dist_all_all = cell(size(info_neu,1),1);

% ni=44;
for ni=1:size(info_neu,1)
  
neuronID = info_neu.neuronID{ni};
spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
% compare between compound syllables and calls
spike_comp = spike_this((ismember(spike_this.category, {'b','x'})) & (spike_this.dur>=0.3),:);
spike_call = spike_this(ismember(spike_this.category, {'v'}),:);
% only look at syllables that have spikes
idx_comp = cellfun(@(x) ~isempty(x), spike_comp.win_loc);
spike_has_comp = spike_comp(idx_comp,:);
idx_call = cellfun(@(x) ~isempty(x), spike_call.win_loc);
spike_has_call = spike_call(idx_call,:);
if isempty(spike_has_comp) || isempty(spike_has_call)
  continue;
end

%% Find the location of the burst
% try first spike of the burst or center of mass
% burst_method = 'first';
burst_method = 'center';
spike_has_comp.burst_loc = cellfun(@(x) ZZfunc_findBurstLoc_v1(x, 10, burst_method), spike_has_comp.win_loc, 'UniformOutput', false);
spike_has_call.burst_loc = cellfun(@(x) ZZfunc_findBurstLoc_v1(x, 10, burst_method), spike_has_call.win_loc, 'UniformOutput', false);

%% 1. Grab the VAE latents with different lags to spikes
% loop through syllables, then loop through spikes, grab the VAE data with different lags
d_vae1 = cell(size(spike_has_comp,1),1);
parfor ri=1:size(spike_has_comp,1)
    % load the vae latent 
    sID = spike_has_comp.syl_ID{ri};
    vae = h5read(fn_vae, ['/' sID]);
    vae = vae';
    % for each spike, determine what window index of vae 
    % assign to the first or last window if out of range
%     win_loc = spike_has.burst_loc{ri} + win_frame;
    win_loc = spike_has_comp.burst_loc{ri} + win_frame/2;
    [res, i_all] = ZZfunc_getVaeWithLag_v2(vae, win_loc, lags);
    d_vae1{ri} = res;
end
win_count1 = cellfun(@(x) size(x,1), d_vae1);
% stack into a large 3d matrix
d_vae1 = cat(1, d_vae1{:});
% get the start/end index of each syllable in this 3d matrix
cum_count1 = cumsum(win_count1);
spike_has_comp.vstart = [1; cum_count1(1:(end-1))+1];
spike_has_comp.vend = cum_count1;


d_vae2 = cell(size(spike_has_call,1),1);
parfor ri=1:size(spike_has_call,1)
    % load the vae latent 
    sID = spike_has_call.syl_ID{ri};
    vae = h5read(fn_vae, ['/' sID]);
    vae = vae';
    % for each spike, determine what window index of vae 
    % assign to the first or last window if out of range
%     win_loc = spike_has.burst_loc{ri} + win_frame;
    win_loc = spike_has_call.burst_loc{ri} + win_frame/2;
    [res, i_all] = ZZfunc_getVaeWithLag_v2(vae, win_loc, lags);
    d_vae2{ri} = res;
end
win_count2 = cellfun(@(x) size(x,1), d_vae2);
% stack into a large 3d matrix
d_vae2 = cat(1, d_vae2{:});
% get the start/end index of each syllable in this 3d matrix
cum_count2 = cumsum(win_count2);
spike_has_call.vstart = [1; cum_count2(1:(end-1))+1];
spike_has_call.vend = cum_count2;


%% 2. Calculate then plot acoustic distance
num_round = 20;  % how many round of sampling to take
dist_all = nan(length(lags), num_round);
for ti=1:length(lags)
  d_this1 = squeeze(d_vae1(:,ti,:));
  d_this2 = squeeze(d_vae2(:,ti,:));
  for round_i=1:num_round
    dist = ZZfunc_pairwiseDist_v2_twoData(d_this1, d_this2, 1000, 'cosine');
    dist_all(ti, round_i) = dist;
  end
end

% plot results
close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 500 500]);
ax = gca; hold(ax, 'on');
plot_mean_with_ci(ax, lags, dist_all', [0 0 0], 2, '-', [0.5 0.5 0.5], 0.25);
xlabel('Time lag from spike (ms)', 'FontSize', 14);
ylabel('Acoustic distance (cosine)', 'FontSize', 14);
fn_fig = fullfile(fd_save, sprintf('%s.%s.comp_dist.%s.pdf', birdID, neuronID, burst_method));
print(fig, fn_fig, '-dpdf', '-painters');

dist_mean = mean(dist_all, 2);
dist_mean_all(ni,:) = dist_mean';

dist_all_all{ni} = dist_all;


end
% save results
fn_dist = fullfile(fd_save, sprintf('%s.comp_dist.dist_all_all.mat', birdID));
save(fn_dist, 'dist_all_all');


%%  Plot all neurons as heatmap
% remove neurons with a lot of nan
sum_nan = sum(isnan(dist_mean_all), 2);
i_keep = find(sum_nan<=(size(dist_mean_all,2)/2));
dist_mean_all = dist_mean_all(i_keep, :);

close all; 
[fig, axs] = generatePanelGrid_v2(2, 1, [0.5 0.35], [0.05], [0.025;0.075], [0.15;0.05], 0.05, [0;0], [10 10 500 800]);
ax1 = axs(1); hold(ax1, 'on');
imagesc(ax1, lags, 1:size(info_neu,1), dist_mean_all, [0.75 1]);
% colormap(ax1, 'p');
xlim(ax1, [lags(1) lags(end)]);
ylim(ax1, [1 size(info_neu,1)]);
ylabel(ax1, 'Neuron ID', 'FontSize', 12);

ax2 = axs(2); cla(ax2); hold(ax2, 'on');
% mean_d = mean(dist_mean_all,1);
plot_mean_with_sem_nan(ax2, lags, dist_mean_all, [0 0 0], 2, '-', [0.5 0.5 0.5], 0.25);
% plot_median_with_sem(ax2, lags, dist_mean_all, [0 0 0], 2, '-', [0.5 0.5 0.5], 0.25);
xlabel('Time lag from spike (ms)', 'FontSize', 12);
ylabel('Acoustic distance (cosine)', 'FontSize', 12);

fn_fig = fullfile(fd_save, sprintf('%s.All.comp_dist.%s.pdf', birdID, burst_method));
print(fig, fn_fig, '-dpdf', '-painters');


%% plot as violin plot
% pre_t = [-75 -50]; pre_idx = find((lags>=pre_t(1)) & (lags<=pre_t(2)));
pre_t = [-50 -25]; pre_idx = find((lags>=pre_t(1)) & (lags<=pre_t(2)));
in_t = [10 35]; in_idx = find((lags>=in_t(1)) & (lags<=in_t(2)));
post_t = [75 100]; post_idx = find((lags>=post_t(1)) & (lags<=post_t(2)));

dist_array = cat(3, dist_all_all{:});
dist_m = squeeze(mean(dist_array, 2));
pre_v = mean(dist_m(pre_idx,:), 1);
in_v = mean(dist_m(in_idx,:), 1);
post_v = mean(dist_m(post_idx,:), 1);

% plot as violin
close all; 
fig = ZZfunc_newFigurePDFsize_v1([50 50 300 400]); 
% boxplot(all_data, group_labels);
d = [pre_v' in_v' post_v'];
violinplot(d);
fn_fig = fullfile(fd_save, sprintf('%s.All.comp_dist.%s.violin.pdf', birdID, burst_method));
print(fig, fn_fig, '-dpdf', '-painters');


%% perfom stats
% paired tests
[~, p1, ~, ~] = ttest(pre_v, in_v);
[~, p2, ~, ~] = ttest(in_v, post_v);

% signed-rank test
[p11,~,~] = signrank(pre_v, in_v);
[p22,~,~] = signrank(in_v, post_v);
end



%% Plot results from all birds
% load all results
res_all = [];
for bi=1:size(birdIDs,2)
  fd_res = fullfile(fd_base, 'Figures', pairIDs{bi}, 'CompoundSyl', birdIDs{bi}, 'embed', [vae_run '_distanceWinCenter_CallComp']);
  fn_res = fullfile(fd_res, sprintf('%s.comp_dist.dist_all_all.mat', birdIDs{bi}));
  a = load(fn_res);
  res_all = [res_all; a.dist_all_all];
end
dist_array = cat(3, res_all{:});
dist_m = squeeze(mean(dist_array, 2));

% plot
dist_mean_all = dist_m';
sum_nan = sum(isnan(dist_mean_all), 2);
i_keep = find(sum_nan<=(size(dist_mean_all,2)/2));
dist_mean_all = dist_mean_all(i_keep, :);

% randomize the order
i_r = randperm(size(dist_mean_all,1));
dist_mean_all = dist_mean_all(i_r,:);

close all; 
[fig, axs] = generatePanelGrid_v2(2, 1, [0.5 0.35], [0.05], [0.025;0.075], [0.15;0.05], 0.05, [0;0], [10 10 320 800]);
ax1 = axs(1); hold(ax1, 'on');
imagesc(ax1, lags, 1:size(dist_mean_all,1), dist_mean_all, [0.8 1.1]);
% colormap(ax1, 'p');
% xlim(ax1, [lags(1) lags(end)]);
xlim(ax1, [-100 100]);
ylim(ax1, [1 size(dist_mean_all,1)]);
ylabel(ax1, 'Neuron ID', 'FontSize', 12);
title(ax1, sprintf('All neurons (n=%d)', size(dist_mean_all,1)), 'FontSize', 12);

ax2 = axs(2); cla(ax2); hold(ax2, 'on');
% mean_d = mean(dist_mean_all,1);
plot_mean_with_sem_nan(ax2, lags, dist_mean_all, [0 0 0], 2, '-', [0.5 0.5 0.5], 0.25);
% plot_median_with_sem(ax2, lags, dist_mean_all, [0 0 0], 2, '-', [0.5 0.5 0.5], 0.25);
xlabel('Time lag from spike (ms)', 'FontSize', 12);
ylabel('Acoustic distance (cosine)', 'FontSize', 12);
xlim(ax2, [-100 100]);

fd_save = fullfile(fd_base, 'Figures', pairIDs{1}, 'CompoundSyl', birdIDs{1}, 'embed', [vae_run '_distanceWinCenter_CallComp']);
fn_fig = fullfile(fd_save, 'AllBirds.comp_distShorterLim.%s.pdf');
print(fig, fn_fig, '-dpdf', '-painters');
  
  
%%  statistical tests
pre_t = [-50 -25]; pre_idx = find((lags>=pre_t(1)) & (lags<=pre_t(2)));
in_t = [10 40]; in_idx = find((lags>=in_t(1)) & (lags<=in_t(2)));
post_t = [75 100]; post_idx = find((lags>=post_t(1)) & (lags<=post_t(2)));

pre_v = mean(dist_m(pre_idx,:), 1);
in_v = mean(dist_m(in_idx,:), 1);
post_v = mean(dist_m(post_idx,:), 1);

% paired t-tests
[~, p1, ~, ~] = ttest(pre_v, in_v);
[~, p2, ~, ~] = ttest(in_v, post_v);
[~, p3, ~, ~] = ttest(pre_v, post_v);

% signed-rank test
[p11,~,~] = signrank(pre_v, in_v);
[p22,~,~] = signrank(in_v, post_v);
[p33,~,~] = signrank(pre_v, post_v);
  

% where is the lowest point
[min_v, min_idx] = min(squeeze(mean(dist_mean_all,1)));
fprintf('Lowest point is at %d ms\n', lags(min_idx));
  
  
  
  
  



