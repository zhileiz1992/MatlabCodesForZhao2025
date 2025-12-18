% calculate the acoutic distance before and after the identified primitive
% differ from v2: calculate for all prims identified from all birds

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% size of sliding window
win_size = 32;  % in unit of frames
sec_per_frame = 0.001;


% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


% bi = 1;
for bi=1:size(birdIDs,2)
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Bird specific settings
fs = 20000;
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre');
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);
% load the re-thresholding results
fn_xprob = fullfile(fd_save, 'compound_call3', 'combined', sprintf('%s.xprob_all.mat', birdID));
load(fn_xprob);

% load information about high-similarity strands
% vi = 3;
for vi=1:size(xprob_all,2)
fd_save_ref = fullfile(fd_save, 'compound_call3');
fn_strand = fullfile(fd_save_ref, sprintf('call_v%d.strands.mat', vi));
load(fn_strand);


%% Re-assess the primitive belongship based on new onset/offset
% determine if the strands belong to which primitive
thre_overlap = 20;
onsets = xprob_all(vi).onsets;
offsets = xprob_all(vi).offsets;
if isempty(onsets); continue; end

for strand_i=1:size(strands,2)
  % loop through identified primitives
  prim_belong = nan;
  for prim_i=1:length(onsets)
    % calculate overlap
    a = [onsets(prim_i) offsets(prim_i)];
    b = strands(strand_i).x_range;
    x_center = strands(strand_i).x_center;
    [overlapRatio, intersectionLength] = calculateOverlapRatio(a, b);
    if (intersectionLength>=thre_overlap) && (x_center>=a(1)) && (x_center<=a(2))
%     if (x_center>=a(1)) && (x_center<=a(2))
      prim_belong = prim_i;
      break
    end
  end
  strands(strand_i).prim_belong = prim_belong;
  % add the information regarding strands onset/offset
  strands(strand_i).onsets = onsets; 
  strands(strand_i).offsets = offsets; 
end


%% Loop through good primitives
% what primitive to look at
num_prim = length(onsets);
for pi=1:num_prim
% pi = 2;
  if xprob_all(vi).is_good(pi)==0
    continue;
  end
strands2 = strands([strands.prim_belong]==pi);
onset = strands2(1).onsets(pi); 
offset = strands2(1).offsets(pi); 
mean_dur = strands2(1).size_d(2);

% where to save results
fd_dist = fullfile(fd_save_ref, 'distance', sprintf('v%d', vi));
if ~exist(fd_dist, 'dir'); mkdir(fd_dist); end


%% 2. Compare acoustic distance before, during and after the primitive
% walk through different time points from the similarity band onset on the compound syllable
for ri=1:size(strands2, 2)
  tar_i = strands2(ri).tar_i;
  strands2(ri).dvae = comp.dvae{tar_i};
end

% how much to calculate before and after the primitive region
% pad = 50; 
pad = 100; 
t_range = -pad:(offset-onset+pad);
dist_all = nan(1, length(t_range));
for ti=1:length(t_range)
  t = t_range(ti);
  d_this = [];
  for ri=1:size(strands2,2)
    t_this = strands2(ri).y_range(1) + t;
    d = strands2(ri).dvae;
    if (t_this>=1) & (t_this<=size(d, 1))
      d_this = [d_this; d(t_this,:)];
    end
  end
  % calcualte pairwise distance
 dist = ZZfunc_pairwiseDist_v1(d_this, 1000, 'cosine');
 dist_all(ti) = dist;
end
% save for later use
dist_struct.t_range = t_range; 
dist_struct.dist_all = dist_all;
dist_struct.onset = onset;
dist_struct.offset = offset;
fn_struct = fullfile(fd_dist, sprintf('%s.distance.v%d.p%d.mat', birdID, vi, pi));
save(fn_struct, 'dist_struct');


%%  plot results
dist_all = dist_struct.dist_all;
t_range = dist_struct.t_range;
close all; 
col_pre = '#377eb8';
col_prim = '#e7298a';
col_post = '#e6ab02';
fig = ZZfunc_newFigurePDFsize_v1([10 10 350 600]);
ax = gca; hold(ax, 'on');
plot(ax, t_range(t_range<0), dist_all(t_range<0), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_pre);
plot(ax, t_range(t_range>=0 & t_range<=(offset-onset)), dist_all(t_range>=0 & t_range<=(offset-onset)), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_prim);
plot(ax, t_range(t_range>(offset-onset)), dist_all(t_range>(offset-onset)), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_post);
xline(0, 'LineStyle', '--', 'LineWidth', 1);
xline(offset-onset, 'LineStyle', '--', 'LineWidth', 1);
xlim(ax, [t_range(1) t_range(end)]);
xlabel('Time from prim. onset (ms)', 'FontSize', 12);
ylabel('Mean cosine distance', 'FontSize', 12);
fn_fig = fullfile(fd_dist, sprintf('%s.distance2.v%d.p%d.pdf', birdID, vi, pi));
print(fig, fn_fig, '-dpdf', '-painters');


end
  
  
end
  
end  

  


%% Combine results from all birds
res_all = [];
count = 0;
for bi=1:size(birdIDs,2)
  fns = dir(fullfile(fd_base, 'Figures', pairIDs{bi}, 'CompoundSyl', birdIDs{bi}, 'pairwiseCorre', 'compound_call3', 'distance', '*', '*.distance.*.mat'));
  for fi=1:size(fns,1)
    count = count+1;
    a = load(fullfile(fns(fi).folder, fns(fi).name));
    res_all(count).t_range = a.dist_struct.t_range;
    res_all(count).dist_all = a.dist_struct.dist_all;
    res_all(count).onset = a.dist_struct.onset;
    res_all(count).offset = a.dist_struct.offset;
    % calculate length of primitive
    res_all(count).dur = res_all(count).offset - res_all(count).onset + 1;
  end
end
  
% save to first bird
fd_save_plot = fullfile(fd_base, 'Figures', pairIDs{1}, 'CompoundSyl', birdIDs{1}, 'pairwiseCorre', 'compound_call3', 'distance');


%% plot raw traces in one figure
col_pre = '#377eb8';
col_prim = '#e7298a';
col_post = '#e6ab02';
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 400 600]);
ax=gca; hold(ax, 'on');
for ri=1:size(res_all,2)
  t_range = res_all(ri).t_range;
  dist_all = res_all(ri).dist_all;
  onset = res_all(ri).onset;
  offset = res_all(ri).offset;
  plot(ax, t_range(t_range<0), dist_all(t_range<0), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_pre);
  plot(ax, t_range(t_range>=0 & t_range<=(offset-onset)), dist_all(t_range>=0 & t_range<=(offset-onset)), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_prim);
  plot(ax, t_range(t_range>(offset-onset)), dist_all(t_range>(offset-onset)), 'LineStyle', '-', 'LineWidth', 2, 'Color', col_post);
end
fn_fig = fullfile(fd_save_plot, 'allBirds.allPrims.rawTrace.pdf');
print(fig, fn_fig, '-dpdf', '-painters');


%% plot two panels, one is aligned to primitive onset, another is align to primitive offset
lens = cellfun(@(x) length(x), {res_all.t_range});
[max_len, max_i] = max(lens);

close all; 
[fig, axs] = generatePanelGrid_v2(1, 2, [0.7], [], [0.1;0.05], [0.15;0.05], 0.05, [0], [10 10 800 500]);
t_cover = 35;   % how much within-primitive to cover for onset/offset aligned plot
% aligned to prim onset
ax1 = axs(1); hold(ax1, 'on'); cla(ax1);
idx = 1:(pad+t_cover);
t_all = cellfun(@(x) x(idx), {res_all.t_range}, 'UniformOutput', false);
rel_t = t_all{1};
v_all = cellfun(@(x) x(idx), {res_all.dist_all}, 'UniformOutput', false);
v_all = cat(1, v_all{:});
plot_mean_with_sem_nan(ax1, rel_t(1:pad), v_all(:, 1:pad), col_pre, 2, '-', hexToRGB(col_pre), 0.5)
plot_mean_with_sem_nan(ax1, rel_t((pad+1):end), v_all(:, (pad+1):end), col_prim, 2, '-', hexToRGB(col_prim), 0.5)
xline(ax1, 0, 'LineStyle', '--', 'LineWidth', 1, 'Color', '#737373');
xlim(ax1, [rel_t(1) rel_t(end)]);
ylim(ax1, [0.38 0.9]);
xlabel(ax1, 'T prim. onset (ms)');
ylabel(ax1, 'Mean acoustic dist.');
title(ax1, sprintf('All prim. (n=%d)', size(res_all,2)));

% aligned to prim offset
ax2 = axs(2); hold(ax2, 'on'); cla(ax2);
t_all2 = cellfun(@(x) x((end-pad-t_cover+1):end), {res_all.t_range}, 'UniformOutput', false);
v_all2 = cellfun(@(x) x((end-pad-t_cover+1):end), {res_all.dist_all}, 'UniformOutput', false);
v_all2 = cat(1, v_all2{:});
rel_t2 = -t_cover:(pad-1);
plot_mean_with_sem_nan(ax2, rel_t2(1:t_cover), v_all2(:, 1:t_cover), col_prim, 2, '-', hexToRGB(col_prim), 0.5)
plot_mean_with_sem_nan(ax2, rel_t2((t_cover+1):end), v_all2(:, (t_cover+1):end), col_post, 2, '-', hexToRGB(col_post), 0.5)
xline(ax2, 0, 'LineStyle', '--', 'LineWidth', 1, 'Color', '#737373');
xlim(ax2, [rel_t2(1) rel_t2(end)]);
ylim(ax2, [0.38 0.9]);
xlabel(ax2, 'T prim. offset (ms)');
% ylabel(ax2, 'Mean cosine distance');
fn_fig = fullfile(fd_save_plot, 'allBirds.allPrims.meanSEM.pdf');
print(fig, fn_fig, '-dpdf', '-painters');


