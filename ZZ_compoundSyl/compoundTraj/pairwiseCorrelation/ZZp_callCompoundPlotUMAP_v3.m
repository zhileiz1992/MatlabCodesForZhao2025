% plot the UMAP trajectories of primitives identified from call-compound comparison
% inputs to UMAP: renditions of calls and compound syllables that have the identified primitives
% differ from v3: for calls, use the time-warpped VAE as inputs for UMAP

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
clims = {[10 21]; [10 21]; [10.5 21.5]; [12 23]};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% size of sliding window
win_size = 32;  % in unit of frames
sec_per_frame = 0.001;


% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
clim = clims{bi};



%% 1. Load compound syllable data
fs = 20000;
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);
% load information about calls
fn_call = fullfile(fd_save, sprintf('%s.call.mat', birdID));
load(fn_call); 
fn_call_mean = fullfile(fd_save, sprintf('%s.call_mean.mat', birdID));
load(fn_call_mean);
% load information about high-similarity strands
vi = 1;
fd_save_ref = fullfile(fd_save, 'compound_call2');
fn_strand = fullfile(fd_save_ref, sprintf('call_v%d.strands.mat', vi));
load(fn_strand);
% what primitive to look at
pi = 1;
strands2 = strands([strands.prim_belong]==pi);
onset = strands2(1).onsets(pi); 
offset = strands2(1).offsets(pi); 
mean_dur = strands2(1).size_d(2);

% where to save results
fd_umap = fullfile(fd_save_ref, 'umap3', sprintf('v%d', vi));
if ~exist(fd_umap, 'dir'); mkdir(fd_umap); end


%% 2. Prepare data for UMAP
% how many call or compound syllables
inc_idx = unique([strands2.tar_i], 'sorted');  % one compound syllable may have several high-similarity strands
num_comp = length(inc_idx);
% grab the VAE data for call
calls_this = calls(strcmp(calls.newCallID, sprintf('v%d', vi)),:);
num_call = 3 * num_comp;   % include more renditions of calls, since compound syllables are generally longer
num_call = min([num_call size(calls_this,1)]);  % sample equal number
rng(1992);
calls_rd = calls_this(randsample(1:size(calls_this,1), num_call), :);
% time-wrapped VAE
mean_len = floor(mean(cellfun(@(x) size(x, 1), calls_rd.dvae)));
temp = cellfun(@(x) ZZfunc_interp_window(x, mean_len), calls_rd.dvae, 'UniformOutput', false);
calls_rd.vae_interp = temp;
% remove the first half_win and last half_win sliding window, which is mostly black
halfw = win_size / 2;
temp = cellfun(@(x) x((halfw+1):(size(x,1)-halfw),:), calls_rd.vae_interp, 'UniformOutput', false);
calls_rd.vae_chop = temp;
d_call = cat(1, calls_rd.vae_chop{:});

% grab the VAE data for compound syllable
comp_this = comp(inc_idx, :);
% remove the first half_win and last half_win sliding window, which is mostly black
temp = cellfun(@(x) x((halfw+1):(size(x,1)-halfw),:), comp_this.dvae, 'UniformOutput', false);
comp_this.vae_chop = temp;
d_comp = cat(1, comp_this.vae_chop{:});

% cat into a large data frame
d_all = [d_call; d_comp];
% d_all = d_call;
fn_d = fullfile(fd_umap, sprintf('%s.call_comp.v%d.p%d.d_all.csv', birdID, vi, pi));
writematrix(d_all, fn_d);


%% 3. Run UMAP
param.n_components = 2;
param.n_neighbors = 25; 
param.min_dist = 0;
% param.metric = 'cosine';
param.metric = 'euclidean';
param.random_state = 1118;
fn_save = fullfile(fd_umap, sprintf('%s.call_comp.v%d.p%d.umap_res.csv', birdID, vi, pi));
% also save the umap model 
fn_p = fullfile(fd_umap, sprintf('%s.call_comp.v%d.p%d.umap_model.p', birdID, vi, pi));
umap_res = ZZfunc_runUMAP_v1(fn_d, fn_save, fn_p, param);

% or read previous UMAP results
umap_res = readmatrix(fn_save);


%% 4. Add umap results to the data table
len_dcall = cellfun(@(x) size(x,1), calls_rd.vae_chop);
len_dcomp = cellfun(@(x) size(x,1), comp_this.vae_chop);
cum_dcall = cumsum(len_dcall);
cum_dcomp = cumsum(len_dcomp);
% loop through rows
for ri=1:size(calls_rd,1)
  if(ri==1); ustart=1; else; ustart=cum_dcall(ri-1)+1; end
  uend = cum_dcall(ri);
  calls_rd.umap_res{ri} = umap_res(ustart:uend,:);
end
for ri=1:size(comp_this,1)
  if(ri==1); ustart=sum(len_dcall)+1; else; ustart=sum(len_dcall)+cum_dcomp(ri-1)+1; end
  uend = sum(len_dcall)+cum_dcomp(ri);
  comp_this.umap_res{ri} = umap_res(ustart:uend,:);
end


%% 5.1 Plot trajectories of equal amount of call and compound syllable renditions
% colors for calls and compound syllables
col_call_bg = '#737373';  % color for the other regions in the calls
col_call_prim = '#66a61e'; % color for the primitive region in the calls
col_comp_bg = '#7570b3'; 
col_comp_prim = '#e7298a'; 

close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
% first plot trajectories of calls as background
to_plot_num = 150;
for ri=1:to_plot_num
% for ri=1:size(calls_rd,1)
% plot trajectory of one call
  u1 = calls_rd.umap_res{ri};
  % highlight area of where the primitive is
%   onset_this = floor(size(u1,1)/mean_dur * onset); 
%   offset_this = floor(size(u1,1)/mean_dur * offset);  
%   offset_this = min([offset_this size(u1,1)]);
  onset_this = onset;
  offset_this = offset;
  other_idx = setdiff(1:size(u1,1), onset_this:offset_this);
  scatter(ax, u1(other_idx,1), u1(other_idx,2), 10, 'filled', 'MarkerFaceColor', col_call_bg, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  scatter(ax, u1(onset_this:offset_this,1), u1(onset_this:offset_this,2), 10, 'filled', 'MarkerFaceColor', col_call_prim, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  
% plot trajectory of one compound syllable 
  u2 = comp_this.umap_res{ri};
  % highlight the area of where the primitive is
  syl_ID = comp_this.syl_ID{ri};
  tar_i = find(strcmp(comp.syl_ID, syl_ID));
  idx_syl = find(([strands2.tar_i]==tar_i) & ([strands2.prim_belong]==pi));
  % the same syllable may have two primitive occurence
  inside_idx = [];
  for ii=1:length(idx_syl)
    onset_this = strands2(idx_syl(ii)).y_range(1);
    offset_this = strands2(idx_syl(ii)).y_range(2);
    inside_idx = [inside_idx onset_this:offset_this];
  end
  other_idx = setdiff(1:size(u2,1), inside_idx);
  scatter(ax, u2(other_idx,1), u2(other_idx,2), 10, 'filled', 'MarkerFaceColor', col_comp_bg, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  scatter(ax, u2(inside_idx,1), u2(inside_idx,2), 10, 'filled', 'MarkerFaceColor', col_comp_prim, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
end
x_lim = [-6 16];
y_lim = [1.5 15];
ax = ZZfunc_rasterizeScatter_v1(ax, x_lim, y_lim);
ax = ZZfunc_addSimpleLegend_v2(ax, {'Call other', 'Call prim.', 'Comp. other', 'Comp. prim'}, {col_call_bg, col_call_prim, col_comp_bg, col_comp_prim}, 14, [-5 -5 -5 -5], linspace(4.5,2.5,4));
xlabel(ax, 'UMAP axis 1', 'FontSize', 12);
ylabel(ax, 'UMAP axis 2', 'FontSize', 12);
fn_fig = fullfile(fd_umap, sprintf('%s.call_compound.v%d.p%d.pdf', birdID, vi, pi));
print(fig, fn_fig, '-dpdf', '-painters');


%% 5.2 Plot as two separate figure: color pre and post primitive differently
col_pre = '#377eb8';
col_prim = '#e7298a';
col_post = '#e6ab02';

to_plot_num = 200;
x_lim = [-6 16];
y_lim = [1.5 15];

close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
for ri=1:to_plot_num
  u1 = calls_rd.umap_res{ri};
  scatter(ax, u1(1:(onset-1),1), u1(1:(onset-1),2), 10, 'filled', 'MarkerFaceColor', col_pre, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  scatter(ax, u1(onset:offset,1), u1(onset:offset,2), 10, 'filled', 'MarkerFaceColor', col_prim, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  scatter(ax, u1((offset+1):size(u1,1),1), u1((offset+1):size(u1,1),2), 10, 'filled', 'MarkerFaceColor', col_post, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
end
ax = ZZfunc_rasterizeScatter_v1(ax, x_lim, y_lim);
ax = ZZfunc_addSimpleLegend_v2(ax, {'Pre-prim.', 'Within-prim.', 'Post-prim.'}, {col_pre, col_prim, col_post}, 14, [-5 -5 -5], linspace(4.5,3,3));
xlabel(ax, 'UMAP axis 1', 'FontSize', 12);
ylabel(ax, 'UMAP axis 2', 'FontSize', 12);
fn_fig = fullfile(fd_umap, sprintf('%s.call_compound.separate1.v%d.p%d.pdf', birdID, vi, pi));
print(fig, fn_fig, '-dpdf', '-painters');


fig2 = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax2 = gca; hold(ax2, 'on');
for ri=1:to_plot_num
% plot trajectory of one compound syllable 
  u2 = comp_this.umap_res{ri};
  % highlight the area of where the primitive is
  syl_ID = comp_this.syl_ID{ri};
  tar_i = find(strcmp(comp.syl_ID, syl_ID));
  idx_syl = find(([strands2.tar_i]==tar_i) & ([strands2.prim_belong]==pi));
  % the same syllable may have multiple primitive occurence
  % sort strands by the location
  onset_all = []; offset_all=[];
  for ii=1:length(idx_syl)
    onset_all = [onset_all strands2(idx_syl(ii)).y_range(1)];
    offset_all = [offset_all strands2(idx_syl(ii)).y_range(2)];
  end
  [~, sort_i] = sort(onset_all);
  onset_sort = onset_all(sort_i);
  offset_sort = offset_all(sort_i);
  
  prev = 1; 
  for ii=1:length(onset_sort)
    onset_this = onset_sort(ii);
    offset_this = offset_sort(ii);
    scatter(ax2, u2(prev:(onset_this-1),1), u2(prev:(onset_this-1),2), 10, 'filled', 'MarkerFaceColor', col_pre, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
    scatter(ax2, u2(onset_this:offset_this,1), u2(onset_this:offset_this,2), 10, 'filled', 'MarkerFaceColor', col_prim, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
    prev = offset_this + 1;
  end
  if prev<size(u2,1)
    scatter(ax2, u2(prev:size(u2,1),1), u2(prev:size(u2,1),2), 10, 'filled', 'MarkerFaceColor', col_post, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  end
end
ax2 = ZZfunc_rasterizeScatter_v1(ax2, x_lim, y_lim);
ax2 = ZZfunc_addSimpleLegend_v2(ax2, {'Pre-prim.', 'Within-prim.', 'Post-prim.'}, {col_pre, col_prim, col_post}, 14, [-5 -5 -5], linspace(4.5,3,3));
xlabel(ax2, 'UMAP axis 1', 'FontSize', 12);
ylabel(ax2, 'UMAP axis 2', 'FontSize', 12);
fn_fig = fullfile(fd_umap, sprintf('%s.call_compound.separate2.v%d.p%d.pdf', birdID, vi, pi));
print(fig2, fn_fig, '-dpdf', '-painters');



%% 5.3 Plot trajectories: calls as background, then a few renditions of compound syllables
col_call_bg = '#969696';  % color for the other regions in the calls
col_call_prim = '#e7298a'; % color for the primitive region in the calls
col_comp_bg = '#7570b3'; 
close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
% first plot trajectories of calls as background
for ri=1:200
% for ri=1:size(calls_rd,1)
% plot trajectory of one call
  u1 = calls_rd.umap_res{ri};
  % highlight area of where the primitive is
  onset_this = floor(size(u1,1)/mean_dur * onset); 
  offset_this = floor(size(u1,1)/mean_dur * offset);  
  offset_this = min([offset_this size(u1,1)]);
  other_idx = setdiff(1:size(u1,1), onset_this:offset_this);
  scatter(ax, u1(other_idx,1), u1(other_idx,2), 10, 'filled', 'MarkerFaceColor', col_call_bg, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  scatter(ax, u1(onset_this:offset_this,1), u1(onset_this:offset_this,2), 10, 'filled', 'MarkerFaceColor', col_call_prim, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1);
  % also mark the start/end of trajectories
  scatter(ax, u1(1,1), u1(1,2), 10, 'filled', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.35);
  scatter(ax, u1(end,1), u1(end,2), 10, 'filled', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.35);
end
x_lim = [-6 16];
y_lim = [2 15.5];
ax = ZZfunc_rasterizeScatter_v1(ax, x_lim, y_lim);
ax = ZZfunc_addSimpleLegend_v2(ax, {'Call other', 'Call prim.'}, {col_call_bg, col_call_prim}, 14, [-5 -5], linspace(3.5,3,2));
xlabel(ax, 'UMAP axis 1', 'FontSize', 12);
ylabel(ax, 'UMAP axis 2', 'FontSize', 12);
% save the intermediate figure
fn_fig = fullfile(fd_umap, sprintf('%s.call_compound.rends_step1.v%d.p%d.pdf', birdID, vi, pi));
print(fig, fn_fig, '-dpdf', '-painters');

% plot trajectory of one compound syllable 
num_plot_comp = 5; 
col_comp = lines(num_plot_comp);
rng(1992);
% idx_plot = randsample(1:size(comp_this,1), num_plot_comp);
idx_plot = [7 20 188];
tar_i_list = {};
for rii=1:length(idx_plot)
  ri = idx_plot(rii);
  u2 = comp_this.umap_res{ri};
  % highlight the area of where the primitive is
  syl_ID = comp_this.syl_ID{ri};
  tar_i = find(strcmp(comp.syl_ID, syl_ID));
  tar_i_list = [tar_i_list; tar_i];
  idx_syl = find(([strands2.tar_i]==tar_i) & ([strands2.prim_belong]==pi));
  % the same syllable may have two primitive occurence
  inside_idx = [];
  for ii=1:length(idx_syl)
    onset_this = strands2(idx_syl(ii)).y_range(1);
    offset_this = strands2(idx_syl(ii)).y_range(2);
    inside_idx = [inside_idx onset_this:offset_this];
  end
  other_idx = setdiff(1:size(u2,1), inside_idx);
%   scatter(ax, u2(other_idx,1), u2(other_idx,2), 10, 'filled', 'MarkerFaceColor', col_comp(ri,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
  scatter(ax, u2(inside_idx,1), u2(inside_idx,2), 50, 'filled', 'MarkerFaceColor', col_comp(rii,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
%   plot(ax, u2(:,1), u2(:,2), 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', col_comp(ri,:), 'MarkerEdgeColor', 'none', 'LineStyle', '-', 'LineWidth', 1, 'Color', col_comp(ri,:));
  plot(ax, u2(:,1), u2(:,2), 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 2, 'Color', col_comp(rii,:));
  % marker the start as triangle
  scatter(ax, u2(1,1), u2(1,2), 80, 'Marker', '^', 'MarkerFaceColor', col_comp(rii,:), 'MarkerEdgeColor', 'none');
end
comp_str = cellfun(@(x) sprintf('comp. %d', x), tar_i_list, 'UniformOutput', false);

fn_fig = fullfile(fd_umap, sprintf('%s.call_compound.rends_step2.v%d.p%d.pdf', birdID, vi, pi));
print(fig, fn_fig, '-dpdf', '-painters');




