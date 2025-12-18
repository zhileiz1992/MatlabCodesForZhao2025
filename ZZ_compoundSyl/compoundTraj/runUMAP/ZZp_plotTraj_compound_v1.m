% plot the result from UMAP

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
% coloring for syllables
syl_full = {'v', 'h', 'e', 'b', 'x'};
col_full = {'#e41a1c', '#984ea3', '#4daf4a', '#377eb8', '#737373'};
col_dict = struct;
for si=1:size(syl_full,2) col_dict.(syl_full{si}) = col_full{si}; end
% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load results from UMAP
fd_umap_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'UMAPcomp');
% umap_run = 'b.x.n-1';
% syls = {'b', 'x'};
% umap_run = 'v.b.x.n-1';
% syls = {'v', 'b', 'x'};
% umap_run = 'v.b.h.e.x';
umap_run = 'v.b.h.e.x.n-1';
syls = {'v', 'b', 'h', 'e', 'x'};
fd_umap = fullfile(fd_umap_base, umap_run);

fn_info = fullfile(fd_umap, sprintf('%s.%s.info.csv', birdID, umap_run));
info = readtable(fn_info, 'Delimiter', ',');
fn_embed = fullfile(fd_umap, sprintf('%s.%s.embedding.csv', birdID, umap_run));
embed = readmatrix(fn_embed);

% where to save plots
fd_save = fullfile(fd_umap, 'trajPlots');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end


%% 2.1 Plot trajectories: all in one plot
close all;
fig_size = [10 10 600 600];
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; 
% determine the xy limits, all subplots use the same
marg = 0.5;
col_alpha = 0.02;  % transparency of points in trajectory, default 0.05
end_alpha = 0.2;  % transparency of start/end points, default 0.2
num_plot = 500;
% num_plot = 10000;
batches = {'batch1'};
r_seed = 1992;
% [img, ax] = ZZfunc_plotTraj_v5(ax, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, [], [], [], false);
% syl_to_plot = {'v', 'b', 'x'};
syl_to_plot = {'v', 'b', 'h', 'e', 'x'};
syl_str = strjoin(syl_to_plot, '');
xx_lim = [0 20]; 
yy_lim = [];
[img, ax, x_lim, y_lim] = ZZfunc_plotTraj_v5(ax, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, xx_lim, yy_lim, syl_to_plot, false);
title(ax, sprintf('%s n=%d', umap_run, num_plot),'FontSize', 14);

% save figure
fn_pdf = fullfile(fd_save, sprintf('%s.%s.n%d.pdf', birdID, syl_str, num_plot));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 2.1 Plot trajectories: one syllable per plot
% determine the xy limits, all subplots use the same
marg = 0.5;
col_alpha = 0.01;  % transparency of points in trajectory, default 0.05
end_alpha = 0.2;  % transparency of start/end points, default 0.2
num_plot = 1000;
% num_plot = 10000;
batches = {'batch1'};
r_seed = 1992;
% [img, ax] = ZZfunc_plotTraj_v5(ax, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, [], [], [], false);
% syl_to_plot = {'b', 'x'};
xx_lim = [0 20];
yy_lim = [];
for si=1:size(syls,2)
  close all;
  [fig, axes] = generatePanelGrid_v2(1, 2, [0.85], [], [0.05;0.05], [0.05;0.05], 0.05, [0], [10 10 1200 600]);
  % first scatter plot
  ax1 = axes(1);
%   [img, ax1, x_lim, y_lim, x_all, y_all] = ZZfunc_plotTraj_v5(ax1, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, [], [], syls(si), false);
  [img, ax1, x_lim, y_lim, x_all, y_all] = ZZfunc_plotTraj_v5(ax1, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, xx_lim, yy_lim, syls(si), false);
  title(ax1, sprintf('%s %s n=%d', umap_run, syls{si}, num_plot),'FontSize', 14);
  
  % then density plot
  ax2 = axes(2); cla(ax2);
  % Create a grid for evaluation
  nbins = 100;
  xx = linspace(x_lim(1), x_lim(2), nbins+1); 
  yy = linspace(y_lim(1), y_lim(2), nbins+1); 
  [N, edgesX, edgesY] = histcounts2(x_all, y_all, xx, yy, 'Normalization', 'pdf');
  imagesc(ax2, xx(1:nbins), yy(1:nbins), N');
  set(ax2, 'YDir', 'normal');
  colormap(custom_map);
  xlabel(ax2, 'UMAP axis 1', 'FontSize', 16);
  ylabel(ax2, 'UMAP axis 2', 'FontSize', 16);
  title(ax2, sprintf('%s %s n=%d', umap_run, syls{si}, num_plot),'FontSize', 14);
  
  % save figure
  fn_pdf = fullfile(fd_save, sprintf('%s.%s.n%d.pdf', birdID, syls{si}, num_plot));
  print(fig, fn_pdf, '-dpdf', '-painters');
end

% combined in one figure, different panels
close all;
[fig, axes] = generatePanelGrid_v2(1, 5, [0.8], [], [0.05;0.05], [0.05;0.05], 0.02, [1], [10 10 2000 400]);
num_plot = 500; 
for si=1:size(syls,2)
  ax1 = axes(si);
%   [img, ax1, x_lim, y_lim, x_all, y_all] = ZZfunc_plotTraj_v5(ax1, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, [], [], syls(si), false);
  [img, ax1, x_lim, y_lim, x_all, y_all] = ZZfunc_plotTraj_v5(ax1, embed, info, num_plot, batches, r_seed, col_dict, marg, col_alpha, end_alpha, xx_lim, yy_lim, syls(si), false);
  title(ax1, sprintf('%s %s n=%d', umap_run, syls{si}, num_plot),'FontSize', 12);
  if si>1; ylabel(ax1, '');end
end
 % save figure
fn_pdf = fullfile(fd_save, sprintf('%s.combined.n%d.pdf', birdID, num_plot));
print(fig, fn_pdf, '-dpdf', '-painters');



















