% Plot the trajectory of calls from the shotgun-VAE/UMAP analysis
% Zhilei, 08/01/2025
% differ from v1: add codes to plot averaged trajectories

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_data = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'vaeWav');
birdID = 'pair5RigCCU29';
% what syllables to analyze
syl_all = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
% col_list = {'#db5f57', '#dbd057', '#75db57', '#57dbaa', '#579bdb', '#8557db', '#db57c0'};  %color used for call clustering
% col_list = {'#a65628','#4daf4a','#984ea3','#e41a1c','#ff7f00','#f781bf','#377eb8','#737373'}; 
col_list = {'#66c2a5','#e5c494','#8da0cb','#e78ac3','#a6d854','#ffd92f','#fc8d62'};
col_dict = struct;
for si=1:size(syl_all,2) col_dict.(syl_all{si}) = col_list{si}; end


%% 1. Plot trajectories: paired UMAP
% syl = {'v4', 'v5'};
syl = {'v1', 'v7'};
syl_str = strjoin(syl, '');
% where is the UMAP embedding located
umap_run = 'traj_chop_32_1_32';
fd_umap = fullfile(fd_data, birdID, 'Traj', 'applySyl5', 'paramSearch5_r2', umap_run);
% fn_umap = 'pair5RigCCU29.paired.v4v5.embedding.csv';
fn_umap = 'pair5RigCCU29.paired.v1v7.embedding.csv';
embed = readtable(fullfile(fd_umap, fn_umap), 'Delimiter', ',');

% where to save plots
fd_save = fullfile(fd_umap, sprintf('plotTraj_%s', syl_str));
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end

% create the scatter plot, then rasterized it to img, then replot with other element editable
close all;
fig_size = [10 10 600 600];
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; 
% determine the xy limits, all subplots use the same
marg = 1;
col_alpha = 0.05;  % transparency of points in trajectory
end_alpha = 0.15;  % transparency of start/end points
[img, ax] = ZZfunc_plotTraj_v1(ax, embed, col_dict, marg, col_alpha, end_alpha);
% save figure
fn_fig = fullfile(fd_save, sprintf('%s.%s.fig', birdID, syl_str));
savefig(fig, fn_fig); 
fn_pdf = fullfile(fd_save, sprintf('%s.%s.pdf', birdID, syl_str));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 2. Plot trajectories: all call UMAP, one call each plot
% Each call subtype one figure, but x/y limits are matched
syl_str = strjoin(syl, '');
% where is the UMAP embedding located
umap_run = 'traj_chop_32_1_32';
fd_umap = fullfile(fd_data, birdID, 'Traj', 'applySyl5', 'paramSearch5_r2', umap_run);
fn_umap = 'pair5RigCCU29.all.v0v1v2v3v4v5v6v7.embedding.csv';
embed = readtable(fullfile(fd_umap, fn_umap), 'Delimiter', ',');

% where to save plots
fd_save = fullfile(fd_umap, 'plotTrajAll');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end

% create the scatter plot, then rasterized it to img, then replot with other element editable
% determine the xy limits, all subplots use the same x/y limits
marg = 1;
col_alpha = 0.025;  % transparency of points in trajectory
end_alpha = 0.1;  % transparency of start/end points
% set a gloabal x/y lim
x_lim = [min(embed.umap1)-marg, max(embed.umap1)+marg];
y_lim = [min(embed.umap2)-marg, max(embed.umap2)+marg];
% loop through call subtypes
for si=1:size(syl_all, 2)
  v = syl_all{si};
  embed_this = embed(strcmp(embed.call_subtype, v), :);
  close all;
  fig_size = [10 10 600 600];
  fig = ZZfunc_newFigurePDFsize_v1(fig_size);
  ax = gca; 
  [img, ax] = ZZfunc_plotTraj_v2(ax, embed_this, col_dict, marg, col_alpha, end_alpha, x_lim, y_lim);
  % save figure
  fn_pdf = fullfile(fd_save, sprintf('%s.%s.pdf', birdID, v));
  print(fig, fn_pdf, '-dpdf', '-painters');
end


%% 3. Plot trajectories: all call UMAP, select calls in one plot
% Each call subtype one figure, but x/y limits are matched
syl_grp = {{'v4', 'v5'}, {'v1', 'v7'}, {'v4', 'v1'}, {'v4', 'v5', 'v1', 'v7'}};
% where is the UMAP embedding located
umap_run = 'traj_chop_32_1_32';
fd_umap = fullfile(fd_data, birdID, 'Traj', 'applySyl5', 'paramSearch5_r2', umap_run);
fn_umap = 'pair5RigCCU29.all.v0v1v2v3v4v5v6v7.embedding.csv';
embed = readtable(fullfile(fd_umap, fn_umap), 'Delimiter', ',');
% where to save plots
fd_save = fullfile(fd_umap, 'plotTrajAll');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end

% create the scatter plot, then rasterized it to img, then replot with other element editable
% determine the xy limits, all subplots use the same x/y limits
marg = 1;
col_alpha = 0.025;  % transparency of points in trajectory
end_alpha = 0.1;  % transparency of start/end points
% set a gloabal x/y lim
x_lim = [min(embed.umap1)-marg, max(embed.umap1)+marg];
y_lim = [min(embed.umap2)-marg, max(embed.umap2)+marg];
% loop through call subtypes
for si=1:size(syl_grp, 2)
  syl = syl_grp{si};
  syl_str = strjoin(syl, '');
  embed_this = embed(ismember(embed.call_subtype, syl), :);
  close all;
  fig_size = [10 10 600 600];
  fig = ZZfunc_newFigurePDFsize_v1(fig_size);
  ax = gca; 
  [img, ax] = ZZfunc_plotTraj_v2(ax, embed_this, col_dict, marg, col_alpha, end_alpha, x_lim, y_lim);
  % save figure
  fn_pdf = fullfile(fd_save, sprintf('%s.group.%s.pdf', birdID, syl_str));
  print(fig, fn_pdf, '-dpdf', '-painters');
end















