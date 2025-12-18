% Plot the trajectory of calls from the shotgun-VAE/UMAP analysis
% Zhilei, 07/28/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_data = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'vaeWav');
birdID = 'pair5RigCCU29';
% what syllables to analyze
syl = {'v4', 'v5'};
% where is the UMAP embedding located
fd_umap = fullfile(fd_data, birdID, 'Traj', 'applySyl1', 'paramSearch2', 'v4v5.nsyl10000');
fn_umap = 'pair5RigCCU29.v4v5.nsyl10000.hop1ms.embedding.csv';
embed = readtable(fullfile(fd_umap, fn_umap), 'Delimiter', ',');


%% 1. Plot trajectories of call subtypes: each syllable is one line
% where to save plots
fd_save = fullfile(fd_umap, 'plotTraj1');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
% different colors for subtypes
col_list = {'#a65628','#4daf4a','#984ea3','#e41a1c','#ff7f00','#f781bf','#377eb8','#737373'};
col_alpha = 0.05;  % transparency of points in trajectory
end_alpha = 0.2;  % transparency of start/end points
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
syl_id = unique(embed.syl_id, 'sorted');  % get all syllable ids
% plot renditions
for si=1:size(syl_id, 1)
  % for si=1:10
  ss = syl_id{si};
  embed_s = embed(strcmp(embed.syl_id, ss), :);
  % sort by the window index
  embed_s = sortrows(embed_s, 'i_i', 'ascend');
  % get color
  v = embed_s.call_subtype{1};
  vi = find(strcmp(syl, v));
  col_hex = col_list{vi};
  col_rgb = sscanf(col_hex(2:end), '%2x%2x%2x', [1 3]) / 255;
  x = embed_s.umap1;
  y = embed_s.umap2;
  % plot as dots on lines
  plot(x, y, '-', 'Color', [col_rgb col_alpha], 'LineWidth', 1); hold on;
  scatter(x, y, 15, 'filled', 'MarkerFaceColor', col_rgb, 'MarkerFaceAlpha', col_alpha, 'MarkerEdgeColor', 'none'); hold on;
  % mark the start and end
  scatter(x(1), y(1), 40, 'Marker', '^', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
end
% fn_fig = fullfile(fd_save, sprintf('%s.%s.allRend.pdf', birdID, v));
% print(fig, fn_fig, '-dpdf', '-painters');
fn_fig = fullfile(fd_save, sprintf('%s.%s.allRend.jpg', birdID, v));
print(fig, fn_fig, '-djpeg', '-r300');


%% 2. Plot trajectories of call subtypes: don't connect windows as lines
% where to save plots
fd_save = fullfile(fd_umap, 'plotTraj1');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
% different colors for subtypes
col_list = {'#a65628','#4daf4a','#984ea3','#e41a1c','#ff7f00','#f781bf','#377eb8','#737373'};
col_alpha = 0.05;  % transparency of points in trajectory
end_alpha = 0.15;  % transparency of start/end points
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
syl_id = unique(embed.syl_id, 'sorted');  % get all syllable ids
% plot renditions
for si=1:size(syl_id, 1)
% for si=1:50
  % for si=1:10
  ss = syl_id{si};
  embed_s = embed(strcmp(embed.syl_id, ss), :);
  % sort by the window index
  embed_s = sortrows(embed_s, 'i_i', 'ascend');
  % get color
  v = embed_s.call_subtype{1};
  vi = find(strcmp(syl, v));
  col_hex = col_list{vi};
  col_rgb = sscanf(col_hex(2:end), '%2x%2x%2x', [1 3]) / 255;
  x = embed_s.umap1;
  y = embed_s.umap2;
  % plot as dots on lines
%   plot(x, y, '-', 'Color', [col_rgb col_alpha], 'LineWidth', 1); hold on;
  scatter(x, y, 15, 'filled', 'MarkerFaceColor', col_rgb, 'MarkerFaceAlpha', col_alpha, 'MarkerEdgeColor', 'none'); hold on;
  % mark the start and end
  scatter(x(1), y(1), 40, 'Marker', '^', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
  scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', end_alpha); hold on;
end
xlim([min(embed.umap1)-0.5 max(embed.umap1)+0.5]);
ylim([min(embed.umap2)-0.5 max(embed.umap2)+0.5]);
% add legend 
h_legend = gobjects(size(syl, 2), 1);
for i = 1:size(syl, 2)
    h_legend(i) = plot(nan, nan, 'o', 'Color', col_list{i}, 'MarkerFaceColor', col_list{i}, 'DisplayName', syl{i}, 'MarkerSize', 40);
end
legend(h_legend, syl, 'Location', 'best', 'FontSize', 12);
% fn_fig = fullfile(fd_save, sprintf('%s.%s.allRend.pdf', birdID, v));
% print(fig, fn_fig, '-dpdf', '-painters');
fn_fig = fullfile(fd_save, sprintf('%s.%s.allRend.noLine.jpg', birdID, v));
print(fig, fn_fig, '-djpeg', '-r300');















