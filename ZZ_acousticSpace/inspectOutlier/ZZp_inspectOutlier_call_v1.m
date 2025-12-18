% Inspect why there are outliers in the call trajectory analysis
% Zhilei, 07/21/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_data = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed', 'vaeWav');
birdID = 'pair5RigCCU29';
% where is the UMAP embedding located
fd_umap = fullfile(fd_data, birdID, 'Traj', 'applySyl1', 'callUMAPone');
% how many syllables used in UMAP
umap_in = 'rd300';
% where is the original VAE latent located
fd_latent = fullfile(fd_data, birdID, 'Traj', 'applySyl1', 'latent.traj_chop_32_2_32');


%% 1. Plot trajectories and identify outliers
v = 'v4';
fn_embed = fullfile(fd_umap, sprintf('%s.%s.%s.embedding.csv', birdID, v, umap_in));
embed = readtable(fn_embed);

% where to save plots
fd_save = fullfile(fd_umap, 'inspectOutlier', v);
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
% plot all syllables
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
syl_id = unique(embed.syl_id, 'sorted');
for si=1:size(syl_id, 1)
% for si=1:10
  ss = syl_id{si};
  embed_s = embed(strcmp(embed.syl_id, ss), :);
  % sort by the window index
  embed_s = sortrows(embed_s, 'i_i', 'ascend');
  x = embed_s.umap1;
  y = embed_s.umap2;
  % plot as dots on lines
  plot(x, y, '-', 'Color', [0.25 0.25 0.25 0.1], 'LineWidth', 1); hold on;
  scatter(x, y, 15, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', 'none'); hold on;
  % mark the start and end
  scatter(x(1), y(1), 40, 'Marker', '^', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none'); hold on;
  scatter(x(end), y(end), 40, 'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none'); hold on;
end
fn_fig = fullfile(fd_save, sprintf('%s.%s.allRend.pdf', birdID, v));
print(fig, fn_fig, '-dpdf', '-painters');

% plot each rendtion in one figure panel 
close all;
b_starts = 1:24:size(syl_id,1);
xlim_all = [min(embed.umap1)-0.5, max(embed.umap1)+0.5];
ylim_all = [min(embed.umap2)-0.5, max(embed.umap2)+0.5];
for b_s=1:(length(b_starts)-1)
  [fig, axes] = generatePanelGrid_v2(4, 6, [0.22;0.22;0.22;0.22], [0.02;0.02;0.02], [0.02;0.02], [0.1;0.05], 0.02, [1;1;1;1], [10 10 1800 900]);
  for ii=0:23
    si = b_starts(b_s)+ii;
    plot_i = idivide(int32(ii), 6)+1;
    plot_j = mod(ii, 6)+1;
    ax = axes(plot_i, plot_j);
    hold(ax, 'on');
    ss = syl_id{si};
    embed_s = embed(strcmp(embed.syl_id, ss), :);
    % sort by the window index
    embed_s = sortrows(embed_s, 'i_i', 'ascend');
    x = embed_s.umap1;
    y = embed_s.umap2;
    % plot as dots on lines
    plot(ax, x, y, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1); 
    scatter(ax, x, y, 20, 'filled', 'MarkerFaceColor', [0.25 0.25 0.25], 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none'); 
    % mark the start and end
    scatter(ax, x(1), y(1), 50, 'Marker', '^', 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', 'none'); 
    scatter(ax, x(end), y(end), 50, 'Marker', 's', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
    xlim(ax, xlim_all);
    ylim(ax, ylim_all);
    set(ax, 'XTickLabel', [], 'YTickLabel', []);
    title(ax, sprintf('%s %d', v, si));
  end
  % save 
  fn_fig = fullfile(fd_save, sprintf('%s.%s.eachRend.%d.pdf', birdID, v, b_s));
  print(fig, fn_fig, '-dpdf', '-painters');
end
    





















