% Plot the mean trajectory of calls from the shotgun-VAE/UMAP analysis
% Zhilei, 08/02/2025

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


%% 1. Plot mean trajectories: paired UMAP
syl = {'v4', 'v5'};
% syl = {'v1', 'v7'};
% syl = syl_all;
syl_str = strjoin(syl, '');
% where is the UMAP embedding located
umap_run = 'traj_chop_32_1_32';
fd_umap = fullfile(fd_data, birdID, 'Traj', 'applySyl5', 'paramSearch5_r2', umap_run);
fn_umap = 'pair5RigCCU29.paired.v4v5.embedding.csv';
% fn_umap = 'pair5RigCCU29.paired.v1v7.embedding.csv';
% fn_umap = 'pair5RigCCU29.all.v0v1v2v3v4v5v6v7.embedding.csv';
embed = readtable(fullfile(fd_umap, fn_umap), 'Delimiter', ',');

% where to save plots
% fd_save = fullfile(fd_umap, sprintf('plotTraj_%s', syl_str));
fd_save = fullfile(fd_umap, 'plotTrajAll');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end

% create the scatter plot, then rasterized it to img, then replot with other element editable
close all;
fig_size = [10 10 600 600];
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca;
n_interp = 200;  % how many points to intepolate to 
% loop through selected syllables
h_legend = gobjects(length(syl), 1);
for vi=1:size(syl,2)
  v = syl{vi};
  embed_this = embed(strcmp(embed.call_subtype, v),:);
  % loop through syllable renditions then plot
  syl_id = unique(embed_this.syl_id, 'sorted');  % get all syllable ids
  interp_d = zeros(size(syl_id,1), n_interp, 2);
  for si=1:size(syl_id,1)
    ss = syl_id{si};
    embed_s = embed(strcmp(embed.syl_id, ss), :);
    % sort by the window index
    embed_s = sortrows(embed_s, 'i_i', 'ascend');
    d = embed_s(:, {'umap1', 'umap2'});
    d = table2array(d);
    % interpolate to standard dims
    d_interp = ZZfunc_interp_window(d, n_interp);
    interp_d(si,:,:) = d_interp;
  end
  % calculate mean and standard deviation
  mean_traj = squeeze(mean(interp_d, 1));
  std_traj = squeeze(std(interp_d, 1));
  % plot mean trajectory
  hold(ax, 'on');
  h_legend(vi) = plot(ax, mean_traj(:,1), mean_traj(:,2), 'Color', col_dict.(v), 'LineWidth', 4, 'DisplayName', syl{vi});
  % marker the start/end
  scatter(ax, mean_traj(1,1), mean_traj(1,2), 80, 'filled', 'MarkerFaceColor', [0 0 1]);
  scatter(ax, mean_traj(end,1), mean_traj(end,2), 80, 'filled', 'MarkerFaceColor', [0 0 0]);
  % mark
%   % convert color code to rgb
%   col_hex = col_dict.(v);
%   col_rgb = reshape(sscanf(col_hex(2:end), '%2x') / 255, 1, 3);
%   ax = GPT_plot_mean_std_3d(ax, interp_d, col_rgb, 3, col_rgb, 0);
end
legend(ax, h_legend, 'Location', 'southeast', 'FontSize', 16);
% determine x and y limits, same as scatter plots
marg = 1;
x_lim = [min(embed.umap1)-marg, max(embed.umap1)+marg];
y_lim = [min(embed.umap2)-marg, max(embed.umap2)+marg];
xlim(ax, x_lim);
ylim(ax, y_lim);
xlabel(ax, 'UMAP axis 1', 'FontSize', 16);
ylabel(ax, 'UMAP axis 2', 'FontSize', 16);
% save figure
fn_pdf = fullfile(fd_save, sprintf('%s.mean.%s.pdf', birdID, syl_str));
print(fig, fn_pdf, '-dpdf', '-painters');
  

















