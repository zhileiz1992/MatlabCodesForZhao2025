% for given syllable types, overlay spikes on the UMAP trajectory
% Zhilei, 08/17/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedPull');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.pull.mat', birdID)));
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% where is the ephys data struct located,
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedPull');
% what syllable types to analyze
syls = {'v4', 'v5'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
umap_run = 'umapAll.v4v5';
fd_embed =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID);
suffix_hahn = 'neuron_orderedPlotted5';
% where to save results and plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_plots']);



%% 1. Load UMAP data, ephys spike location data, and selected sparse neuronID and order
u_all = cell(size(syls, 2), 1);
e_all = cell(size(syls, 2), 1);
n_all = cell(size(syls, 2), 1);
% syl_i = 1;
for syl_i=1:size(syls, 2)
  ss = syls{syl_i};
  fn_u = fullfile(fd_embed, sprintf('%s.%s.%s.pull.mat', birdID, ss, umap_run));
  load(fn_u); u_all{syl_i} = umap;
  fn_e = fullfile(fd_embed, sprintf('%s.%s.sliding_loc.pull.mat', birdID, ss));
  load(fn_e); e_all{syl_i} = spike_embed;
  fn_n = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_hahn));
  load(fn_n); n_all{syl_i} = neuron_ordered;
end


%% 2. Overlay spikes on trajectories: each neuron one plot
% get the batch information
syl_i = 2;
% save results to a subfolder
fd_save_this = fullfile(fd_save, syls{syl_i});
if ~exist(fd_save_this, 'dir')
  mkdir(fd_save_this);
end
umap = u_all{syl_i};
spike_embed = e_all{syl_i};
binfo = cellfun(@(x) x.batch{1}, {umap.umap_meta}, 'UniformOutput', false);
% set the same x/y limit for all neurons
uvs = cat(1, umap.umap);
% loop through neurons, use the same colors as in Hahnloser plots
neuron_ordered = n_all{syl_i};
A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#ffff33'};
N = length(neuron_ordered);
neuron_color = A(mod(0:N-1, numel(A))+1);

% loop through all neurons, saved as png
imgs = cell(size(neuron_ordered));
for ni=1:size(neuron_ordered, 2)
  % for ni=1:5
  close all;
  fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
  ax = gca; hold(ax, 'on');
  % neuronID = '20240912-ch11';
  neuronID = neuron_ordered{ni};
  seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
  % loop through renditions, first plot trajectories
  for sii=1:length(seg_idx)
    si = seg_idx(sii);
    u = umap(si).umap;
    sp = spike_embed(si);
    mat_loc = sp.mat_loc;
    scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', '#737373', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
  end
  % then overlay spikes
  for sii=1:length(seg_idx)
    si = seg_idx(sii);
    u = umap(si).umap;
    sp = spike_embed(si);
    mat_loc = sp.mat_loc;
    scatter(ax, u(mat_loc,1), u(mat_loc,2), 40, 'filled', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'none');
  end
  xlim(ax, [min(uvs(:,1))-1, max(uvs(:,1))+1]);
  ylim(ax, [min(uvs(:,2))-1, max(uvs(:,2))+1]);
  xlabel(ax, 'UMAP axis 1', 'FontSize', 14);
  ylabel(ax, 'UMAP axis 2', 'FontSize', 14);
  title(ax, sprintf('%s n%d %s %d', syls{syl_i}, ni, neuronID, length(seg_idx)), 'FontSize', 14, 'Color', neuron_color{ni});
  fn_fig = fullfile(fd_save_this, sprintf('%s.%s.n%d.%s.%s.png', birdID, syls{syl_i}, ni, neuronID, umap_run));
  % print(fig, fn_pdf, '-dpdf', '-painters');
  exportgraphics(fig, fn_fig, 'Resolution', 300);
  % save imgs to list
  frame = getframe(fig);
  im = frame2im(frame);
  [imind, cm] = rgb2ind(im, 256);
  imgs{ni}.frame = frame;
  imgs{ni}.imind = imind;
  imgs{ni}.cm = cm;
end
% also convert to gif and mp4
fn_gif = fullfile(fd_save_this,  sprintf('%s.%s.%s.gif', birdID, syls{syl_i}, umap_run));
for ni=size(imgs,2):-1:1
  if ni == size(imgs,2)
    imwrite(imgs{ni}.imind, imgs{ni}.cm, fn_gif, 'gif', 'Loopcount', inf, 'DelayTime', 1);
  else
    imwrite(imgs{ni}.imind, imgs{ni}.cm, fn_gif, 'gif', 'WriteMode', 'append', 'DelayTime', 1);
  end
end


%% 3. Overlay spikes: all neurons on the same plot
syl_i = 1;
% save results to a subfolder
fd_save_this = fullfile(fd_save, syls{syl_i});
umap = u_all{syl_i};
spike_embed = e_all{syl_i};
neuron_ordered = n_all{syl_i};
% set the same x/y limit for all neurons
uvs = cat(1, umap.umap);

% determine color to use
% option 1: loop through neurons, use the same colors as in Hahnloser plots
% A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#ffff33'};
% N = length(neuron_ordered);
% neuron_color = A(mod(0:N-1, numel(A))+1);
% option 2: use a rainbow color?
neuron_color = jet(size(neuron_ordered,2));
neuron_color = mat2cell(neuron_color, ones(size(neuron_ordered,2),1), 3);
neuron_color = flipud(neuron_color);

close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
for ni=1:size(neuron_ordered, 2)
  % for ni=1:5
  neuronID = neuron_ordered{ni};
  seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
  % loop through renditions, first plot trajectories
  for sii=1:length(seg_idx)
    si = seg_idx(sii);
    u = umap(si).umap;
    sp = spike_embed(si);
    mat_loc = sp.mat_loc;
    scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', '#737373', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeColor', 'none');
  end
end
for ni=1:size(neuron_ordered, 2)
  % for ni=1:5
  neuronID = neuron_ordered{ni};
  seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
  % then overlay spikes
  for sii=1:length(seg_idx)
    si = seg_idx(sii);
    u = umap(si).umap;
    sp = spike_embed(si);
    mat_loc = sp.mat_loc;
    scatter(ax, u(mat_loc,1), u(mat_loc,2), 20, 'filled', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
  end
end
xlim(ax, [min(uvs(:,1))-1, max(uvs(:,1))+1]);
ylim(ax, [min(uvs(:,2))-1, max(uvs(:,2))+1]);
xlabel(ax, 'UMAP axis 1', 'FontSize', 14);
ylabel(ax, 'UMAP axis 2', 'FontSize', 14);
title(ax, sprintf('%s n=%d', syls{syl_i}, size(neuron_ordered,2)), 'FontSize', 14);
fn_fig = fullfile(fd_save_this, sprintf('%s.%s.nAll.%s.rainbow.png', birdID, syls{syl_i}, umap_run));
% print(fig, fn_pdf, '-dpdf', '-painters');
exportgraphics(fig, fn_fig, 'Resolution', 300);


%% 4. Plot two call subtypes on the same figure: each neuron one plot
fd_save_this = fullfile(fd_save, sprintf('%s_%s_ind', syls{1}, syls{2}));
if ~exist(fd_save_this, 'dir')
  mkdir(fd_save_this);
end
neuron_ordered = unique([n_all{:}], 'stable');
neuron_color = jet(size(neuron_ordered,2));
neuron_color = mat2cell(neuron_color, ones(size(neuron_ordered,2),1), 3);
neuron_color = flipud(neuron_color);
% use different colors for the two call subtypes
col_list = {'#e78ac3', '#a6d854'};
neu_col_list = {'#e41a1c', '#3b8639'};
% the same axis limits for all plots
U = horzcat(u_all{:});            % concatenate all struct arrays
uvs = cat(1, U.umap);
for ni=1:size(neuron_ordered, 2)
  close all;
  fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
  ax = gca; hold(ax, 'on');
  % for ni=1:5
  neuronID = neuron_ordered{ni};
  for syl_i=1:size(syls,2)
    umap = u_all{syl_i};
    spike_embed = e_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % loop through renditions, first plot trajectories
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
      scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', col_list{syl_i}, 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor', 'none');
    end
  end
  
  for syl_i=1:size(syls,2)
    umap = u_all{syl_i};
    spike_embed = e_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % then overlay spikes
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
      scatter(ax, u(mat_loc,1), u(mat_loc,2), 40, '^', 'MarkerFaceColor', neu_col_list{syl_i}, 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
    end
  end
  xlim(ax, [min(uvs(:,1))-1, max(uvs(:,1))+1]);
  ylim(ax, [min(uvs(:,2))-1, max(uvs(:,2))+1]);
  xlabel(ax, 'UMAP axis 1', 'FontSize', 14);
  ylabel(ax, 'UMAP axis 2', 'FontSize', 14);
  title(ax, neuronID, 'FontSize', 14);
  fn_fig = fullfile(fd_save_this, sprintf('%s.n%d.%s.%s.png', birdID, ni, neuronID, umap_run));
  % print(fig, fn_pdf, '-dpdf', '-painters');
  exportgraphics(fig, fn_fig, 'Resolution', 300);
end



%% 5. Jesse's idea: select a few shared/specific neurons
% use outer ring color to show the neuronID, while inner ring color to show syllable type
% Actually, change to use different marker types for the two call subtypes
col_list = {'#e78ac3', '#a6d854'};  % colors of trajectories
% col_list = {'#bdbdbd', '#bdbdbd'};
col_inner = {'#e41a1c', '#3b8639'};
% col_outer = {'#377eb8', '#ff7f00', '#a65628', '#984ea3'};
col_outer = {'#ff7f00', '#a65628', '#377eb8', '#737373'};
neuron_pick = {'20241002-ch7', '20240912-ch8', '20240912-ch1', '20240915-ch4A'};

close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
U = horzcat(u_all{:});            % concatenate all struct arrays
uvs = cat(1, U.umap);
for ni=1:size(neuron_pick, 2)
  neuronID = neuron_pick{ni};
  for syl_i=1:size(syls,2)
    umap = u_all{syl_i};
    spike_embed = e_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % loop through renditions, first plot trajectories
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
%       if syl_i==1
%         scatter(ax, u(:,1), u(:,2), 10,  'o', 'MarkerFaceColor', col_list{syl_i}, 'MarkerFaceAlpha', 0.1, 'MarkerEdgeColor', 'none');
%       else
%         scatter(ax, u(:,1), u(:,2), 30,  '+', 'MarkerEdgeAlpha', 0.1, 'MarkerEdgeColor', col_list{syl_i}, 'LineWidth', 1);
%       end
      scatter(ax, u(:,1), u(:,2), 10,  'o', 'MarkerFaceColor', col_list{syl_i}, 'MarkerFaceAlpha', 0.02, 'MarkerEdgeColor', 'none');
    end
  end
end
x_lim = [min(uvs(:,1))-1, max(uvs(:,1))+1];
y_lim = [min(uvs(:,2))-1, max(uvs(:,2))+1];
ax = ZZfunc_rasterizePlot_v1(ax, x_lim, y_lim);
for ni=1:size(neuron_pick, 2)
  neuronID = neuron_pick{ni};
  for syl_i=1:size(syls,2)
    umap = u_all{syl_i};
    spike_embed = e_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % then overlay spikes
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
      %       scatter(ax, u(mat_loc,1), u(mat_loc,2), 40, 'filled', 'MarkerFaceColor', col_inner{syl_i}, 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', col_outer{ni}, 'LineWidth', 2);
      %       scatter(ax, u(mat_loc,1), u(mat_loc,2), 40, 'filled', 'MarkerFaceColor', col_outer{ni}, 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', col_inner{syl_i} , 'LineWidth', 2);
      if syl_i==1
        scatter(ax, u(mat_loc,1), u(mat_loc,2), 50, 'o', 'MarkerFaceColor', col_outer{ni}, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
      else
        scatter(ax, u(mat_loc,1), u(mat_loc,2), 70, '+', 'MarkerFaceColor', col_outer{ni}, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', col_outer{ni}, 'LineWidth', 2);
      end
    end
  end
end
% add legend
for ni=1:size(neuron_pick,2)
  text(mean(ax.XLim), ax.YLim(1)+ni*0.75, neuron_pick{ni}, 'Color', col_outer{ni}, 'FontSize', 8);
end
xlabel(ax, 'UMAP axis 1', 'FontSize', 14);
ylabel(ax, 'UMAP axis 2', 'FontSize', 14);
title(ax, sprintf('%s vs %s', syls{1}, syls{2}), 'FontSize', 14);
fn_fig = fullfile(fd_save_this, sprintf('%s.JesseSelected.%s.pdf', birdID, umap_run));
print(fig, fn_fig, '-dpdf', '-painters');
% exportgraphics(fig, fn_fig, 'Resolution', 300);















