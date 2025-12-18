% for given syllable types, overlay spikes on the UMAP trajectory
% Zhilei, 09/08/2025
% differ from v2: use the renamed call subtypes and stricter selection of sparse neurons

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.replaced.mat', birdID)));
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% what syllable types to analyze
% syls = {'v4', 'v5'};
syls = {'v1', 'v2'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
umap_run = 'umapAll.v4v5';
fd_embed =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'popRaster');
suffix_hahn = 'neuron_orderedPlotted6';
% where to save results and plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_plots2']);



%% 1. Load UMAP data, ephys spike location data, and selected sparse neuronID and order
u_all = cell(size(syls, 2), 1);
e_all = cell(size(syls, 2), 1);
n_all = cell(size(syls, 2), 1);
% syl_i = 1;
for syl_i=1:size(syls, 2)
  ss = syls{syl_i};
  fn_u = fullfile(fd_embed, sprintf('%s.%s.%s.replaced.mat', birdID, ss, umap_run));
  load(fn_u); u_all{syl_i} = umap;
  fn_e = fullfile(fd_embed, sprintf('%s.%s.sliding_loc.replaced.mat', birdID, ss));
  load(fn_e); e_all{syl_i} = spike_embed;
  fn_n = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_hahn));
  load(fn_n); n_all{syl_i} = neuron_ordered;
end



%% 2. Choose a reference, overlay all neurons, then in a seperate panel, overlay spikes for the same neuron for a different syllable
% use the same color scheme
% in the two panels, could make spikes of one syllable faint, another standout
fd_save_this = fullfile(fd_save, sprintf('%s_%s_ind2', syls{1}, syls{2}));
if ~exist(fd_save_this, 'dir')
  mkdir(fd_save_this);
end
ref_i = 2;
target_i = 1;

% load the neuron order and color
neuron_ordered = n_all{ref_i};
A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#ffff33'};
N = length(neuron_ordered);
neuron_color = A(mod(0:N-1, numel(A))+1);
% set the color for trajectories as background
% traj_color = {'#e78ac3', '#a6d854'};
traj_color = {'#737373', '#737373'};
% set the alpha
alphas = [1 0.25];

% plot background trajectories first
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
U = horzcat(u_all{:});   % concatenate all struct arrays
uvs = cat(1, U.umap);
x_lim = [min(uvs(:,1))-1, max(uvs(:,1))+1];
y_lim = [min(uvs(:,2))-1, max(uvs(:,2))+1];
for ni=1:size(neuron_ordered, 2)
  for syl_i=1:size(syls,2)
    spike_embed = e_all{syl_i};
    umap = u_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuron_ordered{ni}));
    % loop through renditions, first plot trajectories
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', traj_color{syl_i}, 'MarkerFaceAlpha', 0.002, 'MarkerEdgeColor', 'none');
    end
  end
end
xlim(ax, x_lim);
ylim(ax, y_lim);
frame = getframe(ax); img = frame.cdata;
% make a faint image
alpha = 0.7;                                % 0 = original, 1 = white
white_img = 255 * ones(size(img), 'double'); % all white background
faint_img = (1-alpha) * double(img) + alpha * white_img;
faint_img = uint8(faint_img);
close;
% save the img for later use
bg.img=img; bg.faint_img=faint_img; bg.x_lim=x_lim; bg.y_lim=y_lim; bg.fig_size=[10 10 600 600];
fn_bg = fullfile(fd_save_this, sprintf('%s.bg.mat', birdID));
save(fn_bg, 'bg');

% then overlay the spikes
for syl_i=1:size(syls,2)
  umap = u_all{syl_i};
  spike_embed = e_all{syl_i};
  fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
  ax = gca; hold(ax, 'on');
  xlim(ax, x_lim);
  ylim(ax, y_lim);
  % add the scatter background
  image(ax, ax.XLim, ax.YLim, flipud(faint_img));
  for ni=1:size(neuron_ordered, 2)
    neuronID = neuron_ordered{ni};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % then overlay spikes
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
      scatter(ax, u(mat_loc,1), u(mat_loc,2), 20, 'filled', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
    end
  end
  xlim(ax, x_lim);
  ylim(ax, y_lim);
  title(ax, syls{syl_i}, 'FontSize', 14);
  fn_fig = fullfile(fd_save_this, sprintf('%s.overlayAll.%s.pdf', birdID, syls{syl_i}));
  print(fig, fn_fig, '-dpdf', '-painters');
end



%% 3. Differentiate 3 catetories: shared, V1-specific, V2-specific
neuron_both = intersect(n_all{1}, n_all{2});
neuron_spe1 = setdiff(n_all{1}, n_all{2});
neuron_spe2 = setdiff(n_all{2}, n_all{1});
neuron_cat = {neuron_both; neuron_spe1; neuron_spe2};
title_strs = {'Both', 'v1-specific', 'v2-specific'};

close all;
for cat_i=1:size(neuron_cat,1)
neuron_ordered = neuron_cat{cat_i};
A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#ffff33'};
N = length(neuron_ordered);
neuron_color = A(mod(0:N-1, numel(A))+1);
% plot each category in different figures
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
xlim(ax, x_lim);
ylim(ax, y_lim);
% add the scatter background
image(ax, ax.XLim, ax.YLim, flipud(faint_img));
for ni=1:size(neuron_ordered, 2)
  neuronID = neuron_ordered{ni};
  for syl_i=1:size(syls,2)
    spike_embed = e_all{syl_i};
    umap = u_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % then overlay spikes
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
%       scatter(ax, u(mat_loc,1), u(mat_loc,2), 20, 'filled', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
      if syl_i==1
        scatter(ax, u(mat_loc,1), u(mat_loc,2), 20,  'o', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
      else
        scatter(ax, u(mat_loc,1), u(mat_loc,2), 30,  '+', 'MarkerEdgeAlpha', 0.5, 'MarkerEdgeColor', neuron_color{ni}, 'LineWidth', 1);
      end
    end
  end
end
title(ax, title_strs{cat_i}, 'FontSize', 14);
xlim(ax, x_lim); ylim(ax, y_lim);
fn_fig = fullfile(fd_save_this, sprintf('%s.overlayAll.%s.pdf', birdID, title_strs{cat_i}));
print(fig, fn_fig, '-dpdf', '-painters');
end



%% 4. Same as no.3, but sample equal number of renditions
num_rends = 30;
neuron_both = intersect(n_all{1}, n_all{2});
neuron_spe1 = setdiff(n_all{1}, n_all{2});
neuron_spe2 = setdiff(n_all{2}, n_all{1});
neuron_cat = {neuron_both; neuron_spe1; neuron_spe2};
title_strs = {'Both', 'v1-specific', 'v2-specific'};

close all;
for cat_i=1:size(neuron_cat,1)
neuron_ordered = neuron_cat{cat_i};
A = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#f781bf','#377eb8','#ffff33'};
N = length(neuron_ordered);
neuron_color = A(mod(0:N-1, numel(A))+1);
% plot each category in different figures
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
ax = gca; hold(ax, 'on');
xlim(ax, x_lim);
ylim(ax, y_lim);
% add the scatter background
image(ax, ax.XLim, ax.YLim, flipud(faint_img));
for ni=1:size(neuron_ordered, 2)
  neuronID = neuron_ordered{ni};
  for syl_i=1:size(syls,2)
    spike_embed = e_all{syl_i};
    umap = u_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    act_rends = min([length(seg_idx) num_rends]);
    seg_idx = randsample(seg_idx, act_rends);
    % then overlay spikes
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
%       scatter(ax, u(mat_loc,1), u(mat_loc,2), 20, 'filled', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
      if syl_i==1
        scatter(ax, u(mat_loc,1), u(mat_loc,2), 20,  'o', 'MarkerFaceColor', neuron_color{ni}, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
      else
        scatter(ax, u(mat_loc,1), u(mat_loc,2), 30,  '+', 'MarkerEdgeAlpha', 0.5, 'MarkerEdgeColor', neuron_color{ni}, 'LineWidth', 1);
      end
    end
  end
end
title(ax, title_strs{cat_i}, 'FontSize', 14);
xlim(ax, x_lim); ylim(ax, y_lim);
fn_fig = fullfile(fd_save_this, sprintf('%s.overlaySampled%d.%s.pdf', birdID, num_rends, title_strs{cat_i}));
print(fig, fn_fig, '-dpdf', '-painters');
end
















