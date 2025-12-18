% instead of plotting all spikes, calculate and plot the density of spikes
% select the highest density point to represent the neuron
% differ from v1: use the renamed call subtypes and stricter selection of sparse neurons

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced2');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.replaced2.mat', birdID)));
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% what syllable types to analyze
syls = {'v1', 'v2'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
umap_run = 'umapAll.v4v5';
fd_embed =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'popRaster2');
suffix_hahn = 'neuron_orderedPlotted7';
suffix_criteria = 'criteria7'; 
% where to save results and plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_plots2'], sprintf('density_%s_%s', syls{1}, syls{2}));
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
disp(fd_save);


%% 1. Load UMAP data, ephys spike location data, and selected sparse neuronID and order
u_all = cell(size(syls, 2), 1);
e_all = cell(size(syls, 2), 1);
n_all = cell(size(syls, 2), 1);
c_all = cell(size(syls, 2), 1);
% syl_i = 1;
for syl_i=1:size(syls, 2)
  ss = syls{syl_i};
  fn_u = fullfile(fd_embed, sprintf('%s.%s.%s.replaced2.mat', birdID, ss, umap_run));
  load(fn_u); u_all{syl_i} = umap;
  fn_e = fullfile(fd_embed, sprintf('%s.%s.sliding_loc.replaced2.mat', birdID, ss));
  load(fn_e); e_all{syl_i} = spike_embed;
  fn_n = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_hahn));
  load(fn_n); n_all{syl_i} = neuron_ordered;
  fn_c = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_criteria));
  load(fn_c); c_all{syl_i} = criteria;
end


%% 2. Create a background image to represent the trajectory divergence
% avoid plot the same rendition more than once
% plot renditions that happen in all neurons
% neuron_union = union(n_all{:});
neuron_all = cellfun(@(x) {x.neuronID}, e_all, 'UniformOutput', false);
neuron_all = unique(cat(2, neuron_all{:}), 'sorted');

fig_size = [10 10 600 600];
traj_color = {'#737373', '#737373'};
% use the same x/y limits for all plots
U = horzcat(u_all{:});   % concatenate all struct arrays
uvs = cat(1, U.umap);
x_lim = [min(uvs(:,1))-1, max(uvs(:,1))+1];
y_lim = [min(uvs(:,2))-1, max(uvs(:,2))+1];

close all;
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; hold(ax, 'on');
already_plot = {};
for ni=1:size(neuron_all, 2)
  for syl_i=1:size(syls,2)
    spike_embed = e_all{syl_i};
    umap = u_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuron_all{ni}));
    % loop through renditions, first plot trajectories
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      if ~ismember(spike_embed(si).syl_ID, already_plot)
        already_plot = [already_plot; spike_embed(si).syl_ID];
        u = umap(si).umap;
        scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', traj_color{syl_i}, 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none');
      end
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
fn_bg = fullfile(fd_save, sprintf('%s.%s_%s.bg.mat', birdID, syls{1}, syls{2}));
save(fn_bg, 'bg');


%% 3. Loop through neuron, calculate density
spike_col = {'#e78ac3', '#a6d854'};
% a custom colormap for spectrogram
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet
% how many bins when calculate density
num_bin = 100;
fn_bg = fullfile(fd_save, sprintf('%s.%s_%s.bg.mat', birdID, syls{1}, syls{2}));
load(fn_bg);

neuron_union = union(n_all{:}, 'sorted');
% record the max density info
dens_max = zeros(size(neuron_union,2), size(syls,2), 3);
% save density plot in a separate subfolder
fd_save_dens = fullfile(fd_save, 'dens_each'); 
if ~exist(fd_save_dens, 'dir'); mkdir(fd_save_dens); end;
for ni=1:size(neuron_union,2)
  % neuronID = '20240923-ch14';
  neuronID = neuron_union{ni};
  close all;
  [fig, axes] = generatePanelGrid_v2(2, 2, [0.4;0.4], [0.1], [0.05;0.05], [0.05;0.05], 0.05, [1 1], [10 10 800 850]);
  % syl_i = 1;
  for syl_i = 1:size(syls,2)
    umap = u_all{syl_i};
    spike_embed = e_all{syl_i};
    
    ax = axes(syl_i, 1); hold(ax, 'on');
    xlim(ax, bg.x_lim);
    ylim(ax, bg.y_lim);
    % add the scatter background
    image(ax, ax.XLim, ax.YLim, flipud(bg.faint_img));
    % plot where the spikes are
    seg_idx = find(strcmp({spike_embed.neuronID}, neuronID));
    % then overlay spikes
    x = []; y = [];
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      u = umap(si).umap;
      sp = spike_embed(si);
      mat_loc = sp.mat_loc;
      scatter(ax, u(mat_loc,1), u(mat_loc,2), 40, 'filled', 'MarkerFaceColor', spike_col{syl_i}, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none');
      x = [x; u(mat_loc,1)];
      y = [y; u(mat_loc,2)];
    end
    xlim(ax, bg.x_lim);
    ylim(ax, bg.y_lim);
    title(ax, sprintf('%s spikes', syls{syl_i}), 'FontSize', 14);
    
    % calculate the density then plot
    if ~isempty(x)
      ax = axes(syl_i, 2); hold(ax, 'on');
%       [density, x_edges, y_edges] = ZZfunc_calcSpikeDensity_v1(x, y, grid_size, bg.x_lim, bg.y_lim);
      [density, x_edges, y_edges, x_centers, y_centers] = ZZfunc_calcSpikeDensity_count_V1(x, y, num_bin, bg.x_lim, bg.y_lim, 3);
      density = density';
      imagesc(ax, x_centers, y_centers, density, [0, 1]);
      colormap turbo;
      % mark the higheset density point
      [max_val, linear_idx] = max(density(:));
      [row, col] = ind2sub(size(density), linear_idx);
      x_max = x_centers(col);   % because columns correspond to x
      y_max = y_centers(row);   % because rows correspond
      scatter(ax, x_max, y_max, 50, '+', 'MarkerEdgeColor', 'white');
      title(ax, sprintf('%s density max=%.3f', syls{syl_i}, max_val), 'FontSize', 14);

      % save density value
      dens_max(ni, syl_i, 1) = x_max;
      dens_max(ni, syl_i, 2) = y_max;
      dens_max(ni, syl_i, 3) = max_val;
    end
  end
  % export as pdf
  fn_fig = fullfile(fd_save_dens, sprintf('%s.n%d.%s.density.pdf', birdID, ni, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
end

% save the max density results
fn_dens = fullfile(fd_save, sprintf('%s.dens_max.mat', birdID));
save(fn_dens, 'dens_max');



%% 4. Plot peak density location of shared and specific neurons
% use marker types to differentiate call subtype
% add a small line to connect the same neuron
neuron_both = intersect(n_all{1}, n_all{2});
neuron_spe1 = setdiff(n_all{1}, n_all{2});
neuron_spe2 = setdiff(n_all{2}, n_all{1});
neuron_cat = {neuron_both; neuron_spe1; neuron_spe2};
title_strs = {'Both', 'v1-specific', 'v2-specific'};
cat_colors = {'#377eb8', '#e78ac3', '#a6d854'};

cut_off = 0.05;  % not considered as a field if max_dens is smaller
close all;
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; hold(ax, 'on');
xlim(ax, bg.x_lim);
ylim(ax, bg.y_lim);
% add the scatter background
image(ax, ax.XLim, ax.YLim, flipud(bg.faint_img));
% loop through neuron
for cat_i=1:size(neuron_cat, 1)
  neuron_ordered = neuron_cat{cat_i};
  for ni=1:size(neuron_ordered,2)
    n_idx = find(strcmp(neuron_union, neuron_ordered{ni}));
    if cat_i==1 || cat_i==2
      scatter(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), 40,  'o', 'MarkerFaceColor', cat_colors{cat_i}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      text(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), num2str(n_idx));
    end
    if cat_i==1 || cat_i==3
      scatter(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), 50,  '+', 'MarkerEdgeAlpha', 1, 'MarkerEdgeColor', cat_colors{cat_i}, 'LineWidth', 2);
      text(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), num2str(n_idx));
    end
    if cat_i==1
      plot(ax, [dens_max(n_idx, 1, 1) dens_max(n_idx, 2, 1)], [dens_max(n_idx, 1, 2) dens_max(n_idx, 2, 2)], 'LineStyle', '-', 'Color', cat_colors{cat_i});
    end
  end
end
xlim(ax, bg.x_lim);
ylim(ax, bg.y_lim);
% add legend
ax = ZZfunc_addSimpleLegend(ax, title_strs, cat_colors, 12, 0.75);


%% 5. Use stricter criteria to categorize a neuron into syllable-specific or shared 
min_rend = 15; % min number of renditions
c1 = struct2table(c_all{1}); 
c2 = struct2table(c_all{2}); 
cjoin = outerjoin(c1, c2, 'Keys','neuronID', 'MergeKeys',true);
cjoin = cjoin(ismember(cjoin.neuronID, neuron_union),:);


cjoin = cjoin((cjoin.num_rends_c1>=min_rend) & (cjoin.num_rends_c2>=min_rend), :);
% compute a mean time
cjoin.mean_t = nanmean([cjoin.psth_max_smooth_t_c1 cjoin.psth_max_smooth_t_c2], 2);
% sort by mean time
n_sort = sortrows(cjoin, 'mean_t');
neuron_sort = n_sort.neuronID;


%% 6. Change the style of the plot to emphasize co-localization of shared neurons
% use different colors to differentiate neuron categories
neuron_both = intersect(n_all{1}, n_all{2});
neuron_spe1 = setdiff(n_all{1}, n_all{2});
neuron_spe2 = setdiff(n_all{2}, n_all{1});
neuron_cat = {neuron_spe1; neuron_spe2; neuron_both};
title_strs = {'v1-specific', 'v2-specific', 'Shared-v1', 'Shared-v2'};
% cat_colors = {'#e78ac3', '#a6d854', };
% cat_colors = {'#da4ba2', '#8fcb2e', '#e41a1c', '#357933'};
cat_colors = {'#da4ba2', '#8fcb2e', '#ff7f00', '#377eb8'};

cut_off = 0.05;  % not considered as a field if max_dens is smaller
close all;
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; hold(ax, 'on');
xlim(ax, bg.x_lim);
ylim(ax, bg.y_lim);
% add the scatter background
image(ax, ax.XLim, ax.YLim, flipud(bg.faint_img));
% loop through neuron categories
for cat_i=1:size(neuron_cat, 1)
  neuron_ordered = neuron_cat{cat_i};
  for ni=1:size(neuron_ordered,2)
    n_idx = find(strcmp(neuron_union, neuron_ordered{ni}));
    if cat_i==1
      scatter(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), 100,  'o', 'MarkerFaceColor', cat_colors{1}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      text(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), num2str(n_idx));
    elseif cat_i==2
      scatter(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), 100,  'o', 'MarkerFaceColor', cat_colors{2}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      text(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), num2str(n_idx));
    else
      scatter(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), 100,  'o', 'MarkerFaceColor', cat_colors{3}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      scatter(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), 100,  'o', 'MarkerFaceColor', cat_colors{4}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      text(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), num2str(n_idx));
      text(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), num2str(n_idx));
    end
  end
end
xlim(ax, bg.x_lim);
ylim(ax, bg.y_lim);
% add legend
ax = ZZfunc_addSimpleLegend_v2(ax, title_strs, cat_colors, 14, zeros(size(title_strs,2),1)+13, linspace(-5,-8,size(title_strs,2)));
xlabel(ax, 'UMAP axis 1', 'FontSize', 14);
ylabel(ax, 'UMAP axis 2', 'FontSize', 14);

fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.peakLoc.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');




%% 7. Use lines to connect the same neuron
% not showing unconfident neurons
c1 = struct2table(c_all{1}); 
c2 = struct2table(c_all{2}); 
cjoin = outerjoin(c1, c2, 'Keys','neuronID', 'MergeKeys',true);
[~, new_order] = ismember(neuron_union, cjoin.neuronID);
cjoin = cjoin(new_order, :);
cjoin.dens1x = squeeze(dens_max(:, 1, 1)); cjoin.dens1y = squeeze(dens_max(:, 1, 2)); cjoin.dens1v = squeeze(dens_max(:, 1, 3)); 
cjoin.dens2x = squeeze(dens_max(:, 2, 1)); cjoin.dens2y = squeeze(dens_max(:, 2, 2)); cjoin.dens2v = squeeze(dens_max(:, 2, 3)); 

% filter neurons to get confident ones
min_rend = 15; % min number of renditions
min_dens = 0.2;
d = cjoin((cjoin.num_rends_c1>=min_rend) & (cjoin.num_rends_c2>=min_rend), :);

% use different colors to differentiate neuron categories
ni_both = find(d.isPass_c1 & d.isPass_c2);
ni_spe1 = find(d.isPass_c1 & (~d.isPass_c2));
ni_spe2 = find((~d.isPass_c1) & d.isPass_c2);
neuron_cat = {ni_spe1; ni_spe2; ni_both};
title_strs = {'v1-specific', 'v2-specific', 'Shared-v1', 'Shared-v2'};
cat_colors = {'#da4ba2', '#8fcb2e', '#377eb8'};

close all;
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; hold(ax, 'on');
xlim(ax, bg.x_lim);
ylim(ax, bg.y_lim);
% add the scatter background
image(ax, ax.XLim, ax.YLim, flipud(bg.faint_img));
% loop through neuron categories
for cat_i=1:size(neuron_cat, 1)
  ni_this = neuron_cat{cat_i};
  for ni=1:length(ni_this)
    n_idx = ni_this(ni);
    dd = d(n_idx,:);
    if cat_i==1
      scatter(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), 100,  'o', 'MarkerFaceColor', cat_colors{1}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      text(ax, dens_max(n_idx, 1, 1), dens_max(n_idx, 1, 2), num2str(n_idx));
    elseif cat_i==2
      scatter(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), 100,  'o', 'MarkerFaceColor', cat_colors{2}, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
      text(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), num2str(n_idx));
    else
      % only plot if both density pass threshold
      if (dd.dens1v>=min_dens) && (dd.dens2v>=min_dens)
      scatter(ax, dd.dens1x, dd.dens1y, 60,  'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', cat_colors{3}, 'LineWidth',2);
      scatter(ax, dd.dens2x, dd.dens2y, 60,  '+', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', cat_colors{3}, 'LineWidth',2);
      plot(ax, [dd.dens1x dd.dens2x], [dd.dens1y dd.dens2y], 'Color', 'red', 'LineWidth', 2);
      text(ax,  dd.dens1x, dd.dens1y, num2str(n_idx));
%       text(ax, dens_max(n_idx, 2, 1), dens_max(n_idx, 2, 2), num2str(n_idx));
      end
    end
  end
end
xlim(ax, bg.x_lim);
ylim(ax, bg.y_lim);
% add legend
ax = ZZfunc_addSimpleLegend_v2(ax, title_strs, cat_colors, 14, zeros(size(title_strs,2),1)+13, linspace(-5,-8,size(title_strs,2)));
xlabel(ax, 'UMAP axis 1', 'FontSize', 14);
ylabel(ax, 'UMAP axis 2', 'FontSize', 14);

fn_pdf = fullfile(fd_save, sprintf('%s.%s_%s.peakLoc.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');









