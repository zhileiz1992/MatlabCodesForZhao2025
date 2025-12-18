% given call clustering results, plot spectrograms to check clustering quality
% differ from v5: swap the call subtype names; change colors as well

clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
fd_save_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
% what spectrogram input dataset
input_rn = 'spec_goffinet_nn_256_176';
% input_rn = 'spec_goffinet_cutoff_256_176';
% default how to swap names
% bird M1
% oldID = {'v4', 'v5', 'v1', 'v7', 'v3', 'v6', 'v2'};
% newID = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
% % bird M2
% oldID = {'v1', 'v4', 'v3', 'v5', 'v6', 'v2'};
% newID = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'};
% % bird M3
% oldID = {'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v1'};
% newID = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
% bird M4
oldID = {'v1', 'v7', 'v6', 'v2', 'v3', 'v4', 'v5'};
newID = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
col_list = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#ccad25'};
col_rand = {'#252525','#8c510a','#01665e','#252525','#8c510a','#01665e','#252525'};  % use grey for random data
for vi=1:size(newID,2); swapID.(oldID{vi}) = newID{vi}; colDict.(newID{vi})=col_list{vi}; colRand.(newID{vi})=col_rand{vi}; end

syl = 'v';  % what syllable to focus on
% what subfolder has the VAE/UMAP results
% UMAP on VAE
fd_vae = 'UMAPonVAE7'; vae_run = input_rn;


% loop through birds
bi = 4;
bd = birdIDs{bi};
% where to save results
fd_save = fullfile(fd_save_base, bd, fd_vae, syl, vae_run, 'plotsNew');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
disp(fd_save);


%% 1. Plot the embedding data for both real and shuffled datasets
% load the embedding file
fn_embed = fullfile(fd_save_base, bd, fd_vae, syl, input_rn, sprintf('%s.%s.embedding.csv', bd, vae_run));
embed = readtable(fn_embed);
% also read the shuffled dataset
fn_embed_r = fullfile(fd_save_base, bd, fd_vae, syl, input_rn, sprintf('%s.random.embedding.csv', bd));
embed_r = readtable(fn_embed_r);
% locate the spectrogram file
fn_spec = fullfile(fd_save_base, bd, 'Spectrogram2', syl, sprintf('%s.%s.%s.h5', bd, syl, input_rn));
disp(fn_spec);
% read in color txt file
fn_col = fullfile(fd_save_base, bd, fd_vae, syl, input_rn, 'color.txt');
col_list = readtable(fn_col);
% plot embedding to double check, loop through real and shuffled dataset
ds = {embed; embed_r};
suffix = {'real'; 'random'};
for di=1:size(ds,1)
  d = ds{di};
  if di==1; colThis=colDict; else; colThis=colRand;end
  alpha = 0.5;
  x = d.umap1;
  y = d.umap2;
  clusters = d.hdbscan_cluster;
  % Create figure
  close all;
  fig = ZZfunc_newFigurePDFsize_v1([50 50 600 600]);
  hold on;
  % Plot noise cluster (0) in grey
  mask_noise = clusters == 0;
  scatter(x(mask_noise), y(mask_noise), 6, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', alpha);
  % Plot other clusters
  valid_clusters = unique(clusters(clusters > 0));
  for i = 1:numel(valid_clusters)
    cluster_id = valid_clusters(i);
    oldV = sprintf('v%d', cluster_id);
    if di==1; newV = swapID.(oldV); else; newV = oldV; end
    color = colThis.(newV);
    mask = clusters == cluster_id;
    scatter(x(mask), y(mask), 6, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
    % add cluster label
    %   text(mean(x(mask)), mean(y(mask)), num2str(cluster_id), 'FontSize', 16);
  end
  xlim_old = [min(x)-0.5 max(x)+0.5];
  ylim_old = [min(y)-0.5 max(y)+0.5];
  xlim(xlim_old); ylim(ylim_old);
  % rasterized the scatter for fast loading in Illustrator
  % convert to image
  ax = gca;
  frame = getframe(ax);
  img = frame.cdata;
  cla(ax);
  xlim(ax, xlim_old); ylim(ax, ylim_old);
  % add the scatter background
  image(ax, ax.XLim, ax.YLim, flipud(img));
  % add cluster labels
  for i = 1:numel(valid_clusters)
    cluster_id = valid_clusters(i);
    mask = clusters == cluster_id;
    oldV = sprintf('v%d', cluster_id);
     if di==1; newV = swapID.(oldV); else; newV = oldV; end
    text(ax, mean(x(mask)), mean(y(mask)), newV, 'FontSize', 18);
  end
  xlabel('UMAP axis 1');
  ylabel('UMAP axis 2');
  % add the number of vocalizations into the plot
  text(ax, mean(xlim_old), ylim_old(2), sprintf('N=%d calls', size(embed,1)), 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'top', 'FontSize', 16);
  title(bd, 'FontSize', 16);
  fn_fig = fullfile(fd_save, sprintf('%s.%s.embedding.pdf', bd, suffix{di}));
  print(fig, fn_fig, '-dpdf', '-painters');
end


%% 1.2 Color embedding by experimental date
% add a date field to the table
second_last = cellfun(@(s) strsplit(s, '/'), embed.fn_wav, 'UniformOutput', false);
embed.date = cellfun(@(parts) parts{end-1}, second_last, 'UniformOutput', false);
% add a field for relative days
date_dt = datetime(embed.date, 'InputFormat', 'yyyy-MM-dd');
day0 = min(date_dt);
embed.rel_days = days(date_dt - day0)+1;
% plot embedding, colored by relative days
close all;
fig = ZZfunc_newFigurePDFsize_v1([50 50 700 600]);
scatter(embed.umap1, embed.umap2, 6, embed.rel_days, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
colormap jet;
% caxis([0 60]);
cb = colorbar; cb.Label.String = 'Relative days'; cb.Label.FontSize=14;
xlabel('UMAP axis 1'); ylabel('UMAP axis 2');
xlim([min(embed.umap1)-0.5 max(embed.umap1)+0.5]);
ylim([min(embed.umap2)-0.5 max(embed.umap2)+0.5]);
title(bd, 'FontSize', 16);
fn_fig = fullfile(fd_save, sprintf('%s.%s.embeddingByDate.pdf', bd, suffix{di}));
print(fig, fn_fig, '-dpdf');


%% 2. Plot spectrograms
close all;
% a custom colormap for spectrogram
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet
% get the frequency list, assume same across all calls
s = embed{1, 'spec_f'}{1};
s = erase(s, {'[ ', ']'});
s = split(s);
freq = [];
for i=1:(size(s,1)-1)
  freq = [freq str2num(s{i})];
end

rng(1992);
% loop through cluster ids
% cid = 4;
for cid = 1:max(embed.hdbscan_cluster)
  idx = find(embed.hdbscan_cluster==cid);
  % randomly sample
  idx_rd = randsample(idx, 50);
  fig_pos = [10 10 1300 950];
  [fig, axes] = generatePanelGrid_v2(5, 10, [0.15;0.15;0.15;0.15;0.15], [0.02;0.02;0.02;0.02], [0.08;0.02], [0.1;0.05], 0.01, [0;0;0;0;0], fig_pos);
  count = 0;
  for ii=1:size(axes,1)
    for jj=1:size(axes,2)
      count = count+1;
      ax = axes(ii,jj);
      spec = h5read(fn_spec, '/spec_win_all', [1, 1, idx_rd(count)], [Inf, Inf, 1]);
      imagesc(ax, 1:size(spec',1), freq, spec', [0, 1]);
      set(ax, 'YDir', 'Normal');
      colormap(custom_map);
      if jj==1
        set(ax, 'YDir', 'normal');
        set(ax, 'YTick', [1000 3000 5000 7000]);
        set(ax, 'YTickLabel', {'1k', '3k', '5k', '7k'}, 'FontSize', 10);
      else
        set(ax, 'YTick', []);
      end
      set(ax, 'XTick', []);
      % show the HDBSCAN score
      score = embed{idx_rd(count), 'hdbscan_prob'};
      title(ax, sprintf('%d %.3f', idx_rd(count), score), 'FontSize', 8);
      
    end
  end
  sgtitle(sprintf('Cluster %d', cid), 'FontSize', 16);
  % save figure
  fn_pdf = fullfile(fd_save, sprintf('%s.spect.cluster%d.pdf', bd, cid));
  print(fig, fn_pdf, '-dpdf', '-painters');
end


%% 3. Create avereaged spectrogram: with time wrapping
% read sound directly from source nc file
fs = 20000;
% add a little before and after
pad = 200;  % 10 ms
mean_all = {};
t_all = {};
for cid = 1:max(embed.hdbscan_cluster)
  idx = find(embed.hdbscan_cluster==cid);
  % read sound data in
  fns = embed.fn_wav;
  istart = embed.istart;
  iend = embed.iend;
  segs = {};
  parfor ii=1:length(idx)
    si = idx(ii);
    fn = fns{si};
    [sound, fs] = audioread(fn);
    sound = ZZ_butterFilterFunc_v1(sound, fs, 250, 7500, 5);
    i_start = max([1 istart(si)-pad]);
    i_end = min([size(sound,1) iend(si)+pad]);
    seg = sound(i_start:i_end);
    segs{ii} = seg;
  end
  % align sound to the mean duration
  lens = cellfun(@length, segs);
  % Compute mean length
  m = ceil(mean(lens));
  spec_aligned = {};
  parfor ii=1:size(segs,2)
    n = length(segs{ii});
    seg_aligned = interp1(linspace(1, n, n), segs{ii}, linspace(1, n, m), 'linear');
    % calculate spectrogram
    [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(seg_aligned, fs, 256, 256, 236, [250 7500], [10 21]);
    spec_aligned{ii} = powerGrey/1024;  % normalize
  end
  stacked_array = cat(3, spec_aligned{:});
  spec_mean = mean(stacked_array, 3);
  mean_all{cid} = spec_mean;
  % save the time information
  n = length(segs{1});
  seg_aligned = interp1(linspace(1, n, n), segs{1}, linspace(1, n, m), 'linear');
  % calculate spectrogram
  [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(seg_aligned, fs, 256, 256, 236, [250 7500], [10 21]);
  t_all{cid} = t;
end

% plot with different scaling
close all;
fig_pos = [100 100 1600 250];
[fig, axes] = generatePanelGrid_v2(1, 10, [0.75], [], [0.15;0.02], [0.05;0.05], 0.02, [0], fig_pos);% loop through each cluster
for cid = 1:max(embed.hdbscan_cluster)
  oldV = sprintf('v%d', cid);
  newV = swapID.(oldV);
  plot_i = str2num(strrep(newV, 'v', ''));
  ax = axes(plot_i);
  spec_mean = mean_all{cid};
  idx = find(embed.hdbscan_cluster==cid);
  %   imagesc(ax, 1:size(spec_mean,1), freq, spec_mean); colormap jet;
  %   imagesc(ax, 1:size(spec_mean,1), freq, spec_mean, [0.15 max(spec_mean(:))]); colormap(custom_map);
  rel_t = t_all{cid};
  rel_t = rel_t - rel_t(1);
  imagesc(ax, rel_t, freq, spec_mean, [0.15 max(spec_mean(:))]); colormap(custom_map);
  set(ax, 'YDir', 'Normal');
  %   colormap jet;
  set(ax, 'YDir', 'normal');
  set(ax, 'YTick', [1000 3000 5000 7000]);
  set(ax, 'YTickLabel', {'1k', '3k', '5k', '7k'}, 'FontSize', 10);
  %   set(ax, 'XTick', []);
  set(ax, 'XTick', [0 max(rel_t)]);
  %   xticklabels = arrayfun(@(x) sprintf('%.2f', x), [0 max(rel_t)], 'UniformOutput', false);
  title(ax, newV, 'FontSize', 14, 'Color', colDict.(newV));
  xlim_old = ax.XLim;
  ylim_old = ax.YLim;
  % add the number of calls into the plot
  text(ax, mean(xlim_old), ylim_old(2)*0.95, sprintf('N=%d', length(idx)), 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'top', 'FontSize', 12, 'Color',  'white');
end

% save to pdf
fn_pdf = fullfile(fd_save, sprintf('%s.spect.MeanWrap.pdf', bd));
print(fig, fn_pdf, '-dpdf', '-painters');


% plot with padded silence, so duration is comparable
% pad
sizes = cellfun(@size, mean_all, 'UniformOutput', false);
row_sizes = cellfun(@(x) x(1), sizes); % Number of rows
col_sizes = cellfun(@(x) x(2), sizes); % Number of columns
[max_cols, max_i] = max(col_sizes);
% Initialize the output cell array
mean_pad = cell(size(mean_all));
for i = 1:numel(mean_all)
  current_array = mean_all{i};
  current_cols = size(current_array, 2);
  if current_cols < max_cols
    % Calculate padding for left and right
    total_padding = max_cols - current_cols;
    left_pad = floor(total_padding / 2);
    right_pad = total_padding - left_pad; % Ensure total padding matches
    % Pad with zeros on both sides
    mean_pad{i} = [zeros(row_sizes(1), left_pad), current_array, zeros(row_sizes(1), right_pad)];
  else
    % No padding needed if already max_cols
    mean_pad{i} = current_array;
  end
end
close all;
fig_pos = [100 100 1600 250];
[fig, axes] = generatePanelGrid_v2(1, 10, [0.75], [], [0.15;0.02], [0.05;0.05], 0.02, [0], fig_pos);% loop through each cluster
for cid = 1:max(embed.hdbscan_cluster)
  oldV = sprintf('v%d', cid);
  newV = swapID.(oldV);
  plot_i = str2num(strrep(newV, 'v', ''));
  ax = axes(plot_i);
  spec_mean = mean_pad{cid};
  idx = find(embed.hdbscan_cluster==cid);
  %   imagesc(ax, 1:size(spec_mean,1), freq, spec_mean); colormap jet;
  %   imagesc(ax, 1:size(spec_mean,1), freq, spec_mean, [0.15 max(spec_mean(:))]); colormap(custom_map);
  rel_t = t_all{max_i};
  rel_t = rel_t - rel_t(1);
  imagesc(ax, rel_t, freq, spec_mean, [0.15 max(spec_mean(:))]); colormap(custom_map);
  set(ax, 'YDir', 'Normal');
  %   colormap jet;
  set(ax, 'YDir', 'normal');
  set(ax, 'YTick', [1000 3000 5000 7000]);
  set(ax, 'YTickLabel', {'1k', '3k', '5k', '7k'}, 'FontSize', 10);
  %   set(ax, 'XTick', []);
  set(ax, 'XTick', [0 max(rel_t)]);
  
  %   xticklabels = arrayfun(@(x) sprintf('%.2f', x), [0 max(rel_t)], 'UniformOutput', false);
  title(ax, newV, 'FontSize', 14, 'Color', colDict.(newV));
  xlim_old = ax.XLim;
  ylim_old = ax.YLim;
  % add the number of calls into the plot
  text(ax, mean(xlim_old), ylim_old(2)*0.95, sprintf('N=%d', length(idx)), 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'top', 'FontSize', 12, 'Color',  'white');
end

% save to pdf
fn_pdf = fullfile(fd_save, sprintf('%s.spect.MeanWrapPad.pdf', bd));
print(fig, fn_pdf, '-dpdf', '-painters');





