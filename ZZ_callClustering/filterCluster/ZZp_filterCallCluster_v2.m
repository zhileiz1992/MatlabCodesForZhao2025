% derive metrics to quantify the quality of the HDBSCAN clustering output
% filter out low-quality clusters
% rename the call subtype clusters for later use 
% differ from v1: add the option to normalize the amplitude within call subtypes when calculate spectrograms


clear; close all;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
fd_save_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
clims = {[10 21]; [10 21]; [10.5 21.5]; [12 23]};
norms = {0; 0; 0; 1};  % normalize if there are calls that are extra loud or quiet
% what spectrogram input dataset
input_rn = 'spec_goffinet_nn_256_176';
vae_run = input_rn; 
% what UMAP/HDBSCAN run to use
suffix_umap = 'UMAPonVAE7'; 
% color for different call clusters
col_list = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#ccad25'};
col_rand = {'#252525','#8c510a','#01665e','#252525','#8c510a','#01665e','#252525'};  % use grey for random data


bi = 3;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
clim = clims{bi};
norm = norms{bi};
% where to save results
fd_save = fullfile(fd_save_base, birdID, suffix_umap, 'v', vae_run, 'plotsFiltered');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
disp(fd_save);


%% 1. Load HBSCAN results
fd_cluster = fullfile(fd_save_base, birdID, suffix_umap, 'v', input_rn);
fn_embed = fullfile(fd_cluster, sprintf('%s.%s.embedding.csv', birdID, input_rn));
embed = readtable(fn_embed);


%% 2. Calculate the sparseness index of spectrogram features
% Create avereaged spectrogram: with time wrapping
% read sound directly from source nc file
fs = 20000;
% add a little before and after
pad = 200;  % 10 ms
segs_all = {};
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
  segs_all{cid} = segs;
end


%% time-warp the calculate spectrograms then calculate sparseness
mean_all = {};
t_all = {};
sparse_all = {};
for cid = 1:max(embed.hdbscan_cluster)
  segs = segs_all{cid};
  % normalize the amplitude
  if norm==1
    q99 = cellfun(@(x) quantile(x, 0.99), segs);
    segs_raw = segs;
    parfor ii=1:size(segs,2)
      segs{ii} = segs_raw{ii} / q99(ii)*0.25;  % max amplitude is 0.25
    end
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
    [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(seg_aligned, fs, 256, 256, 236, [250 7500], clim);
    spec_aligned{ii} = powerGrey/1024;  % normalize
%     spec_aligned{ii} = powerRaw;
  end
  stacked_array = cat(3, spec_aligned{:});
  spec_mean = mean(stacked_array, 3);
  mean_all{cid} = spec_mean;
  % save the time information
  n = length(segs{1});
  seg_aligned = interp1(linspace(1, n, n), segs{1}, linspace(1, n, m), 'linear');
  [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(seg_aligned, fs, 256, 256, 236, [250 7500], clim);
  t_all{cid} = t;

  % calculate sparseness
  % focus on certain freq range
  f_range = [0 10000]; 
  f_idx = find((f>=f_range(1)) & (f<=f_range(2)));
  sparse_this = zeros(size(stacked_array, 2),1);
  parfor ti=1:size(stacked_array, 2)
    s_this = squeeze(stacked_array(f_idx,ti,:));
    p = mean(s_this, 2)';
    p = p / sum(p);
    % calculate entropy measure of sparseness
    s = p .* log(p);
    sparseness = 1 + nansum(s) / log(size(p,2));
    sparse_this(ti) = sparseness;
  end
  sparse_all{cid} = sparse_this;
  
end
mean_sparse = cellfun(@(x) median(x), sparse_all);



%%  plot mean spectrogram to check
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
close all;
fig_pos = [100 100 1600 250];
[fig, axes] = generatePanelGrid_v2(1, 10, [0.75], [], [0.15;0.02], [0.05;0.05], 0.02, [0], fig_pos);% loop through each cluster
for cid = 1:max(embed.hdbscan_cluster)
  ax = axes(cid);
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
  title(ax, sprintf('V%d %.3f',cid, mean_sparse(cid)), 'FontSize', 10);
  xlim_old = ax.XLim;
  ylim_old = ax.YLim;
  % add the number of calls into the plot
  text(ax, mean(xlim_old), ylim_old(2)*0.95, sprintf('N=%d', length(idx)), 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'top', 'FontSize', 12, 'Color',  'white');
end
% save to pdf
fn_pdf = fullfile(fd_save, sprintf('%s.spect.MeanWrap.pdf', birdID));
print(fig, fn_pdf, '-dpdf', '-painters');



%% plot the spareness of different call subtypes
thre = 0.12;  % a threshold to fail 
[fig, axes] = generatePanelGrid_v2(1, 2, [0.75], [], [0.1;0.05], [0.1;0.05], 0.1, [0], [10 10 1000 500]);% loop through each cluster
% plot how sparseness change with time
ax = axes(1); hold(ax, 'on');
% col_lines = lines(size(sparse_all,2));
for cid=1:size(sparse_all,2)
  plot(ax, sparse_all{cid}, 'DisplayName', sprintf('oldV%d',cid), 'Color', col_list{cid}, 'LineWidth',3); 
end
xlabel(ax, 'Rel. time');
ylabel(ax, 'Frequency sparseness');
legend(ax, 'Location', 'best');
% the plot the average sparseness
ax2 = axes(2);
plot(ax2, 1:length(mean_sparse), mean_sparse, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'LineStyle', 'none', 'MarkerSize', 15);
yline(ax2, thre, 'LineStyle', '--', 'Color', 'red');
xticks(ax2, 1:length(mean_sparse));
xlim(ax2, [0.5 length(mean_sparse)+0.5]);
xticklabels(ax2, 1:length(mean_sparse));
xlabel(ax2, 'Call subtype ID');
ylabel(ax2, 'Median spareness');
fn_fig = fullfile(fd_save, sprintf('%s.filterSubtypeSparseness.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');



%%  set a threshold to remove bad clusters
fail_idx = find(mean_sparse<thre);
disp('Failed to pass threshold:')
disp(fail_idx);

% save the data for later use
sparse.sparse_all = sparse_all; 
sparse.mean_sparse = mean_sparse;
sparse.fail_idx = fail_idx;
fn_mat = fullfile(fd_save, sprintf('%s.filterSubtypeSparseness.mat', birdID));
save(fn_mat, 'sparse');









