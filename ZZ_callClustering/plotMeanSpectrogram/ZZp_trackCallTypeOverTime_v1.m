% check the consistency of call subtypes over time
% Zhilei, 06/19/2025
% 1. check the proportion of each subtype over time
% 2. check the averaged spectrograms over time

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));

%% Folder setting
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% where call embedding results are stored
fd_embed_base = fullfile(fd_home, 'vaeWav');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what clustering run to use
syl = 'v';
input_rn = 'spec_goffinet_cutoff_256_176';
fd_vae = 'UMAPonVAE6'; vae_run = input_rn;
% what call subtype dbase to use
suffix = 'Wsp1Call';
suffix_old = 'Wsp.disVAEcall';
% where to save plots
fd_save_base = fullfile(fd_embed_base, 'combined', 'CallClusteringMetric');

% loop through birds
bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Read the call subtype data from dbase
fn_info = fullfile(fd_home, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
load(fn_info);
% get unique dates
date_unique = unique(info.date_long, 'sorted');
% get all labels from each date
labels_all = {};
parfor di=1:size(date_unique,1)
  data_date = date_unique{di};
  date_short = strrep(data_date,'-','');
  fd_master = fullfile(fd_home, 'DbaseFiles', pairID, data_date, birdID, 'warble');
  % what segmentation dbase to use
  fn_seg = fullfile(fd_master, sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, suffix));
  dbase = load(fn_seg).dbase; 
  % only focus on the production files
  p_idx = find(strcmp(dbase.PropertyNames, 'bProd'));
  pass_idx = find(dbase.Properties(:,p_idx));
  labels = dbase.SegmentTitles;
  labels_all{di} = [labels{pass_idx}];
end


%% 2. Plot the proportion of each type over time
% GPT
% Extract all unique labels
all_labels_flat = [labels_all{:}];
unique_labels = unique(all_labels_flat);
% Count occurrences per label per date
num_dates = numel(labels_all);
num_labels = numel(unique_labels);
label_counts = zeros(num_labels, num_dates);
for d = 1:num_dates
    labels_today = labels_all{d};
    if isempty(labels_today), continue; end
    for l = 1:num_labels
        label_counts(l, d) = sum(strcmp(unique_labels{l}, labels_today));
    end
end
% Convert to proportions
total_counts_per_date = sum(label_counts, 1);
label_proportions = label_counts ./ total_counts_per_date;
% Plot each label in its own subplot
% convert dates to relative days
date_dt = datetime(date_unique, 'InputFormat', 'yyyy-MM-dd');
% Calculate relative days to the first date
relative_days = days(date_dt - date_dt(1))+1;
close all;
fig_pos = [10 10 1400 800];
[fig, axes] = generatePanelGrid_v2(3, 6, [0.22;0.22;0.22], [0.08;0.08], [0.1;0.02], [0.05;0.05], 0.05, [0;0;0], fig_pos);
for l = 1:num_labels
    plot_i = floor((l-1)/6)+1;
    plot_j = mod(l-1, 6)+1; 
    ax = axes(plot_i, plot_j);
    plot(ax, 1:num_dates, label_proportions(l,:), 'LineWidth', 2);
    ylabel(ax, unique_labels{l});
    ylim(ax, [0 max(label_proportions(l,:))*1.1]);
    title(ax, unique_labels{l});
end
sgtitle('Proportion of Each Label Over Time');
% save as pdf
fn_pdf = fullfile(fd_save_base, sprintf('%s.propChange.pdf', birdID));
print(fig, fn_pdf, '-dpdf', '-painters');
  

%% 3. Check the averaged spectrograms overtime
% divide experimental period into equal sized portions
num_portion = 4; 
num_per = floor(size(date_unique,1) / num_portion);
dstart = 1:num_per:(num_portion*num_per);
dend = dstart+num_per-1;
dend(end) = size(date_unique,1);

% grab call subtypes only
is_match = ~cellfun('isempty', regexp(unique_labels, '^v(0*[1-9]\d*)$'));
% Extract matching elements
labels_extract = unique_labels(is_match);

% read embedding table
fn_embed = fullfile(fd_embed_base, birdID, fd_vae, syl, input_rn, sprintf('%s.%s.embedding.csv', birdID, vae_run));
embed = readtable(fn_embed);
% add a date field to the table
second_last = cellfun(@(s) strsplit(s, '/'), embed.fn_wav, 'UniformOutput', false);
embed.date = cellfun(@(parts) parts{end-1}, second_last, 'UniformOutput', false);

% loop through each call subtype, calculate mean spectrogram in each date portion
fs = 20000;
% add a little before and after
pad = 200;  % 10 ms
mean_all = {};
for cid = 1:max(embed.hdbscan_cluster)
  % loop through each date portion
  for dp_i=1:length(dstart)
    date_set = date_unique(dstart(dp_i):dend(dp_i));
    embed_s = embed((embed.hdbscan_cluster==cid) & ismember(embed.date, date_set), :);
  
    % read sound data in
    fns = embed_s.fn_wav;
    istart = embed_s.istart;
    iend = embed_s.iend;
    segs = {};
    parfor si=1:length(fns)
      fn = fns{si};
      [sound, fs] = audioread(fn);
      sound = ZZ_butterFilterFunc_v1(sound, fs, 250, 7500, 5);
      i_start = max([1 istart(si)-pad]);
      i_end = min([size(sound,1) iend(si)+pad]);
      seg = sound(i_start:i_end);
      segs{si} = seg;
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
    mean_all{cid, dp_i} = spec_mean;
  end
end

% plot the averaged spectrogram
s = embed{1, 'spec_f'}{1};
s = erase(s, {'[ ', ']'});
s = split(s);
freq = [];
for i=1:(size(s,1)-1)
  freq = [freq str2num(s{i})];
end
% a custom colormap for spectrogram
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet
% loop through call subtypes
for cid = 1:max(embed.hdbscan_cluster)
  close all;
  fig_pos = [100 100 600 250];
  [fig, axes] = generatePanelGrid_v2(1, num_portion, [0.75], [], [0.15;0.02], [0.05;0.05], 0.05, [0], fig_pos);% loop through each cluster
  for dp_i=1:size(mean_all,2)
    ax = axes(dp_i);
    spec_mean = mean_all{cid, dp_i};
    imagesc(ax, 1:size(spec_mean,1), freq, spec_mean, [0.15 max(spec_mean(:))]); colormap(custom_map);
    set(ax, 'YDir', 'Normal');
    %   colormap jet;
    set(ax, 'YDir', 'normal');
    set(ax, 'YTick', [1000 3000 5000 7000]);
    set(ax, 'YTickLabel', {'1k', '3k', '5k', '7k'}, 'FontSize', 10);
    set(ax, 'XTick', []);
    title(ax, sprintf('V%d part%d', cid, dp_i), 'FontSize', 12);
  end
  % save as pdf
  fn_pdf = fullfile(fd_save_base, sprintf('SpectroChange.%s.v%d.pdf', birdID, cid));
  print(fig, fn_pdf, '-dpdf', '-painters');
end

















