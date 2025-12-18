% a script to examine WhisperSeg results across birds
% Zhilei, 05/27/2025
% 1. Count proportions of syllable categories across birds, show as stacked bar plots
% 2. Export examples of each syllable category to use in supplementary figure
% differ from v1: excluding the auditory files

close all; clear;
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_examineSortedDbase'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_WhisperSegMO'));


%% 0. Inputs
fd_z4 = '/mnt/z4';
fd_home = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
% what birds to extract
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ID = {'M1', 'M2', 'M3', 'M4'};
% what folder has the dbase files
fd_master = fullfile(fd_home, 'DbaseFiles');
% where to save intermediate results
fd_save_itm= fullfile(fd_home, 'Figures');
% what type of WhisperSeg results to use, specify by suffix
suffix = 'Wsp1';


%% 1. Count the number of syllables
syl_type = {'v'; 'h'; 'e'; 'b'; 'x'};
pretty_syl = {'Call', 'Click', 'Squawk', 'Chortle', 'Other'};
syl_count = table(syl_type, 'VariableNames', {'syl_type'});
% loop through birds
% bi = 1;
for bi=1:size(birdIDs,2)
  birdID = birdIDs{bi};
  pairID = pairIDs{bi};
  % load meta info
  fn_meta = fullfile(fd_master, pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
  load(fn_meta);
  % Loop through each date, segment/annotate the warble with WhispserSeg
  % what date have sorted neurons
  date_unique = unique(info.date);
  label_all = {};
  for date_i=1:size(date_unique, 1)
    % load dbase with WhisperSeg results
    date_short = date_unique{date_i};
    data_date = strjoin({date_short(1:4), date_short(5:6), date_short(7:8)}, '-'); % add a dash in between
    fn_base = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, suffix);
    fn_d = fullfile(fd_master, pairID, data_date, birdID, 'warble', fn_base);
    load(fn_d);
    % which field is the bProd, should be the 3rd, but double check 
    p_idx = find(strcmp(dbase.PropertyNames, 'bProd'));
    p = dbase.Properties;
    inc_idx = find(p(:,p_idx)); 
    % get all WhisperSeg labels
    label_this = horzcat(dbase.SegmentTitles{inc_idx});
    label_all = [label_all label_this];
  end
  % sum the counts of each category
  for syl_i=1:size(syl_type,1)
    syl_count.(birdID)(syl_i) = sum(strcmp(label_all, syl_type{syl_i}));
  end
end


%% 2. Plot propertion as stacked barplot
% use default color for each categories
col_list = {'#e41a1c', '#984ea3', '#4daf4a', '#377eb8', '#737373'};
syl_types = syl_count.syl_type;         % 5x1 cell array or categorical
count_data = table2array(syl_count(:, 2:5));  % 5x4 numeric matrix
% Compute proportions per bird
proportion_data = count_data ./ sum(count_data, 1);
% Create stacked bar plot
close all;
fig = ZZfunc_newFigurePDFsize_v1([50 50 350 350]);
h = bar(proportion_data', 'stacked');  % Transpose so birds are on x-axis
% Set color for each syllable type (bar series)
for k = 1:numel(h)
  h(k).FaceColor = col_list{k};
end
% Customize plot
xticks(1:4);
xticklabels(pretty_ID);
xlabel('Bird ID');
ylabel('Proportion of syllable counts');
legend(pretty_syl, 'Location', 'eastoutside');
title('Syllable proportion in counts');
% save as pdf
fn_pdf = fullfile(fd_save_itm, 'CombinedAnalysis', sprintf('SyllablePropChan0.%s.pdf', suffix));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 3. Export example spectrograms
% sample equal number of examples per bird
% save information about syllable, e.g. filename and position, then do sampling, then read data
syl_info = struct();
count = 0;
for bi=1:size(birdIDs,2)
  birdID = birdIDs{bi};
  pairID = pairIDs{bi};
  % load meta info
  fn_meta = fullfile(fd_master, pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
  load(fn_meta);
  date_unique = unique(info.date);
  for date_i=1:size(date_unique, 1)
    % load dbase with WhisperSeg results
    date_short = date_unique{date_i};
    data_date = strjoin({date_short(1:4), date_short(5:6), date_short(7:8)}, '-'); % add a dash in between
    fn_base = sprintf('%s.%s.warble.good.%s.dbase.mat', birdID, date_short, suffix);
    fn_d = fullfile(fd_master, pairID, data_date, birdID, 'warble', fn_base);
    load(fn_d);
    % get all WhisperSeg labels and related info
    sf = dbase.SoundFiles;
    % only include production
    % which field is the bProd, should be the 3rd, but double check 
    p_idx = find(strcmp(dbase.PropertyNames, 'bProd'));
    p = dbase.Properties;
    inc_idx = find(p(:,p_idx)); 
    for sii=1:length(inc_idx)
      si = inc_idx(sii);
      syl_n = dbase.SegmentTitles{si};
      if ~isempty(syl_n)
        % what's the audio file name?
        fn_aud = fullfile(sf(si).folder, sf(si).name);
        % get start and end time as well
        st = dbase.SegmentTimes{si};
        for sj=1:size(syl_n,2)
          count = count+1;
          syl_info(count).birdID = birdID;
          syl_info(count).syl_type = syl_n{sj};
          syl_info(count).start = st(sj,1);
          syl_info(count).end = st(sj,2);
          syl_info(count).fn_aud = fn_aud;
        end
      end
    end
  end
end

% sample certain number of syllables, read in audio data, plot spectrograms
% how much to pad before and after the segmentation
pad = 0.05; 
fs = 20000; 
pad_pt = floor(pad*fs);
% how many renditions to sample
% to_sample = 6;
to_sample_list = [6 12 8 6 12];
% seed = 1992; 
seed = 1118;
rng(seed);
% loop through syllable type
for syl_i=1:size(syl_type,1)
  to_sample = to_sample_list(syl_i);
  % loop through each bird, randomly sample, retrive data, determine figure panel size
  data_lengths = zeros(size(birdIDs,2), to_sample);
  data_heights = zeros(size(data_lengths));
  spec = cell(size(data_lengths));
  rel_t = cell(size(data_lengths));
  for bi=1:size(birdIDs,2)
    % avoid boundary syllables
    idx = find((strcmp({syl_info.syl_type}, syl_type{syl_i})) & (strcmp({syl_info.birdID}, birdIDs{bi})) ...
      & ([syl_info.start]>pad_pt) & ([syl_info.end]<(15*fs-pad_pt)));
    idx_rd = randsample(idx, to_sample);
    for jj=1:length(idx_rd)
%       data_len(syl_i, to_sample*(bi-1)+jj) = syl_info(idx_rd(jj)).end - syl_info(idx_rd(jj)).start + 2*pad_pt;
      % read the data
      d = ncread(syl_info(idx_rd(jj)).fn_aud, 'data');
      % determine start and end
      istart = max([1 syl_info(idx_rd(jj)).start-pad_pt]);
      iend = min([size(d,1) syl_info(idx_rd(jj)).end+pad_pt]);
      d_seg = d(istart:iend);
      % generate spectrogram
      [power, powerRaw, powerGrey, S, f, t] = getAudioSpectrogramZZ_flexible_v1(d_seg, fs, 256, 256, 236, [250 7500], [12.5 25]);
      % what's the time offset between 1st spectrogram column and syllable start
      t_offset = (syl_info(idx_rd(jj)).start - istart) / fs;
      t_spec = t - t_offset;
      data_lengths(bi,jj) = size(power,2);
      data_heights(bi,jj) = size(power,1);
      spec{bi,jj} = power;
      rel_t{bi,jj} = t_spec;
    end
  end
  
  % create figure layout
  xgap_in = 0.05;
  ygap_in = 0.2;
  points_per_inch_x = 800;
  points_per_inch_y = 150;
  left_margin = 0.5;
  bottom_margin = 0.5;
  close all; 
  [fig, axes] = generateFigureByDataLength(data_lengths, data_heights, xgap_in, ygap_in, points_per_inch_x, points_per_inch_y, left_margin, bottom_margin);
  
  % fill in spectrograms
  for bi=1:size(birdIDs,2)
    for jj=1:to_sample
      ax = axes(bi, jj);
      imagesc(ax, rel_t{bi,jj}, f, spec{bi,jj});
      colormap jet;
    %   set(ax2, 'XTickLabel', []);
      set(ax, 'YDir', 'normal');
      set(ax, 'YTick', [1000 3000 5000 7000]);
      set(ax, 'YTickLabel', {'1k', '3k', '5k', '7k'});
      % use the duration of the syllable as xtick labels
      ax.XTick = [rel_t{bi,jj}(end)/2];
      ax.XTickLabel = sprintf('%.2f', rel_t{bi,jj}(end)-pad);
      ax.XTickLabelRotation = 0;
      ax.FontSize = 5;
      if jj>1
        ax.YTick = [];
        ax.YTickLabel = [];
      end
    end
  end
  
  % save figure
  fn_pdf = fullfile(fd_save_itm, 'CombinedAnalysis', sprintf('ExampleSpec.seed%d.%s.pdf', seed, syl_type{syl_i}));
  print(fig, fn_pdf, '-dpdf', '-painters');
end


%% 4. Plot duration distribution
s = syl_info;
unique_types = flipud(syl_types);
col_list2 = fliplr(col_list);
n_types = numel(unique_types);
% Create tiled layout
close all; fig = figure;
t = tiledlayout(n_types, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
xlabel(t, 'Duration (s)');
ylabel(t, 'Syllable count');
all_durations = ([s.end] - [s.start]) / 20000;
xlim_shared = [min(all_durations), max(all_durations)];  % show the full range
disp(xlim_shared);
xlim_shared = [0 1.2];  % or a define range
bin_width = 0.02;  % in seconds
edges = xlim_shared(1):bin_width:xlim_shared(2);
for i = 1:n_types
    syl_type = unique_types{i};
    ax = nexttile;
    % Filter all entries of this syllable type
    mask = strcmp({s.syl_type}, syl_type);
    durations = ([s(mask).end] - [s(mask).start]) / 20000;
    % Plot histogram
    histogram(ax, durations, 'BinEdges', edges, 'FaceColor', col_list2{i}, 'EdgeColor', 'none');
%     title(ax, pretty_syl{i}, 'FontSize', 14, 'Interpreter', 'none');
    xlim(ax, xlim_shared);
    ax.YAxis.Exponent = 3;
    % Share x-axis
    if i < n_types
        ax.XTickLabel = [];
    end
end
set(fig, 'Position', [50 50 400 600]);
% save as PDF
fig.Units = 'inches';
pos = fig.Position;
% Set the paper size to match the figure size
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 pos(3) pos(4)];  % [left bottom width height]
fig.PaperSize = [pos(3) pos(4)];
% save the figure
fn_pdf = fullfile(fd_save_itm, 'CombinedAnalysis', 'Duration_hist.pdf');
print(fig, fn_pdf, '-dpdf', '-painters');


%% 5. Calculate syllable proportion by duration
% Example input:
% data: struct array with fields birdID, syl_type, start, end
% col_list: cell array of color values (e.g., {'#1f77b4','#ff7f0e','#2ca02c'})
bird_list = birdIDs;
syl_list = syl_type;
data = syl_info;
% Initialize duration matrix: rows = birds, columns = syl types
count_data = zeros(length(syl_list), length(bird_list));
% Fill in the duration matrix
for i = 1:length(data)
    bird_idx = find(strcmp(bird_list, data(i).birdID));
    syl_idx = find(strcmp(syl_list, data(i).syl_type));
    dur = (data(i).end - data(i).start) / 20000;  % duration in seconds
    count_data(syl_idx, bird_idx) = count_data(syl_idx, bird_idx) + dur;
end
% Compute proportions per bird
proportion_data = count_data ./ sum(count_data, 1);
% Create stacked bar plot
close all;
fig = ZZfunc_newFigurePDFsize_v1([50 50 350 350]);
h = bar(proportion_data', 'stacked');  % Transpose so birds are on x-axis
% Set color for each syllable type (bar series)
for k = 1:numel(h)
  h(k).FaceColor = col_list{k};
end
% Customize plot
ylim([0 1]);
xticks(1:4);
xticklabels(pretty_ID);
xlabel('Bird ID');
ylabel('Proportion of syllable duration');
legend(pretty_syl, 'Location', 'eastoutside');
title('Syllable proportion in duration');
% save as pdf
fn_pdf = fullfile(fd_save_itm, 'CombinedAnalysis', sprintf('SyllablePropDurationChan0.%s.pdf', suffix));
print(fig, fn_pdf, '-dpdf', '-painters');












  
  
  
  
  
  
  
  
  
  
