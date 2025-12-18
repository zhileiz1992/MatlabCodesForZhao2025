% calculate and plot the firing field of MO neurons in the acoustic space
% step 6: plot the firing fields from assigned sliding windows for spikes
% differ from v2: vary the lag between spike and sliding window start to allow for more premotor lag

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/overlaySpikes'));

%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much extra premotor lag to add
lag_frame = 5;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
% what's the range to include spikes: default only consider spikes within 1-window of syllable boundaries
epad_frame = win_frame;



bi = 2;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
% where is VAE data located
fd_vae = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
% where is spike-window data located
fd_sliding_loc = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_slidingLoc_lag' num2str(lag_frame)]);
fns_loc = dir(fullfile(fd_sliding_loc, sprintf('%s.*.sliding_loc.replaced2.mat', birdID)));
% where is UMAP results located
fd_umap = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_umap']);
% umap_run = 'all_syl';
umap_run = 'all_syl2';
% load UMAP results
fn_info = fullfile(fd_umap, sprintf('%s.%s.info2.mat', birdID, umap_run));
load(fn_info);
tempInfo = info(:, {'syl_ID', 'lens', 'ustart', 'uend'});
fn_umap = fullfile(fd_umap, sprintf('%s.%s.umap_res.csv', birdID, umap_run));
umap_res = readmatrix(fn_umap);

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% where to save results
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField2_lag' num2str(lag_frame)]);
if ~exist(fd_field, 'dir'); mkdir(fd_field); end



%% 1. Load information regarding spikes and associated sliding windows
spike_win = [];
for fi=1:size(fns_loc, 1)
  ss = strsplit(fns_loc(fi).name, '.');
  ss = ss{2};
  fprintf('Loading spike_loc data for %s...\n', ss);
  cate = ss(1);
  fn_loc = fullfile(fns_loc(fi).folder, fns_loc(fi).name);
  load(fn_loc);
  spike_embed = struct2table(spike_embed);
  % add syllable label and category to the table
  spike_embed.label = repelem({ss}, size(spike_embed,1))';
  spike_embed.category = repelem({lower(cate)}, size(spike_embed,1))';
  spike_win = [spike_win; spike_embed];
end
% save for future use
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
save(fn_spike, 'spike_win', '-v7.3');


%% 2. Create background images to show the overall aoucsitc space
% overall space in gray, then category region in brown
% randomly select n number of syllables
num_plot = 300;
rng(1992);
syls_plot = {'v', 'b', 'h', 'e', 'x', 'comp'};
alpha_list = [0.03 0.02 0.03 0.03 0.03 0.02];
i_all = cell(6, 1);
for si=1:5
  idx_this = find(strcmp(info.category, syls_plot{si}));
  i_all{si} = randsample(idx_this, num_plot);
end
i_comb = cat(1, i_all{:});
% then sample index from compound syllables
idx_comp = find((ismember(info.category, {'b','x'})) & (info.dur>=0.3));
i_comp_rd = randsample(idx_comp, num_plot);
i_all{6} = i_comp_rd;

% all figures use the same x/y axis limits
% x_lim = [0.5 19]; y_lim = [-6.5 18.5];  % bird M1
% x_lim = [-8 10]; y_lim = [-10 16];  % bird M2, first half
x_lim = [0 17]; y_lim = [-12 16];  %for bird M2, first half

% loop through syllables, plot embedding
img_all = struct();
for si=1:6
  close all;
  fig_size = [10 10 600 600];
  fig = ZZfunc_newFigurePDFsize_v1(fig_size);
  ax = gca; cla(ax);
  % first all other syllable types
  i_other = setdiff(i_comb, i_all{si});
  ax = ZZfunc_plotEmbedGray_v1(ax, info, umap_res, i_other, x_lim, y_lim, 10, '#d9d9d9', 0.01);
  % then the focal syllable type
  ax = ZZfunc_plotEmbedGray_v1(ax, info, umap_res, i_all{si}, x_lim, y_lim, 10, '#6baed6', alpha_list(si));
  % rasterize to image
  frame = getframe(ax); img = frame.cdata;
  img_all.(syls_plot{si}) = img;
end

% save the image data for re-use
fn_img = fullfile(fd_field, sprintf('%s.background.img_all.mat', birdID));
save(fn_img, 'img_all');



%% 3.2 Plot firing fields and density: restrict b/x to comp (>0.3 sec)
% also calculate and plot the spike density
num_bin = 100;

fd_save_plot2 = fullfile(fd_field, 'plot_comp');
fd_save_overlay = fullfile(fd_save_plot2, 'spike_overlay');
fd_save_density = fullfile(fd_save_plot2, 'spike_density');
if ~exist(fd_save_overlay, 'dir'); mkdir(fd_save_overlay); end
if ~exist(fd_save_density, 'dir'); mkdir(fd_save_density); end
% distinguish different syllable catorgories
cat_plot = {'v', 'b', 'x', 'h', 'e'};
cat_color = {'#e41a1c', '#ff7f00', '#f781bf', '#984ea3', '#4daf4a'};
cat_loc = {[1,1], [1,3], [2,1], [2,2], [2,3]};  % panel location in the figure

% ni = 44;
% for ni=1:size(info_neu,1)
% for ni=1:81
for ni=82:121
  neuronID = info_neu.neuronID{ni};
  %   neuronID = '20240917-ch13';
  % subset the spike_win table
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  spike_this = outerjoin(spike_this, tempInfo, 'Keys', 'syl_ID', 'Type', 'left', 'MergeKeys', true);
  close all;
  % one plot for spike overlay; another for spike density
  [fig, axes] =  generatePanelGrid_v3(2, 3, [0.38;0.38], [0.08], [0.1;0.05], [0.05;0.05], 0.05, [0;0], [10 50 1200 900], false);
  [fig2, axes2] =  generatePanelGrid_v3(2, 3, [0.38;0.38], [0.08], [0.1;0.05], [0.05;0.05], 0.05, [0;0], [10 50 1200 900], false);
  for si=1:size(cat_plot,2)
    cate = cat_plot{si};
    spike = spike_this(strcmp(spike_this.category, cate),:);
    % first plot the spike overlay
    axloc = cat_loc{si};
    ax = axes(axloc(1), axloc(2)); cla(ax); hold(ax, 'on');
    xlim(ax, x_lim); ylim(ax, y_lim);
    % add the scatter background
    img = img_all.(cate);
    image(ax, ax.XLim, ax.YLim, flipud(img));
    % overlay spikes
    [ax, spike_loc] = ZZfunc_plotSpikeOnUmapComp_v2(ax, spike, umap_res, 30, cat_color{si}, 0.2);
    xlim(ax, x_lim); ylim(ax, y_lim);
    title(ax, sprintf('%s (%s)', cate, neuronID), 'FontSize', 14);
    
    % then calculate and plot density
    if ~isempty(spike_loc)
      [density, x_edges, y_edges, x_centers, y_centers] = ZZfunc_calcSpikeDensity_count_V1(spike_loc(:,1), spike_loc(:,2), num_bin, x_lim, y_lim, 3);
      density = density';
      ax2 = axes2(axloc(1), axloc(2)); cla(ax2); hold(ax2, 'on');
      xlim(ax2, x_lim); ylim(ax2, y_lim);
      % add the scatter background
      h1=image(ax2, ax2.XLim, ax2.YLim, flipud(img));
      set(h1, 'AlphaData', 0.75);
      d_qt = quantile(density(:), 0.999);
      d_qt = max([1 d_qt]);
      h2=imagesc(ax2, x_centers, y_centers, density, [0, d_qt]);
      colormap(ax2, flipud(hot));
      set(h2, 'AlphaData', 0.25);
      xlim(ax2, x_lim); ylim(ax2, y_lim);
      title(ax2, sprintf('%s (%s) %.2f', cate, neuronID, d_qt), 'FontSize', 12);
    end
    
  end
  % then plot for the compound syllable
  ax = axes(1,2);  cla(ax); hold(ax, 'on');
  spike = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
  spike = spike(spike.dur>=0.3, :);
  % overlay spikes
  xlim(ax, x_lim); ylim(ax, y_lim);
  % add the scatter background
  img = img_all.comp;
  image(ax, ax.XLim, ax.YLim, flipud(img));
  % overlay spikes
  [ax, spike_loc] = ZZfunc_plotSpikeOnUmapComp_v2(ax, spike, umap_res, 30, '#a65628', 0.2);
  xlim(ax, x_lim); ylim(ax, y_lim);
  title(ax, sprintf('Compound syl. (%s)', neuronID), 'FontSize', 14);
  
  [density, x_edges, y_edges, x_centers, y_centers] = ZZfunc_calcSpikeDensity_count_V1(spike_loc(:,1), spike_loc(:,2), num_bin, x_lim, y_lim, 3);
  density = density';
  ax2 = axes2(1,2); cla(ax2); hold(ax2, 'on');
  xlim(ax2, x_lim); ylim(ax2, y_lim);
  % add the scatter background
  h1=image(ax2, ax2.XLim, ax2.YLim, flipud(img));
  set(h1, 'AlphaData', 0.75);
  d_qt = quantile(density(:), 0.999);
  d_qt = max([1 d_qt]);
  h2=imagesc(ax2, x_centers, y_centers, density, [0, d_qt]);
  colormap(ax2, flipud(hot));
  set(h2, 'AlphaData', 0.25);
  xlim(ax2, x_lim); ylim(ax2, y_lim);
  title(ax2, sprintf('Compound syl. (%s) %.2f', neuronID, d_qt), 'FontSize', 12);
  
  % save results
  fn_fig = fullfile(fd_save_overlay, sprintf('%s.%s.firingFields.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
  
  fn_fig2 = fullfile(fd_save_density, sprintf('%s.%s.firingDensity.pdf', birdID, neuronID));
  print(fig2, fn_fig2, '-dpdf', '-painters');
  
end



%% 4. Plot example spectrograms and spikes
% make sure the same figure width translates into the same temporal duration
% spectrogram to top, spike ticks in middle, UMAP trajectory at bottom
fd_save_spec = fullfile(fd_save_plot2, 'spike_spectrogram');
num_sample = 50;
% for sound how much to pad before syllable onset
fs = 20000;
pad_sound = 0.032; pad_sound_pt = floor(pad_sound*fs);

rng(1992);
% ni = 44;
for ni=26:size(info_neu,1)
  neuronID = info_neu.neuronID{ni};
  fd_save_this = fullfile(fd_save_spec, neuronID);
  if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end
  % subset the spike_win table
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  % merge the umap info table
  spike_this = outerjoin(spike_this, tempInfo, 'Keys', 'syl_ID', 'Type', 'left', 'MergeKeys', true);
  % divide into calls and compound syllables
  calls = spike_this(strcmp(spike_this.category, 'v'),:);
  compounds = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
  compounds = compounds(compounds.dur>=0.3,:);
  
  % embed syllables in a matrix with max visible duration, unit is second
  max_dur = 0.8;
  % get the dims of such a long sound
  soundLong = zeros(floor(max_dur*fs),1);
  [powerLong, ~, ~, ~, ~, tLong] = getAudioSpectrogramZZ_flexible_v1(soundLong, fs,  256, 256, 236, [1000 5000], [12 23]);
  %   powerBlank =
  
  % plot calls, then compound syllables
  to_plot_d = {calls, compounds};
  plot_suffix = {'calls', 'compounds'};
  img_list = {img_all.v, img_all.comp};
  
  for di=1:size(to_plot_d, 2)
    d_plot = to_plot_d{di};
    % read data
    % read dbase
    act_num = min([500 size(d_plot,1)]);
    i_rd = 1:act_num;
    fns_dbase = unique(d_plot.fn_dbase(i_rd), 'sorted');
    dbase_list = cell(size(fns_dbase, 1), 1);
    %     fd_dbase = fullfile(fd_base, 'DbaseFiles', pairID, info_neu.date_long{ni}, birdID, 'warble');
    for fi=1:size(fns_dbase, 1)
      fn = strsplit(fns_dbase{fi}, '/');
      fn_temp = strsplit(fn{end}, '.');
      date_str = fn_temp{2};
      date_long = [date_str(1:4) '-' date_str(5:6) '-' date_str(7:8)];
      fd_dbase = fullfile(fd_base, 'DbaseFiles', pairID, date_long, birdID, 'warble');
      a = load(fullfile(fd_dbase, fn{end}));
      dbase_list{fi} = a.dbase;
    end
    % read sound and spike data
    sound_spike = [];
    for ii=1:length(i_rd)
      idx = i_rd(ii);
      fn = d_plot.fn_dbase{idx};
      fi = find(strcmp(fns_dbase, fn));
      dbase_spike = dbase_list{fi};
      % locate sound files
      sf = dbase_spike.SoundFiles;
      fn_sound = d_plot.fn_audio{idx};
      temp = strsplit(fn_sound, '/');
      sf_i = find(strcmp({sf.name}, temp{end}));
      istart = d_plot.seg_start_ori(idx);
      iend = d_plot.seg_end_ori(idx);
      ch_num = regexp(info_neu.channel{ni}, '\d+', 'match');
      ch_num = str2double(ch_num{1}); % Convert to numeric
      spike_shape = str2double(info_neu.spike_shape{ni});
      % retrieve sound and spike data
      [sound, e_trace, spike_iv] = ZZfunc_getDataFromDbase_v1(dbase_spike, sf_i, istart, iend, pad_sound, fs, ch_num, spike_shape);
      sound_spike(ii).syl_ID = d_plot.syl_ID{idx};
      sound_spike(ii).sound = sound;
      sound_spike(ii).e_trace = e_trace;
      sound_spike(ii).spike_iv = spike_iv;
    end
    
    % how many figures to generate
    num_fig = floor(size(sound_spike,2) / 3);
    
    for fig_i=1:num_fig
      [fig, axs] =  generatePanelGrid_v3(3, 3, [0.2;0.03;0.6], [0.02;0.02], [0.05;0.05], [0.05;0.01], 0.01, [0;0;0], [10 50 1800 900], false);
      for ii=1:3
        idx = (fig_i-1) * 3 + ii;
        % plot spectrogram
        ax1 = axs(1, ii); cla(ax1);
        [power, ~, ~, S, f, t] = getAudioSpectrogramZZ_flexible_v1(sound_spike(idx).sound, fs,  256, 256, 236, [1000 5000], [12 23]);
        powerThis = powerLong;
        powerThis(:,1:size(power,2),:) = power;
        imagesc(ax1, tLong, f, powerThis);
        set(ax1, 'YDir', 'normal');
        % add a title
        i_in_u = find(strcmp(spike_this.syl_ID, sound_spike(idx).syl_ID));
        label_this = spike_this.label{i_in_u};
        title(ax1, sprintf('%d %s %s', idx, label_this, sound_spike(idx).syl_ID), 'FontSize', 6, 'Interpreter', 'none');
        
        % plot spike
        spike_iv = sound_spike(idx).spike_iv;
        eLong = soundLong;
        eLong(1:length(spike_iv)) = spike_iv;
        ax2 = axs(2, ii); cla(ax2); hold(ax2, 'on');
        spike_i = find(eLong==1);
        rel_t = (1:length(eLong)) / fs;
        for iii=1:length(spike_i)
          x = rel_t(spike_i(iii));
          plot(ax2, [x x], [-0.25 0.25], 'LineStyle', '-', 'Color', 'red', 'LineWidth', 1);
        end
        ylim(ax2, [-0.25 0.25]);
        xlim(ax2, [rel_t(1) rel_t(end)]);
        axis(ax2, 'off');
        linkaxes([ax1 ax2], 'x');
        
        % plot trajectories of the syllable; overlay spikes
        ax3 = axs(3, ii); cla(ax3); hold(ax3, 'on');
        % plot the background first
        xlim(ax3, x_lim); ylim(ax3, y_lim);
        % add the scatter background
        img = img_list{di};
        image(ax3, ax3.XLim, ax3.YLim, flipud(img));
        % then plot trajectories
        ustart = spike_this.ustart(i_in_u);
        uend = spike_this.uend(i_in_u);
        u = umap_res(ustart:uend,:);
        plot(ax3, u(:,1), u(:,2), 'LineStyle', '-', 'LineWidth', 2, 'Color', '#a6761d');
        % mark the start and end
        scatter(ax3, u(1,1), u(1,2), 80, 'Marker', '^', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'none');
        scatter(ax3, u(end,1), u(end,2), 80, 'Marker', 's', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none');
        % finally mark the spikes
        mat_loc = spike_this.mat_loc{i_in_u};
        if ~isempty(mat_loc)
          scatter(ax3, u(mat_loc,1), u(mat_loc,2), 80, 'filled', 'MarkerFaceColor', '#e41a1c', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8);
        end
        
      end
      % save figure
      fd_save_call = fullfile(fd_save_this, plot_suffix{di});
      if ~exist(fd_save_call, 'dir'); mkdir(fd_save_call); end
      fn_fig = fullfile(fd_save_call, sprintf('%s.%s.fig%d.pdf', birdID, plot_suffix{di}, fig_i));
      print(fig, fn_fig, '-dpdf', '-painters');
      
      close(fig);
    end
  end
  clear to_plot_d d_plot sound_spike axs power powerThis eLong spike_iv u
  
end





















