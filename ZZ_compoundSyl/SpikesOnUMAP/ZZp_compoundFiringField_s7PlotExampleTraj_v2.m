% calculate and plot the firing field of MO neurons in the acoustic space
% step 7: plot example trajectories for selected neurons to use in the main figure
% differ from v1: plot more examples for neuron 20240917-ch13

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
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
% what's the range to include spikes: default only consider spikes within 1-window of syllable boundaries
epad_frame = win_frame;
% how much extra premotor lag to add
lag_frame = 5;

bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Bird specific folder setting
% where is VAE data located
fd_vae = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
% where is spike-window data located
fd_sliding_loc = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_slidingLoc_lag' num2str(lag_frame)]);
fns_loc = dir(fullfile(fd_sliding_loc, sprintf('%s.*.sliding_loc.replaced2.mat', birdID)));
% where is UMAP results located
fd_umap = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_umap']);
umap_run = 'all_syl';
% load UMAP results
fn_info = fullfile(fd_umap, sprintf('%s.%s.info2.mat', birdID, umap_run));
load(fn_info);
tempInfo = info(:, {'syl_ID', 'lens', 'ustart', 'uend'});
fn_umap = fullfile(fd_umap, sprintf('%s.%s.umap_res.csv', birdID, umap_run));
umap_res = readmatrix(fn_umap);

% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% information regarding spikes and sliding windows
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField_lag' num2str(lag_frame)]);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
load(fn_spike);

% load the background image
fn_img = fullfile(fd_field, sprintf('%s.background.img_all.mat', birdID));
load(fn_img);

% where to save results
fd_traj= fullfile(fd_field, 'plot_comp', 'example_traj2');

% what neurons to plot in the main figure
neu_plot = {'20240830-ch8'};



%% 2. Plot individual spikes and trajectories from selected syllable renditions
% get the call.txt and comp.txt files in each folder
fs = 20000;
pad_sound = 0.032;
% all figures use the same x/y axis limits
x_lim = [0.5 19];
y_lim = [-6.5 18.5];

ni = 1;

% get call and compound data
neuronID = neu_plot{ni};
fd_save_this = fullfile(fd_traj, neuronID);
if ~exist(fd_save_this, 'dir'); mkdir(fd_save_this); end

spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
% merge the umap info table
spike_this = outerjoin(spike_this, tempInfo, 'Keys', 'syl_ID', 'Type', 'left', 'MergeKeys', true);
% divide into calls and compound syllables
calls = spike_this(strcmp(spike_this.category, 'v'),:);
compounds = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
compounds = compounds(compounds.dur>=0.3,:);

% what syllables to plot
% fn_sel_call = fullfile(fd_field, 'plot_comp', 'spike_spectrogram', neuronID, 'call.txt');
fn_sel_call = fullfile(fd_field, 'plot_comp', 'spike_spectrogram', neuronID, 'call3.txt');
idx_call = str2double(strsplit(readlines(fn_sel_call), ','));
fn_sel_comp = fullfile(fd_field, 'plot_comp', 'spike_spectrogram', neuronID, 'comp3.txt');
idx_comp = str2double(strsplit(readlines(fn_sel_comp), ','));
sel_call = calls(idx_call,:);
sel_comp = compounds(idx_comp,:);

% grab sound and ephys data
spike_shape = info_neu{strcmp(info_neu.neuronID, neuronID), 'spike_shape'};
spike_shape = str2double(spike_shape{1});
sound_spike_call = ZZfunc_grabEphysDataForTraj_v1(sel_call, birdID, pairID, spike_shape, pad_sound, fs, fd_base);
sound_spike_comp = ZZfunc_grabEphysDataForTraj_v1(sel_comp, birdID, pairID, spike_shape, pad_sound, fs, fd_base);

% plot
% embed syllables in a matrix with max visible duration, unit is second
max_dur = 1;
% get the dims of such a long sound
soundLong = zeros(floor(max_dur*fs),1);
[powerLong, ~, ~, ~, ~, tLong] = getAudioSpectrogramZZ_flexible_v1(soundLong, fs,  256, 256, 236, [1000 5000], [12 23]);

% plot calls, then compound syllables
to_plot_d = {sound_spike_call, sound_spike_comp};
idx_all = {idx_call, idx_comp};
plot_suffix = {'calls', 'compounds'};
img_list = {img_all.v, img_all.comp};

close all;
for di=1:size(to_plot_d, 2)
  
  % how many figures to generate
  sound_spike = to_plot_d{di};
  idx_this = idx_all{di};
  num_fig = floor(size(sound_spike,2) / 2);
  
  for fig_i=1:num_fig
    [fig, axs] =  generatePanelGrid_v3(3, 2, [0.2;0.03;0.6], [0.02;0.02], [0.05;0.05], [0.05;0.01], 0.03, [0;0;0], [10 50 1700 900], false);
    for ii=1:2
      idx = (fig_i-1) * 2 + ii;
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
      title(ax1, sprintf('%d %s %s', idx_this(idx), label_this, sound_spike(idx).syl_ID), 'FontSize', 6, 'Interpreter', 'none');
      
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
      axis(ax2, 'off');
      
      % set the same x-axis for panel 1 and 2
      xlim(ax1, [0 max_dur]);
      xlim(ax2, [0 max_dur]);
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
    
   

%% 3. Plot trajectories of many renditions in the same figure
% what syllables to plot
fn_sel_call = fullfile(fd_field, 'plot_comp', 'spike_spectrogram', neuronID, 'call2.txt');
idx_call = str2double(strsplit(readlines(fn_sel_call), ','));
fn_sel_comp = fullfile(fd_field, 'plot_comp', 'spike_spectrogram', neuronID, 'comp4.txt');
idx_comp = str2double(strsplit(readlines(fn_sel_comp), ','));
sel_call = calls(idx_call,:);
sel_comp = compounds(idx_comp,:);

% grab sound and ephys data
spike_shape = info_neu{strcmp(info_neu.neuronID, neuronID), 'spike_shape'};
spike_shape = str2double(spike_shape{1});
sound_spike_call = ZZfunc_grabEphysDataForTraj_v1(sel_call, birdID, pairID, spike_shape, pad_sound, fs, fd_base);
sound_spike_comp = ZZfunc_grabEphysDataForTraj_v1(sel_comp, birdID, pairID, spike_shape, pad_sound, fs, fd_base);

% plot
% embed syllables in a matrix with max visible duration, unit is second
max_dur = 1;
% get the dims of such a long sound
soundLong = zeros(floor(max_dur*fs),1);
[powerLong, ~, ~, ~, ~, tLong] = getAudioSpectrogramZZ_flexible_v1(soundLong, fs,  256, 256, 236, [1000 5000], [12 23]);

% plot calls, then compound syllables
to_plot_d = {sound_spike_call, sound_spike_comp};
plot_suffix = {'calls', 'compounds'};
img_list = {img_all.v, img_all.comp};

close all;
[fig, axs] =  generatePanelGrid_v3(1, 2, [0.8], [], [0.1;0.05], [0.05;0.01], 0.03, [0], [10 50 1000 600], true);
for di=1:size(to_plot_d, 2)
  sound_spike = to_plot_d{di};
  ax = axs(di); cla(ax); hold(ax, 'on');
  % first the background embedding
  xlim(ax, x_lim); ylim(ax, y_lim);
  img = img_list{di};
  image(ax, ax.XLim, ax.YLim, flipud(img));
  % then loop through renditions, plot trajectories
  for idx=1:size(sound_spike, 2)
    i_in_u = find(strcmp(spike_this.syl_ID, sound_spike(idx).syl_ID));
    ustart = spike_this.ustart(i_in_u);
    uend = spike_this.uend(i_in_u);
    u = umap_res(ustart:uend,:);
    plot(ax, u(:,1), u(:,2), 'LineStyle', '-', 'LineWidth', 1, 'Color', [0 0.5 0 0.5]);
    % mark the start and end
    scatter(ax, u(1,1), u(1,2), 30, 'Marker', '^', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none');
    scatter(ax, u(end,1), u(end,2), 30, 'Marker', 's', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none');
  end
  % then overlay spikes
  for idx=1:size(sound_spike, 2)
    i_in_u = find(strcmp(spike_this.syl_ID, sound_spike(idx).syl_ID));
    ustart = spike_this.ustart(i_in_u);
    uend = spike_this.uend(i_in_u);
    u = umap_res(ustart:uend,:);
    mat_loc = spike_this.mat_loc{i_in_u};
    if ~isempty(mat_loc)
      scatter(ax, u(mat_loc,1), u(mat_loc,2), 80, 'filled', 'MarkerFaceColor', '#e41a1c', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
    end
  end
  xlim(ax, x_lim); ylim(ax, y_lim);
  title(ax, sprintf('%s (%s)', plot_suffix{di}, neuronID), 'FontSize', 14);
end
% save results
fn_fig = fullfile(fd_save_this, sprintf('%s.%s.exampleTraj.pdf', birdID, neuronID));
print(fig, fn_fig, '-dpdf', '-painters');
















    
    
    
    
    
    
    
    
