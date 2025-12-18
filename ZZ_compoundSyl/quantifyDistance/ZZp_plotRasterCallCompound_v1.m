% for a given neuron, plot the raster aligned to syllable onset for calls and compound syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));


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


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Bird specific folder setting
% load information regarding sparse neurons
fn_info_neu = fullfile(fd_base, 'DbaseFiles', pairID, 'MetaInfo', sprintf('%s_sparseInfo.mat', birdID));
a = load(fn_info_neu); info_neu = a.info;

% information regarding spikes and sliding windows
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField']);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
load(fn_spike);

% where to save results
fd_raster= fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_raster']);
if ~exist(fd_raster, 'dir'); mkdir(fd_raster); end



%% 2. Loop through neurons, plot raster: combine calls and compound syllables
syl_call = {{'v1'}, {'v2'}, {'v3'}, {'v4'}, {'v5'}, {'v6'}, {'v7', 'V'}};  % good call clusters then non-clusters
label_call = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v0'};
% what colors to use
col_call = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#737373'};
rends_plot_call = 50;  % max number of renditions to plot for each category
rends_plot_comp = 200;
fs = 20000;


% ni = 44;
for ni=1:size(info_neu,1)
  
  neuronID = info_neu.neuronID{ni};
  disp(neuronID);
  
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  % divide into calls and compound syllables
  calls = spike_this(strcmp(spike_this.category, 'v'),:);
  compounds = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
  compounds = compounds(compounds.dur>=0.3,:);
  
  close all;
  fig = ZZfunc_newFigurePDFsize_v1([10 10 600 800]);
  ax = gca; hold(ax, 'on');
  set(ax, 'YDir', 'reverse');
  % first plot calls
  y_start = 0;
  x_lim = [-0.05 1];
  rng(1922);
  for si=1:size(syl_call, 2)
    v = syl_call{si};
    spike = calls(ismember(calls.label, v),:);
    act_num = min([rends_plot_call size(spike,1)]);
    i_rd = randsample(1:size(spike,1), act_num);
    d_plot = spike(i_rd,:);
    % plot
    [ax, y_start] = ZZfunc_plotRasterNoWrap_v2(ax, d_plot, y_start, fs, 's', 10, 'black', 1, x_lim, true, col_call{si}, true, label_call{si});
  end
  
  % then plot compound syllable
  act_num = min([rends_plot_comp size(compounds,1)]);
  i_rd = randsample(1:size(compounds,1), act_num);
  d_plot = compounds(i_rd,:);
  [ax, y_start] = ZZfunc_plotRasterNoWrap_v2(ax, d_plot, y_start, fs, 's', 10, 'black', 1, x_lim, true, 'black', true, 'CS');
  
  % save figure
  ylim(ax, [0 y_start+1]);
  yticks(ax, []);
  xlabel(ax, 'Rel. time (sec)', 'FontSize', 12);
  title(ax, neuronID, 'FontSize', 12);
  fn_fig = fullfile(fd_raster, sprintf('%s.%s.compound_raster.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
  
end



%% 3. Separate plots for calls and compound syllables
syl_call = {{'v1'}, {'v2'}, {'v3'}, {'v4'}, {'v5'}, {'v6'}, {'v7', 'V'}};  % good call clusters then non-clusters
label_call = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v0'};
% what colors to use
col_call = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#737373'};
rends_plot_call = 50;  % max number of renditions to plot for each category
rends_plot_comp = 200;
fs = 20000;

% ni = 44;
for ni=1:size(info_neu,1)
  
  neuronID = info_neu.neuronID{ni};
  disp(neuronID);
  
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  % divide into calls and compound syllables
  calls = spike_this(strcmp(spike_this.category, 'v'),:);
  compounds = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
  compounds = compounds(compounds.dur>=0.3,:);
  
  close all;
  % first plot calls
  fig = ZZfunc_newFigurePDFsize_v1([10 10 300 500]);
  ax = gca; hold(ax, 'on');
  set(ax, 'YDir', 'reverse');
  y_start = 0;
  x_lim = [-0.05 0.25];
  rng(1922);
  for si=1:size(syl_call, 2)
    v = syl_call{si};
    spike = calls(ismember(calls.label, v),:);
    act_num = min([rends_plot_call size(spike,1)]);
    i_rd = randsample(1:size(spike,1), act_num);
    d_plot = spike(i_rd,:);
    % plot
    [ax, y_start] = ZZfunc_plotRasterNoWrap_v3(ax, d_plot, y_start, fs, 's', 10, 'red', 1, x_lim, true, col_call{si}, true, label_call{si}, 0.03, '#737373', true);
  end
  % save figure
  ylim(ax, [0 y_start+1]);
  yticks(ax, []);
  xlabel(ax, 'Rel. time (sec)', 'FontSize', 12);
  title(ax, neuronID, 'FontSize', 12);
  fn_fig = fullfile(fd_raster, sprintf('%s.%s.raster1.calls.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
  close(fig);
  
  % then plot compound syllable
  x_lim = [-0.05 1];
  fig = ZZfunc_newFigurePDFsize_v1([10 10 500 500]);
  ax = gca; hold(ax, 'on');
  set(ax, 'YDir', 'reverse');
  y_start = 0;
  act_num = min([rends_plot_comp size(compounds,1)]);
  i_rd = randsample(1:size(compounds,1), act_num);
  d_plot = compounds(i_rd,:);
  [ax, y_start] = ZZfunc_plotRasterNoWrap_v3(ax, d_plot, y_start, fs, 's', 10, 'red', 1, x_lim, true, 'black', true, 'CS', 0.07, '#737373', true);
  
  % save figure
  ylim(ax, [0 y_start+1]);
  yticks(ax, []);
  xlabel(ax, 'Rel. time (sec)', 'FontSize', 12);
  title(ax, neuronID, 'FontSize', 12);
  fn_fig = fullfile(fd_raster, sprintf('%s.%s.raster2.compound.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
  
end



%% 4. Plot relative time
syl_call = {{'v1'}, {'v2'}, {'v3'}, {'v4'}, {'v5'}, {'v6'}, {'v7', 'V'}};  % good call clusters then non-clusters
label_call = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v0'};
% what colors to use
col_call = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#737373'};
rends_plot_call = 50;  % max number of renditions to plot for each category
rends_plot_comp = 200;
fs = 20000;

% ni = 44;
for ni=1:size(info_neu,1)
  
  neuronID = info_neu.neuronID{ni};
  disp(neuronID);
  
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  % divide into calls and compound syllables
  calls = spike_this(strcmp(spike_this.category, 'v'),:);
  compounds = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
  compounds = compounds(compounds.dur>=0.3,:);
  
  close all;
  % first plot calls
  fig = ZZfunc_newFigurePDFsize_v1([10 10 300 500]);
  ax = gca; hold(ax, 'on');
  set(ax, 'YDir', 'reverse');
  y_start = 0;
  x_lim = [-0.05 1.05];
  rng(1922);
  for si=1:size(syl_call, 2)
    v = syl_call{si};
    spike = calls(ismember(calls.label, v),:);
    act_num = min([rends_plot_call size(spike,1)]);
    i_rd = randsample(1:size(spike,1), act_num);
    d_plot = spike(i_rd,:);
    % plot
    [ax, y_start] = ZZfunc_plotRasterRelTime_v1(ax, d_plot, y_start, 's', 10, 'red', 1, x_lim, true, col_call{si}, true, label_call{si}, 0.1, '#737373', true);
  end
  % save figure
  ylim(ax, [0 y_start+1]);
  yticks(ax, []);
  xlabel(ax, 'Rel. time', 'FontSize', 12);
  title(ax, neuronID, 'FontSize', 12);
  fn_fig = fullfile(fd_raster, sprintf('%s.%s.raster1.callsRel.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
  close(fig);
  
  % then plot compound syllable
  x_lim = [-0.05 1];
  fig = ZZfunc_newFigurePDFsize_v1([10 10 500 500]);
  ax = gca; hold(ax, 'on');
  set(ax, 'YDir', 'reverse');
  y_start = 0;
  act_num = min([rends_plot_comp size(compounds,1)]);
  i_rd = randsample(1:size(compounds,1), act_num);
  d_plot = compounds(i_rd,:);
  [ax, y_start] = ZZfunc_plotRasterRelTime_v1(ax, d_plot, y_start, 's', 10, 'red', 1, x_lim, true, 'black', true, 'CS', 0.1, '#737373', true);
  
  % save figure
  ylim(ax, [0 y_start+1]);
  yticks(ax, []);
  xlabel(ax, 'Rel. time (sec)', 'FontSize', 12);
  title(ax, neuronID, 'FontSize', 12);
  fn_fig = fullfile(fd_raster, sprintf('%s.%s.raster2.compoundRel.pdf', birdID, neuronID));
  print(fig, fn_fig, '-dpdf', '-painters');
  
end












