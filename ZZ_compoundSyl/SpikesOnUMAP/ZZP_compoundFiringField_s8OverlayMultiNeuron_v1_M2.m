% overlay spikes of many neurons in the same embedding plot to show that neurons tile the acoustic space
% criteria for selecting neurons to plot:
% 1. basic firing metrics? e.g. sparseness, IRCC
% 2. correlate/encode acoustics? Dip in acoustic variation after bursts; Significant in the ensemble modeling
% for bird M2

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



bi = 2;
birdID = birdIDs{bi};
pairID = pairIDs{bi};

%% Bird specific folder setting
% where is VAE data located
fd_vae = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', vae_run);
% where is spike-window data located
% fd_sliding_loc = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_slidingLoc']);
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

% load pre-calculated sliding window assignment
fd_field = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_firingField_lag' num2str(lag_frame)]);
fn_spike = fullfile(fd_field, sprintf('%s.spike_win.mat', birdID));
a = load(fn_spike); spike_win = a.spike_win;


% where are ensemble modeling results
fd_model = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'ensembleModel');

% where are basic neuron firing metrics stored
fd_hahn = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'popRaster2');

% where to save results
fd_save = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', [vae_run '_multiNeuron']);
if ~exist(fd_save, 'dir'); mkdir(fd_save); end


%% 1. Load the pre-calculated UMAP background
fn_img = fullfile(fd_field, sprintf('%s.background.img_all.mat', birdID));
load(fn_img);


%% 2. Select neurons that burst sparsely in a call subtype
v = 'v1';
% load sparse neurons that burst for this call 
fn_metric = fullfile(fd_hahn, v, sprintf('Hahnloser-%s-chan0.criteria7.mat', v));
a=load(fn_metric);  criteria=a.criteria;
% subset to those that burst
crit = struct2table(criteria([criteria.isPass]==1));

% load ensemble modeling results
for ni=1:size(crit,1)
  neuronID = crit.neuronID{ni};
  fn_model = fullfile(fd_model, 'XGB_res', neuronID, sprintf('%s.%s.metrics_all.mat', birdID, neuronID));
  a=load(fn_model); temp=a.metrics_all;
  fdns = fieldnames(temp);
  if ismember('Compound_Acoustic', fdns); crit.xgb_compound(ni) = median(temp.Compound_Acoustic.improvement_over_null); end
  if ismember('Call_Acoustic', fdns); crit.xgb_call(ni) = median(temp.Call_Acoustic.improvement_over_null); end
  if ismember('v1_Acoustic', fdns); crit.xgb_v1(ni) = median(temp.v1_Acoustic.improvement_over_null); end
end

% what neurons have good modeling for both calls and compound syls
crit2 = crit((crit.xgb_call>=15) & (crit.xgb_compound>=15),:);


% load manually selected neurons
fn_sel = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', sprintf('%s_selected.txt', birdID));
neuron_sel = readtable(fn_sel, 'ReadVariableNames', false);



%% 4. Plot the density instead
% x_lim = [0.5 19]; y_lim = [-6.5 18.5]; % for bird M1
x_lim = [-8 10]; y_lim = [-10 16];  % bird M2
num_bin = 100;
% first calculate density for each selected neurons and color range '20240912-ch10' '20240915-ch4A' '20240917-ch12B' '20240914-ch11A' 
% neuron_sel = {'20240917-ch13', '20240830-ch8', '20240903-ch6', '20240903-ch11', '20240917-ch4', '20240928-ch5', '20240917-ch1', '20240919-ch8'};
% neuron_sel = crit.neuronID';
% neuron_color = {'#e41a1c','#a65628','#4daf4a','#984ea3','#ff7f00','#00aeef','#f781bf','#525252'};
% neuron_color = neuron_color(1:size(neuron_sel,2));
% color_rgb = cellfun(@(x) hexToRGB(x), neuron_color, 'UniformOutput', false);
neuron_sel = neuron_sel{:,1};
neuron_sel = neuron_sel';
cmap_array = {};
colors = turbo(size(neuron_sel, 2));
% colors = flipud(colors);
neuron_color = {};
for ci=1:size(colors, 1)
  x = colors(ci,:);
  a = light_to_color_lab(x, 256);
  cmap_array{ci} = a;
  neuron_color{ci} = x;
end
% col_temp = lines(size(neuron_sel, 2));
% color_rgb = {};
% for ri=1:size(col_temp,1); color_rgb{ri}=col_temp(ri,:); end
% cmap_array = cellfun(@(x) light_to_color_lab(x,256), color_rgb, 'UniformOutput', false);
density_all1 = cell(1, size(neuron_sel,2));
clim_array1 = cell(1, size(neuron_sel,2));
density_all2 = cell(1, size(neuron_sel,2));
clim_array2 = cell(1, size(neuron_sel,2));
% ni =16;
% min_dens = 
for ni=1:size(neuron_sel,2)
%   neuronID = neuron_sel.Var1{ni};
  neuronID = neuron_sel{ni};
  spike_this = spike_win(strcmp(spike_win.neuronID, neuronID),:);
  spike_this = outerjoin(spike_this, tempInfo, 'Keys', 'syl_ID', 'Type', 'left', 'MergeKeys', true);
  
  % plot for calls
  spike1 = spike_this(strcmp(spike_this.category, 'v'),:);
%   spike1 = spike_this(strcmp(spike_this.label, 'v1'),:);
  spike_loc1 = ZZfunc_plotSpikeOnUmapComp_v2_noplot(ax1, spike1, umap_res);
  [density1, x_edges1, y_edges1, x_centers1, y_centers1] = ZZfunc_calcSpikeDensity_count_V1(spike_loc1(:,1), spike_loc1(:,2), num_bin, x_lim, y_lim, 3);
  density1 = density1';
  d_qt1 = quantile(density1(:), 0.99); 
%   d_qt1 = max([1 d_qt1]); 
  density_all1{ni} = density1;
  clim_array1{ni} = [0 d_qt1];
  
  % plot for compound syllables
  spike2 = spike_this(ismember(spike_this.category, {'b', 'x'}),:);
  spike2 = spike2(spike2.dur>=0.3, :);
  spike_loc2 = ZZfunc_plotSpikeOnUmapComp_v2_noplot(ax2, spike2, umap_res);
  [density2, x_edges2, y_edges2, x_centers2, y_centers2] = ZZfunc_calcSpikeDensity_count_V1(spike_loc2(:,1), spike_loc2(:,2), num_bin, x_lim, y_lim, 3);
  density2 = density2';
  d_qt2 = quantile(density2(:), 0.99); 
%   d_qt2 = max([1 d_qt2]); 
  density_all2{ni} = density2;
  clim_array2{ni} = [0 d_qt2];
end

close all;
% add neuron one-by-one to the plot
for ni=1:size(neuron_sel, 2)
  % one plot for calls, another for compound syllables
  [fig, axs] =  generatePanelGrid_v3(1, 2, [0.8], [], [0.1;0.05], [0.05;0.01], 0.03, [0], [10 50 1000 600], true);
  ax1 = axs(1); cla(ax1); hold(ax1, 'on');
  ax2 = axs(2); cla(ax2); hold(ax2, 'on');
  % plot the background first
  xlim(ax1, x_lim); ylim(ax1, y_lim);
  hb1 = image(ax1, ax1.XLim, ax1.YLim, flipud(img_all.v));
  set(hb1, 'AlphaData',0.75);
  xlim(ax2, x_lim); ylim(ax2, y_lim);
  hb2 = image(ax2, ax2.XLim, ax2.YLim, flipud(img_all.comp));
  set(hb2, 'AlphaData',0.75);

  % then overlay the calculate density map
  % overlay_density_images(density_all, cmap_array, 0.25, clim_array);
  [RGB_out1, A_out1] = blend_density_layers(density_all1(1:ni), cmap_array(1:ni), 0.15, clim_array1(1:ni));
  hd1 = image(ax1, ax1.XLim, ax1.YLim, RGB_out1, 'AlphaData', A_out1);
  xlim(ax1, x_lim); ylim(ax1, y_lim);
  [RGB_out2, A_out2] = blend_density_layers(density_all2(1:ni), cmap_array(1:ni), 0.15, clim_array2(1:ni));
  hd2 = image(ax2, ax2.XLim, ax2.YLim, RGB_out2, 'AlphaData', A_out2);
  xlim(ax2, x_lim); ylim(ax2, y_lim);

  % add text
  hold(ax1,'on');
  [ax1] = ZZfunc_addSimpleLegend_v2(ax1, neuron_sel(1:ni), neuron_color(1:ni), 10, zeros(ni,1)+1.5, linspace(10,17,ni));

  fn_fig = fullfile(fd_save, sprintf('%s.firingFields.step%d.pdf', birdID, ni));
  print(fig, fn_fig, '-dpdf', '-painters');
end












