% given identified primitives, plot the time location of primitives as a raster plot


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
clims = {[10 21]; [10 21]; [10.5 21.5]; [12 23]};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% size of sliding window
win_size = 32;  % in unit of frames
sec_per_frame = 0.001;
fs = 20000;

% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
clim = clims{bi};


% only look at well-defined calls
call_syl = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'};
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
fd_save_ref = fullfile(fd_save, 'compound_call3');
% load information about the compound syllable
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);

vi = 1;
% load the strands info
fn_strand = fullfile(fd_save_ref, sprintf('call_v%d.strands.mat', vi));
load(fn_strand);
% add category label 
for ri=1:size(strands, 2)
  tar_i = strands(ri).tar_i;
  strands(ri).category = comp.call_subtype{tar_i};
end


%% Plot the old style raster plot
num_prim = max([strands.prim_belong]);
size_d = strands(1).size_d;
col_prim = lines(num_prim);
close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 250 600]);
ax=gca; hold(ax, 'on');
for prim_i=1:num_prim
  % prim_i = 2;
  strands_plot = strands([strands.prim_belong]==prim_i);
  scatter(ax, [strands_plot.x_center], [strands_plot.y_center], 20, 'Marker', 'o', 'MarkerFaceColor', col_prim(prim_i,:), 'MarkerEdgeColor', 'none');
end
xlim(ax, [1 size_d(2)]);
ylim(ax, [1 800]);
xlabel('Latency in ref. (ms)', 'FontSize', 12);
ylabel('Latency in comp. (ms)', 'FontSize', 12);



%% plot as raster plot
% each occurence of a given primitive in compound syllable is shown as a blue tick
% the end of the syllable as grey tick
% what syllables to check
syls = {'b', 'x'};
for prim_i=1:num_prim
  close all;
  fig = ZZfunc_newFigurePDFsize_v1([10 10 300 600]);
  ax=gca; hold(ax, 'on');
  strands_plot = strands(([strands.prim_belong]==prim_i) & (ismember({strands.category}, syls)));
  
  % sort by the syllable duration
  dur = cellfun(@(x) x(1), {strands_plot.size_d});
  [~, i_sort] = sort(dur, 'descend');
  strands_plot = strands_plot(i_sort);
  
  for ri=1:size(strands_plot,2)
%   for ri=1:30
    scatter(ax, strands_plot(ri).y_center / 1000, ri, 20, 'Marker', 's', 'MarkerFaceColor', '#e7298a', 'MarkerEdgeColor', 'none');
%     plot(ax, strands_plot(ri).y_range / 1000, [ri ri], 'Color', '#e7298a', 'LineStyle', '-', 'LineWidth', 2);
    scatter(ax, strands_plot(ri).size_d(1) / 1000, ri, 15, 'Marker', 's', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none');
  end
  
  xlim(ax, [0 1]);
  ylim(ax, [0 size(strands_plot,2)]);
  xlabel(ax, 'T-onset (sec)', 'FontSize', 12);
  ylabel(ax, 'ID of similarity strands', 'FontSize', 12);
  
  fn_fig = fullfile(fd_save_ref, sprintf('call_v%d.raster.prim%d.pdf', vi, prim_i));
  print(fig, fn_fig, '-dpdf', '-painters');
  
end


%% plot as raster plot along with histrogram
% each occurence of a given primitive in compound syllable is shown as a blue tick
% the end of the syllable as grey tick
% what syllables to check
syls = {'b', 'x'};
for prim_i=1:num_prim
  close all; 
  [fig, axs] = generatePanelGrid_v2(2, 1, [0.75;0.1], [0.02], [0.05;0.05], [0.2;0.05], 0.05, [0;0], [10 10 300 650]);
  strands_plot = strands(([strands.prim_belong]==prim_i) & (ismember({strands.category}, syls)));
  
  % sort by the syllable duration
  dur = cellfun(@(x) x(1), {strands_plot.size_d});
  [~, i_sort] = sort(dur, 'descend');
  strands_plot = strands_plot(i_sort);
  
  ax = axs(1); hold(ax, 'on');
  for ri=1:size(strands_plot,2)
%   for ri=1:30
    scatter(ax, strands_plot(ri).y_center / 1000, ri, 20, 'Marker', 's', 'MarkerFaceColor', '#e7298a', 'MarkerEdgeColor', 'none');
%     plot(ax, strands_plot(ri).y_range / 1000, [ri ri], 'Color', '#e7298a', 'LineStyle', '-', 'LineWidth', 2);
    scatter(ax, strands_plot(ri).size_d(1) / 1000, ri, 15, 'Marker', 's', 'MarkerFaceColor', '#737373', 'MarkerEdgeColor', 'none');
  end
  xlim(ax, [0 1]);
  xticks(ax, []);
  ylim(ax, [0 size(strands_plot,2)]);
  ylabel(ax, 'ID of similarity strands', 'FontSize', 12);
  
  % then plot histogram
  ax2 = axs(2); cla(ax2);
  histogram(ax2, [strands_plot.y_center]/1000, 'BinEdges', 0:0.05:1, 'EdgeColor', 'none', 'FaceColor', '#e7298a');
  xlim(ax2, [0 1]);
  ylabel(ax2, 'Counts', 'FontSize', 12);
  xlabel(ax2, 'T-onset (sec)', 'FontSize', 12);
  
  fn_fig = fullfile(fd_save_ref, sprintf('call_v%d.rasterHisto.prim%d.pdf', vi, prim_i));
  print(fig, fn_fig, '-dpdf', '-painters');
  
end

















