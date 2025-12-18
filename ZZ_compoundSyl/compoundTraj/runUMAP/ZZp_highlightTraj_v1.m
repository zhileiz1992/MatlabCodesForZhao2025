% highlight a few trajectories of compound syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% coloring for syllables
syl_full = {'v', 'h', 'e', 'b', 'x'};
col_full = {'#e41a1c', '#984ea3', '#4daf4a', '#377eb8', '#737373'};
col_dict = struct;
for si=1:size(syl_full,2) col_dict.(syl_full{si}) = col_full{si}; end
% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load results from UMAP
fd_umap_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'UMAPcomp');
umap_run = 'v.b.h.e.x.n-1';
syls = {'v', 'b', 'h', 'e', 'x'};
fd_umap = fullfile(fd_umap_base, umap_run);

fn_info = fullfile(fd_umap, sprintf('%s.%s.info.csv', birdID, umap_run));
info = readtable(fn_info, 'Delimiter', ',');
fn_embed = fullfile(fd_umap, sprintf('%s.%s.embedding.csv', birdID, umap_run));
embed = readmatrix(fn_embed);

% where to save plots
fd_save = fullfile(fd_umap, 'trajHighlight');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end


%% 2. Find regions with high density for compound syllables
% restrict to batch 1 and compound syllable
fs = 20000;
dur_thre = 0.3;  % unit is sec
comp = info(info.duration>=dur_thre & ismember(info.call_subtype, {'b', 'x'}) & ismember(info.batch, {'batch1'}), :);
% comp = info(info.duration>=dur_thre & ismember(info.call_subtype, {'x'}) & ismember(info.batch, {'batch1'}), :);
% get the embedding for these selected syllables; note that python is 0-based, but matlab is 1-based
for ri=1:size(comp, 1)
  i1 = comp.count_start(ri) + 1; 
  i2 = comp.count_end(ri);
  comp.umap1{ri} = embed(i1:i2, 1);
  comp.umap2{ri} = embed(i1:i2, 2);
end

close all;
[fig, axes] = generatePanelGrid_v2(1, 3, [0.8], [], [0.05;0.05], [0.05;0.05], 0.05, [0], [10 10 1900 550]);

% plot the overall distribution
ax1 = axes(1); hold(ax1, 'on');
num_plot = 500; 
x_lim = [0 20];
y_lim = [-6.5366 18.1853];
rng(1992);
idx_r = randsample(1:size(comp,1), num_plot);
x = vertcat(comp.umap1{idx_r});
y = vertcat(comp.umap2{idx_r});
scatter(ax1, x, y, 10, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none'); 
xlim(ax1, x_lim);
ylim(ax1, y_lim);
colorbar(ax1);
title(ax1, sprintf('Compound syllables, %d plotted', num_plot), 'FontSize', 12);

% convert to density
ax2 = axes(2); hold(ax2, 'on'); cla(ax2);
nbins = 100;
xx = linspace(x_lim(1), x_lim(2), nbins+1); 
yy = linspace(y_lim(1), y_lim(2), nbins+1); 
[N, edgesX, edgesY] = histcounts2(x, y, xx, yy);
imagesc(ax2, xx(1:nbins), yy(1:nbins), N');
set(ax2, 'YDir', 'normal');
xlim(ax2, x_lim);
ylim(ax2, y_lim);
colormap(custom_map);
colorbar(ax2);
title(ax2, 'Density plot', 'FontSize', 12);

% find a point to plot trajectories that go through it
ix = 31; 
iy = 33;
% boundaries for x/y
xb = [xx(ix) xx(ix+1)];
yb = [yy(iy) yy(iy+1)];
% what syllables has windows in this range
loc_pass = [];
count = 0;
for ri=1:size(comp, 1)
  temp = [comp.umap1{ri} comp.umap2{ri}];
  ipass = find((temp(:,1)>=xb(1)) & (temp(:,1)<=xb(2)) & (temp(:,2)>=yb(1)) & (temp(:,2)<=yb(2)));
  if ~isempty(ipass)
    count = count+1;
    loc_pass(count).idx = ri;
    loc_pass(count).loc = ipass;
  end
end

plot(ax1, xx(ix), yy(iy), 'Marker', '+', 'MarkerSize', 30, 'Color', 'red');
plot(ax2, xx(ix), yy(iy), 'Marker', '+', 'MarkerSize', 30, 'Color', 'red');
ax3 = axes(3); cla(ax3); hold(ax3, 'on');
% scatter(ax3, x, y, 10, 'Marker', 'o', 'MarkerFaceColor', '#737373', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none'); 
% plot trajectories
to_plot = 5; 
rng(1118);
i_plot = randsample([loc_pass.idx], to_plot);
col_lines = lines(to_plot);
for ii=1:to_plot
  ri = i_plot(ii);
  plot(ax3, comp.umap1{ri}, comp.umap2{ri}, '-', 'Color', [col_lines(ii,:) 0.5], 'LineWidth', 2); 
end
plot(ax3, xx(ix), yy(iy), 'Marker', '+', 'MarkerSize', 30, 'Color', 'red');
xlim(ax3, x_lim);
ylim(ax3, y_lim);
colorbar(ax3);
title(ax3, sprintf('Total %d pass (%.1f%%), only %d plotted', size(loc_pass,2), size(loc_pass,2)/size(comp,1)*100, to_plot), 'FontSize', 12);


% plot the syllable spectrograms and the corresponding window



%% 3. Calculate the latency to syllable onset
for ii=1:size(loc_pass,2)
  % what's the onset/offset of the syllable?
  onset = 32; 
  temp = comp(loc_pass(ii).idx, :);
  offset = temp.num_slides;
  % convert to real time
  sec_per_frame = 0.001;
  rel_t = ((0:(offset-1))-onset) * sec_per_frame; 
  loc_pass(ii).rel_t = rel_t; 
  loc_pass(ii).pass_t = rel_t(loc_pass(ii).loc);
  loc_pass(ii).dur = rel_t(end);
end

% plot
% first sort by duration
[~, sort_i] = sort([loc_pass.dur], 'descend');
loc_pass2 = loc_pass(sort_i);
% black tick at syllable onset, blue tick at syllable offset, red tick at the passed window
blue_loc = [];
red_loc = [];
for ii=1:size(loc_pass2,2)
  rel_t = loc_pass2(ii).rel_t;
  blue_loc = [blue_loc; [rel_t(end) ii]];
  pass_t = loc_pass2(ii).pass_t;
  for jj=1:length(pass_t)
      red_loc = [red_loc; [pass_t(jj) ii]];
  end
end

close all;
[fig, axes] = generatePanelGrid_v2(2, 1, [0.6 0.2], [0.075], [0.05;0.05], [0.1;0.05], 0.05, [0;0], [10 10 600 600]);
ax = axes(1); hold(ax, 'on');
scatter(ax, blue_loc(:,1), blue_loc(:,2), 10, 'Marker', 's', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'none');
scatter(ax, red_loc(:,1), red_loc(:,2), 10, 'Marker', 's', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none');
% xline(ax, 0, 'LineStyle', '-', 'Color', 'black', 'LineWidth', 1);
xlim(ax, [0 max(blue_loc(:,1)+0.05)]);
ylim(ax, [0 size(loc_pass,2)+1]);
xlabel(ax, 'Latency to syllable onset (sec)', 'FontSize', 12);
ylabel(ax, 'Rendition ID', 'FontSize', 12);
title(ax, 'Latency of a tiny region in UMAP', 'FontSize', 12); 
% then a histogram
ax2 = axes(2); cla(ax2);
h = histogram(ax2, red_loc(:,1), 'BinWidth', 0.02);
h.FaceColor = 'red';  % Set bar fill color to blue
h.EdgeColor = 'none'; % Set bar outline color to black
xlabel(ax2, 'Latency to syllable onset (sec)', 'FontSize', 12);
ylabel(ax2, 'Counts', 'FontSize', 12);
linkaxes(axes, 'x');

































