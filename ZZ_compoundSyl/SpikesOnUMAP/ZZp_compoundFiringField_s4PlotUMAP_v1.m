% calculate and plot the firing field of MO neurons in the acoustic space
% step 4: Plot the UMAP embedding without spikes to assess quality

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));

%% General folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what's the sliding window size
win_frame = 32;
% what VAE run to use
vae_run = 'traj_chop_32_1_32';


bi = 2;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird specific folder setting
% load UMAP results
fd_umap =  fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'embed', sprintf('%s_umap', vae_run));
% umap_run = 'all_syl';
umap_run = 'all_syl2';
% umap_run = 'no_eh';
fn_d_info = fullfile(fd_umap, sprintf('%s.%s.info.mat', birdID, umap_run));
load(fn_d_info);
info = struct2table(info);
fn_umap = fullfile(fd_umap, sprintf('%s.%s.umap_res.csv', birdID, umap_run));
umap_res = readmatrix(fn_umap);

% add start/end info
cumlen = cumsum([info.lens]);
info.ustart = [1; cumlen(1:(end-1))+1];
info.uend = cumlen;
% also calculate the duration 
temp = cellfun(@(x) strsplit(x, '_'), info.syl_ID, 'UniformOutput', false);
info.istart = cellfun(@(x) str2double(x{end-1}), temp, 'UniformOutput', false);
info.iend = cellfun(@(x) str2double(x{end}), temp, 'UniformOutput', false);
fs = 20000;
info.dur = (cell2mat(info.iend) - cell2mat(info.istart)) / fs;

fn_d_info = fullfile(fd_umap, sprintf('%s.%s.info2.mat', birdID, umap_run));
save(fn_d_info, 'info', '-v7.3');



%% 1. Plot embedding: all syllables
cat_plot = {'v', 'b', 'x', 'h', 'e'};
cat_color = {'#e41a1c', '#377eb8', '#737373', '#984ea3', '#4daf4a'};
num_plot = 200;
close all;
[fig, axes] =  generatePanelGrid_v2(1, size(cat_plot,2), [0.5], [], [0.1;0.05], [0.05;0.05], 0.02, [0], [10 10 1800 600]);
% marg = 0.25;
% x_lim = [min(umap_res(:,1))-marg max(umap_res(:,1))+marg];
% y_lim = [min(umap_res(:,2))-marg max(umap_res(:,2))+marg];
% x_lim = [0.5 19];
% y_lim = [-6.5 18.5];
% x_lim = [-10 10]; y_lim = [-10 16]; for bird M2, first half
x_lim = [0 17]; y_lim = [-12 16];  %for bird M2, first half
rng(1992);
for si=1:size(cat_plot,2)
  ax = axes(si); hold(ax, 'on');
  idx = find(strcmp(info.category, cat_plot{si}));
  if ~isempty(idx)
    act_num = min([length(idx) num_plot]);
    i_rd = randsample(idx, act_num);
    for ii=1:act_num
      ustart = info.ustart(i_rd(ii));
      uend = info.uend(i_rd(ii));
      u = umap_res(ustart:uend,:);
      scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', cat_color{si}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.05);
    end
    [ax] = ZZfunc_rasterizePlot_v1(ax, x_lim, y_lim);
    title(ax, cat_plot{si}, 'FontSize', 14);
  end
end
% save figures
fn_fig = fullfile(fd_umap, sprintf('%s.%s.syln%d.umap.pdf', birdID, umap_run, num_plot));
print(fig, fn_fig, '-dpdf', '-painters');



%% 2. Plot embedding: calls + compound syllables
% load information regarding compound syllables
fd_comp = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
fn_comp = fullfile(fd_comp, sprintf('%s.comp.mat', birdID));
load(fn_comp);

idx1 = find(strcmp(info.category, 'v'));
idx2 = find(ismember(info.syl_ID, comp.syl_ID));
idx_plot = {idx1; idx2};
titles = {'Calls', 'Compound syllables'};

num_plot = 200;
close all;
[fig, axes] =  generatePanelGrid_v2(1, 2, [0.5], [], [0.1;0.05], [0.05;0.05], 0.05, [0], [10 10 800 700]);
% marg = 0.25;
% x_lim = [min(umap_res(:,1))-marg max(umap_res(:,1))+marg];
% y_lim = [min(umap_res(:,2))-marg max(umap_res(:,2))+marg];
% x_lim = [0.5 19];
% y_lim = [-6.5 18.5];
rng(1992);
for si=1:size(idx_plot, 1)
  ax = axes(si); hold(ax, 'on');
  idx = idx_plot{si};
  act_num = min([length(idx) num_plot]);
  i_rd = randsample(idx, act_num);
  for ii=1:act_num
    ustart = info.ustart(i_rd(ii));
    uend = info.uend(i_rd(ii));
    u = umap_res(ustart:uend,:);
    scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', cat_color{si}, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.05);
  end
  [ax] = ZZfunc_rasterizePlot_v1(ax, x_lim, y_lim);
  title(ax, titles{si}, 'FontSize', 14);
end
% save figures
fn_fig = fullfile(fd_umap, sprintf('%s.%s.syln%d.call_compound.pdf', birdID, umap_run, num_plot));
print(fig, fn_fig, '-dpdf', '-painters');















