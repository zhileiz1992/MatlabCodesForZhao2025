% plot firing density of many neurons in the same UMAP 
% for bird M1, call v1 and v2
% 11/16/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl/SpikesOnUMAP'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
fd_ephys = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'extractedReplaced2');
fns_ephys = dir(fullfile(fd_ephys, sprintf('%s.v*.segments_all.replaced2.mat', birdID)));
% what's the window size
win_frame = 32;
ms_per_frame = 1;
% how much frames when calculating syllable spectrograms
spec_frame = 80;
% what syllable types to analyze
syls = {'v1', 'v2'};
% where is the corresponding VAE/UMAP results located
vae_run = 'traj_chop_32_1_32';
umap_run = 'umapAll.v4v5';
fd_embed =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', vae_run);
% where the ID and order of sparse neurons in Hahnloser plots are located
fd_hahnloser = fullfile(fd_base, 'Figures', pairID, 'HahnloserNew', birdID, 'popRaster2');
suffix_hahn = 'neuron_orderedPlotted7';
suffix_criteria = 'criteria7'; 
% where to save results and plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_plots2'], sprintf('density_%s_%s', syls{1}, syls{2}));
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
disp(fd_save);


%% 1. Load UMAP data, ephys spike location data, and selected sparse neuronID and order
u_all = cell(size(syls, 2), 1);
e_all = cell(size(syls, 2), 1);
n_all = cell(size(syls, 2), 1);
c_all = cell(size(syls, 2), 1);
% syl_i = 1;
for syl_i=1:size(syls, 2)
  ss = syls{syl_i};
  fn_u = fullfile(fd_embed, sprintf('%s.%s.%s.replaced2.mat', birdID, ss, umap_run));
  load(fn_u); u_all{syl_i} = umap;
  fn_e = fullfile(fd_embed, sprintf('%s.%s.sliding_loc.replaced2.mat', birdID, ss));
  load(fn_e); e_all{syl_i} = spike_embed;
  fn_n = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_hahn));
  load(fn_n); n_all{syl_i} = neuron_ordered;
  fn_c = fullfile(fd_hahnloser, ss, sprintf('Hahnloser-%s-chan0.%s.mat', ss, suffix_criteria));
  load(fn_c); c_all{syl_i} = criteria;
end


%% 2. Create a background image to represent the trajectory divergence
% avoid plot the same rendition more than once
% plot renditions that happen in all neurons
% neuron_union = union(n_all{:});
neuron_all = cellfun(@(x) {x.neuronID}, e_all, 'UniformOutput', false);
neuron_all = unique(cat(2, neuron_all{:}), 'sorted');

fig_size = [10 10 600 600];
traj_color = {'#737373', '#737373'};
% use the same x/y limits for all plots
U = horzcat(u_all{:});   % concatenate all struct arrays
uvs = cat(1, U.umap);
x_lim = [min(uvs(:,1))-1, max(uvs(:,1))+1];
y_lim = [min(uvs(:,2))-1, max(uvs(:,2))+1];

close all;
fig = ZZfunc_newFigurePDFsize_v1(fig_size);
ax = gca; hold(ax, 'on');
already_plot = {};
for ni=1:size(neuron_all, 2)
  for syl_i=1:size(syls,2)
    spike_embed = e_all{syl_i};
    umap = u_all{syl_i};
    seg_idx = find(strcmp({spike_embed.neuronID}, neuron_all{ni}));
    % loop through renditions, first plot trajectories
    for sii=1:length(seg_idx)
      si = seg_idx(sii);
      if ~ismember(spike_embed(si).syl_ID, already_plot)
        already_plot = [already_plot; spike_embed(si).syl_ID];
        u = umap(si).umap;
        scatter(ax, u(:,1), u(:,2), 10, 'filled', 'MarkerFaceColor', traj_color{syl_i}, 'MarkerFaceAlpha', 0.01, 'MarkerEdgeColor', 'none');
      end
    end
  end
end
xlim(ax, x_lim);
ylim(ax, y_lim);
frame = getframe(ax); img = frame.cdata;
% make a faint image
alpha = 0.7;                                % 0 = original, 1 = white
white_img = 255 * ones(size(img), 'double'); % all white background
faint_img = (1-alpha) * double(img) + alpha * white_img;
faint_img = uint8(faint_img);
close;
% save the img for later use
bg.img=img; bg.faint_img=faint_img; bg.x_lim=x_lim; bg.y_lim=y_lim; bg.fig_size=[10 10 600 600];
fn_bg = fullfile(fd_save, sprintf('%s.%s_%s.bg.mat', birdID, syls{1}, syls{2}));
save(fn_bg, 'bg');


%% 4. load pre-calculated density
fn_dens = fullfile(fd_save, sprintf('%s.dens_info.mat', birdID));
load(fn_dens); 
dens_info = struct2table(dens_info);
% not showing unconfident neurons
c1 = struct2table(c_all{1}); 
c2 = struct2table(c_all{2}); 
cjoin = outerjoin(c1, c2, 'Keys','neuronID', 'MergeKeys',true);
cjoin = innerjoin(dens_info, cjoin, 'Keys','neuronID');

% identify neurons that have double peaks (if 2nd highest peak is at certain percetage of the highest peak)
per_thre = 0.75; 
for ni=1:size(cjoin, 1)
  v1 = sort(cjoin.peakv1{ni}, 'descend');
  if (length(v1)>1) && (v1(2)>=v1(1)*per_thre); cjoin.has_double1(ni) = 1; else; cjoin.has_double1(ni) = 0; end
  v2 = sort(cjoin.peakv2{ni}, 'descend');
  if (length(v2)>1) && (v2(2)>=v2(1)*per_thre); cjoin.has_double2(ni) = 1; else; cjoin.has_double2(ni) = 0; end
end    
    

% filter neurons to get confident ones
min_rend = 10; % min number of renditions
d = cjoin((cjoin.num_rends_c1>=min_rend) & (cjoin.num_rends_c2>=min_rend), :);


%% 5. Plot density in the same figure
% choose 6 neurons that are shared; then plot 2 specific neurons
% first for call v1
neuron_sel1 = {'20241002-ch7', '20240911-ch12', '20240919-ch8', '20240917-ch3', '20240912-ch11', '20240917-ch13', '20240912-ch10', '20240915-ch11', '20240911-ch3', '20240912-ch1'};
idx1 = cellfun(@(x) find(strcmp(d.neuronID, x)), neuron_sel1);
density_all1 = d.denmap1(idx1);
clim_array1 = cellfun(@(x) [0 quantile(x(:), 0.99)], density_all1, 'UniformOutput', false);

% first for call v2
neuron_sel2 = {'20241002-ch7', '20240911-ch12', '20240919-ch8', '20240917-ch3', '20240912-ch11', '20240917-ch13', '20240912-ch10', '20240915-ch11', '20240917-ch14', '20240914-ch12'};
idx2 = cellfun(@(x) find(strcmp(d.neuronID, x)), neuron_sel2);
density_all2 = d.denmap2(idx2);
clim_array2 = cellfun(@(x) [0 quantile(x(:), 0.99)], density_all2, 'UniformOutput', false);

colors = turbo(size(neuron_sel2, 2));
colors1 = colors;
% replace the last two neuron colors 
colors2 = colors;
colors2(end-1,:) = [0.6 0.3 0.64];
colors2(end,:) = [0.3 0.3 0.3];
cmap_array1 = {};
cmap_array2 = {};
for ci=1:size(colors, 1)
  x = colors1(ci,:);
  a = light_to_color_lab(x, 256);
  cmap_array1{ci} = a;
  x = colors2(ci,:);
  a = light_to_color_lab(x, 256);
  cmap_array2{ci} = a;
end


close all;
% add neuron one-by-one to the plot
% for ni=1:size(neuron_sel, 2)
% one plot for call v1, another for call v2
[fig, axs] =  generatePanelGrid_v3(1, 2, [0.8], [], [0.1;0.05], [0.05;0.01], 0.03, [0], [10 50 1000 600], true);
ax1 = axs(1); cla(ax1); hold(ax1, 'on');
ax2 = axs(2); cla(ax2); hold(ax2, 'on');
% plot the background first
xlim(ax1, x_lim); ylim(ax1, y_lim);
hb1 = image(ax1, ax1.XLim, ax1.YLim, flipud(bg.faint_img));
set(hb1, 'AlphaData',0.75);
xlim(ax2, x_lim); ylim(ax2, y_lim);
hb2 = image(ax2, ax2.XLim, ax2.YLim, flipud(bg.faint_img));
set(hb2, 'AlphaData',0.75);

% then overlay the calculate density map
% overlay_density_images(density_all, cmap_array, 0.25, clim_array);
[RGB_out1, A_out1] = blend_density_layers(density_all1, cmap_array1, 0.08, clim_array1);
hd1 = image(ax1, ax1.XLim, ax1.YLim, RGB_out1, 'AlphaData', A_out1);
xlim(ax1, x_lim); ylim(ax1, y_lim);
[RGB_out2, A_out2] = blend_density_layers(density_all2, cmap_array2, 0.08, clim_array2);
hd2 = image(ax2, ax2.XLim, ax2.YLim, RGB_out2, 'AlphaData', A_out2);
xlim(ax2, x_lim); ylim(ax2, y_lim);

  
fn_fig = fullfile(fd_save, sprintf('%s.firingFields2.v1_v2.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');







