% plot the results from MMD matrix analysis for calls
% reorder based on the latest renaming of call subtypes
% make color consistent

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% what syllable subtypes this bird has
oldID = {'v4', 'v5', 'v1', 'v7', 'v3', 'v6', 'v2'};
newID = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
col_list = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#ccad25'};
% what's the window size
win_frame = 32;
% where is VAE/UMAP embedding located
vae_run = 'traj_chop_32_1_32';
% where is retrieved VAE data stored
fd_home = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed');
fd_vaeRes = fullfile(fd_home, 'sylOnly', 'batch1');
% what .mat to use: time-warped
suffix = 'dvae_warp';
% where to save
fd_save =  fullfile(fd_home, 'MMD', 'matrixCall', 'all');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save MMD results to %s\n', fd_save);
% what sigma to try
sigma_list = [0.5 0.7 0.9 1.1 1.3];
% sigma_list = [1.1 1.3];
% sigma = 0.9;


%% 1. Read data, consider choping some ends off
% load the results and concat into 3d arrays
% how much to chop off at the ends
chopoff = 31;
syls = oldID;
m = cell(size(syls,2), size(syls,2), length(sigma_list));
b = m;  c = m; c_mask=m;
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    for si=1:length(sigma_list)
      fn_res = fullfile(fd_save, 'intermediate', sprintf('%s.%s_%s.%.2f.mmd.mat', birdID, syls{ii}, syls{jj}, sigma_list(si)));
      if exist(fn_res, 'file')
        a = load(fn_res).mmd;
        mmd=a.mmd; B=a.B; C=a.C;
      else
        fn_res = fullfile(fd_save, 'intermediate', sprintf('%s.%s_%s.%.2f.mmd.mat', birdID, syls{jj}, syls{ii}, sigma_list(si)));
        a = load(fn_res).mmd;
        mmd=a.mmd'; B=a.B'; C=a.C';
      end
      % retrieve values, chop if specified
      m{ii,jj,si} = mmd((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      b{ii,jj,si} = B((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      c{ii,jj,si} = C((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      c_mask{ii,jj,si} = C((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      % fill in the other half
      if jj>ii
        m{jj,ii,si} = m{ii,jj,si}'; b{jj,ii,si} = b{ii,jj,si}'; c{jj,ii,si} = c{ii,jj,si}';
        c_mask{jj, ii, si} = NaN(size(c{jj,ii,si}));
      end
    end
  end
end
% concatenate
lens = cellfun(@(x) size(x,1), m(:,1,1));
mc = cell2tiled3d(m); bc = cell2tiled3d(b); cc=cell2tiled3d(c); cc_mask=cell2tiled3d(c_mask);


%% 2. Plot the full matrix
close all;
clims = [-0.05 -0.15 -0.2 -0.25 -0.3];
% where to put the call subtype labels
v_loc = cumsum(lens)-lens/2;
line_loc = cumsum(lens);
for si=1:size(c,3)
  [fig, axes] = generatePanelGrid_v2(1, 3, [0.8], [], [0.125;0.05], [0.05;0.05], 0.05, [0], [10 10 1700 500]);
  ZZfunc_showMMDmatrix_v1(axes(1), squeeze(mc(:,:,si)), [], v_loc, newID, 'Original MMD', true, line_loc);
  ZZfunc_showMMDmatrix_v1(axes(2), squeeze(bc(:,:,si)), [], v_loc, newID, 'Rank-1 MMD', true, line_loc);
  ZZfunc_showMMDmatrix_v1(axes(3), squeeze(cc(:,:,si)), [clims(si),0], v_loc, newID, 'delta MMD', true, line_loc);
  sgtitle(sprintf('Sigma=%.2f', sigma_list(si)), 'FontSize', 16);
  fn_fig = fullfile(fd_save, sprintf('%s.vAllRename.%.2f.MMD.pdf', birdID, sigma_list(si)));
  print(fig, fn_fig, '-dpdf', '-painters');
end
  

%% 3. plot only the upper half
clims1 = [-0.05 -0.15 -0.2 -0.25 -0.3];
clims2 = [0 -0.005 -0.02 -0.02 -0.02];
for si=1:length(sigma_list)
cc_this = squeeze(cc_mask(:, :, si));
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]); ax=gca;
ZZfunc_showMMDmatrix_v2(ax, cc_this, [clims1(si),clims2(si)], v_loc, newID, sprintf('dMMD sigma=%.2f', sigma_list(si)), true, line_loc);
fn_fig = fullfile(fd_save, sprintf('%s.vAllRename.%.2f.MMDhalf.pdf', birdID, sigma_list(si)));
print(fig, fn_fig, '-dpdf', '-painters');
end


%% 4. then plot a zoom version of target syllables
% v_tar = {'v1', 'v2'};
v_tar = {'v3', 'v4'};
si = 3;
v_tar_old = cellfun(@(x) oldID{strcmp(newID, x)}, v_tar, 'UniformOutput', false);
fn_res = fullfile(fd_save, 'intermediate', sprintf('%s.%s_%s.%.2f.mmd.mat', birdID, v_tar_old{1}, v_tar_old{2}, sigma_list(si)));
if exist(fn_res, 'file')
  a = load(fn_res).mmd;
  c = a.C((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
else
  fn_res = fullfile(fd_save, 'intermediate', sprintf('%s.%s_%s.%.2f.mmd.mat', birdID, v_tar_old{2}, v_tar_old{1}, sigma_list(si)));
  a = load(fn_res).mmd; temp=a.C';
  c = temp((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
end
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]); ax=gca;
imagesc(ax, c, [clims1(si) clims2(si)]);
colormap(ax, 'gray'); colorbar(ax);
axis(ax, 'equal');
axis(ax, 'tight');
title(ax, sprintf('%s vs %s', v_tar{1}, v_tar{2}), 'FontSize', 20);
fn_fig = fullfile(fd_save, sprintf('%s.vAllRename.%.2f.MMDzoom.%s_%s.pdf', birdID, sigma_list(si), v_tar{1}, v_tar{2}));
print(fig, fn_fig, '-dpdf', '-painters');
















