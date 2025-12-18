% Given pairs of call subtypes, calculate and plot the MMD matrix
% differ from v1: smaller range of sigma, save shuffled data as well

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% what syllable subtypes this bird has
syls = {'v4', 'v5'};
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
fd_save =  fullfile(fd_home, 'MMD', 'matrixCall2', [syls{1} '_' syls{2}]);
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save MMD results to %s\n', fd_save);


%% 1. Retrieve VAE data
d_all = cell(size(syls,2), 1);
for vi=1:size(syls,2)
  fn = fullfile(fd_vaeRes, sprintf('%s.%s.%s.%s.mat', birdID, vae_run, syls{vi}, suffix));
  a = load(fn);
  d_all{vi} = a.dvae;
end


%% 2. Calculate the MMD matrix
% sample equal number of syllables
num_total = cellfun(@(x) size(x, 2), d_all);
to_sample = min(num_total);
rng(1992);
% sample data then convert to 3d arrays
idx1 = randsample(1:size(d_all{1},2), to_sample);
d1 = cat(3, d_all{1}(idx1).vae_interp);
idx2 = randsample(1:size(d_all{2},2), to_sample);
d2 = cat(3, d_all{2}(idx2).vae_interp);
% split the self population to get null
to_sample_half = floor(to_sample/2);
self1 = sample_subarrays(1:to_sample, [to_sample_half, to_sample_half]);
self2 = sample_subarrays(1:to_sample, [to_sample_half, to_sample_half]);
% calculate MMD matrix, vary sigma
for sigma=0.5:0.2:1.1
  fprintf('Calculating using sigma=%.2f\n', sigma);
  standardize = false;
  mmd_matrix = mmd_matrix_gaussian(d1, d2, sigma, standardize);  
  % decompose MMD as Pearson did
  [B, C, ~, ~] = mmd_matrix_decompose(mmd_matrix);
  mmd.mmd = mmd_matrix;
  mmd.B = B;
  mmd.C = C;
  
  % also calculate MMD on the self-population
  mmd_matrix11 = mmd_matrix_gaussian(d1(:,:,self1{1}), d1(:,:,self1{2}), sigma, standardize); 
  [B, C, ~, ~] = mmd_matrix_decompose(mmd_matrix11);
  mmd.mmd11 = mmd_matrix11;
  mmd.B11 = B; 
  mmd.C11 = C;
  
  mmd_matrix22 = mmd_matrix_gaussian(d2(:,:,self2{1}), d2(:,:,self2{2}), sigma, standardize); 
  [B, C, ~, ~] = mmd_matrix_decompose(mmd_matrix22);
  mmd.mmd22 = mmd_matrix22;
  mmd.B22 = B; 
  mmd.C22 = C;
  
  % save results

  fn_res = fullfile(fd_save, sprintf('%s.%s.%.2f.mmd.mat', birdID, strjoin(syls,''), sigma));
  save(fn_res, 'mmd');
  
  
  % plot
  close all; fig = ZZfunc_newFigurePDFsize_v1([10 10 1100 900]);
  subplot(2,2,1);
  imagesc(mmd_matrix); colormap gray; colorbar;
  title('Original MMD', 'FontSize', 12);
  subplot(2,2,2);
  imagesc(C); colormap gray; colorbar;
  title('delta MMD', 'FontSize', 12);
  % plot; change color range
  subplot(2,2,3);
  imagesc(mmd_matrix, [0 0.1]); colormap gray; colorbar;
  title('Original MMD', 'FontSize', 12);
  subplot(2,2,4);
  imagesc(C, [-0.05 0]); colormap gray; colorbar;
  title('delta MMD', 'FontSize', 12);
  % save
  fn_fig = fullfile(fd_save, sprintf('%s.%s.%.2f.mmd.pdf', birdID, strjoin(syls,''), sigma));
  print(fig, fn_fig, '-dpdf', '-painters');
  
  % check the diagnonal
  diag_v = diag(mmd_matrix);
  fig = ZZfunc_newFigurePDFsize_v1([10 10 600 600]);
  plot(diag_v, 'LineWidth', 3);
  fn_fig = fullfile(fd_save, sprintf('%s.%s.%.2f.diagonal.pdf', birdID, strjoin(syls,''), sigma));
  print(fig, fn_fig, '-dpdf', '-painters');
end











