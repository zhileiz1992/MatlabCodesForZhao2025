% loop all pairs of calls, calculate and plot the MMD matrix
% differ from v1: smaller range of sigma, save shuffled data as well; use the replaced names for call subtypes
clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% what syllable subtypes this bird has
syls = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
% what's the window size
win_frame = 32;
% where is VAE/UMAP embedding located
vae_run = 'traj_chop_32_1_32';
% where is retrieved VAE data stored
fd_home = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed');
fd_vaeRes = fullfile(fd_home, 'sylOnly2', 'batch1');
% what .mat to use: time-warped
suffix = 'dvae_warp';
% where to save
fd_save =  fullfile(fd_home, 'MMD', 'matrixCall2', 'all');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save MMD results to %s\n', fd_save);
% what sigma to try
% sigma_list = [0.5 0.7 0.9];
% sigma_list = [1.1 1.3];
sigma_list = [0.9];


%% 1. Retrieve VAE data
d_all = cell(size(syls,2), 1);
for vi=1:size(syls,2)
  fn = fullfile(fd_vaeRes, sprintf('%s.%s.%s.%s.mat', birdID, vae_run, syls{vi}, suffix));
  a = load(fn);
  d_all{vi} = a.dvae;
end
% what's the size of each calls
lens = cellfun(@(x) size(x,2), d_all);
fprintf('%d ', lens);


%% 2. Calculate the MMD matrix
fd_save_this = fullfile(fd_save, 'intermediate');
if ~exist(fd_save_this, 'dir')
  mkdir(fd_save_this);
end
% loop through pairs, sample equal amount
to_sample = 1500; 
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    dvae1 = d_all{ii};
    dvae2 = d_all{jj};
    
    % sample equal number of syllables
    rng(1992);
    % sample data then convert to 3d arrays
    idx1 = randsample(1:size(dvae1,2), to_sample);
    d1 = cat(3, dvae1(idx1).vae_interp);
    idx2 = randsample(1:size(dvae2,2), to_sample);
    d2 = cat(3, dvae2(idx2).vae_interp);
    
    % split the self population to get null
    to_sample_half = floor(to_sample/2);
    self1 = sample_subarrays(1:to_sample, [to_sample_half, to_sample_half]);
    self2 = sample_subarrays(1:to_sample, [to_sample_half, to_sample_half]);
    
    % calculate MMD matrix, vary sigma
    for sigma_i=1:length(sigma_list)
      sigma = sigma_list(sigma_i);
      fprintf('Comparing %s vs %s using sigma=%.2f\n', syls{ii}, syls{jj}, sigma);
      standardize = false;
      mmd_matrix = mmd_matrix_gaussian(d1, d2, sigma, standardize);
      % decompose MMD as Pearson did
      [B, C, u, v] = mmd_matrix_decompose(mmd_matrix);
      % save results
      mmd.mmd = mmd_matrix;
      mmd.B = B;
      mmd.C = C;
      mmd.u = u;
      mmd.v = v;
      mmd.idx1 = idx1;
      mmd.idx2 = idx2; 
      
      % also perform on the self distribution
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
      
      fn_res = fullfile(fd_save_this, sprintf('%s.%s_%s.%.2f.mmd.mat', birdID, syls{ii}, syls{jj}, sigma));
      save(fn_res, 'mmd');
    end
  end
 
end


%% 3. Plot results
% load the results and concat into 3d arrays
% how much to chop off at the ends
chopoff = 31;
m = cell(size(syls,2), size(syls,2), length(sigma_list));
b = m;  c = m;
lens = zeros(size(syls,2),1);
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    for si=1:length(sigma_list)
      fn_res = fullfile(fd_save, 'intermediate', sprintf('%s.%s_%s.%.2f.mmd.mat', birdID, syls{ii}, syls{jj}, sigma_list(si)));
      a = load(fn_res).mmd;
      % retrieve values, chop if specified
      m{ii,jj,si} = a.mmd((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      b{ii,jj,si} = a.B((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      c{ii,jj,si} = a.C((1+chopoff):(end-chopoff), (1+chopoff):(end-chopoff));
      % fill in the other half
      if jj>ii
        m{jj,ii,si} = m{ii,jj,si}'; b{jj,ii,si} = b{ii,jj,si}'; c{jj,ii,si} = c{ii,jj,si}';
      end       
    end
  end
end
% concatenate
lens = cellfun(@(x) size(x,1), m(:,1,1));
mc = cell2tiled3d(m); bc = cell2tiled3d(b); cc=cell2tiled3d(c); 


% plots
close all;
clims = [-0.05 -0.15 -0.2];
% where to put the call subtype labels
v_loc = cumsum(lens)-lens/2;
line_loc = cumsum(lens);
for si=1:size(c,3)
  [fig, axes] = generatePanelGrid_v2(1, 3, [0.8], [], [0.125;0.05], [0.05;0.05], 0.05, [0], [10 10 1700 500]);
  ZZfunc_showMMDmatrix_v1(axes(1), squeeze(mc(:,:,si)), [], v_loc, syls, 'Original MMD', true, line_loc);
  ZZfunc_showMMDmatrix_v1(axes(2), squeeze(bc(:,:,si)), [], v_loc, syls, 'Rank-1 MMD', true, line_loc);
  ZZfunc_showMMDmatrix_v1(axes(3), squeeze(cc(:,:,si)), [clims(si),0], v_loc, syls, 'delta MMD', true, line_loc);
  sgtitle(sprintf('Sigma=%.2f', sigma_list(si)), 'FontSize', 16);
  fn_fig = fullfile(fd_save, sprintf('%s.vAll.%.2f.MMD.pdf', birdID, sigma_list(si)));
  print(fig, fn_fig, '-dpdf', '-painters');
end
  























