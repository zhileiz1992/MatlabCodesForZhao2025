% using MMD to quantify the overlap and divergence of call trajectories
% 08/06/2025

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/plotTraj'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% what syllable subtypes this bird has
syls = {'v4', 'v5'};
% syls = {'v1', 'v7'};
% syls = {'v5', 'v7'};
% syls = {'v1', 'v4'};
% what's the window size
win_frame = 32;
% where is VAE/UMAP embedding located
vae_run = 'traj_chop_32_1_32';
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'applySylAll', sprintf('latents.%s', vae_run));
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
info = readtable(fn_info, 'Delimiter', ',');
% where to save
fd_save =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD');
if ~exist(fd_save, 'dir')
  mkdir(fd_save);
end
fprintf('Save MMD results to %s\n', fd_save);


%% 1. Retrieve VAE data
% use data from batch1, since batch 2 may have calls from partner bird
batch = {'batch1'};
d_vae = ZZfunc_retrieveVAE_v1(fn_vae, info, batch, syls);
  

%% 2. Align the VAE time-series (32-dim trajectories): interpolate 
% what's the mean length of VAE 
vae_lens1 = cellfun(@(x) size(x,1), {d_vae{1}.vae});
vae_lens2 = cellfun(@(x) size(x,1), {d_vae{2}.vae});
mean_len = floor(mean([vae_lens1 vae_lens2]));
fprintf('Intepolate to %d frames (ms)\n', mean_len);
% loop through syllable types
for si=1:size(syls,2)
  d_syl = d_vae{si};
  parfor ri=1:size(d_syl,2)
    d = d_syl(ri).vae;
    interp_d = ZZfunc_interp_window(d, mean_len);
    d_syl(ri).vae_interp = interp_d;
  end
  d_vae{si} = d_syl;
end
% plot to check
figure; 
subplot(1, 2, 1);
imagesc(d_vae{2}(105).vae); title('Original VAE');
subplot(1, 2, 2);
imagesc(d_vae{2}(105).vae_interp); title('Interpolated VAE');


%% 3.1 Calculate MMD along time series: use all samples
% also calculate between self halves
% mmd_all = [];
% d1 = cat(3, d_vae{1}.vae_interp);
% d2 = cat(3, d_vae{2}.vae_interp);
% rng(1992);
% parfor ti=1:mean_len
%   X = squeeze(d1(ti, :, :));
%   X = X';
%   Y = squeeze(d2(ti, :, :));
%   Y = Y';
%   mmd_all(ti).XY = ZZfunc_mmd_gaussian(X, Y, -1);  % estimate sigma using the median approach
%   
%   idx = randperm(size(X, 1));
%   half = floor(size(X, 1)/2);
%   mmd_all(ti).XX = ZZfunc_mmd_gaussian(X(idx(1:half),:), X(idx((half+1):end),:), -1);  % estimate sigma using the median approach
%   
%   idx = randperm(size(Y, 1));
%   half = floor(size(Y, 1)/2);
%   mmd_all(ti).YY = ZZfunc_mmd_gaussian(Y(idx(1:half),:), Y(idx((half+1):end),:), -1);  % estimate sigma using the median approach 
% end
% mmd_all = struct2table(mmd_all);


%% 3.2 Calculate MMD along time series: partition data, use same sample size
% also calculate between self halves
d1 = cat(3, d_vae{1}.vae_interp);
d2 = cat(3, d_vae{2}.vae_interp);
% partition the data into two halves
part_size = floor(min([size(d1,3)/2, size(d2,3)/2]));
r_seed = 1992;
X_idx = sample_sublists(1:size(d1,3), 2, part_size, r_seed);
X1 = d1(:, :, X_idx(1,:));
X2 = d1(:, :, X_idx(2,:));
Y_idx = sample_sublists(1:size(d2,3), 2, part_size, r_seed);
Y1 = d2(:, :, Y_idx(1,:));
Y2 = d2(:, :, Y_idx(2,:));
% fixed sigma to be the median
sigma_all = [];
parfor ti=1:mean_len
  tempX = squeeze(d1(ti,:,:))';
  tempY = squeeze(d2(ti,:,:))';
  all_data = [tempX; tempY];
  D = pdist2(all_data, all_data);
  sigma = median(D(:));
  sigma_all(ti) = sigma;
end
% add a scale factor: pearson used 0.25 for their Figure 7
factor = 0.25;
sigma_width = mean(sigma_all) * factor;

% for sigma_width=0.3:0.1:0.8
% loop through time points
mmd_all = [];
sigma_width = 0.5; % or used a fixed sigma for all
fprintf('Sigma used %.5f\n', sigma_width);
parfor ti=1:mean_len
  x1 = squeeze(X1(ti, :, :))';
  y1 = squeeze(Y1(ti, :, :))';
  mmd_all(ti).XY = ZZfunc_mmd_gaussian(x1, y1, sigma_width);  % estimate sigma using the median approach
  
  x2 = squeeze(X2(ti, :, :))';
  mmd_all(ti).XX = ZZfunc_mmd_gaussian(x1, x2, sigma_width);  % estimate sigma using the median approach
  
  y2 = squeeze(Y2(ti, :, :))';
  mmd_all(ti).YY = ZZfunc_mmd_gaussian(y1, y2, sigma_width);  % estimate sigma using the median approach 
end
mmd_all = struct2table(mmd_all);
% save results
fn_mmd = fullfile(fd_save, sprintf('%s.%s.%s.mmd_all.mat', birdID, syls{1}, syls{2}));
save(fn_mmd, 'mmd_all');


%% 4. Plot MMD
% close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 600 400]);
ax = gca;
rel_t = (1:mean_len) - win_frame;
hold(ax, 'on');
plot(ax, rel_t, mmd_all.XY, '-r', 'LineWidth', 3, 'DisplayName', sprintf('%s-%s', syls{1}, syls{2})); 
plot(ax, rel_t, mmd_all.XX, '-k', 'LineWidth', 3, 'DisplayName', sprintf('%s-%s', syls{1}, syls{1}));
plot(ax, rel_t, mmd_all.YY, '-b', 'LineWidth', 3, 'DisplayName', sprintf('%s-%s', syls{2}, syls{2}));
% add line for syllable start and end
xline(ax, 0, 'LineStyle', '--', 'HandleVisibility', 'off');
% xline(ax, rel_t(end)-win_frame, 'LineStyle', '--', 'HandleVisibility', 'off');
xlim(ax, [rel_t(1)-10 rel_t(end)+10]);
legend('show', 'Location', 'best', 'FontSize', 12);
xlabel('Relative time (ms)', 'FontSize', 16);
ylabel('Maximum mean discrepancy', 'FontSize', 16);
title(ax, sprintf('Sigma=%.3f', sigma_width), 'FontSize', 14);
% save figure
fn_pdf = fullfile(fd_save, sprintf('%s.%s.%s.mmd.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 5. Does the variance change over time?
var_all = zeros(mean_len, 1);
d_this = d2;
for ti=1:mean_len
  temp = squeeze(d_this(ti, :, :))';
%   v = var(temp, 1);
  cov_m = cov(temp);
  total_v = trace(cov_m);
  var_all(ti) = total_v;
end
figure; 
plot(rel_t, var_all, 'LineWidth', 3);
xlabel('Relative time (ms)', 'FontSize', 16);
ylabel('Total variance', 'FontSize', 16);



