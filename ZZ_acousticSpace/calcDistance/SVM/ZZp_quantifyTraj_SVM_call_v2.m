% build SVM decoders to quantify the overlap and divergence of call trajectories in VAE space
% Zhilei, 08/07/2025
% differ from v1: run many runs to get confidence intervals

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/plotTraj'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdID = 'pair5RigCCU29';
pairID = 'pair5CU29CU55';
% what syllable subtypes this bird has
% syls = {'v4', 'v5'};
syls = {'v1', 'v7'};
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
fd_save =  fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'SVM');
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
% plot an example to check the alignment
figure; 
subplot(1, 2, 1);
imagesc(d_vae{2}(105).vae); title('Original VAE');
subplot(1, 2, 2);
imagesc(d_vae{2}(105).vae_interp); title('Interpolated VAE');


%% 3.2 SVM: PCA first, then control the amount of input info
Xraw1 = cat(3, d_vae{1}.vae_interp);
Xraw2 = cat(3, d_vae{2}.vae_interp);

%% PCA to reduce dimensions
% Sizes
[n1_samp, n_feat, n1_ind] = size(Xraw1);
[n2_samp, ~,    n2_ind]  = size(Xraw2);
% Reshape: merge samples × individuals into one dimension
X1_reshaped = reshape(permute(Xraw1, [1 3 2]), n1_samp * n1_ind, n_feat);
X2_reshaped = reshape(permute(Xraw2, [1 3 2]), n2_samp * n2_ind, n_feat);
% Concatenate and perform PCA
X_all = [X1_reshaped; X2_reshaped];  % size: (n_total × n_feat)
% standarize
X_all = zscore(X_all);
[coeff, score_all, latent] = pca(X_all);
% Separate projections
score1_flat = score_all(1 : n1_samp * n1_ind, :);
score2_flat = score_all(n1_samp * n1_ind + 1 : end, :);
% Number of PCs
n_pcs = size(score_all, 2);
% Reshape back: [samples × PCs × individuals]
X1 = permute(reshape(score1_flat, n1_samp, n1_ind, n_pcs), [1 3 2]);
X2 = permute(reshape(score2_flat, n2_samp, n2_ind, n_pcs), [1 3 2]);
% plot the variance explained 
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 400 400]);
plot(1:size(latent,1), cumsum(latent)/sum(latent), '-bo', 'Marker', '.', 'MarkerSize', 20, 'LineWidth', 2);
ylim([0 1]);
xlabel('No. of PCs', 'FontSize', 14);
ylabel('Variance explained', 'FontSize', 14);
fn_pca = fullfile(fd_save, sprintf('%s.%s.%s.PCA.variance.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pca, '-dpdf', '-painters');


%%  partition into train/test, perform many shuffle runs
Y1 = repmat({syls{1}}, size(X1,3), 1);
Y2 = repmat({syls{2}}, size(X2,3), 1);
% determine the train/test split
train_ratio = 0.5;
num_total = min([size(X1,3), size(X2,3)]);
num_train = floor(train_ratio*num_total);
num_test = num_total - num_train;

num_runs = 50;
pc_use = 1:10;
accu_all = zeros(num_runs, length(pc_use), mean_len, 2);
rng(1992); % set a master seed
% loop through runs
for run_i=1:num_runs
  fprintf('Perform run %d\n', run_i);
  % split data into train and test
  i_X1 = sample_subarrays(1:size(X1,3), [num_train, num_test]);
  i_X2 = sample_subarrays(1:size(X2,3), [num_train, num_test]);
  X_train = cat(3, X1(:,:,i_X1{1}), X2(:,:,i_X2{1}));
  y_train = [Y1(i_X1{1}); Y2(i_X2{1})];
  y_shuffle = y_train(randperm(numel(y_train)));
  X_test = cat(3, X1(:,:,i_X1{2}), X2(:,:,i_X2{2}));
  y_test = [Y1(i_X1{2}); Y2(i_X2{2})];
  
  % loop through number of PCs to use
  for pi=1:length(pc_use)
    pu = pc_use(pi);
    % loop through time points
    accu_this = zeros(mean_len,2);
    parfor ti=1:mean_len
      x_train = X_train(ti, 1:pu, :);
      x_train = reshape(x_train, size(x_train,2), size(x_train,3))';
      x_test = X_test(ti, 1:pu, :);
      x_test = reshape(x_test, size(x_test,2), size(x_test,3))';
      [~, accu_real] = trainSVM(x_train, y_train, x_test, y_test, 'linear');
      % also train on shuffled labels
      [~, accu_shuffle] = trainSVM(x_train, y_shuffle, x_test, y_test, 'linear');
      accu_this(ti,:) = [accu_real accu_shuffle];
    end
    accu_all(run_i, pi, :, :) = accu_this;
  end
end
% save results
fn_accuracy = fullfile(fd_save, sprintf('%s.%s.%s.svm.mat', birdID, syls{1}, syls{2}));
save(fn_accuracy, 'accu_all');


%% 4. Plot SVM accuracy
close all;
t = (1:mean_len) - win_frame;
% remove the baseline on accuracy
baseline = 0.5;
% two colors for real and shuffled
col_list = [0.98 0.55 0.38; 0.5 0.5 0.5];
% fig = ZZfunc_newFigurePDFsize_v1([10 10 1800 800]);
[fig, axes] = generatePanelGrid_v2(2, 5, [0.35;0.35], [0.1], [0.05;0.05], [0.05;0.05], 0.035, [1;1], [10 10 1800 800]);
% loop through number of PCs
for pi=1:size(accu_all,2)
  % what plot panel
  plot_i = idivide((pi-1), int32(5))+1;
  plot_j = mod((pi-1), 5)+1;
  ax = axes(plot_i, plot_j);
  d = squeeze(accu_all(:, pi, :, :));
  % plot mean accuracy with 95% confidence interval
  [n_runs, n_time, n_cond] = size(d);
  hold(ax, 'on');
  for cond = 1:2
    data = squeeze(d(:, :, cond));  % n_runs × n_time
    % Compute mean and 95% confidence interval
    mu = mean(data, 1);  % 1 × n_time
    sem = std(data, 0, 1) / sqrt(n_runs);  % standard error of the mean
    ci95 = 1.96 * sem;
    % Shaded area
    fill(ax, [t, fliplr(t)], [mu + ci95, fliplr(mu - ci95)], col_list(cond,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % Mean line
    plot(ax, t, mu, 'Color', col_list(cond,:), 'LineWidth', 2);
  end
  
  % add line for syllable start
  xline(ax, 0, 'LineStyle', '--', 'HandleVisibility', 'off');
  yline(ax, baseline, 'LineStyle', '--', 'HandleVisibility', 'off');
  % xline(ax, rel_t(end)-win_frame, 'LineStyle', '--', 'HandleVisibility', 'off');
  xlim(ax, [t(1)-10 t(end)+10]);
  ylim(ax, [0.25 1]);
  % legend('show', 'Location', 'best', 'FontSize', 12);
  xlabel(ax, 'Relative time (ms)', 'FontSize', 12);
  ylabel(ax, 'SVM accuracy', 'FontSize', 12);
  title(ax, sprintf('PC used: %d', pc_use(pi)),  'FontSize', 12);
end

% save figure
fn_pdf = fullfile(fd_save, sprintf('%s.%s.%s.SVMonPCA.pdf', birdID, syls{1}, syls{2}));
print(fig, fn_pdf, '-dpdf', '-painters');










