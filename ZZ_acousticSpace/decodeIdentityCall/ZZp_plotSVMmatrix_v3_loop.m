% plot the calculated SVM decoding accuracy matrix for all pairwise comparison between call subtypes
% Zhilei, 09/05/2025
% differ from v2: plot acoustic MMD and neural SVM side-by-side


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
col_list = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#ccad25'};
for si=1:size(syls,2); col_dict.(syls{si})=col_list{si}; end
% where is calculated decoding matrix located
vae_run = 'traj_chop_32_1_32';
fd_data = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix4'], 'intermediate');
% what's the window size and hop size when calculating SVM decoding matrix, unit is seconds
fs = 20000;
pad = 0.08;
bin = 0.01;
hop = 0.002;
% where is the MMD result
fd_mmd = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD', 'MatrixCall2', 'all', 'intermediate');
sigma = 0.9;
win_frame = 32;
hop_frame = 1;
sec_per_frame = 0.001;
% where to save the plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix4']);
fprintf('Save SVM plots to %s\n', fd_save);


%% 1.1 Read neural SVM matrix
neu = cell(size(syls,2), size(syls,2));  % mean accuracy matrix of the real dataset
% also trim to desired time range
rel_t_neu = (0:1000)*hop - pad;  % a long time cooridiante axis for the sliding window start
% rel_t_neu = rel_t_neu + bin/2;  % convert to window center
% tpad_n = [0.032 0];  % how much time before syllable onset and after syllable offset
tpad_n = [0 -0.032]; 
neu2 = cell(size(syls,2), size(syls,2));  % trimmed matrix
neu_shuff = cell(size(syls,2), size(syls,2)); % shuffled data
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    run_name = sprintf('%s_%s', syls{ii}, syls{jj});
    fn_res = fullfile(fd_data, sprintf('%s.%s.accu_allMatrix.mat', birdID, run_name));
    a = load(fn_res); temp = a.accu_all;
    n_this = mean(squeeze(temp(:,:,1,:)), 3);
    neu{ii, jj} = n_this;
    % trim to desired time range
    ty = rel_t_neu(1:size(n_this,1));
    tx = rel_t_neu(1:size(n_this,2));
    tendy = ty(end) - pad + tpad_n(2);
    iy = find((ty>=-tpad_n(1)) & (ty<=tendy));
    tendx = tx(end) - pad + tpad_n(2);
    ix = find((tx>=-tpad_n(1)) & (tx<=tendx));
    neu2{ii,jj} = n_this(iy, ix);
    neu_shuff{ii,jj} = mean(squeeze(temp(:,:,2,:)), 3);
  end
end

% tile into a large 2d matrix
lens_n = cellfun(@(x) size(x,2), neu2(1,:));
neu3 = ZZfunc_tileArrayCell_v1(neu2);
% save for later use
% relative time cooridinates after setting the time range
rel_t_n_chop = rel_t_neu(rel_t_neu>=-tpad_n(1));
d_neu.neu=neu; d_neu.neu2=neu2; d_neu.neu3=neu3; d_neu.lens_n=lens_n; d_neu_shuff=neu_shuff; 
d_neu.rel_t_neu=rel_t_neu; d_neu.rel_t_n_chop=rel_t_n_chop;  
fn_neu = fullfile(fd_save, sprintf('%s.d_neu.mat', birdID));
save(fn_neu, 'd_neu');



%% 1.2 Read acoustic MMD matrix in
acoRaw = cell(size(syls,2), size(syls,2));
aco = cell(size(syls,2), size(syls,2));  % mean accuracy matrix of the real dataset
% also trim to desired time range
rel_t_aco = (0:1000)*hop_frame*sec_per_frame - win_frame*sec_per_frame;  % a long time cooridiante axis for the sliding window start
% rel_t_aco = rel_t_aco + win_frame*sec_per_frame/2;  % convert to window center
% tpad_a = [0.032 0];  % how much time before syllable onset and after syllable offset
tpad_a = [0 -0.032];
aco2 = cell(size(syls,2), size(syls,2));  % time-trimmed matrix
aco_shuff = cell(size(syls,2), size(syls,2)); % shuffled data
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    run_name = sprintf('%s_%s', syls{ii}, syls{jj});
    fn_res = fullfile(fd_mmd, sprintf('%s.%s.%.2f.mmd.mat', birdID, run_name, sigma));
    a = load(fn_res); a_this = a.mmd.C;
    aco{ii, jj} = a_this;
    acoRaw{ii,jj} = a.mmd.mmd;
    % trim to desired time range
    ty = rel_t_aco(1:size(a_this,1));
    tx = rel_t_aco(1:size(a_this,2));
    tendy = ty(end) + tpad_a(2);  % note that for MMD, no more extra slide at the end, sliding window stops at the syllable offset
    iy = find((ty>=-tpad_a(1)) & (ty<=tendy));
    tendx = tx(end) + tpad_a(2);
    ix = find((tx>=-tpad_a(1)) & (tx<=tendx));
    aco2{ii,jj} = a_this(iy, ix);
    aco_shuff{ii,jj} = a.mmd.mmd11; 
  end
end
% tile into a large 2d matrix
lens_a = cellfun(@(x) size(x,2), aco2(1,:));
aco3 = ZZfunc_tileArrayCell_v1(aco2);
% save for later use
rel_t_a_chop = rel_t_aco(rel_t_aco>=-tpad_a(1));
d_aco.aco=aco; d_aco.aco2=aco2; d_aco.aco3=aco3; d_aco.acoRaw=acoRaw; d_aco.lens_a=lens_a; d_aco.aco_shuff=aco_shuff;
d_aco.rel_t_aco=rel_t_aco; d_aco.rel_t_a_chop=rel_t_a_chop;
fn_neu = fullfile(fd_save, sprintf('%s.d_aco.mat', birdID));
save(fn_neu, 'd_aco');


%% 2. Plot the full matrix
% neural SVM
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_n)-lens_n/2;
line_loc = cumsum(lens_n);
clim_n = [0.5 0.9];
% ZZfunc_showMMDmatrix_v3(ax, neu3, [0.5 0.9], v_loc, syls, 'Neural SVM accuracy', true, line_loc);
ZZfunc_showMMDmatrix_v3(ax, neu3, clim_n, v_loc, syls, 'Neural SVM accuracy', true, line_loc);
fn_pdf = fullfile(fd_save, sprintf('%s.neuralSVMmatrix.t%.3f_%.3f.c%.2f_%.2f.pdf', birdID, tpad_n(1), tpad_n(2), clim_n(1), clim_n(2)));
print(fig, fn_pdf, '-dpdf', '-painters');

% acoustic MMD
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_a)-lens_a/2;
line_loc = cumsum(lens_a);
% clim_a = [-0.2 -0.035];
clim_a = [-0.2 -0.03];
% ZZfunc_showMMDmatrix_v3(ax, aco3, [-0.2 -0.02], v_loc, syls, 'Acoustic MMD', true, line_loc);
ZZfunc_showMMDmatrix_v3(ax, aco3, clim_a, v_loc, syls, 'Acoustic MMD', true, line_loc);
fn_pdf = fullfile(fd_save, sprintf('%s.acousticMMDmatrix.t%.3f_%.3f.c%.2f_%.2f.pdf', birdID, tpad_n(1), tpad_n(2), clim_a(1), clim_a(2)));
print(fig, fn_pdf, '-dpdf', '-painters');



%% 3. Plot pairwise: acoustic and neural side-by-side;
% p_all = {{'v1', 'v2'}; {'v3', 'v4'}};
fd_save_pair = fullfile(fd_save, 'pairwise_plots');
if ~exist(fd_save_pair, 'dir'); mkdir(fd_save_pair); end
% loop through all pairs
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    
    close all;
    [fig, axes] = generatePanelGrid_v2(1, 2, [0.75], [], [0.15;0.05], [0.1;0.05], 0.1, [0], [10 10 1200 600]);
    % acoustic first
    ax = axes(1);
    a_this = aco2{ii, jj};
    ty = rel_t_a_chop(1:size(a_this,1));
    tx = rel_t_a_chop(1:size(a_this,2));
    imagesc(ax, tx, ty, a_this, [-0.2 -0.02]);
    %   imagesc(ax, tx, ty, a_this, [0 1]);
    colormap(ax, 'gray'); colorbar(ax);
    set(ax, 'YDir', 'reverse');
    axis(ax, 'equal');
    axis(ax, 'tight');
    ax.XAxisLocation = 'top';
    xlabel(ax, sprintf('%s time (sec)', syls{jj}), 'FontSize', 20, 'Color', col_dict.(syls{jj}));
    ylabel(ax, sprintf('%s time (sec)', syls{ii}), 'FontSize', 20, 'Color', col_dict.(syls{ii}));
    title(ax, 'Acoustic similarity (dMMD)', 'FontSize', 20);
    
    % then neural
    ax = axes(2);
    n_this = neu2{ii, jj};
    ty = rel_t_n_chop(1:size(n_this,1));
    tx = rel_t_n_chop(1:size(n_this,2));
    imagesc(ax, tx, ty, n_this, [0.5 0.9]);
    colormap(ax, 'gray'); colorbar(ax);
    set(ax, 'YDir', 'reverse');
    axis(ax, 'equal');
    axis(ax, 'tight');
    ax.XAxisLocation = 'top';
    xlabel(ax, sprintf('%s time (sec)', syls{jj}), 'FontSize', 20, 'Color', col_dict.(syls{jj}));
    ylabel(ax, sprintf('%s time (sec)', syls{ii}), 'FontSize', 20, 'Color', col_dict.(syls{ii}));
    title(ax, 'Neural similarity (SVM)', 'FontSize', 20);
    
    linkaxes(axes, 'xy');
    
    fn_fig = fullfile(fd_save_pair, sprintf('%s.neuralSVM.%s_%s.pdf', birdID, syls{ii}, syls{jj}));
    print(fig, fn_fig, '-dpdf', '-painters');
  end
end



%% 4. Plot the diagonal: all pairwise
fd_save_diag = fullfile(fd_save, 'diagonal_plots');
if ~exist(fd_save_diag, 'dir'); mkdir(fd_save_diag); end
col_lines = {'#4db34d',  '#f15a29'};
for ii=1:size(syls,2)
  for jj=(ii+1):size(syls,2)
    
    a_this = acoRaw{ii, jj};
    diag_mmd = diag(a_this);
    diag_mmd_shuff = diag(aco_shuff{ii, jj});
    t_mmd = rel_t_aco(1:length(diag_mmd));
    t_mmd_shuff = rel_t_aco(1:length(diag_mmd_shuff));
    
    n_this = neu{ii, jj};
    diag_svm = diag(n_this);
    t_svm = rel_t_neu(1:length(diag_svm));
    diag_svm_shuff = diag(neu_shuff{ii,jj});
    t_svm_shuff = rel_t_neu(1:length(diag_svm_shuff));
    
    xlim1 = max([t_mmd(1) t_svm(1)]);
    xlim2 = min([t_mmd(end) t_svm(end)]);
    
    close all;
    [fig, axes] = generatePanelGrid_v2(1, 2, [0.75], [], [0.15;0.05], [0.1;0.05], 0.15, [0], [10 10 1200 600]);
    ax=axes(1); hold(ax, 'on');
    yyaxis(ax, 'left');
    plot(ax, t_svm, diag_svm, 'Color', col_lines{1}, 'LineStyle', '-', 'LineWidth', 2);
    plot(ax, t_svm_shuff, diag_svm_shuff, 'Color', col_lines{1}, 'LineStyle', ':', 'LineWidth', 2);
    ylabel(ax, 'SVM decoding accuracy', 'FontSize', 14);
    ylim(ax, [0.45 max(diag_svm)*1.05]);
    
    yyaxis(ax, 'right');
    plot(ax, t_mmd, diag_mmd, 'Color', col_lines{2}, 'LineStyle', '-', 'LineWidth', 2);
    plot(ax, t_mmd_shuff, diag_mmd_shuff, 'Color', col_lines{2}, 'LineStyle', ':', 'LineWidth', 2);
    ylabel(ax, 'MMD on VAE latents', 'FontSize', 14);
    xlabel(ax, 'Rel. time (sec)', 'FontSize', 14);
    xlim(ax, [xlim1 xlim2]);
    ylim(ax, [0 max(diag_mmd)*1.05]);
    ax.YAxis(1).Color = col_lines{1};
    ax.YAxis(2).Color = col_lines{2};
    title(ax, sprintf('Neural vs acoustic similarity: %s vs %s', syls{ii}, syls{jj}), 'FontSize', 16);
    
    % then plot the cross-correlation
    ax = axes(2); hold(ax, 'on');
    x = diag_svm - mean(diag_svm);
    y = diag_mmd - mean(diag_mmd);
    [corr, time_lags] = cross_corr_plot(x, t_svm, y, t_mmd);
    
    % also perform on shuffled data
    x = diag_svm_shuff - mean(diag_svm_shuff);
    y = diag_mmd_shuff - mean(diag_mmd_shuff);
    [corr2, time_lags2] = cross_corr_plot(x, t_svm_shuff, y, t_mmd_shuff);
    
    ax = gca; hold(ax, 'on');
    plot(time_lags, corr, 'Color', '#737373', 'LineStyle', '-', 'LineWidth', 3);
    plot(time_lags2, corr2, 'Color', '#737373', 'LineStyle', ':', 'LineWidth', 2);
    xline(0, 'Color', '#737373', 'LineStyle', '--', 'LineWidth', 1);
    % annotate the max point
    [maxv, maxi] = max(corr);
    plot(time_lags(maxi), corr(maxi), 'Marker', '+', 'MarkerSize', 16, 'Color', '#e41a1c', 'LineWidth', 3);
    text(time_lags(maxi), corr(maxi)+0.05, sprintf('%d ms', floor(time_lags(maxi)*1000)), 'FontSize', 14, 'Color', '#e41a1c');
    xlabel('Time lag (sec)', 'FontSize', 14);
    ylabel('Cross correlation', 'FontSize', 14);
    % ylabel('Normalized cross correlation', 'FontSize', 14);
    xlim([-0.1 0.1]);
    ylim([min(corr)-0.1 max(corr)+0.1]);
    title(sprintf('Cross-correlation: %s vs %s', syls{ii}, syls{jj}), 'FontSize', 16);
    
    fn_pdf = fullfile(fd_save_diag, sprintf('%s.neuralAcousticDiag.%s_%s.pdf', birdID,  syls{ii}, syls{jj}));
    print(fig, fn_pdf, '-dpdf', '-painters');
    
  end
end


















