% test the co-occurence between high acoustic similarity and high neural similarity
% Zhilei, 09/05/2025

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
fd_data = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix4']);
% where to save the plots
fd_save = fullfile(fd_data, 'co-occur');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
fprintf('Save SVM plots to %s\n', fd_save);

% load data
fn_neu = fullfile(fd_data, sprintf('%s.d_neu.mat', birdID));  load(fn_neu);
fn_aco = fullfile(fd_data, sprintf('%s.d_aco.mat', birdID));  load(fn_aco);



%% 1. Find regions where acoustic similarity is high
a_thre = -0.035;
close all;
n = size(syls,2);
[fig, axes] = generatePanelGrid_v2(n, n, zeros(7,1)+0.12, zeros(6,1)+0.01, [0.05;0.05], [0.05;0.05], 0.01, zeros(7,1), [10 10 900 900]);
[fig2, axes2] = generatePanelGrid_v2(n, n, zeros(7,1)+0.12, zeros(6,1)+0.01, [0.05;0.05], [0.05;0.05], 0.01, zeros(7,1), [10 10 900 900]);
[fig3, axes3] = generatePanelGrid_v2(n, n, zeros(7,1)+0.12, zeros(6,1)+0.01, [0.05;0.05], [0.05;0.05], 0.01, zeros(7,1), [10 10 900 900]);
[fig4, axes4] = generatePanelGrid_v2(n, n, zeros(7,1)+0.12, zeros(6,1)+0.01, [0.05;0.05], [0.05;0.05], 0.01, zeros(7,1), [10 10 900 900]);
aco2 = d_aco.aco2;
tpad_a = [0 -0.032];  % what time range to check, extra from syllable onset and offset
rel_t_a_chop = d_aco.rel_t_a_chop;
neu2 = d_neu.neu2; 
tpad_n = [0 -0.032];
rel_t_n_chop = d_neu.rel_t_n_chop;
% save SVM accuracy within ROI and randomly choosen region
roi_all = [];
rand_all = [];
rng(1992);

for ii=1:size(syls,2)
  for jj=(ii+1):size(syls,2)
    a_this = aco2{ii,jj};
    tya = rel_t_a_chop(1:size(a_this,1));
    txa = rel_t_a_chop(1:size(a_this,2));
    % restrict to desired time range
    tendy = tya(end) + tpad_a(2);  % note that for MMD, no more extra slide at the end, sliding window stops at the syllable offset
    iy = find((tya>=-tpad_a(1)) & (tya<=tendy));
    tendx = txa(end) + tpad_a(2);
    ix = find((txa>=-tpad_a(1)) & (txa<=tendx));
    a = a_this(iy, ix);
    
    % original MMD
    txaa = txa(ix); tyaa = tya(iy);
    ax = axes(ii, jj); cla(ax); 
    imagesc(ax, txaa, tyaa, a, [-0.2 -0.02]); 
    colormap(ax, 'gray'); 
    set(ax, 'YDir', 'reverse');
    axis(ax, 'off');
    
    % hight acoustically similar regions
    ax2 = axes2(ii, jj); cla(ax2); hold(ax2, 'on');
    imagesc(ax2, txaa, tyaa, a, [-0.2 -0.02]); 
    colormap(ax2, 'gray'); 
    set(ax2, 'YDir', 'reverse');
    % plot points that pass the threshold
    [rows, cols] = find(a <= a_thre);
    plot(ax2, txaa(cols), tyaa(rows), 'r.', 'MarkerSize', 1); % Plot red dots at those points
    axis(ax2, 'off');
    
    % show the neural SVM in the same time range
    ax3 = axes3(ii, jj); cla(ax3); hold(ax3, 'on');
    n_this = neu2{ii,jj};
    tyn = rel_t_n_chop(1:size(n_this,1));
    txn = rel_t_n_chop(1:size(n_this,2));
    % restrict to the same time range as in the acoustic plot
    ix = find((txn>=txaa(1)) & (txn<=txaa(end)));
    iy = find((tyn>=tyaa(1)) & (tyn<=tyaa(end)));
    s = n_this(iy, ix);
    txnn=txn(ix); tynn=tyn(iy);
    imagesc(ax3, txnn, tynn, s, [0.5 0.9]);   
    colormap(ax3, 'gray'); 
    set(ax3, 'YDir', 'reverse');
    axis(ax3, 'off');
    
    % find the acoustically similar region in the SVM plot
    ax4 = axes4(ii, jj); cla(ax4); hold(ax4, 'on');
    imagesc(ax4, txnn, tynn, s, [0.5 0.9]);   
    colormap(ax4, 'gray'); 
    set(ax4, 'YDir', 'reverse');
    % find the nearest point in the SVM matrix
    % Create meshgrids for interpolation
    [Xa, Ya] = meshgrid(txaa, tyaa);
    [Xs, Ys] = meshgrid(txnn, tynn);
    % Create mask for the region of interest in a
    mask_a = a < a_thre;
    % Interpolate the mask to the grid of s, using nearest neighbor and extrapolate outside as 0 (not ROI)
    mask_s_interp = interp2(Xa, Ya, double(mask_a), Xs, Ys, 'nearest', 0);
    % Find indices in s where the interpolated mask indicates the ROI
    [roi_y, roi_x] = find(mask_s_interp == 1);
    % Optional: Plot s with the ROI highlighted, similar to the original
    roi_s = s(sub2ind(size(s), roi_y, roi_x));
    % also sample randomly
    rand_y = randsample(1:size(s,1), length(roi_s), true);
    rand_x = randsample(1:size(s,2), length(roi_s), true);
    rand_s = s(sub2ind(size(s), rand_y', rand_x'));
    
    plot(ax4, txnn(roi_x), tynn(roi_y), 'r.', 'MarkerSize', 1); 
%     plot(ax4, txnn(rand_x), tynn(rand_y), 'g.', 'MarkerSize', 1); 
    axis(ax4, 'off');

    roi_all = [roi_all; roi_s];
    rand_all = [rand_all; rand_s];
  end
end
% save the figures; 
fig_all = {fig; fig2; fig3; fig4};
for fi=1:size(fig_all,1)
  f_this = fig_all{fi};
  fn_pdf = fullfile(fd_save, sprintf('%s.co-occur.thre%.3f.fig%d.pdf', birdID, a_thre, fi));
  print(f_this, fn_pdf, '-dpdf', '-painters');
end



%% 2. Compare the SVM accuracy
close all; 
fig5 = ZZfunc_newFigurePDFsize_v1([50 50 250 400]); 
% boxplot(all_data, group_labels);
all_data = [roi_all; rand_all];
group_labels = [repmat({'1.real'}, numel(roi_all), 1); repmat({'2.random'}, numel(rand_all), 1)];
col_violin = [0.3 0.7 0.3; 0.7 0.7 0.7];
violinplot(all_data, group_labels, 'ViolinColor', col_violin, 'EdgeColor', [1 1 1], 'MedianMarkerSize', 300, 'ShowData', false, 'ViolinAlpha', 0.5);
ylim([0.45 1.05]);
xlim([0.5 2.5]);
ylabel('Neural SVM accuracy', 'FontSize', 14);
fn_pdf = fullfile(fd_save, sprintf('%s.co-occur.thre%.3f.violinplot.pdf', birdID, a_thre));
print(fig5, fn_pdf, '-dpdf', '-painters');

% Perform Wilcoxon rank-sum test
[p, h, stats] = ranksum(roi_all, rand_all);
% Output results
fprintf('Wilcoxon rank-sum test:\n');
fprintf('p-value = %.4f\n', p);
fprintf('Test statistic (z-value) = %.4f\n', stats.zval); % zval is only available for large samples









