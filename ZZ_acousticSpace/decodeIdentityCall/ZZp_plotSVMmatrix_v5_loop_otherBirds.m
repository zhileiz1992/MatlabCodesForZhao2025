% plot the calculated SVM decoding accuracy matrix for all pairwise comparison between call subtypes
% Zhilei, 09/11/2025
% differ from v4: 
% 1. Use the rename where low-quality clusters were identified, exclude those from plotting
% 2. Change metrics to a uniform similarity index (0 to 1)


clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');

birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};

% rename relationship
syls_all_list = {{'v1', 'v2', 'v3', 'v4', 'v5', 'v6'}, ...
             {'v1', 'v2', 'v3', 'v4'}, ...
             {'v1', 'v2', 'v3', 'v4', 'v5'}, ...
             {'v1', 'v2', 'v3', 'v4'}};

           
for bi=1:size(birdIDs,2)
birdID = birdIDs{bi};
pairID = pairIDs{bi};
disp(birdID);
% what good subtypes to analyze
syls = syls_all_list{bi};

col_list = {'#e78ac3', '#a6d854', '#fc8d62', '#8da0cb', '#66c2a5', '#e5c494', '#ccad25'};
for si=1:size(syls,2); col_dict.(syls{si})=col_list{si}; end
% where is calculated decoding matrix located
vae_run = 'traj_chop_32_1_32';
fd_data = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix5'], 'intermediate');
% what's the window size and hop size when calculating SVM decoding matrix, unit is seconds
fs = 20000;
pad = 0.08;
bin = 0.01;
hop = 0.002;
% where is the MMD result
fd_mmd = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', 'MMD', 'MatrixCall3', 'all', 'intermediate');
sigma = 0.9;
win_frame = 32;
hop_frame = 1;
sec_per_frame = 0.001;
% where to save the plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix5'], 'masked');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end
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
    if ~exist(fn_res, 'file'); continue; end
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
% fill nan into empty entry in the upper half
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    if isempty(neu2{ii,jj})
      widths = cellfun(@(x) size(x,1), neu2(ii,:));
      heights = cellfun(@(x) size(x,2), neu2(:,jj));
      temp = nan(mean(widths(widths~=0)), mean(heights(heights~=0)));
      neu2{ii,jj} = temp;
    end
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
fn_aco = fullfile(fd_save, sprintf('%s.d_aco.mat', birdID));
save(fn_aco, 'd_aco');


%% 2.1 Plot the full matrix
% neural SVM
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_n)-lens_n/2;
line_loc = cumsum(lens_n);
clim_n = [0.5 0.85];
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


%% 2.2. Change to a uniform similarity index
% neural SVM
vmin = 0.5; 
vmax = 1;
clim_n2 = [0.35 0.95];
neu4 = (vmax - neu3) / (vmax - vmin);
% neu4 = 2 * (1-neu3);
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_n)-lens_n/2;
line_loc = cumsum(lens_n);
ZZfunc_showMMDmatrix_v4(ax, neu4, clim_n2, v_loc, syls, 'Neural similarity', true, line_loc);
fn_pdf = fullfile(fd_save, sprintf('%s.neuralSVMmatrix.t%.3f_%.3f.c%.2f_%.2f.converted.pdf', birdID, tpad_n(1), tpad_n(2), clim_n2(1), clim_n2(2)));
print(fig, fn_pdf, '-dpdf', '-painters');

% acoustic MMD
vmin = -0.2; 
vmax = -0.03;
aco4 = (vmax - aco3) / (vmax - vmin);
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_a)-lens_a/2;
line_loc = cumsum(lens_a);
% clim_a = [-0.2 -0.035];
ZZfunc_showMMDmatrix_v4(ax, aco4, [0 1], v_loc, syls, 'Acoustic MMD', true, line_loc);
fn_pdf = fullfile(fd_save, sprintf('%s.acousticMMDmatrix.t%.3f_%.3f.c%.2f_%.2f.converted.pdf', birdID, tpad_n(1), tpad_n(2), vmin, vmax));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 2.3 Plot the neural matrix up-side-down
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
v_loc = cumsum(lens_n)-lens_n/2;
line_loc = cumsum(lens_n);
clim_n2 = [0.35 0.95];
ZZfunc_showMMDmatrix_v4_upsideDown(ax, neu4', clim_n2, v_loc, syls, 'Neural similarity', true, line_loc);
fn_pdf = fullfile(fd_save, sprintf('%s.neuralSVMmatrix.t%.3f_%.3f.c%.2f_%.2f.upsidedown.pdf', birdID, tpad_n(1), tpad_n(2), clim_n2(1), clim_n2(2)));
print(fig, fn_pdf, '-dpdf', '-painters');


%% 2.4 Resize the neural matrix to the size of acoustic MMD
neu5 = neu2; 
range_n = [0.5 1];  % what range to normalize to 
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    a = aco2{ii,jj};
    % normalize then interpolate the neural result
    n = neu2{ii,jj}; 
    n = (range_n(2)-n) / (range_n(2)-range_n(1));
    n = imresize(n, size(a), 'Method', 'bicubic');
    neu5{ii,jj} = n;
  end
end
neu6 = ZZfunc_tileArrayCell_v1(neu5);
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]);
ax=gca; hold(ax, 'on');
lens_n2 = cellfun(@(x) size(x,2), neu5(1,:));
v_loc = cumsum(lens_n2)-lens_n2/2;
line_loc = cumsum(lens_n2);
clim_n2 = [0.35 0.95];
ZZfunc_showMMDmatrix_v4_upsideDown(ax, neu6', clim_n2, v_loc, syls, 'Neural similarity', true, line_loc);
fn_pdf = fullfile(fd_save, sprintf('%s.neuralSVMmatrix.t%.3f_%.3f.c%.2f_%.2f.interpolated.pdf', birdID, tpad_n(1), tpad_n(2), clim_n2(1), clim_n2(2)));
print(fig, fn_pdf, '-dpdf', '-painters');



%% 3. Plot pairwise: acoustic and neural side-by-side;
% p_all = {{'v1', 'v2'}; {'v3', 'v4'}};
fd_save_pair = fullfile(fd_save, 'pairwise_plots');
if ~exist(fd_save_pair, 'dir'); mkdir(fd_save_pair); end
% loop through all pairs
range_n = [0.5 1];  % what range to normalize to 
range_a = [-0.2 -0.03];
clims_n2 = [0.35 0.95];
clims_a2 = [0 1];
% acoustic MMD
vmin = -0.2; 
vmax = -0.03;
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    
    close all;
    [fig, axes] = generatePanelGrid_v2(1, 2, [0.75], [], [0.15;0.05], [0.1;0.05], 0.1, [0], [10 10 1200 600]);
    % acoustic first
    ax = axes(1);
    a_this = aco2{ii, jj};
%     a_this = (range_a(2)-a_this) / (range_a(2)-range_a(1));
    a_this = (vmax - a_this) / (vmax - vmin);
    ty = rel_t_a_chop(1:size(a_this,1));
    tx = rel_t_a_chop(1:size(a_this,2));
    imagesc(ax, tx, ty, a_this, clims_a2);
    %   imagesc(ax, tx, ty, a_this, [0 1]);
    colormap(ax, flipud(gray)); colorbar(ax);
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
    n_this = (range_n(2)-n_this) / (range_n(2)-range_n(1));
    ty = rel_t_n_chop(1:size(n_this,1));
    tx = rel_t_n_chop(1:size(n_this,2));
    imagesc(ax, tx, ty, n_this, clims_n2);
    colormap(ax, flipud(gray)); colorbar(ax);
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



%% 3.2 Plot specific pair in desired orientation
p = {'v4', 'v3'}; % first syllable on the y-axis (height)
close all;
ii = find(strcmp(syls, p{1}));
jj = find(strcmp(syls, p{2}));
if ii<=jj
  a_this=aco2{ii,jj}; n_this = neu2{ii, jj};
else
  a_this=aco2{jj,ii}; n_this = neu2{jj, ii};
  a_this=a_this';  n_this=n_this';
end
a_this = (vmax - a_this) / (vmax - vmin);
n_this = (range_n(2)-n_this) / (range_n(2)-range_n(1));
[fig, axes] = generatePanelGrid_v2(1, 2, [0.75], [], [0.15;0.05], [0.1;0.05], 0.1, [0], [10 10 1200 600]);
% acoustic first
ax = axes(1);
ty = rel_t_a_chop(1:size(a_this,1));
tx = rel_t_a_chop(1:size(a_this,2));
imagesc(ax, tx, ty, a_this, clims_a2);
%   imagesc(ax, tx, ty, a_this, [0 1]);
colormap(ax, flipud(gray)); colorbar(ax);
set(ax, 'YDir', 'reverse');
axis(ax, 'equal');
axis(ax, 'tight');
ax.XAxisLocation = 'top';
xlabel(ax, sprintf('%s time (sec)', syls{jj}), 'FontSize', 20, 'Color', col_dict.(syls{jj}));
ylabel(ax, sprintf('%s time (sec)', syls{ii}), 'FontSize', 20, 'Color', col_dict.(syls{ii}));
title(ax, 'Acoustic similarity (dMMD)', 'FontSize', 20);
    
% then neural
ax = axes(2);
ty = rel_t_n_chop(1:size(n_this,1));
tx = rel_t_n_chop(1:size(n_this,2));
imagesc(ax, tx, ty, n_this, clims_n2);
colormap(ax, flipud(gray)); colorbar(ax);
set(ax, 'YDir', 'reverse');
axis(ax, 'equal');
axis(ax, 'tight');
ax.XAxisLocation = 'top';
xlabel(ax, sprintf('%s time (sec)', syls{jj}), 'FontSize', 20, 'Color', col_dict.(syls{jj}));
ylabel(ax, sprintf('%s time (sec)', syls{ii}), 'FontSize', 20, 'Color', col_dict.(syls{ii}));
title(ax, 'Neural similarity (SVM)', 'FontSize', 20);    
linkaxes(axes, 'xy');

fn_fig = fullfile(fd_save_pair, sprintf('%s.neuralSVM.desiredOrient.%s_%s.pdf', birdID, syls{ii}, syls{jj}));
print(fig, fn_fig, '-dpdf', '-painters');
























