% plot the calculated SVM decoding accuracy matrix for all pairwise comparison between call subtypes
% Zhilei, 09/05/2025
% differ from v1: use the latest SVM result where all renditions have the same warped duration

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
% where to save the plots
fd_save = fullfile(fd_base, 'Figures', pairID, 'AcousticSpace', birdID, 'embed', [vae_run '_decodePopMatrix4']);
fprintf('Save SVM plots to %s\n', fd_save);


%% 1. Read all matrix in
d = cell(size(syls,2), size(syls,2));  % full accu_all matrix
m = cell(size(syls,2), size(syls,2));  % mean accuracy matrix of the real dataset
% also trim to desired time range
rel_t = (0:1000)*hop - pad;  % a long time cooridiante axis
tpad = [0.032 -0.025];  % how much time before syllable onset and after syllable offset
m2 = cell(size(syls,2), size(syls,2));  % trimmed matrix
for ii=1:size(syls,2)
  for jj=ii:size(syls,2)
    run_name = sprintf('%s_%s', syls{ii}, syls{jj});
    fn_res = fullfile(fd_data, sprintf('%s.%s.accu_allMatrix.mat', birdID, run_name));
    a = load(fn_res); temp = a.accu_all;
    d{ii, jj} = temp;
    m_this = mean(squeeze(temp(:,:,1,:)), 3);
    m{ii, jj} = m_this;
    % trim
    ty = rel_t(1:size(m_this,1));
    tx = rel_t(1:size(m_this,2));
    tendy = ty(end) - pad + tpad(2);
    iy = find((ty>=-tpad(1)) & (ty<=tendy));
    tendx = tx(end) - pad + tpad(2);
    ix = find((tx>=-tpad(1)) & (tx<=tendx));
    m2{ii,jj} = m_this(iy, ix);
  end
end

% tile into a large 2d matrix
lens = cellfun(@(x) size(x,2), m2(1,:));
cumlens = cumsum(lens);
m3 = nan(sum(lens), sum(lens));
for ii=1:size(m2,1)
  for jj=ii:size(m2,2)
    if ii==1; istart=1; else; istart=cumlens(ii-1)+1;end
    if jj==1; jstart=1; else; jstart=cumlens(jj-1)+1;end
    iend=cumlens(ii); jend=cumlens(jj);
    m3(istart:iend, jstart:jend) = m2{ii, jj};
  end
end
  


%% 2. Plot the full matrix
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]); 
ax=gca; hold(ax, 'on');
% imagesc(ax, m3, [0.5 0.9]); 
% colormap gray; colorbar;
v_loc = cumsum(lens)-lens/2;
line_loc = cumsum(lens);
ZZfunc_showMMDmatrix_v3(ax, m3, [0.5 0.9], v_loc, syls, 'Neural SVM accuracy', true, line_loc);

fn_pdf = fullfile(fd_save, sprintf('%s.neuralSVMmatrix.pdf', birdID));
print(fig, fn_pdf, '-dpdf', '-painters');



%% 3. Plot targeted pairs
% relative time cooridinates after setting the time range
rel_t_chop = rel_t(rel_t>=-tpad(1));
p_all = {{'v1', 'v2'}; {'v3', 'v4'}};
for pi=1:size(p_all, 1)
  p = p_all{pi};
  vi = cellfun(@(x) find(strcmp(syls, x)), p);
  m_this = m2{vi(1), vi(2)};
  ty = rel_t_chop(1:size(m_this,1)); 
  tx = rel_t_chop(1:size(m_this,2)); 
  close all; 
  fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]); 
  ax = gca; hold(ax, 'on');
  imagesc(ax, tx, ty, m_this, [0.5 0.9]);
  colormap(ax, 'gray'); colorbar(ax);
  set(ax, 'YDir', 'reverse');
  axis(ax, 'equal');
  axis(ax, 'tight');
  ax.XAxisLocation = 'top';
  xlabel(ax, sprintf('%s time (sec)', p{2}), 'FontSize', 20, 'Color', col_dict.(p{2}));
  ylabel(ax, sprintf('%s time (sec)', p{1}), 'FontSize', 20, 'Color', col_dict.(p{1}));
  fn_fig = fullfile(fd_save, sprintf('%s.neuralSVM.%s_%s.pdf', birdID, p{1}, p{2}));
  print(fig, fn_fig, '-dpdf', '-painters');
end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  










