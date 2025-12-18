% plot the calculated SVM decoding accuracy matrix for all pairwise comparison between call subtypes
% Zhilei, 09/04/2025

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

% align the matrix along rows and columns; discard small difference in lengths introduced by selecting different neuron
% populations to calculate SVM between different call pairs
row_lens = zeros(size(m2,1),1);
col_lens = zeros(size(m2,2),1);
for ii=1:size(m2,1)
  lens = cellfun(@(x) size(x,1), m2(ii,:));
  row_lens(ii) = min(lens(lens~=0));
  lens2 = cellfun(@(x) size(x,2), m2(:,ii));
  col_lens(ii) = min(lens2(lens2~=0));
end
cum_rows = cumsum(row_lens);
cum_cols = cumsum(col_lens);
m3 = nan(sum(row_lens), sum(col_lens));
for ii=1:size(m2,1)
  for jj=ii:size(m2,2)
    if ii==1; istart=1; else; istart=cum_rows(ii-1)+1;end
    if jj==1; jstart=1; else; jstart=cum_cols(jj-1)+1;end
    iend=cum_rows(ii); jend=cum_cols(jj);
    temp = m2{ii,jj};
    m3(istart:iend, jstart:jend) = temp(1:row_lens(ii), 1:col_lens(jj));
  end
end



%% 2. Plot the matrix
close all;
fig = ZZfunc_newFigurePDFsize_v1([10 10 1000 800]); 
ax=gca; hold(ax, 'on');
% imagesc(ax, m3, [0.5 0.9]); 
% colormap gray; colorbar;
v_loc = cumsum(row_lens)-row_lens/2;
line_loc = cumsum(row_lens);
ZZfunc_showMMDmatrix_v3(ax, m3, [0.5 0.9], v_loc, syls, 'Neural SVM', true, line_loc);

fn_pdf = fullfile(fd_save, sprintf('%s.neuralSVMmatrix.pdf', birdID));
print(fig, fn_pdf, '-dpdf', '-painters');












