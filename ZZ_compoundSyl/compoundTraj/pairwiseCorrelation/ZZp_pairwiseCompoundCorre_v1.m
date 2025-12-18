% generate cross-correlation plot between pairs of compound syllables

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_acousticSpace/calcDistance/MMD'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZBudgieIntanExtract'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% size of sliding window
win_size = 32;  % in unit of frames
sec_per_frame = 0.001; 
% coloring for syllables
syl_full = {'v', 'h', 'e', 'b', 'x'};
col_full = {'#e41a1c', '#984ea3', '#4daf4a', '#377eb8', '#737373'};
col_dict = struct;
for si=1:size(syl_full,2) col_dict.(syl_full{si}) = col_full{si}; end

% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 2;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% 1. Load VAE and UMAP data
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'ApplySylAll', sprintf('latents.%s', vae_run));
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
info_vae = readtable(fn_info, 'Delimiter', ',');
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
% UMAP
fd_umap_base = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'UMAPcomp');
umap_run = 'v.b.h.e.x.n-1';
syls = {'v', 'b', 'h', 'e', 'x'};
fd_umap = fullfile(fd_umap_base, umap_run);
fn_info_umap = fullfile(fd_umap, sprintf('%s.%s.info.csv', birdID, umap_run));
info_umap = readtable(fn_info_umap, 'Delimiter', ',');
fn_embed = fullfile(fd_umap, sprintf('%s.%s.embedding.csv', birdID, umap_run));
embed = readmatrix(fn_embed);

% where to save results
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
if ~exist(fd_save, 'dir'); mkdir(fd_save); end;


%% 2. Focus on compound syllables
fs = 20000;
info_vae.dur = (info_vae.iend - info_vae.istart) / fs;
dur_thre = 0.3;  % unit is sec
comp = info_vae(info_vae.dur>=dur_thre & ismember(info_vae.call_subtype, {'b', 'x'}) & ismember(info_vae.batch, {'batch1'}), :);

% grab the VAE latents
for ii=1:size(comp, 1)
  syl_ID = comp.syl_ID{ii};
  dvae = h5read(fn_vae, ['/' syl_ID]);
  comp.dvae{ii} = dvae';
end
% save for later use
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
save(fn_comp, 'comp', '-v7.3');


%% 3. Calculate pairwise correlation coefficient for a few examples, plot to check
% use distance in the VAE latents
% identify the strands with high similarity
fd_plot = fullfile(fd_save, 'examplePlot');
if ~exist(fd_plot, 'dir'); mkdir(fd_plot); end
rng(1992);
for rand_i=1:50 
% tic;
idx_rd = randsample(1:size(comp,1), 2);

% read wav then plot spectrogram
pad_pre = win_size*sec_per_frame; pad_pre_pt = floor(fs*pad_pre);  % pad equal amount as in VAE sliding
pad_post = win_size*sec_per_frame; pad_post_pt = floor(fs*pad_post);
sounds = cell(2, 1);
for ii=1:length(idx_rd)
  idx = idx_rd(ii);
  [signal, fs] = audioread(comp.fn_wav{idx});
  i_start = max([1 comp.istart(idx)-pad_pre_pt]);
  i_end = min([size(signal,1) comp.iend(idx)+pad_post_pt]);
  sounds{ii} = signal(i_start:i_end);
end
  
close all;
fig_size = [10 10 900 900];
left_margin = 0.05;       % Space from left edge of figure
bottom_margin = 0.05;     % Space from bottom edge of figure
left_col_width = 0.2;    % Width of left column (for B)
% right_col_width = 0.4;   % Width of right column (for A and C; A and C will have same width)
top_row_height = 0.2;    % Height of top row (for A)
% determine bottom height proportionally
[power1, ~, ~, ~, ~, ~] = getAudioSpectrogramZZ_flexible_v1(sounds{1}, fs, 256, 256, 236, [250 7500], [12 23]);
[power2, ~, ~, ~, ~, ~] = getAudioSpectrogramZZ_flexible_v1(sounds{2}, fs, 256, 256, 236, [250 7500], [12 23]);
% bottom_row_height = 0.5; 
if size(power1,2)>size(power2,2)
  right_col_width = 0.5;
  bottom_row_height = right_col_width * size(power2,2) / size(power1,2);
else
  bottom_row_height = 0.5; 
  right_col_width = bottom_row_height * size(power1,2) / size(power2,2);
end

[fig, axes_all] = ZZ_emptyCrossCorrePlot_v1(fig_size, left_margin, bottom_margin, left_col_width, right_col_width, top_row_height, bottom_row_height);
[ax1, ~, ~,  ~, ~,  t1] = showAudioSpectrogramZZ_flexible_v1(sounds{1}, fs, axes_all(1), [250 7500], [12 23], 256, 256, 236);
axis(axes_all(1), 'off');
[ax2, ~, ~,  ~, ~,  t2] = showAudioSpectrogramZZ_flexible_v1_flip(sounds{2}, fs, axes_all(2), [250 7500], [12 23], 256, 256, 236);
axis(axes_all(2), 'off');
  
% computate distance between every two time points of the syllable
d1 = comp.dvae{idx_rd(1)};
d2 = comp.dvae{idx_rd(2)};
% distm = pdist2(d1, d2, 'euclidean');
distm = pdist2(d1, d2, 'cosine');
distm = distm'; 
% decompose 
% [B, C, u, v] = mmd_matrix_decompose(distm);
% imagesc(axes_all(3), t1, t2, distm');
imagesc(axes_all(3), t1, t2, distm, [0 0.5]);
colormap(axes_all(3), gray);
axis(axes_all(3), 'off');

fn_fig = fullfile(fd_plot, sprintf('exampleFig.%d.pdf', rand_i));
print(fig, fn_fig, '-dpdf', '-painters');


% identify the high-similarity strands
% distm2 = smoothdata(distm, 'gaussian', 21);
distm2 = medfilt2(distm, [11 11]);
thre = 0.35; 
binary_mask = distm2 <= thre;
cc = bwconncomp(binary_mask);
close all;
[fig2, ax_all2] = generatePanelGrid_v2(1, 4, [0.7], [], [0.05;0.05], [0.05;0.05], 0.05, [0], [10 10 2000 600]);
% plot the original
imagesc(ax_all2(1), distm, [0 0.5]);
colormap(ax_all2(1), gray); title(ax_all2(1), 'Original matrix'); 
% plot the smooth
imagesc(ax_all2(2), distm2, [0 0.5]);
colormap(ax_all2(2), gray); title(ax_all2(2), 'Median filtered'); 
% plot the hard mask
imagesc(ax_all2(3), binary_mask);
title(ax_all2(3), sprintf('Passed threshold %.2f', thre)); 
% filter out based on min dur
min_dur = 25;  % at least 25 frames on each side, i.e. 25 ms
c_all = cc.PixelIdxList;
pass_i = {};
count = 0;
for ci=1:size(c_all, 2)
  linear_i = c_all{ci};
  % convert to xy index
  [xx, yy] = ind2sub(size(distm2), linear_i);
  x_range = max(xx) - min(xx); 
  y_range = max(yy) - min(yy); 
  if x_range>=min_dur &&  y_range>=min_dur
    count = count + 1;
    pass_i{count} = linear_i;
  end
end
% plot those that pass the size limit
newMask = zeros(size(distm2));
for ci=1:size(pass_i, 2)
  newMask(pass_i{ci}) = 1; 
end
imagesc(ax_all2(4), newMask);
title(ax_all2(4), sprintf('Passed duration limit %d ms', min_dur)); 

fn_fig = fullfile(fd_plot, sprintf('exampleFig.%d.strand.pdf', rand_i));
print(fig2, fn_fig, '-dpdf', '-painters');


end
  

%% 4. Loop over many randomly selected pairs
fd_save_loop = fullfile(fd_save, 'intermediate');
if ~exist(fd_save_loop, 'dir'); mkdir(fd_save_loop); end

num_rd = 1e5;
n = size(comp,1);
rng(1992);
linear_idx = randsample(n*n, num_rd);
[i, j] = ind2sub([n, n], linear_idx);
% Filter out diagonal (i == j)
mask = i ~= j; i = i(mask); j = j(mask);
comp_v = comp.dvae; 
data1 = comp_v(i);
data2 = comp_v(j);

% do it in batches to prevent RAM blowup
batch_size = 1e4; 
num_batch = ceil(num_rd / batch_size);
for bai=1:num_batch
  i_s = (bai-1)*batch_size + 1; 
  i_e = bai * batch_size; 
  res = cell(batch_size, 1);
  parfor idx = (i_s+1):(i_e+1)
      % d1 and d2 are now pulled from SLICED variables
      d1 = data1{idx-1}; 
      d2 = data2{idx-1};  
      distm = pdist2(d1, d2, 'cosine');  
      res{idx - i_s} = distm;
  end
  fn_res = fullfile(fd_save_loop, sprintf('batch%d.res.mat', bai)); 
  save(fn_res, 'res', '-v7.3'); 
  clear res;
end
 


%% 5. Select one syllable as the reference, loop over other syllables
% then for a given position in the reference syllable, find similarity strands to all other target syllables
% overlay the location of these similarity strands onto the reference syllable
ref_i = 8004; 
% ref_i = 50;
fd_save_ref = fullfile(fd_save, 'ref_tar', sprintf('ref_%d', ref_i));
if ~exist(fd_save_ref, 'dir'); mkdir(fd_save_ref); end

% tic;
to_loop = 1000; 
i_all = 1:size(comp,1); 
rng(1992);
rd_list = randsample(i_all(i_all~=ref_i), to_loop);
pass_all = [];
d1 = comp.dvae{ref_i};
d2s = comp.dvae(rd_list);
parfor rdi=1:length(rd_list)
%   idx_rd = [ref_i rd_list(rdi)];
%   idx_rd = [ref_i 3922];
%   close all;
  % calculate then plot the cross-correlation 
%   [fig, axes_all, distm, t1, t2] = ZZfunc_plotCrossCorreTwoSpec_v1(comp, idx_rd);
%   fn_fig = fullfile(fd_save_ref, sprintf('ref%d.tar%d.pdf', ref_i, idx_rd(2)));
%   print(fig, fn_fig, '-dpdf', '-painters');
  % no plotting
%   d1 = comp.dvae{idx_rd(1)};
%   d2 = comp.dvae{idx_rd(2)};
  % distm = pdist2(d1, d2, 'euclidean');
  distm = pdist2(d1, d2s{rdi}, 'cosine');
  distm = distm'; 
  
  % identify the similarity strands
  med_filter = [11 11];
  thre = 0.35; 
  min_dur = 25; 
%   [fig2, ax_all2, pass_i, newMask] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur); 
%   fn_fig = fullfile(fd_save_ref, sprintf('ref%d.tar%d.strand.pdf', ref_i, idx_rd(2)));
%   print(fig2, fn_fig, '-dpdf', '-painters');
  % no plotting
  pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, med_filter, thre, min_dur);

  pass_all(rdi).ref_i = ref_i;
  pass_all(rdi).tar_i = rd_list(rdi);
  pass_all(rdi).distm = distm; 
  pass_all(rdi).pass_i = pass_i;
end
% toc;

% overlay the location of the similarity strands on reference syllable
close all;
[fig3, axes3] = generatePanelGrid_v2(5, 1, [0.12;0.3;0.12;0.12;0.12], [0.02;0.02;0.02;0.02], [0.05;0.05], [0.12;0.05], 0.05, [0;0;0;0;0], [10 10 700 900]);
% first plot spectrogram of reference syllable, include the half-width of the sliding window
pad_pre = win_size*sec_per_frame/2; pad_pre_pt = floor(fs*pad_pre);  % pad equal amount as in VAE sliding
pad_post = win_size*sec_per_frame/2; pad_post_pt = floor(fs*pad_post);
[signal, fs] = audioread(comp.fn_wav{ref_i});
i_start = max([1 comp.istart(ref_i)-pad_pre_pt]);
i_end = min([size(signal,1) comp.iend(ref_i)+pad_post_pt]);
sound = signal(i_start:i_end);
ax1 = axes3(1);
[ax1, ~, ~,  ~, ~,  t1] = showAudioSpectrogramZZ_flexible_v1(sound, fs, ax1, [250 7500], [12 23], 256, 256, 236);
title(ax1, sprintf('Compound syl. %d', ref_i), 'FontSize', 14);
xticks(ax1, []);

% then plot where the strands are
ax2 = axes3(2); cla(ax2); hold(ax2, 'on'); 
count = 0;
xall = {};
yall = {};
for si=1:size(pass_all,2)
  p = pass_all(si);
  pass_i = p.pass_i;
  for pi=1:size(pass_i,2)
    strand = pass_i{pi};
    [yy, xx] = ind2sub(size(p.distm), strand);
    x_range = [min(xx) max(xx)];
    y_range = [min(yy) max(yy)];
    % add the half-width of sliding windows back, so considering the center of sliding window
    x_range = x_range + win_size/2;
    y_range = y_range + win_size/2;
    count = count + 1; 
    plot(ax2, x_range, [count count], 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black');
    % construct indicator varialbe
    x_indi = zeros(1, size(p.distm, 2));
    x_indi(x_range(1):x_range(2)) = 1;
    xall = [xall; x_indi];
    y_indi = zeros(1, size(p.distm, 1));
    y_indi(y_range(1):y_range(2)) = 1;
    yall = [yall; y_indi];
  end   
end
xlim(ax2, [0 size(p.distm, 2)]);
ylim(ax2, [0 count]);
ylabel(ax2, 'Strand ID', 'FontSize', 12); 
% xlabel(ax2, 'Rel. time (ms)', 'FontSize', 12);
xticks(ax2, []);

% then calculate and plot probability
ax3 = axes3(3); 
% truncate what's longer than the syllable VAE dims
n = size(d1, 1);
xall2 = cellfun(@(x) x(1:n), xall, 'UniformOutput', false);
x_cat = vertcat(xall2{:}); 
x_prob = sum(x_cat, 1) / size(x_cat, 1);
plot(ax3, x_prob); 
xlim(ax3, [0 length(x_prob)]);
ylabel(ax3, 'Probability', 'FontSize', 12);     
% xlabel(ax3, 'Rel. time (ms)', 'FontSize', 12);
xticks(ax3, []);

% plot the sound amplitude
win_ms = 10; 
cut_env = 60;
[amp, ~] = ZZfunc_soundAmplitude_v1(sound, fs, [250 5000], win_ms, cut_env);
ax4 = axes3(4);
plot(ax4, amp);
xlim(ax4, [1 length(amp)]);
ylabel(ax4, 'Amplitude', 'FontSize', 12);
xticks(ax4, []);
    
% multiply probablity with sound amplitude
ref_db = 35;
amp2 = amp + ref_db;
amp2(amp2<0) = 0;
% scale between 0 and 1
amp3 = amp2 / ref_db; 
% interpolate to the same length
x_original = 1:length(amp3);
x_new = linspace(1, length(amp3), length(x_prob));
amp_interp = interp1(x_original, amp3, x_new);
% multiply probablity with amplitude
prob_amp = x_prob .* amp_interp;
ax5 = axes3(5); cla(ax5); hold(ax5, 'on');
plot(ax5, prob_amp);
xlim(ax5, [1 length(prob_amp)]);
ylabel(ax5, 'Prob * Amp', 'FontSize', 12);


% find peaks in the prob*amp curve
param.gapSize = 5;  
param.thresholdFlatness = 0.1;  % first threshold to identify peak; initially use 0.03
param.extendFlatness = 0.02; % after thresholding extend to lower threshold, initially 0.01
param.minInterval = 0.01;  % min distance between two adjacent peaks, smaller will be merged
param.minDuration = 0.02;  % minimal duration of the identified primitives
param.maxDuration = 0.2; 

flatness=-prob_amp; dt=0.001; 
[onsets, offsets] = ZZ_GetFlatnessOnsetOffset(flatness, dt, param);

% plot in the curve
% yline(ax5, param.thresholdFlatness, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
% yline(ax5, param.extendFlatness, 'Color', 'green', 'LineStyle', '--', 'LineWidth', 0.5);
if ~isempty(onsets)
  for ii=1:length(onsets)
      xline(ax5, onsets(ii), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 0.5);
      xline(ax5, offsets(ii), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
      % also plot in the spectrogram
      xline(ax1, onsets(ii)/length(prob_amp), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 0.5);
      xline(ax1, offsets(ii)/length(prob_amp), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
  end
end

% save figure
fn_fig = fullfile(fd_save_ref, sprintf('ProbHigh.ref%d.pdf', ref_i));
print(fig3, fn_fig, '-dpdf', '-painters');
  








