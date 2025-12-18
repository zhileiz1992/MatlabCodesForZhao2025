% compare compound syllables to calls, generate the cross-distance matrix, check if they share primitives
% Zhilei, 09/23/2025
% differ from v3: troubleshoot the time scale difference on calls and compound syllables


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
clims = {[10 21]; [10 21]; [10.5 21.5]; [12 23]};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% size of sliding window
win_size = 32;  % in unit of frames
sec_per_frame = 0.001;

% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet


bi = 1;
birdID = birdIDs{bi};
pairID = pairIDs{bi};
clim = clims{bi};


%% 1. Load compound syllable data
fs = 20000;
fd_save = fullfile(fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre'));
fn_comp = fullfile(fd_save, sprintf('%s.comp.mat', birdID));
load(fn_comp);

% grab VAE data for calls
fd_vae = fullfile(fd_base, 'vaeWav', birdID, 'Traj', 'ApplySylAll', sprintf('latents.%s', vae_run));
fn_info = fullfile(fd_vae, sprintf('%s.latents.%s.info.csv', birdID, vae_run));
info_vae = readtable(fn_info, 'Delimiter', ',');
fn_vae = fullfile(fd_vae, sprintf('%s.latents.%s.h5', birdID, vae_run));
% subset the information table: deal with the ID swapping
oldID = {'v4', 'v5', 'v1', 'v7', 'v3', 'v2', 'v6'};
newID = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7'};
for si=1:size(newID,2) ID_dict.(newID{si}) = oldID{si}; end
% only look at well-defined calls
call_syl = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'};
call_syl_old = cellfun(@(x) ID_dict.(x), call_syl, 'UniformOutput', false);
% subset the information table
calls = info_vae(ismember(info_vae.call_subtype, call_syl_old) & ismember(info_vae.batch, {'batch1'}), :);

% grab the VAE latents and sound data (pad half-width of the sliding window)
% pad = 0.016; pad_pt=floor(pad*fs);
pad_sound = 0.032; pad_sound_pt = floor(pad_sound*fs);
for ii=1:size(calls, 1)
  % for ii=1:10
  syl_ID = calls.syl_ID{ii};
  dvae = h5read(fn_vae, ['/' syl_ID]);
  calls.dvae{ii} = dvae';
  % also add the new ID to table
  call_i = find(strcmp(oldID, calls.call_subtype{ii}));
  calls.newCallID{ii} = newID{call_i};
  % then grab the sound data
  [signal, fs] = audioread(calls.fn_wav{ii});
  i_start = max([1 calls.istart(ii)-pad_sound_pt]);
  i_end = min([size(signal,1) calls.iend(ii)+pad_sound_pt]);
  calls.sound{ii} = signal(i_start:i_end);
end
% save for later use
fn_call = fullfile(fd_save, sprintf('%s.call3.mat', birdID));
save(fn_call, 'calls', '-v7.3');



%% 2. Time-warped call renditions and VAE latents; calculate mean
call_mean = [];
for call_i=1:size(call_syl, 2)
  v = call_syl{call_i};
  % warp VAE latents
  rends = calls(strcmp(calls.newCallID, v),:);
  lens = cellfun(@(x) size(x,1), rends.dvae);
  mean_len = floor(mean(lens));
  fprintf('Intepolate %s to %d frames (ms)\n', v, mean_len);
  vae_interp = cellfun(@(x) ZZfunc_interp_window(x, mean_len), rends.dvae, 'UniformOutput', false);
  vae_interp = cat(3, vae_interp{:});
  call_mean(call_i).v = v;
  call_mean(call_i).num_rends = size(rends, 1);
  call_mean(call_i).vae_mean = mean(vae_interp, 3);
  % warp sound
  lens2 = cellfun(@(x) size(x,1), rends.sound);
  mean_len2 = floor(mean(lens2));
  sound_interp = cellfun(@(x) ZZfunc_interpolate1Darray_v1(x, mean_len2), rends.sound, 'UniformOutput', false);
  %   sound_interp = cat(2, sound_interp{:});
  % calculate mean spec
  spec_all = cellfun(@(x) getAudioSpectrogramZZ_flexible_v1_simpleReturn(x, fs, 256, 256, 236, [250 7500], clim), sound_interp, 'UniformOutput', false);
  spec_all = cat(3, spec_all{:});
  spec_mean = mean(spec_all, 3);
  spec_mean = spec_mean/1024;  %normalize to between 0 and 1
  call_mean(call_i).spec_mean = spec_mean;
  call_mean(call_i).mean_len2 = mean_len2; 
  % get the relative time of spectrogram columns
  [~, ~, ~, ~, ~, t_mean_spec] = getAudioSpectrogramZZ_flexible_v1(sound_interp{1}, fs, 256, 256, 236, [250 7500], clim);
  call_mean(call_i).t_mean_spec = t_mean_spec;
end
% save for later use
fn_call_mean = fullfile(fd_save, sprintf('%s.call_mean3.mat', birdID));
save(fn_call_mean, 'call_mean');



%% 3. Loop compound syllable, compare to the mean VAE of a given call
fd_save_ref = fullfile(fd_save, 'compound_call3_temp');

% choose different call subtypes as the reference
vi = 1; 
% vi = 6;
fd_save_this = fullfile(fd_save_ref, sprintf('v%d', vi));
if exist(fd_save_this, 'dir'); rmdir(fd_save_this, 's'); end
mkdir(fd_save_this); 

% loop over all compound syllables
pass_list = [];
% parameters for identifying bands
med_filter = [11 11];
thre = 0.2;
min_dur = 25;
d1 = call_mean(vi).vae_mean;
% remove part of the pad that's silence
half_win = floor(win_size/2);
rm_size = 12; 
d1c = d1((rm_size+1):(size(d1,1)-rm_size), :);
d2s = comp.dvae; 
parfor tar_i=1:size(comp,1)
  % grab VAE latents
  d2 = d2s{tar_i};
  % no plotting
  distm = pdist2(d1c, d2, 'cosine');
  distm = distm';
  % no plotting
  pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, med_filter, thre, min_dur);
  pass_list(tar_i).tar_i = tar_i;
  pass_list(tar_i).size_d = size(distm);
  pass_list(tar_i).pass_i = pass_i;
end
% save result
fd_pass_list = fullfile(fd_save_ref, sprintf('call_v%d.thre%.2f.pass_list.mat', vi, thre));
save(fd_pass_list, 'pass_list');



%% 6. Identify the primitives using density
% a black bar for similarity bands
n = pass_list(1).size_d;
m = cellfun(@(x) size(x,2), {pass_list.pass_i});
xall = zeros(sum(m), n(2));
count = 0;
for ri=1:size(pass_list,2)
  p = pass_list(ri);
  pass_i = p.pass_i;
  for pi=1:size(pass_i,2)
    [yy, xx] = ind2sub(p.size_d, pass_i{pi});
    x_range = [min(xx) max(xx)];
    count = count + 1;
    xall(count, x_range(1):x_range(2)) = 1;
  end   
end

% pad half-width back to the matrix
% xpad = zeros(size(xall,1), half_win);
% xall2 = [xpad xall xpad];
% calculate density
x_prob = sum(xall, 1) / size(pass_list,2);
ref_t = (1:length(x_prob)) - (half_win - rm_size); 


%% 6.2 plot results regarding the primitives
close all;
[fig3, axes3] = generatePanelGrid_v2(3, 1, [0.17;0.3;0.15], [0.02;0.02], [0.05;0.05], [0.25;0.05], 0.05, [0;0;0], [10 10 200 800]);

% plot spectrogram of calls
ax1 = axes3(1);
% get the frequency range
[~, ~, ~, ~, f, ~] = getAudioSpectrogramZZ_flexible_v1(calls.sound{1}, fs, 256, 256, 236, [250 7500], clim);
spec = call_mean(vi).spec_mean;
to_plot_freq = [1000 5000];
freq_i = find((f>=to_plot_freq(1)) & (f<=to_plot_freq(2)));
% also calculate the center of the spectrogramming sliding window
win_spec = 256 / fs * 1000;  % unit is ms
t_spec = (1:size(spec,2)) - 1 + win_spec - pad_sound*1000;
% t_spec = call_mean(vi).t_mean_spec;
% t_spec = t_spec*1000 - pad_sound*1000;
t_idx = find((t_spec>=ref_t(1)) & (t_spec<=ref_t(end)));

spec_plot = spec(freq_i, t_idx);
imagesc(ax1, linspace(0,1,size(spec_plot,2)), f(freq_i), spec_plot, [0.15 max(spec(:))]);
colormap(ax1, custom_map); set(ax1, 'YDir', 'normal');
yticks(ax1, [2000 4000]); 
yticklabels(ax1, {'2k', '4k'});
title(ax1, sprintf('v%d', vi), 'FontSize', 14);
xticks(ax1, []);

% then plot where the strands are
ax2 = axes3(2); cla(ax2); hold(ax2, 'on'); 
for row_i=1:size(xall,1)
  x_strand = find(xall(row_i,:)==1);
  plot(ax2, [x_strand(1) x_strand(end)], [row_i row_i], 'LineStyle', '-', 'LineWidth', 0.5, 'Color', [0.2 0.2 0.2 0.5]);
end
xlim(ax2, [1 size(xall,2)]);
ylim(ax2, [1 size(xall,1)]);
ylabel(ax2, 'Strand ID', 'FontSize', 12); 
% xlabel(ax2, 'Rel. time (ms)', 'FontSize', 12);
xticks(ax2, []);

% then calculate and plot probability
ax3 = axes3(3); cla(ax3); hold(ax3, 'on');
plot(ax3, x_prob); 
xlim(ax3, [1 length(x_prob)]);
ylim(ax3, [0 max(x_prob)*1.1]);
ylabel(ax3, 'Probability', 'FontSize', 12);     
xlabel(ax3, 'Rel. time (ms)', 'FontSize', 12);
% yline(ax3, 0.25, 'LineStyle', '--', 'Color', '#737373', 'LineWidth', 1); 
% xticks(ax3, []);

% find peaks 
param.gapSize = 10;  
param.thresholdFlatness = 0.03;  % first threshold to identify peak, 0.02 for v1
param.extendThreshold = 0.015; % after thresholding extend to lower threshold to certain percentage of the peak height
param.minInterval = 0.01;  % min distance between two adjacent peaks, smaller will be merged
param.minDuration = 0.025;  % minimal duration of the identified primitives
% param.maxDuration = 0.2; 
param.maxDuration = 0.1; 
param.minProminence = 0;  % use in the internal findpeaks function
[onsets, offsets] = ZZfunc_GetProbOnsetOffset_v2(x_prob, sec_per_frame, param);

% indicate where the peaks are
if ~isempty(onsets)
  for ii=1:length(onsets)
    if onsets(ii)==0; onsets(ii)=1; end
      xline(ax3, onsets(ii), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1);
      xline(ax3, offsets(ii), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
      % also plot in the spectrogram
      xline(ax1, onsets(ii)/length(x_prob), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1);
      xline(ax1, offsets(ii)/length(x_prob), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
  end
end

% save figure 
fn_fig = fullfile(fd_save_ref, sprintf('call_compound_density.v%d.pdf', vi));
print(fig3, fn_fig, '-dpdf', '-painters');



%% 7. Plot the location of the primitive in the compound syllable
% get the xy location of all similarity strands 
strands = [];
count = 0;
for si=1:size(pass_list, 2)
    p = pass_list(si);
    pass_i = p.pass_i;
    for strand_i=1:size(pass_i,2)
      [yy, xx] = ind2sub(p.size_d, pass_i{strand_i});
      x_range = [min(xx) max(xx)];
      y_range = [min(yy) max(yy)];
      count = count + 1;
      strands(count).si = si;
      strands(count).tar_i = p.tar_i;
      strands(count).size_d = p.size_d;
      strands(count).strand_i = strand_i;
      strands(count).x_range = x_range;
      strands(count).y_range = y_range;
      strands(count).x_center = mean(x_range);
      strands(count).y_center = mean(y_range);
    end   
end

% determine if the strands belong to which primitive
thre_overlap = 20;
for strand_i=1:size(strands,2)
  % loop through identified primitives
  prim_belong = nan;
  for prim_i=1:length(onsets)
    % calculate overlap
    a = [onsets(prim_i) offsets(prim_i)];
    b = strands(strand_i).x_range;
    x_center = strands(strand_i).x_center;
    [overlapRatio, intersectionLength] = calculateOverlapRatio(a, b);
    if (intersectionLength>=thre_overlap) && (x_center>=a(1)) && (x_center<=a(2))
%     if (x_center>=a(1)) && (x_center<=a(2))
      prim_belong = prim_i;
      break
    end
  end
  strands(strand_i).prim_belong = prim_belong;
  % add the information regarding strands onset/offset
  strands(strand_i).onsets = onsets; 
  strands(strand_i).offsets = offsets; 
end

% save the strands info for later use
fn_strand = fullfile(fd_save_ref, sprintf('call_v%d.strands.mat', vi));
save(fn_strand, 'strands');

col_prim = lines(length(onsets));
close all; 
fig = ZZfunc_newFigurePDFsize_v1([10 10 250 600]);
ax=gca; hold(ax, 'on');
for prim_i=1:length(onsets)
  % prim_i = 2;
  strands_plot = strands([strands.prim_belong]==prim_i);
  scatter(ax, [strands_plot.x_center], [strands_plot.y_center], 20, 'Marker', 'o', 'MarkerFaceColor', col_prim(prim_i,:), 'MarkerEdgeColor', 'none');
end
xlim(ax, [1 p.size_d(2)]);
ylim(ax, [1 800]);
xlabel('Latency in ref. (ms)', 'FontSize', 12);
ylabel('Latency in comp. (ms)', 'FontSize', 12);

fn_fig = fullfile(fd_save_ref, sprintf('call_v%d.target_loc.pdf', vi));
print(fig, fn_fig, '-dpdf', '-painters');




%% 8.1 Plot example cross-distance matrix
% choose compound syllables that have high similarity bands
m = cellfun(@(x) size(x,2), {pass_list.pass_i});
idx_plot = find(m~=0);

% chop the spectrogram of mean call
[~, ~, ~, ~, f, ~] = getAudioSpectrogramZZ_flexible_v1(calls.sound{1}, fs, 256, 256, 236, [250 7500], clim);
spec = call_mean(vi).spec_mean;
to_plot_freq = [1000 5000];
freq_i = find((f>=to_plot_freq(1)) & (f<=to_plot_freq(2)));
% also calculate the center of the spectrogramming sliding window
win_spec = 256 / fs * 1000;  % unit is ms
t_spec = (1:size(spec,2)) - 1 + win_spec - pad_sound*1000;
t_idx = find((t_spec>=ref_t(1)) & (t_spec<=ref_t(end)));
spec_plot = spec(freq_i, t_idx);

% loop through target syllable
for ii=1:50
  tar_i = idx_plot(ii);
  d2 = comp.dvae{tar_i};
  % grab sound
  [signal, fs] = audioread(comp.fn_wav{tar_i});
  i_start = max([1 comp.istart(tar_i)-pad_sound_pt]);
  i_end = min([size(signal,1) comp.iend(tar_i)+pad_sound_pt]);
  sound2 = signal(i_start:i_end);

  close all;
  [fig, axes_all, distm] = ZZfunc_plotCrossCorreTwoSpec_v4(d1c, d2, spec_plot, sound2, [12 23], custom_map, gray, [1000 5000], [0 0.35]);
  fn_fig = fullfile(fd_save_this, sprintf('call_v%d.tar%d.pdf', vi, tar_i));
  print(fig, fn_fig, '-dpdf', '-painters');

  % identify the similarity strands
  [fig2, ax_all2, pass_i, ~] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur);
  fn_fig = fullfile(fd_save_this, sprintf('call_v%d.tar%d.strand.pdf', vi, tar_i));
  print(fig2, fn_fig, '-dpdf', '-painters');
end
  


%% 8.2 Plot example cross-distance matrix: compound syllables on top
% choose compound syllables that have high similarity bands
fd_save_this2 = fullfile(fd_save_ref, sprintf('v%d_rev', vi));
if ~exist(fd_save_this2, 'dir'); mkdir(fd_save_this2); end
m = cellfun(@(x) size(x,2), {pass_list.pass_i});
idx_plot = find(m~=0);

% chop the spectrogram of mean call
[~, ~, ~, ~, f, ~] = getAudioSpectrogramZZ_flexible_v1(calls.sound{1}, fs, 256, 256, 236, [250 7500], clim);
spec = call_mean(vi).spec_mean;
to_plot_freq = [1000 5000];
freq_i = find((f>=to_plot_freq(1)) & (f<=to_plot_freq(2)));
% also calculate the center of the spectrogramming sliding window
win_spec = 256 / fs * 1000;  % unit is ms
t_spec = (1:size(spec,2)) - 1 + win_spec - pad_sound*1000;
t_idx = find((t_spec>=ref_t(1)) & (t_spec<=ref_t(end)));
spec_plot = spec(freq_i, t_idx);

% loop through target syllable
pad_comp = 0.016; pad_comp_pt = floor(pad_comp*fs);
for ii=201:400
  tar_i = idx_plot(ii);
  d2 = comp.dvae{tar_i};
  % grab sound
  [signal, fs] = audioread(comp.fn_wav{tar_i});
  i_start = max([1 comp.istart(tar_i)-pad_sound_pt]);
  i_end = min([size(signal,1) comp.iend(tar_i)+pad_sound_pt]);
  sound2 = signal(i_start:i_end);

  close all;
  [fig, axes_all, distm] = ZZfunc_plotCrossCorreTwoSpec_v4_rev(d1c, d2, spec_plot, sound2, [12 23], custom_map, gray, [1000 5000], [0 0.35]);
  fn_fig = fullfile(fd_save_this2, sprintf('call_v%d.tar%d.pdf', vi, tar_i));
  print(fig, fn_fig, '-dpdf', '-painters');

  % identify the similarity strands
  [fig2, ax_all2, pass_i, ~] = ZZfunc_identifySimilarityStrand_v1(flipud(distm), med_filter, thre, min_dur);
  fn_fig = fullfile(fd_save_this2, sprintf('call_v%d.tar%d.strand.pdf', vi, tar_i));
  print(fig2, fn_fig, '-dpdf', '-painters');
end
  