% compare compound syllables to calls, generate the cross-distance matrix, check if they share primitives
% Zhilei, 09/22/2025


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
pad = 0.016; pad_pt=floor(pad*fs);
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
  i_start = max([1 calls.istart(ii)-pad_pt]);
  i_end = min([size(signal,1) calls.iend(ii)+pad_pt]);
  calls.sound{ii} = signal(i_start:i_end);
end
% save for later use
fn_call = fullfile(fd_save, sprintf('%s.call.mat', birdID));
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
end
% save for later use
fn_call_mean = fullfile(fd_save, sprintf('%s.call_mean.mat', birdID));
save(fn_call_mean, 'call_mean');



%% 3. Concatenate the mean calls into a long syllable
cat_call_vae = cat(1, call_mean(:).vae_mean);
cat_call_spec = cat(2, call_mean(:).spec_mean);
% plot to check
figure; ax=gca; hold(ax, 'on');
imagesc(ax, cat_call_spec, [0.15 max(cat_call_spec(:))]); colormap(custom_map);
set(ax, 'YDir', 'Normal');



%% 5. Loop compound syllable, compare to the mean call VAE
fd_save_ref = fullfile(fd_save, 'compound_call');
if exist(fd_save_ref, 'dir'); rmdir(fd_save_ref, 's'); end
mkdir(fd_save_ref);
% for ref_i=1:size(comp,1)
pad = 0.016; pad_pt=floor(pad*fs);

% loop through call subtypes 
vi = 1; 
d1 = call_mean(vi).vae_mean;
spec1 = call_mean(vi).spec_mean;
% set to defire frequency range
flim = [1000 5000];
[~, ~, ~, ~, f, ~] = getAudioSpectrogramZZ_flexible_v1(calls.sound{1}, fs, 256, 256, 236, [250 7500], clim);
freq_i = find((f>=flim(1)) & (f<=flim(2)));
spec1 = spec1(freq_i,:);

is_plot = true;
% is_plot = false;
pass_list = [];
% for tar_i=1:size(comp,1)
for tar_i=1:100
%   disp(tar_i);
  % grab VAE latents
  d2 = comp.dvae{tar_i};
  % parameters for identifying bands
  med_filter = [11 11];
  thre = 0.35;
  min_dur = 25;
  if is_plot
    % grab sound
    [signal, fs] = audioread(comp.fn_wav{tar_i});
    i_start = max([1 comp.istart(tar_i)-pad_pt]);
    i_end = min([size(signal,1) comp.iend(tar_i)+pad_pt]);
    sound2 = signal(i_start:i_end);
    close all;
    [fig, axes_all, distm] = ZZfunc_plotCrossCorreTwoSpec_v4(d1, d2, spec1, sound2, [12 23], custom_map, gray, [1000 5000]);
    fn_fig = fullfile(fd_save_ref, sprintf('call_v%d.tar%d.pdf', vi, tar_i));
    print(fig, fn_fig, '-dpdf', '-painters');
    
    % identify the similarity strands
    [fig2, ax_all2, pass_i, ~] = ZZfunc_identifySimilarityStrand_v1(distm, med_filter, thre, min_dur);
    fn_fig = fullfile(fd_save_ref, sprintf('call_v%d.tar%d.strand.pdf', vi, tar_i));
    print(fig2, fn_fig, '-dpdf', '-painters');
  else
    % no plotting
    distm = pdist2(d1, d2, 'cosine');
    distm = distm';
    % no plotting
    pass_i = ZZfunc_identifySimilarityStrand_v1_noPlot(distm, med_filter, thre, min_dur);
  end
  pass_list(tar_i).tar_i = tar_i;
  pass_list(tar_i).size_d = size(distm);
  pass_list(tar_i).pass_i = pass_i;
end
% save result
fd_pass_list = fullfile(fd_save_ref, sprintf('call_v%d.pass_list.mat', vi));
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
    % add the half-width of sliding windows back, so considering the center of sliding window
%     x_range = x_range + win_size/2;
    count = count + 1;
%     iend = min([n(2) x_range(2)]);  % ignore those that pass the boundaries
    xall(count, x_range(1):x_range(2)) = 1;
  end   
end
  
% calculate density
% x_prob = sum(x_cat, 1) / size(x_cat, 1);  % divided by the total number of strands
x_prob = sum(xall, 1) / size(pass_list,2);


%% plot results
close all;
[fig3, axes3] = generatePanelGrid_v2(3, 1, [0.15;0.3;0.15], [0.02;0.02], [0.05;0.05], [0.25;0.05], 0.05, [0;0;0], [10 10 200 800]);

% plot spectrogram of calls
ax1 = axes3(1);
% get the frequency range
[~, ~, ~, ~, f, ~] = getAudioSpectrogramZZ_flexible_v1(calls.sound{1}, fs, 256, 256, 236, [250 7500], clim);
spec = call_mean(vi).spec_mean;
to_plot_freq = [1000 5000];
freq_i = find((f>=to_plot_freq(1)) & (f<=to_plot_freq(2)));
imagesc(ax1, 1:size(spec,2), f(freq_i), spec(freq_i,:), [0.15 max(spec(:))]);
colormap(ax1, custom_map); set(ax1, 'YDir', 'normal');
yticks(ax1, [2000 4000]); 
yticklabels(ax1, {'2k', '4k'});
title(ax1, sprintf('v%d', vi), 'FontSize', 14);
xticks(ax1, []);

% then plot where the strands are
ax2 = axes3(2); cla(ax2); hold(ax2, 'on'); 
for row_i=1:5000
  x_strand = find(xall(row_i,:)==1);
  plot(ax2, [x_strand(1) x_strand(end)], [row_i row_i], 'LineStyle', '-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5 0.1]);
end
xlim(ax2, [1 size(xall,2)]);
% ylim(ax2, [0 size(xall,1)]);
ylabel(ax2, 'Strand ID', 'FontSize', 12); 
% xlabel(ax2, 'Rel. time (ms)', 'FontSize', 12);
xticks(ax2, []);

% then calculate and plot probability
ax3 = axes3(3); cla(ax3); hold(ax3, 'on');
rel_t = (1:length(x_prob)) - win_size/2-1;
plot(ax3, rel_t, x_prob); 
xlim(ax3, [rel_t(1) rel_t(end)]);
ylabel(ax3, 'Probability', 'FontSize', 12);     
xlabel(ax3, 'Rel. time (ms)', 'FontSize', 12);
yline(ax3, 0.25, 'LineStyle', '--', 'Color', '#737373', 'LineWidth', 1); 
% xticks(ax3, []);

% find peaks 
param.gapSize = 10;  
param.thresholdFlatness = 0.35;  % first threshold to identify peak
param.extendThreshold = 0.25; % after thresholding extend to lower threshold to certain percentage of the peak height
param.minInterval = 0.01;  % min distance between two adjacent peaks, smaller will be merged
param.minDuration = 0.025;  % minimal duration of the identified primitives
% param.maxDuration = 0.2; 
param.maxDuration = 0.1; 
param.minProminence = 0;  % use in the internal findpeaks function
[onsets, offsets] = ZZfunc_GetProbOnsetOffset_v2(x_prob, sec_per_frame, param);

% indicate where the peaks are
% if ~isempty(onsets)
%   for ii=1:length(onsets)
%     if onsets(ii)==0; onsets(ii)=1; end
%       xline(ax3, rel_t(onsets(ii)), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 2);
%       xline(ax3, rel_t(offsets(ii)), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
%       % also plot in the spectrogram
%       xline(ax1, onsets(ii)/length(prob_amp), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 2);
%       xline(ax1, offsets(ii)/length(prob_amp), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
%   end
% end

% save figure 
fn_fig = fullfile(fd_save_ref, sprintf('call_compound_density.v%d.pdf', vi));
print(fig3, fn_fig, '-dpdf', '-painters');




%% 7. Sort the matrix
dataMatrix = xall(1:1000,:);
Y = pdist(dataMatrix, 'euclidean');
Z = linkage(Y, 'average');
leafOrder = optimalleaforder(Z, Y);
sortedMatrix = dataMatrix(leafOrder, :);


% --- 6. (Optional) Visualize the Result ---
% A great way to check your work is to visualize the
% original and sorted matrices as heatmaps.

figure;

% Plot the original, unsorted data
subplot(2, 1, 1);
imagesc(dataMatrix);
title('Original Matrix (Unsorted)');
ylabel('Samples');
xlabel('Features');

% Plot the new, sorted data
subplot(2, 1, 2);
imagesc(sortedMatrix);
title('Sorted Matrix (Similar Rows are Near Each Other)');
ylabel('Sorted Samples');
xlabel('Features');
colorbar;

















