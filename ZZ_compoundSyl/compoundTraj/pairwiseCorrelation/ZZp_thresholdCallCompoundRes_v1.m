% re-set the threshold when identifying primitives for call-compound syllable comparison
% may need to use different thresholds for different calls
% run this after running ZZp_callCompoundCorre_v3

clear; close all;

addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_ephys_utils'));
addpath(genpath('/home/zz367/ProjectsU/EphysMONAO/Jupyter/MatlabCodes/ZZ_compoundSyl'));


%% Folder setting
fd_z4 = '/mnt/z4';
fd_base = fullfile(fd_z4, 'zz367', 'EphysMONAO', 'Analyzed');
birdIDs = {'pair5RigCCU29', 'pair4RigACU68', 'pair4RigBCU53', 'pair2RigBCU25'};
pairIDs = {'pair5CU29CU55', 'pair4CU68CU53', 'pair4CU68CU53', 'pair2CU20CU25'};
pretty_ids = {'M1', 'M2', 'M3', 'M4'};
% what VAE run to use
vae_run = 'traj_chop_32_1_32';
% what's the pad for sound
fs = 20000;
pad_sound = 0.032; pad_sound_pt = floor(pad_sound*fs);
win_size = 32;
% a custom colormap
n_colors = 256;             % Total number of colors in the colormap
jet_map = jet(n_colors - 1);  % Get standard jet, but with one fewer row
custom_map = [0 0 0; jet_map];  % Prepend black to jet



bi = 4;
birdID = birdIDs{bi};
pairID = pairIDs{bi};


%% Bird-specific folders
% where is results for call-compound analysis stored
fd_base_res = fullfile(fd_base, 'Figures', pairID, 'CompoundSyl', birdID, 'pairwiseCorre');
fd_res = fullfile(fd_base_res, 'compound_call3');
% grab results for all call subtypes
fns_res = dir(fullfile(fd_res, '*.pass_list.mat'));
pass_all = struct();
syls = {};
for fi=1:size(fns_res,1)
  temp = strsplit(fns_res(fi).name, '.');
  v = strrep(temp{1}, 'call_', '');
  a = load(fullfile(fns_res(fi).folder, fns_res(fi).name));
  pass_all.(v) = a.pass_list;
  syls = [syls; v];
end

% also load the mean spectrogram for calls for plotting
fn_call_mean = fullfile(fd_base_res, sprintf('%s.call_mean3.mat', birdID));
load(fn_call_mean);
% also load the call dataset
% fn_call = fullfile(fd_base_res, sprintf('%s.call3.mat', birdID));
% load(fn_call);

% where to save the combined results
fd_save = fullfile(fd_res, 'combined');
if ~exist(fd_save, 'dir'); mkdir(fd_save); end



%% Plot the density of high-similar strands
% get the frequency range from a pseudo-sound
[~, ~, ~, ~, f, ~] = getAudioSpectrogramZZ_flexible_v1(randn(10000,1), fs, 256, 256, 236, [250 7500], [10 21]);
half_win = floor(win_size/2);
rm_size = 12;  % how much to remove from the silent windows before syl onset and after syl offset
sec_per_frame = 0.001;

% save info regarding the prob 
xprob_all = [];
close all; 
[fig, axs] = generatePanelGrid_v2(3, 6, [0.17;0.3;0.15], [0.02;0.02], [0.05;0.05], [0.05;0.05], 0.05, [0;0;0], [10 10 1200 800]);
for vi=1:size(syls,1)
  v = syls{vi};
  pass_list = pass_all.(v);
  
  %% calculate density of high-similar strands on the top of the calls
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
  
  % calculate density
  x_prob = sum(xall, 1) / size(pass_list,2);
  ref_t = (1:length(x_prob)) - (half_win - rm_size); 

  
  %% plot results: spectrogram; strands; density
  ax1 = axs(1, vi);
  spec = call_mean(vi).spec_mean;
  to_plot_freq = [1000 5000];
  freq_i = find((f>=to_plot_freq(1)) & (f<=to_plot_freq(2)));
  % also calculate the center of the spectrogramming sliding window
  win_spec = 256 / fs * 1000;  % unit is ms
  t_spec = (1:size(spec,2)) - 1 + win_spec - pad_sound*1000;
  t_idx = find((t_spec>=ref_t(1)) & (t_spec<=ref_t(end)));
  spec_plot = spec(freq_i, t_idx);
  imagesc(ax1, linspace(0,1,size(spec_plot,2)), f(freq_i), spec_plot, [0.15 max(spec(:))]);
  colormap(ax1, custom_map); set(ax1, 'YDir', 'normal');
  yticks(ax1, [2000 4000]); 
  yticklabels(ax1, {'2k', '4k'});
  title(ax1, sprintf('v%d', vi), 'FontSize', 14);
  xticks(ax1, []);
  
  % then high-similarity strands as bars 
  ax2 = axs(2, vi); cla(ax2); hold(ax2, 'on'); 
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
  ax3 = axs(3, vi); cla(ax3); hold(ax3, 'on');
  plot(ax3, x_prob); 
  xlim(ax3, [1 length(x_prob)]);
  ylim(ax3, [0 max(x_prob)*1.1]);
  ylabel(ax3, 'Probability', 'FontSize', 12);     
  xlabel(ax3, 'Rel. time (ms)', 'FontSize', 12);
  
  xprob_all(vi).v = v;
  xprob_all(vi).xall = xall;
  xprob_all(vi).x_prob = x_prob;
  xprob_all(vi).ref_t = ref_t; 
  xprob_all(vi).spec_plot = spec_plot;
end
% set the same y range for all density plot
dens_max = cellfun(@(x) max(x), {xprob_all.x_prob});
% for vi=1:size(syls,1)
%   ylim(axs(3,vi), [0 max(dens_max)*1.05]);
% end

% finally using simple threshoding to identify peaks
param.gapSize = 10;  
param.minInterval = 0.01;  % min distance between two adjacent peaks, smaller will be merged
param.minDuration = 0.025;  % minimal duration of the identified primitives
% param.maxDuration = 0.2; 
param.maxDuration = 0.1; 
param.minProminence = 0;  % use in the internal findpeaks function
param.smooth_w = 10;
% try different thresholds for different calls
% thre_list = [0.03 0.03  0.015 0.015 0.015 0.015];  % bird M1
% thre_list = [0.015 0.015 0.015 0.015];  % bird M2
thre_list = [0.03 0.03 0.03 0.03]; % bird M4
% ext_list = [0.015 0.015 0.01 0.01 0.01 0.01];
ext_list = thre_list * 0.5;
for vi=1:size(syls,1)
  param.thresholdFlatness = thre_list(vi);
  param.extendThreshold = ext_list(vi);
  x_prob = xprob_all(vi).x_prob;
  [onsets, offsets] = ZZfunc_GetProbOnsetOffset_v3(x_prob, sec_per_frame, param);
%   out = detect_two_peaks(x_prob, 'SmoothW', 7, 'MinProm', 0, 'MinDist', 100);
  % indicate where the peaks are
  if ~isempty(onsets)
    for ii=1:length(onsets)
      if onsets(ii)==0; onsets(ii)=1; end
        xline(axs(3,vi), onsets(ii), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1);
        xline(axs(3,vi), offsets(ii), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
        % also plot in the spectrogram
        xline(axs(1,vi), onsets(ii)/length(x_prob), 'Color', 'green', 'LineStyle', '--', 'LineWidth', 1);
        xline(axs(1,vi), offsets(ii)/length(x_prob), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
    end
  end
  
  xprob_all(vi).onsets = onsets;
  xprob_all(vi).offsets = offsets;
end
 
fn_fig = fullfile(fd_save, sprintf('%s.strand_density.pdf', birdID));
print(fig, fn_fig, '-dpdf', '-painters');

% save info regarding primitives onset/offset as well
% manually assess what primitives are good for this bird
% is_good = {[1 1]; [0]; [1 1]; [0]; [0]; [1]}; % bird M1
% is_good = {[1]; [0]; [1]; [0]};  % bird M2
is_good = {[0]; [0]; [1 1]; [0 0]};  % bird M4
for vi=1:size(xprob_all,2); xprob_all(vi).is_good=is_good{vi}; end
fn_xprob = fullfile(fd_save, sprintf('%s.xprob_all.mat', birdID));
save(fn_xprob, 'xprob_all');
  















